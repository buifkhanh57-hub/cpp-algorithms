// statsengine — outliers.cpp
//
// z-score / IQR-fence / modified-z (MAD) detection and rendering.
// See outliers.hpp for the statistical contracts.
#include "statsengine/outliers.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"
#include "statsengine/table.hpp"

namespace se {

std::string outlierMethodName(OutlierMethod method) {
  switch (method) {
    case OutlierMethod::ZScore:
      return "z-score";
    case OutlierMethod::Iqr:
      return "iqr";
    case OutlierMethod::ModifiedZ:
      break;
  }
  return "modified-z (MAD)";
}

bool outlierMethodFromName(const std::string& name, OutlierMethod& out) {
  const std::string key = toLower(trim(name));
  if (key == "z" || key == "zscore" || key == "z-score") {
    out = OutlierMethod::ZScore;
  } else if (key == "iqr" || key == "tukey") {
    out = OutlierMethod::Iqr;
  } else if (key == "mad" || key == "modified-z" || key == "modifiedz") {
    out = OutlierMethod::ModifiedZ;
  } else {
    return false;
  }
  return true;
}

double medianSorted(const std::vector<double>& sorted) {
  return quantileSorted(sorted, 0.5);
}

double medianAbsoluteDeviation(const std::vector<double>& values) {
  std::vector<double> sorted(values);
  std::sort(sorted.begin(), sorted.end());
  const double median = medianSorted(sorted);
  std::vector<double> deviations;
  deviations.reserve(sorted.size());
  for (double v : sorted) deviations.push_back(std::fabs(v - median));
  std::sort(deviations.begin(), deviations.end());
  return medianSorted(deviations);
}

OutlierReport detectOutliers(const std::vector<double>& values, const OutlierOptions& options) {
  std::vector<std::pair<std::size_t, double>> points;  // (original index, value)
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (std::isfinite(values[i])) points.emplace_back(i, values[i]);
  }
  if (points.size() < 3) {
    throw std::runtime_error("outlier detection needs at least 3 numeric values");
  }

  OutlierReport report;
  report.method = options.method;
  report.n = static_cast<long long>(points.size());
  report.missing = static_cast<long long>(values.size()) - report.n;

  std::vector<double> sorted;
  sorted.reserve(points.size());
  for (const auto& point : points) sorted.push_back(point.second);
  std::sort(sorted.begin(), sorted.end());

  switch (options.method) {
    case OutlierMethod::ZScore: {
      RunningMoments moments;
      for (double v : sorted) moments.push(v);
      const double mean = moments.mean();
      const double sd = moments.sampleStddev();
      report.center = mean;
      report.spread = sd;
      report.usedThreshold = options.threshold > 0.0 ? options.threshold : 3.0;
      if (sd > 0.0) {
        for (const auto& point : points) {
          const double score = (point.second - mean) / sd;
          if (std::fabs(score) > report.usedThreshold) {
            report.flags.push_back({point.first, point.second, score});
          }
        }
      }
      break;
    }
    case OutlierMethod::Iqr: {
      const double q1 = quantileSorted(sorted, 0.25);
      const double q3 = quantileSorted(sorted, 0.75);
      const double iqr = q3 - q1;
      const double k = options.threshold > 0.0 ? options.threshold : 1.5;
      report.center = medianSorted(sorted);
      report.spread = iqr;
      report.usedThreshold = k;
      report.fenceLow = q1 - k * iqr;
      report.fenceHigh = q3 + k * iqr;
      if (iqr > 0.0) {
        for (const auto& point : points) {
          const double v = point.second;
          if (v < report.fenceLow || v > report.fenceHigh) {
            const double beyond =
                v > report.fenceHigh ? v - report.fenceHigh : report.fenceLow - v;
            report.flags.push_back({point.first, v, beyond / iqr});
          }
        }
      }
      break;
    }
    case OutlierMethod::ModifiedZ: {
      const double median = medianSorted(sorted);
      const double mad = medianAbsoluteDeviation(sorted);
      report.center = median;
      report.spread = mad;
      report.usedThreshold = options.threshold > 0.0 ? options.threshold : 3.5;
      for (const auto& point : points) {
        const double deviation = point.second - median;
        if (mad > 0.0) {
          const double score = 0.6745 * deviation / mad;
          if (std::fabs(score) > report.usedThreshold) {
            report.flags.push_back({point.first, point.second, score});
          }
        } else if (deviation != 0.0) {
          // MAD == 0 means more than half the values are identical; any
          // deviation is therefore an outlier with an unbounded score.
          report.flags.push_back(
              {point.first, point.second, deviation > 0.0 ? 1e300 : -1e300});
        }
      }
      break;
    }
  }

  std::sort(report.flags.begin(), report.flags.end(),
            [](const OutlierFlag& a, const OutlierFlag& b) {
              return std::fabs(a.score) > std::fabs(b.score);
            });
  return report;
}

std::string renderOutliers(const OutlierReport& report, const std::string& columnName,
                           bool unicode, bool color, std::size_t maxRows) {
  std::vector<std::pair<std::string, std::string>> items;
  items.emplace_back("column", columnName);
  items.emplace_back("method", outlierMethodName(report.method));
  items.emplace_back("threshold", fmtFixed(report.usedThreshold, 3));
  items.emplace_back("n", std::to_string(report.n) +
                              (report.missing > 0
                                   ? " (excluding " + std::to_string(report.missing) + " missing)"
                                   : ""));
  items.emplace_back("center", fmtNumber(report.center, 6));
  items.emplace_back("spread", fmtNumber(report.spread, 6));
  if (report.method == OutlierMethod::Iqr) {
    items.emplace_back("fences", "[" + fmtFixed(report.fenceLow, 4) + ", " +
                                     fmtFixed(report.fenceHigh, 4) + "]");
  }
  const double flaggedFraction =
      report.n > 0 ? static_cast<double>(report.flags.size()) / static_cast<double>(report.n)
                   : 0.0;
  items.emplace_back("flagged", std::to_string(report.flags.size()) + " (" +
                                    fmtPercent(flaggedFraction, 2) + ")");

  std::ostringstream out;
  out << keyValueBlock(items, color) << "\n";

  if (report.flags.empty()) {
    out << colors::green("no outliers detected") << "\n";
    return out.str();
  }

  Table table;
  table.setTitle("Flagged values (sorted by |score|)");
  table.addColumn("row", Align::Right);
  table.addColumn("value", Align::Right);
  table.addColumn("score", Align::Right);
  const std::size_t shown = std::min(report.flags.size(), maxRows);
  for (std::size_t i = 0; i < shown; ++i) {
    const OutlierFlag& flag = report.flags[i];
    const std::string score = flag.score >= 1e299 || flag.score <= -1e299
                                  ? "inf"
                                  : fmtFixed(flag.score, 4);
    std::string valueText = fmtNumber(flag.value, 6);
    if (color) valueText = colors::red(valueText);
    table.addRow({std::to_string(flag.index + 1), valueText, score});
  }
  if (report.flags.size() > shown) {
    table.setFooter(std::to_string(report.flags.size() - shown) + " more flagged values not shown");
  }
  return out.str() + table.render(unicode);
}

void outliersToJson(json::Writer& writer, const OutlierReport& report,
                    const std::string& columnName) {
  writer.beginObject();
  writer.kv("column", columnName);
  writer.kv("method", outlierMethodName(report.method));
  writer.kv("threshold", report.usedThreshold);
  writer.kv("n", report.n);
  writer.kv("missing", report.missing);
  writer.kv("center", report.center);
  writer.kv("spread", report.spread);
  if (report.method == OutlierMethod::Iqr) {
    writer.key("fences");
    writer.beginArray();
    writer.value(report.fenceLow);
    writer.value(report.fenceHigh);
    writer.endArray();
  }
  writer.key("flagged");
  writer.beginArray();
  for (const OutlierFlag& flag : report.flags) {
    writer.beginObject();
    writer.kv("index", static_cast<long long>(flag.index));
    writer.kv("value", flag.value);
    writer.kv("score", std::fabs(flag.score) >= 1e299 ? 1e299 * (flag.score > 0 ? 1 : -1)
                                                      : flag.score);
    writer.endObject();
  }
  writer.endArray();
  writer.endObject();
}

}  // namespace se
