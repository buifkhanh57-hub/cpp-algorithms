// statsengine — histogram.cpp
//
// Binning rules and ASCII rendering. See histogram.hpp.
#include "statsengine/histogram.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"

namespace se {

std::string binningRuleName(BinningRule rule) {
  switch (rule) {
    case BinningRule::Auto:
      return "auto";
    case BinningRule::Sturges:
      return "sturges";
    case BinningRule::FreedmanDiaconis:
      return "freedman-diaconis";
    case BinningRule::Scott:
      return "scott";
    case BinningRule::Fixed:
      break;
  }
  return "fixed";
}

bool binningRuleFromName(const std::string& name, BinningRule& out) {
  const std::string key = toLower(trim(name));
  if (key == "auto") {
    out = BinningRule::Auto;
  } else if (key == "sturges") {
    out = BinningRule::Sturges;
  } else if (key == "fd" || key == "freedman-diaconis" || key == "freedman") {
    out = BinningRule::FreedmanDiaconis;
  } else if (key == "scott") {
    out = BinningRule::Scott;
  } else if (key == "fixed") {
    out = BinningRule::Fixed;
  } else {
    return false;
  }
  return true;
}

namespace {

struct BinChoice {
  int bins = 1;
  BinningRule rule = BinningRule::Sturges;
};

// Population variance of already-sorted finite values (local helper).
double sortedPopulationVariance(const std::vector<double>& sorted) {
  RunningMoments moments;
  for (double v : sorted) moments.push(v);
  return moments.populationVariance();
}

BinChoice chooseBins(const std::vector<double>& sorted, const HistogramOptions& options) {
  BinChoice choice;
  const double n = static_cast<double>(sorted.size());
  const double min = sorted.front();
  const double max = sorted.back();

  auto sturges = [&]() {
    choice.bins = std::max(1, static_cast<int>(std::ceil(std::log2(n))) + 1);
    choice.rule = BinningRule::Sturges;
  };

  switch (options.rule) {
    case BinningRule::Fixed:
      choice.bins = options.bins;
      choice.rule = BinningRule::Fixed;
      break;
    case BinningRule::Sturges:
      sturges();
      break;
    case BinningRule::Scott: {
      const double sd = std::sqrt(sortedPopulationVariance(sorted));
      const double h = 3.49 * sd / std::cbrt(n);
      if (!(h > 0.0) || !std::isfinite(h) || max <= min) {
        sturges();
      } else {
        choice.bins = std::min(200, std::max(1, static_cast<int>(std::ceil((max - min) / h))));
        choice.rule = BinningRule::Scott;
      }
      break;
    }
    case BinningRule::FreedmanDiaconis:
    case BinningRule::Auto: {
      const double q1 = quantileSorted(sorted, 0.25);
      const double q3 = quantileSorted(sorted, 0.75);
      const double iqr = q3 - q1;
      const double h = 2.0 * iqr / std::cbrt(n);
      if (!(h > 0.0) || !std::isfinite(h) || max <= min) {
        sturges();  // IQR == 0 (or degenerate range) -> Sturges fallback
      } else {
        choice.bins = std::min(200, std::max(1, static_cast<int>(std::ceil((max - min) / h))));
        choice.rule = options.rule;
      }
      break;
    }
  }
  choice.bins = std::min(512, std::max(1, choice.bins));
  return choice;
}

std::string binLabel(const Bin& bin, bool last) {
  const std::string close = last ? "]" : ")";
  return "[" + fmtFixed(bin.lo, 3) + ", " + fmtFixed(bin.hi, 3) + close;
}

std::string barGlyphs(int length, bool unicode) {
  if (unicode) {
    std::string bar;
    for (int i = 0; i < length; ++i) bar += "\xe2\x96\x88";  // U+2588 FULL BLOCK
    return bar;
  }
  return std::string(static_cast<std::size_t>(std::max(0, length)), '#');
}

}  // namespace

HistogramResult computeHistogram(const std::vector<double>& values,
                                 const HistogramOptions& options) {
  std::vector<double> finite;
  finite.reserve(values.size());
  for (double v : values) {
    if (std::isfinite(v)) finite.push_back(v);
  }
  if (finite.empty()) {
    throw std::runtime_error("histogram requires at least one numeric value");
  }
  if (options.rule == BinningRule::Fixed && options.bins <= 0) {
    throw std::invalid_argument("--bins must be a positive integer when rule=fixed");
  }
  if (options.binWidth < 0.0) {
    throw std::invalid_argument("--bin-width must be non-negative");
  }

  std::sort(finite.begin(), finite.end());
  const double min = finite.front();
  const double max = finite.back();

  HistogramResult result;
  result.total = static_cast<long long>(finite.size());
  result.barWidth = options.barWidth;
  result.min = min;
  result.max = max;

  RunningMoments moments;
  for (double v : finite) moments.push(v);
  result.mean = moments.mean();
  result.stddev = moments.sampleStddev();
  result.median = quantileSorted(finite, 0.5);

  int bins = 1;
  if (max > min) {
    if (options.binWidth > 0.0) {
      const double span = max - min;
      bins = std::min(512, std::max(1, static_cast<int>(std::ceil(span / options.binWidth))));
      result.usedRule = BinningRule::Fixed;
      result.usedRuleName = "fixed-width";
    } else {
      const BinChoice choice = chooseBins(finite, options);
      bins = choice.bins;
      result.usedRule = choice.rule;
      result.usedRuleName = binningRuleName(choice.rule);
    }
  } else {
    result.usedRule = BinningRule::Fixed;
    result.usedRuleName = "single-value";
  }

  const double width = (max > min) ? (max - min) / static_cast<double>(bins) : 1.0;
  result.binWidth = width;
  result.bins.assign(static_cast<std::size_t>(bins), Bin{});
  for (int b = 0; b < bins; ++b) {
    result.bins[static_cast<std::size_t>(b)].lo = min + width * static_cast<double>(b);
    result.bins[static_cast<std::size_t>(b)].hi = min + width * static_cast<double>(b + 1);
  }

  for (double v : finite) {
    std::size_t idx = static_cast<std::size_t>((v - min) / width);
    if (idx >= result.bins.size()) idx = result.bins.size() - 1;  // v == max
    result.bins[idx].count += 1;
  }
  return result;
}

std::string renderHistogram(const HistogramResult& result, const std::string& columnName,
                            bool unicode, bool color) {
  (void)color;  // emphasis is applied inline via colors::* helpers
  std::ostringstream out;
  out << colors::bold("Histogram of '" + columnName + "'") << "\n";
  out << colors::dim("n = " + std::to_string(result.total) +
                     " | bins = " + std::to_string(result.bins.size()) +
                     " | rule = " + result.usedRuleName +
                     " | width = " + fmtNumber(result.binWidth, 4))
      << "\n";
  out << colors::dim("min = " + fmtNumber(result.min, 4) + " | max = " + fmtNumber(result.max, 4) +
                     " | mean = " + fmtNumber(result.mean, 4) +
                     " | sd = " + fmtNumber(result.stddev, 4) +
                     " | median = " + fmtNumber(result.median, 4))
      << "\n\n";

  long long maxCount = 0;
  for (const Bin& bin : result.bins) maxCount = std::max(maxCount, bin.count);
  if (maxCount == 0) maxCount = 1;

  std::size_t labelWidth = 0;
  for (std::size_t i = 0; i < result.bins.size(); ++i) {
    labelWidth =
        std::max(labelWidth, binLabel(result.bins[i], i + 1 == result.bins.size()).size());
  }

  long long cumulative = 0;
  auto padBefore = [](std::size_t width, const std::string& text) -> std::string {
    // Two extra columns keep at least one separating space between fields.
    width += 2;
    return text.size() >= width ? std::string() : std::string(width - text.size(), ' ');
  };
  for (std::size_t i = 0; i < result.bins.size(); ++i) {
    const Bin& bin = result.bins[i];
    cumulative += bin.count;
    const std::string label = binLabel(bin, i + 1 == result.bins.size());
    const double fraction =
        static_cast<double>(bin.count) / static_cast<double>(maxCount);
    int barLen = static_cast<int>(fraction * static_cast<double>(result.barWidth) + 0.5);
    if (bin.count > 0 && barLen == 0) barLen = 1;

    const std::string countText = std::to_string(bin.count);
    const std::string pct =
        fmtPercent(static_cast<double>(bin.count) / static_cast<double>(result.total), 1);
    const std::string cum =
        fmtPercent(static_cast<double>(cumulative) / static_cast<double>(result.total), 1);

    out << "  " << colors::gray(label)
        << std::string(labelWidth - colors::visibleWidth(label) + 2, ' ')
        << padBefore(4, countText) << countText
        << padBefore(4, pct) << pct
        << padBefore(4, cum) << colors::dim(cum) << "  "
        << colors::cyan(barGlyphs(barLen, unicode)) << "\n";
  }
  return out.str();
}

void histogramToJson(json::Writer& writer, const HistogramResult& result,
                     const std::string& columnName) {
  writer.beginObject();
  writer.kv("column", columnName);
  writer.kv("n", result.total);
  writer.kvCount("bins", result.bins.size());
  writer.kv("rule", result.usedRuleName);
  writer.kv("bin_width", result.binWidth);
  writer.kv("min", result.min);
  writer.kv("max", result.max);
  writer.kv("mean", result.mean);
  writer.kv("stddev", result.stddev);
  writer.kv("median", result.median);
  writer.key("bins_list");
  writer.beginArray();
  for (const Bin& bin : result.bins) {
    writer.beginObject();
    writer.kv("lo", bin.lo);
    writer.kv("hi", bin.hi);
    writer.kv("count", bin.count);
    writer.endObject();
  }
  writer.endArray();
  writer.endObject();
}

}  // namespace se
