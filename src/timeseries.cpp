// statsengine — timeseries.cpp
//
// SMA / EMA / differencing and the series renderer. See timeseries.hpp.
#include "statsengine/timeseries.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/table.hpp"
#include "statsengine/util.hpp"

namespace se {

std::string maMethodName(MaMethod method) {
  switch (method) {
    case MaMethod::Sma:
      return "sma";
    case MaMethod::Ema:
      return "ema";
    case MaMethod::Diff:
      break;
  }
  return "diff";
}

bool maMethodFromName(const std::string& name, MaMethod& out) {
  const std::string key = toLower(trim(name));
  if (key == "sma" || key == "moving-average" || key == "ma") {
    out = MaMethod::Sma;
  } else if (key == "ema" || key == "ewma") {
    out = MaMethod::Ema;
  } else if (key == "diff" || key == "difference") {
    out = MaMethod::Diff;
  } else {
    return false;
  }
  return true;
}

MaResult movingAverage(const std::vector<double>& x, const MaOptions& options) {
  MaResult result;
  result.method = options.method;

  switch (options.method) {
    case MaMethod::Sma: {
      if (options.window < 1) {
        throw std::invalid_argument("--window must be a positive integer");
      }
      result.param = options.window;
      result.label = "sma(window=" + std::to_string(options.window) + ")";
      result.values.assign(x.size(), kNaN);
      for (std::size_t t = 0; t < x.size(); ++t) {
        const std::size_t begin = t + 1 >= static_cast<std::size_t>(options.window)
                                      ? t + 1 - static_cast<std::size_t>(options.window)
                                      : 0;
        double sum = 0.0;
        std::size_t valid = 0;
        for (std::size_t i = begin; i <= t; ++i) {
          if (std::isfinite(x[i])) {
            sum += x[i];
            ++valid;
          }
        }
        if (valid > 0) result.values[t] = sum / static_cast<double>(valid);
      }
      break;
    }
    case MaMethod::Ema: {
      double alpha = options.alpha;
      int span = options.span;
      if (alpha <= 0.0) {
        if (span > 0) {
          alpha = 2.0 / (static_cast<double>(span) + 1.0);
        } else {
          span = options.window;
          alpha = 2.0 / (static_cast<double>(span) + 1.0);
        }
      }
      if (!(alpha > 0.0) || !(alpha <= 1.0)) {
        throw std::invalid_argument("--alpha must be in (0, 1]");
      }
      result.alpha = alpha;
      result.param = span;
      result.label = "ema(span=" + std::to_string(span) +
                     ", alpha=" + fmtFixed(alpha, 4) + ")";
      result.values.assign(x.size(), kNaN);
      double previous = 0.0;
      bool seeded = false;
      for (std::size_t t = 0; t < x.size(); ++t) {
        if (!std::isfinite(x[t])) {
          if (seeded) result.values[t] = previous;  // carry forward
          continue;
        }
        if (!seeded) {
          previous = x[t];
          seeded = true;
        } else {
          previous = alpha * x[t] + (1.0 - alpha) * previous;
        }
        result.values[t] = previous;
      }
      break;
    }
    case MaMethod::Diff: {
      if (options.lag < 1) {
        throw std::invalid_argument("--lag must be a positive integer");
      }
      result.param = options.lag;
      result.label = "diff(lag=" + std::to_string(options.lag) + ")";
      result.values.assign(x.size(), kNaN);
      for (std::size_t t = 0; t < x.size(); ++t) {
        if (t < static_cast<std::size_t>(options.lag)) continue;
        const double a = x[t];
        const double b = x[t - static_cast<std::size_t>(options.lag)];
        if (std::isfinite(a) && std::isfinite(b)) result.values[t] = a - b;
      }
      break;
    }
  }
  return result;
}

std::string renderSeries(const MaResult& result, const std::vector<double>& original,
                         const std::string& columnName, std::size_t maxRows, bool unicode,
                         bool color) {
  (void)color;
  Table table;
  table.setTitle("Series: " + columnName);
  table.addColumn("#", Align::Right);
  table.addColumn("value", Align::Right);
  table.addColumn(result.label, Align::Right);

  auto cell = [](double v) { return std::isfinite(v) ? fmtNumber(v, 6) : std::string("n/a"); };
  const std::size_t n = result.values.size();
  const bool elide = n > maxRows && maxRows >= 4;
  const std::size_t half = elide ? maxRows / 2 : maxRows;

  for (std::size_t i = 0; i < n; ++i) {
    if (elide && i >= half && i < n - half) {
      if (i == half) {
        table.addRow({"...", std::to_string(n - maxRows) + " more rows", "..."});
      }
      continue;
    }
    const double raw = i < original.size() ? original[i] : kNaN;
    table.addRow({std::to_string(i + 1), cell(raw), cell(result.values[i])});
  }
  table.setFooter("n = " + std::to_string(n) + "; n/a = undefined (insufficient history or "
                                                "missing input)");
  return table.render(unicode);
}

void seriesToJson(json::Writer& writer, const MaResult& result, const std::string& columnName) {
  writer.beginObject();
  writer.kv("column", columnName);
  writer.kv("method", maMethodName(result.method));
  writer.kv("param", result.param);
  writer.kv("alpha", result.alpha);
  writer.kv("label", result.label);
  writer.key("values");
  writer.beginArray();
  for (double v : result.values) writer.value(v);
  writer.endArray();
  writer.endObject();
}

}  // namespace se
