// statsengine — correlation.cpp
//
// Correlation coefficients, rank transform and matrix rendering.
// See correlation.hpp for the statistical contracts.
#include "statsengine/correlation.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"
#include "statsengine/hypothesis.hpp"
#include "statsengine/table.hpp"
#include "statsengine/util.hpp"

namespace se {

std::vector<double> rankAverage(const std::vector<double>& values) {
  const std::size_t n = values.size();
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&values](std::size_t a, std::size_t b) { return values[a] < values[b]; });

  std::vector<double> ranks(n, 0.0);
  std::size_t i = 0;
  while (i < n) {
    std::size_t j = i;
    while (j + 1 < n && values[order[j + 1]] == values[order[i]]) ++j;
    const double average = 0.5 * (static_cast<double>(i) + static_cast<double>(j)) + 1.0;
    for (std::size_t k = i; k <= j; ++k) ranks[order[k]] = average;
    i = j + 1;
  }
  return ranks;
}

double pearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = std::min(x.size(), y.size());
  if (n == 0) return kNaN;
  RunningMoments mx;
  RunningMoments my;
  for (std::size_t i = 0; i < n; ++i) {
    mx.push(x[i]);
    my.push(y[i]);
  }
  const double meanX = mx.mean();
  const double meanY = my.mean();
  double sxy = 0.0;
  double sxx = 0.0;
  double syy = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double dx = x[i] - meanX;
    const double dy = y[i] - meanY;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return kNaN;
  return sxy / std::sqrt(sxx * syy);
}

double spearmanCorrelation(const std::vector<double>& x, const std::vector<double>& y) {
  return pearsonCorrelation(rankAverage(x), rankAverage(y));
}

double kendallTau(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = std::min(x.size(), y.size());
  if (n < 2) return kNaN;
  long long concordant = 0;
  long long discordant = 0;
  long long tiesX = 0;  // Σ t(t-1)/2 over tie groups in x
  long long tiesY = 0;
  for (std::size_t i = 0; i < n; ++i) {
    long long tx = 0;
    long long ty = 0;
    for (std::size_t j = i + 1; j < n; ++j) {
      if (x[i] < x[j] && y[i] < y[j]) ++concordant;
      else if (x[i] > x[j] && y[i] > y[j]) ++concordant;
      else if (x[i] < x[j] && y[i] > y[j]) ++discordant;
      else if (x[i] > x[j] && y[i] < y[j]) ++discordant;
      else {
        if (x[i] == x[j]) ++tx;
        if (y[i] == y[j]) ++ty;
      }
    }
    tiesX += tx * (tx + 1) / 2;
    tiesY += ty * (ty + 1) / 2;
  }
  const double total = 0.5 * static_cast<double>(n) * static_cast<double>(n - 1);
  const double denomX = total - static_cast<double>(tiesX);
  const double denomY = total - static_cast<double>(tiesY);
  if (denomX <= 0.0 || denomY <= 0.0) return kNaN;
  return static_cast<double>(concordant - discordant) / std::sqrt(denomX * denomY);
}

double correlationPValue(double r, long long n) {
  if (n < 3 || !std::isfinite(r)) return kNaN;
  const double absR = std::fabs(r);
  if (absR >= 1.0) return 0.0;
  const double df = static_cast<double>(n - 2);
  const double t = r * std::sqrt(df / (1.0 - r * r));
  return studentTTwoTailedP(t, df);
}

double kendallPValue(double tau, long long n) {
  if (n < 2 || !std::isfinite(tau)) return kNaN;
  const double nn = static_cast<double>(n);
  const double variance = 2.0 * (2.0 * nn + 5.0) / (9.0 * nn * (nn - 1.0));
  if (variance <= 0.0) return kNaN;
  const double z = tau / std::sqrt(variance);
  return 2.0 * (1.0 - normalCdf(std::fabs(z)));
}

CorrelationMatrix computeCorrelationMatrix(const DataFrame& df,
                                           const std::vector<std::string>& columns,
                                           const std::string& method) {
  std::string key = toLower(trim(method));
  if (key.empty()) key = "pearson";
  if (key != "pearson" && key != "spearman" && key != "kendall") {
    throw std::invalid_argument("unknown correlation method '" + method +
                                "' (expected pearson, spearman or kendall)");
  }

  std::vector<std::string> names = columns;
  if (names.empty()) {
    for (std::size_t idx : df.numericColumnIndices()) names.push_back(df.column(idx).name());
  }
  if (names.size() < 2) {
    throw std::runtime_error("correlation needs at least two numeric columns");
  }

  std::vector<const Column*> cols;
  cols.reserve(names.size());
  for (const std::string& name : names) {
    const Column& col = df.column(name);
    if (col.kind() != ColumnKind::Numeric) {
      throw ColumnError(name, "correlation requires numeric columns but '" + name + "' is " +
                                  columnKindName(col.kind()));
    }
    cols.push_back(&col);
  }

  const std::size_t k = cols.size();
  CorrelationMatrix matrix;
  matrix.method = key;
  matrix.names = names;
  matrix.r.assign(k, std::vector<double>(k, kNaN));
  matrix.p.assign(k, std::vector<double>(k, kNaN));
  matrix.n.assign(k, std::vector<long long>(k, 0));

  for (std::size_t i = 0; i < k; ++i) {
    for (std::size_t j = i; j < k; ++j) {
      if (i == j) {
        long long count = 0;
        for (std::size_t row = 0; row < cols[i]->size(); ++row) {
          if (!cols[i]->isMissing(row)) ++count;
        }
        matrix.r[i][j] = count > 0 ? 1.0 : kNaN;
        matrix.p[i][j] = count > 0 ? 0.0 : kNaN;
        matrix.n[i][j] = count;
        continue;
      }
      // Pairwise-complete observations.
      std::vector<double> xs;
      std::vector<double> ys;
      for (std::size_t row = 0; row < cols[i]->size(); ++row) {
        if (cols[i]->isMissing(row) || cols[j]->isMissing(row)) continue;
        xs.push_back(cols[i]->value(row));
        ys.push_back(cols[j]->value(row));
      }
      const long long pairCount = static_cast<long long>(xs.size());
      matrix.n[i][j] = pairCount;
      matrix.n[j][i] = pairCount;
      double r = kNaN;
      if (pairCount >= 2) {
        if (key == "pearson") r = pearsonCorrelation(xs, ys);
        else if (key == "spearman") r = spearmanCorrelation(xs, ys);
        else r = kendallTau(xs, ys);
      }
      double p = kNaN;
      if (std::isfinite(r)) p = key == "kendall" ? kendallPValue(r, pairCount)
                                                 : correlationPValue(r, pairCount);
      matrix.r[i][j] = r;
      matrix.r[j][i] = r;
      matrix.p[i][j] = p;
      matrix.p[j][i] = p;
    }
  }
  return matrix;
}

namespace {

std::string starString(double p, bool color) {
  if (!std::isfinite(p)) return "";
  if (p < 0.001) return color ? colors::dim("***") : "***";
  if (p < 0.01) return color ? colors::dim("**") : "**";
  if (p < 0.05) return color ? colors::dim("*") : "*";
  return "";
}

std::string cellString(double r, double p, bool diagonal, bool color) {
  if (diagonal) return "1.000";
  if (!std::isfinite(r)) return "n/a";
  std::string text = fmtFixed(r, 3);
  if (color && std::isfinite(r) && std::fabs(r) >= 0.7) text = colors::bold(text);
  return text + starString(p, color);
}

}  // namespace

std::string renderCorrelationMatrix(const CorrelationMatrix& matrix, bool unicode, bool color) {
  Table table;
  table.setTitle("Correlation matrix (" + matrix.method + ")");
  table.addColumn("", Align::Left);
  for (const std::string& name : matrix.names) table.addColumn(name, Align::Right);

  const std::size_t k = matrix.names.size();
  for (std::size_t i = 0; i < k; ++i) {
    std::vector<std::string> row{matrix.names[i]};
    for (std::size_t j = 0; j < k; ++j) {
      row.push_back(cellString(matrix.r[i][j], matrix.p[i][j], i == j, color));
    }
    table.addRow(row);
  }
  table.setFooter("pairwise-complete observations; * p<0.05  ** p<0.01  *** p<0.001; "
                  "n/a = undefined (constant column or n<3)");
  return table.render(unicode);
}

void correlationToJson(json::Writer& writer, const CorrelationMatrix& matrix) {
  writer.beginObject();
  writer.kv("method", matrix.method);
  writer.key("columns");
  writer.beginArray();
  for (const std::string& name : matrix.names) writer.value(name);
  writer.endArray();

  writer.key("r");
  writer.beginArray();
  for (const std::vector<double>& row : matrix.r) {
    writer.beginArray();
    for (double v : row) writer.value(v);
    writer.endArray();
  }
  writer.endArray();

  writer.key("p");
  writer.beginArray();
  for (const std::vector<double>& row : matrix.p) {
    writer.beginArray();
    for (double v : row) writer.value(v);
    writer.endArray();
  }
  writer.endArray();

  writer.key("n");
  writer.beginArray();
  for (const std::vector<long long>& row : matrix.n) {
    writer.beginArray();
    for (long long v : row) writer.value(v);
    writer.endArray();
  }
  writer.endArray();
  writer.endObject();
}

}  // namespace se
