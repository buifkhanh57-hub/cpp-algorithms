// statsengine — describe.cpp
//
// Implementation of online moments, quantiles and the describe report.
// See describe.hpp for formulas and contract.
#include "statsengine/describe.hpp"

#include "statsengine/colors.hpp"
#include "statsengine/table.hpp"

#include <algorithm>
#include <cmath>

namespace se {

// --------------------------------------------------------- RunningMoments --

void RunningMoments::push(double x) {
  if (!std::isfinite(x)) return;
  ++n_;
  const double delta = x - mean_;
  const double deltaN = delta / static_cast<double>(n_);
  const double deltaN2 = deltaN * deltaN;
  const double term1 = delta * deltaN * static_cast<double>(n_ - 1);
  const double m2prev = m2_;
  const double m3prev = m3_;
  mean_ += deltaN;
  m4_ += term1 * deltaN2 *
            (static_cast<double>(n_) * static_cast<double>(n_) - 3.0 * static_cast<double>(n_) +
             3.0) +
        6.0 * deltaN2 * m2prev - 4.0 * deltaN * m3prev;
  m3_ += term1 * deltaN * (static_cast<double>(n_) - 2.0) - 3.0 * deltaN * m2prev;
  m2_ += term1;
}

void RunningMoments::merge(const RunningMoments& other) {
  if (other.n_ == 0) return;
  if (n_ == 0) {
    *this = other;
    return;
  }
  const double nA = static_cast<double>(n_);
  const double nB = static_cast<double>(other.n_);
  const double n = nA + nB;
  const double delta = other.mean_ - mean_;
  const double delta2 = delta * delta;
  const double delta3 = delta2 * delta;
  const double delta4 = delta2 * delta2;

  // Pébay / Schubert-Gertz parallel merge. m4 and m3 must be updated before
  // m2 (and everything before mean_), because they use the old accumulators.
  m4_ += other.m4_ + delta4 * nA * nB * (nA * nA - nA * nB + nB * nB) / (n * n * n) +
         6.0 * delta2 / (n * n) * (nA * nA * other.m2_ + nB * nB * m2_) +
         4.0 * delta / n * (nA * other.m3_ - nB * m3_);
  m3_ += other.m3_ + delta3 * nA * nB * (nA - nB) / (n * n) +
         3.0 * delta / n * (nA * other.m2_ - nB * m2_);
  m2_ += other.m2_ + delta2 * nA * nB / n;
  mean_ += delta * nB / n;
  n_ = static_cast<long long>(n);
}

double RunningMoments::populationSkewness() const noexcept {
  if (n_ < 2 || m2_ <= 0.0) return kNaN;
  const double n = static_cast<double>(n_);
  const double m2 = m2_ / n;
  const double m3 = m3_ / n;
  return m3 / (m2 * std::sqrt(m2));
}

double RunningMoments::sampleSkewness() const noexcept {
  if (n_ < 3) return kNaN;
  const double n = static_cast<double>(n_);
  const double g1 = populationSkewness();
  if (!std::isfinite(g1)) return kNaN;
  return std::sqrt(n * (n - 1.0)) / (n - 2.0) * g1;
}

double RunningMoments::populationExcessKurtosis() const noexcept {
  if (n_ < 2 || m2_ <= 0.0) return kNaN;
  const double n = static_cast<double>(n_);
  const double m2 = m2_ / n;
  const double m4 = m4_ / n;
  return m4 / (m2 * m2) - 3.0;
}

double RunningMoments::sampleExcessKurtosis() const noexcept {
  if (n_ < 4) return kNaN;
  const double n = static_cast<double>(n_);
  const double s2 = m2_ / (n - 1.0);  // sample variance
  if (s2 <= 0.0) return kNaN;
  const double sum4 = m4_;  // Σ (x - x̄)⁴
  const double first = n * (n + 1.0) / ((n - 1.0) * (n - 2.0) * (n - 3.0)) * (sum4 / (s2 * s2));
  const double second = 3.0 * (n - 1.0) * (n - 1.0) / ((n - 2.0) * (n - 3.0));
  return first - second;
}

// ------------------------------------------------------------- quantiles --

double quantileSorted(const std::vector<double>& sorted, double q) {
  if (sorted.empty()) return kNaN;
  if (sorted.size() == 1) return sorted.front();
  q = std::min(1.0, std::max(0.0, q));
  const double h = static_cast<double>(sorted.size() - 1) * q;
  const std::size_t lo = static_cast<std::size_t>(std::floor(h));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(h));
  if (lo == hi || hi >= sorted.size()) return sorted[lo];
  return sorted[lo] + (h - static_cast<double>(lo)) * (sorted[hi] - sorted[lo]);
}

// ----------------------------------------------------------------- stats --

DescriptiveStats describeValues(const std::string& name, const std::vector<double>& values) {
  DescriptiveStats stats;
  stats.name = name;

  std::vector<double> sorted;
  sorted.reserve(values.size());
  RunningMoments moments;
  for (double v : values) {
    if (!std::isfinite(v)) continue;
    sorted.push_back(v);
    moments.push(v);
    stats.sum += v;
  }
  stats.missing = values.size() - sorted.size();
  stats.count = static_cast<long long>(sorted.size());
  if (sorted.empty()) return stats;

  std::sort(sorted.begin(), sorted.end());
  stats.mean = moments.mean();
  stats.stddev = moments.populationStddev();
  stats.sampleStddev = moments.sampleStddev();
  stats.min = sorted.front();
  stats.max = sorted.back();
  stats.q1 = quantileSorted(sorted, 0.25);
  stats.median = quantileSorted(sorted, 0.50);
  stats.q3 = quantileSorted(sorted, 0.75);
  stats.iqr = stats.q3 - stats.q1;
  stats.skewness = moments.count() >= 3 ? moments.sampleSkewness() : kNaN;
  stats.kurtosis = moments.count() >= 4 ? moments.sampleExcessKurtosis() : kNaN;
  stats.cv = stats.mean != 0.0 ? stats.sampleStddev / std::fabs(stats.mean) : kNaN;
  stats.valid = true;
  return stats;
}

DescriptiveStats describeColumn(const Column& column) {
  return describeValues(column.name(), column.numeric());
}

namespace {

std::string statCell(double value) {
  return std::isfinite(value) ? fmtNumber(value, 6) : std::string("n/a");
}

}  // namespace

std::string renderDescribeTable(const std::vector<DescriptiveStats>& stats, bool unicode,
                                bool color) {
  (void)color;
  Table table;
  table.setTitle("Descriptive statistics");
  table.addColumn("metric", Align::Left);
  for (const DescriptiveStats& s : stats) table.addColumn(s.name, Align::Right);

  auto numericRow = [&](const char* label, double (*pick)(const DescriptiveStats&)) {
    std::vector<std::string> row{std::string(label)};
    for (const DescriptiveStats& s : stats) row.push_back(statCell(pick(s)));
    table.addRow(row);
  };

  std::vector<std::string> countRow{"count"};
  std::vector<std::string> missingRow{"missing"};
  for (const DescriptiveStats& s : stats) {
    countRow.push_back(std::to_string(s.count));
    missingRow.push_back(std::to_string(s.missing));
  }
  table.addRow(countRow);
  table.addRow(missingRow);
  numericRow("mean", [](const DescriptiveStats& s) { return s.mean; });
  numericRow("stddev", [](const DescriptiveStats& s) { return s.sampleStddev; });
  numericRow("min", [](const DescriptiveStats& s) { return s.min; });
  numericRow("p25", [](const DescriptiveStats& s) { return s.q1; });
  numericRow("p50", [](const DescriptiveStats& s) { return s.median; });
  numericRow("p75", [](const DescriptiveStats& s) { return s.q3; });
  numericRow("max", [](const DescriptiveStats& s) { return s.max; });
  numericRow("iqr", [](const DescriptiveStats& s) { return s.iqr; });
  numericRow("skewness", [](const DescriptiveStats& s) { return s.skewness; });
  numericRow("kurtosis", [](const DescriptiveStats& s) { return s.kurtosis; });
  numericRow("cv", [](const DescriptiveStats& s) { return s.cv; });

  table.setFooter(
      "stddev uses ddof=1; skewness/kurtosis are bias-corrected G1/G2; n/a = undefined");
  return table.render(unicode);
}

void describeToJson(json::Writer& writer, const DescriptiveStats& stats) {
  writer.beginObject();
  writer.kv("column", stats.name);
  writer.kv("count", stats.count);
  writer.kvCount("missing", stats.missing);
  writer.kv("sum", stats.sum);
  writer.kv("mean", stats.mean);
  writer.kv("stddev", stats.sampleStddev);
  writer.kv("min", stats.min);
  writer.kv("q1", stats.q1);
  writer.kv("median", stats.median);
  writer.kv("q3", stats.q3);
  writer.kv("max", stats.max);
  writer.kv("iqr", stats.iqr);
  writer.kv("skewness", stats.skewness);
  writer.kv("kurtosis", stats.kurtosis);
  writer.kv("cv", stats.cv);
  writer.endObject();
}

}  // namespace se
