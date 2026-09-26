// statsengine — anova.cpp
//
// One-way ANOVA computation, F tail probability and report rendering.
// See anova.hpp for the statistical contract.
#include "statsengine/anova.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"  // RunningMoments, quantileSorted
#include "statsengine/hypothesis.hpp"  // incompleteBeta
#include "statsengine/table.hpp"

namespace se {

double fSf(double f, double d1, double d2) {
  if (!(d1 > 0.0) || !(d2 > 0.0)) {
    throw std::invalid_argument("fSf: degrees of freedom must be positive");
  }
  if (std::isnan(f)) return kNaN;
  if (f <= 0.0) return 1.0;
  if (std::isinf(f)) return 0.0;
  // P(F(d1,d2) >= f) = I_{d2/(d2 + d1 f)}(d2/2, d1/2)
  const double x = d2 / (d2 + d1 * f);
  return incompleteBeta(d2 / 2.0, d1 / 2.0, x);
}

AnovaResult oneWayAnova(
    const std::vector<std::pair<std::string, std::vector<double>>>& groups) {
  if (groups.size() < 2) {
    throw std::runtime_error("one-way ANOVA needs at least 2 groups");
  }
  std::set<std::string> seen;
  for (const auto& group : groups) {
    if (!seen.insert(group.first).second) {
      throw std::invalid_argument("duplicate group label '" + group.first + "'");
    }
  }

  AnovaResult result;
  double total = 0.0;
  long long totalN = 0;
  for (const auto& group : groups) {
    RunningMoments moments;
    long long count = 0;
    for (double v : group.second) {
      if (!std::isfinite(v)) continue;
      moments.push(v);
      total += v;
      ++count;
    }
    if (count < 2) {
      throw std::runtime_error("group '" + group.first +
                               "' has fewer than 2 finite values");
    }
    AnovaGroup summary;
    summary.name = group.first;
    summary.n = moments.count();
    summary.mean = moments.mean();
    summary.sampleStddev = moments.sampleStddev();
    result.groups.push_back(summary);
    totalN += moments.count();
  }
  if (result.groups.size() < 2) {
    throw std::runtime_error("one-way ANOVA needs at least 2 groups with >= 2 values");
  }
  if (totalN <= static_cast<long long>(result.groups.size())) {
    throw std::runtime_error("not enough observations (N=" + std::to_string(totalN) +
                             ") for " + std::to_string(result.groups.size()) + " groups");
  }

  result.k = static_cast<long long>(result.groups.size());
  result.n = totalN;
  result.grandMean = total / static_cast<double>(totalN);
  result.dfBetween = result.k - 1;
  result.dfWithin = result.n - result.k;

  // Second pass: sums of squares against group means and the grand mean.
  std::size_t g = 0;
  for (const auto& group : groups) {
    const AnovaGroup& summary = result.groups[g];
    for (double v : group.second) {
      if (!std::isfinite(v)) continue;
      const double devGroup = v - summary.mean;
      result.ssWithin += devGroup * devGroup;
      const double devGrand = v - result.grandMean;
      result.ssTotal += devGrand * devGrand;
    }
    const double meanDev = summary.mean - result.grandMean;
    result.ssBetween += static_cast<double>(summary.n) * meanDev * meanDev;
    ++g;
  }

  result.msBetween = result.ssBetween / static_cast<double>(result.dfBetween);
  result.msWithin = result.ssWithin / static_cast<double>(result.dfWithin);
  result.fStat = result.msWithin > 0.0 ? result.msBetween / result.msWithin : kNaN;
  result.pValue = std::isfinite(result.fStat)
                      ? fSf(result.fStat, static_cast<double>(result.dfBetween),
                            static_cast<double>(result.dfWithin))
                      : kNaN;
  result.etaSquared = result.ssTotal > 0.0 ? result.ssBetween / result.ssTotal : kNaN;
  const double omegaNum = result.ssBetween -
                          static_cast<double>(result.dfBetween) * result.msWithin;
  const double omegaDen = result.ssTotal + result.msWithin;
  result.omegaSquared = omegaDen > 0.0 ? omegaNum / omegaDen : kNaN;
  result.valid = true;
  return result;
}

std::string renderAnova(const AnovaResult& result, double alpha, bool unicode, bool color) {
  (void)unicode;  // tables handle the Unicode/ASCII choice themselves
  std::vector<std::pair<std::string, std::string>> items;
  items.emplace_back("groups", std::to_string(result.k));
  items.emplace_back("observations", std::to_string(result.n));
  items.emplace_back("grand mean", fmtNumber(result.grandMean, 6));
  items.emplace_back("F", fmtFixed(result.fStat, 4));
  items.emplace_back("df", std::to_string(result.dfBetween) + " and " +
                               std::to_string(result.dfWithin));
  items.emplace_back("p-value", fmtFixed(result.pValue, 6));
  items.emplace_back("alpha", fmtFixed(alpha, 3));
  if (std::isfinite(result.etaSquared)) items.emplace_back("eta-squared", fmtFixed(result.etaSquared, 4));
  if (std::isfinite(result.omegaSquared)) {
    items.emplace_back("omega-squared", fmtFixed(result.omegaSquared, 4));
  }

  std::ostringstream out;
  out << keyValueBlock(items, color) << "\n";

  const bool significant = std::isfinite(result.pValue) && result.pValue < alpha;
  const std::string verdict = significant
                                  ? "reject H0 at alpha=" + fmtFixed(alpha, 3) +
                                        " (at least one group mean differs)"
                                  : "fail to reject H0 at alpha=" + fmtFixed(alpha, 3) +
                                        " (no evidence of mean differences)";
  out << (significant ? colors::bold(colors::green("=> " + verdict))
                      : colors::bold("=> " + verdict))
      << "\n\n";

  Table groupTable;
  groupTable.setTitle("Group summaries");
  groupTable.addColumn("group", Align::Left);
  groupTable.addColumn("n", Align::Right);
  groupTable.addColumn("mean", Align::Right);
  groupTable.addColumn("sd", Align::Right);
  for (const AnovaGroup& group : result.groups) {
    groupTable.addRow({group.name, std::to_string(group.n), fmtNumber(group.mean, 6),
                       std::isfinite(group.sampleStddev) ? fmtNumber(group.sampleStddev, 6)
                                                         : std::string("n/a")});
  }
  groupTable.setFooter("H0: all group means are equal (one-way F test)");
  out << groupTable.render(unicode) << "\n";

  Table sourceTable;
  sourceTable.setTitle("ANOVA table");
  sourceTable.addColumn("source", Align::Left);
  sourceTable.addColumn("SS", Align::Right);
  sourceTable.addColumn("df", Align::Right);
  sourceTable.addColumn("MS", Align::Right);
  sourceTable.addColumn("F", Align::Right);
  sourceTable.addColumn("p-value", Align::Right);
  sourceTable.addRow({"between", fmtFixed(result.ssBetween, 4),
                      std::to_string(result.dfBetween), fmtFixed(result.msBetween, 4),
                      std::isfinite(result.fStat) ? fmtFixed(result.fStat, 4)
                                                  : std::string("n/a"),
                      std::isfinite(result.pValue) ? fmtFixed(result.pValue, 6)
                                                   : std::string("n/a")});
  sourceTable.addRow({"within", fmtFixed(result.ssWithin, 4), std::to_string(result.dfWithin),
                      fmtFixed(result.msWithin, 4), "", ""});
  sourceTable.addRow({"total", fmtFixed(result.ssTotal, 4),
                      std::to_string(result.dfBetween + result.dfWithin), "", "", ""});
  out << sourceTable.render(unicode);
  return out.str();
}

void anovaToJson(json::Writer& writer, const AnovaResult& result, double alpha,
                 const std::string& columnName, const std::string& byName) {
  writer.beginObject();
  if (!columnName.empty()) writer.kv("column", columnName);
  if (!byName.empty()) writer.kv("by", byName);
  writer.kv("groups", result.k);
  writer.kv("n", result.n);
  writer.kv("grand_mean", result.grandMean);
  writer.key("group_summaries");
  writer.beginArray();
  for (const AnovaGroup& group : result.groups) {
    writer.beginObject();
    writer.kv("group", group.name);
    writer.kv("n", group.n);
    writer.kv("mean", group.mean);
    writer.kv("stddev", group.sampleStddev);
    writer.endObject();
  }
  writer.endArray();
  writer.kv("ss_between", result.ssBetween);
  writer.kv("ss_within", result.ssWithin);
  writer.kv("ss_total", result.ssTotal);
  writer.kv("df_between", result.dfBetween);
  writer.kv("df_within", result.dfWithin);
  writer.kv("ms_between", result.msBetween);
  writer.kv("ms_within", result.msWithin);
  writer.kv("f", result.fStat);
  writer.kv("p_value", result.pValue);
  writer.kv("alpha", alpha);
  writer.kv("eta_squared", result.etaSquared);
  writer.kv("omega_squared", result.omegaSquared);
  writer.kvBool("significant", std::isfinite(result.pValue) && result.pValue < alpha);
  writer.endObject();
}

}  // namespace se
