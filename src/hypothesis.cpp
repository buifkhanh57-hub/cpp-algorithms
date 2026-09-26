// statsengine — hypothesis.cpp
//
// Distribution functions (incomplete beta / gamma, Student t, chi-square,
// normal) and the t / chi-square test suite. See hypothesis.hpp.
#include "statsengine/hypothesis.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"
#include "statsengine/table.hpp"
#include "statsengine/util.hpp"

namespace se {

namespace {
constexpr int kMaxIterations = 300;
constexpr double kEps = 3e-15;
constexpr double kFpMin = 1e-300;
}  // namespace

// --------------------------------------------------------- distributions --

double logGamma(double x) {
  static const double coefficients[9] = {0.99999999999980993,  676.5203681218851,
                                         -1259.1392167224028,  771.32342877765313,
                                         -176.61502916214059,  12.507343278686905,
                                         -0.13857109526572012, 9.9843695780195716e-6,
                                         1.5056327351493116e-7};
  if (x < 0.5) {
    // Reflection formula for the small-x branch.
    return std::log(M_PI) - std::log(std::fabs(std::sin(M_PI * x))) - logGamma(1.0 - x);
  }
  x -= 1.0;
  double a = coefficients[0];
  for (int i = 1; i < 9; ++i) a += coefficients[i] / (x + static_cast<double>(i));
  const double t = x + 7.5;
  return 0.5 * std::log(2.0 * M_PI) + (x + 0.5) * std::log(t) - t + std::log(a);
}

namespace {
double betaContinuedFraction(double a, double b, double x) {
  const double qab = a + b;
  const double qap = a + 1.0;
  const double qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < kFpMin) d = kFpMin;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= kMaxIterations; ++m) {
    const int m2 = 2 * m;
    double aa = static_cast<double>(m) * (b - static_cast<double>(m)) * x /
                ((qam + static_cast<double>(m2)) * (a + static_cast<double>(m2)));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kFpMin) d = kFpMin;
    c = 1.0 + aa / c;
    if (std::fabs(c) < kFpMin) c = kFpMin;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + static_cast<double>(m)) * (qab + static_cast<double>(m)) * x /
         ((a + static_cast<double>(m2)) * (qap + static_cast<double>(m2)));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kFpMin) d = kFpMin;
    c = 1.0 + aa / c;
    if (std::fabs(c) < kFpMin) c = kFpMin;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < kEps) return h;
  }
  throw std::runtime_error("incomplete beta: continued fraction did not converge");
}
}  // namespace

double incompleteBeta(double a, double b, double x) {
  if (a <= 0.0 || b <= 0.0) throw std::invalid_argument("incompleteBeta: a and b must be > 0");
  if (std::isnan(x)) return kNaN;
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  const double lnBeta =
      logGamma(a + b) - logGamma(a) - logGamma(b) + a * std::log(x) + b * std::log(1.0 - x);
  const double factor = std::exp(lnBeta);
  if (x < (a + 1.0) / (a + b + 2.0)) {
    return factor * betaContinuedFraction(a, b, x) / a;
  }
  return 1.0 - factor * betaContinuedFraction(b, a, 1.0 - x) / b;
}

namespace {
double gammaSeries(double a, double x) {  // regularized P(a, x)
  double sum = 1.0 / a;
  double del = sum;
  double ap = a;
  for (int n = 1; n <= kMaxIterations; ++n) {
    ++ap;
    del *= x / ap;
    sum += del;
    if (std::fabs(del) < std::fabs(sum) * kEps) break;
  }
  return sum * std::exp(-x + a * std::log(x) - logGamma(a));
}

double gammaContinuedFraction(double a, double x) {  // regularized Q(a, x)
  double b = x + 1.0 - a;
  double c = 1.0 / kFpMin;
  double d = 1.0 / b;
  double h = d;
  for (int i = 1; i <= kMaxIterations; ++i) {
    const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < kFpMin) d = kFpMin;
    c = b + an / c;
    if (std::fabs(c) < kFpMin) c = kFpMin;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < kEps) break;
  }
  return std::exp(-x + a * std::log(x) - logGamma(a)) * h;
}
}  // namespace

double chiSquareSf(double x, double df) {
  if (df <= 0.0) throw std::invalid_argument("chiSquareSf: df must be > 0");
  if (std::isnan(x) || x <= 0.0) return 1.0;
  const double a = df / 2.0;
  const double xx = x / 2.0;
  const double q = xx < a + 1.0 ? 1.0 - gammaSeries(a, xx) : gammaContinuedFraction(a, xx);
  return std::min(1.0, std::max(0.0, q));
}

double normalCdf(double z) { return 0.5 * (1.0 + std::erf(z / std::sqrt(2.0))); }

double studentTCdf(double t, double df) {
  if (df <= 0.0) throw std::invalid_argument("studentTCdf: df must be > 0");
  if (t == 0.0) return 0.5;
  const double x = df / (df + t * t);
  const double halfTail = 0.5 * incompleteBeta(df / 2.0, 0.5, x);
  return t > 0.0 ? 1.0 - halfTail : halfTail;
}

double studentTTwoTailedP(double t, double df) {
  if (df <= 0.0) throw std::invalid_argument("studentTTwoTailedP: df must be > 0");
  if (!std::isfinite(t)) return t > 0.0 || t < 0.0 ? 0.0 : kNaN;
  const double x = df / (df + t * t);
  return incompleteBeta(df / 2.0, 0.5, x);
}

double studentTQuantile(double p, double df) {
  if (!(p > 0.0) || !(p < 1.0)) {
    throw std::invalid_argument("studentTQuantile: p must be in (0, 1)");
  }
  double lo = -1000.0;
  double hi = 1000.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (studentTCdf(mid, df) < p) {
      lo = mid;
    } else {
      hi = mid;
    }
    if (hi - lo < 1e-12) break;
  }
  return 0.5 * (lo + hi);
}

// ------------------------------------------------------------------ tests --

namespace {

TestResult baseResult(TestKind kind, const std::string& name, const std::string& hypothesis) {
  TestResult result;
  result.kind = kind;
  result.testName = name;
  result.hypothesis = hypothesis;
  result.effectName = "Cohen's d";
  return result;
}

}  // namespace

double cohensD(const std::vector<double>& a, const std::vector<double>& b) {
  RunningMoments ra;
  RunningMoments rb;
  for (double v : a) ra.push(v);
  for (double v : b) rb.push(v);
  const double na = static_cast<double>(ra.count());
  const double nb = static_cast<double>(rb.count());
  if (na + nb < 3.0) return kNaN;
  const double va = ra.sampleVariance();
  const double vb = rb.sampleVariance();
  const double pooled = std::sqrt(((na - 1.0) * va + (nb - 1.0) * vb) / (na + nb - 2.0));
  if (!(pooled > 0.0)) return kNaN;
  return (ra.mean() - rb.mean()) / pooled;
}

TestResult oneSampleT(const std::vector<double>& x, double mu) {
  std::vector<double> values;
  for (double v : x) {
    if (std::isfinite(v)) values.push_back(v);
  }
  if (values.size() < 2) {
    throw std::runtime_error("one-sample t-test needs at least 2 numeric values");
  }
  RunningMoments moments;
  for (double v : values) moments.push(v);

  const double n = static_cast<double>(values.size());
  const double mean = moments.mean();
  const double sd = moments.sampleStddev();
  const double se = sd / std::sqrt(n);

  TestResult result = baseResult(TestKind::OneSampleT, "One-sample t-test",
                                 "mean(x) = " + fmtNumber(mu, 6));
  result.n1 = moments.count();
  result.meanDiff = mean - mu;
  if (se > 0.0) {
    result.statistic = result.meanDiff / se;
  } else {
    result.statistic = result.meanDiff == 0.0 ? 0.0 : result.meanDiff * 1e300;
  }
  result.df = n - 1.0;
  result.pValue = studentTTwoTailedP(result.statistic, result.df);
  const double critical = studentTQuantile(0.975, result.df);
  result.ciLow = result.meanDiff - critical * se;
  result.ciHigh = result.meanDiff + critical * se;
  result.effectSize = sd > 0.0 ? result.meanDiff / sd : kNaN;
  result.details.push_back("sample mean = " + fmtNumber(mean, 6) +
                           ", sample sd (ddof=1) = " + fmtNumber(sd, 6));
  result.details.push_back("H1: mean(x) != " + fmtNumber(mu, 6) + " (two-tailed)");
  return result;
}

TestResult twoSampleT(const std::vector<double>& a, const std::vector<double>& b, bool welch) {
  std::vector<double> va, vb;
  for (double v : a) if (std::isfinite(v)) va.push_back(v);
  for (double v : b) if (std::isfinite(v)) vb.push_back(v);
  if (va.size() < 2 || vb.size() < 2) {
    throw std::runtime_error("two-sample t-test needs at least 2 numeric values per group");
  }
  RunningMoments ma, mb;
  for (double v : va) ma.push(v);
  for (double v : vb) mb.push(v);
  const double na = static_cast<double>(ma.count());
  const double nb = static_cast<double>(mb.count());
  const double meanA = ma.mean();
  const double meanB = mb.mean();
  const double varA = ma.sampleVariance();
  const double varB = mb.sampleVariance();

  double se = 0.0;
  double df = 0.0;
  if (welch) {
    const double termA = varA / na;
    const double termB = varB / nb;
    se = std::sqrt(termA + termB);
    const double denom = (termA * termA) / (na - 1.0) + (termB * termB) / (nb - 1.0);
    df = denom > 0.0 ? (termA + termB) * (termA + termB) / denom : na + nb - 2.0;
  } else {
    const double pooledVar = ((na - 1.0) * varA + (nb - 1.0) * varB) / (na + nb - 2.0);
    se = std::sqrt(pooledVar * (1.0 / na + 1.0 / nb));
    df = na + nb - 2.0;
  }

  TestResult result = baseResult(TestKind::TwoSampleT,
                                 welch ? "Two-sample t-test (Welch)"
                                       : "Two-sample t-test (pooled)",
                                 "mean(a) = mean(b)");
  result.n1 = ma.count();
  result.n2 = mb.count();
  result.meanDiff = meanA - meanB;
  result.statistic = se > 0.0 ? result.meanDiff / se : (result.meanDiff == 0.0 ? 0.0 : result.meanDiff * 1e300);
  result.df = df;
  result.pValue = studentTTwoTailedP(result.statistic, df);
  const double critical = studentTQuantile(0.975, df);
  result.ciLow = result.meanDiff - critical * se;
  result.ciHigh = result.meanDiff + critical * se;
  result.effectSize = cohensD(va, vb);
  result.details.push_back("mean(a) = " + fmtNumber(meanA, 6) + " (n=" + std::to_string(ma.count()) +
                           "), mean(b) = " + fmtNumber(meanB, 6) + " (n=" + std::to_string(mb.count()) +
                           ")");
  result.details.push_back(std::string(welch ? "Satterthwaite" : "pooled") +
                           " df = " + fmtNumber(df, 2));
  result.details.push_back("H1: mean(a) != mean(b) (two-tailed)");
  return result;
}

TestResult pairedT(const std::vector<double>& a, const std::vector<double>& b) {
  std::vector<double> diffs;
  const std::size_t n = std::min(a.size(), b.size());
  if (a.size() != b.size()) {
    throw std::invalid_argument("paired t-test requires two columns of equal length (" +
                                std::to_string(a.size()) + " vs " + std::to_string(b.size()) + ")");
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(a[i]) && std::isfinite(b[i])) diffs.push_back(a[i] - b[i]);
  }
  if (diffs.size() < 2) {
    throw std::runtime_error("paired t-test needs at least 2 complete pairs");
  }
  TestResult result = oneSampleT(diffs, 0.0);
  result.kind = TestKind::PairedT;
  result.testName = "Paired t-test";
  result.hypothesis = "mean(a - b) = 0";
  result.n1 = static_cast<long long>(diffs.size());
  result.n2 = 0;
  result.details.insert(result.details.begin(),
                        "paired on " + std::to_string(diffs.size()) + " complete row pairs");
  return result;
}

TestResult chiSquareGof(const std::vector<long long>& observed,
                        const std::vector<double>& expectedProbs) {
  const std::size_t k = observed.size();
  if (k < 2) {
    throw std::runtime_error("chi-square goodness-of-fit needs at least 2 categories");
  }
  std::vector<double> probs = expectedProbs;
  if (probs.empty()) {
    probs.assign(k, 1.0 / static_cast<double>(k));
  }
  if (probs.size() != k) {
    throw std::invalid_argument("expected-probability count (" + std::to_string(probs.size()) +
                                ") does not match category count (" + std::to_string(k) + ")");
  }
  long long total = 0;
  for (long long value : observed) {
    if (value < 0) throw std::invalid_argument("observed counts must be non-negative");
    total += value;
  }
  if (total <= 0) throw std::runtime_error("chi-square test needs at least one observation");
  double probSum = 0.0;
  for (double p : probs) {
    if (p < 0.0) throw std::invalid_argument("expected probabilities must be non-negative");
    probSum += p;
  }
  if (std::fabs(probSum - 1.0) > 1e-6) {
    throw std::invalid_argument("expected probabilities must sum to 1 (got " +
                                fmtFixed(probSum, 6) + ")");
  }

  double chi2 = 0.0;
  for (std::size_t i = 0; i < k; ++i) {
    const double expected = static_cast<double>(total) * probs[i];
    if (expected <= 0.0 && observed[i] > 0) {
      throw std::runtime_error("category " + std::to_string(i + 1) +
                               " has zero expected count but " + std::to_string(observed[i]) +
                               " observations");
    }
    if (expected > 0.0) {
      const double diff = static_cast<double>(observed[i]) - expected;
      chi2 += diff * diff / expected;
    }
  }

  TestResult result = baseResult(TestKind::ChiSquareGof, "Chi-square goodness-of-fit",
                                 "categories follow the expected proportions");
  result.kind = TestKind::ChiSquareGof;
  result.effectName = "phi";
  result.statistic = chi2;
  result.df = static_cast<double>(k - 1);
  result.pValue = chiSquareSf(chi2, result.df);
  result.n1 = total;
  result.meanDiff = kNaN;
  result.ciLow = kNaN;
  result.ciHigh = kNaN;
  result.effectSize = std::sqrt(chi2 / static_cast<double>(total));
  result.details.push_back("categories = " + std::to_string(k) + ", total observations = " +
                           std::to_string(total));
  result.details.push_back("H1: at least one proportion differs (right-tailed)");
  return result;
}

// ---------------------------------------------------------------- render --

std::string renderTest(const TestResult& result, double alpha, bool unicode, bool color) {
  (void)unicode;  // the key/value block is plain text in both modes
  std::vector<std::pair<std::string, std::string>> items;
  items.emplace_back("test", result.testName);
  items.emplace_back("H0", result.hypothesis);
  items.emplace_back("statistic", fmtNumber(result.statistic, 6));
  items.emplace_back("df", std::isfinite(result.df) ? fmtNumber(result.df, 2) : std::string("n/a"));
  items.emplace_back("p-value", fmtFixed(result.pValue, 6));
  items.emplace_back("alpha", fmtFixed(alpha, 3));
  if (result.kind == TestKind::ChiSquareGof) {
    items.emplace_back("n", std::to_string(result.n1));
  } else {
    items.emplace_back("n", result.n2 > 0 ? std::to_string(result.n1) + " vs " +
                                                std::to_string(result.n2)
                                          : std::to_string(result.n1));
  }
  if (std::isfinite(result.effectSize)) {
    items.emplace_back(result.effectName, fmtFixed(result.effectSize, 4));
  }
  if (std::isfinite(result.meanDiff)) {
    items.emplace_back("mean diff", fmtNumber(result.meanDiff, 6) + "  [" +
                                        fmtFixed(result.ciLow, 4) + ", " +
                                        fmtFixed(result.ciHigh, 4) + "]");
  }

  std::ostringstream out;
  out << keyValueBlock(items, color) << "\n";
  const bool significant = std::isfinite(result.pValue) && result.pValue < alpha;
  const std::string verdict = significant
                                  ? "reject H0 at alpha=" + fmtFixed(alpha, 3) +
                                        " (statistically significant)"
                                  : "fail to reject H0 at alpha=" + fmtFixed(alpha, 3);
  out << (significant ? colors::bold(colors::green("=> " + verdict))
                      : colors::bold("=> " + verdict))
      << "\n";
  for (const std::string& detail : result.details) {
    out << "  - " << detail << "\n";
  }
  return out.str();
}

void testToJson(json::Writer& writer, const TestResult& result, double alpha) {
  writer.beginObject();
  writer.kv("test", result.testName);
  writer.kv("h0", result.hypothesis);
  writer.kv("statistic", result.statistic);
  writer.kv("df", result.df);
  writer.kv("p_value", result.pValue);
  writer.kv("alpha", alpha);
  writer.kv("n1", result.n1);
  writer.kv("n2", result.n2);
  writer.kv("mean_diff", result.meanDiff);
  writer.kv("ci_low", result.ciLow);
  writer.kv("ci_high", result.ciHigh);
  writer.kv("effect_size", result.effectSize);
  writer.kv("effect_name", result.effectName);
  writer.kvBool("significant", std::isfinite(result.pValue) && result.pValue < alpha);
  writer.endObject();
}

}  // namespace se
