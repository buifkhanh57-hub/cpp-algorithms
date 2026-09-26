// statsengine — tests/test_core.cpp
//
// Assert-based test suite for every statsengine module. The binary is built
// by `make test`, which links every library translation unit except
// main.cpp; this file provides its own main().
//
// Conventions:
//   * CHECK(cond) / CHECK_MSG(cond, msg) record a failure and keep going.
//   * NEAR(a, b, tol) compares doubles with an absolute tolerance.
//   * Each suite prints its own banner; main() prints a summary and exits
//     non-zero when any check failed.
//
// Hand-computed expectations below are derived from first principles
// (closed-form sums, textbook t/F/chi-square critical values, RFC 4180
// examples) so the suite doubles as executable documentation.
#include "statsengine/anova.hpp"
#include "statsengine/args.hpp"
#include "statsengine/colors.hpp"
#include "statsengine/correlation.hpp"
#include "statsengine/csv.hpp"
#include "statsengine/dataframe.hpp"
#include "statsengine/describe.hpp"
#include "statsengine/histogram.hpp"
#include "statsengine/hypothesis.hpp"
#include "statsengine/json.hpp"
#include "statsengine/outliers.hpp"
#include "statsengine/regression.hpp"
#include "statsengine/table.hpp"
#include "statsengine/timeseries.hpp"
#include "statsengine/util.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int gChecks = 0;
int gFailures = 0;
const char* gSuite = "";

void beginSuite(const char* name) {
  gSuite = name;
  std::printf("[ RUN ] %s\n", name);
}

void reportFailure(const std::string& message) {
  ++gFailures;
  std::printf("  FAIL (%s): %s\n", gSuite, message.c_str());
}

#define CHECK(cond)                                                     \
  do {                                                                  \
    ++gChecks;                                                          \
    if (!(cond)) {                                                      \
      reportFailure(std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                    " CHECK(" #cond ")");                             \
    }                                                                   \
  } while (0)

void NEAR(double actual, double expected, double tol, const char* what) {
  ++gChecks;
  if (!(std::fabs(actual - expected) <= tol)) {
    reportFailure(std::string(what) + ": expected " + std::to_string(expected) + " +/- " +
                  std::to_string(tol) + ", got " + std::to_string(actual));
  }
}

// ----------------------------------------------------------------- helpers --

/// Build a DataFrame from a raw CSV string.
se::DataFrame frameFrom(const std::string& text) {
  std::istringstream in(text);
  se::CsvReader reader(in);
  return se::DataFrame::fromCsv(reader.read(), se::CsvOptions{});
}

/// Deterministic pseudo-data for the moments tests.
std::vector<double> sampleValues(std::size_t n, std::uint64_t seed) {
  se::Random rng(seed);
  std::vector<double> values;
  values.reserve(n);
  for (std::size_t i = 0; i < n; ++i) values.push_back(rng.nextDouble() * 200.0 - 100.0);
  return values;
}

// -------------------------------------------------------------------- util --

void testUtil() {
  beginSuite("util");

  CHECK(se::split("a,b,c", ',').size() == 3);
  CHECK(se::split("", ',', false).empty());
  CHECK(se::split(",", ',', false).empty());
  CHECK(se::trim("  hello\t") == "hello");
  CHECK(se::toLower("MiXeD") == "mixed");
  CHECK(se::splitList("a, b,,c") == std::vector<std::string>({"a", "b", "c"}));
  CHECK(se::splitList("").empty());
  CHECK(se::joinList({"x", "y"}) == "x, y");
  CHECK(se::startsWith("abc", "ab"));

  double number = 0.0;
  CHECK(se::parseNumberStrict("3.25", number) && number == 3.25);
  CHECK(se::parseNumberStrict("  -2 ", number) && number == -2.0);
  CHECK(!se::parseNumberStrict("12x", number));
  CHECK(!se::parseNumberStrict("", number));
  long long integer = 0;
  CHECK(se::parseIntStrict("42", integer) && integer == 42);
  CHECK(!se::parseIntStrict("4.5", integer));

  // Dates: Hinnant days_from_civil must round-trip, and 1970-01-02 == 1.
  int y = 0, m = 0, d = 0;
  CHECK(se::parseDate("2025-06-15", y, m, d) && y == 2025 && m == 6 && d == 15);
  CHECK(se::parseDate("2024/02/29", y, m, d));  // leap day
  CHECK(!se::parseDate("2023-02-29", y, m, d));
  CHECK(!se::parseDate("2023-13-01", y, m, d));
  CHECK(!se::parseDate("2023-4-1", y, m, d));  // 1-digit fields rejected
  CHECK(se::dateToSerial(1970, 1, 1) == 0);
  CHECK(se::dateToSerial(1970, 1, 2) == 1);
  CHECK(se::dateToSerial(2025, 1, 1) == 20089);
  CHECK(se::dateToSerial(2025, 1, 31) == 20089 + 30);  // day-of-month matters
  for (long long serial : {0LL, 1LL, 100LL, 20089LL, -719468LL, 738885LL}) {
    const std::string text = se::serialToDate(serial);
    CHECK(se::parseDate(text, y, m, d));
    NEAR(static_cast<double>(se::dateToSerial(y, m, d)), static_cast<double>(serial), 0.0,
         "date round-trip");
  }
  CHECK(se::serialToDate(20089) == "2025-01-01");

  CHECK(se::fmtFixed(1.5, 2) == "1.50");
  CHECK(se::fmtFixed(-0.0001, 3) == "0.000");  // negative zero normalised
  CHECK(se::fmtNumber(3.140000) == "3.14");
  CHECK(se::fmtNumber(42.0) == "42");
  CHECK(se::fmtNumber(se::kNaN) == "nan");
  CHECK(se::fmtPercent(0.1234, 1) == "12.3%");
  CHECK(se::fmtBytes(2048) == "2 KB");

  // PRNG: deterministic for a given seed, shuffle is a permutation.
  se::Random a(42);
  se::Random b(42);
  for (int i = 0; i < 16; ++i) CHECK(a.nextU64() == b.nextU64());
  se::Random rng(7);
  std::vector<std::size_t> items{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  rng.shuffle(items);
  std::vector<bool> seen(items.size(), false);
  for (std::size_t v : items) {
    CHECK(v < 10 && !seen[v]);
    seen[v] = true;
  }
}

// --------------------------------------------------------------------- csv --

void testCsv() {
  beginSuite("csv");

  {
    std::istringstream in("a,b\n1,2\n3,4\n");
    se::CsvReader reader(in);
    se::CsvDocument doc = reader.read();
    CHECK((doc.header == std::vector<std::string>{"a", "b"}));
    CHECK(doc.rows.size() == 2 && doc.rows[0][0] == "1" && doc.rows[1][1] == "4");
  }
  {
    // RFC 4180: quoted fields with delimiter, escaped quotes and newlines.
    std::istringstream in("a,b\n\"1,5\",\"he said \"\"hi\"\"\"\n\"multi\nline\",z\n");
    se::CsvReader reader(in);
    se::CsvDocument doc = reader.read();
    CHECK(doc.rows.size() == 2);
    CHECK(doc.rows[0][0] == "1,5");
    CHECK(doc.rows[0][1] == "he said \"hi\"");
    CHECK(doc.rows[1][0] == "multi\nline");
  }
  {
    // Blank lines are skipped; CRLF normalised; BOM stripped.
    std::istringstream in("\xef\xbb\xbfname,age\r\nbob,30\r\n\r\nann,\r\n");
    se::CsvReader reader(in);
    se::CsvDocument doc = reader.read();
    CHECK(doc.header.front() == "name");
    CHECK(doc.rows.size() == 2);
  }
  {
    // Strict row width errors point at the physical line.
    std::istringstream in("a,b\n1,2\n3,4,5\n");
    se::CsvReader reader(in);
    bool threw = false;
    try {
      (void)reader.read();
    } catch (const se::CsvParseError& e) {
      threw = true;
      CHECK(e.line == 3);
    }
    CHECK(threw);
  }
  {
    // Unterminated quote reports the line where the quote opened.
    std::istringstream in("a\n\"oops\n");
    se::CsvReader reader(in);
    bool threw = false;
    try {
      (void)reader.read();
    } catch (const se::CsvParseError& e) {
      threw = true;
      CHECK(e.line == 2);
    }
    CHECK(threw);
  }
  {
    // Writer: minimal quoting + round-trip through the reader (two records
    // so the first one is not consumed as the header).
    std::ostringstream out;
    se::CsvWriter writer(out);
    writer.row({"plain", "with,comma", "with\"quote", " lead", "multi\nline"});
    writer.row({"plain", "with,comma", "with\"quote", " lead", "multi\nline"});
    const std::string text = out.str();
    // Minimal quoting: plain stays bare, the tricky fields get quoted.
    CHECK(text.find("plain,") == 0);
    CHECK(text.find("\"with,comma\"") != std::string::npos);
    CHECK(text.find("\"with\"\"quote\"") != std::string::npos);
    CHECK(text.find("\" lead\"") != std::string::npos);
    std::istringstream in(text);
    se::CsvReader reader(in);
    se::CsvDocument doc = reader.read();
    CHECK(doc.rows.size() == 1 && doc.rows[0].size() == 5);
    CHECK(doc.header.size() == 5);
    CHECK(doc.rows[0][2] == "with\"quote");
    CHECK(doc.rows[0][4] == "multi\nline");
    CHECK(doc.header[4] == "multi\nline");  // round-trip is lossless
  }
  {
    // Type inference: numbers, dates, mixed with missing tokens.
    CHECK(se::inferColumnKind({"1", "2.5", "-3"}) == se::ColumnKind::Numeric);
    CHECK(se::inferColumnKind({"1", "NA", "null"}, {"NA", "null"}) ==
          se::ColumnKind::Numeric);  // missing tokens carry no signal
    CHECK(se::inferColumnKind({"2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04",
                               "x"}) == se::ColumnKind::Date);  // exactly 80 %
    CHECK(se::inferColumnKind({"2025-01-01", "2025-01-02", "x"}) ==
          se::ColumnKind::Categorical);  // below the 80 % threshold
    CHECK(se::inferColumnKind({"red", "blue"}) == se::ColumnKind::Categorical);
    CHECK(se::inferColumnKind({}) == se::ColumnKind::Categorical);
    CHECK(se::looksLikeNumber("1e3") && !se::looksLikeNumber("1e3.5"));
    CHECK(se::looksLikeDate("2025-01-01") && !se::looksLikeDate("2025-1-1"));
  }
}

// --------------------------------------------------------------- dataframe --

void testDataFrame() {
  beginSuite("dataframe");

  {
    // Header normalisation: trim, synthesise, de-duplicate.
    se::DataFrame df = frameFrom(" a,a,,\n1,2,3,4\n5,6,7,8\n");
    CHECK(df.columnNames() ==
          std::vector<std::string>({"a", "a_2", "c3", "c4"}));
    CHECK(df.rows() == 2 && df.cols() == 4);
  }
  {
    se::DataFrame df = frameFrom("v,w\n1,x\nNA,y\n3,\n");
    const se::Column& v = df.column("v");
    CHECK(v.kind() == se::ColumnKind::Numeric);
    CHECK(v.missingCount() == 1);
    CHECK(v.isMissing(1));
    NEAR(v.value(0), 1.0, 0.0, "v[0]");
    CHECK(v.numericCompact().size() == 2);
    const se::Column& w = df.column("w");
    CHECK(w.kind() == se::ColumnKind::Categorical);
    CHECK(w.categories().size() == 2);  // "x" and "y" (missing "" skipped)
    CHECK(df.missingCells() == 2);
  }
  {
    se::DataFrame df = frameFrom("d\n2025-01-01\n2025-02-01\n");
    const se::Column& d = df.column("d");
    CHECK(d.kind() == se::ColumnKind::Date);
    NEAR(d.value(1) - d.value(0), 31.0, 0.0, "date serial difference");
  }
  {
    se::DataFrame df = frameFrom("a,b\n1,2\n3,4\n5,6\n");
    std::vector<std::size_t> byName = df.resolveColumns("b", false);
    CHECK(byName.size() == 1 && byName[0] == 1);
    std::vector<std::size_t> byIndex = df.resolveColumns("2,1", false);
    CHECK(byIndex.size() == 2 && byIndex[0] == 1 && byIndex[1] == 0);
    std::vector<std::size_t> all = df.resolveColumns("", true);
    CHECK(all.size() == 2);
    bool threw = false;
    try {
      (void)df.resolveColumns("zzz", false);
    } catch (const se::ColumnError&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    se::DataFrame df = frameFrom("a,b\n1,\n2,4\n,6\n3,7\n");
    const std::vector<std::size_t> complete = df.completeRows({"a", "b"});
    CHECK((complete == std::vector<std::size_t>{1, 3}));
    se::DataFrame subset = df.takeRows({2, 0});
    CHECK(subset.rows() == 2);
    CHECK(subset.column("a").rawAt(0) == "" && subset.column("a").rawAt(1) == "1");
    se::DataFrame cols = df.selectColumns({1});
    CHECK(cols.cols() == 1 && cols.column(0).name() == "b");
    CHECK(df.bytesUsed() > 0);
    CHECK(df.rowAt(1) == std::vector<std::string>({"2", "4"}));
  }
  {
    // Integer-like detection and maxDecimals display hint.
    se::DataFrame ints = frameFrom("a\n1\n2\n3\n");
    CHECK(ints.column("a").isIntegerLike());
    se::DataFrame decimals = frameFrom("a\n1.5\n2.25\n");
    CHECK(!decimals.column("a").isIntegerLike());
    CHECK(decimals.column("a").maxDecimals() == 2);
  }
  {
    // No header -> columns named c1..cN (doc.header stays empty).
    se::CsvDocument doc;
    doc.rows = {{"1", "2"}, {"3", "4"}};
    se::CsvOptions options;
    options.hasHeader = false;
    se::DataFrame df = se::DataFrame::fromCsv(doc, options);
    CHECK(df.column(0).name() == "c1" && df.column(1).name() == "c2");
    CHECK(df.rows() == 2);
  }
}

// ---------------------------------------------------------------- describe --

void testDescribe() {
  beginSuite("describe");

  // Welford moments vs naive two-pass computation.
  {
    const std::vector<double> values = sampleValues(500, 12345);
    se::RunningMoments moments;
    for (double v : values) moments.push(v);
    double sum = 0.0;
    for (double v : values) sum += v;
    const double naiveMean = sum / static_cast<double>(values.size());
    double sq = 0.0;
    for (double v : values) sq += (v - naiveMean) * (v - naiveMean);
    const double naiveVar = sq / static_cast<double>(values.size() - 1);
    NEAR(moments.mean(), naiveMean, 1e-9, "welford mean");
    NEAR(moments.sampleVariance(), naiveVar, 1e-6, "welford sample variance");
    CHECK(moments.count() == 500);
  }
  {
    // Parallel merge equals a single pass over the same data.
    const std::vector<double> values = sampleValues(321, 999);
    se::RunningMoments whole;
    for (double v : values) whole.push(v);
    se::RunningMoments left;
    se::RunningMoments right;
    for (std::size_t i = 0; i < values.size() / 2; ++i) left.push(values[i]);
    for (std::size_t i = values.size() / 2; i < values.size(); ++i) right.push(values[i]);
    left.merge(right);
    NEAR(left.mean(), whole.mean(), 1e-12, "merge mean");
    NEAR(left.m2(), whole.m2(), 1e-4, "merge m2");
    NEAR(left.m3(), whole.m3(), 1e-6 + std::fabs(whole.m3()) * 1e-9, "merge m3");
    NEAR(left.m4(), whole.m4(), 1e-2 + std::fabs(whole.m4()) * 1e-9, "merge m4");
    CHECK(left.count() == whole.count());
  }
  {
    // Quantiles with linear interpolation (numpy 'linear').
    std::vector<double> sorted{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    NEAR(se::quantileSorted(sorted, 0.25), 3.25, 1e-12, "q1 of 1..10");
    NEAR(se::quantileSorted(sorted, 0.5), 5.5, 1e-12, "median of 1..10");
    NEAR(se::quantileSorted(sorted, 0.0), 1.0, 1e-12, "min quantile");
    NEAR(se::quantileSorted(sorted, 1.0), 10.0, 1e-12, "max quantile");
    CHECK(std::isnan(se::quantileSorted({}, 0.5)));
  }
  {
    // Known seven-number summary for 1..100.
    std::vector<double> values;
    for (int i = 1; i <= 100; ++i) values.push_back(static_cast<double>(i));
    se::DescriptiveStats stats = se::describeValues("v", values);
    CHECK(stats.valid && stats.count == 100 && stats.missing == 0);
    NEAR(stats.sum, 5050.0, 1e-9, "sum 1..100");
    NEAR(stats.mean, 50.5, 1e-12, "mean 1..100");
    NEAR(stats.sampleStddev, 29.011491975882016, 1e-9, "sample stddev 1..100");
    NEAR(stats.stddev, 28.86607004772212, 1e-9, "population stddev 1..100");
    NEAR(stats.q1, 25.75, 1e-12, "q1 1..100");
    NEAR(stats.median, 50.5, 1e-12, "median 1..100");
    NEAR(stats.q3, 75.25, 1e-12, "q3 1..100");
    NEAR(stats.skewness, 0.0, 1e-12, "symmetric data skewness");
    NEAR(stats.kurtosis, -1.2, 0.01, "uniform excess kurtosis");
    NEAR(stats.cv, 29.011491975882016 / 50.5, 1e-9, "cv 1..100");
  }
  {
    // Non-finite values are treated as missing.
    se::DescriptiveStats stats = se::describeValues("v", {1.0, se::kNaN, 3.0, 5.0});
    CHECK(stats.count == 3 && stats.missing == 1);
    NEAR(stats.mean, 3.0, 1e-12, "mean ignoring NaN");
    CHECK(std::fabs(stats.skewness) < 1e-12);  // symmetric, defined for n = 3
  }
  {
    se::DataFrame df = frameFrom("a,b\n1,x\n2,y\n");
    se::DescriptiveStats catStats = se::describeColumn(df.column("b"));
    CHECK(catStats.count == 0 && !catStats.valid);  // categorical: no values
  }
  {
    const se::DescriptiveStats stats = se::describeValues("v", {1.0, 2.0, 3.0, 4.0});
    const std::string table = se::renderDescribeTable({stats}, false, false);
    CHECK(table.find("mean") != std::string::npos);
    CHECK(table.find("v") != std::string::npos);
  }
}

// --------------------------------------------------------------- regression --

void testRegression() {
  beginSuite("regression");

  {
    // Exact recovery on noise-free y = 2x + 1.
    std::vector<double> y;
    std::vector<std::vector<double>> x;
    for (int i = 0; i < 10; ++i) {
      x.push_back({static_cast<double>(i)});
      y.push_back(2.0 * i + 1.0);
    }
    se::RegressionOptions options;
    se::RegressionResult model = se::fitOls(y, x, {"x"}, options);
    CHECK(model.valid && model.n == 10 && model.dfResidual == 8);
    NEAR(model.coefficients[0], 1.0, 1e-9, "intercept");
    NEAR(model.coefficients[1], 2.0, 1e-9, "slope");
    NEAR(model.r2, 1.0, 1e-12, "R^2 perfect fit");
    NEAR(model.sse, 0.0, 1e-18, "SSE perfect fit");
    NEAR(model.residualStdError, 0.0, 1e-9, "residual SE");
    CHECK(std::isnan(model.durbinWatson));  // documented: NaN on zero residuals
    CHECK(model.pValues[0] == 0.0 || std::isnan(model.pValues[0]));
    const std::vector<double> predictions =
        se::predict(model, {{5.0}, {10.0}});
    NEAR(predictions[0], 11.0, 1e-9, "predict x=5");
    NEAR(predictions[1], 21.0, 1e-9, "predict x=10");
  }
  {
    // Multivariate: y = 3 + 2a - 5b, including a collinearity failure.
    // b must be independent of a and the intercept: with b = i - 3 the design
    // satisfies a = 1.5*(b + 3), making X'X singular (fitOls throws).
    std::vector<double> y;
    std::vector<std::vector<double>> x;
    for (int i = 1; i <= 12; ++i) {
      const double a = i * 1.5;
      const double b = static_cast<double>(i % 3) - 1.0;  // cycles -1, 0, 1
      x.push_back({a, b});
      y.push_back(3.0 + 2.0 * a - 5.0 * b);
    }
    se::RegressionResult model = se::fitOls(y, x, {"a", "b"}, se::RegressionOptions{});
    NEAR(model.coefficients[0], 3.0, 1e-8, "mv intercept");
    NEAR(model.coefficients[1], 2.0, 1e-8, "mv coef a");
    NEAR(model.coefficients[2], -5.0, 1e-8, "mv coef b");
    NEAR(model.r2, 1.0, 1e-12, "mv R^2");
    // Perfect fit => SSE = 0 and F = +inf, which serialises as JSON null;
    // NaN (no signal at all) must never appear here.
    CHECK(!std::isnan(model.fStat) && model.fStat > 1e12);

    std::vector<std::vector<double>> collinear;
    std::vector<double> y2;
    for (int i = 0; i < 8; ++i) {
      collinear.push_back({static_cast<double>(i), 2.0 * i});
      y2.push_back(static_cast<double>(i));
    }
    bool threw = false;
    try {
      (void)se::fitOls(y2, collinear, {"a", "b"}, se::RegressionOptions{});
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    // No-intercept fit recovers a pure proportionality.
    std::vector<double> y;
    std::vector<std::vector<double>> x;
    for (int i = 1; i <= 6; ++i) {
      x.push_back({static_cast<double>(i)});
      y.push_back(3.0 * i);
    }
    se::RegressionOptions options;
    options.hasIntercept = false;
    se::RegressionResult model = se::fitOls(y, x, {"x"}, options);
    CHECK(model.termNames.size() == 1 && model.termNames[0] == "x");
    NEAR(model.coefficients[0], 3.0, 1e-9, "no-intercept slope");
    NEAR(model.r2, 1.0, 1e-12, "uncentered R^2");
  }
  {
    // Statistical sanity on noisy data: SE, t, p are mutually consistent.
    std::vector<double> y{2.9, 4.8, 7.2, 8.8, 11.4, 12.6, 15.1, 16.9};
    std::vector<std::vector<double>> x;
    for (int i = 0; i < 8; ++i) x.push_back({static_cast<double>(i)});
    se::RegressionResult model = se::fitOls(y, x, {"x"}, se::RegressionOptions{});
    CHECK(model.stdErrors[1] > 0.0);
    NEAR(model.tStats[1], model.coefficients[1] / model.stdErrors[1], 1e-9, "t = b/se");
    NEAR(model.pValues[1], se::studentTTwoTailedP(model.tStats[1], model.dfResidual), 1e-12,
         "coefficient p-value");
    NEAR(model.ciLow[1], model.coefficients[1] - model.stdErrors[1] * se::studentTQuantile(0.975, model.dfResidual),
         1e-9, "CI low");
    NEAR(model.ciHigh[1], model.coefficients[1] + model.stdErrors[1] * se::studentTQuantile(0.975, model.dfResidual),
         1e-9, "CI high");
    NEAR(model.rmse * model.rmse, model.sse / 8.0, 1e-12, "RMSE definition");
  }
  {
    bool threw = false;
    try {
      (void)se::fitOls({1.0, 2.0}, {{1.0}, {2.0}, {3.0}}, {"x"}, se::RegressionOptions{});
    } catch (const std::invalid_argument&) {
      threw = true;  // dimension mismatch
    }
    CHECK(threw);
    threw = false;
    try {
      // n = 2 = p (intercept + x): zero residual degrees of freedom make
      // X'X singular, so the fit must throw.
      (void)se::fitOls({1.0, 2.0}, {{1.0}, {2.0}}, {"x"},
                       se::RegressionOptions{});
    } catch (const std::runtime_error&) {
      threw = true;  // n <= parameters
    }
    CHECK(threw);
  }
  {
    const std::string report =
        se::renderRegression(se::fitOls({1.0, 3.0, 2.0, 5.0, 4.0, 6.0},
                                        {{1.0}, {2.0}, {3.0}, {4.0}, {5.0}, {6.0}},
                                        {"x"}, se::RegressionOptions{}),
                             false, false);
    CHECK(report.find("Coefficients") != std::string::npos);
    CHECK(report.find("(intercept)") != std::string::npos);
  }
}

// --------------------------------------------------------------- hypothesis --

void testHypothesis() {
  beginSuite("hypothesis");

  // Distribution functions against textbook critical values.
  NEAR(se::logGamma(0.5), 0.5723649429247001, 1e-12, "logGamma(0.5)");
  NEAR(se::logGamma(5.0), 3.1780538303479458, 1e-12, "logGamma(5)");
  NEAR(se::incompleteBeta(0.5, 0.5, 0.5), 0.5, 1e-12, "I_0.5(0.5,0.5)");
  NEAR(se::incompleteBeta(2.0, 3.0, 0.4), 0.5248, 1e-12, "I_0.4(2,3)");
  NEAR(se::studentTTwoTailedP(2.0, 10.0), 0.0733880347707404, 1e-9, "t p-value");
  NEAR(se::studentTQuantile(0.975, 10.0), 2.2281388519649385, 1e-6, "t quantile");
  NEAR(se::studentTQuantile(0.975, 1.0), 12.706204736432095, 1e-5, "t quantile df=1");
  NEAR(se::studentTCdf(0.0, 5.0), 0.5, 1e-12, "t CDF at 0");
  NEAR(se::chiSquareSf(3.841458820694124, 1.0), 0.05, 1e-6, "chi2 critical df=1");
  NEAR(se::chiSquareSf(5.991464547107979, 2.0), 0.05, 1e-6, "chi2 critical df=2");
  NEAR(se::chiSquareSf(18.307038053275146, 10.0), 0.05, 1e-6, "chi2 critical df=10");
  NEAR(se::normalCdf(1.959963984540054), 0.975, 1e-6, "normal CDF");
  NEAR(se::fSf(4.96460274386841, 1.0, 10.0), 0.05, 1e-6, "F critical (1,10)");
  NEAR(se::fSf(3.3258345301072415, 5.0, 10.0), 0.05, 1e-5, "F critical (5,10)");

  {
    // One-sample t: x = {1..5}, mu = 3 -> t = 0, p = 1, CI symmetric.
    se::TestResult result = se::oneSampleT({1, 2, 3, 4, 5}, 3.0);
    NEAR(result.statistic, 0.0, 1e-12, "t=0");
    NEAR(result.pValue, 1.0, 1e-12, "p=1");
    NEAR(result.meanDiff, 0.0, 1e-12, "mean diff");
    NEAR(result.ciLow, -1.9638, 1e-3, "CI low");
    NEAR(result.ciHigh, 1.9638, 1e-3, "CI high");
    NEAR(result.effectSize, 0.0, 1e-12, "Cohen's d");
  }
  {
    // One-sample against mu = 1: mean diff 2, se = sqrt(2.5)/sqrt(5).
    se::TestResult result = se::oneSampleT({1, 2, 3, 4, 5}, 1.0);
    NEAR(result.statistic, 2.0 / (std::sqrt(2.5) / std::sqrt(5.0)), 1e-12, "t value");
    CHECK(result.pValue < 0.2 && result.pValue > 0.0);
    NEAR(result.df, 4.0, 0.0, "df");
  }
  {
    // Two-sample pooled: a = {1,2,3}, b = {4,5,6} -> t = -3.6742, df = 4.
    se::TestResult pooled = se::twoSampleT({1, 2, 3}, {4, 5, 6}, false);
    NEAR(pooled.statistic, -3.6742346141747673, 1e-9, "pooled t");
    NEAR(pooled.df, 4.0, 1e-12, "pooled df");
    NEAR(pooled.pValue, se::fSf(13.5, 1.0, 4.0), 1e-9, "pooled p (F(1,4) identity)");
    NEAR(pooled.pValue, 0.0213, 1e-4, "pooled p neighbourhood");
    // With equal variances and equal n, Welch df = pooled df = n1+n2-2.
    se::TestResult welch = se::twoSampleT({1, 2, 3}, {4, 5, 6}, true);
    NEAR(welch.df, 4.0, 1e-9, "welch df");
    NEAR(welch.statistic, pooled.statistic, 1e-12, "same t, different df");
    // Unequal variances shrink the Welch df relative to the pooled df,
    // which inflates the two-tailed p-value for the same statistic.
    se::TestResult pooledWide = se::twoSampleT({1, 2, 3}, {4, 5, 9}, false);
    se::TestResult welchWide = se::twoSampleT({1, 2, 3}, {4, 5, 9}, true);
    CHECK(welchWide.df < pooledWide.df);
    CHECK(welchWide.pValue > pooledWide.pValue);
    NEAR(se::cohensD({1, 2, 3}, {4, 5, 6}), -3.0, 1e-9, "Cohen's d = -3");
  }
  {
    // Paired: diffs {2, 2, 3} -> t = 7 on df = 2.
    se::TestResult result = se::pairedT({10, 20, 30}, {8, 18, 27});
    NEAR(result.statistic, 7.0, 1e-9, "paired t");
    NEAR(result.df, 2.0, 1e-12, "paired df");
    CHECK(result.pValue < 0.05);
    CHECK(result.meanDiff > 0.0);
  }
  {
    // Chi-square GOF: uniform -> 0; known off-uniform -> chi2 = 5, df = 1.
    se::TestResult uniform = se::chiSquareGof({10, 10, 10}, {});
    NEAR(uniform.statistic, 0.0, 1e-12, "uniform chi2");
    NEAR(uniform.pValue, 1.0, 1e-12, "uniform p");
    se::TestResult skewed = se::chiSquareGof({15, 5}, {});
    NEAR(skewed.statistic, 5.0, 1e-12, "chi2 = 5");
    NEAR(skewed.df, 1.0, 0.0, "df = 1");
    NEAR(skewed.pValue, se::chiSquareSf(5.0, 1.0), 1e-12, "chi2 p");
    bool threw = false;
    try {
      (void)se::chiSquareGof({10, 10}, {0.4, 0.4});
    } catch (const std::invalid_argument&) {
      threw = true;  // probabilities must sum to 1
    }
    CHECK(threw);
  }
  {
    bool threw = false;
    try {
      (void)se::oneSampleT({1.0}, 0.0);
    } catch (const std::runtime_error&) {
      threw = true;  // needs >= 2 values
    }
    CHECK(threw);
  }
  {
    const se::TestResult result = se::oneSampleT({2.0, 2.2, 1.9, 2.1}, 2.0);
    const std::string report = se::renderTest(result, 0.05, false, false);
    CHECK(report.find("p-value") != std::string::npos);
    CHECK(report.find("=> ") != std::string::npos);
  }
}

// ------------------------------------------------------------- correlation --

void testCorrelation() {
  beginSuite("correlation");

  const std::vector<double> x{1, 2, 3, 4, 5};
  NEAR(se::pearsonCorrelation(x, {2, 4, 6, 8, 10}), 1.0, 1e-12, "perfect positive");
  NEAR(se::pearsonCorrelation(x, {10, 8, 6, 4, 2}), -1.0, 1e-12, "perfect negative");
  CHECK(std::isnan(se::pearsonCorrelation(x, {1, 1, 1, 1, 1})));  // zero variance
  NEAR(se::pearsonCorrelation(x, {1, 3, 2, 5, 4}), 0.8, 1e-12, "known r = 0.8");

  // Spearman is invariant under strictly monotone transforms.
  NEAR(se::spearmanCorrelation(x, {1, 10, 100, 1000, 10000}), 1.0, 1e-12,
       "spearman monotone");
  NEAR(se::spearmanCorrelation(x, {2, 4, 6, 8, 10}), 1.0, 1e-12, "spearman linear");

  // Kendall tau-b on a small untied sample: C = 2, D = 4 -> tau = -1/3.
  NEAR(se::kendallTau({1, 2, 3, 4}, {3, 4, 1, 2}), -1.0 / 3.0, 1e-12, "kendall tau");
  NEAR(se::kendallTau(x, {2, 4, 6, 8, 10}), 1.0, 1e-12, "kendall perfect");

  const std::vector<double> ranks = se::rankAverage({10, 20, 20, 30});
  NEAR(ranks[0], 1.0, 0.0, "rank 1");
  NEAR(ranks[1], 2.5, 0.0, "tied ranks");
  NEAR(ranks[2], 2.5, 0.0, "tied ranks");
  NEAR(ranks[3], 4.0, 0.0, "rank 4");

  NEAR(se::correlationPValue(1.0, 10), 0.0, 0.0, "r=1 -> p=0");
  NEAR(se::correlationPValue(0.9, 5), se::studentTTwoTailedP(0.9 * std::sqrt(3.0 / 0.19), 3),
       1e-12, "p-value formula");
  CHECK(std::isnan(se::correlationPValue(0.5, 2)));

  {
    // Matrix path with a DataFrame, including pairwise-complete n.
    se::DataFrame df = frameFrom("a,b,c\n1,2,5\n2,4,3\n3,6,4\n,8,1\n");
    se::CorrelationMatrix m = se::computeCorrelationMatrix(df, {"a", "b"}, "pearson");
    NEAR(m.r[0][1], 1.0, 1e-12, "matrix r");
    NEAR(m.n[0][1], 3.0, 0.0, "pairwise n");
    CHECK(m.names == std::vector<std::string>({"a", "b"}));
    se::CorrelationMatrix k = se::computeCorrelationMatrix(df, {"a", "c"}, "kendall");
    CHECK(std::isfinite(k.p[0][1]));
    bool threw = false;
    try {
      (void)se::computeCorrelationMatrix(df, {"a", "b"}, "nope");
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    CHECK(threw);
  }
}

// ---------------------------------------------------------------- outliers --

void testOutliers() {
  beginSuite("outliers");

  {
    std::vector<double> data;
    for (int i = 1; i <= 20; ++i) data.push_back(static_cast<double>(i));
    data.push_back(1000.0);
    se::OutlierOptions options;
    se::OutlierReport report = se::detectOutliers(data, options);
    CHECK(report.method == se::OutlierMethod::ZScore);
    CHECK(report.flags.size() == 1 && report.flags[0].value == 1000.0);
    CHECK(report.flags[0].index == 20);
  }
  {
    // Tukey fences on a known sample.
    std::vector<double> data{1, 2, 3, 4, 5, 6, 7, 8, 9, 100};
    se::OutlierOptions options;
    options.method = se::OutlierMethod::Iqr;
    se::OutlierReport report = se::detectOutliers(data, options);
    NEAR(report.center, 5.5, 1e-12, "iqr center");
    NEAR(report.spread, 4.5, 1e-12, "iqr spread");
    NEAR(report.fenceLow, -3.5, 1e-9, "low fence");
    NEAR(report.fenceHigh, 14.5, 1e-9, "high fence");
    CHECK(report.flags.size() == 1);
    NEAR(report.flags[0].score, (100.0 - 14.5) / 4.5, 1e-9, "iqr score");
  }
  {
    // MAD = 0 when more than half the values are identical.
    std::vector<double> data(10, 1.0);
    data.push_back(100.0);
    se::OutlierOptions options;
    options.method = se::OutlierMethod::ModifiedZ;
    se::OutlierReport report = se::detectOutliers(data, options);
    CHECK(report.flags.size() == 1);
    CHECK(report.flags[0].score >= 1e299);
    NEAR(report.center, 1.0, 0.0, "mad center");
    NEAR(report.spread, 0.0, 0.0, "mad spread");
  }
  {
    // Custom threshold for the z method: z(100) ~ 4.37 in this sample.
    std::vector<double> data(20, 0.0);
    data.push_back(100.0);
    se::OutlierOptions options;
    options.threshold = 4.0;
    se::OutlierReport strict = se::detectOutliers(data, options);
    CHECK(strict.flags.size() == 1);
    options.threshold = 5.0;
    se::OutlierReport lax = se::detectOutliers(data, options);
    CHECK(lax.flags.empty());
  }
  {
    // Missing values are excluded from every statistic.
    std::vector<double> data{1.0, se::kNaN, 2.0, se::kNaN, 3.0, 500.0};
    se::OutlierReport report = se::detectOutliers(data, se::OutlierOptions{});
    CHECK(report.n == 4 && report.missing == 2);
  }
  {
    bool threw = false;
    try {
      (void)se::detectOutliers({1.0, 2.0}, se::OutlierOptions{});
    } catch (const std::runtime_error&) {
      threw = true;  // needs >= 3 values
    }
    CHECK(threw);
  }
  {
    const std::vector<double> sorted{1, 2, 3, 4};
    NEAR(se::medianSorted(sorted), 2.5, 0.0, "median even");
    NEAR(se::medianAbsoluteDeviation({1, 1, 2, 2, 2, 3, 3, 3, 3, 9}), 0.5, 1e-12, "MAD");
  }
}

// -------------------------------------------------------------- timeseries --

void testTimeseries() {
  beginSuite("timeseries");

  {
    se::MaOptions options;
    options.method = se::MaMethod::Sma;
    options.window = 3;
    se::MaResult result = se::movingAverage({1, 2, 3, 4, 5}, options);
    NEAR(result.values[0], 1.0, 1e-12, "sma partial window");
    NEAR(result.values[1], 1.5, 1e-12, "sma partial window 2");
    NEAR(result.values[2], 2.0, 1e-12, "sma full window");
    NEAR(result.values[4], 4.0, 1e-12, "sma last");
    CHECK(result.label == "sma(window=3)");
  }
  {
    se::MaOptions options;
    options.method = se::MaMethod::Ema;
    options.span = 2;  // alpha = 2/3
    se::MaResult result = se::movingAverage({1, 2, 3}, options);
    NEAR(result.alpha, 2.0 / 3.0, 1e-12, "ema alpha");
    NEAR(result.values[0], 1.0, 1e-12, "ema seed");
    NEAR(result.values[1], 2.0 / 3.0 * 2.0 + 1.0 / 3.0 * 1.0, 1e-12, "ema step");
    NEAR(result.values[2], 2.0 / 3.0 * 3.0 + 1.0 / 3.0 * (2.0 / 3.0 * 2.0 + 1.0 / 3.0), 1e-12,
         "ema step 2");
  }
  {
    // Missing values carry the previous EMA forward.
    se::MaOptions options;
    options.method = se::MaMethod::Ema;
    options.span = 3;
    se::MaResult result = se::movingAverage({1.0, se::kNaN, 5.0}, options);
    CHECK(!std::isnan(result.values[1]));
    NEAR(result.values[1], result.values[0], 1e-12, "ema carry forward");
  }
  {
    se::MaOptions options;
    options.method = se::MaMethod::Diff;
    options.lag = 2;
    se::MaResult result = se::movingAverage({10, 12, 14, 17}, options);
    CHECK(std::isnan(result.values[0]) && std::isnan(result.values[1]));
    NEAR(result.values[2], 4.0, 1e-12, "diff lag 2");
    NEAR(result.values[3], 5.0, 1e-12, "diff lag 2 last");
  }
  {
    // SMA averages only the valid values inside the window.
    se::MaOptions options;
    options.window = 2;
    se::MaResult result = se::movingAverage({1.0, se::kNaN, 3.0}, options);
    NEAR(result.values[1], 1.0, 1e-12, "sma skips NaN");
    NEAR(result.values[2], 3.0, 1e-12, "sma skips NaN 2");
  }
  {
    se::MaOptions options;
    options.window = 0;
    bool threw = false;
    try {
      (void)se::movingAverage({1.0}, options);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    se::MaOptions options;
    se::MaResult result = se::movingAverage({1, 2, 3, 4}, options);
    const std::string table = se::renderSeries(result, {1, 2, 3, 4}, "v", 40, false, false);
    CHECK(table.find("sma(window=5)") != std::string::npos);
  }
}

// --------------------------------------------------------------- histogram --

void testHistogram() {
  beginSuite("histogram");

  {
    std::vector<double> data;
    for (int i = 0; i < 10; ++i) data.push_back(static_cast<double>(i));
    se::HistogramOptions options;
    options.rule = se::BinningRule::Fixed;
    options.bins = 5;
    se::HistogramResult result = se::computeHistogram(data, options);
    CHECK(result.bins.size() == 5);
    NEAR(result.binWidth, 1.8, 1e-12, "bin width");
    long long sum = 0;
    for (const se::Bin& bin : result.bins) sum += bin.count;
    CHECK(sum == 10);
    for (std::size_t i = 0; i < result.bins.size(); ++i) {
      CHECK(result.bins[i].count == 2);  // 0/1, 2/3, 4/5, 6/7, 8/9
    }
    NEAR(result.bins.back().hi, 9.0, 1e-12, "last bin closed at max");
    CHECK(result.usedRuleName == "fixed");
    NEAR(result.median, 4.5, 1e-12, "histogram median");
  }
  {
    // Constant data -> a single bin.
    se::HistogramResult result =
        se::computeHistogram({7.0, 7.0, 7.0}, se::HistogramOptions{});
    CHECK(result.bins.size() == 1 && result.bins[0].count == 3);
    CHECK(result.usedRuleName == "single-value");
  }
  {
    // Sturges: k = ceil(log2 n) + 1.
    std::vector<double> data;
    for (int i = 0; i < 32; ++i) data.push_back(i * 1.0);
    se::HistogramOptions options;
    options.rule = se::BinningRule::Sturges;
    se::HistogramResult result = se::computeHistogram(data, options);
    CHECK(result.bins.size() == 6 && result.usedRuleName == "sturges");
  }
  {
    // Freedman-Diaconis falls back to Sturges when IQR == 0 while the range
    // stays non-degenerate (30 identical values plus two extremes).
    std::vector<double> data(30, 5.0);
    data.push_back(0.0);
    data.push_back(100.0);
    se::HistogramOptions options;
    options.rule = se::BinningRule::FreedmanDiaconis;
    se::HistogramResult result = se::computeHistogram(data, options);
    CHECK(result.usedRuleName == "sturges");
  }
  {
    // Bin width option overrides the rule.
    std::vector<double> data{0, 0.5, 1.0, 1.5, 2.0};
    se::HistogramOptions options;
    options.binWidth = 1.0;
    se::HistogramResult result = se::computeHistogram(data, options);
    CHECK(result.bins.size() == 2);
    NEAR(result.binWidth, 1.0, 1e-12, "explicit width");
  }
  {
    bool threw = false;
    try {
      se::HistogramOptions options;
      options.rule = se::BinningRule::Fixed;
      options.bins = 0;
      (void)se::computeHistogram({1.0}, options);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    CHECK(threw);
    threw = false;
    try {
      (void)se::computeHistogram({se::kNaN}, se::HistogramOptions{});
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    // Rendering must not depend on colour state and includes every label.
    se::HistogramResult result = se::computeHistogram({1.0, 2.0, 3.0},
                                                      se::HistogramOptions{});
    const std::string text = se::renderHistogram(result, "v", false, false);
    CHECK(text.find("Histogram of 'v'") != std::string::npos);
    CHECK(text.find("100.0%") != std::string::npos);
  }
}

// ------------------------------------------------------------------- anova --

void testAnova() {
  beginSuite("anova");

  {
    // Two groups of three: F = t_pooled^2 = 13.5, eta^2 = 13.5/17.5.
    std::vector<std::pair<std::string, std::vector<double>>> groups{
        {"a", {1, 2, 3}}, {"b", {4, 5, 6}}};
    se::AnovaResult result = se::oneWayAnova(groups);
    CHECK(result.valid && result.k == 2 && result.n == 6);
    NEAR(result.grandMean, 3.5, 1e-12, "grand mean");
    NEAR(result.ssBetween, 13.5, 1e-12, "SS between");
    NEAR(result.ssWithin, 4.0, 1e-12, "SS within");
    NEAR(result.ssTotal, 17.5, 1e-12, "SS total");
    NEAR(result.fStat, 13.5, 1e-12, "F = t^2");
    NEAR(result.dfBetween, 1.0, 0.0, "df1");
    NEAR(result.dfWithin, 4.0, 0.0, "df2");
    NEAR(result.pValue, se::fSf(13.5, 1.0, 4.0), 1e-12, "F p-value = pooled t p-value");
    NEAR(result.pValue, 0.0213, 1e-4, "p in the right neighbourhood");
    NEAR(result.etaSquared, 13.5 / 17.5, 1e-12, "eta squared");
    NEAR(result.msWithin, 1.0, 1e-12, "pooled variance");
  }
  {
    // Three equal groups: everything cancels to F = 0.
    std::vector<std::pair<std::string, std::vector<double>>> groups{
        {"a", {1, 2, 3}}, {"b", {1, 2, 3}}, {"c", {1, 2, 3}}};
    se::AnovaResult result = se::oneWayAnova(groups);
    NEAR(result.fStat, 0.0, 1e-12, "identical groups F = 0");
    NEAR(result.pValue, 1.0, 1e-12, "identical groups p = 1");
    NEAR(result.etaSquared, 0.0, 1e-12, "eta^2 = 0");
  }
  {
    // Non-finite observations are dropped before the decomposition.
    std::vector<std::pair<std::string, std::vector<double>>> groups{
        {"a", {1.0, se::kNaN, 3.0}}, {"b", {4.0, 6.0}}};
    se::AnovaResult result = se::oneWayAnova(groups);
    CHECK(result.n == 4 && result.groups[0].n == 2);
  }
  {
    bool threw = false;
    try {
      (void)se::oneWayAnova({{"a", {1, 2}}, {"a", {3, 4}}});
    } catch (const std::invalid_argument&) {
      threw = true;  // duplicate labels
    }
    CHECK(threw);
    threw = false;
    try {
      (void)se::oneWayAnova({{"a", {1, 2}}});
    } catch (const std::runtime_error&) {
      threw = true;  // single group
    }
    CHECK(threw);
    threw = false;
    try {
      (void)se::oneWayAnova({{"a", {1}}, {"b", {2}}});
    } catch (const std::runtime_error&) {
      threw = true;  // groups need >= 2 values
    }
    CHECK(threw);
  }
  {
    NEAR(se::fSf(1.0, 1.0, 1.0), 0.5, 1e-9, "F(1,1) at 1");
    CHECK(se::fSf(0.0, 3.0, 10.0) == 1.0);
    CHECK(se::fSf(se::kNaN, 3.0, 10.0) != se::fSf(se::kNaN, 3.0, 10.0));
    bool threw = false;
    try {
      (void)se::fSf(1.0, 0.0, 10.0);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    std::vector<std::pair<std::string, std::vector<double>>> groups{
        {"a", {1, 2, 3}}, {"b", {4, 5, 6}}};
    const se::AnovaResult result = se::oneWayAnova(groups);
    const std::string report = se::renderAnova(result, 0.05, false, false);
    CHECK(report.find("ANOVA table") != std::string::npos);
    CHECK(report.find("reject H0") != std::string::npos);
  }
}

// -------------------------------------------------------------------- json --

void testJson() {
  beginSuite("json");

  std::ostringstream out;
  se::json::Writer writer(out);
  writer.beginObject();
  writer.kv("n", 42LL);
  writer.kv("pi", 3.5);
  writer.kv("s", "quote\" back\\slash\n");
  writer.kvBool("yes", true);
  writer.key("list");
  writer.beginArray();
  writer.value(1.0);
  writer.value(se::kNaN);  // non-finite -> null
  writer.null();
  writer.endArray();
  writer.key("nested");
  writer.beginObject();
  writer.kv("deep", std::string("x"));
  writer.endObject();
  writer.endObject();

  const std::string text = out.str();
  CHECK(text ==
        "{\"n\":42,\"pi\":3.5,\"s\":\"quote\\\" back\\\\slash\\n\","
        "\"yes\":true,\"list\":[1,null,null],\"nested\":{\"deep\":\"x\"}}");

  std::ostringstream empty;
  se::json::Writer w2(empty);
  w2.beginArray();
  w2.endArray();
  CHECK(empty.str() == "[]");

  std::ostringstream ctrl;
  se::json::Writer w3(ctrl);
  w3.value(std::string("\x01"));
  CHECK(ctrl.str() == "\"\\u0001\"");
}

// -------------------------------------------------------------------- args --

void testArgs() {
  beginSuite("args");

  {
    char arg0[] = "statsengine";
    char arg1[] = "describe";
    char arg2[] = "--file";
    char arg3[] = "x.csv";
    char arg4[] = "--json";
    char arg5[] = "--mu";
    char arg6[] = "-1.5";
    char* argv[] = {arg0, arg1, arg2, arg3, arg4, arg5, arg6};
    se::Args args(7, argv);
    CHECK(args.command() == "describe");
    CHECK(args.get("file") == "x.csv");
    CHECK(args.jsonMode());
    NEAR(args.getDouble("mu", 0.0), -1.5, 1e-12, "negative option value");
    CHECK(args.has("mu") && !args.has("nope"));
  }
  {
    char arg0[] = "prog";
    char arg1[] = "--flag=value";
    char arg2[] = "--empty=";
    char* argv[] = {arg0, arg1, arg2};
    se::Args args(3, argv);
    CHECK(args.get("flag") == "value");
    CHECK(args.has("empty") && args.get("empty").empty());
  }
  {
    char arg0[] = "prog";
    char arg1[] = "--";
    char arg2[] = "--not-a-flag";
    char* argv[] = {arg0, arg1, arg2};
    se::Args args(3, argv);
    CHECK(args.positionals().size() == 1 && args.positionals()[0] == "--not-a-flag");
  }
  {
    char arg0[] = "prog";
    char arg1[] = "clean";
    char arg2[] = "--bogus";
    char* argv[] = {arg0, arg1, arg2};
    se::Args args(3, argv);
    bool threw = false;
    try {
      args.rejectUnknown({"file"});
    } catch (const se::ArgsError&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    char arg0[] = "prog";
    char arg1[] = "-x";
    char* argv[] = {arg0, arg1};
    bool threw = false;
    try {
      se::Args args(2, argv);
      (void)args;
    } catch (const se::ArgsError&) {
      threw = true;  // single-dash options are rejected
    }
    CHECK(threw);
  }
  {
    char arg0[] = "prog";
    char arg1[] = "cmd";
    char arg2[] = "--bins";
    char arg3[] = "9";
    char* argv[] = {arg0, arg1, arg2, arg3};
    se::Args args(4, argv);
    NEAR(static_cast<double>(args.getInt("bins", 0)), 9.0, 0.0, "getInt");
    args.expectNoPositionals(0);  // only the command itself
    bool threw = false;
    try {
      args.require("missing", "thing");
    } catch (const se::ArgsError&) {
      threw = true;
    }
    CHECK(threw);
  }
}

// ------------------------------------------------------------------ table --

void testTableAndColors() {
  beginSuite("table");

  se::Table table;
  table.setTitle("T");
  table.addColumn("a", se::Align::Left);
  table.addColumn("b", se::Align::Right);
  table.addRow({"xx", "1"});
  table.addRow({"y", "22"});
  table.setFooter("done");
  const std::string ascii = table.render(false);
  CHECK(ascii.find("+----+----+") != std::string::npos);
  CHECK(ascii.find("| y  | 22 |") != std::string::npos);
  CHECK(ascii.find("| xx |  1 |") != std::string::npos);
  CHECK(ascii.find("done") != std::string::npos);
  const std::string unicode = table.render(true);
  CHECK(unicode.find("\xe2\x94\x82") != std::string::npos);  // │

  // ANSI escapes never break the visible width.
  const std::string painted = "\033[31mabc\033[0m";
  CHECK(se::colors::stripAnsi(painted) == "abc");
  CHECK(se::colors::visibleWidth(painted) == 3);
  CHECK(se::colors::paint("x", "1") == "x");  // colours disabled by default
  CHECK(se::colors::paint("x", nullptr) == "x");

  const std::string block =
      se::keyValueBlock({{"long_key", "v"}, {"k", "v"}}, false);
  CHECK(block.find("long_key : v") != std::string::npos);
  CHECK(block.find("k        : v") != std::string::npos);  // padded to 8, then " : v"

  const std::string bar = se::progressBar(0.5, 4, false, false);
  CHECK(bar == "##..");
  CHECK(se::progressBar(1.0, 3, false, false) == "###");
}

}  // namespace

int main() {
  testUtil();
  testCsv();
  testDataFrame();
  testDescribe();
  testRegression();
  testHypothesis();
  testCorrelation();
  testOutliers();
  testTimeseries();
  testHistogram();
  testAnova();
  testJson();
  testArgs();
  testTableAndColors();

  std::printf("\n%d checks, %d failures\n", gChecks, gFailures);
  if (gFailures == 0) {
    std::printf("ALL TESTS PASSED\n");
    return 0;
  }
  std::printf("TESTS FAILED\n");
  return 1;
}
