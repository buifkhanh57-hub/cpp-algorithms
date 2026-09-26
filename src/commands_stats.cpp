// statsengine — commands_stats.cpp
//
// The corr / regress / ma / test / anova subcommands.
#include "statsengine/commands.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/anova.hpp"
#include "statsengine/correlation.hpp"
#include "statsengine/hypothesis.hpp"
#include "statsengine/json.hpp"
#include "statsengine/regression.hpp"
#include "statsengine/table.hpp"
#include "statsengine/timeseries.hpp"
#include "statsengine/util.hpp"

namespace se {

namespace {

// Group a numeric column by the raw values of a grouping column; group
// order follows first appearance; missing values are skipped.
std::vector<std::pair<std::string, std::vector<double>>> groupByRaw(const DataFrame& df,
                                                                    const std::string& valueColumn,
                                                                    const std::string& groupColumn) {
  const Column& values = requireNumericColumn(df, valueColumn);
  const Column& groups = df.column(groupColumn);
  std::map<std::string, std::size_t> order;
  std::vector<std::pair<std::string, std::vector<double>>> out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values.isMissing(i) || groups.isMissing(i)) continue;
    const std::string key = trim(groups.rawAt(i));
    const auto it = order.find(key);
    if (it == order.end()) {
      order[key] = out.size();
      out.emplace_back(key, std::vector<double>());
      out.back().second.push_back(values.value(i));
    } else {
      out[it->second].second.push_back(values.value(i));
    }
  }
  return out;
}

}  // namespace

// ----------------------------------------------------------------- corr --

int cmdCorrelation(const Args& args) {
  args.rejectUnknown({"file", "column", "columns", "method"});
  args.expectNoPositionals(0);
  const DataFrame df = loadFrame(args);

  const std::string method = args.getOr("method", "pearson");
  std::vector<std::string> names = selectedColumnNames(args);
  std::vector<std::string> resolved;
  if (names.empty()) {
    for (std::size_t idx : df.numericColumnIndices()) resolved.push_back(df.column(idx).name());
  } else {
    for (const std::string& name : names) {
      const Column& col = df.column(name);
      if (col.kind() != ColumnKind::Numeric) {
        throw ColumnError(name, "correlation requires numeric columns but '" + name + "' is " +
                                    columnKindName(col.kind()));
      }
      resolved.push_back(col.name());
    }
  }
  if (resolved.size() < 2) {
    throw std::runtime_error("correlation needs at least two numeric columns");
  }

  const CorrelationMatrix matrix = computeCorrelationMatrix(df, resolved, method);
  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    correlationToJson(writer, matrix);
    std::cout << "\n";
    return 0;
  }
  std::cout << renderCorrelationMatrix(matrix, !args.asciiMode(), colors::enabled()) << "\n";
  return 0;
}

// -------------------------------------------------------------- regress --

int cmdRegression(const Args& args) {
  args.rejectUnknown({"file", "target", "column", "columns", "no-intercept", "predict",
                      "predict-file"});
  args.expectNoPositionals(0);
  args.require("target", "numeric response column");

  const DataFrame df = loadFrame(args);
  const std::string targetName = args.get("target");
  const Column& target = requireNumericColumn(df, targetName);

  std::vector<std::string> predictorNames = selectedColumnNames(args);
  if (predictorNames.empty()) {
    for (std::size_t idx : df.numericColumnIndices()) {
      const std::string& name = df.column(idx).name();
      if (name != target.name()) predictorNames.push_back(name);
    }
  } else {
    for (const std::string& name : predictorNames) {
      requireNumericColumn(df, name);  // validates existence + kind
      if (name == target.name()) {
        throw std::invalid_argument("the target column '" + name +
                                    "' cannot also be a predictor");
      }
    }
  }
  if (predictorNames.empty()) {
    throw std::runtime_error("no predictor columns found; use --columns <a,b,c>");
  }

  RegressionOptions options;
  options.hasIntercept = !args.has("no-intercept");
  const double alphaOption = args.getDouble("alpha", 0.05);
  if (!(alphaOption > 0.0) || alphaOption >= 1.0) {
    throw ArgsError("--alpha must be in (0, 1)");
  }
  options.alpha = alphaOption;

  std::vector<std::string> completeNames;
  completeNames.reserve(predictorNames.size() + 1);
  completeNames.push_back(target.name());
  completeNames.insert(completeNames.end(), predictorNames.begin(), predictorNames.end());
  const std::vector<std::size_t> rows = df.completeRows(completeNames);
  std::vector<double> y;
  std::vector<std::vector<double>> x;
  y.reserve(rows.size());
  x.reserve(rows.size());
  for (std::size_t row : rows) {
    y.push_back(target.value(row));
    std::vector<double> rowValues;
    rowValues.reserve(predictorNames.size());
    for (const std::string& name : predictorNames) {
      rowValues.push_back(df.column(name).value(row));
    }
    x.push_back(std::move(rowValues));
  }
  if (rows.size() < predictorNames.size() + 2) {
    throw std::runtime_error("only " + std::to_string(rows.size()) +
                             " complete rows available for " +
                             std::to_string(predictorNames.size() + 1) + " parameters");
  }

  RegressionResult model = fitOls(y, x, predictorNames, options);
  model.targetName = target.name();

  // Optional predictions.
  std::vector<std::vector<double>> predictionRows;
  std::vector<std::string> predictionLabels;
  if (args.has("predict")) {
    std::vector<std::string> tokens = splitList(args.get("predict"));
    if (tokens.size() != predictorNames.size()) {
      throw ArgsError("--predict expects " + std::to_string(predictorNames.size()) +
                      " comma-separated values, got " + std::to_string(tokens.size()));
    }
    std::vector<double> row;
    for (const std::string& token : tokens) {
      double value = 0.0;
      if (!parseNumberStrict(token, value)) {
        throw ArgsError("--predict expects numbers, got '" + token + "'");
      }
      row.push_back(value);
    }
    predictionRows.push_back(row);
    predictionLabels.push_back("manual");
  }
  if (args.has("predict-file")) {
    const std::string path = args.get("predict-file");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open file '" + path + "'");
    CsvReader reader(in);
    CsvDocument doc = reader.read();
    std::size_t firstDataRow = 0;
    if (!doc.rows.empty()) {
      bool headerLike = false;
      for (const std::string& cell : doc.rows[0]) {
        double ignored = 0.0;
        if (!parseNumberStrict(cell, ignored)) headerLike = true;
      }
      if (headerLike) {
        const std::vector<std::string> header = doc.rows[0];
        if (header.size() != predictorNames.size()) {
          throw ArgsError("predict file has " + std::to_string(header.size()) +
                          " columns, expected " + std::to_string(predictorNames.size()));
        }
        for (std::size_t j = 0; j < header.size(); ++j) {
          if (toLower(trim(header[j])) != toLower(predictorNames[j])) {
            throw ArgsError("predict file column '" + header[j] + "' does not match predictor '" +
                            predictorNames[j] + "' (order must match)");
          }
        }
        firstDataRow = 1;
      }
    }
    for (std::size_t i = firstDataRow; i < doc.rows.size(); ++i) {
      const std::vector<std::string>& raw = doc.rows[i];
      if (raw.size() != predictorNames.size()) {
        throw ArgsError("predict file row " + std::to_string(i + 1) + " has " +
                        std::to_string(raw.size()) + " fields, expected " +
                        std::to_string(predictorNames.size()));
      }
      std::vector<double> row;
      for (const std::string& cell : raw) {
        double value = 0.0;
        if (!parseNumberStrict(cell, value)) {
          throw ArgsError("predict file row " + std::to_string(i + 1) + " has invalid value '" +
                          cell + "'");
        }
        row.push_back(value);
      }
      predictionRows.push_back(row);
      predictionLabels.push_back("row " + std::to_string(i - firstDataRow + 1));
    }
  }

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    regressionToJson(writer, model);
    if (!predictionRows.empty()) {
      writer.key("predictions");
      writer.beginArray();
      std::vector<double> values = predict(model, predictionRows);
      for (std::size_t i = 0; i < values.size(); ++i) {
        writer.beginObject();
        writer.kv("source", predictionLabels[i]);
        writer.kv("prediction", values[i]);
        writer.endObject();
      }
      writer.endArray();
    }
    writer.endObject();
    std::cout << "\n";
    return 0;
  }

  std::cout << renderRegression(model, !args.asciiMode(), colors::enabled()) << "\n";
  if (!predictionRows.empty()) {
    const std::vector<double> values = predict(model, predictionRows);
    Table table;
    table.setTitle("Predictions");
    table.addColumn("source", Align::Left);
    table.addColumn("prediction", Align::Right);
    for (std::size_t i = 0; i < values.size(); ++i) {
      table.addRow({predictionLabels[i], fmtFixed(values[i], 6)});
    }
    std::cout << table.render(!args.asciiMode()) << "\n";
  }
  return 0;
}

// ------------------------------------------------------------------- ma --

int cmdMovingAverage(const Args& args) {
  args.rejectUnknown({"file", "column", "method", "window", "span", "alpha", "lag"});
  args.expectNoPositionals(0);
  args.require("column", "numeric column to transform");

  MaOptions options;
  const std::string method = args.getOr("method", "sma");
  if (!maMethodFromName(method, options.method)) {
    throw ArgsError("unknown ma method '" + method + "' (expected sma, ema or diff)");
  }
  options.window = static_cast<int>(args.getInt("window", 5));
  options.span = static_cast<int>(args.getInt("span", 0));
  options.lag = static_cast<int>(args.getInt("lag", 1));
  options.alpha = args.getDouble("alpha", 0.0);

  const DataFrame df = loadFrame(args);
  const std::string columnName = args.get("column");
  const Column& col = requireNumericColumn(df, columnName);
  const MaResult result = movingAverage(col.numeric(), options);

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    seriesToJson(writer, result, col.name());
    std::cout << "\n";
    return 0;
  }
  std::cout << renderSeries(result, col.numeric(), col.name(), 40, !args.asciiMode(),
                            colors::enabled())
            << "\n";
  return 0;
}

// ----------------------------------------------------------------- test --

int cmdTest(const Args& args) {
  args.rejectUnknown({"file", "test", "column", "with", "by", "a", "b", "mu", "alpha", "expected",
                      "no-welch"});
  args.expectNoPositionals(0);
  args.require("column", "column to test");

  const std::string kind = toLower(args.getOr("test", "t1"));
  const double alpha = args.getDouble("alpha", 0.05);
  if (!(alpha > 0.0) || alpha >= 1.0) {
    throw ArgsError("--alpha must be in (0, 1)");
  }

  const DataFrame df = loadFrame(args);
  TestResult result;

  if (kind == "t1") {
    const std::string columnName = args.get("column");
    const Column& col = requireNumericColumn(df, columnName);
    const double mu = args.getDouble("mu", 0.0);
    result = oneSampleT(col.numericCompact(), mu);
  } else if (kind == "t2") {
    const std::string group = args.get("by");
    if (group.empty()) {
      throw ArgsError("two-sample t-test needs --by <group-column>");
    }
    if (!df.hasColumn(group)) {
      throw ColumnError(group, "unknown grouping column '" + group + "'");
    }
    std::vector<std::pair<std::string, std::vector<double>>> groups =
        groupByRaw(df, args.get("column"), group);
    groups.erase(std::remove_if(groups.begin(), groups.end(),
                                [](const std::pair<std::string, std::vector<double>>& g) {
                                  return g.second.size() < 2;
                                }),
                 groups.end());
    if (groups.size() < 2) {
      throw std::runtime_error("--by '" + group +
                               "' does not produce two groups with at least 2 values each");
    }
    if (groups.size() > 2 && !(args.has("a") && args.has("b"))) {
      std::vector<std::string> levelNames;
      for (const auto& g : groups) levelNames.push_back(g.first);
      throw ArgsError("grouping column has more than two levels (" + joinList(levelNames) +
                      "); pick two with --a <G1> --b <G2>");
    }
    const std::string nameA = args.has("a") ? args.get("a") : groups[0].first;
    const std::string nameB = args.has("b") ? args.get("b") : groups[1].first;
    const std::vector<double>* valuesA = nullptr;
    const std::vector<double>* valuesB = nullptr;
    for (const auto& g : groups) {
      if (g.first == nameA) valuesA = &g.second;
      if (g.first == nameB) valuesB = &g.second;
    }
    if (valuesA == nullptr || valuesB == nullptr || nameA == nameB) {
      throw ArgsError("groups '" + nameA + "' / '" + nameB + "' not found in column '" + group +
                      "'");
    }
    result = twoSampleT(*valuesA, *valuesB, !args.has("no-welch"));
    result.details.insert(result.details.begin(),
                          "groups: " + nameA + " vs " + nameB + " (from '" + group + "')");
  } else if (kind == "paired") {
    const std::string second = args.get("with");
    if (second.empty()) {
      throw ArgsError("paired t-test needs --with <second-column>");
    }
    const std::string firstName = args.get("column");
    const Column& a = requireNumericColumn(df, firstName);
    const Column& b = requireNumericColumn(df, second);
    result = pairedT(a.numeric(), b.numeric());
    result.details.insert(result.details.begin(),
                          "columns: " + a.name() + " vs " + b.name());
  } else if (kind == "chisq") {
    const std::string columnName = args.get("column");
    const Column& col = df.column(columnName);
    const std::vector<std::string> categories = col.categories();
    if (categories.size() < 2) {
      throw std::runtime_error("chi-square test needs at least 2 distinct categories in '" +
                               col.name() + "'");
    }
    std::vector<long long> observed;
    observed.reserve(categories.size());
    std::map<std::string, long long> counts;
    for (std::size_t i = 0; i < col.size(); ++i) {
      if (!col.isMissing(i)) ++counts[trim(col.rawAt(i))];
    }
    for (const std::string& category : categories) observed.push_back(counts[category]);

    std::vector<double> expected;
    const std::vector<std::string> expectedTokens = args.getList("expected");
    if (!expectedTokens.empty()) {
      for (const std::string& token : expectedTokens) {
        double p = 0.0;
        if (!parseNumberStrict(token, p)) {
          throw ArgsError("--expected expects numbers, got '" + token + "'");
        }
        expected.push_back(p);
      }
    }
    result = chiSquareGof(observed, expected);
    if (!expected.empty() && expected.size() != categories.size()) {
      // chiSquareGof already validates; keep symmetric behavior here too.
      throw ArgsError("--expected count does not match the number of categories (" +
                      std::to_string(categories.size()) + ")");
    }
    std::string categoryList;
    for (std::size_t i = 0; i < categories.size(); ++i) {
      if (i > 0) categoryList += ", ";
      categoryList += categories[i] + "=" + std::to_string(observed[i]);
    }
    result.details.insert(result.details.begin(), "observed: " + categoryList);
  } else {
    throw ArgsError("unknown test '" + args.get("test") +
                    "' (expected t1, t2, paired or chisq)");
  }

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    testToJson(writer, result, alpha);
    std::cout << "\n";
    return 0;
  }
  std::cout << renderTest(result, alpha, !args.asciiMode(), colors::enabled());
  return 0;
}

// ----------------------------------------------------------------- anova --

const char* kAnovaUsage =
    "Usage:\n"
    "  statsengine anova --file <csv> --column <numeric> --by <group> [--only levels]\n"
    "                    [--alpha A]\n"
    "\n"
    "One-way analysis of variance: is the mean of the numeric column the same\n"
    "across every level of the grouping column? Reports the F statistic, the\n"
    "between/within decomposition, eta-squared and omega-squared.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --column <name>    numeric response column (required)\n"
    "  --by <name>        grouping column (all its levels are used)\n"
    "  --only <list>      restrict the analysis to a comma-separated subset of\n"
    "                     levels, e.g. --only north,south\n"
    "  --alpha <A>        significance level (default 0.05)\n"
    "  --json             machine-readable JSON output\n";

int cmdAnova(const Args& args) {
  args.rejectUnknown({"file", "column", "by", "only", "alpha"});
  args.expectNoPositionals(0);
  args.require("column", "numeric response column");
  const std::string group = args.get("by");
  if (group.empty()) {
    throw ArgsError("anova needs --by <group-column>");
  }
  const double alpha = args.getDouble("alpha", 0.05);
  if (!(alpha > 0.0) || alpha >= 1.0) {
    throw ArgsError("--alpha must be in (0, 1)");
  }

  const DataFrame df = loadFrame(args);
  if (!df.hasColumn(group)) {
    throw ColumnError(group, "unknown grouping column '" + group + "'");
  }
  const std::vector<std::pair<std::string, std::vector<double>>> allGroups =
      groupByRaw(df, args.get("column"), group);

  const std::vector<std::string> only = splitList(args.get("only"));
  std::vector<std::pair<std::string, std::vector<double>>> selected;
  if (only.empty()) {
    selected = allGroups;
  } else {
    std::set<std::string> wanted(only.begin(), only.end());
    for (const auto& entry : allGroups) {
      if (wanted.count(entry.first) > 0) {
        selected.push_back(entry);
        wanted.erase(entry.first);
      }
    }
    if (!wanted.empty()) {
      std::vector<std::string> missing(wanted.begin(), wanted.end());
      std::vector<std::string> available;
      for (const auto& entry : allGroups) available.push_back(entry.first);
      throw ArgsError("--only levels not found in '" + group + "': " + joinList(missing) +
                      " (available: " + joinList(available) + ")");
    }
  }
  if (selected.size() < 2) {
    throw std::runtime_error("anova needs at least 2 groups with observations in '" +
                             group + "'");
  }

  const AnovaResult result = oneWayAnova(selected);
  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    anovaToJson(writer, result, alpha, args.get("column"), group);
    std::cout << "\n";
    return 0;
  }
  std::cout << renderAnova(result, alpha, !args.asciiMode(), colors::enabled()) << "\n";
  return 0;
}

}  // namespace se
