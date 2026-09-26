// statsengine — commands_data.cpp
//
// Shared CLI helpers plus the describe / hist / outliers / clean / sample
// subcommands and the command registry.
#include "statsengine/commands.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/describe.hpp"
#include "statsengine/histogram.hpp"
#include "statsengine/json.hpp"
#include "statsengine/outliers.hpp"
#include "statsengine/table.hpp"
#include "statsengine/util.hpp"

namespace se {

// ------------------------------------------------------- shared helpers --

DataFrame loadFrame(const Args& args) {
  args.require("file", "path to a CSV file ('-' reads stdin)");
  const std::string path = args.get("file");

  if (path == "-") {
    CsvReader reader(std::cin);
    return DataFrame::fromCsv(reader.read(), CsvOptions{});
  }
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open file '" + path + "'");
  }
  CsvReader reader(in);
  return DataFrame::fromCsv(reader.read(), CsvOptions{});
}

const Column& requireColumn(const DataFrame& df, const std::string& name) {
  return df.column(name);
}

const Column& requireNumericColumn(const DataFrame& df, const std::string& name) {
  const Column& col = df.column(name);
  if (col.kind() != ColumnKind::Numeric) {
    throw ColumnError(name, "expected a numeric column but '" + name + "' is " +
                                columnKindName(col.kind()));
  }
  return col;
}

std::vector<std::string> selectedColumnNames(const Args& args) {
  std::vector<std::string> names = args.getList("column");
  std::vector<std::string> extra = args.getList("columns");
  names.insert(names.end(), extra.begin(), extra.end());
  std::vector<std::string> unique;
  std::set<std::string> seen;
  for (const std::string& name : names) {
    if (seen.insert(toLower(name)).second) unique.push_back(name);
  }
  return unique;
}

// ------------------------------------------------------------- registry --

namespace {

const char* kDescribeUsage =
    "Usage:\n"
    "  statsengine describe --file <csv> [--column <name> | --columns <a,b,c> | --all]\n"
    "\n"
    "Summarise each selected numeric column: count, mean, sample stddev (ddof=1),\n"
    "min, quartiles (linear interpolation), max, IQR, skewness and excess kurtosis\n"
    "(bias-corrected G1/G2). Non-numeric columns are rejected; with no selection\n"
    "every numeric column is described.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --column <name>    column to describe (comma list accepted)\n"
    "  --columns <list>   comma-separated column list\n"
    "  --all              describe every numeric column (default)\n"
    "  --json             machine-readable JSON output\n";

const char* kHistogramUsage =
    "Usage:\n"
    "  statsengine hist --file <csv> --column <name> [--rule auto|sturges|fd|scott|fixed]\n"
    "                   [--bins N] [--bin-width W] [--bar-width N] [--cumulative]\n"
    "\n"
    "Binned histogram of a numeric column rendered as an ASCII/Unicode bar chart.\n"
    "Default rule 'auto' uses Freedman-Diaconis with a Sturges fallback.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --column <name>    numeric column to bin (required)\n"
    "  --rule <name>      auto | sturges | fd | scott | fixed\n"
    "  --bins <N>         bin count when --rule fixed\n"
    "  --bin-width <W>    fixed bin width (overrides the rule)\n"
    "  --bar-width <N>    maximum bar length in characters (default 46)\n"
    "  --cumulative       add a cumulative percentage column\n"
    "  --json             machine-readable JSON output\n";

const char* kCorrelationUsage =
    "Usage:\n"
    "  statsengine corr --file <csv> [--columns <a,b,c>] [--method pearson|spearman|kendall]\n"
    "\n"
    "Correlation matrix over numeric columns with pairwise-complete observations,\n"
    "p-values and significance stars.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --columns <list>   columns to correlate (default: all numeric)\n"
    "  --method <name>    pearson (default) | spearman | kendall\n"
    "  --json             machine-readable JSON output\n";

const char* kRegressionUsage =
    "Usage:\n"
    "  statsengine regress --file <csv> --target <y> [--columns <x1,x2,...> | --column <x>]\n"
    "                      [--no-intercept] [--predict \"v1,v2,...\"]\n"
    "                      [--predict-file <csv>]\n"
    "\n"
    "Ordinary least squares via normal equations (Gaussian elimination with\n"
    "pivoting): coefficients with standard errors, t stats, p-values and 95% CIs,\n"
    "R-squared, adjusted R-squared, F test, RMSE and Durbin-Watson.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --target <name>    response variable (numeric, required)\n"
    "  --columns <list>   predictor columns (default: all other numeric)\n"
    "  --no-intercept     suppress the intercept term\n"
    "  --predict <vals>   predict one new row, e.g. --predict \"12.5,3\"\n"
    "  --predict-file <p> predict every row of a CSV (header optional)\n"
    "  --json             machine-readable JSON output\n";

const char* kMovingAverageUsage =
    "Usage:\n"
    "  statsengine ma --file <csv> --column <name> [--method sma|ema|diff]\n"
    "                 [--window N] [--span N] [--alpha A] [--lag N]\n"
    "\n"
    "Smoothing/transform of an ordered column: trailing simple moving average,\n"
    "exponential moving average (seeded at the first value) or lagged differencing.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --column <name>    numeric column to transform (required)\n"
    "  --method <name>    sma (default) | ema | diff\n"
    "  --window <N>       SMA window (default 5); EMA fallback span\n"
    "  --span <N>         EMA span, alpha = 2/(span+1)\n"
    "  --alpha <A>        explicit EMA smoothing factor in (0, 1]\n"
    "  --lag <N>          diff lag (default 1)\n"
    "  --json             machine-readable JSON output\n";

const char* kTestUsage =
    "Usage:\n"
    "  statsengine test --file <csv> --test t1 --column <x> [--mu M]\n"
    "  statsengine test --file <csv> --test t2 --column <x> --by <group> [--a G1 --b G2]\n"
    "  statsengine test --file <csv> --test paired --column <a> --with <b>\n"
    "  statsengine test --file <csv> --test chisq --column <cat> [--expected p1,p2,...]\n"
    "\n"
    "One-sample / two-sample (Welch or pooled) / paired t-tests and chi-square\n"
    "goodness-of-fit. P-values come from the regularized incomplete beta function.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --test <kind>      t1 (default) | t2 | paired | chisq\n"
    "  --column <name>    tested column (numeric except for chisq)\n"
    "  --mu <M>           H0 mean for t1 (default 0)\n"
    "  --by <name>        grouping column for t2\n"
    "  --a G1 --b G2      pick two groups when --by has more than two levels\n"
    "  --with <name>      second column for the paired test\n"
    "  --no-welch         use pooled-variance t instead of Welch\n"
    "  --expected <list>  expected proportions for chisq (default: uniform)\n"
    "  --alpha <A>        significance level (default 0.05)\n"
    "  --json             machine-readable JSON output\n";

const char* kOutliersUsage =
    "Usage:\n"
    "  statsengine outliers --file <csv> --column <name>\n"
    "                       [--method z|iqr|mad] [--threshold T]\n"
    "\n"
    "Flag univariate outliers via z-score (|z|>3), Tukey IQR fences (1.5*IQR) or\n"
    "the MAD-based modified z-score (|0.6745*(x-med)/MAD|>3.5).\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --column <name>    numeric column to scan (required)\n"
    "  --method <name>    z (default) | iqr | mad\n"
    "  --threshold <T>    override the method default cut-off\n"
    "  --json             machine-readable JSON output\n";

const char* kCleanUsage =
    "Usage:\n"
    "  statsengine clean --file <csv>\n"
    "\n"
    "Data-quality report: per-column type, missing counts, unique values and\n"
    "notes (constant columns, high missingness, integer-like), plus file-level\n"
    "duplicate-row and memory statistics.\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --json             machine-readable JSON output\n";

const char* kSampleUsage =
    "Usage:\n"
    "  statsengine sample --file <csv> (--n N | --frac F) [--seed S] [--replace]\n"
    "\n"
    "Emit a random sample of rows as CSV on stdout. Sampling is deterministic\n"
    "for a given --seed (splitmix64; default seed 42).\n"
    "\n"
    "Options:\n"
    "  --file <path>      CSV file to load ('-' reads stdin)\n"
    "  --n <N>            number of rows to sample (without replacement)\n"
    "  --frac <F>         fraction of rows in (0, 1]\n"
    "  --seed <S>         PRNG seed (default 42)\n"
    "  --replace          sample with replacement (requires --n)\n"
    "  --json             machine-readable JSON output\n";

}  // namespace

const std::vector<CommandInfo>& commandRegistry() {
  static const std::vector<CommandInfo> registry = {
      {"describe", "count/mean/stddev/quartiles/skew/kurtosis per column", kDescribeUsage,
       cmdDescribe},
      {"hist", "ASCII histogram of a numeric column", kHistogramUsage, cmdHistogram},
      {"corr", "Pearson/Spearman/Kendall correlation matrix", kCorrelationUsage, cmdCorrelation},
      {"regress", "OLS regression (simple + multivariate) with diagnostics", kRegressionUsage,
       cmdRegression},
      {"ma", "moving average / EMA / differencing on a column", kMovingAverageUsage,
       cmdMovingAverage},
      {"test", "t-tests and chi-square goodness-of-fit with p-values", kTestUsage, cmdTest},
      {"anova", "one-way analysis of variance (F test) across groups", kAnovaUsage, cmdAnova},
      {"outliers", "z-score / IQR / MAD outlier detection", kOutliersUsage, cmdOutliers},
      {"clean", "column types, missing values and data-quality report", kCleanUsage, cmdClean},
      {"sample", "seeded random sample of rows", kSampleUsage, cmdSample},
  };
  return registry;
}

const CommandInfo* findCommand(const std::string& name) {
  for (const CommandInfo& info : commandRegistry()) {
    if (name == info.name) return &info;
  }
  return nullptr;
}

void printHelp(std::ostream& os, bool color) {
  std::size_t nameWidth = 0;
  for (const CommandInfo& info : commandRegistry()) {
    nameWidth = std::max(nameWidth, std::strlen(info.name));
  }

  os << colors::paint("statsengine", color ? "1" : nullptr) << " " << kStatsengineVersion
     << " — CSV statistical analysis from the command line\n\n";
  os << "Usage:\n  statsengine <command> --file data.csv --column x [--options]\n\n";
  os << colors::paint("Commands", color ? "1" : nullptr) << "\n";
  for (const CommandInfo& info : commandRegistry()) {
    std::string name = info.name;
    name += std::string(nameWidth - name.size(), ' ');
    os << "  " << colors::paint(name, color ? "36" : nullptr) << " " << info.summary << "\n";
  }
  os << "\n" << colors::paint("Global options", color ? "1" : nullptr) << "\n";
  os << "  --json        machine-readable JSON output (per command)\n";
  os << "  --ascii       ASCII borders/bars instead of Unicode box drawing\n";
  os << "  --no-color    disable ANSI colours (also honours NO_COLOR)\n";
  os << "  --color       force ANSI colours even when piped\n";
  os << "  --help        show this help, or '<command> --help' for details\n";
  os << "  --version     print version and exit\n\n";
  os << colors::paint("Exit codes", color ? "1" : nullptr) << "\n";
  os << "  0 success   1 runtime error (I/O, CSV, statistics)   2 usage error\n\n";
  os << colors::paint("Examples", color ? "1" : nullptr) << "\n";
  os << colors::paint("  statsengine describe --file examples/sales.csv\n",
                      color ? "2" : nullptr)
     << colors::paint("  statsengine hist --file examples/sales.csv --column revenue --rule fd\n",
                      color ? "2" : nullptr)
     << colors::paint("  statsengine corr --file examples/sales.csv --columns units,revenue\n",
                      color ? "2" : nullptr)
     << colors::paint(
            "  statsengine regress --file examples/sales.csv --target revenue --columns units\n",
            color ? "2" : nullptr)
     << colors::paint(
            "  statsengine test --file examples/sales.csv --test t2 --column revenue --by region\n",
            color ? "2" : nullptr)
     << colors::paint("  statsengine sample --file examples/sales.csv --n 100 --seed 7\n",
                      color ? "2" : nullptr);
}

void printCommandHelp(std::ostream& os, const CommandInfo& command, bool color) {
  os << colors::paint(std::string("statsengine ") + command.name, color ? "1" : nullptr) << " — "
     << command.summary << "\n\n";
  os << command.usage << "\n";
  os << colors::paint("Global flags (--json, --ascii, --no-color, --color) work with every "
                      "command.\n",
                      color ? "2" : nullptr);
}

// ------------------------------------------------------------ describe --

int cmdDescribe(const Args& args) {
  args.rejectUnknown({"file", "column", "columns", "all"});
  args.expectNoPositionals(0);
  const DataFrame df = loadFrame(args);

  std::vector<std::string> names = selectedColumnNames(args);
  if (names.empty() || args.has("all")) {
    for (std::size_t idx : df.numericColumnIndices()) names.push_back(df.column(idx).name());
  }
  if (names.empty()) {
    throw std::runtime_error("no numeric columns found in '" + args.get("file") + "'");
  }

  std::vector<DescriptiveStats> stats;
  for (const std::string& name : names) {
    const Column& col = requireNumericColumn(df, name);
    stats.push_back(describeColumn(col));
  }

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    writer.beginObject();
    writer.kv("command", "describe");
    writer.kv("file", args.getOr("file", "-"));
    writer.key("columns");
    writer.beginArray();
    for (const DescriptiveStats& s : stats) describeToJson(writer, s);
    writer.endArray();
    writer.endObject();
    std::cout << "\n";
    return 0;
  }
  std::cout << renderDescribeTable(stats, !args.asciiMode(), colors::enabled()) << "\n";
  return 0;
}

// ---------------------------------------------------------------- hist --

int cmdHistogram(const Args& args) {
  args.rejectUnknown({"file", "column", "rule", "bins", "bin-width", "bar-width", "cumulative"});
  args.expectNoPositionals(0);
  args.require("column", "numeric column to bin");

  HistogramOptions options;
  const std::string rule = args.getOr("rule", "auto");
  if (!binningRuleFromName(rule, options.rule)) {
    throw ArgsError("unknown binning rule '" + rule +
                    "' (expected auto, sturges, fd, scott or fixed)");
  }
  options.bins = static_cast<int>(args.getInt("bins", 0));
  options.binWidth = args.getDouble("bin-width", 0.0);
  options.barWidth = static_cast<int>(args.getInt("bar-width", 46));
  if (options.barWidth < 8 || options.barWidth > 120) {
    throw ArgsError("--bar-width must be between 8 and 120");
  }
  options.cumulative = args.has("cumulative");

  const DataFrame df = loadFrame(args);
  const std::string columnName = args.get("column");
  const Column& col = requireNumericColumn(df, columnName);
  const HistogramResult result = computeHistogram(col.numeric(), options);

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    histogramToJson(writer, result, col.name());
    std::cout << "\n";
    return 0;
  }
  std::cout << renderHistogram(result, col.name(), !args.asciiMode(), colors::enabled());
  return 0;
}

// ------------------------------------------------------------- outliers --

int cmdOutliers(const Args& args) {
  args.rejectUnknown({"file", "column", "method", "threshold"});
  args.expectNoPositionals(0);
  args.require("column", "numeric column to scan");

  OutlierOptions options;
  const std::string method = args.getOr("method", "z");
  if (!outlierMethodFromName(method, options.method)) {
    throw ArgsError("unknown outlier method '" + method + "' (expected z, iqr or mad)");
  }
  options.threshold = args.getDouble("threshold", 0.0);
  if (options.threshold < 0.0) {
    throw ArgsError("--threshold must be non-negative");
  }

  const DataFrame df = loadFrame(args);
  const std::string columnName = args.get("column");
  const Column& col = requireNumericColumn(df, columnName);
  const OutlierReport report = detectOutliers(col.numeric(), options);

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    outliersToJson(writer, report, col.name());
    std::cout << "\n";
    return 0;
  }
  std::cout << renderOutliers(report, col.name(), !args.asciiMode(), colors::enabled(), 25);
  return 0;
}

// ---------------------------------------------------------------- clean --

int cmdClean(const Args& args) {
  args.rejectUnknown({"file"});
  args.expectNoPositionals(0);
  const DataFrame df = loadFrame(args);

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    writer.beginObject();
    writer.kv("file", args.getOr("file", "-"));
    writer.kvCount("rows", df.rows());
    writer.kvCount("columns", df.cols());
    writer.kvCount("missing_cells", df.missingCells());
    writer.key("columns");
    writer.beginArray();
    for (const Column& col : df.columns()) {
      writer.beginObject();
      writer.kv("name", col.name());
      writer.kv("kind", columnKindName(col.kind()));
      writer.kv("count", static_cast<long long>(col.size()));
      writer.kv("missing", static_cast<long long>(col.missingCount()));
      writer.kv("unique", static_cast<long long>(col.categories().size()));
      writer.kvBool("integer_like", col.isIntegerLike());
      writer.endObject();
    }
    writer.endArray();
    writer.endObject();
    std::cout << "\n";
    return 0;
  }

  std::vector<std::pair<std::string, std::string>> summary;
  summary.emplace_back("file", args.getOr("file", "-"));
  summary.emplace_back("rows", std::to_string(df.rows()));
  summary.emplace_back("columns", std::to_string(df.cols()));
  summary.emplace_back("missing cells", std::to_string(df.missingCells()));
  const double completeness =
      df.rows() * df.cols() > 0
          ? 1.0 - static_cast<double>(df.missingCells()) /
                      static_cast<double>(df.rows() * df.cols())
          : 1.0;
  summary.emplace_back("completeness",
                       progressBar(completeness, 20, !args.asciiMode(), colors::enabled()) +
                           " " + fmtPercent(completeness, 2));
  summary.emplace_back("memory", fmtBytes(df.bytesUsed()));

  // Duplicate rows (exact raw match).
  std::map<std::vector<std::string>, long long> rowCounts;
  for (std::size_t i = 0; i < df.rows(); ++i) ++rowCounts[df.rowAt(i)];
  long long duplicates = 0;
  for (const auto& entry : rowCounts) duplicates += entry.second - 1;
  summary.emplace_back("duplicate rows", std::to_string(duplicates));

  Table table;
  table.setTitle("Column profile");
  table.addColumn("column", Align::Left);
  table.addColumn("type", Align::Left);
  table.addColumn("count", Align::Right);
  table.addColumn("missing", Align::Right);
  table.addColumn("missing %", Align::Right);
  table.addColumn("unique", Align::Right);
  table.addColumn("notes", Align::Left);

  for (const Column& col : df.columns()) {
    std::vector<std::string> notes;
    const std::size_t count = col.size();
    const double missingFraction =
        count > 0 ? static_cast<double>(col.missingCount()) / static_cast<double>(count) : 0.0;
    const std::size_t unique = col.categories().size();
    if (col.missingCount() == count && count > 0) notes.push_back("all values missing");
    else if (missingFraction > 0.30) notes.push_back("high missingness (>30%)");
    if (count > 0 && col.missingCount() < count && unique == 1) notes.push_back("constant column");
    if (col.kind() == ColumnKind::Numeric && col.isIntegerLike()) notes.push_back("integer-like");
    if (col.kind() == ColumnKind::Date) notes.push_back("ISO dates");
    if (notes.empty()) notes.push_back("-");

    table.addRow({col.name(), columnKindName(col.kind()), std::to_string(count),
                  std::to_string(col.missingCount()), fmtPercent(missingFraction, 1),
                  std::to_string(unique), joinList(notes, ", ")});
  }
  table.setFooter("missing tokens: " + joinList(CsvOptions{}.missingTokens, ", "));

  std::cout << keyValueBlock(summary, colors::enabled()) << "\n" << table.render(!args.asciiMode());
  return 0;
}

// --------------------------------------------------------------- sample --

int cmdSample(const Args& args) {
  args.rejectUnknown({"file", "n", "frac", "seed", "replace"});
  args.expectNoPositionals(0);
  const DataFrame df = loadFrame(args);

  const bool hasN = args.has("n");
  const bool hasFrac = args.has("frac");
  if (hasN == hasFrac) {
    throw ArgsError("specify exactly one of --n <rows> or --frac <fraction>");
  }
  const bool withReplacement = args.has("replace");
  if (withReplacement && !hasN) {
    throw ArgsError("--replace requires --n");
  }
  long long n = 0;
  if (hasN) {
    n = args.getInt("n", -1);
    if (n < 0) throw ArgsError("--n must be non-negative");
  } else {
    const double frac = args.getDouble("frac", -1.0);
    if (!(frac > 0.0) || frac > 1.0) throw ArgsError("--frac must be in (0, 1]");
    n = static_cast<long long>(std::llround(frac * static_cast<double>(df.rows())));
  }
  if (!withReplacement && static_cast<std::size_t>(n) > df.rows()) {
    throw ArgsError("cannot sample " + std::to_string(n) + " rows without replacement from a " +
                    std::to_string(df.rows()) + "-row file");
  }
  const std::uint64_t seed = static_cast<std::uint64_t>(args.getInt("seed", 42));
  Random rng(seed);

  std::vector<std::size_t> picked;
  if (withReplacement) {
    picked.reserve(static_cast<std::size_t>(n));
    for (long long i = 0; i < n; ++i) {
      picked.push_back(static_cast<std::size_t>(rng.nextInt(0, static_cast<int>(df.rows()) - 1)));
    }
  } else {
    std::vector<std::size_t> all(df.rows());
    for (std::size_t i = 0; i < df.rows(); ++i) all[i] = i;
    rng.shuffle(all);
    picked.assign(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(n));
  }

  if (args.jsonMode()) {
    json::Writer writer(std::cout);
    writer.beginObject();
    writer.kv("file", args.getOr("file", "-"));
    writer.kv("n", n);
    writer.kv("seed", static_cast<long long>(seed % 9007199254740992LL));
    writer.kvBool("replace", withReplacement);
    writer.key("rows");
    writer.beginArray();
    for (std::size_t idx : picked) {
      writer.beginArray();
      for (const std::string& cell : df.rowAt(idx)) writer.value(cell);
      writer.endArray();
    }
    writer.endArray();
    writer.endObject();
    std::cout << "\n";
    return 0;
  }

  CsvWriter writer(std::cout);
  writer.row(df.columnNames());
  for (std::size_t idx : picked) writer.row(df.rowAt(idx));
  return 0;
}

}  // namespace se
