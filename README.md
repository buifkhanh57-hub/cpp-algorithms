# statsengine

**A dependency-free C++17 command-line toolkit for CSV statistical analysis — descriptive stats, histograms, correlation, OLS regression, hypothesis tests, ANOVA, outliers and data-quality checks, all in one static binary.**

![language](https://img.shields.io/badge/language-C%2B%2B17-00599C?logo=cplusplus&logoColor=white)
![standard](https://img.shields.io/badge/standard-C%2B%2B17-blue)
![build](https://img.shields.io/badge/build-make-success-brightgreen)
![warnings](https://img.shields.io/badge/warnings-0-success)
![dependencies](https://img.shields.io/badge/dependencies-0-8A2BE2)
![tests](https://img.shields.io/badge/tests-361%20checks-passing-brightgreen)
![platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)
![license](https://img.shields.io/badge/license-MIT-green)

---

## Overview

`statsengine` turns the boring half of data work — *look at the CSV before you trust the CSV* —
into one fast command. Point it at any delimited file and it will profile columns, plot ASCII
histograms, correlate variables, fit least-squares models, run t-tests / chi-square / one-way
ANOVA, flag outliers with three different rules, profile missingness and sample rows with a
reproducible seed.

Everything is computed in-process from first principles: the project ships **zero external
dependencies** — no BLAS, no Boost, no CSV library. The special functions (log-gamma,
regularized incomplete beta, Student-t and F distributions) are implemented on top of
Lentz-modified continued fractions, so p-values and confidence intervals come straight out of
the same 4 MB binary you compile with `make`.

Highlights at a glance:

- **10 subcommands** sharing one consistent flag grammar and output style.
- **Human output** (Unicode box tables + ASCII bars) or **`--json`** machine output for every
  single command — same numbers, two renderers.
- **Honest statistics**: ddof=1 standard deviations, bias-corrected skewness/kurtosis (G1/G2),
  Welch-by-default two-sample t-tests, Satterthwaite degrees of freedom, Huber-robust options
  for outlier detection.
- **Fail loudly, exit cleanly**: three exit codes (0/1/2) mapped precisely to runtime vs usage
  errors, Levenshtein "did you mean" hints for typos.

## Features

- **Descriptive statistics** — count/missing/mean/stddev/min/quartiles/max/IQR/skewness/
  kurtosis/CV per numeric column, computed in one streaming pass with running moments.
- **ASCII histograms** — Sturges, Scott, Freedman–Diaconis, Rice, square-root or fixed-width
  binning, with automatic fallbacks when a rule degenerates (IQR = 0) and optional cumulative
  bars.
- **Correlation matrices** — Pearson, Spearman (rank-based) and Kendall tau-a, with pairwise-
  complete observations and significance stars from the exact t approximation.
- **OLS regression** — simple and multivariate via normal equations with Gauss-Jordan
  inversion; full coefficient table (SE, t, p, 95% CI), R²/adjusted R², F-test, RMSE,
  Durbin–Watson, optional no-intercept mode and out-of-sample `--predict`.
- **Time-series transforms** — trailing simple moving average, EMA (span- or alpha-defined),
  and lagged differencing.
- **Hypothesis tests** — one-sample, two-sample (Welch or pooled) and paired t-tests plus
  chi-square goodness-of-fit, all reporting test statistic, df, p-value and effect size.
- **One-way ANOVA** — full ANOVA table (SS/df/MS/F/p), η² and ω² effect sizes, per-group
  summaries.
- **Outlier detection** — classic z-score, Tukey IQR fences and the robust modified z-score
  based on the MAD, each with a ranked flagged-values table.
- **Data quality report** — per-column type inference (numeric / date / categorical / text),
  missing-value tokens, uniqueness, duplicate rows, completeness bar and memory estimate.
- **Seeded sampling** — headless random row samples (fixed `--n` or fraction `--frac`) with an
  optional `--replace` mode; same seed, same sample, every run.

## Requirements

| Dependency | Version | Notes |
|------------|---------|-------|
| C++ compiler | GCC ≥ 8 or Clang ≥ 9 | anything with solid C++17 support |
| `make` | any POSIX make | GNU make recommended |
| OS | Linux, macOS (Windows via MSYS2/WSL) | no platform-specific code |

There are **no third-party libraries**. If `g++ --version` works, `statsengine` builds.

## Installation

```bash
git clone <this repo>
cd statsengine
make            # compiles ./build/statsengine with -std=c++17 -O2 -Wall -Wextra
make test       # builds and runs the assert-based test suite (361 checks)
make run        # build + print the top-level help
make clean      # removes the build/ directory
```

The build is a single `g++` invocation per translation unit — no configure step, no package
manager, no code generation. A full rebuild takes a few seconds.

## Quick Start

```bash
# Profile every numeric column of the bundled example data
./build/statsengine describe --file examples/sales.csv

# ASCII histogram of one column, Freedman–Diaconis binning
./build/statsengine hist --file examples/sales.csv --column revenue --rule fd

# Correlation matrix with significance stars
./build/statsengine corr --file examples/sales.csv --columns units,revenue

# Regress revenue on units (OLS with intercept + full diagnostics)
./build/statsengine regress --file examples/sales.csv --target revenue --columns units

# Welch two-sample t-test between two regions
./build/statsengine test --file examples/sales.csv --test t2 \
    --column revenue --by region --a north --b south

# Same numbers as JSON for the pipeline
./build/statsengine --json describe --file examples/sales.csv --column revenue
```

Recover-the-slope sanity check (exact `y = 2x + 1`, no noise):

```bash
$ printf 'x,y\n' > /tmp/lin.csv
$ for i in $(seq 1 100); do echo "$i,$((2*i+1))" >> /tmp/lin.csv; done
$ ./build/statsengine regress --file /tmp/lin.csv --target y --columns x

  target        : y
  observations  : 100
  parameters    : 2 (with intercept)
  R²           : 1.000000
  adjusted R²  : 1.000000
  residual SE   : 0.000000
  RMSE          : 0.000000
  Durbin-Watson : n/a

Coefficients
┌─────────────┬─────────────┬───────────┬─────┬─────────┬──────────────────┐
 │ term        │ coefficient │ std error │   t │ p-value │           95% CI │
├─────────────┼─────────────┼───────────┼─────┼─────────┼──────────────────┤
 │ (intercept) │    1.000000 │  0.000000 │ n/a │     n/a │ [1.0000, 1.0000] │
 │ x           │    2.000000 │  0.000000 │ n/a │     n/a │ [2.0000, 2.0000] │
└─────────────┴─────────────┴───────────┴─────┴─────────┴──────────────────┘
```

Slope 2, intercept 1, R² = 1 — exactly as designed. (`t`/`p` are `n/a` because a perfect fit
has zero residual variance, which makes the ratio β/SE undefined.)

## Usage

### Command reference

Every command requires `--file` (use `-` to read CSV from stdin) plus the flags below.
Global flags (`--json`, `--ascii`, `--no-color`, `--color`, `--help`, `--version`) work with
every command, **before or after** the command name.

| Command | Purpose | Key options |
|---------|---------|-------------|
| `describe` | Per-column descriptive statistics | `--column <c>` / `--columns a,b` / `--all` |
| `hist` | ASCII histogram of a numeric column | `--rule sturges\|scott\|fd\|rice\|sqrt\|fixed`, `--bins N`, `--bin-width W`, `--bar-width N`, `--cumulative` |
| `corr` | Correlation matrix | `--columns a,b[,c…]`, `--method pearson\|spearman\|kendall` |
| `regress` | OLS regression + diagnostics | `--target <y>`, `--columns x1,x2`, `--no-intercept`, `--predict <"1,2;3,4">`, `--predict-file <csv>` |
| `ma` | Moving average / EMA / differencing | `--method sma\|ema\|diff`, `--window N`, `--span N`, `--alpha A`, `--lag N` |
| `test` | t-tests and chi-square GOF | `--test t1\|t2\|paired\|chisq`, `--mu M`, `--by <g>`, `--a G1 --b G2`, `--with <col>`, `--expected p1,p2`, `--no-welch`, `--alpha A` |
| `anova` | One-way ANOVA (F test) | `--column <x>`, `--by <g>`, `--only g1,g2`, `--alpha A` |
| `outliers` | Outlier detection | `--method z\|iqr\|mad`, `--threshold T` |
| `clean` | Data-quality / missing-value report | — |
| `sample` | Seeded random row sample | `--n N` or `--frac F`, `--seed S`, `--replace` |

### Realistic output samples

`describe` (sales.csv, 500 rows):

```text
Descriptive statistics
┌──────────┬──────────┬────────────┐
 │ metric   │    units │    revenue │
├──────────┼──────────┼────────────┤
 │ count    │      500 │        500 │
 │ missing  │        0 │          0 │
 │ mean     │   22.508 │  422.25176 │
 │ stddev   │ 7.584787 │ 296.920934 │
 │ min      │        1 │      13.49 │
 │ p25      │       18 │     219.25 │
 │ p50      │       22 │    315.295 │
 │ p75      │       27 │    570.085 │
 │ max      │       56 │    1986.72 │
 │ iqr      │        9 │    350.835 │
 │ skewness │ 0.461204 │   1.611164 │
 │ kurtosis │ 1.161369 │   3.205449 │
 │ cv       │ 0.336982 │   0.703185 │
└──────────┴──────────┴────────────┘
```

`corr` (Pearson, with significance stars):

```text
Correlation matrix (pearson)
┌─────────┬──────────┬──────────┐
 │         │    units │  revenue │
├─────────┼──────────┼──────────┤
 │ units   │    1.000 │ 0.509*** │
 │ revenue │ 0.509*** │    1.000 │
└─────────┴──────────┴──────────┘
pairwise-complete observations; * p<0.05  ** p<0.01  *** p<0.001
```

`regress` on real (noisy) data — the diagnostics behave like a textbook printout:

```text
  target        : revenue
  observations  : 500
  parameters    : 2 (with intercept)
  R²           : 0.259315
  adjusted R²  : 0.257828
  F-statistic   : 174.3505 on 1 and 498 DF, p=0.000000
  residual SE   : 255.795691
  RMSE          : 255.283587
  Durbin-Watson : 1.9565

Coefficients
┌─────────────┬─────────────┬───────────┬────────┬──────────┬─────────────────────┐
 │ term        │ coefficient │ std error │      t │  p-value │              95% CI │
├─────────────┼─────────────┼───────────┼────────┼──────────┼─────────────────────┤
 │ (intercept) │  -26.439971 │ 35.854887 │ -0.737 │ 0.461217 │ [-96.8855, 44.0055] │
 │ units       │   19.934767 │  1.509731 │ 13.204 │ 0.000000 │  [16.9685, 22.9010] │
└─────────────┴─────────────┴───────────┴────────┴──────────┴─────────────────────┘
```

`test --test chisq` (goodness-of-fit against uniform region shares):

```text
  test      : Chi-square goodness-of-fit
  H0        : categories follow the expected proportions
  statistic : 51.808
  df        : 3
  p-value   : 0.000000
  alpha     : 0.050
  n         : 500
  phi       : 0.3219

=> reject H0 at alpha=0.050 (statistically significant)
  - observed: east=116, north=180, south=136, west=68
```

`clean` (users.csv, 200 rows with 9 missing `age` cells):

```text
  file           : examples/users.csv
  rows           : 200
  columns        : 6
  missing cells  : 9
  completeness   : ████████████████████ 99.25%
  duplicate rows : 0

Column profile
┌───────────────────┬─────────────┬───────┬─────────┬───────────┬────────┬──────────────┐
 │ column            │ type        │ count │ missing │ missing % │ unique │ notes        │
├───────────────────┼─────────────┼───────┼─────────┼───────────┼────────┼──────────────┤
 │ signup_date       │ date        │   200 │       0 │      0.0% │    180 │ ISO dates    │
 │ plan              │ categorical │   200 │       0 │      0.0% │      3 │ -            │
 │ age               │ numeric     │   200 │       9 │      4.5% │     40 │ -            │
 │ active_days       │ numeric     │   200 │       0 │      0.0% │     66 │ integer-like │
 └───────────────────┴─────────────┴───────┴─────────┴───────────┴────────┴──────────────┘
```

`outliers --method iqr` (Tukey fences, ranked table):

```text
  column    : revenue
  method    : iqr
  threshold : 1.500
  n         : 500
  center    : 315.295
  spread    : 350.835
  fences    : [-307.0025, 1096.3375]
  flagged   : 19 (3.80%)
```

`sample` is deterministic — same `--seed`, same rows:

```text
$ ./build/statsengine sample --file examples/sales.csv --n 5 --seed 7
date,region,product,units,revenue
2025-07-01,north,alpha,35,335.20
2025-10-19,east,gamma,30,247.51
2025-07-11,south,beta,17,419.09
2025-04-15,north,beta,17,477.60
2025-12-30,east,alpha,17,228.63
```

`--json` mode (validated, one line, non-finite values become `null`):

```json
{"command":"describe","file":"examples/sales.csv","columns":[{"column":"revenue","count":500,
"missing":0,"sum":211125.88,"mean":422.25176,"stddev":296.92093375,"min":13.49,"q1":219.25,
"median":315.295,"q3":570.085,"max":1986.72,"iqr":350.835,"skewness":1.61116381488,...}]}
```

## Statistical methods reference

Plain-text formulas, exactly as implemented:

| Quantity | Formula / method |
|----------|------------------|
| Mean | x̄ = (1/n) Σ xᵢ |
| Stddev (sample) | s = sqrt( Σ(xᵢ−x̄)² / (n−1) ) — ddof = 1 |
| Quartiles | linear interpolation on sorted data, q(p) = x⌊(n−1)p⌋ + fractional part |
| Skewness | bias-corrected G1 = [n/((n−1)(n−2))] Σ((xᵢ−x̄)/s)³ |
| Kurtosis | bias-corrected G2 = excess kurtosis adjusted by sample-size factors |
| CV | s / x̄ (undefined for x̄ = 0) |
| Sturges bins | k = ⌈log₂ n⌉ + 1 |
| Scott width | h = 3.49 · σ · n^(−1/3) |
| Freedman–Diaconis | h = 2 · IQR · n^(−1/3), falls back to Sturges when IQR = 0 |
| Pearson r | Σ(xᵢ−x̄)(yᵢ−ȳ) / sqrt(Σ(xᵢ−x̄)² Σ(yᵢ−ȳ)²) |
| Spearman ρ | Pearson r applied to mid-rank-transformed data |
| Kendall τ | concordant − discordant pairs, normal-score p approximation |
| OLS | minimize Σ(yᵢ−xᵢβ)² → β = (XᵀX)⁻¹Xᵀy via Gauss–Jordan with partial pivoting |
| Std errors | Var(β̂) = σ̂² (XᵀX)⁻¹, σ̂² = SSE / (n − p) |
| t statistic | tⱼ = β̂ⱼ / SE(β̂ⱼ), two-tailed p from Student-t with n − p df |
| R² / adj R² | 1 − SSE/SST (centered); adj = 1 − (1−R²)(n−1)/(n−p−1) |
| F test | F = (SSR/dfModel) / (SSE/dfResidual), p from F CDF |
| Durbin–Watson | Σ(eᵢ−eᵢ₋₁)² / Σeᵢ² (undefined → `n/a` on zero residuals) |
| Welch SE | sqrt(s₁²/n₁ + s₂²/n₂), Satterthwaite df = (SE²)² / Σ (sᵢ²/nᵢ)²/(nᵢ−1) |
| Pooled t | s_p² = ((n₁−1)s₁² + (n₂−1)s₂²)/(n₁+n₂−2), df = n₁+n₂−2 |
| Paired t | one-sample t on dᵢ = xᵢ − yᵢ |
| Cohen's d | (x̄₁ − x̄₂) / s_pooled |
| Chi-square GOF | χ² = Σ (Oᵢ−Eᵢ)²/Eᵢ, df = categories − 1; Cramér's φ = sqrt(χ²/n) |
| One-way ANOVA | SSB = Σ nⱼ(x̄ⱼ−x̄)², SSW = Σ Σ(xᵢⱼ−x̄ⱼ)², F = (SSB/dfB)/(SSW/dfW) |
| Effect sizes | η² = SSB/SST, ω² = (SSB − dfB·MSW)/(SST + MSW) |
| z-score outlier | |xᵢ − x̄| / s > threshold |
| IQR outlier | outside [q1 − 1.5·IQR, q3 + 1.5·IQR] (scale via `--threshold`) |
| Modified z | 0.6745(xᵢ − median)/MAD > threshold (default 3.5) |
| EMA | y₀ = x₀; yₜ = α·xₜ + (1−α)·yₜ₋₁, α = 2/(span+1) |
| SMA | trailing mean of the last `window` values (expanding before full) |
| Differencing | Δxₜ = xₜ − xₜ₋lag |
| Special functions | logGamma (Lanczos), regularized incomplete beta (Lentz continued fraction), regularized gamma P/Q (series + continued fraction) |

## Configuration and exit codes

No config files — behaviour is entirely flag-driven, plus standard environment:

| Mechanism | Effect |
|-----------|--------|
| `NO_COLOR` | disables ANSI colours (https://no-color.org convention) |
| `--no-color` / `--color` | explicit override; colours auto-disable when piped |
| `--ascii` | ASCII `+---+` borders and `#` bars instead of Unicode box drawing |
| `--json` | one-line JSON per command; `NaN`/`±inf` serialise as `null` |
| `--file -` | read CSV from stdin |

Exit codes:

| Code | Meaning | Examples |
|------|---------|----------|
| `0` | success | normal output, `--help`, `--version` |
| `1` | runtime error | cannot open file, malformed CSV, singular model matrix |
| `2` | usage error | unknown command/option, wrong flag type, missing `--file`, non-numeric column |

Error messages go to **stderr** and always start with `error:`; typos get a Levenshtein
suggestion (`unknown command 'regrss' (did you mean 'regress'?)`).

## Project Structure

```text
statsengine/
├── Makefile                      build orchestration (build/test/run/clean), 44 lines
├── LICENSE                       MIT, Copyright (c) 2026 Bui Bao Khanh
├── README.md                     this file
├── .gitignore                    build/, objects, binaries, logs
├── examples/
│   ├── sales.csv                 500 rows × 5 cols: date, region, product, units, revenue
│   └── users.csv                 200 rows × 6 cols: signup_date, plan, age, active_days, ...
├── include/statsengine/          public headers (one per module + shared utilities)
│   ├── anova.hpp                 one-way ANOVA API + result types            (92)
│   ├── args.hpp                  argument parser (Args, ArgsError)           (96)
│   ├── colors.hpp                ANSI paint, NO_COLOR, visible width         (96)
│   ├── commands.hpp              command registry + per-command help         (74)
│   ├── correlation.hpp           pearson/spearman/kendall API                (66)
│   ├── csv.hpp                   RFC-4180-ish CSV reader (quotes, CRLF)     (107)
│   ├── dataframe.hpp             typed column store + type inference        (140)
│   ├── describe.hpp              descriptive stats + RunningMoments         (116)
│   ├── histogram.hpp             binning rules + histogram result           (81)
│   ├── hypothesis.hpp            t/chi2/F distributions + tests             (96)
│   ├── json.hpp                  streaming JSON writer (NaN→null)           (195)
│   ├── outliers.hpp              z / IQR / MAD detection API                (76)
│   ├── regression.hpp            OLS fit, predict, render, JSON              (95)
│   ├── table.hpp                 box-drawing table renderer                  (65)
│   ├── timeseries.hpp            SMA/EMA/diff transforms                     (61)
│   └── util.hpp                  string/number helpers, keyValueBlock, etc. (299)
├── src/                          implementations (one .cpp per module)
│   ├── anova.cpp                 ANOVA table, η²/ω², group summaries        (215)
│   ├── args.cpp                  tokenizer: long opts, --key=value, bools   (158)
│   ├── commands_data.cpp         describe/hist/outliers/clean/sample cmds   (564)
│   ├── commands_stats.cpp        corr/regress/ma/test/anova cmds            (497)
│   ├── correlation.cpp           rank transforms, kendall pairs counting    (271)
│   ├── csv.cpp                   streaming parser + CsvParseError           (218)
│   ├── dataframe.cpp             column typing, numeric coercion, missing   (297)
│   ├── describe.cpp              one-pass moments, G1/G2, render + JSON     (216)
│   ├── histogram.cpp             bin rules, bars, cumulative mode           (279)
│   ├── hypothesis.cpp            special functions + t/chi2/paired tests    (441)
│   ├── main.cpp                  entry point, dispatch, exit-code mapping   (112)
│   ├── outliers.cpp              three detection rules + ranked table       (232)
│   ├── regression.cpp            normal equations, invert, diagnostics      (362)
│   ├── table.cpp                 box renderer, keyValueBlock, progressBar   (178)
│   └── timeseries.cpp            sma/ema/diff engines                       (169)
└── tests/
    └── test_core.cpp             361-check assert suite over every module  (1189)
```

Total: ~7,900 lines including headers, tests and build files (≥ 5,800 target), zero
dependencies, zero warnings under `-std=c++17 -Wall -Wextra`.

## Architecture

The codebase is a strict layer cake; data flows upward only:

1. **I/O layer** — `csv.cpp` streams rows out of a file or stdin with quote/CRLF handling;
   `dataframe.cpp` turns raw cells into typed columns (`numeric`, `date`, `categorical`,
   `text`) and records missing tokens (`""`, `NA`, `N/A`, `null`, `NaN`, `?`).
2. **Math layer** — pure functions with no I/O: `describe.cpp` (one-pass running moments),
   `hypothesis.cpp` (log-gamma, incomplete beta/gamma, t/χ²/F distributions), `regression.cpp`
   (normal equations, Gauss–Jordan inversion, diagnostics), `correlation.cpp`, `anova.cpp`,
   `outliers.cpp`, `timeseries.cpp`. These modules never touch `std::cout` for data — they
   return result structs.
3. **Presentation layer** — `table.cpp` (box renderer), `colors.cpp` (ANSI + width-safe
   painting), `json.hpp` (streaming writer that tracks nesting and serialises non-finite
   doubles as `null`). Every result struct has a `render*` and a `*ToJson` pair.
4. **Interface layer** — `args.cpp` parses `--flag value` / `--flag=value` / boolean flags
   (global flags are strictly boolean, so `--json describe ...` works in any position);
   `commands_*.cpp` validate, dispatch and format; `main.cpp` maps exceptions to exit codes
   (`ArgsError`/`ColumnError` → 2, `CsvParseError`/`std::exception` → 1).

Design decisions worth knowing:

- **Two renderers, one truth.** Statistical results are computed once into plain structs;
  text and JSON views are independent projections, so the two can never disagree on numbers.
- **Numerical safety first.** Matrix code pivots and refuses singular systems with a clear
  error instead of emitting garbage betas; degenerate-but-valid results (zero residual
  variance) print `n/a` for the undefined quantities (t, p, Durbin–Watson) rather than fake
  values, and JSON emits `null`.
- **Streaming where it counts.** Moments, ranks and histograms are computed in O(n) memory
  over the CSV stream; only the dataframe commands materialise full columns.
- **P-values from scratch.** The t and F distributions are evaluated through the regularized
  incomplete beta with Lentz continued fractions (verified against textbook critical values
  in the test suite to 1e-9).
- **Registry-driven CLI.** Each subcommand is a `CommandInfo{name, help, run}` entry; adding a
  command is one struct + one function, and help text, suggestions and dispatch follow
  automatically.

## Testing

`make test` builds `tests/test_core.cpp` against every module (excluding `main.cpp`, which
has its own `main`) into `build/test_core` and runs it. The suite is assert-based (no
framework) and currently performs **361 checks, all passing**, covering:

- **util / csv / dataframe** — splitting, strict number parsing, quoted CSV, type inference,
  missing-token handling;
- **describe** — hand-computed moments on toy vectors, G1/G2 against known values;
- **regression** — exact recovery of `y = 2x + 1`, multivariate `y = 3 + 2a − 5b`, no-intercept
  fits, t/SE/CI mutual consistency, collinearity → exception, `n ≤ p` → exception;
- **hypothesis** — textbook critical values (χ² 3.841/5.991/18.307, normal 1.960, F 4.965),
  t = 0 → p = 1, pooled vs Welch df behaviour, paired t;
- **correlation, outliers, timeseries, histogram, anova, json, args, table** — rank
  transforms, fence arithmetic, EMA recursion, bin-rule fallbacks, JSON well-formedness,
  boolean-global-flag parsing, Unicode border rendering and visible-width alignment.

The JSON mode of all ten subcommands is additionally exercised end-to-end (every command's
`--json` output parses with a strict JSON parser), and `regress` recovery of slope 2 /
intercept 1 from noise-free data is part of the regression suite.

## FAQ

**Q: Why not use an existing stats library?**
The point of the project is a small, auditable, dependency-free binary. Every special
function is ~50 lines of readable numerics you can verify against a textbook in minutes.

**Q: How large can my CSV be?**
There is no hard limit; memory scales with columns × rows only for the commands that need
random access (`clean`, `sample`, `corr`, `regress`). `describe` and `hist` stream.

**Q: My file is semicolon- or tab-separated.**
The sniffer accepts `,`, `;` and tab consistently; quoting follows the CSV conventions
(`"hello, world"`, doubled quotes). Convert exotic dialects before ingestion.

**Q: Why does a perfect fit show `n/a` for t and p?**
With SSE = 0 the standard error is 0, so t = β/SE is undefined (0/0 or ∞). statsengine
refuses to print meaningless infinities — the coefficient and its CI are still exact.

**Q: Why does `--json` sometimes show `"f_stat": null`?**
JSON has no `Infinity`/`NaN` literals. A perfect fit makes F = +∞, which serialises as
`null` — the same convention pandas uses.

**Q: Welch by default — how do I get the classic Student pooled t-test?**
Add `--no-welch` to `test --test t2`. The output footer always states which variance
treatment was used.

**Q: Are samples really reproducible?**
Yes. `sample` uses a seeded PRNG; the same `--seed` always yields the same rows for the same
input file, which makes the tool usable in CI snapshots.

**Q: Windows?**
Build under WSL, MSYS2 or any POSIX-compatible toolchain; the source is standard C++17 with
no platform APIs.

## Roadmap

- [ ] Weighted least squares and robust (Huber) regression
- [ ] Two-way ANOVA with interaction terms
- [ ] Mann–Whitney U / Wilcoxon signed-rank nonparametric tests
- [ ] Linear correlation with bootstrap confidence intervals
- [ ] `--out` / `--format md` for Markdown table output
- [ ] Multi-file `describe --compare` side-by-side profiling
- [ ] Optional BLAS backend for large multivariate fits
- [ ] Prebuilt binaries via CI for Linux/macOS

## License

Released under the [MIT License](LICENSE) — free to use, modify and distribute;
attribution required. Copyright (c) 2026 Bui Bao Khanh.

---
**by Bui Bao Khanh**
