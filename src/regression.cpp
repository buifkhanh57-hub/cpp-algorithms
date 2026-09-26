// statsengine — regression.cpp
//
// OLS fitting, linear algebra kernels and the regression report.
// See regression.hpp for the contract and error behavior.
#include "statsengine/regression.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "statsengine/colors.hpp"
#include "statsengine/hypothesis.hpp"
#include "statsengine/table.hpp"
#include "statsengine/util.hpp"

namespace se {

std::vector<double> solveLinearSystem(std::vector<std::vector<double>> a, std::vector<double> b) {
  const std::size_t n = b.size();
  if (a.size() != n) throw std::invalid_argument("solveLinearSystem: dimension mismatch");
  for (const std::vector<double>& row : a) {
    if (row.size() != n) throw std::invalid_argument("solveLinearSystem: matrix not square");
  }

  double scale = 0.0;
  for (const std::vector<double>& row : a) {
    for (double v : row) scale = std::max(scale, std::fabs(v));
  }
  const double eps = 1e-12 * std::max(1.0, scale);

  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    for (std::size_t row = col + 1; row < n; ++row) {
      if (std::fabs(a[row][col]) > std::fabs(a[pivot][col])) pivot = row;
    }
    if (std::fabs(a[pivot][col]) <= eps) {
      throw std::runtime_error("linear system is singular (matrix not invertible)");
    }
    if (pivot != col) {
      std::swap(a[pivot], a[col]);
      std::swap(b[pivot], b[col]);
    }
    for (std::size_t row = col + 1; row < n; ++row) {
      const double factor = a[row][col] / a[col][col];
      if (factor == 0.0) continue;
      for (std::size_t j = col; j < n; ++j) a[row][j] -= factor * a[col][j];
      b[row] -= factor * b[col];
    }
  }

  std::vector<double> x(n, 0.0);
  for (std::size_t i = n; i-- > 0;) {
    double sum = b[i];
    for (std::size_t j = i + 1; j < n; ++j) sum -= a[i][j] * x[j];
    x[i] = sum / a[i][i];
  }
  return x;
}

std::vector<std::vector<double>> invertMatrix(std::vector<std::vector<double>> a) {
  const std::size_t n = a.size();
  if (n == 0) return {};
  for (const std::vector<double>& row : a) {
    if (row.size() != n) throw std::invalid_argument("invertMatrix: matrix not square");
  }

  std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) inv[i][i] = 1.0;

  double scale = 0.0;
  for (const std::vector<double>& row : a) {
    for (double v : row) scale = std::max(scale, std::fabs(v));
  }
  const double eps = 1e-12 * std::max(1.0, scale);

  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    for (std::size_t row = col + 1; row < n; ++row) {
      if (std::fabs(a[row][col]) > std::fabs(a[pivot][col])) pivot = row;
    }
    if (std::fabs(a[pivot][col]) <= eps) {
      throw std::runtime_error("matrix is singular (cannot invert)");
    }
    if (pivot != col) {
      std::swap(a[pivot], a[col]);
      std::swap(inv[pivot], inv[col]);
    }
    const double diag = a[col][col];
    for (std::size_t j = 0; j < n; ++j) {
      a[col][j] /= diag;
      inv[col][j] /= diag;
    }
    for (std::size_t row = 0; row < n; ++row) {
      if (row == col) continue;
      const double factor = a[row][col];
      if (factor == 0.0) continue;
      for (std::size_t j = 0; j < n; ++j) {
        a[row][j] -= factor * a[col][j];
        inv[row][j] -= factor * inv[col][j];
      }
    }
  }
  return inv;
}

RegressionResult fitOls(const std::vector<double>& y,
                        const std::vector<std::vector<double>>& x,
                        const std::vector<std::string>& predictorNames,
                        const RegressionOptions& options) {
  const std::size_t n = y.size();
  if (x.size() != n) {
    throw std::invalid_argument("regression: y has " + std::to_string(n) + " rows but x has " +
                                std::to_string(x.size()));
  }
  const std::size_t k = predictorNames.size();
  for (const std::vector<double>& row : x) {
    if (row.size() != k) {
      throw std::invalid_argument("regression: expected " + std::to_string(k) +
                                  " predictors per row");
    }
  }
  const std::size_t p = k + (options.hasIntercept ? 1u : 0u);
  if (n <= p) {
    throw std::runtime_error("not enough complete observations (n=" + std::to_string(n) +
                             ") to fit " + std::to_string(p) + " parameters");
  }

  // Design matrix A (n x p), normal equations.
  std::vector<std::vector<double>> A(n, std::vector<double>(p, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t col = 0;
    if (options.hasIntercept) A[i][col++] = 1.0;
    for (std::size_t j = 0; j < k; ++j) A[i][col++] = x[i][j];
  }

  std::vector<std::vector<double>> xtx(p, std::vector<double>(p, 0.0));
  std::vector<double> xty(p, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t a = 0; a < p; ++a) {
      xty[a] += A[i][a] * y[i];
      for (std::size_t bidx = 0; bidx <= a; ++bidx) xtx[a][bidx] += A[i][a] * A[i][bidx];
    }
  }
  for (std::size_t a = 0; a < p; ++a) {
    for (std::size_t bidx = 0; bidx < a; ++bidx) xtx[bidx][a] = xtx[a][bidx];
  }

  RegressionResult result;
  result.targetName = "";
  result.predictorNames = predictorNames;
  result.hasIntercept = options.hasIntercept;
  result.n = static_cast<long long>(n);
  result.termNames.clear();
  if (options.hasIntercept) result.termNames.push_back("(intercept)");
  for (const std::string& name : predictorNames) result.termNames.push_back(name);
  result.dfModel = static_cast<long long>(k);
  result.dfResidual = static_cast<long long>(n - p);

  const std::vector<double> beta = solveLinearSystem(xtx, xty);
  const std::vector<std::vector<double>> xtxInv = invertMatrix(xtx);

  result.coefficients = beta;
  result.fitted.resize(n, 0.0);
  result.residuals.resize(n, 0.0);
  double sse = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double fitted = 0.0;
    for (std::size_t j = 0; j < p; ++j) fitted += A[i][j] * beta[j];
    result.fitted[i] = fitted;
    result.residuals[i] = y[i] - fitted;
    sse += result.residuals[i] * result.residuals[i];
  }
  result.sse = sse;

  double meanY = 0.0;
  if (options.hasIntercept) {
    for (double v : y) meanY += v;
    meanY /= static_cast<double>(n);
  }
  double sst = 0.0;
  for (double v : y) {
    const double dev = options.hasIntercept ? v - meanY : v;
    sst += dev * dev;
  }
  result.sst = sst;
  result.ssr = std::max(0.0, sst - sse);
  result.r2 = sst > 0.0 ? 1.0 - sse / sst : kNaN;
  const double nObs = static_cast<double>(n);
  const double pDbl = static_cast<double>(p);
  result.adjR2 = sst > 0.0
                     ? 1.0 - (1.0 - result.r2) * (nObs - (options.hasIntercept ? 1.0 : 0.0)) /
                                  (nObs - pDbl)
                     : kNaN;
  result.rmse = std::sqrt(sse / nObs);
  const double sigma2 = sse / static_cast<double>(n - p);
  result.residualStdError = std::sqrt(sigma2);

  result.stdErrors.assign(p, kNaN);
  result.tStats.assign(p, kNaN);
  result.pValues.assign(p, kNaN);
  result.ciLow.assign(p, kNaN);
  result.ciHigh.assign(p, kNaN);
  if (result.dfResidual > 0) {
    const double critical = studentTQuantile(1.0 - options.alpha / 2.0,
                                             static_cast<double>(result.dfResidual));
    for (std::size_t j = 0; j < p; ++j) {
      const double variance = sigma2 * xtxInv[j][j];
      if (variance >= 0.0) {
        const double se = std::sqrt(variance);
        result.stdErrors[j] = se;
        // With zero residual variance the t statistic is undefined/infinite;
        // report it as n/a instead of an astronomically large finite value.
        result.tStats[j] = se > 0.0 ? beta[j] / se : kNaN;
        result.pValues[j] = studentTTwoTailedP(result.tStats[j], static_cast<double>(result.dfResidual));
        result.ciLow[j] = beta[j] - critical * se;
        result.ciHigh[j] = beta[j] + critical * se;
      }
    }
  }

  if (result.dfModel > 0) {
    const double f = (result.ssr / static_cast<double>(result.dfModel)) / sigma2;
    result.fStat = f;
    result.fPValue = incompleteBeta(static_cast<double>(result.dfResidual) / 2.0,
                                    static_cast<double>(result.dfModel) / 2.0,
                                    static_cast<double>(result.dfResidual) /
                                        (static_cast<double>(result.dfResidual) +
                                         static_cast<double>(result.dfModel) * f));
  } else {
    result.fStat = kNaN;
    result.fPValue = kNaN;
  }

  // Durbin-Watson on the residual sequence.
  if (n >= 3) {
    double num = 0.0;
    double den = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      den += result.residuals[i] * result.residuals[i];
      if (i > 0) {
        const double diff = result.residuals[i] - result.residuals[i - 1];
        num += diff * diff;
      }
    }
    result.durbinWatson = den > 0.0 ? num / den : kNaN;
  } else {
    result.durbinWatson = kNaN;
  }

  result.valid = true;
  return result;
}

std::vector<double> predict(const RegressionResult& model,
                            const std::vector<std::vector<double>>& x) {
  const std::size_t k = model.predictorNames.size();
  std::vector<double> out;
  out.reserve(x.size());
  for (const std::vector<double>& row : x) {
    if (row.size() != k) {
      throw std::invalid_argument("prediction row has " + std::to_string(row.size()) +
                                  " values, expected " + std::to_string(k));
    }
    double value = model.hasIntercept ? model.coefficients[0] : 0.0;
    const std::size_t offset = model.hasIntercept ? 1u : 0u;
    for (std::size_t j = 0; j < k; ++j) value += model.coefficients[offset + j] * row[j];
    out.push_back(value);
  }
  return out;
}

std::string renderRegression(const RegressionResult& result, bool unicode, bool color) {
  std::ostringstream summary;
  std::vector<std::pair<std::string, std::string>> items;
  items.emplace_back("target", result.targetName);
  items.emplace_back("observations", std::to_string(result.n));
  items.emplace_back("parameters", std::to_string(result.termNames.size()) +
                                       (result.hasIntercept ? " (with intercept)" : " (no intercept)"));
  items.emplace_back("R²", fmtFixed(result.r2, 6));
  items.emplace_back("adjusted R²", fmtFixed(result.adjR2, 6));
  if (std::isfinite(result.fStat)) {
    items.emplace_back("F-statistic", fmtFixed(result.fStat, 4) + " on " +
                                          std::to_string(result.dfModel) + " and " +
                                          std::to_string(result.dfResidual) + " DF, p=" +
                                          fmtFixed(result.fPValue, 6));
  }
  items.emplace_back("residual SE", fmtFixed(result.residualStdError, 6));
  items.emplace_back("RMSE", fmtFixed(result.rmse, 6));
  items.emplace_back("Durbin-Watson",
                     std::isfinite(result.durbinWatson) ? fmtFixed(result.durbinWatson, 4)
                                                        : std::string("n/a"));
  summary << keyValueBlock(items, color) << "\n";

  Table table;
  table.setTitle("Coefficients");
  table.addColumn("term", Align::Left);
  table.addColumn("coefficient", Align::Right);
  table.addColumn("std error", Align::Right);
  table.addColumn("t", Align::Right);
  table.addColumn("p-value", Align::Right);
  table.addColumn("95% CI", Align::Right);
  for (std::size_t j = 0; j < result.termNames.size(); ++j) {
    std::vector<std::string> row{result.termNames[j]};
    row.push_back(fmtFixed(result.coefficients[j], 6));
    row.push_back(std::isfinite(result.stdErrors[j]) ? fmtFixed(result.stdErrors[j], 6)
                                                     : std::string("n/a"));
    row.push_back(std::isfinite(result.tStats[j]) ? fmtFixed(result.tStats[j], 3)
                                                  : std::string("n/a"));
    row.push_back(std::isfinite(result.pValues[j]) ? fmtFixed(result.pValues[j], 6)
                                                   : std::string("n/a"));
    if (std::isfinite(result.ciLow[j]) && std::isfinite(result.ciHigh[j])) {
      row.push_back("[" + fmtFixed(result.ciLow[j], 4) + ", " + fmtFixed(result.ciHigh[j], 4) +
                    "]");
    } else {
      row.push_back("n/a");
    }
    table.addRow(row);
  }
  table.setFooter("t and p from two-tailed Student t with " +
                  std::to_string(result.dfResidual) + " residual DF; R² " +
                  (result.hasIntercept ? "uses centered" : "uses uncentered") + " SST");
  return summary.str() + table.render(unicode);
}

void regressionToJson(json::Writer& writer, const RegressionResult& result) {
  writer.beginObject();
  writer.kv("target", result.targetName);
  writer.kv("n", result.n);
  writer.kvBool("has_intercept", result.hasIntercept);
  writer.kv("r2", result.r2);
  writer.kv("adj_r2", result.adjR2);
  writer.kv("f_stat", result.fStat);
  writer.kv("f_pvalue", result.fPValue);
  writer.kv("sse", result.sse);
  writer.kv("rmse", result.rmse);
  writer.kv("residual_se", result.residualStdError);
  writer.kv("durbin_watson", result.durbinWatson);
  writer.kv("df_model", result.dfModel);
  writer.kv("df_residual", result.dfResidual);
  writer.key("coefficients");
  writer.beginArray();
  for (std::size_t j = 0; j < result.termNames.size(); ++j) {
    writer.beginObject();
    writer.kv("term", result.termNames[j]);
    writer.kv("coefficient", result.coefficients[j]);
    writer.kv("std_error", result.stdErrors[j]);
    writer.kv("t", result.tStats[j]);
    writer.kv("p", result.pValues[j]);
    writer.key("ci");
    writer.beginArray();
    writer.value(result.ciLow[j]);
    writer.value(result.ciHigh[j]);
    writer.endArray();
    writer.endObject();
  }
  writer.endArray();
  // NOTE: the top-level object is intentionally left open so callers can
  // append extra keys (e.g. "predictions"); the caller must call endObject().
}

}  // namespace se
