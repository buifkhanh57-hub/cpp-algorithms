// statsengine — dataframe.cpp
//
// Column conversion and DataFrame operations. See dataframe.hpp.
#include "statsengine/dataframe.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <utility>

#include "statsengine/util.hpp"

namespace se {

ColumnError::ColumnError(const std::string& name, const std::string& reason)
    : std::runtime_error(reason.empty() ? "unknown column '" + name + "'"
                                        : reason + " ('" + name + "')") {}

// ----------------------------------------------------------------- Column --

Column::Column(std::string name, std::vector<std::string> raw, ColumnKind kind,
               const std::vector<std::string>& missingTokens)
    : name_(std::move(name)), kind_(kind), raw_(std::move(raw)) {
  const std::size_t n = raw_.size();
  values_.assign(n, kNaN);
  missing_.assign(n, 0);

  bool allInteger = n > 0;
  int maxDecimals = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::string token = trim(raw_[i]);
    if (isMissingToken(token, missingTokens)) {
      missing_[i] = 1;
      ++missingCount_;
      allInteger = false;
      continue;
    }
    switch (kind_) {
      case ColumnKind::Numeric: {
        double parsed = 0.0;
        if (parseNumberStrict(token, parsed)) {
          values_[i] = parsed;
          long long ignored = 0;
          if (!parseIntStrict(token, ignored)) allInteger = false;
          const std::size_t dot = token.find('.');
          if (dot != std::string::npos) {
            int decimals = 0;
            for (std::size_t j = dot + 1; j < token.size(); ++j) {
              if (!std::isdigit(static_cast<unsigned char>(token[j]))) break;
              ++decimals;
            }
            maxDecimals = std::max(maxDecimals, decimals);
          }
        } else {
          missing_[i] = 1;  // unparseable cell in a numeric column
          ++missingCount_;
          allInteger = false;
        }
        break;
      }
      case ColumnKind::Date: {
        int y = 0, m = 0, d = 0;
        if (parseDate(token, y, m, d)) {
          values_[i] = static_cast<double>(dateToSerial(y, m, d));
        } else {
          missing_[i] = 1;
          ++missingCount_;
        }
        break;
      }
      case ColumnKind::Categorical:
        break;  // values_ stays NaN
    }
  }
  integerLike_ = (kind_ == ColumnKind::Numeric) && allInteger;
  maxDecimals_ = std::min(maxDecimals, 10);
}

std::vector<double> Column::numericCompact() const {
  std::vector<double> out;
  out.reserve(values_.size());
  for (std::size_t i = 0; i < values_.size(); ++i) {
    if (missing_[i] == 0) out.push_back(values_[i]);
  }
  return out;
}

std::vector<std::string> Column::categories() const {
  std::set<std::string> unique;
  for (std::size_t i = 0; i < raw_.size(); ++i) {
    if (missing_[i] == 0) unique.insert(trim(raw_[i]));
  }
  return std::vector<std::string>(unique.begin(), unique.end());
}

Column Column::takeRows(const std::vector<std::size_t>& indices) const {
  Column out;
  out.name_ = name_;
  out.kind_ = kind_;
  out.raw_.reserve(indices.size());
  out.values_.reserve(indices.size());
  out.missing_.reserve(indices.size());
  for (std::size_t idx : indices) {
    if (idx >= raw_.size()) throw std::out_of_range("Column::takeRows: row index out of range");
    out.raw_.push_back(raw_[idx]);
    out.values_.push_back(values_[idx]);
    out.missing_.push_back(missing_[idx]);
    if (missing_[idx] != 0) ++out.missingCount_;
  }
  out.integerLike_ = integerLike_ && !indices.empty();
  out.maxDecimals_ = maxDecimals_;
  return out;
}

// -------------------------------------------------------------- DataFrame --

DataFrame DataFrame::fromCsv(const CsvDocument& doc, const CsvOptions& options) {
  DataFrame frame;
  std::size_t width = 0;
  for (const std::vector<std::string>& row : doc.rows) width = std::max(width, row.size());

  std::vector<std::string> names = doc.header;
  if (!names.empty()) {
    width = std::max(width, names.size());
  } else {
    names.reserve(width);
    for (std::size_t j = 0; j < width; ++j) names.push_back("c" + std::to_string(j + 1));
  }
  if (width == 0) return frame;
  names.resize(width);

  // Normalise header names: trim, synthesise empties, de-duplicate.
  std::map<std::string, int> seen;
  for (std::size_t j = 0; j < names.size(); ++j) {
    std::string name = trim(names[j]);
    if (name.empty()) name = "c" + std::to_string(j + 1);
    const auto it = seen.find(name);
    if (it != seen.end()) {
      ++it->second;
      name += "_" + std::to_string(it->second);
    } else {
      seen[name] = 1;
    }
    names[j] = name;
  }

  frame.columns_.reserve(width);
  for (std::size_t j = 0; j < width; ++j) {
    std::vector<std::string> raw;
    raw.reserve(doc.rows.size());
    for (const std::vector<std::string>& row : doc.rows) {
      raw.push_back(j < row.size() ? row[j] : std::string());
    }
    const ColumnKind kind = inferColumnKind(raw, options.missingTokens);
    frame.columns_.emplace_back(names[j], std::move(raw), kind, options.missingTokens);
  }
  frame.rows_ = doc.rows.size();
  return frame;
}

const Column& DataFrame::column(const std::string& name) const {
  for (const Column& col : columns_) {
    if (col.name() == name) return col;
  }
  std::vector<std::string> available;
  available.reserve(columns_.size());
  for (const Column& col : columns_) available.push_back(col.name());
  throw ColumnError(name, "unknown column '" + name + "' (available: " + joinList(available) + ")");
}

bool DataFrame::hasColumn(const std::string& name) const noexcept {
  for (const Column& col : columns_) {
    if (col.name() == name) return true;
  }
  return false;
}

std::vector<std::string> DataFrame::columnNames() const {
  std::vector<std::string> names;
  names.reserve(columns_.size());
  for (const Column& col : columns_) names.push_back(col.name());
  return names;
}

std::vector<std::size_t> DataFrame::numericColumnIndices() const noexcept {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < columns_.size(); ++i) {
    if (columns_[i].kind() == ColumnKind::Numeric) out.push_back(i);
  }
  return out;
}

std::vector<std::size_t> DataFrame::resolveColumns(const std::string& commaList,
                                                   bool allowAll) const {
  const std::vector<std::string> names = splitList(commaList);
  if (names.empty()) {
    if (allowAll) {
      std::vector<std::size_t> all(columns_.size());
      for (std::size_t i = 0; i < columns_.size(); ++i) all[i] = i;
      return all;
    }
    throw ArgsError("no columns selected; use --column <name> or --columns <a,b,c>");
  }
  std::vector<std::size_t> out;
  out.reserve(names.size());
  for (const std::string& rawName : names) {
    const std::string name = trim(rawName);
    bool matched = false;
    for (std::size_t i = 0; i < columns_.size(); ++i) {
      if (columns_[i].name() == name) {
        out.push_back(i);
        matched = true;
        break;
      }
    }
    if (matched) continue;
    // Fall back to 1-based numeric index.
    const bool allDigits = !name.empty() &&
                           name.find_first_not_of("0123456789") == std::string::npos;
    if (allDigits) {
      const long long oneBased = std::stoll(name);
      if (oneBased < 1 || static_cast<std::size_t>(oneBased) > columns_.size()) {
        throw ColumnError(name, "column index out of range (file has " +
                                    std::to_string(columns_.size()) + " columns)");
      }
      out.push_back(static_cast<std::size_t>(oneBased - 1));
      continue;
    }
    std::vector<std::string> available;
    available.reserve(columns_.size());
    for (const Column& col : columns_) available.push_back(col.name());
    throw ColumnError(name,
                      "unknown column '" + name + "' (available: " + joinList(available) + ")");
  }
  return out;
}

DataFrame DataFrame::selectColumns(const std::vector<std::size_t>& indices) const {
  DataFrame out;
  out.rows_ = rows_;
  out.columns_.reserve(indices.size());
  for (std::size_t idx : indices) out.columns_.push_back(columns_.at(idx));
  return out;
}

DataFrame DataFrame::takeRows(const std::vector<std::size_t>& indices) const {
  DataFrame out;
  out.rows_ = indices.size();
  out.columns_.reserve(columns_.size());
  for (const Column& col : columns_) out.columns_.push_back(col.takeRows(indices));
  return out;
}

std::vector<std::size_t> DataFrame::completeRows(const std::vector<std::string>& names) const {
  std::vector<const Column*> cols;
  cols.reserve(names.size());
  for (const std::string& name : names) cols.push_back(&column(name));

  std::vector<std::size_t> keep;
  keep.reserve(rows_);
  for (std::size_t i = 0; i < rows_; ++i) {
    bool complete = true;
    for (const Column* col : cols) {
      if (col->isMissing(i)) {
        complete = false;
        break;
      }
    }
    if (complete) keep.push_back(i);
  }
  return keep;
}

std::size_t DataFrame::missingCells() const noexcept {
  std::size_t total = 0;
  for (const Column& col : columns_) total += col.missingCount();
  return total;
}

std::size_t DataFrame::bytesUsed() const noexcept {
  std::size_t total = 0;
  for (const Column& col : columns_) {
    for (const std::string& cell : col.raw()) total += cell.capacity();
    total += col.size() * (sizeof(double) + 1);
  }
  return total;
}

std::vector<std::string> DataFrame::rowAt(std::size_t i) const {
  std::vector<std::string> cells;
  cells.reserve(columns_.size());
  for (const Column& col : columns_) cells.push_back(col.rawAt(i));
  return cells;
}

}  // namespace se
