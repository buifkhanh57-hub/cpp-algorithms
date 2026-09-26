// statsengine — csv.cpp
//
// RFC 4180 parser and writer implementation. See csv.hpp for the contract.
//
// Parsing notes:
//   * a double quote may only appear at the start of a field (opening a
//     quoted field) or doubled inside a quoted field (escape); anything
//     else is a hard error with the exact physical line number
//   * inside quotes, CRLF and lone CR are normalised to LF
//   * a record is terminated by CRLF, LF or a lone CR
//   * completely blank lines produce no record
#include "statsengine/csv.hpp"

#include <algorithm>
#include <cctype>
#include <utility>

#include "statsengine/util.hpp"

namespace se {

std::string columnKindName(ColumnKind kind) {
  switch (kind) {
    case ColumnKind::Numeric:
      return "numeric";
    case ColumnKind::Date:
      return "date";
    case ColumnKind::Categorical:
      break;
  }
  return "categorical";
}

CsvOptions CsvOptions::withDelimiter(char delimiter) {
  CsvOptions options;
  options.delimiter = delimiter;
  return options;
}

CsvParseError::CsvParseError(int line, const std::string& message)
    : std::runtime_error("CSV parse error at line " + std::to_string(line) + ": " + message),
      line(line) {}

CsvReader::CsvReader(std::istream& in, CsvOptions options)
    : in_(in), options_(std::move(options)) {}

CsvDocument CsvReader::read() {
  std::vector<std::vector<std::string>> records;
  std::vector<int> recordStartLines;
  std::vector<std::string> current;
  std::string field;
  int lineNo = 1;
  int rowStartLine = 1;
  int quoteStartLine = -1;
  bool inQuotes = false;
  bool fieldQuoted = false;
  bool rowOpen = false;

  auto pushField = [&]() {
    if (!rowOpen) rowStartLine = lineNo;
    current.push_back(field);
    field.clear();
    fieldQuoted = false;
    rowOpen = true;
  };

  auto finishRow = [&]() {
    if (!rowOpen && field.empty() && current.empty()) return;  // blank line
    current.push_back(field);
    field.clear();
    fieldQuoted = false;
    records.push_back(current);
    recordStartLines.push_back(rowStartLine);
    current.clear();
    rowOpen = false;
  };

  char c = 0;
  while (in_.get(c)) {
    if (inQuotes) {
      if (c == '"') {
        if (in_.peek() == '"') {  // escaped quote
          field += '"';
          in_.get();
        } else {
          inQuotes = false;  // closing quote
        }
      } else if (c == '\r') {
        if (in_.peek() == '\n') in_.get();
        field += '\n';  // normalise CRLF inside quoted fields
        ++lineNo;
      } else if (c == '\n') {
        field += '\n';
        ++lineNo;
      } else {
        field += c;
      }
      continue;
    }

    if (c == '"') {
      if (field.empty() && !fieldQuoted) {
        if (!rowOpen) rowStartLine = lineNo;
        inQuotes = true;
        fieldQuoted = true;
        quoteStartLine = lineNo;
        rowOpen = true;
      } else {
        throw CsvParseError(lineNo, "unexpected '\"' inside unquoted field");
      }
    } else if (c == options_.delimiter) {
      pushField();
    } else if (c == '\r' || c == '\n') {
      if (c == '\r' && in_.peek() == '\n') in_.get();
      finishRow();
      ++lineNo;
    } else {
      if (!rowOpen) rowStartLine = lineNo;
      field += c;
      rowOpen = true;
    }
  }

  if (inQuotes) {
    throw CsvParseError(quoteStartLine, "unterminated quoted field");
  }
  finishRow();

  if (options_.strictRows && !records.empty()) {
    const std::size_t expected = records.front().size();
    for (std::size_t i = 1; i < records.size(); ++i) {
      if (records[i].size() != expected) {
        throw CsvParseError(recordStartLines[i],
                            "row has " + std::to_string(records[i].size()) +
                                " fields, expected " + std::to_string(expected));
      }
    }
  }

  if (!records.empty() && !records.front().empty()) {
    std::string& firstCell = records.front().front();
    if (firstCell.size() >= 3 && static_cast<unsigned char>(firstCell[0]) == 0xEF &&
        static_cast<unsigned char>(firstCell[1]) == 0xBB &&
        static_cast<unsigned char>(firstCell[2]) == 0xBF) {
      firstCell.erase(0, 3);  // strip UTF-8 BOM
    }
  }

  CsvDocument doc;
  if (options_.hasHeader && !records.empty()) {
    doc.header = records.front();
    doc.rows.assign(records.begin() + 1, records.end());
  } else {
    doc.rows = records;
  }
  return doc;
}

CsvWriter::CsvWriter(std::ostream& out, char delimiter) : out_(out), delimiter_(delimiter) {}

std::string CsvWriter::escapeField(const std::string& field, char delimiter) {
  const bool needsQuotes =
      field.find(delimiter) != std::string::npos ||
      field.find('"') != std::string::npos ||
      field.find('\n') != std::string::npos ||
      field.find('\r') != std::string::npos ||
      (!field.empty() && (field.front() == ' ' || field.back() == ' '));
  if (!needsQuotes) return field;
  std::string out;
  out.reserve(field.size() + 8);
  out.push_back('"');
  for (char c : field) {
    if (c == '"') out.push_back('"');
    out.push_back(c);
  }
  out.push_back('"');
  return out;
}

void CsvWriter::row(const std::vector<std::string>& fields) {
  if (!firstRow_) out_ << '\n';
  firstRow_ = false;
  for (std::size_t i = 0; i < fields.size(); ++i) {
    if (i > 0) out_ << delimiter_;
    out_ << escapeField(fields[i], delimiter_);
  }
}

bool looksLikeNumber(const std::string& token) {
  double ignored = 0.0;
  return parseNumberStrict(token, ignored);
}

bool looksLikeDate(const std::string& token) {
  int y = 0, m = 0, d = 0;
  return parseDate(token, y, m, d);
}

ColumnKind inferColumnKind(const std::vector<std::string>& tokens,
                           const std::vector<std::string>& missingTokens) {
  std::size_t numericCount = 0;
  std::size_t dateCount = 0;
  std::size_t total = 0;
  for (const std::string& token : tokens) {
    const std::string text = trim(token);
    if (text.empty()) continue;
    if (isMissingToken(text, missingTokens)) continue;  // missing cells carry no type signal
    ++total;
    if (looksLikeNumber(text)) ++numericCount;
    if (looksLikeDate(text)) ++dateCount;
  }
  if (total == 0) return ColumnKind::Categorical;
  if (numericCount == total) return ColumnKind::Numeric;
  if (dateCount * 5 >= total * 4) return ColumnKind::Date;  // >= 80 %
  return ColumnKind::Categorical;
}

}  // namespace se
