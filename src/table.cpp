// statsengine — table.cpp
//
// Table rendering with Unicode/ASCII borders and ANSI-aware padding.
// See table.hpp.
#include "statsengine/table.hpp"

#include <cmath>
#include <ostream>
#include <sstream>

#include "statsengine/colors.hpp"

namespace se {

namespace {

std::string repeat(const std::string& glyph, std::size_t count) {
  std::string out;
  out.reserve(glyph.size() * count);
  for (std::size_t i = 0; i < count; ++i) out += glyph;
  return out;
}

}  // namespace

Table& Table::setTitle(const std::string& title) {
  title_ = title;
  return *this;
}

Table& Table::addColumn(const std::string& name, Align align) {
  headers_.push_back(name);
  aligns_.push_back(align);
  return *this;
}

Table& Table::addRow(const std::vector<std::string>& cells) {
  rows_.push_back(cells);
  return *this;
}

Table& Table::setFooter(const std::string& footer) {
  footer_ = footer;
  return *this;
}

Table& Table::setNotes(const std::vector<std::string>& notes) {
  notes_ = notes;
  return *this;
}

std::vector<std::size_t> Table::columnWidths() const {
  std::vector<std::size_t> widths(headers_.size(), 0);
  for (std::size_t i = 0; i < headers_.size(); ++i) {
    widths[i] = colors::visibleWidth(headers_[i]);
  }
  for (const std::vector<std::string>& row : rows_) {
    for (std::size_t i = 0; i < widths.size() && i < row.size(); ++i) {
      widths[i] = std::max(widths[i], colors::visibleWidth(row[i]));
    }
  }
  return widths;
}

std::string Table::padCell(const std::string& cell, std::size_t width, Align align) const {
  const std::size_t textWidth = colors::visibleWidth(cell);
  if (textWidth >= width) return cell;
  const std::size_t padding = width - textWidth;
  switch (align) {
    case Align::Right:
      return std::string(padding, ' ') + cell;
    case Align::Center: {
      const std::size_t left = padding / 2;
      return std::string(left, ' ') + cell + std::string(padding - left, ' ');
    }
    case Align::Left:
      break;
  }
  return cell + std::string(padding, ' ');
}

std::string Table::horizontalLine(const char* left, const char* middle, const char* right,
                                  bool unicode) const {
  const std::vector<std::size_t> widths = columnWidths();
  const std::string glyph = unicode ? "─" : "-";
  std::string out = left;
  for (std::size_t i = 0; i < widths.size(); ++i) {
    if (i > 0) out += middle;
    out += repeat(glyph, widths[i] + 2);
  }
  out += right;
  return out;
}

std::string Table::borderRow(const std::vector<std::string>& cells, bool unicode) const {
  const std::vector<std::size_t> widths = columnWidths();
  const std::string vert = unicode ? "│" : "|";
  std::string out = std::string(" ") + vert;
  for (std::size_t i = 0; i < widths.size(); ++i) {
    out += " ";
    out += padCell(i < cells.size() ? cells[i] : std::string(), widths[i],
                   i < aligns_.size() ? aligns_[i] : Align::Left);
    out += " ";
    out += vert;
  }
  return out;
}

std::string Table::render(bool unicode) const {
  const std::vector<std::size_t> widths = columnWidths();
  if (widths.empty()) return title_.empty() ? std::string() : title_ + "\n";

  const char* topLeft = unicode ? "┌" : "+";
  const char* topMid = unicode ? "┬" : "+";
  const char* topRight = unicode ? "┐" : "+";
  const char* midLeft = unicode ? "├" : "+";
  const char* midMid = unicode ? "┼" : "+";
  const char* midRight = unicode ? "┤" : "+";
  const char* botLeft = unicode ? "└" : "+";
  const char* botMid = unicode ? "┴" : "+";
  const char* botRight = unicode ? "┘" : "+";

  std::ostringstream out;
  if (!title_.empty()) {
    out << colors::bold(title_) << "\n";
  }
  out << colors::dim(horizontalLine(topLeft, topMid, topRight, unicode)) << "\n";
  out << borderRow(headers_, unicode) << "\n";
  out << colors::dim(horizontalLine(midLeft, midMid, midRight, unicode)) << "\n";
  for (const std::vector<std::string>& row : rows_) {
    out << borderRow(row, unicode) << "\n";
  }
  out << colors::dim(horizontalLine(botLeft, botMid, botRight, unicode)) << "\n";
  if (!footer_.empty()) {
    out << colors::dim(footer_) << "\n";
  }
  for (const std::string& note : notes_) {
    out << colors::dim("note: " + note) << "\n";
  }
  return out.str();
}

void Table::print(std::ostream& os, bool unicode) const { os << render(unicode); }

std::string keyValueBlock(const std::vector<std::pair<std::string, std::string>>& items,
                          bool color) {
  std::size_t keyWidth = 0;
  for (const auto& item : items) {
    keyWidth = std::max(keyWidth, item.first.size());
  }
  std::ostringstream out;
  for (const auto& item : items) {
    const std::string key =
        color ? colors::bold(item.first) : item.first;
    out << "  " << key << std::string(keyWidth - item.first.size(), ' ') << " : "
        << item.second << "\n";
  }
  return out.str();
}

std::string progressBar(double fraction, int width, bool unicode, bool color) {
  fraction = std::min(1.0, std::max(0.0, fraction));
  if (width < 1) width = 1;
  const int filled = static_cast<int>(std::round(fraction * static_cast<double>(width)));
  std::string bar;
  for (int i = 0; i < width; ++i) {
    if (unicode) {
      bar += i < filled ? "\xe2\x96\x88" : "\xe2\x96\x91";  // █ / ░
    } else {
      bar += i < filled ? '#' : '.';
    }
  }
  if (!color) return bar;
  const char* code = fraction >= 0.9 ? "32" : (fraction >= 0.5 ? "33" : "31");
  return colors::paint(bar, code);
}

}  // namespace se
