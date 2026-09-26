// statsengine — main.cpp
//
// CLI entry point: global parsing, help/version handling, subcommand
// dispatch and exit-code mapping.
//
// Exit codes:
//   0  success
//   1  runtime error (I/O failure, malformed CSV, statistical failure)
//   2  usage error (unknown command/option, bad or wrong-type column)
#include "statsengine/args.hpp"
#include "statsengine/commands.hpp"
#include "statsengine/colors.hpp"
#include "statsengine/csv.hpp"
#include "statsengine/dataframe.hpp"
#include "statsengine/util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace se {
namespace {

int levenshtein(const std::string& a, const std::string& b) {
  std::vector<std::size_t> prev(b.size() + 1);
  std::vector<std::size_t> curr(b.size() + 1);
  for (std::size_t j = 0; j <= b.size(); ++j) prev[j] = j;
  for (std::size_t i = 1; i <= a.size(); ++i) {
    curr[0] = i;
    for (std::size_t j = 1; j <= b.size(); ++j) {
      const std::size_t substitution = prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
      curr[j] = std::min({prev[j] + 1, curr[j - 1] + 1, substitution});
    }
    std::swap(prev, curr);
  }
  return static_cast<int>(prev[b.size()]);
}

std::string suggestCommand(const std::string& input) {
  const std::string needle = toLower(input);
  const std::vector<CommandInfo>& registry = commandRegistry();
  const CommandInfo* best = nullptr;
  int bestDistance = 1 << 30;
  for (const CommandInfo& info : registry) {
    const int distance = levenshtein(needle, info.name);
    if (distance < bestDistance) {
      bestDistance = distance;
      best = &info;
    }
  }
  if (best == nullptr || bestDistance > 3) return "";
  return best->name;
}

int run(const Args& args) {
  if (args.versionRequested()) {
    std::cout << "statsengine " << kStatsengineVersion << "\n";
    return 0;
  }
  if (args.helpRequested()) {
    const CommandInfo* command = args.command().empty() ? nullptr : findCommand(args.command());
    if (command != nullptr) {
      printCommandHelp(std::cout, *command, colors::enabled());
    } else {
      printHelp(std::cout, colors::enabled());
    }
    return 0;
  }
  if (args.command().empty()) {
    printHelp(std::cout, colors::enabled());
    return 0;
  }

  const CommandInfo* command = findCommand(args.command());
  if (command == nullptr) {
    const std::string suggestion = suggestCommand(args.command());
    std::cerr << colors::red("error: ") << "unknown command '" << args.command() << "'";
    if (!suggestion.empty()) std::cerr << " (did you mean '" << suggestion << "'?)";
    std::cerr << "\nTry 'statsengine --help'.\n";
    return 2;
  }
  return command->run(args);
}

}  // namespace
}  // namespace se

int main(int argc, char** argv) {
  using namespace se;
  try {
    Args args(argc, argv);
    colors::init(args.noColorRequested(), args.forceColor());
    return run(args);
  } catch (const ArgsError& e) {
    std::cerr << colors::red("error: ") << e.what() << "\n";
    return 2;
  } catch (const ColumnError& e) {
    std::cerr << colors::red("error: ") << e.what() << "\n";
    return 2;
  } catch (const std::invalid_argument& e) {
    std::cerr << colors::red("error: ") << e.what() << "\n";
    return 2;
  } catch (const CsvParseError& e) {
    std::cerr << colors::red("error: ") << e.what() << "\n";
    return 1;
  } catch (const std::exception& e) {
    std::cerr << colors::red("error: ") << e.what() << "\n";
    return 1;
  }
}
