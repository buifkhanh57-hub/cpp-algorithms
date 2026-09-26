// statsengine — args.cpp
//
// Argument parser implementation. See args.hpp for the accepted grammar.
#include "statsengine/args.hpp"

#include <cctype>

#include "statsengine/util.hpp"

namespace se {

Args::Args(int argc, char** argv) {
  bool onlyPositionals = false;
  for (int i = 1; i < argc; ++i) {
    const std::string token = argv[i];

    if (onlyPositionals) {
      positionals_.push_back(token);
      continue;
    }
    if (token == "--") {
      onlyPositionals = true;
      continue;
    }
    if (token == "-") {  // stdin by convention
      positionals_.push_back(token);
      continue;
    }

    if (token.size() > 2 && token[0] == '-' && token[1] == '-') {
      std::string body = token.substr(2);
      std::string name = body;
      std::string value;
      bool hasValue = false;
      const std::size_t eq = body.find('=');
      if (eq != std::string::npos) {
        name = body.substr(0, eq);
        value = body.substr(eq + 1);
        hasValue = true;
      }
      if (name.empty()) throw ArgsError("empty option name in '" + token + "'");

      if (name == "help") helpRequested_ = true;
      if (name == "version") versionRequested_ = true;
      if (name == "json") jsonMode_ = true;
      if (name == "ascii") asciiMode_ = true;
      if (name == "no-color") noColorRequested_ = true;
      if (name == "color") forceColor_ = true;

      // Global flags are strictly boolean: they must never swallow the next
      // token, so 'statsengine --json describe ...' keeps 'describe' as the
      // command name instead of turning it into the flag's value.
      const bool isGlobalFlag = name == "help" || name == "version" || name == "json" ||
                                name == "ascii" || name == "no-color" || name == "color";
      if (isGlobalFlag) {
        if (hasValue) {
          throw ArgsError("option --" + name + " does not take a value");
        }
        boolFlags_.insert(name);
        continue;
      }

      if (hasValue) {
        options_[name] = value;
      } else {
        // Consume the next token as the value when it does not look like
        // another option (or when it parses as a number, so negative
        // values such as --mu -1.5 work).
        double numberSink = 0.0;  // parseNumberStrict writes through a reference
        const bool nextIsOption =
            i + 1 < argc && std::string(argv[i + 1]).size() > 1 && argv[i + 1][0] == '-' &&
            !(parseNumberStrict(std::string(argv[i + 1]), numberSink));
        if (i + 1 < argc && !nextIsOption) {
          options_[name] = argv[++i];
        } else {
          boolFlags_.insert(name);
        }
      }
      continue;
    }

    if (token.size() > 1 && token[0] == '-') {
      // Single-dash options are not supported (avoids ambiguity with
      // negative numbers).
      throw ArgsError("unknown option '" + token + "' (options use double dashes, e.g. --file)");
    }

    positionals_.push_back(token);
  }

  command_ = positionals_.empty() ? std::string() : positionals_.front();
}

bool Args::has(const std::string& flag) const {
  return options_.find(flag) != options_.end() || boolFlags_.count(flag) > 0;
}

std::string Args::get(const std::string& flag) const {
  const auto it = options_.find(flag);
  return it == options_.end() ? std::string() : it->second;
}

std::string Args::getOr(const std::string& flag, const std::string& fallback) const {
  const std::string value = get(flag);
  return value.empty() ? fallback : value;
}

long long Args::getInt(const std::string& flag, long long fallback) const {
  if (!has(flag)) return fallback;
  long long parsed = 0;
  if (!parseIntStrict(get(flag), parsed)) {
    throw ArgsError("option --" + flag + " expects an integer, got '" + get(flag) + "'");
  }
  return parsed;
}

double Args::getDouble(const std::string& flag, double fallback) const {
  if (!has(flag)) return fallback;
  double parsed = 0.0;
  if (!parseNumberStrict(get(flag), parsed)) {
    throw ArgsError("option --" + flag + " expects a number, got '" + get(flag) + "'");
  }
  return parsed;
}

std::vector<std::string> Args::getList(const std::string& flag) const {
  return splitList(get(flag));
}

void Args::require(const std::string& flag, const std::string& what) const {
  if (get(flag).empty()) {
    throw ArgsError("missing required option --" + flag + " (" + what + ")");
  }
}

void Args::rejectUnknown(const std::vector<std::string>& allowed) const {
  std::set<std::string> known(allowed.begin(), allowed.end());
  known.insert({"json", "ascii", "no-color", "color", "help", "version"});
  std::vector<std::string> unknown;
  for (const auto& entry : options_) {
    if (known.count(entry.first) == 0) unknown.push_back("--" + entry.first);
  }
  for (const std::string& flag : boolFlags_) {
    if (known.count(flag) == 0) unknown.push_back("--" + flag);
  }
  if (!unknown.empty()) {
    throw ArgsError("unknown option " + joinList(unknown) + " for command '" + command_ +
                    "' (see 'statsengine " + command_ + " --help')");
  }
}

void Args::expectNoPositionals(std::size_t max) const {
  if (positionals_.size() > max + 1) {
    throw ArgsError("unexpected argument '" + positionals_[max + 1] + "'");
  }
}

}  // namespace se
