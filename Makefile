# statsengine — Makefile
#
# Targets:
#   make / make build   compile ./build/statsengine
#   make test           build and run the assert-based test suite
#   make run            build and print the top-level help
#   make clean          remove the build directory

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
CPPFLAGS := -Iinclude

BUILD  := build
BIN    := $(BUILD)/statsengine
TESTBIN := $(BUILD)/test_core

SRCS := $(wildcard src/*.cpp)
LIB_SRCS := $(filter-out src/main.cpp,$(SRCS))
OBJS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRCS))
HDRS := $(wildcard include/statsengine/*.hpp)

.PHONY: all build test run clean

all: build

build: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(BUILD)/%.o: src/%.cpp $(HDRS)
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# The test suite links every module except main.cpp (it has its own main).
test: $(BIN)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) tests/test_core.cpp $(LIB_SRCS) -o $(TESTBIN)
	./$(TESTBIN)

run: $(BIN)
	./$(BIN) --help

clean:
	rm -rf $(BUILD)
