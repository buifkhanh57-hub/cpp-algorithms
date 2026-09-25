CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra
BUILD := build
TARGET := $(BUILD)/algorithms

all: $(TARGET)

$(TARGET): src/main.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(BUILD)

.PHONY: all run clean
