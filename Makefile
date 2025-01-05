CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Werror -Weffc++ -Wconversion             \
           -Wsign-conversion -Wshadow -Wformat=2 -Wunused                     \
           -Woverloaded-virtual -Wnon-virtual-dtor -Wold-style-cast           \
           -Wduplicated-cond -Wduplicated-branches -Wlogical-op               \
           -Wnull-dereference -Wuseless-cast                                  \
           -Iinc                                                              \
           -lhdf5 -lhdf5_cpp                                                  \
					 -fopenmp -O0

SRC_DIR = src
INC_DIR = inc
BIN_DIR = bin

# Source files
SRC = $(SRC_DIR)/maxwell.cpp $(SRC_DIR)/save.cpp
MAIN_SRC = $(SRC_DIR)/main.cpp
PLOT_SRC = $(SRC_DIR)/plot.cpp

# Object files
OBJ_SHARED = $(SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)
OBJ_MAIN = $(MAIN_SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o) $(OBJ_SHARED)
OBJ_PLOT = $(PLOT_SRC:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o) $(OBJ_SHARED)

# Targets
TARGET_MAIN = $(BIN_DIR)/main
TARGET_PLOT = $(BIN_DIR)/plot

all: $(TARGET_MAIN) $(TARGET_PLOT)

$(TARGET_MAIN): $(OBJ_MAIN)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -O3 -o $@ $(OBJ_MAIN)

$(TARGET_PLOT): $(OBJ_PLOT)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_PLOT)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/maxwell.hpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean
