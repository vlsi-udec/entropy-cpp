TARGET   := entropy-ng
TARGET2  := entropy-card
TARGET3  := entropy-mem
CXX      := c++
CXXFLAGS := -Wall -Wextra -fopenmp -O3 -march=native -g
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)
LDFLAGS  := -L /usr/lib/pcapplusplus/ -lstdc++ -lm -lPcap++ -lPacket++ -lCommon++ -lpcap -pthread
INCLUDE  := -Isrc/
SRC      := src/entropy.cpp src/pq_array.cpp src/malloc_count.c src/stack_count.c
SRC2     := src/cardinality.cpp
SRCMEM   := src/entropy-new.cpp src/pq_array.cpp src/malloc_count.c src/stack_count.c src/entropy_p3.cpp


OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS2:= $(SRC2:%.cpp=$(OBJ_DIR)/%.o)
OBJECTSMEM:= $(SRCMEM:%.cpp=$(OBJ_DIR)/%.o)

all: build $(TARGET3)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS) $(LDFLAGS) -o $(BUILD)/$(TARGET)
	
$(TARGET2): $(OBJECTS2)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTS2) $(LDFLAGS) -o $(BUILD)/$(TARGET2)

$(TARGET3): $(OBJECTSMEM)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJECTSMEM) $(LDFLAGS) -o $(BUILD)/$(TARGET3)

.PHONY: all build clean debug release

build:
	@mkdir -p $(BUILD)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(BUILD)/*
