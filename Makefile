CXX      = g++
CXXFLAGS_COMMON = -std=c++11 -Wall -fexceptions

# Change to debug to compile with debugging flags
MODE = release

ifeq ($(MODE),release)
   CXXFLAGS = $(CXXFLAGS_COMMON) -O3 -DNDEBUG
else
   CXXFLAGS = $(CXXFLAGS_COMMON) -g
endif

TARGET = acs

.PHONY: clean all

all: $(TARGET)

$(TARGET): acs.cpp
	$(CXX) $(CXXFLAGS) acs.cpp -o $(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).o
