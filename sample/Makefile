CXX = clang++
CXXFLAGS = -std=c++17 -Wall 
LDFLAGS = 

LIB_DIR = ./lib
SRC_DIR = ./src
OBJ_DIR = ./build
INCLUDE = -I ../include

CXX_SOURCES = $(shell ls $(SRC_DIR)/*.cc)
OBJECTS = $(subst $(SRC_DIR), $(OBJ_DIR), $(CXX_SOURCES:.cc=.o))
SRC_ = $(wildcard $(SRC_DIR)/*.h)
DEPENDS = $(OBJECTS:.o=.d)

TARGET = sample

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) 

$(OBJ_DIR)/%.o: ./$(SRC_DIR)/%.cc $(SRC_H)
	@if [ ! -d $(OBJ_DIR) ]; \
		then echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
		fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

run: $(TARGET)
	./$<

clean:
	$(RM) $(OBJECTS) $(TARGET) $(DEPENDS)

.PHONY: all clean

-include $(DEPENDS)