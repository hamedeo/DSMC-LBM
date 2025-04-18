# Compiler/Flags
SRC = dsmceq.cpp rand.cpp sorter.cpp colider.cpp
OBJ = $(SRC:.cpp=.o)
TAR = dsmcRun
OBJ_DIR = build
$(shell mkdir -p $(OBJ_DIR))

# Rules
all: $(TAR)
$(TAR): $(OBJ:%.o=$(OBJ_DIR)/%.o)
	g++ $(OBJ:%.o=$(OBJ_DIR)/%.o) -o $(TAR)
$(OBJ_DIR)/%.o: %.cpp
	g++ -Wall -g -c $< -o $@
clean:
	rm -f $(OBJ) $(OBJ_DIR)/*.o $(TAR)
.PHONY: all clean
