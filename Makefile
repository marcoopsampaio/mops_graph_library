################
### SETTINGS ###
################

# Compiler
CXX = g++

# Cleanup
RM = rm -f

# Compiling and linking flags
CPPFLAGS = -std=gnu++11
LDFLAGS = 
#LDFLAGS = -DVERBOSE # activate debug messages
#LDFLAGS = -DVERBOSE -DDEBUG # activate verbose and debug messages
LDLIBS =

# Paths
SRC = src
BIN = bin
TESTS = tests
BUILD = build
LIBPATH = lib
$(shell mkdir -p $(BIN))
$(shell mkdir -p $(BUILD))
$(shell mkdir -p $(LIBPATH))

# Sources
SRCS = $(wildcard $(SRC)/*.cpp)
SRCTESTS = $(wildcard $(TESTS)/*.cpp)

# Objects, binaries and library
OBJS = $(subst $(SRC),$(BUILD),$(subst .cpp,.o,$(SRCS)))
BINTESTS = $(subst $(TESTS),$(BIN),$(subst .cpp,,$(SRCTESTS)))
LIBNAME = graph
LIB = $(LIBPATH)/lib$(LIBNAME).a

##############
### RULES ####
##############

all: lib tests

# Rule to compile source
$(BUILD)/%.o: $(SRC)/%.cpp 
	$(CXX) $(CPPFLAGS) -o $@ -c $< $(LDFLAGS)

$(BUILD)/%.o: $(TESTS)/%.cpp
	$(CXX) $(CPPFLAGS) -o $@ -c $< $(LDFLAGS)

# Rule to create library
lib: $(OBJS)
	ar rcs $(LIB) $(OBJS) $(LDLIBS)

# Rules to create tests
tests: $(BINTESTS)

$(BIN)/%: $(BUILD)/%.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(OBJS) $(LDFLAGS) $(LDLIBS)

# Rules to cleanup
clean:
	$(RM) $(OBJS) $(OBJTESTS)

cleanlib: clean 
	$(RM) $(LIB)

cleanall: clean cleanlib
	$(RM) $(BIN)/*
