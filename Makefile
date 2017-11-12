# Compiler
CXX = g++

# Cleanup
RM = rm -f

# Compiling and linking flags
CPPFLAGS = -std=gnu++11
LDFLAGS =
LDLIBS =

# Paths
SRC = src
BIN = bin
TESTS = tests
BUILDSRC = build/src
BUILDTESTS = build/tests
$(shell mkdir -p $(BIN))
$(shell mkdir -p $(BUILDSRC))
$(shell mkdir -p  $(BUILDTESTS))

# Sources
SRCS = $(wildcard $(SRC)/*.cpp)
SRCTESTS = $(wildcard $(TESTS)/*.cpp)
# Objects and library
OBJS = $(subst $(SRC),$(BUILDSRC),$(subst .cpp,.o,$(SRCS)))
OBJTESTS = $(subst $(TESTS),$(BUILDTESTS),$(subst .cpp,.o,$(SRCTESTS)))
LIB = lib/libgraph.a


all: lib tests

# Rule to compile source
$(BUILDSRC)%.o: $(SRC)%.cpp
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

# Rule to create library
lib: $(OBJS)
	ar rcs $(LIB) $(OBJS) $(LDFLAGS) $(LDLIBS)

# Rules to create tests
tests: lib $(OBJTESTS)

$(BUILDTESTS)%.o: $(TESTS)%.cpp
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)
	cp $@ $(subst $(BUILDTESTS),$(BIN),$(subst .o,,$@))

# Rules to cleanup
clean:
	$(RM) $(OBJS) $(OBJTESTS)

cleanlib: clean 
	$(RM) $(LIB)

cleanall: clean cleanlib
	$(RM) $(BIN)/*
