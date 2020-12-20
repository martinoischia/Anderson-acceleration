# Useful variables
SRCS=$(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS=$(wildcard *.hpp)
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
EIGEN_DIR=$(mkEigenInc)
STANDARD=c++17
OPTFLAGS=-g
# DEFINE=-DTraits=VectorTraits
#==========================


LINK.o = $(CXX) $(LDFLAGS)
CPPFLAGS=-I. -I../Include -I$(EIGEN_DIR) $(DEFINE)
CXXFLAGS=-std=$(STANDARD) $(OPTFLAGS)
LDFLAGS=-std=$(STANDARD) $(OPTFLAGS)
LDLIBS=

#==========================
# Targets

.phony= all clean distclean doc run makedep

.DEFAULT_GOAL = all

all: makedep $(EXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS) make.dep doxygen_sqlite3.db; rm -f -r html latex

#distclean: clean

doc:
	doxygen Doxyfile

$(EXEC): $(OBJS)

$(OBJS): $(SRCS)

run:
	./$(EXEC)

makedep: $(SRCS)
	$(RM) make.dep
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> make.dep; \
	done

-include make.dep
