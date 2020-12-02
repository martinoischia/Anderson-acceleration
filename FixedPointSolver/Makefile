SRCS=$(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS=$(wildcard *.hpp)
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
EIGEN_DIR=$(mkEigenInc)
STANDARD=c++17
#==========================

LINK.o = $(CXX) $(LDFLAGS)
CPPFLAGS=-I. -I$(EIGEN_DIR)
CXXFLAGS=-std=$(STANDARD)
LDFLAGS=-std=$(STANDARD)
LDLIBS=

#==========================

.phony= all clean distclean doc

.DEFAULT_GOAL = all

all: make.dep $(EXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS)

#distclean: clean

# doc:
	# doxygen $(DOXYFILE)

$(EXEC): $(OBJS)

$(OBJS): $(SRCS)

make.dep: $(SRCS)
	$(RM) make.dep
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> make.dep; \
	done

-include make.dep
