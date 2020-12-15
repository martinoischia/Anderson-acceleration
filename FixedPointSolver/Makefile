SRCS=$(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS=$(wildcard *.hpp)
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
EIGEN_DIR=$(mkEigenInc)
STANDARD=c++17
OPTFLAGS=-g
DEFINE=-DTraits=VectorTraits
#==========================


LINK.o = $(CXX) $(LDFLAGS)
CPPFLAGS=-I. -I../Include -I$(EIGEN_DIR) 
CXXFLAGS=-std=$(STANDARD) $(OPTFLAGS)
LDFLAGS=-std=$(STANDARD) $(OPTFLAGS)
LDLIBS=

#==========================

.phony= all clean distclean doc

.DEFAULT_GOAL = all

all: makedep $(EXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS) make.dep

#distclean: clean

# doc:
	# doxygen $(DOXYFILE)

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
