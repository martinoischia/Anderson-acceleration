# Useful variables
SRCS=$(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
OBJS1 = $(filter-out main_lin%.o, $(OBJS))
OBJS2 = $(filter-out main_sim%.o, $(OBJS))
HEADERS=$(wildcard *.hpp)
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
exec1=main_simple_problem
exec2=main_linear_problem
EIGEN_DIR=$(mkEigenInc)
STANDARD=c++17
OPTFLAGS=-O3
DEFINE=-DNDEBUG
# OPTFLAGS=-g
# DEFINE=-DTraits=VectorTraits
# DEFINE=EigenSolver
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

all: makedep $(exec1) $(exec2)

clean:
	$(RM) -f $(EXEC) $(OBJS) make.dep doxygen_sqlite3.db; rm -f -r html latex

#distclean: clean

doc:
	doxygen Doxyfile

$(exec1): $(OBJS1)

$(exec2): $(OBJS2)

run:
	./main_linear_problem

makedep: $(SRCS)
	$(RM) make.dep
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> make.dep; \
	done

-include make.dep