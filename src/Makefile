# Useful variables
SRCS=$(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS=$(wildcard *.hpp)
exe_sources=$(filter main%.cpp,$(SRCS))
exe_objects=$(filter-out main%.o,$(OBJS))
EXEC=$(exe_sources:.cpp=)
exec1:=main_simple_problem
exec2:=main_linear_problem
exec3:=main_gauss_seidel
exec4:=main_testing_parameter
exec5:=main_nonlinear
OBJS1 = $(exe_objects) $(exec1:=.o) 
OBJS2 = $(exe_objects) $(exec2:=.o) 
OBJS3 = $(exe_objects) $(exec3:=.o) 
OBJS4 = $(exec4:=.o) 
OBJS5 = $(exe_objects) $(exec5:=.o) 


EIGEN_DIR=$(mkEigenInc)

#==========================
#==========================
#Standard options to be fed into the compilation commands

STANDARD=c++17
# WARNING=-Wall

# OPTFLAGS=-O3
# DEFINE=-DNDEBUG

OPTFLAGS=-g


#============ following block should be all commented or all not commented
DEFINE+=-DTraits=DenseTraits
EXEC:= $(filter-out $(exec2),$(EXEC))
exec2:=
#============

# DEFINE+=-DEigenSolver
LINK.o = $(CXX) $(LDFLAGS)
CPPFLAGS=-I. -I../Include -I$(EIGEN_DIR) $(DEFINE) $(WARNING)
CXXFLAGS=-std=$(STANDARD) $(OPTFLAGS) $(WARNING)
LDFLAGS=-std=$(STANDARD) $(OPTFLAGS) $(WARNING)
LDLIBS=
#==========================
#==========================
# Targets

.phony= all clean distclean doc run makedep parameters $(exec1) $(exec2) $(exec3) $(exec4)

.DEFAULT_GOAL = all

all: makedep $(EXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS) make.dep doxygen_sqlite3.db profiling.txt gmon.out; rm -f -r html latex

# distclean: clean

doc:
	doxygen Doxyfile

$(exec1): $(OBJS1)

$(exec2): $(OBJS2)

$(exec3): $(OBJS3)

$(exec4): $(OBJS4) $(exec3)
	$(LINK.o) $(exec4:=.o) -o $(exec4)

$(exec5): $(OBJS5)

run: $(exec1) $(exec2) $(exec3) $(exec5)
	./$(exec1); ./$(exec2); ./$(exec3); ./$(exec5)
	
parameters: $(exec4) $(exec3)
	./$(exec4) > ../documents/Include/heat_map_data.txt

profiler: OPTFLAGS=-pg
profiler: clean $(exec3)
	./main_gauss_seidel; gprof --graph ./main_gauss_seidel gmon.out > profiling.txt

makedep: $(SRCS)
	$(RM) make.dep
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> make.dep; \
	done

-include make.dep