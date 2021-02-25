## INSTRUCTIONS ##

# A tool for fixed point iterations #

To compile: 
- if you are using the *mk module system*, you have to run *module load eigen*.
- if not, you have to specify the location of the Eigen headers manually in the
apposite variable in the *Makefile* in this directory;

run *make*, then *make run* to run the examples.

Requirements: a compiler supporting the standard C++17, as well as the *doxygen* and *graphviz*
packages for creating the documentation with *make doc*. 
To consult documentation (which is suggested for
approaching the code), open *index.html* in the generated *html* directory.

Many other targets are present in the *Makefile* to build specific examples or provide other tools (check it out if you need).

The example main_linear_problem is not mentioned in the report but contains interesting stuff anyways.

possible TO DO list:
- create an object factory and allow for dynamic loading
- remove the std::deque and rethink the output printing
- generalizing and improving the testing parameter (for example, instead of using a line number for sed,
you should first search the line that matches it, or maybe use *sed* in a different way?). C++ regular expressions
also are ok, but I don't want to parse a file and output a new file each time, better let *sed* do it.
- write in the report the good observation that beta=1 is the Anderson without relaxation (good results in the linear case)
