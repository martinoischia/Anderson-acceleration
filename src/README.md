## INSTRUCTIONS ##

# A tool for fixed point iterations #

To compile: 
- if you are using the *mk module system*, you have to run *module load eigen*.
- if not, you have to specify the location of the Eigen headers manually in the
apposite variable in the *Makefile* in this directory;

run *make*, then *make run* to run the examples.

Other useful targets: *make clean* and *make doc*

Requirements: a compiler supporting the standard C++17, as well as the *doxygen* and *graphviz*
packages for creating the documentation.