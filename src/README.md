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

Many other targets are present in the *Makefile* to build specific examples (check it out if you need).