
# Uncomment the COMPILER and OPTIONS lines corresponding to
# the platform and compiler you plan to use.

# If compiling with the gcc compiler:
COMPILER = gcc
OPTIONS = -O3

# If compiling with the Sun cc compiler on a Sparc machine.
#COMPILER = cc
#OPTIONS = -fast -xO4

# Generic, should work on work on most machines but may not produce
# the most optimal code
#COMPILER = cc
#OPTIONS = -O

# Uncomment the following line to produce code you can debug.
#OPTIONS = -g

# You should not need to change anything below this line
# Type 'make'.  This will build the program 'prest'.
# Type 'make altertest'.  This will build the program 'altertest'.
##################################################################

CC_LIB = -lm

prest:   prest.c 
	$(COMPILER) $(OPTIONS) -o prest prest.c $(CC_LIB)

altertest:   altertest.c 
	$(COMPILER) $(OPTIONS) -o altertest altertest.c $(CC_LIB)

