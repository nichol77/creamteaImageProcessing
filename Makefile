# Makefile 
# will make PROGRAM using PROGRAM_CLASSES and R_CLASSES
# 
# R_DICT and R_LIBRARY are created by root 
#
# - no options will create program then run it
# - clean is as usual
# - profile (in theory) will profile the run program - it makes a special version called a.out
# - test will time and run the program 10 times
# - no run will disable the run feature

# These bits to change - file name macros

PROGRAM    = 3DgradScan
PROGRAM2   = edgeDetection

# Compiler macros

GCC       = g++
CPPFLAGS  = -O2 -g -Wall -fPIC $(INCLUDES)  
INCLUDES  = -I$(ROOTSYS)/include 
RLIBS     = -L$(ROOTSYS)/lib/ $(shell $(ROOTSYS)/bin/root-config --libs) -lMinuit -lTreePlayer
LIBRARIES = -lz -lm $(RLIBS)

.PHONY:all

all: $(PROGRAM)

$(PROGRAM) : $(PROGRAM).o
	@echo "****<LINKING : $@>**** "
	$(GCC) $(CPPFLAGS) $(LIBRARIES) $^ -o  $@

$(PROGRAM2) : $(PROGRAM2).o
	@echo "****<LINKING : $@>**** "
	$(GCC) $(CPPFLAGS) $(LIBRARIES) $^ -o  $@

testHistos: testHistos.o
	@echo "****<LINKING : $@>**** "
	$(GCC) $(CPPFLAGS) $(LIBRARIES) $^ -o  $@

%.o:%.cpp
	@echo "****<COMPILING : $@>****"
	$(GCC) $(CPPFLAGS) -c $^ -o $@

.PHONY:clean

clean:
	rm *.o 

.PHONY:root

root:$(PROGRAM)
	root -l analysed.root 


