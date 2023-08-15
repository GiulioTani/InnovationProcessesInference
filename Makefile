# define the C compiler to use
CC = g++
SH = bash
RM = rm -f
# define any compile-time flags
CFLAGS = -Wall -Ofast -Wextra -std=c++1z # -D_GLIBCXX_DEBUG -D_FORTIFY_SOURCE=2 -pg --coverage -fprofile-abs-path#-fno-ipa-cp-clone -g -Og  #quando faccio il linking di fragprob -lgcov -pg

# define any directories containing header files other than /usr/include
WD=`pwd`
INCLUDES = -I${WD}

#when perfecting the folder structure change this to ../include
HEADERS = lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lm -lpthread
VPATH=bin:src
ODIR=bin

_DEPS = bookprob.hpp bookdkl.hpp supportlib.hpp paramOpt.hpp authorSplitter.hpp datatypes.hpp base_experiment.hpp
DEPS = $(patsubst %,$(HEADERS)/%,$(_DEPS))

_OBJ = bookprob.o bookdkl.o fragprob.o supportlib.o paramOpt.o authorSplitter.o base_experiment.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

CINP = $(wildcard *.c)
COUT = $(patsubst %.c,$(ODIR)/%,$(CINP))
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

all:  structure $(ODIR)/fragprob $(ODIR)/libLZ77dict.so $(ODIR)/libattributor.so
	@echo  All model version compiled

structure:
	@$(SH) structure.sh

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES)
	@echo  $< compiled

$(ODIR)/fragprob: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(INCLUDES)
	@echo fragprob linked

$(ODIR)/libLZ77dict.so: $(ODIR)/LZ77dict_Ctypes.o
	gcc -Wall -Ofast -Wextra -Wno-incompatible-pointer-types -shared -Wl,-soname,LZ77dict.so -o $@ $^ $(INCLUDES)
	@echo $@ shared object created

$(ODIR)/LZ77dict_Ctypes.o: LZ77dict_Ctypes.c
	gcc -c -fPIC -Wall -Ofast -Wextra -Wno-incompatible-pointer-types $^ -o $@ $(INCLUDES)

$(ODIR)/libattributor.so: $(ODIR)/attributor_Ctypes.o
	g++ -Wall -Ofast -Wextra -shared -Wl,-soname,attributor.so -o $@ $^ $(INCLUDES)
	@echo $@ shared object created

$(ODIR)/attributor_Ctypes.o: attributor_Ctypes.cpp
	g++ -c -fPIC -Wall -Ofast -Wextra $^ -o $@ $(INCLUDES)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)


.PHONY: clean

clean:
	$(RM) $(ODIR)/* *~ core $(INCDIR)/*~ 
