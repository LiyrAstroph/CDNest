SHELL=/bin/bash

ifndef __MPICC__  # passed from python setup.py
  $(info "CC not defined, using mpicc")
  CC       = mpicc 
endif

OPTIMIZE = -O2 -Wall -finline-functions -fcommon
#OPTIMIZE += -DDebug

#---------target system
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"
#SYSTEM="TianheII"


ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
OPTIMIZE    += 
endif

ifeq ($(SYSTEM), "Darwin")
NCORE      :=$(shell sysctl machdep.cpu.core_count | awk '{print $2}')
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl)   
OPTIMIZE    += 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/sharefs/mbh/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/sharefs/mbh/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/sharefs/mbh/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/sharefs/mbh/user/liyanrong/soft/mpich3/include
endif

ifeq ($(SYSTEM), "TianheII")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas -lm
MPICHLIB = -lmpich
MPIINCL  =
endif

EXEC     = tests/dnest
LDN     = libdnest.so

.PHONY: all
all: $(EXEC) $(LDN)

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(MPIINCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(MPICHLIB)


SRC      = src/
INCL     = Makefile $(SRC)/dnestvars.h $(SRC)/model1.h $(SRC)/model2.h $(SRC)/model3.h \
           $(SRC)/dnest.h
 
OBJS = $(SRC)/dnest.o $(SRC)/dnestvars.o $(SRC)/dnestpostprocess.o $(SRC)/model1.o \
       $(SRC)/main.o $(SRC)/model2.o $(SRC)/model3.o

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS) $(LIBS) -o $@
	
$(LDN): $(OBJS)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(LIBS) -fPIC -shared -o libdnest.so $(SRC)/dnest.c $(SRC)/dnestvars.c $(SRC)/dnestpostprocess.c
	#ar rcs libdnest.a dnest.o dnestvars.o
	cp $(SRC)/dnest.h .

$(OBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC) libdnest.so *.h
