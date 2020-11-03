SHELL=/bin/bash
CC       ?= mpicc 
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
#LAPACK_INCL = -I /usr/local/share/lapack/include
#LAPACK_LIBS = /usr/local/share/lapack/lib/liblapacke.a -llapack -L/usr/lib64/atlas -lcblas 
#LAPACK_INCL = -I/usr/include/lapacke
#LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
#CBLAS_INCL  = -I/usr/include 
#CBLAS_LIBS  = -L/usr/lib64/atlas -lcblas

MPICHLIB = $(shell pkg-config --libs mpich)
MPIINCL  = $(shell pkg-config --cflags mpich)

#MPICHLIB = -L/usr/local/share/mpich2/lib -lmpich
#MPIINCL  = -I/usr/local/share/mpich2/include

OPTIMIZE    += 
endif

ifeq ($(SYSTEM), "Darwin")
NCORE      :=$(shell sysctl machdep.cpu.core_count | awk '{print $2}')
#GSL_INCL    = -I/opt/local/include
#GSL_LIBS    = -L/opt/local/lib/gsl/lib
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
#LAPACK_INCL = -I /usr/local/share/lapack/include -I/opt/local/include
#LAPACK_LIBS = -framework vecLib -L /usr/local/share/lapack/lib -llapacke -llapack 
#-lcblas 
CBLAS_INCL  =
CBLAS_LIBS  =     
OPTIMIZE    += 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/sharefs/mbh/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/sharefs/mbh/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/sharefs/mbh/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/sharefs/mbh/user/liyanrong/soft/mpich3/include
#LAPACK_INCL = -I/sharefs/mbh/user/liyanrong/soft/lapack/include
#LAPACK_LIBS = -L/sharefs/mbh/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
#CBLAS_INCL  = -I/sharefs/mbh/user/liyanrong/soft/atlas/include
#CBLAS_LIBS  = -L/sharefs/mbh/user/liyanrong/soft/atlas/lib -lcblas
endif

ifeq ($(SYSTEM), "TianheII")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas -lm
MPICHLIB = -lmpich
MPIINCL  =
#LAPACK_INCL = -I/HOME/ihep_yrli_1/BIGDATA/soft/lapack/include
#LAPACK_LIBS = -L/HOME/ihep_yrli_1/BIGDATA/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
endif

EXEC     = tests/dnest
LDN     = libdnest.so

.PHONY: all
all: $(EXEC) $(LDN)

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(MPIINCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(MPICHLIB)


SRC      = src/
INCL     = Makefile $(SRC)/dnestvars.h $(SRC)/model1.h $(SRC)/model2.h $(SRC)/model3.h
 
OBJS = $(SRC)/dnest.o $(SRC)/dnestvars.o $(SRC)/dnestpostprocess.o $(SRC)/model1.o \
       $(SRC)/main.o $(SRC)/model2.o $(SRC)/model3.o

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS) $(LIBS) -o $@
	
$(LDN): $(OBJS)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(LIBS) -fPIC -shared -o libdnest.so $(SRC)/dnest.c $(SRC)/dnestvars.c $(SRC)/dnestpostprocess.c
	#ar rcs libdnest.a dnest.o dnestvars.o
	cp $(SRC)/dnestvars.h .

$(OBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC) libdnest.so *.h
