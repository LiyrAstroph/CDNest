SHELL=/bin/bash
CC       = gcc -O3 -Wall

#---------target system
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
#LAPACK_INCL = -I /usr/local/share/lapack/include
#LAPACK_LIBS = /usr/local/share/lapack/lib/liblapacke.a -llapack -L/usr/lib64/atlas -lcblas 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
#CBLAS_INCL  = -I/usr/include 
#CBLAS_LIBS  = -L/usr/lib64/atlas -lcblas

MPICHLIB = -L/usr/local/share/mpich2/lib -lmpich
MPIINCL  = -I/usr/local/share/mpich2/include

OPTIMIZE    = 
endif

ifeq ($(SYSTEM), "Darwin")
NCORE      :=$(shell sysctl machdep.cpu.core_count | awk '{print $2}')
#GSL_INCL    = -I/opt/local/include
#GSL_LIBS    = -L/opt/local/lib/gsl/lib
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I /usr/local/share/lapack/include -I/opt/local/include
LAPACK_LIBS = -framework vecLib -L /usr/local/share/lapack/lib -llapacke -llapack 
#-lcblas 
CBLAS_INCL  =
CBLAS_LIBS  =     
OPTIMIZE    = 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/mbh/mbhd01/soft/gsl/include
GSL_LIBS = -L/mbh/mbhd01/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/mbh/mbhd01/soft/mpich2/lib -lmpich
MPIINCL  = -I/mbh/mbhd01/soft/mpich2/include
LAPACK_INCL = -I/mbh/mbhd01/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
#CBLAS_INCL  = -I/mbh/mbhd01/user/liyanrong/soft/atlas/include
#CBLAS_LIBS  = -L/mbh/mbhd01/user/liyanrong/soft/atlas/lib -lcblas
endif

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(CBLAS_INCL) $(MPIINCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(CBLAS_LIBS) $(MPICHLIB)

EXEC     = dnest
SRC      = ./
INCL     = Makefile $(SRC)/dnestvars.h $(SRC)/model1.h $(SRC)/model2.h
 
OBJS = $(SRC)/dnest.o $(SRC)/dnestvars.o $(SRC)/dnestpostprocess.o $(SRC)/model1.o \
       $(SRC)/main.o $(SRC)/model2.o

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o $@
	$(CC) $(OPTIMIZE) $(LIBS) -fPIC -shared -o libdnest.so $(SRC)/dnest.c $(SRC)/dnestvars.c $(SRC)/dnestpostprocess.c
	#ar rcs libdnest.a dnest.o dnestvars.o

$(OBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC)
