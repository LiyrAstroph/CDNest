/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 * MPICH version
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

extern void model2(int argc, char **argv);
extern void model2(int argc, char **argv);

int thistask, totaltask;
int namelen;

int main(int argc, char **argv)
{
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);

  MPI_Barrier(MPI_COMM_WORLD);
  model1(argc, argv);
  
  //MPI_Barrier(MPI_COMM_WORLD);
  //model2(argc, argv);

  MPI_Finalize();
  return 0;
}
