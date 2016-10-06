/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "dnestvars.h"
#include "model3.h"

int which_particle_update;
int which_level_update;
int *perturb_accept;
int num_data_points;
int num_params;
int num_particles;

DataType *data;
void *best_model_thismodel, *best_model_std_thismodel;

void model3(int argc, char **argv)
{ 
  /* setup szie of modeltype, which is used for dnest */
  num_params = 2;
  size_of_modeltype = num_params * sizeof(double);
  best_model_thismodel = malloc(size_of_modeltype);
  best_model_std_thismodel = malloc(size_of_modeltype);
  
  /* setup functions used for dnest*/
  from_prior = from_prior_thismodel3;
  data_load = data_load_thismodel3;
  log_likelihoods_cal = log_likelihoods_cal_thismodel3;
  perturb = perturb_thismodel3;
  print_particle = print_particle_thismodel3;
  copy_model = copy_model_thismodel3;
  create_model = create_model_thismodel3;
  get_num_params = get_num_params_thismodel3;
  copy_best_model = copy_best_model_thismodel3;
  
  
  /* run dnest */
  strcpy(options_file, "OPTIONS3");
  
  if(thistask == 0)
  {
    get_num_particles3(options_file);
  }
  MPI_Bcast(&num_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  perturb_accept = malloc(num_particles * sizeof(int));

  dnest(argc, argv);
  dnest_postprocess(1.0);
  if(thistask == 0)
  {
    int j;
    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_thismodel+ j), *((double *)best_model_std_thismodel + j));
  }
  
  /* free memory */
  free(perturb_accept);
}

void get_num_particles3(char *fname)
{
  FILE *fp;
  char buf[BUF_MAX_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, BUF_MAX_LENGTH, fp);
  }
  sscanf(buf, "%d", &num_particles);
  fclose(fp);
}

/*====================================================*/
/* users responsible for following struct definitions */

void from_prior_thismodel3(void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    params[i] = -8.0 + 16.0 * dnest_rand();
    wrap(&params[i], -8.0, 8.0);
  }
}

double log_likelihoods_cal_thismodel3(const void *model)
{
  double *params = (double *)model;
  double logL;
  double logl1;
  double logl2;

  logl1 = -0.5 * pow( (sqrt((params[0] - 4.0)* (params[0] - 4.0) + params[1]*params[1])- 2.0)/0.1,  2.0)
          - 0.5 * log(2.0*M_PI * 0.01);
  logl2 = -0.5 * pow( (sqrt((params[0] + 4.0)* (params[0] + 4.0) + params[1]*params[1])- 2.0)/0.1,  2.0)
          - 0.5 * log(2.0*M_PI * 0.01);
  
  double max = fmax(logl1, logl2);
  
  logL = log( exp(logl1-max) + exp(logl2-max) ) + max;
  
  return logL;
}

double perturb_thismodel3(void *model)
{
  double *params = (double *)model;
  double logH = 0.0, width, limit1, limit2;
  int which = dnest_rand_int(num_params);
  //if(which_level_update != 0 )
  //{
  //  limit1 = limits[(which_level_update-1) * num_params *2 + which *2 ];
  //  limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
  //}
  //else
  {
    limit1 = -8.0;
    limit2 = 8.0;
  }
  width = (limit2 - limit1);
  params[which] += width * dnest_randh();
  wrap(&params[which], -8.0, 8.0);
  return logH;
}
/*=======================================================*/

void print_particle_thismodel3(FILE *fp, const void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", params[i]);
  }
  fprintf(fp, "\n");
}

void data_load_thismodel3()
{
  
}

void copy_best_model_thismodel3(const void *bm, const void *bm_std)
{
  memcpy(best_model_thismodel, bm, size_of_modeltype);
  memcpy(best_model_std_thismodel, bm_std, size_of_modeltype);
}
/*========================================================*/


void copy_model_thismodel3(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void* create_model_thismodel3()
{
  return (void *)malloc(size_of_modeltype);
}

int get_num_params_thismodel3()
{
  return num_params;
}
