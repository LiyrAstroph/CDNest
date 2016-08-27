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
#include "model1.h"

int which_particle_update;
int *perturb_accept;
int num_data_points;
int num_params;
int num_particles;

DataType *data;
void *best_model_thismodel, *best_model_std_thismodel;

void model1(int argc, char **argv)
{ 
  /* setup szie of modeltype, which is used for dnest */
  num_params = 20;
  size_of_modeltype = num_params * sizeof(double);
  best_model_thismodel = malloc(size_of_modeltype);
  best_model_std_thismodel = malloc(size_of_modeltype);
  
  /* setup number of data points and allocate memory */
  num_data_points = 0;
  data = (DataType *)malloc(num_data_points * sizeof(DataType));
  
  /* setup functions used for dnest*/
  from_prior = from_prior_thismodel;
  data_load = data_load_thismodel;
  log_likelihoods_cal = log_likelihoods_cal_thismodel;
  perturb = perturb_thismodel;
  print_particle = print_particle_thismodel;
  copy_model = copy_model_thismodel;
  create_model = create_model_thismodel;
  get_num_params = get_num_params_thismodel;
  copy_best_model = copy_best_model_thismodel;
  
  /* load data */
  if(thistask == 0)
  {
    data_load();
  }
  
  /* run dnest */
  strcpy(options_file, "OPTIONS1");
  if(thistask == 0)
  {
    get_num_particles1(options_file);
  }
  MPI_Bcast(&num_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);

  perturb_accept = malloc(num_particles * sizeof(int));

  dnest(argc, argv);
  
  if(thistask == 0)
  {
    int j;
    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_thismodel + j), *((double *)best_model_std_thismodel+j) );
  }
    
  /* free memory */
  free(data);
  free(perturb_accept);
}

void get_num_particles1(char *fname)
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

void from_prior_thismodel(void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    params[i] = -0.5 + dnest_rand();
  }
}

double log_likelihoods_cal_thismodel(const void *model)
{
  double *params = (double *)model;
  double logL;
  const double u = 0.01;
	const double v = 0.1;
	const double C = log(1.0/sqrt(2*M_PI));

	double logl1 = num_params*(C - log(u));
	double logl2 = num_params*(C - log(v));

  int i;
	for(i=0; i<num_params; i++)
	{
		logl1 += -0.5*pow(( params[i] - 0.031)/u, 2);
		logl2 += -0.5*pow(  params[i]/v, 2);
	}
	logl1 += log(100.);
  
  double max = fmax(logl1, logl2);
  
  logL = log( exp(logl1-max) + exp(logl2-max) ) + max;
  
  return logL;
}

double perturb_thismodel(void *model)
{
  double *params = (double *)model;
  double logH = 0.0;
  int which = dnest_rand_int(num_params);
	params[which] += dnest_randh();
	wrap(&params[which], -0.5, 0.5);
  return logH;
}

void print_particle_thismodel(FILE *fp, const void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", params[i] );
  }
  fprintf(fp, "\n");
}

void data_load_thismodel()
{
  // no data for this model
  return;
}

void copy_best_model_thismodel(const void *bm, const void *bm_std)
{  
  memcpy(best_model_thismodel, bm, size_of_modeltype);
  memcpy(best_model_std_thismodel, bm_std, size_of_modeltype);
}
/*========================================================*/


void copy_model_thismodel(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void* create_model_thismodel()
{
  return (void *)malloc(size_of_modeltype);
}

int get_num_params_thismodel()
{
  return num_params;
}
