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
#include "model2.h"

int which_particle_update;
int num_data_points;
int num_params;

DataType *data;
void *best_model_thismodel, *best_model_std_thismodel;

void model2(int argc, char **argv)
{ 
  /* setup szie of modeltype, which is used for dnest */
  num_params = 3;
  size_of_modeltype = num_params * sizeof(double);
  best_model_thismodel = malloc(size_of_modeltype);
  best_model_std_thismodel = malloc(size_of_modeltype);
  
  /* setup number of data points and allocate memory */
  num_data_points = 30;
  data = (DataType *)malloc(num_data_points * sizeof(DataType));
  
  /* setup functions used for dnest*/
  from_prior = from_prior_thismodel2;
  data_load = data_load_thismodel2;
  log_likelihoods_cal = log_likelihoods_cal_thismodel2;
  perturb = perturb_thismodel2;
  print_particle = print_particle_thismodel2;
  copy_model = copy_model_thismodel2;
  create_model = create_model_thismodel2;
  get_num_params = get_num_params_thismodel2;
  copy_best_model = copy_best_model_thismodel2;
  
  /* load data */
  if(thistask == 0)
  {
    data_load();
  }
  MPI_Bcast(data, num_data_points*sizeof(DataType), MPI_BYTE, 0, MPI_COMM_WORLD);
  
  /* run dnest */
  strcpy(options_file, "OPTIONS2");
  dnest(argc, argv);
  
  if(thistask == 0)
  {
    int j;
    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_thismodel+ j), *((double *)best_model_std_thismodel + j));
  }
  
  /* free memory */
  free(data);
}

/*====================================================*/
/* users responsible for following struct definitions */

void from_prior_thismodel2(void *model)
{
  double *params = (double *)model;

  params[0] = 1E3*dnest_randn();
  params[1] = 1E3*dnest_randn();

	// Log-uniform prior
  params[2] = exp(-10. + 20.*dnest_rand());

}

double log_likelihoods_cal_thismodel2(const void *model)
{
  double *params = (double *)model;
  double logL = 0.0;
  double var = params[2] * params[2];
  
  int i;
  double mu;

	// Conventional gaussian sampling distribution
  for(i=0; i<num_data_points; i++)
  {
    mu = params[0] * data[i].x + params[1];
    logL += -0.5*log(2*M_PI*var) - 0.5*pow(data[i].y - mu, 2.0)/var;
  }
  
  return logL;
}

double perturb_thismodel2(void *model)
{ 
  double *params = (double *)model;
  double logH = 0.0;
  int which = dnest_rand_int(num_params);
  
	if(which == 0)
	{
		logH -= -0.5*pow(params[0]/1E3, 2);
		params[0] += 1E3*dnest_randh();
		logH += -0.5*pow(params[0]/1E3, 2);
	}
	else if(which == 1)
	{
		logH -= -0.5*pow(params[1]/1E3, 2);
		params[1] += 1E3*dnest_randh();
		logH += -0.5*pow(params[1]/1E3, 2);
	}
	else
	{
		// Usual log-uniform prior trick
		params[2] = log(params[2]);
		params[2] += 20.*dnest_randh();
		wrap(&(params[2]), -10., 10.);
		params[2] = exp(params[2]);
	}
  
  return logH;
}
/*=======================================================*/

void print_particle_thismodel2(FILE *fp, const void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", params[i]);
  }
  fprintf(fp, "\n");
}

void data_load_thismodel2()
{
  FILE *fp;
  int i;

  fp = fopen("road.txt", "r");
  if(fp == NULL)
  {
    fprintf(stderr, "ERROR: Cannot open file data.txt.\n");
    exit(0);
  }

  for(i=0; i<num_data_points; i++)
  {
    fscanf(fp, "%lf %lf\n", &data[i].x, &data[i].y);
    //printf("%f %f\n", data[i].x, data[i].y);
  }
  fclose(fp);
}

void copy_best_model_thismodel2(const void *bm, const void *bm_std)
{
  memcpy(best_model_thismodel, bm, size_of_modeltype);
  memcpy(best_model_std_thismodel, bm_std, size_of_modeltype);
}
/*========================================================*/


void copy_model_thismodel2(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void* create_model_thismodel2()
{
  return (void *)malloc(size_of_modeltype);
}

int get_num_params_thismodel2()
{
  return num_params;
}
