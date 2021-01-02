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
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "dnest.h"
#include "model2.h"

int which_level_update;
int num_data_points;
int num_params;
DNestFptrSet *fptrset_thismodel2;

DataType *data;
void *best_model_thismodel, *best_model_std_thismodel;

void model2()
{ 
  int i, argc=0, narg=6;
  char **argv;

  argv = malloc(narg*sizeof(char *));
  for(i=0; i<narg; i++)
  {
    argv[i] = malloc(200*sizeof(char));
  }
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc++], "restart_dnest2.txt");
  strcpy(argv[argc++], "-l"); //level-dependnet sampling
  strcpy(argv[argc++], "-g"); //tag
  strcpy(argv[argc++], "2");

  /* setup szie of modeltype, which is used for dnest */
  num_params = 3;
  
  /* setup number of data points and allocate memory */
  num_data_points = 30;
  data = (DataType *)malloc(num_data_points * sizeof(DataType));
  fptrset_thismodel2 = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_thismodel2->from_prior = from_prior_thismodel2;
  fptrset_thismodel2->log_likelihoods_cal = log_likelihoods_cal_thismodel2;
  fptrset_thismodel2->log_likelihoods_cal_initial = log_likelihoods_cal_thismodel2;
  fptrset_thismodel2->log_likelihoods_cal_restart = log_likelihoods_cal_thismodel2;
  fptrset_thismodel2->perturb = perturb_thismodel2;
  fptrset_thismodel2->print_particle = print_particle_thismodel2;
  fptrset_thismodel2->restart_action = restart_action_model2;
  
  /* load data */
  if(thistask == 0)
  {
    data_load_thismodel2();
  }
  MPI_Bcast(data, num_data_points*sizeof(DataType), MPI_BYTE, 0, MPI_COMM_WORLD);
  
  /* run dnest */
  dnest(argc, argv, fptrset_thismodel2, num_params, NULL, NULL, NULL, "./", "OPTIONS2", NULL);
    
  /* free memory */
  free(data);
  dnest_free_fptrset(fptrset_thismodel2);

  for(i=0; i<narg; i++)
    free(argv[i]);
  free(argv);
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
  double logH = 0.0, width, limit1, limit2;
  int which = dnest_rand_int(num_params), which_level;
  int size_levels;
  
  which_level_update = dnest_get_which_level_update();
  size_levels = dnest_get_size_levels();
  
  which_level = which_level_update > (size_levels - 20)?(size_levels-20):which_level_update;
  if(which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2 ];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = (limit2 - limit1);
  }

	if(which == 0)
	{
    if(which_level_update == 0)
    {
      width = 1E3;
    }
		logH -= -0.5*pow(params[0]/1E3, 2);
		params[0] += width*dnest_randh();
		logH += -0.5*pow(params[0]/1E3, 2);
	}
	else if(which == 1)
	{
    if(which_level_update == 0)
    {
      width = 1E3;
    }

		logH -= -0.5*pow(params[1]/1E3, 2);
		params[1] += width*dnest_randh();
		logH += -0.5*pow(params[1]/1E3, 2);
	}
	else
	{
    if(which_level_update == 0)
    {
      limit1 = -10.0;
      limit2 = 10.0;
      width = 20.0;
    }
    else
    {
      limit1 = log(limit1);
      limit2 = log(limit2);
      width = limit2 - limit1;
    }
    
		// Usual log-uniform prior trick
		params[2] = log(params[2]);
		params[2] += width*dnest_randh();
		dnest_wrap(&(params[2]), -10.0, 10.0);
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
    fprintf(stderr, "ERROR: Cannot open file road.txt.\n");
    exit(0);
  }

  for(i=0; i<num_data_points; i++)
  {
    fscanf(fp, "%lf %lf\n", &data[i].x, &data[i].y);
    //printf("%f %f\n", data[i].x, data[i].y);
  }
  fclose(fp);
}

/*========================================================*/

void restart_action_model2(int iflag)
{
  return;
}