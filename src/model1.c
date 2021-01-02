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
#include "model1.h"

int which_level_update;
int num_params;

DNestFptrSet *fptrset_thismodel;

void model1()
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
  strcpy(argv[argc++], "restart_dnest1.txt");
  strcpy(argv[argc++], "-l"); //level-dependnet sampling
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "1");

  /* setup szie of modeltype, which is used for dnest */
  num_params = 20;
  
  /* allocate memory */
  fptrset_thismodel = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_thismodel->from_prior = from_prior_thismodel;
  fptrset_thismodel->log_likelihoods_cal = log_likelihoods_cal_thismodel;
  fptrset_thismodel->log_likelihoods_cal_initial = log_likelihoods_cal_thismodel;
  fptrset_thismodel->log_likelihoods_cal_restart = log_likelihoods_cal_thismodel;
  fptrset_thismodel->perturb = perturb_thismodel;
  fptrset_thismodel->print_particle = print_particle_thismodel;
  fptrset_thismodel->restart_action = restart_action_model1;
  
  /* run dnest */
  dnest(argc, argv, fptrset_thismodel, num_params, NULL, NULL, NULL, "./", "OPTIONS1", NULL);
    
  /* free memory */
  dnest_free_fptrset(fptrset_thismodel);

  for(i=0; i<narg; i++)
    free(argv[i]);
  free(argv);
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
  }
  else
  {
    limit1 = -0.5;
    limit2 = 0.5;
  }

  width = (limit2 - limit1);
	params[which] += width * dnest_randh();
	dnest_wrap(&params[which], limit1, limit2);
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

/*========================================================*/

void restart_action_model1(int iflag)
{
  return;
}
