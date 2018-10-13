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

int which_level_update;
int num_data_points;
int num_params;
DNestFptrSet *fptrset_thismodel3;

void model3()
{ 
  int i, argc=0, narg=4;
  char **argv;

  argv = malloc(narg*sizeof(char *));
  for(i=0; i<narg; i++)
  {
    argv[i] = malloc(200*sizeof(char));
  }
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc++], "restart_dnest3.txt");
  strcpy(argv[argc++], "-l"); //level-dependnet sampling

  /* setup szie of modeltype, which is used for dnest */
  num_params = 2;
  fptrset_thismodel3 = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_thismodel3->from_prior = from_prior_thismodel3;
  fptrset_thismodel3->log_likelihoods_cal = log_likelihoods_cal_thismodel3;
  fptrset_thismodel3->log_likelihoods_cal_initial = log_likelihoods_cal_thismodel3;
  fptrset_thismodel3->log_likelihoods_cal_restart = log_likelihoods_cal_thismodel3;
  fptrset_thismodel3->perturb = perturb_thismodel3;
  fptrset_thismodel3->print_particle = print_particle_thismodel3;
  fptrset_thismodel3->restart_action = restart_action_model3;
  
  /* run dnest */
  dnest(argc, argv, fptrset_thismodel3, num_params, "OPTIONS3");
  
  /* free memory */
  dnest_free_fptrset(fptrset_thismodel3);

  for(i=0; i<narg; i++)
    free(argv[i]);
  free(argv);

}

/*====================================================*/
/* users responsible for following struct definitions */

void from_prior_thismodel3(void *model)
{
  int i;
  double *params = (double *)model;
  for(i=0; i<num_params; i++)
  {
    params[i] = -6.0 + 12.0 * dnest_rand();
    dnest_wrap(&params[i], -6.0, 6.0);
  }
}

double log_likelihoods_cal_thismodel3(const void *model)
{
  double *params = (double *)model;
  double logL;
  double logl1;
  double logl2;

  logl1 = -0.5 * pow( (sqrt((params[0] - 3.0)* (params[0] - 3.0) + params[1]*params[1])- 2.0)/0.1,  2.0)
          - 0.5 * log(2.0*M_PI * 0.01);
  logl2 = -0.5 * pow( (sqrt((params[0] + 3.0)* (params[0] + 3.0) + params[1]*params[1])- 2.0)/0.1,  2.0)
          - 0.5 * log(2.0*M_PI * 0.01);
  
  double max = fmax(logl1, logl2);
  
  logL = log( exp(logl1-max) + exp(logl2-max) ) + max;
  
  return logL;
}

double perturb_thismodel3(void *model)
{
  double *params = (double *)model;
  double logH = 0.0, width, limit1, limit2;
  int which = dnest_rand_int(num_params), which_level;

  which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels - 20)?(size_levels-20):which_level_update;
  which_level = 0;
  if(which_level > 0 )
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2 ];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
  }
  else
  {
    limit1 = -6.0;
    limit2 = 6.0;
  }
  width = (limit2 - limit1);
  params[which] += width * dnest_randh();
  dnest_wrap(&params[which], -6.0, 6.0);
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

/*========================================================*/

void restart_action_model3(int iflag)
{
  return;
}
