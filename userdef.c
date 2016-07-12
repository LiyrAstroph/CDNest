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
#include <gsl/gsl_rng.h>

#include "userdef.h"
#include "dnestproto.h"

int main(int argc, char **argv)
{
  
  data_load();
  dnest(argc, argv);
  return 1;
}

/*====================================================*/
/* users responsible for following struct definitions */

void from_prior(ModelType *model)
{
  int i;
  for(i=0; i<num_params; i++)
		model->params[i] = -0.5 + dnest_rand();
}

double log_likelihoods_cal(ModelType *model)
{
  double logL;
  const double u = 0.01;
	const double v = 0.1;
	const double C = log(1.0/sqrt(2*M_PI));

	double logl1 = num_params*(C - log(u));
	double logl2 = num_params*(C - log(v));

  int i;
	for(i=0; i<num_params; i++)
	{
		logl1 += -0.5*pow(((*model).params[i] - 0.031)/u, 2);
		logl2 += -0.5*pow((*model).params[i]/v, 2);
	}
	logl1 += log(100.);
  
  double max = fmax(logl1, logl2);
  
  logL = log( exp(logl1-max) + exp(logl2-max) ) + max;
  
  return logL;
}

double perturb(ModelType *model)
{
  double logH = 0.0;
  int which = dnest_rand_int(num_params);
	model->params[which] += dnest_randh();
	wrap(&model->params[which], -0.5, 0.5);
  return logH;
}

void print_particle(FILE *fp, ModelType *model)
{
  int i;
  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", (*model).params[i]);
  }
  fprintf(fp, "\n");
}

void data_load()
{
  FILE *fp;
  int i;

  num_data_points = 100;
  data = (DataType *)malloc(num_data_points * sizeof(DataType));

  fp = fopen("data.txt", "r");
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
