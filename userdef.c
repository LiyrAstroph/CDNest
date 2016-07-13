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

#include "dnestvars.h"
#include "userdef.h"

int num_data_points;
DataType *data;


int main(int argc, char **argv)
{ 
  /* setup szie of modeltype, which is used for dnest */
  size_of_modeltype = sizeof(ModelType);
  
  /* setup number of data points and allocate memory */
  num_data_points = 100;
  data = (DataType *)malloc(num_data_points * sizeof(DataType));
  
  /* load data */
  data_load();
  
  /* run dnest */
  strcpy(options_file, "OPTIONS");
  dnest(argc, argv);
  
  
  strcpy(options_file, "OPTIONS_2D");
  dnest(argc, argv);
  
  /* free memory */
  free(data);
  return 0;
}

/*====================================================*/
/* users responsible for following struct definitions */

void from_prior(const void *model)
{
  int i;
  ModelType *pm = (ModelType *)model;
  for(i=0; i<num_params; i++)
  {
    (*pm).params[i] = -0.5 + dnest_rand();
  }
  
}

double log_likelihoods_cal(const void *model)
{
  ModelType *pm = (ModelType *)model;
  double logL;
  const double u = 0.01;
	const double v = 0.1;
	const double C = log(1.0/sqrt(2*M_PI));

	double logl1 = num_params*(C - log(u));
	double logl2 = num_params*(C - log(v));

  int i;
	for(i=0; i<num_params; i++)
	{
		logl1 += -0.5*pow(((*pm).params[i] - 0.031)/u, 2);
		logl2 += -0.5*pow( (*pm).params[i]/v, 2);
	}
	logl1 += log(100.);
  
  double max = fmax(logl1, logl2);
  
  logL = log( exp(logl1-max) + exp(logl2-max) ) + max;
  
  return logL;
}

double perturb(const void *model)
{
  ModelType *pm = (ModelType *)model;
  double logH = 0.0;
  int which = dnest_rand_int(num_params);
	(*pm).params[which] += dnest_randh();
	wrap(&pm->params[which], -0.5, 0.5);
  return logH;
}

void print_particle(FILE *fp, const void *model)
{
  int i;
  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", ((ModelType *)model)->params[i]);
  }
  fprintf(fp, "\n");
}

void data_load()
{
  FILE *fp;
  int i;

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
/*========================================================*/

void copy_model(const void *dest, const void *src)
{
  memcpy(dest, src, sizeof(ModelType));
}

void *create_model()
{
  return (void *)malloc(sizeof(ModelType));
}

int get_num_params()
{
  return num_params;
}
