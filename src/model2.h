/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _MODEL2_H
#define _MODEL2_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>

/*  data type */
typedef struct 
{
  double x;
  double y;
}DataType;

/* number of model parameters */
extern int num_params;

/* data storage */
extern int num_data_points;
extern DataType *data;

extern int which_particle_update; // which particule to be updated
extern int which_level_update;
extern double *limits;
extern int thistask, totaltask;

extern DNestFptrSet *fptrset_thismodel2;

/* functions */
void from_prior_thismodel2(void *model);
void data_load_thismodel2();
void print_particle_thismodel2(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel2(const void *model);
double perturb_thismodel2(void *model);
void restart_action_model2(int iflag);

#endif
