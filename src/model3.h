/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _MODEL3_H
#define _MODEL3_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>

/* number of model parameters */
extern int num_params;

extern int which_particle_update; // which particule to be updated
extern int which_level_update;
extern double *limits;
extern int thistask, totaltask;

extern DNestFptrSet *fptrset_thismodel3;

/* functions */
void from_prior_thismodel3(void *model);
void print_particle_thismodel3(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel3(const void *model);
double perturb_thismodel3(void *model);
void restart_action_model3(int iflag);

#endif
