/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _MODEL1_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>


/* number of model parameters */
extern int num_params;

extern int which_level_update;
extern double *limits;     // limits from dnest
extern int thistask, totaltask;

extern DNestFptrSet *fptrset_thismodel;

/* functions */
void from_prior_thismodel(void *model);
void print_particle_thismodel(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel(const void *model);
double perturb_thismodel(void *model);
void restart_action_model1(int iflag);
#endif
