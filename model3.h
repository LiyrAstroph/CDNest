/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _MODEL3_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>

/*===================================================*/
// users responsible for following struct definitions

/*  data type */
typedef struct 
{
  double x;
  double y;
}DataType;

/*===================================================*/

/* number of model parameters */
extern int num_params;

/* szie of modeltype, which is used for dnest */
extern int size_of_modeltype;

/* data storage */
extern int num_data_points;
extern DataType *data;

extern int which_particle_update; // which particule to be updated
extern int which_level_update;
extern int *perturb_accept;
extern double *limits;
extern int thistask, totaltask;

/* functions */
void get_num_particles3(char *fname);
void from_prior_thismodel3(void *model);
void data_load_thismodel3();
void print_particle_thismodel3(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel3(const void *model);
double perturb_thismodel3(void *model);
int get_num_params_thismodel3();
void restart_action_model3(int iflag);

void (*data_load)();
void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*perturb)(void *model);
int (*get_num_params)();
void (*restart_action)(int iflag);
#endif
