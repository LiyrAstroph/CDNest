/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _MODEL2_H

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

/* best model */
extern void *best_model_thismodel, *best_model_std_thismodel;

extern int which_particle_update; // which particule to be updated
extern int which_level_update;
extern int *perturb_accept;
extern double *limits;
extern int thisktask, totaltask;

/* functions */
void get_num_particles2(char *fname);
void from_prior_thismodel2(void *model);
void data_load_thismodel2();
void print_particle_thismodel2(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel2(const void *model);
double perturb_thismodel2(void *model);
void copy_model_thismodel2(void *dest, const void *src);
void* create_model_thismodel2();
int get_num_params_thismodel2();
void copy_best_model_thismodel2(const void *bm, const void *bm_std);

void (*data_load)();
void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*perturb)(void *model);
void (*copy_model)(void *dest, const void *src);
void* (*create_model)();
int (*get_num_params)();
void (*copy_best_model)(const void *bm, const void *bm_std);
#endif
