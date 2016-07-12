/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */
#ifndef _DNESTVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "userdef.h"

#define STR_MAX_LENGTH (100)
#define BUF_MAX_LENGTH (200)
#define LEVEL_NUM_MAX (1000)

/* output files */
extern FILE *fsample, *fsample_info;

/* random number generator */
extern gsl_rng_type * dnest_gsl_T;
extern const gsl_rng * dnest_gsl_r;

typedef struct 
{
  double value;
  double tiebreaker;
}LikelihoodType;

typedef struct
{
  LikelihoodType log_likelihood;
  double log_X;
  unsigned long long int visits, exceeds;
  unsigned long long int accepts, tries;
}Level;

// struct for options
typedef struct
{
  unsigned int num_particles;
  unsigned int new_level_interval;
  unsigned int save_interval;
  unsigned int thread_steps;
  unsigned int max_num_levels;
  double lambda, beta;
  unsigned int max_num_saves;

  char sample_file[STR_MAX_LENGTH];
  char sample_info_file[STR_MAX_LENGTH];
  char levels_file[STR_MAX_LENGTH];
  char sampler_state_file[STR_MAX_LENGTH];
}Options;

extern Options options;

// data 
extern DataType *data;
extern int num_data_points;

// sampler
extern bool save_to_disk;
extern unsigned int num_threads;
extern double compression;
extern unsigned int regularisation;

extern ModelType *particles;
extern LikelihoodType *log_likelihoods;
extern unsigned int *level_assignments;

extern int size_levels;  
extern Level *levels;
extern Level *copies_of_levels;
extern LikelihoodType *all_above;
extern unsigned int count_saves;
extern unsigned long long int count_mcmc_steps;
extern LikelihoodType *above;
extern int size_above;

#endif
