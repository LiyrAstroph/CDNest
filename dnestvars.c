/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "dnestvars.h"
 
/* output files */
FILE *fsample, *fsample_info;

/* random number generator */
gsl_rng_type * dnest_gsl_T;
const gsl_rng * dnest_gsl_r;

Options options;

/* data */
DataType *data;
int num_data_points;

// sampler
bool save_to_disk;
unsigned int num_threads;
double compression;
unsigned int regularisation;

ModelType *particles;
LikelihoodType *log_likelihoods;
unsigned int *level_assignments;

int size_levels;  
Level *levels;
Level *copies_of_levels;
LikelihoodType *all_above;
unsigned int count_saves;
unsigned long long int count_mcmc_steps;
LikelihoodType *above;
int size_above;
