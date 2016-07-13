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

#define STR_MAX_LENGTH (100)
#define BUF_MAX_LENGTH (200)
#define LEVEL_NUM_MAX (1000)

/* output files */
extern FILE *fsample, *fsample_info;

/* random number generator */
extern const gsl_rng_type * dnest_gsl_T;
extern gsl_rng * dnest_gsl_r;

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
  char posterior_sample_file[STR_MAX_LENGTH];
}Options;
extern Options options;
extern char options_file[STR_MAX_LENGTH];

extern void *particles;
extern int size_of_modeltype;
extern int particle_offset_size;

// sampler
extern bool save_to_disk;
extern unsigned int num_threads;
extern double compression;
extern unsigned int regularisation;

extern void *particles;
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


//***********************************************
/*                  functions                  */
double mod(double y, double x);
void wrap(double *x, double min, double max);
int mod_int(int y, int x);
int cmp(const void *pa, const void *pb);

void options_load();
void setup(int argc, char** argv);
void finalise();

int dnest(int argc, char **argv);
void run();
void mcmc_run();
void update_particle(unsigned int which);
void update_level_assignment(unsigned int which);
double log_push(unsigned int which_level);
bool enough_levels();
void do_bookkeeping();
void save_levels();
void save_particle();
void kill_lagging_particles();
void renormalise_visits();
void recalculate_log_X();
double dnest_randh();
double dnest_rand();
double dnest_randn();
int dnest_rand_int(int size);
void postprocess();
void initialize_output_file();
void close_output_file();
/*=====================================================*/
// users responsible for following functions
extern void (*data_load)();
extern void (*print_particle)(FILE *fp, const void *model);
extern void (*from_prior)(const void *model);
extern double (*log_likelihoods_cal)(const void *model);
extern double (*perturb)(const void *model);
extern void (*copy_model)(const void *dest, const void *src);
extern void* (*create_model)();
extern int (*get_num_params)();
/*=====================================================*/

#endif
