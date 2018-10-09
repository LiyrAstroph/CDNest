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
  char posterior_sample_info_file[STR_MAX_LENGTH];
  char limits_file[STR_MAX_LENGTH];
}Options;
extern Options options;
extern char options_file[STR_MAX_LENGTH];

extern void *particles;
extern int size_of_modeltype;
extern int particle_offset_size, particle_offset_double;

// sampler
extern bool save_to_disk;
extern unsigned int num_threads;
extern double compression;
extern unsigned int regularisation;

extern void *particles;
extern LikelihoodType *log_likelihoods;
extern unsigned int *level_assignments;

extern int size_levels, size_levels_combine;  
extern Level *levels;
extern Level *copies_of_levels, *levels_combine;
extern LikelihoodType *all_above;
extern unsigned int count_saves, num_saves;
extern unsigned long long int count_mcmc_steps;
extern LikelihoodType *above;
extern unsigned int size_above, size_all_above;

extern int dnest_flag_restart, dnest_flag_postprc, dnest_flag_sample_info, dnest_flag_limits;
extern double dnest_post_temp;
extern char file_restart[STR_MAX_LENGTH], file_save_restart[STR_MAX_LENGTH];

extern double post_logz;

//the limits of parameters for each level;
extern double *limits, *copies_of_limits;

extern int dnest_which_particle_update; // which particle to be updated
extern int dnest_which_level_update;    // which level to be updated;
extern int dnest_thistask, dnest_totaltask;
extern int *dnest_perturb_accept;
extern int dnest_root;
//***********************************************
/*                  functions                  */
double mod(double y, double x);
void wrap(double *x, double min, double max);
void wrap_limit(double *x, double min, double max);
int mod_int(int y, int x);
int dnest_cmp(const void *pa, const void *pb);

void options_load();
void setup(int argc, char** argv, int num_params);
void finalise();

double dnest(int argc, char **argv, int num_params);
void dnest_run();
void dnest_mcmc_run();
void update_particle(unsigned int which);
void update_level_assignment(unsigned int which);
double log_push(unsigned int which_level);
bool enough_levels(Level *l, int size_l);
void do_bookkeeping();
void save_levels();
void save_particle();
void save_limits();
void kill_lagging_particles();
void renormalise_visits();
void recalculate_log_X();
double dnest_randh();
double dnest_rand();
double dnest_randn();
int dnest_rand_int(int size);
void dnest_postprocess(double temperature);
void postprocess(double temperature);
void initialize_output_file();
void close_output_file();
void dnest_save_restart();
void dnest_restart();
int dnest_get_which_level_update();
int dnest_get_which_particle_update();
unsigned int dnest_get_which_num_saves();
unsigned long long int dnest_get_count_mcmc_steps();
/*=====================================================*/
// users responsible for following functions
extern void (*print_particle)(FILE *fp, const void *model);
extern void (*from_prior)(void *model);
extern double (*log_likelihoods_cal)(const void *model);
extern double (*log_likelihoods_cal_initial)(const void *model);
extern double (*log_likelihoods_cal_restart)(const void *model);
extern double (*perturb)(void *model);
extern void (*restart_action)(int iflag);
/*=====================================================*/

#endif
