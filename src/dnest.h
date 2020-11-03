/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */
#ifndef _DNEST_H
#define _DNEST_H

#ifdef __cplusplus
extern "C" {
#endif

#define DNEST_MAJOR_VERSION 0  // Dec 2, 2018
#define DNEST_MINOR_VERSION 1
#define DNEST_PATCH_VERSION 0

typedef struct
{
  void (*from_prior)(void *model);
  double (*log_likelihoods_cal)(const void *model);
  double (*log_likelihoods_cal_initial)(const void *model);
  double (*log_likelihoods_cal_restart)(const void *model);
  double (*perturb)(void *model);
  void (*print_particle)(FILE *fp, const void *model);
  void (*read_particle)(FILE *fp, void *model);
  void (*restart_action)(int iflag);
  void (*accept_action)();
  void (*kill_action)(int i, int i_copy);
}DNestFptrSet;

double dnest_randh();
double dnest_rand();
double dnest_randn();
int dnest_rand_int(int size);
void dnest_wrap(double *x, double min, double max);

int dnest_get_size_levels();
int dnest_get_which_level_update();
int dnest_get_which_particle_update();
void dnest_get_posterior_sample_file(char *fname);
int dnest_check_version(char *verion_str);
unsigned int dnest_get_which_num_saves();
unsigned int dnest_get_count_saves();
unsigned long long int dnest_get_count_mcmc_steps();
void dnest_check_fptrset(DNestFptrSet *fptrset);
DNestFptrSet * dnest_malloc_fptrset();
void dnest_free_fptrset(DNestFptrSet * fptrset);

double dnest(int argc, char **argv, DNestFptrSet *fptrset,  int num_params,  
             double *param_range, int *prior_type, double *prior_info, 
             char *sample_dir, char *optfile, void *args);
             
#ifdef __cplusplus
}
#endif

#endif