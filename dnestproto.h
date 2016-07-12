/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _DNESTPROTO_H

// functions
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
/*=====================================================*/
// users responsible for following functions
void data_load();
void print_particle(FILE *fp, ModelType *model);
void from_prior(ModelType *model);
double log_likelihoods_cal(ModelType *model);
double perturb(ModelType *model);
/*=====================================================*/

#endif
