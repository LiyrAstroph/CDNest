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
#include <time.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "dnestvars.h"
#include "dnestproto.h"

int dnest(int argc, char** argv)
{
  setup(argc, argv);
  //initialize_output_file();
  //run();
  //finalise();
  postprocess();
  
  return 0;
}

void run()
{
  while(true)
  {
    mcmc_run();
    if(options.max_num_saves !=0 &&
        count_saves != 0 && (count_saves%options.max_num_saves == 0))
      break;
      
    count_mcmc_steps += options.thread_steps;

    do_bookkeeping();
  }

  FILE *fp;
  fp = fopen(options.sampler_state_file, "w");
  fprintf(fp, "%d %d\n", size_levels, count_saves);
  fclose(fp);
}

void do_bookkeeping()
{
  int i;
  bool created_level = false;

  if(!enough_levels() && size_above >= options.new_level_interval)
  {
    // in descending order 
    qsort(above, size_above, sizeof(LikelihoodType), cmp);
    int index = (int)( (1.0/compression) * size_above);

    Level level_tmp = {above[index], 0.0, 0, 0, 0, 0};
    levels[size_levels] = level_tmp;
    size_levels++;
    
    printf("# Creating level %d with log likelihood = %e.\n", 
               size_levels-1, levels[size_levels-1].log_likelihood.value);

    // clear out the last index records
    for(i=index; i<size_above; i++)
    {
      above[i].value = 0.0;
      above[i].tiebreaker = 0.0;
    }
    size_above = index;

    if(enough_levels())
    {
      renormalise_visits();
      options.max_num_levels = size_levels;
      printf("# Done creating levles.\n");
    }
    else
    {
      kill_lagging_particles();
    }

    created_level = true;
  }

  recalculate_log_X();

  if(created_level)
    save_levels();

  if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
  {

    save_particle();

    if(!created_level)
      save_levels();
  }
}

void recalculate_log_X()
{
  int i;

  levels[0].log_X = 0.0;
  for(i=1; i<size_levels; i++)
  {
    levels[i].log_X = levels[i-1].log_X 
    + log( (double)( (levels[i-1].exceeds + 1.0/compression * regularisation)
                    /(levels[i-1].visits + regularisation)  ) );
  }
}

void renormalise_visits()
{
  size_t i;

  for(i=0; i<size_levels; i++)
  {
    if(levels[i].tries >= regularisation)
    {
      levels[i].accepts = ((double)(levels[i].accepts+1) / (double)(levels[i].tries+1)) * regularisation;
      levels[i].tries = regularisation;
    }

    if(levels[i].visits >= regularisation)
    {
      levels[i].exceeds = ( (double) (levels[i].exceeds+1) / (double)(levels[i].visits + 1) ) * regularisation;
      levels[i].visits = regularisation;
    }
  }
}

void kill_lagging_particles()
{
  static unsigned int deletions = 0;

  bool *good;
  good = (bool *)malloc(options.num_particles * sizeof(bool));

  double max_log_push = DBL_MAX;

  unsigned int num_bad = 0;
  size_t i;

  for(i=0; i<options.num_particles; i++)good[i] = true;

  for(i=0; i<options.num_particles; i++)
  {
    if( log_push(level_assignments[i]) > max_log_push)
      max_log_push = log_push(level_assignments[i]);

    if( log_push(level_assignments[i]) < -6.0)
    {
      good[i] = false;
      ++num_bad;
    }
  }

  if(num_bad < options.num_particles)
  {
    for(i=0; i< options.num_particles; i++)
    {
      if(!good[i])
      {
        int i_copy;
        do
        {
          i_copy = gsl_rng_uniform_int(dnest_gsl_r, options.num_particles);
        }while(!good[i_copy] || gsl_rng_uniform(dnest_gsl_r) >= exp(log_push(level_assignments[i_copy])));

        particles[i] = particles[i_copy];
        log_likelihoods[i] = log_likelihoods[i_copy];
        level_assignments[i] = level_assignments[i_copy];
        deletions++;

        printf("# Replacing lagging particle.\n");
        printf(" This has happened %d times.\n", deletions);
      }
    }
  }
  else
    printf("# Warning: all particles lagging! Very rare.\n");
}

void save_levels()
{
  if(!save_to_disk)
    return;
  
  int i;
  FILE *fp;

  fp = fopen(options.levels_file, "w");
  fprintf(fp, "# log_X, log_likelihood, tiebreaker, accepts, tries, exceeds, visits\n");
  for(i=0; i<size_levels; i++)
  {
    fprintf(fp, "%14.12g %14.12g %f %llu %llu %llu %llu\n", levels[i].log_X, levels[i].log_likelihood.value, 
      levels[i].log_likelihood.tiebreaker, levels[i].accepts,
      levels[i].tries, levels[i].exceeds, levels[i].visits);
  }
  fclose(fp);
}

void save_particle()
{
  count_saves++;
  
  if(!save_to_disk)
    return;
  
  if(count_saves%1 == 0)
    printf("# Saving particle to disk. N= %d.\n", count_saves);

  int which = gsl_rng_uniform_int(dnest_gsl_r,options.num_particles);


  print_particle(fsample, &particles[which]);

  fprintf(fsample_info, "%d %e %f %d\n", level_assignments[which], 
    log_likelihoods[which].value,
    log_likelihoods[which].tiebreaker,
    which);
}
void mcmc_run()
{
  unsigned int which;
  unsigned int i;
  
  for(i = 0; i<options.thread_steps; i++)
  {
    //if(count_mcmc_steps >= 10000)printf("FFFF %d\n", options.num_particles);
    which = gsl_rng_uniform_int(dnest_gsl_r, options.num_particles);
    //if(count_mcmc_steps >= 10000)printf("FFFF\n");
    //printf("%d\n", which);
    //printf("%f %f %f\n", particles[which].param[0], particles[which].param[1], particles[which].param[2]);
    //printf("level:%d\n", level_assignments[which]);
    //printf("%e\n", log_likelihoods[which].value);
    if(gsl_rng_uniform(dnest_gsl_r) <= 0.5)
    {
      update_particle(which);
      update_level_assignment(which);
    }
    else
    {
      update_level_assignment(which);
      update_particle(which);
    }
    //printf("%f %f %f\n", particles[which].param[0], particles[which].param[1], particles[which].param[2]);
    //printf("level:%d\n", level_assignments[which]);
    //printf("%e\n", log_likelihoods[which].value);

    if( !enough_levels()  && levels[size_levels-1].log_likelihood.value < log_likelihoods[which].value)
    {
      above[size_above] = log_likelihoods[which];
      size_above++;
    }
    //printf("size_above:%d\n", size_above);
  }
}


void update_particle(unsigned int which)
{
  ModelType *particle = &(particles[which]);
  LikelihoodType *logl = &(log_likelihoods[which]);
  
  Level *level = &(levels[level_assignments[which]]);

  ModelType proposal;
  LikelihoodType logl_proposal;
  double log_H;

  //proposal = *particle;
  memcpy(&proposal, particle, sizeof(ModelType));
  
  log_H = perturb(&proposal);
  
  logl_proposal.value = log_likelihoods_cal(&proposal);
  logl_proposal.tiebreaker =  logl_proposal.tiebreaker + gsl_rng_uniform(dnest_gsl_r);
  wrap(&logl_proposal.tiebreaker, 0.0, 1.0);

  if(log_H > 0.0)
    log_H = 0.0;

  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_H) && level->log_likelihood.value < logl_proposal.value)
  {
    //*particle = proposal;
    //*logl = logl_proposal;
    memcpy(particle, &proposal, sizeof(ModelType));
    memcpy(logl, &logl_proposal, sizeof(LikelihoodType));
    level->accepts++;
  }
  level->tries++;

  unsigned int current_level = level_assignments[which];
  for(; current_level < size_levels-1; ++current_level)
  {
    levels[current_level].visits++;
    if(levels[current_level+1].log_likelihood.value < log_likelihoods[which].value)
      levels[current_level].exceeds++;
    else
      break; // exit the loop if it does not satify higher levels
  }
}

void update_level_assignment(unsigned int which)
{
  int proposal = level_assignments[which] 
                 + (int)( pow(10.0, 2*gsl_rng_uniform(dnest_gsl_r))*gsl_ran_gaussian(dnest_gsl_r, 1.0));

  if(proposal == level_assignments[which])
    proposal =  ((gsl_rng_uniform(dnest_gsl_r) < 0.5)?(proposal-1):(proposal+1));

  proposal=mod_int(proposal, size_levels);

  double log_A = -levels[proposal].log_X + levels[level_assignments[which]].log_X;
  
  log_A += log_push(proposal) - log_push(level_assignments[which]);

  if(size_levels == options.max_num_levels)
    log_A += options.beta*log( (double)(levels[level_assignments[which]].tries +1)/ (levels[proposal].tries +1) );

  if(log_A > 0.0)
    log_A = 0.0;

  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_A) && levels[proposal].log_likelihood.value < log_likelihoods[which].value)
  {
    level_assignments[which] = proposal;
  }

}

double log_push(unsigned int which_level)
{
  if(which_level > size_levels)
  {
    printf("level overflow.\n");
    exit(0);
  }
  if(enough_levels())
    return 0.0;

  int i = which_level - (size_levels - 1);
  return i/options.lambda;
}

bool enough_levels()
{
  int i;

  if(options.max_num_levels == 0)
  {
    if(size_levels < 10)
      return false;

    int num_levels_to_check = 20;
    if(size_levels > 80)
      num_levels_to_check = (int)(sqrt(20) * sqrt(0.25*size_levels));

    int k = size_levels - 1;

    for(i= 0; i<num_levels_to_check; i++)
    {
      if(levels[k].log_likelihood.value - levels[k-1].log_likelihood.value
        >=0.8)
        return false;
      k--;
    }
    return true;
  }
  return (size_levels >= options.max_num_levels);
}

void initialize_output_file()
{
  fsample = fopen(options.sample_file, "w");
  if(fsample==NULL)
  {
    fprintf(stderr, "# Cannot open file sample.txt.\n");
    exit(0);
  }
  fprintf(fsample, "# \n");
  fsample_info = fopen(options.sample_info_file, "w");
  if(fsample_info==NULL)
  {
    fprintf(stderr, "# Cannot open file sample_info.txt.\n");
    exit(0);
  }
  fprintf(fsample_info, "# level assignment, log likelihood, tiebreaker, ID.\n");
}

void setup(int argc, char** argv)
{
  int i;

  dnest_gsl_T = (gsl_rng_type *) gsl_rng_default;
  dnest_gsl_r = gsl_rng_alloc (dnest_gsl_T);
  gsl_rng_set(dnest_gsl_r, time(NULL));

  options_load();
  
  num_threads = 1;
  compression = exp(1.0);
  regularisation = options.new_level_interval;
  save_to_disk = true;

// initialise sampler
  all_above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));
  above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  particles = (ModelType *)malloc(options.num_particles*sizeof(ModelType));
  log_likelihoods = (LikelihoodType *)malloc(2*options.num_particles * sizeof(LikelihoodType));
  level_assignments = (unsigned int*)malloc(options.num_particles * sizeof(unsigned int));

  if(options.max_num_levels != 0)
    levels = (Level *)malloc(options.max_num_levels * sizeof(Level));
  else
    levels = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));
  copies_of_levels = (Level *)malloc(options.max_num_levels * sizeof(Level));

  count_mcmc_steps = 0;
  count_saves = 0;

// first level
  size_levels = 0;
  size_above = 0;
  LikelihoodType like_tmp = {-DBL_MAX, gsl_rng_uniform(dnest_gsl_r)};
  Level level_tmp = {like_tmp, 0.0, 0, 0, 0, 0};
  levels[size_levels] = level_tmp;
  size_levels++;
  
  for(i=0; i<options.num_particles; i++)
  {
    from_prior(&particles[i]);
    
    log_likelihoods[i].value = log_likelihoods_cal(&particles[i]);
    log_likelihoods[i].tiebreaker = gsl_rng_uniform(dnest_gsl_r);
    level_assignments[i] = 0;
  }
  
  /*ModelType proposal;
  printf("%f %f %f \n", particles[0].param[0], particles[0].param[1], particles[0].param[2] );
  proposal = particles[0];
  printf("%f %f %f \n", proposal.param[0], proposal.param[1], proposal.param[2] );*/
}

void finalise()
{
  fclose(fsample);
  fclose(fsample_info);
}


void options_load()
{
  FILE *fp;
  char buf[BUF_MAX_LENGTH];

  fp = fopen("OPTIONS", "r");

  if(fp == NULL)
  {
    fprintf(stderr, "ERROR: Cannot open file OPTIONS.\n");
    exit(0);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, BUF_MAX_LENGTH, fp);
  }
  sscanf(buf, "%d", &options.num_particles);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%d", &options.new_level_interval);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%d", &options.save_interval);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%d", &options.thread_steps);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%d", &options.max_num_levels);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%lf", &options.lambda);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%lf", &options.beta);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%d", &options.max_num_saves);

  fclose(fp);

  strcpy(options.sample_file, "sample.txt");
  strcpy(options.sample_info_file, "sample_info.txt");
  strcpy(options.levels_file, "levels.txt");
  strcpy(options.sampler_state_file, "sampler_state.txt");
}


double mod(double y, double x)
{
   if(x <= 0)
     printf("Warning in mod(double, double) (Utils.cpp)");
   return (y/x - floor(y/x))*x;
}

void wrap(double *x, double min, double max)
{
   *x = mod(*x - min, max - min) + min;
}

int mod_int(int y, int x)
{
  if(y >= 0)
    return y - (y/x)*x;
  else
    return (x-1) - mod_int(-y-1, x);
}

double dnest_randh()
{
  return pow(10.0, 1.5 - 3.0*fabs(gsl_ran_tdist(dnest_gsl_r, 2))) * gsl_ran_gaussian(dnest_gsl_r, 1.0);
}

double dnest_rand()
{
  return gsl_rng_uniform(dnest_gsl_r);
}

int dnest_rand_int(int size)
{
  return gsl_rng_uniform_int(dnest_gsl_r, size);
}

double dnest_randn()
{
  return gsl_ran_gaussian(dnest_gsl_r, 1.0);
}

int cmp(const void *pa, const void *pb)
{
  LikelihoodType *a = (LikelihoodType *)pa;
  LikelihoodType *b = (LikelihoodType *)pb;

  // in decesending order
  if(a->value > b->value)
    return false;
  if( a->value == b->value && a->tiebreaker > b->tiebreaker)
    return false;
  
  return true;
}
