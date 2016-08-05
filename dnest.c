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
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "dnestvars.h"

int dnest(int argc, char** argv)
{
  setup(argc, argv);
  initialize_output_file();
  dnest_run();
  close_output_file();
   
  // postprocess, calculate evidence, generate posterior sample.
  double temperature = 1.0;
  if(thistask == root)
  {
    postprocess(temperature);
  }

  finalise();
  return 0;
}

void dnest_run()
{
  int i, j, size_all_above_incr;
  Level *pl, *levels_orig;
  int *buf_size_above, *buf_displs;
  
  if(thistask == root)
  {
    buf_size_above = malloc(totaltask * sizeof(int));  
    buf_displs = malloc(totaltask * sizeof(int));
  }

  while(true)
  {
    dnest_mcmc_run();
    MPI_Barrier(MPI_COMM_WORLD);

    //check for termination
    if(options.max_num_saves !=0 &&
        count_saves != 0 && (count_saves%options.max_num_saves == 0))
      break;
    
    //gather levels
    MPI_Gather(levels, size_levels*sizeof(Level), MPI_BYTE, 
             copies_of_levels, size_levels*sizeof(Level), MPI_BYTE, root, MPI_COMM_WORLD);
    
    //gather size_above 
    MPI_Gather(&size_above, 1, MPI_INT, buf_size_above, 1, MPI_INT, root, MPI_COMM_WORLD);

    // task 0 responsible for updating levels
    if(thistask == root)
    {
      size_all_above_incr = 0;
      for(i = 0; i<totaltask; i++)
      {
        size_all_above_incr += buf_size_above[i];

        buf_size_above[i] *= sizeof(LikelihoodType);
      }

      buf_displs[0] = size_all_above * sizeof(LikelihoodType);
      for(i=1; i<totaltask; i++)
      {
        buf_displs[i] = buf_displs[i-1] + buf_size_above[i-1];
      }
      
      //update size_all_above
      size_all_above += size_all_above_incr; 
    }
    
    // gather above into all_above, stored in task 0, note that its size is different among tasks
    MPI_Gatherv(above, size_above * sizeof(LikelihoodType), MPI_BYTE, 
                all_above, buf_size_above, buf_displs, MPI_BYTE, root, MPI_COMM_WORLD);

    // reset size_above 
    size_above = 0;

    count_mcmc_steps += options.thread_steps * totaltask;

    if(thistask == root)
    {
      //backup levels_combine
      levels_orig = malloc(size_levels_combine * sizeof(Level));
      memcpy(levels_orig, levels_combine, size_levels_combine*sizeof(Level));

      //scan over all copies of levels
      pl = copies_of_levels;
      for(i=0; i< totaltask; i++)
      {
        for(j=0; j<size_levels_combine; j++)
        {
          levels_combine[j].accepts += (pl[j].accepts - levels_orig[j].accepts);
          levels_combine[j].tries += (pl[j].tries - levels_orig[j].tries);
          levels_combine[j].visits += (pl[j].visits - levels_orig[j].visits);
          levels_combine[j].exceeds += (pl[j].exceeds - levels_orig[j].exceeds);
          //printf("%d %d\n", thistask, (pl[j].accepts - levels_orig[j].accepts) );;
        }
        pl += size_levels_combine;
      }

      free(levels_orig);

      do_bookkeeping();

      size_levels = size_levels_combine;
      memcpy(levels, levels_combine, size_levels * sizeof(Level));
    }

    //broadcast levels
    MPI_Bcast(&size_levels, 1, MPI_INT, root, MPI_COMM_WORLD);
    //MPI_Bcast(&count_saves, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(levels, size_levels * sizeof(Level), MPI_BYTE, root,  MPI_COMM_WORLD); 

    if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
    {
      save_particle();
    }
  }

  if(thistask == root)
  {
    /* output state of sampler */
    FILE *fp;
    fp = fopen(options.sampler_state_file, "w");
    fprintf(fp, "%d %d\n", size_levels, count_saves);
    fclose(fp);

    free(buf_size_above);
    free(buf_displs);
  }
}

void do_bookkeeping()
{
  int i;
  bool created_level = false;

  if(!enough_levels(levels_combine, size_levels_combine) && size_all_above >= options.new_level_interval)
  {
    // in descending order 
    qsort(all_above, size_all_above, sizeof(LikelihoodType), cmp);
    int index = (int)( (1.0/compression) * size_all_above);

    Level level_tmp = {all_above[index], 0.0, 0, 0, 0, 0};
    levels_combine[size_levels_combine] = level_tmp;
    size_levels_combine++;
    
    printf("# Creating level %d with log likelihood = %e.\n", 
               size_levels_combine-1, levels_combine[size_levels_combine-1].log_likelihood.value);

    // clear out the last index records
    for(i=index; i<size_all_above; i++)
    {
      all_above[i].value = 0.0;
      all_above[i].tiebreaker = 0.0;
    }
    size_all_above = index;

    if(enough_levels(levels_combine, size_levels_combine))
    {
      renormalise_visits();
      options.max_num_levels = size_levels_combine;
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

    //save_particle();

    if(!created_level)
      save_levels();
  }

}

void recalculate_log_X()
{
  int i;

  levels[0].log_X = 0.0;
  for(i=1; i<size_levels_combine; i++)
  {
    levels_combine[i].log_X = levels_combine[i-1].log_X 
    + log( (double)( (levels_combine[i-1].exceeds + 1.0/compression * regularisation)
                    /(levels_combine[i-1].visits + regularisation)  ) );
  }
}

void renormalise_visits()
{
  size_t i;

  for(i=0; i<size_levels_combine; i++)
  {
    if(levels_combine[i].tries >= regularisation)
    {
      levels_combine[i].accepts = ((double)(levels_combine[i].accepts+1) / (double)(levels_combine[i].tries+1)) * regularisation;
      levels_combine[i].tries = regularisation;
    }

    if(levels_combine[i].visits >= regularisation)
    {
      levels_combine[i].exceeds = ( (double) (levels_combine[i].exceeds+1) / (double)(levels_combine[i].visits + 1) ) * regularisation;
      levels_combine[i].visits = regularisation;
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

        copy_model(particles+i*particle_offset_size, particles + i_copy*particle_offset_size);
        log_likelihoods[i] = log_likelihoods[i_copy];
        level_assignments[i] = level_assignments[i_copy];
        deletions++;

        printf("# Replacing lagging particle.\n");
        printf("# This has happened %d times.\n", deletions);
      }
    }
  }
  else
    printf("# Warning: all particles lagging! Very rare.\n");

  free(good);
}

/* save levels */
void save_levels()
{
  if(!save_to_disk)
    return;
  
  int i;
  FILE *fp;

  fp = fopen(options.levels_file, "w");
  fprintf(fp, "# log_X, log_likelihood, tiebreaker, accepts, tries, exceeds, visits\n");
  for(i=0; i<size_levels_combine; i++)
  {
    fprintf(fp, "%14.12g %14.12g %f %llu %llu %llu %llu\n", levels_combine[i].log_X, levels_combine[i].log_likelihood.value, 
      levels_combine[i].log_likelihood.tiebreaker, levels_combine[i].accepts,
      levels_combine[i].tries, levels[i].exceeds, levels_combine[i].visits);
  }
  fclose(fp);
}

/* save particle */
void save_particle()
{
  count_saves++;
  
  if(!save_to_disk)
    return;
  
  int whichparticle, whichtask;
  void *particle_message;
  
  if(thistask == root)
  {
    if(count_saves%1 == 0)
      printf("# Saving particle to disk. N= %d.\n", count_saves);

    whichtask = gsl_rng_uniform_int(dnest_gsl_r,totaltask);
  }

  MPI_Bcast(&whichtask, 1, MPI_INT, root, MPI_COMM_WORLD);

  if(whichtask != root)
  {
    if(thistask == whichtask)
    {
      int size_message = size_of_modeltype + 2*sizeof(int) + 2*sizeof(double);
      particle_message = (void *)malloc(size_message);
      whichparticle = gsl_rng_uniform_int(dnest_gsl_r,options.num_particles);
      memcpy(particle_message, particles + whichparticle * particle_offset_size, size_of_modeltype);
      memcpy(particle_message + size_of_modeltype, &log_likelihoods[whichparticle].value, 2*sizeof(double));
      memcpy(particle_message + size_of_modeltype + 2*sizeof(double), &level_assignments[whichparticle], sizeof(int));
      memcpy(particle_message + size_of_modeltype + 2*sizeof(double) + sizeof(int), &whichparticle, sizeof(int));

      MPI_Send(particle_message, size_message, MPI_BYTE, root, 1, MPI_COMM_WORLD);

      //printf("%f %f\n", log_likelihoods[whichparticle].value, log_likelihoods[whichparticle].tiebreaker);

      free(particle_message);
    }
    if(thistask == root)
    {
      MPI_Status status;
      int size_message = size_of_modeltype + 2*sizeof(int) + 2*sizeof(double);
      int whichlevel;
      LikelihoodType logl;

      particle_message = (void *)malloc(size_message);

      MPI_Recv(particle_message, size_message, MPI_BYTE, whichtask, 1, MPI_COMM_WORLD, &status);

      memcpy(&logl, particle_message + size_of_modeltype, 2*sizeof(double) );
      memcpy(&whichlevel, particle_message + size_of_modeltype + 2*sizeof(double), sizeof(int) );
      memcpy(&whichparticle, particle_message + size_of_modeltype + 2*sizeof(double) + sizeof(int), sizeof(int) );
      
      //printf("%f %f\n", logl.value, logl.tiebreaker);

      print_particle(fsample, particle_message);

      fprintf(fsample_info, "%d %e %f %d\n", whichlevel, 
        logl.value,
        logl.tiebreaker,
        whichtask * options.num_particles + whichparticle);

      free(particle_message);
    }
  }
  else
  {
    if(thistask == root)
    {
      whichparticle =  gsl_rng_uniform_int(dnest_gsl_r,options.num_particles);

      print_particle(fsample, particles + whichparticle * particle_offset_size);

      fprintf(fsample_info, "%d %e %f %d\n", level_assignments[whichparticle], 
        log_likelihoods[whichparticle].value,
        log_likelihoods[whichparticle].tiebreaker,
        whichtask * options.num_particles + whichparticle);
    }
  }
}

void dnest_mcmc_run()
{
  unsigned int which;
  unsigned int i;
  
  for(i = 0; i<options.thread_steps; i++)
  {

    /* randomly select out one particle to update */
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
    
    if( !enough_levels(levels, size_levels)  && levels[size_levels-1].log_likelihood.value < log_likelihoods[which].value)
    {
      above[size_above] = log_likelihoods[which];
      size_above++;
    }
    
  }
}


void update_particle(unsigned int which)
{
  void *particle = particles+ which*particle_offset_size;
  LikelihoodType *logl = &(log_likelihoods[which]);
  
  Level *level = &(levels[level_assignments[which]]);

  void *proposal = create_model();
  LikelihoodType logl_proposal;
  double log_H;

  copy_model(proposal, particle);
  log_H = perturb(proposal);
  
  logl_proposal.value = log_likelihoods_cal(proposal);
  logl_proposal.tiebreaker =  logl_proposal.tiebreaker + gsl_rng_uniform(dnest_gsl_r);
  wrap(&logl_proposal.tiebreaker, 0.0, 1.0);

  if(log_H > 0.0)
    log_H = 0.0;

  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_H) && level->log_likelihood.value < logl_proposal.value)
  {
    copy_model(particle, proposal);
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
  free(proposal);
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
  if(enough_levels(levels, size_levels))
    return 0.0;

  int i = which_level - (size_levels - 1);
  return i/options.lambda;
}

bool enough_levels(Level *l, int size_l)
{
  int i;

  if(options.max_num_levels == 0)
  {
    if(size_l < 10)
      return false;

    int num_levels_to_check = 20;
    if(size_l > 80)
      num_levels_to_check = (int)(sqrt(20) * sqrt(0.25*size_l));

    int k = size_l - 1;

    for(i= 0; i<num_levels_to_check; i++)
    {
      if(l[k].log_likelihood.value - l[k-1].log_likelihood.value
        >=0.8)
        return false;

      k--;
      if( k < 1 )
        break;
    }
    return true;
  }
  return (size_l >= options.max_num_levels);
}

void initialize_output_file()
{
  if(thistask != root)
    return;

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

void close_output_file()
{
  if(thistask != root )
    return;

  fclose(fsample);
  fclose(fsample_info);
}

void setup(int argc, char** argv)
{
  int i;

  // root task.
  root = 0;
  
  // random number generator
  dnest_gsl_T = (gsl_rng_type *) gsl_rng_default;
  dnest_gsl_r = gsl_rng_alloc (dnest_gsl_T);
  gsl_rng_set(dnest_gsl_r, time(NULL) + thistask);

  // read options
  if(thistask == root)
    options_load();
  MPI_Bcast(&options, sizeof(Options), MPI_BYTE, root, MPI_COMM_WORLD);
  
  compression = exp(1.0);
  regularisation = options.new_level_interval;
  save_to_disk = true;

  // particles
  particle_offset_size = size_of_modeltype/sizeof(void);
  particles = (void *)malloc(options.num_particles*size_of_modeltype);
  
  // initialise sampler
  if(thistask == root)
    all_above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  log_likelihoods = (LikelihoodType *)malloc(2*options.num_particles * sizeof(LikelihoodType));
  level_assignments = (unsigned int*)malloc(options.num_particles * sizeof(unsigned int));

  if(options.max_num_levels != 0)
  {
    levels = (Level *)malloc(options.max_num_levels * sizeof(Level));
    if(thistask == root)
    {
      levels_combine = (Level *)malloc(options.max_num_levels * sizeof(Level));
      copies_of_levels = (Level *)malloc(totaltask * options.max_num_levels * sizeof(Level));
    }
  }
  else
  {
    levels = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));
    if(thistask == root)
    {
      levels_combine = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));
      copies_of_levels = (Level *)malloc(totaltask * LEVEL_NUM_MAX * sizeof(Level));
    }
  }
  

  count_mcmc_steps = 0;
  count_saves = 0;

// first level
  size_levels = 0;
  size_above = 0;
  size_all_above = 0;
  LikelihoodType like_tmp = {-DBL_MAX, gsl_rng_uniform(dnest_gsl_r)};
  Level level_tmp = {like_tmp, 0.0, 0, 0, 0, 0};
  levels[size_levels] = level_tmp;
  size_levels++;

  if(thistask == root)
  {
    size_levels_combine = 0;
    levels_combine[size_levels_combine] = level_tmp;
    size_levels_combine++;
  }
  
  for(i=0; i<options.num_particles; i++)
  {
    from_prior(particles+i*particle_offset_size);
    log_likelihoods[i].value = log_likelihoods_cal(particles+i*particle_offset_size);
    log_likelihoods[i].tiebreaker = dnest_rand();
    level_assignments[i] = 0;
  }
  
  /*ModelType proposal;
  printf("%f %f %f \n", particles[0].param[0], particles[0].param[1], particles[0].param[2] );
  proposal = particles[0];
  printf("%f %f %f \n", proposal.param[0], proposal.param[1], proposal.param[2] );*/
}

void finalise()
{
  free(particles);
  free(above);
  free(log_likelihoods);
  free(level_assignments);
  free(levels);
  if(thistask == root)
  {
    free(all_above);
    free(copies_of_levels);
  }
  gsl_rng_free(dnest_gsl_r);

  if(thistask == root)
    printf("# Finalizing dnest.\n");
}


void options_load()
{
  FILE *fp;
  char buf[BUF_MAX_LENGTH];

  fp = fopen(options_file, "r");

  if(fp == NULL)
  {
    fprintf(stderr, "# ERROR: Cannot open file %s.\n", options_file);
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

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.sample_file);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.sample_info_file);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.levels_file);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.sampler_state_file);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.posterior_sample_file);
  
  fclose(fp);

/*  strcpy(options.sample_file, "sample.txt");
  strcpy(options.sample_info_file, "sample_info.txt");
  strcpy(options.levels_file, "levels.txt");
  strcpy(options.sampler_state_file, "sampler_state.txt");*/
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
