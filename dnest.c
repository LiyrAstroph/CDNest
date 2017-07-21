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
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "dnestvars.h"

int dnest(int argc, char** argv)
{
  int opt;
  
  setup(argc, argv);
  
  // cope with argv
  if(thistask == root )
  {
    post_temp = 1.0;
    flag_restart = 0;
    flag_postprc = 0;
    flag_sample_info = 0;
    
    //int i;
    //for(i=0; i<argc; i++)
    //{
    //  printf("%s\n", argv[i]);
    //}

    opterr = 0;
    optind = 0;
    while( (opt = getopt(argc, argv, "r:s:pt:c")) != -1)
    {
      switch(opt)
      {
        case 'r':
          flag_restart = 1;
          strcpy(file_restart, optarg);
          printf("# Dnest restarts.\n");
          break;
        case 's':
          strcpy(file_save_restart, optarg);
          printf("# Dnest sets restart file %s.\n", file_save_restart);
          break;
        case 'p':
          flag_postprc = 1;
          post_temp = 1.0;
          printf("# Dnest does postprocess.\n");
          break;
        case 't':
          post_temp = atof(optarg);
          printf("# Dnest sets a temperature %f.\n", post_temp);
          if(post_temp == 0.0)
          {
            printf("# Dnest incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(post_temp < 1.0)
          {
            printf("# Dnest temperature should >= 1.0\n");
            exit(0);
          }
          break;
        case 'c':
          flag_sample_info = 1;
          printf("# Dnest recalculates sample information.\n");
          break;
        case '?':
          printf("# Dnest incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;
        default:
          break;
      }
    }
  }

  MPI_Bcast(&flag_restart, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&flag_postprc, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&flag_sample_info, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&post_temp, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if(flag_postprc == 1)
  {
    dnest_postprocess(post_temp);
    return 0;
  }

  if(flag_sample_info == 1)
  {
    dnest_postprocess(post_temp);
    return 0;
  }

  if(flag_restart==1)
    dnest_restart();

  initialize_output_file();
  dnest_run();
  close_output_file();

  dnest_postprocess(post_temp);

  finalise();

  return 0;
}

// postprocess, calculate evidence, generate posterior sample.
void dnest_postprocess(double temperature)
{
  if(thistask == root)
  {
    options_load();
    postprocess(temperature);
  }
}

void dnest_run()
{
  int i, j, k, size_all_above_incr;
  Level *pl, *levels_orig;
  int *buf_size_above, *buf_displs;
  double *plimits;
  
  // used to gather levels' information
  if(thistask == root)
  {
    buf_size_above = malloc(totaltask * sizeof(int));  
    buf_displs = malloc(totaltask * sizeof(int));
  }

  while(true)
  {
    //check for termination
    if(options.max_num_saves !=0 &&
        count_saves != 0 && (count_saves%options.max_num_saves == 0))
      break;

    dnest_mcmc_run();
    MPI_Barrier(MPI_COMM_WORLD);
    
    //gather levels
    MPI_Gather(levels, size_levels*sizeof(Level), MPI_BYTE, 
             copies_of_levels, size_levels*sizeof(Level), MPI_BYTE, root, MPI_COMM_WORLD);

    //gather limits
    MPI_Gather(limits, size_levels*particle_offset_double*2, MPI_DOUBLE, 
               copies_of_limits, size_levels*particle_offset_double*2, MPI_DOUBLE, root, MPI_COMM_WORLD );
    
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

      // store new points following the end of all_above array.
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

    // reset size_above for each task
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

      // scan over all copies of limits
      plimits = copies_of_limits;
      for(i=0; i< totaltask; i++)
      {
        for(j=0; j < size_levels; j++)
        {
          for(k=0; k<particle_offset_double; k++)
          {
            limits[j * particle_offset_double *2 + k*2 ] = fmin(limits[j * particle_offset_double *2 + k*2 ],
              plimits[j * particle_offset_double *2 + k*2]);
            limits[j * particle_offset_double *2 + k*2 +1 ] = fmax(limits[j * particle_offset_double *2 + k*2 +1 ],
              plimits[j * particle_offset_double *2 + k*2 + 1]);
          }
          
        }
        plimits += size_levels * particle_offset_double * 2;
      }

      // limits of smaller levels should be larger than those of higher levels
      for(j=size_levels-2; j >= 0; j--)
        for(k=0; k<particle_offset_double; k++)
        {
          limits[ j * particle_offset_double *2 + k*2 ] = fmin( limits[ j * particle_offset_double *2 + k*2 ],
                    limits[ (j+1) * particle_offset_double *2 + k*2 ] );

          limits[ j * particle_offset_double *2 + k*2 + 1] = fmax( limits[ j * particle_offset_double *2 + k*2 +1 ],
                    limits[ (j+1) * particle_offset_double *2 + k*2 + 1 ] );
        }

      do_bookkeeping();

      size_levels = size_levels_combine;
      memcpy(levels, levels_combine, size_levels * sizeof(Level));
    }

    //broadcast levels
    MPI_Bcast(&size_levels, 1, MPI_INT, root, MPI_COMM_WORLD);
    //MPI_Bcast(&count_saves, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(levels, size_levels * sizeof(Level), MPI_BYTE, root,  MPI_COMM_WORLD); 
    MPI_Bcast(limits, size_levels * particle_offset_double *2, MPI_DOUBLE, root, MPI_COMM_WORLD);

    if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
    {
      save_particle();

      if(thistask == root )
      {
        // save levels, limits, sync samples when running a number of steps
        if( count_saves % (int)(0.2 * options.max_num_saves) == 0 )
        {
          if(size_levels_combine <= options.max_num_levels)
          {
            save_levels();

            printf("# Save levels at N= %d.\n", count_saves);
          }
          save_limits();
          fflush(fsample_info);
          fsync(fileno(fsample_info));
          fflush(fsample);
          fsync(fileno(fsample));
          printf("# Save limits, and sync samples at N= %d.\n", count_saves);
        }
      }

      if( count_saves % (int)(0.2 * options.max_num_saves) == 0 )
      {
        dnest_save_restart();
      }
    }
  }
  
  //dnest_save_restart();

  if(thistask == root)
  {
    //save levels
    save_levels();
    save_limits();

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
  //bool created_level = false;

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

    //created_level = true;
  }
  
  recalculate_log_X();

/*  if(created_level)
    save_levels();

  if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
  {

    //save_particle();

    if(!created_level)
      save_levels();
  }*/

}

void recalculate_log_X()
{
  int i;

  levels_combine[0].log_X = 0.0;
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

  double max_log_push = -DBL_MAX;

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

void save_limits()
{
  int i, j;
  FILE *fp;

  fp = fopen(options.limits_file, "w");
  for(i=0; i<size_levels_combine; i++)
  {
    fprintf(fp, "%d  ", i);
    for(j=0; j<particle_offset_double; j++)
      fprintf(fp, "%f  %f  ", limits[i*2*particle_offset_double+j*2], limits[i*2*particle_offset_double+j*2+1]);

    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* save particle */
void save_particle()
{
  count_saves++;
  
  which_mcmc_steps = count_saves;

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

    which_particle_update = which;
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
  
  which_level_update = level_assignments[which];
  
  log_H = perturb(proposal);
  
  logl_proposal.value = log_likelihoods_cal(proposal);
  logl_proposal.tiebreaker =  (*logl).tiebreaker + gsl_rng_uniform(dnest_gsl_r);
  wrap(&logl_proposal.tiebreaker, 0.0, 1.0);
  
  if(log_H > 0.0)
    log_H = 0.0;

  perturb_accept[which] = 0;
  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_H) && level->log_likelihood.value < logl_proposal.value)
  {
    copy_model(particle, proposal);
    memcpy(logl, &logl_proposal, sizeof(LikelihoodType));
    level->accepts++;

    perturb_accept[which] = 1;
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
  int i;

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

// update the limits of the level
    double *particle = (double *) (particles+ which*particle_offset_size);
    for(i=0; i<particle_offset_double; i++)
    {
      limits[proposal * 2 * particle_offset_double +  i*2] = 
            fmin(limits[proposal * 2* particle_offset_double +  i*2], particle[i]);
      limits[proposal * 2 * particle_offset_double +  i*2+1] = 
            fmax(limits[proposal * 2 * particle_offset_double +  i*2+1], particle[i]);
    }
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
  return ((double)i)/options.lambda;
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

  if(flag_restart !=1)
    fsample = fopen(options.sample_file, "w");
  else
    fsample = fopen(options.sample_file, "a");
  
  if(fsample==NULL)
  {
    fprintf(stderr, "# Cannot open file sample.txt.\n");
    exit(0);
  }
  if(flag_restart != 1)
    fprintf(fsample, "# \n");

  if(flag_restart != 1)
    fsample_info = fopen(options.sample_info_file, "w");
  else
    fsample_info = fopen(options.sample_info_file, "a");

  if(fsample_info==NULL)
  {
    fprintf(stderr, "# Cannot open file sample_info.txt.\n");
    exit(0);
  }
  if(flag_restart != 1)
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
  int i, j;

  // root task.
  root = 0;

  // random number generator
  dnest_gsl_T = (gsl_rng_type *) gsl_rng_default;
  dnest_gsl_r = gsl_rng_alloc (dnest_gsl_T);
#ifndef Debug
  gsl_rng_set(dnest_gsl_r, time(NULL) + thistask);
#else
  gsl_rng_set(dnest_gsl_r, 9999 + thistask);
  printf("# debugging, task %d dnest random seed %d\n", thistask, 9999 + thistask);
#endif  

  // read options
  if(thistask == root)
    options_load();
  MPI_Bcast(&options, sizeof(Options), MPI_BYTE, root, MPI_COMM_WORLD);
  
  post_temp = 1.0;
  compression = exp(1.0);
  regularisation = options.new_level_interval*0.1;
  save_to_disk = true;

  // particles
  particle_offset_size = size_of_modeltype/sizeof(void);
  particle_offset_double = size_of_modeltype/sizeof(double);
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

    limits = malloc(options.max_num_levels * particle_offset_double * 2 * sizeof(double));
    copies_of_limits = malloc( totaltask * options.max_num_levels * particle_offset_double * 2 * sizeof(double));
    for(i=0; i<options.max_num_levels; i++)
    {
      for(j=0; j<particle_offset_double; j++)
      {
        limits[i*2*particle_offset_double+ j*2] = DBL_MAX;
        limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
      }
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

    limits = malloc(LEVEL_NUM_MAX * particle_offset_double * 2 * sizeof(double));
    copies_of_limits = malloc(totaltask * LEVEL_NUM_MAX * particle_offset_double * 2 * sizeof(double));
    for(i=0; i<LEVEL_NUM_MAX; i++)
    {
      for(j=0; j<particle_offset_double; j++)
      {
        limits[i*2*particle_offset_double + j*2] = DBL_MAX;
        limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
      }
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
  
  which_mcmc_steps = 0;
  for(i=0; i<options.num_particles; i++)
  {
    which_particle_update = i;
    which_level_update = 0;
    from_prior(particles+i*particle_offset_size);
    log_likelihoods[i].value = log_likelihoods_cal_initial(particles+i*particle_offset_size);
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
  free(limits);


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
  char buf[BUF_MAX_LENGTH], buf1[BUF_MAX_LENGTH];

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
    if(sscanf(buf, "%s", buf1) < 1) // a blank line
    {
      buf[0] = '#';
    }
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

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.posterior_sample_info_file);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.limits_file);
  
  fclose(fp);

/*  strcpy(options.sample_file, "sample.txt");
  strcpy(options.sample_info_file, "sample_info.txt");
  strcpy(options.levels_file, "levels.txt");
  strcpy(options.sampler_state_file, "sampler_state.txt");*/
}


double mod(double y, double x)
{
  if(x > 0.0)
  {
    return (y/x - floor(y/x))*x;
  }
  else if(x == 0.0)
  {
    return 0.0;
  }
  else
  {
    printf("Warning in mod(double, double) %e\n", x);
    exit(0);
  }
  
}

void wrap(double *x, double min, double max)
{
  *x = mod(*x - min, max - min) + min;
}

void wrap_limit(double *x, double min, double max)
{

  *x = fmax(fmin(*x, max), min);
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

/*!
 *  Save sampler state for later restart. 
 */
void dnest_save_restart()
{
  FILE *fp;
  int i, j;
  void *particles_all;
  LikelihoodType *log_likelihoods_all;
  unsigned int *level_assignments_all;
  char str[200];

  if(thistask == root)
  {
    sprintf(str, "%s_%d", file_save_restart, count_saves);
    fp = fopen(str, "wb");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s. \n", file_save_restart);
      exit(0);
    }

    particles_all = (void *)malloc( options.num_particles *  totaltask * size_of_modeltype );


    log_likelihoods_all = (LikelihoodType *)malloc(totaltask * options.num_particles * sizeof(LikelihoodType));
    level_assignments_all = (unsigned int*)malloc(totaltask * options.num_particles * sizeof(unsigned int));
  }

  MPI_Gather(particles, options.num_particles * size_of_modeltype, MPI_BYTE, 
    particles_all, options.num_particles * size_of_modeltype, MPI_BYTE, root, MPI_COMM_WORLD);

  MPI_Gather(level_assignments, options.num_particles * sizeof(unsigned int), MPI_BYTE, 
    level_assignments_all, options.num_particles * sizeof(unsigned int), MPI_BYTE, root, MPI_COMM_WORLD);

  MPI_Gather(log_likelihoods, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, 
    log_likelihoods_all, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, root, MPI_COMM_WORLD);


  if(thistask == root )
  {
    printf("# Save restart data to file %s.\n", str);

    //fprintf(fp, "%d %d\n", count_saves, count_mcmc_steps);
    //fprintf(fp, "%d\n", size_levels_combine);

    fwrite(&count_saves, sizeof(int), 1, fp);
    fwrite(&count_mcmc_steps, sizeof(int), 1, fp);
    fwrite(&size_levels_combine, sizeof(int), 1, fp);

    for(i=0; i<size_levels_combine; i++)
    {
      //fprintf(fp, "%14.12g %14.12g %f %llu %llu %llu %llu\n", levels_combine[i].log_X, levels_combine[i].log_likelihood.value, 
      //  levels_combine[i].log_likelihood.tiebreaker, levels_combine[i].accepts,
      //  levels_combine[i].tries, levels[i].exceeds, levels_combine[i].visits);

      fwrite(&levels_combine[i], sizeof(Level), 1, fp);
    }

    for(j=0; j<totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        //fprintf(fp, "%d %e %f\n", level_assignments_all[j*options.num_particles + i], 
        //  log_likelihoods_all[j*options.num_particles + i].value,
        //  log_likelihoods_all[j*options.num_particles + i].tiebreaker);  

        fwrite(&level_assignments_all[j*options.num_particles + i], sizeof(int), 1, fp);
        fwrite(&log_likelihoods_all[j*options.num_particles + i], sizeof(LikelihoodType), 1, fp);    
      }
    }
    
    for(i=0; i<size_levels_combine; i++)
    {
      //fprintf(fp, "%d  ", i);
      for(j=0; j<particle_offset_double; j++)
      {
        //fprintf(fp, "%f  %f  ", limits[i*2*particle_offset_double+j*2], limits[i*2*particle_offset_double+j*2+1]);
        fwrite(&limits[i*2*particle_offset_double+j*2], sizeof(double), 2, fp);
      }

      //fprintf(fp, "\n");
    }
    
    for(j=0; j<totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        //print_particle(fp, particles_all + (j * options.num_particles + i) * particle_offset_size);
        fwrite(particles_all + (j * options.num_particles + i) * particle_offset_size, size_of_modeltype, 1, fp);
      } 
    }
    
    fclose(fp);
  }

  restart_clouds(0);
}

void dnest_restart()
{
  FILE *fp;
  int i, j, k, itmp;
  void *particles_all;
  unsigned int *level_assignments_all;
  LikelihoodType *log_likelihoods_all;
  void *particle;
  char  buf[200];

  if(thistask == root)
  {
    fp = fopen(file_restart, "rb");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s. \n", file_restart);
      exit(0);
    }

    printf("# Reading %s\n", file_restart);

    particles_all = (void *)malloc( options.num_particles *  totaltask * size_of_modeltype );
    log_likelihoods_all = (LikelihoodType *)malloc(totaltask * options.num_particles * sizeof(LikelihoodType));
    level_assignments_all = (unsigned int*)malloc(totaltask * options.num_particles * sizeof(unsigned int));

    //fscanf(fp, "%d %lld\n", &count_saves, &count_mcmc_steps);
    // number of levels
    //fscanf(fp, "%d\n", &size_levels_combine);

    fread(&count_saves, sizeof(int), 1, fp);
    fread(&count_mcmc_steps, sizeof(int), 1, fp);
    fread(&size_levels_combine, sizeof(int), 1, fp);

    //printf("%d %d %d\n", count_saves, count_mcmc_steps, size_levels_combine);

    // read levels
    for(i=0; i<size_levels_combine; i++)
    {
      //fscanf(fp, "%lf %lf %lf %lld %lld %lld %lld\n", &(levels_combine[i].log_X), &(levels_combine[i].log_likelihood.value), 
      //  &(levels_combine[i].log_likelihood.tiebreaker), &(levels_combine[i].accepts),
      //  &(levels_combine[i].tries), &(levels_combine[i].exceeds), &(levels_combine[i].visits) );

      fread(&levels_combine[i], sizeof(Level), 1, fp);

      //printf("%d %d %d %d\n", levels_combine[i].accepts, levels_combine[i].tries, levels_combine[i].exceeds, levels_combine[i].visits);
    }
    size_levels = size_levels_combine;
    memcpy(levels, levels_combine, size_levels * sizeof(Level));

    // read level assignment

    for(j=0; j<totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        //fscanf(fp, "%d %lf %lf\n", &(level_assignments_all[j*options.num_particles + i]), 
        //  &(log_likelihoods_all[j*options.num_particles + i].value), 
        //  &(log_likelihoods_all[j*options.num_particles + i].tiebreaker));

        fread(&level_assignments_all[j*options.num_particles + i], sizeof(int), 1, fp);
        fread(&log_likelihoods_all[j*options.num_particles + i], sizeof(LikelihoodType), 1, fp); 
      }
    }

    // read limits
    for(i=0; i<size_levels; i++)
    {
      //fscanf(fp, "%d", &itmp);
      for(j=0; j<particle_offset_double; j++)
      {
        //fscanf(fp, "%lf  %lf", &limits[i*2*particle_offset_double+j*2], &limits[i*2*particle_offset_double+j*2+1]);
        fread(&limits[i*2*particle_offset_double+j*2], sizeof(double), 2, fp);      
      }

      //fscanf(fp, "\n");
    }

    // read particles
    for(j=0; j<totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        particle = (particles_all + (j * options.num_particles + i) * size_of_modeltype);

        //for(k=0; k < particle_offset_double; k++)
        //{
        //  fscanf(fp, "%lf", &particle[k]);
        //}
        //fscanf(fp, "\n");
        fread(particle, size_of_modeltype, 1, fp);
      }
    }

    fclose(fp);
  }

  MPI_Bcast(&count_saves, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&count_mcmc_steps, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&size_levels, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(levels, size_levels * sizeof(Level), MPI_BYTE, root,  MPI_COMM_WORLD); 

  MPI_Bcast(limits, size_levels * particle_offset_double * 2, MPI_DOUBLE, root, MPI_COMM_WORLD);

  MPI_Scatter(level_assignments_all, options.num_particles * sizeof(unsigned int), MPI_BYTE,
      level_assignments, options.num_particles * sizeof(unsigned int), MPI_BYTE, root, MPI_COMM_WORLD);

  MPI_Scatter(log_likelihoods_all, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, 
      log_likelihoods, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, root, MPI_COMM_WORLD);

  MPI_Scatter(particles_all, options.num_particles * size_of_modeltype, MPI_BYTE, 
    particles, options.num_particles * size_of_modeltype, MPI_BYTE, root, MPI_COMM_WORLD);

  
  restart_clouds(1);

  for(i=0; i<options.num_particles; i++)
  {
    which_particle_update = i;
    which_level_update = level_assignments[i];
    //printf("%d %d %f\n", thistask, i, log_likelihoods[i].value);
    log_likelihoods[i].value = log_likelihoods_cal_restart(particles+i*particle_offset_size);
    //printf("%d %d %f\n", thistask, i, log_likelihoods[i].value);
    //due to randomness, the original level assignment may be incorrect. re-asign the level
    while(log_likelihoods[i].value < levels[level_assignments[i]].log_likelihood.value)
    {
      printf("level assignment decrease %d %f %f %d.\n", i, log_likelihoods[i].value, 
        levels[level_assignments[i]].log_likelihood.value, level_assignments[i]);
      level_assignments[i]--;
    }
    
  }
  return;
}


