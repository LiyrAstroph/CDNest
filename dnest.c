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

double dnest(int argc, char** argv, DNestFptrSet *fptrset, int num_params, char *optfile)
{
  int opt;

  MPI_Comm_rank(MPI_COMM_WORLD, &dnest_thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &dnest_totaltask);
  
  dnest_check_fptrset(fptrset);

  // cope with argv
  if(dnest_thistask == dnest_root )
  {
    dnest_post_temp = 1.0;
    dnest_flag_restart = 0;
    dnest_flag_postprc = 0;
    dnest_flag_sample_info = 0;
    dnest_flag_limits = 0;

    strcpy(file_save_restart, "restart_dnest.txt");
    strcpy(dnest_sample_postfix, "\0");
    
    //int i;
    //for(i=0; i<argc; i++)
    //{
    //  printf("%s\n", argv[i]);
    //}

    opterr = 0;
    optind = 0;
    while( (opt = getopt(argc, argv, "r:s:pt:clx:")) != -1)
    {
      switch(opt)
      {
        case 'r':
          dnest_flag_restart = 1;
          strcpy(file_restart, optarg);
          printf("# Dnest restarts.\n");
          break;
        case 's':
          strcpy(file_save_restart, optarg);
          printf("# Dnest sets restart file %s.\n", file_save_restart);
          break;
        case 'p':
          dnest_flag_postprc = 1;
          dnest_post_temp = 1.0;
          printf("# Dnest does postprocess.\n");
          break;
        case 't':
          dnest_post_temp = atof(optarg);
          printf("# Dnest sets a temperature %f.\n", dnest_post_temp);
          if(dnest_post_temp == 0.0)
          {
            printf("# Dnest incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(dnest_post_temp < 1.0)
          {
            printf("# Dnest temperature should >= 1.0\n");
            exit(0);
          }
          break;
        case 'c':
          dnest_flag_sample_info = 1;
          printf("# Dnest recalculates sample information.\n");
          break;
        case 'l':
          dnest_flag_limits = 1;
          printf("# Dnest level-dependent sampling.\n");
          break;
        case 'x':
          strcpy(dnest_sample_postfix, optarg);
          printf("# Dnest sets sample postfix %s.\n", dnest_sample_postfix);
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
  
  MPI_Bcast(&dnest_flag_restart, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&dnest_flag_postprc, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&dnest_flag_sample_info, 1, MPI_INT,dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&dnest_post_temp, 1, MPI_DOUBLE, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&dnest_flag_limits, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);

  setup(argc, argv, fptrset, num_params, optfile);

  if(dnest_flag_postprc == 1)
  {
    dnest_postprocess(dnest_post_temp);
    MPI_Barrier(MPI_COMM_WORLD);
    finalise();
    return post_logz;
  }

  if(dnest_flag_sample_info == 1)
  {
    dnest_postprocess(dnest_post_temp);
    finalise();
    return post_logz;
  }

  if(dnest_flag_restart==1)
    dnest_restart();

  initialize_output_file();
  dnest_run();
  close_output_file();

  dnest_postprocess(dnest_post_temp);

  finalise();
  
  return post_logz;
}

// postprocess, calculate evidence, generate posterior sample.
void dnest_postprocess(double temperature)
{
  if(dnest_thistask == dnest_root)
  {
    options_load();
    postprocess(temperature);
  }
  MPI_Bcast(&post_logz, 1, MPI_DOUBLE, dnest_root, MPI_COMM_WORLD);
}

void dnest_run()
{
  int i, j, k, size_all_above_incr;
  Level *pl, *levels_orig;
  int *buf_size_above, *buf_displs;
  double *plimits;
  
  // used to gather levels' information
  if(dnest_thistask == dnest_root)
  {
    buf_size_above = malloc(dnest_totaltask * sizeof(int));  
    buf_displs = malloc(dnest_totaltask * sizeof(int));
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
             copies_of_levels, size_levels*sizeof(Level), MPI_BYTE, dnest_root, MPI_COMM_WORLD);

    //gather limits
    if(dnest_flag_limits == 1)
    {
      MPI_Gather(limits, size_levels*particle_offset_double*2, MPI_DOUBLE, 
               copies_of_limits, size_levels*particle_offset_double*2, MPI_DOUBLE, dnest_root, MPI_COMM_WORLD );
    }
    
    //gather size_above 
    MPI_Gather(&size_above, 1, MPI_INT, buf_size_above, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);

    // task 0 responsible for updating levels
    if(dnest_thistask == dnest_root)
    {
      size_all_above_incr = 0;
      for(i = 0; i<dnest_totaltask; i++)
      {
        size_all_above_incr += buf_size_above[i];

        buf_size_above[i] *= sizeof(LikelihoodType);
      }

      // store new points following the end of all_above array.
      buf_displs[0] = size_all_above * sizeof(LikelihoodType);
      for(i=1; i<dnest_totaltask; i++)
      {
        buf_displs[i] = buf_displs[i-1] + buf_size_above[i-1];
      }
      
      //update size_all_above
      size_all_above += size_all_above_incr; 

      if(size_all_above > options.new_level_interval*2)
      {
        printf("# Error, all above overflow.\n");
        exit(0);
      }
    }
    

    // gather above into all_above, stored in task 0, note that its size is different among tasks
    MPI_Gatherv(above, size_above * sizeof(LikelihoodType), MPI_BYTE, 
                all_above, buf_size_above, buf_displs, MPI_BYTE, dnest_root, MPI_COMM_WORLD);

    // reset size_above for each task
    size_above = 0;

    count_mcmc_steps += options.thread_steps * dnest_totaltask;

    if(dnest_thistask == dnest_root)
    {
      //backup levels_combine
      levels_orig = malloc(size_levels_combine * sizeof(Level));
      memcpy(levels_orig, levels_combine, size_levels_combine*sizeof(Level));

      //scan over all copies of levels
      pl = copies_of_levels;
      for(i=0; i< dnest_totaltask; i++)
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
      if(dnest_flag_limits == 1)
      {
        plimits = copies_of_limits;
        for(i=0; i< dnest_totaltask; i++)
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
      }

      do_bookkeeping();

      size_levels = size_levels_combine;
      memcpy(levels, levels_combine, size_levels * sizeof(Level));
    }

    //broadcast levels
    MPI_Bcast(&size_levels, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
    //MPI_Bcast(&count_saves, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(levels, size_levels * sizeof(Level), MPI_BYTE, dnest_root,  MPI_COMM_WORLD); 

    if(dnest_flag_limits == 1)
      MPI_Bcast(limits, size_levels * particle_offset_double *2, MPI_DOUBLE, dnest_root, MPI_COMM_WORLD);

    if(count_mcmc_steps >= (count_saves + 1)*options.save_interval)
    {
      save_particle();

      if(dnest_thistask == dnest_root )
      {
        // save levels, limits, sync samples when running a number of steps
        if( count_saves % num_saves == 0 )
        {
          if(size_levels_combine <= options.max_num_levels)
          {
            save_levels();

            printf("# Save levels at N= %d.\n", count_saves);
          }
          if(dnest_flag_limits == 1)
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

  if(dnest_thistask == dnest_root)
  {
    //save levels
    save_levels();
    if(dnest_flag_limits == 1)
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
    qsort(all_above, size_all_above, sizeof(LikelihoodType), dnest_cmp);
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

  double kill_probability = 0.0;
  unsigned int num_bad = 0;
  size_t i;

  for(i=0; i<options.num_particles; i++)good[i] = true;

  for(i=0; i<options.num_particles; i++)
  {
    if( log_push(level_assignments[i]) > max_log_push)
      max_log_push = log_push(level_assignments[i]);

    kill_probability = pow(1.0 - 1.0/(1.0 + exp(-log_push(level_assignments[i]) - 4.0)), 3);
    if(gsl_rng_uniform(dnest_gsl_r) <= kill_probability)
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
        }while(!good[i_copy] || gsl_rng_uniform(dnest_gsl_r) >= exp(log_push(level_assignments[i_copy]) - max_log_push));

        memcpy(particles+i*particle_offset_size, particles + i_copy*particle_offset_size, dnest_size_of_modeltype);
        log_likelihoods[i] = log_likelihoods[i_copy];
        level_assignments[i] = level_assignments[i_copy];
         
        kill_action(i, i_copy);

        deletions++;

        printf("# Replacing lagging particle.\n");
        printf("# This has happened %d times.\n", deletions);
      }
    }
  }
  else
    printf("# Warning: all particles lagging!.\n");

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

  /* update state of sampler */
  fp = fopen(options.sampler_state_file, "w");
  fprintf(fp, "%d %d\n", size_levels, count_saves);
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

  if(!save_to_disk)
    return;
  
  int whichparticle, whichtask;
  void *particle_message;
  
  if(dnest_thistask == dnest_root)
  {
    if(count_saves%1 == 0)
      printf("# Saving particle to disk. N= %d.\n", count_saves);

    whichtask = gsl_rng_uniform_int(dnest_gsl_r,dnest_totaltask);
  }

  MPI_Bcast(&whichtask, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);

  if(whichtask != dnest_root)
  {
    if(dnest_thistask == whichtask)
    {
      int size_message = dnest_size_of_modeltype + 2*sizeof(int) + 2*sizeof(double);
      particle_message = (void *)malloc(size_message);
      whichparticle = gsl_rng_uniform_int(dnest_gsl_r,options.num_particles);
      memcpy(particle_message, particles + whichparticle * particle_offset_size, dnest_size_of_modeltype);
      memcpy(particle_message + dnest_size_of_modeltype, &log_likelihoods[whichparticle].value, 2*sizeof(double));
      memcpy(particle_message + dnest_size_of_modeltype + 2*sizeof(double), &level_assignments[whichparticle], sizeof(int));
      memcpy(particle_message + dnest_size_of_modeltype + 2*sizeof(double) + sizeof(int), &whichparticle, sizeof(int));

      MPI_Send(particle_message, size_message, MPI_BYTE, dnest_root, 1, MPI_COMM_WORLD);

      //printf("%f %f\n", log_likelihoods[whichparticle].value, log_likelihoods[whichparticle].tiebreaker);

      free(particle_message);
    }
    if(dnest_thistask == dnest_root)
    {
      MPI_Status status;
      int size_message = dnest_size_of_modeltype + 2*sizeof(int) + 2*sizeof(double);
      int whichlevel;
      LikelihoodType logl;

      particle_message = (void *)malloc(size_message);

      MPI_Recv(particle_message, size_message, MPI_BYTE, whichtask, 1, MPI_COMM_WORLD, &status);

      memcpy(&logl, particle_message + dnest_size_of_modeltype, 2*sizeof(double) );
      memcpy(&whichlevel, particle_message + dnest_size_of_modeltype + 2*sizeof(double), sizeof(int) );
      memcpy(&whichparticle, particle_message + dnest_size_of_modeltype + 2*sizeof(double) + sizeof(int), sizeof(int) );
      
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
    if(dnest_thistask == dnest_root)
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

    dnest_which_particle_update = which;

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

  void *proposal = (void *)malloc(dnest_size_of_modeltype);
  LikelihoodType logl_proposal;
  double log_H;

  memcpy(proposal, particle, dnest_size_of_modeltype);
  dnest_which_level_update = level_assignments[which];
  
  log_H = perturb(proposal);
  
  logl_proposal.value = log_likelihoods_cal(proposal);
  logl_proposal.tiebreaker =  (*logl).tiebreaker + gsl_rng_uniform(dnest_gsl_r);
  dnest_wrap(&logl_proposal.tiebreaker, 0.0, 1.0);
  
  if(log_H > 0.0)
    log_H = 0.0;

  dnest_perturb_accept[which] = 0;
  if( gsl_rng_uniform(dnest_gsl_r) <= exp(log_H) && level->log_likelihood.value < logl_proposal.value)
  {
    memcpy(particle, proposal, dnest_size_of_modeltype);
    memcpy(logl, &logl_proposal, sizeof(LikelihoodType));
    level->accepts++;

    dnest_perturb_accept[which] = 1;
    accept_action();
    account_unaccepts[which] = 0; /* reset the number of unaccepted perturb */
  }
  else 
  {
    account_unaccepts[which] += 1; /* number of unaccepted perturb */
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
                 + (int)( pow(10.0, 2*gsl_rng_uniform(dnest_gsl_r))*gsl_ran_ugaussian(dnest_gsl_r));

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
    if(dnest_flag_limits == 1)
    {
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

}

double log_push(unsigned int which_level)
{
  if(which_level > size_levels)
  {
    printf("level overflow %d %d.\n", which_level, size_levels);
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
  if(dnest_thistask != dnest_root)
    return;

  if(dnest_flag_restart !=1)
    fsample = fopen(options.sample_file, "w");
  else
    fsample = fopen(options.sample_file, "a");
  
  if(fsample==NULL)
  {
    fprintf(stderr, "# Cannot open file sample.txt.\n");
    exit(0);
  }
  if(dnest_flag_restart != 1)
    fprintf(fsample, "# \n");

  if(dnest_flag_restart != 1)
    fsample_info = fopen(options.sample_info_file, "w");
  else
    fsample_info = fopen(options.sample_info_file, "a");

  if(fsample_info==NULL)
  {
    fprintf(stderr, "# Cannot open file sample_info.txt.\n");
    exit(0);
  }
  if(dnest_flag_restart != 1)
    fprintf(fsample_info, "# level assignment, log likelihood, tiebreaker, ID.\n");
}

void close_output_file()
{
  if(dnest_thistask != dnest_root )
    return;

  fclose(fsample);
  fclose(fsample_info);
}

void setup(int argc, char** argv, DNestFptrSet *fptrset, int num_params, char *optfile)
{
  int i, j;

  // root task.
  dnest_root = 0;

  // setup function pointers
  from_prior = fptrset->from_prior;
  log_likelihoods_cal = fptrset->log_likelihoods_cal;
  log_likelihoods_cal_initial = fptrset->log_likelihoods_cal_initial;
  log_likelihoods_cal_restart = fptrset->log_likelihoods_cal_restart;
  perturb = fptrset->perturb;
  print_particle = fptrset->print_particle;
  read_particle = fptrset->read_particle;
  restart_action = fptrset->restart_action;
  accept_action = fptrset->accept_action;
  kill_action = fptrset->kill_action;
  strcpy(options_file, optfile);

  // random number generator
  dnest_gsl_T = (gsl_rng_type *) gsl_rng_default;
  dnest_gsl_r = gsl_rng_alloc (dnest_gsl_T);
#ifndef Debug
  gsl_rng_set(dnest_gsl_r, time(NULL) + dnest_thistask);
#else
  gsl_rng_set(dnest_gsl_r, 9999 + dnest_thistask);
  printf("# debugging, task %d dnest random seed %d\n", dnest_thistask, 9999 + dnest_thistask);
#endif  
  
  dnest_num_params = num_params;
  dnest_size_of_modeltype = dnest_num_params * sizeof(double);

  // read options
  if(dnest_thistask == dnest_root)
    options_load();
  MPI_Bcast(&options, sizeof(Options), MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  //dnest_post_temp = 1.0;
  compression = exp(1.0);
  regularisation = options.new_level_interval*sqrt(options.lambda);
  save_to_disk = true;

  // particles
  particle_offset_size = dnest_size_of_modeltype/sizeof(void);
  particle_offset_double = dnest_size_of_modeltype/sizeof(double);
  particles = (void *)malloc(options.num_particles*dnest_size_of_modeltype);
  
  // initialise sampler
  if(dnest_thistask == dnest_root)
    all_above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  above = (LikelihoodType *)malloc(2*options.new_level_interval * sizeof(LikelihoodType));

  log_likelihoods = (LikelihoodType *)malloc(2*options.num_particles * sizeof(LikelihoodType));
  level_assignments = (unsigned int*)malloc(options.num_particles * sizeof(unsigned int));

  account_unaccepts = (unsigned int *)malloc(options.num_particles * sizeof(unsigned int));
  for(i=0; i<options.num_particles; i++)
  {
    account_unaccepts[i] = 0;
  }

  if(options.max_num_levels != 0)
  {
    levels = (Level *)malloc(options.max_num_levels * sizeof(Level));
    if(dnest_thistask == dnest_root)
    {
      levels_combine = (Level *)malloc(options.max_num_levels * sizeof(Level));
      copies_of_levels = (Level *)malloc(dnest_totaltask * options.max_num_levels * sizeof(Level));
    }

    if(dnest_flag_limits == 1)
    {
      limits = malloc(options.max_num_levels * particle_offset_double * 2 * sizeof(double));
      for(i=0; i<options.max_num_levels; i++)
      {
        for(j=0; j<particle_offset_double; j++)
        {
          limits[i*2*particle_offset_double+ j*2] = DBL_MAX;
          limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
        }
      }

      if(dnest_thistask == dnest_root)
      {
        copies_of_limits = malloc( dnest_totaltask * options.max_num_levels * particle_offset_double * 2 * sizeof(double));
      }
    }
  }
  else
  {
    levels = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));
    if(dnest_thistask == dnest_root)
    {
      levels_combine = (Level *)malloc(LEVEL_NUM_MAX * sizeof(Level));
      copies_of_levels = (Level *)malloc(dnest_totaltask * LEVEL_NUM_MAX * sizeof(Level));
    }

    if(dnest_flag_limits == 1)
    {
      limits = malloc(LEVEL_NUM_MAX * particle_offset_double * 2 * sizeof(double));
      for(i=0; i<LEVEL_NUM_MAX; i++)
      {
        for(j=0; j<particle_offset_double; j++)
        {
          limits[i*2*particle_offset_double + j*2] = DBL_MAX;
          limits[i*2*particle_offset_double + j*2 + 1] = -DBL_MAX;
        }
      }

      if(dnest_thistask == dnest_root)
      {
        copies_of_limits = malloc(dnest_totaltask * LEVEL_NUM_MAX * particle_offset_double * 2 * sizeof(double));
      }
    }
  }
  
  dnest_perturb_accept = malloc(options.num_particles * sizeof(int));
  for(i=0; i<options.num_particles; i++)
  {
    dnest_perturb_accept[i] = 0;
  }

  count_mcmc_steps = 0;
  count_saves = 0;
  num_saves = (int)fmax(0.02*options.max_num_saves, 1.0);


// first level
  size_levels = 0;
  size_above = 0;
  size_all_above = 0;
  LikelihoodType like_tmp = {-DBL_MAX, gsl_rng_uniform(dnest_gsl_r)};
  Level level_tmp = {like_tmp, 0.0, 0, 0, 0, 0};
  levels[size_levels] = level_tmp;
  size_levels++;

  if(dnest_thistask == dnest_root)
  {
    size_levels_combine = 0;
    levels_combine[size_levels_combine] = level_tmp;
    size_levels_combine++;
  }
  
  for(i=0; i<options.num_particles; i++)
  {
    dnest_which_particle_update = i;
    dnest_which_level_update = 0;
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

  free(account_unaccepts);

  if(dnest_flag_limits == 1)
    free(limits);


  if(dnest_thistask == dnest_root)
  {
    free(all_above);
    free(levels_combine);
    free(copies_of_levels);
    if(dnest_flag_limits == 1)
      free(copies_of_limits);
  }
  gsl_rng_free(dnest_gsl_r);

  free(dnest_perturb_accept);

  if(dnest_thistask == dnest_root)
    printf("# Finalizing dnest.\n");
}


void options_load()
{
  FILE *fp;
  char buf[BUF_MAX_LENGTH], buf1[BUF_MAX_LENGTH];

  fp = fopen(options_file, "r");

  if(fp == NULL)
  {
    fprintf(stderr, "# ERROR: Cannot open options file %s.\n", options_file);
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
  strcat(options.sample_file, dnest_sample_postfix);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.sample_info_file);
  strcat(options.sample_info_file, dnest_sample_postfix);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.levels_file);
  strcat(options.levels_file, dnest_sample_postfix);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.sampler_state_file);
  strcat(options.sampler_state_file, dnest_sample_postfix);
  
  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.posterior_sample_file);
  strcat(options.posterior_sample_file, dnest_sample_postfix);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.posterior_sample_info_file);
  strcat(options.posterior_sample_info_file, dnest_sample_postfix);

  fgets(buf, BUF_MAX_LENGTH, fp);
  sscanf(buf, "%s", options.limits_file);
  strcat(options.limits_file, dnest_sample_postfix);
  
  fclose(fp);

  // check options.
  
  if(options.new_level_interval < dnest_totaltask * options.thread_steps)
  {
    printf("# incorrect options:\n");
    printf("# new level interval should be equal to or larger than"); 
    printf("  totaltask * thread step.\n");
    exit(0);
  }

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

void dnest_wrap(double *x, double min, double max)
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
  return pow(10.0, 1.5 - 3.0*fabs(gsl_ran_tdist(dnest_gsl_r, 2))) * gsl_ran_ugaussian(dnest_gsl_r);
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
  return gsl_ran_ugaussian(dnest_gsl_r);
}

int dnest_cmp(const void *pa, const void *pb)
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


int dnest_get_size_levels()
{
  return size_levels;
}

int dnest_get_which_level_update()
{
  return dnest_which_level_update;
}

int dnest_get_which_particle_update()
{
  return dnest_which_particle_update;
}

unsigned int dnest_get_which_num_saves()
{
  return num_saves;
}
unsigned int dnest_get_count_saves()
{
  return count_saves;
}

unsigned long long int dnest_get_count_mcmc_steps()
{
  return count_mcmc_steps;
}

/* 
 * version check
 * 
 *  1: greater
 *  0: equal
 * -1: lower
 */
int dnest_check_version(char *version_str)
{
  int major, minor, patch;

  sscanf(version_str, "%d.%d.%d", &major, &minor, &patch);
  
  if(major > DNEST_MAJOR_VERSION)
    return 1;
  if(major < DNEST_MAJOR_VERSION)
    return -1;

  if(minor > DNEST_MINOR_VERSION)
    return 1;
  if(minor < DNEST_MINOR_VERSION)
    return -1;

  if(patch > DNEST_PATCH_VERSION)
    return 1;
  if(patch > DNEST_PATCH_VERSION)
    return -1;

  return 0;
}

void dnest_check_fptrset(DNestFptrSet *fptrset)
{
  if(fptrset->from_prior == NULL)
  {
    printf("\"from_prior\" function is not defined at task %d.\n", dnest_thistask);
    exit(0);
  }

  if(fptrset->print_particle == NULL)
  {
    printf("\"print_particle\" function is not defined at task %d. \
      \nSet to be default function in dnest.\n", dnest_thistask);
    fptrset->print_particle = dnest_print_particle;
  }

  if(fptrset->read_particle == NULL)
  {
    printf("\"read_particle\" function is not defined at task %d. \
      \nSet to be default function in dnest.\n", dnest_thistask);
    fptrset->read_particle = dnest_read_particle;
  }

  if(fptrset->log_likelihoods_cal == NULL)
  {
    printf("\"log_likelihoods_cal\" function is not defined at task %d.\n", dnest_thistask);
    exit(0);
  }

  if(fptrset->log_likelihoods_cal_initial == NULL)
  {
    printf("\"log_likelihoods_cal_initial\" function is not defined at task %d. \
      \nSet to the same as \"log_likelihoods_cal\" function.\n", dnest_thistask);
    fptrset->log_likelihoods_cal_initial = fptrset->log_likelihoods_cal;
  }

  if(fptrset->log_likelihoods_cal_restart == NULL)
  {
    printf("\"log_likelihoods_cal_restart\" function is not defined at task %d. \
      \nSet to the same as \"log_likelihoods_cal\" function.\n", dnest_thistask);
    fptrset->log_likelihoods_cal_restart = fptrset->log_likelihoods_cal;
  }

  if(fptrset->perturb == NULL)
  {
    printf("\"perturb\" function is not defined at task %d.\n", dnest_thistask);
    exit(0);
  }

  if(fptrset->restart_action == NULL)
  {
    printf("\"restart_action\" function is not defined at task %d.\
      \nSet to the default function in dnest.\n", dnest_thistask);
    fptrset->restart_action = dnest_restart_action;
  }

  if(fptrset->accept_action == NULL)
  {
    printf("\"accept_action\" function is not defined at task %d.\
      \nSet to the default function in dnest.\n", dnest_thistask);
    fptrset->accept_action = dnest_accept_action;
  }

  if(fptrset->kill_action == NULL)
  {
    printf("\"kill_action\" function is not defined at task %d.\
      \nSet to the default function in dnest.\n", dnest_thistask);
    fptrset->kill_action = dnest_kill_action;
  }

  return;
}

DNestFptrSet * dnest_malloc_fptrset()
{
  DNestFptrSet * fptrset;
  fptrset = (DNestFptrSet *)malloc(sizeof(DNestFptrSet));

  fptrset->from_prior = NULL;
  fptrset->log_likelihoods_cal = NULL;
  fptrset->log_likelihoods_cal_initial = NULL;
  fptrset->log_likelihoods_cal_restart = NULL;
  fptrset->perturb = NULL;
  fptrset->print_particle = NULL;
  fptrset->read_particle = NULL;
  fptrset->restart_action = NULL;
  fptrset->accept_action = NULL;
  fptrset->kill_action = NULL;
  return fptrset;
}

void dnest_free_fptrset(DNestFptrSet * fptrset)
{
  free(fptrset);
  return;
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

  if(dnest_thistask == dnest_root)
  {
    sprintf(str, "%s_%d", file_save_restart, count_saves);
    fp = fopen(str, "wb");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s. \n", file_save_restart);
      exit(0);
    }

    particles_all = (void *)malloc( options.num_particles *  dnest_totaltask * dnest_size_of_modeltype );


    log_likelihoods_all = (LikelihoodType *)malloc(dnest_totaltask * options.num_particles * sizeof(LikelihoodType));
    level_assignments_all = (unsigned int*)malloc(dnest_totaltask * options.num_particles * sizeof(unsigned int));
  }

  MPI_Gather(particles, options.num_particles * dnest_size_of_modeltype, MPI_BYTE, 
    particles_all, options.num_particles * dnest_size_of_modeltype, MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  MPI_Gather(level_assignments, options.num_particles * sizeof(unsigned int), MPI_BYTE, 
    level_assignments_all, options.num_particles * sizeof(unsigned int), MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  MPI_Gather(log_likelihoods, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, 
    log_likelihoods_all, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, dnest_root, MPI_COMM_WORLD);


  if(dnest_thistask == dnest_root )
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

    for(j=0; j<dnest_totaltask; j++)
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
    
    if(dnest_flag_limits == 1)
    {
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
    }
    
    for(j=0; j<dnest_totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        //print_particle(fp, particles_all + (j * options.num_particles + i) * particle_offset_size);
        fwrite(particles_all + (j * options.num_particles + i) * particle_offset_size, dnest_size_of_modeltype, 1, fp);
      } 
    }
    
    fclose(fp);
    free(particles_all);
    free(log_likelihoods_all);
    free(level_assignments_all);
  }

  restart_action(0);
}

void dnest_restart()
{
  FILE *fp;
  int i, j;
  void *particles_all;
  unsigned int *level_assignments_all;
  LikelihoodType *log_likelihoods_all;
  void *particle;

  if(dnest_thistask == dnest_root)
  {
    fp = fopen(file_restart, "rb");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s. \n", file_restart);
      exit(0);
    }

    printf("# Reading %s\n", file_restart);

    particles_all = (void *)malloc( options.num_particles *  dnest_totaltask * dnest_size_of_modeltype );
    log_likelihoods_all = (LikelihoodType *)malloc(dnest_totaltask * options.num_particles * sizeof(LikelihoodType));
    level_assignments_all = (unsigned int*)malloc(dnest_totaltask * options.num_particles * sizeof(unsigned int));

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

    for(j=0; j<dnest_totaltask; j++)
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
    if(dnest_flag_limits == 1)
    {
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
    }

    // read particles
    for(j=0; j<dnest_totaltask; j++)
    {
      for(i=0; i<options.num_particles; i++)
      {
        particle = (particles_all + (j * options.num_particles + i) * dnest_size_of_modeltype);

        //for(k=0; k < particle_offset_double; k++)
        //{
        //  fscanf(fp, "%lf", &particle[k]);
        //}
        //fscanf(fp, "\n");
        fread(particle, dnest_size_of_modeltype, 1, fp);
      }
    }

    fclose(fp);
  }

  MPI_Bcast(&count_saves, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&count_mcmc_steps, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(&size_levels, 1, MPI_INT, dnest_root, MPI_COMM_WORLD);
  MPI_Bcast(levels, size_levels * sizeof(Level), MPI_BYTE, dnest_root,  MPI_COMM_WORLD); 

  if(dnest_flag_limits == 1)
    MPI_Bcast(limits, size_levels * particle_offset_double * 2, MPI_DOUBLE, dnest_root, MPI_COMM_WORLD);

  MPI_Scatter(level_assignments_all, options.num_particles * sizeof(unsigned int), MPI_BYTE,
      level_assignments, options.num_particles * sizeof(unsigned int), MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  MPI_Scatter(log_likelihoods_all, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, 
      log_likelihoods, options.num_particles * sizeof(LikelihoodType), MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  MPI_Scatter(particles_all, options.num_particles * dnest_size_of_modeltype, MPI_BYTE, 
    particles, options.num_particles * dnest_size_of_modeltype, MPI_BYTE, dnest_root, MPI_COMM_WORLD);

  
  restart_action(1);

  for(i=0; i<options.num_particles; i++)
  {
    dnest_which_particle_update = i;
    dnest_which_level_update = level_assignments[i];
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

  if(dnest_thistask == dnest_root)
  {
    free(particles_all);
    free(log_likelihoods_all);
    free(level_assignments_all);
  }
  return;
}

void dnest_print_particle(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<dnest_num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

void dnest_read_particle(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < dnest_num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file %s.\n", options.sample_file);
      exit(0);
    }
  }
  return;
}

void dnest_restart_action(int iflag)
{
  return;
}

void dnest_accept_action()
{
  return;
}

void dnest_kill_action(int i, int i_copy)
{
  return;
}



