/*
 * C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * Jun 30, 2016
 *
 */

#ifndef _USERDEF_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>

/* random number generator */
extern const gsl_rng * gsl_r;

/*===========================================*/
// users responsible for following struct definitions
// data 
typedef struct 
{
  double x;
  double y;
}DataType;

// model
#define num_params (20)
typedef struct
{
  double params[num_params];
}ModelType;

/*==========================================*/

/* data */
extern DataType *data;
extern int num_data_points;

#endif
