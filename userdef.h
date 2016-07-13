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

/*===================================================*/
// users responsible for following struct definitions

/*  data type */
typedef struct 
{
  double x;
  double y;
}DataType;

/* model type */
#define num_params (20)
typedef struct
{
  double params[num_params];
}ModelType;
/*===================================================*/

/* szie of modeltype, which is used for dnest */
extern int size_of_modeltype;

/* data storage */
extern int num_data_points;
extern DataType *data;


#endif
