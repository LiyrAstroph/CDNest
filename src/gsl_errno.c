/* 
 * code for gsl error functions
 * extracted from GSL package (https://www.gnu.org/software/gsl/)
 * Yan-Rong Li, Thu Dec 4, 2025
 * liyanrong@mail.ihep.ac.cn
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "gsl_errno.h"

gsl_error_handler_t * gsl_error_handler = NULL;

static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);

void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler) 
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}

static void
no_error_handler (const char *reason, const char *file, int line, int gsl_errno)
{
  /* do nothing */
  reason = 0;
  file = 0;
  line = 0;
  gsl_errno = 0;
  return;
}

FILE * gsl_stream = NULL ;
gsl_stream_handler_t * gsl_stream_handler = NULL;

void
gsl_stream_printf (const char *label, const char *file, int line, 
                   const char *reason)
{
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}
