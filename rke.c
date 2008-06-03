/* Ordinary differential equation solver, Runge-Kutta-England technique.
   Copyright © 1988 Free Software Foundation, Inc.
   François Pinard <pinard@iro.umontreal.ca>, 1988.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/



/* Global declarations, initialization and termination. */


#include "rke.h"


extern double fabs ();
extern char *malloc ();
extern char *alloca ();



/* Initialize a new system of equations. */

rke_variables
rke_init (number, routine)	/* Newly allocated reentrancy block */
     int number;		/* Number of simultaneous equations */
     int (*routine) ();		/* Problem function routine */
{
  rke_variables var;

  var = (rke_variables) malloc (sizeof (struct struct_rke_variables));

  var->n_equations = number;
  var->eval_routine = routine;
  var->minimum_step = 0.000001;
  var->maximum_step = 1000000.0;
  var->current_step = 1.0;
  var->error_slope = 0.0000001;
  var->error_biais = 0.00000001;
  var->accepted_steps = 0;
  var->rejected_steps = 0;

  return var;
}



/* Terminate a set of equations. */

void
rke_term (var)
     rke_variables var;		/* Reentrency block */
{
  free (var);
}


/* Main routine of the module, ODE solver. */



/* Perform a consistent move of time in the system. */

int				/* !0 if success */
rke_solve (var, time, variables, aimed_time)
     rke_variables var;		/* Reentrency block */
     double *time;		/* Current value of time */
     double variables[];	/* Current variables */
     double aimed_time;		/* Value of time which is aimed for */
{
  double whole_step;		/* Signed integration step size */
  double quarter_step;		/* 0.25 * whole_step */
  double half_step;		/* 0.50 * whole_step */
  double three_quarter_step;	/* 0.75 * whole_step */
  double estimated_error;	/* Error as estimated by England method */
  double allowable_error;	/* Maximum error that user permits */
  int within_tolerance;		/* Allowable error has not been passed */
  int all_errors_small;		/* All errors within 2% of tolerances */
  int length_of_array;		/* Length of temporary arrays, is bytes */
  int k;			/* Index in various arrays */

  double *dp, *vt, *v, *d;
  double *a1, *a2, *a3, *a4, *a5, *a6, *a7;

  /* Allocate the work arrays. */

  length_of_array = var->n_equations * sizeof (double);

  dp = (double *) alloca (length_of_array);
  vt = (double *) alloca (length_of_array);
  v  = (double *) alloca (length_of_array);
  d  = (double *) alloca (length_of_array);
  a1 = (double *) alloca (length_of_array);
  a2 = (double *) alloca (length_of_array);
  a3 = (double *) alloca (length_of_array);
  a4 = (double *) alloca (length_of_array);
  a5 = (double *) alloca (length_of_array);
  a6 = (double *) alloca (length_of_array);
  a7 = (double *) alloca (length_of_array);

  /* The integration will continue if a minimum step could bring the
     system closer to the time that is aimed for, even if we have to
     overshoot it a little. */

  while (2 * fabs (aimed_time - *time) > var->minimum_step)
    {

      /* Evaluate initial step size and direction. */

      if ((whole_step = aimed_time - *time) > 0.0)
	{
	  if (whole_step > var->current_step)
	    whole_step = var->current_step;
	}
      else
	{
	  if (whole_step < - var->current_step)
	    whole_step = - var->current_step;
	}

      /*  Evaluate initial differentials. */

      if (! (*var->eval_routine) (*time, variables, dp))
	return 0;

      do

	/* Loop integrating at this time point until integration error is
	   within tolerances.  In any case, adjust integration step size. */

	{
	  /* Calculate various step sizes. */

	  quarter_step = 0.25 * whole_step;
	  half_step = quarter_step + quarter_step;
	  three_quarter_step = half_step + quarter_step;

	  /* Perform a partial computation for one step of Runge-Kutta
	     4th order integration, as far as necessary to chain it to
	     England method for estimating integration errors. */

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a1[k] = half_step * dp[k];
	      v[k] = variables[k]
		+ 0.5*a1[k];
	    }

	  if (! (*var->eval_routine) (*time + quarter_step, v, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a2[k] = half_step * d[k];
	      v[k] = variables[k]
		+ 0.25 * (a1[k] + a2[k]);
	    }

	  if (! (*var->eval_routine) (*time + quarter_step, v, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a3[k] = half_step * d[k];
	      v[k] = variables[k]
		+ (-a2[k] + a3[k] + a3[k]);
	    }

	  if (! (*var->eval_routine) (*time + half_step, v, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a4[k] = half_step * d[k];
	      vt[k] = variables[k]
		+ (a1[k] + 4.0*a3[k] + a4[k]) / 6.0;
	    }

	  if (! (*var->eval_routine) (*time + half_step, vt, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a5[k] = half_step * d[k];
	      v[k] = vt[k]
		+ 0.5*a5[k];
	    }

	  if (! (*var->eval_routine) (*time + three_quarter_step, v, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a6[k] = half_step * d[k];
	      v[k] = vt[k]
		+ 0.25*(a5[k] + a6[k]);
	    }

	  if (! (*var->eval_routine) (*time + three_quarter_step, v, d))
	    return 0;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      a7[k] = half_step * d[k];
	      v[k] = variables[k]
		+ (-a1[k] - 96.0*a2[k] + 92.0*a3[k] - 121.0*a4[k]
		   + 144.0*a5[k] + 6.0*a6[k] - 12.0*a7[k]) / 6.0;
	    }

	  /* Perform England error analysis on partial integration. */

	  if (! (*var->eval_routine) (*time + whole_step, v, d))
	    return 0;

	  within_tolerance = 1;
	  all_errors_small = 1;

	  for (k = 0; k < var->n_equations; ++k)
	    {
	      estimated_error
		= fabs ((-a1[k] + 4.0*a3[k] + 17.0*a4[k]
			  - 23.0*a5[k] + 4.0*a7[k] - half_step*d[k])
		        / 90.0);
	      allowable_error = fabs (whole_step)
		* (var->error_slope*fabs (vt[k]) + var->error_biais);
	      if (estimated_error > allowable_error)
		{
		  within_tolerance = 0;
		  break;
		}
	      else if (estimated_error > 0.02 * allowable_error)
		all_errors_small = 0;
	    }
	  if (within_tolerance)
	    {
	      ++var->accepted_steps;

	      /* Complete the Runge-Kutta step and return values. */

	      for (k = 0; k < var->n_equations; ++k)
		v[k] = vt[k] + (-a6[k] + a7[k] + a7[k]);

	      if (! (*var->eval_routine) (*time + whole_step, v, d))
		return 0;

	      *time += whole_step;

	      for (k = 0; k < var->n_equations; ++k)
		variables[k] = vt[k]
		  + (a5[k] + 4.0*a7[k] + half_step*d[k]) / 6.0;

	      /* Increment step size if desirable. */

	      if (all_errors_small
		  && fabs (whole_step) == var->current_step)
		if (2 * var->current_step > var->maximum_step)
		  var->current_step = var->maximum_step;
		else
		  var->current_step *= 2;
	    }
	  else
	    {

	      ++var->rejected_steps;

	      /* Decrement step size if possible. */

	      if (fabs (whole_step) > var->minimum_step)
		{
		  if (var->current_step < 2 * var->minimum_step)
		    var->current_step = var->minimum_step;
		  else
		    var->current_step *= 0.5;
		  if (aimed_time > *time)
		    whole_step = var->current_step;
		  else
		    whole_step = - var->current_step;
		}
	      else
		return 0;	/* Convergence failed */
	    }
	}

      while (!within_tolerance);
    }
  return 1;			/* Convergence succeeded */
}
