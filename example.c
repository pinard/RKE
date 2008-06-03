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


#include "rke.h"


/* Check how close we can get back to our initial conditions. */

print_return (back, initial)
     double back;		/* Value from solving forth and back */
     double initial;		/* True value before the whole process */
{
  printf ("  returning to %12.6lf, got %12.6lf\n", initial, back);
}



/* Print statistics about number of steps. */

print_steps (var)
     rke_variables var;
{
  printf ("    using %3d accepted and %3d rejected steps\n",
	 var->accepted_steps, var->rejected_steps);
}


/* Integration under a normal curve. */

#include <math.h>

static double example_1_const;	/* 1.0 / sqrt (2 * pi) */


static int
problem_function_1 (t, v, d)
     double t;			/* Abcissa point */
     double v[1];		/* v[0] is the surface under the curve */
     double d[1];		/* d[0] is the normal curve itself */
{
  d[0] = example_1_const * exp (-0.5 * t * t);
  return 1;
}


static
example_1 ()
{
  rke_variables p;
  double t;
  double v[1];

  example_1_const = 1.0 / sqrt (2 * 3.1415926);

  p = rke_init (1, problem_function_1);

  t = -1.0;			/* Start at -1.0 */
  v[0] = 0.0;			/* Surface is 0.0 at this point */

  /* Now, simply move to +1.0, and collect the answer. */

  if (rke_solve (p, &t, v, 1.0))
    printf ("\nProbability	= %12.6lf.\n", v[0]);
  else
    printf ("\nProbability not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, -1.0))
      print_return (v[0], 0.0);
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
}



/* Rediscovering cos and sin. */

static int
problem_function_2 (t, v, d)
     double t;
     double v[2];		/* v[0] is cos t */
				/* v[1] is sin t */
     double d[2];		/* d[0] is d cos t / dt == - sin t */
				/* d[1] is d sin t / dt == cos t */
{
  d[0] = -v[1];
  d[1] = v[0];
  return 1;
}


static
example_2 ()
{
  rke_variables p;
  double t;
  double v[2];

  p = rke_init (2, problem_function_2);

  t = 0.0;			/* Start where we know the values */
  v[0] = 1.0;			/* cos 0 = 1.0 */
  v[1] = 0.0;			/* sin 0 = 0.0 */

  /* Now, simply move to 1.5, and collect the answer. */

  if (rke_solve (p, &t, v, 1.5))
    printf ("\ncos (1.5)	= %12.6lf.\n", v[0]);
  else
    printf ("\ncos (1.5) not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, 0.0))
    {
      print_return (v[0], 1.0);
      print_return (v[1], 0.0);
    }
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
}


/* Box slowing by friction in air. */

static int
problem_function_3 (t, v, d)
     double t;			/* Time so far */
     double v[2];		/* v[0] is the distance so far */
				/* v[1] is the current speed */
     double d[2];		/* d[0] is also the current speed */
				/* d[1] is the current acceleration */
{
  d[0] = v[1];
  d[1] = -0.01 * v[1] * v[1];
  return 1;
}


static
example_3 ()
{
  rke_variables p;
  double t;
  double v[2];

  p = rke_init (2, problem_function_3);

  t = 0.0;			/* Start the clock... */
  v[0] = 0.0;			/* ... with no distance so far */
  v[1] = 100.0;		/* ... but some initial speed */

  /* Now, simply ask the clock to be 5.0, and collect the answer. */

  if (rke_solve (p, &t, v, 5.0))
    printf ("\nDistance	= %12.6lf.\n", v[0]);
  else
    printf ("\nDistance not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, 0.0))
    {
      print_return (v[0], 0.0);
      print_return (v[1], 100.0);
    }
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
}


/* Main program. */

main ()
{
  example_1 ();
  example_2 ();
  example_3 ();
}
