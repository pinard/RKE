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



typedef struct struct_rke_variables
  {

    /* The following are saved from rke_init call arguments. */

    int n_equations;		/* Number of simultaneous equations */
    int (*eval_routine) ();	/* Routine to compute derivatives */

    /* These may be changed by the user between two solve calls. */

    double minimum_step;	/* Minimum allowable step size */
    double maximum_step;	/* Maximum allowable step size */
    double current_step;	/* Current integration step size */
    double error_slope;		/* Slope of maximum error per time unit */
    double error_biais;		/* Biais of maximum error per time unit */
    int accepted_steps;		/* Accumulated number of accepted steps */
    int rejected_steps;		/* Accumulated number of rejected steps */
  }
*rke_variables;


extern rke_variables rke_init (int, int (*) ());
extern void rke_term (rke_variables);
extern int rke_solve (rke_variables, double *, double[], double);
