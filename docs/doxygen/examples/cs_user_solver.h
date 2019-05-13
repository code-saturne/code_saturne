/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*-----------------------------------------------------------------------------*/

/*!
  \page user_solver Solving a heat equation with a user solver (cs_user_solver.c)

  \section cs_user_solver_h_intro Introduction

  C user functions for the setting of a user solver. The \ref cs_user_solver
  function allows one to replace the Code_Saturne solver by a solver of his
  own.

  The following example shows the setting of a user solver which solves
  a heat equation with the finite volume method.

  \section cs_user_solver_h_cs_user_solver_set Activating the user solver

  If the function \ref cs_user_solver_set returns 1, the Code_Saturne
  solver is replaced by the solver defined in \ref cs_user_solver .

  \snippet cs_user_solver-heat-equation.c set_solver

  \section cs_user_solver_h_cs_user_solver User solver implementation

  \subsection cs_user_solver_h_variable Variable declaration

  The following variables are declared locally and will be used
  later in the code. Three pointers to \ref cs_real_t are declared
  and will stand for the temperature at the time steps n and n-1, and
  for the analytical solution of the heat equation;

  \snippet cs_user_solver-heat-equation.c local_variables

  \subsection cs_user_solver_h_initialization Initialization

  Three arrays of size \c n are created. Several constants are initialized,
  and the array t_old is initialized as follows :

  \f[
    t_{old}(x = i*dx) = sin\frac{pi(0.5+i)}{n};
  \f]

  \snippet cs_user_solver-heat-equation.c initialization

  \subsection cs_user_solver_h_restart Restart creation

  One can define a restart computation, as for the real Code_Saturne solver.

  \snippet cs_user_solver-heat-equation.c restart

  \subsection cs_user_solver_h_time_monitoring Time monitoring

  Monitoring points can equally be created in order to plot the evolution
  of the temperature over the time.

  \snippet cs_user_solver-heat-equation.c time_monitoring

  \subsection cs_user_solver_h_calculation Heat equation solving

  At each iteration, the new temperature is computed using the Finite Volume
  Method.

  \f[
     t(x = i \cdot dx) = t_{old}(x = i \cdot dx) + r\cdot(t_{old}(x = (i+1) \cdot dx) - 2 \cdot t_{old}(x = i \cdot dx) + t_{old}(x = (i-1) \cdot dx))
  \f]

  The boundary conditions at \f$ x = 0 \f$ and \f$ x = L \f$ are :

  \f[
     t(x = 0) = t_{old}(x = 0) + r \cdot (t_{old}(x = dx) - 3 \cdot t_{old}(x = 0) + 2 \cdot t_0)
  \f]

  \f[
     t(x = (n-1) \cdot dx) = t_{old}(x = (n-1) \cdot dx) + r \cdot (2 \cdot t_L - 3 \cdot t_{old}(x = (n-1) \cdot dx) + t_{old}(x = (n-2) \cdot dx))
  \f]

  t_old is then updated and t_sol computed. Finally, the temperature at \f$ x = \frac{L}{2} \f$ and
  at the current iteration is plotted.

  \snippet cs_user_solver-heat-equation.c calculation

  \subsection cs_user_solver_h_checkpoint Checkpoint creation

  A checkpoint can be created in order to restart the computation later on.

  \snippet cs_user_solver-heat-equation.c checkpoint

  \subsection cs_user_solver_h_post_processing Post-processing

  Finally, t and t_sol are postprocessed in order to be compared.

  \snippet cs_user_solver-heat-equation.c post_processing

  \subsection cs_user_solver_h_post_finalization Finalization

  \snippet cs_user_solver-heat-equation.c finalization

*/
