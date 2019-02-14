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
  \page cs_user_extra_operations_examples cs_user_extra_operations.f90


  \section cs_user_extra_operations_examples_intro Introduction

  This page provides several examples of code blocks that may be used
  to perform data extraction or modify values
  in \ref cs_user_extra_operations.

  \section cs_user_extra_operations_examples_cs_user_extra_op_examples Extra operations examples
  Here is the list of examples dedicated to different physics:

  - \subpage cs_user_extra_operations_examples_energy_balance_p
  - \subpage cs_user_extra_operations_examples_balance_by_zone_p
  - \subpage cs_user_extra_operations_examples_force_temperature_p
  - \subpage cs_user_extra_operations_examples_global_efforts_p
  - \subpage cs_user_extra_operations_examples_parallel_operations_p
  - \subpage cs_user_extra_operations_examples_stopping_criterion_p

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_energy_balance_p Energy balance

  \section cs_user_extra_operations_examples_energy_balance Energy balance

  \subsection cs_user_extra_operations_examples_loc_var_eb Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_extra_operations-energy_balance.f90 loc_var_dec

  \subsection cs_user_extra_operations_examples_init Initialization and finalization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_extra_operations-energy_balance.f90 init

  Ad the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_extra_operations-energy_balance.f90 finalize

  In theory Fortran 95 deallocates locally-allocated arrays automatically,
  but deallocating arrays in a symmetric manner to their allocation is good
  practice, and it avoids using a different logic for C and Fortran.

  \subsection cs_user_extra_operations_examples_eb_body Body

  This example computes energy balance relative to temperature
  We assume that we want to compute balances  (convective and diffusive)
  at the boundaries of the calculation domain represented below
  (with boundaries marked by colors).

  The scalar considered if the temperature. We will also use the
  specific heat (to obtain balances in Joules)


  Domain and associated boundary colors:
  - 2, 4, 7 : adiabatic walls
  - 6       : wall with fixed temperature
  - 3       : inlet
  - 5       : outlet
  - 1       : symmetry


  To ensure calculations have physical meaning, it is best to use
  a spatially uniform time step (\ref optcal::idtvar "idtvar" = 0 or 1).
  In addition, when restarting a calculation, the balance may be
  incorrect at the first time step.

  Temperature variable
  - ivar = \ref isca(\ref optcal::iscalt "iscalt")

  Boundary coefficients coefap/coefbp are those of
  \ref numvar::ivarfl "ivarfl"(ivar).


  The balance at time step n is equal to:

  \f[
  \begin{array}{r c l}
  Balance^n &=& \displaystyle
               \sum_{\celli=1}^{\ncell}
                  \norm{\vol{\celli}} C_p \rho_\celli
                  \left(T_\celli^{n-1} -T_\celli^n \right)  \\
           &+& \displaystyle
               \sum_{\fib}
                  C_p \Delta t_\celli \norm{\vect{S}_\ib}
                  \left(A_\ib^f + B_\ib^f T_\celli^n \right) \\
           &+& \displaystyle
               \sum_{\fib}
                  C_p \Delta t_\celli \dot{m}_\ib
                  \left(A_\ib^g + B_\ib^g T_\celli^n \right)
  \end{array}
  \f]

  The first term is negative if the amount of energy in the volume
  has increased (it is 0 in a steady regime).

  The other terms (convection, diffusion) are positive if the amount
  of energy in the volume has increased due to boundary conditions.

  In a steady regime, a positive balance thus indicates an energy gain.


  With \f$ \rho \f$ (\c rom) calculated using the density law from the
  \ref usphyv subroutine, for example:

  \f[
  \rho^{n-1}_\celli = P_0 / \left( R T_\celli^{n-1} + T_0 \right)
  \f]
  where \f$ R\f$ is \c cs_physical_constants_r and \f$ T_0 \f$ is \c tkelv.


  \f$ C_p \f$ and \f$ \lambda/C_p \f$ may vary.


  Here is the corresponding code:

  \snippet cs_user_extra_operations-energy_balance.f90 example_1

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_balance_by_zone_p Scalar and head loss balance by zone

  \section cs_user_extra_operations_examples_scalar_balance_by_zone Scalar balance by zone

  This is an example of \ref cs_user_extra_operations which performs scalar
  balances on specified zones.

  \subsection cs_user_extra_operations_examples_bz_body body

  The algorithm implemented in the subroutine balance_by_zone adds up
  contributions of fluxes on the boundary of the sub-domain defined by the user.
  The different contributions are classified according to their nature (inlet,
  outlet, wall, symmetry...) if they are boundary faces of the whole domain or
  according to the sign of the mass flux if they are boundary faces of the
  sub-domain but internal faces of the whole domain.

  To ensure calculations have physical meaning, it is best to use
  a spatially uniform time step (\ref optcal::idtvar "idtvar" = 0 or 1).

  The balance at time step n over a subdomain \f$ \Omega \f$ of boundary
  \f$ \partial \Omega \f$ is equal to:

  \f[
  \begin{array}{r c l}
  Balance^n &=& \displaystyle
                \sum_{\Omega_i \in \Omega}
                  \norm{\vol{\celli}} \rho_\celli
                  \left(\varia_\celli^{n-1} -\varia_\celli^n \right)  \\
           &+& \displaystyle
               \sum_{\Omega_i \in \Omega}
               \sum_{\face \in \Face{\celli}}
                  \Delta t_\celli
                  \varia_\celli^n \left(\rho \vect{u}\right)_\face^n \cdot \vect{S}_{\iface} \\
           &-& \displaystyle
               \sum_{\face \in \partial \Omega}
                  \Delta t_\celli
                  \varia_\face^n \left(\rho \vect{u}\right)_\face^n \cdot \vect{S}_{\iface} \\
           &+& \displaystyle
               \sum_{\face \in \partial \Omega}
                   \Delta t_\celli
                   K_\face \grad_\face \varia^n \cdot \vect{S}_{\iface} \\
           &-& \displaystyle
               \sum_{\fib \in \partial \Omega}
                  \Delta t_\celli \dot{m}_\ib
                  \left(A_\ib^g + B_\ib^g \varia_\celli^n \right) \\
           &+& \displaystyle
               \sum_{\fib \in \partial \Omega}
                  \Delta t_\celli \norm{\vect{S}_\ib}
                  \left(A_\ib^f + B_\ib^f \varia_\celli^n \right)
  \end{array}
  \f]

  The first term is negative if the amount of scalar in the volume
  has increased (it is 0 in a steady regime).

  The terms of convection and diffusion (at internal or boundary faces) are positive
  if the amount of scalar in the volume has increased.

  In a steady regime, a positive balance thus indicates a scalar gain.

  \subsection cs_user_extra_operations_examples_example_1 Example 1

  This example computes energy balance relative to temperature.
  We assume that we want to compute balances  (convective and diffusive)
  at the boundaries of the calculation domain.

  The scalar considered is the temperature, nevertheless it is multiplied by the
  specific heat at each cell so that the computed balance is on energy, hence in Joules.

  Here is the corresponding code:

  \snippet cs_user_extra_operations-balance_by_zone.c example_1

  \subsection cs_user_extra_operations_examples_example_2 Example 2

  This example computes the balance relative to a scalar named "scalar1".
  We assume that we want to compute balances (convective and diffusive)
  on a box defined by two diagonally opposite points (the extrema in terms
  of coordinates).

  The box criterion can be used as follows:
  box[\f$ x_{min}\f$, \f$ y_{min}\f$, \f$ z_{min}\f$, \f$ x_{max}\f$, \f$ y_{max}\f$, \f$ z_{max}\f$].

  Here is the corresponding code:

  \snippet cs_user_extra_operations-balance_by_zone.c example_2

  \subsection cs_user_extra_operations_examples_example_3 Scalar balance through a surface

  This example computes the balance relative to a scalar named "scalar1".
  Here is the corresponding code:

  \snippet cs_user_extra_operations-balance_by_zone.c example_3

  \subsection cs_user_extra_operations_examples_example_4 Specific terms of a scalar balance

  Instead of simply logging the various balance terms, it is possible to access
  them using the lower level functions, and the \ref cs_balance_term_t
  components of the computed balance.

  The following exemple shows how to access for example the mass flow
  components of the scalar balance:
  \snippet cs_user_extra_operations-balance_by_zone.c example_4

  \section cs_user_extra_operations_examples_head_balance_by_zone Head loss balance by zone

  This example computes the head balance for a volume zone.
  Here is the corresponding code:

  \snippet cs_user_extra_operations-balance_by_zone.c example_5

  \subsection cs_user_extra_operations_examples_example_6 Specific terms of a head balance

  Instead of simply logging the various balance terms, it is possible to access
  them using the lower level functions, and the \ref cs_balance_p_term_t
  components of the computed balance.

  The following exemple shows how to access for example the mass flow
  components of the pressure drop computation:
  \snippet cs_user_extra_operations-balance_by_zone.c example_6

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_force_temperature_p Force temperature in a given region

  \section cs_user_extra_operations_examples_force_temperature Force temperature in a given region

  This is an example of \ref cs_user_extra_operations
  which sets temperature to 20 in a given region starting at t = 12s

  \subsection cs_user_extra_operations_examples_loc_var_ft Local variables to be added

  \snippet cs_user_extra_operations-force_temperature.f90 loc_var_dec

  \subsection cs_user_extra_operations_examples_ft_body Body


  Do this with precaution...
  The user is responsible for the validity of results.

  Here is the corresponding code:

  \snippet cs_user_extra_operations-force_temperature.f90 example_1

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_global_efforts_p Global efforts

  \section cs_user_extra_operations_examples_global_efforts Global efforts

  This is an example of \ref cs_user_extra_operations which computes global efforts

  \subsection cs_user_extra_operations_examples_loc_var_geff Local variables to be added

  \snippet cs_user_extra_operations-global_efforts.f90 loc_var_dec

  \subsection cs_user_extra_operations_examples_geff_body Body

  Example: compute global efforts on a subset of faces.
  If boundary stresses have been calculated correctly:

  \snippet cs_user_extra_operations-global_efforts.f90 example_1

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_parallel_operations_p Parallel operations

  \section cs_user_extra_operations_examples_parallel_operations Parallel operations

  This is an example of \ref cs_user_extra_operations which performs parallel operations.

  \subsection cs_user_extra_operations_examples_loc_var_po Local variables to be added

  \snippet cs_user_extra_operations-parallel_operations.f90 loc_var_dec

  \subsection cs_user_extra_operations_examples_example_1_po Example 1

  Sum of an integer counter 'ii', here the number of cells.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_1

  \subsection cs_user_extra_operations_examples_example_2_po Example 2

  Maximum of an integer counter 'ii', here the number of cells.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_2

  \subsection cs_user_extra_operations_examples_example_3_po Example 3

  Sum of a real 'rrr', here the volume.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_3

  \subsection cs_user_extra_operations_examples_example_4_po Example 4

  Minimum of a real 'rrr', here the volume.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_4

  \subsection cs_user_extra_operations_examples_example_5_po Example 5

  Maximum of a real 'rrr', here the volume.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_5

  \subsection cs_user_extra_operations_examples_example_6_po Example 6

  Maximum of a real and associated real values;
  here the volume and its location (3 coordinates).

  \snippet cs_user_extra_operations-parallel_operations.f90 example_6

  \subsection cs_user_extra_operations_examples_example_7_po Example 7

  Minimum of a real and associated real values;
  here the volume and its location (3 coordinates).

  \snippet cs_user_extra_operations-parallel_operations.f90 example_7

  \subsection cs_user_extra_operations_examples_example_8_po Example 8

  Sum of an array of integers;
  here, the number of cells, faces, and boundary faces.

  local values; note that to avoid counting interior faces on
  parallel boundaries twice, we check if 'ifacel(1,ifac) .le. ncel',
  as on a parallel boundary, this is always true for one domain
  and false for the other.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_8

  \subsection cs_user_extra_operations_examples_example_9_po Example 9

  Maxima from an array of integers;
  here, the number of cells, faces, and boundary faces.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_9


  \subsection cs_user_extra_operations_examples_example_10_po Example 10

  Minima from an array of integers;
  here, the number of cells, faces, and boundary faces.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_10

  \subsection cs_user_extra_operations_examples_example_11_po Example 11

  Sum of an array of reals;
  here, the 3 velocity components (so as to compute a mean for example).

  \snippet cs_user_extra_operations-parallel_operations.f90 example_11

  \subsection cs_user_extra_operations_examples_example_12_po Example 12

  Maximum of an array of reals;
  here, the 3 velocity components.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_12

  \subsection cs_user_extra_operations_examples_example_13_po Example 13

  Maximum of an array of reals;
  here, the 3 velocity components.

  \snippet cs_user_extra_operations-parallel_operations.f90 example_13

  \subsection cs_user_extra_operations_examples_example_14_po Example 14

  Broadcast an array of local integers to other ranks;
  in this example, we use the number of cells, interior faces, and boundary
  faces from process rank 0 (irangv).

  \snippet cs_user_extra_operations-parallel_operations.f90 example_14

  \subsection cs_user_extra_operations_examples_example_15_po Example 15

  Broadcast an array of local reals to other ranks;
  in this example, we use 3 velocity values from process rank 0 (irangv).

  \snippet cs_user_extra_operations-parallel_operations.f90 example_15

*/
// __________________________________________________________________________________
/*!

  \page cs_user_extra_operations_examples_stopping_criterion_p Stopping criterion based on L2 time residuals

  \section cs_user_extra_operations_examples_stopping_criterion Stopping criterion based on L2 time residuals

  This is an example of \ref cs_user_extra_operations allowing to properly stop a computation
  when the L2 time residuals (displayed in the run_solver.log file) of all
  solved variables have decreased below a value of 1e-3.

  L2 time residuals of a variable at a given time step are a relative measure of
  the unsteady term of its transport equation:

  \f[
  \sqrt{\int_\Omega \left| \der{\varia}{t} \right|^2 \dd \Omega / \int_\Omega \left| \varia \right|^2 \dd \Omega}
  \f]

  \snippet cs_user_extra_operations-stopping_criterion.c extra_stopping_criterion

*/
