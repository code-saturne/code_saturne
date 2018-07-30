/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
  \page cs_user_boundary_conditions_examples cs_user_boundary_conditions.f90


  \section intro_bc Introduction

  This page provides several examples of code blocks that may be used
  to define boundary conditions in \ref cs_user_boundary_conditions.

  \section cs_user_bc_examples Boundary condition definition examples
  Here is the list of examples dedicated to different physics:

  - \subpage base_examples
  - \subpage channel_inlet_mapped
  - \subpage channel_inlet_auto
  - \subpage advanced_examples
  - \subpage atmospheric_examples
  - \subpage compressible_examples
  - \subpage electric_arcs_examples
  - \subpage electric_arcs_ieljou_3_or_4_examples
  - \subpage fuel_examples
  - \subpage gas_3ptchem_examples
  - \subpage gas_ebu_examples
  - \subpage gas_libby_williams_examples
  - \subpage pulverized_coal1
  - \subpage pulverized_coal2

*/
// __________________________________________________________________________________
/*!

  \page base_examples Basic examples

  \section base_loc_var_bc Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_boundary_conditions-advanced.f90 loc_var_dec

  \section base_init_ex Initialization and finalization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_boundary_conditions-base.f90 init

  Ad the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_boundary_conditions-base.f90 finalize

  In theory Fortran 95 deallocates locally-allocated arrays automatically,
  but deallocating arrays in a symmetric manner to their allocation is good
  practice, and it avoids using a different logic for C and Fortran.

  \section base_body_bc Body

  In the body, we may define several boundary conditions. Here are a few examples.

  \section base_example_1_bc Inlet example with hydraulic diameter

  Assign an inlet to boundary faces of group '2' and x < 0.01.

  \warning the <, <=, >, and >= operators may only be used with variables x, y, and z.
  This syntax is not a full equation interpreter, so formulas involving x, y, or z
  are not allowed.

  Set a a Dirichlet value on the three components of \f$ \vect{u} \f$
  on the faces with the selection criterion '2 and x < 0.01' and set a Dirichlet
  to all the scalars \f$ \varia \f$.

  Turbulence example computed using equations valid for a pipe.

  We will be careful to specify a hydraulic diameter adapted
  to the current inlet.

  We will also be careful if necessary to use a more precise
  formula for the dynamic viscosity use in the calculation of
  the Reynolds number (especially if it is variable, it may be
  useful to take the law from \ref usphyv. Here, we use by default
  the 'viscl0" value.

  Regarding the density, we have access to its value at boundary
  faces (romb) so this value is the one used here (specifically,
  it is consistent with the processing in \ref usphyv, in case of
  variable density).

  \snippet cs_user_boundary_conditions-base.f90 example_1

  \section base_example_2_bc Inlet example with turbulence intensity

  Assign an inlet to boundary faces of group '3'.

  Set a a Dirichlet value on the three components of \f$ \vect{u} \f$
  on the faces with the selection criterion '3' and set a Dirichlet
  to all the scalars \f$ \varia \f$.

  Turbulence example computed using turbulence intensity data.

  We will be careful to specify a hydraulic diameter adapted
  to the current inlet.

  Calculation of \f$ k \f$ and \f$ \varepsilon \f$
  at the inlet (xkent and xeent) using
  the turbulence intensity and standard laws for a circular pipe
  (their initialization is not needed here but is good practice)

  \snippet cs_user_boundary_conditions-base.f90 example_2

  \section base_example_3_bc Assign an outlet to boundary faces of group 'outlet'

  Outlet:
  - zero flux for velocity and temperature, prescribed pressure
  - Note that the pressure will be set to P0 at the first
  - free outlet face (isolib)

  \snippet cs_user_boundary_conditions-base.f90 example_3

  \section base_example_4_bc Wall example

  Assign a wall to boundary faces of group '5'.

  Wall:
  - zero flow (zero flux for pressure)
  - friction for velocities (+ turbulent variables)
  - zero flux for scalars

  \snippet cs_user_boundary_conditions-base.f90 example_4

  \section base_example_5_bc Rough wall example

  Assign a rough wall to boundary faces of group '7'.

  Wall:
  - zero flow (zero flux for pressure)
  - rough friction for velocities (+ turbulent variables)
  - zero flux for scalars

  \snippet cs_user_boundary_conditions-base.f90 example_5

  \section base_example_6_bc Symmetry example

  Assign a symmetry condition to boundary faces of group '4'

  \snippet cs_user_boundary_conditions-base.f90 example_6

*/
// __________________________________________________________________________________
/*!

  \page channel_inlet_mapped Infinite channel inlet

  A method of defining an inlet condition which converges towards an infinite
  channel profile is based simply on feedback from values computed downstream
  from the inlet.

  Here, we define an inlet boundary condition for a very long channel or duct
  with a section matching the boundary faces of group 'INLET'.

  \section channel_inlet_mapped_loc_var_bc Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_boundary_conditions-mapped_inlet.f90 loc_var_dec

  Note that the \c inlet_l pointer has the \c save attribute, as the
  matching structure usually needs to be built only once, then reused at
  each time step.

  \section channel_inlet_mapped_init_bc Initialization and finalization

  Initialization and finalization is similar to that of the base examples,
  with the addition of the mapping object, described in a specific section
  hereafter

  \section channel_inlet_mapped_example_1_base Base definition

  \warning We assume other boundary conditions are defined before this one
  (ideally, using the GUI).

  The base definition given here ensures initialization of a (flat) inlet profile
  with the required mean velocity.

  \note We may skip the addition of the following block in the user
  subroutine if we define an equivalent inlet condition using the GUI,
  which will ensure the appropriate initialization before entering the
  user subroutine.

  \snippet cs_user_boundary_conditions-mapped_inlet.f90 example_1_base

  \section channel_inlet_mapped_example_1_map_init Mapping definition

  To define the appropriate mapping object, the
  \ref cs_c_bindings::boundary_conditions_map "boundary_conditions_map"
  function is used. In this example, coordinates of the selected inlet face
  are shifted by a constant distance (5.95 meters along the \em x axis),
  and mapped to the corresponding mesh cells. Here, all cells are selected
  as point location candidates, but optionally, a restricted list of cells
  may be provided, which may accelerate location (or ensure it is done on a
  selected subset). In most cases, as in this example, a constant coordinate
  shift is used for all inlet points, but it is possible to define a
  specific shift for each face (defining a \c stride argument of 1 instead of 0
  and \c coord_shift(:,ifac) instead of \c coord_shift(:,1)).

  \snippet cs_user_boundary_conditions-mapped_inlet.f90 example_1_map_init

  In general, when defining a pseudo-periodic boundary condition, it is assumed
  that the mesh is not moving, so the mapping may be defined once and for
  all. Hence the test on \ref optcal::ntcabs "ntcabs" and
  \ref optcal::ntpabs "ntpabs" and the \c save attribute
  for the \c inlet_1 pointer.

  \section channel_inlet_mapped_example_1_map_apply Applying the map

  At all time steps after the first (possibly even the first if the flow
  at the mapping locations is initialized to nonzero values), the values at
  points mapped to inlet face centers are applied to the \c rcodcl(ifac,1)
  values of the corresponding faces, using the
  \ref cs_c_bindings::boundary_conditions_mapped_set
  "boundary_conditions_mapped_set" subroutine.
  Optionally, a normalization by be applied,
  ensuring the mean values initially defined are preserved. Normalization
  is recommended for the velocity, and possibly scalars, but not for
  turbulent quantities, which should adapt to the velocity.

  \snippet cs_user_boundary_conditions-mapped_inlet.f90 example_1_map_apply

  \section channel_inlet_mapped_example_1_map_clean Finalization

  At the end of the computation, it is good practice to free the mapping
  structure:

  \snippet cs_user_boundary_conditions-mapped_inlet.f90 example_1_map_free

*/
// __________________________________________________________________________________
/*!

  \page channel_inlet_auto Infinite channel inlet with near-inlet feedback

  A second method of defining an inlet condition which converges towards
  an infinite channel profile is based simply on feedback from values
  computed at inlet cell centers (combining the inlet boundary conditions
  and the effect of wall friction on the inlet walls). It assumes the mesh
  is very regular and orthogonal, at least on the first two inlet cell layers
  (using a gradient correction in the case of a less regular mesh might work,
  but has never been tested.

  \section channel_inlet_auto_loc_var_bc Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_boundary_conditions-auto_inlet_profile.f90 loc_var_dec

  \section channel_inlet_auto_init_bc Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section channel_inlet_auto_example_1_bc Body

  Here, we define an inlet boundary condition for a very long channel or duct
  with a section matching the boundary faces of group 'INLET'.

  We fix a mean inlet velocity, and use a feedback loop assuming a fixed-point type
  behavior will allow us to reach a state matching that of a very long inlet channel.

  \warning We assume other boundary conditions are defined before this one
  (ideally, using the GUI).

  \warning We also assume that the mesh is orthogonal at the inlet, and we are using a RANS
  (not LES) computation.
  to the current inlet.

  For EBRSM of V2f models, to avoid laminarization at the inlet, the initial velocity
  (at the first time step) is divided by 10 on inlet faces adjacent to the boundary, so
  as to ensure a velocity gradient and turbulent production. Otherwise, the initialization
  may fail.

  \snippet cs_user_boundary_conditions-auto_inlet_profile.f90 example_1

*/
// __________________________________________________________________________________
/*!

  \page advanced_examples Advanced examples

  \section advanced_loc_var_be Local variables to be added

  \snippet cs_user_boundary_conditions-advanced.f90 loc_var_dec

  \section base_init_be Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section advanced_ex_1 Example 1

  Example of specific boundary conditions fully defined by the user,
  on the basis of wall conditions selection
  (mass flow computation, specific logging, ...)

  We prescribe for group '1234' a wall, with in addition:
  - a Dirichlet condition on velocity (sliding wall with no-slip condition)
  - a Dirichlet condition on the first scalar.

  \snippet cs_user_boundary_conditions-advanced.f90 example_1

  \section advanced_ex_2 Example 2

  Example of specific boundary conditions fully defined by the user,
  with no definition of a specific type.

  We prescribe at group '5678' a homogeneous Neumann condition for
  all variables.

  \snippet cs_user_boundary_conditions-advanced.f90 example_2

  \section advanced_ex_3 Example 3

  Example of specific boundary conditions fully defined by the user,
  with the definition of a specific type, for example for future
  selection (mass flow computation, specific logging, ...)

  We prescribe for group '6789' a homogeneous Neumann condition for
  all variables, except for the first
  scalar, for which we select a homogeneous Dirichlet.

  \snippet cs_user_boundary_conditions-advanced.f90 example_3

  \section advanced_ex_4 Example 4

  Example of wall boundary condition with automatic continuous switch
  between rough and smooth.

  \snippet cs_user_boundary_conditions-advanced.f90 example_4

*/
// __________________________________________________________________________________
/*!

  \page atmospheric_examples Atmospheric examples

  \section advanced_loc_var_atm Local variables to be added

  \snippet cs_user_boundary_conditions-atmospheric.f90 loc_var_dec

  \section base_init_atm Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section atmospheric_ex_1 Example 1

  For boundary faces of color 11,
  assign an inlet boundary condition prescribed from the
  meteo profile with automatic choice between inlet/ outlet according to
  the meteo profile.

  \snippet cs_user_boundary_conditions-atmospheric.f90 example_1

  \section atmospheric_ex_2 Example 2

  For boundary faces of color 21,
  assign an inlet boundary condition prescribed from the
  meteo profile.

  \snippet cs_user_boundary_conditions-atmospheric.f90 example_2

  \section atmospheric_ex_3 Example 3

  For boundary faces of color 31,
  assign an inlet boundary condition prescribed from the
  meteo profile except for dynamical variables which are prescribed with
  a rough log law.

  \snippet cs_user_boundary_conditions-atmospheric.f90 example_3

  \section atmospheric_ex_4 Example 4

   Prescribe at boundary faces of color '12' an outlet.

  \snippet cs_user_boundary_conditions-atmospheric.f90 example_4

  \section atmospheric_ex_5 Example 5

  Prescribe at boundary faces of color 4 a symmetry.

  \snippet cs_user_boundary_conditions-atmospheric.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page compressible_examples Compressible examples

  \section advanced_loc_var_ce Local variables to be added

  \snippet cs_user_boundary_conditions-compressible.f90 loc_var_dec

  \section base_init_ce Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section compressible_ex_1 Example 1

   Example of input / output for which everything is known.

   Without presuming subsonic or supersonic character,
   the user wishes to impose all the flow characteristics.
   A supersonic inlet is a special case.

   If the speed is outgoing, an homogenous Neumann is imposed
   on turbulence and user scalars.

  \snippet cs_user_boundary_conditions-compressible.f90 example_1

  \section compressible_ex_2 Example 2

  Example supersonic output

   All features are output.
   Internal values are used to calculate the flow edge, they should not be imposed.


   For turbulence and scalar, if RCODCL values are provided,
   they will be used in Dirichlet if the mass flow is incoming,
   otherwise a null flow is imposed (flow outgoing mass or RCODCL informed here).
   Note that for turbulence RCODCL must be defined for all turbulent variables.
   Otherwise a null  flow is applied).

TODO : verifier la traduction.
  Exemple de sortie supersonique

     toutes les caracteristiques sortent
     on ne doit rien imposer (ce sont les valeurs internes qui sont
       utilisees pour le calcul des flux de bord)

     pour la turbulence et les scalaires, si on fournit ici des
       valeurs de RCODCL, on les impose en Dirichlet si le flux
       de masse est entrant ; sinon, on impose un flux nul (flux de
       masse sortant ou RCODCL renseigne ici).
       Noter que pour la turbulence, il faut renseigner RCODCL pour
       toutes les variables turbulentes (sinon, on applique un flux
      nul).

  \snippet cs_user_boundary_conditions-compressible.f90 example_2

  \section compressible_ex_3 Example 3

  Subsonic input example (flow, flow calorimetry)

  Two of three characteristics are incoming: two informations must be provided,
  third is deduced by a scenario of 2-contact and 3-relaxation in the field.
  Here we choose to give (rho * (U.n) * rho (U.n) * H)
   - where    H = 1/2 * G G + P / E + rho
   - and n is the unit normal incoming

  \warning DENSITIES debit (per unit area) are provided.

TODO : Verifier la traduction

 Exemple d'entree subsonique (debit, debit enthalpique)

     2 caracteristiques sur 3 entrent : il faut donner 2 informations
       la troisieme est deduite par un scenario de 2-contact et
       3-detente dans le domaine
     ici on choisit de donner (rho*(U.n), rho*(U.n)*H)
       avec H = 1/2 U*U + P/rho + e
            n la normale unitaire entrante

     ATTENTION, on donne des DENSITES de debit (par unite de surface)


  \snippet cs_user_boundary_conditions-compressible.f90 example_3

  \section compressible_ex_4 Example 4

  Subsonic input example with density and velocity.

  Two of three characteristics are incoming: two informations must be provided,
  third is deduced by a scenario of 2-contact and 3-relaxation in the field.
  Here we choose to give (rho, U).

TODO : Verifier la traduction
 Exemple d'entree subsonique (masse volumique, vitesse)

 2 caracteristiques sur 3 entrent : il faut donner 2 informations
 la troisieme est deduite par un scenario de 2-contact et
  3-detente dans le domaine
   ici on choisit de donner (rho, U)


  \snippet cs_user_boundary_conditions-compressible.f90 example_4

  \section compressible_ex_5 Example 5

   Subsonic outlet example

   1 characteristic out of 3 exits: 1 information must be given
   the 2 others are deduced by a 2-contact and 3-relaxation in the domain.
   Here we choose to definer P.

   Turbulence and user scalars take a zero flux.

  \snippet cs_user_boundary_conditions-compressible.f90 example_5

  \section compressible_ex_6 Example 6

  Wall example

  \snippet cs_user_boundary_conditions-compressible.f90 example_6

  \section compressible_ex_7 Example 7

  Symmetry example

  \snippet cs_user_boundary_conditions-compressible.f90 example_7

*/
// __________________________________________________________________________________
/*!

  \page electric_arcs_examples Electric arcs examples

  \section elec_loc_var Local variables to be added

  \snippet cs_user_boundary_conditions-electric_arcs.f90 loc_var_dec

  \section elec_init Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section electric_arcs_ex_1 Example 1

  For boundary faces of color 1 assign an inlet
  and assign a cathode for "electric" variables.

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_1

  \section electric_arcs_ex_2 Example 2

  For boundary faces of color 5 assign an free outlet
  and example of electrode for Joule Effect by direct conduction.

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_2

  \section electric_arcs_ex_3 Example 3

   For boundary faces of color 2 assign a free outlet
   and example of anode for electric arc.

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_3

  \section electric_arcs_ex_4 Example 4

  For boundary faces of color 3 assign a wall
  and example of potential vector Dirichlet condition

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_4

  \section electric_arcs_ex_5 Example 5

   For boundary faces of color 51 assign a wall
   and restriking model for electric arc (anode boundaray condition).

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_5

  \section electric_arcs_ex_6 Example 6

  For boundary faces of color 4 assign a symmetry.

  \snippet cs_user_boundary_conditions-electric_arcs.f90 example_6

*/
// __________________________________________________________________________________
/*!

  \page electric_arcs_ieljou_3_or_4_examples Electric arcs Joule examples

  \section elec2_loc_var Local variables to be added

  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c loc_var_dec

  \section elec2_init Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c init

  \section elec2_step1 Computation of intensity (A/m2) for each electrode

  Pre initialisation

  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c pre_init

  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_1

  \section elec2_step2 Definition of Voltage on each termin of transformers

  Computation of Intensity on each termin of transformers:
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_2_1

  RVoltage on each termin:
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_2_2

  Total intensity for a transformer (zero valued WHEN Offset established):
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_2_3

  Take in account of Offset:
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_2_4

  Take in account of Boundary Conditions
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_2_5

  \section elec2_step3 Finalization step

  Test, if not any reference transformer a piece of wall may be at ground:
  \snippet cs_user_boundary_conditions-electric_arcs_ieljou_3_or_4.c step_3

*/
// __________________________________________________________________________________
/*!

  \page fuel_examples Fuel examples

  \section advanced_loc_var_fe Local variables to be added

  \snippet cs_user_boundary_conditions-fuel.f90 loc_var_dec

  \section base_init_fe Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section fuel_ex_1 Example 1

  The 12 color is a pure air inlet

  \snippet cs_user_boundary_conditions-fuel.f90 example_1

  \section fuel_ex_2 Example 2

  Inlet of both primary Air and Fuel

  \snippet cs_user_boundary_conditions-fuel.f90 example_2

  \section fuel_ex_3 Example 3

  Color 15 is a wall

  \snippet cs_user_boundary_conditions-fuel.f90 example_3

  \section fuel_ex_4 Example 4

  Color 19 is an outlet

  \snippet cs_user_boundary_conditions-fuel.f90 example_4

  \section fuel_ex_5 Example 5

  14 and 4 are symmetry

  \snippet cs_user_boundary_conditions-fuel.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page gas_3ptchem_examples Gas 3 PTCHEM examples

  \section advanced_loc_var_ptchem Local variables to be added

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 loc_var_dec

  \section base_init_ptchem Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section gas_3ptchem_ex_1 Example 1

  Definition of a fuel flow inlet for each face of colour 11

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 example_1

  \section gas_3ptchem_ex_2 Example 2

  Definition of an air flow inlet for each face of colour 21

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 example_2

  \section gas_3ptchem_ex_3 Example 3

  Definition of a wall for each face of colour 51 up to 59

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 example_3

  \section gas_3ptchem_ex_4 Example 4

  Definition of an exit for each face of colour 91

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 example_4

  \section gas_3ptchem_ex_5 Example 5

  Definition of symmetric boundary conditions for each face of colour 41 and 4.

  \snippet cs_user_boundary_conditions-gas_3ptchem.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page gas_ebu_examples Gas EBU examples

  \section advanced_loc_var_ebu Local variables to be added

  \snippet cs_user_boundary_conditions-gas_ebu.f90 loc_var_dec

  \section base_init_ebu Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section gas_ebu_ex_1 Example 1

  Definition of a burned gas inlet (pilot flame) for each face of colour 11

  \snippet cs_user_boundary_conditions-gas_ebu.f90 example_1

  \section gas_ebu_ex_2 Example 2

  Definition of an unburned gas inlet for each face of colour 12

  \snippet cs_user_boundary_conditions-gas_ebu.f90 example_2

  \section gas_ebu_ex_3 Example 3

  Definition of a wall for each face of colour 51 and 5

  \snippet cs_user_boundary_conditions-gas_ebu.f90 example_3

  \section gas_ebu_ex_4 Example 4

  Definition of an exit for each face of colour 91 and 9

  \snippet cs_user_boundary_conditions-gas_ebu.f90 example_4

  \section gas_ebu_ex_5 Example 5

  Definition of symmetric boundary conditions for each face of colour 41 and 4.

  \snippet cs_user_boundary_conditions-gas_ebu.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page gas_libby_williams_examples Gas Libby-Williams examples

  \section advanced_loc_var_lw Local variables to be added

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 loc_var_dec

  \section base_init_lw Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section gas_libby_williams_ex_1 Example 1

  Definition of a burned gas inlet (pilot flame) for each face of colour 11

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 example_1

  \section gas_libby_williams_ex_2 Example 2

  Definition of an unburned gas inlet for each face of colour 12

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 example_2

  \section gas_libby_williams_ex_3 Example 3

  Definition of a wall for each face of colour 51 and 5

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 example_3

  \section gas_libby_williams_ex_4 Example 4

  Definition of an exit for each face of colour 91 and 9

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 example_4

  \section gas_libby_williams_ex_5 Example 5

  Definition of symmetric boundary conditions for each face of colour 41 and 4.

  \snippet cs_user_boundary_conditions-gas_libby_williams.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page pulverized_coal1 Pulverized coal

  \section advanced_loc_var_pc Local variables to be added

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 loc_var_dec

  \section base_init_pc Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  \section pulverized_coal1_ex_1 Example 1

  BOUNDARY FACE corresponding to AIR INLET, e.g. : secondary or tertiary air.

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 example_1

  \section pulverized_coal1_ex_2 Example 2

  BOUNDARY FACE for pulverised COAL & primary air INLET.

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 example_2

  \section pulverized_coal1_ex_3 Example 3

  The color 15 becomes a WALL.

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 example_3

  \section pulverized_coal1_ex_4 Example 4

  The color 19 becomes an OUTLET.

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 example_4

  \section pulverized_coal1_ex_5 Example 5

  The color 14 becomes a symmetry plane.

  \snippet cs_user_boundary_conditions-pulverized_coal.f90 example_5

*/
// __________________________________________________________________________________
/*!

  \page pulverized_coal2 Pulverized coal lagrangian

  \section advanced_loc_var_pcl Local variables to be added

  \snippet cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 loc_var_dec

  \section base_init_pcl Initialization and finalization

  Initialization and finalization is similar to that of the base examples

  Assign boundary conditions to boundary faces here
  For each subset:
   - use selection criteria to filter boundary faces of a given subset
   - loop on faces from a subset
   - set the boundary condition for each face
  The color 15 become a WALL.
  \section pulverized_coal2_ex_1 Example 1

  BOUNDARY FACE corresponding to AIR INLET, e.g.: secondary or tertiary air.

  \snippet cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 example_1

  \section pulverized_coal2_ex_2 Example 2

  The color 15 becomes a WALL.

  \snippet cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 example_2

  \section pulverized_coal2_ex_3 Example 3

  The color 19 becomes an OUTLET.

  \snippet cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 example_3

  \section pulverized_coal2_ex_4 Example 4

  The colors 14 and 4 become a symmetry plane.

  \snippet cs_user_boundary_conditions-pulverized_coal_lagrangian.f90 example_4

*/
