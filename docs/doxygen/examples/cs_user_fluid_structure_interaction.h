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

/*----------------------------------------------------------------------------*/

/*!
  \page cs_user_fluid_structure_interaction Fluid Structure Interaction (cs_user_fluid_structure_interaction.f90)

  \section cs_user_fluid_structure_interaction_h_intro Introduction

  This page provides code blocks of several examples that may be used
  to perform fluid structure interaction computations.

  \section cs_user_fluid_structure_interaction_h_usstr1  Internal structures and corresponding initial conditions

  \subsection cs_user_fluid_structure_interaction_h_usstr1_init Initialization

  The \c lstelt array is allocated.

  \snippet cs_user_fluid_structure_interaction.f90 usstr1_init

  \subsection cs_user_fluid_structure_interaction_h_usstr1_example Example

  The following code block provides an example of definition of two internal
  structures.

  Here one fills array \c idfstr(\ref nfabor). For each boundary face \c ifac,
  \c idfstr(ifac) is the number of the structure the face belongs to
  (if \c idfstr(ifac) = 0, the face ifac doesn't belong to any structure.
  When using internal coupling, structure number necessarily needs to be positive
  (as shown in following examples).

  The number of "internal" structures is automatically defined with the
  maximum value of \ref idfstr table, meaning that internal structure numbers must
  be defined sequentially with positive values, beginning with integer value 1.

  In the following example, boundary faces with color 4 belong to internal structure 1.
  Boundary faces with color 2 belong to internal structure 2.
  The total number of internal structures equals 2.
  The boundary faces are identified using the \ref getfbr subroutine.

  \snippet cs_user_fluid_structure_interaction.f90 usstr1_example_a

  For each internal structure one can here define :
  - an initial velocity \c vstr0
  - an initial displacement \c xstr0 (i.e. \c xstr0 is the value of the
    displacement \ref xstr compared to the initial mesh at time t=0)
  - a displacement compared to equilibrium \c \ref xstreq.  \c \ref xstreq is the
    initial displacement of the internal structure compared to its position
    at equilibrium; at each time step t and for a displacement \ref xstr(t),
    the associated internal structure will be subjected to a force due to the spring:
    \f[
       -k (xstr(t)+xstreq).
    \f]

  Note that \c xstr0, \c \ref xstreq and \c vstr0 arrays are initialized at the beginning
  of the calculations to the value of 0.

  When starting a calculation using ALE, or re-starting a calculation with
  ALE basing on a first calculation without ALE, an initial iteration 0 is
  automatically calculated in order to take initial arrays \c xstr0, \c vstr0 and
  \c \ref xstreq into account. In another case, set the following option
  \code{.f90}
  italin=1
  \endcode
  in subroutine \ref usipsu, so that the code can deal with arrays
  \c xstr0, \c vstr0 or \c \ref xstreq.

  In the following example :
  - internal structure 1 has got an initial displacement \c xstr0 of 2 meters in
    the y direction and a displacement compared to equilibrium \c \ref xstreq of
    1 meter in the y direction.
  - Initial velocity \c vstr0 in the z direction of structure 2 equals -0.5 m/s.

  \snippet cs_user_fluid_structure_interaction.f90 usstr1_example_b

  Here one can modify the values of the prediction coefficients for
  displacements anf fluid forces in internal FSI coupled algorithm.

  The use of these coefficients is reminded here :
  - predicted displacement
    \f[
    \vect{X_{n+1}} = \vect{X_n} + aexxst \Delta t \vect{X^{\prime}_n} +
    bexxst \Delta t (\vect{X^{\prime}_n}-\vect{X^{\prime}_{n-1}})
    \f]
    with \f$ \vect{X_n} \f$ standing for the displacement at iteration \f$ n \f$,
    \f$ \vect{X^{\prime}_n} \f$ and \f$ \vect{X^{\prime}_{n-1}} \f$ representing the internal
    structure velocity respectively at iteration \f$ n \f$ and \f$ n-1 \f$.
  - fluid force sent to structure internal solver
    \f[
    \vect{F_{n+1}} = cfopre \vect{F_n} + (1-cfopre) \vect{F_{n-1}}
    \f]
    \f$ \vect{F_n} \f$ and \f$ \vect{F_{n-1}} \f$ standing for the fluid force acting on the
    structure respectively at iteration \f$ n \f$ and \f$ n-1 \f$.

  \snippet cs_user_fluid_structure_interaction.f90 usstr1_example_c

  Activation of structural history output (i.e. displacement, structural
  velocity, structural acceleration and fluid force)
  (\c ihistr=0, disabled ; \c ihistr=1, enabled)
  The value of structural history output step is the same as the one for
  standard variables (\ref nthist).

  \snippet cs_user_fluid_structure_interaction.f90 usstr1_example_d

  \section cs_user_fluid_structure_interaction_h_usstr2  Structural parameters in case of Fluid Structure internal coupling

  Note that the subroutine \ref usstr2 is called at each time step of the calculation.

  For each internal structure one defines here :
  - its mass \f$ \tens{M} \f$ (\c xmstru)
  - its friction coefficient \f$ \tens{C} \f$ (\c xcstru)
  - its stiffness \f$ \tens{K} \f$ (\c xkstru)

  \c forstr array gives fluid stresses acting on each internal structure.
  Moreover it's possible to take external forces (gravity for example) into
  account, too.

  \c \ref xstr array indicates the displacement of the structure compared to its
  position in the initial mesh.

  \c xstr0 array gives the displacement of the structures in initial mesh
  compared to structural equilibrium.

  \c \ref vstr array stands for structural velocity.

  \c \ref xstr, \c xstr0, and \c \ref vstr arrays are data tables that can be used to define
  arrays mass, friction and stiffness. THOSE ARE NOT TO BE MODIFIED.

  The 3D structural equation that is solved is the following one :
  \f[
    \tens{M} \vect{X^{\prime\prime}} + \tens{C} \vect{X^{\prime}} + \tens{K} (\vect{X}+\vect{X_0}) = \vect{F}
  \f]
  where \f$\vect{X}\f$ stands for the structural displacement compared to initial mesh
      position (\c \ref xstr) and \f$\vect{X_0}\f$ represents the displacement of the structure in initial mesh compared to equilibrium.

    Note that \f$\tens{M}\f$, \f$\tens{C}\f$ and \f$\tens{K}\f$ are 3x3 matrices.

    This equation is solved using a Newmark HHT algorithm.
    Note that the time step used to solve this equation (\c dtstr) can be
    different from the one of fluid calculations. user is free to define
    \c dtstr array. At the beginning of the calculation \c dtstr is
    initialized to the value of \c \ref dt (fluid time step).

  \subsection cs_user_fluid_structure_interaction_h_usstr2_init Initialization

  The matrices \c xmstru, \c xcstru and \c xkstru are set to zero.

  \snippet cs_user_fluid_structure_interaction.f90 usstr2_init

  \subsection cs_user_fluid_structure_interaction_h_usstr2_example_1 Example 1

  Two examples of definition of the mass, the friction coefficient and the stiffness
  of an internal structure are provided hereafter.

  \snippet cs_user_fluid_structure_interaction.f90 usstr2_example_1

  \subsection cs_user_fluid_structure_interaction_h_usstr2_example_2 Example 2

  \snippet cs_user_fluid_structure_interaction.f90 usstr2_example_2

  \section cs_user_fluid_structure_interaction_h_usaste Fluid Structure external coupling with Code_Aster

  \subsection cs_user_fluid_structure_interaction_h_usaste1_init Initialization

  The \c lstelt array is allocated.

  \snippet cs_user_fluid_structure_interaction.f90 usaste_init

  \subsection cs_user_fluid_structure_interaction_h_usaste_example Example

  In the following example, two sets of boundary faces that will be coupled with
  Code_Aster are defined as well as the fluid force components which are given
  to structural calculations.

  Here one fills array \c idfstr(\ref nfabor)
  For each boundary face \c ifac, \c idfstr(ifac) is the number of the structure
  the face belongs to (if \c idfstr(ifac) = 0, the face \c ifac doesn't
  belong to any structure.)
  When using external coupling with Code_Aster, structure number necessarily
  needs to be negative (as shown in following examples).

  The number of "external" structures with Code_Aster is automatically
  defined with the maximum absolute value of \c idfstr table, meaning that
  external structure numbers must be defined sequentially with negative
  values beginning with integer value -1.

  In following example, boundary faces with color 2 and which
  abscissa \f$ x < 2 \f$ belong to external structure -1.
  Boundary faces with color 2 and which abscissa \f$ x > 2 \f$ belong to external
  structure -2. The total number of external structures coupled with
  Code_Aster equals 2.

  \snippet cs_user_fluid_structure_interaction.f90 usaste_example

*/
