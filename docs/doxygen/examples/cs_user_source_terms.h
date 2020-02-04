/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

  \page user_source_terms Examples of data settings for source terms (cs_user_source_terms.c)

   See also the \ref cs_user_source_terms reference documentation.


  \brief Additional right-hand side source terms


   - \subpage subroutine_ustsnv
   - \subpage function_cs_user_source_terms_scalar
   - \subpage function_cs_user_turbulence_source_terms
   - \subpage cs_user_source_terms-scalar_in_a_channel


*/
//________________________________________________________________________________________________
/*!

   \page subroutine_ustsnv For velocity components equation (Navier-Stokes): ustsnv subroutine

 \brief Additional right-hand side source terms for velocity components equation
 (Navier-Stokes)

  \section use_src1  Usage

  The additional source term is decomposed into an explicit part \f$(\vect{crvexp})\f$ and
  an implicit part \f$(\tens{crvimp})\f$ that must be provided here.
  The resulting equation solved by the code for a velocity is:
  \f[
  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
   = \tens{crvimp}\cdot  \vect{u} + \vect{crvexp}
  \f]

  Note that \f$\vect{crvexp}\f$ and \f$\tens{crvimp}\f$ are defined after the Finite Volume integration
  over the cells, so they include the \f$\norm{\vol{\celli}}\f$ term. More precisely:
   - \f$\vect{crvexp}\f$ is expressed in \f$ kg\cdot ms^{-2} \f$
   - \f$\tens{crvimp}\f$ is expressed in \f$ kg\cdot s^{-1} \f$

  The \f$\vect{crvexp}\f$ and \f$\tens{crvimp}\f$ arrays are already initialized to 0
  before entering the
  the routine. It is not needed to do it in the routine (waste of CPU time).

  For stability reasons, \c Code_Saturne will not add \c -crvimp directly to the
  diagonal of the matrix, but \c Max(-crvimp,0). This way, the \f$ \tens{crvimp}\f$ term  is
  treated implicitely only if it strengthens the diagonal of the matrix.
  However, when using the second-order in time scheme, this limitation cannot
  be done anymore and \c -crvimp is added directly. The user should therefore test
  the negativity of \c crvimp by himself.

  When using the second-order in time scheme, one should supply:
   - \f$\vect{crvexp}\f$ at time n
   - \f$\tens{crvimp}\f$ at time n+1/2

  The selection of cells where to apply the source terms is based on a
  \ref getcel command. For more info on the syntax of the \ref getcel command,
  refer to the user manual or to the comments on the similar command
  \ref getfbr in the routine \ref cs_user_boundary_conditions.

   \section loc_var1 Local variables

   \snippet cs_user_source_terms.f90 loc_var_dec_1

   \section init_and_finit_1 Initialization and finalization

The following initialization block needs to be added for the following examples:
   \snippet cs_user_source_terms.f90 allocate_1

At the end of the subroutine, it is recommended to deallocate the work array:
   \snippet cs_user_source_terms.f90 deallocate_1

In theory Fortran 95 deallocates locally-allocated arrays automatically, but deallocating arrays in a symetric manner to their alloacation is good pratice, and avoids using a different logic C and Fortran.

   \section example_source_terms_1  Example
Example of arbitrary source term for component \f$\vect{u}\f$:

   \f$ \vect{S} = \tens{A} \cdot \vect{u} + \vect{B} \f$ appearing in the equation under the form:

   \f$ \rho \dfrac{d\vect{u}}{dt} = \vect{S} \: (+ \text{standard Navier-Stokes terms})\f$


In the following example:
  \f[  \tens{A} = -\rho \cdot \tens{CKP} \f]
  \f[  \vect{B} =  \vect{XMMT}  \f]

with:
 - <tt> CKP = 1.0 </tt>(in \f$ s^{-1}\f$) (return term on velocity) \n
 - <tt> MMT = 100.0 </tt>(in \f$kg \cdot m^{-2} \cdot s^{-2}\f$) (momentum production by volume and time unit)

which yields:
 - <tt>  crvimp(1, 1, iel) = volume(iel)* A = - volume(iel)*(rho*CKP )</tt> \n
 - <tt>  crvexp(1, iel) = volume(iel)* B = volume(iel)*(XMMT) </tt>

 \section body1 Body

   \snippet cs_user_source_terms.f90 remaining_1

   \section example_source_terms_2 Example of a boussinesq momentum source term
   Example to add Boussinesq source to the z component of \f$\vect{u}\f$:

 \section body2 Body

   \snippet cs_user_source_terms.f90 boussinesq_st

*/
//_________________________________________________________________________________________________
/*!
   \page function_cs_user_source_terms_scalar Transported scalar source terms

Source terms for transported scalars may be defined using the \ref cs_user_source_terms
user-defined function.

 \subsection field_meta Field access and information

The following initialization block or portions thereof needs to be added for the following examples:
   \snippet cs_user_source_terms-base.c st_meta

 Indicator of variance scalars:
To determine whether a scalar is a variance, the following info can be accessed:
   \snippet  cs_user_source_terms-base.c field_is_variance

 - If <tt> var_f_id == -1 </tt>,
   the scalar is not a variance
 - If <tt> var_f_id >= 0 </tt>,
   the field is the variance of the scalar with field id \c var_f_id

 Density

  \snippet cs_user_source_terms-base.c density_2

  \section examplesource_2_1 Example 1

Example of arbitrary source term for the scalar f, named "scalar_2" in the calculation.

  \f$ S=A \cdot f+B \f$

appearing in the equation under the form

\f$ \rho \dfrac{df}{dt}=S \: \text{(+ regular other terms in the equation)} \f$
In the following example:

\f[A=-\frac{\rho}{\tau_f} \f]
\f[B=\rho \cdot prod_f \f]

with:
 - tauf = 10.0 (in \f$ s \f$) (dissipation time for \f$f\f$)
 - prodf = 100.0 (in \f$ [f]\cdot s^{-1} \f$) (production of \f$f\f$ by unit of time)

which yields:
 - <tt> st_imp[i] = volume[i]*A   = -volume[i]*rho/tauf </tt>
 - <tt> st_exp[i] = volume[i]*B =  volume[i]*rho*prod_f </tt>

 \subsection bodysource2 Body

 <b> Source term applied to second scalar</b>

  \snippet cs_user_source_terms-base.c src_term_applied

  \section examplesource2_2 Example 2

Example of arbitrary volumic heat term in the equation for enthalpy h.


 In the considered example, a uniform volumic source of heating is imposed
 in the cells with coordinate X in [0;1.2] and Y in [3.1;4].

 The global heating power if \c Pwatt (in \f$W\f$) and the total volume of the
 selected cells is \c volf (in \f$m^3\f$).

 This yields:
    - <tt> st_imp[i] = 0 </tt>
    - <tt> st_exp[i] = volume[i]* pwatt/volf </tt>


  \subsection end2 Body

   \warning It is assumed here that the thermal scalar is an enthalpy. If the scalar is a  temperature. PWatt does not need to  be divided by \f$ C_p \f$ because \f$C_p\f$ is put outside the diffusion term and multiplied in the temperature equation as follows:
 \f[ \rho C_p \norm{\vol{\celli}} \frac{dT}{dt} + ... = \norm{\vol{\celli}[i]} \frac{pwatt}{voltf} \f]

with  <tt> pwatt = 100.0 </tt>

   \subsection cs_user_st_3_apply Apply source term

   \snippet cs_user_source_terms-base.c ex_3_apply

*/
//________________________________________________________________________________________________
/*!

  \page function_cs_user_turbulence_source_terms Turbulence model source terms.

  Turbulence source terms may be modified using the \ref cs_user_source_terms
  user-defined function.

  \brief Additional right_hand side source terms for turbulence models

  \section loc_var3 Local variables

  \snippet cs_user_source_terms-turbulence.c st_meta

  <b> Remaining initialization</b>

  Get the density array in \c cpro_rom

  \snippet cs_user_source_terms-turbulence.c dens_array_3

  Get the array of the current turbulent variable and its name

  \snippet cs_user_source_terms-turbulence.c  current_turb_3

  \section example3_1 Example

  Example of arbitrary additional source term for turbulence models (Source term on the TKE "k" here).

  Source term for \f$\varia:\f$
   \f[ \rho \norm{\vol{\celli}} \frac{d(\varia)}{dt} = ... - \rho \norm{\vol{\celli}} \cdot ff - \rho \frac{ \varia}{\tau}\f]
  with \f$ ff \f$ = <tt>  3.0 </tt> an \f$ \tau \f$ = <tt> 4.0 </tt>


  \section body3 Body

  \note The turbulence variable names are:
   - 'k' and 'epsilon' for the k-epsilon models
   - 'rij' and 'epsilon' for the Rij-epsilon LRR and S SG
   - 'rij', 'epsilon' and 'alpha' for the EBRSM
   - 'k', 'epsilon', 'phi' and 'f_bar' for the phi-model
   - 'k', 'epsilon', 'phi' and 'alpha' for the Bl-v2-k model
   - 'k' and 'omega' for the k-omega turbulence model
   - 'nu_tilda' for the Spalart Allmaras model


  \subsection cal_exp_imp3 Calculation of the explicit and implicit source terms

  \snippet cs_user_source_terms-turbulence.c rem_code_3
*/
