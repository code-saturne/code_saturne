/*============================================================================
 * code_saturne documentation page
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

   - \subpage function_cs_user_source_terms_momentum
   - \subpage function_cs_user_source_terms_scalar
   - \subpage function_cs_user_turbulence_source_terms
   - \subpage cs_user_source_terms-scalar_in_a_channel
*/

//________________________________________________________________________________________________
/*!

 \page function_cs_user_source_terms_momentum For velocity components equation (Navier-Stokes)

 \brief Additional right-hand side source terms for velocity components equation
 (Navier-Stokes)

   \section loc_var1 Local variables and initialization

   \snippet cs_user_source_terms-momentum.c st_meta

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
 - <tt>  st_imp[i][0][0] = volume[i] * A = - volume[i]*(rho*CKP)</tt> \n
 - <tt>  st_exp[i][0]    = volume[i] * B =   volume[i]*(XMMT) </tt>

 \section body1 Body

 \snippet cs_user_source_terms-momentum.c st_momentum_e_1

 \section example_source_terms_2 Example of a boussinesq momentum source term

 Example to add Boussinesq source to the z component of \f$\vect{u}\f$:

 \section body2 Body

 \snippet cs_user_source_terms-momentum.c boussinesq_st
*/

//_________________________________________________________________________________________________
/*!
   \page function_cs_user_source_terms_scalar Transported scalar source terms

Source terms for transported scalars may be defined using the \ref cs_user_source_terms
user-defined function.

 \subsection field_meta Field access and information

The following initialization block or portions thereof needs to be added for the
following examples:
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
