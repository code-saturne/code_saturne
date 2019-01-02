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

  \page cs_user_source_terms-scalar_in_a_channel Examples of data settings for source terms with scalar in a channel (cs_user_source_terms-scalar_in_a_channel.f90)

\brief Additional right-hand side source terms for scalar equations (user
 scalars and specific physics scalars) with \ref ustssc subroutine.

\section usa_scal_channel Usage

The routine is called for each scalar, user or specific physisc. It is
therefore necessary to test the value of the scalar number iscal to separate
the treatments of the different scalars <tt>(if (iscal.eq.p) then ....)</tt>.

The additional source term is decomposed into an explicit part \f$(crvexp)\f$ and
an implicit part \f$(crvimp)\f$ that must be provided here.
The resulting equation solved by the code for a scalar \f$f\f$ is:

 \f[\rho \norm{\vol{\celli}}  \frac{df}{dt} + .... = crvimp \cdot f + crvexp\f]


\note \f$ crvexp \f$ and \f$ crvimp\f$ are defined after the Finite Volume integration
over the cells, so they include the \f$\norm{\vol{\celli}}\f$ term. More precisely:
  - \f$ crvexp\f$ is expressed in \f$ kg\cdot [scal]\cdot s^{-1}\f$, where \f$ [scal]\f$ is the unit of the scalar
  - \f$crvimp\f$ is expressed in \f$kg \cdot s^{-1} \f$


The \f$ crvexp \f$  and \f$ crvimp \f$ arrays are already initialized to 0 before entering the
the routine. It is not needed to do it in the routine (waste of CPU time).

For stability reasons, \c Code_Saturne will not add \c -crvimp directly to the
diagonal of the matrix, but \c Max(-crvimp,0). This way, the \c crvimp term is
treated implicitely only if it strengthens the diagonal of the matrix.
However, when using the second-order in time scheme, this limitation cannot
be done anymore and \c -crvimp is added directly. The user should therefore test
the negativity of \c crvimp by himself.

When using the second-order in time scheme, one should supply:
  - \f$ crvexp \f$ at time n
  - \f$ crvimp \f$ at time n+1/2


The selection of cells where to apply the source terms is based on a getcel
command. For more info on the syntax of the getcel command, refer to the
user manual or to the comments on the similar command \ref getfbr in the routine
\ref cs_user_boundary_conditions.

\warning If scalar is the temperature, the resulting equation
         solved by the code is:

 \f[\rho C_p \norm{\vol{\celli}}  \frac{dT}{dt} + .... = crvimp \cdot T + crvexp\f]


\note \f$ crvexp \f$ and \f$ crvimp \f$ are defined after the Finite Volume integration
over the cells, so they include the \f$\norm{\vol{\celli}}\f$ term. More precisely:
  - \f$ crvexp \f$ is expressed in \f$W\f$
  - \f$ crvimp \f$ is expressed in \f$ W \cdot K^{-1}\f$



\section steep_scalar_channel Steep source terms

In case of a complex, non-linear source term, say \f$F(f)\f$, for scalar \f$f\f$, the
easiest method is to implement the source term explicitely.

  \f[\frac{df}{dt} = .... + F(f^{(n)})\f]
  where \f$f^{(n)}\f$ is the value of \f$f\f$ at time \f$t^n\f$, the beginning of the time step.

This yields :
 - <tt>crvexp = volume*F(f(n))</tt>
 - <tt>crvimp = 0</tt>

However, if the source term is potentially steep, this fully explicit
method will probably generate instabilities. It is therefore wiser to
partially implicit the term by writing:

  \f[\frac{df}{dt} = .... + \frac{dF}{df} f^{(n+1)} - \frac{dF}{df} f^{(n)} + F(f^{(n)})\f]

This yields:
 - <tt>crvexp = volume*(F(f(n)) - dF/df*f(n)) </tt>
 - <tt>crvimp = volume*dF/df</tt>

\section loc_var_scal Local variables

\snippet cs_user_source_terms-scalar_in_a_channel.f90 loc_var_scal


\section body_scal Body
Initialization
\snippet cs_user_source_terms-scalar_in_a_channel.f90 init_scal

Map field arrays
\snippet cs_user_source_terms-scalar_in_a_channel.f90 map_field_scal

Parallelism
 \snippet cs_user_source_terms-scalar_in_a_channel.f90 parallelism_scal

\subsection format_scal Formats
 \snippet cs_user_source_terms-scalar_in_a_channel.f90 format_scal

\subsection end_scal End
 \snippet cs_user_source_terms-scalar_in_a_channel.f90 end_scal


*/
