/*============================================================================
 * This function is called each time step to define physical properties
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms.c
 *
 * \brief Additional right-hand side source terms for variable equations
 *   (momentum, user scalars and specific physics scalars, turbulence...).
 *
 *  Usage
 *  -----
 *  The routine is called for each variable. It is
 *  therefore necessary to test the value of the field id to separate
 *  the treatments of the different variables (if (f_id.eq.CS_F(p)->id) then ....).
 *
 *  The additional source term is decomposed into an explicit part (st_exp) and
 *  an implicit part (st_imp) that must be provided here.
 *  The resulting equation solved by the code for a scalar f is:
 *
 *    \f[ \rho*volume*\frac{df}{dt} + .... = st\_imp*f + st\_exp \f]
 *
 *
 *  Note that st_exp and st_imp are defined after the Finite Volume integration
 *  over the cells, so they include the "volume" term. More precisely:
 *    - st_exp is expressed in kg.[var]/s, where [var] is the unit of the variable.
 *      Its dimension is the one of the variable (3 for vectors)
 *    - st_imp is expressed in kg/s.
 *      Its dimension is 1 for scalars, 3x3 for vectors.
 *
 *  The st_exp and st_imp arrays are already initialized to 0 before entering the
 *  the routine. It is not needed to do it in the routine (waste of CPU time).
 *
 *  For stability reasons, Code_Saturne will not add -st_imp directly to the
 *  diagonal of the matrix, but Max(-st_imp,0). This way, the st_imp term is
 *  treated implicitely only if it strengthens the diagonal of the matrix.
 *  However, when using the second-order in time scheme, this limitation cannot
 *  be done anymore and -st_imp is added directly. The user should therefore test
 *  the negativity of st_imp by himself.
 *
 *  When using the second-order in time scheme, one should supply:
 *    - st_exp at time n
 *    - st_imp at time n+1/2
 *
 *
 *  The selection of cells where to apply the source terms is based on a getcel
 *  command. For more info on the syntax of the getcel command, refer to the
 *  user manual or to the comments on the similar command \ref getfbr in the routine
 *  \ref cs_user_boundary_conditions.
 *
 *  WARNING: If variable is the temperature, the resulting equation
 *           solved by the code is:
 *
 *   rho*Cp*volume*dT/dt + .... = st_imp*T + st_exp
 *
 *
 *  Note that st_exp and st_imp are defined after the Finite Volume integration
 *  over the cells, so they include the "volume" term. More precisely:
 *    - st_exp is expressed in W
 *    - st_imp is expressed in W/K
 *
 *
 *  STEP SOURCE TERMS
 * ===================
 *  In case of a complex, non-linear source term, say F(f), for variable f, the
 *  easiest method is to implement the source term explicitely.
 *
 *    df/dt = .... + F(f(n))
 *    where f(n) is the value of f at time tn, the beginning of the time step.
 *
 *  This yields :
 *    st_exp = volume*F(f(n))
 *    st_imp = 0
 *
 *  However, if the source term is potentially steep, this fully explicit
 *  method will probably generate instabilities. It is therefore wiser to
 *  partially implicit the term by writing:
 *
 *    df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
 *
 *  This yields:
 *    st_exp = volume*( F(f(n)) - dF/df*f(n) )
 *    st_imp = volume*dF/df
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t *domain,
                     int         f_id,
                     cs_real_t   *st_exp,
                     cs_real_t   *st_imp)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
