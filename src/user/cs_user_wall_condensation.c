/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*
===============================================================================
 Function:
 ---------

> \file cs_user_wall_condensation.c
>
> \brief Source terms associated at the boundary faces and the neighboring
> cells with surface condensation.
>
> This subroutine fills the condensation source terms for each variable at
> the cell center associated to the boundary faces identifed in the mesh.
> The fluid exchange coefficient is computed with a empiric law to be
> imposed at the boundary face where the condensation phenomenon occurs.
>
> This user subroutine is called which allows the setting of
> \f$ \gamma_{\mbox{cond}} \f$ the condensation source term.
>
> This function fills the condensation source term array gamma_cond adding
> to the following equations:
>
> - The equation for mass conservation:
> \f[ D\frac{rho}{dt} + divs \left( \rho \vect{u}^n\right) = \Gamma _{cond}
> \f]
>
> - The equation for a variable \f$\Phi \f$:
> \f[ D\frac{\phi}{dt} = ... + \Gamma _{cond}*(\Phi _i - \Phi)
> \f]
>
> discretized as below:
>
> \f[ \rho*\dfrac{\Phi^{n+1}-\Phi^{n}}/dt = ...
>                            + \Gamma _{cond}*(\Phi _i - \Phi^{n+1})
> \f]
>
> \remarks
>  - \f$ \Phi _i \f$ is the value of \f$ \Phi \f$ associated to the
>    injected condensation rate.
>
>    With 2 options are available:
>       - the condensation rate is injected with the local value
>         of variable \f$ \Phi = \Phi ^{n+1}\f$
>         in this case the \f$ \Phi \f$ variable is not modified.
>
>       - the condensation rate is injected with a specific value
>         for \f$ \Phi = \Phi _i \f$ the specified value given by the
>         user.
>
> \section use Usage
>
> The three stages in the code where this User subroutine
> is called (with \code iappel = 1, 2 and 3\endcode)
>
> \code iappel = 1 \endcode
>  - Calculation of the number of cells where a mass source term is
>    imposed: ncesmp
>    Called once at the beginning of the calculation
>
> \code iappel = 2 \endcode
>   - Identification of the cells where a mass source term is imposed:
>     array icesmp(ncesmp)
>     Called once at the beginning of the calculation
>
> \code iappel = 3 \endcode
>   - Calculation of the values of the mass source term
>     Called at each time step
>
> \section the specific variables to define with is user subroutine
>
>  - ifbpcd(ieltcd): identification of the faces where a condensation
>                    source term is imposed.
>
>  - itypcd(ieltcd,ivar): type of treatment for variable ivar in the
>                       ieltcd cell with condensation source term.
>                     - itypcd = 0 --> injection of ivar at local value
>                     - itypcd = 1 --> injection of ivar at user
>                                      specified value.
>
>  - spcond(ielscd,ipr): value of the injection condensation rate
>                       gamma_cond (kg/m3/s) in the ieltcd cell
>                       with condensation source term.
>
>  - spcond(ieltcd,ivar): specified value for variable ivar associated
>                        to the injected condensation in the ieltcd
>                        cell with a condensation source term except
>                        for ivar=ipr.
>
> \remarks
>  - For each face where a condensation source terms is imposed ielscd
>    in [1;nfbpcd]), ifbpcd(ielscd) is the global index number of the
>    corresponding face (ifbpcd(ieltcd) in [1;ncel]).
>  - if itypcd(ieltcd,ivar)=0, spcond(ielpcd,ivar) is not used.
>  - if spcond(ieltcd,ipr)<0, mass is removed from the system,
>     therefore Code_Saturna automatically considers f_i=f^(n+1),
>     whatever the values of itypcd or smacel specified by the user
>
>   \par Examples of settings for boundary condensation mass source terms
>        Examples are available
>        \ref condens_h_boundary "here".
>
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
 Arguments
______________________________________________________________________________.
  mode           name          role                                           !
______________________________________________________________________________!
> \param[in]     nvar          total number of variables
> \param[in]     nscal         total number of scalars
> \param[in]     iappel        indicates which at which stage the routine is
_______________________________________________________________________________
*/

#pragma weak cs_user_wall_condensation
void
cs_user_wall_condensation(int nvar, int nscal, int iappel)
{
  CS_UNUSED(nvar);
  CS_UNUSED(nscal);
  CS_UNUSED(iappel);
};

/*----------------------------------------------------------------------------*/

END_C_DECLS
