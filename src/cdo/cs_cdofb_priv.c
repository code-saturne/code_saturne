/*============================================================================
 * Common functions between scalar-valued and vector-valued CDO face-based
 * schemes
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_PRIV_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the advection-related parameters in the context structure of
 *         CDO face-based schemes
 *
 * \param[in]      eqp    pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb    pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] eqc    pointer to a \ref cs_cdofb_priv_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_set_advection_function(const cs_equation_param_t   *eqp,
                                cs_equation_builder_t       *eqb,
                                cs_cdofb_priv_t             *eqc)
{
  if (eqc == NULL || eqb == NULL)
    return;

  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  eqc->adv_func = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    const cs_xdef_t *adv_def = eqp->adv_field->definition;
    if (adv_def != NULL) { /* If linked to a NS equation, it might be null */
      if (adv_def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) {
        /* Required by cs_advection_field_cw_face_flux */
        eqb->msh_flag |= CS_FLAG_COMP_FEQ;
        eqb->msh_flag |= cs_quadrature_get_flag(adv_def->qtype,
                                                cs_flag_primal_face);
      }
    }

    /* Boundary conditions for advection */
    eqb->bd_msh_flag |= CS_FLAG_COMP_PFQ;

    switch (eqp->adv_formulation) {

    case CS_PARAM_ADVECTION_FORM_CONSERV:
      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        eqc->adv_func = cs_cdo_advection_fb_upwcsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_CENTERED:

        if (cs_equation_param_has_diffusion(eqp))
          eqc->adv_func = cs_cdo_advection_fb_cencsv;

        else {
          if (cs_equation_param_has_time(eqp) == false) {
            /* Remark 5 about static condensation of paper (DiPietro, Droniou,
             * Ern, 2015) */
            bft_error(__FILE__, __LINE__, 0,
                      " %s: Centered advection scheme is not valid option for"
                      " face-based discretization and without diffusion.",
                      __func__);
          }
          else
            eqc->adv_func = cs_cdo_advection_fb_cencsv;

        } /* else has diffusion */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for face-based discretization",
                  __func__);

      } /* Scheme */
      break; /* Conservative formulation */

    case CS_PARAM_ADVECTION_FORM_NONCONS:
      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        eqc->adv_func = cs_cdo_advection_fb_upwnoc;
        break;

      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        if (cs_equation_param_has_diffusion(eqp))
          eqc->adv_func = cs_cdo_advection_fb_cencsv;

        else {
          if (cs_equation_param_has_time(eqp) == false) {
            /* Remark 5 about static condensation of paper (DiPietro, Droniou,
             * Ern, 2015) */
            bft_error(__FILE__, __LINE__, 0,
                      " %s: Centered advection scheme is not valid option for"
                      " face-based discretization and without diffusion.",
                      __func__);
          }
          else
            eqc->adv_func = cs_cdo_advection_fb_cencsv;
        } /* else has diffusion */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for face-based discretization",
                  __func__);

      } /* Scheme */
      break; /* Non-conservative formulation */

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of formulation for the advection term",
                __func__);

    } /* Switch on the formulation */

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
