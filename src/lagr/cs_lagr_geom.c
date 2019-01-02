/*============================================================================
 * Methods for particle parameters
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

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_geom.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build geometric information needed by the deposition model.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_geom(void)
{
  cs_mesh_t            *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq  = cs_glob_mesh_quantities;

  BFT_REALLOC(cs_glob_lagr_b_u_normal, mesh->n_b_faces, cs_real_4_t);
  BFT_REALLOC(cs_glob_lagr_b_face_proj, mesh->n_b_faces, cs_real_33_t);

/* ==============================================================================*/
/* 1. Boundary faces planes computation     */
/* A X + B Y + C Z + D = 0   */
/* ==============================================================================*/

  for (cs_lnum_t ifac = 0; ifac < mesh->n_b_faces; ifac++) {

    /* Recover the first face nodes */
    cs_lnum_t v_id0  = mesh->b_face_vtx_lst[mesh->b_face_vtx_idx[ifac]];
    cs_lnum_t v_id1  = mesh->b_face_vtx_lst[mesh->b_face_vtx_idx[ifac] + 1];

    cs_real_t xs1     = mesh->vtx_coord[v_id0*3];
    cs_real_t ys1     = mesh->vtx_coord[v_id0*3+1];
    cs_real_t zs1     = mesh->vtx_coord[v_id0*3+2];

    cs_real_t xs2     = mesh->vtx_coord[v_id1*3];
    cs_real_t ys2     = mesh->vtx_coord[v_id1*3+1];
    cs_real_t zs2     = mesh->vtx_coord[v_id1*3+2];

    /* Face plane equation  */

    cs_real_t xnor = cs_math_3_norm(fvq->b_face_normal + ifac*3);

    cs_real_t xn   = fvq->b_face_normal[ifac*3]   / xnor;
    cs_real_t yn   = fvq->b_face_normal[ifac*3+1] / xnor;
    cs_real_t zn   = fvq->b_face_normal[ifac*3+2] / xnor;

    cs_glob_lagr_b_u_normal[ifac][0] = xn;
    cs_glob_lagr_b_u_normal[ifac][1] = yn;
    cs_glob_lagr_b_u_normal[ifac][2] = zn;

    cs_glob_lagr_b_u_normal[ifac][3] = -(xn*xs1 + yn*ys1 + zn*zs1);

    /* Matrix of Reference frame change    */

    xnor = sqrt (  (xs2 - xs1) * (xs2 - xs1)
                 + (ys2 - ys1) * (ys2 - ys1)
                 + (zs2 - zs1) * (zs2 - zs1));

    cs_real_t xt   = (xs2 - xs1) / xnor;
    cs_real_t yt   = (ys2 - ys1) / xnor;
    cs_real_t zt   = (zs2 - zs1) / xnor;

    cs_real_t xtt  = yn * zt - zn * yt;
    cs_real_t ytt  = zn * xt - xn * zt;
    cs_real_t ztt  = xn * yt - yn * xt;

    xnor = sqrt(xtt*xtt + ytt*ytt + ztt*ztt);

    xtt  = xtt / xnor;
    ytt  = ytt / xnor;
    ztt  = ztt / xnor;

    cs_glob_lagr_b_face_proj[ifac][0][0] = xn;
    cs_glob_lagr_b_face_proj[ifac][0][1] = yn;
    cs_glob_lagr_b_face_proj[ifac][0][2] = zn;
    cs_glob_lagr_b_face_proj[ifac][1][0] = xt;
    cs_glob_lagr_b_face_proj[ifac][1][1] = yt;
    cs_glob_lagr_b_face_proj[ifac][1][2] = zt;
    cs_glob_lagr_b_face_proj[ifac][2][0] = xtt;
    cs_glob_lagr_b_face_proj[ifac][2][1] = ytt;
    cs_glob_lagr_b_face_proj[ifac][2][2] = ztt;

  }

}
