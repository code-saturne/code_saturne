/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_MESH_WARPING_H__
#define __CS_MESH_WARPING_H__

/*============================================================================
 * Cut warped faces in serial or parallel with/without periodicity.
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set the threshold to cut warped faces.
 *
 * Fortran interface :
 *
 * subroutine setcwf (cwfthr)
 * *****************
 *
 * integer          cwfpst      : <-> : if 1, activate postprocessing when
 *                                      cutting warped faces (default 0)
 * double precision cwfthr      : <-> : threshold angle (in degrees) if
 *                                      positive, do not cut warped faces
 *                                      if negative (default -1)
 *----------------------------------------------------------------------------*/

void
CS_PROCF (setcwf, SETCWF) (const cs_int_t   *cwfpst,
                           const cs_real_t  *cwfthr);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Cut warped faces.
 *
 * Updates border face connectivity and associated mesh quantities.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure.
 *   max_warp_angle <-- criterion to know which face to cut
 *   post_flag      <-- 1 if we have to post-process cut faces, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_cut_faces(cs_mesh_t    *mesh,
                          double        max_warp_angle,
                          cs_bool_t     post_flag);

/*----------------------------------------------------------------------------
 * Set defaults for cutting of warped faces.
 *
 * parameters:
 *   max_warp_angle <-- maximum warp angle (in degrees) over which faces will
 *                      be cut; negative (-1) if faces should not be cut
 *   postprocess    <-- 1 if postprocessing should be activated when cutting
 *                      warped faces, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_set_defaults(double  max_warp_angle,
                             int     postprocess);

/*----------------------------------------------------------------------------
 * Get defaults for cutting of warped faces.
 *
 * parameters:
 *   max_warp_angle --> if non NULL, returns maximum warp angle (in degrees)
 *                      over which faces will be cut, or -1 if faces should
 *                      not be cut
 *   postprocess    --> if non NULL, returns 1 if postprocessing should be
 *                      activated when cutting warped faces, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_warping_get_defaults(double  *max_warp_angle,
                             int     *postprocess);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_WARPING_H__ */
