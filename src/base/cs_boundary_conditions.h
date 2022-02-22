#ifndef __CS_BOUNDARY_CONDITIONS_H__
#define __CS_BOUNDARY_CONDITIONS_H__

/*============================================================================
 * Boundary condition handling.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

#include "fvm_nodal.h"
#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*! Boundary condition type (code) associated with each boundary face */

extern const int  *cs_glob_bc_type;

/*! boundary zone number associated with each boundary face
  (specific physical models)*/

extern const int  *cs_glob_bc_face_zone;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Handling of boundary condition definition errors and associated output.
 *
 * This function checks for errors, and simply returns if no error is
 * encountered. In case of error, it outputs helpful information so as to
 * make it easier to locate the matching faces.
 *
 * For each boundary face, bc_type defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 * An optional label may be used if the error is related to another
 * attribute than the boundary type, for appropriate error reporting.
 *
 * parameters:
 *   bc_flag   <-- array of BC type ids
 *   type_name <-- name of attribute in error, or NULL
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int   *bc_type,
                             const char  *type_name);

/*----------------------------------------------------------------------------
 * Locate shifted boundary face coordinates on possibly filtered
 * cells or boundary faces for later interpolation.
 *
 * parameters:
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS
 *                        or CS_MESH_LOCATION_BOUNDARY_FACES)
 *   n_location_elts <-- number of selected location elements
 *   n_faces         <-- number of selected boundary faces
 *   location_elts   <-- list of selected location elements (0 to n-1),
 *                       or NULL if no indirection is needed
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or NULL if no indirection is needed
 *   coord_shift     <-- array of coordinates shift relative to selected
 *                       boundary faces
 *   coord_stride    <-- access stride in coord_shift: 0 for uniform
 *                       shift, 1 for "per face" shift.
 *   tolerance       <-- relative tolerance for point location.
 *
 * returns:
 *   associated locator structure
 *----------------------------------------------------------------------------*/

ple_locator_t *
cs_boundary_conditions_map(cs_mesh_location_type_t    location_type,
                           cs_lnum_t                  n_location_elts,
                           cs_lnum_t                  n_faces,
                           const cs_lnum_t           *location_elts,
                           const cs_lnum_t           *faces,
                           cs_real_3_t               *coord_shift,
                           int                        coord_stride,
                           double                     tolerance);

/*----------------------------------------------------------------------------
 * Set mapped boundary conditions for a given field and mapping locator.
 *
 * parameters:
 *   field           <-- field whose boundary conditions are set
 *   locator         <-- associated mapping locator, as returned
 *                       by cs_boundary_conditions_map().
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS or
 *                       CS_MESH_LOCATION_BOUNDARY_FACES)
 *   normalize       <-- normalization option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   interpolate     <-- interpolation option:
 *                         0: values are simply based on matching
 *                            cell or face center values
 *                         1: values are based on matching cell
 *                            or face center values, corrected
 *                            by gradient interpolation
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
 *   nvar            <-- number of variables requiring BC's
 *   rcodcl          <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(const cs_field_t          *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w,
                                  int                        nvar,
                                  cs_real_t                  rcodcl[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add location of locate shifted boundary face coordinates on
 *        cells or boundary faces for automatic interpolation.
 *
 * \note
 * This function is currently restricted to mapping of boundary face
 * locations (usually from boundary zones) to cell of boundary face
 * locations, but could be extended to other location types in the future.
 *
 * \param[in]  bc_location_id      id of selected boundary mesh location;
 *                                 currently restricted to subsets of
 *                                 boundary faces (i.e. boundary zone
 *                                 location ids).
 * \param[in]  source_location_id  id of selected location  mesh location
 *                                 (usually CS_MESH_LOCATION_CELLS but can be
 *                                 a more restricted cell or boundary face zone
 *                                 location location id).
 * \param[in]  coord_shift      coordinates shift relative to selected
 *                              boundary faces
 * \param[in]  tolerance        relative tolerance for point location.
 *
 * \return  id of added map
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_conditions_add_map(int         bc_location_id,
                               int         source_location_id,
                               cs_real_t   coord_shift[3],
                               double      tolerance);

/*----------------------------------------------------------------------------
 * Create the boundary conditions face type and face zone arrays
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create(void);

/*----------------------------------------------------------------------------
 * Free the boundary conditions face type and face zone arrays.
 *
 * This also frees boundary condition mappings which may have been defined.
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_free(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Neumann BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   qimp   flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_scalar(cs_real_t  *a,
                                          cs_real_t  *af,
                                          cs_real_t  *b,
                                          cs_real_t  *bf,
                                          cs_real_t   qimp,
                                          cs_real_t   hint)
{
  /* Gradient BCs */
  *a = -qimp/hint;
  *b = 1.;

  /* Flux BCs */
  *af = qimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Neumann BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   qimpv  flux value to impose
 * \param[in]   hint   internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector(cs_real_t        a[3],
                                          cs_real_t        af[3],
                                          cs_real_t        b[3][3],
                                          cs_real_t        bf[3][3],
                                          const cs_real_t  qimpv[3],
                                          cs_real_t        hint)
{
  /* Gradient BCs */

  for (size_t i = 0; i < 3; i++) {
    a[i] = -qimpv[i] / fmax(hint, 1.e-300);
  }

  b[0][0] = 1., b[0][1] = 0., b[0][2] = 0.;
  b[1][0] = 0., b[1][1] = 1., b[1][2] = 0.;
  b[2][0] = 0., b[2][1] = 0., b[2][2] = 1.;

  /* Flux BCs */

  for (size_t i = 0; i < 3; i++) {
    af[i] = qimpv[i];

    for (size_t j = 0; j < 3; j++)
      bf[i][j] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a scalar for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimp   dirichlet value to impose
 * \param[in]   hint   internal exchange coefficient
 * \param[in]   hext   external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_scalar(cs_real_t  *a,
                                            cs_real_t  *af,
                                            cs_real_t  *b,
                                            cs_real_t  *bf,
                                            cs_real_t   pimp,
                                            cs_real_t   hint,
                                            cs_real_t   hext)
{
  if (hext < 0.) {

    /* Gradient BCs */
    *a = pimp;
    *b = 0.;

    /* Flux BCs */
    *af = -hint*pimp;
    *bf =  hint;

  }
  else {

    /* Gradient BCs */
    *a = hext*pimp/(hint + hext);
    *b = hint     /(hint + hext);

    /* Flux BCs */
    cs_real_t heq = hint*hext/(hint + hext);
    *af = -heq*pimp;
    *bf =  heq;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a vector for a given face.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimpv  dirichlet value to impose
 * \param[in]   hint   internal exchange coefficient
 * \param[in]   hextv  external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector(cs_real_3_t    a,
                                            cs_real_3_t    af,
                                            cs_real_33_t   b,
                                            cs_real_33_t   bf,
                                            cs_real_3_t    pimpv,
                                            cs_real_t      hint,
                                            cs_real_3_t    hextv)
{
  for (int isou = 0 ; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[isou] = pimpv[isou];
      for (int jsou = 0 ; jsou < 3 ; jsou++)
        b[isou][jsou] = 0.;

      /* Flux BCs */
      af[isou] = -hint*pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++) {
        if (jsou == isou)
          bf[isou][jsou] = hint;
        else
          bf[isou][jsou] = 0.;
      }
    } else {

      cs_real_t heq = hint*hextv[isou]/(hint + hextv[isou]);

      /* Gradient BCs */
      a[isou] = hextv[isou]*pimpv[isou]/(hint + hextv[isou]);
      for (int jsou = 0 ; jsou < 3 ; jsou++) {
        if (jsou == isou)
          b[isou][jsou] = hint/(hint + hextv[isou]);
        else
          b[isou][jsou] = 0.;
      }

      /* Flux BCs */
      af[isou] = -heq*pimpv[isou];
      for (int jsou = 0 ; jsou < 3 ; jsou++) {
        if (jsou == isou)
          bf[isou][jsou] = heq;
        else
          bf[isou][jsou] = 0.;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective oulet boundary condition for a scalar
 *
 * Parameters:
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimp          Flux value to impose
 * \param[in]     cfl           Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t *coefa,
                                                    cs_real_t *cofaf,
                                                    cs_real_t *coefb,
                                                    cs_real_t *cofbf,
                                                    cs_real_t  pimp,
                                                    cs_real_t  cfl,
                                                    cs_real_t  hint);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set Dirichlet BC for a vector for a given face with left anisotropic
 *        diffusion.
 *
 * \param[out]  a      explicit BC coefficient for gradients
 * \param[out]  af     explicit BC coefficient for diffusive flux
 * \param[out]  b      implicit BC coefficient for gradients
 * \param[out]  bf     implicit BC coefficient for diffusive flux
 * \param[in]   pimpv  dirichlet value to impose
 * \param[in]   hintt  internal exchange coefficient
 * \param[in]   hextv  external exchange coefficient
 *                     (assumed infinite/ignored if < 0)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_vector_aniso(cs_real_3_t    a,
                                                  cs_real_3_t    af,
                                                  cs_real_33_t   b,
                                                  cs_real_33_t   bf,
                                                  cs_real_3_t    pimpv,
                                                  cs_real_6_t    hintt,
                                                  cs_real_3_t    hextv)
{
  for (int isou = 0 ; isou < 3 ; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[isou] = pimpv[isou];
      for (int jsou = 0 ; jsou < 3 ; jsou++)
        b[isou][jsou] = 0.;

    } else {

      cs_exit(1);

    }
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, pimpv, af);
  for (int isou = 0 ; isou < 3 ; isou++)
    af[isou] = -af[isou];

  bf[0][0] = hintt[0];
  bf[1][1] = hintt[1];
  bf[2][2] = hintt[2];
  bf[0][1] = hintt[3];
  bf[1][0] = hintt[3];
  bf[1][2] = hintt[4];
  bf[2][1] = hintt[4];
  bf[0][2] = hintt[5];
  bf[2][0] = hintt[5];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update per variable boundary condition codes.
 *
 * \param[in]       nvar             number of variables requiring BC's
 * \param[in]       itypfb           type of boundary for each face
 * \param[in, out]  icodcl           boundary condition codes
 * \param[in, out]  rcodcl           boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_compute(int     nvar,
                               int     itypfb[],
                               int     icodcl[],
                               double  rcodcl[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic adjustments for boundary condition codes.
 *
 * Currently handles mapped inlets, after the call to \ref stdtcl.
 * As portions of stdtcl are migrated to C, they should be called here,
 * before mapped inlets.
 *
 * \param[in]       nvar             number of variables requiring BC's
 * \param[in]       itypfb           type of boundary for each face
 * \param[in, out]  icodcl           boundary condition codes
 * \param[in, out]  rcodcl           boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_complete(int     nvar,
                                int     itypfb[],
                                int     icodcl[],
                                double  rcodcl[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_H__ */
