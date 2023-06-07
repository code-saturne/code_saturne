#ifndef __CS_BOUNDARY_CONDITIONS_H__
#define __CS_BOUNDARY_CONDITIONS_H__

/*============================================================================
 * Boundary condition handling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_base.h"
#include "cs_field.h"
#include "cs_function.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_time_control.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*! Maximum number of physical model zones
 *
 * \deprecated This is used for Fortran compatibilty (and maps to the Fortran
 * \ref ppppar:nozppm  "nozppm" parameter). In C, we should move to high
 * level boundary condition definitions not requiring indexing by legacy
 * zone numbers indexed by \ref.
*/
#define  CS_MAX_BC_PM_ZONE_NUM  2000

/*============================================================================
 * Local type definitions
 *============================================================================*/

typedef cs_real_t  cs_real_5_t[5];          /* Vector of 5 real values */

typedef cs_real_t  cs_real_5_20_t[5][20];   /* Matrix of 5x20 real values */

/*=============================================================================
 * Global variables
 *============================================================================*/

/*! Boundary condition type (code) associated with each boundary face */

extern const int  *cs_glob_bc_type;

/*! boundary zone number associated with each boundary face
 *  (specific physical models)
 *
 * \deprecated This is used for \ref cs_boundary_condition_pm_info_t only.
*/

/*----------------------------------------------------------------------------*/

/*! Legacy physical model boundary conditions.
 *
 * \remark The amppings of member arrays of this structure are shifted
 * by 1 when mapped to Fortran, so that both in Fortran and C, we can use
 * the natural indexing (1-based and 0-based respectively) without needing
 * to shift the zone id/number.
 *
 * \deprecated This should be used for Fortran compatibilty and migration
 * to C only. In C, we should then move to high level boundary condition
 * definitions not requiring indexing by legacy zone numbers indexed by
 * \ref cs_glob_bc_face_zone. */

typedef struct {

  /*! Legacy physical model zone id per boundary face */

  int  *izfppp;

  /*! indirection array allowing to sort the boundary faces
   *  according to their boundary condition type \c bc_type */
  int  *itrifb;

  /*! Imposed flow zone indicator (for inlet zones).
   * If the mass flow is imposed (\c iqimp(z_id) = 1), the matching
   * \c qimp value must be set, and the defined velocity boundary condition
   * will be rescaled so as to match the given mass flow (i.e. only its original
   * direction is used. Otherwise, the given velocity boundary condition
   * given by \c rcodcl1 is unchanged. */
  int   iqimp[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Turbulence inlet type:
    * - 0: given by the user
    * - 1: automatic, from hydraulic diameter and input velocity performed.
    * - 2: automatic, from turbulent intensity and input velocity performed.
    */
  int   icalke[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Imposed flow value (for inlet zones).
   * If the mass flow is imposed (\c iqimp(z_num - 1) = 1), the matching \c qimp
   * value must be set, and the defined velocity boundary condition will be
   * rescaled so as to match the given mass flow (i.e. only its original
   * direction is used. Otherwise, the given velocity boundary condition
   * given by \c rcodcl1 is unchanged. */
  cs_real_t  qimp[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! hydraulic diameter */
  cs_real_t  dh[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! turbulent intensity */
  cs_real_t  xintur[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Air temperature in K (per zone) for pulverized coal combustion. */
  cs_real_t  *timpat;

  /*! Coal mass flow per coal for inlet zones
    (qimpcp[zone_num][coal_id],  with 5 coals max) */

  cs_real_5_t  *qimpcp;

  /*! Coal temperature (in K) per coal for inlet zones
    (timpcp[zone_num][coal_id], with 5 coals max) */

  cs_real_5_t  *timpcp;

  /*! Coal class mass distribution ratio (in percentage) for inlet zones
   *  (distch[zone_num][coal_id][class_id], with 4 coals max and
      20 coal classes par coal max) */

  cs_real_5_20_t  *distch;

  /*! Air indicator by input facet type */
  int *ientat;

  /*! Cp indicator by input facet type */
  int *ientcp;

  /*! coal: number of oxydant for the current inlet */
  int *inmoxy;

  /*! gas combustion (cogz) */
  int ientfu[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for fuel flow inlet (gas combustion)

  int ientox[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for an air fow inlet (gas combustion)

  int ientgb[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for burned gas inlet (gas combustion)

  int ientgf[CS_MAX_BC_PM_ZONE_NUM+1]; // <-- 1 for unburned gas inlet (gas combustion)

  /*! inlet temperature (gas combustion) */
  double tkent[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! Mean Mixture Fraction at Inlet (gas combustion) */
  double fment[CS_MAX_BC_PM_ZONE_NUM+1];

  /*! atmo */
  /* atmospheric flows: auto inlet/outlet flag */
  cs_lnum_t *iautom;

  /* atmospheric flows: on/off for profile from data */
  int iprofm[CS_MAX_BC_PM_ZONE_NUM+1];

} cs_boundary_condition_pm_info_t;

extern cs_boundary_condition_pm_info_t  *cs_glob_bc_pm_info;

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
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(const cs_field_t          *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w);

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
  *a = -qimp/cs_math_fmax(hint, 1.e-300);
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
      bf[i][j] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set neumann BC for an anisotropic vector for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     qimpv         Flux value to impose
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_vector_aniso(cs_real_t        a[3],
                                                cs_real_t        af[3],
                                                cs_real_t        b[3][3],
                                                cs_real_t        bf[3][3],
                                                const cs_real_t  qimpv[3],
                                                const cs_real_t  hint[6])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};
  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  cs_real_t invdet = 1./(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  /* Gradient BCs */
  cs_math_sym_33_3_product(invh, qimpv, a);
  for (int isou = 0; isou < 3; isou++)
    a[isou] = -a[isou];

  b[0][0] = 1.0, b[0][1] = 0.0, b[0][2] = 0.0;
  b[1][0] = 0.0, b[1][1] = 1.0, b[1][2] = 0.0;
  b[2][0] = 0.0, b[2][1] = 0.0, b[2][2] = 1.0;

  for (int isou = 0; isou < 3; isou++) {

    /* Flux BCs */
    af[isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief  Set Neumann boundary conditions for a tensor for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     qimpts        Flux value to impose
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_tensor(cs_real_t        a[6],
                                          cs_real_t        af[6],
                                          cs_real_t        b[6][6],
                                          cs_real_t        bf[6][6],
                                          const cs_real_t  qimpts[6],
                                          cs_real_t        hint)
{
  for (int isou = 0; isou < 6; isou++) {

    /* Gradient BC */
    a[isou] = -qimpts[isou]/cs_math_fmax(hint, 1.e-300);
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        b[isou][jsou] = 1.0;
      else
        b[isou][jsou] = 0.0;
    }

    /* Flux BCs */
    af[isou] = qimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      bf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Dirichlet BC for a scalar for a given face.
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
  //if (fabs(hext) > cs_math_infinite_r*0.5) {
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
 * \brief  Set Dirichlet BC for a vector for a given face.
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
cs_boundary_conditions_set_dirichlet_vector(cs_real_t        a[3],
                                            cs_real_t        af[3],
                                            cs_real_t        b[3][3],
                                            cs_real_t        bf[3][3],
                                            const cs_real_t  pimpv[3],
                                            cs_real_t        hint,
                                            const cs_real_t  hextv[3])
{
  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {

      /* Gradient BCs */
      a[isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[isou][jsou] = 0.;

      /* Flux BCs */
      af[isou] = -hint*pimpv[isou];

      bf[0][0] = hint, bf[0][1] = 0.,   bf[0][2] = 0.;
      bf[1][0] = 0.,   bf[1][1] = hint, bf[1][2] = 0.;
      bf[2][0] = 0.,   bf[2][1] = 0.,   bf[2][2] = hint;

    }
    else {

      const cs_real_t val = hint/(hint + hextv[isou]);
      const cs_real_t heq = hextv[isou]*val;

      /* Gradient BCs */
      a[isou] = hextv[isou]*pimpv[isou]/(hint + hextv[isou]);

      b[0][0] = val, b[0][1] = 0.,  b[0][2] = 0.;
      b[1][0] = 0.,  b[1][1] = val, b[1][2] = 0.;
      b[2][0] = 0.,  b[2][1] = 0.,  b[2][2] = val;

      /* Flux BCs */
      af[isou] = -heq*pimpv[isou];

      bf[0][0] = heq, bf[0][1] = 0.,  bf[0][2] = 0.;
      bf[1][0] = 0.,  bf[1][1] = heq, bf[1][2] = 0.;
      bf[2][0] = 0.,  bf[2][1] = 0.,  bf[2][2] = heq;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective oulet BC for a scalar for a given face.
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
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t  *a,
                                                    cs_real_t  *af,
                                                    cs_real_t  *b,
                                                    cs_real_t  *bf,
                                                    cs_real_t   pimp,
                                                    cs_real_t   cfl,
                                                    cs_real_t   hint);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Acess the time control structure of an inlet.
 *
 * This allows modifying that structure, for example updating the inlet
 * velocity values only in a certain time range, and avoiding
 * uneeded recomputations outside that range.
 *
 * \param[in]  zone  pointer to associated zone
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t *
cs_boundary_conditions_get_inlet_time_control(const  cs_zone_t  *zone);

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
cs_boundary_conditions_set_dirichlet_vector_aniso(cs_real_t        a[3],
                                                  cs_real_t        af[3],
                                                  cs_real_t        b[3][3],
                                                  cs_real_t        bf[3][3],
                                                  const cs_real_t  pimpv[3],
                                                  const cs_real_t  hintt[6],
                                                  const cs_real_t  hextv[3])
{
  /* Gradient BCs */
  for (int isou = 0; isou < 3; isou++) {
    if (fabs(hextv[isou]) > 0.5*cs_math_infinite_r) {
      a[isou] = pimpv[isou];
      for (int jsou = 0; jsou < 3; jsou++)
        b[isou][jsou] = 0.;
    }
    else {
      /* FIXME: at least log error message */
      cs_exit(1);
    }
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, pimpv, af);
  for (int isou = 0; isou < 3; isou++)
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
 * \brief  Set Dirichlet BC for a tensor for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpts        Dirichlet value to impose
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     hextts        External exchange coefficient (10^30 by default)
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_tensor(cs_real_t        coefa[6],
                                            cs_real_t        cofaf[6],
                                            cs_real_t        coefb[6][6],
                                            cs_real_t        cofbf[6][6],
                                            const cs_real_t  pimpts[6],
                                            cs_real_t        hint,
                                            const cs_real_t  hextts[6])
{
  for (int isou = 0; isou < 6; isou++) {

    if (fabs(hextts[isou]) > 0.5*cs_math_infinite_r) {
      /* Gradient BCs */
      coefa[isou] = pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++)
        coefb[isou][jsou] = 0.;

      /* Flux BCs */
      cofaf[isou] = -hint * pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          cofbf[isou][jsou] = hint;
        else
          cofbf[isou][jsou] = 0.;
      }
    }
    else {

      const cs_real_t heq = hint * hextts[isou] / (hint + hextts[isou]);

      /* Gradient BCs */
      coefa[isou] = hextts[isou] * pimpts[isou] / (hint + hextts[isou]);
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          coefb[isou][jsou] = hint / (hint + hextts[isou]);
        else
          coefb[isou][jsou] = 0.;
      }

      /* Flux BCs */
      cofaf[isou] = -heq * pimpts[isou];
      for (int jsou = 0; jsou < 6; jsou++) {
        if (jsou == isou)
          cofbf[isou][jsou] = heq;
        else
          cofbf[isou][jsou] = 0.;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for a symmetric vector for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose on the normal
 *                              component
 * \param[in]     qimpv         Flux value to impose on the
 *                              tangential components
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_sym_vector(cs_real_t        coefa[3],
                                                  cs_real_t        cofaf[3],
                                                  cs_real_t        coefb[3][3],
                                                  cs_real_t        cofbf[3][3],
                                                  const cs_real_t  pimpv[3],
                                                  const cs_real_t  qimpv[3],
                                                  cs_real_t        hint,
                                                  const cs_real_t  normal[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    coefa[isou] = - qimpv[isou]/cs_math_fmax(hint, 1.e-300);
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++) {

      coefa[isou] = coefa[isou] + normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      if (jsou == isou)
        coefb[isou][jsou] = 1.0 - normal[isou] * normal[jsou];
      else
        coefb[isou][jsou] = - normal[isou] * normal[jsou];
    }

    /* Flux BCs */
    cofaf[isou] = qimpv[isou];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++){

      cofaf[isou] = cofaf[isou] - normal[isou]*normal[jsou]
                  * (hint * pimpv[jsou] + qimpv[jsou]);

      cofbf[isou][jsou] = hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose on the normal
 *                              component
 * \param[in]     qimpv         Flux value to impose on the
 *                              tangential components
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for a vector for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose on the tangential
 *                              components
 * \param[in]     qimpv         Flux value to impose on the
 *                              normal component
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_generalized_dirichlet_vector
  (cs_real_t         coefa[3],
   cs_real_t         cofaf[3],
   cs_real_t         coefb[3][3],
   cs_real_t         cofbf[3][3],
   const cs_real_t   pimpv[3],
   const cs_real_t   qimpv[3],
   cs_real_t         hint,
   const cs_real_t   normal[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BC*/
    /* "[1 -n(x)n] Pimp" is divided into two */
    coefa[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      coefa[isou] = coefa[isou] - normal[isou]*normal[jsou]
        * (pimpv[jsou] + qimpv[jsou] / cs_math_fmax(hint, 1.e-300));

      coefb[isou][jsou] = normal[isou] * normal[jsou];
    }

    /* Flux BC */
    /* "[1 -n(x)n] Pimp" is divided into two */
    cofaf[isou] = -hint*pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      cofaf[isou] = cofaf[isou] + normal[isou]*normal[jsou]
        * (qimpv[jsou] + pimpv[jsou] * hint);

      if (jsou == isou)
        cofbf[isou][jsou] = hint * (1.0 - normal[isou] * normal[jsou]);
      else
        cofbf[isou][jsou] = - hint * normal[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 * \param[in]     pimpv         Dirichlet value to impose on the tangential
 *                              components
 * \param[in]     qimpv         Flux value to impose on the
 *                              normal component
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for a vector for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose
 * \param[in]     cflv          Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector(cs_real_t       coefa[3],
                                                    cs_real_t       cofaf[3],
                                                    cs_real_t       coefb[3][3],
                                                    cs_real_t       cofbf[3][3],
                                                    const cs_real_t   pimpv[3],
                                                    const cs_real_t   cflv[3],
                                                    cs_real_t         hint)
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou ++) {
      if (jsou == isou)
        coefb[isou][jsou] = cflv[isou] / (1.0 + cflv[isou]);
      else
        coefb[isou][jsou] = 0.0;
    }
    coefa[isou] = pimpv[isou] * (1.0 - coefb[isou][isou]);

    /* Flux BCs */
    cofaf[isou] = -hint * coefa[isou];
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        cofbf[isou][jsou] = hint * (1.0 - coefb[isou][jsou]);
    else
      cofbf[isou][jsou] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective outlet BC for a tensor for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpts         Dirichlet value to impose
 * \param[in]     cflts          Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_tensor(cs_real_t        coefa[6],
                                                    cs_real_t        cofaf[6],
                                                    cs_real_t        coefb[6][6],
                                                    cs_real_t        cofbf[6][6],
                                                    const cs_real_t  pimpts[6],
                                                    const cs_real_t  cflts[6],
                                                    cs_real_t        hint)
{
  for (int isou = 0; isou < 6; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        coefb[isou][jsou] = cflts[isou] / (1.0 + cflts[isou]);
      else
        coefb[isou][jsou] = 0.0;
    }
    coefa[isou] = (1.0 - coefb[isou][isou]) * pimpts[isou];

    /* Flux BCs */
    cofaf[isou] = -hint*coefa[isou];
    for (int jsou = 0; jsou < 6; jsou++) {
      if (jsou == isou)
        cofbf[isou][jsou] = hint * (1.0 - coefb[isou][jsou]);
      else
        cofbf[isou][jsou] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set convective outlet BC for an anisotropic vector for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose
 * \param[in]     cflv          Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_convective_outlet_vector_aniso
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  cflv[3],
   const cs_real_t  hintt[6])
{
  for(int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    for (int jsou = 0; jsou < 3; jsou++) {
      if (jsou == isou)
        coefb[isou][jsou] = cflv[isou]/(1.0+cflv[isou]);
      else
        coefb[isou][jsou] = 0.0;
    }
  coefa[isou] = (1.0-coefb[isou][isou])*pimpv[isou];
  }

  /* Flux BCs */
  cs_math_sym_33_3_product(hintt, coefa, cofaf);
  for (int isou = 0; isou < 3; isou++)
    cofaf[isou] = -cofaf[isou];

  cofbf[0][0] = hintt[0]*(1.0 - coefb[0][0]);
  cofbf[1][1] = hintt[1]*(1.0 - coefb[1][1]);
  cofbf[2][2] = hintt[2]*(1.0 - coefb[2][2]);
  cofbf[0][1] = hintt[3]*(1.0 - coefb[0][0]);
  cofbf[1][0] = hintt[3]*(1.0 - coefb[0][0]);
  cofbf[1][2] = hintt[4]*(1.0 - coefb[1][1]);
  cofbf[2][1] = hintt[4]*(1.0 - coefb[1][1]);
  cofbf[0][2] = hintt[5]*(1.0 - coefb[2][2]);
  cofbf[2][0] = hintt[5]*(1.0 - coefb[2][2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set BC for an affine scalar function for a given face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pinf          affine part
 * \param[in]     ratio         linear part
 * \param[in]     hint          internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_scalar(cs_real_t  *a,
                                                  cs_real_t  *af,
                                                  cs_real_t  *b,
                                                  cs_real_t  *bf,
                                                  cs_real_t   pinf,
                                                  cs_real_t   ratio,
                                                  cs_real_t   hint)
{
  /* Gradient BCs */
  *b = ratio;
  *a = pinf;

  /* Flux BCs */
  *af = -hint * *a;
  *bf =  hint * (1. - *b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, zero flux for diffusion.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     dimp          Flux value to impose
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_neumann_conv_h_neumann_diff_scalar(cs_real_t  *a,
                                                              cs_real_t  *af,
                                                              cs_real_t  *b,
                                                              cs_real_t  *bf,
                                                              cs_real_t   dimp,
                                                              cs_real_t   hint)

{
  /* Gradient BCs */
  cs_boundary_conditions_set_neumann_scalar(&*a,
                                            &*af,
                                            &*b,
                                            &*bf,
                                            dimp,
                                            hint);

  /* Flux BCs */
  *af = 0.;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Neumann BC for the convection operator, imposed flux for
 *         diffusion.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pinf          affine part
 * \param[in]     ratio         linear part
 * \param[in]     dimp          Flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_affine_function_conv_neumann_diff_scalar
  (cs_real_t  *a,
   cs_real_t  *af,
   cs_real_t  *b,
   cs_real_t  *bf,
   cs_real_t   pinf,
   cs_real_t   ratio,
   cs_real_t   dimp)
{
  /* Gradient BCs */
  *b = ratio;
  *a = pinf;

  /* Flux BCs */
  *af = dimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set total flux as a Robin condition.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     hext          convective flux to be imposed
 * \param[in]     dimp          Flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_total_flux(cs_real_t *a,
                                      cs_real_t *af,
                                      cs_real_t *b,
                                      cs_real_t *bf,
                                      cs_real_t hext,
                                      cs_real_t dimp)
{
  /* Gradients BCs */
  *a = 0.;
  *b = 1.;

  /* Flux BCs */
  *af = dimp;
  *bf = hext;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a scalar.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimp          Dirichlet value to impose
 * \param[in]     dimp          Flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_scalar(cs_real_t *a,
                                                              cs_real_t *af,
                                                              cs_real_t *b,
                                                              cs_real_t *bf,
                                                              cs_real_t pimp,
                                                              cs_real_t dimp)
{
  /* Gradients BC */
  *a = pimp;
  *b = 0.;

  /* Flux BCs */
  *af = dimp;
  *bf = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a vector.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose
 * \param[in]     qimpv         Flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_vector
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3])
{
  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    coefa[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      coefb[isou][jsou] = 0.0;


    /* Flux BCs */
    cofaf[isou] = qimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++)
      cofbf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Imposed value for the convection operator, imposed flux for
 *         diffusion, for a tensor
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpts         Dirichlet value to impose
 * \param[in]     qimpts         Flux value to impose
 */
/*----------------------------------------------------------------------------*/

inline static void
cs_boundary_conditions_set_dirichlet_conv_neumann_diff_tensor
  (cs_real_t        coefa[6],
   cs_real_t        cofaf[6],
   cs_real_t        coefb[6][6],
   cs_real_t        cofbf[6][6],
   const cs_real_t  pimpts[6],
   const cs_real_t  qimpts[6])
{
  for (int isou = 0; isou < 6; isou++) {

    /* BS test sur hextv ? if (abs(hextv[isou])>rinfin*0.5d0) then */

    /* Gradient BCs */
    coefa[isou] = pimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      coefb[isou][jsou] = 0.0;

    /* Flux BCs */
    cofaf[isou] = qimpts[isou];
    for (int jsou = 0; jsou < 6; jsou++)
      cofbf[isou][jsou] = 0.0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare (reset) condition coefficients for all variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_reset(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update per variable boundary condition codes.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_compute(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic adjustments for boundary condition codes.
 *
 * Currently handles mapped inlets, after the call to \ref stdtcl.
 * As portions of stdtcl are migrated to C, they should be called here,
 * before mapped inlets.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_complete(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant velocity to an open (inlet/outlet) boundary.
 *
 * This function may also be used to define the flow direction if called
 * before one of the \c cs_boundary_conditions_open_set_mass_flow_rate
 * or \c cs_boundary_conditions_set_volume_flow_rate functions.
 *
 * \param[in]  z       pointer to associated zone
 * \param[in]  u_norm  associated constant normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_value(const cs_zone_t  *z,
                                                  const cs_real_t   u[3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant velocity normal to an inlet.
 *
 * \param[in]  z       pointer to associated zone
 * \param[in]  u_norm  associated constant normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_normal_value(const  cs_zone_t  *z,
                                                         cs_real_t     u_norm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a normal velocity to an inlet using a provided function.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated velocity vector evaluation function
 *                    at zone faces
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_func(const  cs_zone_t       *z,
                                                 cs_eval_at_location_t  *func,
                                                 void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant mass flow rate to an inlet.
 *
 * By default, the flow direction is considered normal to the boundary.
 * The flow direction may be specified by calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func
 * for the appropriate zone before calling this function.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_value(const  cs_zone_t  *z,
                                                        cs_real_t          q);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a mass flow rate to an inlet based on provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (mass flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a constant volume flow rate to an inlet.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant volume flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_value(const  cs_zone_t  *z,
                                                          cs_real_t          q);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a volume flow rate to an inlet based on provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (volume flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign a flow rate to an inlet based directly on the velocity.
 *
 * This is the default, so this function is useful only if need to restore
 * that behavior if needed after calling one of the
 * \c cs_boundary_conditions_set_volume_volume_flow_rate
 * \c cs_boundary_conditions_open_set_volume_flow_rate functions.
 *
 * \param[in]  z      pointer to associated zone
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_flow_rate_by_velocity(const  cs_zone_t  *z);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_H__ */
