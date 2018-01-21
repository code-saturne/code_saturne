#ifndef __CS_MESH_QUANTITIES_H__
#define __CS_MESH_QUANTITIES_H__

/*============================================================================
 * Management of mesh quantities
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_internal_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @defgroup bad_cells_flags Flags specifying bad cells treatment
 *
 * @{
 */

/*
 * Field property type
 */

/*! Correct bad cells warping for gradients */
#define CS_BAD_CELLS_WARPED_CORRECTION (1 << 0)

/*! Regularise on bad cells */
#define CS_BAD_CELLS_REGULARISATION (1 << 1)

/*! Recompute face centers */
#define CS_CELL_FACE_CENTER_CORRECTION (1 << 2)

/*! Recompute cell centers */
#define CS_CELL_CENTER_CORRECTION (1 << 3)

/*! Clip face distance when negative or too small */
#define CS_FACE_DISTANCE_CLIP (1 << 4)

/*! Clip geometrical quantities used in flux reconstuction */
#define CS_FACE_RECONSTRUCTION_CLIP (1 << 5)

/*! Limit cells volumle ratio */
#define CS_CELL_VOLUME_RATIO_CORRECTION (1 << 6)

/*! @} */

/*============================================================================
 * Type definition
 *============================================================================*/

/* Structure associated to mesh quantities management */

typedef struct {

  cs_real_t     *cell_cen;       /* Cell center coordinates  */
  cs_real_t     *cell_vol;       /* Cell volume */
  cs_real_t     *cell_f_vol;     /* Cell fluid volume */

  cs_real_t     *i_face_normal;  /* Surface normal of interior faces.
                                    (L2 norm equals area of the face) */
  cs_real_t     *b_face_normal;  /* Surface normal of border faces.
                                    (L2 norm equals area of the face) */
  cs_real_t     *i_f_face_normal;/* Fluid Surface normal of interior faces.
                                    (L2 norm equals area of the face) */
  cs_real_t     *b_f_face_normal;/* Fluid Surface normal of border faces.
                                    (L2 norm equals area of the face) */
  cs_real_t     *i_face_cog;     /* Center of gravity of interior faces */
  cs_real_t     *b_face_cog;     /* Center of gravity of border faces */

  cs_real_t     *i_face_surf;    /* Surface of interior faces. */
  cs_real_t     *b_face_surf;    /* Surface of boundary faces. */

  cs_real_t     *i_f_face_surf;  /* Fluid surface of interior faces. */
  cs_real_t     *b_f_face_surf;  /* Fluid surface of boundary faces. */

  cs_real_2_t   *i_f_face_factor;/* Fluid surface factor of interior faces. */
  cs_real_t     *b_f_face_factor;/* Fluid surface factor of boundary faces. */

  cs_real_t     *dijpf;          /* Vector I'J' for interior faces */
  cs_real_t     *diipb;          /* Vector II'  for border faces */
  cs_real_t     *dofij;          /* Vector OF   for interior faces */
  cs_real_t     *diipf;          /* Vector II'  for interior faces */
  cs_real_t     *djjpf;          /* Vector JJ'  for interior faces */

  cs_real_t     *i_dist;         /* Distance between the cell center and
                                    the center of gravity of interior faces */
  cs_real_t     *b_dist;         /* Distance between the cell center and
                                    the center of gravity of border faces */

  cs_real_t     *weight;         /* Interior faces weighting factor */

  cs_real_t      min_vol;        /* Minimum cell volume */
  cs_real_t      max_vol;        /* Maximum cell volume */
  cs_real_t      tot_vol;        /* Total volume */

  cs_real_t      min_f_vol;        /* Minimum cell volume */
  cs_real_t      max_f_vol;        /* Maximum cell volume */
  cs_real_t      tot_f_vol;        /* Total volume */

  cs_real_33_t  *cocgb_s_it;     /* coupling of gradient components for
                                    iterative reconstruction at boundary */
  cs_real_33_t  *cocg_s_it;      /* coupling of gradient components for
                                    iterative reconstruction */
  cs_real_33_t  *cocgb_s_lsq;    /* coupling of gradient components for
                                    least-square reconstruction at boundary */

  cs_real_33_t  *cocg_it;        /* Interleaved cocg matrix
                                    for iterative gradients */
  cs_real_33_t  *cocg_lsq;       /* Interleaved cocg matrix
                                    for least square gradients */

  cs_real_t     *corr_grad_lin_det;  /* Determinant of geometrical matrix
                                        linear gradient correction */
  cs_real_33_t  *corr_grad_lin;      /* Geometrical matrix
                                        linear gradient correction */

  cs_int_t      *b_sym_flag;     /* Symmetry flag for boundary faces */
  cs_int_t      *c_solid_flag;   /* Is the fluid volume 0 flag */
  unsigned      *bad_cell_flag;  /* Flag (mask) for bad cells detected */

} cs_mesh_quantities_t ;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to mesh quantities structure associated to the main mesh */

extern cs_mesh_quantities_t  *cs_glob_mesh_quantities;

/* Flag (mask) to activate bad cells correction */
extern unsigned cs_glob_mesh_quantities_flag;

/* Choice of the porous model */
extern int cs_glob_porous_model;

/*============================================================================
 * Public function prototypes for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for computing cell centers.
 *
 * This function returns 1 or 2 according to the selected algorithm.
 *
 * Fortran interface :
 *
 * SUBROUTINE ALGCEN (IOPT)
 * *****************
 *
 * INTEGER          IOPT        : <-> : Choice of the algorithm
 *                                      < 0 : query
 *                                        0 : computation based
 *                                            on faces (default choice)
 *                                        1 : computation based
 *                                            on vertices
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algcen, ALGCEN) (cs_int_t  *const iopt);

/*----------------------------------------------------------------------------
 * Set behavior for computing the cocg matrixes for the iterative algo
 * and for the Least square method for scalar and vector gradients.
 *
 * Fortran interface :
 *
 * subroutine comcoc (imrgra)
 * *****************
 *
 * integer          imrgra        : <-- : gradient reconstruction option
 *----------------------------------------------------------------------------*/

void
CS_PROCF (comcoc, COMCOC) (const cs_int_t  *const imrgra);

/*----------------------------------------------------------------------------
 * Set porous model
 *
 * Fortran interface :
 *
 * subroutine compor (iporos)
 * *****************
 *
 * integer          iporos        : <-- : porous model
 *----------------------------------------------------------------------------*/

void
CS_PROCF (compor, COMPOR) (const cs_int_t  *const iporos);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for computing cell centers.
 *
 *  < 0 : query
 *    0 : computation based on faces (default choice)
 *    1 : computation based on vertices
 *
 * algo_choice  <--  choice of algorithm to compute cell centers.
 *
 * returns:
 *  1 or 2 according to the selected algorithm.
 *----------------------------------------------------------------------------*/

int
cs_mesh_quantities_cell_cen_choice(const int algo_choice);

/*----------------------------------------------------------------------------
 * Compute cocg for iterative gradient reconstruction for scalars.
 *
 * parameters:
 *   gradient_option <-- gradient option (Fortran IMRGRA)
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_set_cocg_options(int  gradient_option);

/*----------------------------------------------------------------------------
 * Compute Fluid volumes and fluid surface in addition to cell volume and surfaces.
 *
 * parameters:
 *   porous_model <-- gradient option (Fortran iporos)
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_set_porous_model(int  porous_model);

/*----------------------------------------------------------------------------
 * Create a mesh quantities structure.
 *
 * returns:
 *   pointer to created cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

cs_mesh_quantities_t  *
cs_mesh_quantities_create(void);

/*----------------------------------------------------------------------------
 * Destroy a mesh quantities structure
 *
 * parameters:
 *   mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *
 * returns:
 *  NULL
 *----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_destroy(cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a mesh quantities structure to its empty initial state.
 *
 * \param[in]   mq           pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_free_all(cs_mesh_quantities_t  *mq);

/*----------------------------------------------------------------------------
 * Compute mesh quantities needed fo preprocessing
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute_preprocess(const cs_mesh_t       *mesh,
                                      cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute mesh quantities
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute(const cs_mesh_t       *mesh,
                           cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute fluid mesh quantities
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_fluid_compute(const cs_mesh_t       *mesh,
                                 cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute fluid section mesh quantities at the initial step
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_fluid_sections(const cs_mesh_t       *mesh,
                            cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute mesh quantities
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_sup_vectors(const cs_mesh_t       *mesh,
                               cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute internal and border face normal.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_normal <-> pointer to the internal face normal array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_face_normal(const cs_mesh_t   *mesh,
                               cs_real_t         *p_i_face_normal[],
                               cs_real_t         *p_b_face_normal[]);

/*----------------------------------------------------------------------------
 * Compute interior face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_cog    <-> pointer to the interior face center array
 *   p_i_face_normal <-> pointer to the interior face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_i_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_i_face_cog[],
                           cs_real_t         *p_i_face_normal[]);

/*----------------------------------------------------------------------------
 * Compute border face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_b_face_cog    <-> pointer to the border face center array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_b_face_cog[],
                           cs_real_t         *p_b_face_normal[]);

/*----------------------------------------------------------------------------
 * Compute cell centers.
 *
 * The corresponding array is allocated by this function, and it is the
 * caller's responsibility to free it when they are no longer needed.
 *
 * parameters:
 *   mesh       <-- pointer to a cs_mesh_t structure
 *   p_cell_cen <-> pointer to the cell centers array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_cell_cen(const cs_mesh_t  *mesh,
                            cs_real_t        *cell_cen[]);

/*----------------------------------------------------------------------------
 * Check that no negative volumes are present, and exit on error otherwise.
 *
 * parameters:
 *   mesh            <-- pointer to mesh structure
 *   mesh_quantities <-- pointer to mesh quantities structure
 *   allow_error     <-- 1 if errors are allowed, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_check_vol(const cs_mesh_t             *mesh,
                             const cs_mesh_quantities_t  *mesh_quantities,
                             int                          allow_error);

/*----------------------------------------------------------------------------
 * Update mesh quantities relative to extended ghost cells when the
 * neighborhood is reduced.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_reduce_extended(const cs_mesh_t       *mesh,
                                   cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Return the number of times mesh quantities have been computed.
 *
 * returns:
 *   number of times mesh quantities have been computed
 *----------------------------------------------------------------------------*/

int
cs_mesh_quantities_compute_count(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each vertex.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of smoothing passes
 * \param[out]  b_thickness  thickness for each mesh vertex
 *                           (0 at non-boundary vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_v(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each boundary face.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of optional smoothing passes
 * \param[out]  b_thickness  thickness for each mesh boundary face
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_f(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[]);

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_quantities_t structure
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_dump(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient least squares algorithm
 * adapted for internal coupling.
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *   ce   <->  coupling
 *----------------------------------------------------------------------------*/

void
cs_compute_cell_cocg_lsq_coupling(const cs_mesh_t         *m,
                                  cs_mesh_quantities_t    *fvq,
                                  cs_internal_coupling_t  *ce);

/*----------------------------------------------------------------------------
 * Compute 3x3 matrix cocg for the scalar gradient iterative algorithm
 * adapted for internal coupling.
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *   ce   <->  coupling
 *----------------------------------------------------------------------------*/

void
cs_compute_cell_cocg_it_coupling(const cs_mesh_t         *m,
                                 cs_mesh_quantities_t    *fvq,
                                 cs_internal_coupling_t  *ce);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_QUANTITIES_H__ */
