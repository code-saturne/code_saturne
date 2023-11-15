#ifndef __CS_FLAG_H__
#define __CS_FLAG_H__

/*============================================================================
 * Manage the definition/setting of a computation
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*!
 * @name Flags specifying the type of cell
 * @{
 *
 * \def CS_FLAG_BOUNDARY_CELL_BY_FACE
 * \brief (= 1) boundary cell with at least one border face
 *
 * \def CS_FLAG_BOUNDARY_CELL_BY_VERTEX
 * \brief (= 2) boundary cell with at least one border vertex
 *
 * \def CS_FLAG_BOUNDARY_CELL_BY_EDGE
 * \brief (= 4) boundary cell with at least one border edge
 *
 * \def CS_FLAG_SOLID_CELL
 * \brief (= 8) cell attached to a solid zone during the calculation
 *
 */

#define CS_FLAG_BOUNDARY_CELL_BY_FACE     (1 << 0)
#define CS_FLAG_BOUNDARY_CELL_BY_VERTEX   (1 << 1)
#define CS_FLAG_BOUNDARY_CELL_BY_EDGE     (1 << 2)
#define CS_FLAG_SOLID_CELL                (1 << 3)

/*!
 * @name Flags specifying the type of update to perform
 * @{
 *
 * \def CS_FLAG_CURRENT_TO_PREVIOUS
 * \brief (value = 1) perform a previous to current operation on fields during
 *        the update step
 *
 * \def CS_FLAG_INITIALIZATION
 * \brief (value = 2) This is the first time that the update is called so it is
 *        also an initialization step
 */

#define CS_FLAG_CURRENT_TO_PREVIOUS       (1 << 0)
#define CS_FLAG_INITIALIZATION            (1 << 1)

/*!
 * @}
 * @name Flag related to the way a CDO system is built
 * @{
 */

#define CS_FLAG_SYS_MASS_MATRIX  (1 << 0) /*!<  1: build a mass matrix */
#define CS_FLAG_SYS_SYM          (1 << 1) /*!<  2: system matrix is symmetric */
#define CS_FLAG_SYS_TIME_DIAG    (1 << 2) /*!<  4: lump the time term */
#define CS_FLAG_SYS_REAC_DIAG    (1 << 3) /*!<  8: lump the reaction term */
#define CS_FLAG_SYS_SOURCES_HLOC (1 << 4) /*!< 16: Hodge op. for source terms */
#define CS_FLAG_SYS_VECTOR       (1 << 5) /*!< 32: vector-valued variable */

/*!
 * @}
 * @name Flags use to identify the nature/status of an object
 * Apply this flag on a variable, a property or an advection field
 * @{
 */

#define CS_FLAG_STATE_UNIFORM      (1 << 0) /*!<   1: uniform (in space) */
#define CS_FLAG_STATE_CELLWISE     (1 << 1) /*!<   2: cellwise uniform */
#define CS_FLAG_STATE_FACEWISE     (1 << 2) /*!<   4: uniform on each face */
#define CS_FLAG_STATE_STEADY       (1 << 3) /*!<   8: steady */
#define CS_FLAG_STATE_POTENTIAL    (1 << 4) /*!<  16: potential */
#define CS_FLAG_STATE_CIRCULATION  (1 << 5) /*!<  32: circulation */
#define CS_FLAG_STATE_FLUX         (1 << 6) /*!<  64: flux */
#define CS_FLAG_STATE_DENSITY      (1 << 7) /*!< 128: density */

/*!
 * @}
 * @name Flags use to identify metadata for a quantity
 * where is located this quantity and how to access to its values
 * @{
 */

#define CS_FLAG_FULL_LOC  (1 <<  0) /*!<    1: defined on the whole location */
#define CS_FLAG_SCALAR    (1 <<  1) /*!<    2: scalar-valued (stride = 1) */
#define CS_FLAG_VECTOR    (1 <<  2) /*!<    4: vector-valued (stride = 3) */
#define CS_FLAG_TENSOR    (1 <<  3) /*!<    8: tensor-valued (stride = 9) */
#define CS_FLAG_VERTEX    (1 <<  4) /*!<   16: on vertices */
#define CS_FLAG_EDGE      (1 <<  5) /*!<   32: on edges */
#define CS_FLAG_FACE      (1 <<  6) /*!<   64: on faces */
#define CS_FLAG_CELL      (1 <<  7) /*!<  128: on cells */
#define CS_FLAG_PRIMAL    (1 <<  8) /*!<  256: on primal mesh */
#define CS_FLAG_DUAL      (1 <<  9) /*!<  512: on dual mesh */
#define CS_FLAG_BORDER    (1 << 10) /*!< 1024: located on the boundary */
#define CS_FLAG_BY_CELL   (1 << 11) /*!< 2048: by cell (c2e, c2f, c2v) */
#define CS_FLAG_BY_FACE   (1 << 12) /*!< 4096: by face (bf2v) */

/*!
 * @}
 * @name Flags use to identify the type of equations requested.

 * There exists one flag for each type of space discretization.
 * If this flag is activated, then at least one equation or system of this
 * kind is discretized with this space discretization scheme.
 * @{
 */

#define CS_FLAG_SCHEME_SCALAR  (1 << 0) /*!<  1: scheme for scalar eq. */
#define CS_FLAG_SCHEME_VECTOR  (1 << 1) /*!<  2: scheme for a vector eq. */
#define CS_FLAG_SCHEME_NAVSTO  (1 << 2) /*!<  4: Navier-Stokes system */
#define CS_FLAG_SCHEME_POLY0   (1 << 3) /*!<  8: lowest-order scheme */
#define CS_FLAG_SCHEME_POLY1   (1 << 4) /*!< 16: Linear gradient approx. */
#define CS_FLAG_SCHEME_POLY2   (1 << 5) /*!< 32: Quadratic gradient approx. */

/*!
 * @}
 */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef unsigned char       cs_mask_t;   /* Elementary flag */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Default locations */

extern const cs_flag_t  cs_flag_primal_vtx;
extern const cs_flag_t  cs_flag_primal_edge;
extern const cs_flag_t  cs_flag_primal_face;
extern const cs_flag_t  cs_flag_primal_cell;
extern const cs_flag_t  cs_flag_vertex;      /* equal to cs_flag_primal_vtx */
extern const cs_flag_t  cs_flag_cell;        /* equal to cs_flag_primal_cell */
extern const cs_flag_t  cs_flag_boundary_face;

extern const cs_flag_t  cs_flag_dual_vtx;
extern const cs_flag_t  cs_flag_dual_face;
extern const cs_flag_t  cs_flag_dual_cell;

extern const cs_flag_t  cs_flag_primal_edge_byc;
extern const cs_flag_t  cs_flag_dual_cell_byc;
extern const cs_flag_t  cs_flag_dual_face_byc;

/* Part of dual cell closure belonging to a boundary primal face */

extern const cs_flag_t  cs_flag_dual_closure_byf;

/* According to the extended flag defined below one can identify which set of
 * quantities or connectivities have to be built on-the-fly and stored in a
 * local structure possibly owned by each thread and with a cellwise scope */

typedef unsigned int  cs_eflag_t;

/* Store predefined flags */

extern const cs_eflag_t  cs_flag_need_v;
extern const cs_eflag_t  cs_flag_need_e;
extern const cs_eflag_t  cs_flag_need_f;
extern const cs_eflag_t  cs_flag_need_fe;
extern const cs_eflag_t  cs_flag_need_ef;
extern const cs_eflag_t  cs_flag_need_peq;
extern const cs_eflag_t  cs_flag_need_dfq;
extern const cs_eflag_t  cs_flag_need_pfq;
extern const cs_eflag_t  cs_flag_need_deq;
extern const cs_eflag_t  cs_flag_need_pfc;

typedef enum {

  /* Compute simple and cellwise information for vertices */

  CS_FLAG_COMP_PV      = 1 << 0,   /* = 1 */

  /* Compute cellwise quantities for vertices */

  CS_FLAG_COMP_PVQ     = 1 << 1,   /* = 2 */

  /* Compute simple and cellwise information for edges */

  CS_FLAG_COMP_PE      = 1 << 2,   /* = 4 */

  /* Compute cellwise quantities for edges */

  CS_FLAG_COMP_PEQ     = 1 << 3,   /* = 8 */

  /* Compute cellwise quantities for dual faces (associated to edges) */

  CS_FLAG_COMP_DFQ     = 1 << 4,   /* = 16 */

  /* Compute simple and cellwise information for faces */

  CS_FLAG_COMP_PF      = 1 << 5,   /* = 32 */

  /* Compute cellwise quantities for faces */

  CS_FLAG_COMP_PFQ     = 1 << 6,   /* = 64 */

  /* Compute cellwise quantities for dual edges (associated to faces) */

  CS_FLAG_COMP_DEQ     = 1 << 7,   /* = 128 */

  /* Compute the cellwise connectivity edge to vertices */

  CS_FLAG_COMP_EV      = 1 << 8,   /* = 256 */

  /* Compute cellwise connectivity face to edges */

  CS_FLAG_COMP_FE      = 1 << 9,   /* = 512 */

  /* Compute cellwise quantities associated to the couple (face, edge) */

  CS_FLAG_COMP_FEQ     = 1 << 10,  /* = 1024 */

  /* Compute cellwise connectivity face to vertices */

  CS_FLAG_COMP_FV      = 1 << 11,  /* = 2048 */

  /* Compute cellwise connectivity edge to faces */

  CS_FLAG_COMP_EF      = 1 << 12,  /* = 4096 */

  /* Compute elemental portion of dual faces associated to the couple
     (edge, face) */

  CS_FLAG_COMP_SEF     = 1 << 13,  /* = 8192 */

  /* Compute cellwise quantities related to the height of the pyramid with
     basis spanned by a face and with apex the cell center */

  CS_FLAG_COMP_HFQ     = 1 << 14,  /* = 16384 */

  /* Compute cellwise orientation of oriented edges belonging to a face */

  CS_FLAG_COMP_FES     = 1 << 15,  /* = 32768 */

  /* Compute cellwise quantities related to the volume of the pyramid with
     basis spanned by a face and with apex the cell center */

  CS_FLAG_COMP_PFC     = 1 << 16,  /* = 65536 */

  /* Compute cellwise quantities related to the volume surrounding an edge */

  CS_FLAG_COMP_PEC     = 1 << 17,  /* = 131072 */

  /* Compute cellwise diameters */

  CS_FLAG_COMP_DIAM    = 1 << 18,  /* = 262144 */

} cs_flag_comp_bits_t;

typedef enum {

  CS_FLAG_LOCATION_PRIMAL_VTX = CS_FLAG_PRIMAL | CS_FLAG_VERTEX,
  CS_FLAG_LOCATION_PRIMAL_EDGE = CS_FLAG_PRIMAL | CS_FLAG_EDGE,
  CS_FLAG_LOCATION_PRIMAL_FACE = CS_FLAG_PRIMAL | CS_FLAG_FACE,
  CS_FLAG_LOCATION_PRIMAL_CELL = CS_FLAG_PRIMAL | CS_FLAG_CELL,
  CS_FLAG_LOCATION_DUAL_VTX  = CS_FLAG_DUAL | CS_FLAG_VERTEX,
  CS_FLAG_LOCATION_DUAL_EDGE = CS_FLAG_DUAL | CS_FLAG_EDGE,
  CS_FLAG_LOCATION_DUAL_FACE = CS_FLAG_DUAL | CS_FLAG_FACE,
  CS_FLAG_LOCATION_DUAL_CELL = CS_FLAG_DUAL | CS_FLAG_CELL,

  CS_FLAG_N_LOCATIONS

} cs_flag_location_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a two flag share the same pattern
 *         Return true if the flag to check has at least the pattern of the
 *         reference flag.
 *
 * \param[in]  flag_to_check   flag corresponding to the location to check
 * \param[in]  reference       flag corresponding to the referenced support
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_flag_test(cs_flag_t    flag_to_check,
             cs_flag_t    reference)
{
  if ((flag_to_check & reference) == reference)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a flag has in common at least one of the given masks
 *         Return true if the test is satisfied.
 *
 * \param[in]  flag_to_check   flag corresponding to the location to check
 * \param[in]  n_masks         number of masks to check
 * \param[in]  masks           array of masks
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_flag_at_least(cs_flag_t    flag_to_check,
                 int          n_masks,
                 cs_flag_t    masks[])
{
  for (int i = 0; i < n_masks; i++)
    if ((flag_to_check & masks[i]) == masks[i])
      return true;
  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a two compute flag share the same pattern
 *         Return true if the computed flag to check has at least the pattern
 *         of the reference compute flag.
 *         Case of extended flags.
 *
 * \param[in]  flag_to_check   flag corresponding to the location to check
 * \param[in]  reference       flag corresponding to the referenced support
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_eflag_test(cs_eflag_t    flag_to_check,
              cs_eflag_t    reference)
{
  if ((flag_to_check & reference) == reference)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the label associated to a location flag
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_flag_str_location(cs_flag_t  loc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FLAG_H__ */
