#ifndef __CS_CDO_H__
#define __CS_CDO_H__

/*============================================================================
 * General functions or variables for the INNOV module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Flag related to the activation (or not) of the CDO schemes */
#define CS_CDO_OFF      -1 // CDO schemes are not used (no activation)
#define CS_CDO_WITH_FV   0 // CDO schemes are used as well as finite volume
#define CS_CDO_ONLY      1 // CDO schemes are exclusively used

/* Flag associated to each cell */
#define CS_FLAG_BOUNDARY     (1 << 0)  //  1: cell with at least one border face

/* Flag related to the way a CDO system is built */
#define CS_FLAG_SYS_HLOC_CONF    (1 << 0) //   1: build conforming Hodge op.
#define CS_FLAG_SYS_SYM          (1 << 1) //   2: system matrix is symmetric
#define CS_FLAG_SYS_TIME_DIAG    (1 << 2) //   4: lumping/diag by construction
#define CS_FLAG_SYS_SOURCES_HLOC (1 << 3) //   8: source terms need a hodge op.
#define CS_FLAG_SYS_VECTOR       (1 << 4) //  16: vector-valued variable
#define CS_FLAG_SYS_DEBUG        (1 << 5) //  32: activate debug mode

/* Flags use to identify the nature/status of an object (variable, property) */
#define CS_FLAG_STATE_UNIFORM     (1 << 0) //    1: uniform (in space)
#define CS_FLAG_STATE_CELLWISE    (1 << 1) //    2: cellwise uniform
#define CS_FLAG_STATE_FACEWISE    (1 << 2) //    4: uniform on each face
#define CS_FLAG_STATE_STEADY      (1 << 3) //    8: steady
#define CS_FLAG_STATE_POTENTIAL   (1 << 4) //   16: potential
#define CS_FLAG_STATE_CIRCULATION (1 << 5) //   32: circulation
#define CS_FLAG_STATE_FLUX        (1 << 6) //   64: flux
#define CS_FLAG_STATE_DENSITY     (1 << 7) //  128: density
#define CS_FLAG_STATE_OWNER       (1 << 8) //  256: owner

/* Flags use to identify where is located a variable and how to access to
   its values */
#define CS_FLAG_PRIMAL    (1 <<  0) //    1: on primal mesh
#define CS_FLAG_DUAL      (1 <<  1) //    2: on dual mesh
#define CS_FLAG_VERTEX    (1 <<  2) //    4: on vertices
#define CS_FLAG_EDGE      (1 <<  3) //    8: on edges
#define CS_FLAG_FACE      (1 <<  4) //   16: on faces
#define CS_FLAG_CELL      (1 <<  5) //   32: on cells
#define CS_FLAG_BORDER    (1 <<  6) //   64: located on the boundary
#define CS_FLAG_SCALAR    (1 <<  7) //  128: scalar-valued (stride = 1)
#define CS_FLAG_VECTOR    (1 <<  8) //  256: vector-valued (stride = 3)
#define CS_FLAG_TENSOR    (1 <<  9) //  512: tensor-valued (stride = 9)
#define CS_FLAG_BY_CELL   (1 << 10) // 1024: by cell (c2e, c2f, c2v)
#define CS_FLAG_FULL_LOC  (1 << 11) // 2048: defined on the whole location

/* Flags use to identify the type of numerical schemes requested for computing
   the different equations attached to the computational domain. If flag is
   activated, then at least one equation solved is discretized using thiq kind
   of numerical scheme. */
#define CS_SCHEME_FLAG_SCALAR   (1 << 0) //  1: scheme for scalar eq.
#define CS_SCHEME_FLAG_VECTOR   (1 << 1) //  2: scheme for a vector eq.
#define CS_SCHEME_FLAG_POLY0    (1 << 2) //  4: lowest-order scheme
#define CS_SCHEME_FLAG_POLY1    (1 << 3) //  8: approx. with linear polynomials
#define CS_SCHEME_FLAG_POLY2    (1 << 4) // 16: approx. with quadratic poly.

/* Size of the buffer used to collect global ids for rows and columns
   when assembling the values in the global matrix from the local cellwise
   matrices */
#define CS_CDO_ASSEMBLE_BUF_SIZE  99

/* The following limitation only results from an optimization in the size of
   the bit mask (can be changed if needed by changing the definition of
   the type cs_mask_t) */
#define CS_CDO_N_MAX_REACTIONS  8 // Max number of reaction terms in an equation

/* Specifications for open mp loops */
#define CS_CDO_OMP_CHUNK_SIZE  128
#define CS_CDO_OMP_SCHEDULE  schedule(static, CS_CDO_OMP_CHUNK_SIZE)

/* Number of DoFs on faces and cells according to the polynomial space */
#define CS_N_FACE_DOFS_0TH  1
#define CS_N_FACE_DOFS_1ST  3
#define CS_N_FACE_DOFS_2ND  6

#define CS_N_CELL_DOFS_0TH  1
#define CS_N_CELL_DOFS_1ST  4
#define CS_N_CELL_DOFS_2ND  10

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef  unsigned char  cs_mask_t;

/* Vector-valued quantity stored using its measure (i.e. length) and
   its direction given by a unitary vector */
typedef struct {

  double  meas;
  double  unitv[3];

} cs_nvec3_t;

/* Type of numerical scheme for the discretization in space */
typedef enum {

  CS_SPACE_SCHEME_CDOVB,  /* CDO scheme with vertex-based positionning */
  CS_SPACE_SCHEME_CDOVCB, /* CDO scheme with vertex+cell-based positionning */
  CS_SPACE_SCHEME_CDOFB,  /* CDO cell-based scheme with hybridization */
  CS_SPACE_SCHEME_HHO_P0, /* Hybrid High Order scheme; P0 approx. of gradient */
  CS_SPACE_SCHEME_HHO_P1, /* Hybrid High Order scheme; P1 approx. of gradient */
  CS_SPACE_SCHEME_HHO_P2, /* Hybrid High Order scheme; P2 approx. of gradient */

  CS_SPACE_N_SCHEMES

} cs_space_scheme_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for an analytic function
 *         elt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *
 * \param[in]      time     when ?
 * \param[in]      n_elts   number of elements to consider
 * \param[in]      elt_ids  list of elements ids (to access coords and fill)
 * \param[in]      coords   where ?
 * \param[in]      compact  true:no indirection, false:indirection for filling
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] retval   result of the function
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_analytic_func_t) (cs_real_t            time,
                      cs_lnum_t            n_elts,
                      const cs_lnum_t     *elt_ids,
                      const cs_real_t     *coords,
                      bool                 compact,
                      void                *input,
                      cs_real_t           *retval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function which defines the time step according to the number of
 *         iteration already done, the current time and any structure given as
 *         a parameter
 *
 * \param[in]   time_iter   current number of iterations
 * \param[in]   time        value of the time at the end of the last iteration
 * \param[in]   input       pointer to a structure cast on-the-fly
 *
 * \return the value of the time step
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_timestep_func_t) (int       time_iter,
                      double    time,
                      void     *input);

/*============================================================================
 * Global variables
 *============================================================================*/

/* Activation of the CDO/HHO module */
extern int  cs_cdo_activation_mode;

/* Separation lines: long, medium, short */
extern const char lsepline[80];
extern const char msepline[60];
extern const char ssepline[40];

/* Default locations */
extern const cs_flag_t  cs_cdo_primal_vtx;
extern const cs_flag_t  cs_cdo_primal_face;
extern const cs_flag_t  cs_cdo_primal_cell;
extern const cs_flag_t  cs_cdo_dual_vtx;
extern const cs_flag_t  cs_cdo_dual_face;
extern const cs_flag_t  cs_cdo_dual_cell;
extern const cs_flag_t  cs_cdo_dual_face_byc;

/*============================================================================
 * Static inline function prototypes
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
cs_test_flag(cs_flag_t    flag_to_check,
             cs_flag_t    reference)
{
  if ((flag_to_check & reference) == reference)
    return true;
  else
    return false;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string "true" or "false" according to the boolean
 *
 * \param[in]  boolean  bool type
 *
 * \return a string "true" or "false"
 */
/*----------------------------------------------------------------------------*/

static inline const char *
cs_base_strtf(bool  boolean)
{
  if (boolean)
    return "true";
  else
    return "false";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_nvec3_t structure from a cs_real_3_t
 *
 * \param[in]  v     vector of size 3
 * \param[out] qv    pointer to a cs_nvec3_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_nvec3(const cs_real_3_t    v,
         cs_nvec3_t          *qv)
{
  cs_real_t  magnitude = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  qv->meas = magnitude;
  if (fabs(magnitude) > cs_math_zero_threshold) {

    const cs_real_t  inv = 1/magnitude;
    qv->unitv[0] = inv * v[0];
    qv->unitv[1] = inv * v[1];
    qv->unitv[2] = inv * v[2];

  }
  else
    qv->unitv[0] = qv->unitv[1] = qv->unitv[2] = 0;

}

#if defined(DEBUG) && !defined(NDEBUG)
/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a cs_sdm_t structure which is defined by block
 *          Print into the file f if given otherwise open a new file named
 *          fname if given otherwise print into the standard output
 *          The usage of threshold allows one to compare more easier matrices
 *          without taking into account numerical roundoff.
 *
 * \param[in]  fp         pointer to a file structure or NULL
 * \param[in]  fname      filename or NULL
 * \param[in]  thd        threshold (below this value --> set 0)
 * \param[in]  n_elts     size of the array
 * \param[in]  array      list of values to dump
 * \param[in]  n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_array_fprintf(FILE             *fp,
                     const char       *fname,
                     cs_real_t         thd,
                     cs_lnum_t         n_elts,
                     const cs_real_t   array[],
                     int               n_cols);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In debug mode, dump an array of double into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_darray_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_real_t    array[],
                         int                n_cols);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  In debug mode, dump an array of integer into the listing
 *
 * \param[in] header     header message to write
 * \param[in] size       number of elements in array
 * \param[in] array      pointer to the array of values
 * \param[in] n_cols     print array with n_cols columns
 */
/*----------------------------------------------------------------------------*/

void
cs_dbg_iarray_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_lnum_t    array[],
                         int                n_cols);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_H__ */
