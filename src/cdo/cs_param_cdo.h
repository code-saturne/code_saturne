#ifndef __CS_PARAM_CDO_H__
#define __CS_PARAM_CDO_H__

/*============================================================================
 * Manage the definition/setting of a computation
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

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

#define CS_ALL_FACES   0        /* All faces: interior + border */
#define CS_BND_FACES   1        /* Boundary faces */
#define CS_INT_FACES   2        /* Interior faces */

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

/* DISCRETE HODGE OPERATORS */
/* ======================== */

typedef enum {

  CS_PARAM_HODGE_TYPE_VPCD, // from primal vertices to dual cells
  CS_PARAM_HODGE_TYPE_EPFD, // from primal edges to dual faces
  CS_PARAM_HODGE_TYPE_FPED, // from primal faces to dual edges
  CS_PARAM_HODGE_TYPE_EDFP, // from dual edges to primal faces
  CS_PARAM_HODGE_TYPE_CPVD, // from primal cells to dual vertices
  CS_PARAM_HODGE_TYPE_VC,   // primal vertices + primal cells
  CS_PARAM_N_HODGE_TYPES

} cs_param_hodge_type_t;

typedef enum {

  CS_PARAM_HODGE_ALGO_VORONOI, // Under othogonality condition gives a diag. op.
  CS_PARAM_HODGE_ALGO_WBS,     // WBS: Whitney Barycentric Subdivision
  CS_PARAM_HODGE_ALGO_COST,    // COST: COnsistency & STabilization splitting
  CS_PARAM_HODGE_ALGO_AUTO,    /* Switch between the previous algo. according to
                                  the type of cell and its property */
  CS_PARAM_N_HODGE_ALGOS

} cs_param_hodge_algo_t;

typedef struct {

  bool   is_unity; /* No associated property. Property is equalt to the unity */
  bool   is_iso;   /* Associated property is isotroopic ? */
  bool   inv_pty;  /* Definition based on the property or its inverse */

  cs_param_hodge_type_t   type;  /* type of discrete Hodge operator */
  cs_param_hodge_algo_t   algo;  /* type of algorithm used to build this op. */
  double                  coef;  /* Value of the stabilization parameter
                                    if the COST algo. is used, otherwise 0. */

} cs_param_hodge_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: long, medium, short */
extern const char lsepline[80];
extern const char msepline[60];
extern const char ssepline[40];

/* Activation of the CDO/HHO module */
extern int  cs_param_cdo_mode;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of algorithm related to a discrete Hdoge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the algorithm
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_algo_name(const cs_param_hodge_t   h_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the type of discrete Hodge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_type_name(const cs_param_hodge_t   h_info);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_CDO_H__ */
