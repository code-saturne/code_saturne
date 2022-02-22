#ifndef __CS_BOUNDARY_H__
#define __CS_BOUNDARY_H__

/*============================================================================
 * Handle the boundaries of a computational domain
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Rescaling or not of the pressure */
#define CS_BOUNDARY_PRESSURE_NO_RESCALING   1
#define CS_BOUNDARY_PRESSURE_RESCALING      0

/* Name of the boundary zone gathering all domain boundary walls */
#define CS_BOUNDARY_WALLS_NAME   "auto:wall"

#define CS_BOUNDARY_UNDEFINED   0

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Boundary categories */

typedef enum {

  CS_BOUNDARY_CATEGORY_FLOW,       /*< flow related boundaries */
  CS_BOUNDARY_CATEGORY_ALE,        /*< ALE related boundaries */
  CS_BOUNDARY_CATEGORY_RADIATIVE   /*< Radiative boundaries */

} cs_boundary_category_t;

/*! Flag values defining boundary condition subtypes (0 for none) */

typedef int  cs_boundary_type_t;

/* Bit values for flow boundaries
   ------------------------------ */

typedef enum {

  /* Main types
     ---------- */

  /*! wall */
  CS_BOUNDARY_WALL                = 1<<0,

  /*! inlet */
  CS_BOUNDARY_INLET               = 1<<1,

  /*! outlet */
  CS_BOUNDARY_OUTLET              = 1<<2,

  /*! symmetry */
  CS_BOUNDARY_SYMMETRY            = 1<<3,

  /* Additional flags
     ---------------- */

  /*! rough wall */
  CS_BOUNDARY_ROUGH_WALL          = 1<<4,

  /*! sliding wall */
  CS_BOUNDARY_SLIDING_WALL        = 1<<5,

  /*! imposed velocity */
  CS_BOUNDARY_IMPOSED_VEL         = 1<<6,

  /*! imposed pressure*/
  CS_BOUNDARY_IMPOSED_P           = 1<<7,

  /*! free inlet-outlet */
  CS_BOUNDARY_FREE_INLET_OUTLET   = 1<<8,

  /*! convective inlet */
  CS_BOUNDARY_CONVECTIVE_INLET    = 1<<9,

  /*! compressible inlet, imposed flux and enthalpy */
  CS_BOUNDARY_INLET_QH            = 1<<10,

  /*! compressible (subsonic) inlet, imposed pressure and enthalpy */
  CS_BOUNDARY_INLET_SUBSONIC_PH   = 1<<11,

  /*! compressible subsonic */
  CS_BOUNDARY_SUBSONIC             = 1<<12,

  /*! compressible supersonic */
  CS_BOUNDARY_SUPERSONIC           = 1<<13,

  /*! free surface */
  CS_BOUNDARY_FREE_SURFACE        = 1<<14,

  /*! coupled */
  CS_BOUNDARY_COUPLED             = 1<<15,

  /*! coupled with decentered flux */
  CS_BOUNDARY_COUPLED_DF          = 1<<16

} cs_boundary_flow_subtype_bits_t;

/* Bit values for ALE boundaries
   ----------------------------- */

typedef enum {

  CS_BOUNDARY_ALE_FIXED               = 1<<0, /*!< fixed */
  CS_BOUNDARY_ALE_SLIDING             = 1<<1, /*!< sliding */
  CS_BOUNDARY_ALE_IMPOSED_VEL         = 1<<2, /*!< imposed velocity */
  CS_BOUNDARY_ALE_IMPOSED_DISP        = 1<<3, /*!< imposed displacement */
  CS_BOUNDARY_ALE_INTERNAL_COUPLING   = 1<<4, /*!< internal coupling */
  CS_BOUNDARY_ALE_EXTERNAL_COUPLING   = 1<<5, /*!< external coupling */
  CS_BOUNDARY_ALE_FREE_SURFACE        = 1<<6  /*!< free surface */

} cs_boundary_ale_subtype_bits_t;

/*! \struct cs_boundary_t
 *  \brief Structure storing information related to the "physical"
 *  boundaries associated with the computational domain
 */

typedef struct {

  cs_boundary_category_t  category;      /*!< boundary category */
  cs_boundary_type_t      default_type;  /*!< default boundary */

  int                     n_boundaries;  /*!< number of boundaries */
  cs_boundary_type_t     *types;         /*!< type of each boundary */
  int                    *zone_ids;      /*!< associated zone ids */

} cs_boundary_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_boundary_t  *cs_glob_boundaries; /* Pointer to the shared boundaries
                                            * on the computational domain */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a boundary with a given flag is present.
 *
 * \param[in]  boundaries   pointer to a cs_boundary_t structure
 * \param[in]  type_flag    boundary type flag
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_boundary_has_type(const cs_boundary_t  *boundaries,
                     int                   type_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the related id associated to a boundary from its zone id
 *
 * \param[in] boundaries       pointer to a cs_boundary_t structure
 * \param[in] z_id             id of the related zone
 *
 * \return the associated boundary id in the boundary list
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_id_by_zone_id(const cs_boundary_t  *boundaries,
                          int                   z_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to the given \ref cs_boundary_t
 *         structure
 *
 * \param[in, out]   boundaries   pointer to a structure storing boundary info
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_set_default(cs_boundary_t        *boundaries,
                        cs_boundary_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a default boundary structure for the computational domain
 *
 * \param[in]  category       default type of boundary to set
 * \param[in]  default_type   default type of boundary to set
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_boundary_t *
cs_boundary_create(cs_boundary_category_t  category,
                   cs_boundary_type_t      default_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all metadate related to the domain boundaries
 *
 * \param[in, out]   p_boundaries   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_free(cs_boundary_t   **p_boundaries);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new boundary type for a given boundary zone
 *
 * \param[in, out] bdy          pointer to a structure storing boundary info
 * \param[in]      type         type of boundary to set
 * \param[in]      zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_add(cs_boundary_t        *bdy,
                cs_boundary_type_t    type,
                const char           *zone_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build an array on boundary faces which specifies the boundary type
 *         for each face.
 *
 * \param[in]       boundaries    pointer to the domain boundaries
 * \param[in]       n_b_faces     number of boundaries faces
 * \param[in, out]  bf_type       boundary type flag
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_build_type_array(const cs_boundary_t   *boundaries,
                             cs_lnum_t              n_b_faces,
                             cs_boundary_type_t     bf_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new zone gathering all CS_BOUNDARY_WALL type zones
 *
 * \param[in, out]  boundaries    pointer to the domain boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_def_wall_zones(cs_boundary_t   *boundaries);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to perform a pressure rescaling (in case of a
 *         Dirichlet on the velocity for the whole boundary)
 *         Use in CDO schemes for Navier--Stokes
 *
 * \param[in] n_b_faces    number of border faces
 * \param[in] bf_type      array of types of boundary for each boundary face
 *
 * \return 1 if a pressure rescaling is needed otherwise 0
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_need_pressure_rescaling(cs_lnum_t                  n_b_faces,
                                    const cs_boundary_type_t   bf_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a boundary type description
 *
 * \param[in]   bdy            pointer to a structure storing boundary info
 * \param[in]   b_type         type flag
 * \param[in]   descr_len_max  maximum name length
 * \param[out]  descr          subtype name
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_get_type_descr(const cs_boundary_t  *bdy,
                           cs_boundary_type_t    b_type,
                           int                   descr_len_max,
                           char                  descr[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the setup of the boundary of the computational domain
 *
 * \param[in] bdy          pointer to a structure storing boundary info
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_log_setup(const cs_boundary_t     *bdy);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_H__ */
