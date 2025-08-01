#ifndef __CS_POROSITY_FROM_SCAN_H__
#define __CS_POROSITY_FROM_SCAN_H__

/*============================================================================
 * Main for cooling towers related functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm/fvm_nodal.h"

#include "base/cs_base.h"
#include "base/cs_halo.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------
 * Porosity from scan model options descriptor
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_COG_FROM_FLUID_FACES = 0,
  CS_COG_FROM_PYRAMID = 1,
  CS_COG_WITHOUT_RECONSTRUCTION_FOR_IBM_PLANE = 2

} cs_ibm_cog_location_t;

typedef enum {

  CS_FILL_RADIAL = 0,
  CS_FILL_DIRECTION = 1,

} cs_fill_type_t;

typedef struct {
  bool  compute_porosity_from_scan;
  char *file_names;
  int n_headers;
  int *header_type;
  char **headers;
  char *output_name;
  bool  postprocess_points;
  /*! Matrix of associated transformation
     (3x4 matrix, 3 first rows of a homogeneous
     coordinates transformation matrix,
     with last row = [0 0 0 1]) */
  cs_real_34_t transformation_matrix;
  /*! To define the direction of fill,
      by default it is filled from top (z-direction).*/
  cs_real_3_t direction_vector;
  cs_fill_type_t type_fill;
  int   nb_sources;
  cs_real_3_t *sources;
  cs_lnum_t *source_c_ids;
  cs_lnum_t threshold;
  cs_lnum_t n_agglomeration;
  cs_real_t porosity_threshold;
  cs_real_t convection_porosity_threshold;
  bool      use_staircase;
  cs_real_t eigenvalue_criteria;
  int       use_restart;
  cs_ibm_cog_location_t cog_location;
  cs_real_33_t *mom_mat;
  bool  has_classification;
  float *classification_values;
  bool  *class_used;
  int   n_classifications;
} cs_porosity_from_scan_opt_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to options structure */
extern cs_porosity_from_scan_opt_t *cs_glob_porosity_from_scan_opt;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the file name of points for the computation of the
 *        porosity from scan.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_file_name(const char  *file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the output name for the FVM writer of scan points.
 *
 * \param[in] output_name  name of the output (a suffix will be added)
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_output_name(const char  *output_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a scanner source point.
 *
 * \param[in] source     source vector
 * \param[in] transform  flag to apply the transformation matrix to the source
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_add_source(const cs_real_t  source[3],
                                 bool             transform);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the scanner sources from csv file to fill fluid space.
 *
 * \param[in] csv file containing the (x,y,z) coordinates of each scanner
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_sources_by_file_name(const char *file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the porosity which is equal to one from
 *        a source, radiating sphericaly, and is 0 when touching points
 *        of the scan.
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{e}_r \right)
 *      - \divs \left( \vect{e}_r \right) \varia = 0
 *  \f]
 *  where \f$ \vect{e}_r = \dfrac{\vect{x} - \vect{x}_0}{\norm{\vect{x} - \vect{x}_0}} \f$
 *  is the radial direction from the source \f$\vect{x}_0 \f$.
 *
 *  The boundary conditions on \f$ \varia \f$ is an homogeneous Neumann, and
 *  a penalisation term is impose in the cell of center \f$ \vect{x}_0\f$.
 *  \f[
 *   \dfrac{\partial \varia}{\partial n} = 0 \textrm{everywhere}
 *  \f]
 *
 *  Remarks:
 *  - a steady state is looked for.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_porosity_from_scan(void);

/*----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*
 * \brief Write the restart file of the ibm module
 */
/*--------------------------------------------------------------------*/

void
cs_porous_model_restart_write(void);

/*--------------------------------------------------------------------*/
/*
 * \brief Read the restart file of the ibm module
 */
/*--------------------------------------------------------------------*/

void
cs_porous_model_restart_read(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POROSITY_FROM_SCAN_H__ */
