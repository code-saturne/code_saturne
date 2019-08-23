#ifndef __CS_POROSITY_FROM_SCAN_H__
#define __CS_POROSITY_FROM_SCAN_H__

/*============================================================================
 * Main for cooling towers related functions
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

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

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

typedef struct {
  bool  compute_porosity_from_scan;
  char *file_name;
  bool  postprocess_points;
  /*! Matrix of associated transformation
     (3x4 matrix, 3 first rows of a homogeneous
     coordinates transformation matrix,
     with last row = [0 0 0 1]) */
  cs_real_34_t transformation_matrix;
  int   nb_sources;
  cs_real_3_t *sources;
  cs_lnum_t *source_c_ids;
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
 * \brief This function set the file name of points for the computation of the
 * porosity from scan.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_set_file_name(const char *file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function add a scanner source point
 *
 * \param[in] source     source vector
 * \param[in] transform  flag to apply the transformation matrix to the source
 */
/*----------------------------------------------------------------------------*/

void
cs_porosity_from_scan_add_source(const cs_real_3_t source,
                                 const bool transform);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes the porosity which is equal to one from
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

END_C_DECLS

#endif /* __CS_POROSITY_FROM_SCAN_H__ */
