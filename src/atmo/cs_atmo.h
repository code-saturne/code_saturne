#ifndef __CS_ATMO_H__
#define __CS_ATMO_H__

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
 *  Local headers
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
 * Atmospheric model options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {
  /* Space and tim reference of the run */
  /*! Starting year */
  int syear;
  /*! Starting quantile */
  int squant;
  /*! Starting hour */
  int shour;
  /*! Starting minute */
  int smin;
  /*! Starting second */
  cs_real_t ssec;
  /*! longitude of the domain origin */
  cs_real_t longitute;
  /*! latitude of the domain origin */
  cs_real_t latitude;

  /* Model options */
  bool compute_z_ground;
} cs_atmo_option_t;

/*----------------------------------------------------------------------------
 * Atmospheric chemistry options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {
  /*! Choice of chemistry resolution scheme
       - 0 --> no atmospheric chemistry
       - 1 --> quasi steady equilibrium NOx scheme with 4 species and 5 reactions
       - 2 --> scheme with 20 species and 34 reactions
       - 3 --> scheme CB05 with 52 species and 155 reactions
       - 4 --> user defined schema from SPACK */
  int model;
  int n_species;
  int n_reactions;
  char *spack_file_name;
  int *species_to_scalar_id; // used only in fortran
  int *species_to_field_id;
  /*! Molar mass of the chemical species (g/mol) */
  cs_real_t *molar_mass;
  int *chempoint;
} cs_atmo_chemistry_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to atmo options structure */
extern cs_atmo_option_t        *cs_glob_atmo_option;

/* Pointer to atmo chemistry structure */
extern cs_atmo_chemistry_t *cs_glob_atmo_chemistry;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the SPACK file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_spack_file_name(const char *file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes the ground elevation
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{g} \right)
 *      - \divs \left( \vect{V} \right) \varia = 0
 *  \f]
 *  where \f$ \vect{g} \f$ is the gravity field
 *
 *  The boundary conditions on \f$ \varia \f$ read:
 *  \f[
 *   \varia = z \textrm{ on walls}
 *  \f]
 *  \f[
 *   \dfrac{\partial \varia}{\partial n} = 0 \textrm{ elsewhere}
 *  \f]
 *
 *  Remarks:
 *  - a steady state is looked for.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_z_ground_compute(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function declare additional transported variables for
 *        atmospheric module  for the chemistry defined from SPACK.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_declare_chem_from_spack(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_H__ */
