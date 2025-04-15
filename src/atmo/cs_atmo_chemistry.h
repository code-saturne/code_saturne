#ifndef __CS_ATMO_CHEMISTRY_H__
#define __CS_ATMO_CHEMISTRY_H__

/*============================================================================
 * Main for atmospheric chemistry related functions
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Atmospheric aerosol external library
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_ATMO_AEROSOL_OFF = 0,
  CS_ATMO_AEROSOL_SSH = 1

} cs_atmo_aerosol_type_t;

/*----------------------------------------------------------------------------
 * Atmospheric chemistry options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  /*! Choice of chemistry resolution scheme
    - 0 --> no atmospheric chemistry
    - 1 --> quasi steady equilibrium NOx scheme with 4 species and 5 reactions
    - 2 --> scheme with 20 species and 34 reactions
    - 3 --> CB05 scheme with 52 species and 155 reactions
    - 4 --> user defined schema from SPACK */

  int model;
  int n_species;
  int n_reactions;

  /*! split (=1) or semi-coupled (=2, pu-sun) resolution of chemistry */
  int chemistry_sep_mode;

  /* Flag to deactivate photolysis */
  bool chemistry_with_photolysis;

  /*! Choice of the aerosol model
       - CS_ATMO_AEROSOL_OFF ---> no aerosol model
       - CS_ATMO_AEROSOL_SSH ---> external library SSH-aerosol */
  cs_atmo_aerosol_type_t aerosol_model;

  /*! Flag to deactivate gaseous chemistry when using aerosols */
  bool frozen_gas_chem;
  /*! Flag to initialize gaseous species with the aerosol library */
  bool init_gas_with_lib;
  /*! Flag to initialize aerosols species with the aerosol library */
  bool init_aero_with_lib;
  /*! Number of layers within each aerosol */
  int n_layer;
  /*! Number of aerosols */
  int n_size;
  char *spack_file_name;
  int *species_to_scalar_id; // used only in Fortran
  int *species_to_field_id;
  int *species_profiles_to_field_id;
  /*! Molar mass of the chemical species (g/mol) */
  cs_real_t *molar_mass;
  int *chempoint;
  /*! conversion factors for reaction rates Jacobian matrix */
  cs_real_t *conv_factor_jac;
  /*! kinetics constants */
  cs_real_t *reacnum;
  /*! Initial gaseous and particulate concentrations
    and aerosol number read in file */
  cs_real_t *dlconc0;
  /*! Name of the file used to initialize the aerosol shared library */
  char *aero_file_name;
  /*! Name of the file used to initialize and to apply boundary
   *  conditions on chemical species */
  char *chem_conc_file_name;
  /*! Name of the file used to initialize and to apply boundary
   *  conditions on aerosol species */
  char *aero_conc_file_name;

  // option for chemestry profiles file

  /*! number of time steps for the concentration profiles file */
  int nt_step_profiles;
  /*! number of altitudes for the concentration profiles file */
  int n_z_profiles;
  /*! number of initialized chemical species
      in the concentration profiles file */
  int n_species_profiles;

  /*! concentration profiles */
  cs_real_t *conc_profiles;
  /*! altitudes of the concentration profiles*/
  cs_real_t *z_conc_profiles;
  /*! time steps of the concentration profiles */
  cs_real_t *t_conc_profiles;
  /*! X coordinates of concentration profiles */
  cs_real_t *x_conc_profiles;
  /*! Y coordinates of concentration profiles */
  cs_real_t *y_conc_profiles;

} cs_atmo_chemistry_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to atmo chemistry structure */
extern cs_atmo_chemistry_t *cs_glob_atmo_chemistry;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 *
 * \param[out]  array  gas concentrations
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t  *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute aerosol cloud droplets nucleation when using the atmospheric
 * humid model using a microphysical model.
 *
 * It is taken into account as an additional step split from advection-diffusion
 * equation, hence the droplet number is first clipped if necessary.
 *
 * \param[out]  nc      droplet number (scalar) in 1/cm**3
 * \param[in]   rom     density of air in kg/m**3
 * \param[in]   qldia   mass fraction of liquid water in kg/kg
 * \param[in]   pphy    true pressure in pascals
 * \param[in]   refrad  radiative cooling
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_nuclea(cs_real_t         *nc,
                       const cs_real_t   *rom,
                       const cs_real_t   *qldia,
                       const cs_real_t   *pphy,
                       const cs_real_t   *refrad);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate chemistry initialization procedure
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_initialization_deactivate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize gaseous and particulate concentrations and aerosol number
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chem_initialize_dlconc0(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric chemistry options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric aerosols options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the chemistry module needs initialization
 *
 * \return int value : 1 if needed, 0 if not
 */
/*----------------------------------------------------------------------------*/

int
cs_atmo_chemistry_need_initialization(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function sets the file name to initialize the aerosol library.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_aerosol_file_name(const char *file_name);

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
 * \brief This function declare additional transported variables for
 *        atmospheric module  for the chemistry defined from SPACK.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_declare_chem_from_spack(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize chemistry array.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_init_chemistry(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reads initial aerosol concentration and number
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_read_aerosol(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reads the chemistry profile data for the atmospheric chemistry
 *
 * \param[in] mode    if false reading of dimensions only else reading of data
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_read_chemistry_profile(int mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the aerosol concentration file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_aero_conc_file_name(const char *file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the chemistry concentration file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_chem_conc_file_name(const char *file_name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_CHEMISTRY_H__ */
