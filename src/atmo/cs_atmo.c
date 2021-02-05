/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_atmo_profile_std.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* global atmo options structure */
static cs_atmo_option_t  _atmo_option = {
  .syear = -1,
  .squant = -1,
  .shour = -1,
  .smin = -1,
  .ssec = -1.,
  .longitude = 1e12, // TODO use cs_math_big_r
  .latitude = 1e12,
  .domain_orientation = 0.,
  .compute_z_ground = false,
  .open_bcs_treatment = 0,
  .sedimentation_model = 0,
  .deposition_model = 0,
  .nucleation_model = 0,
  .subgrid_model = 0,
  .distribution_model = 1, /* all or nothing */
  .meteo_profile = 0, /* no meteo profile */
  .meteo_file_name = NULL,
  .meteo_dlmo = 0.,
  .meteo_z0 = -1.,
  .meteo_zref = -1.,
  .meteo_zi = -1.,
  .meteo_zu1 = -1.,
  .meteo_zu2 = -1.,
  .meteo_zt1 = -1.,
  .meteo_zt2 = -1.,
  .meteo_uref = -1.,
  .meteo_u1 = -1.,
  .meteo_u2 = -1.,
  .meteo_ustar0 = -1.,
  .meteo_wstar0 = -1.,
  .meteo_angle = -1.,
  .meteo_t0 = 284.15, /* 11 deg C */
  .meteo_t1 = 0.,
  .meteo_t2 = 0.,
  .meteo_tstar = 0.,
  .meteo_psea = 101325.,
  .nbmetd = 0,
  .nbmett = 0,
  .nbmetm = 0,
  .nbmaxt = 0,
  .z_temp_met = NULL,
  .time_met   = NULL,
  .hyd_p_met  = NULL
};

/* global atmo constants structure */
static cs_atmo_constants_t _atmo_constants = {
  .ps = 1.e5
};

/* atmo chemistry options structure */
static cs_atmo_chemistry_t _atmo_chem = {
  .model = 0,
  .n_species = 0,
  .n_reactions = 0,
  .chemistry_with_photolysis = true,
  .aerosol_model = CS_ATMO_AEROSOL_OFF,
  .frozen_gas_chem = false,
  .init_gas_with_lib = false,
  .init_aero_with_lib = false,
  .n_layer = 0,
  .n_size = 0,
  .spack_file_name = NULL,
  .species_to_scalar_id = NULL,
  .species_to_field_id = NULL,
  .molar_mass = NULL,
  .chempoint = NULL,
  .aero_file_name = NULL
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_option_t *cs_glob_atmo_option = &_atmo_option;

cs_atmo_constants_t *cs_glob_atmo_constants = &_atmo_constants;

cs_atmo_chemistry_t *cs_glob_atmo_chemistry = &_atmo_chem;

static const char *cs_atmo_aerosol_type_enum_name[]
  = {"CS_ATMO_AEROSOL_OFF",
     "CS_ATMO_AEROSOL_SSH"};

static const char *cs_atmo_aerosol_type_name[]
  = {N_("No atmospheric aerosol"),
     N_("Atmospheric aerosol using external code SSH-aerosol")};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

cs_real_t
cs_mo_phim(cs_real_t              z,
           cs_real_t              dlmo);

cs_real_t
cs_mo_phih(cs_real_t              z,
           cs_real_t              dlmo);

cs_real_t
cs_mo_psim(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo);

cs_real_t
cs_mo_psih(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo);

void
cs_atmo_init_meteo_profiles(void);

void
cs_atmo_compute_meteo_profiles(void);

void
cs_f_atmo_get_meteo_file_name(int           name_max,
                              const char  **name,
                              int          *name_len);

void
cs_f_atmo_get_pointers(cs_real_t              **ps,
                       int                    **syear,
                       int                    **squant,
                       int                    **shour,
                       int                    **smin,
                       cs_real_t              **ssec,
                       cs_real_t              **longitude,
                       cs_real_t              **latitude,
                       bool                   **compute_z_ground,
                       int                    **open_bcs_treatment,
                       int                    **sedimentation_model,
                       int                    **deposition_model,
                       int                    **nucleation_model,
                       int                    **subgrid_model,
                       int                    **distribution_model,
                       int                    **model,
                       int                    **n_species,
                       int                    **n_reactions,
                       bool                   **chemistry_with_photolysis,
                       cs_atmo_aerosol_type_t **aerosol_model,
                       bool                   **frozen_gas_chem,
                       bool                   **init_gas_with_lib,
                       bool                   **init_aero_with_lib,
                       int                    **n_layer,
                       int                    **n_size,
                       int                    **meteo_profile,
                       int                    **nbmetd,
                       int                    **nbmett,
                       int                    **nbmetm,
                       int                    **nbmaxt);

void
cs_f_atmo_arrays_get_pointers(cs_real_t **z_temp_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              int         dim_hyd_p_met[2]);

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint);

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid);

void
cs_f_atmo_chem_finalize(void);

void
cs_f_atmo_finalize(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the name meteo file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_meteo_file_name(int           name_max,
                              const char  **name,
                              int          *name_len)
{
  *name = _atmo_option.meteo_file_name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving meteo file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       _atmo_option.meteo_file_name, name_max, *name, *name_len);
  }
}

/*----------------------------------------------------------------------------
 * Get pointer
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_pointers(cs_real_t              **ps,
                       int                    **syear,
                       int                    **squant,
                       int                    **shour,
                       int                    **smin,
                       cs_real_t              **ssec,
                       cs_real_t              **longitude,
                       cs_real_t              **latitude,
                       bool                   **compute_z_ground,
                       int                    **open_bcs_treatment,
                       int                    **sedimentation_model,
                       int                    **deposition_model,
                       int                    **nucleation_model,
                       int                    **subgrid_model,
                       int                    **distribution_model,
                       int                    **model,
                       int                    **n_species,
                       int                    **n_reactions,
                       bool                   **chemistry_with_photolysis,
                       cs_atmo_aerosol_type_t **aerosol_model,
                       bool                   **frozen_gas_chem,
                       bool                   **init_gas_with_lib,
                       bool                   **init_aero_with_lib,
                       int                    **n_layer,
                       int                    **n_size,
                       int                    **meteo_profile,
                       int                    **nbmetd,
                       int                    **nbmett,
                       int                    **nbmetm,
                       int                    **nbmaxt)
{
  *ps        = &(_atmo_constants.ps);
  *syear     = &(_atmo_option.syear);
  *squant    = &(_atmo_option.squant);
  *shour     = &(_atmo_option.shour);
  *smin      = &(_atmo_option.smin);
  *ssec      = &(_atmo_option.ssec);
  *longitude = &(_atmo_option.longitude);
  *latitude  = &(_atmo_option.latitude);
  *compute_z_ground = &(_atmo_option.compute_z_ground);
  *open_bcs_treatment = &(_atmo_option.open_bcs_treatment);
  *sedimentation_model = &(_atmo_option.sedimentation_model);
  *deposition_model = &(_atmo_option.deposition_model);
  *nucleation_model = &(_atmo_option.nucleation_model);
  *subgrid_model = &(_atmo_option.subgrid_model);
  *distribution_model = &(_atmo_option.distribution_model);
  *meteo_profile = &(_atmo_option.meteo_profile);
  *nbmetd     = &(_atmo_option.nbmetd);
  *nbmett     = &(_atmo_option.nbmett);
  *nbmetm     = &(_atmo_option.nbmetm);
  *nbmaxt     = &(_atmo_option.nbmaxt);
  *model = &(_atmo_chem.model);
  *n_species = &(_atmo_chem.n_species);
  *n_reactions = &(_atmo_chem.n_reactions);
  *chemistry_with_photolysis = &(_atmo_chem.chemistry_with_photolysis);
  *aerosol_model = &(_atmo_chem.aerosol_model);
  *frozen_gas_chem = &(_atmo_chem.frozen_gas_chem);
  *init_gas_with_lib = &(_atmo_chem.init_gas_with_lib);
  *init_aero_with_lib = &(_atmo_chem.init_aero_with_lib);
  *n_layer = &(_atmo_chem.n_layer);
  *n_size = &(_atmo_chem.n_size);
}

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint)
{
  if (_atmo_chem.species_to_scalar_id == NULL)
    BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);
  if (_atmo_chem.species_to_field_id == NULL)
    BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
  if (_atmo_chem.molar_mass == NULL)
    BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);
  if (_atmo_chem.chempoint == NULL)
    BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);

  *species_to_scalar_id = (_atmo_chem.species_to_scalar_id);
  *molar_mass = (_atmo_chem.molar_mass);
  *chempoint = (_atmo_chem.chempoint);
}

void
cs_f_atmo_arrays_get_pointers(cs_real_t **z_temp_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              int         dim_hyd_p_met[2])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  if (_atmo_option.z_temp_met == NULL)
    BFT_MALLOC(_atmo_option.z_temp_met, _atmo_option.nbmaxt, cs_real_t);
  if (_atmo_option.time_met == NULL)
    BFT_MALLOC(_atmo_option.time_met, _atmo_option.nbmetm, cs_real_t);
  if (_atmo_option.hyd_p_met == NULL)
    BFT_MALLOC(_atmo_option.hyd_p_met,
               _atmo_option.nbmetm*_atmo_option.nbmaxt, cs_real_t);

  *hyd_p_met       = _atmo_option.hyd_p_met;
  dim_hyd_p_met[0] = _atmo_option.nbmaxt;
  dim_hyd_p_met[1] = _atmo_option.nbmetm;

  *z_temp_met = _atmo_option.z_temp_met;
  *time_met   = _atmo_option.time_met;
}

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid)
{
  assert(species_fid != NULL);
  assert(_atmo_chem.species_to_field_id != NULL);

  for (int i = 0; i < _atmo_chem.n_species; i++)
    _atmo_chem.species_to_field_id[i] = species_fid[i];
}

void
cs_f_atmo_chem_finalize(void)
{
  if (_atmo_chem.aerosol_model != CS_ATMO_AEROSOL_OFF)
    cs_atmo_aerosol_finalize();

  BFT_FREE(_atmo_chem.species_to_scalar_id);
  BFT_FREE(_atmo_chem.species_to_field_id);
  BFT_FREE(_atmo_chem.molar_mass);
  BFT_FREE(_atmo_chem.chempoint);
  BFT_FREE(_atmo_chem.spack_file_name);
  BFT_FREE(_atmo_chem.aero_file_name);
}


void
cs_f_atmo_finalize(void)
{
  BFT_FREE(_atmo_option.meteo_file_name);
  BFT_FREE(_atmo_option.z_temp_met);
  BFT_FREE(_atmo_option.time_met);
  BFT_FREE(_atmo_option.hyd_p_met);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert string to lower case
 *----------------------------------------------------------------------------*/

static void
_strtolower(char        *dest,
            const char  *src)
{
  char *_dest = dest;
  while (*src) {
    *_dest = tolower(*src);
    src++;
    _dest++;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize meteo profiles if no meteo file is given
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_init_meteo_profiles(void)
{
  /* Some turbulence constants */
  cs_real_t kappa = cs_turb_xkappa;

  cs_atmo_option_t *aopt = &_atmo_option;
  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  /* potential temp at ref */
  cs_real_t pref = cs_glob_atmo_constants->ps;
  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t cp0 = phys_pro->cp0;
  cs_real_t rscp = rair/cp0;
  cs_real_t g = cs_math_3_norm(cs_glob_physical_constants->gravity);
  if (g <= 0.)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo meteo profiles: gravity must not be 0.\n"));

  cs_real_t theta0 = aopt->meteo_t0 * pow(pref/ aopt->meteo_psea, rscp);

  /* Reference fluid properties set from meteo values */
  phys_pro->p0 = aopt->meteo_psea;
  phys_pro->t0 = theta0; /* ref potential temp theta0*/

  cs_real_t z0 = aopt->meteo_z0;
  cs_real_t zref = aopt->meteo_zref;
  if (aopt->meteo_ustar0 < 0. && aopt->meteo_uref < 0.)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo meteo profiles: meteo_ustar0 or meteo_uref.\n"));

  /* Recompute LMO inverse */
  if (aopt->meteo_ustar0 >= 0. && aopt->meteo_uref >= 0.) {
    //TODO iterative process
  }

  /* Compute ground friction velocity from dlmo and uref */
  if (aopt->meteo_ustar0 < 0.)
    aopt->meteo_ustar0 =
      aopt->meteo_uref * kappa
      / cs_mo_psim(zref + z0,
                   z0,
                   aopt->meteo_dlmo);

  /* Compute uref from ground friction velocity and dlmo */
  if (aopt->meteo_uref < 0.)
    aopt->meteo_uref =
      aopt->meteo_ustar0 / kappa
      * cs_mo_psim(zref + z0,
                   z0,
                   aopt->meteo_dlmo);

  /* LMO inverse, ustar at ground */
  cs_real_t dlmo = aopt->meteo_dlmo;
  cs_real_t ustar0 = aopt->meteo_ustar0;

  /* Friction temperature */
  aopt->meteo_tstar = cs_math_pow2(ustar0) * theta0 * dlmo / (kappa * g);
  cs_real_t tstar = aopt->meteo_tstar;

  /* BL height according to Marht 1982 formula */
  /* value of C=0.2, 0.185, 0.06, 0.14, 0.07, 0.04 */
  cs_real_t zi_coef = 0.2;
  cs_real_t corio_f = 4. * cs_math_pi / 86164.
    * sin(aopt->latitude * cs_math_pi / 180.);
  aopt->meteo_zi = zi_coef * ustar0 / corio_f;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute meteo profiles if no meteo file is given
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_compute_meteo_profiles(void)
{
  cs_domain_t *domain = cs_glob_domain;
  cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;
  /* Get fields */
  cs_real_t *cpro_met_potemp = cs_field_by_name("meteo_pot_temperature")->val;
  cs_real_3_t *cpro_met_vel =
    (cs_real_3_t *) (cs_field_by_name("meteo_velocity")->val);
  cs_real_t *cpro_met_k = cs_field_by_name("meteo_tke")->val;
  cs_real_t *cpro_met_eps = cs_field_by_name("meteo_eps")->val;
  cs_real_t *cpro_met_p = cs_field_by_name("meteo_pressure")->val;
  cs_real_t *cpro_met_rho = cs_field_by_name("meteo_density")->val;
  cs_real_t *cpro_met_t = cs_field_by_name("meteo_temperature")->val;

  /* Some turbulence constants */
  cs_real_t kappa = cs_turb_xkappa;
  cs_real_t cmu = cs_turb_cmu;

  cs_atmo_option_t *aopt = &_atmo_option;
  cs_real_t z0 = aopt->meteo_z0;

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  /* potential temp at ref */
  cs_real_t pref = cs_glob_atmo_constants->ps;
  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t cp0 = phys_pro->cp0;
  cs_real_t rscp = rair/cp0;
  cs_real_t theta0 = aopt->meteo_t0 * pow(pref/ aopt->meteo_psea, rscp);

  /* LMO inverse, ustar at ground */
  cs_real_t dlmo = aopt->meteo_dlmo;
  cs_real_t ustar0 = aopt->meteo_ustar0;
  cs_real_t angle = aopt->meteo_angle;

  /* Friction temperature */
  cs_real_t tstar = aopt->meteo_tstar;

  /* Variables used for clipping */
  cs_real_t ri_max = cs_math_big_r;
  cs_real_t *dlmo_var = NULL;
  cs_real_t z_min = cs_math_big_r;
  cs_real_t u_met_min;
  cs_real_t theta_met_min;

  if (aopt->compute_z_ground == true){
    cs_atmo_z_ground_compute();
  }
  cs_real_t *z_ground = cs_field_by_name_try("z_ground")->val;

  BFT_MALLOC(dlmo_var, m->n_cells, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    dlmo_var[cell_id]=0.0;
  }

  if (dlmo>0) {
    ri_max = 0.75; // Value chosen to limit buoyancy vs shear production
  }

  /* Profiles */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

    //TODO reference altitude or use z_ground?
    cs_real_t z = cell_cen[cell_id][2] - z_ground[cell_id];

    /* Velocity profile */
    cs_real_t u_norm = ustar0 / kappa * cs_mo_psim(z+z0, z0, dlmo);

    cpro_met_vel[cell_id][0] = - sin(angle * cs_math_pi/180.) * u_norm;
    cpro_met_vel[cell_id][1] = - cos(angle * cs_math_pi/180.) * u_norm;

    /* Potential temperature profile
     * Note: same roughness as dynamics */
    cpro_met_potemp[cell_id] = theta0 + tstar / kappa * cs_mo_psih(z+z0, z0, dlmo);

    /* Richardson flux number profile */
    // Note : ri_f = z/(Pr_t L) * phih/phim^2 = z/Lmo * phim
    cs_real_t ri_f = (z+z0) * dlmo / cs_mo_phim(z+z0, dlmo);

    /* TKE profile */
    cpro_met_k[cell_id] = cs_math_pow2(ustar0) / sqrt(cmu)
      * sqrt(1. - CS_MIN(ri_f, 1.));

    /* epsilon profile */
    cpro_met_eps[cell_id] = cs_math_pow3(ustar0) / (kappa * (z+z0))
      * cs_mo_phim(z+z0, dlmo)*(1.-ri_f); //FIXME min (1, ri_f) ?

    /* TODO compute hydrostatic pressure and density profiles
     * with Laplace integration  (see atlecm.f90) */
    cs_atmo_profile_std((z+z0),
                        &(cpro_met_p[cell_id]),
                        &(cpro_met_t[cell_id]),
                        &(cpro_met_rho[cell_id]));


    /* Very stable cases */
    if (ri_f > ri_max) {
      if (z < z_min) {
        //Ri_f is an increasing monotonic function, so the lowest value of
        //z for which Ri_f>Ri_max is needed
        z_min = z;
        u_met_min=u_norm;
        theta_met_min=cpro_met_potemp[cell_id];
      }
    }
  }

  /* Very stable cases, corresponding to mode 0 in the Python prepro */
  if (z_min < cs_math_big_r) { // Clipping only if there are cells to be clipped
    bft_printf("Switching to very stable clipping for meteo profile.\n");
    bft_printf("All altitudes above %f have been modified by clipping.\n",z_min);
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
      cs_real_t z = cell_cen[cell_id][2] - z_ground[cell_id];
      if (z >= z_min) {
         /* mode = 0 is ustar=cst */
        dlmo_var[cell_id] = dlmo * (z_min + z0) / (z + z0);

        /* Velocity profile */
        cs_real_t u_norm = u_met_min + ustar0 / kappa * cs_mo_phim(z_min + z0, dlmo) * log((z+z0) / (z_min+z0));

        cpro_met_vel[cell_id][0] = - sin(angle * cs_math_pi/180.) * u_norm;
        cpro_met_vel[cell_id][1] = - cos(angle * cs_math_pi/180.) * u_norm;

       /* Potential temperature profile
        * Note: same roughness as dynamics */
        cpro_met_potemp[cell_id] = theta_met_min
          + tstar * (z_min+z0) / kappa * cs_mo_phih(z_min+z0, dlmo) * (-1./(z+z0) + 1./(z_min+z0)) ;
       /* TKE profile
          ri_max is necessarily lower than 1, but CS_MIN might be useful if
          that changes in the future */
        cpro_met_k[cell_id] = cs_math_pow2(ustar0) / sqrt(cmu)
          * sqrt(1. - CS_MIN(ri_max, 1.));

        /* epsilon profile */
        cpro_met_eps[cell_id] = cs_math_pow3(ustar0) / kappa  / dlmo_var[cell_id]
         * (1- CS_MIN(ri_max, 1.)) / CS_MIN(ri_max, 1.);
      }
    }
  }
  BFT_FREE(dlmo_var);
}

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
cs_atmo_z_ground_compute(void)
{
  /* Initialization
   *===============*/

  if (!_atmo_option.compute_z_ground)
    return;

  const cs_domain_t *domain = cs_glob_domain;
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_real_3_t *restrict i_face_normal =
     (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal =
     (const cs_real_3_t *restrict)mq->b_face_normal;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;

  const int *bc_type = cs_glob_bc_type;

  /* Pointer to z_ground field */
  cs_field_t *f = cs_field_by_name_try("z_ground");

  cs_real_t *restrict i_massflux
    = cs_field_by_id
        (cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id")))->val;
  cs_real_t *restrict b_massflux
    = cs_field_by_id
        (cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id")))->val;

  cs_var_cal_opt_t vcopt;
  cs_field_get_key_struct(f, cs_field_key_id("var_cal_opt"), &vcopt);

  cs_real_t normal[3];
  /* Normal direction is given by the gravity */
  cs_math_3_normalise((const cs_real_t *)(cs_glob_physical_constants->gravity),
                      normal);

  for (int i = 0; i < 3; i++)
    normal[i] *= -1.;

  /* Compute the mass flux due to V = - g / ||g||
   *=============================================*/

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
    i_massflux[face_id] = cs_math_3_dot_product(normal, i_face_normal[face_id]);

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++)
    b_massflux[face_id] = cs_math_3_dot_product(normal, b_face_normal[face_id]);

  /* Boundary conditions
   *====================*/

  cs_real_t norm = 0.;
  cs_real_t ground_surf = 0.;

  /* Dirichlet at walls, homogeneous Neumann elsewhere */

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    /* Dirichlet BCs */
    if ((bc_type[face_id] == CS_SMOOTHWALL || bc_type[face_id] == CS_ROUGHWALL)
        && b_massflux[face_id] <= 0.) {

      vcopt.ndircl = 1;
      cs_real_t hint = 1. / mq->b_dist[face_id];
      cs_real_t pimp = cs_math_3_dot_product(b_face_cog[face_id], normal);

      cs_boundary_conditions_set_dirichlet_scalar(&(f->bc_coeffs->a[face_id]),
                                                  &(f->bc_coeffs->af[face_id]),
                                                  &(f->bc_coeffs->b[face_id]),
                                                  &(f->bc_coeffs->bf[face_id]),
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);
      norm += cs_math_pow2(f->bc_coeffs->a[face_id]) * mq->b_face_surf[face_id];
      ground_surf += mq->b_face_surf[face_id];
    }
    /* Neumann Boundary Conditions */
    else {

      cs_real_t hint = 1. / mq->b_dist[face_id];
      cs_real_t qimp = 0.;

      cs_boundary_conditions_set_neumann_scalar(&(f->bc_coeffs->a[face_id]),
                                                &(f->bc_coeffs->af[face_id]),
                                                &(f->bc_coeffs->b[face_id]),
                                                &(f->bc_coeffs->bf[face_id]),
                                                qimp,
                                                hint);

    }
  }

  cs_parall_max(1, CS_INT_TYPE, &(vcopt.ndircl));

  /* Matrix
   *=======*/

  cs_real_t *rovsdt, *dpvar;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    rovsdt[cell_id] = 0.;

  /* Right hand side
   *================*/

  cs_real_t *rhs;
  BFT_MALLOC(rhs, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id< m->n_cells_with_ghosts; cell_id++)
    rhs[cell_id] = 0.;

  /* Norm
   *======*/

  cs_parall_sum(1, CS_REAL_TYPE, &norm);
  cs_parall_sum(1, CS_REAL_TYPE, &ground_surf);

  if (ground_surf > 0.)
    norm = sqrt(norm / ground_surf) * mq->tot_vol;
  else {
    bft_printf("No ground BC or no gravity: no computation of ground elevation.\n");
    return;
  }

  /* Solving
   *=========*/

  /* In case of a theta-scheme, set theta = 1;
     no relaxation in steady case either */

  cs_equation_iterative_solve_scalar(0,   /* idtvar: no steady state algo */
                                     -1,  /* no over loops */
                                     f->id,
                                     f->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     norm,
                                     &vcopt,
                                     f->val_pre,
                                     f->val,
                                     f->bc_coeffs->a,
                                     f->bc_coeffs->b,
                                     f->bc_coeffs->af,
                                     f->bc_coeffs->bf,
                                     i_massflux,
                                     b_massflux,
                                     i_massflux, /* viscosity, not used */
                                     b_massflux, /* viscosity, not used */
                                     i_massflux, /* viscosity, not used */
                                     b_massflux, /* viscosity, not used */
                                     NULL,
                                     NULL,
                                     NULL,
                                     0, /* icvflb (upwind) */
                                     NULL,
                                     rovsdt,
                                     rhs,
                                     f->val,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Free memory */

  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes hydrostatic profiles of density and pressure
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \divs \left( \grad \varia \right)
 *      = \divs \left( \dfrac{\vect{g}}{c_p \theta} \right) \varia 0
 *  \f]
 *  where \f$ \vect{g} \f$ is the gravity field and \f$ \theta \f$
 *  is the potential temperature.
 *
 *  The boundary conditions on \f$ \varia \f$ read:
 *  \f[
 *   \varia = \left(\dfrac{P_{sea}}{p_s}\right)^{R/C_p} \textrm{on the ground}
 *  \f]
 *  and Neumann elsewhere.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_hydrostatic_profiles_compute(void)
{
 //TODO
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the meteo file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_meteo_file_name(const char *file_name)
{
  if (file_name == NULL) {
    return;
  }

  BFT_MALLOC(_atmo_option.meteo_file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_atmo_option.meteo_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the SPACK file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_spack_file_name(const char *file_name)
{
  if (file_name == NULL) {
    _atmo_chem.model = 0;
    return;
  }

  _atmo_chem.model = 4;

  BFT_MALLOC(_atmo_chem.spack_file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_atmo_chem.spack_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function sets the file name to initialize the aerosol library.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_aerosol_file_name(const char *file_name)
{
  if (file_name == NULL) {
    _atmo_chem.aerosol_model = CS_ATMO_AEROSOL_OFF;
    return;
  }

  _atmo_chem.aerosol_model = CS_ATMO_AEROSOL_SSH;

  BFT_MALLOC(_atmo_chem.aero_file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_atmo_chem.aero_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function declares additional transported variables for
 *        atmospheric module for the chemistry defined from SPACK.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_declare_chem_from_spack(void)
{
  assert(_atmo_chem.model == 4);

  if (_atmo_chem.spack_file_name == NULL)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo chemistry from SPACK file: requires a SPACK file."));

  char line[512] = "";

  /* Open file */
  bft_printf("SPACK file for atmo chemistry:\n    %s \n",
             _atmo_chem.spack_file_name);

  FILE *file = fopen(_atmo_chem.spack_file_name, "rt");
  if (file == NULL)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo chemistry from SPACK file: Could not open file."));

  int line_num = 0;

  /* Read "[species]" */
  while (true) {
    line_num++;
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not skip header."));

    if (strncmp("[species]", line, 512) == 0)
      break;
  }

  /* Read SPACK: first loop count the number of species */
  for (int i = 1; true; i++ ) {
    /* Read species */
    line_num++;
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not read line %d."),
                line_num);

    /* When reach [molecular_waight]: break */
    if (strncmp("[molecular_weight]", line, 512) == 0)
      break;
    else {
      _atmo_chem.n_species = i;
    }
  }

  /* Now allocate arrays */
  BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
  BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);
  BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);
  BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);

  /* Read SPACK: second loop Create variables and read molar mass */
  for (int i = 0; i < _atmo_chem.n_species; i++ ) {
    char name[512] = "";
    char label[512] = "";

    /* Read species name and molar mass */
    line_num++;
    if (fscanf(file, "%s %lf\n", line, &(_atmo_chem.molar_mass[i])) != 2)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmospheric chemistry from SPACK file:\n"
                  "  could not read species name and molar mass, line %d."),
                line_num);

    /* The order is already ok */
    _atmo_chem.chempoint[i] = i+1; //FIXME ?

    /* Build name of the field:
     * species_name in lower case */
    strcpy(name, "species_");
    _strtolower(label, line);
    strcat(name, label);

    /* Field of dimension 1 */
    /* Give the original name as label */
    _atmo_chem.species_to_field_id[i]
      = cs_variable_field_create(name, line, CS_MESH_LOCATION_CELLS, 1);

    /* Scalar field, store in isca_chem/species_to_scalar_id (FORTRAN/C) array */
    _atmo_chem.species_to_scalar_id[i]
      = cs_add_model_field_indexes(_atmo_chem.species_to_field_id[i]);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief 1D Radiative scheme - Solar data + zenithal angle)
 *
 * Compute:
 *   - zenithal angle
 *   - solar contant (with correction for earth - solar length)
 *   - albedo if above the sea
 *   (Use analytical formulae of Paltrige and Platt
 *              dev.in atm. science no 5)
 * \param[in]   latitude    latitude
 * \param[in]   longitude   longitude
 * \param[in]   squant      start day in the year
 * \param[in]   utc         Universal time (hour)
 * \param[in]   sea_id      sea index
 * \param[out]  albedo      albedo
 * \param[out]  muzero      cosin of zenithal angle
 * \param[out]  omega       solar azimut angle
 * \param[out]  fo          solar constant
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_compute_solar_angles(cs_real_t latitude,
                             cs_real_t longitude,
                             cs_real_t squant,
                             cs_real_t utc,
                             int       sea_id,
                             cs_real_t *albedo,
                             cs_real_t *muzero,
                             cs_real_t *omega,
                             cs_real_t *fo)
{
  /* 1 - initialisations */
  *fo = 1370.;

  /* conversions sexagesimal-decimal */

  cs_real_t flat = latitude  *cs_math_pi/180.;
  cs_real_t flong = longitude * 4. / 60.;

  cs_real_t t00 = 2. * cs_math_pi * squant/365.;

  /* 2 - compute declinaison (maximum error < 3 mn) */

  cs_real_t decl
    =   0.006918 - 0.399912*cos(t00) + 0.070257*sin(t00)
      - 0.006758*cos(2.*t00) + 0.000907*sin(2.*t00) - 0.002697*cos(3.*t00)
      + 0.001480*sin(3.*t00);

  /* 3 - compute local solar hour
   * equation du temps     erreur maxi    35 secondes
   */

  cs_real_t eqt
    = (  0.000075 + 0.001868*cos(t00) - 0.032077*sin(t00)
       - 0.014615*cos(2.*t00) - 0.040849*sin(2.*t00)) * 12./cs_math_pi;

  cs_real_t local_time = utc + flong + eqt;

  /* Transformation local_time-radians */

  /* On retire cs_math_pi et on prend le modulo 2pi du resultat */
  cs_real_t hr = (local_time - 12.)*cs_math_pi/12.;
  if (local_time < 12.)
    hr = (local_time + 12.)*cs_math_pi/12.;

  /* 4 - compute of cosinus of the zenitghal angle */

  *muzero = sin(decl)*sin(flat) + cos(decl)*cos(flat)*cos(hr);

  cs_real_t za = acos(*muzero);

  /* 5 - compute solar azimut */
  *omega = 0.;
  if (CS_ABS(sin(za)) > cs_math_epzero) {
    /* Cosinus of the zimut angle */
    cs_real_t co = (sin(decl)*cos(flat)-cos(decl)*sin(flat)*cos(hr))/sin(za);
    *omega = acos(co);
    if (local_time > 12.)
      *omega = 2. * cs_math_pi - acos(co);
  }
  *omega -= cs_glob_atmo_option->domain_orientation * cs_math_pi / 180.;

  /* 5 - compute albedo at sea which depends on the zenithal angle */

  if (sea_id == 1) {
    cs_real_t ho = acos(*muzero);
    ho = 180.*(cs_math_pi/2. - ho)/cs_math_pi;
    if (ho < 8.5)
      ho = 8.5;
    if (ho > 60.)
      ho = 60.;
    *albedo = 3./ho;
  }

 /* 6 - Compute solar constant
    distance correction earth-sun
    corfo=(r0/r)**2
    precision better than e-04 */

  cs_real_t corfo
    =   1.000110 + 0.034221*cos(t00) + 0.001280*sin(t00)
      + 0.000719*cos(2.*t00) + 0.000077*sin(2.*t00);
  *fo *= corfo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_OFF)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Atmospheric module options\n"
                  "--------------------------\n\n"));

  switch(cs_glob_physical_model_flag[CS_ATMOSPHERIC]) {
    case CS_ATMO_CONSTANT_DENSITY:
      /* Constant density */
      cs_log_printf(CS_LOG_SETUP,
                  _("  Constant density\n\n"));
      break;
    case CS_ATMO_DRY:
      /* Dry */
      cs_log_printf(CS_LOG_SETUP,
                  _("  Dry atmosphere\n\n"));
      break;
    case CS_ATMO_HUMID:
      /* Humid */
      cs_log_printf(CS_LOG_SETUP,
                  _("  Humid atmosphere\n\n"));
      break;
    default:
      break;
  }

  if (cs_glob_atmo_option->compute_z_ground > 0)
    cs_log_printf(CS_LOG_SETUP,
        _("  Compute ground elevation\n\n"));

  if (cs_glob_atmo_option->open_bcs_treatment > 0)
    cs_log_printf(CS_LOG_SETUP,
        _("  Impose open BCs with momentum source terms\n"));

  if (cs_glob_atmo_option->open_bcs_treatment == 2)
    cs_log_printf(CS_LOG_SETUP,
        _("  and impose profiles at ingoing faces\n\n"));

  /* CUT */
  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "  Starting Coordinated Universal Time (CUT):\n"
       "    Year:      %4d\n"
       "    Quant:     %4d\n"
       "    Hour:      %4d\n"
       "    Min:       %4d\n"
       "    Sec:       %4f\n\n"),
     cs_glob_atmo_option->syear,
     cs_glob_atmo_option->squant,
     cs_glob_atmo_option->shour,
     cs_glob_atmo_option->smin,
     cs_glob_atmo_option->ssec);

  /* Centre of the domain latitude */
  cs_log_printf
    (CS_LOG_SETUP,
     _("  Domain center:\n"
       "    latitude:  %6f\n"
       "    longitude: %6f\n\n"),
     cs_glob_atmo_option->latitude,
     cs_glob_atmo_option->longitude);

  if (cs_glob_atmo_option->meteo_profile == 1) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Large scale Meteo file: %s\n\n"),
       cs_glob_atmo_option->meteo_file_name);
  }


  if (cs_glob_atmo_option->meteo_profile == 2) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Large scale Meteo profile info:\n"
         "    roughness: %12f [m]\n"
         "    LMO inv:   %12f [m^-1]\n"
         "    ustar0:    %12f [m/s]\n"
         "    uref:      %12f [m/s]\n"
         "    zref:      %12f [m]\n"
         "    angle:     %12f [°]\n"
         "    P sea:     %12f [Pa]\n"
         "    T0:        %12f [K]\n"
         "    Tstar:     %12f [K]\n"
         "    BL height: %12f [m]\n\n"),
       cs_glob_atmo_option->meteo_z0,
       cs_glob_atmo_option->meteo_dlmo,
       cs_glob_atmo_option->meteo_ustar0,
       cs_glob_atmo_option->meteo_uref,
       cs_glob_atmo_option->meteo_zref,
       cs_glob_atmo_option->meteo_angle,
       cs_glob_atmo_option->meteo_psea,
       cs_glob_atmo_option->meteo_t0,
       cs_glob_atmo_option->meteo_tstar,
       cs_glob_atmo_option->meteo_zi);
  }

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric chemistry options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 0)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Atmospheric chemistry options\n"
                  "-----------------------------\n\n"));

  if (cs_glob_atmo_chemistry->model == 0) {

    /* No atmospheric chemistry */
    cs_log_printf(CS_LOG_SETUP,
                  _("  No atmospheric chemistry\n\n"));

  }
  else if (   cs_glob_atmo_chemistry->model == 1
           || cs_glob_atmo_chemistry->model == 2
           || cs_glob_atmo_chemistry->model == 3) {

    /* Pre-defined schemes */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric chemistry activated\n\n"
         "    Pre-defined scheme %12d\n\n"
         "      n_species: %18d (Number of species)\n"
         "      n_reactions: %16d (Number of reactions)\n"
         "      photolysis: %17s\n"
         "      frozen_gas_chem: %12s\n\n"),
       cs_glob_atmo_chemistry->model,
       cs_glob_atmo_chemistry->n_species,
       cs_glob_atmo_chemistry->n_reactions,
       cs_glob_atmo_chemistry->chemistry_with_photolysis ? "Yes": "No",
       cs_glob_atmo_chemistry->frozen_gas_chem ? "Yes": "No");

  }
  else if (cs_glob_atmo_chemistry->model == 4) {

    /* User defined SPACK chemistry */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric chemistry activated\n\n"
         "    User-defined SPACK scheme\n\n"
         "      n_species: %18d (Number of species)\n"
         "      n_reactions: %16d (Number of reactions)\n"
         "      photolysis: %17s\n"
         "      frozen_gas_chem: %12s\n"
         "      Spack file: %17s\n"),
       cs_glob_atmo_chemistry->n_species,
       cs_glob_atmo_chemistry->n_reactions,
       cs_glob_atmo_chemistry->chemistry_with_photolysis ? "Yes": "No",
       cs_glob_atmo_chemistry->frozen_gas_chem ? "Yes": "No",
       cs_glob_atmo_chemistry->spack_file_name);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric aerosols options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 0)
    return;

  cs_log_printf
    (CS_LOG_SETUP,
     _("\nAtmospheric aerosols options\n"
       "----------------------------\n\n"));

  cs_atmo_aerosol_type_t atm_aer_type = cs_glob_atmo_chemistry->aerosol_model;

  if (atm_aer_type == CS_ATMO_AEROSOL_OFF)
    cs_log_printf(CS_LOG_SETUP,_("  %s\n"),
                  cs_atmo_aerosol_type_name[atm_aer_type]);

  else if (atm_aer_type == CS_ATMO_AEROSOL_SSH) {

    /* Atmospheric chemistry with SSH */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric aerosols activated\n\n"
         "    Global parameters\n\n"
         "      model: %22s (%s)\n"
         "      n_layer: %20d (Number of layers inside aerosols)\n"
         "      n_size:  %20d (Number of resolved aerosols)\n"
         "      init_gas_with_lib: %10s\n"
         "      init_aero_with_lib: %9s\n"
         "      namelist: %s\n\n"),
       cs_atmo_aerosol_type_enum_name[atm_aer_type],
       cs_atmo_aerosol_type_name[atm_aer_type],
       cs_glob_atmo_chemistry->n_layer,
       cs_glob_atmo_chemistry->n_size,
       cs_glob_atmo_chemistry->init_gas_with_lib ? "Yes": "No",
       cs_glob_atmo_chemistry->init_aero_with_lib ? "Yes": "No",
       cs_glob_atmo_chemistry->aero_file_name);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
