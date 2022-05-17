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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_atmo_profile_std.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
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
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_volume_zone.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_convection_diffusion.h"
#include "cs_parameters.h"
#include "cs_porous_model.h"
#include "cs_timer.h"
#include "cs_matrix_building.h"
#include "cs_sles.h"
#include "cs_sles_default.h"
#include "cs_face_viscosity.h"
#include "cs_divergence.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_velocity_pressure.h"

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
  .x_l93= 1e12,
  .y_l93 = 1e12,
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
  .z_dyn_met  = NULL,
  .z_temp_met = NULL,
  .u_met      = NULL,
  .v_met      = NULL,
  .time_met   = NULL,
  .hyd_p_met  = NULL,
  .pot_t_met  = NULL
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
  .aero_file_name = NULL,
  .chem_conc_file_name = NULL,
  .aero_conc_file_name = NULL
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

void
cs_f_atmo_get_meteo_file_name(int           name_max,
                              const char  **name,
                              int          *name_len);

void
cs_f_atmo_get_chem_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len);

void
cs_f_atmo_get_aero_conc_file_name(int           name_max,
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
                       cs_real_t              **x_l93,
                       cs_real_t              **y_l93,
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
                       int                    **nbmaxt,
                       cs_real_t              **meteo_zi);

void
cs_f_atmo_arrays_get_pointers(cs_real_t **z_dyn_met,
                              cs_real_t **z_temp_met,
                              cs_real_t **u_met,
                              cs_real_t **v_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              cs_real_t **pot_t_met,
                              int         dim_u_met[2],
                              int         dim_hyd_p_met[2],
                              int         dim_pot_t_met[2]);

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint);

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid);

void
cs_f_atmo_chem_finalize(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the name of the meteo file
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
 * Return the name of the chemistry concentration file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_chem_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len)
{
  *name = _atmo_chem.chem_conc_file_name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving chemistry concentration file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       _atmo_chem.chem_conc_file_name, name_max, *name, *name_len);
  }
}

/*----------------------------------------------------------------------------
 * Return the name of the aerosol concentration file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_aero_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len)
{
  *name = _atmo_chem.aero_conc_file_name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving chemistry concentration file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       _atmo_chem.aero_conc_file_name, name_max, *name, *name_len);
  }
}

/*----------------------------------------------------------------------------
 * Access pointers for Fortran mapping.
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
                       cs_real_t              **x_l93,
                       cs_real_t              **y_l93,
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
                       int                    **nbmaxt,
                       cs_real_t              **meteo_zi)
{
  *ps        = &(_atmo_constants.ps);
  *syear     = &(_atmo_option.syear);
  *squant    = &(_atmo_option.squant);
  *shour     = &(_atmo_option.shour);
  *smin      = &(_atmo_option.smin);
  *ssec      = &(_atmo_option.ssec);
  *longitude = &(_atmo_option.longitude);
  *latitude  = &(_atmo_option.latitude);
  *x_l93 = &(_atmo_option.x_l93);
  *y_l93 = &(_atmo_option.y_l93);
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
  *meteo_zi   = &(_atmo_option.meteo_zi);
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
cs_f_atmo_arrays_get_pointers(cs_real_t **z_dyn_met,
                              cs_real_t **z_temp_met,
                              cs_real_t **u_met,
                              cs_real_t **v_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              cs_real_t **pot_t_met,
                              int         dim_u_met[2],
                              int         dim_hyd_p_met[2],
                              int         dim_pot_t_met[2])
{
  int n_level = 0, n_level_t = 0;
  int n_times = 0;
  if (_atmo_option.meteo_profile) {
    n_level = CS_MAX(1, _atmo_option.nbmetd);
    n_level_t = CS_MAX(1, _atmo_option.nbmaxt);
    n_times = CS_MAX(1, _atmo_option.nbmetm);
  }

  if (_atmo_option.z_dyn_met == NULL)
    BFT_MALLOC(_atmo_option.z_dyn_met, n_level, cs_real_t);
  if (_atmo_option.z_temp_met == NULL)
    BFT_MALLOC(_atmo_option.z_temp_met, _atmo_option.nbmaxt, cs_real_t);
  if (_atmo_option.u_met == NULL)
    BFT_MALLOC(_atmo_option.u_met, n_level*n_times, cs_real_t);
  if (_atmo_option.v_met == NULL)
    BFT_MALLOC(_atmo_option.v_met, n_level*n_times, cs_real_t);
  if (_atmo_option.time_met == NULL)
    BFT_MALLOC(_atmo_option.time_met, _atmo_option.nbmetm, cs_real_t);
  if (_atmo_option.hyd_p_met == NULL)
    BFT_MALLOC(_atmo_option.hyd_p_met,
               _atmo_option.nbmetm*_atmo_option.nbmaxt, cs_real_t);
  if (_atmo_option.pot_t_met == NULL)
    BFT_MALLOC(_atmo_option.pot_t_met, n_level_t*n_times, cs_real_t);

  *u_met           = _atmo_option.u_met;
  *v_met           = _atmo_option.v_met;
  *hyd_p_met       = _atmo_option.hyd_p_met;
  *pot_t_met       = _atmo_option.pot_t_met;
  dim_u_met[0]     = _atmo_option.nbmetd;
  dim_u_met[1]     = _atmo_option.nbmetm;
  dim_hyd_p_met[0] = _atmo_option.nbmaxt;
  dim_hyd_p_met[1] = _atmo_option.nbmetm;
  dim_pot_t_met[0] = _atmo_option.nbmaxt;
  dim_pot_t_met[1] = _atmo_option.nbmetm;

  *z_dyn_met  = _atmo_option.z_dyn_met;
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
  BFT_FREE(_atmo_chem.chem_conc_file_name);
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief This auxiliary function solve a Poisson equation for hydrostatic
 *  pressure :
 *  \f$ \divs ( \grad P ) = \divs ( f ) \f$
 * \param[in]      f_ext       external forcing
 * \param[in]      dfext       external forcing increment
 * \param[in]      name        field name
 * \param[in]      min_z       Minimum altitude of the domain
 * \param[in]      p_ground    Pressure at the minimum altitude
 * \param[in,out]  i_massflux  Internal mass flux
 * \param[in,out]  b_massflux  Boundary mass flux
 * \param[in,out]  i_viscm     Internal face viscosity
 * \param[in,out]  i_viscm     Boundary face viscosity
 * \param[in,out]  dam         Working array
 * \param[in,out]  xam         Working array
 * \param[in,out]  dpvar       Pressure increment
 * \param[in,out]  rhs         Working array
 */
/*----------------------------------------------------------------------------*/

static void
_hydrostatic_pressure_compute(cs_real_3_t  f_ext[],
                              cs_real_3_t  dfext[],
                              cs_real_t    pvar[],
                              const char  *name,
                              cs_real_t    min_z,
                              cs_real_t    p_ground,
                              cs_real_t    i_massflux[],
                              cs_real_t    b_massflux[],
                              cs_real_t    i_viscm[],
                              cs_real_t    b_viscm[],
                              cs_real_t    dam[],
                              cs_real_t    xam[],
                              cs_real_t    dpvar[],
                              cs_real_t    rhs[])

{
  /* Local variables */
  cs_domain_t *domain = cs_glob_domain;
  cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  cs_field_t *f = cs_field_by_name(name);
  cs_equation_param_t *eqp_p = cs_field_get_equation_param(f);
  int f_id = f->id;
  int niterf;

  /*==========================================================================
   * 0.  Initialization
   *==========================================================================*/

  /* solving info */
  cs_solving_info_t sinfo;
  int key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    cs_field_get_key_struct(f, key_sinfo_id, &sinfo);
  }

  /* Symmetric matrix, except if advection */
  int isym = 1;
  bool symmetric = true;

  cs_real_3_t *next_fext;
  BFT_MALLOC(next_fext, m->n_cells_with_ghosts, cs_real_3_t);
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    next_fext[cell_id][0] = f_ext[cell_id][0] + dfext[cell_id][0];
    next_fext[cell_id][1] = f_ext[cell_id][1] + dfext[cell_id][1];
    next_fext[cell_id][2] = f_ext[cell_id][2] + dfext[cell_id][2];
  }

  /* --> Handle parallelism and periodicity */
  if (cs_glob_rank_id  >= 0 || cs_glob_mesh->n_init_perio > 0)
    cs_mesh_sync_var_vect((cs_real_t *)next_fext);

  /* Boundary conditions
   *====================*/

  eqp_p->ndircl = 0;

  /* To solve hydrostatic pressure:
   * p_ground on the lowest face, homogeneous Neumann everywhere else */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    cs_real_t hint = 1. / mq->b_dist[face_id];
    cs_real_t qimp = 0.;
    if ((b_face_cog[face_id][2] - min_z) < 0.005) {//FIXME dimensionless value
      cs_real_t pimp = p_ground;
      cs_boundary_conditions_set_dirichlet_scalar(&(f->bc_coeffs->a[face_id]),
                                                  &(f->bc_coeffs->af[face_id]),
                                                  &(f->bc_coeffs->b[face_id]),
                                                  &(f->bc_coeffs->bf[face_id]),
                                                  pimp,
                                                  hint,
                                                  cs_math_big_r);
      eqp_p->ndircl = 1;
    }
    else {
      cs_boundary_conditions_set_neumann_scalar(&(f->bc_coeffs->a[face_id]),
                                                &(f->bc_coeffs->af[face_id]),
                                                &(f->bc_coeffs->b[face_id]),
                                                &(f->bc_coeffs->bf[face_id]),
                                                qimp,
                                                hint);
    }
  }

  cs_parall_max(1, CS_INT_TYPE, &(eqp_p->ndircl));
  cs_real_t *rovsdt;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    rovsdt[cell_id] = 0.;

  /* Faces viscosity */
  cs_real_t *c_visc;
  BFT_MALLOC(c_visc, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    c_visc[cell_id] = 1.;

  cs_face_viscosity(m, mq, eqp_p->imvisf, c_visc, i_viscm, b_viscm);

  cs_matrix_wrapper_scalar(eqp_p->iconv,
                           eqp_p->idiff,
                           eqp_p->ndircl,
                           isym,
                           eqp_p->thetav,
                           0,
                           f->bc_coeffs->b,
                           f->bc_coeffs->bf,
                           rovsdt,
                           i_massflux,
                           b_massflux,
                           i_viscm,
                           b_viscm,
                           NULL,
                           dam,
                           xam);

  /* Right hand side
   *================*/

  cs_ext_force_flux(m,
                    mq,
                    1, /* init */
                    eqp_p->nswrgr,
                    next_fext,
                    f->bc_coeffs->bf,
                    i_massflux,
                    b_massflux,
                    i_viscm,
                    b_viscm,
                    c_visc,
                    c_visc,
                    c_visc);

  cs_real_t *divergfext;
  BFT_MALLOC(divergfext, m->n_cells_with_ghosts, cs_real_t);

  cs_divergence(m,
                1, /* init */
                i_massflux,
                b_massflux,
                divergfext);

  /* --- Right hand side residual */
  cs_real_t rnorm = sqrt(cs_gdot(m->n_cells, divergfext, divergfext));
  cs_real_t residu = rnorm;

  /* Initial Right-Hand-Side */
  cs_diffusion_potential(f_id,
                         m,
                         mq,
                         1, /* init */
                         1, /* inc */
                         eqp_p->imrgra,
                         1, /* iccocg */
                         eqp_p->nswrgr,
                         eqp_p->imligr,
                         1, /* iphydp */
                         eqp_p->iwgrec,
                         eqp_p->iwarni,
                         eqp_p->epsrgr,
                         eqp_p->climgr,
                         next_fext,
                         pvar,
                         f->bc_coeffs->a,
                         f->bc_coeffs->b,
                         f->bc_coeffs->af,
                         f->bc_coeffs->bf,
                         i_viscm,
                         b_viscm,
                         c_visc,
                         rhs);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
    rhs[cell_id] = - divergfext[cell_id] - rhs[cell_id];

  cs_real_t ressol = 0;
  for (int sweep = 0;
       sweep < eqp_p->nswrsm && residu > (rnorm * eqp_p->epsrsm);
       sweep++) {

    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      dpvar[cell_id] = 0.;

    ressol = residu;
    cs_sles_solve_native(f_id,
                         "",
                         symmetric,
                         1, /* db_size */
                         1, /* eb_size */
                         dam,
                         xam,
                         eqp_p->epsilo,
                         rnorm,
                         &niterf,
                         &ressol,
                         rhs,
                         dpvar);

    /* Update variable and right-hand-side */
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      pvar[cell_id] += dpvar[cell_id];

    cs_diffusion_potential(f_id,
                           m,
                           mq,
                           1, /* init */
                           1, /* inc */
                           eqp_p->imrgra,
                           1, /* iccocg */
                           eqp_p->nswrgr,
                           eqp_p->imligr,
                           1, /* iphydp */
                           eqp_p->iwgrec,
                           eqp_p->iwarni,
                           eqp_p->epsrgr,
                           eqp_p->climgr,
                           next_fext,
                           pvar,
                           f->bc_coeffs->a,
                           f->bc_coeffs->b,
                           f->bc_coeffs->af,
                           f->bc_coeffs->bf,
                           i_viscm,
                           b_viscm,
                           c_visc,
                           rhs);

    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      rhs[cell_id] = - divergfext[cell_id] - rhs[cell_id];

    /* --- Convergence test */
    residu = sqrt(cs_gdot(m->n_cells, rhs, rhs));

    /* Writing */
    if (eqp_p->iwarni >= 2) {
      bft_printf("%s: CV_DIF_TS, IT: %d, Res: %12.5e, Norm: %12.5e\n",
          name, sweep, residu, rnorm);
      bft_printf("%s: Current reconstruction sweep: %d, "
          "Iterations for solver: %d\n", name, sweep, niterf);
    }

  }

  cs_sles_free_native(f_id, "");

  BFT_FREE(divergfext);
  BFT_FREE(next_fext);
  BFT_FREE(rovsdt);
  BFT_FREE(c_visc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert coordinates from Lambert-93 to WGS84 inside atmo options
 */
/*----------------------------------------------------------------------------*/

static void
_convert_from_l93_to_wgs84(void)
{
  //computation from https://georezo.net/forum/viewtopic.php?id=94465
  cs_real_t c = 11754255.426096; // projection constant
  cs_real_t e = 0.0818191910428158; // ellipsoid excentricity
  cs_real_t n = 0.725607765053267; // projection exponent
  cs_real_t xs = 700000; // projection's pole x-coordinate
  cs_real_t ys = 12655612.049876;  // projection's pole y-coordinate
  cs_real_t a = (log(c/(sqrt(  cs_math_pow2(cs_glob_atmo_option->x_l93-xs)
                             + cs_math_pow2(cs_glob_atmo_option->y_l93-ys))))/n);

  double t1 = a + e*atanh(e*(tanh(a+e*atanh(e*sin(1)))));
  double t2 = e*tanh(a+e*atanh(e*(tanh(t1))));
  double t3 = a+e*atanh(e*tanh(a+e*atanh(e*tanh(a+e*atanh(t2)))));

  cs_glob_atmo_option->longitude = ((atan(-(cs_glob_atmo_option->x_l93-xs)
                                          /(cs_glob_atmo_option->y_l93-ys)))/n
                                    + 3*cs_math_pi/180)/cs_math_pi*180;
  cs_glob_atmo_option->latitude
    = asin(tanh((log(c/sqrt(  cs_math_pow2(cs_glob_atmo_option->x_l93-xs)
                            + cs_math_pow2(cs_glob_atmo_option->y_l93-ys)))/n)
                +e*atanh(e*tanh(t3))))/cs_math_pi*180;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert coordinates from WGS84 to Lambert-93 inside atmo options
 */
/*----------------------------------------------------------------------------*/

static void
_convert_from_wgs84_to_l93(void)
{
  //computation from https://georezo.net/forum/viewtopic.php?id=94465
  cs_real_t  c = 11754255.426096; // projection constant
  cs_real_t  e = 0.0818191910428158; // ellipsoid excentricity
  cs_real_t  n = 0.725607765053267; // projection exponent
  cs_real_t  xs = 700000; // projection's pole x-coordinate
  cs_real_t  ys = 12655612.049876;  // projection's pole y-coordinate

  cs_real_t lat_rad= cs_glob_atmo_option->latitude*cs_math_pi/180; //latitude in rad
  cs_real_t lat_iso= atanh(sin(lat_rad))-e*atanh(e*sin(lat_rad)); // isometric latitude

  cs_glob_atmo_option->x_l93= ((c*exp(-n*(lat_iso)))
                               *sin(n*(cs_glob_atmo_option->longitude-3)
                                    *cs_math_pi/180)+xs);
  cs_glob_atmo_option->y_l93= (ys-(c*exp(-n*(lat_iso)))
                               *cos(n*(cs_glob_atmo_option->longitude-3)
                                    *cs_math_pi/180));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize meteo profiles if no meteo file is given.
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
  phys_pro->t0 = aopt->meteo_t0; /* ref temp T0 */
  phys_pro->ro0 = phys_pro->p0/(rair * aopt->meteo_t0); /* ref density T0 */

  cs_real_t z0 = aopt->meteo_z0;
  cs_real_t zref = aopt->meteo_zref;
  if (aopt->meteo_ustar0 < 0. && aopt->meteo_uref < 0.)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo meteo profiles: meteo_ustar0 or meteo_uref.\n"));

  /* Recompute LMO inverse */
  if (aopt->meteo_ustar0 >= 0. && aopt->meteo_uref >= 0.) {

    /* U+ */
    cs_real_t up = aopt->meteo_uref / aopt->meteo_ustar0;

    /* U+ if neutral, dlmo = 0 */
    cs_real_t up_l = cs_mo_psim(zref + z0,
                                z0,
                                0.) / kappa;

    cs_real_t dlmo = 0.;
    cs_real_t error = up_l - up;

    /* Dichotomy */
    cs_real_t dl_min = -1.e6;
    cs_real_t dl_max =  1.e6;
    cs_real_t tol = 1e-6;
    int it;
    int it_max = 1000;
    for (it = 0;
         it < it_max && CS_ABS(error) > tol && 0.5*(dl_max - dl_min) > tol;
         it++) {
      cs_real_t dl_mid = 0.5 * (dl_min + dl_max);

      cs_real_t error_min = cs_mo_psim(zref + z0,
                                       z0,
                                       dl_min) / kappa - up;
      cs_real_t error_mid = cs_mo_psim(zref + z0,
                                       z0,
                                       dl_mid) / kappa - up;

      /* The solution is between min and mid */
      if (error_min * error_mid < 0) {
        dl_max = dl_mid;
        if (CS_ABS(error_min) < CS_ABS(error_mid)) {
          dlmo = dl_min;
          error = error_min;
        }
        else {
          dlmo = dl_mid;
          error = error_mid;
        }
      }
      /* The solution is between mid and max */
      else {
        cs_real_t error_max = cs_mo_psim(zref + z0,
                                         z0,
                                         dl_max) / kappa - up;
        dl_min = dl_mid;
        if (CS_ABS(error_mid) < CS_ABS(error_max)) {
          dlmo = dl_mid;
          error = error_mid;
        }
        else {
          dlmo = dl_max;
          error = error_max;
        }
      }
#if 0
      bft_printf("IT %d: dlmo = %f, error = %f\n",it, dlmo, error);
#endif
    }

    if (it == it_max)
      bft_printf("Warning: meteo preprocessor did not converge to find inverse\n"
                 " of LMO length, current value is %f.\n"
                 "Please, check reference velocity, reference altitude and ustar\n",
                 dlmo);

    aopt->meteo_dlmo = dlmo;
  }

  /* Compute ground friction velocity from dlmo and uref */
  if (aopt->meteo_ustar0 < 0.)
    aopt->meteo_ustar0 =
      aopt->meteo_uref * kappa
      / cs_mo_psim(zref + z0,
                   z0,
                   aopt->meteo_dlmo);

  /* Compute uref from ground friction velocity and dlmo */
  if (aopt->meteo_uref < 0. && zref > 0.)
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

  /* Center of the domain */
  /* if neither latitude/longitude nor lambert coordinates are given */
  if (((aopt->latitude > 0.5 * cs_math_big_r || aopt->longitude > 0.5 * cs_math_big_r)
      && (aopt->x_l93 > 0.5 * cs_math_big_r || aopt->y_l93 > 0.5 * cs_math_big_r))) {
    bft_printf("Neither latitude nor center in Lambert-93 was given \n");
    bft_printf("It is set to Paris values \n");
    aopt->latitude = 45.44;
    aopt->longitude = 4.39;
    _convert_from_wgs84_to_l93();
  }
  /* else if latitude/longitude are given */
  else if (aopt->x_l93 > 0.5 * cs_math_big_r || aopt->y_l93 > 0.5 * cs_math_big_r) {
    bft_printf("Latitude and longitude were given, Lambert center's coordinates"
               "are automatically computed\n");
    _convert_from_wgs84_to_l93();
  }
  /* All other cases */
  else{
    bft_printf("Lambert coordinates were given, latitude"
                "and longitude are automatically computed\n");
    _convert_from_l93_to_wgs84();
  }

  /* BL height according to Marht 1982 formula */
  /* value of C=0.2, 0.185, 0.06, 0.14, 0.07, 0.04 */
  cs_real_t zi_coef = 0.2;
  cs_real_t corio_f = 4. * cs_math_pi / 86164.
                         * sin(aopt->latitude * cs_math_pi / 180.);
  aopt->meteo_zi = zi_coef * ustar0 / fabs(corio_f);

  /* Force the computation of z_ground */
  aopt->compute_z_ground = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute meteo profiles if no meteo file is given.
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

  /* In the log */
  bft_printf(" Computing meteo profiles from CS\n\n");

  /* Get fields */
  cs_real_t *cpro_met_potemp = cs_field_by_name("meteo_pot_temperature")->val;
  cs_real_3_t *cpro_met_vel =
    (cs_real_3_t *) (cs_field_by_name("meteo_velocity")->val);
  cs_real_t *cpro_met_k = cs_field_by_name("meteo_tke")->val;
  cs_real_t *cpro_met_eps = cs_field_by_name("meteo_eps")->val;

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
  cs_real_t z_lim = cs_math_big_r;
  cs_real_t u_met_min= cs_math_big_r;
  cs_real_t theta_met_min= cs_math_big_r;

  cs_real_t *z_ground = NULL;
  if (aopt->compute_z_ground == true) {

    cs_field_t *f_z_ground = cs_field_by_name("z_ground");

    /* Do not recompute in case of restart */
    int has_restart = cs_restart_present();
    if (has_restart == 1) {
      cs_restart_t *rp = cs_restart_create("main.csc",
                                           NULL,
                                           CS_RESTART_MODE_READ);

      int retval = cs_restart_read_field_vals(rp,
                                              f_z_ground->id,
                                              0);    /* current value */
      if (retval != CS_RESTART_SUCCESS)
        has_restart = 0;
    }

    /* z_ground needs to be computed? */
    if (has_restart == 0)
      cs_atmo_z_ground_compute();

    z_ground = f_z_ground->val;
  }

  BFT_MALLOC(dlmo_var, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++) {
    dlmo_var[cell_id] = 0.0;
  }
  /* For DRSM models store Rxz/k */
  cs_field_t *f_axz = cs_field_by_name_try("meteo_shear_anisotropy");

  if (dlmo > 0)
    ri_max = 0.75; // Value chosen to limit buoyancy vs shear production

  /* Profiles */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

    cs_real_t z_grd = 0.;
    if (z_ground != NULL)
      z_grd = z_ground[cell_id];

    /* Local elevation */
    cs_real_t z = cell_cen[cell_id][2] - z_grd;

    /* Velocity profile */
    cs_real_t u_norm = ustar0 / kappa * cs_mo_psim(z+z0, z0, dlmo);

    cpro_met_vel[cell_id][0] = - sin(angle * cs_math_pi/180.) * u_norm;
    cpro_met_vel[cell_id][1] = - cos(angle * cs_math_pi/180.) * u_norm;

    /* Potential temperature profile
     * Note: same roughness as dynamics */
    cpro_met_potemp[cell_id] = theta0 + tstar / kappa * cs_mo_psih(z+z0, z0, dlmo);

    /* Richardson flux number profile */
    // Note : ri_f = z/(Pr_t L) * phih/phim^2 = z/Lmo / phim
    cs_real_t ri_f = (z+z0) * dlmo / cs_mo_phim(z+z0, dlmo);

    /* TKE profile */
    cpro_met_k[cell_id] = cs_math_pow2(ustar0) / sqrt(cmu)
      * sqrt(1. - CS_MIN(ri_f, 1.));
    if (f_axz != NULL)
      f_axz->val[cell_id] = -sqrt(cmu / (1. - CS_MIN(ri_f, ri_max)));

    /* epsilon profile */
    cpro_met_eps[cell_id] = cs_math_pow3(ustar0) / (kappa * (z+z0))
       * cs_mo_phim(z+z0, dlmo)*(1.-CS_MIN(ri_f, 1.));

    /* Very stable cases */
    if (ri_f > ri_max) {
      if (z < z_lim) {
        //Ri_f is an increasing monotonic function, so the lowest value of
        //z for which Ri_f>Ri_max is needed
        z_lim = z;
        u_met_min=u_norm;
        theta_met_min=cpro_met_potemp[cell_id];
      }
    }
  }

  cs_parall_min(1,CS_REAL_TYPE, &z_lim);
  cs_parall_min(1,CS_REAL_TYPE, &u_met_min);
  cs_parall_min(1,CS_REAL_TYPE, &theta_met_min);

  /* Very stable cases, corresponding to mode 0 in the Python prepro */
  if (z_lim < 0.5*cs_math_big_r) { // Clipping only if there are cells to be clipped
    bft_printf("Switching to very stable clipping for meteo profile.\n");
    bft_printf("All altitudes above %f have been modified by clipping.\n",z_lim);
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

      cs_real_t z_grd = 0.;
      if (z_ground != NULL)
        z_grd = z_ground[cell_id];

      cs_real_t z = cell_cen[cell_id][2] - z_grd;
      if (z >= z_lim) {
         /* mode = 0 is ustar=cst */
        dlmo_var[cell_id] = dlmo * (z_lim + z0) / (z + z0);

        /* Velocity profile */
        cs_real_t u_norm =   u_met_min + ustar0 / kappa
                           * cs_mo_phim(z_lim + z0, dlmo)
                           * log((z+z0) / (z_lim+z0));

        cpro_met_vel[cell_id][0] = - sin(angle * cs_math_pi/180.) * u_norm;
        cpro_met_vel[cell_id][1] = - cos(angle * cs_math_pi/180.) * u_norm;

       /* Potential temperature profile
        * Note: same roughness as dynamics */
        cpro_met_potemp[cell_id] =   theta_met_min
                                   + tstar * (z_lim+z0) / kappa
                                     * cs_mo_phih(z_lim+z0, dlmo)
                                     * (-1./(z+z0) + 1./(z_lim+z0)) ;
       /* TKE profile
          ri_max is necessarily lower than 1, but CS_MIN might be useful if
          that changes in the future */
        cpro_met_k[cell_id] = cs_math_pow2(ustar0) / sqrt(cmu)
          * sqrt(1. - CS_MIN(ri_max, 1.));

        if (f_axz != NULL)
          f_axz->val[cell_id] = -sqrt(cmu / (1. - CS_MIN(ri_max, 1.)));

        /* epsilon profile */
        cpro_met_eps[cell_id] = cs_math_pow3(ustar0) / kappa  * dlmo_var[cell_id]
         * (1- CS_MIN(ri_max, 1.)) / CS_MIN(ri_max, 1.);
      }
    }
  }

  cs_atmo_hydrostatic_profiles_compute();
  BFT_FREE(dlmo_var);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the ground elevation.
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

  cs_equation_param_t *eqp_p = cs_field_get_equation_param(f);

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

      eqp_p->ndircl = 1;
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

  cs_parall_max(1, CS_INT_TYPE, &(eqp_p->ndircl));

  /* Matrix
   *=======*/

  cs_real_t *rovsdt, *dpvar;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++) {
    rovsdt[cell_id] = 0.;
    dpvar[cell_id] = 0.;
  }
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
  cs_real_t inf_norm = 1.;

  /* Overall loop in order to ensure convergence */
  for (int sweep = 0; sweep < eqp_p->nswrsm && inf_norm > eqp_p->epsrsm; sweep++) {

    cs_equation_iterative_solve_scalar(0,   /* idtvar: no steady state algo */
                                       -1,  /* no over loops */
                                       f->id,
                                       f->name,
                                       0,   /* iescap */
                                       0,   /* imucpp */
                                       norm,
                                       eqp_p,
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

    /* Compute the L_infinity norm */
    inf_norm = 0.;
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
      //FIXME make this dimensionless
      inf_norm = fmax(inf_norm, fabs(f->val[cell_id]-f->val_pre[cell_id]));

      /* Current to previous */
      f->val_pre[cell_id] = f->val[cell_id];
    }

    cs_parall_max(1, CS_REAL_TYPE, &inf_norm);
  }

  /* Free memory */
  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute hydrostatic profiles of density and pressure.
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \divs \left( \grad \varia \right)
 *      = \divs \left( \dfrac{\vect{g}}{c_p \theta} \right)
 *  \f]
 *  where \f$ \vect{g} \f$ is the gravity field and \f$ \theta \f$
 *  is the potential temperature.
 *
 *  The boundary conditions on \f$ \varia \f$ read:
 *  \f[
 *   \varia = \left(\dfrac{P_{sea}}{p_s}\right)^{R/C_p} \textrm{on the ground}
 *  \f]
 *  and Neumann elsewhere.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_hydrostatic_profiles_compute(void)
{
   /* Initialization
   *===============*/

  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;

  cs_physical_constants_t *phys_cst = cs_get_glob_physical_constants();
  cs_field_t *f = cs_field_by_name("meteo_pressure");
  cs_field_t *potemp = cs_field_by_name("meteo_pot_temperature");
  cs_field_t *density = cs_field_by_name("meteo_density");
  cs_field_t *temp = cs_field_by_name("meteo_temperature");
  cs_equation_param_t *eqp_p = cs_field_get_equation_param(f);
  cs_atmo_option_t *aopt = &_atmo_option;
  cs_real_t g = cs_math_3_norm(phys_cst->gravity);

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  /* potential temp at ref */
  cs_real_t pref = cs_glob_atmo_constants->ps;
  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t cp0 = phys_pro->cp0;
  cs_real_t rscp = rair/cp0; /* Around 2/7 */

  const cs_velocity_pressure_model_t  *vp_model
    = cs_glob_velocity_pressure_model;
  const int idilat = vp_model->idilat;

  /* Get the lowest altitude (should also be minimum of z_ground)
   *=============================================================*/

  cs_real_t z_min = cs_math_big_r;

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id ++) {
    if (b_face_cog[face_id][2] < z_min)
      z_min = b_face_cog[face_id][2];
  }

  cs_parall_min(1, CS_REAL_TYPE, &z_min);

  /* p_ground is pressure at the lowest level */

  cs_real_t p_ground = aopt->meteo_psea;

  /* Check if restart is read and if meteo_pressure is present */
  int has_restart = cs_restart_present();
  if (has_restart == 1) {
    cs_restart_t *rp = cs_restart_create("main.csc",
                                         NULL,
                                         CS_RESTART_MODE_READ);

    int retval = cs_restart_read_field_vals(rp,
                                            f->id, /* meteo_pressure */
                                            0);    /* current value */
    if (retval != CS_RESTART_SUCCESS)
      has_restart = 0;

  }

  /* Initialize temperature, pressure and density from neutral conditions
   *=====================================================================*/
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    cs_real_t z  = cell_cen[cell_id][2] - z_min;
    cs_real_t zt = fmin(z, 11000.);
    cs_real_t factor = fmax(1. - g * zt / (cp0 * aopt->meteo_t0), 0.);
    temp->val[cell_id] = aopt->meteo_t0 * factor;

    /* Do not overwrite pressure in case of restart */
    if (has_restart == 0)
      f->val[cell_id] = p_ground * pow(factor, rscp)
                      /* correction factor for z > 11000m */
                      * exp(- g/(rair*temp->val[cell_id]) * (z - zt));

    if (idilat > 0)
      temp->val[cell_id] = potemp->val[cell_id]
                         * pow((f->val[cell_id]/pref), rscp);
    density->val[cell_id] = f->val[cell_id] / (rair * temp->val[cell_id]);
  }

  if (has_restart == 1)
    return;

  /* Boussinesq hypothesis */
  if (idilat == 0) {
    bft_printf(
        "Meteo profiles are computed according to Boussinesq approximation.\n"
        "Using adiabatic profiles for temperature and pressure."
        "Density is computed accordingly.\n");
  }

  cs_real_t *i_massflux = NULL;
  cs_real_t *b_massflux = NULL;
  BFT_MALLOC(i_massflux, m->n_i_faces, cs_real_t);
  BFT_MALLOC(b_massflux, m->n_b_faces, cs_real_t);

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
    i_massflux[face_id] = 0.;

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++)
    b_massflux[face_id] = 0.;

  cs_real_t *i_viscm = NULL;
  BFT_MALLOC(i_viscm, m->n_i_faces, cs_real_t);

  cs_real_t *b_viscm = NULL;
  BFT_MALLOC(b_viscm, m->n_b_faces, cs_real_t);

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
    i_viscm[face_id] = 0.;

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++)
    b_viscm[face_id] = 0.;

  cs_real_3_t *f_ext, *dfext;
  BFT_MALLOC(f_ext, m->n_cells_with_ghosts, cs_real_3_t);
  BFT_MALLOC(dfext, m->n_cells_with_ghosts, cs_real_3_t);

  /* dfext is actually a dummy used to copy calhyd */
  /* f_ext is initialized with an initial density */

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    f_ext[cell_id][0] = density->val[cell_id] * phys_cst->gravity[0];
    dfext[cell_id][0] = 0.;
    f_ext[cell_id][1] = density->val[cell_id] * phys_cst->gravity[1];
    dfext[cell_id][1] = 0.;
    f_ext[cell_id][2] = density->val[cell_id] * phys_cst->gravity[2];
    dfext[cell_id][2] = 0;
  }

  /* Solving
  *=========*/

  cs_real_t *dam = NULL;
  BFT_MALLOC(dam, m->n_cells_with_ghosts, cs_real_t);

  cs_real_t *xam = NULL;
  BFT_MALLOC(xam, m->n_i_faces, cs_real_t);

  cs_real_t *rhs = NULL;
  BFT_MALLOC(rhs, m->n_cells_with_ghosts, cs_real_t);

  cs_real_t *dpvar = NULL;
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    dpvar[cell_id] = 0.;
    dam[cell_id] = 0.;
    rhs[cell_id] = 0.;
  }

  cs_real_t inf_norm = 1.;

  /* Loop to compute pressure profile */
  for (int sweep = 0; sweep < eqp_p->nswrsm && inf_norm > eqp_p->epsrsm; sweep++) {
    //FIXME 100 or nswrsm

    /* Update previous values of pressure for the convergence test */
    cs_field_current_to_previous(f);

    _hydrostatic_pressure_compute(f_ext,
                                  dfext,
                                  f->val,
                                  f->name,
                                  z_min,
                                  p_ground,
                                  i_massflux,
                                  b_massflux,
                                  i_viscm,
                                  b_viscm,
                                  dam,
                                  xam,
                                  dpvar,
                                  rhs);

    /* L infinity residual computation and forcing update */
    inf_norm = 0.;
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
      inf_norm = fmax(fabs(f->val[cell_id] - f->val_pre[cell_id])/pref, inf_norm);

      /* Boussinesq hypothesis: do not update adiabatic temperature profile */
      if (idilat > 0) {
        temp->val[cell_id] = potemp->val[cell_id]
                           * pow((f->val[cell_id]/pref), rscp);
      }

      /* f_ext = rho^k * g */
      cs_real_t rho_k = f->val[cell_id] / (rair * temp->val[cell_id]);
      density->val[cell_id] = rho_k;

      f_ext[cell_id][0] = rho_k * phys_cst->gravity[0];
      f_ext[cell_id][1] = rho_k * phys_cst->gravity[1];
      f_ext[cell_id][2] = rho_k * phys_cst->gravity[2];
    }
    cs_parall_max(1, CS_REAL_TYPE ,&inf_norm);

    if (cs_log_default_is_active())
      bft_printf
        (_("Meteo profiles: iterative process to compute hydrostatic pressure\n"
           "  sweep %d, L infinity norm (delta p) / ps =%e\n"), sweep, inf_norm);
  }

  /* Free memory */
  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(i_viscm);
  BFT_FREE(b_viscm);
  BFT_FREE(xam);
  BFT_FREE(dam);
  BFT_FREE(f_ext);
  BFT_FREE(dfext);
  BFT_FREE(b_massflux);
  BFT_FREE(i_massflux);
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

  if (_atmo_option.meteo_file_name == NULL) {
    BFT_MALLOC(_atmo_option.meteo_file_name,
               strlen(file_name) + 1,
               char);
  }
  else {
    BFT_REALLOC(_atmo_option.meteo_file_name,
                strlen(file_name) + 1,
                char);
  }

  sprintf(_atmo_option.meteo_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the chemistry file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_chem_conc_file_name(const char *file_name)
{
  if (file_name == NULL) {
    return;
  }

  if (_atmo_chem.chem_conc_file_name == NULL) {
    BFT_MALLOC(_atmo_chem.chem_conc_file_name,
               strlen(file_name) + 1,
               char);
  }
  else {
    BFT_REALLOC(_atmo_chem.chem_conc_file_name,
                strlen(file_name) + 1,
                char);
  }

  sprintf(_atmo_chem.chem_conc_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the aerosol file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_aero_conc_file_name(const char *file_name)
{
  if (file_name == NULL) {
    return;
  }
  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_OFF) {
    return;
  }

  if (_atmo_chem.aero_conc_file_name == NULL) {
    BFT_MALLOC(_atmo_chem.aero_conc_file_name,
               strlen(file_name) + 1,
               char);
  }
  else {
    BFT_REALLOC(_atmo_chem.aero_conc_file_name,
                strlen(file_name) + 1,
                char);
  }

  sprintf(_atmo_chem.aero_conc_file_name, "%s", file_name);
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
       "    longitude: %6f\n"
       "    x center (in Lambert-93) : %6f\n"
       "    y center (in Lmabert-93) : %6f\n\n"),
     cs_glob_atmo_option->latitude,
     cs_glob_atmo_option->longitude,
     cs_glob_atmo_option->x_l93,
     cs_glob_atmo_option->y_l93);

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
/*!
 * \brief Deallocate arrays for atmo module
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_finalize(void)
{
  BFT_FREE(_atmo_option.meteo_file_name);
  BFT_FREE(_atmo_option.z_dyn_met);
  BFT_FREE(_atmo_option.z_temp_met);
  BFT_FREE(_atmo_option.u_met);
  BFT_FREE(_atmo_option.v_met);
  BFT_FREE(_atmo_option.time_met);
  BFT_FREE(_atmo_option.hyd_p_met);
  BFT_FREE(_atmo_option.pot_t_met);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
