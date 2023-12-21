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

#include "cs_air_props.h"
#include "cs_array.h"
#include "cs_atmo_profile_std.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
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
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
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
#include "cs_intprf.h"

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
  .met_1d_nlevels_d = 0,
  .met_1d_nlevels_t = 0,
  .met_1d_ntimes = 0,
  .met_1d_nlevels_max_t = 0,
  .radiative_model_1d = 0,
  .rad_1d_nvert = 1,
  .rad_1d_nlevels = 20,
  .rad_1d_nlevels_max = 0,
  .rad_1d_xy = NULL,
  .rad_1d_z = NULL,
  .rad_1d_acinfe = NULL,
  .rad_1d_dacinfe = NULL,
  .rad_1d_aco2 = NULL,
  .rad_1d_aco2s = NULL,
  .rad_1d_daco2 = NULL,
  .rad_1d_daco2s = NULL,
  .rad_1d_acsup = NULL,
  .rad_1d_acsups = NULL,
  .rad_1d_dacsup = NULL,
  .rad_1d_dacsups = NULL,
  .rad_1d_tauzq = NULL,
  .rad_1d_tauz = NULL,
  .rad_1d_zq = NULL,
  .rad_1d_zray = NULL,
  .rad_1d_ir_div = NULL,
  .rad_1d_sol_div = NULL,
  .rad_1d_iru = NULL,
  .rad_1d_ird = NULL,
  .rad_1d_solu = NULL,
  .rad_1d_sold = NULL,
  .rad_1d_qw = NULL,
  .rad_1d_ql = NULL,
  .rad_1d_qv = NULL,
  .rad_1d_nc = NULL,
  .rad_1d_fn = NULL,
  .rad_1d_aerosols = NULL,
  .domain_orientation = 0.,
  .compute_z_ground = false,
  .open_bcs_treatment = 0,
  .theo_interp = 0,
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
  .meteo_qw0 = 0.,
  .meteo_qwstar = DBL_MAX,
  .meteo_qw1 = DBL_MAX,
  .meteo_qw2 = DBL_MAX,
  .meteo_ql0 = 0.,
  .meteo_evapor = DBL_MAX,
  .meteo_sensi = DBL_MAX,
  .meteo_psea = 101325.,
  .meteo_phim_s = 0, /* Cheng 2005 by default */
  .meteo_phih_s = 0, /* Cheng 2005 by default */
  .meteo_phim_u = 1, /* Hogstrom 1988 by default */
  .meteo_phih_u = 1, /* Hogstrom 1988 by default */
  .u_met      = NULL,
  .v_met      = NULL,
  .w_met      = NULL,
  .ek_met     = NULL,
  .ep_met     = NULL,
  .z_dyn_met  = NULL,
  .z_temp_met = NULL,
  .time_met   = NULL,
  .hyd_p_met  = NULL,
  .pot_t_met  = NULL,
  .soil_model = 0, /* off or user defined */
  .soil_cat = 0, /* CS_ATMO_SOIL_5_CAT */
  .soil_zone_id = -1,
  .soil_meb_model = 0
};

static const char *_univ_fn_name[] = {N_("Cheng 2005"),
                                      N_("Hogstrom 1988"),
                                      N_("Businger 1971"),
                                      N_("Hartogensis 2007"),
                                      N_("Carl 1973")};

/* global atmo constants structure */
static cs_atmo_constants_t _atmo_constants = {
  .ps = 1.e5
};

/* atmo chemistry options structure */
static cs_atmo_chemistry_t _atmo_chem = {
  .model = 0,
  .n_species = 0,
  .n_reactions = 0,
  .chemistry_sep_mode = 2,
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
  .reacnum = NULL,
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
                       int                    **theo_interp,
                       int                    **sedimentation_model,
                       int                    **deposition_model,
                       int                    **nucleation_model,
                       int                    **subgrid_model,
                       int                    **distribution_model,
                       int                    **model,
                       int                    **isepchemistry,
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
                       int                    **radiative_model_1d,
                       int                    **nbmaxt,
                       cs_real_t              **meteo_zi,
                       int                    **soil_model,
                       int                    **nvert,
                       int                    **kvert,
                       int                    **kmx);

void
cs_f_atmo_arrays_get_pointers(cs_real_t **z_dyn_met,
                              cs_real_t **z_temp_met,
                              cs_real_t **u_met,
                              cs_real_t **v_met,
                              cs_real_t **w_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              cs_real_t **pot_t_met,
                              cs_real_t **ek_met,
                              cs_real_t **ep_met,
                              cs_real_t **xyvert,
                              cs_real_t **zvert,
                              cs_real_t **acinfe,
                              cs_real_t **dacinfe,
                              cs_real_t **aco2,
                              cs_real_t **aco2s,
                              cs_real_t **daco2,
                              cs_real_t **daco2s,
                              cs_real_t **acsup,
                              cs_real_t **acsups,
                              cs_real_t **dacsup,
                              cs_real_t **dacsups,
                              cs_real_t **tauzq,
                              cs_real_t **tauz,
                              cs_real_t **zq,
                              cs_real_t **zray,
                              cs_real_t **rayi,
                              cs_real_t **rayst,
                              cs_real_t **iru,
                              cs_real_t **ird,
                              cs_real_t **solu,
                              cs_real_t **sold,
                              int         dim_u_met[2],
                              int         dim_hyd_p_met[2],
                              int         dim_pot_t_met[2],
                              int         dim_ek_met[2],
                              int         dim_ep_met[2],
                              int         dim_xyvert[2],
                              int         dim_kmx2[2],
                              int         dim_kmx_nvert[2]);

void
cs_f_atmo_rad_1d_arrays_get_pointers(cs_real_t **qwvert,
                                     cs_real_t **qlvert,
                                     cs_real_t **qvvert,
                                     cs_real_t **ncvert,
                                     cs_real_t **fnvert,
                                     cs_real_t **aevert);

void
cs_f_atmo_get_soil_zone(cs_lnum_t         *n_elts,
                        int               *n_soil_cat,
                        const cs_lnum_t  **elt_ids);

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint);

void
cs_f_atmo_chem_initialize_reacnum(cs_real_t **reacnum);

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid);

void
cs_f_atmo_chem_finalize(void);

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
  cs_solving_info_t *sinfo = NULL;
  int key_sinfo_id = cs_field_key_id("solving_info");
  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    sinfo = cs_field_get_key_struct_ptr(f, key_sinfo_id);
    sinfo->n_it = 0;
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

  /* Log */
  if (sinfo != NULL)
    sinfo->rhs_norm = residu;

  /* Initial Right-Hand-Side */
  cs_diffusion_potential(f_id,
                         m,
                         mq,
                         1, /* init */
                         1, /* inc */
                         eqp_p->imrgra,
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
                         NULL,
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

    /* Log */
    if (sinfo != NULL)
      sinfo->n_it += niterf;

    /* Update variable and right-hand-side */
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
      pvar[cell_id] += dpvar[cell_id];

    cs_diffusion_potential(f_id,
                           m,
                           mq,
                           1, /* init */
                           1, /* inc */
                           eqp_p->imrgra,
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

  /* For log */
  if (sinfo != NULL) {
    if (rnorm > 0.)
      sinfo->res_norm = residu/rnorm;
    else
      sinfo->res_norm = 0.;
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
                                    + 3.*cs_math_pi/180.)/cs_math_pi*180.;
  cs_glob_atmo_option->latitude
    = asin(tanh((log(c/sqrt(  cs_math_pow2(cs_glob_atmo_option->x_l93-xs)
                            + cs_math_pow2(cs_glob_atmo_option->y_l93-ys)))/n)
                +e*atanh(e*tanh(t3))))/cs_math_pi*180.;
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions, for neutral
 *        (derivative function Phi_m and Phi_h)
 *
 * \param[in]  z             altitude
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_phim_n(cs_real_t              z,
           cs_real_t              dlmo)
{
  CS_UNUSED(z);
  CS_UNUSED(dlmo);

  return 1.0;
}

static cs_real_t
_mo_phih_n(cs_real_t              z,
           cs_real_t              dlmo)
{
  CS_UNUSED(z);
  CS_UNUSED(dlmo);

  return 1.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions, for neutral
 *        (Integrated version from z0 to z)
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_psim_n(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  CS_UNUSED(dlmo);

  return log(z/z0);
}

static cs_real_t
_mo_psih_n(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  CS_UNUSED(dlmo);

  return log(z/z0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions for stable
 *        (derivative function Phi_m and Phi_h)
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_phim_s(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;

  switch(cs_glob_atmo_option->meteo_phim_s) {

  case CS_ATMO_UNIV_FN_CHENG:
    {
      cs_real_t a = 6.1;
      cs_real_t b = 2.5;

      return 1. + a*(x+(pow(x,b))*( pow(1.+pow(x,b),(1.-b)/b) ))
             / (x+ pow(1.+ pow(x,b),1./b));
    }

  case CS_ATMO_UNIV_FN_HOGSTROM:
    {
      if (x < 0.5) {
        cs_real_t b = 4.8;

        return 1. + b*x;
      }
      else if (x < 10.) {
        cs_real_t a = 7.9;

        return a - 4.25/x + 1./pow(x,2.);
      }
      else {
        cs_real_t a = 0.7485;

        return a*x;
      }
    }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      if (x < 0.5) {
        cs_real_t b = 4.7;

        return 1. + b*x;
      }
      else if (x < 10.) {
        cs_real_t a = 7.85;

        return a - 4.25/x + 1./pow(x,2.);
      }
      else {
        cs_real_t a = 0.7435;

        return a*x;
      }
    }

  case CS_ATMO_UNIV_FN_HARTOGENSIS:
    {
      cs_real_t a = 1.;
      cs_real_t b = 2./3.;
      cs_real_t c = 5.;
      cs_real_t d = 0.35;

      return 1. + x*(a + b*exp(-d*x) - b*d*(x - c/d)*exp(-d*x));
    }

  default:
    assert(0);
    return -1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions for unstable
 *        (derivative function)
 *
 * \param[in]  z             altitude
 * \param[in]  dlmo          inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_phih_s(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;

  switch(cs_glob_atmo_option->meteo_phih_s) {

    case CS_ATMO_UNIV_FN_CHENG:
      {
        cs_real_t a = 5.3;
        cs_real_t b = 1.1;

        return 1.+a*(x+(pow(x,b))*(pow((1.+pow(x,b)), ((1.-b)/b))))
          / (x + pow((1.+pow(x,b)),1./b));
      }

  case CS_ATMO_UNIV_FN_HOGSTROM:
    {
      cs_real_t a = 0.95;
      cs_real_t b = 7.8;

      return a + b*x;
    }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      cs_real_t a = 0.74;
      cs_real_t b = 4.7;

      return a + b*x;
    }

  case CS_ATMO_UNIV_FN_HARTOGENSIS:
    {
      cs_real_t a = 1.;
      cs_real_t b = 2./3.;
      cs_real_t c = 5.;
      cs_real_t d = 0.35;

      return 1. + x*(a*pow((1. + 2./3. * a * x),0.5)
             + b*exp(-d*x) - b*d*(x - c/d)*exp(-d*x));
    }

  default:
    assert(0);
    return -1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions for unstable
 *        (derivative function)
 *
 * \param[in]  z             altitude
 * \param[in]  dlmo          inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_phim_u(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;

  switch (cs_glob_atmo_option->meteo_phim_u) {

  case  CS_ATMO_UNIV_FN_HOGSTROM:
    {
      cs_real_t a = 1.;
      cs_real_t b = 19.3;
      cs_real_t e = -0.25;

      return a*pow((1.-b*x),e);
    }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      cs_real_t a = 1.;
      cs_real_t b = 15.;
      cs_real_t e = -0.25;

      return a*pow((1.-b*x),e);
    }

  case CS_ATMO_UNIV_FN_CARL:
    {
      cs_real_t a = 1.;
      cs_real_t b = 16.;
      cs_real_t e = -1./3.;

      return a*pow((1.-b*x),e);
    }

  default:
    assert(0);
    return -1;
  }
}

static cs_real_t
_mo_phih_u(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;

  switch(cs_glob_atmo_option->meteo_phih_u) {

  case CS_ATMO_UNIV_FN_HOGSTROM:
    {
      cs_real_t a = 0.95;
      cs_real_t b = 11.6;
      cs_real_t e = -0.5;

      return a*pow(1.-b*x, e);
    }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      cs_real_t a = 0.74;
      cs_real_t b = 9.;
      cs_real_t e = -0.5;

      return a*pow(1.-b*x, e);
    }

  case CS_ATMO_UNIV_FN_CARL:
    {
      cs_real_t a = 0.74;
      cs_real_t b = 16.;
      cs_real_t e = -0.5;

      return a*pow(1.-b*x, e);
    }

  default:
    assert(0);
    return -1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions for stable
 *        (integral functions Psi_m and Psi_h, integrated version from z0 to z)
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_psim_s(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;
  cs_real_t x0 = z0 * dlmo;

  switch(cs_glob_atmo_option->meteo_phim_s) {

  case CS_ATMO_UNIV_FN_CHENG:
    {
      cs_real_t a = 6.1;
      cs_real_t b = 2.5;

      return log(z/z0) + a*log(x + pow((1. + pow(x,b)),1./b))
         - a*log(x0 + pow((1. + pow(x0,b)),1./b));
    }

    case CS_ATMO_UNIV_FN_HOGSTROM:
      {
        if (x < 0.5) {
          cs_real_t b = 4.8;

          return log(z/z0) + b*(x - x0);
        }
        else if (x < 10.) {
          cs_real_t a = 7.9;
          cs_real_t b = 4.8;
          cs_real_t c = 4.1;

          return a*log(2.*x) + 4.25/x - 0.5/pow(x,2.) - log(2.*x0) - b*x0 - c;
        }
        else {
          cs_real_t a = 0.7485;
          cs_real_t b = 7.9;
          cs_real_t c = 4.8;

          return a*x + b*log(20.) - 11.165 - log(2.*x0) - c*x0;
        }
      }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      if (x < 0.5) {
        cs_real_t b = 4.7;

        return log(z/z0) + b*(x - x0);
      }
      else if (x < 10.) {
        cs_real_t a = 7.85;
        cs_real_t b = 4.7;
        cs_real_t c = 4.15;

        return a*log(2.*x) + 4.25/x - 0.5/pow(x,2.) - log(2.*x0) - b*x0 - c;
      }
      else {
        cs_real_t a = 0.7435;
        cs_real_t b = 7.85;
        cs_real_t c = 4.7;

        return a*x + b*log(20.) - 11.165 - log(2.*x0) - c*x0;
      }
    }

  case CS_ATMO_UNIV_FN_HARTOGENSIS:
    {
      cs_real_t a = 1.;
      cs_real_t b = 2./3.;
      cs_real_t c = 5.;
      cs_real_t d = 0.35;

      return log(z/z0) + a*(x - x0)
             + b*(x - c/d)*exp(-d*x) - b*(x0 - c/d)*exp(-d*x0);
    }

  default:
    assert(0);
    return -1;
  }
}

static cs_real_t
_mo_psih_s(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  cs_real_t x = z * dlmo;
  cs_real_t x0 = z0 * dlmo;

  switch(cs_glob_atmo_option->meteo_phih_s) {

  case CS_ATMO_UNIV_FN_CHENG:
    {
      cs_real_t a = 5.3;
      cs_real_t b = 1.1;

      return log(z/z0) + a*log(x + pow((1. + pow(x,b)),1./b))
        - a*log(x0 + pow((1. + pow(x0,b)),1./b));
    }

  case CS_ATMO_UNIV_FN_HOGSTROM:
    {
      cs_real_t a = 0.95;
      cs_real_t b = 7.8;

      return a*log(z/z0) + b*(x - x0);
    }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      cs_real_t a = 0.74;
      cs_real_t b = 4.7;

      return a*log(z/z0) + b*(x - x0);
    }

  case CS_ATMO_UNIV_FN_HARTOGENSIS:
    {
      cs_real_t a = 1.;
      cs_real_t b = 2./3.;
      cs_real_t c = 5.;
      cs_real_t d = 0.35;

      return log(z/z0) + pow((1. + 2./3. * a * x),3./2.)
        + b*(x - c/d)*exp(-d*x) - pow((1. + 2./3. * a * x0),3./2.)
        - b*(x0 - c/d)*exp(-d*x0);
    }

  default:
    assert(0);
    return -1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal functions for unstable
 *        (integral functions Psi_m and Psi_h, integrated version from z0 to z)
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          inverse Monin Obukhov length
 *
 * \return coef
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_mo_psim_u(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  switch (cs_glob_atmo_option->meteo_phim_u) {

    case CS_ATMO_UNIV_FN_HOGSTROM:
      {
        cs_real_t b = 19.3;
        cs_real_t e = 0.25;
        cs_real_t x = pow((1. - b*z*dlmo), e);
        cs_real_t x0 = pow((1. - b*z0*dlmo), e);

        return log(z/z0) - 2.*log((1. + x)/(1. + x0))
          - log((1. + pow(x,2.))/(1. + pow(x0,2.))) + 2.*atan(x) - 2.*atan(x0);
      }

  case CS_ATMO_UNIV_FN_BUSINGER:
    {
      cs_real_t b = 15.;
      cs_real_t e = 0.25;
      cs_real_t x = pow((1. - b*z*dlmo), e);
      cs_real_t x0 = pow((1. - b*z0*dlmo), e);

      return log(z/z0) - 2.*log((1. + x)/(1. + x0))
        - log((1. + pow(x,2.))/(1. + pow(x0,2.))) + 2.*atan(x) - 2.*atan(x0);
    }

  case CS_ATMO_UNIV_FN_CARL:
    {
      cs_real_t b = 16.;
      cs_real_t e = 1./3.;
      cs_real_t x = pow((1. - b*z*dlmo), e);
      cs_real_t x0 = pow((1. - b*z0*dlmo), e);

      return log(z/z0) - 1.5*log((1. + x + pow(x,2.))/3.)
        + pow(3., 0.5)*atan((1. + 2.*x)/pow(3., 0.5))
        + 1.5*log((1. + x0 + pow(x0,2.))/3.)
        - pow(3., 0.5)*atan((1. + 2.*x0)/pow(3., 0.5));
    }

  default:
    assert(0);
    return -1;
  }
}

static cs_real_t
_mo_psih_u(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
   switch (cs_glob_atmo_option->meteo_phih_u) {

     case CS_ATMO_UNIV_FN_HOGSTROM:
       {
         cs_real_t a = 0.95;
         cs_real_t b = 11.6;
         cs_real_t e = 0.5;
         cs_real_t x = pow((1. - b*z*dlmo), e);
         cs_real_t x0 = pow((1. - b*z0*dlmo), e);

         return a*(log(z/z0) - 2.*log((1. + x)/(1. + x0)));
       }

   case  CS_ATMO_UNIV_FN_BUSINGER:
     {
       cs_real_t a = 0.74;
       cs_real_t b = 9.;
       cs_real_t e = 0.5;
       cs_real_t x = pow((1. - b*z*dlmo), e);
       cs_real_t x0 = pow((1. - b*z0*dlmo), e);

       return a*(log(z/z0) - 2.*log((1. + x)/(1. + x0)));
     }

   case  CS_ATMO_UNIV_FN_CARL:
     {
       cs_real_t a = 0.74;
       cs_real_t b = 16.;
       cs_real_t e = 0.5;
       cs_real_t x = pow((1. - b*z*dlmo), e);
       cs_real_t x0 = pow((1. - b*z0*dlmo), e);

       return a*(log(z/z0) - 2.*log((1. + x)/(1. + x0)));
     }

   default:
     assert(0);
     return -1;
   }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal function phim for neutral, stable and unstable
 *
 * \param[in]  z             altitude
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return                   factor
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_mo_phim(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t dlmoneutral = 1.e-12;
  cs_real_t coef;

  if (CS_ABS(dlmo) < dlmoneutral)
    coef = _mo_phim_n(z,dlmo);
  else if (dlmo >= 0.)
    coef = _mo_phim_s(z,dlmo);
  else
    coef = _mo_phim_u(z,dlmo);

  return coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal function phih for neutral, stable and unstable
 *
 * \param[in]  z             altitude
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return                   factor
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_mo_phih(cs_real_t              z,
           cs_real_t              dlmo)
{
  cs_real_t dlmoneutral = 1.e-12;
  cs_real_t coef;

  if (CS_ABS(dlmo) < dlmoneutral)
    coef = _mo_phih_n(z,dlmo);
  else if (dlmo >= 0.)
    coef = _mo_phih_s(z,dlmo);
  else
    coef = _mo_phih_u(z,dlmo);

  return coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal function psim for neutral, stable and unstable
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return                   factor
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_mo_psim(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  cs_real_t dlmoneutral = 1.e-12;
  cs_real_t coef;

  if (CS_ABS(dlmo) < dlmoneutral)
    coef = _mo_psim_n(z, z0, dlmo);
  else if (dlmo >= 0.)
    coef = _mo_psim_s(z, z0, dlmo);
  else
    coef = _mo_psim_u(z, z0, dlmo);

  return coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Universal function psih for neutral, stable and unstable
 *
 * \param[in]  z             altitude
 * \param[in]  z0            altitude of the starting point integration
 * \param[in]  dlmo          Inverse Monin Obukhov length
 *
 * \return                   factor
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_mo_psih(cs_real_t              z,
           cs_real_t              z0,
           cs_real_t              dlmo)
{
  cs_real_t dlmoneutral = 1.e-12;
  cs_real_t coef;

  if (CS_ABS(dlmo) < dlmoneutral)
    coef = _mo_psih_n(z, z0, dlmo);
  else if (dlmo >= 0.)
    coef = _mo_psih_s(z, z0, dlmo);
  else
    coef = _mo_psih_u(z, z0, dlmo);

  return coef;
}

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
                       int                    **theo_interp,
                       int                    **sedimentation_model,
                       int                    **deposition_model,
                       int                    **nucleation_model,
                       int                    **subgrid_model,
                       int                    **distribution_model,
                       int                    **model,
                       int                    **isepchemistry,
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
                       int                    **radiative_model_1d,
                       int                    **nbmaxt,
                       cs_real_t              **meteo_zi,
                       int                    **soil_model,
                       int                    **nvert,
                       int                    **kvert,
                       int                    **kmx)
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
  *theo_interp = &(_atmo_option.theo_interp);
  *sedimentation_model = &(_atmo_option.sedimentation_model);
  *deposition_model = &(_atmo_option.deposition_model);
  *nucleation_model = &(_atmo_option.nucleation_model);
  *subgrid_model = &(_atmo_option.subgrid_model);
  *distribution_model = &(_atmo_option.distribution_model);
  *meteo_profile = &(_atmo_option.meteo_profile);
  *nbmetd     = &(_atmo_option.met_1d_nlevels_d);
  *nbmett     = &(_atmo_option.met_1d_nlevels_t);
  *nbmetm     = &(_atmo_option.met_1d_ntimes);
  *nbmaxt     = &(_atmo_option.met_1d_nlevels_max_t);
  *meteo_zi   = &(_atmo_option.meteo_zi);
  *soil_model = &(_atmo_option.soil_model);
  *model = &(_atmo_chem.model);
  *isepchemistry = &(_atmo_chem.chemistry_sep_mode);
  *n_species = &(_atmo_chem.n_species);
  *n_reactions = &(_atmo_chem.n_reactions);
  *chemistry_with_photolysis = &(_atmo_chem.chemistry_with_photolysis);
  *aerosol_model = &(_atmo_chem.aerosol_model);
  *frozen_gas_chem = &(_atmo_chem.frozen_gas_chem);
  *init_gas_with_lib = &(_atmo_chem.init_gas_with_lib);
  *init_aero_with_lib = &(_atmo_chem.init_aero_with_lib);
  *n_layer = &(_atmo_chem.n_layer);
  *n_size = &(_atmo_chem.n_size);
  *radiative_model_1d = &(_atmo_option.radiative_model_1d);
  *nvert = &(_atmo_option.rad_1d_nvert);
  *kvert = &(_atmo_option.rad_1d_nlevels);
  *kmx = &(_atmo_option.rad_1d_nlevels_max);
}

void
cs_f_atmo_get_soil_zone(cs_lnum_t         *n_elts,
                        int               *n_soil_cat,
                        const cs_lnum_t  **elt_ids)
{
  *n_elts = 0;
  *elt_ids = NULL;

  /* Not defined */
  if (cs_glob_atmo_option->soil_zone_id < 0) {
    *n_soil_cat = 0;
    return;
  }

  const cs_zone_t *z = cs_boundary_zone_by_id(cs_glob_atmo_option->soil_zone_id);
  *elt_ids = z->elt_ids;
  *n_elts = z->n_elts;
  switch (cs_glob_atmo_option->soil_cat) {
    case CS_ATMO_SOIL_5_CAT:
      *n_soil_cat = 5;
      break;
    case CS_ATMO_SOIL_7_CAT:
      *n_soil_cat = 7;
      break;
    case CS_ATMO_SOIL_23_CAT:
      *n_soil_cat = 23;
      break;
    default:
      *n_soil_cat = 0;
      break;
  }

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
                              cs_real_t **w_met,
                              cs_real_t **time_met,
                              cs_real_t **hyd_p_met,
                              cs_real_t **pot_t_met,
                              cs_real_t **ek_met,
                              cs_real_t **ep_met,
                              cs_real_t **xyvert,
                              cs_real_t **zvert,
                              cs_real_t **acinfe,
                              cs_real_t **dacinfe,
                              cs_real_t **aco2,
                              cs_real_t **aco2s,
                              cs_real_t **daco2,
                              cs_real_t **daco2s,
                              cs_real_t **acsup,
                              cs_real_t **acsups,
                              cs_real_t **dacsup,
                              cs_real_t **dacsups,
                              cs_real_t **tauzq,
                              cs_real_t **tauz,
                              cs_real_t **zq,
                              cs_real_t **zray,
                              cs_real_t **rayi,
                              cs_real_t **rayst,
                              cs_real_t **iru,
                              cs_real_t **ird,
                              cs_real_t **solu,
                              cs_real_t **sold,
                              int         dim_u_met[2],
                              int         dim_hyd_p_met[2],
                              int         dim_pot_t_met[2],
                              int         dim_ek_met[2],
                              int         dim_ep_met[2],
                              int         dim_xyvert[2],
                              int         dim_kmx2[2],
                              int         dim_kmx_nvert[2])
{
  int n_level = 0, n_level_t = 0;
  int n_times = 0;
  if (_atmo_option.meteo_profile) {
    n_level = CS_MAX(1, _atmo_option.met_1d_nlevels_d);
    n_level_t = CS_MAX(1, _atmo_option.met_1d_nlevels_max_t);
    n_times = CS_MAX(1, _atmo_option.met_1d_ntimes);
  }

  if (_atmo_option.z_dyn_met == NULL)
    BFT_MALLOC(_atmo_option.z_dyn_met, n_level, cs_real_t);
  if (_atmo_option.z_temp_met == NULL)
    BFT_MALLOC(_atmo_option.z_temp_met, n_level_t, cs_real_t);
  if (_atmo_option.u_met == NULL)
    BFT_MALLOC(_atmo_option.u_met, n_level*n_times, cs_real_t);
  if (_atmo_option.v_met == NULL)
    BFT_MALLOC(_atmo_option.v_met, n_level*n_times, cs_real_t);
  if (_atmo_option.w_met == NULL)
    BFT_MALLOC(_atmo_option.w_met, n_level*n_times, cs_real_t);
  if (_atmo_option.time_met == NULL)
    BFT_MALLOC(_atmo_option.time_met, n_times, cs_real_t);
  if (_atmo_option.hyd_p_met == NULL)
    BFT_MALLOC(_atmo_option.hyd_p_met,
               n_times * n_level_t, cs_real_t);
  if (_atmo_option.pot_t_met == NULL)
    BFT_MALLOC(_atmo_option.pot_t_met, n_level_t*n_times, cs_real_t);
  if (_atmo_option.ek_met == NULL)
    BFT_MALLOC(_atmo_option.ek_met, n_level*n_times, cs_real_t);
  if (_atmo_option.ep_met == NULL)
    BFT_MALLOC(_atmo_option.ep_met, n_level*n_times, cs_real_t);

  *u_met     = _atmo_option.u_met;
  *v_met     = _atmo_option.v_met;
  *w_met     = _atmo_option.w_met;
  *hyd_p_met = _atmo_option.hyd_p_met;
  *pot_t_met = _atmo_option.pot_t_met;
  *ek_met    = _atmo_option.ek_met;
  *ep_met    = _atmo_option.ep_met;

  *z_dyn_met  = _atmo_option.z_dyn_met;
  *z_temp_met = _atmo_option.z_temp_met;
  *time_met   = _atmo_option.time_met;

  n_level = 0;
  int n_vert = 0;
  if (_atmo_option.radiative_model_1d != 0) {
    n_level = CS_MAX(1, _atmo_option.rad_1d_nlevels_max);
    n_vert = CS_MAX(1, _atmo_option.rad_1d_nvert);
  }

  if (         _atmo_option.rad_1d_xy == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_xy , 3*n_vert, cs_real_t);
  if (         _atmo_option.rad_1d_z == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_z , n_level, cs_real_t);
  if (         _atmo_option.rad_1d_acinfe == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_acinfe , n_level, cs_real_t);
  if (         _atmo_option.rad_1d_dacinfe == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_dacinfe , n_level, cs_real_t);
  if (         _atmo_option.rad_1d_aco2 == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_aco2, n_level*n_level, cs_real_t);
  if (         _atmo_option.rad_1d_aco2s == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_aco2s, n_level*n_level, cs_real_t);
  if (         _atmo_option.rad_1d_daco2 == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_daco2, n_level*n_level, cs_real_t);
  if (         _atmo_option.rad_1d_daco2s == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_daco2s, n_level*n_level, cs_real_t);
  if (         _atmo_option.rad_1d_acsup == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_acsup, n_level, cs_real_t);
  if (         _atmo_option.rad_1d_acsups == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_acsups, n_level, cs_real_t);
  if (         _atmo_option.rad_1d_dacsup == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_dacsup, n_level, cs_real_t);
  if (         _atmo_option.rad_1d_dacsups == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_dacsups, n_level, cs_real_t);
  if (         _atmo_option.rad_1d_tauzq == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_tauzq, n_level+1, cs_real_t);
  if (         _atmo_option.rad_1d_tauz == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_tauz, n_level+1, cs_real_t);
  if (         _atmo_option.rad_1d_zq == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_zq, n_level+1, cs_real_t);
  if (         _atmo_option.rad_1d_zray == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_zray, n_level, cs_real_t);
  if (         _atmo_option.rad_1d_ir_div == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_ir_div, n_level * n_vert, cs_real_t);
  if (         _atmo_option.rad_1d_sol_div == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_sol_div, n_level * n_vert, cs_real_t);
  if (         _atmo_option.rad_1d_iru == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_iru, n_level * n_vert, cs_real_t);
  if (         _atmo_option.rad_1d_ird == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_ird, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_ird);
  }
  if (         _atmo_option.rad_1d_solu == NULL)
    BFT_MALLOC(_atmo_option.rad_1d_solu, n_level * n_vert, cs_real_t);
  if (         _atmo_option.rad_1d_sold == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_sold, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_sold);
  }
  if (         _atmo_option.rad_1d_qw == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_qw, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_qw);
  }
  if (         _atmo_option.rad_1d_ql == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_ql, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_ql);
  }
  if (         _atmo_option.rad_1d_qv == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_qv, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_qv);
  }
  if (         _atmo_option.rad_1d_nc == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_nc, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_nc);
  }
  if (         _atmo_option.rad_1d_fn == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_fn, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_fn);
  }
  if (         _atmo_option.rad_1d_aerosols == NULL) {
    BFT_MALLOC(_atmo_option.rad_1d_aerosols, n_level * n_vert, cs_real_t);
    cs_array_real_fill_zero(n_level * n_vert, _atmo_option.rad_1d_aerosols);
  }

  *xyvert = _atmo_option.rad_1d_xy;
  *zvert  = _atmo_option.rad_1d_z;
  *acinfe = _atmo_option.rad_1d_acinfe;
  *dacinfe = _atmo_option.rad_1d_dacinfe;
  *aco2   = _atmo_option.rad_1d_aco2  ;
  *aco2s  = _atmo_option.rad_1d_aco2s ;
  *daco2  = _atmo_option.rad_1d_daco2 ;
  *daco2s = _atmo_option.rad_1d_daco2s;
  *acsup  = _atmo_option.rad_1d_acsup ;
  *acsups = _atmo_option.rad_1d_acsups;
  *dacsup = _atmo_option.rad_1d_dacsup;
  *dacsups = _atmo_option.rad_1d_dacsups;
  *tauzq  = _atmo_option.rad_1d_tauzq ;
  *tauz   = _atmo_option.rad_1d_tauz  ;
  *zq     = _atmo_option.rad_1d_zq    ;
  *zray   = _atmo_option.rad_1d_zray  ;
  *rayi   = _atmo_option.rad_1d_ir_div  ;
  *rayst  = _atmo_option.rad_1d_sol_div ;
  *iru    = _atmo_option.rad_1d_iru   ;
  *ird    = _atmo_option.rad_1d_ird   ;
  *solu   = _atmo_option.rad_1d_solu  ;
  *sold   = _atmo_option.rad_1d_sold  ;

  dim_u_met[0]     = _atmo_option.met_1d_nlevels_d;
  dim_u_met[1]     = _atmo_option.met_1d_ntimes;
  dim_hyd_p_met[0] = _atmo_option.met_1d_nlevels_max_t;
  dim_hyd_p_met[1] = _atmo_option.met_1d_ntimes;
  dim_pot_t_met[0] = _atmo_option.met_1d_nlevels_max_t;
  dim_pot_t_met[1] = _atmo_option.met_1d_ntimes;
  dim_ek_met[0]    = _atmo_option.met_1d_nlevels_d;
  dim_ek_met[1]    = _atmo_option.met_1d_ntimes;
  dim_ep_met[0]    = _atmo_option.met_1d_nlevels_d;
  dim_ep_met[1]    = _atmo_option.met_1d_ntimes;
  dim_xyvert[0]    = _atmo_option.rad_1d_nvert;
  dim_xyvert[1]    = 3;
  dim_kmx2[0]    = _atmo_option.rad_1d_nlevels_max;
  dim_kmx2[1]    = _atmo_option.rad_1d_nlevels_max;
  dim_kmx_nvert[0] = _atmo_option.rad_1d_nlevels_max;
  dim_kmx_nvert[1] = _atmo_option.rad_1d_nvert;

}

void
cs_f_atmo_rad_1d_arrays_get_pointers(cs_real_t **qwvert,
                                     cs_real_t **qlvert,
                                     cs_real_t **qvvert,
                                     cs_real_t **ncvert,
                                     cs_real_t **fnvert,
                                     cs_real_t **aevert)
{
  *qwvert = _atmo_option.rad_1d_qw;
  *qlvert = _atmo_option.rad_1d_ql;
  *qvvert = _atmo_option.rad_1d_qv;
  *ncvert = _atmo_option.rad_1d_nc;
  *fnvert = _atmo_option.rad_1d_fn;
  *aevert = _atmo_option.rad_1d_aerosols;
}

void
cs_f_atmo_chem_initialize_reacnum(cs_real_t **reacnum)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  if (_atmo_chem.reacnum == NULL)
    BFT_MALLOC(_atmo_chem.reacnum, _atmo_chem.n_reactions*n_cells, cs_real_t);

  *reacnum = _atmo_chem.reacnum;
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

  BFT_FREE(_atmo_chem.reacnum);
  BFT_FREE(_atmo_chem.species_to_scalar_id);
  BFT_FREE(_atmo_chem.species_to_field_id);
  BFT_FREE(_atmo_chem.molar_mass);
  BFT_FREE(_atmo_chem.chempoint);
  BFT_FREE(_atmo_chem.spack_file_name);
  BFT_FREE(_atmo_chem.aero_file_name);
  BFT_FREE(_atmo_chem.chem_conc_file_name);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * add properties fields
 *----------------------------------------------------------------------------*/

void
cs_atmo_add_property_fields(void)
{
  cs_field_t *f;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  int z_id = cs_glob_atmo_option->soil_zone_id;

  if (z_id > -1) {
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);
    int soil_num = 5;
    switch (cs_glob_atmo_option->soil_cat) {
      case CS_ATMO_SOIL_5_CAT:
        soil_num = 5;
        break;
      case CS_ATMO_SOIL_7_CAT:
        soil_num = 7;
        break;
      case CS_ATMO_SOIL_23_CAT:
        soil_num = 23;
        break;
    }

    /* Note the size is soil_num+1, first value is undef.
     * */
    f = cs_field_create("atmo_soil_percentages",
                        field_type,
                        z->location_id,
                        soil_num + 1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil percentages");

    /* Boundary variable fields for the soil */
    /*---------------------------------------*/

    /* Soil surface temperature in Celcius */
    f = cs_field_create("soil_temperature",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil T");

    /* Soil surface potential temperature (K) */
    f = cs_field_create("soil_pot_temperature",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        true); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil theta");

    /* Soil total water content */
    f = cs_field_create("soil_total_water",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        true); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil qw");

    /* ratio of the shallow reservoir water content to its maximum capacity */
    f = cs_field_create("soil_w1",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil w1");

    /* ratio of the deep reservoir water content to its maximum capacity */
    f = cs_field_create("soil_w2",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil w2");

    /* Incident solar radiative flux */
    f = cs_field_create("soil_solar_incident_flux",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil solar incid flux");

    /* Incident solar radiative flux */
    f = cs_field_create("soil_infrared_incident_flux",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil IR incid flux");

    /* Boundary parameters fields characterizing soil */
    /*------------------------------------------------*/

    f = cs_field_create("boundary_albedo",
                        field_type,
                        /* Note: as for boundary_roughness,
                         * location can be reduced in the future */
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Albedo");

    f = cs_field_by_name_try("emissivity");
    if (f == NULL)
      f = cs_field_create("emissivity",
                          field_type,
                          /* Note: as for boundary_roughness,
                           * location can be reduced in the future */
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1, /* dim */
                          false); /* has_previous */

    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Emissivity");

    f = cs_field_create("boundary_vegetation",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Vegetation");

    /* maximum water capacity of shallow reservoir */
    f = cs_field_create("soil_water_capacity",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil c1w");

    /* ratio of the maximum water capacity of the shallow reservoir to the deep
     * reservoir [0,1] */
    f = cs_field_create("soil_water_ratio",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil c2w");

    /* Thermal inertia of the soil */
    f = cs_field_create("soil_thermal_capacity",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil Cp");

    f = cs_field_create("soil_r1",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil r1");

    f = cs_field_create("soil_r2",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil r2");

    /* Deep soil temperature (in Celsius)
     * FIXME potential value? */
    f = cs_field_create("soil_temperature_deep",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil deep T");

    /* Fields usefull for heat budget plot on the soil boundary
     * */
    f = cs_field_create("soil_sensible_heat",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil sensible heat");

    f = cs_field_create("soil_latent_heat",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil latent heat");

    f = cs_field_create("soil_thermal_rad_upward",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil thermal radiation upward");

    f = cs_field_create("soil_thermal_rad_downward",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil thermal radiation downward");

    f = cs_field_create("soil_visible_rad_absorbed",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil visible radiation absorbed");

    if (cs_glob_atmo_option->soil_meb_model > CS_ATMO_SOIL_GENUINE) {
      f = cs_field_create("cover_geometry_ratio",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover geometry ratio");

      f = cs_field_create("cover_reflectivity",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover reflectivity");

      f = cs_field_create("cover_temperature_radiative",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover temperature radiation");
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deardorff force restore model
 */
/*----------------------------------------------------------------------------*/

void
cs_soil_model(void)
{
  int z_id = cs_glob_atmo_option->soil_zone_id;
  if (z_id > -1) {
    int micro_scale_option = cs_glob_atmo_option->soil_meb_model;
    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(CS_F_(t));

    const cs_domain_t *domain = cs_glob_domain;
    cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

    const cs_real_t stephn = cs_physical_constants_stephan;
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);
    const cs_lnum_t *elt_ids = z->elt_ids;
    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
    const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;
    cs_real_t *dt = cs_field_by_name("dt")->val;
    /* Post treatment fields */
    cs_field_t *soil_sensible_heat = cs_field_by_name_try("soil_sensible_heat");
    cs_field_t *soil_latent_heat = cs_field_by_name_try("soil_latent_heat");
    cs_field_t *soil_thermal_rad_upward = cs_field_by_name_try("soil_thermal_rad_upward");
    cs_field_t *soil_thermal_rad_downward = cs_field_by_name_try("soil_thermal_rad_downward");
    cs_field_t *soil_visible_rad_absorbed = cs_field_by_name_try("soil_visible_rad_absorbed");

    /* Soil related fields */
    cs_field_t *soil_temperature = cs_field_by_name("soil_temperature");
    cs_field_t *soil_pot_temperature = cs_field_by_name("soil_pot_temperature");
    cs_field_t *soil_total_water = cs_field_by_name("soil_total_water");
    cs_field_t *soil_w1 = cs_field_by_name("soil_w1");
    cs_field_t *soil_w2 = cs_field_by_name("soil_w2");
    cs_real_t *f_fos = cs_field_by_name("soil_solar_incident_flux")->val;
    cs_real_t *f_foir = cs_field_by_name("soil_infrared_incident_flux")->val;
    cs_field_t *soil_temperature_deep = cs_field_by_name("soil_temperature_deep");
    cs_field_t *soil_r1 = cs_field_by_name("soil_r1");
    cs_field_t *soil_r2 = cs_field_by_name("soil_r2");
    cs_field_t *soil_water_capacity = cs_field_by_name("soil_water_capacity");
    cs_field_t *soil_water_ratio = cs_field_by_name("soil_water_ratio");
    cs_field_t *soil_thermal_capacity = cs_field_by_name("soil_thermal_capacity");
    cs_field_t *soil_percentages = cs_field_by_name("atmo_soil_percentages");
    cs_field_t *boundary_vegetation = cs_field_by_name("boundary_vegetation");
    /* Fields related to all faces */
    cs_field_t *boundary_albedo = cs_field_by_name_try("boundary_albedo");
    cs_field_t *emissivity = cs_field_by_name_try("emissivity");
    /* Cell fields  */
    cs_field_t *atm_total_water = cs_field_by_name_try("ym_water");
    cs_field_t *atm_temp = CS_F_(t);
    cs_field_t *meteo_pressure = cs_field_by_name_try("meteo_pressure");
    /* Radiative tables */
    cs_real_t *sold = (cs_real_t *)  cs_glob_atmo_option->rad_1d_sold;
    cs_real_t *ird = (cs_real_t *) cs_glob_atmo_option->rad_1d_ird;
    /* Pointer to the spectral flux density field */
    cs_field_t *f_qinspe = NULL;
    if (rt_params->atmo_model != CS_RAD_ATMO_3D_NONE)
      f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

    /* Exchange coefficients*/
    cs_real_t *h_t = CS_F_(t)->bc_coeffs->bf;
    cs_real_t *h_q = atm_total_water->bc_coeffs->bf;

    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    /* Deardorff parameterisation */
    const cs_real_t tau_1 = 86400.;

    cs_air_fluid_props_t *ct_prop = cs_glob_air_props;

    /* Update previous values */
    cs_field_current_to_previous(cs_field_by_name("soil_pot_temperature"));
    cs_field_current_to_previous(cs_field_by_name("soil_total_water"));

    /* In case of multi energy balance (MEB) models including PV */
    cs_field_t *cover_geometry_ratio = cs_field_by_name_try("cover_geometry_ratio");
    cs_field_t *cover_reflectivity = cs_field_by_name_try("cover_reflectivity");
    cs_field_t *cover_temperature_radiative = cs_field_by_name_try("cover_temperature_radiative");

    for (cs_lnum_t soil_id = 0; soil_id < z->n_elts; soil_id++) {
      cs_real_t ray2 = 0;
      cs_real_t chas2 = 0;
      cs_real_t chal2 = 0;
      cs_real_t rapp2 = 0;
      cs_real_t secmem = 0.;
      cs_real_t premem = 0.;
      cs_real_t pphy = 0.;
      cs_real_t dum = 0.;
      cs_real_t qvs_plus = 0.;
      cs_real_t w1_min = 0.;
      cs_real_t w1_max = 1.;
      cs_real_t w2_min = 0.;
      cs_real_t w2_max = 1.;
      cs_real_t precip = 0.;
      cs_real_t tseuil = 16. + cs_physical_constants_celsius_to_kelvin;
      cs_real_t ts_plus = 0.;
      cs_real_t cpvcpa = ct_prop->cp_v / ct_prop->cp_a;
      cs_real_t rvsra = phys_pro->rvsra;
      cs_real_t clatev = phys_pro->clatev;
      cs_real_t cp0 = phys_pro->cp0;
      cs_real_t rair = phys_pro->r_pg_cnst;
      cs_real_t ps = cs_glob_atmo_constants->ps;
      cs_real_t esat = 0.;

      cs_lnum_t face_id = elt_ids[soil_id];
      cs_lnum_t cell_id = b_face_cells[face_id];
      cs_real_t foir = 0.;
      cs_real_t fos = 0.;

      cs_real_t gcr = 0.;
      cs_real_t refl = 0.;
      if (cover_geometry_ratio != NULL)
        gcr = cover_geometry_ratio->val[soil_id];
      if (cover_reflectivity != NULL)
        refl = cover_reflectivity->val[soil_id];

      /* Infrared and Solar radiative fluxes
       * Warning: should be adapted for many verticales */
      if (cs_glob_atmo_option->radiative_model_1d == 1
         && rt_params->atmo_model == CS_RAD_ATMO_3D_NONE) {
        foir = ird[0];
        fos = sold[0];
      }
      /* In case of 3D, take the incident flux */
      else if (rt_params->atmo_model != CS_RAD_ATMO_3D_NONE) {

        /* Number of bands (stride) */
        cs_lnum_t stride = rt_params->nwsgg;

        /* Compute fos */
        int gg_id = rt_params->atmo_dr_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_dr_o3_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_df_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_df_o3_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        /* Compute foir */
        gg_id = rt_params->atmo_ir_id;
        if (gg_id > -1)
          foir += f_qinspe->val[gg_id + face_id * stride];
      }

      /* f_fos and f_foir store the radiations coming ontop of
       * the boundary cell */
      f_fos[soil_id] = fos;
      f_foir[soil_id] = foir;

      if (cs_glob_atmo_option->meteo_profile == 0) {
        cs_atmo_profile_std(cell_cen[cell_id][2], &pphy, &dum, &dum);
      }
      else if (cs_glob_atmo_option->meteo_profile == 1) {
        int nbmett = cs_glob_atmo_option->met_1d_nlevels_t;
        int nbmetm = cs_glob_atmo_option->met_1d_ntimes;
        pphy = cs_intprf(nbmett,
            nbmetm,
            cs_glob_atmo_option->z_temp_met,
            cs_glob_atmo_option->time_met,
            cs_glob_atmo_option->hyd_p_met,
            cell_cen[cell_id][2],
            cs_glob_time_step->t_cur);
      }
      else {
        pphy = meteo_pressure->val[cell_id];
      }

      /* ====================================
       * Specific case for water faces (sea/lake)
       * ==================================== */

      /* Water is second component of soil_percentages */
      cs_lnum_t cat_id = 1 + soil_percentages->dim * soil_id;
      if (soil_percentages->val[cat_id] > 50.) {
        /* NB: soil_temperature in °C */
        esat = cs_air_pwv_sat(soil_temperature->val[soil_id]);
        qvs_plus = esat / ( rvsra * pphy
          + esat * (1.0 - rvsra));
        /* Constant temperature at the moment
         * TODO integrate the lake model GLM to match with LNHE models */
        ts_plus = soil_temperature->val[soil_id]
          + cs_physical_constants_celsius_to_kelvin;
      }
      else {

      cs_real_t rscp1 = (rair / cp0) * (1. + (rvsra - cpvcpa)
          * soil_total_water->val[soil_id] );

      /* ===============================
       * Compute coefficients for heat and latent heat fluxes
       * =============================== */

      /* ratio specific heat of humide air/ specidif heat of dry air
       * Cph/Cpd */
      cs_real_t cph_dcpd = (1. + (cpvcpa - 1.)
        * soil_total_water->val[soil_id] );
      /* Conversion theta -> T */
      cs_real_t cht = h_t[face_id] * pow(ps / pphy, rscp1) * cph_dcpd;

      cs_real_t chq = h_q[face_id]
        * (clatev - 2370.* (soil_temperature->val[soil_id]) );

      /* ===============================
       * Compute reservoirs water (shallow and deep)
       * =============================== */

      cs_real_t evapor = - h_q[face_id] * (atm_total_water->val[cell_id]
        - soil_total_water->val[soil_id]);
      cs_real_t dtref  = dt[cell_id];

      cs_real_t w1_num = soil_w1->val[soil_id]
        + dtref * (precip - evapor) / soil_water_capacity->val[soil_id]
        + soil_w2->val[soil_id] * dtref
          / (tau_1 + soil_water_ratio->val[soil_id] * dtref);
      cs_real_t w1_den = 1.
        + 1. / (tau_1/dtref + soil_water_ratio->val[soil_id]);
      cs_real_t w1_plus = w1_num / w1_den;
      w1_plus = CS_MAX(w1_plus, w1_min);
      w1_plus = CS_MIN(w1_plus, w1_max);
      cs_real_t w2_num = soil_w2->val[soil_id] * tau_1
        + w1_plus * dtref * soil_water_ratio->val[soil_id];
      cs_real_t w2_den = tau_1 + dtref * soil_water_ratio->val[soil_id];
      cs_real_t w2_plus = w2_num / w2_den;
      w2_plus = CS_MAX(w2_plus, w2_min);
      w2_plus = CS_MIN(w2_plus, w2_max);

      soil_w1->val[soil_id] = w1_plus;
      soil_w2->val[soil_id] = w2_plus;

      cs_real_t hu = 0.5 * (1. - cos(cs_math_pi * w1_plus));

      /* ============================
       * Compute saturated pressure and DL1
       * ============================ */

      esat = cs_air_pwv_sat(soil_temperature->val[soil_id]);
      cs_real_t rapsat = rvsra * pphy + esat * (1. - rvsra);
      cs_real_t qsat = esat / rapsat;
      cs_real_t cstder = 17.438 * 239.78; //#TODO transfer in cs_air_props.c
      cs_real_t dqsat = pphy * rvsra * cstder * esat
        / cs_math_pow2(rapsat * (soil_temperature->val[soil_id] + 239.78));

      /* Compute equivalent emissivity and reflexivity du to pannels/plants */
      cs_real_t c_refl = 1. - refl;
      cs_real_t emi = emissivity->val[face_id];
      cs_real_t emi_eq = gcr * emi / (emi * refl + c_refl);
      cs_real_t c_gcr = 1. - gcr;
      /* ============================
       * Compute the first member of soil_temperature equation
       * ============================ */

      cs_lnum_t iseuil = 0;
      if (soil_temperature->val[soil_id] < tseuil)
        iseuil = 1;

      cs_lnum_t ichal = 1;

      cs_real_t ray1 = 4. * c_gcr * emissivity->val[face_id] * stephn
        * cs_math_pow3(soil_temperature->val[soil_id]
        + cs_physical_constants_celsius_to_kelvin);
      cs_real_t chas1 = cht;
      cs_real_t chal1 = chq * hu * dqsat;
      cs_real_t rapp1 = 2.* cs_math_pi / tau_1;

      premem = soil_thermal_capacity->val[soil_id] * ( ray1
        + chas1 + chal1 * ichal + soil_r2->val[soil_id] * iseuil ) + rapp1;

      /* ============================
       * Compute the second member of soil_temperature equation
       * ============================ */

      ray2 = c_gcr * fos*(1. - boundary_albedo->val[face_id])
        + emi * c_gcr * foir
        + 3. * (emi_eq * c_refl  + c_gcr) * emi * stephn
        * cs_math_pow4(soil_temperature->val[soil_id]
        + cs_physical_constants_celsius_to_kelvin);

      if ((micro_scale_option == CS_ATMO_SOIL_PHOTOVOLTAICS) ||
          (micro_scale_option == CS_ATMO_SOIL_VEGETATION)) {
        ray2 += cs_math_pow4(cover_temperature_radiative->val[soil_id]
            + cs_physical_constants_celsius_to_kelvin) * stephn
          * emi_eq * c_refl * refl;
      }
      chas2 = cht * atm_temp->val[cell_id]
        * pow(pphy / ps , rscp1);
      chal2 = chq * (atm_total_water->val[cell_id]
        * (1. - boundary_vegetation->val[soil_id] * ( 1. - hu))
        - hu * (qsat - (soil_temperature->val[soil_id]
              + cs_physical_constants_celsius_to_kelvin)
              * dqsat));
      rapp2 = 2.*cs_math_pi * (soil_temperature_deep->val[soil_id]
        + cs_physical_constants_celsius_to_kelvin) / tau_1;

      secmem = soil_thermal_capacity->val[soil_id] * ( ray2
        + chas2 + chal2 * ichal + soil_r1->val[soil_id]
        + tseuil * soil_r2->val[soil_id] * iseuil ) + rapp2;

      /* ============================
       * Compute new soil variables
       * ============================ */

      ts_plus = (soil_temperature->val[soil_id]
        + cs_physical_constants_celsius_to_kelvin + dtref * secmem )
        / ( 1. + dtref * premem );

      qvs_plus = hu * ( qsat + dqsat * ( ts_plus
        - cs_physical_constants_celsius_to_kelvin
        - soil_temperature->val[soil_id] ))
        + boundary_vegetation->val[soil_id] * atm_total_water->val[cell_id]
        * ( 1. - hu );

      /* TODO filling soil fields should be done one for loop below
       * by computing cht and chq for the ""lake model""
       * At the moment the allocation is performed ONLY for ""soil model""  */

      if (soil_latent_heat != NULL)
        soil_latent_heat->val[soil_id] = chq * (qvs_plus
            - atm_total_water->val[cell_id]);
      if (soil_sensible_heat != NULL)
        soil_sensible_heat->val[soil_id] = cht * ((ts_plus
              * (pow(ps/pphy,(rair/cp0)
                  * (1. + (rvsra - cpvcpa) * qvs_plus ))))
            - atm_temp->val[cell_id]);
      if (soil_thermal_rad_upward != NULL)
        soil_thermal_rad_upward->val[soil_id] = stephn * emi
          * cs_math_pow4(ts_plus);
      if (soil_thermal_rad_downward != NULL)
        soil_thermal_rad_downward->val[soil_id] = emi * foir;
      if (soil_visible_rad_absorbed != NULL)
        soil_visible_rad_absorbed->val[soil_id] =
          (1. - boundary_albedo->val[cell_id]) * fos;

      }

      /* ============================
       * Update new soil variables
       * ============================ */

      soil_temperature->val[soil_id] = ts_plus
        - cs_physical_constants_celsius_to_kelvin;

      soil_pot_temperature->val[soil_id] = ts_plus
        * (pow(ps/pphy,(rair/cp0)
        * (1. + (rvsra - cpvcpa) * qvs_plus )));
      soil_total_water->val[soil_id] = qvs_plus;
    }
  }
  else {
    bft_error(__FILE__,__LINE__, 0,
        _("cs_glob_atmo_option->soil_zone_id is missing \n"));
  }
}

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

  /* Reference fluid properties set from meteo values */
  phys_pro->p0 = aopt->meteo_psea;
  phys_pro->t0 = aopt->meteo_t0; /* ref temp T0 */

  /* Compute reference q_l, theta_liq and rho */
  cs_real_t t_c = aopt->meteo_t0 - cs_physical_constants_celsius_to_kelvin;
  cs_real_t q_sat = cs_air_yw_sat(t_c, aopt->meteo_psea);
  aopt->meteo_ql0 = CS_MAX(aopt->meteo_qw0 - q_sat, 0.);
  cs_real_t rvsra = phys_pro->rvsra;
  cs_real_t rhum = rair*(1. + (rvsra - 1.)*(aopt->meteo_qw0 - aopt->meteo_ql0)
                        - aopt->meteo_ql0);
  phys_pro->ro0 = phys_pro->p0/(rhum * aopt->meteo_t0); /* ref density T0 */
  cs_real_t clatev = phys_pro->clatev;
  cs_real_t theta0 = (aopt->meteo_t0 - clatev/cp0 * aopt->meteo_ql0)
                   * pow(pref/ aopt->meteo_psea, rscp);


  cs_real_t z0 = aopt->meteo_z0;
  cs_real_t zref = aopt->meteo_zref;
  bool is_humid = (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2 );
  aopt->meteo_qwstar = is_humid ? aopt->meteo_qwstar : 0;
  if (aopt->meteo_ustar0 <= 0. && aopt->meteo_uref <= 0. && aopt->meteo_u1 <= 0.
      && aopt->meteo_u2 <= 0.)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Meteo preprocessor: meteo_ustar0 or meteo_uref"
                " or velocity measurements.\n"));

  if (aopt->meteo_ustar0 > 0. && aopt->meteo_uref > 0.) {
  /* Recompute LMO inverse */

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
      bft_printf("IT %d: dlmo = %f, error = %f\n", it, dlmo, error);
#endif
    }

    if (it == it_max) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("Meteo preprocessor did not converge to find inverse\n"
                   "of LMO length, current value is %f.\n"
                   "Please, check reference velocity, reference altitude and ustar\n"),
                 dlmo);
    }

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
    aopt->meteo_uref =   aopt->meteo_ustar0 / kappa
                       * cs_mo_psim(zref + z0, z0, aopt->meteo_dlmo);

  /* LMO inverse, ustar at ground */
  cs_real_t dlmo = aopt->meteo_dlmo;
  cs_real_t ustar0 = aopt->meteo_ustar0;
  if (is_humid && (aopt->meteo_qwstar > 0.5*DBL_MAX)
      && (aopt->meteo_evapor < 0.5*DBL_MAX)) {
    aopt->meteo_qwstar = aopt->meteo_evapor / (ustar0 * phys_pro->ro0);
  }
  aopt->meteo_tstar = cs_math_pow2(ustar0) * theta0 * dlmo / (kappa * g)
                      - 0.61 * theta0 * aopt->meteo_qwstar;
    if ((aopt->meteo_zu1 < 0. && aopt->meteo_u1 > 0.)
      || (aopt->meteo_zu2 < 0. && aopt->meteo_u2 > 0.))
    bft_error(__FILE__,
              __LINE__,
              0,
              _("When computing idealised atmospheric profiles,\n"
                "No heights are provided for velocities u1 and u2.\n"));
  if ((aopt->meteo_zu1 > 0. && aopt->meteo_t1 < 0.)
      || (aopt->meteo_zu2 > 0. && aopt->meteo_t2 < 0.))
    bft_error(__FILE__,
              __LINE__,
              0,
              _("When computing idealised atmospheric profiles,\n"
                "No temperature measurements are provided.\n"));
  if ((aopt->meteo_zu1 > 0. && aopt->meteo_qw1 > 0.5*DBL_MAX && is_humid)
      || (aopt->meteo_zu2 > 0. && aopt->meteo_qw2 > 0.5*DBL_MAX && is_humid))
    bft_error(__FILE__,
              __LINE__,
              0,
              _("When computing idealised atmospheric profiles,\n"
                "No humidity measurements are provided, though humid "
                "atmosphere is selected.\n"));
  cs_real_t u1  = aopt->meteo_u1;
  cs_real_t u2  = aopt->meteo_u2;
  cs_real_t z1  = aopt->meteo_zu1;
  cs_real_t z2  = aopt->meteo_zu2;
  cs_real_t t1  = aopt->meteo_t1;
  cs_real_t t2  = aopt->meteo_t2;
  cs_real_t qw1 = aopt->meteo_qw1;
  cs_real_t qw2 = aopt->meteo_qw2;
  if ((u1 > 0.) || (u2 > 0.)) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Meteo preprocessor uses u1=%17.9e, u2=%17.9e\n"
                 "collected at z1=%17.9e [m], z2=%17.9e [m]\n"),
               u1,
               u2,
               z1,
               z2);

    cs_real_t du1u2 = u2 - u1;
    cs_real_t dt1t2 = t2 - t1; /* Should be negative for unstable case */
    cs_real_t tmoy  = (t1 + t2) * 0.5; /* Proxy for Lmo */

    cs_real_t dqw1qw2 = qw2 - qw1; /* is zero by default */

    cs_real_t tol = 1e-12;
    int it;
    int it_max = 1000;

    /* Neutral Assumption */
    cs_real_t dlmos = 0.; /* dlmos = dlmo starts */
    cs_real_t dlmou = 0.; /* dlmou = dlmo updated */
    cs_real_t err   = 1.; /* in order to enter the loop */
    for (it = 0; it < it_max && CS_ABS(err) > tol; it++) {
      if (it != 0) {
        dlmos = dlmou;
      }
      cs_real_t ustaru  = kappa * du1u2 / (cs_mo_psim(z2, z1, dlmos));
      cs_real_t tstaru  = kappa * dt1t2 / (cs_mo_psih(z2, z1, dlmos));
      cs_real_t qwstaru = kappa * dqw1qw2 / (cs_mo_psih(z2, z1, dlmos));
      dlmou = kappa * (g / tmoy) * (tstaru + 0.61 * tmoy * qwstaru)
              / (ustaru * ustaru);
      err = dlmou - dlmos;
    }
    if (it == it_max) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("Meteo preprocessor did not converge to find inverse\n"
                   "of LMO length, current value is %f.\n"),
                 dlmou);
    }

    /* Update ustar and tstar */
    aopt->meteo_dlmo   = dlmou;
    aopt->meteo_ustar0 = kappa * du1u2 / (cs_mo_psim(z2, z1, dlmou));
    aopt->meteo_tstar  = kappa * dt1t2 / (cs_mo_psih(z2, z1, dlmou));
    aopt->meteo_qwstar = kappa * dqw1qw2 / (cs_mo_psih(z2, z1, dlmou));
    aopt->meteo_qw0    = qw1
      - aopt->meteo_qwstar * cs_mo_psih(z1 + z0, z0, dlmo) / kappa;
    aopt->meteo_t0     = t1
      - aopt->meteo_tstar * cs_mo_psih(z1 + z0, z0, dlmo) / kappa;//FIXME conversion theta->T

  }

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
               " and longitude are automatically computed\n");
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

  if (is_humid && aopt->meteo_qw0 < 0. && aopt->meteo_qw1 < 0
      && aopt->meteo_qw2 < 0)
    bft_error(__FILE__,
              __LINE__,
              0,
              _("Meteo preprocessor: at least one information is required"
                " for humidity field preprocessing (qw0 || qw1 || qw2).\n"));


  bft_printf("\n Meteo preprocessing values for computation:\n"
             "  dlmo=%17.9e\n ustar=%17.9e\n tstar=%17.9e\n qwstar=%17.9e\n"
             " t0=%17.9e\n qw0=%17.9e\n ql0=%17.9e\n",
             aopt->meteo_dlmo,
             aopt->meteo_ustar0,
             aopt->meteo_tstar,
             aopt->meteo_qwstar,
             aopt->meteo_t0,
             aopt->meteo_qw0,
             aopt->meteo_ql0);
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
  bft_printf(" Computing meteo profiles from large scale meteo data\n\n");

  /* Get fields */
  cs_real_t *cpro_met_potemp = cs_field_by_name("meteo_pot_temperature")->val;
  cs_real_3_t *cpro_met_vel
    = (cs_real_3_t *) (cs_field_by_name("meteo_velocity")->val);
  cs_real_t *cpro_met_k = cs_field_by_name("meteo_tke")->val;
  cs_real_t *cpro_met_eps = cs_field_by_name("meteo_eps")->val;
  cs_field_t *f_met_qw = cs_field_by_name_try("meteo_humidity");

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

  /* Humidity field */

  cs_real_t qw0    = aopt->meteo_qw0;
  cs_real_t qwstar = aopt->meteo_qwstar;

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
    cpro_met_potemp[cell_id] =   theta0
                               + tstar / kappa * cs_mo_psih(z+z0, z0, dlmo);

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

  cs_parall_min(1, CS_REAL_TYPE, &z_lim);
  cs_parall_min(1, CS_REAL_TYPE, &u_met_min);
  cs_parall_min(1, CS_REAL_TYPE, &theta_met_min);
  if (f_met_qw != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
      cs_real_t z_grd = 0.;
      if (z_ground != NULL)
        z_grd = z_ground[cell_id];
      cs_real_t z = cell_cen[cell_id][2] - z_grd;
      f_met_qw->val[cell_id]
        = qw0 + qwstar / kappa * cs_mo_psih(z + z0, z0, dlmo);
    }
  }

  /* Very stable cases, corresponding to mode 0 in the Python prepro */
  if (z_lim < 0.5*cs_math_big_r) {
    // Clipping only if there are cells to be clipped
    bft_printf("Switching to very stable clipping for meteo profile.\n");
    bft_printf("All altitudes above %f have been modified by clipping.\n",
               z_lim);
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
                                     * (-1./(z+z0) + 1./(z_lim+z0));
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
   * ============== */

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
  cs_math_3_normalize((const cs_real_t *)(cs_glob_physical_constants->gravity),
                      normal);

  for (int i = 0; i < 3; i++)
    normal[i] *= -1.;

  /* Compute the mass flux due to V = - g / ||g||
   * ============================================ */

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
    i_massflux[face_id] = cs_math_3_dot_product(normal, i_face_normal[face_id]);

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++)
    b_massflux[face_id] = cs_math_3_dot_product(normal, b_face_normal[face_id]);

  /* Boundary conditions
   * =================== */

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
   * ====== */

  cs_real_t *rovsdt = NULL, *dpvar = NULL;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++) {
    rovsdt[cell_id] = 0.;
    dpvar[cell_id] = 0.;
  }

  /* Right hand side
   * =============== */

  cs_real_t *rhs = NULL;
  BFT_MALLOC(rhs, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    rhs[cell_id] = 0.;

  /* Norm
   * ==== */

  cs_parall_sum(1, CS_REAL_TYPE, &norm);
  cs_parall_sum(1, CS_REAL_TYPE, &ground_surf);

  if (ground_surf > 0.)
    norm = sqrt(norm / ground_surf) * mq->tot_vol;
  else {
    bft_printf("No ground BC or no gravity:"
               " no computation of ground elevation.\n");
    return;
  }

  /* Solving
   * ======= */

  /* In case of a theta-scheme, set theta = 1;
     no relaxation in steady case either */
  cs_real_t inf_norm = 1.;

  /* Overall loop in order to ensure convergence */
  for (int sweep = 0; sweep < eqp_p->nswrsm && inf_norm > eqp_p->epsrsm;
       sweep++) {

    cs_equation_iterative_solve_scalar(0,   /* idtvar: no steady state algo */
                                       -1,  /* no over loops */
                                       f->id,
                                       NULL,
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
   * ============== */

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
   * ============================================================ */

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
   *
   * p(z) = p0 (1 - g z / (Cp T0))^(Cp/R)
   * Note: Cp/R = gamma/(gamma-1)
   *=====================================================================*/

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    cs_real_t z  = cell_cen[cell_id][2] - z_min;
    cs_real_t zt = fmin(z, 11000.);
    cs_real_t factor = fmax(1. - g * zt / (cp0 * aopt->meteo_t0), 0.);
    temp->val[cell_id] = aopt->meteo_t0 * factor;

    /* Do not overwrite pressure in case of restart */
    if (has_restart == 0)
      f->val[cell_id] =   p_ground * pow(factor, 1./rscp)
                          /* correction factor for z > 11000m */
                        * exp(- g/(rair*temp->val[cell_id]) * (z - zt));

    if (idilat > 0)
      temp->val[cell_id] =   potemp->val[cell_id]
                           * pow((f->val[cell_id]/pref), rscp);
    density->val[cell_id] = f->val[cell_id] / (rair * temp->val[cell_id]);
  }

  if (has_restart == 1)
    return;

  /* Boussinesq hypothesis */
  if (idilat == 0) {
    bft_printf
      ("Meteo profiles are computed according to Boussinesq approximation.\n"
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

  /* dfext is actually a dummy used to copy _hydrostatic_pressure_compute */
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
   * ======= */

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
        temp->val[cell_id] =   potemp->val[cell_id]
                             * pow((f->val[cell_id]/pref), rscp);
      }

      /* f_ext = rho^k * g */
      cs_real_t rho_k = f->val[cell_id] / (rair * temp->val[cell_id]);
      density->val[cell_id] = rho_k;

      f_ext[cell_id][0] = rho_k * phys_cst->gravity[0];
      f_ext[cell_id][1] = rho_k * phys_cst->gravity[1];
      f_ext[cell_id][2] = rho_k * phys_cst->gravity[2];
    }

    if (cs_log_default_is_active()) {
      cs_parall_max(1, CS_REAL_TYPE, &inf_norm);
      bft_printf
        (_("Meteo profiles: iterative process to compute hydrostatic pressure\n"
           "  sweep %d, L infinity norm (delta p) / ps =%e\n"), sweep, inf_norm);
    }
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
 * \param[out]  za          zenithal angle
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
                             cs_real_t *za,
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

  /* 4 - compute of cosinus of the zenithal angle */

  *muzero = sin(decl)*sin(flat) + cos(decl)*cos(flat)*cos(hr);

  *za = acos(*muzero);

  /* 5 - compute solar azimut */
  *omega = 0.;
  if (CS_ABS(sin(*za)) > cs_math_epzero) {
    /* Cosinus of the zimut angle */
    cs_real_t co = (sin(decl)*cos(flat)-cos(decl)*sin(flat)*cos(hr))/sin(*za);
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

   /* Correction for very low zenithal angle */
   /* Optical air mass
    * cf. Kasten, F., Young, A.T., 1989. Revised optical air mass tables and
    * approximation formula.
    *
    * Note: old formula (LH74)
    * m = 35.d0/sqrt(1224.d0*muzero*muzero + 1.d0) */
#if 1
   if (*muzero > 0)
     *muzero += 0.50572 * pow(96.07995-180./cs_math_pi * *za, -1.6364);
#endif
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
       "    y center (in Lambert-93) : %6f\n\n"),
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
         "    BL height: %12f [m]\n"
         "    phim_s:    %s\n"
         "    phih_s:    %s\n"
         "    phim_u:    %s\n"
         "    phih_u:    %s\n\n"),
       cs_glob_atmo_option->meteo_z0,
       cs_glob_atmo_option->meteo_dlmo,
       cs_glob_atmo_option->meteo_ustar0,
       cs_glob_atmo_option->meteo_uref,
       cs_glob_atmo_option->meteo_zref,
       cs_glob_atmo_option->meteo_angle,
       cs_glob_atmo_option->meteo_psea,
       cs_glob_atmo_option->meteo_t0,
       cs_glob_atmo_option->meteo_tstar,
       cs_glob_atmo_option->meteo_zi,
       _univ_fn_name[cs_glob_atmo_option->meteo_phim_s],
       _univ_fn_name[cs_glob_atmo_option->meteo_phih_s],
       _univ_fn_name[cs_glob_atmo_option->meteo_phim_u],
       _univ_fn_name[cs_glob_atmo_option->meteo_phih_u]);
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
    cs_log_printf(CS_LOG_SETUP, _("  %s\n"),
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
  BFT_FREE(_atmo_option.w_met);
  BFT_FREE(_atmo_option.time_met);
  BFT_FREE(_atmo_option.hyd_p_met);
  BFT_FREE(_atmo_option.pot_t_met);
  BFT_FREE(_atmo_option.ek_met);
  BFT_FREE(_atmo_option.ep_met);

  BFT_FREE(_atmo_option.rad_1d_xy);
  BFT_FREE(_atmo_option.rad_1d_z);
  BFT_FREE(_atmo_option.rad_1d_acinfe);
  BFT_FREE(_atmo_option.rad_1d_dacinfe);
  BFT_FREE(_atmo_option.rad_1d_aco2);
  BFT_FREE(_atmo_option.rad_1d_aco2s);
  BFT_FREE(_atmo_option.rad_1d_daco2);
  BFT_FREE(_atmo_option.rad_1d_daco2s);
  BFT_FREE(_atmo_option.rad_1d_acsup);
  BFT_FREE(_atmo_option.rad_1d_acsups);
  BFT_FREE(_atmo_option.rad_1d_dacsup);
  BFT_FREE(_atmo_option.rad_1d_dacsups);
  BFT_FREE(_atmo_option.rad_1d_tauzq);
  BFT_FREE(_atmo_option.rad_1d_tauz);
  BFT_FREE(_atmo_option.rad_1d_zq);
  BFT_FREE(_atmo_option.rad_1d_zray);
  BFT_FREE(_atmo_option.rad_1d_ir_div);
  BFT_FREE(_atmo_option.rad_1d_sol_div);
  BFT_FREE(_atmo_option.rad_1d_iru);
  BFT_FREE(_atmo_option.rad_1d_ird);
  BFT_FREE(_atmo_option.rad_1d_solu);
  BFT_FREE(_atmo_option.rad_1d_sold);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
