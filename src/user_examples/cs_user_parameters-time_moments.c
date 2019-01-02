/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_convection_diffusion.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gui_util.h"
#include "cs_grid.h"
#include "cs_internal_coupling.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_model.h"
#include "cs_selector.h"
#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-time_moments.c
 *
 * \brief Time moments example
 *
 * See \subpage parameters for examples.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * User function example with simple data computation.
 *
 * This function computes a sum of 2 specific user scalars defined on cells.
 *
 * parameters:
 *   input <-- pointer to simple data array (here, containing a single
 *             character ket)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *             radial velocity for input 0, tangential for input 1, and
 *             axial for input 2
 *----------------------------------------------------------------------------*/

/*! [tmom_simple_sum_data] */
static void
_simple_data_sum(const void  *input,
                 cs_real_t   *vals)
{
  const int location_id = CS_MESH_LOCATION_CELLS;
  const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];

  const cs_real_t *s1 = cs_field_by_name("species_1")->val;
  const cs_real_t *s2 = cs_field_by_name("species_2")->val;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    vals[i] = s1[i] + s2[i];
  }
}
/*! [tmom_simple_sum_data] */

/*----------------------------------------------------------------------------
 * User function which computes the thermal flux through a boundary.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *             radial velocity for input 0, tangential for input 1, and
 *             axial for input 2
 *----------------------------------------------------------------------------*/

/*! [tmom_b_thermal_flux_data] */
static void
_boundary_thermal_flux(const void  *input,
                       cs_real_t   *vals)
{
  const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
  const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];

  cs_field_t *f = cs_thermal_model_field();

  cs_post_boundary_flux(f->name, n_elts, NULL, vals);
}

/*! [tmom_b_thermal_flux_data] */

/*----------------------------------------------------------------------------
 * User function for velocity  values for moments computation.
 *
 * With a rotating frame of reference, the velocity is separated into
 * radial, tangential, and axial components.
 *
 * parameters:
 *   input <-- pointer to simple data array (here, containing a single
 *             character ket)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *             radial velocity for input 0, tangential for input 1, and
 *             axial for input 2
 *----------------------------------------------------------------------------*/

/*! [tmom_velocity_rotation_data] */
static void
_velocity_moment_data(const void  *input,
                      cs_real_t   *vals)
{
  const char key = *((const char *)input);

  const int location_id = CS_MESH_LOCATION_CELLS;

  const cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];
  const cs_real_3_t *vel = (const cs_real_3_t *)(CS_F_(vel)->val);

  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const cs_rotation_t *rot = cs_glob_rotation;

  double omgnrm = fabs(rot->omega);

  /* Axial, tangential and radial unit vectors */

  cs_real_3_t e_ax = {rot->axis[0], rot->axis[1], rot->axis[2]};

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    cs_real_3_t e_th;
    cs_rotation_velocity(rot, cell_cen[i], e_th);

    double xnrm = sqrt(cs_math_3_square_norm(e_th));

    e_th[0] /= xnrm;
    e_th[1] /= xnrm;
    e_th[2] /= xnrm;

    cs_real_3_t e_r;
    cs_rotation_coriolis_v(rot, -1., e_th, e_r);

    xnrm = sqrt(cs_math_3_square_norm(e_r));

    e_r[0] /= xnrm;
    e_r[1] /= xnrm;
    e_r[2] /= xnrm;

    /* Radius */
    cs_real_t xr = cs_math_3_dot_product(cell_cen[i], e_r);

    /* Axial, tangential and radial components of velocity */
    cs_real_t xva = vel[i][0]*e_ax[0] + vel[i][1]*e_ax[1] + vel[i][2]*e_ax[2];
    cs_real_t xvt = vel[i][0]*e_th[0] + vel[i][1]*e_th[1] + vel[i][2]*e_th[2];
    cs_real_t xvr = vel[i][0]*e_r[0]  + vel[i][1]*e_r[1]  + vel[i][2]*e_r[2];

    /* Entrainment velocity is removed */
    xvt -= omgnrm*xr;

    /* Store value */

    if (key == 'r')
      vals[i] = xvr;
    else if (key == 't')
      vals[i] = xvt;
    else if (key == 'a')
      vals[i] = xva;

  }
}
/*! [tmom_velocity_rotation_data] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time moments.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined, and before fine control of field output options
 * is defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_moments(void)
{
  /*
   * We compute temporal means of the type <f1*f2*f3*...*fn>
   * The fi's are variables defined by fields of a same location
   * (usually cells or boundary faces)

   * The parameters for time_moment_define_by_field_ids are:
   *   name         <--  name of associated moment
   *   n_fields     <--  number of associated fields
   *   field_id     <--  ids of associated fields
   *   component_id <--  ids of matching field components (-1 for all)
   *   type         <--  moment type (CS_TIME_MOMENT_MEAN
   *                     or CS_TIME_MOMENT_VARIANCE)
   *   nt_start     <--  starting time step (or -1 to use t_start)
   *   t_start      <--  starting time
   *   restart_mode <--  behavior in case or restart:
   *                     CS_TIME_MOMENT_RESTART_RESET,
   *                     CS_TIME_MOMENT_RESTART_AUTO, or
   *                     CS_TIME_MOMENT_RESTART_EXACT
   *   restart_name <--  name in previous run, NULL for default
   */

  {
    /* Moment <U> calculated starting from time step 1000. */

    /*! [tmom_u] */
    int moment_f_id[] = {CS_F_(vel)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("U_mean",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1000, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
    /*! [tmom_u] */
  }

  /*! [tmom_rho_u] */
  {
    /* Moment <U> calculated starting from time step 1000. */

    int moment_f_id[] = {CS_F_(rho)->id, CS_F_(vel)->id};
    int moment_c_id[] = {-1, -1};
    int n_fields = 2;
    cs_time_moment_define_by_field_ids("U_mean",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1000, /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }
  /*! [tmom_rho_u] */

  {
    /* Moment <u v> is calculated from physical time 20 s
       (reinitialized at each restart). */

    /*! [tmom_rho_u_v] */
    int moment_f_id[] = {CS_F_(rho)->id, CS_F_(vel)->id, CS_F_(vel)->id};
    int moment_c_id[] = {-1, 0, 1};
    int n_fields = 3;
    cs_time_moment_define_by_field_ids("rho_u_v_mean",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       -1,     /* nt_start */
                                       20.0,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_RESET,
                                       NULL);
    /*! [tmom_rho_u_v] */
  }

  /* Moments for sum of user scalars "species_1" and "species_2". */

  {
    /*! [tmom_simple_sum] */
    const char *sum_comp_name[] = {"species_sum_mean", "species_sum_variance"};
    cs_time_moment_type_t m_type[] = {CS_TIME_MOMENT_MEAN,
                                      CS_TIME_MOMENT_VARIANCE};

    for (int i = 0; i < 2; i++) {
      cs_time_moment_define_by_func(sum_comp_name[i],
                                    CS_MESH_LOCATION_CELLS,
                                    1,                      /* field dimension */
                                    _simple_data_sum,       /* data_func */
                                    NULL,                   /* data_input */
                                    NULL,                   /* w_data_func */
                                    NULL,                   /* w_data_input */
                                    m_type[i],
                                    1000,                   /* nt_start */
                                    -1,                     /* t_start */
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }
    /*! [tmom_simple_sum] */
  }

  /* Moments for the boundary thermal flux. */

  {
    /*! [tmom_b_thermal_flux] */
    const char *sum_comp_name[] = {"thermal_flux_mean", "thermal_flux_variance"};
    cs_time_moment_type_t m_type[] = {CS_TIME_MOMENT_MEAN,
                                      CS_TIME_MOMENT_VARIANCE};

    for (int i = 0; i < 2; i++) {
      cs_time_moment_define_by_func(sum_comp_name[i],
                                    CS_MESH_LOCATION_BOUNDARY_FACES,
                                    1,                      /* field dimension */
                                    _boundary_thermal_flux, /* data_func */
                                    NULL,                   /* data_input */
                                    NULL,                   /* w_data_func */
                                    NULL,                   /* w_data_input */
                                    m_type[i],
                                    1,                      /* nt_start */
                                    -1,                     /* t_start */
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }
    /*! [tmom_b_thermal_flux] */
  }

  /* Moments for radial, tangential, and axial velocity components
     require extracting those components first, so a more advanced
     function is needed. */

  {
    /*! [tmom_velocity_rotation] */
    const char *vel_comp_name[] = {"Wr_moy", "Wt,moy", "Wa_moy"};

    /* Data input must be "static" so it can be used in later calls */
    static char vel_comp_input[3] = {'r', 't', 'a'};

    for (int comp_id = 0; comp_id < 3; comp_id++)  {
      cs_time_moment_define_by_func(vel_comp_name[comp_id],
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    _velocity_moment_data,       /* data_func */
                                    &(vel_comp_input[comp_id]),  /* data_input */
                                    NULL,                        /* w_data_func */
                                    NULL,                        /* w_data_input */
                                    CS_TIME_MOMENT_MEAN,
                                    74000,                       /* nt_start */
                                    -1,                          /* t_start */
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }
    /*! [tmom_velocity_rotation] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
