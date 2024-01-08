/*============================================================================
 * User function. Define immersed boundaries in time and space.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_stl.h"
#include "cs_ibm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_ibm.c
 *
 * \brief User function. Define immersed boundaries in time and space.
 */
/*----------------------------------------------------------------------------*/

/* ===========================================================================
 * Definitions of cut-cell functions for different objets need to be placed
 * here!
 * ===========================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut-cell user function template. Define immersed boundaries in time
 *         and space (solid(s) interior part).
 *
 *  This function is called several times during each time step.
 *
 *  Ipenal: 1 means only solid and 0 only fluid.
 *
 *  Warning, porosity values have to be 0 or 1.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 *
 * \returns    ipenal       indicator for cut cells algo (int)
 */
/*----------------------------------------------------------------------------*/

/* -------------------------------------------------------------------------- */
/* Solid object for x < 0.5 */
/* -------------------------------------------------------------------------- */

static int
cutcell_func1(const cs_lnum_t c_id,
              const cs_real_3_t xyz,
              const cs_real_t t,
              const int num_object)
{

  /* Function which defines a solid for x < 0.5 */

  int ipenal = 0;

  if  (xyz[0] < 0.5)
    ipenal = 1;

  return ipenal;
}

/* -------------------------------------------------------------------------- */
/* Dam break gate opening between t = 0 and t = 0.127 s */
/* -------------------------------------------------------------------------- */
static int
cutcell_dambreak_func(const cs_lnum_t c_id,
                      const cs_real_3_t xyz,
                      const cs_real_t t,
                      const int num_object)
{

  int ipenal = 0;

  if (t <= 0.127) {
    cs_real_t gate_pos = -2.8440e2*cs_math_pow3(t)
                       + 7.2640e1*cs_math_pow2(t)
                       + 8.8722e-2*t
                       -9.5122e-4;

    if (xyz[0] < 0.21 && xyz[0] > 0.2 && xyz[1] > gate_pos)
      ipenal = 1;
  }

  return ipenal;
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function in which the user defines the objects to model.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_define_objects(void)
{
  /*!< [add_cutcell_object] */
  /* Add object using a cutcell method */
  cs_ibm_add_object_from_func("cutcell_object1",  /* Name of the object */
                              cutcell_func1, /* Pointer to the cutcell function of the object */
                              false,  /* Object used in the FSI resolution if true */
                              0); /* Number of nodes if the object is deformable (0 otherwise) */
  /*!< [add_cutcell_object] */

  /*!< [add_stl_object] */
  /* Add object from an STL file */
  cs_ibm_add_object_from_file("stl_object1", /* Name of the object */
                              CS_IBM_ALGO_STL, /* Porosity computation method */
                              "toto.stl", /* File name */
                              false); /* Object used in the FSI resolution if true */
  int n_pts = 2;
  cs_real_t pts_coords[6] = { 0.9, 0.9, 1.9,
                             -0.9,-0.9,-1.9};

  cs_ibm_stl_define_ext_points("stl_object1", /* Name of the object */
                               n_pts,
                               pts_coords);
  /*!< [add_stl_object] */

  /*!< [add_med_object] */
  /* Add object from a MED file */
  cs_ibm_add_object_from_file("cylinder", /* Name of the object */
                              CS_IBM_ALGO_MEDCOUPLING, /* Porosity computation method */
                              "cylinder.med", /* File name */
                              true); /* Object used in the FSI resolution if true */
  /*!< [add_med_object] */

  /* Eventually move the object at initialisation */

  /* Rotation */
  cs_real_t theta     = cs_math_pi / 2.0;  // Angle of rotation
  cs_real_t axis[3]   = {0.0, 1.0, 0.0}; // Axis of rotation
  cs_real_t center[3] = {0.0, 0.0, 0.1}; // Center of rotation

  cs_ibm_object_rotate("object_to_rotate",
                                     theta,
                                     axis,
                                     center);

  /* Translation */
  cs_real_t trans[3]  = {0.0, 0.0, 0.1}; // Translation vector
  cs_ibm_object_translate("object_to_translate",
                                        trans);

  /* Scaling */
  cs_real_t scaling_factor = 1.2;
  cs_ibm_object_scale("object_to_scale",
                                    scaling_factor);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function to set global parameters for the immersed boundary
 *         module.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_parameters(void)
{
  /*!< [set_dim] */
  /* Set the problem dimension type to optimize the solid boundary
   * detection algorithm:
  CS_IBM_3D    ->  3D computation
  CS_IBM_2D_X  ->  2D computation with symmetry in X direction
  CS_IBM_2D_Y  ->  2D computation with symmetry in Y direction
  CS_IBM_2D_Z  ->  2D computation with symmetry in Z direction */
  cs_ibm->prob_dim = CS_IBM_2D_Z;
  /*!< [set_dim] */

  /*!< [algo_choice] */
  /* Algorithm choice for porosity computation:
     CS_IBM_ALGO_CUT_CELLS   ->  Cut-cells: optimised cutting (default)
     CS_IBM_ALGO_MEDCOUPLING -> MEDCoupling volume intersection
     CS_IBM_ALGO_STL -> computation from binary STL file
  */
  cs_ibm->algo_choice = CS_IBM_ALGO_CUT_CELLS;
  /*!< [algo_choice] */

  /* Number of sub-cut for cells and faces */

  /*!< [n_sub_cut] */
  cs_ibm->nb_cut_cells = 1;
  cs_ibm->nb_cut_faces = 1;
  /*!< [n_sub_cut] */

  /*!< [vel_bc_choice] */
  /* Choice for velocity B.C. at walls
  CS_IBM_SLIP_WALL_CONDITION     ->  Slip wall b.c.
  CS_IBM_NO_SLIP_WALL_CONDITION  ->  No-Slip wall b.c.
  CS_IBM_WALL_LAW_WALL_CONDITION ->  Wall law b.c. */
  cs_ibm->wall_condition = CS_IBM_NO_SLIP_WALL_CONDITION;
  /*!< [vel_bc_choice] */

  /* Default: no user source term for porosity (abrasion for example) */

  /*!< [ts_poros] */
  cs_ibm->porosity_user_source_term_modification = false;
  /*!< [ts_poros] */

  /* Ensure the same volume for porous object at each iteration
   * (default: none) */
  cs_ibm->ensure_isovol = false;
  /* Cell porosity based on nodes porosity (smoothing effect)
   * (default: none) */
  cs_ibm->porosity_from_nodes = false;

  /* Bounding box limiting the varaible porosity computation
   * (default: whole domain) */
  /*!< [bounding_box] */
  cs_ibm->xyzmin_moving_porosity[0] = -1.e20;
  cs_ibm->xyzmin_moving_porosity[1] = -1.e20;
  cs_ibm->xyzmin_moving_porosity[2] = -1.e20;
  cs_ibm->xyzmax_moving_porosity[0] = +1.e20;
  cs_ibm->xyzmax_moving_porosity[1] = +1.e20;
  cs_ibm->xyzmax_moving_porosity[2] = +1.e20;
  /*!< [bounding_box] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function where to apply predefined transformations to med/stl
 *         based objects.
 *
 * \param[in]  t            time value for the current time step
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_object_transformations(const cs_real_t time)
{
  CS_UNUSED(time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function which allows the definition of a 'porous' object.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_solid_por(const cs_lnum_t    c_id,
                      const cs_real_3_t  xyz,
                      const cs_real_t    t,
                      const int          num_object)
{
  /*!< [solid_intern_poro] */
  cs_ibm->solid_porosity[c_id] = 0.;
  /*!< [solid_intern_poro] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
