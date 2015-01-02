/*============================================================================
 * Fan modeling through velocity source terms.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_fan.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_fan.c
        Fan modeling through velocity source terms.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure associated to a fan */

struct _cs_fan_t {

  int            num;                    /* Fan number */
  int            model_dim;              /* 1D, 2D, or 3D modeling */
  int            fan_dim;                /* 2D or 3D geometry */

  double         inlet_axis_coords[3];   /* Axis point coordinates of the
                                            inlet face */
  double         outlet_axis_coords[3];  /* Axis point coordinates of the
                                            outlet face */
  double         axis_dir[3];            /* Unit vector of the axis
                                            (inlet to outlet) */
  double         thickness;              /* Fan thickness */
  double         surface;                /* Fan total surface */

  double         fan_radius;             /* Fan radius */
  double         blades_radius;          /* Blades radius */
  double         hub_radius;             /* Hub radius */
  double         curve_coeffs[3];        /* Coefficients of the terms of
                                            degree 0, 1 and 2 of the
                                            characteristic curve */
  double         axial_torque;           /* Fan axial torque */

  cs_lnum_t      n_cells;                /* Number of cells */

  cs_lnum_t     *cell_list;              /* List of the cells belonging
                                            to the fan */

  double         in_flow;                /* Current inlet flow */
  double         out_flow;               /* Current outlet flow */

};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Fans array */

static cs_lnum_t    cs_glob_n_fans_max = 0;

static cs_lnum_t    cs_glob_n_fans = 0;
static cs_fan_t  ** cs_glob_fans = NULL;

/*============================================================================
 * Macro definitions
 *============================================================================*/

enum {X, Y, Z};

#define CS_LOC_CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define CS_LOC_DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define CS_LOC_MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Flag the cells belonging to the different fans
 * (by the fan number, 0 otherwise)
 *
 * parameters:
 *   mesh        <-- associated mesh structure
 *   cell_fan_id --> indicator by cell
 *----------------------------------------------------------------------------*/

static void
_flag_fan_cells(const cs_mesh_t  *mesh,
                cs_lnum_t         cell_fan_id[])
{
  cs_lnum_t   cell_id;
  cs_lnum_t   fan_id;

  cs_fan_t  *fan;

  const cs_lnum_t  n_ext_cells = mesh->n_cells_with_ghosts;

  /* Flag the cells */

  for (cell_id = 0; cell_id < n_ext_cells; cell_id++)
    cell_fan_id[cell_id] = -1;

  for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {

    fan = cs_glob_fans[fan_id];

    for (cs_lnum_t i = 0; i < fan->n_cells; i++) {
      cell_id = fan->cell_list[i] - 1;
      cell_fan_id[cell_id] = fan_id;
    }

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of fans.
 *
 * Fortran interface:
 *
 * subroutine tstvtl
 * *****************
 *
 * integer          nbrvtl         : --> : number of fans
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstvtl, TSTVTL)
(
 cs_lnum_t  *const nbrvtl
)
{
  *nbrvtl = cs_glob_n_fans;
}

/*----------------------------------------------------------------------------
 * Adds a fan.
 *
 * Fortran interface:
 *
 * subroutine defvtl
 * *****************
 *
 * integer          dimmod     : <-- : fan model dimension:
 *                             :     : 1: constant_f; 2: force_profile;
 *                             :     : 3: force_profile + tangential couple
 *                  dimvtl     : <-- : fan dimension:
 *                             :     : 2: pseudo-2d (extruded mesh)
 *                             :     : 3: 3d (standard)
 * double precision xyzvt1(3)  : <-- : coo. of the axis point in inlet face
 * double precision xyzvt2(3)  : <-- : coo. of the axis point in outlet face
 * double precision rvvt       : <-- : fan radius
 * double precision rpvt       : <-- : blades radius
 * double precision rmvt       : <-- : hub radius
 * double precision ccarac(3)  : <-- : coefficients of degre 0, 1 and 2
 *                             :     : of the characteristic curve
 * double precision tauvt      : <-- : Fan axial couple
 *----------------------------------------------------------------------------*/

void CS_PROCF (defvtl, DEFVTL)
(
 const cs_int_t   *dimmod,
 const cs_int_t   *dimvtl,
 const cs_real_t   xyzvt1[3],
 const cs_real_t   xyzvt2[3],
 const cs_real_t  *rvvt,
 const cs_real_t  *rpvt,
 const cs_real_t  *rmvt,
 const cs_real_t   ccarac[3],
 const cs_real_t  *tauvt
)
{
  cs_fan_define(*dimmod,
                *dimvtl,
                xyzvt1,
                xyzvt2,
                *rvvt,
                *rpvt,
                *rmvt,
                ccarac,
                *tauvt);
}

/*----------------------------------------------------------------------------
 * Build the list of cells associated to the fans
 *
 * Fortran interface:
 *
 * subroutine inivtl
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (inivtl, INIVTL)
(
 void
)
{
  cs_fan_build_all(cs_glob_mesh,
                   cs_glob_mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Flag the fans and associate the fan number to the cells belonging to
 * thus fan, 0 otherwise.
 *
 * Fortran interface:
 *
 * subroutine numvtl (indic)
 * *****************
 *
 * integer indic(ncelet)       : --> : Fan number (0 if outside the fan)
 *----------------------------------------------------------------------------*/

void CS_PROCF (numvtl, NUMVTL)
(
 cs_int_t  indic[]
)
{
  _flag_fan_cells(cs_glob_mesh, indic);
}

/*----------------------------------------------------------------------------
 * Compute the flows through the fans
 *
 * Fortran interface:
 *
 * subroutine debvtl
 * *****************
 *
 * double precision flumas(*)      : <-- : interior faces mass flux
 * double precision flumab(*)      : <-- : boundary faces mass flux
 * double precision rhofac(*)      : <-- : density at cells
 * double precision rhofab(*)      : <-- : density at boundary faces
 * double precision debent(nbrvtl) : --> : inlet flow through the fan
 * double precision debsor(nbrvtl) : --> : Outlet flow through the fan
 *----------------------------------------------------------------------------*/

void CS_PROCF (debvtl, DEBVTL)
(
 cs_real_t  flumas[],
 cs_real_t  flumab[],
 cs_real_t  rho[],
 cs_real_t  rhofab[],
 cs_real_t  debent[],
 cs_real_t  debsor[]
)
{
  cs_fan_compute_flows(cs_glob_mesh,
                       cs_glob_mesh_quantities,
                       flumas,
                       flumab,
                       rho,
                       rhofab);

  for (int i = 0; i < cs_glob_n_fans; i++) {
    debent[i] = cs_glob_fans[i]->in_flow;
    debsor[i] = cs_glob_fans[i]->out_flow;
  }
}

/*----------------------------------------------------------------------------
 * Compute the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array crvxep (which can have other
 * contributions).
 *
 * Fortran interface:
 *
 * subroutine tsvvtl
 * *****************
 *
 * parameters:
 *  idimts         <-- Dimension associated to the source
 *                     term of velocity (1: X; 2: Y; 3: Z)
 *  crvexp         <-> Explicit source term (velocity)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tsvvtl, TSVVTL)
(
 cs_int_t  *idimts,
 cs_real_t  crvexp[]
)
{
  cs_fan_compute_force(cs_glob_mesh_quantities,
                       (*idimts) - 1,
                       crvexp);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fan definition (added to the ones previously defined)
 *
 * \param[in]    model_dim           fan model dimension:
 *                                     1: constant_f
 *                                     2: force_profile
 *                                     3: force_profile + tangential couple
 * \param[in]    fan_dim             fan dimension:
 *                                     2: pseudo-2D (extruded mesh)
 *                                     3: 3D (standard)
 * \param[in]    inlet_axis_coords   intersection coords. of axis and inlet face
 * \param[in]    outlet_axis_coords  intersection coords. od axis and outlet face
 * \param[in]    fan_radius          fan radius
 * \param[in]    blades_radius       blades radius
 * \param[in]    hub_radius          hub radius
 * \param[in]    curve_coeffs        coefficients of degre 0, 1 and 2 of
                                     the characteristic curve
 * \param[in]    axial_torque        fan axial torque
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_define(int              model_dim,
              int              fan_dim,
              const cs_real_t  inlet_axis_coords[3],
              const cs_real_t  outlet_axis_coords[3],
              cs_real_t        fan_radius,
              cs_real_t        blades_radius,
              cs_real_t        hub_radius,
              const cs_real_t  curve_coeffs[3],
              cs_real_t        axial_torque)
{
  cs_fan_t  *fan = NULL;

  /* Define a new fan */

  BFT_MALLOC(fan, 1, cs_fan_t);

  fan->num = cs_glob_n_fans + 1;

  fan->model_dim = model_dim;
  fan->fan_dim = fan_dim;

  for (int i = 0; i < 3; i++) {
    fan->inlet_axis_coords[i] = inlet_axis_coords[i];
    fan->outlet_axis_coords[i] = outlet_axis_coords[i];
  }

  fan->fan_radius = fan_radius;
  fan->blades_radius  = blades_radius;
  fan->hub_radius  = hub_radius;

  for (int i = 0; i < 3; i++)
    fan->curve_coeffs[i] = curve_coeffs[i];
  fan->axial_torque = axial_torque;

  fan->n_cells = 0;
  fan->cell_list = NULL;

  /* Compute the axis vector */

  fan->thickness = 0.0;

  for (int i = 0; i < 3; i++) {
    fan->axis_dir[i] = outlet_axis_coords[i] - inlet_axis_coords[i];
    fan->thickness += (fan->axis_dir[i] * fan->axis_dir[i]);
  }
  fan->thickness = sqrt(fan->thickness);

  for (int i = 0; i < 3; i++)
    fan->axis_dir[i] /= fan->thickness;

  /* Surface initialized to 0, will be set by cs_fan_cree_listes */

  fan->surface = 0.0;

  /* Flows initialized to 0 */

  fan->in_flow = 0.0;
  fan->out_flow = 0.0;

  /* Increase the fans array if necessary */

  if (cs_glob_n_fans == cs_glob_n_fans_max) {
    cs_glob_n_fans_max = (cs_glob_n_fans_max + 1) * 2;
    BFT_REALLOC(cs_glob_fans, cs_glob_n_fans_max, cs_fan_t *);
  }

  /* Adds in the fans array */

  cs_glob_fans[cs_glob_n_fans] = fan;
  cs_glob_n_fans += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy the structures associated with fans.
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_destroy_all(void)
{
  int i;

  cs_fan_t  *fan = NULL;

  for (i = 0; i < cs_glob_n_fans; i++) {

    fan = cs_glob_fans[i];

    BFT_FREE(fan->cell_list);

    BFT_FREE(fan);

  }

  cs_glob_n_fans_max = 0;
  cs_glob_n_fans = 0;
  BFT_FREE(cs_glob_fans);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different fans.
 *
 * \param[in]   mesh             associated mesh structure
 * \param[in]   mesh_quantities  mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_build_all(const cs_mesh_t              *mesh,
                 const cs_mesh_quantities_t   *mesh_quantities)
{
  cs_lnum_t  cell_id, cell_id_1, cell_id_2;
  cs_lnum_t  face_id;
  cs_lnum_t  fan_id;
  cs_lnum_t  coo_id;

  cs_real_t  coo_axe;
  cs_real_t  d_2_axe;
  cs_real_t  d_cel_axe[3];
  cs_real_t  l_surf;

  cs_fan_t  *fan = NULL;
  cs_lnum_t  *cpt_cel_vtl = NULL;
  cs_lnum_t  *cell_fan_id = NULL;

  const cs_lnum_t  n_ext_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t  *)(mesh->i_face_cells);
  const cs_lnum_t  *b_face_cells = mesh->b_face_cells;
  const cs_real_t  *coo_cen  = mesh_quantities->cell_cen;
  const cs_real_t  *surf_fac = mesh_quantities->i_face_normal;
  const cs_real_t  *surf_fbr = mesh_quantities->b_face_normal;

  /* Create an array for cells flaging */
  /*-----------------------------------*/

  BFT_MALLOC(cell_fan_id, n_ext_cells, cs_lnum_t);

  for (cell_id = 0; cell_id < n_ext_cells; cell_id++)
    cell_fan_id[cell_id] = -1;

  /* Main loop on cells */

  for (cell_id = 0; cell_id < n_ext_cells; cell_id++) {

    /* Loop on fans */

    for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {

      fan = cs_glob_fans[fan_id];

      /* Vector from the outlet face axis point to the cell centre */

      for (coo_id = 0; coo_id < 3; coo_id++) {
        d_cel_axe[coo_id] =   (coo_cen[cell_id*3 + coo_id])
                          - fan->inlet_axis_coords[coo_id];
      }

      /* Dot product with the axis vector */

      coo_axe = (  d_cel_axe[0] * fan->axis_dir[0]
                 + d_cel_axe[1] * fan->axis_dir[1]
                 + d_cel_axe[2] * fan->axis_dir[2]);

      /* Cell potentially in the fan if its centre projection on the axis
         is within the thickness */

      if (coo_axe >= 0.0 && coo_axe <= fan->thickness) {

        /* Projection of the vector from the outlet face axis point
           to the cell centre in the fan plane */

        for (coo_id = 0; coo_id < 3; coo_id++)
          d_cel_axe[coo_id] -= coo_axe * fan->axis_dir[coo_id];

        /* Square distance to the axis */

        d_2_axe = (  d_cel_axe[0] * d_cel_axe[0]
                   + d_cel_axe[1] * d_cel_axe[1]
                   + d_cel_axe[2] * d_cel_axe[2]);

        /* If the cell is in the fan */

        if (d_2_axe <= fan->fan_radius * fan->fan_radius) {

          cell_fan_id[cell_id] = fan_id;
          fan->n_cells += 1;
          break;

        }

      }

    } /* End of loop on fans */

  } /* End of main loop on cells */

  /* Create the lists of cells belonging to each fan */
  /*-------------------------------------------------*/

  BFT_MALLOC(cpt_cel_vtl, cs_glob_n_fans, cs_lnum_t);

  for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {

    fan = cs_glob_fans[fan_id];
    BFT_MALLOC(fan->cell_list, fan->n_cells, cs_lnum_t);

    cpt_cel_vtl[fan_id] = 0;
  }

  for (cell_id = 0; cell_id < n_ext_cells; cell_id++) {

    if (cell_fan_id[cell_id] > -1) {
      fan_id = cell_fan_id[cell_id];
      fan = cs_glob_fans[fan_id];
      fan->cell_list[cpt_cel_vtl[fan_id]] = cell_id + 1;
      cpt_cel_vtl[fan_id] += 1;
    }

  }

#if defined(DEBUG) && !defined(NDEBUG)
  for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {
    fan = cs_glob_fans[fan_id];
    assert(cpt_cel_vtl[fan_id] == fan->n_cells);
  }
#endif

  /* Compute each fan surface */
  /*--------------------------*/

  /* Contribution to the domain interior */

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cell_id_1 = i_face_cells[face_id][0];
    cell_id_2 = i_face_cells[face_id][1];

    if (   cell_id_1 < mesh->n_cells /* ensure the contrib is from one domain */
        && cell_fan_id[cell_id_1] != cell_fan_id[cell_id_2]) {

      l_surf = CS_LOC_MODULE((surf_fac + 3*face_id));
      if (cell_fan_id[cell_id_1] > -1) {
        fan_id = cell_fan_id[cell_id_1];
        fan = cs_glob_fans[fan_id];
        fan->surface += l_surf;
      }
      if (cell_fan_id[cell_id_2] > -1) {
        fan_id = cell_fan_id[cell_id_2];
        fan = cs_glob_fans[fan_id];
        fan->surface += l_surf;
      }
    }

  }

  /* Contribution to the domain boundary */

  for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {

    if (cell_fan_id[b_face_cells[face_id]] > -1) {
      l_surf = CS_LOC_MODULE((surf_fbr + 3*face_id));
      fan_id = cell_fan_id[b_face_cells[face_id]];
      fan = cs_glob_fans[fan_id];
      fan->surface += l_surf;
    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {
      cs_real_t g_surf;
      l_surf = (cs_glob_fans[fan_id])->surface;
      MPI_Allreduce (&l_surf, &g_surf, 1, CS_MPI_REAL, MPI_SUM,
                     cs_glob_mpi_comm);
      (cs_glob_fans[fan_id])->surface = g_surf;
    }

  }
#endif

  /* Free memory */

  BFT_FREE(cpt_cel_vtl);
  BFT_FREE(cell_fan_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the flows through the fans.
 *
 * \param[in]  mesh             mesh structure
 * \param[in]  mesh_quantities  mesh quantities
 * \param[in]  i_mass_flux      interior faces mass flux
 * \param[in]  b_mass_flux      boundary faces mass flux
 * \param[in]  c_rho            density at cells
 * \param[in]  b_rho            density at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_compute_flows(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     const cs_real_t              i_mass_flux[],
                     const cs_real_t              b_mass_flux[],
                     const cs_real_t              c_rho[],
                     const cs_real_t              b_rho[])
{
  cs_lnum_t   cell_id, cell_id_1, cell_id_2;
  cs_lnum_t   face_id;
  cs_lnum_t   fan_id;
  cs_lnum_t   coo_id;
  cs_lnum_t   i, sens;

  cs_real_t  flow;
  cs_real_t  orient[3];

  cs_fan_t  *fan = NULL;
  cs_lnum_t  *cell_fan_id = NULL;

  const cs_lnum_t  n_ext_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_t  nbr_fac = mesh->n_i_faces;
  const cs_lnum_t  nbr_fbr = mesh->n_b_faces;
  const cs_real_t  *coo_cen = mesh_quantities->cell_cen;
  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const cs_lnum_t   *b_face_cells = mesh->b_face_cells;

  /* Flag the cells */

  BFT_MALLOC(cell_fan_id, n_ext_cells, cs_lnum_t);

  _flag_fan_cells(mesh, cell_fan_id);

  /* Set the fans flows to zero */

  for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {
    fan = cs_glob_fans[fan_id];
    fan->in_flow = 0.0;
    fan->out_flow = 0.0;
  }

  /* Contribution to the domain interior */

  for (face_id = 0; face_id < nbr_fac; face_id++) {

    cell_id_1 = i_face_cells[face_id][0];
    cell_id_2 = i_face_cells[face_id][1];

    if (   cell_id_1 < mesh->n_cells /* Make sure the contrib is from one domain */
        && cell_fan_id[cell_id_1] != cell_fan_id[cell_id_2]) {

      for (coo_id = 0; coo_id < 3; coo_id++)
        orient[coo_id] =   coo_cen[cell_id_2*3 + coo_id]
                         - coo_cen[cell_id_1*3 + coo_id];

      for (i = 0; i < 2; i++) {

        cell_id = i_face_cells[face_id][i];
        fan_id = cell_fan_id[cell_id];

        if (fan_id > -1) {
          fan = cs_glob_fans[fan_id];
          flow = i_mass_flux[face_id]/c_rho[cell_id];
          sens = (i == 0 ? 1 : - 1);
          if (CS_LOC_DOT_PRODUCT(fan->axis_dir, orient) * sens > 0.0)
            fan->out_flow += flow;
          else
            fan->in_flow += flow;
        }

      }

    }

  }

  /* Contribution to the domain boundary */

  for (face_id = 0; face_id < nbr_fbr; face_id++) {

    fan_id = cell_fan_id[b_face_cells[face_id]];

    if (fan_id > -1) {

      fan = cs_glob_fans[fan_id];

      for (coo_id = 0; coo_id < 3; coo_id++)
        orient[coo_id] = mesh_quantities->b_face_normal[face_id * 3 + coo_id];

      flow = b_mass_flux[face_id]/b_rho[face_id];
      if (CS_LOC_DOT_PRODUCT(fan->axis_dir, orient) > 0.0)
        fan->out_flow += flow;
      else
        fan->in_flow += flow;

    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {

      cs_real_t flow_glob[2];
      cs_real_t flow_loc[2];

      fan = cs_glob_fans[fan_id];

      flow_loc[0] = fan->out_flow;
      flow_loc[1] = fan->in_flow;

      MPI_Allreduce (flow_loc, flow_glob, 2, CS_MPI_REAL, MPI_SUM,
                     cs_glob_mpi_comm);

      fan->out_flow = flow_glob[0];
      fan->in_flow = flow_glob[1];

    }
  }
#endif

  /* In 2D, the flow is normalized by the surface */

  if (fan->fan_dim == 2) {
    cs_real_t  surf_2d;
    surf_2d =   (0.5*fan->surface - 2*fan->fan_radius*fan->thickness)
              /                       (2*fan->fan_radius+fan->thickness);
    fan->out_flow = fan->out_flow / surf_2d;
    fan->in_flow = fan->in_flow / surf_2d;
  }

  /* Free memory */

  BFT_FREE(cell_fan_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the force induced by the fans
 *        (needs a previous calculation of the flows through each fan).
 *
 * \param[in]  mesh_quantities  mesh quantities
 * \param[in]  source_coo_id    coordinate associated to the source term
 *                              of velocity (0: X; 0: Y; 0: Z)
 * \param[in]  source_t         explicit source term for the velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_fan_compute_force(const cs_mesh_quantities_t  *mesh_quantities,
                     int                          source_coo_id,
                     cs_real_t                    source_t[])
{
  cs_lnum_t  cell_id;
  cs_lnum_t  fan_id;
  int  coo_id;

  cs_real_t  f_z, f_theta;
  cs_real_t  f_rot[3];

  const cs_real_t  *coo_cen = mesh_quantities->cell_cen;
  const cs_real_t  pi = 3.14159265358979323846;

  /* Compute the force induced by fans */

  /* Loop on fans */
  /*--------------*/

  for (fan_id = 0; fan_id < cs_glob_n_fans; fan_id++) {

    const cs_fan_t  *fan = cs_glob_fans[fan_id];

    const cs_real_t  hub_radius  = fan->hub_radius;
    const cs_real_t  blades_radius  = fan->blades_radius;
    const cs_real_t  fan_radius = fan->fan_radius;

    const cs_real_t  mean_flow = 0.5 * (  fan->out_flow
                                        - fan->in_flow);

    const cs_real_t  delta_p = - (fan->curve_coeffs[2] * mean_flow*mean_flow)
                               + (fan->curve_coeffs[1] * mean_flow)
                               + (fan->curve_coeffs[0]);

    /* Loop on fan cells */
    /*-------------------*/

    for (cs_lnum_t i = 0; i < fan->n_cells; i++) {

      cell_id = fan->cell_list[i] - 1;

      f_z = 0.0;
      f_theta = 0.0;
      f_rot[0] = 0.0, f_rot[1] = 0.0, f_rot[2] = 0.0;

      if (blades_radius < 1.0e-12 && hub_radius < 1.0e-12) {

        f_z = delta_p / fan->thickness;
        f_theta = 0.0;

      }
      else if (hub_radius < blades_radius) {

        cs_real_t  r_1, r_2, aux, aux_1, aux_2, coo_axe, d_axe, d_cel_axe[3];

        r_1 = 0.7  * fan->blades_radius;
        r_2 = 0.85 * fan->blades_radius;

        if (fan->fan_dim == 2) {
          aux_1 =   (delta_p * 2.0 * fan_radius)
                  / (fan->thickness * (1.15*blades_radius - hub_radius));
          aux_2 = 0.0;
        }
        else {
          cs_real_t f_base;
          const cs_real_t hub_radius3
            = hub_radius * hub_radius * hub_radius;
          const cs_real_t blades_radius3
            = blades_radius * blades_radius * blades_radius;
          const cs_real_t blades_radius2 = blades_radius * blades_radius;
          const cs_real_t fan_radius2 = fan_radius * fan_radius;
          f_base =   (0.7*blades_radius - hub_radius)
                   / (1.0470*fan->thickness * (  hub_radius3
                                               + 1.4560*blades_radius3
                                               - 2.570*blades_radius2*hub_radius));
          aux_1 = f_base * delta_p * pi * fan_radius2;
          aux_2 = f_base * fan->axial_torque;
        }

        /* Vector from the outlet face axis point to the cell centre */

        for (coo_id = 0; coo_id < 3; coo_id++) {
          d_cel_axe[coo_id] =   (coo_cen[cell_id*3 + coo_id])
                            - fan->inlet_axis_coords[coo_id];
        }

        /* Projection of the cell centre on the fan axis */

        coo_axe = (  d_cel_axe[0] * fan->axis_dir[0]
                   + d_cel_axe[1] * fan->axis_dir[1]
                   + d_cel_axe[2] * fan->axis_dir[2]);

        /* Projection of the vector from the outlet face axis point
           to the cell centre in the fan plane */

        for (coo_id = 0; coo_id < 3; coo_id++)
          d_cel_axe[coo_id] -= coo_axe * fan->axis_dir[coo_id];

        d_axe = CS_LOC_MODULE(d_cel_axe); /* Distance to the axis */

        CS_LOC_CROSS_PRODUCT(f_rot, fan->axis_dir, d_cel_axe);

        aux = CS_LOC_MODULE(f_rot);
        for (coo_id = 0; coo_id < 3; coo_id++)
          f_rot[coo_id] /= aux;

        if (d_axe < hub_radius) {
          f_z     = 0.0;
          f_theta = 0.0;
        }
        else if (d_axe < r_1) {
          f_z     = aux_1 * (d_axe - hub_radius) / (r_1 - hub_radius);
          f_theta = aux_2 * (d_axe - hub_radius) / (r_1 - hub_radius);
        }
        else if (d_axe < r_2) {
          f_z     = aux_1;
          f_theta = aux_2;
        }
        else if (d_axe < blades_radius) {
          f_z     = aux_1 * (blades_radius - d_axe) / (blades_radius - r_2);
          f_theta = aux_2 * (blades_radius - d_axe) / (blades_radius - r_2);
        }
        else {
          f_z     = 0.0;
          f_theta = 0.0;
        }

      }

      source_t[cell_id] +=   (f_z * fan->axis_dir[source_coo_id])
                           + (f_theta * f_rot[source_coo_id]);

    }  /* End of loop on fan cells */

  } /* End of loop on fans */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
