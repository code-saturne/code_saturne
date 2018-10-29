/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_interface.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_convection_diffusion.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell and face centers of gravity, cell volumes
 *         and update bad cells.
 *
 * \param[out]       min_vol        Minimum cell volume
 * \param[out]       max_vol        Maximum cell volume
 * \param[out]       tot_vol        Total cell volume
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh_quantities(cs_real_t  *min_vol,
                              cs_real_t  *max_vol,
                              cs_real_t  *tot_vol)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_mesh_quantities_compute(m, mq);
  cs_mesh_bad_cells_detect(m, mq);

  *min_vol = mq->min_vol;
  *max_vol = mq->max_vol;
  *tot_vol = mq->tot_vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Project the displacement on mesh vertices (solved on cell center).
 *
 * \param[in]       ialtyb        Type of boundary for ALE
 * \param[in]       meshv         Mesh velocity
 * \param[in]       gradm         Mesh velocity gradient
 *                                (du_i/dx_j : gradv[][i][j])
 * \param[in]       claale        Boundary conditions A
 * \param[in]       clbale        Boundary conditions B
 * \param[in]       dt            Time step
 * \param[out]      disp_proj     Displacement projected on vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_project_displacement(const int           ialtyb[],
                            const cs_real_3_t  *meshv,
                            const cs_real_33_t  gradm[],
                            const cs_real_3_t  *claale,
                            const cs_real_33_t *clbale,
                            const cs_real_t    *dt,
                            cs_real_3_t        *disp_proj)
{
  int  j, face_id, vtx_id, cell_id, cell_id1, cell_id2;
  bool *vtx_interior_indicator = NULL;
  cs_real_t *vtx_counter = NULL;
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const int n_vertices = m->n_vertices;
  const int n_cells = m->n_cells;
  const int n_b_faces = m->n_b_faces;
  const int n_i_faces = m->n_i_faces;
  const int dim = m->dim;
  const cs_real_3_t *restrict vtx_coord
    = (const cs_real_3_t *restrict)m->vtx_coord;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict face_cen
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  BFT_MALLOC(vtx_counter, n_vertices, cs_real_t);
  BFT_MALLOC(vtx_interior_indicator, n_vertices, bool);

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++) {

    vtx_counter[vtx_id] = 0.;
    vtx_interior_indicator[vtx_id] = true;

    for (int i = 0; i < dim; i++)
      disp_proj[vtx_id][i] = 0.;

  }

  /* All nodes wich belongs to a boundary face where the
     displacement is imposed (that is all faces except sliding BCs)
     are boundary nodes, the others are interior nodes. */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    if (ialtyb[face_id] != 2) {

      for (j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1];
           j++) {


        vtx_id = m->b_face_vtx_lst[j];
        vtx_interior_indicator[vtx_id] = false;

      } /* End of loop on vertices of the face */

    }

  } /* End of loop on border faces */


  /* Interior face and nodes treatment */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell_id1 = m->i_face_cells[face_id][0];
    cell_id2 = m->i_face_cells[face_id][1];

    cs_real_t dvol1 = 1./fvq->cell_vol[cell_id1];
    cs_real_t dvol2 = 1./fvq->cell_vol[cell_id2];

    if (cell_id1 < n_cells) { /* Test to take into account face only once */

      for (j = m->i_face_vtx_idx[face_id];
           j < m->i_face_vtx_idx[face_id+1];
           j++) {

        /* Get the vertex number */

        vtx_id = m->i_face_vtx_lst[j];

        if (vtx_interior_indicator[vtx_id]) {

          /* Get the vector from the cell center to the node */

          cs_real_3_t cen1_node;
          cs_real_3_t cen2_node;
          for (int i = 0; i < 3; i++) {
            cen1_node[i] = vtx_coord[vtx_id][i]-cell_cen[cell_id1][i];
            cen2_node[i] = vtx_coord[vtx_id][i]-cell_cen[cell_id2][i];
          }

          for (int i = 0; i < 3; i++) {
            disp_proj[vtx_id][i] +=
              dvol1*(meshv[cell_id1][i] + gradm[cell_id1][i][0]*cen1_node[0]
                                        + gradm[cell_id1][i][1]*cen1_node[1]
                                        + gradm[cell_id1][i][2]*cen1_node[2])
              *dt[cell_id1]
             +dvol2*(meshv[cell_id2][i] + gradm[cell_id2][i][0]*cen2_node[0]
                                        + gradm[cell_id2][i][1]*cen2_node[1]
                                        + gradm[cell_id2][i][2]*cen2_node[2])
              *dt[cell_id2];
          }

          vtx_counter[vtx_id] += dvol1+dvol2;

        } /* End of Interior nodes */

      }

    }

  } /* End of loop on internal faces */

  /* Border face treatment.
     only border face contribution */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = m->b_face_cells[face_id];

    for (j = m->b_face_vtx_idx[face_id];
         j < m->b_face_vtx_idx[face_id+1];
         j++) {

      vtx_id = m->b_face_vtx_lst[j];

      if (!vtx_interior_indicator[vtx_id]) {

        /* Get the vector from the face center to the node*/

        cs_real_3_t face_node;
        for (int i = 0; i<3; i++)
          face_node[i] = -face_cen[face_id][i] + vtx_coord[vtx_id][i];

        /* 1st order extrapolation of the mesh velocity at the face center
         * to the node */

        cs_real_3_t vel_node;
        for (int i = 0; i<3; i++)
          vel_node[i] = claale[face_id][i]
                      + gradm[cell_id][i][0]*face_node[0]
                      + gradm[cell_id][i][1]*face_node[1]
                      + gradm[cell_id][i][2]*face_node[2];

        cs_real_t dsurf = 1./fvq->b_face_surf[face_id];

        for (int i = 0; i<3; i++)
          disp_proj[vtx_id][i] += dsurf*dt[cell_id]*
            (vel_node[i] + clbale[face_id][i][0]*meshv[cell_id][0]
                         + clbale[face_id][i][1]*meshv[cell_id][1]
                         + clbale[face_id][i][2]*meshv[cell_id][2]);

        vtx_counter[vtx_id] += dsurf;

      } /* End of boundary nodes */

    } /* End of loop on vertices of the face */

  } /* End of loop on border faces */


  /* If the boundary face IS a sliding face.
     We project the displacment paralelly to the face. */

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    if (ialtyb[face_id] == 2) {

      for (j = m->b_face_vtx_idx[face_id];
           j < m->b_face_vtx_idx[face_id+1];
           j++) {


        vtx_id = m->b_face_vtx_lst[j];
        disp_proj[vtx_id][0] = clbale[face_id][0][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][0][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][0][2]*disp_proj[vtx_id][2];
        disp_proj[vtx_id][1] = clbale[face_id][1][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][1][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][1][2]*disp_proj[vtx_id][2];
        disp_proj[vtx_id][2] = clbale[face_id][2][0]*disp_proj[vtx_id][0]
                             + clbale[face_id][2][1]*disp_proj[vtx_id][1]
                             + clbale[face_id][2][2]*disp_proj[vtx_id][2];

      } /* End of loop on vertices of the face */

    }

  } /* End of loop on border faces */

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         3,
                         true,
                         CS_REAL_TYPE,
                         disp_proj);
    cs_interface_set_sum(m->vtx_interfaces,
                         n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         vtx_counter);
  }

  for (vtx_id = 0; vtx_id < n_vertices; vtx_id++)
    for (int i = 0; i < dim; i++)
      disp_proj[vtx_id][i] /= vtx_counter[vtx_id];

  BFT_FREE(vtx_counter);
  BFT_FREE(vtx_interior_indicator);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update mesh in the ALE framework.
 *
 * \param[in]       itrale        number of the current ALE iteration
 * \param[in]       xyzno0        nodes coordinates of the initial mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_update_mesh(const int           itrale,
                   const cs_real_3_t  *xyzno0)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  cs_real_3_t *vtx_coord = (cs_real_3_t *)m->vtx_coord;
  cs_real_3_t *mshvel, *mshvela, *disale, *disala;
  cs_var_cal_opt_t var_cal_opt;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const cs_lnum_t n_vertices = cs_glob_mesh->n_vertices;
  const int ndim = cs_glob_mesh->dim;
  cs_time_step_t *ts = cs_get_glob_time_step();

  /* Initialization */
  cs_field_get_key_struct(CS_F_(mesh_u), key_cal_opt_id, &var_cal_opt);

  if (var_cal_opt.iwarni >= 1) {
    bft_printf("\n ------------------------"
               "---------------------------"
               "---------\n\n\n"
               "  Update mesh (ALE)\n"
               "  =================\n\n");
  }

  /* Retrieving fields */
  mshvel = (cs_real_3_t *)CS_F_(mesh_u)->val;
  mshvela = (cs_real_3_t *)CS_F_(mesh_u)->val_pre;

  disale = (cs_real_3_t *)cs_field_by_name("disale")->val;
  disala = (cs_real_3_t *)cs_field_by_name("disale")->val_pre;

  /* Update geometry */
  for (int inod = 0 ; inod < n_vertices ; inod++) {
    for (int idim = 0 ; idim < ndim ; idim++) {
      vtx_coord[inod][idim] = xyzno0[inod][idim] + disale[inod][idim];
      disala[inod][idim] = vtx_coord[inod][idim] - xyzno0[inod][idim];
    }
  }

  cs_ale_update_mesh_quantities(&(fvq->min_vol),
                                &(fvq->max_vol),
                                &(fvq->tot_vol));

  /* Abort at the end of the current time-step if there is a negative volume */
  if (fvq->min_vol <= 0.) {
    ts->nt_max = ts->nt_cur;
  }

  /* The mesh velocity is reverted to its initial value if the current time step is
     the initialization time step */
  if (itrale == 0) {
    for (cs_lnum_t cell_id = 0 ; cell_id < n_cells_ext ; cell_id++) {
      for (int idim = 0 ; idim < ndim ; idim++) {
        mshvel[cell_id][idim] = mshvela[cell_id][idim];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve a Poisson equation on the mesh velocity in ALE framework.
 *
 * It also updates the mesh displacement
 * so that it can be used to update mass fluxes (due to mesh displacement).
 *
 * \param[in]       iterns        Navier-Stokes iteration number
 * \param[in]       ndircl        Number of Dirichlet BCs for mesh velocity
 * \param[in]       impale        Indicator for fixed node displacement
 * \param[in]       ialtyb        Type of boundary for ALE
 */
/*----------------------------------------------------------------------------*/

void
cs_ale_solve_mesh_velocity(const int   iterns,
                           const int   ndircl,
                           const int  *impale,
                           const int  *ialtyb)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = (const cs_lnum_t *)m->b_face_cells;
  const cs_real_t *b_dist = (const cs_real_t *)fvq->b_dist;
  const cs_real_t *b_face_surf = (const cs_real_t *)fvq->b_face_surf;
  const cs_real_3_t *b_face_normal = (const cs_real_3_t *)fvq->b_face_normal;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  /* The mass flux is necessary to call cs_equation_iterative_solve_vector
     but not used (iconv = 0), except for the free surface, where it is used
     as a boundary condition */

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *i_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(u), kimasf))->val;
  const cs_real_t *b_massflux
    = cs_field_by_id(cs_field_get_key_int(CS_F_(u), kbmasf))->val;

  /* 1. Initialization */

  cs_real_3_t rinfiv =
  { cs_math_infinite_r,
    cs_math_infinite_r,
    cs_math_infinite_r};

  cs_real_3_t *smbr;
  cs_real_33_t *fimp;
  BFT_MALLOC(smbr, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(fimp, n_cells_ext, cs_real_33_t);

  cs_real_3_t *mshvel = (cs_real_3_t *)CS_F_(mesh_u)->val;
  cs_real_3_t *mshvela = (cs_real_3_t *)CS_F_(mesh_u)->val_pre;

  cs_real_3_t *disale = (cs_real_3_t *)cs_field_by_name("disale")->val;
  cs_real_3_t *disala = (cs_real_3_t *)cs_field_by_name("disale")->val_pre;

  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(CS_F_(mesh_u), key_cal_opt_id, &var_cal_opt);

  if (var_cal_opt.iwarni >= 1) {
    bft_printf("\n   ** SOLVING MESH VELOCITY\n"
               "      ---------------------\n");
  }

  /* We compute the boundary condition on the mesh velocity at the free surface
   * from the new mass flux. */

  /* Density at the boundary */
  cs_real_t *brom = CS_F_(rho_b)->val;

  cs_field_bc_coeffs_t *bc_coeffs = CS_F_(mesh_u)->bc_coeffs;

  cs_real_3_t  *bc_a   = (cs_real_3_t  *)bc_coeffs->a;
  cs_real_3_t  *bc_af  = (cs_real_3_t  *)bc_coeffs->af;
  cs_real_33_t *bc_b   = (cs_real_33_t *)bc_coeffs->b;
  cs_real_33_t *bc_bf  = (cs_real_33_t *)bc_coeffs->bf;

  int idftnp = var_cal_opt.idften;

  /* The mesh moves in the direction of the gravity in case of free-surface */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    if (ialtyb[face_id] == CS_FREE_SURFACE) {
      cs_lnum_t cell_id = b_face_cells[face_id];
      cs_real_t distbf = b_dist[face_id];

      cs_real_6_t hintt = {0., 0., 0., 0., 0., 0.};
      if (idftnp & CS_ISOTROPIC_DIFFUSION) {
        for (int isou = 0; isou < 3; isou++)
          hintt[isou] = CS_F_(vism)->val[cell_id] / distbf;
      } else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
          for (int isou = 0; isou < 6; isou++)
            hintt[isou] = CS_F_(vism)->val[6*cell_id+isou] / distbf;
     }

     cs_real_t prosrf = cs_math_3_dot_product(grav, b_face_normal[face_id]);

     cs_real_3_t pimpv;
     for (int i = 0 ; i < 3 ; i++)
       pimpv[i] = grav[i]*b_massflux[face_id]/(brom[face_id]*prosrf);

     cs_boundary_conditions_set_dirichlet_vector_aniso((bc_a[face_id]),
                                                       (bc_af[face_id]),
                                                       (bc_b[face_id]),
                                                       (bc_bf[face_id]),
                                                       pimpv,
                                                       hintt,
                                                       rinfiv);

    }
  }

  /* 2. Solving of the mesh velocity equation */

  if (var_cal_opt.iwarni >= 1)
    bft_printf("\n\n           SOLVING VARIABLE %s\n\n",
               CS_F_(mesh_u)->name);

  for (cs_lnum_t cell_id = 0 ; cell_id < n_cells_ext ; cell_id++) {
    for (int isou = 0 ; isou < 3 ; isou++) {
      smbr[cell_id][isou] = 0.;
      for (int jsou = 0 ; jsou < 3 ; jsou++)
        fimp[cell_id][jsou][isou] = 0.;
    }
  }

  cs_real_t *i_visc, *b_visc;

  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  if (idftnp & CS_ISOTROPIC_DIFFUSION) {
    BFT_MALLOC(i_visc, n_i_faces, cs_real_t);

    cs_face_viscosity(m,
                      fvq,
                      cs_glob_space_disc->imvisf,
                      CS_F_(vism)->val,
                      i_visc,
                      b_visc);

  } else if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
    BFT_MALLOC(i_visc, 9*n_i_faces, cs_real_t);

    cs_face_anisotropic_viscosity_vector(m,
                                         fvq,
                                         cs_glob_space_disc->imvisf,
                                         (cs_real_6_t *)CS_F_(vism)->val,
                                         (cs_real_33_t *)i_visc,
                                         b_visc);
  }

  var_cal_opt.relaxv = 1.;
  var_cal_opt.thetav = 1.;
  var_cal_opt.istat  = -1;
  var_cal_opt.idifft = -1;

  cs_equation_iterative_solve_vector(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     CS_F_(mesh_u)->id,
                                     CS_F_(mesh_u)->name,
                                     ndircl,
                                     0, /* ivisep */
                                     0, /* iescap */
                                     &var_cal_opt,
                                     (const cs_real_3_t *)mshvela,
                                     (const cs_real_3_t *)mshvela,
                                     (const cs_real_3_t *)bc_coeffs->a,
                                     (const cs_real_33_t *)bc_coeffs->b,
                                     (const cs_real_3_t *)bc_coeffs->af,
                                     (const cs_real_33_t *)bc_coeffs->bf,
                                     i_massflux,
                                     b_massflux,
                                     i_visc,
                                     b_visc,
                                     i_visc,
                                     b_visc,
                                     NULL, /* i_secvis */
                                     NULL, /* b_secvis */
                                     NULL, /* viscel */
                                     NULL, /* weighf */
                                     NULL, /* weighb */
                                     0,    /* icvflv */
                                     NULL, /* icvfli */
                                     (const cs_real_33_t *)fimp,
                                     smbr,
                                     mshvel,
                                     NULL); /* eswork */

  /* Free memory */
  BFT_FREE(smbr);
  BFT_FREE(fimp);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);

  /* 3. Update nodes displacement */

  cs_real_3_t *dproj;
  cs_real_33_t *gradm;

  /* Allocate a temporary array */
  BFT_MALLOC(dproj, n_vertices, cs_real_3_t);
  BFT_MALLOC(gradm, n_cells_ext, cs_real_33_t);

  bool use_previous_t = false;
  int inc = 1;

  cs_field_gradient_vector(CS_F_(mesh_u),
                           use_previous_t,
                           inc,
                           gradm);

  cs_ale_project_displacement(ialtyb,
                              (const cs_real_3_t *)mshvel,
                              (const cs_real_33_t *)gradm,
                              (const cs_real_3_t *)bc_coeffs->a,
                              (const cs_real_33_t *)bc_coeffs->b,
                              (const cs_real_t *)CS_F_(dt)->val,
                              dproj);

  /* FIXME : warning if nterup > 1, use itrale ? */
  /* Update mesh displacement only where it is not
     imposed by the user (ie when impale <> 1) */
  for (cs_lnum_t inod = 0 ; inod < n_vertices ; inod++) {
    if (impale[inod] == 0) {
      for (int isou = 0 ; isou < 3 ; isou++)
        disale[inod][isou] = disala[inod][isou] + dproj[inod][isou];
    }
  }

  /* Free memory */
  BFT_FREE(dproj);
  BFT_FREE(gradm);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
