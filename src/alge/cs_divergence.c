/*============================================================================
 * Divergence operators.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_divergence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_divergence.c

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_mass_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (inimav, INIMAV)
(
 const cs_int_t   *const  f_id,
 const cs_int_t   *const  itypfl,
 const cs_int_t   *const  iflmb0,
 const cs_int_t   *const  init,
 const cs_int_t   *const  inc,
 const cs_int_t   *const  imrgra,
 const cs_int_t   *const  nswrgu,
 const cs_int_t   *const  imligu,
 const cs_int_t   *const  iwarnu,
 const cs_real_t  *const  epsrgu,
 const cs_real_t  *const  climgu,
 const cs_real_t          rom[],
 const cs_real_t          romb[],
 const cs_real_3_t        vel[],
 const cs_real_3_t        coefav[],
 const cs_real_33_t       coefbv[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_mass_flux(m,
               fvq,
               *f_id,
               *itypfl,
               *iflmb0,
               *init,
               *inc,
               *imrgra,
               *nswrgu,
               *imligu,
               *iwarnu,
               *epsrgu,
               *climgu,
               rom,
               romb,
               vel,
               coefav,
               coefbv,
               i_massflux,
               b_massflux);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_divergence
 *----------------------------------------------------------------------------*/

void CS_PROCF (divmas, DIVMAS)
(
 const cs_int_t  *const   init,
 const cs_real_t          i_massflux[],
 const cs_real_t          b_massflux[],
 cs_real_t                diverg[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  cs_divergence(m,
                *init,
                i_massflux,
                b_massflux,
                diverg);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_tensor_divergence
 *----------------------------------------------------------------------------*/

void CS_PROCF (divmat, DIVMAT)
(
 const cs_int_t    *const   init,
 const cs_real_3_t          i_massflux[],
 const cs_real_3_t          b_massflux[],
 cs_real_3_t                diverg[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  cs_tensor_divergence(m,
                       *init,
                       i_massflux,
                       b_massflux,
                       diverg);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_ext_force_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (projts, PROJTS)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_t          viselx[],
 const cs_real_t          visely[],
 const cs_real_t          viselz[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_ext_force_flux(m,
                    fvq,
                    *init,
                    *nswrgu,
                    frcxt,
                    cofbfp,
                    i_massflux,
                    b_massflux,
                    i_visc,
                    b_visc,
                    viselx,
                    visely,
                    viselz);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_ext_force_anisotropic_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (projtv, PROJTV)
(
 const cs_int_t  *const   init,
 const cs_int_t  *const   nswrgu,
 const cs_int_t  *const   ircflp,
 const cs_real_3_t        frcxt[],
 const cs_real_t          cofbfp[],
 const cs_real_t          i_visc[],
 const cs_real_t          b_visc[],
 const cs_real_6_t        viscel[],
 const cs_real_2_t        weighf[],
 cs_real_t                i_massflux[],
 cs_real_t                b_massflux[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_ext_force_anisotropic_flux(m,
                                fvq,
                                *init,
                                *nswrgu,
                                *ircflp,
                                frcxt,
                                cofbfp,
                                i_visc,
                                b_visc,
                                viscel,
                                weighf,
                                i_massflux,
                                b_massflux);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_tensor_face_flux
 *----------------------------------------------------------------------------*/

void CS_PROCF (divrij, DIVRIJ)
(
 const cs_int_t   *const  f_id,
 const cs_int_t   *const  itypfl,
 const cs_int_t   *const  iflmb0,
 const cs_int_t   *const  init,
 const cs_int_t   *const  inc,
 const cs_int_t   *const  imrgra,
 const cs_int_t   *const  nswrgu,
 const cs_int_t   *const  imligu,
 const cs_int_t   *const  iwarnu,
 const cs_real_t  *const  epsrgu,
 const cs_real_t  *const  climgu,
 const cs_real_t          rom[],
 const cs_real_t          romb[],
 const cs_real_6_t        tensorvel[],
 const cs_real_6_t        coefat[],
 const cs_real_66_t       coefbt[],
 cs_real_3_t              i_massflux[],
 cs_real_3_t              b_massflux[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_tensor_face_flux(m,
                      fvq,
                      *f_id,
                      *itypfl,
                      *iflmb0,
                      *init,
                      *inc,
                      *imrgra,
                      *nswrgu,
                      *imligu,
                      *iwarnu,
                      *epsrgu,
                      *climgu,
                      rom,
                      romb,
                      tensorvel,
                      coefat,
                      coefbt,
                      i_massflux,
                      b_massflux);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \vect{u} \cdot \vect{s}_\ij\f$ to
 * the mass flux \f$ \dot{m}_\ij \f$.
 *
 * For the reconstruction, \f$ \gradt \left(\rho \vect{u} \right) \f$ is
 * computed with the following approximated boundary conditions:
 *  - \f$ \vect{a}_{\rho u} = \rho_\fib \vect{a}_u \f$
 *  - \f$ \tens{b}_{\rho u} = \tens{b}_u \f$
 *
 * For the mass flux at the boundary we have:
 * \f[
 * \dot{m}_\ib = \left[ \rho_\fib \vect{a}_u  + \rho_\fib \tens{b}_u \vect{u}
 * + \tens{b}_u \left(\gradt \vect{u} \cdot \vect{\centi \centip}\right)\right]
 * \cdot \vect{s}_\ij
 * \f]
 * The last equation uses some approximations detailed in the theory guide.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f_id          field id (or -1)
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on walls and
 *                               symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     rom           cell density
 * \param[in]     romb          density at boundary faces
 * \param[in]     vel           vector variable
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part - vector array )
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part - 3x3 tensor array)
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_mass_flux(const cs_mesh_t          *m,
             cs_mesh_quantities_t     *fvq,
             int                       f_id,
             int                       itypfl,
             int                       iflmb0,
             int                       init,
             int                       inc,
             int                       imrgra,
             int                       nswrgu,
             int                       imligu,
             int                       iwarnu,
             double                    epsrgu,
             double                    climgu,
             const cs_real_t           rom[],
             const cs_real_t           romb[],
             const cs_real_3_t         vel[],
             const cs_real_3_t         coefav[],
             const cs_real_33_t        coefbv[],
             cs_real_t       *restrict i_massflux,
             cs_real_t       *restrict b_massflux)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  /* Local variables */

  char var_name[32];

  cs_real_3_t *qdm, *f_momentum, *coefaq;
  cs_real_33_t *grdqdm;

  cs_field_t *f;

  BFT_MALLOC(qdm, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(f_momentum, m->n_b_faces, cs_real_3_t);
  BFT_MALLOC(coefaq, m->n_b_faces, cs_real_3_t);

  /*==========================================================================
    1.  Initialization
    ==========================================================================*/

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else {
    strcpy(var_name, "Work array"); var_name[31] = '\0';
  }

  /* ---> Momentum computation */

  if (init == 1) {
#   pragma omp parallel for
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = rom[cell_id]*vel[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = rom[cell_id]*vel[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        qdm[cell_id][0] = ( porosf[cell_id][0]*vel[cell_id][0]
                          + porosf[cell_id][3]*vel[cell_id][1]
                          + porosf[cell_id][5]*vel[cell_id][2] )
                        * rom[cell_id];
        qdm[cell_id][1] = ( porosf[cell_id][3]*vel[cell_id][0]
                          + porosf[cell_id][1]*vel[cell_id][1]
                          + porosf[cell_id][4]*vel[cell_id][2] )
                        * rom[cell_id];
        qdm[cell_id][2] = ( porosf[cell_id][5]*vel[cell_id][0]
                          + porosf[cell_id][4]*vel[cell_id][1]
                          + porosf[cell_id][2]*vel[cell_id][2] )
                        * rom[cell_id];
      }
    }

    /* Velocity flux */
  } else {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = vel[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 3; isou++) {
          qdm[cell_id][isou] = vel[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        qdm[cell_id][0] = porosf[cell_id][0]*vel[cell_id][0]
                        + porosf[cell_id][3]*vel[cell_id][1]
                        + porosf[cell_id][5]*vel[cell_id][2];
        qdm[cell_id][1] = porosf[cell_id][3]*vel[cell_id][0]
                        + porosf[cell_id][1]*vel[cell_id][1]
                        + porosf[cell_id][4]*vel[cell_id][2];
        qdm[cell_id][2] = porosf[cell_id][5]*vel[cell_id][0]
                        + porosf[cell_id][4]*vel[cell_id][1]
                        + porosf[cell_id][2]*vel[cell_id][2];
      }
    }
  }

  /* ---> Periodicity and parallelism treatment */

  if (halo != NULL) {
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)qdm, 3);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_vect(halo, halo_type, (cs_real_t *)qdm, 3);
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = romb[face_id]*coefav[face_id][isou];
          f_momentum[face_id][isou] = romb[face_id]*vel[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = romb[face_id]
                                 *coefav[face_id][isou]*porosi[cell_id];
          f_momentum[face_id][isou] = romb[face_id]*vel[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        coefaq[face_id][0] = ( porosf[cell_id][0]*coefav[face_id][0]
                             + porosf[cell_id][3]*coefav[face_id][1]
                             + porosf[cell_id][5]*coefav[face_id][2] )
                           * romb[face_id];
        coefaq[face_id][1] = ( porosf[cell_id][3]*coefav[face_id][0]
                             + porosf[cell_id][1]*coefav[face_id][1]
                             + porosf[cell_id][4]*coefav[face_id][2] )
                           * romb[face_id];
        coefaq[face_id][2] = ( porosf[cell_id][5]*coefav[face_id][0]
                             + porosf[cell_id][4]*coefav[face_id][1]
                             + porosf[cell_id][2]*coefav[face_id][2] )
                           * romb[face_id];
        f_momentum[face_id][0] = ( porosf[cell_id][0]*vel[cell_id][0]
                             + porosf[cell_id][3]*vel[cell_id][1]
                             + porosf[cell_id][5]*vel[cell_id][2] )
                           * romb[face_id];
        f_momentum[face_id][1] = ( porosf[cell_id][3]*vel[cell_id][0]
                             + porosf[cell_id][1]*vel[cell_id][1]
                             + porosf[cell_id][4]*vel[cell_id][2] )
                           * romb[face_id];
        f_momentum[face_id][2] = ( porosf[cell_id][5]*vel[cell_id][0]
                             + porosf[cell_id][4]*vel[cell_id][1]
                             + porosf[cell_id][2]*vel[cell_id][2] )
                           * romb[face_id];
      }
    }

    /* Velocity flux */
  } else {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou];
          f_momentum[face_id][isou] = vel[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 3; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou]*porosi[cell_id];
          f_momentum[face_id][isou] = vel[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        coefaq[face_id][0] = porosf[cell_id][0]*coefav[face_id][0]
                           + porosf[cell_id][3]*coefav[face_id][1]
                           + porosf[cell_id][5]*coefav[face_id][2];
        coefaq[face_id][1] = porosf[cell_id][3]*coefav[face_id][0]
                           + porosf[cell_id][1]*coefav[face_id][1]
                           + porosf[cell_id][4]*coefav[face_id][2];
        coefaq[face_id][2] = porosf[cell_id][5]*coefav[face_id][0]
                           + porosf[cell_id][4]*coefav[face_id][1]
                           + porosf[cell_id][2]*coefav[face_id][2];
        f_momentum[face_id][0] = ( porosf[cell_id][0]*vel[cell_id][0]
                             + porosf[cell_id][3]*vel[cell_id][1]
                             + porosf[cell_id][5]*vel[cell_id][2] );
        f_momentum[face_id][1] = ( porosf[cell_id][3]*vel[cell_id][0]
                             + porosf[cell_id][1]*vel[cell_id][1]
                             + porosf[cell_id][4]*vel[cell_id][2] );
        f_momentum[face_id][2] = ( porosf[cell_id][5]*vel[cell_id][0]
                             + porosf[cell_id][4]*vel[cell_id][1]
                             + porosf[cell_id][2]*vel[cell_id][2] );
      }
    }

  }

  /*==========================================================================
    2. Compute mass flux without recontructions
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];
          double pnd = weight[face_id];
          /* u, v, w Components */
          for (int isou = 0; isou < 3; isou++) {
            i_massflux[face_id] += (pnd*qdm[ii][isou]+(1.-pnd)*qdm[jj][isou])
                                  * i_f_face_normal[face_id][isou];
          }

        }
      }
    }

    /* Boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
            face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
            face_id++) {

          /* u, v, w Components */
          for (int isou = 0; isou < 3; isou++) {
            double pfac = inc*coefaq[face_id][isou];

            /* coefbv is a matrix */
            for (int jsou = 0; jsou < 3; jsou++) {
              pfac += coefbv[face_id][jsou][isou]*f_momentum[face_id][jsou];
            }

            b_massflux[face_id] += pfac*b_f_face_normal[face_id][isou];
          }

        }
      }
    }

  }

  /*==========================================================================
    4. Compute mass flux with reconstruction technics if the mesh is
       non orthogonal
    ==========================================================================*/

  if (nswrgu > 1) {

    BFT_MALLOC(grdqdm, n_cells_ext, cs_real_33_t);


    /* Computation of qdm gradient
       (vectorial gradient, the periodicity has already been treated) */

    cs_gradient_vector(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       nswrgu,
                       iwarnu,
                       imligu,
                       epsrgu,
                       climgu,
                       (const cs_real_3_t*)coefaq,
                       coefbv,
                       qdm,
                       NULL, /* weighted gradient */
                       NULL, /* cpl */
                       grdqdm);

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double pnd = weight[face_id];

          double dofx = dofij[face_id][0];
          double dofy = dofij[face_id][1];
          double dofz = dofij[face_id][2];

          /* Terms along U, V, W */
          for (int isou = 0; isou < 3; isou++) {

            i_massflux[face_id] = i_massflux[face_id]
              /* Non-reconstructed term */
              + (pnd*qdm[ii][isou] + (1.-pnd)*qdm[jj][isou]

                 /*  --->     ->    -->      ->
                     (Grad(rho U ) . OFij ) . Sij */
                 + 0.5*(grdqdm[ii][isou][0] +grdqdm[jj][isou][0])*dofx
                 + 0.5*(grdqdm[ii][isou][1] +grdqdm[jj][isou][1])*dofy
                 + 0.5*(grdqdm[ii][isou][2] +grdqdm[jj][isou][2])*dofz
                 )*i_f_face_normal[face_id][isou];
          }

        }
      }

    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
            face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
            face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];
          double diipbx = diipb[face_id][0];
          double diipby = diipb[face_id][1];
          double diipbz = diipb[face_id][2];

          /* Terms along U, V, W */
          for (int isou = 0; isou < 3; isou++) {

            double pfac = inc*coefaq[face_id][isou];

            /* coefu is a matrix */
            for (int jsou = 0; jsou < 3; jsou++) {

              double pip = f_momentum[face_id][jsou]
                + grdqdm[ii][jsou][0]*diipbx
                + grdqdm[ii][jsou][1]*diipby
                + grdqdm[ii][jsou][2]*diipbz;

              pfac += coefbv[face_id][jsou][isou]*pip;

            }

            b_massflux[face_id] += pfac*b_f_face_normal[face_id][isou];

          }

        }
      }
    }


    /* Deallocation */
    BFT_FREE(grdqdm);

  }

  BFT_FREE(qdm);
  BFT_FREE(coefaq);
  BFT_FREE(f_momentum);

  /*==========================================================================
    6. Here, we make sure that the mass flux is null at the boundary faces of
       type symmetry and coupled walls.
    ==========================================================================*/

  if (iflmb0 == 1) {
    /* Force flumab to 0 for velocity */
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      if (fvq->b_sym_flag[face_id] == 0) {
        b_massflux[face_id] = 0.;
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add \f$ \rho \tens{r} \vect{s}_\ij\f$ to a flux.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     f_id          field id (or -1)
 * \param[in]     itypfl        indicator (take rho into account or not)
 *                               - 1 compute \f$ \rho\vect{u}\cdot\vect{s} \f$
 *                               - 0 compute \f$ \vect{u}\cdot\vect{s} \f$
 * \param[in]     iflmb0        the mass flux is set to 0 on walls and
 *                               symmetries if = 1
 * \param[in]     init          the mass flux is initialized to 0 if > 0
 * \param[in]     inc           indicator
 *                               - 0 solve an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least square gradient
 * \param[in]     nswrgu        number of sweeps for the reconstruction
 *                               of the gradients
 * \param[in]     imligu        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thanks to neighbooring gradients
 *                               - = 1 thanks to the mean gradient
 * \param[in]     iwarnu        verbosity
 * \param[in]     epsrgu        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgu        clipping coefficient for the computation of
 *                               the gradient
 * \param[in]     c_rho         cell density
 * \param[in]     b_rho         density at boundary faces
 * \param[in]     c_var         variable
 * \param[in]     coefav        boundary condition array for the variable
 *                               (explicit part - symmetric tensor array)
 * \param[in]     coefbv        boundary condition array for the variable
 *                               (implicit part - 6x6 symmetric tensor array)
 * \param[in,out] i_massflux    mass flux at interior faces \f$ \dot{m}_\fij \f$
 * \param[in,out] b_massflux    mass flux at boundary faces \f$ \dot{m}_\fib \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_face_flux(const cs_mesh_t          *m,
                    cs_mesh_quantities_t     *fvq,
                    int                       f_id,
                    int                       itypfl,
                    int                       iflmb0,
                    int                       init,
                    int                       inc,
                    int                       imrgra,
                    int                       nswrgu,
                    int                       imligu,
                    int                       iwarnu,
                    double                    epsrgu,
                    double                    climgu,
                    const cs_real_t           c_rho[],
                    const cs_real_t           b_rho[],
                    const cs_real_6_t         c_var[],
                    const cs_real_6_t         coefav[],
                    const cs_real_66_t        coefbv[],
                    cs_real_3_t     *restrict i_massflux,
                    cs_real_3_t     *restrict b_massflux)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;
  const cs_real_3_t *restrict dofij
    = (const cs_real_3_t *restrict)fvq->dofij;

  /* Local variables */

  char var_name[32];

  cs_real_6_t *c_mass_var, *b_mass_var, *coefaq;

  cs_field_t *f;

  BFT_MALLOC(c_mass_var, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(b_mass_var, m->n_b_faces, cs_real_6_t);
  BFT_MALLOC(coefaq, m->n_b_faces, cs_real_6_t);

  /*==========================================================================
    1.  Initialization
    ==========================================================================*/

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 31, "%s", f->name); var_name[31] = '\0';
  }
  else {
    strcpy(var_name, "Work array"); var_name[31] = '\0';
  }

  /* ---> Momentum computation */

  if (init == 1) {
#   pragma omp parallel for
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      for (int i = 0; i < 3; i++)
      i_massflux[face_id][i] = 0.;
    }
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      for (int i = 0; i < 3; i++)
        b_massflux[face_id][i] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_rho[cell_id]*c_var[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_rho[cell_id]*c_var[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               c_mass_var[cell_id]);

        for (int isou = 0; isou < 6; isou++)
          c_mass_var[cell_id][isou] *= c_rho[cell_id];
      }
    }

    /* Velocity flux */
  } else {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_var[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int isou = 0; isou < 6; isou++) {
          c_mass_var[cell_id][isou] = c_var[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               c_mass_var[cell_id]);
      }
    }
  }

  /* ---> Periodicity and parallelism treatment */

  if (halo != NULL) {
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)c_mass_var, 6);
    if (cs_glob_mesh->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens(halo, halo_type, (cs_real_t *)c_mass_var);
  }

  /* Standard mass flux */
  if (itypfl == 1) {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = b_rho[face_id]*coefav[face_id][isou];
          b_mass_var[face_id][isou] = b_rho[face_id]*c_var[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = b_rho[face_id]
                                 *coefav[face_id][isou]*porosi[cell_id];
          b_mass_var[face_id][isou] = b_rho[face_id]*c_var[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               coefav[face_id],
                               coefaq[face_id]);

        for (int isou = 0; isou < 6; isou++)
          coefaq[face_id][isou] *= b_rho[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               b_mass_var[face_id]);

        for (int isou = 0; isou < 6; isou++)
          b_mass_var[face_id][isou] *= b_rho[face_id];

      }
    }

    /* Velocity flux */
  } else {

    /* Without porosity */
    if (porosi == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou];
          b_mass_var[face_id][isou] = c_var[cell_id][isou];
        }
      }
      /* With porosity */
    } else if (porosi != NULL && porosf == NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];
        for (int isou = 0; isou < 6; isou++) {
          coefaq[face_id][isou] = coefav[face_id][isou]*porosi[cell_id];
          b_mass_var[face_id][isou] = c_var[cell_id][isou]*porosi[cell_id];
        }
      }
      /* With anisotropic porosity */
    } else if (porosi != NULL && porosf != NULL) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_math_sym_33_product(porosf[cell_id],
                               coefav[face_id],
                               coefaq[face_id]);

        cs_math_sym_33_product(porosf[cell_id],
                               c_var[cell_id],
                               b_mass_var[face_id]);
      }
    }

  }

  /*==========================================================================
    2. Compute mass flux without recontructions
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];
          double pnd = weight[face_id];

          cs_real_6_t f_mass_var;

          for (int isou = 0; isou < 6; isou++)
            f_mass_var[isou] = pnd*c_mass_var[ii][isou]+(1.-pnd)*c_mass_var[jj][isou];

          cs_math_sym_33_3_product_add(f_mass_var,
                                       i_f_face_normal[face_id],
                                       i_massflux[face_id]);

        }
      }
    }

    /* Boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
            face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
            face_id++) {

          cs_real_6_t f_mass_var;

          /* var_f = a + b * var_i */
          for (int isou = 0; isou < 6; isou++)
            f_mass_var[isou] = inc*coefaq[face_id][isou];

          cs_math_66_6_product_add(coefbv[face_id],
                                   b_mass_var[face_id],
                                   f_mass_var);

          cs_math_sym_33_3_product_add(f_mass_var,
                                       b_f_face_normal[face_id],
                                       b_massflux[face_id]);

        }
      }
    }

  }

  /*==========================================================================
    4. Compute mass flux with reconstruction technics if the mesh is
       non orthogonal
    ==========================================================================*/

  if (nswrgu > 1) {

    cs_real_63_t *c_grad_mvar;
    BFT_MALLOC(c_grad_mvar, n_cells_ext, cs_real_63_t);

    /* Computation of c_mass_var gradient
       (vectorial gradient, the periodicity has already been treated) */

    cs_gradient_tensor(var_name,
                       gradient_type,
                       halo_type,
                       inc,
                       nswrgu,
                       iwarnu,
                       imligu,
                       epsrgu,
                       climgu,
                       coefaq,
                       coefbv,
                       c_mass_var,
                       c_grad_mvar);

    /* Mass flow through interior faces */

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double pnd = weight[face_id];

          cs_real_t f_mass_var[6];

          for (cs_lnum_t isou = 0; isou < 6; isou++)
            /* Non-reconstructed face value */
            f_mass_var[isou] = pnd*c_mass_var[ii][isou]+(1.-pnd)*c_mass_var[jj][isou]
              /* Reconstruction: face gradient times OF */
              + 0.5*(c_grad_mvar[ii][isou][0] +c_grad_mvar[jj][isou][0])*dofij[face_id][0]
              + 0.5*(c_grad_mvar[ii][isou][1] +c_grad_mvar[jj][isou][1])*dofij[face_id][1]
              + 0.5*(c_grad_mvar[ii][isou][2] +c_grad_mvar[jj][isou][2])*dofij[face_id][2];

          cs_math_sym_33_3_product_add(f_mass_var,
                                       i_f_face_normal[face_id],
                                       i_massflux[face_id]);

        }
      }

    }

    /* Mass flow through boundary faces */

    for (int g_id = 0; g_id < n_b_groups; g_id++) {
#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
      for (int t_id = 0; t_id < n_b_threads; t_id++) {
        for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
            face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
            face_id++) {

          cs_lnum_t ii = b_face_cells[face_id];

          cs_real_6_t f_mass_var;

          /* var_f = a + b * var_I' */
          for (int isou = 0; isou < 6; isou++)
            f_mass_var[isou] = inc*coefaq[face_id][isou];


          /* Add the reconstruction to get value in I' */
          for (int jsou = 0; jsou < 6; jsou++)
              b_mass_var[face_id][jsou] += c_grad_mvar[ii][jsou][0]*diipb[face_id][0]
                + c_grad_mvar[ii][jsou][1]*diipb[face_id][1]
                + c_grad_mvar[ii][jsou][2]*diipb[face_id][2];

          cs_math_66_6_product_add(coefbv[face_id],
                                   b_mass_var[face_id],
                                   f_mass_var);

          cs_math_sym_33_3_product_add(f_mass_var,
                                       b_f_face_normal[face_id],
                                       b_massflux[face_id]);

        }
      }
    }


    /* Deallocation */
    BFT_FREE(c_grad_mvar);

  }

  BFT_FREE(c_mass_var);
  BFT_FREE(coefaq);
  BFT_FREE(b_mass_var);

  /*==========================================================================
    6. Here, we make sure that the mass flux is null at the boundary faces of
       type symmetry and coupled walls.
    ==========================================================================*/

  if (iflmb0 == 1) {
    /* Force flumab to 0 for velocity */
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      if (fvq->b_sym_flag[face_id] == 0) {
        for (int isou = 0; isou < 3; isou++)
          b_massflux[face_id][isou] = 0.;
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_divergence(const cs_mesh_t          *m,
              int                       init,
              const cs_real_t           i_massflux[],
              const cs_real_t           b_massflux[],
              cs_real_t       *restrict diverg)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      diverg[cell_id] = 0.;
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }


  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        diverg[ii] += i_massflux[face_id];
        diverg[jj] -= i_massflux[face_id];

      }
    }
  }


  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];
        diverg[ii] += b_massflux[face_id];

      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the integrated mass flux on the cells for a tensor variable.
 *
 * \f[
 * \dot{m}_i = \dot{m}_i + \sum_{\fij \in \Facei{\celli}} \dot{m}_\ij
 * \f]
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     init          indicator
 *                               - 1 initialize the divergence to 0
 *                               - 0 otherwise
 * \param[in]     i_massflux    mass flux vector at interior faces
 * \param[in]     b_massflux    mass flux vector at boundary faces
 * \param[in,out] diverg        mass flux divergence vector
 */
/*----------------------------------------------------------------------------*/

void
cs_tensor_divergence(const cs_mesh_t            *m,
                     int                         init,
                     const cs_real_3_t           i_massflux[],
                     const cs_real_3_t           b_massflux[],
                     cs_real_3_t       *restrict diverg)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_groups = m->b_face_numbering->n_groups;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        diverg[cell_id][isou] = 0.;
      }
    }
  } else if (init == 0 && n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells+0; cell_id < n_cells_ext; cell_id++) {
      for (int isou = 0; isou < 3; isou++) {
        diverg[cell_id][isou] = 0.;
      }
    }
  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }


  /*==========================================================================
    2. Integration on internal faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_i_groups; g_id++) {
#   pragma omp parallel for
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
           face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        for (int isou = 0; isou < 3; isou++) {
          diverg[ii][isou] += i_massflux[face_id][isou];
          diverg[jj][isou] -= i_massflux[face_id][isou];
        }

      }
    }
  }


  /*==========================================================================
    3. Integration on border faces
    ==========================================================================*/

  for (int g_id = 0; g_id < n_b_groups; g_id++) {
#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[(t_id*n_b_groups + g_id)*2];
           face_id < b_group_index[(t_id*n_b_groups + g_id)*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];
        for (int isou = 0; isou < 3; isou++) {
          diverg[ii][isou] += b_massflux[face_id][isou];
        }

      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_diffusion_scalar for the improved hydrostatic pressure algorithm
 * (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgu        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viselx        viscosity by cell, dir x
 * \param[in]     visely        viscosity by cell, dir y
 * \param[in]     viselz        viscosity by cell, dir z
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_force_flux(const cs_mesh_t          *m,
                  cs_mesh_quantities_t     *fvq,
                  int                       init,
                  int                       nswrgu,
                  const cs_real_3_t         frcxt[],
                  const cs_real_t           cofbfp[],
                  cs_real_t       *restrict i_massflux,
                  cs_real_t       *restrict b_massflux,
                  const cs_real_t           i_visc[],
                  const cs_real_t           b_visc[],
                  const cs_real_t           viselx[],
                  const cs_real_t           visely[],
                  const cs_real_t           viselz[])
{
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_f_face_surf;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgu <= 1) {

    /* Mass flow through interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      i_massflux[face_id] =  i_massflux[face_id]
                           + i_visc[face_id]*(
                                              ( i_face_cog[face_id][0]
                                               -cell_cen[ii][0] )*frcxt[ii][0]
                                             +( i_face_cog[face_id][1]
                                               -cell_cen[ii][1] )*frcxt[ii][1]
                                             +( i_face_cog[face_id][2]
                                               -cell_cen[ii][2] )*frcxt[ii][2]
                                             -( i_face_cog[face_id][0]
                                               -cell_cen[jj][0] )*frcxt[jj][0]
                                             -( i_face_cog[face_id][1]
                                               -cell_cen[jj][1] )*frcxt[jj][1]
                                             -( i_face_cog[face_id][2]
                                               -cell_cen[jj][2] )*frcxt[jj][2]
                                             );

    }

    /* Mass flux through boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];
      double surfn = b_f_face_surf[face_id];
      double distbf = b_dist[face_id];

      b_massflux[face_id] =  b_massflux[face_id]
                           +  b_visc[face_id]*distbf/surfn
                             *cofbfp[face_id]
                             *( frcxt[ii][0]*b_f_face_normal[face_id][0]
                               +frcxt[ii][1]*b_f_face_normal[face_id][1]
                               +frcxt[ii][2]*b_f_face_normal[face_id][2] );

    }

  /*==========================================================================
    3. Update mass flux with reconstruction technics
    ==========================================================================*/

  } else {


    /* Mass flux through interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      double pnd = weight[face_id];

      double dijpfx = dijpf[face_id][0];
      double dijpfy = dijpf[face_id][1];
      double dijpfz = dijpf[face_id][2];

      double surfn = i_f_face_surf[face_id];

      /* Recompute II' and JJ' at this level */

      double diipx = i_face_cog[face_id][0]-cell_cen[ii][0]-(1.-pnd)*dijpfx;
      double diipy = i_face_cog[face_id][1]-cell_cen[ii][1]-(1.-pnd)*dijpfy;
      double diipz = i_face_cog[face_id][2]-cell_cen[ii][2]-(1.-pnd)*dijpfz;
      double djjpx = i_face_cog[face_id][0]-cell_cen[jj][0]+pnd*dijpfx;
      double djjpy = i_face_cog[face_id][1]-cell_cen[jj][1]+pnd*dijpfy;
      double djjpz = i_face_cog[face_id][2]-cell_cen[jj][2]+pnd*dijpfz;

      i_massflux[face_id] =  i_massflux[face_id]
                           + i_visc[face_id]*(
                                               ( i_face_cog[face_id][0]
                                                -cell_cen[ii][0] )*frcxt[ii][0]
                                              +( i_face_cog[face_id][1]
                                                -cell_cen[ii][1] )*frcxt[ii][1]
                                              +( i_face_cog[face_id][2]
                                                -cell_cen[ii][2] )*frcxt[ii][2]
                                              -( i_face_cog[face_id][0]
                                                -cell_cen[jj][0] )*frcxt[jj][0]
                                              -( i_face_cog[face_id][1]
                                                -cell_cen[jj][1] )*frcxt[jj][1]
                                              -( i_face_cog[face_id][2]
                                                -cell_cen[jj][2] )*frcxt[jj][2]
                                              )
                            + surfn/i_dist[face_id]*0.5
                             *( (djjpx-diipx)*( viselx[ii]*frcxt[ii][0]
                                               +viselx[jj]*frcxt[jj][0] )
                               +(djjpy-diipy)*( visely[ii]*frcxt[ii][1]
                                               +visely[jj]*frcxt[jj][1] )
                               +(djjpz-diipz)*( viselz[ii]*frcxt[ii][2]
                                               +viselz[jj]*frcxt[jj][2] )
                              );

    }


    /* Mass flux through boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];
      double surfn = b_f_face_surf[face_id];
      double distbf = b_dist[face_id];

      b_massflux[face_id] = b_massflux[face_id]
                           + b_visc[face_id]*distbf/surfn
                            *cofbfp[face_id]
                            *( frcxt[ii][0]*b_f_face_normal[face_id][0]
                              +frcxt[ii][1]*b_f_face_normal[face_id][1]
                              +frcxt[ii][2]*b_f_face_normal[face_id][2] );

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the external source terms to the faces in coherence with
 * cs_face_anisotropic_diffusion_scalar for the improved hydrostatic pressure
 * algorithm (iphydr=1).
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     cofbfp        boundary condition array for the diffusion
 *                               of the variable (implicit part)
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_force_anisotropic_flux(const cs_mesh_t          *m,
                              cs_mesh_quantities_t     *fvq,
                              int                       init,
                              int                       nswrgp,
                              int                       ircflp,
                              const cs_real_3_t         frcxt[],
                              const cs_real_t           cofbfp[],
                              const cs_real_t           i_visc[],
                              const cs_real_t           b_visc[],
                              const cs_real_6_t         viscel[],
                              const cs_real_2_t         weighf[],
                              cs_real_t       *restrict i_massflux,
                              cs_real_t       *restrict b_massflux)
{
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_f_face_normal
    = (const cs_real_3_t *restrict)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_f_face_normal
    = (const cs_real_3_t *restrict)fvq->b_f_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;

  /* Local variables */

  double diippf[3], djjppf[3];
  double visci[3][3], viscj[3][3];

  /*==========================================================================*/

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init == 1) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      i_massflux[face_id] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
      b_massflux[face_id] = 0.;
    }

  } else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* ---> Contribution from interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      i_massflux[face_id] =  i_massflux[face_id]+ i_visc[face_id]*(
                                               ( i_face_cog[face_id][0]
                                                -cell_cen[ii][0])*frcxt[ii][0]
                                              +( i_face_cog[face_id][1]
                                                -cell_cen[ii][1])*frcxt[ii][1]
                                              +( i_face_cog[face_id][2]
                                                -cell_cen[ii][2])*frcxt[ii][2]
                                              -( i_face_cog[face_id][0]
                                                -cell_cen[jj][0])*frcxt[jj][0]
                                              -( i_face_cog[face_id][1]
                                                -cell_cen[jj][1])*frcxt[jj][1]
                                              -( i_face_cog[face_id][2]
                                                -cell_cen[jj][2])*frcxt[jj][2]
                                              );

    }

    /* ---> Contribution from boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];
      double surfn = b_f_face_surf[face_id];
      double distbf = b_dist[face_id];

      b_massflux[face_id] =  b_massflux[face_id]
                           + b_visc[face_id]*distbf/surfn
                            *cofbfp[face_id]
                            *( frcxt[ii][0]*b_f_face_normal[face_id][0]
                              +frcxt[ii][1]*b_f_face_normal[face_id][1]
                              +frcxt[ii][2]*b_f_face_normal[face_id][2] );

    }

    /*========================================================================
      3. Update mass flux with reconstruction technics
      ========================================================================*/

  } else {

    /* ---> Contribution from interior faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      /* Recompute II' and JJ' at this level */

      visci[0][0] = viscel[ii][0];
      visci[1][1] = viscel[ii][1];
      visci[2][2] = viscel[ii][2];
      visci[1][0] = viscel[ii][3];
      visci[0][1] = viscel[ii][3];
      visci[2][1] = viscel[ii][4];
      visci[1][2] = viscel[ii][4];
      visci[2][0] = viscel[ii][5];
      visci[0][2] = viscel[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      double fikdvi = weighf[face_id][0];

      /* II" = IF + FI" */
      for (int i = 0; i < 3; i++) {
        diippf[i] =  i_face_cog[face_id][i]-cell_cen[ii][i]
                   - fikdvi*(  visci[0][i]*i_f_face_normal[face_id][0]
                             + visci[1][i]*i_f_face_normal[face_id][1]
                             + visci[2][i]*i_f_face_normal[face_id][2] );
      }

      viscj[0][0] = viscel[jj][0];
      viscj[1][1] = viscel[jj][1];
      viscj[2][2] = viscel[jj][2];
      viscj[1][0] = viscel[jj][3];
      viscj[0][1] = viscel[jj][3];
      viscj[2][1] = viscel[jj][4];
      viscj[1][2] = viscel[jj][4];
      viscj[2][0] = viscel[jj][5];
      viscj[0][2] = viscel[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      double fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (int i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*(  viscj[0][i]*i_f_face_normal[face_id][0]
                              + viscj[1][i]*i_f_face_normal[face_id][1]
                              + viscj[2][i]*i_f_face_normal[face_id][2] );
      }

      i_massflux[face_id] =  i_massflux[face_id]
                            + i_visc[face_id]
                             *(  frcxt[ii][0]*( i_face_cog[face_id][0]
                                               -cell_cen[ii][0] )
                               + frcxt[ii][1]*( i_face_cog[face_id][1]
                                               -cell_cen[ii][1] )
                               + frcxt[ii][2]*( i_face_cog[face_id][2]
                                               -cell_cen[ii][2] )
                               - frcxt[jj][0]*( i_face_cog[face_id][0]
                                               -cell_cen[jj][0] )
                               - frcxt[jj][1]*( i_face_cog[face_id][1]
                                               -cell_cen[jj][1] )
                               - frcxt[jj][2]*( i_face_cog[face_id][2]
                                               -cell_cen[jj][2] )
                              )
                            + i_visc[face_id]*ircflp
                             *(- frcxt[ii][0]*diippf[0]
                               - frcxt[ii][1]*diippf[1]
                               - frcxt[ii][2]*diippf[2]
                               + frcxt[jj][0]*djjppf[0]
                               + frcxt[jj][1]*djjppf[1]
                               + frcxt[jj][2]*djjppf[2]
                              );

    }

    /* ---> Contribution from boundary faces */

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      double surfn = b_f_face_surf[face_id]; //FIXME when 0
      double distbf = b_dist[face_id];

      /* FIXME: wrong if dirichlet and viscel is really a tensor */
      b_massflux[face_id] =  b_massflux[face_id]
                            + b_visc[face_id]*distbf/surfn*cofbfp[face_id]
                             *(  frcxt[ii][0]*b_f_face_normal[face_id][0]
                               + frcxt[ii][1]*b_f_face_normal[face_id][1]
                               + frcxt[ii][2]*b_f_face_normal[face_id][2] );

    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
