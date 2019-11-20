/*============================================================================
 * LES Balance
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

/*============================================================================
 * Functions dealing with LES balance
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_divergence.h"
#include "cs_face_viscosity.h"
#include "cs_convection_diffusion.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_geom.h"
#include "cs_gradient.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_properties.h"
#include "cs_prototypes.h"
#include "cs_restart.h"
#include "cs_time_moment.h"
#include "cs_timer_stats.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_les_balance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_les_balance.c
        LES balance computation and related data.
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_les_balance_rij_t

  \brief Reynolds tensor (Rij) LES balance descriptor.

  Members of this LES balance model are publicly accessible, to allow for
  concise syntax, as it is expected to be used in many places.

  \var  cs_les_balance_rij_t::pp2
        Pressure variance \f$ \overline{{p^{\prime}}^2} \f$.
  \var  cs_les_balance_rij_t::smagp2
        Variance of Smagorinsky constant.

  \var  cs_les_balance_rij_t::prodij
  \var  cs_les_balance_rij_t::phiij
  \var  cs_les_balance_rij_t::epsij
  \var  cs_les_balance_rij_t::difftij
  \var  cs_les_balance_rij_t::difftpij
  \var  cs_les_balance_rij_t::unstij
  \var  cs_les_balance_rij_t::convij
  \var  cs_les_balance_rij_t::diffflamij
  \var  cs_les_balance_rij_t::budsgsij
  \var  cs_les_balance_rij_t::budsgsfullij
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_les_balance_tui_t

  \brief Turbulent thermal flux vector (Tui) LES balance descriptor.

  Members of this LES balance model are publicly accessible, to allow for
  concise syntax, as it is expected to be used in many places.

  \var  cs_les_balance_tui_t::unstvar
  \var  cs_les_balance_tui_t::tptp
  \var  cs_les_balance_tui_t::prodvar
  \var  cs_les_balance_tui_t::epsvar
  \var  cs_les_balance_tui_t::difftvar
  \var  cs_les_balance_tui_t::convvar
  \var  cs_les_balance_tui_t::difflamvar
  \var  cs_les_balance_tui_t::budsgsvar
  \var  cs_les_balance_tui_t::tpuip
  \var  cs_les_balance_tui_t::unstti
  \var  cs_les_balance_tui_t::prodtUi
  \var  cs_les_balance_tui_t::prodtTi
  \var  cs_les_balance_tui_t::phiti
  \var  cs_les_balance_tui_t::epsti
  \var  cs_les_balance_tui_t::difftti
  \var  cs_les_balance_tui_t::diffttpi
  \var  cs_les_balance_tui_t::convti
  \var  cs_les_balance_tui_t::difflamti
  \var  cs_les_balance_tui_t::budsgstui
  \var  cs_les_balance_tui_t::budsgsvarfull
  \var  cs_les_balance_tui_t::budsgstuifull
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_les_balance_t

  \brief LES balance general options descriptor.


  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_les_balance_t::brij
  \var  cs_les_balance_t::btui
  \var  cs_les_balance_t::i_les_balance
  \var  cs_les_balance_t::type
  \var  cs_les_balance_t::frequency_n
*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define DEBUG_LES 0

/*============================================================================
 * Static variables
 *============================================================================*/

/* Tensor directions */

static const int idirtens[6][2] = {{0,0},
                                   {1,1},
                                   {2,2},
                                   {0,1},
                                   {0,2},
                                   {1,2}};

static const int idirtens3[10][3] = {{0,1,2},
                                     {0,0,0},
                                     {0,0,1},
                                     {0,0,2},
                                     {1,1,0},
                                     {1,1,1},
                                     {1,1,2},
                                     {2,2,0},
                                     {2,2,1},
                                     {2,2,2}};

static int ipdirtens[3][3];

static int ipdirtens3[3][3][3];

/* Total number of scalars */
static int nscal = 0;

/* Global static structure _les_balance */
static cs_les_balance_t  _les_balance = {NULL,
                                         NULL,
                                         0,
                                         0,
                                         -1};

static cs_field_t *_gradv   = NULL;
static cs_field_t *_gradnut = NULL;
static cs_field_t **_gradt  = NULL;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Global public structure pointer cs_glob_les_balance */
cs_les_balance_t  *cs_glob_les_balance = &_les_balance;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * Get a time moment field by name.
 *
 * parameters:
 *   name <--  time moment name
 *
 * returns  pointer to a field or NULL.
 *----------------------------------------------------------------------------*/

static cs_field_t *
_les_balance_get_tm_by_name(const char *name)
{
  int n_moments = cs_time_moment_n_moments();
  cs_field_t *f = NULL;

  for (int imom = 0 ; imom < n_moments ; imom++) {
    f = cs_time_moment_get_field(imom);
    if (strcmp(f->name, name) == 0)
      return f;
  }

  return NULL;
}

/*----------------------------------------------------------------------------*
 * Get a time moment label for a given scalar and a name.
 *
 * parameters:
 *   isca   <--  scalar number
 *   name   <--  time moment prefix
 *   buffer <--> time moment label
 *
 * returns a time moment label including the scalar number.
 *----------------------------------------------------------------------------*/

static void
_les_balance_get_tm_label(int         isca,
                          const char *name,
                          char       *buffer)
{
  /* Initializing an empty string with name */
  char cSca[5];

  memset(buffer, '\0', 32);
  memset(cSca, '\0', sizeof(cSca));
  strcpy(buffer, name);

  /* Add the scalar number to this name */
  sprintf(cSca, "_%d", isca);
  strcat(buffer, cSca);
}

/*----------------------------------------------------------------------------*
 * Get a time moment field by scalar id.
 *
 * parameters:
 *   scalar_id <-- scalar id
 *   name      <-- time moment name
 *
 * returns  pointer to a field or NULL.
 *----------------------------------------------------------------------------*/

static cs_field_t *
_les_balance_get_tm_by_scalar_id(int scalar_id,
                                 const char *name)
{
  cs_field_t *f = NULL;
  char *buffer;

  BFT_MALLOC(buffer, 32, char);

  _les_balance_get_tm_label(scalar_id, name, buffer);
  f = _les_balance_get_tm_by_name((const char *)buffer);

  BFT_FREE(buffer);

  return f;
}

/*----------------------------------------------------------------------------*
 * Compute the Laplacian of a scalar.
 *
 * parameters:
 *   wa   <--  scalar array
 *   res  <->  Laplacian of wa
 *   type <--  called for Rij (0) or Tui (1) LES balance
 *----------------------------------------------------------------------------*/

static void
_les_balance_laplacian(cs_real_t   *wa,
                       cs_real_t   *res,
                       int         type)
{
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const int *bc_type = cs_glob_bc_type;

  cs_real_t *coefaf, *coefbf;
  BFT_MALLOC(coefaf, n_b_faces, cs_real_t);
  BFT_MALLOC(coefbf, n_b_faces, cs_real_t);

  const cs_real_t visc = 1., pimp = 0., qimp = 0., hext = -1;
  cs_real_t a, b;

  for (cs_lnum_t face_id = 0 ; face_id < n_b_faces ; face_id++) {
    cs_real_t hint = visc / fvq->b_dist[face_id];

    if (   type == 0
        && (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) )
      cs_boundary_conditions_set_dirichlet_scalar(&a,
                                                  &coefaf[face_id],
                                                  &b,
                                                  &coefbf[face_id],
                                                  pimp,
                                                  hint,
                                                  hext);
    else
      cs_boundary_conditions_set_neumann_scalar(&a,
                                                &coefaf[face_id],
                                                &b,
                                                &coefbf[face_id],
                                                qimp,
                                                hint);
  }

  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);
  var_cal_opt.iconv = 0; /* only diffusion */
  var_cal_opt.thetav = 1.;

  cs_real_t *c_visc, *i_visc, *b_visc;
  BFT_MALLOC(c_visc, n_cells_ext, cs_real_t);
  for (cs_lnum_t c_id = 0 ; c_id < n_cells ; c_id++)
    c_visc[c_id] = visc;

  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);
  cs_face_viscosity(m,
                    fvq,
                    0,      /* mean type */
                    c_visc,
                    i_visc,
                    b_visc);
  BFT_FREE(c_visc);

  cs_convection_diffusion_scalar(0,           /* idtvar */
                                 -1,          /* f_id */
                                 var_cal_opt,
                                 0,           /* icvflb (not used) */
                                 1,           /* inc */
                                 true,        /* recompute cocg */
                                 1,           /* imasac (not used) */
                                 wa,          /* pvar */
                                 NULL,        /* pvara (not used) */
                                 0,           /* icvfli (not used) */
                                 coefaf,      /* coefa (not used) */
                                 coefbf,      /* coefb (not used) */
                                 coefaf,
                                 coefbf,
                                 i_visc,      /* mass flux (not used) */
                                 b_visc,      /* mass flux (not used) */
                                 i_visc,
                                 b_visc,
                                 res);

  BFT_FREE(coefaf);
  BFT_FREE(coefbf);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);

  for (cs_lnum_t c_id = 0 ; c_id < n_cells ; c_id++)
    res[c_id] /= cs_glob_mesh_quantities->cell_f_vol[c_id];
}

/*----------------------------------------------------------------------------
 * Compute the divergence of a cell-based vector.
 *
 * Computation options are the ones of the velocity
 *
 * parameters:
 *   wa   <--  vector array
 *   res  <->  divergence of wa
 *
 *----------------------------------------------------------------------------*/

 static void
 _les_balance_divergence_vector(cs_real_3_t  *wa,
                                cs_real_t    *res)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const int *bc_type = cs_glob_bc_type;
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  int f_id, itypfl, iflmb0, init, inc;

  cs_real_t *i_massflux, *b_massflux;
  cs_real_3_t *coefav;
  cs_real_33_t *coefbv;

  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);

  BFT_MALLOC(coefbv, n_b_faces, cs_real_33_t);
  BFT_MALLOC(coefav, n_b_faces, cs_real_3_t);
  BFT_MALLOC(i_massflux, n_i_faces, cs_real_t);
  BFT_MALLOC(b_massflux, n_b_faces, cs_real_t);

  f_id = -1;
  itypfl = 0;
  iflmb0 = 1;
  init = 1;
  inc = 1;

  /* Bc coeffs */
  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
    for (int ii = 0 ; ii < 3 ; ii++) {
        coefav[ifac][ii] = 0.;
        for (int pp = 0 ; pp < 3 ; pp++) {
          if (bc_type[ifac] == CS_SMOOTHWALL
           || bc_type[ifac] == CS_ROUGHWALL)
            coefbv[ifac][ii][pp] = 0.;
          else
            coefbv[ifac][ii][pp] = 1.;
      }
    }
  }

  cs_mass_flux(m,
               mq,
               f_id,
               itypfl,
               iflmb0,
               init,
               inc,
               var_cal_opt.imrgra,
               var_cal_opt.nswrgr,
               var_cal_opt.imligr,
               var_cal_opt.iwarni,
               var_cal_opt.epsrgr,
               var_cal_opt.climgr,
               NULL,
               NULL,
               (const cs_real_3_t *)wa,
               (const cs_real_3_t *)coefav,
               (const cs_real_33_t *)coefbv,
               i_massflux,
               b_massflux);

  cs_divergence(m,
                init,
                i_massflux,
                b_massflux,
                res);

  BFT_FREE(coefav);
  BFT_FREE(coefbv);
  BFT_FREE(i_massflux);
  BFT_FREE(b_massflux);
}

/*----------------------------------------------------------------------------
 *
 *  Compute the most needed gradients at each iteration.
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_gradients(void)
{
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const int *bc_type = cs_glob_bc_type;

  cs_real_t *coefas, *coefbs;

  cs_var_cal_opt_t var_cal_opt;
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;

  /* Computation of the velocity gradient */
  cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  bool use_previous_t = false;
  int inc = 1;

  bft_printf_flush();
  cs_field_gradient_vector(CS_F_(vel),
                           use_previous_t,
                           inc,
                           (cs_real_33_t *)_gradv->val);

  /* Computation of the nu_t gradient */
  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL ||
      _les_balance.type & CS_LES_BALANCE_TUI_FULL) {

    BFT_MALLOC(coefas, n_b_faces, cs_real_t);
    BFT_MALLOC(coefbs, n_b_faces, cs_real_t);

    /* Bc coeffs */
    for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
      coefas[ifac] = 0.;
      if (bc_type[ifac] == CS_SMOOTHWALL
       || bc_type[ifac] == CS_ROUGHWALL)
        coefbs[ifac] = 0.;
      else
        coefbs[ifac] = 1.;
    }

    cs_gradient_scalar("nu_t",
                       gradient_type,
                       halo_type,
                       inc,
                       true,
                       var_cal_opt.nswrgr,
                       0,
                       0,
                       1,
                       var_cal_opt.iwarni,
                       var_cal_opt.imligr,
                       var_cal_opt.epsrgr,
                       var_cal_opt.extrag,
                       var_cal_opt.climgr,
                       NULL,
                       coefas,
                       coefbs,
                       CS_F_(mu_t)->val,
                       NULL,
                       NULL,
                       (cs_real_3_t *)_gradnut->val);

    BFT_FREE(coefas);
    BFT_FREE(coefbs);
  }

  /* Computation of scalar gradients */
  if (_les_balance.type & CS_LES_BALANCE_TUI) {

    const int keysca = cs_field_key_id("scalar_id");
    int iii = 0;

    for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int isca = cs_field_get_key_int(f, keysca);
      if (isca > 0) {
        cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

        cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                                   &gradient_type,
                                   &halo_type);

        cs_field_gradient_scalar(f,
                                 false, /* use_previous_t */
                                 inc,
                                 true, /* _recompute_cocg */
                                 (cs_real_3_t *)_gradt[iii]->val);
         iii++;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes the pressure times the velocity gradient.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_pdjuisym(const void   *input,
                              cs_real_t    *vals)

{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    cs_real_t pre = CS_F_(p)->val[iel];
    for (int ii = 0 ; ii < 6 ; ii++) {
      int i = idirtens[ii][0];
      int j = idirtens[ii][1];
      vals[6*iel + ii] = pre*(grdv[iel][i][j] + grdv[iel][j][i]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes smag
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

inline static void
_les_balance_compute_smag(const void   *input,
                          cs_real_t    *vals)

{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_t *cpro_smago = cs_field_by_name("smagorinsky_constant^2")->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
    vals[iel] = cs_math_sq(cpro_smago[iel]);
}

/*----------------------------------------------------------------------------
 * Function which computes dkui+dkuj.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_dkuidkuj(const void   *input,
                              cs_real_t    *vals)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {

    for (int ii = 0 ; ii < 6 ; ii++)
      vals[6*iel + ii] = 0.;

    for (int ii = 0 ; ii < 6 ; ii++) {
      int i = idirtens[ii][0];
      int j = idirtens[ii][1];
      for (int k = 0 ; k < 3 ; k++)
        vals[6*iel + ii] += grdv[iel][i][k] + grdv[iel][j][k];
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes uidtaujkdxk.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_uidktaujk(const void   *input,
                                 cs_real_t    *vals)
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_t *diverg;
  cs_real_3_t *vel, *velocity;

  velocity = (cs_real_3_t *)CS_F_(vel)->val;

  BFT_MALLOC(diverg, n_cells_ext, cs_real_t);
  BFT_MALLOC(vel, n_cells, cs_real_3_t);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (int i = 0 ; i < 3 ; i++) {

    for (int j = 0 ; j < 3 ; j++) {

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int k = 0 ; k < 3 ; k++)
          vel[iel][k] = -CS_F_(mu_t)->val[iel]*( grdv[iel][j][k]
                                                +grdv[iel][k][j]);

      _les_balance_divergence_vector(vel, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        vals[9*iel+i*3+j] = velocity[iel][i]*diverg[iel];
    }
  }

  BFT_FREE(diverg);
  BFT_FREE(vel);
}

/*----------------------------------------------------------------------------
 * Function which computes nutdkuidkuj.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_nutdkuidkuj(const void   *input,
                                 cs_real_t    *vals)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {

    for (int ii = 0 ; ii < 6 ; ii++)
      vals[6*iel + ii] = 0.;

    for (int ii = 0 ; ii < 6 ; ii++) {
      int i = idirtens[ii][0];
      int j = idirtens[ii][1];
      for (int k = 0 ; k < 3 ; k++)
        vals[6*iel + ii] +=
          CS_F_(mu_t)->val[iel]*grdv[iel][i][k]*grdv[iel][j][k];
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes dknutuidjuksym.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_dknutuidjuksym(const void   *input,
                                    cs_real_t    *vals)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  int i, j;

  CS_UNUSED(input);

  cs_real_3_t *vel;
  vel = (cs_real_3_t *)CS_F_(vel)->val;

  cs_real_33_t *grdv  = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {

    for (int ii = 0 ; ii < 6 ; ii++)
      vals[6*iel + ii] = 0.;

    for (int ii = 0 ; ii < 6 ; ii++) {
      i = idirtens[ii][0];
      j = idirtens[ii][1];
      for (int k = 0 ; k < 3 ; k++)
        vals[6*iel + ii] += grdnu[iel][i]
                            *( vel[iel][i]*grdv[iel][k][j]
                              +vel[iel][j]*grdv[iel][k][i]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes nutdkuiuj.
 *
 * parameters:
 *   input <-- pointer to simple data array
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_nutdkuiuj(const void   *input,
                               cs_real_t    *vals)
{
  const int *k = (const int *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  int i,j ;

  cs_real_3_t *vel   = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_6_t *tens  = (cs_real_6_t *)vals;
  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int ii = 0 ; ii < 6 ; ii++) {
      i = idirtens[ii][0];
      j = idirtens[ii][1];
      tens[iel][ii] = CS_F_(mu_t)->val[iel]
                      *( vel[iel][i]*grdv[iel][j][*k]
                        +vel[iel][j]*grdv[iel][i][*k]);
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes dknutdiuk.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_dknutdiuk(const void   *input,
                               cs_real_t    *vals)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {

    for (int i = 0 ; i < 3 ; i++)
      vals[3*iel + i] = 0.;

    for (int i = 0 ; i < 3 ; i++)
      for (int j = 0 ; j < 3 ; j++)
        vals[3*iel + i] += grdnu[iel][i]*grdv[iel][i][j];
  }
}

/*----------------------------------------------------------------------------
 * Function which computes uidjnut.
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_uidjnut(const void   *input,
                               cs_real_t    *vals)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  CS_UNUSED(input);

  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;
  cs_real_3_t *vel   = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_33_t *tens = (cs_real_33_t *)vals;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int i = 0 ; i < 3 ; i++)
      for (int j = 0 ; j < 3 ; j++)
        tens[iel][i][j] = vel[iel][i]*grdnu[iel][j];
  }
}

/*----------------------------------------------------------------------------
 * Function which computes djtdjui
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_djtdjui(const void   *input,
                             cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const int keysca = cs_field_key_id("scalar_id");
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t dtdxjduidxj;
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int i = 0 ; i < 3 ; i++) {
      dtdxjduidxj = 0.;
      for (int kk = 0 ; kk < 3 ; kk++)
        dtdxjduidxj += grdt[iel][kk]*grdv[iel][i][kk];

      vals[3*iel+i] = dtdxjduidxj;
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes tuiuj
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_tuiuj(const void   *input,
                           cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  int i, j;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int ii = 0 ; ii < 6 ; ii++) {
      i = idirtens[ii][0];
      j = idirtens[ii][1];
      vals[6*iel+ii] = sca->val[iel]*vel[iel][i]*vel[iel][j];
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes uidjt
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_uidjt(const void   *input,
                           cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_33_t *tens = (cs_real_33_t *)vals;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) 
    for (int i = 0 ; i < 3 ; i++) 
      for (int j = 0 ; j < 3 ; j++) 
        tens[iel][i][j] = vel[iel][i]*grdt[iel][j];
}

/*----------------------------------------------------------------------------
 * Function which computes ditdit
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_ditdit(const void   *input,
                                cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_id");
  cs_real_t dtdxidtdxi;
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    dtdxidtdxi = 0.;
    for (int i = 0 ; i < 3 ; i++)
      dtdxidtdxi += grdt[iel][i]; // FIXME: missing SQUARE??
    vals[iel] = dtdxidtdxi;
  }
}

/*----------------------------------------------------------------------------
 * Function which computes tdjtauij
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_tdjtauij(const void   *input,
                              cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const int keysca = cs_field_key_id("scalar_id");
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_t *diverg;
  cs_real_3_t *w1;

  BFT_MALLOC(diverg, n_cells_ext, cs_real_t);
  BFT_MALLOC(w1, n_cells, cs_real_3_t);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (int i = 0 ; i < 3 ; i++) {
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int k = 0 ; k < 3 ; k++)
        w1[iel][k] = -CS_F_(mu_t)->val[iel]*( grdv[iel][i][k]
                                             +grdv[iel][k][i]);

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      vals[3*iel+i] = sca->val[iel]*diverg[iel];
  }

  BFT_FREE(diverg);
  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------
 * Function which computes uidivturflux
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_uidivturflux(const void   *input,
                                  cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const int keysca = cs_field_key_id("scalar_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_t sigmas, *diverg;
  cs_real_3_t *w1, *vel;

  vel = (cs_real_3_t *)CS_F_(vel)->val;
  sigmas = cs_field_get_key_double(sca, ksigmas);

  BFT_MALLOC(diverg, n_cells_ext, cs_real_t);
  BFT_MALLOC(w1, n_cells, cs_real_3_t);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (int i = 0 ; i < 3 ; i++) {
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int k = 0 ; k < 3 ; k++)
        w1[iel][k] =  cs_math_sq(CS_F_(mu_t)->val[iel])/sigmas
                     *(grdv[iel][i][k]+grdv[iel][k][i]);

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      vals[3*iel+i] = vel[iel][i]*diverg[iel];
  }

  BFT_FREE(diverg);
  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------
 * Function which computes tdivturflux
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_tdivturflux(const void   *input,
                                 cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const int keysca = cs_field_key_id("scalar_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_t sigmas, *diverg;
  cs_real_3_t *w1;

  sigmas = cs_field_get_key_double(sca, ksigmas);

  BFT_MALLOC(diverg, n_cells_ext, cs_real_t);
  BFT_MALLOC(w1, n_cells, cs_real_3_t);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  /* TODO : bug dans le fortran, boucle sur ii ? */
  for (int i = 0 ; i < 3 ; i++)
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int k = 0 ; k < 3 ; k++)
        w1[iel][k] =  cs_math_sq(CS_F_(mu_t)->val[iel])/sigmas
                     *(grdv[iel][i][k]+grdv[iel][k][i]);

  _les_balance_divergence_vector(w1, diverg);

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
    vals[iel] = sca->val[iel]*diverg[iel];

  BFT_FREE(diverg);
  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------
 * Function which computes nutdtdxidtdxi
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_nutditdit(const void   *input,
                               cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_id");
  cs_real_t nutditdit;
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    nutditdit = 0.;
    for (int i = 0 ; i < 3 ; i++)
      nutditdit +=  CS_F_(mu_t)->val[iel]*sca->val[iel]
                       *cs_math_sq(grdt[iel][i]);

    vals[iel] = nutditdit;
  }
}

/*----------------------------------------------------------------------------
 * Function which computes nutuidjt
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_nutuidjt(const void   *input,
                              cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_33_t *tens = (cs_real_33_t *)vals;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) 
    for (int i = 0 ; i < 3 ; i++) 
      for (int j = 0 ; j < 3 ; j++) 
        tens[iel][i][j] = CS_F_(mu_t)->val[iel]*vel[iel][i]*grdt[iel][j];
}

/*----------------------------------------------------------------------------
 * Function which computes nutdjuidjt
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_nutdjuidjt(const void   *input,
                                cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_id");
  cs_real_t nutdjuidjt;
  int isca = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
      isca++;
    }
  }

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int i = 0 ; i < 3 ; i++) {
      nutdjuidjt = 0.;
      for (int kk = 0 ; kk < 3 ; kk++)
        nutdjuidjt += CS_F_(mu_t)->val[iel]*grdv[iel][i][kk]*grdt[iel][kk];

      vals[3*iel+i] = nutdjuidjt;
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes djnutdiuj
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_djnutdiuj(const void   *input,
                               cs_real_t    *vals)
{
  CS_UNUSED(input);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t djnutdiuj;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int i = 0 ; i < 3 ; i++) {
      djnutdiuj = 0.;
      for (int kk = 0 ; kk < 3 ; kk++)
        djnutdiuj += grdnu[iel][kk]*grdv[iel][kk][i];

      vals[3*iel+i] = djnutdiuj;
    }
  }
}

/*----------------------------------------------------------------------------
 * Function which computes djnuttdiuj
 *
 * parameters:
 *   input <-- pointer to simple data array (ignored here)
 *   vals  --> pointer to values (size: n_local elements*dimension)
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_djnuttdiuj(const void   *input,
                                cs_real_t    *vals)
{
  const cs_field_t *sca = (const cs_field_t *)input;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t djnuttdiuj;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int i = 0 ; i < 3 ; i++) {
      djnuttdiuj = 0.;
      for (int kk = 0 ; kk < 3 ; kk++)
        djnuttdiuj += grdnu[iel][kk]*sca->val[iel]*grdv[iel][kk][i];

      vals[3*iel+i] = djnuttdiuj;
    }
  }
}

/*----------------------------------------------------------------------------
 *
 * Declare generic time moments for either the Rij or the Tui LES balance.
 * Time moments are defined either by field ids or by function.
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_time_moment(void)
{
  /* Define time moments for Rij balance */
  {
    /* ui */
    int moment_f_id[] = {CS_F_(vel)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("ui_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* p */
    int moment_f_id[] = {CS_F_(p)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("p_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* djui */
    int moment_f_id[] = {_gradv->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;

    cs_time_moment_define_by_field_ids("djui_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    int moment_f_id[] = {CS_F_(vel)->id, CS_F_(vel)->id};
    int moment_c_id[] = {-1, -1};
    int n_fields = 2;

    /* uiuj mean */
    cs_time_moment_define_by_field_ids("uiuj_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,
                                       -1,
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
   
  }

  {
    int moment_f_id[] = {CS_F_(vel)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    
    /* vect(u) variance */
    cs_time_moment_define_by_field_ids("u_v",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_VARIANCE,
                                       1,
                                       -1,
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* nut mean */
    int moment_f_id[] = {CS_F_(mu_t)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("nut_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* nutdjui mean */
    int moment_f_id[] = {CS_F_(mu_t)->id, _gradv->id};
    int moment_c_id[] = {-1, -1};
    int n_fields = 2;
    cs_time_moment_define_by_field_ids("nutdjui_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  if (   _les_balance.type & CS_LES_BALANCE_RIJ_FULL
      || _les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    {
      /* nutui mean */
      int moment_f_id[] = {CS_F_(mu_t)->id, CS_F_(vel)->id};
      int moment_c_id[] = {-1, -1};
      int n_fields = 2;
      cs_time_moment_define_by_field_ids("nutui_m",
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id,
                                         CS_TIME_MOMENT_MEAN,
                                         1,    /* nt_start */
                                         -1,   /* t_start */
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
    }
    {
      /* dinut mean */
      int moment_f_id[] = {_gradnut->id};
      int moment_c_id[] = {-1};
      int n_fields = 1;
      cs_time_moment_define_by_field_ids("dinut_m",
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id,
                                         CS_TIME_MOMENT_MEAN,
                                         1,    /* nt_start */
                                         -1,   /* t_start */
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
    }
  }
}

/*----------------------------------------------------------------------------
 *
 * Declare time moments that can be define either by field ids or by function.
 * for the Rij LES balance
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_time_moment_rij(void)
{

  {
    /* p variance */
    int moment_f_id[] = {CS_F_(p)->id};
    int moment_c_id[] = {-1};
    int n_fields = 1;
    cs_time_moment_define_by_field_ids("p_v",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_VARIANCE,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* pui mean */
    int moment_f_id[] = {CS_F_(p)->id, CS_F_(vel)->id};
    int moment_c_id[] = {-1, -1};
    int n_fields = 2;
    cs_time_moment_define_by_field_ids("pu_m",
                                       n_fields,
                                       moment_f_id,
                                       moment_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       1,    /* nt_start */
                                       -1,   /* t_start */
                                       CS_TIME_MOMENT_RESTART_AUTO,
                                       NULL);
  }

  {
    /* uiujuk mean */
    int moment_f_id[] = {CS_F_(vel)->id, CS_F_(vel)->id, CS_F_(vel)->id};
    int n_fields = 3;

    const char *name[] = { "u1u2u3_m",
                           "u1u1u1_m",
                           "u1u1u2_m",
                           "u1u1u3_m",
                           "u2u2u1_m",
                           "u2u2u2_m",
                           "u2u2u3_m",
                           "u3u3u1_m",
                           "u3u3u2_m",
                           "u3u3u3_m" };

    int moment_c_id[][3] = { {0, 1, 2},
                             {0, 0, 0},
                             {0, 0, 1},
                             {0, 0, 2},
                             {1, 1, 0},
                             {1, 1, 1},
                             {1, 1, 2},
                             {2, 2, 0},
                             {2, 2, 1},
                             {2, 2, 2} };

    for (int ii = 0; ii < 10; ii++) {
      cs_time_moment_define_by_field_ids(name[ii],
                                         n_fields,
                                         moment_f_id,
                                         moment_c_id[ii],
                                         CS_TIME_MOMENT_MEAN,
                                         1,    /* nt_start */
                                         -1,   /* t_start */
                                         CS_TIME_MOMENT_RESTART_AUTO,
                                         NULL);
    }
  }

  /* p(djui+diuj) mean */
  cs_time_moment_define_by_func("pdjuisym_m",
                                CS_MESH_LOCATION_CELLS,
                                6,
                                _les_balance_compute_pdjuisym,
                                NULL,
                                NULL,
                                NULL,
                                CS_TIME_MOMENT_MEAN,
                                1,
                                -1,
                                CS_TIME_MOMENT_RESTART_AUTO,
                                NULL);

  /* dk(ui+uj) mean */
  cs_time_moment_define_by_func("dkuidkuj_m",
                                CS_MESH_LOCATION_CELLS,
                                6,
                                _les_balance_compute_dkuidkuj,
                                NULL,
                                NULL,
                                NULL,
                                CS_TIME_MOMENT_MEAN,
                                1,
                                -1,
                                CS_TIME_MOMENT_RESTART_AUTO,
                                NULL);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {

    if (cs_glob_turb_model->iturb == 41) {

      /* smagorinsky variance */
      cs_time_moment_define_by_func("smag_v",
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    _les_balance_compute_smag,
                                    NULL,
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_VARIANCE,
                                    1,
                                      -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }

    /* uidktaujk */
    cs_time_moment_define_by_func("uidktaujk_m",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  _les_balance_compute_uidktaujk,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {

    {
      /* uidkuj mean */
      int moment_f_id[] = {CS_F_(vel)->id, _gradv->id};
      int n_fields = 2;

      const char *name[] = { "u1dkuj_m",
                             "u2dkuj_m",
                             "u3dkuj_m" };

      int moment_c_id[][2] = { {0, -1},
                               {1, -1},
                               {2, -1} };

      for (int ii = 0; ii < 3; ii++) {
        cs_time_moment_define_by_field_ids(name[ii],
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id[ii],
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }
    }

    /* nutdkuidkuj mean */
    cs_time_moment_define_by_func("nutdkuidkuj_m",
                                  CS_MESH_LOCATION_CELLS,
                                  6,
                                  _les_balance_compute_nutdkuidkuj,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);

    /* dknutuidjuksym mean */
    cs_time_moment_define_by_func("dknutuidjuksym_m",
                                  CS_MESH_LOCATION_CELLS,
                                  6,
                                  _les_balance_compute_dknutuidjuksym,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);

    /* nutdkuiuj mean */
    {
      const char *name[] = { "nutd1uiuj_m",
                             "nutd2uiuj_m",
                             "nutd3uiuj_m" };

      static int index[3] = {0, 1, 2};

      for (int k = 0; k < 3; k++) {
        cs_time_moment_define_by_func(name[k],
                                      CS_MESH_LOCATION_CELLS,
                                      6,
                                      _les_balance_compute_nutdkuiuj,
                                      &(index[k]),
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);
      }
    }

    /* dknutdiuk mean */
    cs_time_moment_define_by_func("dknutdiuk_m",
                                  CS_MESH_LOCATION_CELLS,
                                  3,
                                  _les_balance_compute_dknutdiuk,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);

    {
      cs_time_moment_define_by_func("uidjnut_m",
                                    CS_MESH_LOCATION_CELLS,
                                    9,
                                    _les_balance_compute_uidjnut,
                                    NULL,
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
    }
  }
}

/*----------------------------------------------------------------------------
 *
 * Declare time moments that can be define either by field ids or by function.
 * for the Tui LES balance
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_time_moment_tui(void)
{
  const int keysca = cs_field_key_id("scalar_id");
  char *buffer;
  int isca = 0;

  BFT_MALLOC(buffer, 32, char);

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    /* _djnutdiuj */
    cs_time_moment_define_by_func("djnutdiuj_m",
                                  CS_MESH_LOCATION_CELLS,
                                  3,
                                  _les_balance_compute_djnutdiuj,
                                  NULL,
                                  NULL,
                                  NULL,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  NULL);
  }

  /* Define time moments for Tui balance */
  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int iscal = cs_field_get_key_int(f, keysca)-1;
    if (iscal > -1) {
      {
        /* t mean */
        int moment_f_id[] = {f_id};
        int moment_c_id[] = {-1};
        int n_fields = 1;
        _les_balance_get_tm_label(isca, "t_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* t variance */
        int moment_f_id[] = {f_id};
        int moment_c_id[] = {-1};
        int n_fields = 1;
        _les_balance_get_tm_label(isca, "t_v", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_VARIANCE,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* dtdxi mean */
        int moment_f_id[] = {_gradt[iscal]->id};
        int moment_c_id[] = {-1};
        int n_fields = 1;
        _les_balance_get_tm_label(isca, "dtdxi_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* pdtdxi mean */
        int moment_f_id[] = {_gradt[iscal]->id, CS_F_(p)->id};
        int moment_c_id[] = {-1, -1};
        int n_fields = 2;
        _les_balance_get_tm_label(isca, "pdtdxi_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* tui mean */
        int moment_f_id[] = {f_id, CS_F_(vel)->id};
        int moment_c_id[] = {-1, -1};
        int n_fields = 2;
        _les_balance_get_tm_label(isca, "tui_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* tui variance */
        int moment_f_id[] = {f_id, CS_F_(vel)->id};
        int moment_c_id[] = {-1, -1};
        int n_fields = 2;
        _les_balance_get_tm_label(isca, "tui_v", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_VARIANCE,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      /* djtdjui mean */
      _les_balance_get_tm_label(isca, "djtdjui_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    _les_balance_compute_djtdjui,
                                    f,
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
      /* tuiuj mean */
      _les_balance_get_tm_label(isca, "tuiuj_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    6,
                                    _les_balance_compute_tuiuj,
                                    f,
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);

      {
        /* tp mean */
        int moment_f_id[] = {f_id, CS_F_(p)->id};
        int moment_c_id[] = {-1, -1};
        int n_fields = 2;
        _les_balance_get_tm_label(isca, "tp_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        /* tdjui mean */
        int moment_f_id[] = {f_id, _gradv->id};
        int moment_c_id[] = {-1, -1};
        int n_fields = 2;
        _les_balance_get_tm_label(isca, "tdjui_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      {
        _les_balance_get_tm_label(isca, "uidjt_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      9,
                                      _les_balance_compute_uidjt,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);
      }

      /* ditdit mean */
      _les_balance_get_tm_label(isca, "ditdit_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    _les_balance_compute_ditdit,
                                    f,
                                    NULL,
                                    NULL,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    NULL);
      {
        /* ttui mean */
        int moment_f_id[] = {f_id, f_id, CS_F_(vel)->id};
        int moment_c_id[] = {-1, -1, -1};
        int n_fields = 3;
        _les_balance_get_tm_label(isca, "ttui_m", buffer);
        cs_time_moment_define_by_field_ids(buffer,
                                           n_fields,
                                           moment_f_id,
                                           moment_c_id,
                                           CS_TIME_MOMENT_MEAN,
                                           1,    /* nt_start */
                                           -1,   /* t_start */
                                           CS_TIME_MOMENT_RESTART_AUTO,
                                           NULL);
      }

      if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
        /* _tdjtauij mean */
        _les_balance_get_tm_label(isca, "tdjtauij_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      _les_balance_compute_tdjtauij,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);

        /* _uidivturflux mean */
        _les_balance_get_tm_label(isca, "uidivturflux_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      _les_balance_compute_uidivturflux,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);

        /* _tdivturflux mean */
        _les_balance_get_tm_label(isca, "tdivturflux_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      _les_balance_compute_tdivturflux,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);
      }

      if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
        {
          /* nutdit mean */
          int moment_f_id[] = {CS_F_(mu_t)->id, _gradt[iscal]->id};
          int moment_c_id[] = {-1, -1};
          int n_fields = 2;
          _les_balance_get_tm_label(isca, "nutdit_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

        /* nutditdit mean */
        _les_balance_get_tm_label(isca, "nutditdit_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      _les_balance_compute_nutditdit,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);

        {
          /* nuttdjui mean */
          int moment_f_id[] = {CS_F_(mu_t)->id, f_id, _gradv->id};
          int moment_c_id[] = {-1, -1, -1};
          int n_fields = 3;
          _les_balance_get_tm_label(isca, "nuttdjui_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

        {
          /* nutuidjt mean */
          _les_balance_get_tm_label(isca, "nutuidjt_m", buffer);
          cs_time_moment_define_by_func(buffer,
                                        CS_MESH_LOCATION_CELLS,
                                        9,
                                        _les_balance_compute_nutuidjt,
                                        f,
                                        NULL,
                                        NULL,
                                        CS_TIME_MOMENT_MEAN,
                                        1,
                                        -1,
                                        CS_TIME_MOMENT_RESTART_AUTO,
                                        NULL);
        }

        /* nutdjuidjt mean */
        _les_balance_get_tm_label(isca, "nutdjuidjt_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      _les_balance_compute_nutdjuidjt,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);

        {
          /* tdit mean */
          int moment_f_id[] = {f_id, _gradt[iscal]->id};
          int moment_c_id[] = {-1, -1};
          int n_fields = 2;
          _les_balance_get_tm_label(isca, "tdit_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

        {
          /* nuttdit mean */
          int moment_f_id[] = {CS_F_(mu_t)->id, f_id, _gradt[iscal]->id};
          int moment_c_id[] = {-1, -1, -1};
          int n_fields = 3;
          _les_balance_get_tm_label(isca, "nuttdit_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

        {
          /* nutt */
          int moment_f_id[] = {CS_F_(mu_t)->id, f_id};
          int moment_c_id[] = {-1, -1};
          int n_fields = 2;
          _les_balance_get_tm_label(isca, "nutt_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

        /* djnuttdiuj mean */
        _les_balance_get_tm_label(isca, "djnuttdiuj_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      _les_balance_compute_djnuttdiuj,
                                      f,
                                      NULL,
                                      NULL,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      NULL);
        {
          /* tdinut mean */
          int moment_f_id[] = {f_id, _gradnut->id};
          int moment_c_id[] = {-1, -1};
          int n_fields = 2;
          _les_balance_get_tm_label(isca, "tdinut_m", buffer);
          cs_time_moment_define_by_field_ids(buffer,
                                             n_fields,
                                             moment_f_id,
                                             moment_c_id,
                                             CS_TIME_MOMENT_MEAN,
                                             1,    /* nt_start */
                                             -1,   /* t_start */
                                             CS_TIME_MOMENT_RESTART_AUTO,
                                             NULL);
        }

      }
      isca++;
    }
  }

  BFT_FREE(buffer);

}

/*----------------------------------------------------------------------------
 * Create and allocate a cs_les_balance_rij_t structure
 *
 * returns a pointer to the created structure
 *----------------------------------------------------------------------------*/

static cs_les_balance_rij_t *
_les_balance_allocate_rij(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  /* Initialization */
  cs_les_balance_rij_t *brij = NULL;

  BFT_MALLOC(brij, 1, cs_les_balance_rij_t);

  /* Allocation of the working arrays  that cannot
     be declared using time moments */

  BFT_MALLOC(brij->prodij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->epsij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->phiij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->difftij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->difftpij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->convij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->difflamij, n_cells, cs_real_6_t);
  BFT_MALLOC(brij->unstij, n_cells, cs_real_6_t);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) 
    BFT_MALLOC(brij->budsgsij, n_cells, cs_real_6_t);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) 
    BFT_MALLOC(brij->budsgsfullij, n_cells, cs_real_96_t);

  return brij;
}

/*----------------------------------------------------------------------------
 * Create and allocate a cs_les_balance_tui_t structure
 *
 * returns a pointer to the created structure
 *----------------------------------------------------------------------------*/

static cs_les_balance_tui_t *
_les_balance_allocate_tui(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  /* Initialization */
  cs_les_balance_tui_t *btui = NULL;

  BFT_MALLOC(btui, 1, cs_les_balance_tui_t);

  /* Allocation of the working arrays  that cannot
     be declared using time moments */

  /* Working arrays */
  BFT_MALLOC(btui->unstvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->tptp, n_cells, cs_real_t);
  BFT_MALLOC(btui->prodvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->epsvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->difftvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->convvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->difflamvar, n_cells, cs_real_t);
  BFT_MALLOC(btui->tpuip, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->unstti, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->prodtUi, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->prodtTi, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->phiti, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->epsti, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->difftti, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->diffttpi, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->convti, n_cells, cs_real_3_t);
  BFT_MALLOC(btui->difflamti, n_cells, cs_real_3_t);

  if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
    BFT_MALLOC(btui->budsgsvar, n_cells, cs_real_t);
    BFT_MALLOC(btui->budsgstui, n_cells, cs_real_3_t);
  }

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    BFT_MALLOC(btui->budsgsvarfull, n_cells, cs_real_6_t);
    BFT_MALLOC(btui->budsgstuifull, 10, cs_real_3_t *);
    for (int ii = 0 ; ii < 10 ; ii++)
      BFT_MALLOC(btui->budsgstuifull[ii], n_cells, cs_real_3_t);
  }

  return btui;
}

/*----------------------------------------------------------------------------
 *
 * Initialize the brij structure of _les_balance
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_initialize_rij(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_les_balance_rij_t *brij = _les_balance.brij;

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {

    for (int i = 0 ; i < 6 ; i++) {
      brij->prodij[iel][i] = 0.;
      brij->epsij[iel][i] = 0.;
      brij->phiij[iel][i] = 0.;
      brij->difftij[iel][i] = 0.;
      brij->difftpij[iel][i] = 0.;
      brij->convij[iel][i] = 0.;
      brij->difflamij[iel][i] = 0.;
      brij->unstij[iel][i] = 0.;

      if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE)
        brij->budsgsij[iel][i] = 0.;
    }

    if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL)
      for (int i = 0 ; i < 9 ; i++)
        for (int j = 0 ; j < 6 ; j++)
          brij->budsgsfullij[iel][i][j] = 0.;

  }
}

/*----------------------------------------------------------------------------
 *
 * Check that in case of a restart calculation, the current LES balance type
 * is the same as in the previous calculation.
 *
 *----------------------------------------------------------------------------*/

static void
_check_restart_type(void)
{
  int ierror;

  char const *ficsui = "les_balance";
  cs_restart_t *suite = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_READ);

  {
    int itysup = 0;
    int type;

    ierror = cs_restart_read_section(suite,
                                     "les_balance_type",
                                     itysup,
                                     1,
                                     CS_TYPE_cs_int_t,
                                     &type);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Abort while opening the LES balance restart file: %s\n"
                  "This file does not seem to be a LES balance checkpoint file."),
                ficsui);

    if (!(type & _les_balance.type))
      bft_error(__FILE__, __LINE__, 0,
                _("Abort while reading the LES balance restart file: %s\n"
                  "The previous balance type is different from the current\n"
                  "balance type:\n"
                  "  previous type: %d\n"
                  "  current type:  %d\n"),
                  ficsui, type, _les_balance.type);
  }
}

/*----------------------------------------------------------------------------
 *
 * Initialize a btui structure of _les_balance
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_initialize_tui(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_id");
  int iscal = 0;

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      _les_balance.btui[iscal]->f_id = f_id;
      iscal++;
    }
  }

  /* Since every time averages for Tui were created using
     time moments, there is no need to initialize them by
     reading the les_balance restart */
  for (int isca = 0 ; isca < nscal ; isca++) {
    cs_les_balance_tui_t *btui = _les_balance.btui[isca];

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      btui->unstvar[iel] = 0.;
      btui->tptp[iel] = 0.;
      btui->prodvar[iel] = 0.;
      btui->epsvar[iel] = 0.;
      btui->difftvar[iel] = 0.;
      btui->convvar[iel] = 0.;
      btui->difflamvar[iel] = 0.;

      for (int ii = 0 ; ii < 3 ; ii++) {
        btui->unstti[iel][ii] = 0.;
        btui->tpuip[iel][ii] = 0.;
        btui->prodtUi[iel][ii] = 0.;
        btui->prodtTi[iel][ii] = 0.;
        btui->phiti[iel][ii] = 0.;
        btui->epsti[iel][ii] = 0.;
        btui->difftti[iel][ii] = 0.;
        btui->diffttpi[iel][ii] = 0.;
        btui->convti[iel][ii] = 0.;
        btui->difflamti[iel][ii] = 0.;
      }

      if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
        btui->budsgsvar[iel] = 0.;

        for (int ii = 0 ; ii < 3 ; ii++)
          btui->budsgstui[iel][ii] = 0.;
      }

      if (_les_balance.type & CS_LES_BALANCE_TUI_FULL)
        for (int ii = 0 ; ii < 10 ; ii++)
          for (int jj = 0 ; jj < 3 ; jj++)
            btui->budsgstuifull[ii][iel][jj] = 0.;
    }
  }
}

/*----------------------------------------------------------------------------
 * Destroy a given cs_les_balance_rij_t structure pointer
 *
 * returns NULL
 *----------------------------------------------------------------------------*/

static cs_les_balance_rij_t *
_les_balance_destroy_rij(cs_les_balance_rij_t *brij)
{
  if (brij == NULL)
    return brij;

  BFT_FREE(brij->prodij);
  BFT_FREE(brij->epsij);
  BFT_FREE(brij->phiij);
  BFT_FREE(brij->difftij);
  BFT_FREE(brij->difftpij);
  BFT_FREE(brij->convij);
  BFT_FREE(brij->difflamij);

  BFT_FREE(brij->unstij);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {
    BFT_FREE(brij->budsgsij);
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    BFT_FREE(brij->budsgsfullij);
  }

  BFT_FREE(brij);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Destroy an array of cs_les_balance_tui_t structure pointers
 *
 * returns NULL
 *----------------------------------------------------------------------------*/

static cs_les_balance_tui_t **
_les_balance_destroy_tui(cs_les_balance_tui_t **btui)
{
  if (btui == NULL)
    return btui;

  for (int isca = 0 ; isca < nscal ; isca++) {

    BFT_FREE(btui[isca]->unstvar);
    BFT_FREE(btui[isca]->tptp);
    BFT_FREE(btui[isca]->tpuip);
    BFT_FREE(btui[isca]->prodvar);
    BFT_FREE(btui[isca]->epsvar);
    BFT_FREE(btui[isca]->difftvar);
    BFT_FREE(btui[isca]->convvar);
    BFT_FREE(btui[isca]->difflamvar);
    BFT_FREE(btui[isca]->unstti);
    BFT_FREE(btui[isca]->prodtUi);
    BFT_FREE(btui[isca]->prodtTi);
    BFT_FREE(btui[isca]->phiti);
    BFT_FREE(btui[isca]->epsti);
    BFT_FREE(btui[isca]->difftti);
    BFT_FREE(btui[isca]->diffttpi);
    BFT_FREE(btui[isca]->convti);
    BFT_FREE(btui[isca]->difflamti);

    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
      BFT_FREE(btui[isca]->budsgsvar);
      BFT_FREE(btui[isca]->budsgstui);
    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
      for (int ii = 0 ; ii < 10 ; ii++)
        BFT_FREE(btui[isca]->budsgstuifull[ii]);
      BFT_FREE(btui[isca]->budsgstuifull);
      BFT_FREE(btui[isca]->budsgsvarfull);
    }

    BFT_FREE(btui[isca]);
  }

  BFT_FREE(btui);

  return btui;
}

/*----------------------------------------------------------------------------
 *
 * Create, allocate and initialize a the brij structure of _les_balance.
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_create_rij(void)
{
  /* Creation and allocation of the structure containing
     other time averages that can't be time moments
     and other working arrays */
  _les_balance.brij = _les_balance_allocate_rij();

  /* Initialization of working arrays and reading
     restart file in case of restart */
  _les_balance_initialize_rij();
}

/*----------------------------------------------------------------------------
 *
 * Create, allocate and initialize btui structures of _les_balance.
 *
 *----------------------------------------------------------------------------*/

static void
_les_balance_create_tui(void)
{
  BFT_MALLOC(_les_balance.btui, nscal, cs_les_balance_tui_t *);

  /* Creation and allocation of the structure containing
     working arrays */
  for (int isca = 0 ; isca < nscal ; isca++)
    _les_balance.btui[isca] = _les_balance_allocate_tui();

  /* Initialization of working arrays */
  _les_balance_initialize_tui();
}

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_les_balance_get_pointer(int **i_les_balance);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointer to member of the global les balance structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   i_les_balance --> pointer to cs_glob_les_balance->i_les_balance
 *----------------------------------------------------------------------------*/

void
cs_f_les_balance_get_pointer(int **i_les_balance)
{
  *i_les_balance = &(_les_balance.i_les_balance);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *! \brief Create fields used in LES balance computation
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_create_fields(void)
{
  /* Velocity gradient */
  {
    const char *name = "vel_grad";
    int dim = 9;
    _gradv = cs_field_create(name,
                             CS_FIELD_PROPERTY,
                             CS_MESH_LOCATION_CELLS,
                             dim,
                             false);
  }

  /* Nu_t gradient */
  if (   _les_balance.type & CS_LES_BALANCE_RIJ_FULL
      || _les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    const char *name = "nut_grad";
    int dim = 3;
    _gradnut = cs_field_create(name,
                               CS_FIELD_PROPERTY,
                               CS_MESH_LOCATION_CELLS,
                               dim,
                               false);
  }

  if (_les_balance.type & CS_LES_BALANCE_TUI) {
    const int k_sca = cs_field_key_id("scalar_id");

    /* First we allocate _gradt */
    int n_scal = 0;

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f_sca = cs_field_by_id(f_id);
      int i_sca = cs_field_get_key_int(f_sca, k_sca);
      if (i_sca > 0) n_scal ++;
    }

    BFT_MALLOC(_gradt, n_scal, cs_field_t*);

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f_sca = cs_field_by_id(f_id);

      int i_sca = cs_field_get_key_int(f_sca, k_sca)-1;
      if (i_sca > -1) {
        char *name;
        int len = strlen(f_sca->name)+6;
        BFT_MALLOC(name, len, char);
        snprintf(name, len, "%s_grad", f_sca->name);

        int dim = 3;
        _gradt[i_sca] = cs_field_create(name,
                                        CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_CELLS,
                                        dim,
                                        false);
        BFT_FREE(name);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 *! \brief Provide access to cs_glob_les_balance
 *
 * \return pointer to LES balance global structure
 */
/*----------------------------------------------------------------------------*/

cs_les_balance_t *
cs_get_glob_les_balance(void)
{
  return &_les_balance;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a LES balance descriptor.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_create(void)
{
  /* Rij active if Rij basic or full is active */
  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE
   || _les_balance.type & CS_LES_BALANCE_RIJ_FULL)
    _les_balance.type |= CS_LES_BALANCE_RIJ;

  /* Tui active if Tui basic or full is active */
  if (_les_balance.type & CS_LES_BALANCE_TUI_BASE
   || _les_balance.type & CS_LES_BALANCE_TUI_FULL)
    _les_balance.type |= CS_LES_BALANCE_TUI;

  /* Count the number of scalars nscal */
  const int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0 ; f_id < cs_field_n_fields() ; f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int isca = cs_field_get_key_int(f, keysca);
    if (isca > 0) nscal++;
  }

  /* Remplissage des tableaux de directions de tenseurs */
  for (int ii = 0 ; ii < 3 ; ii++)
    for (int jj = 0 ; jj < 3 ; jj++)
      for (int iii = 0 ; iii < 6 ; iii++)
        if (ii*jj == idirtens[iii][0]*idirtens[iii][1])
          ipdirtens[ii][jj] = iii;

  for (int ii = 0 ; ii < 3 ; ii++)
    for (int jj = 0 ; jj < 3 ; jj++)
      for (int kk = 0 ; kk < 3 ; kk++)
        for (int iii = 0 ; iii < 10 ; iii++)
          if (ii*jj*kk == idirtens3[iii][0]*idirtens3[iii][1]*idirtens3[iii][2])
            ipdirtens3[ii][jj][kk] = iii;

  if (cs_restart_present())
    _check_restart_type();

  cs_les_balance_create_fields();

  /* Creation of the generic time moments used for both Rij
     and Tui LES balance */
  _les_balance_time_moment();

  /* If Rij balance is active */
  if (_les_balance.type & CS_LES_BALANCE_RIJ) {

    /* Creation of the time moments for Rij */
    _les_balance_time_moment_rij();

    /* Creation of the brij structure inside _les_balance */
    _les_balance_create_rij();
  }

  /* If Tui balance is active */
  if (_les_balance.type & CS_LES_BALANCE_TUI) {

    /* Creation of the time moments for Tui */
    _les_balance_time_moment_tui();

    /* Creation of the btui structure inside _les_balance */
    _les_balance_create_tui();
  }

  /* Add time moments log in the listing in DEBUG mode */
#if DEBUG_LES == 1
  const int log_key_id = cs_field_key_id("log");
  const int n_moments = cs_time_moment_n_moments();

  for (int m_id = 0; m_id < n_moments; m_id++) {
    cs_field_t *f = cs_time_moment_get_field(m_id);
    if (f != NULL && !cs_field_is_key_set(f, log_key_id))
      cs_field_set_key_int(f, log_key_id, 1);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update gradients needed in LES balance
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_update_gradients(void)
{
  _les_balance_compute_gradients();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Rij LES balance.
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_compute_rij(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);

  /* Rij balance structure */
  cs_les_balance_rij_t *brij = _les_balance.brij;

  cs_real_t *p, *nut;
  cs_real_3_t *ui, *pu, *nutui, *dnutdxkdukdxi, *dnutdxi;
  cs_real_6_t *uiuj, *rij, *pduidxj, *duidxkdujdxk, *nutduidxkdujdxk, *dnutdxkuidukdxjsym;
  cs_real_33_t *nutduidxj, *uidtaujkdxk, *duidxj, *uidnutdxj;

  /* Time moments retrieved by name */
  p            =                 _les_balance_get_tm_by_name("p_m")->val;
  nut          =                 _les_balance_get_tm_by_name("nut_m")->val;
  ui           = (cs_real_3_t  *)_les_balance_get_tm_by_name("ui_m")->val;
  pu           = (cs_real_3_t  *)_les_balance_get_tm_by_name("pu_m")->val;
  uiuj         = (cs_real_6_t  *)_les_balance_get_tm_by_name("uiuj_m")->val;
  rij          = (cs_real_6_t  *)_les_balance_get_tm_by_name("u_v")->val;
  nutduidxj    = (cs_real_33_t *)_les_balance_get_tm_by_name("nutdjui_m")->val;
  pduidxj      = (cs_real_6_t  *)_les_balance_get_tm_by_name("pdjuisym_m")->val;
  duidxj       = (cs_real_33_t *)_les_balance_get_tm_by_name("djui_m")->val;
  duidxkdujdxk = (cs_real_6_t  *)_les_balance_get_tm_by_name("dkuidkuj_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE)
    uidtaujkdxk  = (cs_real_33_t *)_les_balance_get_tm_by_name("uidktaujk_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    nutui              = (cs_real_3_t  *)_les_balance_get_tm_by_name("nutui_m")->val;
    dnutdxkdukdxi      = (cs_real_3_t  *)_les_balance_get_tm_by_name("dknutdiuk_m")->val;
    dnutdxi            = (cs_real_3_t  *)_les_balance_get_tm_by_name("dinut_m")->val;
    nutduidxkdujdxk    = (cs_real_6_t  *)_les_balance_get_tm_by_name("nutdkuidkuj_m")->val;
    dnutdxkuidukdxjsym = (cs_real_6_t  *)_les_balance_get_tm_by_name("dknutuidjuksym_m")->val;
    uidnutdxj          = (cs_real_33_t *)_les_balance_get_tm_by_name("uidjnut_m")->val;
  }

  /* Get the triple corrleations mean UiUjUk*/
  cs_real_t** uiujuk;
  BFT_MALLOC(uiujuk, 10, cs_real_t*);

  uiujuk[0] = cs_field_by_name("u1u2u3_m")->val;
  uiujuk[1] = cs_field_by_name("u1u1u1_m")->val;
  uiujuk[2] = cs_field_by_name("u1u1u2_m")->val;
  uiujuk[3] = cs_field_by_name("u1u1u3_m")->val;
  uiujuk[3] = cs_field_by_name("u2u2u1_m")->val;
  uiujuk[5] = cs_field_by_name("u2u2u2_m")->val;
  uiujuk[6] = cs_field_by_name("u2u2u3_m")->val;
  uiujuk[7] = cs_field_by_name("u3u3u1_m")->val;
  uiujuk[8] = cs_field_by_name("u3u3u2_m")->val;
  uiujuk[9] = cs_field_by_name("u3u3u3_m")->val;

  /* Get additional averaged fields */
  cs_real_6_t** nutdkuiuj;
  cs_real_33_t** uidujdxk;
  BFT_MALLOC(nutdkuiuj, 3, cs_real_6_t*);
  BFT_MALLOC(uidujdxk, 3, cs_real_33_t*);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    nutdkuiuj[0] = cs_field_by_name("nutd1uiuj_m")->val;
    nutdkuiuj[1] = cs_field_by_name("nutd2uiuj_m")->val;
    nutdkuiuj[2] = cs_field_by_name("nutd3uiuj_m")->val;
    
    uidujdxk[0] = cs_field_by_name("u1dkuj_m")->val;
    uidujdxk[1] = cs_field_by_name("u2dkuj_m")->val;
    uidujdxk[2] = cs_field_by_name("u3dkuj_m")->val;
  }

  /* Working arrays */
  cs_real_t *diverg, *w2, *lapl;
  cs_real_t *coefas, *coefbs;
  cs_real_3_t *w1, *w3;

  cs_real_t dtref = cs_glob_time_step->dt_ref;
  cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
  int i, j, jj, ll, jjj, kkk, lll, inc;

  const int *bc_type = cs_glob_bc_type;

  inc = 1;

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  /* Working arrays memory allocation */
  BFT_MALLOC(w1, n_cells, cs_real_3_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(w3, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(coefbs, n_b_faces, cs_real_t);
  BFT_MALLOC(coefas, n_b_faces, cs_real_t);
  BFT_MALLOC(diverg, n_cells_ext, cs_real_t);
  BFT_MALLOC(lapl, n_cells_ext, cs_real_t);

  for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
    for (int iii = 0 ; iii < 6 ; iii++) {
      brij->prodij[iel][iii]    = 0.;
      brij->epsij[iel][iii]     = 0.;
      brij->phiij[iel][iii]     = 0.;
      brij->difftij[iel][iii]   = 0.;
      brij->difftpij[iel][iii]  = 0.;
      brij->convij[iel][iii]    = 0.;
      brij->difflamij[iel][iii] = 0.;
    }
  }

  /* unstij, epsij, prodij, phiij */
  for (int ii = 0 ; ii < 6 ; ii++) {
    i = idirtens[ii][0];
    j = idirtens[ii][1];

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      brij->unstij[iel][ii] = (rij[iel][ii] - brij->unstij[iel][ii])/dtref;
      brij->epsij[iel][ii] = duidxkdujdxk[iel][ii];

      for (int kk = 0 ; kk < 3 ; kk++) {
        jj = ipdirtens[i][kk];
        ll = ipdirtens[j][kk];
        brij->prodij[iel][ii] -= rij[iel][ll]*duidxj[iel][i][kk]
                               + rij[iel][jj]*duidxj[iel][j][kk];
        brij->epsij[iel][ii] -= duidxj[iel][i][kk]*duidxj[iel][j][kk];
      }

      brij->phiij[iel][ii] = (pduidxj[iel][ii] - p[iel]*(duidxj[iel][i][j]+duidxj[iel][j][i]))/ro0;
      brij->epsij[iel][ii] *= -2.*viscl0/ro0;
    }
  }

  /* convij */
  for (int iii = 0 ; iii < 6 ; iii++) {
    i = idirtens[iii][0];
    j = idirtens[iii][1];

    /* convij */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int kk = 0 ; kk < 3 ; kk++)
        w1[iel][kk] = ui[iel][kk]*rij[iel][iii];

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      brij->convij[iel][iii] = diverg[iel];

    /* difftij */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      for (int kk = 0 ; kk < 3 ; kk++) {
        jjj = ipdirtens[i][kk];
        kkk = ipdirtens[j][kk];
        lll = ipdirtens3[i][j][kk];
        cs_real_t *triple_corr = uiujuk[lll];

        w1[iel][kk] = - triple_corr[iel] - 2.*ui[iel][i]*ui[iel][j]*ui[iel][kk]
                                         + ui[iel][i]*uiuj[iel][kkk]
                                         + ui[iel][j]*uiuj[iel][jjj]
                                         + ui[iel][kk]*uiuj[iel][iii];
      }
    }

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      brij->difftij[iel][iii] = diverg[iel];

    /* difftpij */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int kk = 0 ; kk < 3 ; kk++)
        w1[iel][kk] = 0.;

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      for (int kk = 0 ; kk < 3 ; kk++) {
        if (kk == i)
          w1[iel][kk] -= pu[iel][j] - p[iel]*ui[iel][j];
        else if (kk == j)
          w1[iel][kk] -= pu[iel][i] - p[iel]*ui[iel][i];
      }
    }

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      brij->difftpij[iel][iii] = diverg[iel];

    /* Laminar diffusion */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      w2[iel] = rij[iel][iii];

    _les_balance_laplacian(w2, lapl, 0);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      brij->difflamij[iel][iii] = viscl0*lapl[iel]/ro0;
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {

    for (int iii = 0 ; iii < 6 ; iii++) {
      i = idirtens[iii][0];
      j = idirtens[iii][1];

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = -nutduidxj[iel][j][kk] - nutduidxj[iel][kk][j];

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        brij->budsgsij[iel][iii] = -(uidtaujkdxk[iel][i][j]-ui[iel][i]*diverg[iel]);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = -nutduidxj[iel][i][kk] - nutduidxj[iel][kk][i];

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsij[iel][iii] -= (uidtaujkdxk[iel][j][i]-ui[iel][j]*diverg[iel]);
        brij->budsgsij[iel][iii] /= ro0;
      }
    }
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {

    for (int iii = 0 ; iii < 6 ; iii++) {
      i = idirtens[iii][0];
      j = idirtens[iii][1];
      cs_real_33_t *uidujdxk_ii = (cs_real_33_t*)uidujdxk[i];
      cs_real_33_t *uidujdxk_jj = (cs_real_33_t*)uidujdxk[j];

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = nut[iel]*((uidujdxk_ii[iel][j][kk] - ui[iel][i]*duidxj[iel][j][kk])
                                + (uidujdxk_jj[iel][i][kk] - ui[iel][j]*duidxj[iel][i][kk]));

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsfullij[iel][iii][0] = diverg[iel]/ro0;
        brij->budsgsfullij[iel][iii][1] = nut[iel]/viscl0*brij->epsij[iel][iii];
      }

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        for (int kk = 0 ; kk < 3 ; kk++) {
          cs_real_6_t *nutdkuiuj_loc = (cs_real_6_t*)nutdkuiuj[kk];

          w1[iel][kk] = nutdkuiuj_loc[iel][iii]
                      + 2.*nut[iel]*(ui[iel][i]*duidxj[iel][j][kk]
                                   + ui[iel][j]*duidxj[iel][i][kk])
                      - nut[iel]*(uidujdxk_ii[iel][j][kk]+ uidujdxk_jj[iel][i][kk])
                      - ui[iel][i]*nutduidxj[iel][j][kk] - ui[iel][j]*nutduidxj[iel][i][kk]
                      - duidxj[iel][j][kk]*nutui[iel][i] - duidxj[iel][i][kk]*nutui[iel][j];
        }
      }

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        brij->budsgsfullij[iel][iii][2] = diverg[iel]/ro0;

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsfullij[iel][iii][3] = nutduidxkdujdxk[iel][iii] - nut[iel]*duidxkdujdxk[iel][iii];

        cs_real_t xx = 0.;
        for (int kk = 0 ; kk < 3 ; kk++)
          xx += 2.*nut[iel]*duidxj[iel][i][kk]*duidxj[iel][j][kk]
              - duidxj[iel][i][kk]*nutduidxj[iel][j][kk]
              - duidxj[iel][j][kk]*nutduidxj[iel][i][kk];

        brij->budsgsfullij[iel][iii][3] += xx;
        brij->budsgsfullij[iel][iii][3] = -2./ro0*brij->budsgsfullij[iel][iii][3];
      }
    }

    /* Bc coeffs */
    for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
      coefas[ifac] = 0.;
      if (bc_type[ifac] == CS_SMOOTHWALL
       || bc_type[ifac] == CS_ROUGHWALL)
        coefbs[ifac] = 0.;
      else
        coefbs[ifac] = 1.;
    }

    cs_gradient_scalar("nu_t",
                       gradient_type,
                       halo_type,
                       inc,
                       true,
                       var_cal_opt.nswrgr,
                       0,
                       0,
                       1,
                       var_cal_opt.iwarni,
                       var_cal_opt.imligr,
                       var_cal_opt.epsrgr,
                       var_cal_opt.extrag,
                       var_cal_opt.climgr,
                       NULL,
                       coefas,
                       coefbs,
                       nut,
                       NULL,
                       NULL,
                       w3);

    for (int iii = 0 ; iii < 6 ; iii++) {
      i = idirtens[iii][0];
      j = idirtens[iii][1];
      cs_real_33_t *uidujdxk_ii = (cs_real_33_t*)uidujdxk[i];
      cs_real_33_t *uidujdxk_jj = (cs_real_33_t*)uidujdxk[j];

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsfullij[iel][iii][4] = 0.;
        for (int kk = 0 ; kk < 3 ; kk++)
          brij->budsgsfullij[iel][iii][4] += w3[iel][kk]*
                                           ( uidujdxk_ii[iel][kk][j]-ui[iel][i]*duidxj[iel][kk][j]
                                           + uidujdxk_jj[iel][kk][i]-ui[iel][j]*duidxj[iel][kk][i] );
        brij->budsgsfullij[iel][iii][4] /= ro0;
      }
    }

    for (int iii = 0 ; iii < 6 ; iii++) {
      i = idirtens[iii][0];
      j = idirtens[iii][1];
      cs_real_33_t *uidujdxk_ii = (cs_real_33_t*)uidujdxk[i];
      cs_real_33_t *uidujdxk_jj = (cs_real_33_t*)uidujdxk[j];

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
       brij->budsgsfullij[iel][iii][5] = dnutdxkuidukdxjsym[iel][iii]-ui[iel][i]*dnutdxkdukdxi[iel][j]
                                                                     -ui[iel][j]*dnutdxkdukdxi[iel][i];
      }

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        cs_real_t xx = 0.;
        for (int kk = 0 ; kk < 3 ; kk++) {
          xx += 2.*dnutdxi[iel][kk]*(ui[iel][i]*duidxj[iel][kk][j]+ui[iel][j]*duidxj[iel][kk][i])
                 - dnutdxi[iel][kk]*(uidujdxk_ii[iel][kk][j] + uidujdxk_jj[iel][kk][i])
                 - duidxj[iel][kk][j]*uidnutdxj[iel][i][kk] - duidxj[iel][kk][i]*uidnutdxj[iel][j][kk];
        }

        brij->budsgsfullij[iel][iii][5] += xx;
        brij->budsgsfullij[iel][iii][5] /= ro0;
      }

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = (nutui[iel][j]-nut[iel]*ui[iel][j])*duidxj[iel][i][kk]
                      + (nutui[iel][i]-nut[iel]*ui[iel][i])*duidxj[iel][j][kk];

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        brij->budsgsfullij[iel][iii][6] = diverg[iel]/ro0;


      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsfullij[iel][iii][7] = 0.;
        for (int kk = 0  ; kk < 3 ; kk++)
          brij->budsgsfullij[iel][iii][7] -= (nutduidxj[iel][j][kk]-nut[iel]*duidxj[iel][j][kk])
                                             *duidxj[iel][i][kk]
                                           + (nutduidxj[iel][i][kk]-nut[iel]*duidxj[iel][i][kk])
                                             *duidxj[iel][j][kk];

        brij->budsgsfullij[iel][iii][7] /= ro0;
      }

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        brij->budsgsfullij[iel][iii][8] = 0.;
        for (int kk = 0  ; kk < 3 ; kk++)
          brij->budsgsfullij[iel][iii][8] += (uidnutdxj[iel][i][kk]-ui[iel][i]*dnutdxi[iel][kk])
                                             *duidxj[iel][kk][j]
                                           + (uidnutdxj[iel][j][kk]-ui[iel][j]*dnutdxi[iel][kk])
                                             *duidxj[iel][kk][i];

        brij->budsgsfullij[iel][iii][8] /= ro0;
      }

    }
  }

  BFT_FREE(uiujuk);
  BFT_FREE(nutdkuiuj);
  BFT_FREE(uidujdxk);

  BFT_FREE(w1);
  BFT_FREE(w2);
  BFT_FREE(w3);
  BFT_FREE(coefas);
  BFT_FREE(coefbs);
  BFT_FREE(diverg);
  BFT_FREE(lapl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Tui LES balance.
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_compute_tui(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const int kvisls0 = cs_field_key_id("diffusivity_ref");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  const int *bc_type = cs_glob_bc_type;

  cs_var_cal_opt_t var_cal_opt;
  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  cs_real_t *tdivturflux, *nut, *nutt, *nutdtdxidtdxi;
  cs_real_3_t *nutdtdxi, *uidivturflux, *tdtauijdxj;
  cs_real_3_t *nutui, *nutduidxjdtdxj, *dnutdxjtdujdxi, *nuttdtdxi;
  cs_real_3_t *dnutdxjdujdxi, *dnutdxi, *tdnutdxi, *tdtdxi;
  cs_real_33_t *nuttduidxj, *nutuidtdxj;

  cs_real_t    *p         = (cs_real_t    *)_les_balance_get_tm_by_name("p_m")->val;
  cs_real_3_t  *ui        = (cs_real_3_t  *)_les_balance_get_tm_by_name("ui_m")->val;
  cs_real_33_t *duidxj    = (cs_real_33_t *)_les_balance_get_tm_by_name("djui_m")->val;
  cs_real_6_t  *uiuj      = (cs_real_6_t  *)_les_balance_get_tm_by_name("uiuj_m")->val;
  cs_real_33_t *nutduidxj = (cs_real_33_t *)_les_balance_get_tm_by_name("nutdjui_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    nut           =                _les_balance_get_tm_by_name("nut_m")->val;
    nutui         = (cs_real_3_t *)_les_balance_get_tm_by_name("nutui_m")->val;
    dnutdxi       = (cs_real_3_t *)_les_balance_get_tm_by_name("dinut_m")->val;
    dnutdxjdujdxi = (cs_real_3_t *)_les_balance_get_tm_by_name("djnutdiuj_m")->val;
  }

  cs_real_t dtref  = cs_glob_time_step->dt_ref;
  cs_real_t ro0    = cs_glob_fluid_properties->ro0;
  cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  cs_real_t wvar, xx;
  int iii, jjj;

  /* Working local arrays */
  cs_real_t *diverg, *lapl, *w2;
  cs_real_t *coefas, *coefbs;
  cs_real_3_t *w1;

  BFT_MALLOC(w1    , n_cells_ext, cs_real_3_t);
  BFT_MALLOC(w2    , n_cells_ext, cs_real_t  );
  BFT_MALLOC(diverg, n_cells_ext, cs_real_t  );
  BFT_MALLOC(lapl  , n_cells_ext, cs_real_t  );
  BFT_MALLOC(coefbs, n_b_faces  , cs_real_t  );
  BFT_MALLOC(coefas, n_b_faces  , cs_real_t  );

  /* For each scalar */
  for (int isca = 0 ; isca < nscal ; isca++) {

    cs_les_balance_tui_t *b_sca = _les_balance.btui[isca];

    cs_field_t *sca = cs_field_by_id(b_sca->f_id);
    cs_real_t sigmas = cs_field_get_key_double(sca, ksigmas);
    cs_real_t visls0 = cs_field_get_key_double(sca, kvisls0);
    cs_real_t xvistot = visls0 + viscl0;

    /* Time moments retrieved by name */
    cs_real_t    *dtdxidtdxi = _les_balance_get_tm_by_scalar_id(isca, "ditdit_m")->val;
    cs_real_t    *t          = _les_balance_get_tm_by_scalar_id(isca, "t_m")->val;
    cs_real_t    *tp         = _les_balance_get_tm_by_scalar_id(isca, "tp_m")->val;
    cs_real_t    *t2         = _les_balance_get_tm_by_scalar_id(isca, "t_v")->val;

    cs_real_3_t  *tui         = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "tui_m")->val;
    cs_real_3_t  *ttui        = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "ttui_m")->val;
    cs_real_3_t  *dtdxi       = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "dtdxi_m")->val;
    cs_real_3_t  *pdtdxi      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "pdtdxi_m")->val;
    cs_real_3_t  *dtdxjduidxj = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "djtdjui_m")->val;

    cs_real_6_t  *tuiuj = (cs_real_6_t *)_les_balance_get_tm_by_scalar_id(isca, "tuiuj_m")->val;

    cs_real_33_t *tduidxj = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "tdjui_m")->val;
    cs_real_33_t *uidtdxj = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "uidjt_m")->val;

    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
      tdivturflux  =                _les_balance_get_tm_by_scalar_id(isca, "tdivturflux_m")->val;
      tdtauijdxj   = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "tdjtauij_m")->val;
      uidivturflux = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "uidivturflux_m")->val;
    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
      nutt           =                 _les_balance_get_tm_by_scalar_id(isca, "nutt_m")->val;
      nutdtdxidtdxi  =                 _les_balance_get_tm_by_scalar_id(isca, "nutditdit_m")->val;
      nutduidxjdtdxj = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nutdjuidjt_m")->val;
      dnutdxjtdujdxi = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "djnuttdiuj_m")->val;
      tdnutdxi       = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "tdinut_m")->val;
      tdtdxi         = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "tdit_m")->val;
      nuttdtdxi      = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nuttdit_m")->val;
      nuttduidxj     = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "nuttdjui_m")->val;
      nutuidtdxj     = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "nutuidjt_m")->val;
    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL ||
        _les_balance.type & CS_LES_BALANCE_TUI_BASE   )
      nutdtdxi       = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nutdit_m")->val;

    /* Initialization */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      for (int ii = 0 ; ii < 3 ; ii++) {
        b_sca->prodtUi[iel][ii]   = 0.;
        b_sca->prodtTi[iel][ii]   = 0.;
        b_sca->phiti[iel][ii]     = 0.;
        b_sca->epsti[iel][ii]     = 0.;
        b_sca->difftti[iel][ii]   = 0.;
        b_sca->diffttpi[iel][ii]  = 0.;
        b_sca->convti[iel][ii]    = 0.;
        b_sca->difflamti[iel][ii] = 0.;
      }
    }

    /* tptp */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      b_sca->tptp[iel] = t2[iel] - cs_math_sq(t[iel]);
      for (int ii = 0 ; ii < 3 ; ii++)
        b_sca->tpuip[iel][ii] = tui[iel][ii] - t[iel]*ui[iel][ii];
    }

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      for (int ii = 0 ; ii < 3 ; ii++) {
        b_sca->unstti[iel][ii] = (b_sca->tpuip[iel][ii]-b_sca->unstti[iel][ii])/dtref;
        b_sca->prodtUi[iel][ii] = 0.;
        b_sca->prodtTi[iel][ii] = 0.;
        b_sca->epsti[iel][ii] = dtdxjduidxj[iel][ii];

        for (int kk = 0 ; kk < 3 ; kk++) {
          iii = ipdirtens[ii][kk];
          b_sca->prodtUi[iel][ii] -= b_sca->tpuip[iel][kk]*duidxj[iel][ii][kk];
          wvar = uiuj[iel][iii] - ui[iel][ii]*ui[iel][kk];
          b_sca->prodtTi[iel][ii] -= wvar*dtdxi[iel][kk];
          b_sca->epsti[iel][ii] -= dtdxi[iel][kk]*duidxj[iel][ii][kk];
        }

        b_sca->epsti[iel][ii] *= -xvistot/ro0;
        b_sca->phiti[iel][ii] = pdtdxi[iel][ii]-p[iel]*dtdxi[iel][ii];
      }
    }

    /* convti */
    for (int ii = 0 ; ii < 3 ; ii++) {
      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = ui[iel][kk]*b_sca->tpuip[iel][ii];

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        b_sca->convti[iel][ii] = diverg[iel];

      /* difftti */
      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
        for (int kk = 0 ; kk < 3 ; kk++) {
          jjj = ipdirtens[ii][kk];
          w1[iel][kk] = - tuiuj[iel][jjj] - 2.*t[iel]*ui[iel][ii]*ui[iel][kk]
                                             + ui[iel][ii]*tui[iel][kk]
                                             + ui[iel][kk]*tui[iel][ii]
                                             + t[iel]*uiuj[iel][jjj];
        }
      }

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        b_sca->difftti[iel][ii] = diverg[iel];


      /* diffttpi */
      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = tp[iel] - t[iel]*p[iel];

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        b_sca->diffttpi[iel][iii] = diverg[iel]/ro0;


      /* Laminar diffusion */
      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = viscl0*(tduidxj[iel][ii][kk]-t[iel]*duidxj[iel][ii][kk])
                      + visls0*(uidtdxj[iel][ii][kk]-ui[iel][ii]*dtdxi[iel][kk]);

      _les_balance_divergence_vector(w1, diverg);

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        b_sca->difflamti[iel][ii] = diverg[iel]/ro0;
    }

    /* Variance budgets */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      b_sca->prodvar[iel]    = 0.;
      b_sca->epsvar[iel]     = 0.;
      b_sca->difftvar[iel]   = 0.;
      b_sca->convvar[iel]    = 0.;
      b_sca->difflamvar[iel] = 0.;
    }

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      b_sca->unstvar[iel] = (b_sca->tptp[iel] - b_sca->unstvar[iel])/dtref;
      b_sca->prodvar[iel] = 0.;
      b_sca->epsvar[iel] = dtdxidtdxi[iel];
      for (int kk = 0 ; kk < 3 ; kk++) {
        b_sca->prodvar[iel] += b_sca->tpuip[iel][kk]*dtdxi[iel][kk];
        b_sca->epsvar[iel] -= cs_math_sq(dtdxi[iel][kk]);
      }
      b_sca->prodvar[iel] *= -2.;
      b_sca->epsvar[iel] *= -2.*visls0/ro0;
    }

    /* convvar */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      for (int kk = 0 ; kk < 3 ; kk++)
        w1[iel][kk] = ui[iel][kk]*b_sca->tptp[iel];

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      b_sca->convvar[iel] = diverg[iel];

    /* difftti */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      for (int kk = 0 ; kk < 3 ; kk++) {
        w1[iel][kk] = ttui[iel][kk] + 2.*cs_math_sq(t[iel])*ui[iel][kk]
                                    - ui[iel][kk]*t2[iel]
                                    - 2.*t[iel]*tui[iel][kk];
      }
    }

    _les_balance_divergence_vector(w1, diverg);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      b_sca->difftvar[iel] = diverg[iel];

    /* Laminar diffusion */
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      w2[iel] = b_sca->tptp[iel];

    _les_balance_laplacian(w2, lapl, 1);

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
      b_sca->difflamvar[iel] = visls0*lapl[iel]/ro0;


    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {

      for (int ii = 0 ; ii < 3 ; ii++) {
        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = -nutduidxj[iel][ii][kk] - nutduidxj[iel][kk][ii];

        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstui[iel][ii] = -tdtauijdxj[iel][ii]-t[iel]*diverg[iel];
      }

      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        for (int kk = 0 ; kk < 3 ; kk++)
          w1[iel][kk] = -nutdtdxi[iel][kk];


      _les_balance_divergence_vector(w1, diverg);


      for (int ii = 0 ; ii < 3 ; ii++) {
        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstui[iel][ii] -= (uidivturflux[iel][ii]-ui[iel][ii]*diverg[iel]/sigmas);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstui[iel][ii] /= ro0;
      }

      /* Total SGS contribution for the variance of a scalar */
      for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
        b_sca->budsgsvar[iel] = -2.*(tdivturflux[iel]-t[iel]*diverg[iel]/sigmas)/ro0;

    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {

      for (int ii = 0 ; ii < 3 ; ii++) {
        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = nut[iel]*(tduidxj[iel][ii][kk]-t[iel]*duidxj[iel][ii][kk])
                        + nut[iel]*(uidtdxj[iel][ii][kk]-ui[iel][ii]*dtdxi[iel][kk])/sigmas;


        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstuifull[0][iel][ii] = diverg[iel]/ro0;


        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstuifull[1][iel][ii] = nut[iel]*(1.+1./sigmas)*b_sca->epsti[iel][ii]/xvistot;


        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          for (int kk = 0 ; kk < 3 ; kk++) {
            w1[iel][kk] = nuttduidxj[iel][ii][kk]
                        + 2.*nut[iel]*t[iel]*duidxj[iel][ii][kk]
                        - nut[iel]*tduidxj[iel][ii][kk]
                        - t[iel]*nutduidxj[iel][ii][kk]
                        - duidxj[iel][ii][kk]*nutt[iel];
            xx          = nutuidtdxj[iel][ii][kk]
                        + 2.*nut[iel]*uidtdxj[iel][ii][kk]
                        - ui[iel][ii]*nutdtdxi[iel][kk]
                        - dtdxi[iel][kk]*nutui[iel][ii];
            w1[iel][kk] += xx/sigmas;
          }
        }

        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstuifull[2][iel][ii] = diverg[iel]/ro0;

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          b_sca->budsgstuifull[3][iel][ii] = nutduidxjdtdxj[iel][ii] - nut[iel]*dtdxjduidxj[iel][ii];

          xx = 0.;
          for (int kk = 0 ; kk < 3 ; kk++)
            xx = xx + 2.*nut[iel]*duidxj[iel][ii][kk]*dtdxi[iel][kk]
                    - duidxj[iel][ii][kk]*nutdtdxi[iel][kk]
                    - dtdxi[iel][kk]*nutduidxj[iel][ii][kk];

          b_sca->budsgstuifull[3][iel][ii] += xx;
          b_sca->budsgstuifull[3][iel][ii] *= -(1.+1./sigmas)/ro0;
        }
      }

      /* Bc coeffs */
      for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
        coefas[ifac] = 0.;
        if (bc_type[ifac] == CS_SMOOTHWALL
         || bc_type[ifac] == CS_ROUGHWALL)
          coefbs[ifac] = 0.;
        else
          coefbs[ifac] = 1.;
      }

      cs_gradient_scalar("nu_t",
                         gradient_type,
                         halo_type,
                         1,
                         true,
                         var_cal_opt.nswrgr,
                         0,
                         0,
                         1,
                         var_cal_opt.iwarni,
                         var_cal_opt.imligr,
                         var_cal_opt.epsrgr,
                         var_cal_opt.extrag,
                         var_cal_opt.climgr,
                         NULL,
                         coefas,
                         coefbs,
                         nut,
                         NULL,
                         NULL,
                         w1);


      for (int ii = 0 ; ii < 3 ; ii++) {
        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          b_sca->budsgstuifull[4][iel][ii] = 0.;
          for (int kk = 0 ; kk < 3 ; kk++)
            b_sca->budsgstuifull[4][iel][ii] += w1[iel][kk]*(tduidxj[iel][ii][kk]-t[iel]*duidxj[iel][ii][kk]);
          b_sca->budsgstuifull[4][iel][ii] /= ro0;

          b_sca->budsgstuifull[5][iel][ii] = dnutdxjtdujdxi[iel][ii] - t[iel]*dnutdxjdujdxi[iel][ii];
          xx = 0.;
          for (int kk = 0 ; kk < 3 ; kk++)
            xx = xx + 2.*t[iel]*dnutdxi[iel][kk]*duidxj[iel][kk][ii]
                    - duidxj[iel][kk][ii]*tdnutdxi[iel][kk]
                    - dnutdxi[iel][kk]*tduidxj[iel][kk][ii];

          b_sca->budsgstuifull[5][iel][ii] += xx;
          b_sca->budsgstuifull[5][iel][ii] *= viscl0/ro0;
        }

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = duidxj[iel][ii][kk]*(nutt[iel]-nut[iel]*t[iel])
                        + dtdxi[iel][kk]*(nutui[iel][ii]-nut[iel]*ui[iel][ii])/sigmas;

        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgstuifull[6][iel][ii] = diverg[iel]/ro0;

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          b_sca->budsgstuifull[7][iel][ii] = 0.;
          b_sca->budsgstuifull[8][iel][ii] = 0.;
          for (int kk = 0 ; kk < 3 ; kk++) {
            b_sca->budsgstuifull[7][iel][ii] -= duidxj[iel][ii][kk]*(nutdtdxi[iel][kk]-nut[iel]*dtdxi[iel][kk])
                                              + dtdxi[iel][kk]*(nutduidxj[iel][ii][kk]-nut[iel]*duidxj[iel][ii][kk])/sigmas;
            b_sca->budsgstuifull[8][iel][ii] += duidxj[iel][kk][ii]*(tdnutdxi[iel][kk]-t[iel]*dnutdxi[iel][kk]);
          }
          b_sca->budsgstuifull[7][iel][ii] /= ro0;
          b_sca->budsgstuifull[8][iel][ii] /= ro0;
        }


        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = 2.*nut[iel]/sigmas*(tdtdxi[iel][kk]-t[iel]*dtdxi[iel][kk]);


        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgsvarfull[iel][0] = diverg[iel]/ro0;

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgsvarfull[iel][1] = b_sca->epsvar[iel]*nut[iel]/viscl0;

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = nuttdtdxi[iel][kk]
                        + 2.*nut[iel]*t[iel]*dtdxi[iel][kk]
                        - nut[iel]*tdtdxi[iel][kk]
                        - t[iel]*nutdtdxi[iel][kk]
                        - dtdxi[iel][kk]*nutt[iel];

        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgsvarfull[iel][2] = 2.*diverg[iel]/(ro0*sigmas);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          b_sca->budsgsvarfull[iel][3] = nutdtdxidtdxi[iel]-nut[iel]*dtdxidtdxi[iel];
          xx = 0.;
          for (int kk = 0 ; kk < 3 ; kk++)
            xx += 2.*nut[iel]*t[iel]*dtdxi[iel][kk]-t[iel]*nutdtdxi[iel][kk]-dtdxi[iel][kk]*nutt[iel];

          b_sca->budsgsvarfull[iel][3] += xx;
          b_sca->budsgsvarfull[iel][3] *= -2./(sigmas*ro0);
        }

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          for (int kk = 0 ; kk < 3 ; kk++)
            w1[iel][kk] = dtdxi[iel][kk]*(nutt[iel]-nut[iel]*t[iel]);

        _les_balance_divergence_vector(w1, diverg);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++)
          b_sca->budsgsvarfull[iel][4] = 2.*diverg[iel]/(sigmas*ro0);

        for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
          b_sca->budsgsvarfull[iel][5] = 0.;
          for (int kk = 0 ; kk < 3 ; kk++)
            b_sca->budsgsvarfull[iel][5] += dtdxi[iel][kk]*(nutdtdxi[iel][kk]-nut[iel]*dtdxi[iel][kk]);
          b_sca->budsgsvarfull[iel][5] *= -2./(sigmas*ro0);
        }
      }
    }
  }

  BFT_FREE(w1);
  BFT_FREE(w2);
  BFT_FREE(diverg);
  BFT_FREE(coefas);
  BFT_FREE(coefbs);
  BFT_FREE(lapl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write the LES balance structure in the LES balance restart file.
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_write_restart(void)
{
  char const *ficsui = "les_balance";
  cs_restart_t *suite = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_WRITE);

  if (suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the auxiliary restart "
                "file in write mode for the LES balance module.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              ficsui);

  { /* Write the header */
    char nomrub[] = "les_balance_type";

    int support = CS_MESH_LOCATION_NONE;

    cs_restart_write_section(suite,
                             nomrub,
                             support,
                             1,
                             CS_TYPE_cs_int_t,
                             &(_les_balance.type));

  }

  cs_restart_destroy(&suite);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy the LES balance structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_finalize(void)
{
  /* Freeing of the brij structure */
  _les_balance.brij = _les_balance_destroy_rij(_les_balance.brij);

  /* Freeing of the btui structure */
  _les_balance.btui = _les_balance_destroy_tui(_les_balance.btui);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Active the LES balance module.
 *
 * \param[in]  type_flag     mask of LES balance type
 * \param[in]  frequency_n   balance computing frequency in time-steps
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_activate(int     type_flag,
                        int     frequency_n)
{
  _les_balance.i_les_balance = 1;

  /* Filling _les_balance */
  _les_balance.type |= type_flag;
  _les_balance.frequency_n = frequency_n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the LES balance for Tui or Rij
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_compute(void)
{
  int nt_cur = cs_glob_time_step->nt_cur;

  if ((_les_balance.frequency_n > 0
   && nt_cur % _les_balance.frequency_n == 0) ||
      nt_cur == cs_glob_time_step->nt_max) {
    if (_les_balance.type & CS_LES_BALANCE_RIJ)
      cs_les_balance_compute_rij();

    if (_les_balance.type & CS_LES_BALANCE_TUI)
      cs_les_balance_compute_tui();
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
