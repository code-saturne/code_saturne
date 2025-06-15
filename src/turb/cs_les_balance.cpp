/*============================================================================
 * LES Balance
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

/*============================================================================
 * Functions dealing with LES balance
 *============================================================================*/

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_file.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "alge/cs_divergence.h"
#include "alge/cs_face_viscosity.h"
#include "alge/cs_convection_diffusion.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_mem.h"
#include "mesh/cs_geom.h"
#include "alge/cs_gradient.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_physical_properties.h"
#include "base/cs_prototypes.h"
#include "base/cs_restart.h"
#include "base/cs_time_moment.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_les_balance.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_les_balance.cpp
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
  \var  cs_les_balance_rij_t::difflamij
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

#undef _IJV2T
#undef _IJV2T3
#undef _PIJV2T
#undef _PIJV2T3
#define _IJV2T { \
    {0, 0}, \
    {1, 1}, \
    {2, 2}, \
    {0, 1}, \
    {1, 2}, \
    {0, 2}  \
};

#define _PIJV2T { \
    {5, 5, 5}, \
    {5, 1, 4}, \
    {5, 4, 2}  \
};

#define _IJV2T3 { \
    {0, 1, 2}, \
    {0, 0, 0}, \
    {0, 0, 1}, \
    {0, 0, 2}, \
    {1, 1, 0}, \
    {1, 1, 1}, \
    {1, 1, 2}, \
    {2, 2, 0}, \
    {2, 2, 1}, \
    {2, 2, 2}  \
};

#define _PIJV2T3 { \
    { \
        {7, 7, 7}, \
        {7, 7, 7}, \
        {7, 7, 7}  \
    }, \
    { \
        {7, 7, 7}, \
        {7, 5, 6}, \
        {7, 6, 8}  \
    }, \
    { \
        {7, 7, 7}, \
        {7, 6, 8}, \
        {7, 8, 9}  \
    }  \
};

/* Total number of scalars */
static int nscal = 0;

/* Global static structure _les_balance */
static cs_les_balance_t  _les_balance = {.brij = nullptr,
                                         .btui = nullptr,
                                         .i_les_balance = 0,
                                         .type = 0,
                                         .frequency_n = -1};

static cs_field_t *_gradv   = nullptr;
static cs_field_t *_gradnut = nullptr;
static cs_field_t **_gradt  = nullptr;

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
 * returns  pointer to a field or nullptr.
 *----------------------------------------------------------------------------*/

static cs_field_t *
_les_balance_get_tm_by_name(const char *name)
{
  cs_field_t *f = cs_field_by_name_try(name);
  if (! (f->type & CS_FIELD_ACCUMULATOR))
    f = nullptr;

  return f;
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
  char csca[12];

  memset(buffer, '\0', 32);
  strcpy(buffer, name);

  /* Add the scalar number to this name */
  snprintf(csca, 12, "_%d", isca);
  csca[11] = '\0';
  strcat(buffer, csca);
}

/*----------------------------------------------------------------------------*
 * Get a time moment field by scalar id.
 *
 * parameters:
 *   scalar_id <-- scalar id
 *   name      <-- time moment name
 *
 * returns  pointer to a field or nullptr.
 *----------------------------------------------------------------------------*/

static cs_field_t *
_les_balance_get_tm_by_scalar_id(int scalar_id,
                                 const char *name)
{
  cs_field_t *f = nullptr;
  char *buffer;

  CS_MALLOC(buffer, 32, char);

  _les_balance_get_tm_label(scalar_id, name, buffer);
  f = _les_balance_get_tm_by_name((const char *)buffer);

  CS_FREE(buffer);

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
                       int          type)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_real_t *b_dist = fvq->b_dist;
  const cs_real_t *cell_vol = fvq->cell_vol;

  const int *bc_type = cs_glob_bc_type;

  cs_dispatch_context ctx;

  cs_field_bc_coeffs_t bc_coeffs_loc; // a voir ci dessous (a=af,b=bf)
  cs_field_bc_coeffs_init(&bc_coeffs_loc);

  CS_MALLOC_HD(bc_coeffs_loc.a,  n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_loc.b,  n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_loc.af, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_loc.bf, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *coefa = bc_coeffs_loc.a;
  cs_real_t *coefb = bc_coeffs_loc.b;
  cs_real_t *coefaf = bc_coeffs_loc.af;
  cs_real_t *coefbf = bc_coeffs_loc.bf;

  const cs_real_t visc = 1., pimp = 0., qimp = 0.; // hext = -1;
  //cs_real_t a, b;

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    cs_real_t hint = visc / b_dist[face_id];

    if (   type == 0
        && (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL) ) {
      /*cs_boundary_conditions_set_dirichlet_scalar(&a,
                                                  &coefaf[face_id],
                                                  &b,
                                                  &coefbf[face_id],
                                                  pimp,
                                                  hint,
                                                  hext);
      */

      coefaf[face_id] = -hint*pimp;
      coefbf[face_id] =  hint;

      coefa[face_id] = coefaf[face_id];
      coefb[face_id] = coefbf[face_id];

    }
    else {
      /*
      cs_boundary_conditions_set_neumann_scalar(&a,
                                                &coefaf[face_id],
                                                &b,
                                                &coefbf[face_id],
                                                qimp,
                                                hint);
      */

      coefaf[face_id] = qimp;
      coefbf[face_id] = 0.;

      coefa[face_id] = coefaf[face_id];
      coefb[face_id] = coefbf[face_id];

    }
  });

  ctx.wait();

  cs_real_t *c_visc, *i_visc, *b_visc;
  CS_MALLOC_HD(c_visc, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_arrays_set_value<cs_real_t, 1>(n_cells, visc, c_visc);

  cs_face_viscosity(m,
                    fvq,
                    0,      /* mean type */
                    c_visc,
                    i_visc,
                    b_visc);
  CS_FREE(c_visc);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(CS_F_(vel));
  cs_equation_param_t _eqp = *eqp;
  _eqp.iconv = 0; /* only diffusion */
  _eqp.theta = 1.;

  cs_convection_diffusion_scalar(0,              /* idtvar */
                                 -1,             /* f_id */
                                 _eqp,
                                 0,              /* icvflb (not used) */
                                 1,              /* inc */
                                 1,              /* imasac (not used) */
                                 wa,             /* pvar */
                                 nullptr,           /* pvara (not used) */
                                 0,              /* icvfli (not used) */
                                 &bc_coeffs_loc, /* coefa & b not used */
                                 i_visc,         /* mass flux (not used) */
                                 b_visc,         /* mass flux (not used) */
                                 i_visc,
                                 b_visc,
                                 res,
                                 nullptr, nullptr);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    res[c_id] /= cell_vol[c_id];
  });

  ctx.wait();

  CS_FREE(coefa);
  CS_FREE(coefb);
  CS_FREE(coefaf);
  CS_FREE(coefbf);

  CS_FREE(i_visc);
  CS_FREE(b_visc);
}

/*----------------------------------------------------------------------------
 * Compute the divergence of a cell-based vector.
 *
 * Computation options are the ones of the velocity
 *
 * parameters:
 *   wa   <--  vector array
 *   res  <->  divergence of wa
 *----------------------------------------------------------------------------*/

 static void
 _les_balance_divergence_vector(cs_real_3_t  *wa,
                                cs_real_t    *res)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const int *bc_type = cs_glob_bc_type;

  cs_dispatch_context ctx;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(CS_F_(vel));

  cs_field_bc_coeffs_t bc_coeffs_v_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_v_loc);
  CS_MALLOC_HD(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_3_t  *coefav = (cs_real_3_t  *)bc_coeffs_v_loc.a;
  cs_real_33_t *coefbv = (cs_real_33_t *)bc_coeffs_v_loc.b;

  cs_real_t *i_massflux, *b_massflux;
  CS_MALLOC_HD(i_massflux, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_massflux, n_b_faces, cs_real_t, cs_alloc_mode);

  /* Bc coeffs */
  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
        coefav[face_id][ii] = 0.;

        for (cs_lnum_t jj = 0; jj < 3; jj++)
          coefbv[face_id][ii][jj] = 0.;

        if (!(   bc_type[face_id] == CS_SMOOTHWALL
              || bc_type[face_id] == CS_ROUGHWALL))
            coefbv[face_id][ii][ii] = 1.;
    }
  });

  ctx.wait();

  int f_id = -1;
  int itypfl = 0;
  int iflmb0 = 1;
  int init = 1;
  int inc = 1;

  cs_mass_flux(m,
               mq,
               f_id,
               itypfl,
               iflmb0,
               init,
               inc,
               eqp->imrgra,
               eqp->nswrgr,
               static_cast<cs_gradient_limit_t>(eqp->imligr),
               eqp->verbosity,
               eqp->epsrgr,
               eqp->climgr,
               nullptr,
               nullptr,
               (const cs_real_3_t *)wa,
               &bc_coeffs_v_loc,
               i_massflux,
               b_massflux);

  cs_divergence(m,
                init,
                i_massflux,
                b_massflux,
                res);

  CS_FREE(coefav);
  CS_FREE(coefbv);

  CS_FREE(i_massflux);
  CS_FREE(b_massflux);
}

/*----------------------------------------------------------------------------
 *  Compute the most needed gradients at each iteration.
 *----------------------------------------------------------------------------*/

static void
_les_balance_compute_gradients(void)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const int *bc_type = cs_glob_bc_type;
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(CS_F_(vel));

  cs_dispatch_context ctx;

  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;

  /* Computation of the velocity gradient */
  cs_gradient_type_by_imrgra(eqp->imrgra,
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
  if (   _les_balance.type & CS_LES_BALANCE_RIJ_FULL
      || _les_balance.type & CS_LES_BALANCE_TUI_FULL) {

    cs_field_bc_coeffs_t bc_coeffs_loc;
    CS_MALLOC_HD(bc_coeffs_loc.a, n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_loc.b, n_b_faces, cs_real_t, cs_alloc_mode);

    cs_real_t *coefas = bc_coeffs_loc.a;
    cs_real_t *coefbs = bc_coeffs_loc.b;

    /* Bc coeffs */
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      coefas[face_id] = 0.;

      if (   bc_type[face_id] == CS_SMOOTHWALL
          || bc_type[face_id] == CS_ROUGHWALL)
        coefbs[face_id] = 0.;
      else
        coefbs[face_id] = 1.;
    });

    ctx.wait();

    cs_gradient_scalar("nu_t",
                       gradient_type,
                       halo_type,
                       inc,
                       eqp->nswrgr,
                       0,
                       1,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       nullptr,
                       &bc_coeffs_loc,
                       CS_F_(mu_t)->val,
                       nullptr,
                       nullptr,
                       (cs_real_3_t *)_gradnut->val);

    CS_FREE(coefas);
    CS_FREE(coefbs);
  }

  /* Computation of scalar gradients */
  if (_les_balance.type & CS_LES_BALANCE_TUI) {

    const int keysca = cs_field_key_id("scalar_id");
    int iii = 0;

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int isca = cs_field_get_key_int(f, keysca);
      if (isca > 0) {
        const cs_equation_param_t *eqps
          = cs_field_get_equation_param_const(f);

        cs_gradient_type_by_imrgra(eqps->imrgra,
                                   &gradient_type,
                                   &halo_type);

        cs_field_gradient_scalar(f,
                                 false, /* use_previous_t */
                                 inc,
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

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_t *pressure = CS_F_(p)->val;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t pre = pressure[c_id];

    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      cs_lnum_t i = idirtens[ii][0];
      cs_lnum_t j = idirtens[ii][1];

      const cs_lnum_t id = 6*c_id + ii;
      vals[id] = pre*(grdv[c_id][i][j] + grdv[c_id][j][i]);
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_t *cpro_smago = cs_field_by_name("smagorinsky_constant^2")->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    vals[c_id] = cs_math_sq(cpro_smago[c_id]);
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      cs_lnum_t i = idirtens[ii][0];
      cs_lnum_t j = idirtens[ii][1];
      const cs_lnum_t id = 6*c_id + ii;
      vals[id] = 0.;

      for (cs_lnum_t k = 0; k < 3; k++)
        vals[id] += grdv[c_id][i][k] + grdv[c_id][j][k];
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_3_t *velocity = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  cs_real_t *diverg;
  cs_real_3_t *vel;
  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(vel, n_cells, cs_real_3_t, cs_alloc_mode);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  for (cs_lnum_t i = 0; i < 3; i++) {

    for (cs_lnum_t j = 0; j < 3; j++) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t k = 0; k < 3; k++)
          vel[c_id][k] = -mu_t[c_id]*(  grdv[c_id][j][k]
                                      + grdv[c_id][k][j]);
      });

      ctx.wait();

      _les_balance_divergence_vector(vel, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        const cs_lnum_t id = 9*c_id+i*3+j;
        vals[id] = velocity[c_id][i]*diverg[c_id];
      });

    }
  }

  ctx.wait();

  CS_FREE(diverg);
  CS_FREE(vel);
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

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      int i = idirtens[ii][0];
      int j = idirtens[ii][1];
      const cs_lnum_t id = 6*c_id + ii;

      vals[id] = 0.;
      for (cs_lnum_t k = 0; k < 3; k++)
        vals[id] += mu_t[c_id]*grdv[c_id][i][k]*grdv[c_id][j][k];
    }
  });

  ctx.wait();
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

  CS_UNUSED(input);

  cs_dispatch_context ctx;

  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      cs_lnum_t i = idirtens[ii][0];
      cs_lnum_t j = idirtens[ii][1];

      const cs_lnum_t id = 6*c_id + ii;
      vals[id] = 0.;
      for (cs_lnum_t k = 0; k < 3; k++)
        vals[id] += grdnu[c_id][i]
                      *(  vel[c_id][i]*grdv[c_id][k][j]
                        + vel[c_id][j]*grdv[c_id][k][i]);
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_3_t *vel   = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t   *mu_t  = CS_F_(mu_t)->val;
  cs_real_6_t *tens  = (cs_real_6_t *)vals;
  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      cs_lnum_t i = idirtens[ii][0];
      cs_lnum_t j = idirtens[ii][1];
      tens[c_id][ii] = mu_t[c_id]
                      *( vel[c_id][i]*grdv[c_id][j][*k]
                        +vel[c_id][j]*grdv[c_id][i][*k]);
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      const cs_lnum_t id = 3*c_id + i;
      vals[id] = 0.;
      for (cs_lnum_t j = 0; j < 3; j++)
        vals[id] += grdnu[c_id][i]*grdv[c_id][i][j];
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;
  cs_real_3_t *vel   = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_33_t *tens = (cs_real_33_t *)vals;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        tens[c_id][i][j] = vel[c_id][i]*grdnu[c_id][j];
    }
  });

  ctx.wait();
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
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      const cs_lnum_t id = 3*c_id+i;
      cs_real_t dtdxjduidxj = 0.;

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        dtdxjduidxj += grdt[c_id][kk]*grdv[c_id][i][kk];

      vals[id] = dtdxjduidxj;
    }
  });

  ctx.wait();
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
  cs_real_t *sca_val = sca->val;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

  cs_dispatch_context ctx;

  const int idirtens[6][2] = _IJV2T;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      cs_lnum_t i = idirtens[ii][0];
      cs_lnum_t j = idirtens[ii][1];

      const cs_lnum_t id = 6*c_id+ii;
      vals[id] = sca_val[c_id]*vel[c_id][i]*vel[c_id][j];
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
        tens[c_id][i][j] = vel[c_id][i]*grdt[c_id][j];
  });

  ctx.wait();
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
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_dispatch_context ctx;

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t dtdxidtdxi = 0.;
    for (cs_lnum_t i = 0; i < 3; i++)
      dtdxidtdxi += grdt[c_id][i]; // FIXME: missing SQUARE??

    vals[c_id] = dtdxidtdxi;
  });

  ctx.wait();
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
  cs_real_t *sca_val = sca->val;
  const int keysca = cs_field_key_id("scalar_id");
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_dispatch_context ctx;

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
    }
  }

  cs_real_t *diverg;
  cs_real_3_t *w1;
  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w1, n_cells, cs_real_3_t, cs_alloc_mode);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  for (cs_lnum_t i = 0; i < 3; i++) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t k = 0; k < 3; k++)
        w1[c_id][k] = -mu_t[c_id]*(  grdv[c_id][i][k]
                                   + grdv[c_id][k][i]);
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_lnum_t id = 3*c_id+i;
      vals[id] = sca_val[c_id]*diverg[c_id];
    });
  }

  ctx.wait();

  CS_FREE(diverg);
  CS_FREE(w1);
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

  cs_dispatch_context ctx;

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      if (f_id == sca->id)
        break;
    }
  }

  cs_real_t sigmas, *diverg;
  cs_real_3_t *w1, *vel;

  vel = (cs_real_3_t *)CS_F_(vel)->val;
  sigmas = cs_field_get_key_double(sca, ksigmas);

  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w1, n_cells, cs_real_3_t, cs_alloc_mode);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  for (cs_lnum_t i = 0; i < 3; i++) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t k = 0; k < 3; k++)
        w1[c_id][k] =  cs_math_sq(mu_t[c_id])/sigmas
                     *(grdv[c_id][i][k]+grdv[c_id][k][i]);
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_lnum_t id = 3*c_id+i;
      vals[id] = vel[c_id][i]*diverg[c_id];
    });
  }

  ctx.wait();

  CS_FREE(diverg);
  CS_FREE(w1);
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
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_dispatch_context ctx;

  cs_real_t sigmas, *diverg;
  cs_real_3_t *w1;

  sigmas = cs_field_get_key_double(sca, ksigmas);

  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w1, n_cells, cs_real_3_t, cs_alloc_mode);

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t k = 0; k < 3; k++)
      w1[c_id][k] = mu_t[c_id]/sigmas * grdt[c_id][k];
  });

  ctx.wait();

  _les_balance_divergence_vector(w1, diverg);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    vals[c_id] = sca->val[c_id]*diverg[c_id];
  });

  ctx.wait();

  CS_FREE(diverg);
  CS_FREE(w1);
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
  cs_real_t *sca_val = sca->val;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_dispatch_context ctx;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    cs_real_t nutditdit = 0.;
    for (cs_lnum_t i = 0; i < 3; i++)
      nutditdit += mu_t[c_id]*sca_val[c_id]*cs_math_sq(grdt[c_id][i]);

    vals[c_id] = nutditdit;
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++)
      for (cs_lnum_t j = 0; j < 3; j++)
        tens[c_id][i][j] = mu_t[c_id]*vel[c_id][i]*grdt[c_id][j];
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  const int keysca = cs_field_key_id("scalar_id");
  int isca = cs_field_get_key_int(sca, keysca) - 1;
  assert(isca > -1);

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdt = (cs_real_3_t *)_gradt[isca]->val;
  cs_real_t *mu_t = CS_F_(mu_t)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      cs_real_t nutdjuidjt = 0.;

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        nutdjuidjt += mu_t[c_id]*grdv[c_id][i][kk]*grdt[c_id][kk];

      const cs_lnum_t id = 3*c_id+i;
      vals[id] = nutdjuidjt;
    }
  });

  ctx.wait();
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

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      cs_real_t djnutdiuj = 0.;
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        djnutdiuj += grdnu[c_id][kk]*grdv[c_id][kk][i];

      const cs_lnum_t id = 3*c_id+i;
      vals[id] = djnutdiuj;
    }
  });

  ctx.wait();
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
  cs_real_t *sca_val = sca->val;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_dispatch_context ctx;

  cs_real_33_t *grdv = (cs_real_33_t *)_gradv->val;
  cs_real_3_t *grdnu = (cs_real_3_t *)_gradnut->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      cs_real_t djnuttdiuj = 0.;

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        djnuttdiuj += grdnu[c_id][kk]*sca_val[c_id]*grdv[c_id][kk][i];

      const cs_lnum_t id = 3*c_id+i;
      vals[id] = djnuttdiuj;
    }
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Declare generic time moments for either the Rij or the Tui LES balance.
 * Time moments are defined either by field ids or by function.
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
                                       nullptr);
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
                                       nullptr);
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
                                       nullptr);
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
                                       nullptr);

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
                                       nullptr);
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
                                       nullptr);
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
                                       nullptr);
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
                                         nullptr);
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
                                         nullptr);
    }
  }
}

/*----------------------------------------------------------------------------
 * Declare time moments that can be define either by field ids or by function.
 * for the Rij LES balance
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
                                       nullptr);
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
                                       nullptr);
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
                                         nullptr);
    }
  }

  /* p(djui+diuj) mean */
  cs_time_moment_define_by_func("pdjuisym_m",
                                CS_MESH_LOCATION_CELLS,
                                6,
                                true, /* intensive*/
                                _les_balance_compute_pdjuisym,
                                nullptr,
                                nullptr,
                                nullptr,
                                CS_TIME_MOMENT_MEAN,
                                1,
                                -1,
                                CS_TIME_MOMENT_RESTART_AUTO,
                                nullptr);

  /* dk(ui+uj) mean */
  cs_time_moment_define_by_func("dkuidkuj_m",
                                CS_MESH_LOCATION_CELLS,
                                6,
                                true, /* intensive*/
                                _les_balance_compute_dkuidkuj,
                                nullptr,
                                nullptr,
                                nullptr,
                                CS_TIME_MOMENT_MEAN,
                                1,
                                -1,
                                CS_TIME_MOMENT_RESTART_AUTO,
                                nullptr);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {

    if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {

      /* smagorinsky variance */
      cs_time_moment_define_by_func("smag_v",
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    true, /* intensive*/
                                    _les_balance_compute_smag,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    CS_TIME_MOMENT_VARIANCE,
                                    1,
                                      -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    nullptr);
    }

    /* uidktaujk */
    cs_time_moment_define_by_func("uidktaujk_m",
                                  CS_MESH_LOCATION_CELLS,
                                  9,
                                  true, /* intensive*/
                                  _les_balance_compute_uidktaujk,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  nullptr);
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
                                           nullptr);
      }
    }

    /* nutdkuidkuj mean */
    cs_time_moment_define_by_func("nutdkuidkuj_m",
                                  CS_MESH_LOCATION_CELLS,
                                  6,
                                  true, /* intensive*/
                                  _les_balance_compute_nutdkuidkuj,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  nullptr);

    /* dknutuidjuksym mean */
    cs_time_moment_define_by_func("dknutuidjuksym_m",
                                  CS_MESH_LOCATION_CELLS,
                                  6,
                                  true, /* intensive*/
                                  _les_balance_compute_dknutuidjuksym,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  nullptr);

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
                                      true, /* intensive*/
                                      _les_balance_compute_nutdkuiuj,
                                      &(index[k]),
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);
      }
    }

    /* dknutdiuk mean */
    cs_time_moment_define_by_func("dknutdiuk_m",
                                  CS_MESH_LOCATION_CELLS,
                                  3,
                                  true, /* intensive*/
                                  _les_balance_compute_dknutdiuk,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  nullptr);

    {
      cs_time_moment_define_by_func("uidjnut_m",
                                    CS_MESH_LOCATION_CELLS,
                                    9,
                                    true, /* intensive*/
                                    _les_balance_compute_uidjnut,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    nullptr);
    }
  }
}

/*----------------------------------------------------------------------------
 * Declare time moments that can be define either by field ids or by function.
 * for the Tui LES balance
 *----------------------------------------------------------------------------*/

static void
_les_balance_time_moment_tui(void)
{
  const int keysca = cs_field_key_id("scalar_id");
  char buffer[32];
  int isca = 0;

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    /* _djnutdiuj */
    cs_time_moment_define_by_func("djnutdiuj_m",
                                  CS_MESH_LOCATION_CELLS,
                                  3,
                                  true, /* intensive*/
                                  _les_balance_compute_djnutdiuj,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  CS_TIME_MOMENT_MEAN,
                                  1,
                                  -1,
                                  CS_TIME_MOMENT_RESTART_AUTO,
                                  nullptr);
  }

  /* Define time moments for T.ui balance */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
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
                                           nullptr);
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
                                           nullptr);
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
                                           nullptr);
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
                                           nullptr);
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
                                           nullptr);
      }

      {
        /* T.ui variance */
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
                                           nullptr);
      }

      /* djtdjui mean */
      _les_balance_get_tm_label(isca, "djtdjui_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    3,
                                    true, /* intensive*/
                                    _les_balance_compute_djtdjui,
                                    f,
                                    nullptr,
                                    nullptr,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    nullptr);
      /* tuiuj mean */
      _les_balance_get_tm_label(isca, "tuiuj_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    6,
                                    true, /* intensive*/
                                    _les_balance_compute_tuiuj,
                                    f,
                                    nullptr,
                                    nullptr,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    nullptr);

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
                                           nullptr);
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
                                           nullptr);
      }

      {
        _les_balance_get_tm_label(isca, "uidjt_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      9,
                                      true, /* intensive*/
                                      _les_balance_compute_uidjt,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);
      }

      /* ditdit mean */
      _les_balance_get_tm_label(isca, "ditdit_m", buffer);
      cs_time_moment_define_by_func(buffer,
                                    CS_MESH_LOCATION_CELLS,
                                    1,
                                    true, /* intensive*/
                                    _les_balance_compute_ditdit,
                                    f,
                                    nullptr,
                                    nullptr,
                                    CS_TIME_MOMENT_MEAN,
                                    1,
                                    -1,
                                    CS_TIME_MOMENT_RESTART_AUTO,
                                    nullptr);
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
                                           nullptr);
      }

      if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
        /* mean of T.d(tau_ij)/dxj :  _tdjtauij mean */
        _les_balance_get_tm_label(isca, "tdjtauij_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      true, /* intensive*/
                                      _les_balance_compute_tdjtauij,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);

        /* mean of U_i. d(U'_j. T')/dxi : _uidivturflux mean */
        _les_balance_get_tm_label(isca, "uidivturflux_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      true, /* intensive*/
                                      _les_balance_compute_uidivturflux,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);

        /* mean of T.div(U'T') : _tdivturflux mean */
        _les_balance_get_tm_label(isca, "tdivturflux_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      true, /* intensive*/
                                      _les_balance_compute_tdivturflux,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);
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
                                             nullptr);
        }

        /* nutditdit mean */
        _les_balance_get_tm_label(isca, "nutditdit_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      true, /* intensive*/
                                      _les_balance_compute_nutditdit,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);

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
                                             nullptr);
        }

        {
          /* nutuidjt mean */
          _les_balance_get_tm_label(isca, "nutuidjt_m", buffer);
          cs_time_moment_define_by_func(buffer,
                                        CS_MESH_LOCATION_CELLS,
                                        9,
                                        true, /* intensive*/
                                        _les_balance_compute_nutuidjt,
                                        f,
                                        nullptr,
                                        nullptr,
                                        CS_TIME_MOMENT_MEAN,
                                        1,
                                        -1,
                                        CS_TIME_MOMENT_RESTART_AUTO,
                                        nullptr);
        }

        /* nutdjuidjt mean */
        _les_balance_get_tm_label(isca, "nutdjuidjt_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      true, /* intensive*/
                                      _les_balance_compute_nutdjuidjt,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);

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
                                             nullptr);
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
                                             nullptr);
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
                                             nullptr);
        }

        /* djnuttdiuj mean */
        _les_balance_get_tm_label(isca, "djnuttdiuj_m", buffer);
        cs_time_moment_define_by_func(buffer,
                                      CS_MESH_LOCATION_CELLS,
                                      3,
                                      true, /* intensive*/
                                      _les_balance_compute_djnuttdiuj,
                                      f,
                                      nullptr,
                                      nullptr,
                                      CS_TIME_MOMENT_MEAN,
                                      1,
                                      -1,
                                      CS_TIME_MOMENT_RESTART_AUTO,
                                      nullptr);
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
                                             nullptr);
        }

      }
      isca++;
    }
  }
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
  cs_les_balance_rij_t *brij = nullptr;

  CS_MALLOC(brij, 1, cs_les_balance_rij_t);

  /* Allocation of the working arrays  that cannot
     be declared using time moments */

  CS_MALLOC_HD(brij->prodij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->epsij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->phiij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->difftij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->difftpij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->convij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->difflamij, n_cells, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(brij->unstij, n_cells, cs_real_6_t, cs_alloc_mode);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE)
    CS_MALLOC_HD(brij->budsgsij, n_cells, cs_real_6_t, cs_alloc_mode);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL)
    CS_MALLOC_HD(brij->budsgsfullij, n_cells, cs_real_69_t, cs_alloc_mode);

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
  cs_les_balance_tui_t *btui = nullptr;

  CS_MALLOC(btui, 1, cs_les_balance_tui_t);

  /* Allocation of the working arrays  that cannot
     be declared using time moments */

  /* Working arrays */
  CS_MALLOC_HD(btui->unstvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->tptp, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->prodvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->epsvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->difftvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->convvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->difflamvar, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->tpuip, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->unstti, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->prodtUi, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->prodtTi, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->phiti, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->epsti, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->difftti, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->diffttpi, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->convti, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(btui->difflamti, n_cells, cs_real_3_t, cs_alloc_mode);

  if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
    CS_MALLOC_HD(btui->budsgsvar, n_cells, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(btui->budsgstui, n_cells, cs_real_3_t, cs_alloc_mode);
  }

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    CS_MALLOC_HD(btui->budsgsvarfull, n_cells, cs_real_6_t, cs_alloc_mode);
    CS_MALLOC_HD(btui->budsgstuifull, 10, cs_real_3_t *, cs_alloc_mode);
    for (int ii = 0; ii < 10; ii++)
      CS_MALLOC_HD(btui->budsgstuifull[ii], n_cells, cs_real_3_t, cs_alloc_mode);
  }

  return btui;
}

/*----------------------------------------------------------------------------
 * Initialize the brij structure of _les_balance
 *----------------------------------------------------------------------------*/

static void
_les_balance_initialize_rij(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_les_balance_rij_t *brij = _les_balance.brij;

  cs_dispatch_context ctx;

  cs_real_6_t *prodij = brij->prodij;
  cs_real_6_t *epsij = brij->epsij;
  cs_real_6_t *phiij = brij->phiij;
  cs_real_6_t *difftij = brij->difftij;
  cs_real_6_t *difftpij = brij->difftpij;
  cs_real_6_t *convij = brij->convij;
  cs_real_6_t *difflamij = brij->difflamij;
  cs_real_6_t *unstij = brij->unstij;
  cs_real_6_t *budsgsij = brij->budsgsij;

  cs_real_69_t *budsgsfullij = brij->budsgsfullij;

  const bool is_rij_base = (_les_balance.type & CS_LES_BALANCE_RIJ_BASE);
  const bool is_rij_full = (_les_balance.type & CS_LES_BALANCE_RIJ_FULL);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    for (cs_lnum_t i = 0; i < 6; i++) {
      prodij[c_id][i] = 0.;
      epsij[c_id][i] = 0.;
      phiij[c_id][i] = 0.;
      difftij[c_id][i] = 0.;
      difftpij[c_id][i] = 0.;
      convij[c_id][i] = 0.;
      difflamij[c_id][i] = 0.;
      unstij[c_id][i] = 0.;

      if (is_rij_base)
        budsgsij[c_id][i] = 0.;
    }

    if (is_rij_full)
      for (cs_lnum_t i = 0; i < 9; i++)
        for (cs_lnum_t j = 0; j < 6; j++)
          budsgsfullij[c_id][i][j] = 0.;

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Check that in case of a restart calculation, the current LES balance type
 * is the same as in the previous calculation.
 *----------------------------------------------------------------------------*/

static void
_check_restart_type(void)
{
  int b_type = -1;

  cs_restart_t *rp = nullptr;
  const char default_path[] = "restart/les_balance.csc";

  int is_reg = cs_file_isreg(default_path);
  if (is_reg)
    rp = cs_restart_create("les_balance.csc",
                           nullptr,
                           CS_RESTART_MODE_READ);

  else {
    is_reg = cs_file_isreg("restart/les_balance");
    if (is_reg)
      rp = cs_restart_create("les_balance",
                             nullptr,
                             CS_RESTART_MODE_READ);
  }

  if (rp == nullptr)
    bft_error
      (__FILE__, __LINE__, 0,
       _("LES balance restart file not present: %s."),
       default_path);

  if (rp != nullptr)  {
    cs_restart_read_section(rp,
                            "les_balance_type",
                            CS_MESH_LOCATION_NONE,
                            1,
                            CS_TYPE_int,
                            &b_type);

    if (!(b_type & _les_balance.type))
      bft_error(__FILE__, __LINE__, 0,
                _("Abort while reading the LES balance restart file: %s\n"
                  "The previous balance type is different from the current\n"
                  "balance type:\n"
                  "  previous type: %d\n"
                  "  current type:  %d\n"),
                cs_restart_get_name(rp), b_type, _les_balance.type);
  }

  cs_restart_destroy(&rp);
}

/*----------------------------------------------------------------------------
 * Initialize a btui structure of _les_balance
 *----------------------------------------------------------------------------*/

static void
_les_balance_initialize_tui(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int keysca = cs_field_key_id("scalar_id");
  cs_dispatch_context ctx;

  int iscal = 0;

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      _les_balance.btui[iscal]->f_id = f_id;
      iscal++;
    }
  }

  /* Since every time averages for Tui were created using
     time moments, there is no need to initialize them by
     reading the les_balance restart */

  for (int isca = 0; isca < nscal; isca++) {
    cs_les_balance_tui_t *btui = _les_balance.btui[isca];

    cs_real_t *unstvar = btui->unstvar;
    cs_real_t *tptp = btui->tptp;
    cs_real_t *prodvar = btui->prodvar;
    cs_real_t *epsvar = btui->epsvar;
    cs_real_t *difftvar = btui->difftvar;
    cs_real_t *convvar = btui->convvar;
    cs_real_t *difflamvar = btui->difflamvar;

    cs_real_3_t *unstti = btui->unstti;
    cs_real_3_t *tpuip = btui->tpuip;
    cs_real_3_t *prodtUi = btui->prodtUi;
    cs_real_3_t *prodtTi = btui->prodtTi;
    cs_real_3_t *phiti = btui->phiti;
    cs_real_3_t *epsti = btui->epsti;
    cs_real_3_t *difftti = btui->difftti;
    cs_real_3_t *diffttpi = btui->diffttpi;
    cs_real_3_t *convti = btui->convti;
    cs_real_3_t *difflamti = btui->difflamti;

    cs_real_t *budsgsvar = btui->budsgsvar;
    cs_real_3_t *budsgstui = btui->budsgstui;
    cs_real_3_t **budsgstuifull = btui->budsgstuifull;

    const bool is_tui_base = (_les_balance.type & CS_LES_BALANCE_TUI_BASE);
    const bool is_tui_full = (_les_balance.type & CS_LES_BALANCE_TUI_FULL);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      unstvar[c_id] = 0.;
      tptp[c_id] = 0.;
      prodvar[c_id] = 0.;
      epsvar[c_id] = 0.;
      difftvar[c_id] = 0.;
      convvar[c_id] = 0.;
      difflamvar[c_id] = 0.;

      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        unstti[c_id][ii] = 0.;
        tpuip[c_id][ii] = 0.;
        prodtUi[c_id][ii] = 0.;
        prodtTi[c_id][ii] = 0.;
        phiti[c_id][ii] = 0.;
        epsti[c_id][ii] = 0.;
        difftti[c_id][ii] = 0.;
        diffttpi[c_id][ii] = 0.;
        convti[c_id][ii] = 0.;
        difflamti[c_id][ii] = 0.;
      }

      if (is_tui_base) {
        budsgsvar[c_id] = 0.;

        for (cs_lnum_t ii = 0; ii < 3; ii++)
          budsgstui[c_id][ii] = 0.;
      }

      if (is_tui_full)
        for (cs_lnum_t ii = 0; ii < 10; ii++)
          for (cs_lnum_t jj = 0; jj < 3; jj++)
            budsgstuifull[ii][c_id][jj] = 0.;
    });
  } /* End loop on scalars */

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Destroy a given cs_les_balance_rij_t structure pointer
 *
 * returns nullptr
 *----------------------------------------------------------------------------*/

static cs_les_balance_rij_t *
_les_balance_destroy_rij(cs_les_balance_rij_t *brij)
{
  if (brij == nullptr)
    return brij;

  CS_FREE(brij->prodij);
  CS_FREE(brij->epsij);
  CS_FREE(brij->phiij);
  CS_FREE(brij->difftij);
  CS_FREE(brij->difftpij);
  CS_FREE(brij->convij);
  CS_FREE(brij->difflamij);

  CS_FREE(brij->unstij);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {
    CS_FREE(brij->budsgsij);
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    CS_FREE(brij->budsgsfullij);
  }

  CS_FREE(brij);

  return nullptr;
}

/*----------------------------------------------------------------------------
 * Destroy an array of cs_les_balance_tui_t structure pointers
 *
 * returns nullptr
 *----------------------------------------------------------------------------*/

static cs_les_balance_tui_t **
_les_balance_destroy_tui(cs_les_balance_tui_t **btui)
{
  if (btui == nullptr)
    return btui;

  for (int isca = 0; isca < nscal; isca++) {

    CS_FREE(btui[isca]->unstvar);
    CS_FREE(btui[isca]->tptp);
    CS_FREE(btui[isca]->tpuip);
    CS_FREE(btui[isca]->prodvar);
    CS_FREE(btui[isca]->epsvar);
    CS_FREE(btui[isca]->difftvar);
    CS_FREE(btui[isca]->convvar);
    CS_FREE(btui[isca]->difflamvar);
    CS_FREE(btui[isca]->unstti);
    CS_FREE(btui[isca]->prodtUi);
    CS_FREE(btui[isca]->prodtTi);
    CS_FREE(btui[isca]->phiti);
    CS_FREE(btui[isca]->epsti);
    CS_FREE(btui[isca]->difftti);
    CS_FREE(btui[isca]->diffttpi);
    CS_FREE(btui[isca]->convti);
    CS_FREE(btui[isca]->difflamti);

    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
      CS_FREE(btui[isca]->budsgsvar);
      CS_FREE(btui[isca]->budsgstui);
    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
      for (int ii = 0; ii < 10; ii++)
        CS_FREE(btui[isca]->budsgstuifull[ii]);
      CS_FREE(btui[isca]->budsgstuifull);
      CS_FREE(btui[isca]->budsgsvarfull);
    }

    CS_FREE(btui[isca]);
  }

  CS_FREE(btui);

  return btui;
}

/*----------------------------------------------------------------------------
 * Create, allocate and initialize a the brij structure of _les_balance.
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
 * Create, allocate and initialize btui structures of _les_balance.
 *----------------------------------------------------------------------------*/

static void
_les_balance_create_tui(void)
{
  CS_MALLOC(_les_balance.btui, nscal, cs_les_balance_tui_t *);

  /* Creation and allocation of the structure containing
     working arrays */
  for (int isca = 0; isca < nscal; isca++)
    _les_balance.btui[isca] = _les_balance_allocate_tui();

  /* Initialization of working arrays */
  _les_balance_initialize_tui();
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

    CS_MALLOC(_gradt, n_scal, cs_field_t*);

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f_sca = cs_field_by_id(f_id);

      int i_sca = cs_field_get_key_int(f_sca, k_sca)-1;
      if (i_sca > -1) {
        char *name;
        int len = strlen(f_sca->name)+6;
        CS_MALLOC(name, len, char);
        snprintf(name, len, "%s_grad", f_sca->name);

        int dim = 3;
        _gradt[i_sca] = cs_field_create(name,
                                        CS_FIELD_PROPERTY,
                                        CS_MESH_LOCATION_CELLS,
                                        dim,
                                        false);
        CS_FREE(name);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 *\brief Provide access to cs_glob_les_balance
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
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_create(void)
{
  /* Rij active if Rij basic or full is active */
  if (   _les_balance.type & CS_LES_BALANCE_RIJ_BASE
      || _les_balance.type & CS_LES_BALANCE_RIJ_FULL)
    _les_balance.type |= CS_LES_BALANCE_RIJ;

  /* Tui active if Tui basic or full is active */
  if (   _les_balance.type & CS_LES_BALANCE_TUI_BASE
      || _les_balance.type & CS_LES_BALANCE_TUI_FULL)
    _les_balance.type |= CS_LES_BALANCE_TUI;

  /* Count the number of scalars nscal */
  const int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int isca = cs_field_get_key_int(f, keysca);
    if (isca > 0) nscal++;
  }

  /* I keep this comment to undestand how _PIJV2T3 is computed */
  /* Remplissage des tableaux de directions de tenseurs */
  /*for (cs_lnum_t ii = 0; ii < 3; ii++)
    for (cs_lnum_t jj = 0; jj < 3; jj++)
      for (cs_lnum_t iii = 0; iii < 6; iii++)
        if (ii*jj == idirtens[iii][0]*idirtens[iii][1])
          ipdirtens[ii][jj] = iii;

  for (cs_lnum_t ii = 0; ii < 3; ii++)
    for (cs_lnum_t jj = 0; jj < 3; jj++)
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        for (cs_lnum_t iii = 0; iii < 10; iii++)
          if (ii*jj*kk == idirtens3[iii][0]*idirtens3[iii][1]*idirtens3[iii][2])
            ipdirtens3[ii][jj][kk] = iii;*/

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
    if (f != nullptr && !cs_field_is_key_set(f, log_key_id))
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

  cs_dispatch_context ctx;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(CS_F_(vel));

  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;

  /* Rij balance structure */
  cs_les_balance_rij_t *brij = _les_balance.brij;

  cs_real_3_t *nutui = nullptr, *dnutdxkdukdxi = nullptr, *dnutdxi = nullptr;
  cs_real_6_t *nutduidxkdujdxk = nullptr, *dnutdxkuidukdxjsym = nullptr;
  cs_real_33_t *uidtaujkdxk = nullptr, *uidnutdxj = nullptr;

  /* Time moments retrieved by name */
  cs_real_t *p      =                 _les_balance_get_tm_by_name("p_m")->val;
  cs_real_t *nut    =                 _les_balance_get_tm_by_name("nut_m")->val;
  cs_real_3_t *ui   = (cs_real_3_t *)_les_balance_get_tm_by_name("ui_m")->val;
  cs_real_3_t *pu   = (cs_real_3_t *)_les_balance_get_tm_by_name("pu_m")->val;
  cs_real_6_t *uiuj = (cs_real_6_t  *)_les_balance_get_tm_by_name("uiuj_m")->val;
  cs_real_6_t *rij  = (cs_real_6_t  *)_les_balance_get_tm_by_name("u_v")->val;
  cs_real_33_t *nutduidxj
    = (cs_real_33_t *)_les_balance_get_tm_by_name("nutdjui_m")->val;
  cs_real_6_t * pduidxj
    = (cs_real_6_t  *)_les_balance_get_tm_by_name("pdjuisym_m")->val;
  cs_real_33_t * duidxj
    = (cs_real_33_t *)_les_balance_get_tm_by_name("djui_m")->val;
  cs_real_6_t *duidxkdujdxk
    = (cs_real_6_t  *)_les_balance_get_tm_by_name("dkuidkuj_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE)
    uidtaujkdxk  = (cs_real_33_t *)_les_balance_get_tm_by_name("uidktaujk_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    nutui
      = (cs_real_3_t *)_les_balance_get_tm_by_name("nutui_m")->val;
    dnutdxkdukdxi
      = (cs_real_3_t *)_les_balance_get_tm_by_name("dknutdiuk_m")->val;
    dnutdxi
      = (cs_real_3_t *)_les_balance_get_tm_by_name("dinut_m")->val;
    nutduidxkdujdxk
      = (cs_real_6_t  *)_les_balance_get_tm_by_name("nutdkuidkuj_m")->val;
    dnutdxkuidukdxjsym
      = (cs_real_6_t *)_les_balance_get_tm_by_name("dknutuidjuksym_m")->val;
    uidnutdxj
      = (cs_real_33_t *)_les_balance_get_tm_by_name("uidjnut_m")->val;
  }

  /* Get the triple corrleations mean UiUjUk */
  cs_real_t **uiujuk;
  CS_MALLOC(uiujuk, 10, cs_real_t*);

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
  cs_real_6_t  **nutdkuiuj;
  cs_real_33_t **uidujdxk;
  CS_MALLOC(nutdkuiuj, 3, cs_real_6_t*);
  CS_MALLOC(uidujdxk, 3, cs_real_33_t*);

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {
    nutdkuiuj[0] = (cs_real_6_t *)cs_field_by_name("nutd1uiuj_m")->val;
    nutdkuiuj[1] = (cs_real_6_t *)cs_field_by_name("nutd2uiuj_m")->val;
    nutdkuiuj[2] = (cs_real_6_t *)cs_field_by_name("nutd3uiuj_m")->val;

    uidujdxk[0] = (cs_real_33_t *)cs_field_by_name("u1dkuj_m")->val;
    uidujdxk[1] = (cs_real_33_t *)cs_field_by_name("u2dkuj_m")->val;
    uidujdxk[2] = (cs_real_33_t *)cs_field_by_name("u3dkuj_m")->val;
  }

  cs_field_bc_coeffs_t bc_coeffs;
  cs_field_bc_coeffs_init(&bc_coeffs);
  CS_MALLOC_HD(bc_coeffs.a, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs.b, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *coefas = bc_coeffs.a, *coefbs = bc_coeffs.b;

  cs_real_t dtref = cs_glob_time_step->dt_ref;
  cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  const int *bc_type = cs_glob_bc_type;

  cs_gradient_type_by_imrgra(eqp->imrgra,
                             &gradient_type,
                             &halo_type);

  /* Working arrays memory allocation */
  cs_real_t *diverg, *w2, *lapl;
  cs_real_3_t *w1, *w3;
  CS_MALLOC_HD(w1, n_cells, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(w2, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w3, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(lapl, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_6_t *prodij = brij->prodij;
  cs_real_6_t *epsij = brij->epsij;
  cs_real_6_t *phiij = brij->phiij;
  cs_real_6_t *difftij = brij->difftij;
  cs_real_6_t *difftpij = brij->difftpij;
  cs_real_6_t *convij = brij->convij;
  cs_real_6_t *difflamij = brij->difflamij;
  cs_real_6_t *unstij = brij->unstij;
  cs_real_6_t *budsgsij = brij->budsgsij;
  cs_real_69_t *budsgsfullij = brij->budsgsfullij;

  const int idirtens[6][2] = _IJV2T;
  const int ipdirtens[3][3] = _PIJV2T;
  const int ipdirtens3[3][3][3] = _PIJV2T3;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    /* unstij, epsij, prodij, phiij */

    for (cs_lnum_t ii = 0; ii < 6; ii++) {
      prodij[c_id][ii] = 0.;
      epsij[c_id][ii] = 0.;
      phiij[c_id][ii] = 0.;
      difftij[c_id][ii] = 0.;
      difftpij[c_id][ii] = 0.;
      convij[c_id][ii] = 0.;
      difflamij[c_id][ii] = 0.;
      unstij[c_id][ii] = (rij[c_id][ii] - unstij[c_id][ii]) / dtref;
      epsij[c_id][ii] = duidxkdujdxk[c_id][ii];

      const cs_lnum_t i = idirtens[ii][0];
      const cs_lnum_t j = idirtens[ii][1];

      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        const cs_lnum_t jj = ipdirtens[i][kk];
        const cs_lnum_t ll = ipdirtens[j][kk];

        prodij[c_id][ii] -=   rij[c_id][ll]*duidxj[c_id][i][kk]
                            + rij[c_id][jj]*duidxj[c_id][j][kk];

        epsij[c_id][ii] -= duidxj[c_id][i][kk]*duidxj[c_id][j][kk];
      }

      phiij[c_id][ii] = (  pduidxj[c_id][ii]
                         - p[c_id]*(  duidxj[c_id][i][j]
                                    + duidxj[c_id][j][i])) / ro0;

      epsij[c_id][ii] *= -2.*viscl0/ro0;
    }
  });

  /* convij */
  for (cs_lnum_t iii = 0; iii < 6; iii++) {
    const cs_lnum_t i = idirtens[iii][0];
    const cs_lnum_t j = idirtens[iii][1];

    /* convij */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        w1[c_id][kk] = ui[c_id][kk]*rij[c_id][iii];
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      convij[c_id][iii] = diverg[c_id];

      /* difftij */
      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        const cs_lnum_t jjj = ipdirtens[i][kk];
        const cs_lnum_t kkk = ipdirtens[j][kk];
        const cs_lnum_t lll = ipdirtens3[i][j][kk];
        cs_real_t *triple_corr = uiujuk[lll];

        w1[c_id][kk] = - triple_corr[c_id] - 2.*ui[c_id][i]*ui[c_id][j]*ui[c_id][kk]
                                         + ui[c_id][i]*uiuj[c_id][kkk]
                                         + ui[c_id][j]*uiuj[c_id][jjj]
                                         + ui[c_id][kk]*uiuj[c_id][iii];
      }
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      difftij[c_id][iii] = diverg[c_id];

      /* difftpij */
      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        w1[c_id][kk] = 0.;

        if (kk == i)
          w1[c_id][kk] -= pu[c_id][j] - p[c_id]*ui[c_id][j];
        else if (kk == j)
          w1[c_id][kk] -= pu[c_id][i] - p[c_id]*ui[c_id][i];
      }
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      difftpij[c_id][iii] = diverg[c_id];

      /* Laminar diffusion */
      w2[c_id] = rij[c_id][iii];
    });

    ctx.wait();

    _les_balance_laplacian(w2, lapl, 0);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      difflamij[c_id][iii] = viscl0*lapl[c_id]/ro0;
    });

    ctx.wait();
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_BASE) {

    for (cs_lnum_t iii = 0; iii < 6; iii++) {
      const cs_lnum_t i = idirtens[iii][0];
      const cs_lnum_t j = idirtens[iii][1];

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = -nutduidxj[c_id][j][kk] - nutduidxj[c_id][kk][j];
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsij[c_id][iii]
          = -(uidtaujkdxk[c_id][i][j]-ui[c_id][i]*diverg[c_id]);

        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = -nutduidxj[c_id][i][kk] - nutduidxj[c_id][kk][i];
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsij[c_id][iii] -= (  uidtaujkdxk[c_id][j][i]
                                     - ui[c_id][j]*diverg[c_id]);
        budsgsij[c_id][iii] /= ro0;
      });

      ctx.wait();
    }
  }

  if (_les_balance.type & CS_LES_BALANCE_RIJ_FULL) {

    for (cs_lnum_t iii = 0; iii < 6; iii++) {
      const cs_lnum_t i = idirtens[iii][0];
      const cs_lnum_t j = idirtens[iii][1];
      cs_real_33_t *uidujdxk_ii = (cs_real_33_t *)uidujdxk[i];
      cs_real_33_t *uidujdxk_jj = (cs_real_33_t *)uidujdxk[j];

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = nut[c_id]*( (  uidujdxk_ii[c_id][j][kk]
                                      - ui[c_id][i]*duidxj[c_id][j][kk])
                                    +(  uidujdxk_jj[c_id][i][kk]
                                      - ui[c_id][j]*duidxj[c_id][i][kk]));
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsfullij[c_id][iii][0] = diverg[c_id]/ro0;
        budsgsfullij[c_id][iii][1] = nut[c_id]/viscl0*epsij[c_id][iii];

        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          cs_real_6_t *nutdkuiuj_loc = (cs_real_6_t*)nutdkuiuj[kk];

          w1[c_id][kk] = nutdkuiuj_loc[c_id][iii]
                      + 2.*nut[c_id]*(  ui[c_id][i]*duidxj[c_id][j][kk]
                                      + ui[c_id][j]*duidxj[c_id][i][kk])
                      - nut[c_id]*(  uidujdxk_ii[c_id][j][kk]
                                   + uidujdxk_jj[c_id][i][kk])
                      - ui[c_id][i]*nutduidxj[c_id][j][kk]
                      - ui[c_id][j]*nutduidxj[c_id][i][kk]
                      - duidxj[c_id][j][kk]*nutui[c_id][i]
                      - duidxj[c_id][i][kk]*nutui[c_id][j];
        }
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsfullij[c_id][iii][2] = diverg[c_id]/ro0;

        budsgsfullij[c_id][iii][3] =  nutduidxkdujdxk[c_id][iii]
                                    - nut[c_id]*duidxkdujdxk[c_id][iii];

        cs_real_t xx = 0.;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          xx += 2.*nut[c_id]*duidxj[c_id][i][kk]*duidxj[c_id][j][kk]
              - duidxj[c_id][i][kk]*nutduidxj[c_id][j][kk]
              - duidxj[c_id][j][kk]*nutduidxj[c_id][i][kk];

        budsgsfullij[c_id][iii][3] += xx;
        budsgsfullij[c_id][iii][3] = -2./ro0*budsgsfullij[c_id][iii][3];
      });
    }

    /* Bc coeffs */
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      coefas[face_id] = 0.;
      if (   bc_type[face_id] == CS_SMOOTHWALL
          || bc_type[face_id] == CS_ROUGHWALL)
        coefbs[face_id] = 0.;
      else
        coefbs[face_id] = 1.;
    });

    ctx.wait();

    cs_gradient_scalar("nu_t",
                       gradient_type,
                       halo_type,
                       1, // inc
                       eqp->nswrgr,
                       0,
                       1,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       nullptr,
                       &bc_coeffs,
                       nut,
                       nullptr,
                       nullptr,
                       w3);

    for (cs_lnum_t iii = 0; iii < 6; iii++) {
      const cs_lnum_t i = idirtens[iii][0];
      const cs_lnum_t j = idirtens[iii][1];
      cs_real_33_t *uidujdxk_ii = (cs_real_33_t*)uidujdxk[i];
      cs_real_33_t *uidujdxk_jj = (cs_real_33_t*)uidujdxk[j];

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsfullij[c_id][iii][4] = 0.;

        budsgsfullij[c_id][iii][5] =  dnutdxkuidukdxjsym[c_id][iii]
                                    - ui[c_id][i]*dnutdxkdukdxi[c_id][j]
                                    - ui[c_id][j]*dnutdxkdukdxi[c_id][i];

        cs_real_t xx = 0.;
        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          budsgsfullij[c_id][iii][4]
            += w3[c_id][kk]
               *(  uidujdxk_ii[c_id][kk][j]-ui[c_id][i]*duidxj[c_id][kk][j]
                 + uidujdxk_jj[c_id][kk][i]-ui[c_id][j]*duidxj[c_id][kk][i] );

          xx += 2.*dnutdxi[c_id][kk]*(  ui[c_id][i]*duidxj[c_id][kk][j]
                                      + ui[c_id][j]*duidxj[c_id][kk][i])
                 - dnutdxi[c_id][kk]*(  uidujdxk_ii[c_id][kk][j]
                                      + uidujdxk_jj[c_id][kk][i])
                 - duidxj[c_id][kk][j]*uidnutdxj[c_id][i][kk]
                 - duidxj[c_id][kk][i]*uidnutdxj[c_id][j][kk];

          w1[c_id][kk]
            =  (nutui[c_id][j]-nut[c_id]*ui[c_id][j])*duidxj[c_id][i][kk]
             + (nutui[c_id][i]-nut[c_id]*ui[c_id][i])*duidxj[c_id][j][kk];
        }

        budsgsfullij[c_id][iii][4] /= ro0;

        budsgsfullij[c_id][iii][5] += xx;
        budsgsfullij[c_id][iii][5] /= ro0;
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        budsgsfullij[c_id][iii][6] = diverg[c_id]/ro0;
        budsgsfullij[c_id][iii][7] = 0.;
        budsgsfullij[c_id][iii][8] = 0.;

        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          budsgsfullij[c_id][iii][7]
            -=    (  nutduidxj[c_id][j][kk]
                   - nut[c_id]*duidxj[c_id][j][kk]) * duidxj[c_id][i][kk]
                + (  nutduidxj[c_id][i][kk]
                   - nut[c_id]*duidxj[c_id][i][kk]) * duidxj[c_id][j][kk];

          budsgsfullij[c_id][iii][8]
            +=      (uidnutdxj[c_id][i][kk]-ui[c_id][i]*dnutdxi[c_id][kk])
                  * duidxj[c_id][kk][j]
                +   (uidnutdxj[c_id][j][kk]-ui[c_id][j]*dnutdxi[c_id][kk])
                  * duidxj[c_id][kk][i];
        }

        budsgsfullij[c_id][iii][7] /= ro0;
        budsgsfullij[c_id][iii][8] /= ro0;
      });

      ctx.wait();

    } /* End loop on iii */
  } /* end test on CS_LES_BALANCE_RIJ_FULL */

  CS_FREE(uiujuk);
  CS_FREE(nutdkuiuj);
  CS_FREE(uidujdxk);

  CS_FREE(w1);
  CS_FREE(w2);
  CS_FREE(w3);
  CS_FREE(diverg);
  CS_FREE(lapl);

  CS_FREE(coefas);
  CS_FREE(coefbs);
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
  const int *bc_type = cs_glob_bc_type;

  cs_dispatch_context ctx;

  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;
  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(CS_F_(vel));

  cs_gradient_type_by_imrgra(eqp->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_real_t *tdivturflux = nullptr, *nut = nullptr, *nutt = nullptr;
  cs_real_t *nutdtdxidtdxi = nullptr;
  cs_real_3_t *nutdtdxi = nullptr;
  cs_real_3_t *uidivturflux = nullptr, *tdtauijdxj = nullptr;
  cs_real_3_t *nutui = nullptr, *nutduidxjdtdxj = nullptr;
  cs_real_3_t *dnutdxjtdujdxi = nullptr, *nuttdtdxi = nullptr;
  cs_real_3_t *dnutdxjdujdxi = nullptr, *dnutdxi = nullptr;
  cs_real_3_t *tdnutdxi = nullptr, *tdtdxi = nullptr;
  cs_real_33_t *nuttduidxj = nullptr, *nutuidtdxj = nullptr;

  cs_real_t    *p
    = (cs_real_t    *)_les_balance_get_tm_by_name("p_m")->val;
  cs_real_3_t  *ui
    = (cs_real_3_t  *)_les_balance_get_tm_by_name("ui_m")->val;
  cs_real_33_t *duidxj
    = (cs_real_33_t *)_les_balance_get_tm_by_name("djui_m")->val;
  cs_real_6_t  *uiuj
    = (cs_real_6_t  *)_les_balance_get_tm_by_name("uiuj_m")->val;
  cs_real_33_t *nutduidxj
    = (cs_real_33_t *)_les_balance_get_tm_by_name("nutdjui_m")->val;

  if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
    nut           =                _les_balance_get_tm_by_name("nut_m")->val;
    nutui         = (cs_real_3_t *)_les_balance_get_tm_by_name("nutui_m")->val;
    dnutdxi       = (cs_real_3_t *)_les_balance_get_tm_by_name("dinut_m")->val;
    dnutdxjdujdxi = (cs_real_3_t *)_les_balance_get_tm_by_name("djnutdiuj_m")->val;
  }

  cs_real_t dtref  = cs_glob_time_step->dt_ref;
  cs_real_t ro0    = cs_glob_fluid_properties->ro0;
  cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  /* Working local arrays */
  cs_real_t *diverg, *lapl, *w2;
  cs_real_3_t *w1;
  CS_MALLOC_HD(w1    , n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(w2    , n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(diverg, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(lapl  , n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_field_bc_coeffs_t bc_coeffs;
  cs_field_bc_coeffs_init(&bc_coeffs);
  CS_MALLOC_HD(bc_coeffs.a, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs.b, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *coefas = bc_coeffs.a, *coefbs = bc_coeffs.b;

  /* For each scalar */
  for (int isca = 0; isca < nscal; isca++) {

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

    cs_real_3_t  *tui
      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "tui_m")->val;
    cs_real_3_t  *ttui
      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "ttui_m")->val;
    cs_real_3_t  *dtdxi
      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "dtdxi_m")->val;
    cs_real_3_t  *pdtdxi
      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "pdtdxi_m")->val;
    cs_real_3_t  *dtdxjduidxj
      = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "djtdjui_m")->val;

    cs_real_6_t  *tuiuj
      = (cs_real_6_t *)_les_balance_get_tm_by_scalar_id(isca, "tuiuj_m")->val;

    cs_real_33_t *tduidxj
      = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "tdjui_m")->val;
    cs_real_33_t *uidtdxj
      = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "uidjt_m")->val;

    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {
      tdivturflux
        = _les_balance_get_tm_by_scalar_id(isca, "tdivturflux_m")->val;
      tdtauijdxj
        = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "tdjtauij_m")->val;
      uidivturflux
        = (cs_real_3_t *)_les_balance_get_tm_by_scalar_id(isca, "uidivturflux_m")->val;
    }

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {
      nutt
        =                 _les_balance_get_tm_by_scalar_id(isca, "nutt_m")->val;
      nutdtdxidtdxi
        =                 _les_balance_get_tm_by_scalar_id(isca, "nutditdit_m")->val;
      nutduidxjdtdxj
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nutdjuidjt_m")->val;
      dnutdxjtdujdxi
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "djnuttdiuj_m")->val;
      tdnutdxi
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "tdinut_m")->val;
      tdtdxi
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "tdit_m")->val;
      nuttdtdxi
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nuttdit_m")->val;
      nuttduidxj
        = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "nuttdjui_m")->val;
      nutuidtdxj
        = (cs_real_33_t *)_les_balance_get_tm_by_scalar_id(isca, "nutuidjt_m")->val;
    }

    if (   _les_balance.type & CS_LES_BALANCE_TUI_FULL
        || _les_balance.type & CS_LES_BALANCE_TUI_BASE   )
      nutdtdxi
        = (cs_real_3_t  *)_les_balance_get_tm_by_scalar_id(isca, "nutdit_m")->val;

    cs_real_t *unstvar = b_sca->unstvar;
    cs_real_t *tptp = b_sca->tptp;
    cs_real_t *prodvar = b_sca->prodvar;
    cs_real_t *epsvar = b_sca->epsvar;
    cs_real_t *difftvar = b_sca->difftvar;
    cs_real_t *convvar = b_sca->convvar;
    cs_real_t *difflamvar = b_sca->difflamvar;

    cs_real_3_t *unstti = b_sca->unstti;
    cs_real_3_t *tpuip = b_sca->tpuip;
    cs_real_3_t *prodtUi = b_sca->prodtUi;
    cs_real_3_t *prodtTi = b_sca->prodtTi;
    cs_real_3_t *phiti = b_sca->phiti;
    cs_real_3_t *epsti = b_sca->epsti;
    cs_real_3_t *difftti = b_sca->difftti;
    cs_real_3_t *diffttpi = b_sca->diffttpi;
    cs_real_3_t *convti = b_sca->convti;
    cs_real_3_t *difflamti = b_sca->difflamti;

    cs_real_t *budsgsvar = b_sca->budsgsvar;
    cs_real_3_t *budsgstui = b_sca->budsgstui;
    cs_real_6_t *budsgsvarfull = b_sca->budsgsvarfull;
    cs_real_3_t **budsgstuifull = b_sca->budsgstuifull;

    const int ipdirtens[3][3] = _PIJV2T;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      tptp[c_id] = t2[c_id] - cs_math_sq(t[c_id]);

      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        tpuip[c_id][ii] = tui[c_id][ii] - t[c_id]*ui[c_id][ii];
        unstti[c_id][ii] = (tpuip[c_id][ii] - unstti[c_id][ii]) / dtref;
      }

      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        prodtUi[c_id][ii]   = 0.;
        prodtTi[c_id][ii]   = 0.;
        phiti[c_id][ii]     = 0.;
        epsti[c_id][ii]     = dtdxjduidxj[c_id][ii];
        difftti[c_id][ii]   = 0.;
        diffttpi[c_id][ii]  = 0.;
        convti[c_id][ii]    = 0.;
        difflamti[c_id][ii] = 0.;

        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          cs_lnum_t iii = ipdirtens[ii][kk];

          prodtUi[c_id][ii] -= tpuip[c_id][kk]*duidxj[c_id][ii][kk];
          cs_real_t wvar = uiuj[c_id][iii] - ui[c_id][ii]*ui[c_id][kk];
          prodtTi[c_id][ii] -= wvar*dtdxi[c_id][kk];
          epsti[c_id][ii] -= dtdxi[c_id][kk]*duidxj[c_id][ii][kk];
        }

        epsti[c_id][ii] *= -xvistot/ro0;
        phiti[c_id][ii] = pdtdxi[c_id][ii]-p[c_id]*dtdxi[c_id][ii];
      }
    });

    /* convti */
    for (cs_lnum_t ii = 0; ii < 3; ii++) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = ui[c_id][kk]*tpuip[c_id][ii];
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        convti[c_id][ii] = diverg[c_id];

        /* difftti */
        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          cs_lnum_t jjj = ipdirtens[ii][kk];
          w1[c_id][kk] = - tuiuj[c_id][jjj] - 2.*t[c_id]*ui[c_id][ii]*ui[c_id][kk]
                                             + ui[c_id][ii]*tui[c_id][kk]
                                             + ui[c_id][kk]*tui[c_id][ii]
                                             + t[c_id]*uiuj[c_id][jjj];
        }
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        difftti[c_id][ii] = diverg[c_id];

        /* diffttpi */
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = tp[c_id] - t[c_id]*p[c_id];
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        diffttpi[c_id][ii] = diverg[c_id]/ro0;

        /* Laminar diffusion */
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk]
            =  viscl0*(tduidxj[c_id][ii][kk]-t[c_id]*duidxj[c_id][ii][kk])
             + visls0*(uidtdxj[c_id][ii][kk]-ui[c_id][ii]*dtdxi[c_id][kk]);
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        difflamti[c_id][ii] = diverg[c_id]/ro0;
      });
    } /* End loop on ii */

    /* Variance budgets */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      prodvar[c_id]    = 0.;
      epsvar[c_id]     = 0.;
      difftvar[c_id]   = 0.;
      convvar[c_id]    = 0.;
      difflamvar[c_id] = 0.;

      unstvar[c_id] = (tptp[c_id] - unstvar[c_id])/dtref;
      prodvar[c_id] = 0.;
      epsvar[c_id] = dtdxidtdxi[c_id];

      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        prodvar[c_id] += tpuip[c_id][kk]*dtdxi[c_id][kk];
        epsvar[c_id] -= cs_math_sq(dtdxi[c_id][kk]);
        w1[c_id][kk] = ui[c_id][kk]*tptp[c_id];
      }

      prodvar[c_id] *= -2.;
      epsvar[c_id] *= -2.*visls0/ro0;
    });

    ctx.wait();

    /* convvar */

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      convvar[c_id] = diverg[c_id];

      /* difftti */
      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        w1[c_id][kk] = ttui[c_id][kk] + 2.*cs_math_sq(t[c_id])*ui[c_id][kk]
                                    - ui[c_id][kk]*t2[c_id]
                                    - 2.*t[c_id]*tui[c_id][kk];
      }
    });

    ctx.wait();

    _les_balance_divergence_vector(w1, diverg);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      difftvar[c_id] = diverg[c_id];

      /* Laminar diffusion */
      w2[c_id] = tptp[c_id];
    });

    ctx.wait();

    _les_balance_laplacian(w2, lapl, 1);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      difflamvar[c_id] = visls0*lapl[c_id]/ro0;
    });

    if (_les_balance.type & CS_LES_BALANCE_TUI_BASE) {

      for (cs_lnum_t ii = 0; ii < 3; ii++) {

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          for (cs_lnum_t kk = 0; kk < 3; kk++)
            w1[c_id][kk] = -nutduidxj[c_id][ii][kk] - nutduidxj[c_id][kk][ii];
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgstui[c_id][ii] = -tdtauijdxj[c_id][ii]-t[c_id]*diverg[c_id];
        });
      }

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w1[c_id][kk] = -nutdtdxi[c_id][kk];
      });

      ctx.wait();

      _les_balance_divergence_vector(w1, diverg);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        /* Total SGS contribution for the variance of a scalar */
        budsgsvar[c_id]
          = -2.*(tdivturflux[c_id]-t[c_id]*diverg[c_id]/sigmas)/ro0;

        for (cs_lnum_t ii = 0; ii < 3; ii++) {
          budsgstui[c_id][ii]
            -= (uidivturflux[c_id][ii]-ui[c_id][ii]*diverg[c_id]/sigmas);

          budsgstui[c_id][ii] /= ro0;
        }
      });

    } /* End test on CS_LES_BALANCE_TUI_BASE */

    if (_les_balance.type & CS_LES_BALANCE_TUI_FULL) {

      for (cs_lnum_t ii = 0; ii < 3; ii++) {

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

          for (cs_lnum_t kk = 0; kk < 3; kk++)
            w1[c_id][kk] = nut[c_id]*( tduidxj[c_id][ii][kk]
                                    -t[c_id]*duidxj[c_id][ii][kk])
                        + nut[c_id]*( uidtdxj[c_id][ii][kk]
                                    -ui[c_id][ii]*dtdxi[c_id][kk])/sigmas;
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgstuifull[0][c_id][ii] = diverg[c_id]/ro0;
          budsgstuifull[1][c_id][ii]
            = nut[c_id]*(1.+1./sigmas)*epsti[c_id][ii]/xvistot;

          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            w1[c_id][kk] = nuttduidxj[c_id][ii][kk]
                        + 2.*nut[c_id]*t[c_id]*duidxj[c_id][ii][kk]
                        - nut[c_id]*tduidxj[c_id][ii][kk]
                        - t[c_id]*nutduidxj[c_id][ii][kk]
                        - duidxj[c_id][ii][kk]*nutt[c_id];

            cs_real_t xx =   nutuidtdxj[c_id][ii][kk]
                           + 2.*nut[c_id]*uidtdxj[c_id][ii][kk]
                           - ui[c_id][ii]*nutdtdxi[c_id][kk]
                           - dtdxi[c_id][kk]*nutui[c_id][ii];

            w1[c_id][kk] += xx/sigmas;
          }
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgstuifull[2][c_id][ii] = diverg[c_id]/ro0;

          budsgstuifull[3][c_id][ii]
            = nutduidxjdtdxj[c_id][ii] - nut[c_id]*dtdxjduidxj[c_id][ii];

          cs_real_t xx = 0.;
          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            xx = xx + 2.*nut[c_id]*duidxj[c_id][ii][kk]*dtdxi[c_id][kk]
                    - duidxj[c_id][ii][kk]*nutdtdxi[c_id][kk]
                    - dtdxi[c_id][kk]*nutduidxj[c_id][ii][kk];
          }

          budsgstuifull[3][c_id][ii] += xx;
          budsgstuifull[3][c_id][ii] *= -(1.+1./sigmas)/ro0;
        });
      } /* End loop on ii */

      /* Bc coeffs */
      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        coefas[face_id] = 0.;
        if (   bc_type[face_id] == CS_SMOOTHWALL
            || bc_type[face_id] == CS_ROUGHWALL)
          coefbs[face_id] = 0.;
        else
          coefbs[face_id] = 1.;
      });

      ctx.wait();

      cs_gradient_scalar("nu_t",
                         gradient_type,
                         halo_type,
                         1,
                         eqp->nswrgr,
                         0,
                         1,
                         eqp->verbosity,
                         static_cast<cs_gradient_limit_t>(eqp->imligr),
                         eqp->epsrgr,
                         eqp->climgr,
                         nullptr,
                         &bc_coeffs,
                         nut,
                         nullptr,
                         nullptr,
                         w1);

      for (cs_lnum_t ii = 0; ii < 3; ii++) {

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          cs_real_t xx = 0.;
          budsgstuifull[4][c_id][ii] = 0.;
          budsgstuifull[5][c_id][ii] =   dnutdxjtdujdxi[c_id][ii]
                                       - t[c_id]*dnutdxjdujdxi[c_id][ii];

          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            budsgstuifull[4][c_id][ii]
              += w1[c_id][kk]*(  tduidxj[c_id][ii][kk]
                               -t[c_id]*duidxj[c_id][ii][kk]);

            w1[c_id][kk] =     duidxj[c_id][ii][kk]
                            * (nutt[c_id]-nut[c_id]*t[c_id])
                          +   dtdxi[c_id][kk]
                            * (  nutui[c_id][ii]
                               - nut[c_id]*ui[c_id][ii]) / sigmas;

            xx = xx + 2.*t[c_id]*dnutdxi[c_id][kk]*duidxj[c_id][kk][ii]
                    - duidxj[c_id][kk][ii]*tdnutdxi[c_id][kk]
                    - dnutdxi[c_id][kk]*tduidxj[c_id][kk][ii];
          }

          budsgstuifull[4][c_id][ii] /= ro0;
          budsgstuifull[5][c_id][ii] += xx;
          budsgstuifull[5][c_id][ii] *= viscl0/ro0;
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgstuifull[6][c_id][ii] = diverg[c_id]/ro0;
          budsgstuifull[7][c_id][ii] = 0.;
          budsgstuifull[8][c_id][ii] = 0.;

          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            budsgstuifull[7][c_id][ii]
              -=   duidxj[c_id][ii][kk]*(nutdtdxi[c_id][kk]
                 - nut[c_id]*dtdxi[c_id][kk])
                 + dtdxi[c_id][kk]*(nutduidxj[c_id][ii][kk]
                 - nut[c_id]*duidxj[c_id][ii][kk])/sigmas;

            budsgstuifull[8][c_id][ii]
              += duidxj[c_id][kk][ii]*(  tdnutdxi[c_id][kk]
                                       - t[c_id]*dnutdxi[c_id][kk]);

            w1[c_id][kk] =   2.*nut[c_id]/sigmas*(tdtdxi[c_id][kk]
                           - t[c_id]*dtdxi[c_id][kk]);
          }
          budsgstuifull[7][c_id][ii] /= ro0;
          budsgstuifull[8][c_id][ii] /= ro0;
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgsvarfull[c_id][0] = diverg[c_id]/ro0;
          budsgsvarfull[c_id][1] = epsvar[c_id]*nut[c_id]/viscl0;

          for (cs_lnum_t kk = 0; kk < 3; kk++)
            w1[c_id][kk] = nuttdtdxi[c_id][kk]
                        + 2.*nut[c_id]*t[c_id]*dtdxi[c_id][kk]
                        - nut[c_id]*tdtdxi[c_id][kk]
                        - t[c_id]*nutdtdxi[c_id][kk]
                        - dtdxi[c_id][kk]*nutt[c_id];
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgsvarfull[c_id][2] = 2.*diverg[c_id]/(ro0*sigmas);
          budsgsvarfull[c_id][3]
            = nutdtdxidtdxi[c_id]-nut[c_id]*dtdxidtdxi[c_id];

          cs_real_t xx = 0.;
          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            xx += 2.*nut[c_id]*t[c_id]*dtdxi[c_id][kk]
                  -t[c_id]*nutdtdxi[c_id][kk]-dtdxi[c_id][kk]*nutt[c_id];

            w1[c_id][kk] = dtdxi[c_id][kk]*(nutt[c_id]-nut[c_id]*t[c_id]);
          }

          budsgsvarfull[c_id][3] += xx;
          budsgsvarfull[c_id][3] *= -2./(sigmas*ro0);
        });

        ctx.wait();

        _les_balance_divergence_vector(w1, diverg);

        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          budsgsvarfull[c_id][4] = 2.*diverg[c_id]/(sigmas*ro0);
          budsgsvarfull[c_id][5] = 0.;

          for (cs_lnum_t kk = 0; kk < 3; kk++) {
            budsgsvarfull[c_id][5]
              += dtdxi[c_id][kk]*(nutdtdxi[c_id][kk]-nut[c_id]*dtdxi[c_id][kk]);
          }
          budsgsvarfull[c_id][5] *= -2./(sigmas*ro0);
        });
      } /* End loop on ii */
    } /* End test on CS_LES_BALANCE_TUI_FULL */
  } /* End loop on scalars */

  ctx.wait();

  CS_FREE(w1);
  CS_FREE(w2);
  CS_FREE(diverg);
  CS_FREE(coefas);
  CS_FREE(coefbs);
  CS_FREE(lapl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write the LES balance structure in the LES balance restart file.
 */
/*----------------------------------------------------------------------------*/

void
cs_les_balance_write_restart(void)
{
  const char  restart_name[] = "les_balance.csc";
  cs_restart_t *rp = cs_restart_create(restart_name,
                                       nullptr,
                                       CS_RESTART_MODE_WRITE);

  if (rp == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the auxiliary restart "
                "file in write mode for the LES balance module.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              restart_name);

  /* Write the header */
  cs_restart_write_section(rp,
                           "les_balance_type",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_int,
                           &(_les_balance.type));

  cs_restart_destroy(&rp);
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
