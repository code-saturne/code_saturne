/*============================================================================
 * Scalar balance on zones.
 *============================================================================*/

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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_convection_diffusion.h"
#include "cs_convection_diffusion_priv.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"
#include "cs_face_viscosity.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_volume_mass_injection.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_balance_by_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_balance_by_zone.c

*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Compute convection and diffusion contributions to the flux of a scalar at
 * a boundary face.
 *
 * parameters:
 *   icvflf        -->  imposed convective flux (1: yes, 0: upwind flux)
 *   idtvar        -->  indicator of the temporal scheme
 *   iconvp        -->  convection flag
 *   idiffp        -->  diffusion flag
 *   ircflp        -->  recontruction flag
 *   relaxp        -->  relaxation coefficient
 *   diipb         -->  distance I'I'
 *   gradi         -->  gradient at boundary cell i
 *   pi            -->  value at cell i
 *   pia           -->  old value at cell i
 *   bc_type       -->  type of boundary face
 *   b_visc        -->  boundary face surface
 *   a_F           -->  explicit boundary coefficient for convection operator
 *   b_F           -->  implicit boundary coefficient for convection operator
 *   a_F           -->  explicit boundary coefficient for diffusion operator
 *   b_F           -->  implicit boundary coefficient for diffusion operator
 *   ac_F          -->  explicit imposed convective flux value (0 otherwise).
 *   bc_F          -->  implicit part of imp. conv. flux value
 *   b_mass_flux   -->  boundary mass flux
 *   xcppi         -->  specific heat value if the scalar is the temperature,
 *                      1 otherwise at cell i
 *   term_balance  <->  flux contribution
 *----------------------------------------------------------------------------*/

inline static void
_balance_boundary_faces(const int          icvflf,
                        const int          idtvar,
                        const int          iconvp,
                        const int          idiffp,
                        const int          ircflp,
                        const cs_real_t    relaxp,
                        const cs_real_3_t  diipb,
                        const cs_real_3_t  gradi,
                        const cs_real_t    pi,
                        const cs_real_t    pia,
                        const int          bc_type,
                        const cs_real_t    b_visc,
                        const cs_real_t    a_F,
                        const cs_real_t    b_F,
                        const cs_real_t    af_F,
                        const cs_real_t    bf_F,
                        const cs_real_t    ac_F,
                        const cs_real_t    bc_F,
                        const cs_real_t    b_mass_flux,
                        const cs_real_t    xcppi,
                        cs_real_t         *term_balance)
{
  /* Steady */
  if (idtvar < 0) {

    cs_real_t pir, pipr;

    cs_b_cd_steady(ircflp,
                   relaxp,
                   diipb,
                   gradi,
                   pi,
                   pia,
                   &pir,
                   &pipr);

    cs_b_imposed_conv_flux(iconvp,
                           1.,/* thetap */
                           0, /* Conservative formulation,
                                 no mass accumulation */
                           1.,
                           bc_type,
                           icvflf,
                           pi,
                           pir,
                           pipr,
                           a_F,
                           b_F,
                           ac_F,
                           bc_F,
                           b_mass_flux,
                           xcppi,
                           term_balance);

    cs_b_diff_flux(idiffp,
                   1., /* thetap */
                   1.,
                   pipr,
                   af_F,
                   bf_F,
                   b_visc,
                   term_balance);

    /* Unsteady */
  } else {

    cs_real_t pip;

    cs_b_cd_unsteady(ircflp,
                     diipb,
                     gradi,
                     pi,
                     &pip);

    cs_b_imposed_conv_flux(iconvp,
                           1.,/* thetap */
                           0, /* Conservative formulation,
                                 no mass accumulation */
                           1.,
                           bc_type,
                           icvflf,
                           pi,
                           pi, /* no relaxation */
                           pip,
                           a_F,
                           b_F,
                           ac_F,
                           bc_F,
                           b_mass_flux,
                           xcppi,
                           term_balance);

    cs_b_diff_flux(idiffp,
                   1., /* thetap */
                   1.,
                   pip,
                   af_F,
                   bf_F,
                   b_visc,
                   term_balance);
  }
}

/*----------------------------------------------------------------------------
 * Compute convection and diffusion contributions to the flux of a scalar at
 * an internal face.
 *
 * parameters:
 *   iupwin          -->  upwind scheme enabled (1: yes, 0: no)
 *   idtvar          -->  indicator of the temporal scheme
 *   iconvp          -->  convection flag
 *   idiffp          -->  diffusion flag
 *   ircflp          -->  recontruction flag
 *   ischcp          -->  second order convection scheme flag
 *   isstpp          -->  slope test flag
 *   limiter_choice  -->  choice of limiter
 *   relaxp          -->  relaxation coefficient
 *   blencp          -->  proportion of centered or SOLU scheme,
 *                        (1-blencp) is the proportion of upwind.
 *   blend_st        -->  proportion of centered or SOLU scheme,
 *                        after slope test
 *                        (1-blend_st) is the proportion of upwind.
 *   weight          -->  geometrical weight
 *   i_dist          -->  distance IJ.Nij
 *   i_face_surf     -->  face surface
 *   cell_ceni       -->  center of gravity coordinates of cell i
 *   cell_cenj       -->  center of gravity coordinates of cell j
 *   cell_cenc       -->  center of gravity coordinates of central cell
 *   cell_cend       -->  center of gravity coordinates of downwind cell
 *   i_face_u_normal -->  face unit normal
 *   i_face_cog      -->  center of gravity coordinates of face ij
 *   hybrid_blend_i  -->  blending factor between SOLU and centered
 *   hybrid_blend_j  -->  blending factor between SOLU and centered
 *   diipf           -->  distance I'I'
 *   djjpf           -->  distance J'J'
 *   gradi           -->  gradient at cell i
 *   gradj           -->  gradient at cell j
 *   gradc           -->  gradient at central cell
 *   gradupi         -->  upwind gradient at cell i
 *   gradupj         -->  upwind gradient at cell j
 *   gradsti         -->  slope test gradient at cell i
 *   gradstj         -->  slope test gradient at cell j
 *   pi              -->  value at cell i
 *   pj              -->  value at cell j
 *   pc              -->  value at central cell
 *   pd              -->  value at downwind cell
 *   pia             -->  old value at cell i
 *   pja             -->  old value at cell j
 *   i_visc          -->  diffusion coefficient (divided by IJ) at face ij
 *   i_mass_flux     -->  mass flux at face ij
 *   xcppi           -->  specific heat value if the scalar is the temperature,
 *                        1 otherwise at cell i
 *   xcppj           -->  specific heat value if the scalar is the temperature,
 *                       1 otherwise at cell j
 *   local_max       -->  local maximum of variable
 *   local_min       -->  local minimum of variable
 *   courant_c       -->  central cell courant number
 *   bi_bterms       <->  flux contribution
 *----------------------------------------------------------------------------*/

inline static void
_balance_internal_faces(int              iupwin,
                        int              idtvar,
                        int              iconvp,
                        int              idiffp,
                        int              ircflp,
                        int              ischcp,
                        int              isstpp,
                        cs_nvd_type_t    limiter_choice,
                        cs_real_t        relaxp,
                        cs_real_t        blencp,
                        cs_real_t        blend_st,
                        cs_real_t        weight,
                        cs_real_t        i_dist,
                        const cs_real_t  cell_ceni[3],
                        const cs_real_t  cell_cenj[3],
                        const cs_real_t  cell_cenc[3],
                        const cs_real_t  cell_cend[3],
                        const cs_real_t  i_face_u_normal[3],
                        const cs_real_t  i_face_cog[3],
                        cs_real_t        hybrid_blend_i,
                        cs_real_t        hybrid_blend_j,
                        const cs_real_t  diipf[3],
                        const cs_real_t  djjpf[3],
                        const cs_real_t  gradi[3],
                        const cs_real_t  gradj[3],
                        const cs_real_t  gradc[3],
                        const cs_real_t  gradupi[3],
                        const cs_real_t  gradupj[3],
                        const cs_real_t  gradsti[3],
                        const cs_real_t  gradstj[3],
                        cs_real_t        pi,
                        cs_real_t        pj,
                        cs_real_t        pc,
                        cs_real_t        pd,
                        cs_real_t        pia,
                        cs_real_t        pja,
                        cs_real_t        i_visc,
                        cs_real_t        i_mass_flux,
                        cs_real_t        xcppi,
                        cs_real_t        xcppj,
                        cs_real_t        local_max,
                        cs_real_t        local_min,
                        cs_real_t        courant_c,
                        cs_real_t        bi_bterms[2])
{
  if (iupwin == 1) {

    /* Upwind
       ====== */

    /* Steady */
    if (idtvar < 0) {

      cs_real_t pip, pjp, pipr, pjpr;
      cs_real_t pifri, pjfri, pifrj, pjfrj;

      cs_i_cd_steady_upwind(ircflp,
                            relaxp,
                            diipf,
                            djjpf,
                            gradi,
                            gradj,
                            pi,
                            pj,
                            pia,
                            pja,
                            &pifri,
                            &pifrj,
                            &pjfri,
                            &pjfrj,
                            &pip,
                            &pjp,
                            &pipr,
                            &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     0, /* Conservative formulation, no mass accumulation */
                     pi,
                     pj,
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_mass_flux,
                     xcppi,
                     xcppj,
                     bi_bterms);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc,
                     bi_bterms);

      /* Unsteady */
    } else {

      cs_real_t pip, pjp;
      cs_real_t pif, pjf;

      cs_i_cd_unsteady_upwind(ircflp,
                              diipf,
                              djjpf,
                              gradi,
                              gradj,
                              pi,
                              pj,
                              &pif,
                              &pjf,
                              &pip,
                              &pjp);

      cs_i_conv_flux(iconvp,
                     1.,
                     0, /* Conservative formulation, no mass accumulation */
                     pi,
                     pj,
                     pif,
                     pif, /* no relaxation */
                     pjf,
                     pjf, /* no relaxation */
                     i_mass_flux,
                     xcppi,
                     xcppj,
                     bi_bterms);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pip, /* no relaxation */
                     pjp, /* no relaxation */
                     i_visc,
                     bi_bterms);

    }

    /* Flux with no slope test
       ======================= */

  } else if (isstpp == 1 || isstpp == 2) {

    /* Steady */
    if (idtvar < 0) {

      cs_real_t pip, pjp, pipr, pjpr;
      cs_real_t pifri, pjfri, pifrj, pjfrj;

      cs_i_cd_steady(ircflp,
                     ischcp,
                     relaxp,
                     blencp,
                     weight,
                     cell_ceni,
                     cell_cenj,
                     i_face_cog,
                     diipf,
                     djjpf,
                     gradi,
                     gradj,
                     gradupi,
                     gradupj,
                     pi,
                     pj,
                     pia,
                     pja,
                     &pifri,
                     &pifrj,
                     &pjfri,
                     &pjfrj,
                     &pip,
                     &pjp,
                     &pipr,
                     &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     0, /* Conservative formulation, no mass accumulation */
                     pi,
                     pj,
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_mass_flux,
                     xcppi,
                     xcppj,
                     bi_bterms);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc,
                     bi_bterms);

      /* Unsteady */
    } else {

      cs_real_t pip, pjp;
      cs_real_t pif, pjf;

      if (ischcp == 4) {
        cs_i_cd_unsteady_nvd(limiter_choice,
                             blencp,
                             cell_cenc,
                             cell_cend,
                             i_face_u_normal,
                             i_face_cog,
                             gradc,
                             pc,
                             pd,
                             local_max,
                             local_min,
                             courant_c,
                             &pif,
                             &pjf);

        cs_i_conv_flux(iconvp,
                       1.,
                       0,
                       pi,
                       pj,
                       pif,
                       pif, /* no relaxation */
                       pjf,
                       pjf, /* no relaxation */
                       i_mass_flux,
                       xcppi,
                       xcppj,
                       bi_bterms);

        /* Compute required quantities for diffusive flux */
        cs_real_t recoi, recoj;

        cs_i_compute_quantities(ircflp,
                                diipf,
                                djjpf,
                                gradi,
                                gradj,
                                pi,
                                pj,
                                &recoi,
                                &recoj,
                                &pip,
                                &pjp);
      } else {
        cs_i_cd_unsteady(ircflp,
                         ischcp,
                         blencp,
                         weight,
                         cell_ceni,
                         cell_cenj,
                         i_face_cog,
                         hybrid_blend_i,
                         hybrid_blend_j,
                         diipf,
                         djjpf,
                         gradi,
                         gradj,
                         gradupi,
                         gradupj,
                         pi,
                         pj,
                         &pif,
                         &pjf,
                         &pip,
                         &pjp);

        cs_i_conv_flux(iconvp,
                       1.,
                       0, /* Conservative formulation, no mass accumulation */
                       pi,
                       pj,
                       pif,
                       pif, /* no relaxation */
                       pjf,
                       pjf, /* no relaxation */
                       i_mass_flux,
                       xcppi,
                       xcppj,
                       bi_bterms);
      }

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pip, /* no relaxation */
                     pjp, /* no relaxation */
                     i_visc,
                     bi_bterms);

    }

    /* --> Flux with slope test
       ======================== */

  } else { /* isstpp = 0 */

    /* Steady */
    if (idtvar < 0) {

      cs_real_t pip, pjp, pipr, pjpr;
      cs_real_t pifri, pjfri, pifrj, pjfrj;

      /* Upwind indicator (useless here) */
      bool indic = false;

      cs_i_cd_steady_slope_test(&indic,
                                iconvp,
                                ircflp,
                                ischcp,
                                relaxp,
                                blencp,
                                blend_st,
                                weight,
                                i_dist,
                                cell_ceni,
                                cell_cenj,
                                i_face_u_normal,
                                i_face_cog,
                                diipf,
                                djjpf,
                                i_mass_flux,
                                gradi,
                                gradj,
                                gradupi,
                                gradupj,
                                gradsti,
                                gradstj,
                                pi,
                                pj,
                                pia,
                                pja,
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

      cs_i_conv_flux(iconvp,
                     1.,
                     0, /* Conservative formulation, no mass accumulation */
                     pi,
                     pj,
                     pifri,
                     pifrj,
                     pjfri,
                     pjfrj,
                     i_mass_flux,
                     xcppi,
                     xcppj,
                     bi_bterms);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pipr,
                     pjpr,
                     i_visc,
                     bi_bterms);

      /* Unsteady */
    } else {

      cs_real_t pip, pjp;
      cs_real_t pif, pjf;

      /* Upwind indicator (useless here) */
      bool indic = false;

      cs_i_cd_unsteady_slope_test(&indic,
                                  iconvp,
                                  ircflp,
                                  ischcp,
                                  blencp,
                                  blend_st,
                                  weight,
                                  i_dist,
                                  cell_ceni,
                                  cell_cenj,
                                  i_face_u_normal,
                                  i_face_cog,
                                  diipf,
                                  djjpf,
                                  i_mass_flux,
                                  gradi,
                                  gradj,
                                  gradupi,
                                  gradupj,
                                  gradsti,
                                  gradstj,
                                  pi,
                                  pj,
                                  &pif,
                                  &pjf,
                                  &pip,
                                  &pjp);

      cs_i_conv_flux(iconvp,
                     1.,
                     0, /* Conservative formulation, no mass accumulation */
                     pi,
                     pj,
                     pif,
                     pif, /* no relaxation */
                     pjf,
                     pjf, /* no relaxation */
                     i_mass_flux,
                     xcppi,
                     xcppj,
                     bi_bterms);

      cs_i_diff_flux(idiffp,
                     1.,
                     pip,
                     pjp,
                     pip, /* no relaxation */
                     pjp, /* no relaxation */
                     i_visc,
                     bi_bterms);

    }
  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the different terms of the balance of a given scalar,
 *        on a volume zone defined by selected cell ids/
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     scalar_name         scalar name
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone_compute(const char      *scalar_name,
                           cs_lnum_t        n_cells_sel,
                           const cs_lnum_t  cell_sel_ids[],
                           cs_real_t        balance[CS_BALANCE_N_TERMS])
{
  int idtvar = cs_glob_time_step_options->idtvar;

  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_u_normal
    = (const cs_real_3_t *restrict)fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict diipf
    = (const cs_real_3_t *restrict)fvq->diipf;
  const cs_real_3_t *restrict djjpf
    = (const cs_real_3_t *restrict)fvq->djjpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  /* initialize output */

  for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
    balance[i] = 0;

  /* all boundary convective fluxes are upwind */
  int icvflb = 0; // TODO handle total energy balance
  int icvflf = 0;

  /* Get physical fields */
  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *rho = CS_F_(rho)->val;
  const cs_field_t *f = cs_field_by_name_try(scalar_name);
  const int field_id = cs_field_id_by_name(scalar_name);

  /* If the requested scalar field is not computed, return */
  if (field_id == -1) {
    bft_printf("Scalar field does not exist. Balance will not be computed.\n");
    return;
  }

  /* Internal coupling variables */
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  /* Get the calculation option from the field */
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  cs_real_t *pvar_local = NULL;
  cs_real_t *pvar_distant = NULL;
  cs_real_t  hint, rcodcl2, heq;

  const cs_lnum_t *faces_local = NULL;
  cs_lnum_t  n_local = 0;
  cs_lnum_t  n_distant = 0;
  const cs_lnum_t *faces_distant = NULL;
  cs_internal_coupling_t *cpl = NULL;

  /* Temperature indicator.
     Will multiply by CP in order to have energy. */
  const int itemperature
    = cs_field_get_key_int(f, cs_field_key_id("is_temperature"));

  /* Specific heat (CP) */
  cs_real_t *cpro_cp = NULL;
  const int icp = cs_field_id_by_name("specific_heat");
  if (itemperature) {
    if (icp != -1)
      cpro_cp = CS_F_(cp)->val;
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_cp[c_id] = cp0;
      }
    }
  }
  else {
    BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cpro_cp[c_id] = 1.;
    }
  }

  /* Internal coupling initialization*/
  if (var_cal_opt.icoupl > 0) {
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    const int coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* Zone cells selection variables*/
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = NULL;
  cs_lnum_t n_bb_faces_sel = 0;
  cs_lnum_t *bb_face_sel_ids = NULL;
  cs_lnum_t n_bi_faces_sel = 0;
  cs_lnum_t *bi_face_sel_ids = NULL;
  cs_lnum_2_t *bi_face_cells = NULL;
  cs_lnum_t *cells_tag_ids = NULL;

  cs_real_t *local_min = NULL;
  cs_real_t *local_max = NULL;
  cs_real_t *courant = NULL;

  cs_real_t *cv_limiter = NULL;
  cs_real_t *df_limiter = NULL;

  const int key_lim_choice = cs_field_key_id("limiter_choice");

  /* Initialize balance contributions
    ---------------------------------

    vol_balance   : volume contribution of unsteady terms
    div_balance   : volume contribution due to to term in div(rho u)
    mass_i_balance: contribution from mass injections
    mass_o_balance: contribution from mass suctions
    bi_i_balance  : contribution from inlet boundary faces of the selected zone
                    which are internal in the total mesh
    bi_o_balance  : contribution from outlet boundary faces of the selected zone
                    which are internal in the total mesh
    in_balance    : contribution from inlets
    out_balance   : contribution from outlets
    sym_balance   : contribution from symmetry boundaries
    s_wall_balance: contribution from smooth walls
    r_wall_balance: contribution from rough walls
    cpl_balance   : contribution from coupled faces
    i_cpl_balance : contribution from internal coupled faces
    ndef_balance  : contribution from undefined faces
    tot_balance   : total balance */

  double vol_balance = 0.;
  double tot_vol_balance2 = 0.;
  double div_balance = 0.;
  double mass_i_balance = 0.;
  double mass_o_balance = 0.;
  double bi_i_balance = 0.;
  double bi_o_balance = 0.;
  double in_balance = 0.;
  double out_balance = 0.;
  double sym_balance = 0.;
  double s_wall_balance = 0.;
  double r_wall_balance = 0.;
  double cpl_balance = 0.;
  double i_cpl_balance = 0.;
  double ndef_balance = 0.;

  /* Boundary condition coefficient for h */
  const cs_real_t *a_F = f->bc_coeffs->a;
  const cs_real_t *b_F = f->bc_coeffs->b;
  const cs_real_t *af_F = f->bc_coeffs->af;
  const cs_real_t *bf_F = f->bc_coeffs->bf;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  const int imrgra = var_cal_opt.imrgra;
  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  /* Limiters */

  int limiter_choice = -1;
  int ischcp = var_cal_opt.ischcv;
  if (field_id != -1) {
    f = cs_field_by_id(field_id);

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      limiter_choice = cs_field_get_key_int(f, key_lim_choice);
      BFT_MALLOC(local_max, n_cells_ext, cs_real_t);
      BFT_MALLOC(local_min, n_cells_ext, cs_real_t);
      cs_field_local_extrema_scalar(field_id,
                                    halo_type,
                                    local_max,
                                    local_min);
      if (limiter_choice >= CS_NVD_VOF_HRIC) {
        BFT_MALLOC(courant, n_cells_ext, cs_real_t);
        cs_cell_courant_number(field_id, courant);
      }
    }

    int cv_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("convection_limiter_id"));
    if (cv_limiter_id > -1)
      cv_limiter = cs_field_by_id(cv_limiter_id)->val;

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;
  }

  /* Allocate temporary array */
  cs_real_t *f_reconstructed;
  BFT_MALLOC(f_reconstructed, n_b_faces, cs_real_t);

  /* Reconstructed value */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  halo_type = CS_HALO_STANDARD;
  cs_field_gradient_scalar(f,
                           true, /* use_previous_t */
                           1, /* inc */
                           grad);

  if (false) { //FIXME
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      /* Associated boundary cell */
      cs_lnum_t c_id = b_face_cells[f_id];
      f_reconstructed[f_id] = f->val[c_id]
                               + grad[c_id][0]*diipb[f_id][0]
                               + grad[c_id][1]*diipb[f_id][1]
                               + grad[c_id][2]*diipb[f_id][2];
    }

  /* Non-reconstructed value */
  } else {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      /* Associated boundary cell */
      cs_lnum_t c_id = b_face_cells[f_id];
      f_reconstructed[f_id] = f->val[c_id];
    }
  }

  int inc = 1;

  /* Compute the gradient for convective scheme (the slope test, limiter, SOLU, etc) */
  cs_real_3_t *gradup = NULL;
  cs_real_3_t *gradst = NULL;
  if (var_cal_opt.blencv > 0 && var_cal_opt.isstpc == 0) {
    BFT_MALLOC(gradst, n_cells_ext, cs_real_3_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      gradst[c_id][0] = 0.;
      gradst[c_id][1] = 0.;
      gradst[c_id][2] = 0.;
    }
    /* Slope test gradient */
    if (var_cal_opt.iconv > 0)
      cs_slope_test_gradient(field_id,
                             inc,
                             halo_type,
                             (const cs_real_3_t *)grad,
                             gradst,
                             f->val,
                             a_F,
                             b_F,
                             i_mass_flux);

  }
  /* Pure SOLU scheme without using gradient_slope_test function
     or Roe and Sweby limiters */
  if (var_cal_opt.blencv > 0
      && (var_cal_opt.ischcv==2 || var_cal_opt.ischcv==4)) {
    BFT_MALLOC(gradup, n_cells_ext, cs_real_3_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      gradup[c_id][0] = 0.;
      gradup[c_id][1] = 0.;
      gradup[c_id][2] = 0.;
    }

    if (var_cal_opt.iconv > 0)
      cs_upwind_gradient(field_id,
                         inc,
                         halo_type,
                         a_F,
                         b_F,
                         i_mass_flux,
                         b_mass_flux,
                         f->val,
                         gradup);

  }

  /* Face viscosity */
  int imvisf = var_cal_opt.imvisf;
  cs_real_t *i_visc;
  cs_real_t *b_visc;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  cs_real_t *c_visc = NULL;
  BFT_MALLOC(c_visc, n_cells_ext, cs_real_t);
  const int kivisl
    = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
  if (kivisl != -1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] = cs_field_by_id(kivisl)->val[c_id];
  }
  else {
    const double visls0
      = cs_field_get_key_double(f, cs_field_key_id("diffusivity_ref"));
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      c_visc[c_id] = visls0;
    }
  }

  /* Turbulent part */
  cs_real_t *c_visct = cs_field_by_name("turbulent_viscosity")->val;

  if (var_cal_opt.idifft == 1) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] += cpro_cp[c_id] * c_visct[c_id]/turb_schmidt;
  }

  cs_face_viscosity(m, fvq, imvisf, c_visc, i_visc, b_visc);

  /* Get user-selected zone
     ====================== */

  /* Initialize arrays */

  /* Internal faces of the selected zone */
  BFT_MALLOC(i_face_sel_ids, n_i_faces, cs_lnum_t);
  /* Boundary faces of the selected zone,
     which are internal faces of the global mesh.
     Faces -> cells connectivity */
  BFT_MALLOC(bi_face_sel_ids, n_i_faces, cs_lnum_t);
  BFT_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_face_sel_ids[f_id] = -1;
    bi_face_sel_ids[f_id] = -1;
    bi_face_cells[f_id][0] = -999;
    bi_face_cells[f_id][1] = -999;
  }

  /* Boundary faces of the selected zone,
     which are also boundary faces of the global mesh */
  BFT_MALLOC(bb_face_sel_ids, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bb_face_sel_ids[f_id] = -1;
  }

  /* Synchronization for parallelism */
  BFT_MALLOC(cells_tag_ids, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    cells_tag_ids[c_id] = 0;
  }
  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {
    cs_lnum_t c_id_sel = cell_sel_ids[c_id];
    cells_tag_ids[c_id_sel] = 1;
  }
  if (halo != NULL) {
    cs_halo_sync_num(halo, halo_type, cells_tag_ids);
  }

  /* Classify mesh faces with respect to the selected zone */

  /* Check boundary faces:
     if they are in the selected zone, they are boundary as well */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id = b_face_cells[f_id];

    if (cells_tag_ids[c_id] == 1) {
      n_bb_faces_sel++;
      bb_face_sel_ids[n_bb_faces_sel-1] = f_id;
    }
  }

  /* Check internal faces:
     if they are in the selected zone, they can be either
     internal or boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    bool indic1 = false;
    bool indic2 = false;

    if (cells_tag_ids[c_id1] == 1)
      indic1 = true;
    if (cells_tag_ids[c_id2] == 1)
      indic2 = true;

    if (indic1 && indic2) {
      n_i_faces_sel++;
      i_face_sel_ids[n_i_faces_sel-1] = f_id;
    }
    else if (indic1 || indic2) {
      bi_face_sel_ids[n_bi_faces_sel] = f_id;
      n_bi_faces_sel++;
      /* Build the faces -> cells connectivity as done in
         i_face_cells */
      if (indic1)
        bi_face_cells[f_id][0] = c_id1;
      else
        bi_face_cells[f_id][1] = c_id2;
    }

  }

  /* Compute the balance at time step n
    ===================================

    --> Balance on interior volumes and
        total quantity on interior volumes
        ---------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {

    cs_lnum_t c_id_sel = cell_sel_ids[c_id];

    vol_balance += cell_vol[c_id_sel] * rho[c_id_sel]
                 * cpro_cp[c_id_sel]
                 * (f->val_pre[c_id_sel] - f->val[c_id_sel]);

    cs_real_t rho_y_dt =  rho[c_id_sel] * cpro_cp[c_id_sel]
                        * f->val_pre[c_id_sel] * dt[c_id_sel];
    tot_vol_balance2 += cell_vol[c_id_sel] * rho_y_dt * rho_y_dt;
  }

  /* Balance on all faces (interior and boundary), for div(rho u)
     ------------------------------------------------------------ */

  /* Interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
    /* Associated internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    /* Contribution to flux from the two cells of the current face
      (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */

    if (c_id1 < n_cells)
      div_balance += i_mass_flux[f_id_sel] * dt[c_id1] * f->val[c_id1]
                   * cpro_cp[c_id1];

    if (c_id2 < n_cells)
      div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                   * cpro_cp[c_id2];

  }

  /* Boundary faces which are internal in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = bi_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = bi_face_cells[f_id_sel][1];

    /* Contribution to flux from the only cell of the current face
       lying inside the selected zone
      (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */

    if (c_id1 >= 0) {

      if (c_id1 < n_cells)
        div_balance += i_mass_flux[f_id_sel] * dt[c_id1] * f->val[c_id1]
                     * cpro_cp[c_id1];
    }

    else {

      if (c_id2 < n_cells)
        div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                     * cpro_cp[c_id2];
    }

  }

  /* Boundary faces which are also boundary in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    /* Contribution to flux from the current face */
      div_balance += b_mass_flux[f_id_sel] * dt[c_id] * f->val[c_id]
                   * cpro_cp[c_id];

  }

  /* Mass source terms and mass accumulation term.
     In case of a mass source term, add contribution from Gamma*Tn+1 */

  cs_lnum_t ncesmp = 0;
  cs_lnum_t *icetsm = NULL;
  int *itpsmp = NULL;
  cs_real_t *smcelp, *gamma = NULL;

  cs_volume_mass_injection_get_arrays(f, &ncesmp, &icetsm, &itpsmp,
                                      &smcelp, &gamma);

  if (ncesmp > 0) {

    const cs_real_t *cell_f_vol = fvq->cell_f_vol;
    const double cp0 = cs_glob_fluid_properties->cp0;

    for (cs_lnum_t c_idx = 0; c_idx < ncesmp; c_idx++) {
      cs_lnum_t c_id_sel = icetsm[c_idx] - 1;

      if (cells_tag_ids[c_id_sel]) {

        cs_real_t vg = gamma[c_idx];
        cs_real_t v;
        if (itpsmp[c_idx] == 0 || vg < 0)
          v = f->val[c_id_sel];
        else
          v = smcelp[c_idx];

        cs_real_t c_st = cell_f_vol[c_id_sel] * dt[c_id_sel]* vg * v;

        if (itemperature) {
          if (icp >= 0)
            c_st *= cpro_cp[c_id_sel];
          else
            c_st *= cp0;
        }

        if (vg < 0)
          mass_o_balance += c_st;
        else
          mass_i_balance += c_st;
      }

    }
  }

  int iconvp = var_cal_opt.iconv;
  int idiffp = var_cal_opt.idiff;
  int ircflp = var_cal_opt.ircflu;
  double relaxp = var_cal_opt.relaxv;

  /* Balance on boundary faces
     -------------------------

     We handle different types of boundary faces separately to better
     analyze the information, but this is not mandatory. */

  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    cs_real_t term_balance = 0.;

    cs_real_t ac_F = 0.;
    cs_real_t bc_F = 0.;

    if (icvflb == 1) {
      ac_F = f->bc_coeffs->ac[f_id_sel];
      bc_F = f->bc_coeffs->bc[f_id_sel];
      icvflf = 0; /* = icvfli */
    }

    _balance_boundary_faces(icvflf,
                            idtvar,
                            iconvp,
                            idiffp,
                            ircflp,
                            relaxp,
                            diipb[f_id_sel],
                            grad[c_id],
                            f->val[c_id],
                            f->val_pre[c_id],
                            bc_type[f_id_sel],
                            b_visc[f_id_sel],
                            a_F[f_id_sel],
                            b_F[f_id_sel],
                            af_F[f_id_sel],
                            bf_F[f_id_sel],
                            ac_F,
                            bc_F,
                            b_mass_flux[f_id_sel],
                            cpro_cp[c_id],
                            &term_balance);

    if (bc_type[f_id_sel] == CS_INLET ||
        bc_type[f_id_sel] == CS_FREE_INLET ||
        bc_type[f_id_sel] == CS_ESICF ||
        bc_type[f_id_sel] == CS_EPHCF)
      in_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_OUTLET ||
             bc_type[f_id_sel] == CS_SSPCF ||
             bc_type[f_id_sel] == CS_SOPCF)
      out_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_SYMMETRY)
      sym_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_SMOOTHWALL)
      s_wall_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_ROUGHWALL)
      r_wall_balance -= term_balance*dt[c_id];
    else if (   bc_type[f_id_sel] == CS_COUPLED
             || bc_type[f_id_sel] == CS_COUPLED_FD)
      cpl_balance -= term_balance*dt[c_id];
    else
      ndef_balance -= term_balance*dt[c_id];

  }

  /* Balance on coupled faces
     ------------------------

     We handle different types of boundary faces separately to better
     analyze the information, but this is not mandatory. */

  if (var_cal_opt.icoupl > 0) {

    /* Prepare data for sending from distant */
    BFT_MALLOC(pvar_distant, n_distant, cs_real_t);

    for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
      cs_lnum_t f_id = faces_distant[ii];
      cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t pip;
      cs_b_cd_unsteady(ircflp,
                       diipb[f_id],
                       grad[c_id],
                       f->val[c_id],
                       &pip);
      pvar_distant[ii] = pip;
    }

    /* Receive data */
    BFT_MALLOC(pvar_local, n_local, cs_real_t);
    cs_internal_coupling_exchange_var(cpl,
                                      1, /* Dimension */
                                      pvar_distant,
                                      pvar_local);

    /* flux contribution */
    for (cs_lnum_t ii = 0; ii < n_local; ii++) {
      cs_lnum_t f_id = faces_local[ii];
      cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t surf = b_face_surf[f_id];

      if (cells_tag_ids[c_id] == 1) {
        cs_real_t pip, pjp;
        cs_real_t term_balance = 0.;

        cs_b_cd_unsteady(ircflp,
                         diipb[f_id],
                         grad[c_id],
                         f->val[c_id],
                         &pip);

        pjp = pvar_local[ii];

        hint = f->bc_coeffs->hint[f_id];
        rcodcl2 = f->bc_coeffs->rcodcl2[f_id];
        heq = surf * hint * rcodcl2 / (hint + rcodcl2);

        cs_b_diff_flux_coupling(idiffp,
                                pip,
                                pjp,
                                heq,
                                &term_balance);

        i_cpl_balance -= term_balance*dt[c_id];
      }
    }

    BFT_FREE(pvar_local);
    BFT_FREE(pvar_distant);

  }

  /* Balance on boundary faces of the selected zone
     that are interior to the total mesh
     ---------------------------------------------- */

  int isstpp = var_cal_opt.isstpc;
  double blencp = var_cal_opt.blencv;
  double blend_st = var_cal_opt.blend_st;
  int iupwin = (blencp > 0.) ? 0 : 1;

  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    cs_real_t beta = blencp;
    /* Beta blending coefficient ensuring positivity of the scalar */
    if (isstpp == 2) {
      beta = CS_MAX(CS_MIN(cv_limiter[c_id1], cv_limiter[c_id2]), 0.);
    }

    int bldfrp = ircflp;
    /* Local limitation of the reconstruction */
    if (df_limiter != NULL && ircflp > 0) {
      cs_real_t _bldfrp = fmax(fmin(df_limiter[c_id1],
                                    df_limiter[c_id2]), 0.);
      bldfrp = (int)_bldfrp;
    }

    cs_real_t hybrid_coef_ii, hybrid_coef_jj;
    cs_lnum_t ic = -1, id = -1;
    cs_real_t courant_c = -1., _local_max = 0., _local_min = 0.;
    if (ischcp == 3) {
      hybrid_coef_ii = CS_F_(hybrid_blend)->val[c_id1];
      hybrid_coef_jj = CS_F_(hybrid_blend)->val[c_id2];
    }
    else if (ischcp == 4) {
      hybrid_coef_ii = 0.;
      hybrid_coef_jj = 0.;
      /* Determine central and downwind sides w.r.t. current face */
      cs_central_downwind_cells(c_id1,
                                c_id2,
                                i_mass_flux[f_id_sel],
                                &ic,  /* central cell id */
                                &id); /* downwind cell id */

      if (courant != NULL)
        courant_c = courant[ic];

      if (local_max != NULL) {
        _local_max = local_max[ic];
        _local_min = local_min[ic];
      }
    }
    else {
      hybrid_coef_ii = 0.;
      hybrid_coef_jj = 0.;
    }

    cs_real_2_t bi_bterms = {0., 0.};

    _balance_internal_faces(iupwin,
                            idtvar,
                            iconvp,
                            idiffp,
                            bldfrp,
                            ischcp,
                            isstpp,
                            (cs_nvd_type_t)limiter_choice,
                            relaxp,
                            beta,
                            blend_st,
                            weight[f_id_sel],
                            i_dist[f_id_sel],
                            cell_cen[c_id1],
                            cell_cen[c_id2],
                            cell_cen[ic],
                            cell_cen[id],
                            i_face_u_normal[f_id_sel],
                            i_face_cog[f_id_sel],
                            hybrid_coef_ii,
                            hybrid_coef_jj,
                            diipf[f_id_sel],
                            djjpf[f_id_sel],
                            grad[c_id1],
                            grad[c_id2],
                            grad[ic],
                            gradup[c_id1],
                            gradup[c_id2],
                            gradst[c_id1],
                            gradst[c_id2],
                            f->val[c_id1],
                            f->val[c_id2],
                            f->val[ic],
                            f->val[id],
                            f->val_pre[c_id1],
                            f->val_pre[c_id2],
                            i_visc[f_id_sel],
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            _local_max,
                            _local_min,
                            courant_c,
                            bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0)
          bi_o_balance -= bi_bterms[0]*dt[c_id1];
        else
          bi_i_balance -= bi_bterms[0]*dt[c_id1];
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0)
          bi_i_balance += bi_bterms[1]*dt[c_id2];
        else
          bi_o_balance += bi_bterms[1]*dt[c_id2];
      }
    }

  }

  /* Free memory */

  BFT_FREE(grad);
  BFT_FREE(gradup);
  BFT_FREE(gradst);
  BFT_FREE(f_reconstructed);
  BFT_FREE(local_max);
  BFT_FREE(local_min);
  BFT_FREE(courant);

  if (!itemperature || icp == -1)
    BFT_FREE(cpro_cp);
  BFT_FREE(c_visc);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);

  BFT_FREE(cells_tag_ids);
  BFT_FREE(bi_face_cells);
  BFT_FREE(i_face_sel_ids);
  BFT_FREE(bb_face_sel_ids);
  BFT_FREE(bi_face_sel_ids);

  /* Sum of values on all ranks (parallel calculations) */

  balance[CS_BALANCE_TOTAL_NORMALIZED] = tot_vol_balance2; /* temporary */

  balance[CS_BALANCE_VOLUME] = vol_balance;
  balance[CS_BALANCE_DIV] = div_balance;
  balance[CS_BALANCE_UNSTEADY] = vol_balance + div_balance;
  balance[CS_BALANCE_MASS] = mass_i_balance + mass_o_balance;
  balance[CS_BALANCE_MASS_IN] = mass_i_balance;
  balance[CS_BALANCE_MASS_OUT] = mass_o_balance;
  balance[CS_BALANCE_INTERIOR_IN] = bi_i_balance;
  balance[CS_BALANCE_INTERIOR_OUT] = bi_o_balance;
  balance[CS_BALANCE_BOUNDARY_IN] = in_balance;
  balance[CS_BALANCE_BOUNDARY_OUT] = out_balance;
  balance[CS_BALANCE_BOUNDARY_SYM] = sym_balance;
  balance[CS_BALANCE_BOUNDARY_WALL] = s_wall_balance + r_wall_balance;
  balance[CS_BALANCE_BOUNDARY_WALL_S] = s_wall_balance;
  balance[CS_BALANCE_BOUNDARY_WALL_R] = r_wall_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED] = cpl_balance + i_cpl_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED_E] = cpl_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED_I] = i_cpl_balance;
  balance[CS_BALANCE_BOUNDARY_OTHER] = ndef_balance;

  cs_parall_sum(CS_BALANCE_N_TERMS, CS_REAL_TYPE, balance);

  /* Total balance: add the different contributions calculated above */

  balance[CS_BALANCE_TOTAL]
    =   balance[CS_BALANCE_UNSTEADY] + balance[CS_BALANCE_MASS]
      + balance[CS_BALANCE_INTERIOR_IN] + balance[CS_BALANCE_INTERIOR_OUT]
      + balance[CS_BALANCE_BOUNDARY_IN] + balance[CS_BALANCE_BOUNDARY_OUT]
      + balance[CS_BALANCE_BOUNDARY_SYM] + balance[CS_BALANCE_BOUNDARY_WALL]
      + balance[CS_BALANCE_BOUNDARY_COUPLED]
      + balance[CS_BALANCE_BOUNDARY_OTHER];

  tot_vol_balance2 = balance[CS_BALANCE_TOTAL_NORMALIZED]; /* from temporary above */
  balance[CS_BALANCE_TOTAL_NORMALIZED] = balance[CS_BALANCE_TOTAL];

  if (tot_vol_balance2 > 0.)
    balance[CS_BALANCE_TOTAL_NORMALIZED] /= sqrt(tot_vol_balance2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and log the different terms of the balance of a given scalar,
 *        on a volumic zone defined by selection criteria.
 *        The different contributions to the balance are printed in the
 *        run_solver.log.
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * The scalar and the zone are selected at the top of the routine
 * by the user.
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const char  *selection_crit,
                   const char  *scalar_name)
{
  cs_real_t balance[CS_BALANCE_N_TERMS];

  const cs_mesh_t *m = cs_glob_mesh;
  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Select cells */

  cs_lnum_t n_cells_sel = 0;
  cs_lnum_t *cells_sel_ids = NULL;

  BFT_MALLOC(cells_sel_ids, m->n_cells, cs_lnum_t);
  cs_selector_get_cell_list(selection_crit, &n_cells_sel, cells_sel_ids);

  /* Compute balance */

  cs_balance_by_zone_compute(scalar_name,
                             n_cells_sel,
                             cells_sel_ids,
                             balance);

  BFT_FREE(cells_sel_ids);

  /* Log results at time step n */

  bft_printf
    (_("   ** SCALAR BALANCE BY ZONE at iteration %6i\n"
       "   ---------------------------------------------\n"
       "------------------------------------------------------------\n"
       "   SCALAR: %s\n"
       "   ZONE SELECTION CRITERIA: \"%s\"\n"
       "------------------------------------------------------------\n"
       "   Unst. term   Inj. Mass.   Suc. Mass.\n"
       "  %12.4e %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   IB inlet     IB outlet\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Inlet        Outlet\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Sym.         Smooth W.    Rough W.\n"
       "  %12.4e %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Coupled      Int. Coupling    Undef. BC\n"
       "  %12.4e %12.4e     %12.4e\n"
       "------------------------------------------------------------\n"
       "   Total        Instant. norm. total\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n\n"),
     nt_cur, scalar_name, selection_crit,
     balance[CS_BALANCE_UNSTEADY],
     balance[CS_BALANCE_MASS_IN], balance[CS_BALANCE_MASS_OUT],
     balance[CS_BALANCE_INTERIOR_IN], balance[CS_BALANCE_INTERIOR_OUT],
     balance[CS_BALANCE_BOUNDARY_IN], balance[CS_BALANCE_BOUNDARY_OUT],
     balance[CS_BALANCE_BOUNDARY_SYM],
     balance[CS_BALANCE_BOUNDARY_WALL_S], balance[CS_BALANCE_BOUNDARY_WALL_R],
     balance[CS_BALANCE_BOUNDARY_COUPLED_E],
     balance[CS_BALANCE_BOUNDARY_COUPLED_I],
     balance[CS_BALANCE_BOUNDARY_OTHER],
     balance[CS_BALANCE_TOTAL], balance[CS_BALANCE_TOTAL_NORMALIZED]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 *        on a volume zone defined by selected cell ids/
 *
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_p_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone_compute(cs_lnum_t        n_cells_sel,
                                 const cs_lnum_t  cell_sel_ids[],
                                 cs_real_t        balance[CS_BALANCE_P_N_TERMS])
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict diipf
    = (const cs_real_3_t *restrict)fvq->diipf;
  const cs_real_3_t *restrict djjpf
    = (const cs_real_3_t *restrict)fvq->djjpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  /* initialize output */

  for (int i = 0; i < CS_BALANCE_P_N_TERMS; i++)
    balance[i] = 0;

  /* Get physical fields */
  const cs_real_t *rho = CS_F_(rho)->val;
  const cs_field_t *f_pres = CS_F_(p);
  const cs_real_t *pressure = f_pres->val;
  const cs_field_t *f_vel = CS_F_(vel);
  const cs_real_3_t *velocity =  (const cs_real_3_t *)f_vel->val;
  const cs_real_t gravity[3] = {cs_glob_physical_constants->gravity[0],
                                cs_glob_physical_constants->gravity[1],
                                cs_glob_physical_constants->gravity[2]};

  /* Zone cells selection variables*/
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = NULL;
  cs_lnum_t n_bb_faces_sel = 0;
  cs_lnum_t *bb_face_sel_ids = NULL;
  cs_lnum_t n_bi_faces_sel = 0;
  cs_lnum_t *bi_face_sel_ids = NULL;
  cs_lnum_2_t *bi_face_cells = NULL;
  cs_lnum_t *cells_tag_ids = NULL;

  /* Initialization of balance contributions
     ---------------------------------------

    in_pressure   : contribution from inlets
    out_pressure  : contribution from outlets
    in_u2         : contribution from inlets
    out_u2        : contribution from outlets
    in_rhogx      : contribution from inlets
    out_rhogx     : contribution from outlets
    in_debit      : debit from inlets
    out_debit     : debit from outlets
    in_m_debit    : mass flow from inlets
    out_m_debit   : mass flow from outlets

  */

  double in_pressure= 0.;
  double out_pressure= 0.;
  double in_u2 = 0.;
  double out_u2 = 0.;
  double in_rhogx = 0.;
  double out_rhogx = 0.;
  double in_debit = 0.;
  double out_debit = 0.;
  double in_m_debit = 0.;
  double out_m_debit = 0.;

  /* Boundary condition coefficient for p */
  const cs_real_t *a_p = f_pres->bc_coeffs->a;
  const cs_real_t *b_p = f_pres->bc_coeffs->b;

  /* Boundary condition coefficient for u */
  const cs_real_3_t *a_u = (const cs_real_3_t *)f_vel->bc_coeffs->a;
  const cs_real_33_t *b_u = (const cs_real_33_t *)f_vel->bc_coeffs->b;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = cs_field_get_key_int(f_pres, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = cs_field_get_key_int(f_pres, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  int inc = 1;

  /* Get user-selected zone
     ====================== */

  /* Initialize arrays */

  /* Internal faces of the selected zone */
  BFT_MALLOC(i_face_sel_ids, n_i_faces, cs_lnum_t);
  /* Boundary faces of the selected zone,
     which are internal faces of the global mesh.
     Faces -> cells connectivity */
  BFT_MALLOC(bi_face_sel_ids, n_i_faces, cs_lnum_t);
  BFT_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_face_sel_ids[f_id] = -1;
    bi_face_sel_ids[f_id] = -1;
    bi_face_cells[f_id][0] = -999;
    bi_face_cells[f_id][1] = -999;
  }

  /* Boundary faces of the selected zone,
     which are also boundary faces of the global mesh */
  BFT_MALLOC(bb_face_sel_ids, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bb_face_sel_ids[f_id] = -1;
  }


  /* Synchronization for parallelism */
  BFT_MALLOC(cells_tag_ids, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    cells_tag_ids[c_id] = 0;
  }
  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {
    cs_lnum_t c_id_sel = cell_sel_ids[c_id];
    cells_tag_ids[c_id_sel] = 1;
  }
  if (halo != NULL) {
    cs_halo_sync_num(halo, CS_HALO_STANDARD, cells_tag_ids);
  }

  /* Classify mesh faces with respect to the selected zone */

  /* Check boundary faces:
     if they are in the selected zone, they are boundary as well */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id = b_face_cells[f_id];

    if (cells_tag_ids[c_id] == 1) {
      n_bb_faces_sel++;
      bb_face_sel_ids[n_bb_faces_sel-1] = f_id;
    }
  }

  /* Check internal faces:
     if they are in the selected zone, they can be either
     internal or boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    bool indic1 = false;
    bool indic2 = false;

    if (cells_tag_ids[c_id1] == 1)
      indic1 = true;
    if (cells_tag_ids[c_id2] == 1)
      indic2 = true;

    if (indic1 && indic2) {
      n_i_faces_sel++;
      i_face_sel_ids[n_i_faces_sel-1] = f_id;
    }
    else if (indic1 || indic2) {
      n_bi_faces_sel++;
      bi_face_sel_ids[n_bi_faces_sel-1] = f_id;
      /* Build the faces -> cells connectivity as done in
         i_face_cells */
      if (indic1)
        bi_face_cells[f_id][0] = c_id1;
      else
        bi_face_cells[f_id][1] = c_id2;
    }

  }

  /* Balance computation
     =================== */

  /* Compute the balance at time step n */

  int iconvp = 1;
  int ircflp = 0; /* No reconstruction */

  /* Balance on boundary faces
     -------------------------

     We handle different types of boundary faces separately to better
     analyze the information, but this is not mandatory. */

  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    cs_real_t pip;

    /* Pressure term FIXME rho0*gravity*(X-X0) should be added */
    cs_real_t p_rho = pressure[c_id] / rho[c_id];
    cs_real_t a_p_rho = a_p[f_id_sel] / rho[c_id];
    cs_real_t b_p_rho = b_p[f_id_sel];

    cs_real_3_t grad = {0, 0, 0};

    cs_b_cd_unsteady(ircflp,
                     diipb[f_id_sel],
                     grad,
                     p_rho,
                     &pip);

    cs_real_t term_balance = 0.;

    cs_b_upwind_flux(iconvp,
                     1., /* thetap */
                     0, /* Conservative formulation, no mass accumulation */
                     inc,
                     bc_type[f_id_sel],
                     p_rho,
                     p_rho, /* no relaxation */
                     pip,
                     a_p_rho,
                     b_p_rho,
                     b_mass_flux[f_id_sel],
                     1.,
                     &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_debit += b_mass_flux[f_id_sel]/rho[c_id];
      out_m_debit += b_mass_flux[f_id_sel];
      out_pressure += term_balance;
    } else {
      in_debit += b_mass_flux[f_id_sel]/rho[c_id];
      in_m_debit += b_mass_flux[f_id_sel];
      in_pressure += term_balance;
    }

    /* Kinematic term */
    cs_real_t u2 = 0.5 * cs_math_3_square_norm(velocity[c_id]);
    cs_real_t a_u2 = 0.5 * cs_math_3_square_norm(a_u[f_id_sel]);
    /* Approximation of u^2 BC */
    cs_real_t b_u2 = 1./6.*( b_u[f_id_sel][0][0] * b_u[f_id_sel][0][0]
                           + b_u[f_id_sel][1][1] * b_u[f_id_sel][1][1]
                           + b_u[f_id_sel][2][2] * b_u[f_id_sel][2][2]);

    cs_b_cd_unsteady(ircflp,
                     diipb[f_id_sel],
                     grad,
                     u2,
                     &pip);

    term_balance = 0.;

    cs_b_upwind_flux(iconvp,
                     1., /* thetap */
                     0, /* Conservative formulation, no mass accumulation */
                     inc,
                     bc_type[f_id_sel],
                     u2,
                     u2, /* no relaxation */
                     pip,
                     a_u2,
                     b_u2,
                     b_mass_flux[f_id_sel],
                     1.,
                     &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_u2 += term_balance;
    } else {
      in_u2 += term_balance;
    }

    /* Gravity term */
    cs_real_t gx = - cs_math_3_dot_product(gravity, b_face_cog[f_id_sel]);
    /* Trivial BCs */
    cs_real_t a_gx = gx;
    cs_real_t b_gx = 0.;

    cs_b_cd_unsteady(ircflp,
                     diipb[f_id_sel],
                     grad,
                     gx,
                     &pip);

    term_balance = 0.;

    cs_b_upwind_flux(iconvp,
                     1., /* thetap */
                     0, /* Conservative formulation, no mass accumulation */
                     inc,
                     bc_type[f_id_sel],
                     gx,
                     gx, /* no relaxation */
                     pip,
                     a_gx,
                     b_gx,
                     b_mass_flux[f_id_sel],
                     1.,
                     &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_rhogx += term_balance;
    } else {
      in_rhogx += term_balance;
    }

  }

  /* Balance on boundary faces of the selected zone
     that are internal of the total mesh
     ---------------------------------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    cs_real_2_t bi_bterms = {0.,0.};
    cs_real_3_t grad = {0, 0, 0};

    cs_real_t pip, pjp;
    cs_real_t pif, pjf;

    /* Pressure term */
    cs_real_t p_rho_id1 = pressure[c_id1] / rho[c_id1];
    cs_real_t p_rho_id2 = pressure[c_id2] / rho[c_id2];

    cs_i_cd_unsteady_upwind(ircflp,
                            diipf[f_id_sel],
                            djjpf[f_id_sel],
                            grad,
                            grad,
                            p_rho_id1,
                            p_rho_id2,
                            &pif,
                            &pjf,
                            &pip,
                            &pjp);

    cs_i_conv_flux(iconvp,
                   1.,
                   0, /* Conservative formulation, no mass accumulation */
                   p_rho_id1,
                   p_rho_id2,
                   pif,
                   pif, /* no relaxation */
                   pjf,
                   pjf, /* no relaxation */
                   i_mass_flux[f_id_sel],
                   1.,
                   1.,
                   bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_pressure += bi_bterms[0];
          out_debit += i_mass_flux[f_id_sel] / rho[c_id1];
          out_m_debit += i_mass_flux[f_id_sel];
        } else {
          in_pressure += bi_bterms[0];
          in_debit += i_mass_flux[f_id_sel] / rho[c_id1];
          in_m_debit += i_mass_flux[f_id_sel];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_pressure -= bi_bterms[1];
          in_debit -= i_mass_flux[f_id_sel] / rho[c_id2];
          in_m_debit -= i_mass_flux[f_id_sel];
        } else {
          out_pressure -= bi_bterms[1];
          out_debit -= i_mass_flux[f_id_sel] / rho[c_id2];
          out_m_debit -= i_mass_flux[f_id_sel];
        }
      }
    }

    /* Kinematic term */
    bi_bterms[0] = 0.;
    bi_bterms[1] = 0.;

    cs_real_t u2_id1 = 0.5 * cs_math_3_square_norm(velocity[c_id1]);
    cs_real_t u2_id2 = 0.5 * cs_math_3_square_norm(velocity[c_id2]);

    cs_i_cd_unsteady_upwind(ircflp,
                            diipf[f_id_sel],
                            djjpf[f_id_sel],
                            grad,
                            grad,
                            u2_id1,
                            u2_id2,
                            &pif,
                            &pjf,
                            &pip,
                            &pjp);

    cs_i_conv_flux(iconvp,
                   1.,
                   0, /* Conservative formulation, no mass accumulation */
                   u2_id1,
                   u2_id2,
                   pif,
                   pif, /* no relaxation */
                   pjf,
                   pjf, /* no relaxation */
                   i_mass_flux[f_id_sel],
                   1.,
                   1.,
                   bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_u2 += bi_bterms[0];
        } else {
          in_u2 += bi_bterms[0];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_u2 -= bi_bterms[1];
        } else {
          out_u2 -= bi_bterms[1];
        }
      }
    }

    /* Gravity term */
    bi_bterms[0] = 0.;
    bi_bterms[1] = 0.;

    cs_real_t gx_id1 = - cs_math_3_dot_product(gravity, i_face_cog[f_id_sel]);
    cs_real_t gx_id2 = - cs_math_3_dot_product(gravity, i_face_cog[f_id_sel]);

    cs_i_cd_unsteady_upwind(ircflp,
                            diipf[f_id_sel],
                            djjpf[f_id_sel],
                            grad,
                            grad,
                            gx_id1,
                            gx_id2,
                            &pif,
                            &pjf,
                            &pip,
                            &pjp);

    cs_i_conv_flux(iconvp,
                   1.,
                   0, /* Conservative formulation, no mass accumulation */
                   gx_id1,
                   gx_id2,
                   pif,
                   pif, /* no relaxation */
                   pjf,
                   pjf, /* no relaxation */
                   i_mass_flux[f_id_sel],
                   1.,
                   1.,
                   bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_rhogx += bi_bterms[0];
        } else {
          in_rhogx += bi_bterms[0];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_rhogx -= bi_bterms[1];
        } else {
          out_rhogx -= bi_bterms[1];
        }
      }
    }

  }

  /* Free memory */

  BFT_FREE(cells_tag_ids);
  BFT_FREE(bi_face_cells);
  BFT_FREE(i_face_sel_ids);
  BFT_FREE(bb_face_sel_ids);
  BFT_FREE(bi_face_sel_ids);

  /* Sum of values on all ranks (parallel calculations) */

  balance[CS_BALANCE_P_IN] = in_pressure;
  balance[CS_BALANCE_P_OUT] = out_pressure;
  balance[CS_BALANCE_P_U2_IN] = in_u2;
  balance[CS_BALANCE_P_U2_OUT] = out_u2;
  balance[CS_BALANCE_P_RHOGX_IN] = in_rhogx;
  balance[CS_BALANCE_P_RHOGX_OUT] = out_rhogx;
  balance[CS_BALANCE_P_U_IN] = in_debit;
  balance[CS_BALANCE_P_U_OUT] = out_debit;
  balance[CS_BALANCE_P_RHOU_IN] = in_m_debit;
  balance[CS_BALANCE_P_RHOU_OUT] = out_m_debit;

  cs_parall_sum(CS_BALANCE_P_N_TERMS, CS_REAL_TYPE, balance);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 * volumic zone defined by the criterion also given as argument.
 * The different contributions are printed in the run_solver.log.
 *
 * \param[in]     selection_crit      zone selection criterion
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone(const char * selection_crit)
{
  cs_real_t balance[CS_BALANCE_P_N_TERMS];

  const cs_mesh_t *m = cs_glob_mesh;
  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Select cells */

  cs_lnum_t n_cells_sel = 0;
  cs_lnum_t *cells_sel_ids = NULL;

  BFT_MALLOC(cells_sel_ids, m->n_cells, cs_lnum_t);
  cs_selector_get_cell_list(selection_crit, &n_cells_sel, cells_sel_ids);

  /* Compute pressure drop terms */

  cs_pressure_drop_by_zone_compute(n_cells_sel,
                                   cells_sel_ids,
                                   balance);

  BFT_FREE(cells_sel_ids);

  /* Log results at time step n */

  bft_printf(_("   ** PRESSURE DROP BY ZONE at iteration %6i\n"
               "   ---------------------------------------------\n"
               "------------------------------------------------------------\n"
               "   ZONE SELECTION CRITERIA: \"%s\"\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | p u . dS        | p u . dS\n"
               "  |   -    -        |   -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | u^2/2 rho u . dS| u^2/2 rho u . dS\n"
               "  | -         -    -| -         -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  |-rho(g . x)u . dS|-rho(g . x)u . dS\n"
               "  |     -   - -    -|     -   - -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | u . dS          | u . dS\n"
               "  | -    -          | -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | rho u . dS      | rho u . dS\n"
               "  |     -    -      |     -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n\n"),
             nt_cur, selection_crit,
             balance[CS_BALANCE_P_IN], balance[CS_BALANCE_P_OUT],
             balance[CS_BALANCE_P_U2_IN], balance[CS_BALANCE_P_U2_OUT],
             balance[CS_BALANCE_P_RHOGX_IN], balance[CS_BALANCE_P_RHOGX_OUT],
             balance[CS_BALANCE_P_U_IN], balance[CS_BALANCE_P_U_OUT],
             balance[CS_BALANCE_P_RHOU_IN], balance[CS_BALANCE_P_RHOU_OUT]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface balance of a given scalar.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 * \param[in]     normal              outwards normal direction
 */
/*----------------------------------------------------------------------------*/

void
cs_surface_balance(const char       *selection_crit,
                   const char       *scalar_name,
                   const cs_real_t   normal[3])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;

  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Faces selection */

  cs_lnum_t n_b_faces_sel = 0;
  cs_lnum_t *b_face_sel_ids = NULL;
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = NULL;

  BFT_MALLOC(i_face_sel_ids, m->n_i_faces, cs_lnum_t);
  BFT_MALLOC(b_face_sel_ids, m->n_b_faces, cs_lnum_t);

  cs_selector_get_i_face_list(selection_crit, &n_i_faces_sel, i_face_sel_ids);
  cs_selector_get_b_face_list(selection_crit, &n_b_faces_sel, b_face_sel_ids);

  /* Balance on selected faces */

  cs_real_t  balance[CS_BALANCE_N_TERMS];

  cs_flux_through_surface(scalar_name,
                          normal,
                          n_b_faces_sel,
                          n_i_faces_sel,
                          b_face_sel_ids,
                          i_face_sel_ids,
                          balance,
                          NULL,   /* flux_b_faces */
                          NULL);  /* flux_i_faces */

  /* Recount selected interior faces (parallel test) */

  cs_gnum_t n_sel[2] = {(cs_gnum_t)n_b_faces_sel, 0};

  for (cs_lnum_t i = 0; i < n_i_faces_sel; i++) {
    cs_lnum_t f_id = i_face_sel_ids[i];
    if (i_face_cells[f_id][0] < n_cells)
      n_sel[1] += 1;
  }

  cs_parall_sum(2, CS_GNUM_TYPE, n_sel);

  /* Free memory */

  BFT_FREE(i_face_sel_ids);
  BFT_FREE(b_face_sel_ids);

  /* Compute some sums */

  cs_real_t flux_b_faces
    = balance[CS_BALANCE_BOUNDARY_IN] + balance[CS_BALANCE_BOUNDARY_OUT]
    + balance[CS_BALANCE_BOUNDARY_SYM] + balance[CS_BALANCE_BOUNDARY_WALL]
    + balance[CS_BALANCE_BOUNDARY_COUPLED_E]
    + balance[CS_BALANCE_BOUNDARY_OTHER];

  cs_real_t flux_i_faces
    = balance[CS_BALANCE_INTERIOR_IN] + balance[CS_BALANCE_INTERIOR_OUT];

  /* Log balance */

  bft_printf
    (_("\n   ** SURFACE BALANCE at iteration %6i\n"
       "     ------------------------------------\n"
       "------------------------------------------------------------\n"
       "   SCALAR: %s\n"
       "   ZONE SELECTION CRITERIA: \"%s\"\n"
       "   OUTGOING NORMAL: [%.2e, %.2e, %.2e] \n"
       "------------------------------------------------------------\n"
       "   Interior faces selected: %llu of %llu \n"
       "   Boundary faces selected: %llu of %llu \n"
       "------------------------------------------------------------\n"
       "    Boundary faces:        %12.4e \n"
       "    Int. Coupling faces:   %12.4e \n"
       "    Interior faces:        \n"
       "      In:                  %12.4e \n"
       "      Out:                 %12.4e \n"
       "      Balance:             %12.4e \n"
       "------------------------------------------------------------\n"),
     nt_cur, scalar_name, selection_crit,
     normal[0], normal[1], normal[2],
     (unsigned long long)n_sel[1], (unsigned long long)(m->n_g_i_faces),
     (unsigned long long)n_sel[0], (unsigned long long)(m->n_g_b_faces),
     flux_b_faces, balance[CS_BALANCE_BOUNDARY_COUPLED_E],
     balance[CS_BALANCE_INTERIOR_IN], balance[CS_BALANCE_INTERIOR_OUT],
     flux_i_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the face by face surface flux of a given scalar, through a
 *        surface area defined by the given face ids.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]   scalar_name       scalar name
 * \param[in]   normal            outwards normal direction
 * \param[in]   n_b_faces_sel     number of selected boundary faces
 * \param[in]   n_i_faces_sel     number of selected internal faces
 * \param[in]   b_face_sel_ids    ids of selected boundary faces
 * \param[in]   i_face_sel_ids    ids of selected internal faces
 * \param[out]  balance           optional array of computed balance terms
 *                                (see \ref cs_balance_term_t), of
 *                                size CS_BALANCE_N_TERMS, or NULL
 * \param[out]  flux_b_faces      optional surface flux through selected
 *                                boundary faces, or NULL
 * \param[out]  flux_i_faces      optional surface flux through selected
 *                                interior faces, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_flux_through_surface(const char         *scalar_name,
                        const cs_real_t     normal[3],
                        cs_lnum_t           n_b_faces_sel,
                        cs_lnum_t           n_i_faces_sel,
                        const cs_lnum_t     b_face_sel_ids[],
                        const cs_lnum_t     i_face_sel_ids[],
                        cs_real_t          *balance,
                        cs_real_t          *flux_b_faces,
                        cs_real_t          *flux_i_faces)
{
  int idtvar = cs_glob_time_step_options->idtvar;

  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_real_t *local_min = NULL;
  cs_real_t *local_max = NULL;
  cs_real_t *courant = NULL;

  cs_real_t *cv_limiter = NULL;
  cs_real_t *df_limiter = NULL;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict diipf
    = (const cs_real_3_t *restrict)fvq->diipf;
  const cs_real_3_t *restrict djjpf
    = (const cs_real_3_t *restrict)fvq->djjpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  const cs_field_t *f = cs_field_by_name_try(scalar_name);
  const int field_id = cs_field_id_by_name(scalar_name);

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  const int key_lim_choice = cs_field_key_id("limiter_choice");

  /* initialize output */

  cs_real_t  _balance[CS_BALANCE_N_TERMS];
  for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
    _balance[i] = 0;

  /* all boundary convective fluxes are upwind */
  int icvflb = 0; // TODO handle total energy balance
  int icvflf = 0;

  /* Internal cuplin varibale initialization*/
  cs_real_t *pvar_local = NULL;
  cs_real_t *pvar_distant = NULL;
  const cs_lnum_t *faces_local = NULL;
  cs_lnum_t n_local = 0;
  cs_lnum_t n_distant = 0;
  const cs_lnum_t *faces_distant = NULL;
  cs_internal_coupling_t *cpl = NULL;

 /* Physical properties
    ------------------- */

  /* Temperature indicator.
     Will multiply by CP in order to have energy. */
  bool itemperature = false;
  if (f == cs_thermal_model_field()) {
    if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE)
      itemperature = true;
  }

  /* Specific heat (CP) */
  cs_real_t *cpro_cp = NULL;
  const int icp = cs_field_id_by_name("specific_heat");
  if (itemperature) {
    if (icp != -1)
      cpro_cp = CS_F_(cp)->val;
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_cp[c_id] = cp0;
      }
    }
  }
  else {
    BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cpro_cp[c_id] = 1.;
    }
  }

  /* Boundary condition coefficient for h */
  const cs_real_t *a_F = f->bc_coeffs->a;
  const cs_real_t *b_F = f->bc_coeffs->b;
  const cs_real_t *af_F = f->bc_coeffs->af;
  const cs_real_t *bf_F = f->bc_coeffs->bf;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Face viscosity */
  int imvisf = var_cal_opt.imvisf;
  cs_real_t *i_visc;
  cs_real_t *b_visc;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  cs_real_t *c_visc = NULL;
  BFT_MALLOC(c_visc, n_cells_ext, cs_real_t);
  const int kivisl
    = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
  if (kivisl != -1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] = cs_field_by_id(kivisl)->val[c_id];
  }
  else {
    const double visls0
      = cs_field_get_key_double(f, cs_field_key_id("diffusivity_ref"));
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      c_visc[c_id] = visls0;
    }
  }

  /* Turbulent part */
  cs_real_t *c_visct = cs_field_by_name("turbulent_viscosity")->val;

  if (var_cal_opt.idifft == 1) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] += cpro_cp[c_id] * c_visct[c_id]/turb_schmidt;

  }

  cs_face_viscosity(m, fvq, imvisf, c_visc, i_visc, b_visc);

  /* Internal coupling*/

  if (var_cal_opt.icoupl > 0) {
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    const int coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  const int imrgra = var_cal_opt.imrgra;
  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  /* Limiters */

  int limiter_choice = -1;
  int ischcp = var_cal_opt.ischcv;
  if (field_id != -1) {
    f = cs_field_by_id(field_id);

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      limiter_choice = cs_field_get_key_int(f, key_lim_choice);
      BFT_MALLOC(local_max, n_cells_ext, cs_real_t);
      BFT_MALLOC(local_min, n_cells_ext, cs_real_t);
      cs_field_local_extrema_scalar(field_id,
                                    halo_type,
                                    local_max,
                                    local_min);
      if (limiter_choice >= CS_NVD_VOF_HRIC) {
        BFT_MALLOC(courant, n_cells_ext, cs_real_t);
        cs_cell_courant_number(field_id, courant);
      }
    }

    int cv_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("convection_limiter_id"));
    if (cv_limiter_id > -1)
      cv_limiter = cs_field_by_id(cv_limiter_id)->val;

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;
  }

  /* Gradient calculation
     -------------------- */

  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_field_gradient_scalar(f,
                           true, /* use_previous_t */
                           1, /* inc */
                           grad);

  /* Compute the gradient for convective scheme
     (the slope test, limiter, SOLU, etc) */
  cs_real_3_t *gradup = NULL;
  cs_real_3_t *gradst = NULL;
  if (var_cal_opt.blencv > 0 && var_cal_opt.isstpc == 0) {
    BFT_MALLOC(gradst, n_cells_ext, cs_real_3_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      gradst[c_id][0] = 0.;
      gradst[c_id][1] = 0.;
      gradst[c_id][2] = 0.;
    }
    /* Slope test gradient */
    if (var_cal_opt.iconv > 0)
      cs_slope_test_gradient(f->id,
                             1, /* inc */
                             CS_HALO_STANDARD,
                             (const cs_real_3_t *)grad,
                             gradst,
                             f->val,
                             a_F,
                             b_F,
                             i_mass_flux);
  }

  /* Pure SOLU scheme without using gradient_slope_test function
     or Roe and Sweby limiters */
  if (var_cal_opt.blencv > 0
      && (var_cal_opt.ischcv==2 || var_cal_opt.ischcv==4)) {
    BFT_MALLOC(gradup, n_cells_ext, cs_real_3_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      gradup[c_id][0] = 0.;
      gradup[c_id][1] = 0.;
      gradup[c_id][2] = 0.;
    }

    if (var_cal_opt.iconv > 0)
      cs_upwind_gradient(f->id,
                         1, /* inc */
                         CS_HALO_STANDARD,
                         a_F,
                         b_F,
                         i_mass_flux,
                         b_mass_flux,
                         f->val,
                         gradup);
  }

  /* Faces selection
     --------------- */

  cs_lnum_2_t *bi_face_cells = NULL;

  if (n_i_faces_sel > 0) {

    BFT_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      bi_face_cells[f_id][0] = -999;
      bi_face_cells[f_id][1] = -999;
    }

    for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {
      cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
      cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
      cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

      cs_real_t dot_pro = cs_math_3_dot_product(normal, i_face_normal[f_id_sel]);
      if (fabs(dot_pro) < 1.0e-14)//FIXME
        dot_pro = 0;
      if(dot_pro > 0.)
        bi_face_cells[f_id_sel][0] = c_id1;
      else if (dot_pro < 0.)
        bi_face_cells[f_id_sel][1] = c_id2;
    }

    if (flux_i_faces != NULL) {
      for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++)
        flux_i_faces[f_id] = 0.;
    }

    if (flux_b_faces != NULL) {
      for (cs_lnum_t f_id = 0; f_id < n_b_faces_sel; f_id++)
        flux_b_faces[f_id] = 0.;
    }

  }

  /* Boundary faces contribution
     --------------------------- */

  int iconvp = var_cal_opt.iconv;
  int idiffp = var_cal_opt.idiff;
  int ircflp = var_cal_opt.ircflu;
  double relaxp = var_cal_opt.relaxv;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = (b_face_sel_ids != NULL) ? b_face_sel_ids[f_id] : f_id;

    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    cs_real_t term_balance = 0.;

    cs_real_t ac_F = 0.;
    cs_real_t bc_F = 0.;

    if (icvflb == 1) {
      ac_F = f->bc_coeffs->ac[f_id_sel];
      bc_F = f->bc_coeffs->bc[f_id_sel];
      icvflf = 0; /* = icvfli */
    }

    _balance_boundary_faces(icvflf,
                            idtvar,
                            iconvp,
                            idiffp,
                            ircflp,
                            relaxp,
                            diipb[f_id_sel],
                            grad[c_id],
                            f->val[c_id],
                            f->val_pre[c_id],
                            bc_type[f_id_sel],
                            b_visc[f_id_sel],
                            a_F[f_id_sel],
                            b_F[f_id_sel],
                            af_F[f_id_sel],
                            bf_F[f_id_sel],
                            ac_F,
                            bc_F,
                            b_mass_flux[f_id_sel],
                            cpro_cp[c_id],
                            &term_balance);

    if (flux_b_faces != NULL)
      flux_b_faces[f_id] = term_balance;

    if (bc_type[f_id_sel] == CS_INLET ||
        bc_type[f_id_sel] == CS_FREE_INLET ||
        bc_type[f_id_sel] == CS_ESICF ||
        bc_type[f_id_sel] == CS_EPHCF)
      _balance[CS_BALANCE_BOUNDARY_IN] -= term_balance;
    else if (bc_type[f_id_sel] == CS_OUTLET ||
             bc_type[f_id_sel] == CS_SSPCF ||
             bc_type[f_id_sel] == CS_SOPCF)
      _balance[CS_BALANCE_BOUNDARY_OUT] -= term_balance;
    else if (bc_type[f_id_sel] == CS_SYMMETRY)
      _balance[CS_BALANCE_BOUNDARY_SYM] -= term_balance;
    else if (bc_type[f_id_sel] == CS_SMOOTHWALL)
      _balance[CS_BALANCE_BOUNDARY_WALL_S] -= term_balance;
    else if (bc_type[f_id_sel] == CS_ROUGHWALL)
      _balance[CS_BALANCE_BOUNDARY_WALL_R] -= term_balance;
    else if (   bc_type[f_id_sel] == CS_COUPLED
             || bc_type[f_id_sel] == CS_COUPLED_FD)
      _balance[CS_BALANCE_BOUNDARY_COUPLED_E] -= term_balance;
    else
      _balance[CS_BALANCE_BOUNDARY_OTHER] -= term_balance;

  }

  /* Balance on coupled faces
     ------------------------

    We handle different types of boundary faces separately to better
    analyze the information, but this is not mandatory. */

  if (var_cal_opt.icoupl > 0) {

    cs_lnum_t *inv_b_face_sel_ids = NULL;

    BFT_MALLOC(inv_b_face_sel_ids, n_b_faces, cs_lnum_t);
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
      inv_b_face_sel_ids[f_id] = -1;

    if (b_face_sel_ids != NULL) {
      for (cs_lnum_t f_id = 0; f_id < n_b_faces_sel; f_id++) {
        cs_lnum_t f_id_sel = b_face_sel_ids[f_id];
        inv_b_face_sel_ids[f_id_sel] = f_id;
      }
    }
    else {
      for (cs_lnum_t f_id_sel = 0; f_id_sel < n_b_faces_sel; f_id_sel++)
        inv_b_face_sel_ids[f_id_sel] = f_id_sel;
    }

    /* Prepare data for sending from distant */
    BFT_MALLOC(pvar_distant, n_distant, cs_real_t);

    for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
      cs_lnum_t f_id = faces_distant[ii];
      cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t pip;
      cs_b_cd_unsteady(ircflp,
                       diipb[f_id],
                       grad[c_id],
                       f->val[c_id],
                       &pip);
      pvar_distant[ii] = pip;
    }

    /* Receive data */
    BFT_MALLOC(pvar_local, n_local, cs_real_t);
    cs_internal_coupling_exchange_var(cpl,
                                      1, /* Dimension */
                                      pvar_distant,
                                      pvar_local);

    /* Flux contribution */

    for (cs_lnum_t ii = 0; ii < n_local; ii++) {

      cs_lnum_t f_id = faces_local[ii];
      cs_lnum_t sel_f_id = inv_b_face_sel_ids[f_id];

      if (sel_f_id < 0)
        continue;

      cs_lnum_t c_id = b_face_cells[f_id];
      cs_real_t surf = b_face_surf[f_id];

      cs_real_t pip, pjp;
      cs_real_t term_balance = 0.;

      cs_b_cd_unsteady(ircflp,
                       diipb[f_id],
                       grad[c_id],
                       f->val[c_id],
                       &pip);

      pjp = pvar_local[ii];

      cs_real_t hint = f->bc_coeffs->hint[f_id];
      cs_real_t rcodcl2 = f->bc_coeffs->rcodcl2[f_id];
      cs_real_t heq = 0;
      if (fabs(hint + rcodcl2) > 0)
        heq = surf * hint * rcodcl2 / (hint + rcodcl2);

      cs_b_diff_flux_coupling(idiffp,
                              pip,
                              pjp,
                              heq,
                              &term_balance);

      if (flux_b_faces != NULL)
        flux_b_faces[inv_b_face_sel_ids[f_id]] = term_balance;

      _balance[CS_BALANCE_BOUNDARY_COUPLED_I] -= term_balance;

    }

    BFT_FREE(pvar_local);
    BFT_FREE(pvar_distant);

    BFT_FREE(inv_b_face_sel_ids);
  }

  /* Balance on selected interior faces */

  int isstpp = var_cal_opt.isstpc;
  double blencp = var_cal_opt.blencv;
  double blend_st = var_cal_opt.blend_st;
  int iupwin = (blencp > 0.) ? 0 : 1;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    cs_real_t beta = blencp;
    /* Beta blending coefficient ensuring positivity of the scalar */
    if (isstpp == 2) {
      beta = CS_MAX(CS_MIN(cv_limiter[c_id1], cv_limiter[c_id2]), 0.);
    }

    int bldfrp = ircflp;
    /* Local limitation of the reconstruction */
    if (df_limiter != NULL && ircflp > 0) {
      cs_real_t _bldfrp = fmax(fmin(df_limiter[c_id1],
                                    df_limiter[c_id2]), 0.);
      bldfrp = (int)_bldfrp;
    }

    cs_real_t hybrid_coef_ii, hybrid_coef_jj;
    cs_lnum_t ic = -1, id = -1;
    cs_real_t courant_c = -1., _local_max = 0., _local_min = 0.;
    if (ischcp == 3) {
      hybrid_coef_ii = CS_F_(hybrid_blend)->val[c_id1];
      hybrid_coef_jj = CS_F_(hybrid_blend)->val[c_id2];
    }
    else if (ischcp == 4) {
      hybrid_coef_ii = 0.;
      hybrid_coef_jj = 0.;
      /* Determine central and downwind sides w.r.t. current face */
      cs_central_downwind_cells(c_id1,
                                c_id2,
                                i_mass_flux[f_id_sel],
                                &ic,  /* central cell id */
                                &id); /* downwind cell id */

      if (courant != NULL)
        courant_c = courant[ic];

      if (local_max != NULL) {
        _local_max = local_max[ic];
        _local_min = local_min[ic];
      }
    }
    else {
      hybrid_coef_ii = 0.;
      hybrid_coef_jj = 0.;
    }

    cs_real_2_t bi_bterms = {0., 0.};

    _balance_internal_faces(iupwin,
                            idtvar,
                            iconvp,
                            idiffp,
                            bldfrp,
                            ischcp,
                            isstpp,
                            (cs_nvd_type_t)limiter_choice,
                            relaxp,
                            beta,
                            blend_st,
                            weight[f_id_sel],
                            i_dist[f_id_sel],
                            cell_cen[c_id1],
                            cell_cen[c_id2],
                            cell_cen[ic],
                            cell_cen[id],
                            i_face_normal[f_id_sel],
                            i_face_cog[f_id_sel],
                            hybrid_coef_ii,
                            hybrid_coef_jj,
                            diipf[f_id_sel],
                            djjpf[f_id_sel],
                            grad[c_id1],
                            grad[c_id2],
                            grad[ic],
                            gradup[c_id1],
                            gradup[c_id2],
                            gradst[c_id1],
                            gradst[c_id2],
                            f->val[c_id1],
                            f->val[c_id2],
                            f->val[ic],
                            f->val[id],
                            f->val_pre[c_id1],
                            f->val_pre[c_id2],
                            i_visc[f_id_sel],
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            _local_max,
                            _local_min,
                            courant_c,
                            bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check i_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (flux_i_faces != NULL)
          flux_i_faces[f_id] -= bi_bterms[0];
        if (i_mass_flux[f_id_sel] > 0)
          _balance[CS_BALANCE_INTERIOR_IN] -= bi_bterms[0];
        else
          _balance[CS_BALANCE_INTERIOR_OUT] -= bi_bterms[0];
      }
    }
    /* Face normal direction reversed */
    else if (bi_face_cells[f_id_sel][1] >= 0) {
      if (c_id2 < n_cells) {
        if (flux_i_faces != NULL)
          flux_i_faces[f_id] += bi_bterms[1];
        if (i_mass_flux[f_id_sel] > 0)
          _balance[CS_BALANCE_INTERIOR_IN] += bi_bterms[1];
        else
          _balance[CS_BALANCE_INTERIOR_OUT] += bi_bterms[1];
      }
    }
  }

  if (balance != NULL) {

    _balance[CS_BALANCE_BOUNDARY_WALL] =   _balance[CS_BALANCE_BOUNDARY_WALL_S]
                                         + _balance[CS_BALANCE_BOUNDARY_WALL_R];
    _balance[CS_BALANCE_BOUNDARY_COUPLED]
      =   _balance[CS_BALANCE_BOUNDARY_COUPLED_E]
        + _balance[CS_BALANCE_BOUNDARY_COUPLED_I];

    for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
      balance[i] = _balance[i];
    cs_parall_sum(CS_BALANCE_N_TERMS, CS_REAL_TYPE, balance);
  }

  /* Free memory */

  BFT_FREE(bi_face_cells);
  BFT_FREE(grad);
  BFT_FREE(gradup);
  BFT_FREE(gradst);
  BFT_FREE(local_max);
  BFT_FREE(local_min);
  BFT_FREE(courant);

  if (!itemperature || icp == -1)
    BFT_FREE(cpro_cp);
  BFT_FREE(c_visc);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
