/*============================================================================
 * Low-level functions and structures used to build the algebraic system with
 * a cellwise process when Hybrid High Order schemes are set for the space
 * discretization
 *============================================================================*/

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
#include <float.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_cdo_bc.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_property.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_macfb_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_MAC_BUILDER_DBG 0

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointer of pointers to global structures */

cs_macfb_builder_t **cs_mac_builders = nullptr;

/*============================================================================
 * Local static variables
 *============================================================================*/

static int cs_macfb_builder_n_structures = 0;

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate global structures used for MAC builder
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_initialize(void)
{
  assert(cs_glob_n_threads > 0);

  int nthr = cs_glob_n_threads;

  cs_macfb_builder_n_structures = nthr;
  BFT_MALLOC(cs_mac_builders, nthr, cs_macfb_builder_t *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_mac_builders[t_id] = cs_macfb_builder_create();
  }
#else

  assert(cs_glob_n_threads == 1);

  cs_mac_builders[0] = cs_macfb_builder_create();

#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to \ref cs_macfb_builder_t
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_finalize(void)
{
  if (cs_macfb_builder_n_structures < 1)
    return;

  assert(cs_macfb_builder_n_structures == cs_glob_n_threads);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_macfb_builder_free(&(cs_mac_builders[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_macfb_builder_free(&(cs_mac_builders[0]));
#endif /* openMP */

  BFT_FREE(cs_mac_builders);
  cs_macfb_builder_n_structures = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate  a cs_macfb_builder_t structure
 *
 * \return a pointer to a new allocated cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_builder_t *
cs_macfb_builder_create(void)
{
  cs_macfb_builder_t *b = nullptr;

  BFT_MALLOC(b, 1, cs_macfb_builder_t);

  cs_macfb_builder_reset(b);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_macfb_builder_t structure
 *
 * \param[in, out] p_builder  pointer of pointer on a cs_macfb_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_free(cs_macfb_builder_t **p_builder)
{
  if (p_builder == nullptr)
    return;

  cs_macfb_builder_t *b = *p_builder;

  BFT_FREE(b);

  *p_builder = nullptr;
  p_builder  = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize to invalid values a cs_macfb_builder_t structure
 *
 * \param[in, out]  macb       pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_reset(cs_macfb_builder_t *macb)
{
  macb->n_fc   = 0;
  macb->n_dofs = 0;

  /* Face information */

  for (short int f = 0; f < 6; f++) {
    macb->f_axis[f]     = CS_FLAG_N_AXIS;
    macb->f_sgn_axis[f] = 0;
    macb->f_vol_cv[f]   = -DBL_MAX;
    macb->f_h_cv[f]     = -DBL_MAX;
    macb->f_opp_idx[f]  = -1;
  }

  for (short int f = 0; f < 24; f++) {
    macb->f2f_idx[f]       = -1;
    macb->f2f_ids[f]       = -1;
    macb->f2f_h[f]         = -DBL_MAX;
    macb->f2f_surf_cv_c[f] = -DBL_MAX;
    macb->dir_values[f]    = -DBL_MAX;
  }

  for (short int f = 0; f < 30; f++) {
    macb->f_ids[f]   = -1;
    macb->dof_ids[f] = -1;
  }

  memset(macb->f2fo_idx, -1, 48 * sizeof(short int));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_macfb_builder_t structure corresponding to id
 *
 * \param[in]   id   id in the array of pointer to cs_macfb_builder_t struct.
 *
 * \return a pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_builder_t *
cs_macfb_get_builder(int id)
{
  if (id < 0 || id >= cs_glob_n_threads)
    return nullptr;

  cs_macfb_builder_t *macb = cs_mac_builders[id];

  return macb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set-up face informations needed to build operators
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  macb       pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_cellwise_setup(const cs_cell_mesh_t      *cm,
                                const cs_cdo_connect_t    *connect,
                                const cs_cdo_quantities_t *quant,
                                cs_macfb_builder_t        *macb)
{
  if (macb == nullptr)
    return;

  /* Sanity checks */
  assert(cm != nullptr && quant != nullptr && connect != nullptr);

  cs_macfb_builder_reset(macb);

  /* cell id */
  macb->c_id = cm->c_id;

  /* number of face - updated later */
  macb->n_fc = 6;

  /* number of dofs - updated later */
  macb->n_dofs = 6;

  /* maximum number of dofs - 6 faces + 4 faces by face = 6 + 4*4 =30*/
  macb->n_max_dofs = 30;

  const cs_adjacency_t *f2f_ed = connect->f2f_ed;
  assert(f2f_ed != nullptr); /* Has to be build before */

  cs_adjacency_t *f2e = connect->f2e;
  assert(f2e != nullptr); /* Has to be build before */

  /* loop on internal faces */
  for (short int f = 0; f < 6; f++) {
    const cs_lnum_t f_id = cm->f_ids[f];

    macb->f_ids[f]  = f_id;
    macb->f_axis[f] = quant->face_axis[f_id];

    const cs_real_t *f_c = cs_quant_get_face_center(f_id, quant);

    if (cm->f_sgn[f]
          * cs_quant_get_face_vector_area(f_id, quant)[macb->f_axis[f]]
        > 0) {
      macb->f_sgn_axis[f] = 1;
    }
    else {
      macb->f_sgn_axis[f] = -1;
    }

    /* add dof id*/
    macb->dof_ids[f] = f_id;

    /* normal length of the control volume for a face */
    macb->f_h_cv[f] = (f_id < cm->bface_shift)
                        ? quant->i_dist[f_id]
                        : quant->b_dist[f_id - cm->bface_shift];

    /* total volume of the control volume */
    macb->f_vol_cv[f] = macb->f_h_cv[f] * cm->face[f].meas;

    /* volume of the control volume restricted to the cell */
    const cs_real_t f_vol_cv_c = cm->hfc[f] * cm->face[f].meas;

    /* number of connected faces to the current face throught edge and
     * same direction */
    const short int nb_f_outer = f2f_ed->idx[f_id + 1] - f2f_ed->idx[f_id];

    /* Physical face */
    for (short int fj = 0; fj < nb_f_outer; fj++) {
      const short int shift_j = 4 * f + fj;

      const cs_lnum_t  fj_id = f2f_ed->ids[f2f_ed->idx[f_id] + fj];
      const cs_real_t *fj_c  = cs_quant_get_face_center(fj_id, quant);

      macb->f_ids[macb->n_fc]     = fj_id;
      macb->f_axis[macb->n_fc]    = quant->face_axis[fj_id];
      macb->dof_ids[macb->n_dofs] = fj_id;

      macb->f2f_idx[shift_j] = macb->n_dofs;
      macb->f2f_ids[shift_j] = fj_id;
      macb->f2f_h[shift_j]   = cs_math_3_distance(f_c, fj_c);

      /* Find common edge */
      for (short int ej = 0; ej < 4; ej++) {
        const cs_lnum_t ej_id = f2e->ids[f2e->idx[fj_id] + ej];
        for (short int ei = 0; ei < 4; ei++) {
          const cs_lnum_t ei_id = f2e->ids[f2e->idx[f_id] + ei];
          if (ei_id == ej_id) {
            macb->f2e_ids[shift_j] = ei_id;
            break;
          }
        }
      }

      const cs_quant_t ei_q
        = cs_quant_get_edge_center(macb->f2e_ids[shift_j], connect, quant);
      const cs_real_t length       = 2.0 * cs_math_3_distance(f_c, ei_q.center);
      macb->f2f_surf_cv_c[shift_j] = f_vol_cv_c / length;

      macb->n_dofs++;
      macb->n_fc++;
    }

    /* Virtual Face - use edge instead */
    /* Find bounday edges */
    if (nb_f_outer < 4) {
      /* Boundary edges */
      assert(f2e->idx[f_id + 1] - f2e->idx[f_id] == 4);

      cs_lnum_t fe_id[4] = { f2e->ids[f2e->idx[f_id] + 0],
                             f2e->ids[f2e->idx[f_id] + 1],
                             f2e->ids[f2e->idx[f_id] + 2],
                             f2e->ids[f2e->idx[f_id] + 3] };

      /* Loop on neighbourhood faces */
      for (short int fj = 0; fj < nb_f_outer; fj++) {
        const cs_lnum_t fj_id = f2f_ed->ids[f2f_ed->idx[f_id] + fj];

        /* Loop on neighbourhood faces' edges */
        for (short int ej = 0; ej < 4; ej++) {
          const cs_lnum_t ej_id = f2e->ids[f2e->idx[fj_id] + ej];

          /* Loop on faces' edges */
          for (short int ei = 0; ei < 4; ei++) {
            if (ej_id == fe_id[ei]) {
              /* Same id -> not a boundary edge */
              fe_id[ei] = -1;
              break;
            }
          }
        }
      }

      short int n_ed_find = 0;
      for (short int ei = 0; ei < 4; ei++) {
        if (fe_id[ei] >= 0) {
          n_ed_find++;
          const cs_lnum_t ei_id = fe_id[ei];

          const short int shift_j = 4 * f + nb_f_outer - 1 + n_ed_find;

          const cs_quant_t ei_q
            = cs_quant_get_edge_center(ei_id, connect, quant);

          macb->f2f_idx[shift_j] = -1;
          macb->f2f_ids[shift_j] = -(ei_id + 1);
          macb->f2e_ids[shift_j] = ei_id;

          macb->f2f_h[shift_j] = cs_math_3_distance(f_c, ei_q.center);
          macb->f2f_surf_cv_c[shift_j]
            = f_vol_cv_c / (2.0 * macb->f2f_h[shift_j]);
        }
      }
      assert(n_ed_find == 4 - nb_f_outer);
    }
  }

  assert(macb->n_dofs <= macb->n_max_dofs);
  assert(macb->n_dofs == macb->n_fc);

  /* loop on internal faces to find orthogonal faces */
  for (short int f = 0; f < 6; f++) {

    /* Find opposite face */
    for (short int fj = 0; fj < 6; fj++) {
      if (f != fj && macb->f_axis[f] == macb->f_axis[fj]) {
        macb->f_opp_idx[f] = fj;
        break;
      }
    }

    /* Loop on outer faces */
    for (short int fj = 0; fj < 4; fj++) {
      const short int shift_j = 4 * f + fj;

      const cs_lnum_t ej_id = macb->f2e_ids[shift_j];

      short int n_e_find = 0;
      for (short int fk = 0; fk < macb->n_fc; fk++) {
        /* Orthogonal direction */
        if (macb->f_axis[f] != macb->f_axis[fk]) {
          const cs_lnum_t fk_id = macb->f_ids[fk];

          /* Boundary edges */
          assert(f2e->idx[fk_id + 1] - f2e->idx[fk_id] == 4);

          for (short int ek = 0; ek < 4; ek++) {
            cs_lnum_t ek_id = f2e->ids[f2e->idx[fk_id] + ek];
            if (ej_id == ek_id) {
              macb->f2fo_idx[2 * shift_j + n_e_find] = fk;
              n_e_find++;
              break;
            }
          }
        }
        if (n_e_find == 2) {
          break;
        }
      }
      assert(n_e_find > 0);
    }
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_MAC_BUILDER_DBG > 2
  if (macb->c_id == 0)
    cs_macfb_builder_dump(macb);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_macfb_builder_t structure
 *
 * \param[in]    macb    pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_builder_dump(const cs_macfb_builder_t *macb)
{

  if (macb == nullptr) {
    bft_printf("\n>> Dump cs_macfb_builder_t %p\n", (const void *)macb);
    return;
  }

  bft_printf("\n>> Dump cs_macfb_builder_t of cell id: %d\n", macb->c_id);

  /* Information related to primal faces */

  bft_printf(" %s | %6s | %5s | %5s | %8s  | %8s | %8s \n",
             "f",
             "id",
             "dof id",
             "axis",
             "s_axis",
             "vol_cv",
             "h_cv");

  for (short int f = 0; f < 6; f++) {
    bft_printf("%2d | %6ld | %6ld | %4d  | %5d  | %.3e | %.3e \n",
               f,
               (long)macb->f_ids[f],
               (long)macb->dof_ids[f],
               macb->f_axis[f],
               macb->f_sgn_axis[f],
               macb->f_vol_cv[f],
               macb->f_h_cv[f]);
  }

  for (short int f = 0; f < 6; f++) {

    bft_printf("face %d with id %d : \n", f, macb->f_ids[f]);
    bft_printf(" %s | %6s | %6s | %6s | %6s    | %8s \n",
               "f",
               "f id",
               "dof id",
               "e id",
               "h",
               "surf_cv");

    for (short int fc = 0; fc < 4; fc++) {
      bft_printf("%2d | %6ld | %6ld | %6ld | %.3e | %.3e \n",
                 fc,
                 (long)macb->f2f_ids[4 * f + fc],
                 macb->f2f_idx[4 * f + fc] >= 0
                   ? (long)macb->dof_ids[macb->f2f_idx[4 * f + fc]]
                   : macb->f2f_idx[4 * f + fc],
                 (long)macb->f2e_ids[4 * f + fc],
                 macb->f2f_h[4 * f + fc],
                 macb->f2f_surf_cv_c[4 * f + fc]);
    }
  }
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
