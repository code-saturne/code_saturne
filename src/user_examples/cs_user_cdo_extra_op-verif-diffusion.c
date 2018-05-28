/*============================================================================
 * Additional post-processing functions defined by user related to CDO schemes
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

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_array_reduce.h"
#include "cs_blas.h"
#include "cs_domain.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_field.h"
#include "cs_sdm.h"
#include "cs_cdofb_scaleq.h"
#include "cs_equation.h"
#include "cs_equation_param.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_cdo_extra_op-verif-diffusion.c
 *
 * \brief Additional user-defined post-processing and analysis functions
*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

static FILE  *resume = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for an analytic function
 *         elt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *
 * \param[in]      time     when ?
 * \param[in]      n_elts   number of elements to consider
 * \param[in]      elt_ids  list of elements ids (to access coords and fill)
 * \param[in]      coords   where ?
 * \param[in]      compact  true:no indirection, false:indirection for filling
 * \param[in]      input    pointer to a structure cast on-the-fly (may be NULL)
 * \param[in, out] retval   result of the function
 */
/*----------------------------------------------------------------------------*/

static inline void
_get_sol(cs_real_t          time,
         cs_lnum_t          n_pts,
         const cs_lnum_t    pt_ids[],
         const cs_real_t   *xyz,
         bool               compact,
         void              *input,
         cs_real_t         *retval)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const double  pi = 4.0*atan(1.0);
  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {

      const cs_lnum_t  id = pt_ids[p];
      const cs_real_t  *_xyz = xyz + 3*id;
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      retval[id] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));

    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_real_t  *_xyz = xyz + 3*pt_ids[p];
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      retval[p] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));
    }

  }
  else {

    assert(pt_ids == NULL);
    for (cs_lnum_t p = 0; p < n_pts; p++) {
      const cs_real_t  *_xyz = xyz + 3*p;
      const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

      retval[p] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));
    }

  }

}

static cs_analytic_func_t *get_sol = _get_sol;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a time step structure
 * \param[in]  eq         pointer to a cs_equation_t structure
 * \param[in]  anacomp    do an analytic comparison or not
 */
/*----------------------------------------------------------------------------*/

static void
_cdovb_post(const cs_cdo_connect_t     *connect,
            const cs_cdo_quantities_t  *cdoq,
            const cs_time_step_t       *time_step,
            const cs_equation_t        *eq,
            bool                        anacomp)
{
  int  len;

  char  *postlabel = NULL;
  double  *ddip = NULL, *rpex = NULL;

  const double  tcur = time_step->t_cur;
  const cs_lnum_t  n_vertices = cdoq->n_vertices;
  const cs_field_t  *field = cs_equation_get_field(eq);
  const cs_real_t  *pdi = field->val;
  const cs_adjacency_t  *c2v = connect->c2v;

  /* Analyze the discrete solution */
  cs_real_t  pdi_min, pdi_max, pdi_wsum, pdi_asum, pdi_ssum;
  cs_array_reduce_minmax_l(cdoq->n_vertices, 1, NULL, pdi,
                           &pdi_min, &pdi_max);

  cs_array_scatter_reduce_norms_l(cdoq->n_cells, c2v->idx, c2v->ids,
                                  NULL, // filter list
                                  1,    // dim
                                  cdoq->n_vertices,
                                  pdi,
                                  cdoq->dcell_vol,
                                  &pdi_wsum,
                                  &pdi_asum,
                                  &pdi_ssum);
  /* Parallel treatment */
  if (cs_glob_n_ranks > 1) {

    cs_real_t  minmax[2] = {-pdi_min, pdi_max};
    cs_parall_max(2, CS_REAL_TYPE, minmax);
    pdi_min = -minmax[0], pdi_max = minmax[1];

    cs_real_t  sums[3] = {pdi_wsum, pdi_asum, pdi_ssum};
    cs_parall_sum(3, CS_REAL_TYPE, sums);
    pdi_wsum = sums[0], pdi_asum = sums[1], pdi_ssum = sums[2];

  }

  if (cs_glob_rank_id < 1) { /* Only the first rank write something */

    const double  inv = 1/cdoq->vol_tot;
    const cs_real_t  pdi_mean = pdi_wsum*inv;

    fprintf(resume, " -bnd- Scal.Min   % 10.6e\n", pdi_min);
    fprintf(resume, " -bnd- Scal.Max   % 10.6e\n", pdi_max);
    fprintf(resume, " -bnd- Scal.Mean  % 10.6e\n", pdi_mean);
    fprintf(resume, " -bnd- Scal.Varia % 10.6e\n",
            pdi_ssum*inv - pdi_mean*pdi_mean);
    fprintf(resume, "%s", msepline);

  }

  if (anacomp) { /* Comparison with an analytical solution */

    /* pex = exact potential
       pdi = discrete potential (solution of the discrete system)
       rpex = Red(vtxsp)(Pex) = reduction on vertices of the exact potential
       ddip = rpex - pdi
    */

    BFT_MALLOC(rpex, n_vertices, double);
    BFT_MALLOC(ddip, n_vertices, double);
    get_sol(tcur, n_vertices, NULL, cdoq->vtx_coord, true, NULL, rpex);
    for (int i = 0; i < n_vertices; i++)
      ddip[i] = rpex[i] - pdi[i];

    len = strlen(field->name) + 7 + 1;
    BFT_MALLOC(postlabel, len, char);
    sprintf(postlabel, "%s.Error", field->name);

    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_ALL_ASSOCIATED,
                             postlabel,
                             1,           /* dim */
                             false,       /* interlace */
                             true,        /* parent mesh */
                             CS_POST_TYPE_cs_real_t,
                             ddip,        /* values on vertices */
                             time_step);  /* time step structure */

    sprintf(postlabel, "%s.RefSol", field->name);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_ALL_ASSOCIATED,
                             postlabel,
                             1,           /* dim */
                             false,       /* interlace */
                             true,        /* parent mesh */
                             CS_POST_TYPE_cs_real_t,
                             rpex,        /* values on vertices */
                             time_step);  /* time step structure */

    /* Free */
    BFT_FREE(postlabel);
    BFT_FREE(ddip);
    BFT_FREE(rpex);

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional operations on results produced by CDO schemes.
 *         Define advanced post-processing and/or analysis for instance.
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(cs_domain_t          *domain)
{
  const cs_cdo_connect_t  *connect = domain->connect;
  const cs_cdo_quantities_t  *cdoq = domain->cdo_quantities;
  const cs_time_step_t  *time_step = domain->time_step;

  cs_equation_t  *eq = cs_equation_by_name("FVCA6.1");
  const char *eqname = cs_equation_get_name(eq);
  const cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid equation name. Stop extra operations.");

  /* Open a file */
  char  *filename = NULL;
  int len = strlen("Resume-.log")+strlen(eqname)+1;

  if (eqp->flag & CS_EQUATION_UNSTEADY) {
    if (time_step->nt_cur == 0)
      return;
    if (time_step->nt_cur % domain->output_nt > 0)
      return;

    len += 9;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Resume-%s-t%.f.log", eqname, time_step->t_cur);
  }
  else {
    if (time_step->nt_cur > 0)
      return;

    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Resume-%s.log", eqname);
  }

  resume = fopen(filename, "w");

  bft_printf("\n%s", msepline);
  bft_printf("    Extra operations\n");
  bft_printf("%s", msepline);

  /* Extra-operation depends on the numerical scheme */
  cs_param_space_scheme_t  space_scheme = cs_equation_get_space_scheme(eq);

  switch (space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    _cdovb_post(connect, cdoq, time_step, eq, true);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid space scheme. Stop post-processing.\n"));
  }

  bft_printf("\n");
  bft_printf(" >> Equation %s (done)\n", eqname);
  printf("\n >> Extra operation for equation: %s\n", eqname);

  /* Free */
  BFT_FREE(filename);
  fclose(resume);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
