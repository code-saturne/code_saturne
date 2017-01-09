/*============================================================================
 * Additional post-processing functions defined by user related to CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_blas.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"
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

/*! \endcond (end ignore by Doxygen) */

static FILE  *resume = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* ---------------------------------------------------------------------------
 * Retrieve the analytical solution for a given problem
 * ---------------------------------------------------------------------------*/

static inline void
_get_sol(cs_real_t          time,
         const cs_real_3_t  xyz,
         cs_get_t          *retval)
{
  CS_UNUSED(time);

  const double  x = xyz[0], y = xyz[1], z = xyz[2];
  const double  pi = 4.0*atan(1.0);

  double  solution = 1+sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));

  (*retval).val = solution;
}

static cs_analytic_func_t *get_sol = _get_sol;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump information into resume file to specify the parameters of the
 *         current simulation
 *
 * \param[in]  cdoq          pointer to a cs_cdo_quantities_t structure
 * \param[in]  space_scheme  scheme for the discretization in space
 * \param[in]  eqp           pointer to the setting structure related to an
 *                           equation
 */
/*----------------------------------------------------------------------------*/

static void
_dump_info(const cs_cdo_quantities_t   *cdoq,
           const cs_space_scheme_t      space_scheme,
           const cs_equation_param_t   *eqp)

{
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

  fprintf(resume," -dim- n_vertices  %d\n", cdoq->n_vertices);
  fprintf(resume," -dim- n_edges     %d\n", cdoq->n_edges);
  fprintf(resume," -dim- n_faces     %d\n", cdoq->n_faces);
  fprintf(resume," -dim- n_cells     %d\n", cdoq->n_cells);
  fprintf(resume, "\n%s", msepline);

  if (eqp->flag & CS_EQUATION_DIFFUSION) {

    if (space_scheme == CS_SPACE_SCHEME_CDOFB) {
      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
        fprintf(resume, " -hdg- Hodge.Op         EDFP:COST\n");
        break;
      case CS_PARAM_HODGE_ALGO_VORONOI:
        fprintf(resume, " -hdg- Hodge.Op         EDFP:VORONOI\n");
        break;
      default:
        fprintf(resume, " -hdg- Hodge.Op         EDFP:Unknown\n");
        break;
      } // end of switch
    }
    else if (space_scheme == CS_SPACE_SCHEME_CDOVB) {
      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
        fprintf(resume, " -hdg- Hodge.Op         EPFD:COST\n");
        break;
      case CS_PARAM_HODGE_ALGO_VORONOI:
        fprintf(resume, " -hdg- Hodge.Op         EPFD:VORONOI\n");
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        fprintf(resume, " -hdg- Hodge.Op         EPFD:WHITNEY_BARY\n");
        break;
      default:
        fprintf(resume, " -hdg- Hodge.Op         EPFD:Unkown\n");
        break;

      } // end of switch
    }

    if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
      fprintf(resume, " -hdg- Beta.Coef        %5.3e\n", h_info.coef);

  } // Diffusion term is activated

  fprintf(resume," -bc-  Enforcement    %s",
          cs_param_get_bc_enforcement_name(eqp->bc->enforcement));

  fprintf(resume, "\n%s", msepline);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the L2 norm for a potential in vertex-based scheme from
 *         a conforming reconstruction
 *
 * \param[in]   topo       pointer to the connectivity struct.
 * \param[in]   geom       pointer to the additional quantities struct.
 * \param[in]   time_step  pointer to a time step structure
 * \param[in]   pdi        pointer to the field of vtx-based DoFs
 * \param[out]  num        int_Omega (pex - Reco{vtxsp}(pdi))^2
 * \param[out]  denum      int_Omega pex^2
 */
/*----------------------------------------------------------------------------*/

static void
_compute_vb_l2pot(const cs_mesh_t             *m,
                  const cs_cdo_connect_t      *topo,
                  const cs_cdo_quantities_t   *geom,
                  const cs_time_step_t        *time_step,
                  const double                *pdi,
                  double                      *num,
                  double                      *denum)
{
  int  i, j, k, p, c_id;
  double  ddip_gpts, n_add, d_add;
  cs_real_3_t  xc, gpts[5];
  double  rpex_gpts[5], weights[5], pdi_gpts[5];

  double  _num = 0.0, _denum = 0.0;
  double  *pdi_recc = NULL, *pdi_recf = NULL;

  const double  tcur = time_step->t_cur;
  const double  *xyz = m->vtx_coord;

  /* Reconstruct potentials at face centers and cell centers */
  cs_reco_conf_vtx_dofs(topo, geom, pdi, &pdi_recc, &pdi_recf);

  /* Compute conforming norm */
  for (c_id = 0; c_id < geom->n_cells; c_id++) {

    cs_get_t  get;
    double  num_cell = 0.0, denum_cell = 0.0;

    const double  pc = pdi_recc[c_id];

    /* Get cell center */
    for (k = 0; k < 3; k++)
      xc[k] = geom->cell_centers[3*c_id+k];

    for (i = topo->c2f->idx[c_id]; i < topo->c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = topo->c2f->col_id[i];
      const cs_quant_t  fq = geom->face[f_id]; // Face quantities
      const double  pf = pdi_recf[f_id];

      for (j = topo->f2e->idx[f_id]; j < topo->f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  e_id = topo->f2e->col_id[j];
        const cs_lnum_t  v_id1 = topo->e2v->col_id[topo->e2v->idx[e_id]  ];
        const cs_lnum_t  v_id2 = topo->e2v->col_id[topo->e2v->idx[e_id]+1];
        const double  pv1 = pdi[v_id1], pv2 = pdi[v_id2];

        const double  voltet = cs_math_voltet(&(xyz[3*v_id1]),
                                              &(xyz[3*v_id2]),
                                              fq.center,
                                              xc);

        /* Analytical function and integration with the highest
           available quadrature */
        cs_quadrature_tet_5pts(&(xyz[3*v_id1]), &(xyz[3*v_id2]), fq.center, xc,
                               voltet, gpts, weights);

        for (p = 0; p < 5; p++) {
          get_sol(tcur, gpts[p], &get);
          rpex_gpts[p] = get.val;
        }

        /* Should keep the same order */
        pdi_gpts[0] = cs_math_onesix *(pv1 + pv2 + pf) + 0.5*pc;
        pdi_gpts[1] = cs_math_onesix *(pv2 + pf + pc)  + 0.5*pv1;
        pdi_gpts[2] = cs_math_onesix *(pf + pc + pv1)  + 0.5*pv2;
        pdi_gpts[3] = cs_math_onesix *(pc + pv1 + pv2) + 0.5*pf;
        pdi_gpts[4] = 0.25*(pv1 + pv2 + pf + pc);

        n_add = 0, d_add = 0;
        for (p = 0; p < 5; p++) { /* Loop on Gauss points */
          ddip_gpts = rpex_gpts[p] - pdi_gpts[p];
          n_add += weights[p] * ddip_gpts*ddip_gpts;
          d_add += weights[p] * rpex_gpts[p]*rpex_gpts[p];
        }
        num_cell += n_add;
        denum_cell += d_add;

      } /* Loop on face edges */

    } /* Loop on cell faces */

    _denum += denum_cell;
    _num += num_cell;

  } /* Loop on cells */

  *num = _num;
  *denum = _denum;

  BFT_FREE(pdi_recc);
  BFT_FREE(pdi_recf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]  m          pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a time step structure
 * \param[in]  eq         pointer to a cs_equation_t structure
 * \param[in]  anacomp    do an analytic comparison or not
 */
/*----------------------------------------------------------------------------*/

static void
_cdovb_post(const cs_mesh_t            *m,
            const cs_cdo_connect_t     *connect,
            const cs_cdo_quantities_t  *cdoq,
            const cs_time_step_t       *time_step,
            const cs_equation_t        *eq,
            bool                        anacomp)
{
  cs_data_info_t  dinfo;
  int  i, len, work_size;
  double  num, denum, l2pot, l2dpot, enerd, l2dgrd;
  cs_get_t  get;

  char  *postlabel = NULL;
  double  *ddip = NULL, *ddig = NULL, *rpex = NULL, *gdi = NULL, *rgex = NULL;
  double  *pvol = NULL, *work = NULL;

  const double  tcur = time_step->t_cur;
  const cs_lnum_t  n_vertices = cdoq->n_vertices, n_edges = cdoq->n_edges;
  const cs_equation_param_t  *eqp = cs_equation_get_param(eq);
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_field_t  *field = cs_equation_get_field(eq);
  const cs_real_t  *pdi = field->val;

  /* Output a summary of results */
  dinfo = cs_analysis_data(n_vertices,  // n_elts
                           1,           // stride
                           CS_DOUBLE,   // cs_datatype_t
                           pdi,         // data
                           false);      // compute with absolute values

  fprintf(resume, " -bnd- Scal.Min   % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Scal.Max   % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -bnd- Scal.Mean  % 10.6e\n", dinfo.mean);
  fprintf(resume, " -bnd- Scal.Sigma % 10.6e\n", dinfo.sigma);
  fprintf(resume, "%s", msepline);

  if (anacomp) { /* Comparison with an analytical solution */

    /* Work buffer */
    work_size = CS_MAX(3*cdoq->n_cells, n_vertices);
    work_size = CS_MAX(n_edges, work_size);
    BFT_MALLOC(work, work_size, double);

    /* pex = exact potential
       pdi = discrete potential (solution of the discrete system)
       rpex = Red(vtxsp)(Pex) = reduction on vertices of the exact potential
       ddip = rpex - pdi
    */

    BFT_MALLOC(rpex, n_vertices, double);
    BFT_MALLOC(ddip, n_vertices, double);
    for (i = 0; i < n_vertices; i++) {
      get_sol(tcur, &(m->vtx_coord[3*i]), &get);
      rpex[i] = get.val;
      ddip[i] = rpex[i] - pdi[i];
    }

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

    /* Analyse the exact solution */
    dinfo = cs_analysis_data(n_vertices, // n_elts
                             1,                // stride
                             CS_DOUBLE,        // cs_datatype_t
                             rpex,             // data
                             false);           // compute with absolute values

    fprintf(resume, " -bnd- Ref.Min    % 10.6e\n", dinfo.min.value);
    fprintf(resume, " -bnd- Ref.Max    % 10.6e\n", dinfo.max.value);
    fprintf(resume, " -bnd- Ref.Mean   % 10.6e\n", dinfo.mean);
    fprintf(resume, " -bnd- Ref.Sigma  % 10.6e\n", dinfo.sigma);
    fprintf(resume, "%s", msepline);

    for (i = 0; i < n_vertices; i++)
      work[i] = fabs(ddip[i]);

    sprintf(postlabel, "%s.AbsErr", field->name);
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_ALL_ASSOCIATED,
                             postlabel,
                             1,           /* dim */
                             false,       /* interlace */
                             true,        /* parent mesh */
                             CS_POST_TYPE_cs_real_t,
                             work,        /* values on vertices */
                             time_step);  /* time step structure */

    dinfo = cs_analysis_data(n_vertices, // n_elts
                             1,                // stride
                             CS_DOUBLE,        // cs_datatype_t
                             work,             // data
                             false);           // compute with absolute values

    fprintf(resume, " -cvg- ErrAbs.Scal.Max  % 10.6e\n", dinfo.max.value);
    fprintf(resume, " -cvg- ErrAbs.Scal.Min  % 10.6e\n", dinfo.min.value);

    /* Compute pvol related to vertices */
    BFT_MALLOC(pvol, CS_MAX(n_vertices, n_edges), double);
    cs_compute_pvol_vtx(connect, cdoq, &pvol);

    /* Compute discrete L2 error norm on the potential */
    num = cs_weighted_sum_square(n_vertices, ddip, pvol);
    denum = cs_weighted_sum_square(n_vertices, rpex, pvol);
    if (fabs(denum) > cs_math_zero_threshold)
      l2dpot = sqrt(num/denum);
    else
      l2dpot = sqrt(num);

    /* Continuous L2 norm based on a conforming reconstruction of the potential
       ||*||_2^2 = int_Omega (pex - Reco{vtxsp}(pdi))^2

       Sudivide each cell into tetrehadron (xv, xe, xf, xc) and use Gauss
       quadrature to compute the elementary integrals.
    */

    _compute_vb_l2pot(m, connect, cdoq, time_step, pdi, &num, &denum);
    if (fabs(denum) > cs_math_zero_threshold)
      l2pot = sqrt(num/denum);
    else
      l2pot = sqrt(num);

    fprintf(resume, " -cvg- l2dpot           % 10.6e\n", l2dpot);
    fprintf(resume, " -cvg- l2pot            % 10.6e\n", l2pot);
    fprintf(resume, "%s", msepline);

    printf(" >> l2dpot = %7.4e\n", l2dpot);
    printf(" >> l2pot  = %7.4e\n", l2pot);

    /* Compute norm related to the gradient of the error
       - Discrete L2 norm on the discrete gradient
       - Energy norm

       gdi = the discrete gradient related to pdi
       rgex = the reduction of the exact gradient on each edge
       ddig = rgex - gdi
    */

    cs_sla_matvec(connect->e2v, pdi, &gdi, true);
    cs_sla_matvec(connect->e2v, rpex, &rgex, true);

    BFT_MALLOC(ddig, n_edges, double);
    for (i = 0; i < n_edges; i++)
      ddig[i] = rgex[i] - gdi[i];

    dinfo = cs_analysis_data(n_edges,      // n_elts
                             1,            // stride
                             CS_DOUBLE,    // cs_datatype_t
                             ddig,         // data
                             true);        // compute with absolute values

    /* Discrete L2 norm for the error on the gradient
       ||| * |||_{edgesp,2}^2 = Sum_e |pvol_e| (*_e/|e|)^2
       |pvol_e| volume related to each edge -> store in pvol  */
    cs_compute_pvol_edge(connect, cdoq, &pvol);
    for (i = 0; i < n_edges; i++)
      work[i] = ddig[i]/cdoq->edge[i].meas;
    num = cs_weighted_sum_square(n_edges, work, pvol);
    for (i = 0; i < n_edges; i++)
      work[i] = rgex[i]/cdoq->edge[i].meas;
    denum = cs_weighted_sum_square(n_edges, work, pvol);
    if (fabs(denum) > cs_math_zero_threshold)
      l2dgrd = sqrt(num/denum);
    else
      l2dgrd = sqrt(num);
    printf(" >> l2dgrd = %7.4e\n", l2dgrd);

    /* Energy norm^2 = [ ddig, Hodge*ddig ]_EpFd / [ rgex, Hodge*rgex]_EpFd */
    if (h_info.algo == CS_PARAM_HODGE_ALGO_WBS) {

      // TODO
    }
    else {

      cs_hodge_matvec(connect, cdoq, h_info, eqp->diffusion_property,
                      ddig, work);
      num = cs_dot(n_edges, ddig, work);
      cs_hodge_matvec(connect, cdoq, h_info, eqp->diffusion_property,
                      rgex, work);
      denum = cs_dot(n_edges, rgex, work);

    }

    if (fabs(denum) > cs_math_zero_threshold)
      enerd = sqrt(num/denum);
    else
      enerd = sqrt(num);

    printf(" >> enerd  = %7.4e\n", enerd);

    /* Output results */
    fprintf(resume, " -cvg- ErrAbs.Grd.Max  % 10.6e\n", dinfo.max.value);
    fprintf(resume, " -cvg- ErrAbs.Grd.Min  % 10.6e\n", dinfo.min.value);
    fprintf(resume, " -cvg- l2dgrd          % 10.6e\n", l2dgrd);
    fprintf(resume, " -cvg- enerd           % 10.6e\n", enerd);
    fprintf(resume, "%s", msepline);

    /* Post-processing of reconstructed vector fields at each cell center */
    cs_reco_ccen_edge_dofs(connect, cdoq, gdi, &work);

    sprintf(postlabel, "%s.GrdRec", field->name);
    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      postlabel,
                      3,               // dim
                      true,            // interlace
                      true,            // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      work,            // values on cells
                      NULL,            // values at internal faces
                      NULL,            // values at border faces
                      time_step);      // time step management structure

    cs_reco_ccen_edge_dofs(connect, cdoq, ddig, &work);

    sprintf(postlabel, "%s.ErrGrd", field->name);
    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      postlabel,
                      3,               // dim
                      true,            // interlace
                      true,            // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      work,            // values on cells
                      NULL,            // values at internal faces
                      NULL,            // values at border faces
                      time_step);      // time step management structure

    dinfo = cs_analysis_data(cdoq->n_cells,  // n_elts
                             3,               // stride
                             CS_DOUBLE,       // cs_datatype_t
                             work,            // data
                             true);           // compute with absolute values

    fprintf(resume, " -cvg- ErrAbs.GrdReco.Max  % 10.6e\n", dinfo.max.value);
    fprintf(resume, " -cvg- ErrAbs.GrdReco.Min  % 10.6e\n", dinfo.min.value);

    /* Free */
    BFT_FREE(postlabel);
    BFT_FREE(pvol);
    BFT_FREE(work);
    BFT_FREE(ddip);
    BFT_FREE(rpex);
    BFT_FREE(ddig);
    BFT_FREE(rgex);
    BFT_FREE(gdi);

  }

}

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
cs_user_cdo_extra_op(const cs_domain_t          *domain)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  const cs_mesh_t  *m = domain->mesh;
  const cs_cdo_connect_t  *connect = domain->connect;
  const cs_cdo_quantities_t  *cdoq = domain->cdo_quantities;
  const cs_time_step_t  *time_step = domain->time_step;

  cs_equation_t  *eq = cs_domain_get_equation(domain, "FVCA6.1");
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
  cs_space_scheme_t  space_scheme = cs_equation_get_space_scheme(eq);

  _dump_info(cdoq, space_scheme, eqp);

  switch (space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    _cdovb_post(m, connect, cdoq, time_step, eq, true);
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
