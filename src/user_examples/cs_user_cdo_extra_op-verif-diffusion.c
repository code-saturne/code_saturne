/*============================================================================
 * Additional post-processing functions defined by user related to CDO schemes
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

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_field.h"
#include "cs_cdo.h"
#include "cs_reco.h"
#include "cs_quadrature.h"
#include "cs_cdo_toolbox.h"
#include "cs_param.h"
#include "cs_param_eq.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_cdofb_codits.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_user_cdo_post.c

  \brief Additional user-defined post-processing and analysis functions
*/

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
_get_sol(cs_real_t    time,
         cs_real_3_t  xyz,
         cs_get_t    *retval)
{
  double  solution;

  const double  x = xyz[0], y = xyz[1], z = xyz[2];
  const double  pi = 4.0*atan(1.0);

  solution = 1+sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+1./3.));

  (*retval).val = solution;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the discrete L2 norm and energy norm of the error on the
 *         gradient in face-based scheme
 *
 *
 * \param[in] topo       pointer to the connectivity struct.
 * \param[in] cdoq       pointer to the additional quantities struct.
 * \param[in] h_info     information about the discrete Hodge op.
 * \param[in] cell_rpex  reduction of the exact solution at cell centers
 * \param[in] face_rpex  reduction of the exact solution at face centers
 * \param[in] cell_pdi   computed solution at cell centers
 * \param[in] face_pdi   computed solution at face centers
 */
/*----------------------------------------------------------------------------*/

static void
_compute_fb_errgrd(const cs_cdo_connect_t      *topo,
                   const cs_cdo_quantities_t   *cdoq,
                   const cs_param_hodge_t       h_info,
                   const double                 cell_rpex[],
                   const double                 face_rpex[],
                   const double                 cell_pdi[],
                   const double                 face_pdi[])
{
  short int sgn;
  int  i, j, k, ij, c_id, f_id;
  double  gcontrib, pvol, l2dgrd, enerd, gdic;
  cs_real_3_t  xc;
  cs_quant_t  fq;

  double  num_l2d = 0, denum_l2d = 0, num_end = 0, denum_end = 0;
  double  *gexc = NULL, *dgc = NULL;

  cs_toolbox_locmat_t  *_h = cs_toolbox_locmat_create(topo->n_max_fbyc);
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(topo->n_max_fbyc);

  const double  _over3 = 1./3.;

  /* Initialize local buffers */
  BFT_MALLOC(gexc, 2*topo->n_max_fbyc, double);
  dgc = gexc + topo->n_max_fbyc;
  for (i = 0; i < 2*topo->n_max_fbyc; i++) gexc[i] = 0;

  /* Loop on cells */
  for (c_id = 0; c_id < cdoq->n_cells; c_id++) {

    double  _nl2c = 0, _dl2c = 0, _denc = 0, _nenc = 0;

    /* Build a local discrete Hodge operator */
    if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
      cs_hodge_cost_build_local(c_id, topo, cdoq, h_info, _h, hb);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Incompatible choice of discrete Hodge operator\n"));

    for (k = 0; k < 3; k++)
      xc[k] = cdoq->cell_centers[3*c_id+k];

    for (i = 0; i < _h->n_ent; i++) {

      cs_lnum_t  shift = topo->c2f->idx[c_id] + i;

      const double  ed_meas =  cdoq->dedge[4*shift];  /* Dual edge quantities */
      const double  *ed_uvec = &(cdoq->dedge[4*shift+1]);

      f_id = _h->ids[i];
      assert(f_id == topo->c2f->col[shift]-1);

      fq = cdoq->face[f_id];
      sgn = topo->c2f->sgn[shift];

      /* Compute pvol_{f,c} */
      pvol = _over3 * ed_meas * fq.meas * _dp3(ed_uvec, fq.unitv);
      gcontrib = pvol/(ed_meas*ed_meas);

      gexc[i] = sgn*(face_rpex[f_id] - cell_rpex[c_id]);
      gdic = sgn*(face_pdi[f_id] - cell_pdi[c_id]);
      dgc[i] = gexc[i] - gdic;

      _nl2c += gcontrib * dgc[i]*dgc[i];
      _dl2c += gcontrib * gexc[i]*gexc[i];

    } // Loop on cell faces

    for (i = 0; i < _h->n_ent; i++) {
      for (j = 0; j < _h->n_ent; j++) {
        ij = i*_h->n_ent+j;
        _nenc += dgc[i]*_h->mat[ij]*dgc[j];
        _denc += gexc[i]*_h->mat[ij]*gexc[j];
      }
    }

    num_l2d += _nl2c;
    denum_l2d += _dl2c;
    num_end += _nenc;
    denum_end += _denc;

  } /* Loop on cells */

  /* Compute the L2 discrete norm on gradient */
  if (fabs(denum_l2d) > cs_base_zthreshold)
    l2dgrd = sqrt(num_l2d/denum_l2d);
  else
    l2dgrd = sqrt(num_l2d);
  printf(" >> l2dgrd       % 10.6e\n", l2dgrd);

  /* Compute the discrete energy norm */
  if (fabs(denum_end) > cs_base_zthreshold)
    enerd = sqrt(num_end/denum_end);
  else
    enerd = sqrt(num_end);
  printf(" >> enerd        % 10.6e\n", enerd);

  /* Output results */
  fprintf(resume, " -cvg- l2dgrd              % 10.6e\n", l2dgrd);
  fprintf(resume, " -cvg- enerd               % 10.6e\n", enerd);

  /* Free buffers */
  _h = cs_toolbox_locmat_free(_h);
  hb = cs_hodge_builder_free(hb);
  BFT_FREE(gexc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the L2 norm for a potential in vertex-based scheme from
 *         a conforming reconstruction
 *
 * \param[in]   topo     pointer to the connectivity struct.
 * \param[in]   geom     pointer to the additional quantities struct.
 * \param[in]   pdi      pointer to the field of vtx-based DoFs
 * \param[out]  num      int_Omega (pex - Reco{vtxsp}(pdi))^2
 * \param[out]  denum    int_Omega pex^2
 */
/*----------------------------------------------------------------------------*/

static void
_compute_vb_l2pot(const cs_mesh_t             *m,
                  const cs_cdo_connect_t      *topo,
                  const cs_cdo_quantities_t   *geom,
                  const double                *pdi,
                  double                      *num,
                  double                      *denum)
{
  int  i, j, k, p, c_id;
  double  ddip_gpts, voltet, n_add, d_add;
  cs_real_3_t  xc, gpts[5];
  double  rpex_gpts[5], weights[5], pdi_gpts[5];

  double  _num = 0.0, _denum = 0.0, tcur = 0.0;
  double  *pdi_recc = NULL, *pdi_recf = NULL;

  const double  _over_six = 1./6.;
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

      const cs_lnum_t  f_id = topo->c2f->col[i]-1;
      const cs_quant_t  fq = geom->face[f_id]; // Face quantities
      const double  pf = pdi_recf[f_id];

      for (j = topo->f2e->idx[f_id]; j < topo->f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  e_id = topo->f2e->col[j]-1;
        const cs_lnum_t  v_id1 = topo->e2v->col[topo->e2v->idx[e_id]  ] - 1;
        const cs_lnum_t  v_id2 = topo->e2v->col[topo->e2v->idx[e_id]+1] - 1;
        const double  pv1 = pdi[v_id1], pv2 = pdi[v_id2];

        voltet = cs_voltet(&(xyz[3*v_id1]), &(xyz[3*v_id2]), fq.center, xc);

        /* Analytical function and integration with the highest
           available quadrature */
        cs_quadrature_tet_5pts(&(xyz[3*v_id1]), &(xyz[3*v_id2]), fq.center, xc,
                               voltet, gpts, weights);

        for (p = 0; p < 5; p++) {
          _get_sol(tcur, gpts[p], &get);
          rpex_gpts[p] = get.val;
        }

        /* Should keep the same order */
        pdi_gpts[0] = _over_six *(pv1 + pv2 + pf) + 0.5*pc;
        pdi_gpts[1] = _over_six *(pv2 + pf + pc)  + 0.5*pv1;
        pdi_gpts[2] = _over_six *(pf + pc + pv1)  + 0.5*pv2;
        pdi_gpts[3] = _over_six *(pc + pv1 + pv2) + 0.5*pf;
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
 * \param[in]  m        pointer to a cs_mesh_t struct.
 * \param[in]  connect  pointer to a cs_cdo_connect_t struct.
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t struct.
 * \param[in]  eq_id    id of the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

static void
_cdovb_post(const cs_mesh_t            *m,
            const cs_cdo_connect_t     *connect,
            const cs_cdo_quantities_t  *cdoq,
            int                         eq_id)
{
  cs_data_info_t  dinfo;
  int  i, len, work_size;
  double  num, denum, l2pot, l2dpot, enerd, l2dgrd;
  cs_get_t  get;

  char  *postlabel = NULL;
  double  tcur = 0.0;
  double  *ddip = NULL, *ddig = NULL, *rpex = NULL, *gdi = NULL, *rgex = NULL;
  double  *pvol = NULL, *work = NULL;

  const cs_lnum_t  n_vertices = cdoq->n_vertices, n_edges = cdoq->n_edges;
  const cs_param_eq_t  *eq = cs_param_eq_get_by_id(eq_id);
  cs_sla_matrix_t  *H = NULL;

  cs_field_t  *fld = cs_field_by_id(eq->field_id);
  double  *pdi = fld->val;

  /* Work buffer */
  work_size = CS_MAX(3*cdoq->n_cells, n_vertices);
  work_size = CS_MAX(n_edges, work_size);
  BFT_MALLOC(work, work_size, double);

  /* Output a summary of results */
  assert(eq->diffusion_hodge.type == CS_PARAM_HODGE_TYPE_EPFD);

  switch (eq->diffusion_hodge.algo) {
  case CS_PARAM_HODGE_ALGO_COST:
    fprintf(resume, " -hdg- Hodge.Op         EPFD:COST\n");
    H = cs_hodge_cost_build(connect, cdoq, eq->diffusion_hodge);
    break;
  case CS_PARAM_HODGE_ALGO_VORONOI:
    fprintf(resume, " -hdg- Hodge.Op         EPFD:VORONOI\n");
    H = cs_hodge_voronoi_build(connect, cdoq, eq->diffusion_hodge);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid algorithm to define a discrete Hodge operator."));
  } // end of switch

  if (eq->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_COST)
    fprintf(resume, " -hdg- Beta.Coef        %5.3e\n",
            eq->diffusion_hodge.coef);

  fprintf(resume, lsepline);
  fprintf(resume, " -dim- n_vertices       %d\n", n_vertices);
  fprintf(resume, " -dim- n_edges          %d\n", n_edges);
  fprintf(resume, " -dim- n_faces          %d\n", cdoq->n_faces);
  fprintf(resume, " -dim- n_cells          %d\n", cdoq->n_cells);
  fprintf(resume, lsepline);

  cs_post_write_vertex_var(-1,              // id du maillage de post
                           eq->varname,
                           1,               // dim
                           false,           // interlace
                           true,            // true = original mesh
                           CS_POST_TYPE_cs_real_t,
                           pdi,             // values on vertices
                           NULL);           // time step management structure

  dinfo = cs_analysis_data(n_vertices, // n_elts
                           1,                // stride
                           CS_DOUBLE,        // cs_datatype_t
                           pdi,              // data
                           false);           // compute with absolute values

  fprintf(resume, " -bnd- Scal.Min   % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Scal.Max   % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -bnd- Scal.Mean  % 10.6e\n", dinfo.mean);
  fprintf(resume, " -bnd- Scal.Sigma % 10.6e\n", dinfo.sigma);
  fprintf(resume, lsepline);

  /* Comparison with an analytical solution */
  /* ====================================== */

  /* pex = exact potential
     pdi = discrete potential (solution of the discrete system)
     rpex = Red(vtxsp)(Pex) = reduction on vertices of the exact potential
     ddip = rpex - pdi
  */

  BFT_MALLOC(rpex, n_vertices, double);
  BFT_MALLOC(ddip, n_vertices, double);
  for (i = 0; i < n_vertices; i++) {
    _get_sol(tcur, &(m->vtx_coord[3*i]), &get);
    rpex[i] = get.val;
    ddip[i] = rpex[i] - pdi[i];
  }

  len = strlen(eq->varname) + 7 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Error", eq->varname);

  cs_post_write_vertex_var(-1,          /* id du maillage de post */
                           postlabel,
                           1,           /* dim */
                           false,       /* interlace */
                           true,        /* parent mesh */
                           CS_POST_TYPE_cs_real_t,
                           ddip,        /* values on vertices */
                           NULL);       /* time step structure */

  sprintf(postlabel, "%s.RefSol", eq->varname);
  cs_post_write_vertex_var(-1,          /* id du maillage de post */
                           postlabel,
                           1,           /* dim */
                           false,       /* interlace */
                           true,        /* parent mesh */
                           CS_POST_TYPE_cs_real_t,
                           rpex,        /* values on vertices */
                           NULL);       /* time step structure */

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
  fprintf(resume, lsepline);

  for (i = 0; i < n_vertices; i++)
    work[i] = fabs(ddip[i]);

  sprintf(postlabel, "%s.AbsErr", eq->varname);
  cs_post_write_vertex_var(-1,          /* id du maillage de post */
                           postlabel,
                           1,           /* dim */
                           false,       /* interlace */
                           true,        /* parent mesh */
                           CS_POST_TYPE_cs_real_t,
                           work,        /* values on vertices */
                           NULL);       /* time step structure */

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
  num = cs_sum(n_vertices, ddip, pvol, CS_TOOLBOX_WSUM2);
  denum = cs_sum(n_vertices, rpex, pvol, CS_TOOLBOX_WSUM2);
  if (fabs(denum) > cs_base_zthreshold)
    l2dpot = sqrt(num/denum);
  else
    l2dpot = sqrt(num);

  /* Continuous L2 norm based on a conforming reconstruction of the potential
     ||*||_2^2 = int_Omega (pex - Reco{vtxsp}(pdi))^2

     Sudivide each cell into tetrehadron (xv, xe, xf, xc) and use Gauss
     quadrature to compute the elementary integrals.
  */

  _compute_vb_l2pot(m, connect, cdoq, pdi, &num, &denum);
  if (fabs(denum) > cs_base_zthreshold)
    l2pot = sqrt(num/denum);
  else
    l2pot = sqrt(num);

  fprintf(resume, " -cvg- l2dpot           % 10.6e\n", l2dpot);
  fprintf(resume, " -cvg- l2pot            % 10.6e\n", l2pot);
  fprintf(resume, lsepline);

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
  num = cs_sum(n_edges, work, pvol, CS_TOOLBOX_WSUM2);
  for (i = 0; i < n_edges; i++)
    work[i] = rgex[i]/cdoq->edge[i].meas;
  denum = cs_sum(cdoq->n_edges, work, pvol, CS_TOOLBOX_WSUM2);
  if (fabs(denum) > cs_base_zthreshold)
    l2dgrd = sqrt(num/denum);
  else
    l2dgrd = sqrt(num);
  printf(" >> l2dgrd = %7.4e\n", l2dgrd);

  /* Energy norm^2 = [ ddig, Hodge*ddig ]_EpFd / [ rgex, Hodge*rgex]_EpFd */
  cs_sla_matvec(H, ddig, &work, true);
  num = cs_dp(n_edges, ddig, work);
  cs_sla_matvec(H, rgex, &work, true);
  denum = cs_dp(n_edges, rgex, work);
  if (fabs(denum) > cs_base_zthreshold)
    enerd = sqrt(num/denum);
  else
    enerd = sqrt(num);
  printf(" >> enerd  = %7.4e\n", enerd);

  /* Output results */
  fprintf(resume, " -cvg- ErrAbs.Grd.Max  % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrAbs.Grd.Min  % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -cvg- l2dgrd          % 10.6e\n", l2dgrd);
  fprintf(resume, " -cvg- enerd           % 10.6e\n", enerd);
  fprintf(resume, lsepline);

  /* Post-processing of reconstructed vector fields at each cell center */
  cs_reco_ccen_edge_dofs(connect, cdoq, gdi, &work);

  sprintf(postlabel, "%s.GrdRec", eq->varname);
  cs_post_write_var(-1,              // id du maillage de post
                    postlabel,
                    3,               // dim
                    true,            // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    work,            // values on cells
                    NULL,            // values at internal faces
                    NULL,            // values at border faces
                    NULL);           // time step management structure

  cs_reco_ccen_edge_dofs(connect, cdoq, ddig, &work);

  sprintf(postlabel, "%s.ErrGrd", eq->varname);
  cs_post_write_var(-1,              // id du maillage de post
                    postlabel,
                    3,               // dim
                    true,            // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    work,            // values on cells
                    NULL,            // values at internal faces
                    NULL,            // values at border faces
                    NULL);           // time step management structure

  dinfo = cs_analysis_data(cdoq->n_cells,  // n_elts
                           3,               // stride
                           CS_DOUBLE,       // cs_datatype_t
                           work,            // data
                           true);           // compute with absolute values

  fprintf(resume, " -cvg- ErrAbs.GrdReco.Max  % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrAbs.GrdReco.Min  % 10.6e\n", dinfo.min.value);

  /* Free */
  H = cs_sla_matrix_free(H);

  BFT_FREE(postlabel);
  BFT_FREE(pvol);
  BFT_FREE(work);
  BFT_FREE(ddip);
  BFT_FREE(rpex);
  BFT_FREE(ddig);
  BFT_FREE(rgex);
  BFT_FREE(gdi);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO face-based scheme.
 *
 * \param[in]  m        pointer to a cs_mesh_t struct.
 * \param[in]  connect  pointer to a cs_cdo_connect_t struct.
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t struct.
 * \param[in]  eq_id    id of the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

static void
_cdofb_post(const cs_mesh_t            *m,
            const cs_cdo_connect_t     *connect,
            const cs_cdo_quantities_t  *cdoq,
            int                         eq_id)
{
  cs_data_info_t  dinfo;
  int  i, len, work_size;
  double  num, denum, l2dpotc, l2dpotf, l2r0potc;
  cs_get_t  get;

  double  tcur = 0.0;
  char  *postlabel = NULL;
  double  *cell_dpdi = NULL, *face_dpdi = NULL;
  double  *cell_rpex = NULL, *face_rpex = NULL;
  double  *pvol = NULL, *work = NULL;

  const cs_lnum_t  n_cells = cdoq->n_cells;
  const cs_lnum_t  n_faces = cdoq->n_faces;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_param_eq_t  *eq = cs_param_eq_get_by_id(eq_id);

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  const double  *cell_pdi = fld->val;
  const double  *face_pdi = cs_cdofb_codits_get_face_values(eq_id);

  /* Output a summary of results */
  assert(eq->diffusion_hodge.type == CS_PARAM_HODGE_TYPE_EDFP);

  switch (eq->diffusion_hodge.algo) {
  case CS_PARAM_HODGE_ALGO_COST:
    fprintf(resume, " -hdg- Hodge.Op         EDFP:COST\n");
    break;
  case CS_PARAM_HODGE_ALGO_VORONOI:
    fprintf(resume, " -hdg- Hodge.Op         EDFP:VORONOI\n");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid algorithm to define a discrete Hodge operator."));
  } // end of switch

  if (eq->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_COST)
    fprintf(resume, " -hdg- Beta.Coef        %5.3e\n",
            eq->diffusion_hodge.coef);

  fprintf(resume, lsepline);
  fprintf(resume, " -dim- n_vertices       %d\n", cdoq->n_vertices);
  fprintf(resume, " -dim- n_edges          %d\n", cdoq->n_edges);
  fprintf(resume, " -dim- n_faces          %d\n", n_faces);
  fprintf(resume, " -dim- n_cells          %d\n", cdoq->n_cells);
  fprintf(resume, lsepline);

  cs_post_write_var(-1,              // id du maillage de post
                    eq->varname,
                    1,               // dim
                    false,           // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    cell_pdi,        // values on cells
                    NULL,            // values at internal faces
                    NULL,            // values at border faces
                    NULL);           // time step management structure

  dinfo = cs_analysis_data(n_cells,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           cell_pdi,   // data
                           false);     // compute with absolute values

  fprintf(resume, " -bnd- Scal.Cell.Min   % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Scal.Cell.Max   % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -bnd- Scal.Cell.Mean  % 10.6e\n", dinfo.mean);
  fprintf(resume, " -bnd- Scal.Cell.Sigma % 10.6e\n", dinfo.sigma);
  fprintf(resume, lsepline);

  dinfo = cs_analysis_data(n_faces,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           face_pdi,   // data
                           false);     // compute with absolute values

  fprintf(resume, " -bnd- Scal.Face.Min   % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Scal.Face.Max   % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -bnd- Scal.Cell.Mean  % 10.6e\n", dinfo.mean);
  fprintf(resume, " -bnd- Scal.Cell.Sigma % 10.6e\n", dinfo.sigma);
  fprintf(resume, lsepline);

  /* Comparison with an analytical solution */
  /* ====================================== */

  /* pex = exact potential
     pdi = discrete potential (solution of the discrete system)
     rpex = Red(vtxsp)(Pex) = reduction on vertices of the exact potential
     dpdi = rpex - pdi
  */

  BFT_MALLOC(cell_rpex, n_cells + n_faces, double);
  BFT_MALLOC(cell_dpdi, n_cells + n_faces, double);
  face_rpex = cell_rpex + n_cells;
  face_dpdi = cell_dpdi + n_cells;

  for (i = 0; i < n_cells; i++) {
    _get_sol(tcur, &(cdoq->cell_centers[3*i]), &get);
    cell_rpex[i] = get.val;
    cell_dpdi[i] = cell_rpex[i] - cell_pdi[i];
  }

  /* Analyse the exact solution */
  dinfo = cs_analysis_data(n_cells,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           cell_rpex,  // data
                           false);     // compute with absolute values

  fprintf(resume, " -bnd- Ref.Cell.Min    % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Ref.Cell.Max    % 10.6e\n", dinfo.max.value);

  for (i = 0; i < n_faces; i++) {

    cs_quant_t  fq = cdoq->face[i];

    _get_sol(tcur, fq.center, &get);
    face_rpex[i] = get.val;
    face_dpdi[i] = face_rpex[i] - face_pdi[i];
  }

  len = strlen(eq->varname) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.RefSolB", eq->varname);

  cs_post_write_var(-2,                    // id du maillage de post
                    postlabel,
                    1,                     // dim
                    false,                 // interlace
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_rpex + n_i_faces, // values at border faces
                    NULL);                 // time step management structure

  sprintf(postlabel, "%s.Border", eq->varname);
  cs_post_write_var(-2,                    // id du maillage de post
                    "Pot.Border",
                    1,                     // dim
                    false,                 // interlace
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_pdi + n_i_faces,  // values at border faces
                    NULL);                 // time step management structure

  sprintf(postlabel, "%s.ErrBord", eq->varname);
  cs_post_write_var(-2,                    // id du maillage de post
                    postlabel,
                    1,                     // dim
                    false,                 // interlace
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_dpdi + n_i_faces, // values at border faces
                    NULL);                 // time step management structure

  /* Analyse the exact solution */
  dinfo = cs_analysis_data(n_cells,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           cell_rpex,  // data
                           false);     // compute with absolute values

  fprintf(resume, " -bnd- Ref.Face.Min    % 10.6e\n", dinfo.min.value);
  fprintf(resume, " -bnd- Ref.Face.Max    % 10.6e\n", dinfo.max.value);
  fprintf(resume, lsepline);

  /* Post-processing of the error */
  sprintf(postlabel, "%s.Error", eq->varname);
  cs_post_write_var(-1,              // id du maillage de post
                    postlabel,
                    1,               // dim
                    false,           // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    cell_dpdi,       // values on cells
                    NULL,            // values at interior faces
                    NULL,            // values at border faces
                    NULL);           // time step structure

  /* Compute an approximation of int_Omega (delta_p)^2 / int_Omega pex^2
     which approximates the normalized L2 norm */
  num = cs_sum(n_cells, cell_dpdi, cdoq->cell_vol, CS_TOOLBOX_WSUM2);
  denum = cs_sum(n_cells, cell_rpex, cdoq->cell_vol, CS_TOOLBOX_WSUM2);
  if (fabs(denum) > cs_base_zthreshold)
    l2dpotc = sqrt(num/denum);
  else
    l2dpotc = sqrt(num);
  printf(" >> l2dpotc      % 10.6e\n", l2dpotc);

  dinfo = cs_analysis_data(n_cells,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           cell_dpdi,  // data
                           false);     // compute with absolute values

  /* Output results to resume (for analysis) */
  fprintf(resume, " -cvg- ErrPot.Cell.Max     % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrPot.Cell.Min     % 10.6e\n", dinfo.min.value);

  /* Work buffer */
  work_size = CS_MAX(3*n_cells, n_cells + n_faces);
  BFT_MALLOC(work, work_size, double);
  for (i = 0; i < n_cells + n_faces; i++)
    work[i] = fabs(cell_dpdi[i]);

  sprintf(postlabel, "%s.AbsErr", eq->varname);
  cs_post_write_var(-1,              // id du maillage de post
                    postlabel,
                    1,               // dim
                    false,           // interlace
                    true,            // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    work,            // values on cells
                    NULL,            // values at interior faces
                    NULL,            // values at border faces
                    NULL);           // time step structure

  dinfo = cs_analysis_data(n_cells,   // n_elts
                           1,         // stride
                           CS_DOUBLE, // cs_datatype_t
                           work,      // data
                           false);    // compute with absolute values

  fprintf(resume, " -cvg- ErrAbsPot.Cell.Max  % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrAbsPot.Cell.Min  % 10.6e\n", dinfo.min.value);

  /* Compute pvol related to faces */
  BFT_MALLOC(pvol, n_faces, double);
  cs_compute_pvol_face(connect, cdoq, &pvol);

  /* Compute discrete L2 error norm on the potential */
  num = cs_sum(n_faces, face_dpdi, pvol, CS_TOOLBOX_WSUM2);
  denum = cs_sum(n_faces, face_rpex, pvol, CS_TOOLBOX_WSUM2);
  if (fabs(denum) > cs_base_zthreshold)
    l2dpotf = sqrt(num/denum);
  else
    l2dpotf = sqrt(num);
  printf(" >> l2dpotf      % 10.6e\n", l2dpotf);

  dinfo = cs_analysis_data(n_faces,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           face_dpdi,  // data
                           false);     // compute with absolute values

  /* Output results to resume (for analysis) */
  fprintf(resume, lsepline);
  fprintf(resume, " -cvg- ErrPot.Face.Max     % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrPot.Face.Min     % 10.6e\n", dinfo.min.value);

  dinfo = cs_analysis_data(n_faces,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           work + n_cells,  // data
                           false);     // compute with absolute values

  fprintf(resume, " -cvg- ErrAbsPot.Face.Max  % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrAbsPot.Face.Min  % 10.6e\n", dinfo.min.value);
  fprintf(resume, lsepline);

  /* Compute the error in the norm that should be of order 2:
      -> compute the mean volume of the integral over each cell of the
         exact solution: L_cell*R_cell(solu) (using a piecewise constant
         reconstruction)
      -> dpdi = LcRc(solu) - L_vtxd^0(pdi) where L_vtxd^0(pdi) is a piecewise
         constant reconstruction over (primal) cells
      -> Compute the L^2 norm of dpdi
  */

  cs_def_t  def;
  def.analytic = _get_sol;

  cs_flag_t  flag = CS_PARAM_FLAG_CELL|CS_PARAM_FLAG_PRIMAL|CS_PARAM_FLAG_SCAL;

  cs_evaluate(m, cdoq, connect,
              0.0,   //tcur
              flag,  // DoF flag (where to compute the evaluation)
              cs_mesh_location_get_id_by_name("cells"),
              CS_PARAM_DEF_BY_ANALYTIC_FUNCTION,
              CS_QUADRATURE_BARY,
              def,
              &work); // work --> int_cell solu(x,y,z)

  for (i = 0; i < n_cells; i++) {
    work[i] /= cdoq->cell_vol[i];
    cell_dpdi[i] = work[i] - cell_pdi[i];
  }

  num = cs_sum(n_cells, cell_dpdi, cdoq->cell_vol, CS_TOOLBOX_WSUM2);
  denum = cs_sum(n_cells, work, cdoq->cell_vol, CS_TOOLBOX_WSUM2);
  if (fabs(denum) > cs_base_zthreshold)
    l2r0potc = sqrt(num/denum);
  else
    l2r0potc = sqrt(num);
  printf(" >> l2r0potc     % 10.6e\n", l2r0potc);

  dinfo = cs_analysis_data(n_cells,    // n_elts
                           1,          // stride
                           CS_DOUBLE,  // cs_datatype_t
                           cell_dpdi,  // data
                           true);      // compute with absolute values

  /* Output results to resume (for analysis) */
  fprintf(resume, " -cvg- ErrAbs.Pot.Lc0.Max  % 10.6e\n", dinfo.max.value);
  fprintf(resume, " -cvg- ErrAbs.Pot.Lc0.Min  % 10.6e\n", dinfo.min.value);
  fprintf(resume, lsepline);
  fprintf(resume, " -cvg- l2dpotc             % 10.6e\n", l2dpotc);
  fprintf(resume, " -cvg- l2dpotf             % 10.6e\n", l2dpotf);
  fprintf(resume, " -cvg- l2r0potc            % 10.6e\n", l2r0potc);

  /* Compute the L2 discrete norm on the gradient error */
  _compute_fb_errgrd(connect, cdoq, eq->diffusion_hodge,
                     cell_rpex, face_rpex, cell_pdi, face_pdi);

  fprintf(resume, lsepline);

  /* Free */
  BFT_FREE(postlabel);
  BFT_FREE(pvol);
  BFT_FREE(work);
  BFT_FREE(cell_dpdi);
  BFT_FREE(cell_rpex);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional operations on results produced by CDO schemes.
 *         Define advanced post-processing and/or analysis for instance.
 *
 * \param[in]  m        pointer to a cs_mesh_t struct.
 * \param[in]  mq       pointer to a cs_quantities_t struct.
 * \param[in]  connect  pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]  tcur     current physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(const cs_mesh_t             *m,
                     const cs_mesh_quantities_t  *mq,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *cdoq,
                     double                       tcur)
{
  _Bool do_navsto;
  int  eq_id, len;

  char  *filename = NULL;
  int  n_scal_equations = 0;

  /* Retrieve general information */
  cs_param_eq_get_info(&do_navsto, &n_scal_equations);

  /* Loop on equations */
  for (eq_id = 0; eq_id < n_scal_equations; eq_id++) {

    cs_space_scheme_t  space_scheme = cs_param_eq_get_space_scheme(eq_id);

    bft_printf("\n");
    bft_printf(" >> Extra operation for %s equation\n",
               cs_param_eq_get_name(eq_id));

    /* Up to now only this type of equation is handled */
    assert(cs_param_eq_get_type(eq_id) == CS_PARAM_EQ_TYPE_SCAL);

    /* Open a file */
    len = strlen("Resume-.log")+strlen(cs_param_eq_get_name(eq_id))+1;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Resume-%s.log", cs_param_eq_get_name(eq_id));
    resume = fopen(filename, "w");

    /* Build algebraic system */
    switch (space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      _cdovb_post(m, connect, cdoq, eq_id);
      break;

    case CS_SPACE_SCHEME_CDOFB:
      _cdofb_post(m, connect, cdoq, eq_id);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid space scheme. Stop post-processing.\n"));

    } /* space_scheme */

    /* Free */
    BFT_FREE(filename);
    fclose(resume);

  } /* Loop on equations */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
