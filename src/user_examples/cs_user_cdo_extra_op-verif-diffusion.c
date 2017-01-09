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
         cs_lnum_t          n_pts,
         const cs_real_t   *xyz,
         cs_real_t         *retval)
{
  CS_UNUSED(time);

  const double  pi = 4.0*atan(1.0);

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_real_t  *_xyz = xyz + 3*p;
    const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

    retval[p] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+cs_math_onethird));

  }

}

static cs_analytic_func_t *get_sol = _get_sol;

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
  CS_UNUSED(connect);

  cs_data_info_t  dinfo;
  int  i, len;
  cs_get_t  get;

  char  *postlabel = NULL;
  double  *ddip = NULL, *rpex = NULL;

  const double  tcur = time_step->t_cur;
  const cs_lnum_t  n_vertices = cdoq->n_vertices;
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

    /* pex = exact potential
       pdi = discrete potential (solution of the discrete system)
       rpex = Red(vtxsp)(Pex) = reduction on vertices of the exact potential
       ddip = rpex - pdi
    */

    BFT_MALLOC(rpex, n_vertices, double);
    BFT_MALLOC(ddip, n_vertices, double);
    get_sol(tcur, n_vertices, cdoq->vtx_coord, rpex);
    for (i = 0; i < n_vertices; i++)
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

    /* Free */
    BFT_FREE(postlabel);
    BFT_FREE(ddip);
    BFT_FREE(rpex);

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
