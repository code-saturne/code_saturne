/*============================================================================
 * Coal combustion model.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal.c

  \brief Coal combustion model.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Coal combustion model parameters structure */

cs_coal_model_t  *cs_glob_coal_model = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cpincl_coal_get_pointers(int     **ncharb,
                              int     **nclacp,
                              int     **nclpch,
                              int     **idrift,
                              int     **ich,
                              int     **ick,
                              int     **iash,
                              int     **iwat,
                              double  **ehsoli,
                              double  **wmols,
                              double  **eh0sol,
                              int     **ichcor,
                              double  **diam20,
                              double  **dia2mn,
                              double  **rho20,
                              double  **rho2mn,
                              double  **xmp0,
                              double  **xmasch);

void
cs_f_cpincl_get_pointers_1(double  **cch,
                           double  **hch ,
                           double  **och,
                           double  **sch,
                           double  **nch,
                           double  **alpha,
                           double  **beta,
                           double  **teta,
                           double  **omega,
                           double  **pcich,
                           double  **rho0ch,
                           double  **thcdch,
                           double  **cck,
                           double  **hck,
                           double  **ock,
                           double  **sck,
                           double  **nck,
                           double  **gamma,
                           double  **delta,
                           double  **kappa,
                           double  **zeta,
                           double  **rhock,
                           double  **pcick,
                           double  **cpashc,
                           double  **h0ashc,
                           double  **h02ch,
                           double  **cp2wat,
                           double  **crepn1,
                           double  **crepn2,
                           double  **cp2ch,
                           double  **xashsec,
                           double  **xashch,
                           double  **xwatch);

void
cs_f_cpincl_get_pointers_2(int     **iy1ch,
                           int     **iy2ch,
                           int     **iochet,
                           int     **ioetc2,
                           int     **ioetwt,
                           double  **y1ch,
                           double  **a1ch,
                           double  **e1ch,
                           double  **y2ch,
                           double  **a2ch,
                           double  **e2ch,
                           double  **ahetch,
                           double  **ehetch,
                           double  **ahetc2,
                           double  **ehetc2,
                           double  **ahetwt,
                           double  **ehetwt);

void
cs_f_cpincl_get_pointers_3(int     **ico,
                           int     **ico2,
                           int     **ih2o,
                           int     **io2,
                           int     **in2,
                           int     **ichx1c,
                           int     **ichx2c,
                           int     **ichx1,
                           int     **ichx2,
                           double  **chx1,
                           double  **chx2,
                           double  **a1,
                           double  **b1,
                           double  **c1,
                           double  **d1,
                           double  **e1,
                           double  **f1,
                           double  **a2,
                           double  **b2,
                           double  **c2,
                           double  **d2,
                           double  **e2,
                           double  **f2,
                           double  **thc,
                           int     **npoc);

void
cs_f_coal_incl_get_pointers(int     **ihth2o,
                            int     **ighh2o,
                            int     **ipci,
                            double  **qpr,
                            double  **fn,
                            double  **yhcnle,
                            double  **yhcnlo,
                            double  **ynh3le,
                            double  **ynh3lo,
                            double  **repnle,
                            double  **repnlo,
                            double  **repnck,
                            double  **yhcnc1,
                            double  **ynoch1,
                            double  **yhcnc2,
                            double  **ynoch2,
                            double  **wmchx1,
                            double  **wmchx2);

void
cs_f_co_models_init(void);

void
cs_f_cp_model_map_coal(void);

void
cs_f_coal_incl_init(void);

void
cs_f_ppcpfu_models_init(void);

void
cs_f_thch_models_init(void);

void
cs_f_coal_radst(int         id,
                cs_real_t  *smbrs,
                cs_real_t  *rovsdt);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coal model.
 *
 * \pram[in, out]  cm  pointer to coal model pointer to destroy.
 */
/*----------------------------------------------------------------------------*/

static void
_coal_model_finalize(void)
{
  BFT_FREE(cs_glob_coal_model);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_coal_get_pointers(int     **ncharb,
                              int     **nclacp,
                              int     **nclpch,
                              int     **idrift,
                              int     **ich,
                              int     **ick,
                              int     **iash,
                              int     **iwat,
                              double  **ehsoli,
                              double  **wmols,
                              double  **eh0sol,
                              int     **ichcor,
                              double  **diam20,
                              double  **dia2mn,
                              double  **rho20,
                              double  **rho2mn,
                              double  **xmp0,
                              double  **xmasch)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *ncharb = &(cm->n_coals);
  *nclacp = &(cm->nclacp);
  *nclpch = cm->n_classes_per_coal;
  *idrift = &(cm->idrift);

  *ich    = cm->ich;
  *ick    = cm->ick;
  *iash   = cm->iash;
  *iwat   = cm->iwat;
  *ehsoli = (double *)cm->ehsoli;
  *wmols  = cm->wmols;
  *eh0sol = cm->eh0sol;

  *ichcor = cm->ichcor;
  *diam20 = cm->diam20;
  *dia2mn = cm->dia2mn;
  *rho20  = cm->rho20;
  *rho2mn = cm->rho2mn;
  *xmp0   = cm->xmp0;
  *xmasch = cm->xmasch;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers_1(double  **cch,
                           double  **hch ,
                           double  **och,
                           double  **sch,
                           double  **nch,
                           double  **alpha,
                           double  **beta,
                           double  **teta,
                           double  **omega,
                           double  **pcich,
                           double  **rho0ch,
                           double  **thcdch,
                           double  **cck,
                           double  **hck,
                           double  **ock,
                           double  **sck,
                           double  **nck,
                           double  **gamma,
                           double  **delta,
                           double  **kappa,
                           double  **zeta,
                           double  **rhock,
                           double  **pcick,
                           double  **cpashc,
                           double  **h0ashc,
                           double  **h02ch,
                           double  **cp2wat,
                           double  **crepn1,
                           double  **crepn2,
                           double  **cp2ch,
                           double  **xashsec,
                           double  **xashch,
                           double  **xwatch)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *cch = cm->cch;
  *hch  = cm->hch;
  *och = cm->och;
  *sch = cm->sch;
  *nch = cm->nch;
  *alpha = cm->alpha;
  *beta = cm->beta;
  *teta = cm->teta;
  *omega = cm->omega;
  *pcich = cm->pcich;
  *rho0ch = cm->rho0ch;
  *thcdch = cm->thcdch;
  *cck = cm->cck;
  *hck = cm->hck;
  *ock = cm->ock;
  *sck = cm->sck;
  *nck = cm->nck;
  *gamma = cm->gamma;
  *delta = cm->delta;
  *kappa = cm->kappa;
  *zeta = cm->zeta;
  *rhock = cm->rhock;
  *pcick = cm->pcick;
  *cpashc = cm->cpashc;
  *h0ashc = cm->h0ashc;
  *h02ch = cm->h02ch;
  *cp2wat = cm->cp2wat;
  *crepn1 = (double *)cm->crepn1;
  *crepn2 = (double *)cm->crepn2;
  *cp2ch  = cm->cp2ch;
  *xashsec = cm->xashsec;
  *xashch = cm->xashch;
  *xwatch = cm->xwatch;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers_2(int     **iy1ch,
                           int     **iy2ch,
                           int     **iochet,
                           int     **ioetc2,
                           int     **ioetwt,
                           double  **y1ch,
                           double  **a1ch,
                           double  **e1ch,
                           double  **y2ch,
                           double  **a2ch,
                           double  **e2ch,
                           double  **ahetch,
                           double  **ehetch,
                           double  **ahetc2,
                           double  **ehetc2,
                           double  **ahetwt,
                           double  **ehetwt)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *iy1ch = cm->iy1ch;
  *iy2ch = cm->iy2ch;
  *iochet = cm->iochet;
  *ioetc2 = cm->ioetc2;
  *ioetwt = cm->ioetwt;
  *y1ch = cm->y1ch;
  *a1ch = cm->a1ch;
  *e1ch = cm->e1ch;
  *y2ch = cm->y2ch;
  *a2ch = cm->a2ch;
  *e2ch = cm->e2ch;
  *ahetch = cm->ahetch;
  *ehetch = cm->ehetch;
  *ahetc2 = cm->ahetc2;
  *ehetc2 = cm->ehetc2;
  *ahetwt = cm->ahetwt;
  *ehetwt = cm->ehetwt;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ico    --> pointer to cm->ico
 *   ico2   --> pointer to cm->ico2
 *   ih2o   --> pointer to cm->ih2o
 *   io2    --> pointer to cm->io2
 *   in2    --> pointer to cm->in2
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers_3(int     **ico,
                           int     **ico2,
                           int     **ih2o,
                           int     **io2,
                           int     **in2,
                           int     **ichx1c,
                           int     **ichx2c,
                           int     **ichx1,
                           int     **ichx2,
                           double  **chx1,
                           double  **chx2,
                           double  **a1,
                           double  **b1,
                           double  **c1,
                           double  **d1,
                           double  **e1,
                           double  **f1,
                           double  **a2,
                           double  **b2,
                           double  **c2,
                           double  **d2,
                           double  **e2,
                           double  **f2,
                           double  **thc,
                           int     **npoc)
{
  if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *ico   = &(cm->ico);
    *io2   = &(cm->io2);
    *ico2   = &(cm->ico2);
    *ih2o   = &(cm->ih2o);
    *in2   = &(cm->in2);

    *ichx1c = cm->ichx1c;
    *ichx2c = cm->ichx2c;
    *ichx1 = &(cm->ichx1);
    *ichx2 = &(cm->ichx2);
    *chx1 = cm->chx1;
    *chx2 = cm->chx2;
    *a1 = cm->a1;
    *b1 = cm->b1;
    *c1 = cm->c1;
    *d1 = cm->d1;
    *e1 = cm->e1;
    *f1 = cm->f1;
    *a2 = cm->a2;
    *b2 = cm->b2;
    *c2 = cm->c2;
    *d2 = cm->d2;
    *e2 = cm->e2;
    *f2 = cm->f2;

    *thc = cm->thc;
    *npoc = &(cm->npoc);

  }
}


/*----------------------------------------------------------------------------
 * Get pointers to members of the global compbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_coal_incl_get_pointers(int     **ihth2o,
                            int     **ighh2o,
                            int     **ipci,
                            double  **qpr,
                            double  **fn,
                            double  **yhcnle,
                            double  **yhcnlo,
                            double  **ynh3le,
                            double  **ynh3lo,
                            double  **repnle,
                            double  **repnlo,
                            double  **repnck,
                            double  **yhcnc1,
                            double  **ynoch1,
                            double  **yhcnc2,
                            double  **ynoch2,
                            double  **wmchx1,
                            double  **wmchx2)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *ihth2o = &(cm->ihth2o);
  *ighh2o = cm->ighh2o;
  *ipci = cm->ipci;

  *qpr = cm->qpr;
  *fn = cm->fn;
  *yhcnle = cm->yhcnle;
  *yhcnlo = cm->yhcnlo;
  *ynh3le = cm->ynh3le;
  *ynh3lo = cm->ynh3lo;
  *repnle = cm->repnle;
  *repnlo = cm->repnlo;
  *repnck = cm->repnck;
  *yhcnc1 = cm->yhcnc1;
  *ynoch1 = cm->ynoch1;
  *yhcnc2 = cm->yhcnc2;
  *ynoch2 = cm->ynoch2;
  *wmchx1 = &(cm->wmchx1);
  *wmchx2 = &(cm->wmchx2);
}

/*----------------------------------------------------------------------------
 * Take in account the radiative source terms in the particle equation
 * of a given class for pulverized coal flame.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id      <-- field id
 *   smbrs   <-- right side of the system
 *   rovsdt  <-- system diagonal
 *----------------------------------------------------------------------------*/

void
cs_f_coal_radst(int         id,
                cs_real_t  *smbrs,
                cs_real_t  *rovsdt)
{
  const cs_field_t *f = cs_field_by_id(id);

  cs_coal_rad_transfer_st(f, smbrs, rovsdt);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return coal combustion model type.
 *
 * \return type of active coal combustion model
 *         (CS_COMBUSTION_COAL_NONE if model not active)
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_type_t
cs_coal_model_get_type(void)
{
  if (cs_glob_coal_model == NULL)
   return CS_COMBUSTION_COAL_NONE;
  else
    return cs_glob_coal_model->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate coal combustion model.
 *
 * \return  pointer to coal combustion model structure.
 *
 * \param[in]  type  coal combustion model type
 */
/*----------------------------------------------------------------------------*/

cs_coal_model_t *
cs_coal_model_set_model(cs_coal_model_type_t  type)
{
  cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = type;

  if (type == CS_COMBUSTION_COAL_NONE) {
    BFT_FREE(cs_glob_coal_model);
    return NULL;
  }
  else if (cs_glob_coal_model != NULL) {
    cs_glob_coal_model->type = type;
    return cs_glob_coal_model;
  }

  /* Create and initialize model structure */

  cs_coal_model_t *cm = NULL;

  BFT_MALLOC(cm, 1, cs_coal_model_t);
  memset(cm, 0, sizeof(cs_coal_model_t));

  cs_glob_coal_model = cm;

  /* Members whose initial value is not 0 */

  cm->type = type;

  cm->ieqnox = 1;
  cm->ieqco2 = 1;

  cm->ico = -1;
  cm->io2 = -1;
  cm->in2 = -1;
  cm->ico2 = -1;
  cm->ih2o = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++)
    cm->ighh2o[i] = -1;

  /* Set finalization callback */

  cs_base_at_finalize(_coal_model_finalize);

  /* Set mappings with Fortran */

  cs_f_co_models_init();
  cs_f_cp_model_map_coal();
  cs_f_coal_incl_init();
  cs_f_ppcpfu_models_init();
  cs_f_thch_models_init();

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take in account the radiative source terms in the particle equation
 *        of a given class for pulverized coal flame.
 *
 * \param[in]      f       pointer to scalar field
 * \param[in, out] smbrs   right and side (explicit ST contribution)
 * \param[in, out] rovsdt  system diagonal (implicit ST contribution)
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_rad_transfer_st(const cs_field_t  *f,
                        cs_real_t         *smbrs,
                        cs_real_t         *rovsdt)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  /* Initialization
   * -------------- */

  const int keyccl = cs_field_key_id("scalar_class");
  const int numcla = cs_field_get_key_int(f, keyccl);
  const int ipcl   = 1 + numcla;

  /* Radiative source terms
   * ---------------------- */

  char f_name[64];

  snprintf(f_name, 63, "rad_st_implicit_%02d", ipcl);
  cs_real_t *cpro_tsri = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "rad_st_%02d", ipcl);
  cs_real_t *cpro_tsre = cs_field_by_name(f_name)->val;

  snprintf(f_name, 63, "x_p_%02d", numcla);
  const cs_real_t *cval_x_p = cs_field_by_name(f_name)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cpro_tsri[c_id] = cs_math_fmax(-cpro_tsri[c_id], 0.);

    if (cval_x_p[c_id] > cs_math_epzero) {
      /* Explicit part */
      smbrs[c_id] += cpro_tsre[c_id]*cell_f_vol[c_id]*cval_x_p[c_id];

      /* Implicit part */
      rovsdt[c_id] += cpro_tsri[c_id]*cell_f_vol[c_id];
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
