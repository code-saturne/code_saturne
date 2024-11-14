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

/*!>  molar volume under normal pressure and temperature conditions
   (1 atmosphere, 0 degres C) in m-3 */

/*! reference temperature for molar volume */
const double  cs_coal_trefth = 25. + 273.15;

/*! reference pressure for molar volume */
const double  cs_coal_prefth = 1.01325e5;

/*! molar volume under normal pressure and temperature conditions
  (1 atmosphere, 0 \f$\text{\degresC}\f$) in \f$m^{-3}\f$ */
const double  cs_coal_volmol = 22.41e-3;

/* ids for atom types in wmolat */
const int  cs_coal_atom_id_c = 0;  /*!< id for C in wmolat */
const int  cs_coal_atom_id_h = 1;  /*!< id for H in wmolat */
const int  cs_coal_atom_id_o = 2;  /*!< id for O in wmolat */
const int  cs_coal_atom_id_n = 3;  /*!< id for N in wmolat */
const int  cs_coal_atom_id_s = 4;  /*!< id for S in wmolat */

/* precision for tests */
const double cs_coal_epsilon = 1.e-8;

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
                              int     **nsolid,
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
                              double  **xmash);

void
cs_f_cpincl_get_pointers_1(double  **cch,
                           double  **hch ,
                           double  **och,
                           double  **sch,
                           double  **nch,
                           double  **pcich,
                           double  **rho0ch,
                           double  **thcdch,
                           double  **cck,
                           double  **hck,
                           double  **ock,
                           double  **sck,
                           double  **nck,
                           double  **rhock,
                           double  **pcick,
                           double  **cpashc,
                           double  **h0ashc,
                           double  **h02ch,
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
                           double  **ehetwt,
                           double  **ehgaze);

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
                           double  **f2);

void
cs_f_cpincl_get_pointers_4(int  **ihgas,
                           int  **if1m,
                           int  **if2m,
                           int  **if4m,
                           int  **if5m,
                           int  **if6m,
                           int  **if7m,
                           int  **if8m,
                           int  **if9m,
                           int  **ifvp2m,
                           int  **ixck,
                           int  **ixch,
                           int  **inp,
                           int  **ih2,
                           int  **ixwt,
                           int  **iym1,
                           int  **irom1,
                           int  **immel,
                           int  **itemp2,
                           int  **irom2,
                           int  **idiam2,
                           int  **ix2,
                           int  **igmdch,
                           int  **igmhet,
                           int  **igmtr,
                           int  **ighco2,
                           int  **igmdv1,
                           int  **igmdv2,
                           int  **igmsec,
                           int  **ibcarbone,
                           int  **iboxygen,
                           int  **ibhydrogen);

void
cs_f_cpincl_get_pointers_5(double  **af3,
                           double  **af4,
                           double  **af5,
                           double  **af6,
                           double  **af7,
                           double  **af8,
                           double  **af9,
                           int     **ihy,
                           int     **ih2s,
                           int     **iso2,
                           int     **ihcn,
                           int     **inh3,
                           int     **ieqco2,
                           int     **iyco2,
                           int     **ihtco2,
                           int     **ieqnox,
                           int     **imdnox,
                           int     **irb,
                           int     **iyhcn,
                           int     **iyno,
                           int     **iynh3,
                           int     **ihox,
                           int     **igrb,
                           int     **noxyd,
                           int     **ighcn1,
                           int     **ighcn2,
                           int     **ignoth,
                           int     **ignh31,
                           int     **ignh32,
                           int     **ifhcnd,
                           int     **ifhcnc,
                           int     **ifnh3d,
                           int     **ifnh3c,
                           int     **ifnohc,
                           int     **ifnonh,
                           int     **ifnoch,
                           int     **ifnoth,
                           int     **ifhcnr,
                           int     **icnohc,
                           int     **icnonh,
                           int     **icnorb);

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
cs_f_ppincl_combustion_init(void);

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
                              int     **nsolid,
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
                              double  **xmash)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *ncharb = &(cm->n_coals);
  *nclacp = &(cm->nclacp);
  *nclpch = cm->n_classes_per_coal;
  *idrift = &(cm->idrift);
  *nsolid = &(cm->nsolid);

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
  *xmash  = cm->xmash;
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
                           double  **pcich,
                           double  **rho0ch,
                           double  **thcdch,
                           double  **cck,
                           double  **hck,
                           double  **ock,
                           double  **sck,
                           double  **nck,
                           double  **rhock,
                           double  **pcick,
                           double  **cpashc,
                           double  **h0ashc,
                           double  **h02ch,
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
  *pcich = cm->pcich;
  *rho0ch = cm->rho0ch;
  *thcdch = cm->thcdch;
  *cck = cm->cck;
  *hck = cm->hck;
  *ock = cm->ock;
  *sck = cm->sck;
  *nck = cm->nck;
  *rhock = cm->rhock;
  *pcick = cm->pcick;
  *cpashc = cm->cpashc;
  *h0ashc = cm->h0ashc;
  *h02ch = cm->h02ch;
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
                           double  **ehetwt,
                           double  **ehgaze)
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
  *ehgaze = (double *)cm->ehgaze;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global combbustion model (cpincl).
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
                           double  **f2)
{
  if (cs_glob_coal_model == NULL)
    return;

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
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global combbustion model (cpincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers_4(int  **ihgas,
                           int  **if1m,
                           int  **if2m,
                           int  **if4m,
                           int  **if5m,
                           int  **if6m,
                           int  **if7m,
                           int  **if8m,
                           int  **if9m,
                           int  **ifvp2m,
                           int  **ixck,
                           int  **ixch,
                           int  **inp,
                           int  **ih2,
                           int  **ixwt,
                           int  **iym1,
                           int  **irom1,
                           int  **immel,
                           int  **itemp2,
                           int  **irom2,
                           int  **idiam2,
                           int  **ix2,
                           int  **igmdch,
                           int  **igmhet,
                           int  **igmtr,
                           int  **ighco2,
                           int  **igmdv1,
                           int  **igmdv2,
                           int  **igmsec,
                           int  **ibcarbone,
                           int  **iboxygen,
                           int  **ibhydrogen)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *ihgas = &(cm->ihgas);
  *if1m = cm->if1m;
  *if2m = cm->if2m;
  *if4m = &(cm->if4m);
  *if5m = &(cm->if5m);
  *if6m = &(cm->if6m);
  *if7m = &(cm->if7m);
  *if8m = &(cm->if8m);
  *if9m = &(cm->if9m);
  *ifvp2m = &(cm->ifvp2m);
  *ixck = cm->ixck;
  *ixch = cm->ixch;
  *inp = cm->inp;
  *ih2 = cm->ih2;
  *ixwt = cm->ixwt;
  *iym1 = cm->iym1;
  *irom1 = &(cm->irom1);
  *immel = &(cm->immel);
  *itemp2 = cm->itemp2;
  *irom2 = cm->irom2;
  *idiam2 = cm->idiam2;
  *ix2 = cm->ix2;
  *igmdch = cm->igmdch;
  *igmhet = cm->igmhet;
  *igmtr = cm->igmtr;
  *ighco2 = cm->ighco2;
  *igmdv1 = cm->igmdv1;
  *igmdv2 = cm->igmdv2;
  *igmsec = cm->igmsec;
  *ibcarbone = &(cm->ibcarbone);
  *iboxygen = &(cm->iboxygen);
  *ibhydrogen = &(cm->ibhydrogen);
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
 * Get pointers to members of combustion model (ppcpfu).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_cpincl_get_pointers_5(double  **af3,
                           double  **af4,
                           double  **af5,
                           double  **af6,
                           double  **af7,
                           double  **af8,
                           double  **af9,
                           int     **ihy,
                           int     **ih2s,
                           int     **iso2,
                           int     **ihcn,
                           int     **inh3,
                           int     **ihtco2,
                           int     **ieqco2,
                           int     **iyco2,
                           int     **ieqnox,
                           int     **imdnox,
                           int     **irb,
                           int     **iyhcn,
                           int     **iyno,
                           int     **iynh3,
                           int     **ihox,
                           int     **igrb,
                           int     **noxyd,
                           int     **ighcn1,
                           int     **ighcn2,
                           int     **ignoth,
                           int     **ignh31,
                           int     **ignh32,
                           int     **ifhcnd,
                           int     **ifhcnc,
                           int     **ifnh3d,
                           int     **ifnh3c,
                           int     **ifnohc,
                           int     **ifnonh,
                           int     **ifnoch,
                           int     **ifnoth,
                           int     **ifhcnr,
                           int     **icnohc,
                           int     **icnonh,
                           int     **icnorb)
{
  if (cs_glob_coal_model == NULL)
    return;

  cs_coal_model_t  *cm = cs_glob_coal_model;

  *af3 = cm->af3;
  *af4 = cm->af4;
  *af5 = cm->af5;
  *af6 = cm->af6;
  *af7 = cm->af7;
  *af8 = cm->af8;
  *af9 = cm->af9;

  *ihy = &(cm->ihy);
  *ih2s = &(cm->ih2s);
  *iso2 = &(cm->iso2);
  *ihcn = &(cm->ihcn);
  *inh3 = &(cm->inh3);
  *ihtco2 = &(cm->ihtco2);
  *ieqco2 = &(cm->ieqco2);
  *iyco2 = &(cm->iyco2);
  *ieqnox = &(cm->ieqnox);
  *imdnox = &(cm->imdnox);
  *irb = &(cm->irb);

  *iyhcn = &(cm->iyhcn);
  *iyno = &(cm->iyno);
  *iynh3 = &(cm->iynh3);
  *ihox = &(cm->ihox);
  *igrb = &(cm->igrb);
  *noxyd = &(cm->noxyd);

  *ighcn1 = &(cm->ighcn1);
  *ighcn2 = &(cm->ighcn2);
  *ignoth = &(cm->ignoth);
  *ignh31 = &(cm->ignh31);
  *ignh32 = &(cm->ignh32);
  *ifhcnd = &(cm->ifhcnd);
  *ifhcnc = &(cm->ifhcnc);
  *ifnh3d = &(cm->ifnh3d);
  *ifnh3c = &(cm->ifnh3c);
  *ifnohc = &(cm->ifnohc);
  *ifnonh = &(cm->ifnonh);
  *ifnoch = &(cm->ifnoch);
  *ifnoth = &(cm->ifnoth);
  *ifhcnr = &(cm->ifhcnr);
  *icnohc = &(cm->icnohc);
  *icnonh = &(cm->icnonh);
  *icnorb = &(cm->icnorb);
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

  cm->ihgas = -1;
  cm->iyco2 = -1;
  cm->iyhcn = -1;
  cm->iynh3 = -1;
  cm->iyno  = -1;

  cm->ihox  = -1;

  for (int i = 0; i < CS_COMBUSTION_MAX_COALS; i++) {
    cm->if1m[i] = -1;
    cm->if2m[i] = -1;
  }

  cm->if4m = -1;
  cm->if5m = -1;
  cm->if6m = -1;
  cm->if7m = -1;
  cm->if8m = -1;
  cm->if9m = -1;
  cm->ifvp2m = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++) {
    cm->ixck[i] = -1;
    cm->ixch[i] = -1;
    cm->inp[i] = -1;
    cm->ih2[i] = -1;
    cm->ixwt[i] = -1;
  }

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS; i++)
    cm->iym1[i] = -1;

  cm->irom1 = -1;
  cm->immel = -1;

  for (int i = 0; i < CS_COMBUSTION_COAL_MAX_CLASSES; i++) {
    cm->itemp2[i] = -1;
    cm->irom2[i] = -1;
    cm->idiam2[i] = -1;
    cm->ix2[i] = -1;
    cm->igmdch[i] = -1;
    cm->igmhet[i] = -1;
    cm->igmtr[i] = -1;
    cm->ighco2[i] = -1;
    cm->igmdv1[i] = -1;
    cm->igmdv2[i] = -1;
    cm->igmsec[i] = -1;
  }

  cm->ibcarbone = -1;
  cm->iboxygen = -1;
  cm->ibhydrogen = -1;

  cm->ighcn1 = -1;
  cm->ighcn2 = -1;
  cm->ignoth = -1;
  cm->ignh31 = -1;
  cm->ignh32 = -1;
  cm->ifhcnd = -1;
  cm->ifhcnc = -1;
  cm->ifnh3d = -1;
  cm->ifnh3c = -1;
  cm->ifnohc = -1;
  cm->ifnonh = -1;
  cm->ifnoch = -1;
  cm->ifnoth = -1;
  cm->icnohc = -1;
  cm->icnonh = -1;
  cm->ifhcnr = -1;
  cm->icnorb = -1;

  cm->srrom = 0.95;

  /* Set finalization callback */

  cs_base_at_finalize(_coal_model_finalize);

  /* Set mappings with Fortran */

  cs_f_ppincl_combustion_init();
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
