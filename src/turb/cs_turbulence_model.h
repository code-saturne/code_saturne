#ifndef __CS_TURBULENCE_MODEL_H__
#define __CS_TURBULENCE_MODEL_H__

/*============================================================================
 * Base turbulence model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* turbulence model general options descriptor */
/*---------------------------------------------*/

typedef struct {

  int           iturb;        /* turbulence model
                                 - 0: no turbulence model (laminar flow)
                                 - 10: mixing length model
                                 - 20: standard k-epsilon model
                                 - 21: k-epsilon model with Linear
                                       Production (LP) correction
                                 - 30: Rij-epsilon (LRR)
                                 - 31: Rij-epsilon (SSG)
                                 - 32: Rij-epsilon (EBRSM)
                                 - 40: LES (constant Smagorinsky model)
                                 - 41: LES ("classical" dynamic Smagorisky
                                       model)
                                 - 42: LES (WALE)
                                 - 50: v2f phi-model
                                 - 51: v2f BL-v2-k
                                 - 60: k-omega SST
                                 - 70: Spalart-Allmaras model */
  int           itytur;       /* class of turbulence model (integer value
                                 iturb/10) */
  int           nvarcl;       /* number of variable plus number of turbulent
                                 fluxes (used by the boundary conditions) */

} cs_turb_model_t;

/* rans turbulence model descriptor */
/*----------------------------------*/

typedef struct {

  int           irccor;       /* activation of rotation/curvature correction for
                                 an eddy viscosity turbulence models
                                 - 0: false
                                 - 1: true */
  int           itycor;       /* type of rotation/curvature correction for an
                                 eddy viscosity turbulence models
                                 - 1: Cazalbou correction (default when irccor=1
                                      and itytur=2 or 5)
                                 - 2: Spalart-Shur correction (default when
                                      irccor=1 and iturb=60 or 70) */
  int           idirsm;       /* turbulent diffusion model for second moment
                                 closure
                                 - 0: scalar diffusivity (Shir model)
                                 - 1: tensorial diffusivity (Daly and Harlow
                                      model, default model) */
  int           iclkep;       /* clipping of k and epsilon
                                 - 0: absolute value clipping
                                 - 1: coupled clipping based on physical
                                      relationships */
  int           igrhok;       /* take (2/3 rho grad k) in the momentum
                                 equation
                                 - 1: true
                                 - 0: false (default) */
  int           igrake;       /* buoyant term in k-epsilon
                                 - 1: true (default if rho is variable)
                                 - 0: false */
  int           igrari;       /* buoyant term in Rij-epsilon
                                 - 1: true (default if rho is variable)
                                 - 0: false */
  int           ikecou;       /* partially coupled version of
                                 k-epsilon (only for iturb=20)
                                 - 1: true (default)
                                 - 0: false */
  int           irijnu;       /* pseudo eddy viscosity in the matrix of momentum
                                 equation to partially implicit div( rho R )
                                 - 1: true
                                 - 0: false (default) */
  int           irijrb;       /* accurate treatment of R at the boundary (see
                                 \ref condli)
                                 - 1: true
                                 - 0: false (default) */
  int           irijec;       /* wall echo term of R
                                 - 1: true
                                 - 0: false (default) */
  int           idifre;       /* whole treatment of the diagonal part of the
                                 diffusion tensor of R and epsilon
                                 - 1: true (default)
                                 - 0: simplified treatment */
  int           iclsyr;       /* partial implicitation of symmetry BCs of R
                                 - 1: true (default)
                                 - 0: false */
  int           iclptr;       /* partial implicitation of wall BCs of R
                                 - 1: true
                                 - 0: false (default) */
  double        almax;        /* characteristic macroscopic length of the
                                 domain */
  double        uref;         /* characteristic flow velocity */
  double        xlomlg;       /* mixing length */

} cs_turb_rans_model_t;

/* les turbulence model descriptor */
/*---------------------------------*/

typedef struct {

  int           idries;       /* Van Driest smoothing at the wall (only for
                                 itytur=4)
                                 - 1: true
                                 - 0: false */
  int           ivrtex;       /* vortex method (in LES)
                                 - 1: true
                                 - 0: false (default) */

} cs_turb_les_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main turbulence model descriptor structure */

extern const cs_turb_model_t         *cs_glob_turb_model;

/* Pointer to RANS turbulence model descriptor structure */

extern const cs_turb_rans_model_t    *cs_glob_turb_rans_model;

/* Pointer to LES turbulence model descriptor structure */

extern const cs_turb_les_model_t     *cs_glob_turb_les_model;

/* Constant for turbulence models */

extern const double xkappa;
extern const double cstlog;
extern const double apow;
extern const double bpow;
extern double dpow;
extern const double cmu;
extern double cmu025;
extern const double ce1;
extern const double ce2;
extern const double ce4;
extern const double sigmak;
extern double sigmae;
extern const double crij1;
extern const double crij2;
extern const double crij3;
extern const double crijp1;
extern const double crijp2;
extern const double cssge2;
extern const double cssgs1;
extern const double cssgs2;
extern const double cssgr1;
extern const double cssgr2;
extern const double cssgr3;
extern const double cssgr4;
extern const double cssgr5;
extern const double cebms1;
extern const double cebms2;
extern const double cebmr1, cebmr2, cebmr3, cebmr4, cebmr5, cebmr6;
extern double csrij;
extern const double cebme2;
extern const double cebmmu;
extern const double xcl;
extern const double xa1;
extern const double xct;
extern const double xceta;
extern const double cpale1;
extern const double cpale2;
extern const double cpale3;
extern const double cpale4;
extern const double cpalse;
extern const double cpalmu;
extern const double cpalc1;
extern const double cpalc2;
extern const double cpalct;
extern const double cpalcl;
extern const double cpalet;
extern const double ckwsk1;
extern const double ckwsk2;
extern const double ckwsw1;
extern const double ckwsw2;
extern const double ckwbt1;
extern const double ckwbt2;
extern double ckwgm1;
extern double ckwgm2;
extern const double ckwa1;
extern const double ckwc1;
extern const double csab1;
extern const double csab2;
extern const double csasig;
extern const double csav1;
extern double csaw1;
extern const double csaw2;
extern const double csaw3;
extern const double cssr1;
extern const double cssr2;
extern const double cssr3;
extern const double ccaze2;
extern const double ccazsc;
extern const double ccaza;
extern const double ccazb;
extern const double ccazc;
extern const double ccazd;
extern const double xlesfl;
extern const double ales;
extern const double bles;
extern const double csmago;
extern const double xlesfd;
extern double smagmx;
extern const double cdries;
extern const double cv2fa1;
extern const double cv2fe2;
extern const double cv2fmu;
extern const double cv2fc1;
extern const double cv2fc2;
extern const double cv2fct;
extern const double cv2fcl;
extern const double cv2fet;
extern const double cwale;
extern const double xiafm;
extern const double etaafm;
extern const double c1trit;
extern const double c2trit;
extern const double c3trit;
extern const double c4trit;
extern const double cthafm;
extern const double cthdfm;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_MODEL_H__ */
