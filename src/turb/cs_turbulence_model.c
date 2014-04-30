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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_model.c
        Base turbulence model data.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_model_t

  \brief Turbulence model general options descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_model_t::iturb
        turbulence model
        - 0: no turbulence model (laminar flow)
        - 10: mixing length model
        - 20: standard \f$ k-\varepsilon \f$ model
        - 21: \f$ k-\varepsilon \f$ model with Linear Production (LP) correction
        - 30: \f$ R_{ij}-\epsilon \f$ (LRR)
        - 31: \f$ R_{ij}-\epsilon \f$ (SSG)
        - 32: \f$ R_{ij}-\epsilon \f$ (EBRSM)
        - 40: LES (constant Smagorinsky model)
        - 41: LES ("classical" dynamic Smagorisky model)
        - 42: LES (WALE)
        - 50: v2f phi-model
        - 51: v2f \f$ BL-v^2-k \f$
        - 60: \f$ k-\omega \f$ SST
        - 70: Spalart-Allmaras model
  \var  cs_turb_model_t::itytur
        class of turbulence model (integer value iturb/10)
  \var  cs_turb_model_t::nvarcl
        number of variable plus number of turbulent fluxes (used by the
        boundary conditions)
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_rans_model_t

  \brief RANS turbulence model descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_rans_model_t::irccor
        activation of rotation/curvature correction for an eddy viscosity
        turbulence models
        - 0: false
        - 1: true
  \var  cs_turb_rans_model_t::itycor
        type of rotation/curvature correction for an eddy viscosity turbulence
        models
        - 1: Cazalbou correction (default when irccor=1 and itytur=2 or 5)
        - 2: Spalart-Shur correction (default when irccor=1 and iturb=60 or 70)
  \var  cs_turb_rans_model_t::idirsm
        turbulent diffusion model for second moment closure
        - 0: scalar diffusivity (Shir model)
        - 1: tensorial diffusivity (Daly and Harlow model, default model)
  \var  cs_turb_rans_model_t::iclkep
        clipping of k and epsilon
        - 0: absolute value clipping
        - 1: coupled clipping based on physical relationships
  \var  cs_turb_rans_model_t::igrhok
        take \f$ 2/3 \rho \grad k \f$ in the momentum equation
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::igrake
        buoyant term in \f$ k- \varepsilon \f$
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
  \var  cs_turb_rans_model_t::igrari
        buoyant term in \f$ R_{ij}- \varepsilon \f$
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
  \var  cs_turb_rans_model_t::ikecou
        partially coupled version of \f$ k-\varepsilon \f$ (only for iturb=20)
        - 1: true (default)
        - 0: false
  \var  cs_turb_rans_model_t::irijnu
        pseudo eddy viscosity in the matrix of momentum equation to partially
        implicit \f$ \divv \left( \rho \tens{R} \right) \f$
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijrb
        accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijec
        wall echo term of \f$ \tens{R} \f$
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::idifre
        whole treatment of the diagonal part of the diffusion tensor of \f$
        \tens{R} \f$ and \f$ \varepsilon \f$
        - 1: true (default)
        - 0: simplified treatment
  \var  cs_turb_rans_model_t::iclsyr
        partial implicitation of symmetry BCs of \f$ \tens{R} \f$
        - 1: true (default)
        - 0: false
  \var  cs_turb_rans_model_t::iclptr
        partial implicitation of wall BCs of \f$ \tens{R} \f$
        - 1: true
        - 0: false (default)
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_les_model_t

  \brief LES turbulence model descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_les_model_t::idries
        Van Driest smoothing at the wall (only for itytur=4)
        - 1: true
        - 0: false
  \var  cs_turb_les_model_t::ivrtex
        vortex method (in LES)
        - 1: true
        - 0: false (default)
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* main turbulence model structure and associated pointer */

static cs_turb_model_t  _turb_model =
  {-999, -999, 0};

const cs_turb_model_t  *cs_glob_turb_model = &_turb_model;

/* RANS turbulence model structure and associated pointer */

static cs_turb_rans_model_t  _turb_rans_model =
  {0, -999, 0, 0, 1, 1, -999, 0, 0, 0, 1, 1, 0};

const cs_turb_rans_model_t  *cs_glob_turb_rans_model = &_turb_rans_model;

/* LES turbulence model structure and associated pointer */

static cs_turb_les_model_t  _turb_les_model =
  {-1, 0};

const cs_turb_les_model_t  *cs_glob_turb_les_model = &_turb_les_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **nvarcl);

void
cs_f_turb_rans_model_get_pointers(int     **irccor,
                                  int     **itycor,
                                  int     **idirsm,
                                  int     **iclkep,
                                  int     **igrhok,
                                  int     **igrake,
                                  int     **igrari,
                                  int     **ikecou,
                                  int     **irijnu,
                                  int     **irijrb,
                                  int     **irijec,
                                  int     **idifre,
                                  int     **iclsyr,
                                  int     **iclptr);

void
cs_f_turb_les_model_get_pointers(int     **idries,
                                 int     **ivrtex);

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to members of the global turbulence model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iturb  --> pointer to cs_glob_turb_model->iturb
 *   itytur --> pointer to cs_glob_turb_model->itytur
 *   nvarcl --> pointer to cs_glob_turb_model->nvarcl
 *----------------------------------------------------------------------------*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **nvarcl)
{
  *iturb  = &(_turb_model.iturb);
  *itytur = &(_turb_model.itytur);
  *nvarcl = &(_turb_model.nvarcl);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the RANS turbulence functions structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   irccor --> pointer to cs_glob_turb_rans_model->irccor
 *   itycor --> pointer to cs_glob_turb_rans_model->itycor
 *   idirsm --> pointer to cs_glob_turb_rans_model->idirsm
 *   iclkep --> pointer to cs_glob_turb_rans_model->iclkep
 *   igrhok --> pointer to cs_glob_turb_rans_model->igrhok
 *   igrake --> pointer to cs_glob_turb_rans_model->igrake
 *   igrari --> pointer to cs_glob_turb_rans_model->igrari
 *   ikecou --> pointer to cs_glob_turb_rans_model->ikecou
 *   irijnu --> pointer to cs_glob_turb_rans_model->irijnu
 *   irijrb --> pointer to cs_glob_turb_rans_model->irijrb
 *   irijec --> pointer to cs_glob_turb_rans_model->irijec
 *   idifre --> pointer to cs_glob_turb_rans_model->idifre
 *   iclsyr --> pointer to cs_glob_turb_rans_model->iclsyr
 *   iclptr --> pointer to cs_glob_turb_rans_model->iclptr
 *----------------------------------------------------------------------------*/

void
cs_f_turb_rans_model_get_pointers(int     **irccor,
                                  int     **itycor,
                                  int     **idirsm,
                                  int     **iclkep,
                                  int     **igrhok,
                                  int     **igrake,
                                  int     **igrari,
                                  int     **ikecou,
                                  int     **irijnu,
                                  int     **irijrb,
                                  int     **irijec,
                                  int     **idifre,
                                  int     **iclsyr,
                                  int     **iclptr)
{
  *irccor = &(_turb_rans_model.irccor);
  *itycor = &(_turb_rans_model.itycor);
  *idirsm = &(_turb_rans_model.idirsm);
  *iclkep = &(_turb_rans_model.iclkep);
  *igrhok = &(_turb_rans_model.igrhok);
  *igrake = &(_turb_rans_model.igrake);
  *igrari = &(_turb_rans_model.igrari);
  *ikecou = &(_turb_rans_model.ikecou);
  *irijnu = &(_turb_rans_model.irijnu);
  *irijrb = &(_turb_rans_model.irijrb);
  *irijec = &(_turb_rans_model.irijec);
  *idifre = &(_turb_rans_model.idifre);
  *iclsyr = &(_turb_rans_model.iclsyr);
  *iclptr = &(_turb_rans_model.iclptr);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the LES turbulence model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   idries --> pointer to cs_glob_turb_les_model->idries
 *   ivrtex --> pointer to cs_glob_turb_les_model->ivrtex
 *----------------------------------------------------------------------------*/

void
cs_f_turb_les_model_get_pointers(int     **idries,
                                 int     **ivrtex)
{
  *idries = &(_turb_les_model.idries);
  *ivrtex = &(_turb_les_model.ivrtex);
}

/*! \endcond (end ignore by Doxygen) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
