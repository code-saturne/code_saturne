/*============================================================================
 * Stokes equation model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_stokes_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_stokes_model.c
        Stokes equation model data.

  \struct cs_stokes_model_t

  \brief Stokes equation model descriptor.

  Members of these Stokes equation model descriptor are publicly accessible, to
  allow for concise syntax, as it is expected to be used in many places.

  \var  cs_stokes_model_t::ivisse
        take \f$ \divs \left( \mu \transpose{\gradt \, \vect{u}} - 2/3 \mu
        \trace{\gradt \, \vect{u}} \right) \f$
        into account in the momentum equation
  \var  cs_stokes_model_t::irevmc
        reconstruction of the velocity field with the updated pressure option
        - 0: default
  \var  cs_stokes_model_t::iprco
        compute the pressure step thanks to the continuity equation
        - 1: true (default)
        - 0: false
  \var  cs_stokes_model_t::irnpnw
        compute the normed residual for the pressure step in the prediction step
        - 1: true (default)
        - 0: false
  \var  cs_stokes_model_t::rnormp
        normed residual for the pressure step
  \var  cs_stokes_model_t::arak
        Arakawa multiplicator for the Rhie and Chow filter (1 by default)
  \var  cs_stokes_model_t::ipucou
        pseudo coupled pressure-velocity solver
        - 1: true (default)
        - 0: false
  \var  cs_stokes_model_t::iccvfg
        calculation with a fixed velocity field
        - 1: true
        - 0: false (default)
  \var  cs_stokes_model_t::idilat
        algorithm to take into account the density variation in time
        - 1: dilatable steady algorithm (default)
        - 2: dilatable unsteady algorithm
        - 3: low-Mach algorithm
        - 4: algorithm for fire
  \var  cs_stokes_model_t::epsdp
        parameter of diagonal pressure strengthening
  \var  cs_stokes_model_t::itbrrb
        accurate treatment of the wall temperature
        (reconstruction of wall temperature)
        - 1: true
        - 0: false (default)
        (see \ref condli, useful in case of coupling with syrthes)
  \var  cs_stokes_model_t::iphydr
        improve static pressure algorithm
        Take into account the balance or imbalance between the pressure
        gradient and source terms (as gravity and head losses)
        - 1: impose the equilibrium of the hydrostaic part of the pressure with
          any external force, even head losses
        - 0: no treatment (default)
        - 2: hydrostatic pressure computation with a apriori momentum equation
             to obtain a hydrostatic pressure taking into account the imbalance
             between the pressure gradient and the gravity source term
  \var  cs_stokes_model_t::igprij
        improve static pressure algorithm
        - 1: take -div(rho R) in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)
  \var  cs_stokes_model_t::igpust
        improve static pressure algorithm
        - 1: take user momentum source terms in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)
  \var  cs_stokes_model_t::iifren
        indicates the presence of a Bernoulli boundary face (automatically
        computed)
        - 0: no face
        - 1: at least one face
  \var  cs_stokes_model_t::icalhy
        compute the hydrostatic pressure in order to compute the Dirichlet
        conditions on the pressure at outlets
        - 1: true
        - 0: false (default)
  \var  cs_stokes_model_t::irecmf
        use interpolated face diffusion coefficient instead of cell diffusion coefficient
        for the mass flux reconstruction for the non-orthogonalities
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

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main Stokes equation model descriptor structure and associated pointer:
 * Default Options */

static cs_stokes_model_t  _stokes_model = {
  .ivisse = 1,
  .irevmc = 0,
  .iprco  = 1,
  .irnpnw = 1,
  .rnormp = 0,
  .arak   = 1.0,
  .ipucou = 0,
  .iccvfg = 0,
  .idilat = 1,
  .epsdp  = 1.e-12,
  .itbrrb = 0,
  .iphydr = 0,
  .igprij = 0,
  .igpust = 1,
  .iifren = 0,
  .icalhy = -1,
  .irecmf = 0};

const cs_stokes_model_t  *cs_glob_stokes_model = &_stokes_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_stokes_options_get_pointers(int     **ivisse,
                                 int     **irevmc,
                                 int     **iprco,
                                 int     **irnpnw,
                                 double  **rnormp,
                                 double  **arak,
                                 int     **ipucou,
                                 int     **iccvfg,
                                 int     **idilat,
                                 double  **epsdp,
                                 int     **itbrrb,
                                 int     **iphydr,
                                 int     **igprij,
                                 int     **igpust,
                                 int     **iifren,
                                 int     **icalhy,
                                 int     **irecmf);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global Stokes model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ivisse  --> pointer to cs_glob_stokes_model->ivisse
 *   irevmc  --> pointer to cs_glob_stokes_model->irevmc
 *   iprco   --> pointer to cs_glob_stokes_model->iprco
 *   irnpnw  --> pointer to cs_glob_stokes_model->irnpnw
 *   rnormp  --> pointer to cs_glob_stokes_model->rnormp
 *   arak    --> pointer to cs_glob_stokes_model->arak
 *   ipucou  --> pointer to cs_glob_stokes_model->ipucou
 *   iccvfg  --> pointer to cs_glob_stokes_model->iccvfg
 *   idilat  --> pointer to cs_glob_stokes_model->idilat
 *   epsdp   --> pointer to cs_glob_stokes_model->epsdp
 *   itbrrb  --> pointer to cs_glob_stokes_model->itbrrb
 *   iphydr  --> pointer to cs_glob_stokes_model->iphydr
 *   igprij  --> pointer to cs_glob_stokes_model->igprij
 *   igpust  --> pointer to cs_glob_stokes_model->igpust
 *   iifren  --> pointer to cs_glob_stokes_model->iifren
 *   icalhy  --> pointer to cs_glob_stokes_model->icalhy
 *   irecmf  --> pointer to cs_glob_stokes_model->irecmf
 *----------------------------------------------------------------------------*/

void
cs_f_stokes_options_get_pointers(int     **ivisse,
                                 int     **irevmc,
                                 int     **iprco,
                                 int     **irnpnw,
                                 double  **rnormp,
                                 double  **arak,
                                 int     **ipucou,
                                 int     **iccvfg,
                                 int     **idilat,
                                 double  **epsdp,
                                 int     **itbrrb,
                                 int     **iphydr,
                                 int     **igprij,
                                 int     **igpust,
                                 int     **iifren,
                                 int     **icalhy,
                                 int     **irecmf)
{
  *ivisse = &(_stokes_model.ivisse);
  *irevmc = &(_stokes_model.irevmc);
  *iprco  = &(_stokes_model.iprco);
  *irnpnw = &(_stokes_model.irnpnw);
  *rnormp = &(_stokes_model.rnormp);
  *arak   = &(_stokes_model.arak);
  *ipucou = &(_stokes_model.ipucou);
  *iccvfg = &(_stokes_model.iccvfg);
  *idilat = &(_stokes_model.idilat);
  *epsdp  = &(_stokes_model.epsdp );
  *itbrrb = &(_stokes_model.itbrrb);
  *iphydr = &(_stokes_model.iphydr);
  *igprij = &(_stokes_model.igprij);
  *igpust = &(_stokes_model.igpust);
  *iifren = &(_stokes_model.iifren);
  *icalhy = &(_stokes_model.icalhy);
  *irecmf = &(_stokes_model.irecmf);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_stokes_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_stokes_model_t *
cs_get_glob_stokes_model(void)
{
  return &_stokes_model;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
