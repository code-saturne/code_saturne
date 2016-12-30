#ifndef __CS_CTWR_H__
#define __CS_CTWR_H__

/*============================================================================
 * Main for cooling towers related functions
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

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_ctwr_zone_t cs_ctwr_zone_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int        num;             /* Exchange zone number */
  char      *ze_name;         /* Exchange zone elements name */
  int        imctch;          /* 0: None; 1: Poppe's model; 2: Merkel's model */
  int        ct_type;         /* 1: Counter currents; 2: Crossed-currents;
                                 3: Rain zone */

  cs_real_t  hmin;            /* Minimum vertical height of exchange zone */
  cs_real_t  hmax;            /* Maximum height of exchange zone */
  cs_real_t  delta_t;         /* Temperature delta required for exchange zone
                                 if positive */
  cs_real_t  relax;           /* Relaxation of the imposed temperature */

  cs_real_t  t_l_bc;          /* Water entry temperature */
  cs_real_t  q_l_bc;          /* Water flow */
  cs_real_t  y_l_bc;          /* Mass fraction of water */

  cs_real_t  xap;             /* Exchange law lambda coefficient */
  cs_real_t  xnp;             /* Exchange law n exponent */

  cs_real_t  surface_in;      /* Water inlet surface */
  cs_real_t  surface_out;     /* Water outlet surface */
  cs_real_t  surface;         /* Total surface */

  cs_int_t   n_cells;         /* Number of air cells belonging to the zone */

  cs_int_t   up_ct_id;        /* Id of upstream exchange zone (if any) */

  cs_int_t   *ze_cell_list;   /* List of cells of ct criteria */
  cs_int_t   *mark_ze;        /* Cell marker for ct */

  cs_lnum_t n_inlet_faces;    /* Number of inlet faces */
  cs_lnum_t n_outlet_faces;   /* Number of outlet faces */
  cs_lnum_t *inlet_faces_list; /* List of inlet faces */
  cs_lnum_t *outlet_faces_list; /* List of outlet faces */

  cs_real_t  q_l_in;          /* Water entry flow */
  cs_real_t  q_l_out;         /* Water exit flow */
  cs_real_t  t_l_in;          /* Mean water entry temperature */
  cs_real_t  t_l_out;         /* Mean water exit temperature */
  cs_real_t  h_l_in;          /* Mean water entry enthalpy */
  cs_real_t  h_l_out;         /* Mean water exit enthalpy */
  cs_real_t  t_h_in;          /* Mean air entry temperature */
  cs_real_t  t_h_out;         /* Mean air exit temperature */
  cs_real_t  xair_e;          /* Mean air entry humidity */
  cs_real_t  xair_s;          /* Mean air exit humidity */
  cs_real_t  h_h_in;          /* Mean air entry enthalpy */
  cs_real_t  h_h_out;         /* Mean air exit enthalpy */
  cs_real_t  q_h_in;          /* Air entry flow */
  cs_real_t  q_h_out;         /* Air exit flow */

  cs_real_t  droplet_diam;    /* Drop diameter for rain zones */

};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* array of exchanges area */

extern cs_int_t            cs_glob_ct_nbr_max;
extern cs_int_t            cs_glob_ct_nbr;
extern cs_ctwr_zone_t     **cs_glob_ct_tab;

/* array containing the stacking of the exchange area*/
extern cs_int_t  *  cs_stack_ct;

/* array containing the treatment order of the exchanges areas */
extern cs_int_t  *  cs_chain_ct;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

void
cs_ctwr_field_pointer_map(void);

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]   zone_criterion  Zone name for selction
 * \param[in]   imctch          model : 1 Poppe
 *                                      2 Merkel
 *                                      0 None
 * \param[in]   ct_type         type  : 1 Counter current
 *                                      2 Crossed current
 *                                      3 Rain zone
 * \param[in]   delta_t         Imposed delta temperature delta between inlet
 *                              and oulet of the zone
 * \param[in]   relax           Relaxation of the imposed delta temperature
 * \param[in]   t_l_bc          Liquid water temperature at the inlet
 * \param[in]   q_l_bc          Flow at the inlet
 * \param[in]   xap             Beta_x_0 of the exchange law
 * \param[in]   xnp             Exponent n of the exchange law
 * \param[in]   surface         Total Surface of ingoing water
 * \param[in]   droplet_diam    Droplet diameter in rain zones
 */
/*----------------------------------------------------------------------------*/

void cs_ctwr_define(const char       zone_criterion[],
                    const int        imctch,
                    const int        ct_type,
                    const cs_real_t  delta_t,
                    const cs_real_t  relax,
                    const cs_real_t  t_l_bc,
                    const cs_real_t  q_l_bc,
                    const cs_real_t  xap,
                    const cs_real_t  xnp,
                    const cs_real_t  surface,
                    const cs_real_t  droplet_diam);

/*----------------------------------------------------------------------------
* Function cs_ctwr_source_term
* Phase change source terms - Exchange terms between the injected liquid
* and the water vapour phase in the bulk, humid air
*----------------------------------------------------------------------------*/
void cs_ctwr_source_term
(
   const int          f_id,             /* field id */
   const cs_real_t    p0,               /* Reference pressure */
   const cs_real_t    molmassrat,       /* dry air to water vapour molecular mass ratio */
   cs_real_t          exp_st[],         /* Explicit source term */
   cs_real_t          imp_st[]          /* Implicit source term */
 );

/*----------------------------------------------------------------------------
* Function cs_ctwr_init_flow
* Initialise the flow variables relevant to the cooling tower scalars
* inside the packing zones
*----------------------------------------------------------------------------*/
void cs_ctwr_init_flow_vars(cs_real_t liq_mass_flow[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Phase change mass source term from the evaporating liquid to the bulk, humid air
* Careful, this is different from an injection source term, which would normally
* be handled with 'cs_user_mass_source_term'
 *
 * \param[in]   iappel          Calling sequence flag
 * \param[in]   p0              Reference pressure
 * \param[in]   molmassrat      Dry air to water vapour molecular mass ratio
 * \param[in]   n_tot           Pointer to the total number of cells in the packing zones
 * \param[in]   packing_cell    Packing cell ids
 * \param[in]   mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void cs_ctwr_bulk_mass_source_term(const int       iappel,
                                   const cs_real_t p0,
                                   const cs_real_t molmassrat,
                                   int             *n_tot,
                                   cs_lnum_t       packing_cell[],
                                   cs_real_t       mass_source[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 *
 * \param[in]   mesh             associated mesh structure
 * \param[in]   mesh_quantities  mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(const cs_mesh_t              *mesh,
                  const cs_mesh_quantities_t   *mesh_quantities);

/*----------------------------------------------------------------------------
 * Destruction des structures associees aux ventilateurs
 *----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log Packing zone definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_setup(void);

/*----------------------------------------------------------------------------
* Function cs_ctwr_init_field_vars
* Initialise the field variables
*----------------------------------------------------------------------------*/
void cs_ctwr_init_field_vars
(
 const cs_real_t r0,         /* Reference density of humid air*/
 const cs_real_t t0,         /* Reference temperature of humid air*/
 const cs_real_t p0,         /* Reference pressure */
 const cs_real_t molmassrat  /* Dry air to water vapour molecular mass ratio */
 );

/*----------------------------------------------------------------------------
* Function cs_ctwr_phyvar_update
* Update the thermo physical properties fields for the humid air and the liquid
*----------------------------------------------------------------------------*/
void cs_ctwr_phyvar_update
(
 const cs_real_t r0,         /* Reference density of humid air*/
 const cs_real_t t0,         /* Reference temperature of humid air*/
 const cs_real_t p0,         /* Reference pressure */
 const cs_real_t molmassrat  /* Dry air to water vapour molecular mass ratio */
);

/*----------------------------------------------------------------------------
* Convert the injected liquid scalars from and to their transported form
*
* iflag = 1 : Convert transported variables to physical variables
* iflag = 2 : Convert physical variables to transported variables
*----------------------------------------------------------------------------*/
void cs_ctwr_transport_vars
(
 const int iflag    /* Operation flag */
);

/*----------------------------------------------------------------------------
 * Calcul des PdC induites dans les zones de pluie
 *----------------------------------------------------------------------------*/

void cs_ctwr_aetsvi
(
  const int         idim,
  const cs_real_t   rho[],       /* masse volumique air */
  const cs_real_t   vitx[],      /* vitesse air suivant x */
  const cs_real_t   vity[],      /* vitesse air suivant y */
  const cs_real_t   vitz[],      /* vitesse air suivant z */
  const cs_real_t   xair[],             /* humidite de l'air */
  cs_real_t   utsex[]            /* terme source explicite */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform balances in packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_balance(void);

/*----------------------------------------------------------------------------
 * Get pointer to exchange area.
 *
 * parameters:
 *   ct_id  <--  Id (0 to n-1) of exchange area
 *
 * returns:
 *   pointer to exchange area structure
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t *
cs_ctwr_by_id(int ct_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_H__ */
