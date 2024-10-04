#ifndef __CS_CTWR_H__
#define __CS_CTWR_H__

/*============================================================================
 * Cooling towers related functions
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Type of cooling tower model */

typedef enum {

  CS_CTWR_NONE = 0,              /*!< no cooling tower model */
  CS_CTWR_POPPE = 1,             /*!< Poppe's model */
  CS_CTWR_MERKEL = 2             /*!< Merkel's model */

} cs_ctwr_model_t;

/*! Type of cooling tower exchange zone */

typedef enum {

  CS_CTWR_COUNTER_CURRENT = 1,   /*!< counter-current zone */
  CS_CTWR_CROSS_CURRENT = 2,     /*!< cross-current zone */
  CS_CTWR_INJECTION = 3,         /*!< water injection zone */

} cs_ctwr_zone_type_t;

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int                  num;        /* Exchange zone number */
  char                *criteria;   /* Exchange zone selection criteria */
  int                  z_id;       /* id of the volume zone */
  char                *name;       /* Exchange zone name */
  char                *file_name;  /* Exchange zone budget file name */
  cs_ctwr_zone_type_t  type;       /* Zone type */

  cs_real_t  hmin;               /* Minimum vertical height of exchange zone */
  cs_real_t  hmax;               /* Maximum height of exchange zone */
  cs_real_t  delta_t;            /* Temperature delta required for exchange zone
                                    if positive */
  cs_real_t  relax;              /* Relaxation of the imposed temperature */

  cs_real_t  t_l_bc;             /* Water entry temperature */
  cs_real_t  q_l_bc;             /* Water flow */

  cs_real_t  xap;                /* Exchange law a_0 coefficient */
  cs_real_t  xnp;                /* Exchange law n exponent */

  cs_real_t  surface_in;         /* Water inlet surface */
  cs_real_t  surface_out;        /* Water outlet surface */
  cs_real_t  surface;            /* Total surface */

  cs_real_t  xleak_fac;          /* Leakage factor (ratio of outlet/inlet
                                    flow rate) */
  cs_real_t  v_liq_pack;         /* Vertical liquid film velocity in packing */

  cs_lnum_t  n_cells;            /* Number of air cells belonging to the zone */
  cs_real_t  vol_f;              /* Cooling tower zone total volume */

  int        up_ct_id;           /* Id of upstream exchange zone (if any) */

  cs_lnum_t  n_inlet_faces;      /* Number of inlet faces */
  cs_lnum_t  n_outlet_faces;     /* Number of outlet faces */
  cs_lnum_t *inlet_faces_ids;    /* List of inlet faces */
  cs_lnum_t *outlet_faces_ids;   /* List of outlet faces */

  cs_lnum_t  n_outlet_cells;     /* Number of outlet cells */
  cs_lnum_t *outlet_cells_ids;   /* List of outlet cells */

  cs_real_t  p_in;            /* Average inlet pressure */
  cs_real_t  p_out;           /* Average outlet pressure */
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

};

typedef struct _cs_ctwr_zone_t cs_ctwr_zone_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------
 * Cooling Tower model options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_ctwr_model_t  evap_model;
  bool             has_rain;
  bool             mixture_model; /*!< Activate multiphase mixture model
                                       for rain + air (not functioning yet) >!*/
  bool             solve_rain_velocity; /*!< Activate drift velocity
                                          resolution */
  bool             air_rain_friction; /*!< Activate interfacial friction */
  bool             rain_to_packing; /*!< Activate liquid water transfer
                                         from rain to packing */
} cs_ctwr_option_t;


/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to cooling tower model options structure */
extern const cs_ctwr_option_t  *cs_glob_ctwr_option;

/* Make number of cooling towers zones accessible */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_ctwr_option
 *----------------------------------------------------------------------------*/

cs_ctwr_option_t *
cs_get_glob_ctwr_option(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_ctwr_zone
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t **
cs_get_glob_ctwr_zone(void);

/*----------------------------------------------------------------------------
 * Provide access to number of ct zones
 *----------------------------------------------------------------------------*/

int *
cs_get_glob_ctwr_n_zones(void);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]  zone_criteria  zone selection criteria
 * \param[in]  z_id           z_id if zone already created (-1 otherwise)
 * \param[in]  zone_type      exchange zone type
 * \param[in]  delta_t        imposed delta temperature delta between inlet
 *                            and oulet of the zone
 * \param[in]  relax          relaxation of the imposed delta temperature
 * \param[in]  t_l_bc         liquid water temperature at the inlet
 * \param[in]  q_l_bc         mass flow rate at the inlet
 * \param[in]  xap            beta_x_0 of the exchange law
 * \param[in]  xnp            exponent n of the exchange law
 * \param[in]  surface        total Surface of ingoing water
 * \param[in]  xleak_fac      leakage factor (ratio of outlet/inlet flow rate)
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define(const char           zone_criteria[],
               int                  z_id,
               cs_ctwr_zone_type_t  zone_type,
               cs_real_t            delta_t,
               cs_real_t            relax,
               cs_real_t            t_l_bc,
               cs_real_t            q_l_bc,
               cs_real_t            xap,
               cs_real_t            xnp,
               cs_real_t            surface,
               cs_real_t            xleak_fac);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fields used by the cooling tower module to pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_field_pointer_map(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy cs_ctwr_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log Packing zone definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform balances in packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_balance(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Convert injected liquid scalars from and to their transported form.
 *
 * \param[in]   iflag     1: Convert transported variables to physical variables
 *                        2: Convert physical variables to
 *                           transported variables
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_transport_vars(int  iflag);

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
