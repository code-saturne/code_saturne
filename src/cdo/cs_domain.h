#ifndef __CS_DOMAIN_H__
#define __CS_DOMAIN_H__

/*============================================================================
 * Manage a computational domain
 *  - Physical boundary conditions attached to a domain
 *  - Equations
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_time_step.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Store information about the kind of boundary attached to the computational
   domain: inlet, outlet, wall, symmetry... */

typedef struct _cs_domain_boundary_t cs_domain_boundary_t;

typedef struct {

  /* Code_Saturne mesh and mesh quantities structures already computed */
  const  cs_mesh_t              *mesh;
  const  cs_mesh_quantities_t   *mesh_quantities;

  /* CDO structures:
     - cs_cdo_connect_t contains additional information about connectivity
     - cs_cdo_quantities_t contains additional information on mesh quantities
  */

  cs_cdo_connect_t              *connect;
  cs_cdo_quantities_t           *cdo_quantities;

  /* Physical boundary conditions on the computational domain */
  cs_domain_boundary_t          *boundaries;

  /* Time step management */
  double                   dt_cur;             // current time step
  cs_param_def_type_t      time_step_def_type; // Way of defining the time step
  cs_def_t                 time_step_def;      // Definition of the time_step
  cs_time_step_t          *time_step;          // time step descriptor
  cs_time_step_options_t   time_options;       // time step options

  /* Overview of pre-defined equations to solve
      - Navier-Stokes equations (named NavierStokes)
      - Wall distance (named WallDistance)
  */
  bool             do_navsto;

  // TODO: add a specific equation for solving Navier-Stokes

  /* Number of equations defined on this domain splitted into
     predefined equations and user equations.
     Predefined equations are stored first.
  */
  int              n_equations;
  int              n_predef_equations;
  int              n_user_equations;
  cs_equation_t  **equations;

  bool             only_steady;

  /* Predefined equations
     If xxxxx_eq_id = -1, then this equation is not activated */

  int              wall_distance_eq_id;

  /* Output options */
  int              output_freq;

} cs_domain_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize of cs_domain_t structure
 *
 * \param[in]   mesh              pointer to a cs_mesh_t struct.
 * \param[in]   mesh_quantities   pointer to a cs_mesh_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_init(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Proceed to the last settings of a cs_domain_t structure
 *
 * \param[in, out]  domain    pointer to the cs_domain_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_last_setup(cs_domain_t    *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_domain_t structure
 *
 * \param[in, out]   domain    pointer to the cs_domain_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_free(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_domain_t structure
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_summary(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary type by default
 *
 * \param[in, out]   domain        pointer to a cs_domain_t structure
 * \param[in]        bdy_name      key name of the default boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_default_boundary(cs_domain_t       *domain,
                               const char        *bdy_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a boundary type defined on a mesh location
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        ml_name      mesh location name
 * \param[in]        bdy_name     key name of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_boundary(cs_domain_t               *domain,
                       const char                *ml_name,
                       const char                *bdy_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the frequency at which output is done in listing
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        freq     each freq iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_output_freq(cs_domain_t   *domain,
                          int            freq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the time step structure related to a domain
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        t_end     final physical time
 * \param[in]        nt_max    max. number of temporal iterations
 * \param[in]        defkey    way of defining the time step
 * \param[in]        defval    definition of the time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_time_step(cs_domain_t   *domain,
                        double         t_end,
                        int            nt_max,
                        const char    *defkey,
                        void          *defval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure whith name eqname
 *         Return NULL if not find
 *
 * \param[in]  domain    pointer to a cs_domain_t structure
 * \param[in]  eqname    name of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_domain_get_equation(const cs_domain_t  *domain,
                       const char         *eqname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the computation of the wall distance
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_activate_wall_distance(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup predefined equations which are activated
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_predefined_equations(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user equation to a domain
 *
 * \param[in, out] domain         pointer to a cs_domain_t structure
 * \param[in]      eqname         name of the equation
 * \param[in]      varname        name of the related variable
 * \param[in]      key_type       type of equation: "scalar", "vector", "tensor"
 * \param[in]      is_steady      add an unsteady term or not
 * \param[in]      do_convection  add a convection term or not
 * \param[in]      do_diffusion   add a diffusion term or not
 * \param[in]      key_bc         type of boundary condition set by default
 *                                "zero_value" or "zero_flux"
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_user_equation(cs_domain_t         *domain,
                            const char          *eqname,
                            const char          *varname,
                            const char          *key_type,
                            bool                 is_steady,
                            bool                 do_convection,
                            bool                 do_diffusion,
                            const char          *key_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_field_t structure for each equation defined in the
 *         domain
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_create_fields(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to continue iterations in time
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_iterate(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_solve(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_H__ */
