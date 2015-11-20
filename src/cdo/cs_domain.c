/*============================================================================
 * Manage a computational domain within the CDO framework
 *  - Physical boundary conditions attached to a domain
 *============================================================================*/

/* VERS */

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
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_evaluate.h"
#include "cs_walldistance.h"

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

typedef struct {

  int     id;          // mesh location id (-1 if not used)
  int     n_sub_ids;   // number of related (sub) mesh location attached
  int    *sub_ids;     // list of mesh location ids

} cs_domain_bdy_ml_t;

struct _cs_domain_boundary_t {

  cs_param_boundary_type_t    default_boundary; // boundary set by default

  cs_lnum_t                   n_b_faces;        // number of boundary faces
  cs_param_boundary_type_t   *types;            // type of each boundary face

  /* Number of border faces related to each type of boundary */
  cs_lnum_t                   n_elts[CS_PARAM_N_BOUNDARY_TYPES];

  /* Id related to specific mesh location structures corresponding to
     different types of boundaries
     >> "domain_walls"
     >> "domain_inlets"
     >> "domain_outlets"
     >> "domain_symmetries"
  */
  cs_domain_bdy_ml_t           bdy_ml[CS_PARAM_N_BOUNDARY_TYPES];

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

static const char
_domain_boundary_ml_name[CS_PARAM_N_BOUNDARY_TYPES][CS_CDO_LEN_NAME] =
  { N_("domain_walls"),
    N_("domain_inlets"),
    N_("domain_outlets"),
    N_("domain_symmetries") };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Postprocessing of quantities attached to the current domain s.t.
 *         material properties and/or advection fields
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  time_step  pointer to a time_step structure
 */
/*----------------------------------------------------------------------------*/

static void
_domain_post(const cs_mesh_t       *mesh,
             const cs_time_step_t  *time_step)
{
  cs_lnum_t  v_id, id;

  const double  t_cur = time_step->t_cur;
  const int  nt_cur = time_step->nt_cur;
  const int  n_adv_fields = cs_param_get_n_adv_fields();

  /* Up to now, only advection field are post-processed */
  for (id = 0; id < n_adv_fields; id++) {

    cs_param_adv_field_t  *adv = cs_param_adv_field_get(id);

    if (adv->post_freq > -1) { // Post-processing is requested

      bool  do_post = true;

      if (nt_cur == 0) {
        if (adv->post_freq > 0)
          do_post = false;
      }
      else { /* nt_cur > 0 */
        if (adv->post_freq == 0)
          do_post = false;
        else if (nt_cur % adv->post_freq > 0)
          do_post = false;
      }

      if (do_post) {

        cs_lnum_t  shft;
        cs_nvec3_t  v;
        cs_real_3_t  advect;

        cs_field_t  *fld = cs_field_by_id(adv->field_id);
        cs_real_t  *val = fld->val;
        cs_real_t  *unitv = NULL, *magnitude = NULL;

        const  cs_real_t  *xyz = mesh->vtx_coord;

        bft_printf(" <post/advection_field> %s\n", adv->name);

        BFT_MALLOC(unitv, 3*mesh->n_vertices, cs_real_t);
        BFT_MALLOC(magnitude, mesh->n_vertices, cs_real_t);

        /* Evaluate the value of the advection field at each vertex */
        for (v_id = 0; v_id < mesh->n_vertices; v_id++) {

          cs_evaluate_adv_field(id,
                                t_cur,
                                xyz + 3*v_id,
                                &advect);
          cs_nvec3(advect, &v);
          magnitude[v_id] = v.meas;

          shft = 3*v_id;
          for (int k = 0; k < 3; k++) {
            val[shft+k] = advect[k];
            unitv[shft+k] = v.unitv[k];
          }

        } // Loop on vertices

        cs_post_write_vertex_var(-1,              // id du maillage de post
                                 fld->name,
                                 3,               // dim
                                 true,            // interlace
                                 true,            // true = original mesh
                                 CS_POST_TYPE_cs_real_t,
                                 val,             // values on vertices
                                 time_step);      // time step structure

        char  *label = NULL;
        int  len = strlen(fld->name) + 1 + 10;

        BFT_MALLOC(label, len, char);
        sprintf(label, "%s.Magnitude", fld->name);

        cs_post_write_vertex_var(-1,              // id du maillage de post
                                 label,
                                 1,               // dim
                                 true,            // interlace
                                 true,            // true = original mesh
                                 CS_POST_TYPE_cs_real_t,
                                 magnitude,       // values on vertices
                                 time_step);      // time step structure

        sprintf(label, "%s.Unit", fld->name);
        cs_post_write_vertex_var(-1,              // id du maillage de post
                                 label,
                                 3,               // dim
                                 true,            // interlace
                                 true,            // true = original mesh
                                 CS_POST_TYPE_cs_real_t,
                                 unitv,           // values on vertices
                                 time_step);      // time step structure

        /* Free temporay buffers */
        BFT_FREE(label);
        BFT_FREE(unitv);
        BFT_FREE(magnitude);

      } // Do post

    } // Post is potentially requested

  } // Loop on advection fields

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the setup of boundary is reliable
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_check_boundary_setup(cs_domain_t   *domain)
{
  int  error_count = 0;

  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_domain_t structure is not allocated."));

  cs_domain_boundary_t  *bcs = domain->boundaries;

  for (unsigned int k = 0; k < CS_PARAM_N_BOUNDARY_TYPES; k++)
    bcs->n_elts[k] = 0;

  /* Sanity check */
  assert(bcs !=NULL);

  for (int i = 0; i < bcs->n_b_faces; i++) {
    if (bcs->types[i] == CS_PARAM_N_BOUNDARY_TYPES)
      error_count++;
    else
      bcs->n_elts[bcs->types[i]] += 1;
  }

  if (error_count > 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Problem detected during the setup.\n"
                " %d boundary faces have no boundary type."), error_count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the setup is reliable
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_mesh_locations(cs_domain_t   *domain)
{
  cs_domain_boundary_t  *bcs = domain->boundaries;
  cs_param_boundary_type_t  default_type = bcs->default_boundary;

  /* Store all mesh locations not associated to the default type */
  int  n_sub_ids = 0;
  int  *sub_ids = NULL;

  for (unsigned int type = 0; type < CS_PARAM_N_BOUNDARY_TYPES; type++) {

    int _n_sub_ids = bcs->bdy_ml[type].n_sub_ids;
    int  *_sub_ids = bcs->bdy_ml[type].sub_ids;

    if (_n_sub_ids > 0 && type != default_type) {
      // There is at least one mesh location

      bcs->bdy_ml[type].id  =
        cs_mesh_location_add_by_union(_domain_boundary_ml_name[type],
                                      CS_MESH_LOCATION_BOUNDARY_FACES,
                                      _n_sub_ids,
                                      _sub_ids,
                                      false);  // complement ?

      /* Incremant the list of mesh locations */
      BFT_REALLOC(sub_ids, n_sub_ids + _n_sub_ids, int);
      for (int i = 0; i < _n_sub_ids; i++)
        sub_ids[n_sub_ids + i] = _sub_ids[i];
      n_sub_ids += _n_sub_ids;

    }

  } /* Loop on boundary types */

  /* Treatment of the default boundary type */
  bcs->bdy_ml[default_type].id  =
    cs_mesh_location_add_by_union(_domain_boundary_ml_name[default_type],
                                  CS_MESH_LOCATION_BOUNDARY_FACES,
                                  n_sub_ids,
                                  sub_ids,
                                  true);

  BFT_FREE(sub_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_domain_boundary_t structure
 *
 * \param[in]   n_b_faces   number of boundary faces
 *
 * \return a pointer to a new allocated cs_domain_boundary_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_boundary_t *
_init_domain_boundaries(cs_lnum_t     n_b_faces)
{
  cs_domain_boundary_t  *bcs = NULL;

  BFT_MALLOC(bcs, 1, cs_domain_boundary_t);

  bcs->n_b_faces = n_b_faces;
  bcs->default_boundary = CS_PARAM_N_BOUNDARY_TYPES;

  BFT_MALLOC(bcs->types, n_b_faces, cs_param_boundary_type_t);
  for (int i = 0; i < n_b_faces; i++)
    bcs->types[i] = CS_PARAM_N_BOUNDARY_TYPES;

  for (int i = 0; i < CS_PARAM_N_BOUNDARY_TYPES; i++) {
    bcs->bdy_ml[i].id = -1;
    bcs->bdy_ml[i].n_sub_ids = 0;
    bcs->bdy_ml[i].sub_ids = NULL;
    bcs->n_elts[i] = 0;
  }

  return bcs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_domain_boundary_t structure
 *
 * \param[in,out]  bcs       pointer to the cs_domain_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static cs_domain_boundary_t *
_free_domain_boundaries(cs_domain_boundary_t   *bcs)
{
  if (bcs == NULL)
    return bcs;

  BFT_FREE(bcs->types);

  for (int type = 0; type < CS_PARAM_N_BOUNDARY_TYPES; type++)
    BFT_FREE(bcs->bdy_ml[type].sub_ids);

  BFT_FREE(bcs);

  return NULL;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a cs_domain_t structure
 *
 * \param[in]   mesh              pointer to a cs_mesh_t struct.
 * \param[in]   mesh_quantities   pointer to a cs_mesh_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_init(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_real_t  default_time_step = -1e13;
  cs_domain_t  *domain = NULL;

  BFT_MALLOC(domain, 1, cs_domain_t);

  domain->mesh = mesh;
  domain->mesh_quantities = mesh_quantities;

  /* Build additional connectivity structures */
  domain->connect = cs_cdo_connect_build(mesh);

  /* Build additional mesh quantities in a seperate structure */
  domain->cdo_quantities =  cs_cdo_quantities_build(mesh,
                                                    mesh_quantities,
                                                    domain->connect);

  /* Default initialization of the time step */
  domain->dt_cur = default_time_step;
  domain->time_step_def_type = CS_PARAM_DEF_BY_VALUE;
  domain->time_step_def.get.val = default_time_step;

  BFT_MALLOC(domain->time_step, 1, cs_time_step_t);
  domain->time_step->is_variable = 0;
  domain->time_step->is_local = 0; // In CDO, this is always equal to 0
  domain->time_step->nt_prev = 0; // Changed in case of restart
  domain->time_step->nt_cur =  0; // Do not modify this value
  domain->time_step->nt_max = 0;
  domain->time_step->nt_ini = 2;  // Not useful in CDO module
  domain->time_step->t_prev = 0.;
  domain->time_step->t_cur = 0.;  // Assume initial time is 0.0
  domain->time_step->t_max = -1.;

  domain->time_options.inpdt0 = 0; // standard calculation
  domain->time_options.iptlro = 0;
  domain->time_options.idtvar = 0; // constant time step by default
  domain->time_options.dtref = default_time_step;
  domain->time_options.coumax = 1.;
  domain->time_options.cflmmx = 0.99;
  domain->time_options.foumax = 10.;
  domain->time_options.varrdt = 0.1;
  domain->time_options.dtmin = default_time_step;
  domain->time_options.dtmax = default_time_step;
  domain->time_options.relxst = 0.7; // Not useful in CDO schemes

  /* Equations */
  domain->do_navsto = false;
  domain->n_equations = 0;
  domain->n_predef_equations = 0;
  domain->n_user_equations = 0;
  domain->equations = NULL;
  domain->only_steady = true;

  /* Other options */
  domain->output_freq = 10;

  /* Predefined equations */
  domain->wall_distance_eq_id = -1;

  /* Specify the "physical" domain boundaries. Define a mesh location for
     each boundary type */
  domain->boundaries = _init_domain_boundaries(mesh->n_b_faces);

  /* User-defined settings for this domain */
  cs_user_cdo_setup_domain(domain);

  /* Update mesh locations */
  _add_mesh_locations(domain);
  _check_boundary_setup(domain);

  /* Sanity check */
  assert(domain->n_equations ==
         domain->n_user_equations + domain->n_predef_equations);

  return domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Proceed to the last settings of a cs_domain_t structure
 *
 * \param[in, out]  domain    pointer to the cs_domain_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_last_setup(cs_domain_t    *domain)
{
  /* Proceed to the last settings of a cs_equation_t structure
     - Assign to a cs_equation_t structure a list of function to manage this
       structure during the computation.
     - The set of functions chosen for each equation depends on the parameters
       specifying the cs_equation_t structure
     - Setup the structure related to cs_sles_*
  */

  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

    cs_equation_t  *eq = domain->equations[eq_id];

    cs_equation_last_setup(eq);
    if (!cs_equation_is_steady(eq))
      domain->only_steady = false;

  } // Loop on equations

}

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
cs_domain_free(cs_domain_t   *domain)
{
  if (domain == NULL)
    return domain;

  /* cs_mesh_t and cs_mesh_quantities_t structure are not freed since they
     are only shared */
  domain->mesh = NULL;
  domain->mesh_quantities = NULL;

  domain->cdo_quantities = cs_cdo_quantities_free(domain->cdo_quantities);
  domain->connect = cs_cdo_connect_free(domain->connect);

  domain->boundaries = _free_domain_boundaries(domain->boundaries);

  BFT_FREE(domain->time_step);

  /* Free memory related to equations */
  for (int i = 0; i < domain->n_equations; i++)
    domain->equations[i] = cs_equation_free(domain->equations[i]);
  BFT_FREE(domain->equations);

  BFT_FREE(domain);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_domain_t structure
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_summary(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return;

  /* Output information */
  bft_printf("\n%s", lsepline);
  bft_printf("\tSummary of domain settings\n");
  bft_printf("%s", lsepline);

  bft_printf(" -msg- n_cdo_equations          %d\n", domain->n_equations);
  bft_printf(" -msg- n_predefined_equations   %d\n",
             domain->n_predef_equations);
  bft_printf(" -msg- n_user_equations         %d\n", domain->n_user_equations);

  cs_domain_boundary_t  *bdy = domain->boundaries;

  /* Default boundary */
  bft_printf("\n  Domain boundary by default: ");
  switch (bdy->default_boundary) {
  case CS_PARAM_BOUNDARY_WALL:
    bft_printf(" wall\n");
    break;
  case CS_PARAM_BOUNDARY_SYMMETRY:
    bft_printf(" symmetry\n");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid boundary by default.\n"
                " Please modify your settings."));
  }

  /* Number of border faces for each type of boundary */
  bft_printf("  >> Number of faces with a wall boundary:      %d\n",
             bdy->n_elts[CS_PARAM_BOUNDARY_WALL]);
  bft_printf("  >> Number of faces with a inlet boundary:     %d\n",
             bdy->n_elts[CS_PARAM_BOUNDARY_INLET]);
  bft_printf("  >> Number of faces with a outlet boundary:    %d\n",
             bdy->n_elts[CS_PARAM_BOUNDARY_OUTLET]);
  bft_printf("  >> Number of faces with a symmetry boundary:  %d\n",
             bdy->n_elts[CS_PARAM_BOUNDARY_SYMMETRY]);

  /* Time step summary */
  bft_printf("\n  Time step information\n");
  if (domain->only_steady)
    bft_printf("  >> Steady-state computation");

  else { /* Time information */

    bft_printf("  >> Time step status:");
    if (domain->time_options.idtvar == 0)
      bft_printf("  constant\n");
    else if (domain->time_options.idtvar == 1)
      bft_printf("  variable in time\n");
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid idtvar value for the CDO module.\n"));
    bft_printf("  >> Type of definition: %s",
               cs_param_get_def_type_name(domain->time_step_def_type));
    if (domain->time_step_def_type == CS_PARAM_DEF_BY_VALUE)
      bft_printf(" => %5.3e\n", domain->dt_cur);
    else
      bft_printf("\n");
  }
  bft_printf("\n");

  /* Summary for each equation */
  for (int  eq_id = 0; eq_id < domain->n_equations; eq_id++)
    cs_equation_summary(domain->equations[eq_id]);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary type by default
 *
 * \param[in, out]   domain        pointer to a cs_domain_t structure
 * \param[in]        bdy_name      key name of the default boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_default_boundary(cs_domain_t     *domain,
                               const char      *bdy_name)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_domain_t structure is not allocated."));

  cs_domain_boundary_t  *bcs = domain->boundaries;
  cs_param_boundary_type_t  type = CS_PARAM_N_BOUNDARY_TYPES;

  /* Sanity check */
  assert(bcs != NULL);

  /* Find boundary type associated to bdy_name */
  if (strcmp(bdy_name, "wall") == 0)
    type = CS_PARAM_BOUNDARY_WALL;
  else if (strcmp(bdy_name, "symmetry") == 0)
    type = CS_PARAM_BOUNDARY_SYMMETRY;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key name %s for setting a boundary by default.\n"
                " Available choices are: wall or symmetry."),
              bdy_name);

  bcs->default_boundary = type;
  for (int i = 0; i < bcs->n_b_faces; i++)
    bcs->types[i] = type;
}

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
                       const char                *bdy_name)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_domain_t structure is not allocated."));

  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  /* Sanity check */
  assert(cs_mesh_location_get_type(ml_id) == CS_MESH_LOCATION_BOUNDARY_FACES);

  /* Add this mesh location id to a list of mesh location ids of
     the same type */
  cs_param_boundary_type_t  type = CS_PARAM_N_BOUNDARY_TYPES;
  cs_domain_boundary_t  *bcs = domain->boundaries;

  /* Find boundary type associated to bdy_name */
  if (strcmp(bdy_name, "wall") == 0)
    type = CS_PARAM_BOUNDARY_WALL;
  else if (strcmp(bdy_name, "inlet") == 0)
    type = CS_PARAM_BOUNDARY_INLET;
  else if (strcmp(bdy_name, "outlet") == 0)
    type = CS_PARAM_BOUNDARY_OUTLET;
  else if (strcmp(bdy_name, "symmetry") == 0)
    type = CS_PARAM_BOUNDARY_SYMMETRY;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key name %s for setting a boundary.\n"
                " Available choices are: wall, inlet, outlet or symmetry."),
              bdy_name);

  /* Number of mesh locations already defined for this type of boundary */
  int  n_ml_ids = bcs->bdy_ml[type].n_sub_ids;

  /* Add a new mesh location for this type of boundary */
  BFT_REALLOC(bcs->bdy_ml[type].sub_ids, n_ml_ids + 1, int);
  bcs->bdy_ml[type].sub_ids[n_ml_ids] = ml_id;
  bcs->bdy_ml[type].n_sub_ids += 1;

  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);

  if (elt_ids == NULL)
    for (int i = 0; i < n_elts[0]; i++)
      bcs->types[i] = type;
  else
    for (int i = 0; i < n_elts[0]; i++)
      bcs->types[elt_ids[i]] = type;
}

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
                          int            freq)
{
  if (domain == NULL)
    return;

  if (freq == 0)
    freq = -1;
  domain->output_freq = freq;
}

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
                        void          *defval)
{
  assert(domain != NULL);

  domain->time_step->nt_max = nt_max;
  domain->time_step->t_max = t_end;

  /* Get the type of definition */
  if (strcmp(defkey, "value") == 0) {
    domain->time_step_def_type = CS_PARAM_DEF_BY_VALUE;
    domain->time_step->is_variable = 0; // constant time step
    domain->time_options.idtvar = 0;    // constant time step by default
  }
  else {

    domain->time_step->is_variable = 1; // not constant time step
    domain->time_options.idtvar = 1;    /* uniform in space but can change
                                          from one time step to the other */

    if (strcmp(defkey, "time_func") == 0)
      domain->time_step_def_type = CS_PARAM_DEF_BY_TIME_FUNCTION;
    else if (strcmp(defkey, "analytic") == 0)
      domain->time_step_def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
    else if (strcmp(defkey, "user") == 0)
      domain->time_step_def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key for setting the type of definition.\n"
                  " Given key: %s\n"
                  " Choice among value, field, evaluator, analytic, user, law"
                  " or file\n"
                  " Please modify your settings."), defkey);


  } // type of definition

  cs_param_set_def(domain->time_step_def_type,
                   CS_PARAM_VAR_SCAL,
                   defval,
                   &(domain->time_step_def));

  if (domain->time_step_def_type == CS_PARAM_DEF_BY_VALUE) {

    domain->dt_cur = domain->time_step_def.get.val;
    domain->time_options.dtref = domain->dt_cur;
    domain->time_options.dtmin = domain->dt_cur;
    domain->time_options.dtmax = domain->dt_cur;

  }
  else { /* Default initialization.
            To be changed at first call to cs_domain_time_step_increment() */

    domain->dt_cur = t_end;;
    domain->time_options.dtref = t_end;
    domain->time_options.dtmin = t_end;
    domain->time_options.dtmax = 0.;

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain)
{
  if (domain->only_steady)
    return;
  if (domain->time_step_def_type == CS_PARAM_DEF_BY_VALUE)
    return; // job has already been done during the initial setting

  int  nt_cur = domain->time_step->nt_cur;

  if (nt_cur > 0) {

    cs_param_def_type_t  def_type = domain->time_step_def_type;

    if (def_type == CS_PARAM_DEF_BY_TIME_FUNCTION)
      domain->dt_cur = domain->time_step_def.time_func(nt_cur);

    else if (def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION) {

      cs_real_3_t  xyz = {0, 0, 0}; // not useful
      cs_get_t  get;

      domain->time_step_def.analytic(domain->time_step->t_cur, xyz, &get);
      domain->dt_cur = get.val;

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid way of defining the current time step.\n"
                " Please modify your settings.");

  } /* nt_cur > 0 */

}

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
                       const char         *eqname)
{
  cs_equation_t  *eq = NULL;

  for (int i = 0; i < domain->n_equations; i++) {

    cs_equation_t  *_eq = domain->equations[i];
    if (strcmp(eqname, cs_equation_get_name(_eq)) == 0) {
      eq = _eq;
      break;
    }

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the computation of the wall distance
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_activate_wall_distance(cs_domain_t   *domain)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_domain_t structure is not allocated."));

  domain->wall_distance_eq_id = domain->n_equations;

  domain->n_predef_equations += 1;
  domain->n_equations += 1;
  BFT_REALLOC(domain->equations, domain->n_equations, cs_equation_t *);

  domain->equations[domain->wall_distance_eq_id] =
    cs_equation_create("WallDistance",           // equation name
                       "WallDistance",           // variable name
                       CS_EQUATION_PREDEFINED,   // status
                       CS_EQUATION_TYPE_SCAL,    // type of equation
                       true,                     // steady ?
                       false,                    // convection term ?
                       true,                     // diffusion term ?
                       CS_PARAM_BC_HMG_NEUMANN); // default BC

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup predefined equations which are activated
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_predefined_equations(cs_domain_t   *domain)
{
  assert(domain != NULL);

  cs_domain_boundary_t  *bdy = domain->boundaries;
  cs_equation_t  *eq = NULL;

  if (domain->wall_distance_eq_id > -1) {
    eq = domain->equations[domain->wall_distance_eq_id];

    cs_domain_bdy_ml_t  wall_bdy_ml = bdy->bdy_ml[CS_PARAM_BOUNDARY_WALL];

    cs_walldistance_setup(eq, wall_bdy_ml.id);
  }

}

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
                            const char          *key_bc)
{
  cs_equation_type_t  type = CS_EQUATION_N_TYPES;
  cs_param_bc_type_t  default_bc = CS_PARAM_N_BC_TYPES;

  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_domain_t structure is not allocated."));

  BFT_REALLOC(domain->equations, domain->n_equations + 1, cs_equation_t *);

  /* Define the type of equation */
  if (strcmp(key_type, "scalar") == 0)
    type = CS_EQUATION_TYPE_SCAL;
  else if (strcmp(key_type, "vector") == 0)
    type = CS_EQUATION_TYPE_VECT;
  else if (strcmp(key_type, "tensor") == 0)
    type = CS_EQUATION_TYPE_TENS;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation: %s\n"
                " Choices are scalar, vector or tensor."), key_type);

  /* Define a boundary condition by default */
  if (strcmp(key_bc, "zero_value") == 0)
    default_bc = CS_PARAM_BC_HMG_DIRICHLET;
  else if (strcmp(key_bc, "zero_flux") == 0)
    default_bc = CS_PARAM_BC_HMG_NEUMANN;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of boundary condition by default: %s\n"
                " Choices are zero_value or zero_flux."), key_bc);

  domain->equations[domain->n_equations] =
    cs_equation_create(eqname,           // equation name
                       varname,          // variable name
                       CS_EQUATION_USER, // status
                       type,             // type of equation
                       is_steady,        // steady ?
                       do_convection,    // convection term ?
                       do_diffusion,     // diffusion term ?
                       default_bc);      // default BC

  domain->n_user_equations += 1;
  domain->n_equations += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure whith name eqname
 *         Return NULL if not find
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_create_fields(cs_domain_t  *domain)
{
  /* Loop on user-defined equations */
  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++)
    cs_equation_create_field(domain->equations[eq_id]);
}

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
cs_domain_needs_iterate(cs_domain_t  *domain)
{
  bool  one_more_iter = true;

  cs_time_step_t  *ts = domain->time_step;

  if (ts->nt_cur > ts->nt_max)
    one_more_iter = false;

  if (ts->nt_cur > 0) {

    if (domain->only_steady)
      one_more_iter = false;
    if (ts->t_cur > ts->t_max) /* Be careful to avoid wrong stop due to
                                  an accumulation of small truncation error */
      if (fabs(ts->t_cur - ts->t_max) > 0.1*domain->dt_cur)
        one_more_iter = false;

  } /* nt_cur > 0 */

  return one_more_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain)
{
  cs_time_step_t  *ts = domain->time_step;

  /* Increment time iteration */
  ts->nt_cur++;
  ts->t_prev = ts->t_cur;
  ts->t_cur += domain->dt_cur;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_solve(cs_domain_t  *domain)
{
  int  eq_id;

  int  nt_cur = domain->time_step->nt_cur;

  if (nt_cur == 0) { // Time loop has not yet begun

    /* Output information */
    bft_printf("\n%s", lsepline);
    bft_printf("    Initialize system and/or solve steady-state problem\n");
    bft_printf("%s", lsepline);

    /* Initialize system before resolution for all equations
       - create system builder
       - initialize field according to initial conditions
       - initialize source term
     */
    for (eq_id = 0; eq_id < domain->n_equations; eq_id++) {

      cs_equation_t  *eq = domain->equations[eq_id];

      cs_equation_init_system(domain->mesh,
                              domain->connect,
                              domain->cdo_quantities,
                              domain->time_step,
                              eq);

      if (cs_equation_is_steady(eq)) {

        /* Define the algebraic system */
        cs_equation_build_system(domain->mesh,
                                 domain->connect,
                                 domain->cdo_quantities,
                                 domain->time_step,
                                 domain->dt_cur,
                                 eq);

        /* Solve the algebraic system */
        cs_equation_solve(domain->connect,
                          domain->cdo_quantities,
                          domain->time_step,
                          eq);

      } /* for steady-state equations */

    } // Loop on equations

    /* Specific treatment for wall distance */
    if (domain->wall_distance_eq_id > -1) {
      eq_id = domain->wall_distance_eq_id;
      cs_walldistance_compute(domain->connect,
                              domain->cdo_quantities,
                              domain->equations[eq_id]);
    }

  }
  else { /* nt_cur > 0 */

    /* Output information */
    if (nt_cur % domain->output_freq == 0) {
      bft_printf("\n%s", lsepline);
      bft_printf("    Solve domain for iteration %5d (time = %5.3e s)\n",
                 nt_cur, domain->time_step->t_cur);
      bft_printf("%s", lsepline);
    }

    for (eq_id = 0; eq_id < domain->n_equations; eq_id++) {

      cs_equation_t  *eq = domain->equations[eq_id];

      if (!cs_equation_is_steady(eq)) { /* Unsteady equation */

        /* Define the algebraic system */
        if (cs_equation_needs_build(eq))
          cs_equation_build_system(domain->mesh,
                                   domain->connect,
                                   domain->cdo_quantities,
                                   domain->time_step,
                                   domain->dt_cur,
                                   eq);

        /* Solve domain */
        cs_equation_solve(domain->connect,
                          domain->cdo_quantities,
                          domain->time_step,
                          eq);

      } // unsteady eq.

    } // Loop on equations

  } /* test on nt_cur */

  /* Generic post-processing */
  _domain_post(domain->mesh, domain->time_step);

  /* Generic post-processing related to the equation structure */
  for (eq_id = 0; eq_id < domain->n_equations; eq_id++)
    cs_equation_post(domain->mesh,
                     domain->cdo_quantities,
                     domain->time_step,
                     domain->equations[eq_id]);

  /* User-defined extra operations */
  cs_user_cdo_extra_op(domain);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
