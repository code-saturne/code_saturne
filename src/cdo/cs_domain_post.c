/*============================================================================
 * Manage specific post-processing related to a computational domain
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_post.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain_post.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

typedef struct { /* Only shared with a cs_domain_t structure */

  double                       dt_cur;

  const cs_cdo_quantities_t   *quant;

  /* Properties attached to the computational domain */
  int                          n_properties;
  cs_property_t              **properties;

  /* Advection fields attached to the computational domain */
  int                          n_adv_fields;
  cs_adv_field_t             **adv_fields;

  /* Number of equations defined on this domain */
  int                          n_equations;
  cs_equation_t              **equations;

} cs_domain_post_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_domain_post_t  *domain_post = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for advection fields.
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 * \param[in]  dt_cur      value of the current time step
 */
/*----------------------------------------------------------------------------*/

static void
_post_advection_field(const cs_adv_field_t       *adv,
                      const cs_cdo_quantities_t  *quant,
                      const cs_time_step_t       *time_step,
                      double                      dt_cur)
{
  if (adv == NULL)
    return;
  if (adv->post_flag == 0)
    return;

  const bool post_courant =
    (adv->post_flag & CS_ADVECTION_FIELD_POST_COURANT) ? true : false;

  if (post_courant) { /* Compute and postprocess the Courant number */

    double  hc;
    cs_nvec3_t  adv_c;

    char  *label = NULL;
    double  *courant = NULL;
    int  len = strlen(adv->name) + 8 + 1;

    BFT_MALLOC(courant, quant->n_cells, double);
    BFT_MALLOC(label, len, char);
    sprintf(label, "%s.Courant", adv->name);

    if (adv->cell_field_id > -1) { /* field is defined at cell centers */

      cs_field_t  *fld = cs_field_by_id(adv->cell_field_id);

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_nvec3(fld->val + 3*c_id, &adv_c);
        hc = pow(quant->cell_vol[c_id], cs_math_onethird);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }

    }
    else {

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_advection_field_get_cell_vector(c_id, adv, &adv_c);
        hc = pow(quant->cell_vol[c_id], cs_math_onethird);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }

    }

    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      label,
                      1,
                      true,           // interlace
                      true,           // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      courant,        // values on cells
                      NULL,           // values at internal faces
                      NULL,           // values at border faces
                      time_step);     // time step management struct.

    BFT_FREE(label);
    BFT_FREE(courant);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for property structures
 *
 * \param[in]  pty         pointer to the property structure
 * \param[in]  quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t structure
 * \param[in]  dt_cur      value of the current time step
 */
/*----------------------------------------------------------------------------*/

static void
_post_property(const cs_property_t        *pty,
               const cs_cdo_quantities_t  *quant,
               const cs_time_step_t       *time_step,
               double                      dt_cur)
{
  if (pty == NULL)
    return;
  if (pty->post_flag == 0)
    return;

  if (pty->post_flag & CS_PROPERTY_POST_FOURIER) {

    char  *label = NULL;
    double  *fourier = NULL;
    int len = strlen(pty->name) + 8 + 1;

    BFT_MALLOC(fourier, quant->n_cells, double);
    BFT_MALLOC(label, len, char);
    sprintf(label, "%s.Fourier", pty->name);

    /* Compute Fourier number */
    cs_property_get_fourier(pty, dt_cur, fourier);

    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      label,
                      1,
                      true,           // interlace
                      true,           // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      fourier,        // values on cells
                      NULL,           // values at internal faces
                      NULL,           // values at border faces
                      time_step);     // time step management struct.

    BFT_FREE(label);
    BFT_FREE(fourier);

  } // Post Fourier number

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the computational domain
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_domain_post(void                      *input,
             int                        mesh_id,
             int                        cat_id,
             int                        ent_flag[5],
             cs_lnum_t                  n_cells,
             cs_lnum_t                  n_i_faces,
             cs_lnum_t                  n_b_faces,
             const cs_lnum_t            cell_ids[],
             const cs_lnum_t            i_face_ids[],
             const cs_lnum_t            b_face_ids[],
             const cs_time_step_t      *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == NULL)
    return;

  if (mesh_id != -1) /* Post-processing only on the generic volume mesh */
    return;

  cs_domain_post_t  *dp = (cs_domain_post_t *)input;

  /* Post-processing related to advection fields */
  for (int adv_id = 0; adv_id < dp->n_adv_fields; adv_id++)
    _post_advection_field(dp->adv_fields[adv_id],
                          dp->quant, time_step, dp->dt_cur);

  /* Post-processing related to properties */
  for (int pty_id = 0; pty_id < dp->n_properties; pty_id++)
    _post_property(dp->properties[pty_id],
                   dp->quant, time_step, dp->dt_cur);

  /* Post-processing related to equations */
  for (int eq_id = 0; eq_id < dp->n_equations; eq_id++)
    cs_equation_extra_post(dp->equations[eq_id], time_step, dp->dt_cur);

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the generic post-processing related to a domain
 *
 * \param[in]  dt             reference time step value
 * \param[in]  n_adv_fields   number of advection fields
 * \param[in]  adv_fields     pointer on advection field pointers
 * \param[in]  n_properties   number of properties
 * \param(in]  properties     pointer on properties pointers
 * \param[in]  n_equations    number of equations
 * \param[in]  equations      pointer on equations pointers
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_init(double                dt,
                    cs_cdo_quantities_t  *quant,
                    int                   n_adv_fields,
                    cs_adv_field_t      **adv_fields,
                    int                   n_properties,
                    cs_property_t       **properties,
                    int                   n_equations,
                    cs_equation_t       **equations)
{
  BFT_MALLOC(domain_post, 1, cs_domain_post_t);

  domain_post->dt_cur = dt;

  /* Shared */
  domain_post->quant = quant;
  domain_post->n_adv_fields = n_adv_fields;
  domain_post->adv_fields = adv_fields;
  domain_post->n_properties = n_properties;
  domain_post->properties = properties;
  domain_post->n_equations = n_equations;
  domain_post->equations = equations;

  /* Set pointers of function if additional postprocessing is requested */
  cs_post_add_time_mesh_dep_output(_domain_post, domain_post);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the hidden view of the domain dedicated for post-processing
 *
 * \param[in]    dt      current value of the time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_update(double    dt)
{
  assert(domain_post != NULL);
  domain_post->dt_cur = dt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate writers and output meshes if needed
 *
 * \param[in]  time_step    pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_activate(cs_time_step_t    *time_step)
{
  cs_post_time_step_begin(time_step);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the hidden view of the domain dedicated for post-processing
 *
 * \param[in]  time_step    pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post(cs_time_step_t    *time_step)
{
  assert(domain_post != NULL);

  /* Predefined extra-operations related to
      - the domain (advection fields and properties),
      - equations
      - groundwater flows
     are also handled during the call of this function thanks to
     cs_post_add_time_mesh_dep_output() function pointer
  */
  cs_post_time_step_output(time_step);

  cs_post_time_step_end();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize post-processing related to the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_finalize(void)
{
  BFT_FREE(domain_post);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
