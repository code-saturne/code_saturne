#ifndef __CS_PROBE_H__
#define __CS_PROBE_H__

/*============================================================================
 * Set of structures and functions to handle probes and profiles
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_probe_set_t  cs_probe_set_t;

typedef enum {

  CS_PROBE_MODE_EXACT,               /* Get the value exactly where it is
                                        requested. This mode may require an
                                        interpolation step. */
  CS_PROBE_MODE_NEAREST_CELL_CENTER, /* Get the value at the nearest cell
                                        center */
  CS_PROBE_MODE_NEAREST_VERTEX,      /* Get the value at the nearest vertex */
  CS_PROBE_N_MODES

} cs_probe_mode_t;

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of probe sets defined
 *
 * \return the number of probe sets defined
 */
/*----------------------------------------------------------------------------*/

int
cs_probe_get_n_sets(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_probe_set_t structure
 *
 * \param[in]   name        name of the set of probes to find
 *
 * \return a pointer to a cs_probes_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_get(const char    *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_probe_set_t structure from its id
 *
 * \param[in]   pset_id       id related to the set of probes to find
 *
 * \return a pointer to a cs_probes_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_get_by_id(int   pset_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name related to a cs_probe_set_t structure
 *
 * \param[in]   pset       pointer to a cs_probe_set_t structure
 *
 * \return the name of the cs_probes_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const char *
cs_probe_set_get_name(cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve information useful for the postprocessing step
 *
 * \param[in]  pset          pointer to a cs_probe_set_t structure
 * \param[out] time_varying  true if probe coords may change during computation
 * \param[out] is_profile    true if probe set is related to a profile
 * \param[out] on_boundary   true if probes are located on boundary
 * \param[out] is_automatic  true if the set of variables to post is predefined
 * \param[out] n_writers     number of associated  user-defined writers
 * \param[out] writer_ids    list of writer ids
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_get_post_info(const cs_probe_set_t   *pset,
                           bool                   *time_varying,
                           bool                   *is_profile,
                           bool                   *on_boundary,
                           bool                   *is_automatic,
                           int                    *n_writers,
                           int                    *writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a set of monitoring probes has been defined among all the
 *         probe sets
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_probe_set_have_monitoring(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new set of probes
 *
 * \param[in]   name        name of the set of probes
 *
 * \return a pointer to a new allocated cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_create(const char    *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new probe to an existing set of probes
 *
 * \param[in, out]  pset    set of probes
 * \param[in]       xyz     coordinates of the point to add
 * \param[in]       label   NULL or the name of the point (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_add_probe(cs_probe_set_t     *pset,
                       const cs_real_t    *xyz,
                       const char         *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new set of probes from an array of coordinates
 *
 * \param[in]   name      name of the set of probes
 * \param[in]   n_probes  number of probes in coords and labels
 * \param[in]   coords    list of coordinates related to each probe
 * \param[in]   labels    list of label related to each probe (optional)
 *
 * \return a pointer to a new allocated cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_create_from_array(const char         *name,
                               int                 n_probes,
                               const cs_real_t    *coords,
                               const char        **labels);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new set of probes from the segment spanned by two points
 *
 * \param[in]  name          name of the set of probes
 * \param[in]  n_probes      number of probes
 * \param[in]  start_coords  coordinates of the starting point
 * \param[in]  end_coords    coordinates of the ending point
 *
 * \return a pointer to a new allocated cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_create_from_segment(const char        *name,
                                 int                n_probes,
                                 const cs_real_t    start_coords[3],
                                 const cs_real_t    end_coords[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a list of writers to a probe set
 *
 * \param[in, out] pset        pointer to a cs_probe_set_t structure to set
 * \param[in]      n_writers   number of writers assocuated to this probe set
 * \param[in]      writer_ids  list of writer ids
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_associate_writers(cs_probe_set_t   *pset,
                               int               n_writers,
                               const int        *writer_ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to the management of a set of probes
 *
 * \param[in, out] pset     pointer to a cs_probe_set_t structure to set
 * \param[in]      keyname  name of the keyword related to the parameter to set
 * \param[in]      keyval   value of the keyword to set
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_option(cs_probe_set_t   *pset,
                    const char       *keyname,
                    const char       *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Try to locate each probe and define the coordinate really use for
 *         the postprocessing step
 *
 * \param[in, out]  pset    pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_locate(cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a fvm_nodal_t structure from the set of probes
 *
 * \param[in, out]  pset        pointer to a cs_probe_set_t structure
 * \param[in]       mesh_name   name of the mesh to export
 *
 * \return a pointer to a fvm_nodal_t structure
 */
/*----------------------------------------------------------------------------*/

fvm_nodal_t *
cs_probe_set_export_mesh(cs_probe_set_t   *pset,
                         const char       *mesh_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to a set of probes
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_probe_set_t structure
 *
 * \param[in]  pset    pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_dump(const cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main members of a cs_probe_set_t structure
 *
 * \param[in]      pset      pointer to a cs_probe_set_t structure
 * \param[in, out] mode      mode of location
 * \param[in, out] n_probes  number of probes
 * \param[in, out] coords    probe coordinates
 * \param[in, out] ent_num   entity numbers (-1 if not located on this rank)
 * \param[in, out] distances distance of each probe from its related cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_get_members(const cs_probe_set_t   *pset,
                         int                    *mode,
                         int                    *n_probes,
                         cs_real_t              *coords[],
                         cs_lnum_t              *ent_num[],
                         float                  *distances[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROBE_H__ */
