#ifndef __CS_PROBE_H__
#define __CS_PROBE_H__

/*============================================================================
 * Set of structures and functions to handle probes and profiles
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "fvm_nodal.h"

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

  CS_PROBE_SNAP_NONE,            /*!< No position change */
  CS_PROBE_SNAP_ELT_CENTER,      /*!< snap to nearest cell or face center */
  CS_PROBE_SNAP_VERTEX           /*!> snap to nearest vertex */

} cs_probe_snap_t;

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer to definition of probes based on rank-local points.
 *
 * If non-empty and not containing all elements, a list of coordinates
 * as well as a list of curvilinear coordinates should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * Those list's lifecycle is then managed automatically by the
 * probe set object.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * \param[in, out]  input   pointer to optional (untyped) value or structure
 * \param[out]      n_elts  number of selected coordinates
 * \param[out]      coords  coordinates of selected elements.
 * \param[out]      s       curvilinear coordinates of selected elements
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_probe_set_define_local_t) (void          *input,
                               cs_lnum_t     *n_elts,
                               cs_real_3_t  **coords,
                               cs_real_t    **s);

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to a set of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of probe sets defined.
 *
 * \return the number of probe sets defined
 */
/*----------------------------------------------------------------------------*/

int
cs_probe_get_n_sets(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a cs_probe_set_t structure.
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
 * \brief  Retrieve a cs_probe_set_t structure from its id.
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
 * \brief  Retrieve the name related to a cs_probe_set_t structure.
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
 * \brief  Retrieve information useful for the postprocessing step.
 *
 * Output arguments may be set to NULL if we do not need to query them.
 *
 * \param[in]  pset            pointer to a cs_probe_set_t structure
 * \param[out] time_varying    true if probe locations may change with time
 * \param[out] on_boundary     true if probes are located on boundary
 * \param[out] on_curve        true if the probe set has cuvilinear coordinates
 * \param[out] auto_variables  true if set of variables to output is predefined
 * \param[out] n_writers       number of associated  user-defined writers,
 *                             or -1 if default unchanged
 * \param[out] writer_ids      pointer to a list of writer ids
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_get_post_info(const cs_probe_set_t   *pset,
                           bool                   *time_varying,
                           bool                   *on_boundary,
                           bool                   *on_curve,
                           bool                   *auto_variables,
                           int                    *n_writers,
                           int                    *writer_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the location filter selection criteria string for a
 *         given probe set
 *
 * \param[in]   pset       pointer to a cs_probe_set_t structure
 *
 * \return selection criteria string, or NULL if no filter defined
 */
/*----------------------------------------------------------------------------*/

const char *
cs_probe_set_get_location_criteria(cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new set of probes.
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
 * \brief  Add a new probe to an existing set of probes.
 *
 * \param[in, out]  pset    set of probes
 * \param[in]       x       x coordinate  of the point to add
 * \param[in]       y       y coordinate  of the point to add
 * \param[in]       z       z coordinate  of the point to add
 * \param[in]       label   NULL or the name of the point (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_add_probe(cs_probe_set_t     *pset,
                       cs_real_t           x,
                       cs_real_t           y,
                       cs_real_t           z,
                       const char         *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new set of probes from an array of coordinates.
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
cs_probe_set_create_from_array(const char          *name,
                               int                  n_probes,
                               const cs_real_3_t   *coords,
                               const char         **labels);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new set of probes from the segment spanned by two points.
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
 * \brief Define a new set of probes from rank-local definition function.
 *
 * The local definition function given by the \ref p_define_func pointer
 * is called just before locating probes on the parent mesh, so this allows
 * building probe sets based on subsets of the computational mesh.
 *
 * Note: if the p_define_input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * \param[in]  name           name of the set of probes
 * \param[in]  p_define_func  function used for local definition
 * \param[in]  p_define_input optional input for local definition function
 *
 * \return a pointer to a new allocated cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_create_from_local(const char                   *name,
                               cs_probe_set_define_local_t  *p_define_func,
                               void                         *p_define_input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  allow overwriting the definition of a given probe set.
 *
 * If no a probe set of the given name exists, the operation is ignored.
 *
 * \param[in]  name  name of the probe set
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_allow_overwrite(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a list of writers to a probe set.
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
 * \brief  Set to true or false the automatic post-processing of variables
 *
 * \param[in, out] pset     pointer to a cs_probe_set_t structure
 * \param[in]      mode     true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_auto_var(cs_probe_set_t   *pset,
                      bool              mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set snap mode related to the management of a set of probes.
 *
 * \param[in, out] pset        pointer to a cs_probe_set_t structure
 * \param[in]      snap_mode   snap mode to set
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_snap_mode(cs_probe_set_t   *pset,
                       cs_probe_snap_t   snap_mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to the management of a set of probes
 *
 * Available option key names accepting \c true or \c false:
 *
 * - \c \b transient_location if \c true, relocate probes relative to
 *         deforming or moving mesh (default: \c false)
 * - \c \b boundary  if \ c true, locate on boundary mesh; if \c false,
 *         locate on volume mesh (default)
 *
 * Other options:
 *
 * - \c \b selection_criteria where keyval is selection criteria string
 * - \c \b tolerance  where keyval is for instance "0.05" (default "0.10")
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
 * \brief  Try to locate each probe and define the coordinate really used for
 *         the postprocessing step.
 *
 * For better performance when using multiple probe sets, a pointer to
 * an existing location mesh may be passed to this function. The caller is
 * responsible for ensuring this mesh matches selection criteria for the
 * probe set.
 *
 * \param[in, out]  pset           pointer to a cs_probe_set_t structure
 * \param[in]       location_mesh  optional pointer to mesh relative to which
 *                                 probe set should be located, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_locate(cs_probe_set_t     *pset,
                    const fvm_nodal_t  *location_mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a fvm_nodal_t structure from the set of probes.
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
 * \brief  Define a fvm_nodal_t structure from the set of unlocated probes.
 *
 * \param[in, out]  pset        pointer to a cs_probe_set_t structure
 * \param[in]       mesh_name   name of the mesh to export
 *
 * \return a pointer to a fvm_nodal_t structure
 */
/*----------------------------------------------------------------------------*/

fvm_nodal_t *
cs_probe_set_unlocated_export_mesh(cs_probe_set_t   *pset,
                                   const char       *mesh_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_probe_set_t structure.
 *
 * \param[in]  pset    pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_dump(const cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the main members of a cs_probe_set_t structure.
 *
 * \param[in]       pset       pointer to a cs_probe_set_t structure
 * \param[in, out]  snap_mode  mode of location
 * \param[in, out]  n_probes   number of probes
 * \param[in, out]  coords     probe coordinates
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_get_members(const cs_probe_set_t   *pset,
                         cs_probe_snap_t        *snap_mode,
                         int                    *n_probes,
                         cs_real_3_t            *coords[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the number probes in the local domain.
 *
 * \param[in]       pset       pointer to a cs_probe_set_t structure
 *
 * \return  number of probes in local domain
 */
/*----------------------------------------------------------------------------*/

int
cs_probe_set_get_n_local(const cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the list of curvilinear abscissa for the given probe set
 *
 * \param[in]  pset              pointer to a cs_probe_set_t structure
 *
 * \return NULL or the pointer to the array of abscissa
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_probe_set_get_curvilinear_abscissa(const cs_probe_set_t   *pset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the ids of a probe set's local matching elements, relative
 *         to a given mesh location.
 *
 * The mesh_location id must match one of \ref CS_MESH_LOCATION_CELLS,
 * \ref CS_MESH_LOCATION_BOUNDARY_FACES, or \ref CS_MESH_LOCATION_VERTICES.
 *
 * \param[in]  pset              pointer to a cs_probe_set_t structure
 * \param[in]  mesh_location_id  id of parent mesh location
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_probe_set_get_elt_ids(const cs_probe_set_t  *pset,
                         int                    mesh_location_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROBE_H__ */
