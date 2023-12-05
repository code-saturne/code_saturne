/*============================================================================
 * Volume zones handling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_flag_check.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_volume_zone.c
        Volume zone handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Zone descriptor allocation block size */

#define _CS_ZONE_S_ALLOC_SIZE       16

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*! Volume zone ids */

typedef struct {
  cs_lnum_t         n_elts;   /*!< local number of associated elements */
  const cs_lnum_t  *max_id;   /*!< local cell ids, or NULL if trivial */
} cs_volume_zone_id_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Zone definitions */

static int  _n_zones = 0;
static int  _n_zones_max = 0;
static cs_zone_t     **_zones = NULL;
static cs_map_name_to_id_t   *_zone_map = NULL;

/* Volume zone id associated with cells */

static int *_zone_id = NULL;

/* Names for logging */

static const int _n_type_flags = 5;

static const int _type_flag_mask[] = {CS_VOLUME_ZONE_INITIALIZATION,
                                      CS_VOLUME_ZONE_POROSITY,
                                      CS_VOLUME_ZONE_HEAD_LOSS,
                                      CS_VOLUME_ZONE_SOURCE_TERM,
                                      CS_VOLUME_ZONE_MASS_SOURCE_TERM,
                                      CS_VOLUME_ZONE_GWF_SOIL,
                                      CS_VOLUME_ZONE_SOLID};

static const char *_type_flag_name[] = {N_("initialization"),
                                        N_("porosity"),
                                        N_("head loss"),
                                        N_("source term"),
                                        N_("mass source term"),
                                        N_("groundwater soil"),
                                        N_("solid")};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its name if present.
 *
 * If no volume zone of the given name is defined, NULL is returned.
 *
 * \param[in]  name  volume zone name
 *
 * \return  pointer to (read only) zone structure, or NULL
 */
/*----------------------------------------------------------------------------*/

static cs_zone_t  *
_zone_by_name_try(const char  *name)
{
  cs_zone_t *z = NULL;
  int id = cs_map_name_to_id_try(_zone_map, name);

  if (id > -1)
    z = _zones[id];

  return z;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a volume zone.
 *
 * \param[in]  name  zone name
 *
 * \return   pointer to new zone.
 */
/*----------------------------------------------------------------------------*/

static cs_zone_t *
_zone_define(const char  *name)
{
  int zone_id = -1;

  cs_zone_t *z = _zone_by_name_try(name);

  /* Check this name was not already used */

  if (z != NULL)
    return z;

  /* Initialize if necessary */

  if (_zone_map == NULL)
    _zone_map = cs_map_name_to_id_create();

  size_t l = 0;
  if (name != NULL)
    l = strlen(name);
  if (l == 0)
    bft_error(__FILE__, __LINE__, 0, _("Defining a zone requires a name."));

  /* Insert entry in map */

  zone_id = cs_map_name_to_id(_zone_map, name);

  if (zone_id == _n_zones)
    _n_zones = zone_id + 1;

  /* Reallocate zones pointer if necessary */

  if (_n_zones > _n_zones_max) {
    if (_n_zones_max == 0)
      _n_zones_max = 8;
    else
      _n_zones_max *= 2;
    BFT_REALLOC(_zones, _n_zones_max, cs_zone_t *);
  }

  /* Allocate zones descriptor block if necessary
     (to reduce fragmentation and improve locality of zone
     descriptors, they are allocated in blocks) */

  int shift_in_alloc_block = zone_id % _CS_ZONE_S_ALLOC_SIZE;
  if (shift_in_alloc_block == 0)
    BFT_MALLOC(_zones[zone_id], _CS_ZONE_S_ALLOC_SIZE, cs_zone_t);
  else
    _zones[zone_id] =   _zones[zone_id - shift_in_alloc_block]
                      + shift_in_alloc_block;

  /* Assign zone */

  z = _zones[zone_id];

  z->name = cs_map_name_to_id_reverse(_zone_map, zone_id);

  z->id = zone_id;
  z->type = 0;

  z->location_id = 0;

  z->n_elts = 0;
  z->elt_ids = NULL;

  z->time_varying = false;
  z->allow_overlay = true;

  z->n_g_elts = 0;

  z->measure = -1.;
  z->boundary_measure = -1.;

  for (int idim = 0; idim < 3; idim++)
    z->cog[idim] = 0.;

  return z;
}

/*----------------------------------------------------------------------------
 * Add type flag info to the current position in the setup log.
 *
 * parameters:
 *   type <-- type flag
 *----------------------------------------------------------------------------*/

static inline void
_log_type(int type)
{
  if (type == 0)
    return;

  int n_loc_flags = 0;

  cs_log_printf(CS_LOG_SETUP,
                _("    type:                       %d"), type);

  for (int i = 0; i < _n_type_flags; i++) {
    if (type & _type_flag_mask[i]) {
      if (n_loc_flags == 0)
        cs_log_printf(CS_LOG_SETUP, " (%s", _(_type_flag_name[i]));
      else
        cs_log_printf(CS_LOG_SETUP, ", %s", _(_type_flag_name[i]));
      n_loc_flags++;
    }
  }

  if (n_loc_flags > 0)
    cs_log_printf(CS_LOG_SETUP, ")\n");
  else
    cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute geometrical measure of a volume zone (volume and surface)
 *
 * For time-varying zones, the associated mesh location is updated.
 *
 * \param[in]  mesh_modified indicate if mesh has been modified
 * \param[in]  z             zone for which measures need to be computed
 */
/*----------------------------------------------------------------------------*/

static void
_volume_zone_compute_metadata(bool       mesh_modified,
                              cs_zone_t *z)
{
  /* We recompute values only if mesh is modified or zone is time varying.
   * FIXME: For the moment, the boundary measure is not computed, but set to -1.
   * to be improved in the future.
   */
  if (z->time_varying || mesh_modified) {

    cs_real_t *cell_vol   = cs_glob_mesh_quantities->cell_vol;
    cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
    cs_real_3_t *cell_cen = (cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

    z->measure = 0.;
    z->f_measure = 0.;
    z->boundary_measure = -1.;
    z->f_boundary_measure = -1.;
    for (int idim = 0; idim < 3; idim++)
      z->cog[idim] = 0.;

    for (cs_lnum_t e_id = 0; e_id < z->n_elts; e_id++) {
      cs_lnum_t c_id = z->elt_ids[e_id];
      z->measure   += cell_vol[c_id];
      z->f_measure += cell_f_vol[c_id];
      for (int idim = 0; idim < 3; idim++)
        z->cog[idim] += cell_cen[c_id][idim] * cell_vol[c_id];
    }

    cs_real_t measures[7] = {z->measure, z->f_measure,
                             z->boundary_measure, z->f_boundary_measure,
                             z->cog[0], z->cog[1], z->cog[2]};

    cs_parall_sum(7, CS_REAL_TYPE, measures);

    cs_gnum_t n_g_elts = z->n_elts;
    cs_parall_sum(1, CS_GNUM_TYPE, &n_g_elts);

    z->n_g_elts = n_g_elts;

    z->measure = measures[0];
    z->f_measure = measures[1];
    z->boundary_measure = measures[2];
    z->f_boundary_measure = measures[3];

    /* Avoid a SIGFPE error */
    if (fabs(measures[0]) > DBL_MIN)
      for (int idim = 0; idim < 3; idim++)
        z->cog[idim] = measures[4+idim] / measures[0];

  } /* Need to compute metadata */

}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize volume zone structures.
 *
 * This defines a default volume zone. This is the first function of
 * the volume zone handling functions which should be called, and it should
 * only be called after \ref cs_mesh_location_initialize.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_initialize(void)
{
  assert(_n_zones == 0);

  cs_mesh_location_set_explicit_ids(CS_MESH_LOCATION_CELLS, true);

  const char *name = cs_mesh_location_get_name(CS_MESH_LOCATION_CELLS);

  cs_zone_t *z = _zone_define(name);

  z->location_id = CS_MESH_LOCATION_CELLS;

  z->type = 0;

  z->allow_overlay = true;

  assert(z->id == 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all volume zone structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_finalize(void)
{
  BFT_FREE(_zone_id);

  for (int i = 0; i < _n_zones; i++) {
    if (i % _CS_ZONE_S_ALLOC_SIZE == 0)
      BFT_FREE(_zones[i]);
  }

  BFT_FREE(_zones);

  cs_map_name_to_id_destroy(&_zone_map);

  _n_zones = 0;
  _n_zones_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones defined.
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_n_zones(void)
{
  return _n_zones;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones which may vary in time.
 *
 * \return  number of zones which may vary in time
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_n_zones_time_varying(void)
{
  int count = 0;

  for (int i = 0; i < _n_zones; i++) {
    if (_zones[i]->time_varying)
      count += 1;
  }

  return count;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update association of volume zones with a mesh.
 *
 * For time-varying zones, the associated mesh location is updated.
 *
 * \param[in]  mesh_modified  indicate if mesh has been modified
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_build_all(bool  mesh_modified)
{
  cs_mesh_t  *m = cs_glob_mesh;
  bool has_time_varying = false;

  /* update zone lists */

  for (int i = 0; i < _n_zones; i++) {
    cs_zone_t *z = _zones[i];
    if (z->time_varying) {
      cs_mesh_location_build(m, z->location_id);
      has_time_varying = true;
    }
    z->n_elts = cs_mesh_location_get_n_elts(z->location_id)[0];
    z->elt_ids = cs_mesh_location_get_elt_ids(z->location_id);
  }

  /* Assign maximum zone id and check for overlap errors
     (start with zone 1, as 0 is default) */

  if (mesh_modified)
    BFT_REALLOC(_zone_id, m->n_cells_with_ghosts, int);

  if (mesh_modified || has_time_varying) {

    cs_lnum_t n_cells = m->n_cells_with_ghosts;

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i <n_cells; i++)
      _zone_id[i] = 0;

    int overlap_error[2] = {_n_zones, _n_zones};

    for (int i = 1; i < _n_zones; i++) {
      cs_zone_t *z = _zones[i];
      for (cs_lnum_t j = 0; j < z->n_elts; j++) {
        cs_lnum_t c_id = z->elt_ids[j];
        int z_id_prev = _zone_id[c_id];
        if (z_id_prev == 0)
          _zone_id[c_id] = z->id;
        else if (_zones[z_id_prev]->allow_overlay)
          _zone_id[c_id] = z->id;
        else if (overlap_error[0] == _n_zones) {
          overlap_error[0] = z_id_prev;
          overlap_error[1] = z->id;
          break;
        }
      }
    }

    cs_parall_min(2, CS_INT_TYPE, overlap_error);

    if (overlap_error[0] < _n_zones) {

      for (int i = 1; i < _n_zones; i++) {
        cs_zone_t *z = _zones[i];
        for (cs_lnum_t j = 0; j < z->n_elts; j++) {
          cs_lnum_t c_id = z->elt_ids[j];
          int z_id_prev = CS_ABS(_zone_id[c_id]);
          if (z_id_prev == 0)
            _zone_id[c_id] = z->id;
          else if (   _zones[z_id_prev]->allow_overlay
                   && _zone_id[c_id] > 0)
            _zone_id[c_id] = z->id;
          else
            _zone_id[c_id] = -z->id;
        }
      }

      cs_flag_check_error_info(_("cell with forbidden zone overlap"),
                               _("zone id"),
                               _("zone_id"),
                               _("Cells with zone error"),
                               _("Cells with valid zones"),
                               CS_MESH_LOCATION_CELLS,
                               0, /* min_flag */
                               _zone_id);

      int i0 = overlap_error[0], i1 = overlap_error[1];

      bft_error(__FILE__, __LINE__, 0,
                _("Volume zone %i (\"%s\") contains at least\n"
                  "one cell already marked with zone id %d (\"%s\").\n\n"
                  "Check definitions or allow overlays for this zone."),
                i1, _zones[i1]->name, i0, _zones[i0]->name);

    }

    /* Compute or update zone geometrical measures */
    for (int i = 0; i < _n_zones; i++) {
      cs_zone_t *z = _zones[i];
      _volume_zone_compute_metadata(mesh_modified, z);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new volume zone using a selection criteria string.
 *
 * \param[in]  name       name of location to define
 * \param[in]  criteria   selection criteria for associated elements
 * \param[in]  type_flag  mask of zone category values
 *
 * \return  id of newly defined volume zone
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_define(const char  *name,
                      const char  *criteria,
                      int          type_flag)
{
  if (criteria == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: selection criteria string must be non-null."),
              __func__);

  cs_zone_t *z = _zone_define(name);

  if (strcmp(criteria, "all[]"))
    z->location_id = cs_mesh_location_add(name,
                                          CS_MESH_LOCATION_CELLS,
                                          criteria);
  else
    z->location_id = CS_MESH_LOCATION_CELLS;

  z->type = type_flag;

  return z->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new mesh location with an associated selection function.
 *
 * So as to define a subset of mesh entities of a given type, a pointer
 * to a selection function may be given.
 *
 * This requires more programming but allows finer control than selection
 * criteria, as the function has access to the complete mesh structure.
 *
 * \param[in]  name        name of location to define
 * \param[in]  func        pointer to selection function for associated elements
 * \param[in, out]  input  pointer to optional (untyped) value
 *                         or structure.
 * \param[in]  type_flag   mask of zone category values
 *
 * \return  id of newly defined created mesh location
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_define_by_func(const char                 *name,
                              cs_mesh_location_select_t  *func,
                              void                       *input,
                              int                         type_flag)
{
  if (func == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: selection function pointer must be non-null."),
              __func__);

  cs_zone_t *z = _zone_define(name);

  z->location_id = cs_mesh_location_add_by_func(name,
                                                CS_MESH_LOCATION_CELLS,
                                                func,
                                                input);

  z->type = type_flag;

  return z->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its id.
 *
 * This function requires that a volume zone of the given id is defined.
 *
 * \param[in]  id   zone id
 *
 * \return  pointer to the volume zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_zone_t  *
cs_volume_zone_by_id(int  id)
{
  if (id > -1 && id < _n_zones)
    return _zones[id];
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Volume zone with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its name if present.
 *
 * This function requires that a volume zone of the given name is defined.
 *
 * \param[in]  name  volume zone name
 *
 * \return  pointer to (read-only) zone structure
 */
/*----------------------------------------------------------------------------*/

const cs_zone_t  *
cs_volume_zone_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_zone_map, name);

  if (id > -1)
    return _zones[id];
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Volume zone \"%s\" is not defined."), name);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a volume zone based on its name if present.
 *
 * If no volume zone of the given name is defined, NULL is returned.
 *
 * \param[in]  name  volume zone name
 *
 * \return  pointer to (read only) zone structure, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_zone_t  *
cs_volume_zone_by_name_try(const char  *name)
{
  const cs_zone_t *z = NULL;
  int id = cs_map_name_to_id_try(_zone_map, name);

  if (id > -1)
    z = _zones[id];

  return z;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the volume zone id from its zone name.
 *         If the zone name is equal to NULL or has an empty length, then
 *         the default zone id (=0) corresponding to all entities is returned
 *
 * \param[in] z_name        name of the zone or NULL or ""
 *
 * \return the id of the volume zone
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_id_by_name(const char   *z_name)
{
  int  id = 0; /* Return the default zone correspoinding to all cells by
                  default */

  if (z_name != NULL) {
    if (strlen(z_name) > 0) {

      id = cs_map_name_to_id_try(_zone_map, z_name);
      if (id < 0)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Volume zone \"%s\" not found.", __func__, z_name);

    }
  }

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set type flag for a given volume zone.
 *
 * \param[in]  id         volume zone id
 * \param[in]  type_flag  volume zone type flag
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_type(int   id,
                        int   type_flag)
{
  const cs_zone_t *z0 = cs_volume_zone_by_id(id);

  _zones[z0->id]->type |= type_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set time varying behavior for a given volume zone.
 *
 * \param[in]  id            volume zone id
 * \param[in]  time_varying  true if the zone's definition varies in time
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_time_varying(int   id,
                                bool  time_varying)
{
  const cs_zone_t *z0 = cs_volume_zone_by_id(id);

  _zones[z0->id]->time_varying = time_varying;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set overlay behavior for a given volume zone.
 *
 * \param[in]  id             volume zone id
 * \param[in]  allow_overlay  true if the zone may be overlayed by another
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_set_overlay(int   id,
                           bool  allow_overlay)
{
  const cs_zone_t *z0 = cs_volume_zone_by_id(id);

  _zones[z0->id]->allow_overlay = allow_overlay;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to zone id associated with each cell.
 *
 * In case of overlayed zones, the highest zone id associated with
 * a given cell is given.
 */
/*----------------------------------------------------------------------------*/

const int *
cs_volume_zone_cell_zone_id(void)
{
  return (const int *)_zone_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given volume zone to log file.
 *
 * \param[in]  z   pointer to volume zone structure
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_log_info(const cs_zone_t  *z)
{
  if (z == NULL)
    return;

  /* Global indicators */
  /*-------------------*/

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "  Zone: \"%s\"\n"
                  "    id:                         %d\n"),
                z->name, z->id);

  _log_type(z->type);

  cs_log_printf(CS_LOG_SETUP,
                _("    location_id:                %d\n"),
                z->location_id);

  if (z->time_varying)
    cs_log_printf(CS_LOG_SETUP, _("    time varying\n"));
  if (z->allow_overlay)
    cs_log_printf(CS_LOG_SETUP, _("    allow overlay\n"));

  cs_mesh_location_def_t ml_def =
    cs_mesh_location_get_definition_method(z->location_id);

  if (ml_def == CS_MESH_LOCATION_DEF_SELECTION_STR) {
    const char *sel_str = cs_mesh_location_get_selection_string(z->location_id);
    cs_log_printf(CS_LOG_SETUP,
                  _("    selection criteria:         \"%s\"\n"),
                  sel_str);
  }
  else if (ml_def == CS_MESH_LOCATION_DEF_SELECTION_FUNC) {
    cs_mesh_location_select_t *sel_fp
      = cs_mesh_location_get_selection_function(z->location_id);

    cs_log_printf(CS_LOG_SETUP,
                  _("    selection function:         %p\n"),
                  (void *)sel_fp);
  }
  else if (ml_def == CS_MESH_LOCATION_DEF_UNION) {
    int n_sub_ids = cs_mesh_location_get_n_sub_ids(z->location_id);
    int *sub_ids  = cs_mesh_location_get_sub_ids(z->location_id);

    bool is_complement = cs_mesh_location_is_complement(z->location_id);
    if (!is_complement)
      cs_log_printf(CS_LOG_SETUP,
                    _("    Union of %d mesh locations:\n"),
                    n_sub_ids);
    else
      cs_log_printf(CS_LOG_SETUP,
                    _("    Complement of %d mesh locations:\n"),
                    n_sub_ids);
    for (int iid = 0; iid < n_sub_ids; iid++) {
      cs_log_printf(CS_LOG_SETUP,
                    _("      sub-location %d/%d\n"),
                    iid+1,
                    n_sub_ids);

      int sub_ml_id = sub_ids[iid];

      cs_log_printf(CS_LOG_SETUP,
                    _("        location_id:            %d\n"),
                    sub_ml_id);

      cs_mesh_location_def_t sub_ml_def =
        cs_mesh_location_get_definition_method(sub_ml_id);

      if (sub_ml_def == CS_MESH_LOCATION_DEF_SELECTION_STR) {
        const char *sub_sel_str =
          cs_mesh_location_get_selection_string(sub_ml_id);

        cs_log_printf(CS_LOG_SETUP,
                      _("        selection criteria:     \"%s\"\n"),
                      sub_sel_str);
      }
      else if (sub_ml_def == CS_MESH_LOCATION_DEF_SELECTION_FUNC) {
        cs_mesh_location_select_t *sub_select_fp
          = cs_mesh_location_get_selection_function(sub_ml_id);

        cs_log_printf(CS_LOG_SETUP,
                      _("        selection function:     %p\n"),
                      (void *)sub_select_fp);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log setup information relative to defined volume zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_log_setup(void)
{
  if (_n_zones == 0)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Volume zones\n"
                  "------------\n"));

  for (int i = 0; i < _n_zones; i++)
    cs_volume_zone_log_info(_zones[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zones associated with a
 *        given zone flag.
 *
 * \param[in]  type_flag  flag to compare to zone type
 *
 * \return  number of zones matching the given type flag
 */
/*----------------------------------------------------------------------------*/

int
cs_volume_zone_n_type_zones(int  type_flag)
{
  int count = 0;

  for (int i = 0; i < _n_zones; i++) {
    if (_zones[i]->type & type_flag)
      count += 1;
  }

  return count;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of volume zone cells associated with a
 *        given zone flag.
 *
 * Note that in the case of overlapping zones, a cell may be accounted
 * for multiple times.
 *
 * \param[in]  type_flag  flag to compare to zone type
 *
 * \return  number of cells in zones matching the given type flag
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_volume_zone_n_type_cells(int  type_flag)
{
  cs_lnum_t count = 0;

  for (int i = 0; i < _n_zones; i++) {
    if (_zones[i]->type & type_flag)
      count += _zones[i]->n_elts;
  }

  return count;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells associated with volume zones of a given type.
 *
 * Note that in the case of overlapping zones, a cell may be accounted
 * for multiple times.
 *
 * \param[in]   type_flag   flag to compare to zone type
 * \param[out]  cell_ids    ids of selected cells (size: given by
 *                          \ref cs_volume_zone_n_type_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_select_type_cells(int        type_flag,
                                 cs_lnum_t  cell_ids[])
{
  cs_lnum_t count = 0;

  for (int i = 0; i < _n_zones; i++) {
    const cs_zone_t *z = _zones[i];
    if (z->type & type_flag) {
      const cs_lnum_t _n_cells = z->n_elts;
      const cs_lnum_t *_cell_ids = z->elt_ids;
      if (_cell_ids != NULL) {
        for (cs_lnum_t j = 0; j < _n_cells; j++) {
          cell_ids[count] = _cell_ids[j];
          count++;
        }
      }
      else {
        for (cs_lnum_t j = 0; j < _n_cells; j++) {
          cell_ids[count] = j;
          count++;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tag cells of a given zone type
 *
 * The tag array should be initialized. The given tag_value is associted
 * to cells of the given zone type using a logical "or", so multiple flag
 * bits can be handled using the same array if necessary.
 *
 * \param[in]       zone_type_flag  zone types to tag
 * \param[in]       tag_value  tag value to add to cells of matching zones
 * \param[in, out]  tag        tag value for each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_tag_cell_type(int  zone_type_flag,
                             int  tag_value,
                             int  tag[])
{
  for (int i = 0; i < _n_zones; i++) {
    const cs_zone_t *z = _zones[i];
    if (z->type & zone_type_flag) {

      const cs_lnum_t _n_cells = z->n_elts;
      const cs_lnum_t *_cell_ids = z->elt_ids;
      if (_cell_ids != NULL) {
        for (cs_lnum_t j = 0; j < _n_cells; j++)
          tag[_cell_ids[j]] |= tag_value;
      }
      else {
        for (cs_lnum_t j = 0; j < _n_cells; j++)
          tag[j] = tag_value;
      }

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print volume zones information to listing file
 */
/*----------------------------------------------------------------------------*/

void
cs_volume_zone_print_info(void)
{
  bft_printf("\n");
  bft_printf(" --- Information on volume zones\n");

  /* Cell volumes and fluid cell volume */
  cs_real_t *cell_vol   = cs_glob_mesh_quantities->cell_vol;
  cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  cs_real_t *b_face_surf   = cs_glob_mesh_quantities->b_face_surf;
  cs_real_t *b_f_face_surf = cs_glob_mesh_quantities->b_f_face_surf;

  for (int i = 0; i < _n_zones; i++) {
    cs_zone_t *z = _zones[i];
    bft_printf(_("  Volume zone \"%s\"\n"
                 "    id              = %d\n"
                 "    Number of cells = %llu\n"
                 "    Volume          = %1.5g\n"
                 "    Center of gravity = (%1.5g, %1.5g, %1.5g)\n"),
               z->name,
               z->id,
               (unsigned long long)z->n_g_elts,
               z->measure,
               z->cog[0], z->cog[1], z->cog[2]);
    /* Only log fluid volumes when different to volumes */
    if (cell_f_vol != cell_vol && cell_f_vol != NULL)
      bft_printf(_("    Fluid volume    = %1.5g\n"),
                 z->f_measure);

    if (z->boundary_measure >= 0.) {
      bft_printf(_("    Surface         = %1.5g\n"),
                 z->f_boundary_measure);
      /* Only log fluid fluid when different to surface */
      if (b_f_face_surf != b_face_surf && b_f_face_surf != NULL)
        bft_printf(_("    Fluid surface   = %1.5g\n"),
                   z->f_boundary_measure);
    }
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
