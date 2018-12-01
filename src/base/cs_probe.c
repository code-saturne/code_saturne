/*============================================================================
 * Set of structures and functions to handle probes and profiles
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_point_location.h"

#include "cs_base.h"
#include "cs_map.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_selector.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_probe.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_probe.c

  \brief Probes and profiles management.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Predefined masks to set a flag for a probe set structure */

#define CS_PROBE_TRANSIENT   (1 << 0) //   1: locations of probes may change
#define CS_PROBE_BOUNDARY    (1 << 1) //   2: locations only on the border mesh
#define CS_PROBE_ON_CURVE    (1 << 2) //   4: locations are on a curve
#define CS_PROBE_AUTO_VAR    (1 << 3) //   8: automatic output of variables
#define CS_PROBE_OVERWRITE   (1 << 4) //  16: allow re-creation of probe set

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure to handle a set of probes */

struct _cs_probe_set_t {

  char            *name;          /* Associated name */

  int32_t          flags;         /* Metadata related to a set of probes */

  char            *sel_criter;    /* Selection criterion to filter entities
                                     before the location step */
  double           tolerance;     /* Criterion to define a threshold during
                                     the location step. This is a relative
                                     criterion. */
  cs_probe_snap_t  snap_mode;     /* Indicate how the positions of the
                                     probes are computed */

  cs_lnum_t     n_max_probes;   /* Number of probes initially requested */
  cs_lnum_t     n_probes;       /* Number of probes really used */
  cs_lnum_t     n_loc_probes;   /* Number of probes located on this rank */

  cs_real_3_t  *coords;         /* Coordinates of the set of probes
                                   Initially allocated to n_max_probes. */
  cs_real_t    *s_coords;       /* NULL or curvilinear coordinates
                                   for a profile */

  char        **labels;         /* List of labels for each probe (optional)
                                   By default, this is set to NULL */

  cs_probe_set_define_local_t  *p_define_func;   /* Advanced local definition
                                                    function */
  void                         *p_define_input;  /* Advanced local definition
                                                    input */

  int          *loc_id;         /* ids of probes located on local domain */

  cs_lnum_t    *elt_id;         /* Element ids where the probes have been
                                   located (size: n_loc_probes); -1 for
                                   unlocated probes assigned to local domain*/
  cs_lnum_t    *vtx_id;         /* Vertex ids closest to probes; -1 for
                                   unlocated probes assigned to local domain*/

  char         *located;        /* 1 for located probes, 0 for unlocated */

  /* User-defined writers associated to this set of probes */

  int           n_writers;      /* Number of writers (-1 if unset) */
  int          *writer_ids;     /* List of writer ids */

};

/* List of available keys for setting a set of probes */

typedef enum {

  PSETKEY_TRANSIENT_LOC,
  PSETKEY_BOUNDARY,
  PSETKEY_SELECT_CRIT,
  PSETKEY_TOLERANCE,
  PSETKEY_ERROR

} psetkey_t;

/*============================================================================
 *  Global static variables
 *============================================================================*/

/* Each structure stores a set of probes with some metadata */
static int  _n_probe_sets = 0;
static cs_probe_set_t  **_probe_set_array = NULL;

static const char _err_empty_pset[]
  = N_(" Stop execution since the given cs_probe_set_t structure is empty.\n"
       " Please check your settings.\n");
static const char _err_truefalse_key[]
  = N_(" Invalid value %s for setting key %s\n"
       " Valid choices are true or false.\n"
       " Please modify your setting.\n");

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_psetkey(psetkey_t  key)
{
  switch (key) {

  case PSETKEY_TRANSIENT_LOC:
    return "transient_location";
  case PSETKEY_BOUNDARY:
    return "boundary";
  case PSETKEY_SELECT_CRIT:
    return "selection_criteria";
  case PSETKEY_TOLERANCE:
    return "tolerance";
  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of a key.
 *         If not found, return a key error.
 *
 * \param[in] keyname    name of the key
 *
 * \return a psetkey_t
 */
/*----------------------------------------------------------------------------*/

static psetkey_t
_get_psetkey(const char  *keyname)
{
  psetkey_t  key = PSETKEY_ERROR;

  if (strcmp(keyname, "transient_location") == 0)
    key = PSETKEY_TRANSIENT_LOC;
  else if (strcmp(keyname, "boundary") == 0)
    key = PSETKEY_BOUNDARY;
  else if (strcmp(keyname, "selection_criteria") == 0)
    key = PSETKEY_SELECT_CRIT;
  else if (strcmp(keyname, "tolerance") == 0)
    key = PSETKEY_TOLERANCE;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy a label to a new string
 *
 * \param[in, out]  label    label to set
 *
 * \return   copy of label, or NULL
 */
/*----------------------------------------------------------------------------*/

inline static char *
_copy_label(const char  *name)
{
  char *label = NULL;

  if (name) {
    size_t  len = strlen(name) + 1;
    BFT_MALLOC(label, len, char);
    strcpy(label, name);
  }

  return label;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_probe_set_t structure
 *
 * \param[in, out]  pset          pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_probe_set_free(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return;

  BFT_FREE(pset->name);
  BFT_FREE(pset->coords);
  BFT_FREE(pset->sel_criter);
  BFT_FREE(pset->loc_id);
  BFT_FREE(pset->elt_id);
  BFT_FREE(pset->vtx_id);
  BFT_FREE(pset->located);

  if (pset->labels != NULL) {
    for (int i = 0; i < pset->n_probes; i++)
      BFT_FREE(pset->labels[i]);
    BFT_FREE(pset->labels);
  }

  if (pset->s_coords != NULL)
    BFT_FREE(pset->s_coords);

  if (pset->n_writers > 0)
    BFT_FREE(pset->writer_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize by default a new set of probes
 *
 * \param[in]  name          name of the probe set
 * \param[in]  n_max_probes  maximal number of probes
 *
 * \return a pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_probe_set_t *
_probe_set_create(const char    *name,
                  cs_lnum_t      n_max_probes)
{
  cs_probe_set_t *pset = cs_probe_set_get(name);

  /* Already defined */

  if (pset != NULL) {
    if (pset->flags & CS_PROBE_OVERWRITE)
      _probe_set_free(pset);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Error adding a new set of probes.\n"
                  " %s is already used as a name for a set of probes.\n"
                  " Please check your settings."), name);
  }

  /* Add a new set of probes */

  else {
    int pset_id = _n_probe_sets;

    _n_probe_sets++;
    BFT_REALLOC(_probe_set_array, _n_probe_sets, cs_probe_set_t *);
    BFT_MALLOC(pset, 1, cs_probe_set_t);
    _probe_set_array[pset_id] = pset;
  }

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pset->name, len, char);
  strncpy(pset->name, name, len);

  pset->flags = CS_PROBE_AUTO_VAR;
  pset->tolerance = 0.1;
  pset->sel_criter = NULL;

  pset->snap_mode = CS_PROBE_SNAP_NONE;

  pset->n_max_probes = n_max_probes;
  pset->n_probes = 0;
  pset->n_loc_probes = 0;

  BFT_MALLOC(pset->coords, n_max_probes, cs_real_3_t);
  pset->s_coords = NULL;
  pset->labels = NULL;

  pset->p_define_func = NULL;
  pset->p_define_input = NULL;

  pset->loc_id = NULL;
  pset->elt_id = NULL;
  pset->vtx_id = NULL;
  pset->located = NULL;

  pset->n_writers = -1;
  pset->writer_ids = NULL;

  return pset;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given name has the same name as the given set of probes
 *
 * \param[in]  pset        pointer to a cs_probe_set_t structure to test
 * \param[in]  ref_name    name to check
 *
 * \return true if the name of the set is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_check_probe_set_name(const cs_probe_set_t   *pset,
                      const char             *ref_name)
{
  if (pset == NULL)
    return false;

  int  reflen = strlen(ref_name);
  int  len = strlen(pset->name);

  if (reflen == len) {
    if (strcmp(ref_name, pset->name) == 0)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build probes based on a function-definition
 *
 * Note: if the p_define_input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * \param[in, out]  pset  pointer to a cs_probe_set_t structure to update
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_probe_set(cs_probe_set_t  *pset)
{
  assert(pset->p_define_func != NULL);

  pset->n_max_probes = 0;
  pset->n_probes = 0;
  pset->n_loc_probes = 0;

  BFT_FREE(pset->coords);
  BFT_FREE(pset->s_coords);

  cs_lnum_t     n_elts = 0;
  cs_real_3_t  *coords = NULL;
  cs_real_t    *s = NULL;

  pset->p_define_func(pset->p_define_input,
                      &n_elts,
                      &coords,
                      &s);

  pset->n_probes = n_elts;
  pset->coords = coords;
  pset->s_coords = s;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to a set of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_finalize(void)
{
  for (int i = 0; i < _n_probe_sets; i++) {
    cs_probe_set_t *pset = _probe_set_array[i];
    _probe_set_free(pset);
    BFT_FREE(pset);
  }

  _n_probe_sets = 0;
  BFT_FREE(_probe_set_array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of probe sets defined.
 *
 * \return the number of probe sets defined
 */
/*----------------------------------------------------------------------------*/

int
cs_probe_get_n_sets(void)
{
  return _n_probe_sets;
}

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
cs_probe_set_get(const char  *name)
{
  cs_probe_set_t  *pset = NULL;

  if (name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given name for this set of probes is empty."));

  /* Check if the given name is already used */
  for (int pset_id = 0; pset_id < _n_probe_sets; pset_id++) {
    if (_check_probe_set_name(_probe_set_array[pset_id], name)) {
      pset = _probe_set_array[pset_id];
      break;
    }
  }

  return pset;
}

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
cs_probe_set_get_by_id(int   pset_id)
{
  /* Check if the given name is already used */
  if (pset_id < 0 || pset_id >= _n_probe_sets)
    return  NULL;

  return _probe_set_array[pset_id];
}

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
cs_probe_set_get_name(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return  NULL;

  return pset->name;
}

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
                           int                    *writer_ids[])
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_pset);

  if (time_varying != NULL)
    *time_varying = (pset->flags & CS_PROBE_TRANSIENT) ? true : false;
  if (auto_variables != NULL)
    *auto_variables = (pset->flags & CS_PROBE_AUTO_VAR) ? true: false;
  if (on_curve != NULL)
    *on_curve = (pset->flags & CS_PROBE_ON_CURVE) ? true : false;
  if (on_boundary != NULL)
    *on_boundary = (pset->flags & CS_PROBE_BOUNDARY) ? true : false;

  if (n_writers != NULL)
    *n_writers = pset->n_writers;
  if (writer_ids != NULL)
    *writer_ids = pset->writer_ids;
}

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
cs_probe_set_get_location_criteria(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return  NULL;

  return pset->sel_criter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new set of probes.
 *
 * \param[in]   name   name of the set of probes
 *
 * \return a pointer to a new allocated cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_probe_set_t *
cs_probe_set_create(const char  *name)
{
  /* Default initialization of a set of probes (max number is set to 4 by
     default but a realloc is available if the max. number of probes is
     reached) */

  cs_probe_set_t *pset = _probe_set_create(name, 4);

  return  pset;
}

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
                       const char         *label)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_pset);

  cs_lnum_t  point_id = pset->n_probes;

  pset->n_probes++;

  if (point_id >= pset->n_max_probes) { /* Reallocate arrays */
    pset->n_max_probes *= 2;
    BFT_REALLOC(pset->coords, pset->n_max_probes, cs_real_3_t);
    if (pset->labels != NULL)
      BFT_REALLOC(pset->labels, pset->n_max_probes, char *);
  }

  /* Set coordinates */
  pset->coords[point_id][0] = x;
  pset->coords[point_id][1] = y;
  pset->coords[point_id][2] = z;

  if (label != NULL) { /* Manage the label */
    if (pset->labels == NULL)
      BFT_MALLOC(pset->labels, pset->n_max_probes, char *);

    /* Copy label */
    pset->labels[point_id] = _copy_label(label);
  }
}

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
                               const char         **labels)
{
  cs_probe_set_t  *pset = _probe_set_create(name, n_probes);

  pset->n_probes = n_probes;

  /* Coordinates */
  for (int i = 0; i < n_probes; i++) {
    pset->coords[i][0] = coords[i][0];
    pset->coords[i][1] = coords[i][1];
    pset->coords[i][2] = coords[i][2];
  }

  /* Copy labels */
  if (labels != NULL) {
    BFT_MALLOC(pset->labels, n_probes, char *);
    for (int i = 0; i < n_probes; i++)
      pset->labels[i] = _copy_label(labels[i]);
  }

  return  pset;
}

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
                                 const cs_real_t    end_coords[3])
{
  cs_probe_set_t  *pset = _probe_set_create(name, n_probes);

  pset->n_probes = n_probes;
  pset->flags |= CS_PROBE_ON_CURVE;
  if (pset->flags & CS_PROBE_AUTO_VAR)
    pset->flags -= CS_PROBE_AUTO_VAR;

  BFT_MALLOC(pset->s_coords, n_probes, cs_real_t);

  /* 2 probes are already defined (the starting and ending points)
     Define the additional probes */
  cs_real_t  distance;
  cs_real_3_t  unitv, delta_vect;

  cs_math_3_length_unitv(start_coords, end_coords, &distance, unitv);

  const double  delta = distance / (n_probes - 1);
  for (int k = 0; k < 3; k++)
    delta_vect[k] = delta*unitv[k];

  /* Set the starting probe */
  pset->s_coords[0] = 0;
  for (int k = 0; k < 3; k++)
    pset->coords[0][k] = start_coords[k];

  /* Set additional probes */
  for (int i = 1; i < n_probes - 1; i++) {

    pset->s_coords[i] = pset->s_coords[i-1] + delta;
    for (int k = 0; k < 3; k++)
      pset->coords[i][k] = pset->coords[i-1][k] + delta_vect[k];

  }

  /* Set the ending probe */
  pset->s_coords[n_probes-1] = distance;
  for (int k = 0; k < 3; k++)
    pset->coords[n_probes-1][k] = end_coords[k];

  return  pset;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new set of probes from rank-local definition function.
 *
 * The local definition function given by the p_define_func pointer
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
                               void                         *p_define_input)
{
  cs_probe_set_t  *pset = _probe_set_create(name, 0);

  pset->flags |= CS_PROBE_ON_CURVE;
  if (pset->flags & CS_PROBE_AUTO_VAR)
    pset->flags -= CS_PROBE_AUTO_VAR;

  pset->p_define_func = p_define_func;
  pset->p_define_input = p_define_input;

  return  pset;
}

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
cs_probe_set_allow_overwrite(const char  *name)
{
  cs_probe_set_t *pset = cs_probe_set_get(name);

  if (pset != NULL)
    pset->flags = pset->flags | CS_PROBE_OVERWRITE;
}

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
                               const int        *writer_ids)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pset));

  if (pset->n_writers < 0)
    pset->n_writers = 0;

  int  n_init_writers = pset->n_writers;
  pset->n_writers += n_writers;
  BFT_REALLOC(pset->writer_ids, pset->n_writers, int);

  for (int i = n_init_writers, j = 0; i < pset->n_writers; i++, j++)
    pset->writer_ids[i] = writer_ids[j];
}

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
                      bool              mode)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pset));

  if (mode == false) {
    if (pset->flags & CS_PROBE_AUTO_VAR)
      pset->flags -= CS_PROBE_AUTO_VAR;
  }
  else
    pset->flags |= CS_PROBE_AUTO_VAR;
}

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
                       cs_probe_snap_t   snap_mode)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pset));

  pset->snap_mode = snap_mode;
}

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
                    const char       *keyval)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pset));

  psetkey_t  key = _get_psetkey(keyname);

  if (key == PSETKEY_ERROR) {

    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < PSETKEY_ERROR; i++) {
      bft_printf("%s ", _print_psetkey(i));
      if (i > 0 && i % 4 == 0) {
        bft_printf("\n");
        if (i + 1 < PSETKEY_ERROR)
          bft_printf("\t");
      }
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for probe options %s.\n"
                " Please read run_solver.log for more details and"
                " modify your settings."), pset->name);

  } /* Error message */

  switch(key) {

  case PSETKEY_BOUNDARY:
    if (strcmp(keyval, "true") == 0)
      pset->flags |= CS_PROBE_BOUNDARY;
    else if (strcmp(keyval, "false") == 0) { // remove the flags if it is set
      if (pset->flags & CS_PROBE_BOUNDARY)
        pset->flags ^= CS_PROBE_BOUNDARY;
    }
    else
      bft_error(__FILE__, __LINE__, 0, _err_truefalse_key, keyval, keyname);
    break;

  case PSETKEY_SELECT_CRIT:
    {
      int  len = strlen(keyval) + 1;
      BFT_MALLOC(pset->sel_criter, len, char);
      strncpy(pset->sel_criter, keyval, len);
    }
    break;

  case PSETKEY_TRANSIENT_LOC:
    if (strcmp(keyval, "true") == 0)
      pset->flags |= CS_PROBE_TRANSIENT;
    else if (strcmp(keyval, "false") == 0) { // remove the flag if it is set
      if (pset->flags & CS_PROBE_TRANSIENT)
        pset->flags ^= CS_PROBE_TRANSIENT;
    }
    else
      bft_error(__FILE__, __LINE__, 0, _err_truefalse_key, keyval, keyname);
    break;

  case PSETKEY_TOLERANCE:
    pset->tolerance = atof(keyval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);
    break;

  } /* Switch on keys */

}

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
                    const fvm_nodal_t  *location_mesh)
{
  if (pset == NULL)
    return;

  bool first_location = false;
  const double  tolerance_base = 0.;
  const cs_mesh_t  *mesh = cs_glob_mesh;
  fvm_nodal_t *_location_mesh = NULL;

  const bool  on_boundary = (pset->flags & CS_PROBE_BOUNDARY) ? true : false;

  /* Build in local case */

  if (pset->p_define_func != NULL)
    _build_local_probe_set(pset);

  /* Allocate on first pass */

  if (pset->located == NULL) {
    BFT_MALLOC(pset->located, pset->n_probes, char);
    first_location = true;
  }

  /* Reallocate on all passes, in case local sizes change */

  BFT_REALLOC(pset->loc_id, pset->n_probes, int);
  BFT_REALLOC(pset->elt_id, pset->n_probes, cs_lnum_t);
  BFT_FREE(pset->vtx_id);

  if (location_mesh == NULL) {

    cs_lnum_t  n_select_elements = 0;
    cs_lnum_t  *selected_elements = NULL;

    if (on_boundary) { /* Deal with the surface mesh related to the boundary */

      n_select_elements = mesh->n_b_faces;
      if (pset->sel_criter != NULL) {
        if (strcmp(pset->sel_criter, "all[]")) {
          BFT_MALLOC(selected_elements, mesh->n_b_faces, cs_lnum_t);
          cs_selector_get_b_face_num_list(pset->sel_criter,
                                          &n_select_elements, selected_elements);
        }
      } /* Need to define a list of faces ? */

      _location_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                      "probe_location_mesh",
                                                      false, // no family info
                                                      0,     // interior faces
                                                      n_select_elements,
                                                      NULL,  // interior faces
                                                      selected_elements);

    }
    else { /* Deal with the volumic mesh */

      n_select_elements = mesh->n_cells;
      if (pset->sel_criter != NULL) {
        if (strcmp(pset->sel_criter, "all[]")) {
          BFT_MALLOC(selected_elements, mesh->n_cells, cs_lnum_t);
          cs_selector_get_cell_num_list(pset->sel_criter,
                                        &n_select_elements, selected_elements);
        }
      } /* Need to define a list of cells ? */

      _location_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                                      "probe_location_mesh",
                                                      false, // no family info
                                                      n_select_elements,
                                                      selected_elements);

    } /* volumic or surfacic mesh */

    if (selected_elements != NULL)
      BFT_FREE(selected_elements);

    location_mesh = _location_mesh;
  }

  /* Locate probes on this location mesh */

  float *distance;

  BFT_MALLOC(distance, pset->n_probes, float);

  for (int i = 0; i < pset->n_probes; i++) {
    pset->elt_id[i] = -1;
    distance[i] = -1.0;
  }

  fvm_point_location_nodal(location_mesh,
                           tolerance_base,
                           pset->tolerance,
                           0, /* locate_on_parents */
                           pset->n_probes,
                           NULL, /* point_tag */
                           (const cs_coord_t *)(pset->coords),
                           pset->elt_id,
                           distance);

  for (int i = 0; i < pset->n_probes; i++) {
    if (pset->elt_id[i] < 0) /* Not found */
      distance[i] = HUGE_VAL;
  }

  /* Warning if points have not beeen located */
  cs_gnum_t  n_unlocated_probes = 0;

  int n_loc_probes = 0;

  if (cs_glob_n_ranks == 1 || pset->p_define_func != NULL) {
    for (int i = 0; i < pset->n_probes; i++) {
      if (distance[i] >= 0.5*HUGE_VAL) {
        pset->located[i] = 0;
        n_unlocated_probes++;
      }
      else {
        pset->loc_id[n_loc_probes] = i;
        pset->elt_id[n_loc_probes] = pset->elt_id[i];
        pset->located[i] = 1;
        n_loc_probes += 1;
      }
    }
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1 && pset->p_define_func == NULL) {

    cs_double_int_t  *gmin_loc = NULL, *loc = NULL;

    BFT_MALLOC(gmin_loc, pset->n_probes, cs_double_int_t);
    BFT_MALLOC(loc, pset->n_probes, cs_double_int_t);

    for (int i = 0; i < pset->n_probes; i++) {
      gmin_loc[i].id = loc[i].id = cs_glob_rank_id;
      gmin_loc[i].val = loc[i].val = distance[i];
    }

    MPI_Allreduce(loc, gmin_loc, pset->n_probes, MPI_DOUBLE_INT, MPI_MINLOC,
                  cs_glob_mpi_comm);

    for (int i = 0; i < pset->n_probes; i++) {
      if (gmin_loc[i].val >= 0.5*HUGE_VAL) {
        pset->located[i] = 0;
        n_unlocated_probes++;
      }
      else {
        pset->located[i] = 1;
        if (gmin_loc[i].id == cs_glob_rank_id) {
          pset->loc_id[n_loc_probes] = i;
          pset->elt_id[n_loc_probes] = pset->elt_id[i];
          n_loc_probes += 1;
        }
      }
    }

    BFT_FREE(gmin_loc);
    BFT_FREE(loc);
  }
#endif

  BFT_FREE(distance);

  if (n_unlocated_probes > 0 && first_location) {
    bft_printf(_("\n Warning: probe set \"%s\"\n"
                 "   %lu (of %d) probes are not located"
                 " on the associated mesh:\n"),
               pset->name, (unsigned long)n_unlocated_probes,
               pset->n_probes);
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->located[i] == 0) {
        if (pset->labels == NULL)
          bft_printf(_("    %2d ([%8.3e, %8.3e, %8.3e])\n"),
                     i+1,
                     pset->coords[i][0],
                     pset->coords[i][1],
                     pset->coords[i][2]);
        else
          bft_printf(_("    %s ([%8.3e, %8.3e, %8.3e])\n"),
                     pset->labels[i],
                     pset->coords[i][0],
                     pset->coords[i][1],
                     pset->coords[i][2]);
      }
    }
  }

  pset->n_loc_probes = n_loc_probes;

  if (n_unlocated_probes) {

    /* For a transient mesh, non-located probes may vary, so we need to maintain
       them in the structure to avoid having a varying number of columns in a
       plot. We add those locations to the last rank */

    if (pset->flags & CS_PROBE_TRANSIENT) {
      if (cs_glob_rank_id == cs_glob_n_ranks - 1 || cs_glob_n_ranks == 1)
        pset->n_loc_probes += n_unlocated_probes;
    }

    /* For a fixed mesh, if we do not have labels, we add labels based on the
       probe global ids so as to maintain probe ids in the output metadata. */

    else if (   pset->labels == NULL
             && (pset->flags & CS_PROBE_ON_CURVE) == false) {
      BFT_MALLOC(pset->labels, pset->n_probes, char *);
      char label[16];
      for (int i = 0; i < pset->n_probes; i++) {
        snprintf(label, 15, "%d", i+1); label[15] = '\0';
        pset->labels[i] = _copy_label(label);
      }
    }

  }

  BFT_REALLOC(pset->loc_id, pset->n_loc_probes, cs_lnum_t);
  BFT_REALLOC(pset->elt_id, pset->n_loc_probes, cs_lnum_t);
  BFT_MALLOC(pset->vtx_id, pset->n_loc_probes, cs_lnum_t);

  /* Now also locate relative to closest vertices and update
     element num relative to parents */

  cs_coord_3_t  *probe_coords = NULL;
  BFT_MALLOC(probe_coords, pset->n_loc_probes, cs_coord_3_t);
  for (int i = 0; i < n_loc_probes; i++) {
    int j = pset->loc_id[i];
    for (int k = 0; k < 3; k++)
      probe_coords[i][k] = pset->coords[j][k];
  }

  fvm_point_location_closest_vertex(location_mesh,
                                    1, /* locate on parents */
                                    n_loc_probes,
                                    (const cs_coord_t *)(probe_coords),
                                    pset->elt_id,
                                    pset->vtx_id);

  BFT_FREE(probe_coords);

  /* Now switch to 0-based location */

  for (int i = 0; i < n_loc_probes; i++) {
    if (pset->elt_id[i] > -1) {
      pset->elt_id[i] -= 1;
      pset->vtx_id[i] -= 1;
    }
  }

  if (_location_mesh != NULL)
    _location_mesh = fvm_nodal_destroy(_location_mesh);

  /* Finally, effectively add non-located probe info when required */

  if (pset->n_loc_probes > n_loc_probes) {
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->located[i] == 0) {
        pset->loc_id[n_loc_probes] = i;
        pset->elt_id[n_loc_probes] = -1;
        pset->vtx_id[n_loc_probes] = -1;
        n_loc_probes++;
      }
    }
  }

  assert(n_loc_probes == pset->n_loc_probes);
}

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
                         const char       *mesh_name)
{
  if (pset == NULL)
    return  NULL;

  cs_gnum_t  *global_num = NULL;
  cs_coord_3_t  *probe_coords = NULL;
  fvm_nodal_t  *exp_mesh = fvm_nodal_create(mesh_name, 3);

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;

  const cs_real_t *centers = mq->cell_cen;

  if (pset->flags & CS_PROBE_BOUNDARY) centers = mq->b_face_cog;

  BFT_MALLOC(probe_coords, pset->n_loc_probes, cs_coord_3_t);
  BFT_MALLOC(global_num, pset->n_loc_probes, cs_gnum_t);

  /* Build the final list of probe coordinates */

  cs_real_t  max_distance = 0.;

  for (int i = 0; i < pset->n_loc_probes; i++) {
    int j = pset->loc_id[i];

    for (int k = 0; k < 3; k++)
      probe_coords[i][k] = pset->coords[j][k];

    global_num[i] = j+1;

    if (pset->elt_id[i] > -1) {
      const cs_real_t  *elt_coords = centers + pset->elt_id[i]*3;
      cs_real_t v[3];
      for (int k = 0; k < 3; k++)
        v[k] = elt_coords[k] - pset->coords[j][k];
      max_distance = fmax(max_distance, cs_math_3_square_norm(v));
    }
  }

  cs_real_t gmax_distance = max_distance;

  /* Handle snap mode if active */

  if (pset->snap_mode == CS_PROBE_SNAP_ELT_CENTER) {
    for (int i = 0; i < pset->n_loc_probes; i++) {
      int j = pset->loc_id[i];
      if (pset->elt_id[i] > -1) {
        const cs_real_t  *elt_coords = centers + pset->elt_id[i]*3;
        for (int k = 0; k < 3; k++)
          pset->coords[j][k] = elt_coords[k];
      }
    }
  }
  else if (pset->snap_mode == CS_PROBE_SNAP_VERTEX) {
    for (int i = 0; i < pset->n_loc_probes; i++) {
      int j = pset->loc_id[i];
      if (pset->vtx_id[i] > -1) {
        const cs_real_t  *vtx_coords = m->vtx_coord + pset->vtx_id[i]*3;
        for (int k = 0; k < 3; k++)
          pset->coords[j][k] = vtx_coords[k];
      }
    }
  }

  /* Update the probe set structure */

  fvm_nodal_define_vertex_list(exp_mesh, pset->n_loc_probes, NULL);
  fvm_nodal_transfer_vertices(exp_mesh, (cs_coord_t *)probe_coords);

  /* Set a global numbering if needed */

  if (pset->p_define_func != NULL) {
    cs_real_t *s;
    BFT_MALLOC(s, pset->n_loc_probes, cs_real_t);
    for (int i = 0; i < pset->n_loc_probes; i++) {
      int j = pset->loc_id[i];
      s[i] = pset->s_coords[j];
    }
    fvm_io_num_t *vtx_io_num
      = fvm_io_num_create_from_real(pset->s_coords, pset->n_loc_probes);
    BFT_FREE(s);
    fvm_nodal_transfer_vertex_io_num(exp_mesh, &vtx_io_num);
  }
  else if (cs_glob_n_ranks > 1) {
    fvm_nodal_init_io_num(exp_mesh, global_num, 0); // 0 = vertices

#if defined(HAVE_MPI)
    MPI_Reduce(&max_distance, &gmax_distance, 1, MPI_DOUBLE, MPI_MAX, 0,
               cs_glob_mpi_comm);
#endif
  }

  if (! (   pset->flags & CS_PROBE_ON_CURVE
         || pset->flags & CS_PROBE_TRANSIENT))
    bft_printf(_("\n Probe set: \"%s\":\n"
                 "   maximum distance between cell centers and"
                 " requested coordinates:"
                 " %5.3e\n"), pset->name, gmax_distance);

  BFT_FREE(global_num);

  /* Add global labels */

  if (pset->labels != NULL) {
    char **g_labels;
    int ngl = fvm_nodal_get_n_g_vertices(exp_mesh);
    BFT_MALLOC(g_labels, ngl, char *);

    int j = 0;
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->located[i] != 0)
        g_labels[j++] = _copy_label(pset->labels[i]);
    }
    assert(j == ngl);
    fvm_nodal_transfer_global_vertex_labels(exp_mesh, g_labels);

  }

  /* probe_coords is managed by exp_mesh */

  return exp_mesh;
}

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
                                   const char       *mesh_name)
{
  if (pset == NULL)
    return  NULL;

  int  n_exp_probes = 0;
  cs_gnum_t  *global_num = NULL;
  cs_coord_3_t  *probe_coords = NULL;
  fvm_nodal_t  *exp_mesh = fvm_nodal_create(mesh_name, 3);

  BFT_MALLOC(probe_coords, pset->n_probes, cs_coord_3_t);
  BFT_MALLOC(global_num, pset->n_loc_probes, cs_gnum_t);

  /* Build the final list of probe coordinates */

  for (int i = 0; i < pset->n_probes; i++) {
    if (pset->located[i] == 0) {
      for (int k = 0; k < 3; k++)
        probe_coords[n_exp_probes][k] = pset->coords[i][k];
      global_num[n_exp_probes] = i+1;
      n_exp_probes++;
    }
  }

  /* Update the probe set structure */

  fvm_nodal_define_vertex_list(exp_mesh, n_exp_probes, NULL);
  fvm_nodal_transfer_vertices(exp_mesh, (cs_coord_t *)probe_coords);

  /* Set a global numbering if needed */

  if (pset->p_define_func != NULL) {
    cs_real_t *s;
    BFT_MALLOC(s, pset->n_probes, cs_real_t);
    int j = 0;
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->located[i] == 0)
        s[j++] = pset->s_coords[i];
    }
    fvm_io_num_t *vtx_io_num
      = fvm_io_num_create_from_real(pset->s_coords, j);
    BFT_FREE(s);
    fvm_nodal_transfer_vertex_io_num(exp_mesh, &vtx_io_num);
  }
  else if (cs_glob_n_ranks > 1)
    fvm_nodal_init_io_num(exp_mesh, global_num, 0); // 0 = vertices

  BFT_FREE(global_num);

  /* Add global labels */

  if (pset->labels != NULL) {
    char **g_labels;
    int ngl = fvm_nodal_get_n_g_vertices(exp_mesh);
    BFT_MALLOC(g_labels, ngl, char *);

    int j = 0;
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->located[i] == 0)
        g_labels[j++] = _copy_label(pset->labels[i]);
    }
    assert(j == ngl);
    fvm_nodal_transfer_global_vertex_labels(exp_mesh, g_labels);
  }

  /* probe_coords is managed by exp_mesh */

  return exp_mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_probe_set_t structure.
 *
 * \param[in]  pset    pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_dump(const cs_probe_set_t   *pset)
{
  bft_printf("\n\n Dump cs_probe_set_t structure %p\n", (const void *)pset);

  if (pset == NULL)
    return;

  bft_printf(" name:                %s\n"
             " flags:               %d\n"
             " location criteria:   %s\n"
             " tolerance:           %5.3e\n",
             pset->name, pset->flags, pset->sel_criter,
             pset->tolerance);

  if (pset->sel_criter != NULL)
    bft_printf(" selection:  %s\n", pset->sel_criter);

  bft_printf(" n_probes:   %d; %d; %d (locally located; defined; max.)\n",
             pset->n_loc_probes, pset->n_probes, pset->n_max_probes);

  for (int i = 0; i < pset->n_probes; i++) {

    bft_printf(" %4d | %-5.3e %-5.3e %-5.3e |", i,
               pset->coords[i][0], pset->coords[i][1], pset->coords[i][2]);

    if (pset->s_coords != NULL)
      bft_printf(" %5.3e |", pset->s_coords[i]);
    if (pset->elt_id != NULL && pset->located != NULL)
      bft_printf(" %6d | %c |", pset->elt_id[i], pset->located[i]);
    if (pset->labels != NULL)
      if (pset->labels[i] != NULL)
        bft_printf(" %s", pset->labels[i]);
    bft_printf("\n");

  }

}

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
                         cs_real_3_t            *coords[])
{
  if (pset == NULL)
    return;

  /* Return pointers */

  if (snap_mode != NULL)
    *snap_mode = pset->snap_mode;

  if (n_probes != NULL)
    *n_probes = pset->n_probes;

  if (coords != NULL)
    *coords = pset->coords;
}

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
cs_probe_set_get_n_local(const cs_probe_set_t   *pset)
{
  int retval = 0;

  if (pset != NULL)
    retval = pset->n_loc_probes;

  return retval;
}

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
cs_probe_set_get_curvilinear_abscissa(const cs_probe_set_t   *pset)
{
  const cs_real_t *retval = NULL;

  if (pset == NULL)
    return retval;
  else
    return pset->s_coords;
}

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
                         int                    mesh_location_id)
{
  const cs_lnum_t *retval = NULL;

  if (pset == NULL)
    return retval;

  bool on_boundary = (pset->flags & CS_PROBE_BOUNDARY) ? true : false;

  if (mesh_location_id == CS_MESH_LOCATION_CELLS && on_boundary == false)
    retval = pset->elt_id;
  else if (mesh_location_id == CS_MESH_LOCATION_BOUNDARY_FACES && on_boundary)
    retval = pset->elt_id;

  else if (CS_MESH_LOCATION_VERTICES)
    retval = pset->vtx_id;

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
