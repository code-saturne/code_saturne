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
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
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

#define CS_PROBE_ACTIVATED   (1 << 0) //   1: state = activated
#define CS_PROBE_PREDEFINED  (1 << 1) //   2: predefined (ie not user)
#define CS_PROBE_TRANSIENT   (1 << 2) //   4: locations of probes may change
#define CS_PROBE_INTERPOLATE (1 << 3) //   8: switch on the interpolation step
#define CS_PROBE_BOUNDARY    (1 << 4) //  16: locations only on the border mesh
#define CS_PROBE_PROFILE     (1 << 5) //  32: probe set is a profile
#define CS_PROBE_POST_DIST   (1 << 6) //  64: post distance to exact location
#define CS_PROBE_USER_VAR    (1 << 7) // 128: post user-defined set of variables

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure to handle a set of probes */

struct _cs_probe_set_t {

  char            *name;          /* Associated name */

  cs_flag_t        flag;          /* Metadata related to a set of probes */

  char            *sel_criter;    /* Selection criterion to filter entities
                                     before the location step */
  double           tolerance;     /* Criterion to define a threshold during
                                     the location step. This is a relative
                                     criterion. */
  cs_probe_mode_t  mode;          /* Indicate how the value of the probes are
                                     computed */

  cs_lnum_t     n_max_probes;   /* Number of probes initially requested */
  cs_lnum_t     n_probes;       /* Number of probes really used */
  cs_lnum_t     n_loc_probes;   /* Number of probes located on this rank */

  cs_real_t    *s_coords;       /* NULL or curvilinear abscissa for a profile */
  cs_real_t    *coords;         /* Coordinates of the set of probes
                                   Initially allocated to n_max_probes. */
  char        **labels;         /* List of labels for each probe (optional)
                                   By default, this is set to NULL */
  cs_lnum_t    *entity_num;     /* Ids related to the entity where the probe
                                   has been located.
                                   -1 if not found on this rank */
  float        *distances;      /* Distance between the probe coordinates and
                                   the entity if found */

  fvm_nodal_t  *location_mesh;  /* Nodal mesh used for locating probes */

  /* User-defined writers associated to this set of probes */

  int           n_writers;      /* Number of writers */
  int          *writer_ids;     /* List of writer ids */

};

/* List of available keys for setting a set of probes */

typedef enum {

  PSETKEY_ACTIVATED,
  PSETKEY_BOUNDARY,
  PSETKEY_MODE,
  //PSETKEY_OUTPUT_DISTANCE, // TODO
  PSETKEY_PROFILE,
  PSETKEY_SELECT_CRIT,
  PSETKEY_TOLERANCE,
  PSETKEY_TRANSIENT_LOC,
  //PSETKEY_USER_VAR, //TODO
  PSETKEY_ERROR

} psetkey_t;

/*============================================================================
 *  Global static variables
 *============================================================================*/

/* Each structure stores a set of probes with some metadata */
static int  _n_probe_sets = 0;
static cs_probe_set_t   *_probe_set_array = NULL;

static const char _err_empty_pset[] =
  N_(" Stop execution since the given cs_probe_set_t structure is empty.\n"
     " Please check your settings.\n");
static const char _err_truefalse_key[] =
  N_(" Invalid value %s for setting key %s\n"
     " Valid choices are true or false.\n"
     " Please modify your setting.\n");

static const char
cs_probe_mode_name[CS_PROBE_N_MODES][CS_BASE_STRING_LEN]=
  { N_("exact"),
    N_("nearest cell center"),
    N_("nearest vertex") };

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

  case PSETKEY_ACTIVATED:
    return "activated";
  case PSETKEY_BOUNDARY:
    return "boundary";
  case PSETKEY_MODE:
    return "mode";
  case PSETKEY_PROFILE:
    return "profile";
  case PSETKEY_SELECT_CRIT:
    return "selection_criterion";
  case PSETKEY_TOLERANCE:
    return "tolerance";
  case PSETKEY_TRANSIENT_LOC:
    return "moving_probes";
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

  if (strcmp(keyname, "activation") == 0)
    key = PSETKEY_ACTIVATED;
  else if (strcmp(keyname, "boundary") == 0)
    key = PSETKEY_BOUNDARY;
  else if (strcmp(keyname, "mode") == 0)
    key = PSETKEY_MODE;
  else if (strcmp(keyname, "profile") == 0)
    key = PSETKEY_PROFILE;
  else if (strcmp(keyname, "selection_criterion") == 0)
    key = PSETKEY_SELECT_CRIT;
  else if (strcmp(keyname, "tolerance") == 0)
    key = PSETKEY_TOLERANCE;
  else if (strcmp(keyname, "moving_probes") == 0)
    key = PSETKEY_TRANSIENT_LOC;

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
  /* Add a new set of probes */
  int  pset_id = _n_probe_sets;

  _n_probe_sets++;
  BFT_REALLOC(_probe_set_array, _n_probe_sets, cs_probe_set_t);

  cs_probe_set_t  *pset = _probe_set_array + pset_id;

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pset->name, len, char);
  strncpy(pset->name, name, len);

  pset->flag = CS_PROBE_ACTIVATED;
  pset->mode = CS_PROBE_MODE_NEAREST_CELL_CENTER;
  pset->tolerance = 0.1;
  pset->sel_criter = NULL;

  pset->n_max_probes = n_max_probes;
  pset->n_probes = 0;
  pset->n_loc_probes = 0;

  BFT_MALLOC(pset->coords, 3*n_max_probes, cs_real_t);

  pset->labels = NULL;
  pset->s_coords = NULL;
  pset->entity_num = NULL;
  pset->distances = NULL;

  pset->location_mesh = NULL;

  pset->n_writers = 0;
  pset->writer_ids = NULL;

  return pset;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_probe_set_t structure
 *
 * \param[in, out]  pset          pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_probe_set(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return;

  BFT_FREE(pset->name);
  BFT_FREE(pset->coords);
  BFT_FREE(pset->sel_criter);
  BFT_FREE(pset->entity_num);
  BFT_FREE(pset->distances);

  if (pset->labels != NULL) {
    for (int i = 0; i < pset->n_probes; i++)
      BFT_FREE(pset->labels[i]);
    BFT_FREE(pset->labels);
  }

  if (pset->s_coords != NULL)
    BFT_FREE(pset->s_coords);

  pset->location_mesh = fvm_nodal_destroy(pset->location_mesh);

  if (pset->n_writers > 0)
    BFT_FREE(pset->writer_ids);
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of probe sets defined
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
 * \brief  Retrieve a cs_probe_set_t structure
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
    if (_check_probe_set_name(_probe_set_array + pset_id, name)) {
      pset = _probe_set_array + pset_id;
      break;
    }
  }

  return pset;
}

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
cs_probe_set_get_by_id(int   pset_id)
{
  /* Check if the given name is already used */
  if (pset_id < 0 || pset_id >= _n_probe_sets)
    return  NULL;

  return _probe_set_array + pset_id;
}

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
cs_probe_set_get_name(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return  NULL;

  return pset->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve information useful for the postprocessing step
 *
 * \param[in]  pset          pointer to a cs_probe_set_t structure
 * \param[out] time_varying  true if probe coords may change during computation
 * \param[out] is_profile    true if probe set is related to a profile
 * \param[out] on_boundary   true if probes are located on boundary
 * \param[out] is_automatic  true if set of variables to output is predefined
 * \param[out] n_writers     number of associated  user-defined writers
 * \param[out] writer_ids    pointer to a list of writer ids
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_get_post_info(const cs_probe_set_t   *pset,
                           bool                   *time_varying,
                           bool                   *is_profile,
                           bool                   *on_boundary,
                           bool                   *is_automatic,
                           int                    *n_writers,
                           int                    *writer_ids[])
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_pset);

  *time_varying = (pset->flag & CS_PROBE_TRANSIENT) ? true : false;
  *is_automatic = (pset->flag & CS_PROBE_USER_VAR) ? false : true;
  *is_profile = (pset->flag & CS_PROBE_PROFILE) ? true : false;
  *on_boundary = (pset->flag & CS_PROBE_BOUNDARY) ? true : false;

  *n_writers = pset->n_writers;
  *writer_ids = pset->writer_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a set of monitoring probes has been defined among all the
 *         probe sets
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_probe_set_have_monitoring(void)
{
  bool  have_monitoring = false;

  for (int i = 0; i < _n_probe_sets; i++) {
    cs_probe_set_t  *pset = _probe_set_array + i;
    if (!(pset->flag & CS_PROBE_PROFILE))
      have_monitoring = true;
  }

  return have_monitoring;
}

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
cs_probe_set_create(const char  *name)
{
  cs_probe_set_t *pset = cs_probe_set_get(name);

  if (pset != NULL) /* Already defined */
    bft_error(__FILE__, __LINE__, 0,
              _(" Error adding a new set of probes.\n"
                " %s is already used as a name for a set of probes.\n"
                " Please check your settings."), name);

  /* Default initialization of a set of probes (max number is set to 4 by
     default but a realloc is available is the max. number of probes is
     reached) */

  pset = _probe_set_create(name, 4);

  return  pset;
}

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
                       const char         *label)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_pset);

  cs_lnum_t  point_id = pset->n_probes;

  pset->n_probes++;

  if (point_id >= pset->n_max_probes) { /* Reallocate arrays */
    pset->n_max_probes *= 2;
    BFT_REALLOC(pset->coords, 3*pset->n_max_probes, cs_real_t);
    if (pset->labels != NULL)
      BFT_REALLOC(pset->labels, pset->n_max_probes, char *);
  }

  /* Set coordinates */
  pset->coords[3*point_id]     = xyz[0];
  pset->coords[3*point_id + 1] = xyz[1];
  pset->coords[3*point_id + 2] = xyz[2];

  if (label != NULL) { /* Manage the label */
    if (pset->labels == NULL)
      BFT_MALLOC(pset->labels, pset->n_max_probes, char *);

    /* Copy label */
    pset->labels[point_id] = _copy_label(label);
  }
}

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
                               const char        **labels)
{
  cs_probe_set_t  *pset = cs_probe_set_get(name);

  if (pset != NULL) /* Already defined */
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop adding a new set of probes.\n"
                " %s is already used as a name for a set of probes.\n"
                " Please check your settings."), name);

  pset = _probe_set_create(name, n_probes);
  pset->n_probes = n_probes;

  /* Coordinates */
  for (int i = 0; i < n_probes; i++) {
    const cs_lnum_t shift = 3*i;
    pset->coords[shift  ] = coords[shift];
    pset->coords[shift+1] = coords[shift+1];
    pset->coords[shift+2] = coords[shift+2];
  }

  /* Copy labels */
  if (labels != NULL) {

    BFT_MALLOC(pset->labels, n_probes, char *);
    for (int i = 0; i < n_probes; i++)
      pset->labels[i] = _copy_label(labels[i]);

  } // labels != NULL

  return  pset;
}

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
                                 const cs_real_t    end_coords[3])
{
  cs_probe_set_t  *pset = cs_probe_set_get(name);

  if (pset != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop adding a new set of probes.\n"
                " %s is already used as a name for a set of probes.\n"
                " Please check your settings."), name);

  pset = _probe_set_create(name, n_probes);
  pset->n_probes = n_probes;
  pset->flag |= CS_PROBE_PROFILE;

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
    pset->coords[k] = start_coords[k];

  /* Set additional probes */
  for (int i = 1; i < n_probes - 1; i++) {

    pset->s_coords[i] = pset->s_coords[i-1] + delta;
    for (int k = 0; k < 3; k++)
      pset->coords[3*i+k] = pset->coords[3*(i-1)+k] + delta_vect[k];

  }

  /* Set the ending probe */
  pset->s_coords[n_probes-1] = distance;
  for (int k = 0; k < 3; k++)
    pset->coords[3*(n_probes-1)+k] = end_coords[k];

  return  pset;
}

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
                               const int        *writer_ids)
{
  if (pset == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pset));

  int  n_init_writers = pset->n_writers;
  pset->n_writers += n_writers;
  BFT_REALLOC(pset->writer_ids, pset->n_writers, int);

  for (int i = n_init_writers, j = 0; i < pset->n_writers; i++, j++)
    pset->writer_ids[i] = writer_ids[j];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to the management of a set of probes
 *
 * Available option key names are the following:
 *
 * \li \c \b activated where keyval is either \c true or \c false (default)
 * \li \c \b boundary  where keyval is either \c true or \c false (default)
 * \li \c \b mode      where keyval is \c exact, \c nearest_vertex or
 *                   \c nearest_center (default)
 * \li \c \b profile   where keyval is either \c true or \c false
 * \li \c \b selection_criteria where keyval is selection criteria string
 * \li \c \b tolerance  where keyval is for instance "0.05" (default "0.10")
 * \li \c \b moving_probes  where keyval is either \c true or \c false (default)
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
              _(" Invalid key for setting the set of probes %s.\n"
                " Please read listing for more details and"
                " modify your settings."), pset->name);

  } /* Error message */

  switch(key) {

  case PSETKEY_ACTIVATED:
    if (strcmp(keyval, "true") == 0)
      pset->flag |= CS_PROBE_ACTIVATED;
    else if (strcmp(keyval, "false") == 0) { // remove the flag if it is set
      if (pset->flag & CS_PROBE_ACTIVATED)
        pset->flag ^= CS_PROBE_ACTIVATED;
    }
    else
      bft_error(__FILE__, __LINE__, 0, _err_truefalse_key, keyval, keyname);
    break;

  case PSETKEY_BOUNDARY:
    if (strcmp(keyval, "true") == 0)
      pset->flag |= CS_PROBE_BOUNDARY;
    else if (strcmp(keyval, "false") == 0) { // remove the flag if it is set
      if (pset->flag & CS_PROBE_BOUNDARY)
        pset->flag ^= CS_PROBE_BOUNDARY;
    }
    else
      bft_error(__FILE__, __LINE__, 0, _err_truefalse_key, keyval, keyname);
    break;

  case PSETKEY_MODE:
    if (strcmp(keyval, "exact") == 0)
      pset->mode = CS_PROBE_MODE_EXACT;
    else if (strcmp(keyval, "nearest_cell_center") == 0)
      pset->mode = CS_PROBE_MODE_NEAREST_CELL_CENTER;
    else if (strcmp(keyval, "nearest_vertex") == 0)
      pset->mode = CS_PROBE_MODE_NEAREST_VERTEX;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value %s for setting key %s.\n"
                  " Valid choices are exact, nearest_center or"
                  " nearest_vertex\n"
                  " Please check your settings."), keyval, keyname);
    break;

  case PSETKEY_PROFILE:
    if (strcmp(keyval, "true") == 0)
      pset->flag |= CS_PROBE_PROFILE;
    else if (strcmp(keyval, "false") == 0) { // remove the flag if it is set
      if (pset->flag & CS_PROBE_PROFILE)
        pset->flag ^= CS_PROBE_PROFILE;
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
      pset->flag |= CS_PROBE_TRANSIENT;
    else if (strcmp(keyval, "false") == 0) { // remove the flag if it is set
      if (pset->flag & CS_PROBE_TRANSIENT)
        pset->flag ^= CS_PROBE_TRANSIENT;
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
 * \brief  Try to locate each probe and define the coordinate really use for
 *         the postprocessing step
 *
 * \param[in, out]  pset    pointer to a cs_probe_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_set_locate(cs_probe_set_t   *pset)
{
  if (pset == NULL)
    return;

  const double  tolerance_base = 0.;
  const cs_mesh_t  *mesh = cs_glob_mesh;

  cs_lnum_t  n_select_elements = 0;
  cs_lnum_t  *selected_elements = NULL;

  const bool  on_boundary = (pset->flag & CS_PROBE_BOUNDARY) ? true : false;

  /* Temporary name for the location mesh */
  char  *tmp_name = NULL;

  BFT_MALLOC(tmp_name, strlen(pset->name) + strlen("_tmp") + 1, char);
  sprintf(tmp_name, "%s_tmp", pset->name);

  if (on_boundary) { /* Deal with the surfacic mesh related to the boundary */

    n_select_elements = mesh->n_b_faces;
    if (pset->sel_criter != NULL) {
      if (strcmp(pset->sel_criter, "all[]")) {
        BFT_MALLOC(selected_elements, mesh->n_b_faces, cs_lnum_t);
        cs_selector_get_b_face_num_list(pset->sel_criter,
                                        &n_select_elements, selected_elements);
      }
    } /* Need to define a list of faces ? */

    pset->location_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                         tmp_name,
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

    pset->location_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                                         tmp_name,
                                                         false, // no family info
                                                         n_select_elements,
                                                         selected_elements);

  } /* volumic or surfacic mesh */

  /* Locate probes on this location mesh */
  BFT_MALLOC(pset->entity_num, pset->n_probes, cs_lnum_t);
  BFT_MALLOC(pset->distances, pset->n_probes, float);

  for (int i = 0; i < pset->n_probes; i++) {
    pset->entity_num[i] = -1;
    pset->distances[i] = -1.0;
  }

  int  locate_on_parents = 1; // true by default
  if (pset->mode == CS_PROBE_MODE_NEAREST_VERTEX)
    /* One needs a second step which is easier if parent_num is not used */
    locate_on_parents = 0;

  fvm_point_location_nodal(pset->location_mesh,
                           tolerance_base,
                           pset->tolerance,
                           locate_on_parents,
                           pset->n_probes,
                           NULL,             // point_tag (not useful here)
                           pset->coords,
                           pset->entity_num,
                           pset->distances);

  for (int i = 0; i < pset->n_probes; i++)
    if (pset->entity_num[i] < 0) // Not found
      pset->distances[i] = HUGE_VAL;

  /* Free memory */
  BFT_FREE(tmp_name);
  if (selected_elements != NULL)
    BFT_FREE(selected_elements);

  /* Warning if points have not beeen located */
  cs_gnum_t  n_unlocated_probes = 0;

  if (cs_glob_n_ranks == 1) {
    for (int i = 0; i < pset->n_probes; i++)
      if (pset->distances[i] > 0.99*HUGE_VAL)
        n_unlocated_probes++;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_double_int_t  *gmin_loc = NULL, *loc = NULL;

    BFT_MALLOC(gmin_loc, pset->n_probes, cs_double_int_t);
    BFT_MALLOC(loc, pset->n_probes, cs_double_int_t);

    for (int i = 0; i < pset->n_probes; i++) {
      gmin_loc[i].id = loc[i].id = cs_glob_rank_id;
      gmin_loc[i].val = loc[i].val = pset->distances[i];
    }

    MPI_Allreduce(loc, gmin_loc, pset->n_probes, MPI_DOUBLE_INT, MPI_MINLOC,
                  cs_glob_mpi_comm);

    for (int i = 0; i < pset->n_probes; i++) {
      if (gmin_loc[i].id != cs_glob_rank_id)
        pset->entity_num[i] = -1;
      if (gmin_loc[i].val > 0.99*HUGE_VAL)
        n_unlocated_probes++;
    }

    BFT_FREE(gmin_loc);
    BFT_FREE(loc);
  }
#endif

  if (n_unlocated_probes > 0) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("\nWarning: probe set \"%s\"\n"
                 "         %lu probes have not been located"
                 "         on the associated mesh."),
               pset->name, (unsigned long)n_unlocated_probes);
  }
}

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
                         const char       *mesh_name)
{
  if (pset == NULL)
    return  NULL;

  int  n_exp_probes = 0;
  cs_gnum_t  *global_num = NULL;
  cs_coord_3_t  *probe_coords = NULL;
  fvm_nodal_t  *exp_mesh = fvm_nodal_create(mesh_name, 3);

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *quant = cs_glob_mesh_quantities;

  BFT_MALLOC(probe_coords, pset->n_probes, cs_coord_3_t);
  BFT_MALLOC(global_num, pset->n_probes, cs_gnum_t);

  /* Build the final list of probe coordinates */
  switch (pset->mode) {

  case CS_PROBE_MODE_NEAREST_CELL_CENTER:
    {
      bool  *cell_tag = NULL;

      BFT_MALLOC(cell_tag, mesh->n_cells, bool);
      for (int i = 0; i < mesh->n_cells; i++)
        cell_tag[i] = false;

      for (int i = 0; i < pset->n_probes; i++) {
        if (pset->entity_num[i] > -1) { /* Located on this rank */

          const cs_lnum_t  cell_id = pset->entity_num[i] - 1;

          if (cell_tag[cell_id] == false) {
            cell_tag[cell_id] = true;

            const cs_real_t  *xp = pset->coords + 3*i;
            const cs_real_t  *xc = quant->cell_cen + 3*cell_id;
            for (int k = 0; k < 3; k++)
              probe_coords[n_exp_probes][k] = xc[k];
            pset->distances[i] = cs_math_3_length(xp, xc);
            global_num[n_exp_probes] = i+1;
            n_exp_probes++;
          }

        }
      } // Loop on probes

      BFT_FREE(cell_tag);

    } /* Nearest cell center */
    break;

  case CS_PROBE_MODE_NEAREST_VERTEX:
    {
      fvm_point_location_closest_vertex(pset->location_mesh,
                                        1, // locate on parents
                                        pset->n_probes,
                                        pset->coords,
                                        pset->entity_num,
                                        pset->distances);

      bool  *vtx_tag = NULL;

      BFT_MALLOC(vtx_tag, mesh->n_vertices, bool);
      for (int i = 0; i < mesh->n_vertices; i++)
        vtx_tag[i] = false;

      for (int i = 0; i < pset->n_probes; i++) {
        if (pset->entity_num[i] > -1) { /* Located on this rank */

          const cs_lnum_t  vtx_id = pset->entity_num[i] - 1;

          if (vtx_tag[vtx_id] == false) {
            vtx_tag[vtx_id] = true;
            for (int k = 0; k < 3; k++)
              probe_coords[n_exp_probes][k] = mesh->vtx_coord[3*vtx_id+k];
            global_num[n_exp_probes] = i+1;
            n_exp_probes++;
          }

        }
      } // Loop on probes

      BFT_FREE(vtx_tag);

    } /* Nearest vertex */
    break;

  case CS_PROBE_MODE_EXACT:

    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->entity_num[i] > -1) { /* Located on this rank */
        for (int k = 0; k < 3; k++)
          probe_coords[n_exp_probes][k] = pset->coords[3*i+k];
        global_num[n_exp_probes] = i+1;
        n_exp_probes++;
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              N_(" This mode is not yet implemented to handle probe set.\n"
                 " Please modify your settings."));

  } // Switch on pset->mode

  /* Update the probe set structure */
  pset->n_loc_probes = n_exp_probes;
  BFT_REALLOC(probe_coords, n_exp_probes, cs_coord_3_t);

  fvm_nodal_define_vertex_list(exp_mesh, n_exp_probes, NULL);
  fvm_nodal_transfer_vertices(exp_mesh, (cs_coord_t *)probe_coords);

  cs_real_t  max_distance = 0., gmax_distance = 0.;
  for (int i = 0; i < pset->n_probes; i++)
    if (pset->entity_num[i] > -1)  /* Located on this rank */
      max_distance = fmax(max_distance, pset->distances[i]);

  if (cs_glob_n_ranks > 1) {

    /* Set a global numbering */
    fvm_nodal_init_io_num(exp_mesh, global_num, 0); // 0 = vertices

#if defined(HAVE_MPI)
    MPI_Reduce(&max_distance, &gmax_distance, 1, MPI_DOUBLE, MPI_MAX, 0,
               cs_glob_mpi_comm);
#endif
  }
  else
    gmax_distance = max_distance;

  bft_printf("\n Probe set: \"%s\":\n"
             "   maximum distance between real and requested coordinates:"
             " %5.3e\n", pset->name, gmax_distance);

  BFT_FREE(global_num);

  /* Add global labels */

  if (pset->labels != NULL) {
    char **g_labels;
    int ngl = fvm_nodal_get_n_g_vertices(exp_mesh);
    BFT_MALLOC(g_labels, ngl, char *);

    int j = 0;
    for (int i = 0; i < pset->n_probes; i++) {
      if (pset->distances[i] <= 0.99*HUGE_VAL)
        g_labels[j++] = _copy_label(pset->labels[i]);
    }
    assert(j == ngl);
    fvm_nodal_transfer_global_vertex_labels(exp_mesh, g_labels);

  }

  // probe_coords is managed by exp_mesh

  return exp_mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to a set of probes
 */
/*----------------------------------------------------------------------------*/

void
cs_probe_finalize(void)
{
  for (int i = 0; i < _n_probe_sets; i++)
    _free_probe_set(_probe_set_array + i);

  _n_probe_sets = 0;
  BFT_FREE(_probe_set_array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_probe_set_t structure
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

  bft_printf(" name:       %s\n"
             " flag:       %d\n"
             " mode:       %s\n"
             " tolerance:  %5.3e\n",
             pset->name, pset->flag, cs_probe_mode_name[pset->mode],
             pset->tolerance);

  if (pset->sel_criter != NULL)
    bft_printf(" selection:  %s\n", pset->sel_criter);

  bft_printf(" n_probes:   %d; %d; %d (locally located; defined; max.)\n",
             pset->n_loc_probes, pset->n_probes, pset->n_max_probes);
  bft_printf(" nodal mesh: %p\n\n", (void *)pset->location_mesh);

  for (int i = 0; i < pset->n_probes; i++) {

    bft_printf(" %4d | %-5.3e %-5.3e %-5.3e |", i,
               pset->coords[3*i], pset->coords[3*i+1], pset->coords[3*i+2]);

    if (pset->s_coords != NULL)
      bft_printf(" %5.3e |", pset->s_coords[i]);
    if (pset->entity_num != NULL && pset->distances != NULL)
      bft_printf(" %6d | %5.3e |", pset->entity_num[i], pset->distances[i]);
    if (pset->labels != NULL)
      if (pset->labels[i] != NULL)
        bft_printf(" %s", pset->labels[i]);
    bft_printf("\n");

  } // Loop on probes

  if (true && pset->location_mesh != NULL)
    fvm_nodal_dump(pset->location_mesh);
}

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
                         float                  *distances[])
{
  if (pset == NULL)
    return;

  switch (pset->mode) {

  case CS_PROBE_MODE_NEAREST_CELL_CENTER:
    *mode = 0;
    break;
  case CS_PROBE_MODE_NEAREST_VERTEX:
    *mode = 1;
    break;
  case CS_PROBE_MODE_EXACT:
    *mode = 2;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              N_(" This mode is not yet implemented to handle probe set.\n"
                 " Please modify your settings."));

  } // Switch on pset->mode

  /* Return pointers */
  *n_probes = pset->n_probes;
  *coords = pset->coords;
  *ent_num = pset->entity_num;
  *distances = pset->distances;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
