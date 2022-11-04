/*============================================================================
 * Function objects management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_function.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_function.c

  \brief Function objects management.

  Function objects can have various roles. Their main use is to unify
  handling of expression evaluations for mesh-location based data that
  can be re-evaluated on the fly (rather than requiring persistent storage
  such as field data), mostly for logging and post-processing.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* Function descriptor allocation block size */

#define _CS_FUNCTION_S_ALLOC_SIZE       16

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static int  _n_functions = 0;
static int  _n_functions_max = 0;

static cs_function_t  **_functions = NULL;
static cs_map_name_to_id_t  *_function_map = NULL;

/* Names for logging */

static const int _n_type_flags = 3;
static const int _type_flag_mask[] = {CS_FUNCTION_INTENSIVE,
                                      CS_FUNCTION_EXTENSIVE,
                                      CS_FUNCTION_USER};
static const char *_type_flag_name[] = {N_("intensive"),
                                        N_("extensive"),
                                        N_("user")};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a function descriptor.
 *
 * parameters:
 *   name         <-- function name
 *   is_intensive <-- are the function values intensive or not (extensive) ?
 *   location_id  <-- id of associated location
 *   dim          <-- function dimension (number of components)
 *   datatype     <-- associated data type
 *
 * returns:
 *   pointer to new function.
 *----------------------------------------------------------------------------*/

static cs_function_t *
_function_create(const char    *name,
                 bool           is_intensive,
                 int            location_id,
                 int            dim,
                 cs_datatype_t  datatype)
{
  int function_id = -1;
  size_t l = strlen(name);

  cs_function_t *f = cs_function_by_name_try(name);

  /* Check this name was not already used */

  if (f != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating function:\n"
                "  name:        \"%s\"\n"
                "  location_id: %d\n"
                "  dimension:   %d\n\n"
                "A function with that name has already been defined:\n"
                "  id:          %d\n"
                "  location_id: %d\n"
                "  dimension:   %d"),
              name, location_id, dim, f->id, f->location_id, f->dim);

  /* Initialize if necessary */

  if (_function_map == NULL)
    _function_map = cs_map_name_to_id_create();

  if (l == 0)
    bft_error(__FILE__, __LINE__, 0, _("Defining a function requires a name."));

  for (size_t i = 0; i < l; i++) {
    if (name[i] == '[' || name[i] == ']')
      bft_error(__FILE__, __LINE__, 0,
                _("Function \"%s\" is not allowed,\n"
                  "as \'[\' and \']\' are reserved for component access."),
                name);
  }

  /* Insert entry in map */

  function_id = cs_map_name_to_id(_function_map, name);

  if (function_id == _n_functions)
    _n_functions = function_id + 1;

  /* Reallocate functions pointer if necessary */

  if (_n_functions > _n_functions_max) {
    if (_n_functions_max == 0)
      _n_functions_max = 8;
    else
      _n_functions_max *= 2;
    BFT_REALLOC(_functions, _n_functions_max, cs_function_t *);
  }

  /* Allocate functions descriptor block if necessary
     (to reduce fragmentation and improve locality of function
     descriptors, they are allocated in blocks) */

  int shift_in_alloc_block = function_id % _CS_FUNCTION_S_ALLOC_SIZE;
  if (shift_in_alloc_block == 0)
    BFT_MALLOC(_functions[function_id], _CS_FUNCTION_S_ALLOC_SIZE,
               cs_function_t);
  else
    _functions[function_id] = _functions[function_id - shift_in_alloc_block]
                                                     + shift_in_alloc_block;

  f = _functions[function_id];

  /* Check type flags and location id */

  if (location_id < 0 || location_id >= cs_mesh_location_n_locations())
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh location %d associated with function \"%s\"\n"
                " has not been defined yet."),
              location_id, name);

  /* Assign function */

  cs_function_t f_ini
    = {.name = cs_map_name_to_id_reverse(_function_map, function_id),
       .label = NULL,
       .id = function_id,
       .type = 0,
       .dim = dim,
       .location_id = location_id,
       .datatype = datatype,
       .post_vis = 0,
       .log = 0,
       .time_stamp = -1,
       .restart_file = CS_RESTART_DISABLED,
       .eval_func = NULL,
       .analytic_func = NULL,
       .dof_func = NULL,
       .func_input = NULL};

  f = _functions[function_id];
  memcpy(f, &f_ini, sizeof(cs_function_t));

  f->type = (is_intensive) ? CS_FUNCTION_INTENSIVE : CS_FUNCTION_EXTENSIVE;
  f->type |= CS_FUNCTION_USER;  /* to be unset explicitely for predefined
                                   functions */

  return f;
}

/*----------------------------------------------------------------------------
 * Add type flag info to the current position in the setup log.
 *
 * parameters:
 *   type <-- type flag
 *----------------------------------------------------------------------------*/

static inline void
_log_add_type_flag(int type)
{
  int n_loc_flags = 0;

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
    cs_log_printf(CS_LOG_SETUP, ")");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided evaluation function.
 *
 * If of dimension > 1, the evaluated values are always interleaved.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  datatype      associated data values type
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_func(const char             *name,
                           int                     location_id,
                           int                     dim,
                           bool                    is_intensive,
                           cs_datatype_t           datatype,
                           cs_eval_at_location_t  *data_func,
                           void                   *data_input)
{
  cs_function_t *f = _function_create(name,
                                      is_intensive,
                                      location_id,
                                      dim,
                                      datatype);

  f->eval_func = data_func;
  f->func_input = data_input;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided "degree of freedom" type evaluation function.
 *
 * The provided function and optional associated input is of the same
 * form as and may be shared with some boundary condition or property
 * definitions.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_analytic_func(const char          *name,
                                    int                  location_id,
                                    int                  dim,
                                    bool                 is_intensive,
                                    cs_analytic_func_t  *data_func,
                                    void                *data_input)
{
  cs_function_t *f = _function_create(name,
                                      is_intensive,
                                      location_id,
                                      dim,
                                      CS_REAL_TYPE);

  f->analytic_func = data_func;
  f->func_input = data_input;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a function whose data values will be computed using the
 *        provided "degree of freedom" type evaluation function.
 *
 * The provided function and optional associated input is of the same
 * form as and may be shared with some boundary condition or property
 * definitions.
 *
 * \param[in]  name          name of associated function
 * \param[in]  location_id   id of associated mesh location
 * \param[in]  dim           dimension associated with element data
 * \param[in]  is_intensive  is the function intensive?
 * \param[in]  data_func     function used to define data values
 * \param[in]  data_input    pointer to optional (untyped) value or structure
 *                           to be used by data_func

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_by_dof_func(const char      *name,
                               int              location_id,
                               int              dim,
                               bool             is_intensive,
                               cs_dof_func_t   *data_func,
                               void            *data_input)
{
  cs_function_t *f = _function_create(name,
                                      is_intensive,
                                      location_id,
                                      dim,
                                      CS_REAL_TYPE);

  f->dof_func = data_func;
  f->func_input = data_input;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all functions management metadata.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_destroy_all(void)
{
  for (int i = 0; i < _n_functions; i++) {
    cs_function_t  *f = _functions[i];
    BFT_FREE(f->label);
  }

  for (int i = 0; i < _n_functions; i++) {
    if (i % _CS_FUNCTION_S_ALLOC_SIZE == 0)
      BFT_FREE(_functions[i]);
  }

  BFT_FREE(_functions);

  cs_map_name_to_id_destroy(&_function_map);

  _n_functions = 0;
  _n_functions_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined functions.
 *
 * \return  number of defined functions
 */
/*----------------------------------------------------------------------------*/

int
cs_function_n_functions(void)
{
  return _n_functions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its id.
 *
 * This function requires that a function of the given id is defined.
 *
 * \param[in]  id   function id
 *
 * \return  pointer to the function structure
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_id(int  id)
{
  if (id > -1 && id < _n_functions)
    return _functions[id];
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Function with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its name.
 *
 * This function requires that a function of the given name is defined.
 *
 * \param[in]  name  function name
 *
 * \return  pointer to the function structure
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_function_map, name);

  if (id > -1)
    return _functions[id];
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Function \"%s\" is not defined."), name);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a function object based on its name if present.
 *
 * If no function of the given name is defined, NULL is returned.
 *
 * \param[in]  name  function name
 *
 * \return  pointer to the function structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_function_t  *
cs_function_by_name_try(const char  *name)
{
  int id = cs_map_name_to_id_try(_function_map, name);

  if (id > -1)
    return _functions[id];
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assig a label to a function object.
 *
 * \param[in, out]  f      pointer to associated function handle
 * \param[in]       label  associated label
 */
/*----------------------------------------------------------------------------*/

void
cs_function_set_label(cs_function_t   *f,
                      const char      *label)
{
  BFT_REALLOC(f->label, strlen(label) + 1, char);
  strcpy(f->label, label);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log function definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_log_defs(void)
{
  if (_n_functions == 0)
    return;

  int cat_id_max = -1;
  int *in_cat_id;
  BFT_MALLOC(in_cat_id, _n_functions, int);
  for (int i = 0; i < _n_functions; i++)
    in_cat_id[i] = -1;

  const char *cat_id_name[] = {N_("mesh-independent"),
                               N_("mesh-independent, user"),
                               N_("mesh-based"),
                               N_("mesh-based, user")};

  /* First loop to determine category id */

  for (int i = 0; i < _n_functions; i++) {
    const cs_function_t *f = _functions[i];

    int cat_id = 0;
    if (f->location_id > 0)
      cat_id += 2;
    if (f->type & CS_FUNCTION_USER)
      cat_id += 1;

    in_cat_id[i] = cat_id;
    if (cat_id > cat_id_max)
      cat_id_max = cat_id;
  }

  /* Functions by category */

  for (int cat_id = 0; cat_id < cat_id_max + 1; cat_id++) {

    size_t name_width = 24;

    /* First loop to determine name width */

    int n_cat_functions = 0;

    for (int i = 0; i < _n_functions; i++) {

      const cs_function_t *f = _functions[i];

      if (in_cat_id[i] != cat_id)
        continue;

      size_t l = strlen(f->name);
      if (l > name_width)
        name_width = l;
    }

    if (name_width > 63)
      name_width = 63;

    /* Main loop */

    for (int i = 0; i < _n_functions; i++) {

      char ilv_c = ' ';

      const cs_function_t *f = _functions[i];

      if (in_cat_id[i] != cat_id)
        continue;

      char tmp_s[6][64] =  {"", "", "", "", "", ""};

      /* Print header for first function of each category */

      if (n_cat_functions == 0) {

        cs_log_strpad(tmp_s[0], _("Function"), name_width, 64);
        cs_log_strpad(tmp_s[1], _("Dim."), 4, 64);
        cs_log_strpad(tmp_s[2], _("Location"), 20, 64);
        cs_log_strpad(tmp_s[3], _("Id"), 4, 64);
        cs_log_strpad(tmp_s[4], _("Type"), 6, 64);

        /* Print logging header */

        cs_log_printf(CS_LOG_SETUP,
                      _("\n"
                        "Functions of type: %s\n"
                        "-----------------\n\n"), _(cat_id_name[cat_id]));

        cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s %s Type flag\n"),
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

        for (int j = 0; j < 5; j++)
          memset(tmp_s[j], '-', 64);

        tmp_s[0][name_width] = '\0';
        tmp_s[1][4] = '\0';
        tmp_s[2][20] = '\0';
        tmp_s[3][4] = '\0';
        tmp_s[4][6] = '\0';

        cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s %s ---------\n"),
                      tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3], tmp_s[4]);

      }

      /* Print function info */

      cs_log_strpad(tmp_s[0], f->name, name_width, 64);

      cs_log_strpad(tmp_s[1],
                    _(cs_mesh_location_get_name(f->location_id)),
                    20,
                    64);

      cs_log_strpad(tmp_s[4],
                    cs_datatype_name[f->datatype],
                    6,
                    64);

      cs_log_printf(CS_LOG_SETUP,
                    "  %s %d %c  %s %-4d %s ",
                    tmp_s[0], f->dim, ilv_c,
                    tmp_s[1], f->id, tmp_s[4]);

      if (f->type != 0) {
        cs_log_printf(CS_LOG_SETUP, "%-4d", f->type);
        _log_add_type_flag(f->type);
        cs_log_printf(CS_LOG_SETUP, "\n");
      }
      else
        cs_log_printf(CS_LOG_SETUP, "0\n");

      n_cat_functions++;

    } /* End of loop on functions */

  } /* End fo loop on categories */

  BFT_FREE(in_cat_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all given function object settings
 *        to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_log_all_settings(void)
{
  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Settings per function:\n"
                  "---------------------\n"));

  /* First loop to determine field width */

  size_t name_width = 24;

  for (int i = 0; i < _n_functions; i++) {
    const cs_function_t *f = _functions[i];
    size_t l = strlen(f->name);
    if (l > name_width)
      name_width = l;
  }
  if (name_width > 63)
    name_width = 63;

  /* Now loop on settings */

  const char null_str[] = "(null)";

  const char *setting[] = {"label", "post_vis", "log", "restart_file"};

  for (int i = 0; i < 4; i++) {

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "  Member: \"%s\", values per function object\n"
                    "  -------\n"),
                  setting[i]);

    for (int j = 0; j < _n_functions; j++) {

      const cs_function_t *f = _functions[j];
      char name_s[64] =  "";
      cs_log_strpad(name_s, f->name, name_width, 64);

      switch (i) {
      case 0:
        {
          const char *s = f->label;
          if (s == NULL)
            s = null_str;
          cs_log_printf(CS_LOG_SETUP, _("    %s %s\n"), name_s, s);
        }
        break;
      case 1:
        cs_log_printf(CS_LOG_SETUP, "    %s %d\n",
                      name_s, f->post_vis);
        break;
      case 2:
        cs_log_printf(CS_LOG_SETUP, "    %s %d\n",
                      name_s, f->log);
        break;
      case 4:
        if (f->restart_file < 0)
          cs_log_printf(CS_LOG_SETUP, _("    %s <disabled>\n"),
                        name_s);
        else
          cs_log_printf(CS_LOG_SETUP, "    %s %d\n",
                        name_s, (int)(f->restart_file));
        break;
      default:
        break;
      }

    } /* End of loop on function objects */

  } /* End of loop on members */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate function values.
 *
 * If the matching values are multidimensional, they must be interleaved.
 * The output values are assumed to use a dense storage (i.e. of size
 * \c n_elts * <value dimension> for the associated data type, in the same
 * order as \c elt_ids if present.)
 *
 * \param[in]       f            pointer to associated function handle
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_function_evaluate(const cs_function_t   *f,
                     const cs_time_step_t  *ts,
                     int                    location_id,
                     cs_lnum_t              n_elts,
                     const cs_lnum_t       *elt_ids,
                     void                  *vals)
{
  if (f->eval_func != NULL)
    f->eval_func(location_id,
                 n_elts,
                 elt_ids,
                 f->func_input,
                 vals);

  else if (f->analytic_func != NULL) {

    const double t_cur = (ts != NULL) ? ts->t_cur : 0.;

    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
    const cs_mesh_location_type_t loc_type
      = cs_mesh_location_get_type(location_id);
    const cs_real_t *base_coords = NULL;
    if (loc_type == CS_MESH_LOCATION_CELLS)
      base_coords = (const cs_real_t *)(mq->cell_cen);
    else if (loc_type == CS_MESH_LOCATION_INTERIOR_FACES)
      base_coords = (const cs_real_t *)(mq->i_face_cog);
    else if (loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
      base_coords = (const cs_real_t *)(mq->b_face_cog);
    else if (loc_type == CS_MESH_LOCATION_VERTICES)
      base_coords = (const cs_real_t *)(cs_glob_mesh->vtx_coord);

    f->analytic_func(t_cur,
                     n_elts,
                     elt_ids,
                     base_coords,
                     true,
                     f->func_input,
                     vals);
  }

  else if (f->dof_func != NULL) {
    f->dof_func(n_elts,
                elt_ids,
                true,
                f->func_input,
                vals);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
