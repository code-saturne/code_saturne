/*============================================================================
 * Mesh element flag checking and error handling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_flag_check.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_flag_check.c
        Mesh element flag checking and error handling..
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Mappings to MPI datatypes */

#if defined(HAVE_MPI)

typedef struct
{
  int val;
  int rank;
} _mpi_int_int_t;

#endif /* defined(HAVE_MPI) */

/* Face marker structure for selecting error postprocessing output faces. */

typedef struct {

  cs_lnum_t   n_elts;    /* Number of elements faces */
  int         min_flag;  /* Minimum valid flag value */
  const int  *flag;      /* >= min_flag for valid elements,
                            < min_flag for marked faces */

} _error_elt_marker_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transfer info on element with lowest global number to rank 0.
 *
 * \param[in]  elt_gnum  pointer to global element number
 * \param[in]  elt_type  pointer to element type, or NULL
 * \param[in]  elt_coo   pointer toelement coordinates
 */
/*----------------------------------------------------------------------------*/

static void
_min_gnum_elt(cs_gnum_t  *elt_gnum,
              int        *elt_type,
              double      elt_coo[3])
{
#if defined(HAVE_MPI)

  /* local variables */

  cs_gnum_t  min_elt_gnum;
  _mpi_int_int_t  val_in, val_min;

  /* Return immediately if not running under MPI */

  if (cs_glob_n_ranks < 2)
    return;

  int _elt_type = 0;
  if (elt_type != NULL)
    _elt_type = *elt_type;

  /* Obtain the lowest global elt number with an error; use minloc
     with a marker, rather than with a global number directly, in
     case cs_gnum_t is larger than an integer) */

  MPI_Allreduce(elt_gnum, &min_elt_gnum, 1, CS_MPI_GNUM, MPI_MIN,
                cs_glob_mpi_comm);

  if (*elt_gnum == min_elt_gnum)
    val_in.val = 0;
  else
    val_in.val = 1;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_min, 1, MPI_2INT, MPI_MINLOC,
                cs_glob_mpi_comm);

  /* Now exchange values */

  if (val_min.rank > 0) {

    if (val_min.rank == cs_glob_rank_id) {
      MPI_Send(elt_gnum, 1, CS_MPI_GNUM, 0, 1, cs_glob_mpi_comm);
      MPI_Send(&_elt_type, 1, MPI_INT, 0, 2, cs_glob_mpi_comm);
      MPI_Send(elt_coo, 3, MPI_DOUBLE, 0, 3, cs_glob_mpi_comm);
    }
    else if (cs_glob_rank_id == 0) {
      MPI_Status status;
      MPI_Recv(elt_gnum, 1, CS_MPI_GNUM, val_min.rank, 1,
               cs_glob_mpi_comm, &status);
      MPI_Recv(&_elt_type, 1, MPI_INT, val_min.rank, 2,
               cs_glob_mpi_comm, &status);
      MPI_Recv(elt_coo, 3, MPI_DOUBLE, val_min.rank, 3,
               cs_glob_mpi_comm, &status);
    }

  }

  if (elt_type != NULL)
    *elt_type = _elt_type;

#endif
}

/*----------------------------------------------------------------------------
 * Function for selection of elements with flag errors.
 *
 * parameters:
 *   input    <-- pointer to input (elt_marker structure)
 *   n_eltd   --> number of selected elements
 *   elt_ids  --> array of selected element ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_post_error_elt_select(void         *input,
                       cs_lnum_t    *n_elts,
                       cs_lnum_t   **elt_ids)
{
  cs_lnum_t elt_id;

  cs_lnum_t _n_elts = 0;
  cs_lnum_t *_elt_ids = NULL;

  const _error_elt_marker_t  *marker = input;

  BFT_MALLOC(_elt_ids, marker->n_elts, cs_lnum_t);

  for (elt_id = 0; elt_id < marker->n_elts; elt_id++) {
    if (marker->flag[elt_id] < marker->min_flag)
      _elt_ids[_n_elts++] = elt_id;
  }

  *n_elts = _n_elts;
  *elt_ids = _elt_ids;
}

/*----------------------------------------------------------------------------
 * Function for selection of elements with valid flags.
 *
 * parameters:
 *   input   <-- pointer to input (elt_marker structure)
 *   n_elts  --> number of selected elts
 *   elt_ids --> array of selected elt ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_post_valid_elt_select(void         *input,
                       cs_lnum_t    *n_elts,
                       cs_lnum_t   **elt_ids)
{
  cs_lnum_t elt_id;

  cs_lnum_t _n_elts = 0;
  cs_lnum_t *_elt_ids = NULL;

  const _error_elt_marker_t  *marker = input;

  BFT_MALLOC(_elt_ids, marker->n_elts, cs_lnum_t);

  for (elt_id = 0; elt_id < marker->n_elts; elt_id++) {
    if (marker->flag[elt_id] >= marker->min_flag)
      _elt_ids[_n_elts++] = elt_id;
  }

  *n_elts = _n_elts;
  *elt_ids = _elt_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle boundary face postprocessing output based on flags.
 *
 * It is assumed that element flags are usually positive integers, and that
 * in case of a detected error, their signs have been set to a negative value.
 *
 * A minimum allowed value may be specified, so for example 0 may be
 * considered a valid or invalid flag depending on that minimum.
 *
 * \param[in]  error_mesh_name  mesh name for elements with error
 * \param[in]  valid_mesh_name  mesh name for elements without
 * \param[in]  flag_label       field name for flag
 * \param[in]  n_elts           number of elements
 * \param[in]  min_flag         minimum allowed flag
 * \param[in]  flag             current element flag
 */
/*----------------------------------------------------------------------------*/

static void
_postprocess(const char   *error_mesh_name,
             const char   *valid_mesh_name,
             const char   *flag_label,
             int           location_id,
             int           min_flag,
             const int     flag[])
{
  cs_lnum_t n_elts = 0;
  const cs_mesh_t *m = cs_glob_mesh;

  switch(location_id) {
  case CS_MESH_LOCATION_CELLS:
    n_elts = m->n_cells;
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    n_elts = m->n_b_faces;
    break;
  default:
    assert(0);
    return;
  }

  _error_elt_marker_t  marker = {.n_elts = n_elts,
                                 .min_flag = min_flag,
                                 .flag = flag};

  cs_gnum_t n_g_valid_elts = 0;
  int mesh_id[2] = {0, 0};

  const int writer_id = -2;
  const int writer_ids[] = {writer_id};

  cs_post_init_error_writer();

  /* Mesh for invalid elements */

  mesh_id[0] = cs_post_get_free_mesh_id();

  switch(location_id) {
  case CS_MESH_LOCATION_CELLS:
    cs_post_define_volume_mesh_by_func(mesh_id[0],
                                       error_mesh_name,
                                       _post_error_elt_select,
                                       &marker,
                                       false, /* time varying */
                                       true,  /* add groups if present */
                                       false, /* auto variables */
                                       1,
                                       writer_ids);
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    cs_post_define_surface_mesh_by_func(mesh_id[0],
                                        error_mesh_name,
                                        NULL,
                                        _post_error_elt_select,
                                        NULL,
                                        &marker,
                                        false, /* time varying */
                                        true,  /* add groups if present */
                                        false, /* auto variables */
                                        1,
                                        writer_ids);
    break;
  default:
    break;
  }

  /* Mesh for valid faces */

  for (cs_lnum_t elt_id = 0; elt_id < n_elts; elt_id++) {
    if (flag[elt_id] >= min_flag)
      n_g_valid_elts += 1;
  }

  cs_parall_counter(&n_g_valid_elts, 1);

  if (n_g_valid_elts > 0) {

    mesh_id[1] = cs_post_get_free_mesh_id();

    switch(location_id) {
    case CS_MESH_LOCATION_CELLS:
      cs_post_define_volume_mesh_by_func(mesh_id[1],
                                         valid_mesh_name,
                                         _post_valid_elt_select,
                                         &marker,
                                         false, /* time varying */
                                         true,  /* add groups if present */
                                         false, /* auto variables */
                                         1,
                                         writer_ids);
      break;
    case CS_MESH_LOCATION_BOUNDARY_FACES:
      cs_post_define_surface_mesh_by_func(mesh_id[1],
                                          valid_mesh_name,
                                          NULL,
                                          _post_valid_elt_select,
                                          NULL,
                                          &marker,
                                          false, /* time varying */
                                          true,  /* add groups if present */
                                          false, /* auto variables */
                                          1,
                                          writer_ids);
      break;
    default:
      break;
    }

  }

  cs_post_activate_writer(writer_id, 1);

  cs_post_write_meshes(NULL);

  {
    char var_name[32];
    strncpy(var_name, flag_label, 31);
    var_name[31] = '\0';

    int *_flag;
    BFT_MALLOC(_flag, n_elts, int);
    for (cs_lnum_t i = 0; i < n_elts; i++)
      _flag[i] = abs(flag[i]);

    for (int ii = 0; ii < 2; ii++) {
      if (mesh_id[ii] != 0)
        cs_post_write_var(mesh_id[ii],
                          writer_id,
                          var_name,
                          1,
                          false, /* no interlace */
                          true,  /* use parents */
                          CS_POST_TYPE_cs_int_t,
                          NULL,
                          NULL,
                          _flag,
                          NULL);

    }

    BFT_FREE(_flag);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for and handle errors with an associated element flag
 *
 * It is assumed that element flags are usually positive integers, and that
 * in case of a detected error, their signs have been set to a negative value.
 *
 * A minimum allowed value may be specified, so for example 0 may be
 * considered a valid or invalid flag depending on that minimum.
 *
 * This function exits silently if no such marked elements are present in the
 * computational domain.
 *
 * Otherwise, it logs information on the first detected error location, and
 * outputs postprocessing visualization information to assist debugging.
 *
 * If the error status (i.e. negative flag) is known locally but not
 * globally, use \ref cs_flag_check.
 *
 * Currently supported locations are CS_MESH_LOCATION_CELLS and
 * CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * \param[in]  err_elt_descr    description fro first element with error
 * \param[in]  flag_descr       flag type description
 * \param[in]  flag_label       field label for flag postprocessing
 * \param[in]  error_mesh_name  postprocessing mesh name for elements with error
 * \param[in]  valid_mesh_name  postprocessing mesh name for valid elements
 * \param[in]  location_id      associated mesh location
 * \param[in]  min_flag         minimum allowed flag
 * \param[in]  elt_flag         current element flag
 */
/*----------------------------------------------------------------------------*/

int
cs_flag_check(const char   *err_elt_descr,
              const char   *flag_descr,
              const char   *flag_label,
              const char   *error_mesh_name,
              const char   *valid_mesh_name,
              int           location_id,
              int           min_flag,
              const int     elt_flag[])
{
  cs_lnum_t n_elts = cs_mesh_location_get_n_elts(location_id)[0];

  /* Check for error */

  int error_flag = 0;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    if (elt_flag[i] < min_flag) {
      error_flag = 1;
      break;
    }
  }

  cs_parall_max(1, CS_INT_TYPE, &error_flag);

  /* Handle error */

  if (error_flag)
    cs_flag_check_error_info(err_elt_descr,
                             flag_descr,
                             flag_label,
                             error_mesh_name,
                             valid_mesh_name,
                             location_id,
                             min_flag,
                             elt_flag);

  return error_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle an error with an associated element flag
 *
 * This function logs information on the first detected error location, and
 * outputs postprocessing visualization information to assist debugging.
 *
 * It is assumed that element flags are usually positive integers, and that
 * in case of a detected error, their signs have been set to a negative value.
 *
 * A minimum allowed value may be specified, so for example 0 may be
 * considered a valid or invalid flag depending on that minimum.
 *
 * This function should be called when the error status has been previously
 * checked, and all ranks know that an error is present.
 *
 * If the error status (i.e. negative flag) is known locally but not
 * globally, use \ref cs_flag_check.
 *
 * Currently supported locations are CS_MESH_LOCATION_CELLS and
 * CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * \param[in]  err_elt_descr    description fro first element with error
 * \param[in]  flag_descr       flag type description
 * \param[in]  flag_label       field label for flag postprocessing
 * \param[in]  error_mesh_name  postprocessing mesh name for elements with error
 * \param[in]  valid_mesh_name  postprocessing mesh name for valid elements
 * \param[in]  location_id      associated mesh location
 * \param[in]  min_flag         minimum allowed flag
 * \param[in]  elt_flag         current element flag
 */
/*----------------------------------------------------------------------------*/

void
cs_flag_check_error_info(const char   *err_elt_descr,
                         const char   *flag_descr,
                         const char   *flag_label,
                         const char   *error_mesh_name,
                         const char   *valid_mesh_name,
                         int           location_id,
                         int           min_flag,
                         const int     elt_flag[])
{
  cs_gnum_t   n_g_errors = 0;

  cs_lnum_t         n_elts = 0;
  const cs_gnum_t  *g_elt_num = NULL;
  const cs_real_t  *elt_coo = NULL;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  switch(location_id) {
  case CS_MESH_LOCATION_CELLS:
    n_elts = m->n_cells;
    g_elt_num = m->global_cell_num;
    elt_coo = mq->cell_cen;
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    n_elts = m->n_b_faces;
    g_elt_num = m->global_b_face_num;
    elt_coo = mq->b_face_cog;
    break;
  default:
    assert(0);
    return;
  }

  /* Count and mark cells with problems */

  cs_real_t  err_elt_coo[3] = {0, 0, 0};
  int        err_flag = 0;
  cs_gnum_t  err_elt_gnum = 0;

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    if (elt_flag[i] < min_flag) {

      cs_gnum_t elt_gnum;

      if (g_elt_num != NULL)
        elt_gnum = g_elt_num[i];
      else
        elt_gnum = i + 1;

      if (err_elt_gnum == 0 || elt_gnum < err_elt_gnum) {
        err_flag = elt_flag[i];
        for (cs_lnum_t j = 0; j < 3; j++)
          err_elt_coo[j] = elt_coo[i*3 + j];
      }

      n_g_errors += 1;

    }
  }

  /* Obtain the lowest global elt number with an error,
     and print associated info */

  _min_gnum_elt(&err_elt_gnum, &err_flag, err_elt_coo);
  cs_parall_counter(&n_g_errors, 1);

  if (cs_glob_rank_id < 1)
    bft_printf(_("\nFirst %s\n"
                 "  (out of %llu)\n"
                 "  has %s %d, center (%g, %g, %g)\n\n"),
               err_elt_descr,
               (unsigned long long)n_g_errors,
               flag_descr, abs(err_flag),
               err_elt_coo[0], err_elt_coo[1], err_elt_coo[2]);

  /* Activate postprocessing */

  _postprocess(error_mesh_name,
               valid_mesh_name,
               flag_label,
               location_id,
               min_flag,
               elt_flag);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
