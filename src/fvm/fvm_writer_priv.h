#ifndef __FVM_WRITER_PRIV_H__
#define __FVM_WRITER_PRIV_H__

/*============================================================================
 * Private types for mesh and field writers
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_timer.h"

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Writer format implementation and functionality info
 */

#define FVM_WRITER_FORMAT_USE_EXTERNAL      (1 << 0)

#define FVM_WRITER_FORMAT_HAS_POLYGON       (1 << 1)
#define FVM_WRITER_FORMAT_HAS_POLYHEDRON    (1 << 2)

#define FVM_WRITER_FORMAT_SEPARATE_MESHES   (1 << 3)

#define FVM_WRITER_FORMAT_NAME_IS_OPTIONAL  (1 << 4)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int
(fvm_writer_n_version_strings_t) (void);

typedef const char *
(fvm_writer_version_string_t)(int string_index,
                              int compile_time_version);

#if defined(HAVE_MPI)

typedef void *
(fvm_writer_init_t) (const char             *name,
                     const char             *path,
                     const char             *options,
                     fvm_writer_time_dep_t   time_dependency,
                     MPI_Comm                comm);

#else

typedef void *
(fvm_writer_init_t) (const char             *name,
                     const char             *path,
                     const char             *options,
                     fvm_writer_time_dep_t   time_dependency);

#endif /* defined(HAVE_MPI) */

typedef void *
(fvm_writer_finalize_t) (void  *this_writer);

typedef void
(fvm_writer_set_mesh_time_t) (void    *this_writer,
                              int      time_step,
                              double   time_value);

typedef int
(fvm_writer_needs_tesselation_t) (void               *this_writer,
                                  const fvm_nodal_t  *mesh,
                                  fvm_element_t       element_type);

typedef void
(fvm_writer_export_nodal_t) (void               *this_writer,
                             const fvm_nodal_t  *mesh);

typedef void
(fvm_writer_export_field_t) (void                   *this_writer,
                             const fvm_nodal_t      *mesh,
                             const char             *name,
                             fvm_writer_var_loc_t    location,
                             int                     dimension,
                             cs_interlace_t          interlace,
                             int                     n_parent_lists,
                             const cs_lnum_t         parent_num_shift[],
                             cs_datatype_t           datatype,
                             int                     time_step,
                             double                  time_value,
                             const void       *const field_values[]);

typedef void
(fvm_writer_flush_t) (void  *this_writer);

/*----------------------------------------------------------------------------
 * Format information structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                     name[32];     /* Format name */
  char                     version[16];  /* Format version (if defined) */
  int                      info_mask;    /* Additional format info */
  fvm_writer_time_dep_t    max_time_dep; /* Maximum time dependency level
                                            possible with this format */

  int                      dl_count;     /* Number of writers using the
                                            dynamically loadable library
                                            for this format, if relevant */
  void                    *dl_lib;       /* Pointer to dynamically loadable
                                            library, if used */
  const char              *dl_name;      /* Prefix for name of dynamically
                                            loadable library, or NULL */
  const char              *dl_prefix;    /* Prefix for exported symbols of
                                            dynamically loadable library,
                                            or NULL */

  fvm_writer_n_version_strings_t  *n_version_strings_func;
  fvm_writer_version_string_t     *version_string_func;
  fvm_writer_init_t               *init_func;
  fvm_writer_finalize_t           *finalize_func;
  fvm_writer_set_mesh_time_t      *set_mesh_time_func;
  fvm_writer_needs_tesselation_t  *needs_tesselation_func;
  fvm_writer_export_nodal_t       *export_nodal_func;
  fvm_writer_export_field_t       *export_field_func;
  fvm_writer_flush_t              *flush_func;

} fvm_writer_format_t;

/*----------------------------------------------------------------------------
 * Structure defining a writer definition
 *----------------------------------------------------------------------------*/

struct _fvm_writer_t {

  char                   *name;              /* Writer name */
  fvm_writer_format_t    *format;            /* Output format */
  char                   *options;           /* Output options */
  char                   *path;              /* Output path */
  fvm_writer_time_dep_t   time_dep;          /* Geometry time dependency */
  int                     n_format_writers;  /* number of actual writers */
  void                  **format_writer;     /* Format-specific writers */
  char                  **mesh_names;        /* List of mesh names if one
                                                writer per mesh is required */

  cs_timer_counter_t      mesh_time;         /* Meshes output timer */
  cs_timer_counter_t      field_time;        /* Fields output timer */
  cs_timer_counter_t      flush_time;        /* output "completion" timer */

};

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_WRITER_PRIV_H__ */
