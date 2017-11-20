#ifndef __FVM_TO_HISTOGRAM_H__
#define __FVM_TO_HISTOGRAM_H__

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to histogram (whitepace-separated or csv) files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of histogram file format */

typedef enum {
  CS_HISTOGRAM_TXT,  /* .txt file */
  CS_HISTOGRAM_TEX,  /* .tex file */
  CS_HISTOGRAM_PNG   /* .png file */
} cs_histogram_format_t;

/*----------------------------------------------------------------------------
 * Histogram writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;               /* Writer name */
  char        *path;               /* Path prefix */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  cs_histogram_format_t  format;   /* Histogram format */

  int               nt;            /* Time step */
  double            t;             /* Time value */

  cs_real_t        *buffer;        /* Values buffer */

  char             *file_name;     /* File name */
  FILE             *f;             /* Associated histograms */

  int               n_sub;         /* Number of subdivisions */

#if defined(HAVE_MPI)
  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

} fvm_to_histogram_writer_t;

/*----------------------------------------------------------------------------*/
/*!
 * Function pointer to display histograms.
 *
 * parameters:
 *  var_min  <-- minimum variable value
 *  var_max  <-- maximum variable value
 *  count    <-- count for each histogram slice
 *  f        <-- pointer to the histogram file
 */
/*----------------------------------------------------------------------------*/

typedef void
(fvm_to_histogram_display_t) (cs_real_t                   var_min,
                              cs_real_t                   var_max,
                              cs_gnum_t                   count[],
                              fvm_to_histogram_writer_t  *w,
                              char                       *var_name);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to histogram file writer.
 *
 * Options are:
 *   txt                 output txt (space-separated) files
 *   tex                 output TeX (TixZ) files
 *   png                 output PNG files
 *   [n_sub]             number of subdivisions
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque histogram writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void *
fvm_to_histogram_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t   time_dependency,
                             MPI_Comm                comm);

#else

void *
fvm_to_histogram_init_writer(const char             *name,
                             const char             *path,
                             const char             *options,
                             fvm_writer_time_dep_t  time_dependency);

#endif

/*----------------------------------------------------------------------------
 * Finalize FVM to histogram file writer.
 *
 * parameters:
 *   writer <-- pointer to opaque histogram writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_histogram_finalize_writer(void  *writer);

/*----------------------------------------------------------------------------
 * Associate new time step with a histogram geometry.
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_histogram_set_mesh_time(void          *writer,
                               const int      time_step,
                               const double   time_value);

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a histogram file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   writer           <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_histogram_export_field(void                  *writer,
                              const fvm_nodal_t     *mesh,
                              const char            *name,
                              fvm_writer_var_loc_t   location,
                              int                    dimension,
                              cs_interlace_t         interlace,
                              int                    n_parent_lists,
                              const cs_lnum_t        parent_num_shift[],
                              cs_datatype_t          datatype,
                              int                    time_step,
                              double                 time_value,
                              const void      *const field_values[]);

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * In this case, the effective writing to file is done.
 *
 * parameters:
 *   writer <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_histogram_flush(void  *writer);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TO_HISTOGRAM_H__ */
