#ifndef __CS_TIME_PLOT_H__
#define __CS_TIME_PLOT_H__

/*============================================================================
 * Time_Plot helper structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

typedef struct _cs_time_plot_t  cs_time_plot_t;

/*============================================================================
 * Local type definitions
 *============================================================================*/

/* Type of 1D plot file format */

typedef enum {
  CS_TIME_PLOT_DAT,  /* .dat file (usable by Qtplot or Grace) */
  CS_TIME_PLOT_CSV   /* .csv file (readable by ParaView or spreadsheat) */
} cs_time_plot_format_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a plot file writer for probe-type plots
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   plot_name        <-- plot (variable) name
 *   file_prefix      <-- file name prefix
 *   format           <-- associated file format
 *   use_iteration    <-- should we use the iteration number instead of the
 *                        physical time ?
 *   flush_wtime      <-- elapsed time interval between file flushes
 *                        (if < 0, no forced flush)
 *   n_buffer_steps   <-- number of time steps in output buffer if
 *                        file is not to be kept open
 *   n_probes         <-- number of probes associated with this plot
 *   probe_list       <-- numbers (1 to n) of probes if filtered, or NULL
 *   probe_coords     <-- probe coordinates, or NULL
 *   probe_names      <-- probe names, or NULL
 *
 * returns:
 *   pointer to new time plot writer
 *----------------------------------------------------------------------------*/

cs_time_plot_t *
cs_time_plot_init_probe(const char             *plot_name,
                        const char             *file_prefix,
                        cs_time_plot_format_t   format,
                        bool                    use_iteration,
                        double                  flush_wtime,
                        int                     n_buffer_steps,
                        int                     n_probes,
                        const int              *probe_list,
                        const cs_real_t         probe_coords[],
                        const char             *probe_names[]);

/*----------------------------------------------------------------------------
 * Initialize a plot file writer for structure-type plots
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   plot_name          <-- plot (variable) name
 *   file_prefix        <-- file name prefix
 *   format             <-- associated file format
 *   use_iteration      <-- should we use the iteration number instead of the
 *                          physical time ?
 *   flush_wtime        <-- elapsed time interval between file flushes
 *                          (if < 0, no forced flush)
 *   n_buffer_steps     <-- number of time steps in output buffer if
 *                          file is not to be kept open
 *   n_structures       <-- number of structures associated with this plot
 *   mass_matrixes      <-- mass matrix coefficients (3x3 blocks)
 *   damping_matrixes   <-- damping matrix coefficients (3x3 blocks)
 *   stiffness_matrixes <-- stiffness matrix coefficients (3x3 blocks)
 *
 * returns:
 *   pointer to new time plot writer
 *----------------------------------------------------------------------------*/

cs_time_plot_t *
cs_time_plot_init_struct(const char             *plot_name,
                         const char             *file_prefix,
                         cs_time_plot_format_t   format,
                         bool                    use_iteration,
                         double                  flush_wtime,
                         int                     n_buffer_steps,
                         int                     n_structures,
                         const cs_real_t         mass_matrixes[],
                         const cs_real_t         damping_matrixes[],
                         const cs_real_t         stiffness_matrixes[]);

/*----------------------------------------------------------------------------
 * Finalize time plot writer for a given variable
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   p <-> time plot values file handler
 *----------------------------------------------------------------------------*/

void
cs_time_plot_finalize(cs_time_plot_t  **p);

/*----------------------------------------------------------------------------
 * Write time plot values
 *
 * This function should only be called by one rank for a given data series.
 *
 * parameters:
 *   p      <-- pointer to associated plot structure
 *   tn     <-- associated time step number
 *   t      <-- associated time value
 *   n_vals <-- number of associated time values
 *   vals   <-- associated time values
 *----------------------------------------------------------------------------*/

void
cs_time_plot_vals_write(cs_time_plot_t  *p,
                        int              tn,
                        double           t,
                        int              n_vals,
                        const cs_real_t  vals[]);

/*----------------------------------------------------------------------------
 * Flush buffered values to file if applicable
 *
 * parameters:
 *   p <-> time plot values file handler
 *----------------------------------------------------------------------------*/

void
cs_time_plot_flush(cs_time_plot_t  *p);

/*----------------------------------------------------------------------------
 * flush all time plots
 *----------------------------------------------------------------------------*/

void
cs_time_plot_flush_all(void);

/*----------------------------------------------------------------------------
 * Set time plot file writer flush behavior defaults.
 *
 * parameters:
 *   flush_wtime     <-- elapsed time interval between file flushes;
 *                       if < 0, no forced flush
 *   n_buffer_steps  <-- number of time steps in output buffer if
 *                       file is not to be kept open
 *----------------------------------------------------------------------------*/

void
cs_time_plot_set_flush_default(float  flush_wtime,
                               int    n_buffer_steps);

/*----------------------------------------------------------------------------
 * Return time plot file writer flush behavior defaults.
 *
 * parameters:
 *   flush_wtime     --> elapsed time interval between file flushes;
 *                       if < 0, no forced flush (NULL if not queried)
 *   n_buffer_steps  <-- number of time steps in output buffer if
 *                       file is not to be kept open (NULL if not queried)
 *----------------------------------------------------------------------------*/

void
cs_time_plot_get_flush_default(float  *flush_wtime,
                               int    *n_buffer_steps);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROBE_H__ */
