/*============================================================================
 * User definition of time tables
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_time_table.c
 *
 * \brief User definitions of time tables
 *
 * See \ref cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time tables using C API.
 * This function is called at the begin of the simulation only.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_table()
{
  /*!< [define_time_table1] */
  {
    // Define Table and headers
    cs_time_table_t *table = cs_time_table_from_csv_file_simple("table1",
                                                                "data.csv",
                                                                ",");

    const char *_headers[] = {"time", "flowrate"};
    cs_time_table_set_headers(table, 2, _headers);
  }
  /*!< [define_time_table1] */

  /*!< [define_time_table2] */
  {
    // Define Table and skip 2 first lines of file
    cs_time_table_t *table =
      cs_time_table_from_csv_file_simple_headers("table2",    /* table name */
                                                 "file2.csv", /* data file */
                                                 ";",         /* Separator */
                                                 2);          /* Number of header lines to skip */

    // Optionnal: add headers
    const char *_headers[] = {"time", "temperature"};
    cs_time_table_set_headers(table, 2, _headers);
  }
  /*!< [define_time_table2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
