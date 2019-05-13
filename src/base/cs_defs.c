/*============================================================================
 * Base macro and typedef definitions for system portability
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

#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_defs.c
        Base macro and typedef definitions for system portability.

  \typedef cs_gnum_t
           \brief global mesh entity number
           \details Global mesh-entity numbers are strictly positive
           (<em>1 to n</em> based) integers, so they are declared as a form of
           unsigned integer. Such a number is unique across MPI ranks;
           2 mesh elements on different ranks share the same number
           <em>if and only if</em> they are local instances of the same global
           element (such as shared faces or vertices on rank boundaries).
           A value of 0 is commonly used to mark undefined (or not yet
           defined) element ids in various pre or post-processing stages.

  \typedef cs_lnum_t
           \brief local mesh entity id
           \details Local mesh-entity ids are signed integers, and be either
           <em>0 to n-1</em> or <em>1 to n</em> based. When 0-based,
           the \e id prefix or postfix is preferred
           for associated variable names, while \e num is preferred when
           1-based.
           In C code, using this type is recommended over using simple
           \c int integers, as 64-bit variants could be used in the future
           for shared-memory machines with large memory. This type should
           \b not be used to declare identifiers which are not mesh entities,
           such as groups, fields, or any other entity whose count does
           not depend on mesh size, so as not to pollute the readability
           and comprehensibility of the code.

  \typedef cs_int_t
            \brief Fortran-compatible integer
            \deprecated
            Currently, this integer type is necessarily of the same
            type as \ref cs_lnum_t, but it should only be used in Fortran
            wrapper function definitions. Moving to ISO_C_BINDINGS,
            and converting more code to C, this type should eventually
            disappear.

  \typedef cs_lnum_2_t
           \brief vector of 2 local mesh-entity ids

  \typedef cs_real_t
           \brief Floating-point value

  \typedef cs_real_2_t
           \brief vector of 2 floating-point values

  \typedef cs_real_3_t
           \brief vector of 3 floating-point values

  \typedef cs_real_4_t
           \brief vector of 4 floating-point values

  \typedef cs_real_6_t
           \brief vector of 6 floating-point values

  \typedef cs_real_33_t
           \brief 3x3 matrix of floating-point values

  \typedef cs_real_66_t
           \brief 6x6 matrix of floating-point values

  \typedef cs_real_332_t
           \brief vector of 2 3x3 matrices of floating-point values

*/

/*! \fn inline static cs_lnum_t cs_align(cs_lnum_t i, cs_lnum_t m)
 *
 * \brief Given a base index i, return the next index aligned with a size m.
 *
 * This index is computed as follows:
 *
 * if i > 0:
 *   ((i - 1) / m + 1) * m
 * if i = 0:
 *   0
 *
 * \param[in]  i  base index
 * \param[in]  m  block size to align with
 *
 * \return  aligned index
 */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Sizes associated with datatypes */

const size_t  cs_datatype_size[] = {0,
                                    1,
                                    sizeof(float),
                                    sizeof(double),
                                    sizeof(unsigned short int),
                                    4,
                                    8,
                                    4,
                                    8};

const char   *cs_datatype_name[] = {"",
                                    "char",
                                    "float",
                                    "double",
                                    "flag",
                                    "int32",
                                    "int64",
                                    "uint32",
                                    "uint64"};

#if defined(HAVE_MPI)

/* MPI Datatypes associated with Code_Saturne datatypes */

MPI_Datatype  cs_datatype_to_mpi[] = {MPI_DATATYPE_NULL,
                                      MPI_CHAR,
                                      MPI_FLOAT,
                                      MPI_DOUBLE,
                                      MPI_UNSIGNED_SHORT,
                                      MPI_INT,            /* CS_INT32 */
                                      MPI_LONG_INT,       /* CS_INT64 */
                                      MPI_UNSIGNED,       /* CS_UINT32 */
                                      MPI_UNSIGNED_LONG}; /* CS_UINT64 */


#endif

/* Global variables indicationg task state */

int  cs_glob_n_threads = 1;    /* Number of threads */

int  cs_glob_rank_id = -1;     /* Rank of process in communicator */
int  cs_glob_n_ranks =  1;     /* Number of processes in communicator */

#if defined(HAVE_MPI)

MPI_Comm  cs_glob_mpi_comm = MPI_COMM_NULL;   /* Main MPI intra-communicator */

#endif

/*=============================================================================
 * Local static variables
 *============================================================================*/

/*=============================================================================
 * Public functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS
