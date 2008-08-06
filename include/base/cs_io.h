/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_IO_H__
#define __CS_IO_H__

/*============================================================================
 *  Low level file I/O utility functions for Preprocessor and restart files
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_IO_NAME_LEN   32    /* Section header nam length */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*  Input or output mode */

  typedef enum {

  CS_IO_MODE_READ,
  CS_IO_MODE_WRITE

} cs_io_mode_t;

/* Structure associated with opaque pre-processing structure object */

typedef struct _cs_io_t cs_io_t;

/* Structure used to save section header data, so as to simplify
   passing this data to various functions */

typedef struct {

  char            sec_name[CS_IO_NAME_LEN + 1];     /* Section name */
  char            type_read_name[3];                /* Name of type in file */
  fvm_gnum_t      n_elts;                           /* Number of elements */
  fvm_datatype_t  elt_type;                         /* Type if n_elts > 0 */
  fvm_datatype_t  type_read;                        /* Type in file */

} cs_io_sec_header_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Global structure associated with pre-processor output (kernel input) file */

extern cs_io_t  *cs_glob_pp_io;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a preprocessor output file structure.
 *
 * parameters:
 *   name         --> file name
 *   magic_string --> magic string associated with file type
 *   mode         --> read or write
 *   echo         --> echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *
 * returns:
 *   pointer to preprocessor IO structure
 *----------------------------------------------------------------------------*/

cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 size_t         echo);

/*----------------------------------------------------------------------------
 * Free a preprocessor output file structure, closing the associated file.
 *
 * parameters:
 *   pp_io <-> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_finalize(cs_io_t **pp_io);

/*----------------------------------------------------------------------------
 * Return a pointer to a preprocessor IO structure's name.
 *
 * parameters:
 *   pp_io --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_name(const cs_io_t  *pp_io);

/*----------------------------------------------------------------------------
 * Return a preprocessor IO structure's echo (verbosity) level.
 *
 * parameters:
 *   pp_io --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_echo(const cs_io_t  *pp_io);

/*----------------------------------------------------------------------------
 * Read a message header.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_read_header(cs_io_sec_header_t  *header,
                  cs_io_t             *pp_io);

/*----------------------------------------------------------------------------
 * Set a message's final data type to fvm_lnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_lnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *pp_io);

/*----------------------------------------------------------------------------
 * Set a message's final data type to fvm_gnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_gnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *pp_io);

/*----------------------------------------------------------------------------
 * Check that a message's final data type corresponds to cs_real_t.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_assert_cs_real(const cs_io_sec_header_t  *header,
                     const cs_io_t             *pp_io);

/*----------------------------------------------------------------------------
 * Read a message body and replicate it to all processors.
 *
 * If global_num_start and global_num_end are > 0, a different block is
 * assigned to each processor. Otherwise, the full data is replicated
 * for each processor.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   elts             <-> pointer to data array, or NULL
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_global(const cs_io_sec_header_t  *header,
                  void                         *elts,
                  cs_io_t                   *pp_io);

/*----------------------------------------------------------------------------
 * Read a message body, assigning a different block to each processor.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elts             <-> pointer to data array, or NULL
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_block(const cs_io_sec_header_t  *header,
                 fvm_gnum_t                 global_num_start,
                 fvm_gnum_t                 global_num_end,
                 void                      *elts,
                 cs_io_t                   *pp_io);

/*----------------------------------------------------------------------------
 * Read a message body, assigning a different block to each processor,
 * when the body corresponds to an index.
 *
 * In serial mode, this function behaves just like cs_io_read_block(),
 * except that it allows only unsigned integer values (fvm_gnum_t).
 *
 * In parallel mode, global_num_end should be set to the past-the-end value
 * of the base data block, the same as for regular data (and not increased
 * by 1 for the last rank, as this will be handled internally).
 * On each rank, the buffer size should be:
 * global_num_end - global_num_start + 1, as the past-the end index
 * for the local block is added automatically.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elts             <-> pointer to data array, or NULL
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_index_block(cs_io_sec_header_t  *header,
                       fvm_gnum_t           global_num_start,
                       fvm_gnum_t           global_num_end,
                       fvm_gnum_t          *elts,
                       cs_io_t             *pp_io);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_IO_H__ */
