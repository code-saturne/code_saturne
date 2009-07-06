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

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_file.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_IO_NAME_LEN   32    /* Section header name length */

#define CS_IO_ECHO_NONE -2        /* No verbosity at all */
#define CS_IO_ECHO_OPEN_CLOSE -1  /* Echo open or close operations */
#define CS_IO_ECHO_HEADERS 0      /* Echo headers */

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

  const char     *sec_name;           /* Pointer to section name */
  fvm_file_off_t  n_vals;             /* Number of associated values */
  size_t          location_id;        /* Id of associated location, or 0 */
  size_t          index_id;           /* Id of associated index, or 0 */
  size_t          n_location_vals;    /* Number of values per location */
  fvm_datatype_t  elt_type;           /* Type if n_elts > 0 */
  fvm_datatype_t  type_read;          /* Type in file */

} cs_io_sec_header_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Default hints for files using this API (for MPI-IO) */

extern int       cs_glob_io_hints;

/* Global structure associated with pre-processor output (kernel input) file */

extern cs_io_t  *cs_glob_pp_io;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set the default kernel IO hints to the specified value.
 *
 * Fortran interface :
 *
 * SUBROUTINE IOHINT (IHINT)
 * *****************
 *
 * INTEGER          IHINT       : <-> : IO hints (bit mask)
 *                                        0: default
 *                                        1: disable MPI IO
 *                                        4: MPI IO uses explicit offsets
 *                                        8: MPI IO uses individual pointers
 *----------------------------------------------------------------------------*/

void
CS_PROCF (iohint, IOHINT) (const cs_int_t  *iopt);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a kernel IO file structure.
 *
 * The magic string may be NULL only in read mode;
 *
 * If the position of section bodies is already known (after initial
 * analysis for example), the file may be opened for reading section bodies
 * only by using "seek_read_section_bodies_only" as a magic string. This may
 * be used to map another type of file to kernel io files, if header data is
 * different but body data is similar (binary, using the same datatypes).
 *
 * parameters:
 *   name         <-- file name
 *   magic_string <-- magic string associated with file type, or NULL
 *   mode         <-- read or write
 *   hints        <-- optional flags for file access method (see fvm_file.h)
 *   echo         <-- echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)

cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 int            hints,
                 long           echo,
                 MPI_Comm       comm);

#else

cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 int            hints,
                 long           echo);

#endif /* FVM_HAVE_MPI */

/*----------------------------------------------------------------------------
 * Initialize a kernel IO file structure in read mode, building an index.
 *
 * The magic string may be NULL, if we choose to ignore it.
 *
 * parameters:
 *   name         <-- file name
 *   magic_string <-- magic string associated with file type
 *   hints        <-- optional flags for file access method (see fvm_file.h)
 *   echo         <-- echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)

cs_io_t *
cs_io_initialize_with_index(const char    *file_name,
                            const char    *magic_string,
                            int            hints,
                            long           echo,
                            MPI_Comm       comm);

#else

cs_io_t *
cs_io_initialize_with_index(const char    *file_name,
                            const char    *magic_string,
                            int            hints,
                            long           echo);

#endif /* FVM_HAVE_MPI */

/*----------------------------------------------------------------------------
 * Free a preprocessor output file structure, closing the associated file.
 *
 * parameters:
 *   pp_io <-> kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_finalize(cs_io_t **pp_io);

/*----------------------------------------------------------------------------
 * Return a pointer to a preprocessor IO structure's name.
 *
 * parameters:
 *   pp_io <-- kernel IO structure
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_name(const cs_io_t  *pp_io);

/*----------------------------------------------------------------------------
 * Return the number of indexed entries in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *
 * returns:
 *   size of index if present, 0 otherwise,
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_index_size(const cs_io_t  *inp);

/*----------------------------------------------------------------------------
 * Return the name of an indexed section in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *   id  <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   pointer to section name if id in index range, NULL otherwise
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_indexed_sec_name(const cs_io_t  *inp,
                           size_t          id);

/*----------------------------------------------------------------------------
 * Return header data for an indexed section in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *   id  <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   section header data (if id not in index range, fields set to zero)
 *----------------------------------------------------------------------------*/

cs_io_sec_header_t
cs_io_get_indexed_sec_header(const cs_io_t  *inp,
                             size_t          id);

/*----------------------------------------------------------------------------
 * Return a kernel IO structure's echo (verbosity) level.
 *
 * parameters:
 *   pp_io <-- kernel IO structure
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_echo(const cs_io_t  *pp_io);

/*----------------------------------------------------------------------------
 * Read a message header.
 *
 * parameters:
 *   pp_io  <-- kernel IO structure
 *   header --> header structure
 *
 * returns:
 *   0 if a header was read, 1 in case of error or end-of-file
 *----------------------------------------------------------------------------*/

int
cs_io_read_header(cs_io_t             *inp,
                  cs_io_sec_header_t  *header);

/*----------------------------------------------------------------------------
 * Set a kernel IO's state so as to be ready to read an indexed section.
 *
 * The header values and position in the file are set so as to be equivalent
 * to those they would have if the corresponding header had just been read.
 *
 * parameters:
 *   inp    <-> input kernel IO structure
 *   header --> associated header
 *   id     <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   0 in case of success, 1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_io_set_indexed_position(cs_io_t             *inp,
                           cs_io_sec_header_t  *header,
                           size_t               id);

/*----------------------------------------------------------------------------
 * Set a message's final data type to fvm_lnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> kernel IO structure
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
 *   header <-> header structure
 *   pp_io  <-- kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_gnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *pp_io);

/*----------------------------------------------------------------------------
 * Check that a message's final data type corresponds to cs_real_t.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  <-- kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_assert_cs_real(const cs_io_sec_header_t  *header,
                     const cs_io_t             *pp_io);

/*----------------------------------------------------------------------------
 * Read a message body and replicate it to all processors.
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
 *   pp_io            --> kernel IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_global(const cs_io_sec_header_t  *header,
                  void                      *elts,
                  cs_io_t                   *pp_io);

/*----------------------------------------------------------------------------
 * Read a message body, assigning a different block to each processor.
 *
 * If location_id > 0 and header->n_location_vals > 1, then
 * global_num_start and global_num_end will be based on location element
 * numbers, so the total number of values read equals
 * (global_num_end - global_num_start) * header->n_location_vals.
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
 *   pp_io            --> kernel IO structure
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
 *   pp_io            --> kernel IO structure
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

/*----------------------------------------------------------------------------
 * Write a global section.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The section data on other ranks is ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   section_name     <-- section name
 *   n_vals           <-- total number of values
 *   location_id      <-- id of associated location, or 0
 *   index_id         <-- id of associated index, or 0
 *   n_location_vals  <-- number of values per location
 *   elt_type         <-- element type
 *   elts             <-- pointer to element data
 *   outp             <-> output kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_write_global(const char      *sec_name,
                   fvm_gnum_t       n_vals,
                   size_t           location_id,
                   size_t           index_id,
                   size_t           n_location_vals,
                   fvm_datatype_t   elt_type,
                   const void      *elts,
                   cs_io_t         *outp);

/*----------------------------------------------------------------------------
 * Write a section to file, each associated process providing a contiguous
 * of the section's body.
 *
 * Each process should provide a (possibly empty) block of the body,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * If location_id > 0 and n_location_vals > 1, then global_num_start
 * and global_num_end will be based on location element numbers, so the
 * total number of values read equals
 * (global_num_end - global_num_start) * header->n_location_vals.
 *
 * This function is intended to be used mainly data that is already of
 * copy of original data (such as data that has been redistributed across
 * processors just for the sake of output), or that is to be deleted after
 * writing, so it may modify the values in its input buffer (notably to
 * convert from little-endian to big-endian of vice-versa if necessary).
 *
 * parameters:
 *   section_name     <-- section name
 *   n_g_elts         <-- number of global elements (locations)
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *   location_id      <-- id of associated location, or 0
 *   index_id         <-- id of associated index, or 0
 *   n_location_vals  <-- number of values per location
 *   elt_type         <-- element type
 *                        (1 to n numbering)
 *   elts             <-- pointer to element data
 *   outp             <-> output kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_write_block_buffer(const char      *sec_name,
                         fvm_gnum_t       n_g_elts,
                         fvm_gnum_t       global_num_start,
                         fvm_gnum_t       global_num_end,
                         size_t           location_id,
                         size_t           index_id,
                         size_t           n_location_vals,
                         fvm_datatype_t   elt_type,
                         void            *elts,
                         cs_io_t         *outp);

/*----------------------------------------------------------------------------
 * Dump a kernel IO file handle's metadata.
 *
 * parameters:
 *   cs_io  <-- kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_dump(const cs_io_t  *cs_io);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_IO_H__ */
