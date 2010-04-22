#ifndef __FVM_FILE_H__
#define __FVM_FILE_H__

/*============================================================================
 * Parallel file I/O
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2007-2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * File hints and semantics
 */

#define FVM_FILE_NO_MPI_IO            (1 << 0)
#define FVM_FILE_NO_PREDISTRIBUTE     (1 << 1)

/* MPI-IO positioning semantics */

#define FVM_FILE_EXPLICIT_OFFSETS     (1 << 2)
#define FVM_FILE_INDIVIDUAL_POINTERS  (1 << 3)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* FVM file descriptor */

typedef struct _fvm_file_t  fvm_file_t;

/* FVM file modes */

typedef enum {

  FVM_FILE_MODE_READ,   /* Read mode */
  FVM_FILE_MODE_WRITE,  /* Write mode */
  FVM_FILE_MODE_APPEND  /* Append mode */

} fvm_file_mode_t;

/* Hints for file management */

typedef unsigned int fvm_file_hints_t;

/* Offset for FVM file position indicator (int64_t in C99) */

#if defined(SIZEOF_LONG_LONG)
typedef long long fvm_file_off_t;
#else
typedef long fvm_file_off_t;
#endif

/* Possibilities for the third argument of fvm_file_seek() */

typedef enum {

  FVM_FILE_SEEK_SET,   /* Seek from beginning of file */
  FVM_FILE_SEEK_CUR,   /* Seek from current position */
  FVM_FILE_SEEK_END    /* Seek from end of file */

} fvm_file_seek_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a file descriptor and open the associated file.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by fvm_file_set_swap_endian().
 *
 * parameters:
 *   name  <-- file name
 *   mode  <-- file acces mode: read, write, or append
 *   hints <-- file I/O hints (for MPI and MPI I/O behavior)
 *
 * returns:
 *   pointer to fvm_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

fvm_file_t *
fvm_file_open(const char         *name,
              fvm_file_mode_t     mode,
              fvm_file_hints_t    hints,
              MPI_Comm            comm);

#else

fvm_file_t *
fvm_file_open(const char         *name,
              fvm_file_mode_t     mode,
              fvm_file_hints_t    hints);

#endif

/*----------------------------------------------------------------------------
 * Destroy a file descriptor and close the associated file.
 *
 * parameters:
 *   f <-> file descriptor to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvm_file_t *
fvm_file_free(fvm_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's name.
 *
 * parameters:
 *   f <-- fvm_file_t descriptor
 *
 * returns:
 *   pointer to the file's name.
 *----------------------------------------------------------------------------*/

const char *
fvm_file_get_name(const fvm_file_t  *f);

/*----------------------------------------------------------------------------
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * parameters:
 *   f <-> fvm_file_t descriptor
 *----------------------------------------------------------------------------*/

void
fvm_file_set_big_endian(fvm_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's byte-swapping behavior.
 *
 * parameters:
 *   f <-- fvm_file_t descriptor
 *
 * returns:
 *   0 if file's endianness is the same as the system's, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_file_get_swap_endian(const fvm_file_t  *f);

/*----------------------------------------------------------------------------
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f    <-> fvm_file_t descriptor
 *   swap --> 1 if bytes must be swapped, 0 otherwise
 *----------------------------------------------------------------------------*/

void
fvm_file_set_swap_endian(fvm_file_t  *f,
                         int          swap);

/*----------------------------------------------------------------------------
 * Read data to a buffer, distributing it to all processes associated
 * with a file.
 *
 * parameters:
 *   f    <-- fvm_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_read_global(fvm_file_t  *f,
                     void        *buf,
                     size_t       size,
                     size_t       ni);

/*----------------------------------------------------------------------------
 * Write global data to a file.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The buffers on other ranks are ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   f    <-- fvm_file_t descriptor
 *   buf  <-- pointer to location containing data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_write_global(fvm_file_t  *f,
                      const void  *buf,
                      size_t       size,
                      size_t       ni);

/*----------------------------------------------------------------------------
 * Read data to a buffer, distributing a contiguous part of it to each
 * process associated with a file.
 *
 * Each process should receive a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_read_block(fvm_file_t  *f,
                    void        *buf,
                    size_t       size,
                    size_t       stride,
                    fvm_gnum_t   global_num_start,
                    fvm_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Write data to a file, each associated process providing a contiguous part
 * of this data.
 *
 * Each process should provide a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * This function may require an internal copy of the data to ensure that
 * the buffer contents are not modified, so if the buffer contents are
 * temporary values, to be deleted after writing, using
 * fvm_file_write_block_buffer() instead may be used to avoid an unneeded
 * memory allocation and copy.
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_write_block(fvm_file_t  *f,
                     void        *buf,
                     size_t       size,
                     size_t       stride,
                     fvm_gnum_t   global_num_start,
                     fvm_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Write data to a file, each associated process providing a contiguous part
 * of this data.
 *
 * Each process should provide a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * This function is intended to be used mainly data that is already of
 * copy of original data (such as data that has been redistributed across
 * processors just for the sake of output), or that is to be deleted after
 * writing, so it may modify the values in its input buffer (notably to
 * convert from little-endian to big-endian of vice-versa if necessary).
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_write_block_buffer(fvm_file_t  *f,
                            void        *buf,
                            size_t       size,
                            size_t       stride,
                            fvm_gnum_t   global_num_start,
                            fvm_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Update the file pointer according to whence.
 *
 * parameters:
 *   f      <-- fvm_file_t descriptor.
 *   offset <-- add to position specified to whence to obtain new position,
 *              measured in characters from the beginning of the file.
 *   whence <-- beginning if FVM_FILE_SEEK_SET, current if FVM_FILE_SEEK_CUR,
 *               or end-of-file if FVM_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise; currently, errors are fatal.
 *----------------------------------------------------------------------------*/

int
fvm_file_seek(fvm_file_t       *f,
              fvm_file_off_t    offset,
              fvm_file_seek_t   whence);

/*----------------------------------------------------------------------------
 * Return the position of the file pointer.
 *
 * When using MPI-IO with individual file pointers, we consider the file
 * pointer to be equal to the highest value of then individual file pointers.
 *
 * parameters:
 *   f <-- fvm_file_t descriptor
 *
 * returns:
 *   current position of the file pointer
 *----------------------------------------------------------------------------*/

fvm_file_off_t
fvm_file_tell(fvm_file_t  *f);

/*----------------------------------------------------------------------------
 * Get the default semantics for file access.
 *
 * returns:
 *   current default semantics for file access
 *----------------------------------------------------------------------------*/

fvm_file_hints_t
fvm_file_get_default_semantics(void);

/*----------------------------------------------------------------------------
 * Set the default semantics for file access.
 *
 * This may fail if semantics given contain incompatible values,
 * such as (FVM_FILE_EXPLICIT_OFFSETS | FVM_FILE_INDIVIDUAL_POINTERS),
 * or when setting MPI-IO access semantics when MPI-IO is not available.
 *
 * returns:
 *   0 if the semantics were valid, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
fvm_file_set_default_semantics(fvm_file_hints_t  hints);

/*----------------------------------------------------------------------------
 * Dump the metadata of a file structure in human readable form
 *
 * parameters:
 *   f <-- pointer to file
 *----------------------------------------------------------------------------*/

void
fvm_file_dump(const fvm_file_t  *f);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_FILE_H__ */
