#ifndef __BFT_FILE_H__
#define __BFT_FILE_H__

/*============================================================================
 * Base file wrapper type and associated functions
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2010  EDF

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

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if defined(__STDC_VERSION__)
#  if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdio.h>
#  endif
#else
#  include <stdio.h>
#endif

/* BFT headers */

#include "bft_error.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/* BFT file descriptor */

typedef struct _bft_file_t bft_file_t;

/* BFT file types */

typedef enum {

  BFT_FILE_TYPE_TEXT,          /* Text file */
  BFT_FILE_TYPE_BINARY,        /* Simple C binary file */
  BFT_FILE_TYPE_FORTRAN_BINARY /* Common Fortran binary file */

} bft_file_type_t;

/* BFT file modes */

typedef enum {

  BFT_FILE_MODE_READ,   /* Read mode */
  BFT_FILE_MODE_WRITE,  /* Write mode */
  BFT_FILE_MODE_APPEND  /* Append mode */

} bft_file_mode_t;

/* Offset for BFT file position indicator */

#if defined(BFT_SIZEOF_LONG_LONG)
typedef long long bft_file_off_t;
#else
typedef long bft_file_off_t;
#endif

/* Possibilities for the third argument of bft_file_seek() */

typedef enum {

  BFT_FILE_SEEK_SET,   /* Seek from beginning of file */
  BFT_FILE_SEEK_CUR,   /* Seek from current position */
  BFT_FILE_SEEK_END    /* Seek from end of file */

} bft_file_seek_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Create a `bft_file_t' file descriptor and open the associated file.
 *
 * The associated file is also opened. By default, data is written
 * or read as big-endian data. This behavior may be modified by
 * bft_file_set_swap_endian().
 *
 * parameters:
 *   name: <-- file name.
 *   mode: <-- file acces mode: read, write, or append.
 *   type: <-- file type: text, binary, or Fortran binary.
 *
 * returns:
 *   pointer to bft_file_t file descriptor (NULL in case of failure).
 */

bft_file_t *
bft_file_open(const char             *const name,
              const bft_file_mode_t         mode,
              const bft_file_type_t         type);

/*
 * Destroy a `bft_file_t' descriptor and close the associated file.
 *
 * The descriptor may only be destroyed if the file was successfully
 * closed. To force destruction of a bft_file_t descriptor even
 * if the associated file was not closed, use (bft_file_free_force()).
 *
 * The associated file is only closed if this was not already the case.
 *
 * parameters:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   pointer to bft_file_t file descriptor (NULL in case of success,
 *   f in case of failure).
 */

bft_file_t *
bft_file_free(bft_file_t  *f);

/*
 * Destroy a `bft_file_t' descriptor without closing its associated file.
 *
 * parameters:
 *   f:  -> bft_file_t descriptor.
 *
 * returns:
 *   NULL pointer.
 */

bft_file_t *
bft_file_free_descriptor(bft_file_t *f);

/*
 * Open `bft_file_t' descriptor's associated file.
 *
 * If the file is already open, this function does nothing.
 *
 * parameters:
 *  f:    <-- bft_file_t descriptor.
 *  mode: <-- file acces mode: read, write, or append.
 *
 * returns:
 *   0 in case of success, system error code in case of failure.
 */

int
bft_file_open_stream(bft_file_t       *const f,
                     bft_file_mode_t   const mode);

/*
 * Close a bft_file_t file descriptor's associated file.
 *
 * If the file is already closed, this function does nothing.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   0 in case of success, system error code in case of failure.
 */

int
bft_file_close_stream(bft_file_t  *const f);

/*
 * Test the end-of-file indicator for a given file.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   0 if the end-of-file has not been reached, or non-zero
 *   (1 or feof() return value) otherwise.
 */

int
bft_file_eof(const bft_file_t  *const f);

/*
 * Force write of all user-space buffered data for a given file.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   0 upon successful completion, system error code otherwise.
 */

int
bft_file_flush(bft_file_t  *const f);

/*
 * Obtain the current value of a file's position indicator.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   current value of the file's position indicator, or -1 in case of failure.
 */

bft_file_off_t
bft_file_tell(bft_file_t  *const f);

/*
 * Sets the file position indicator to the beginning of the file.
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 */

void
bft_file_rewind(bft_file_t  *const f);

/*
 * This function may call the libc's fseek() function. The C 99 standard draft
 * specifies that for a text file, the offset argument to fseek() should be
 * zero or a value returned by an earlier successful call to ftell()
 * (here bft_file_ftell()) on a stream (here a bft_file_t structure).
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * parameters:
 *   f:      <-- bft_file_t descriptor.
 *   offset: <-- add to position specified to whence to obtain new position,
 *               measured in characters from the beginning of the file.
 *   whence: <-- beginning if BFT_FILE_SEEK_SET, current if BFT_FILE_SEEK_CUR,
 *               or end-of-file if BFT_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise.
 */

int
bft_file_seek(bft_file_t             *const f,
              const bft_file_off_t          offset,
              const bft_file_seek_t         whence);

/*
 * Return a file's name.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   pointer to file's name.
 */

const char *
bft_file_get_name(const bft_file_t  *const f);

/*
 * Return a file's type.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   file's type.
 */

bft_file_type_t
bft_file_get_type(const bft_file_t  *const f);

/*
 * Change a file's type.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f:    <-> bft_file_t descriptor.
 *   type: <-- text, binary, or Fortran binary type descriptor.
 */

void
bft_file_set_type(bft_file_t             *const f,
                  const bft_file_type_t         type);

/*
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * By default, data is written or read in native format (as regards
 * big-endian or little-endian)..
 *
 * parameter:
 *   f: <-> bft_file_t descriptor.
 */

void
bft_file_set_big_endian(bft_file_t  *const f);

/*
 * Return a file's byte-swapping behavior.
 *
 * parameter:
 *   f: <-- bft_file_t descriptor.
 *
 * returns:
 *   0 if file's endianness is the same as the system's, 1 otherwise.
 */

int
bft_file_get_swap_endian(const bft_file_t  *const f);

/*
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f:    <-> bft_file_t descriptor.
 *   swap: <-- 1 if bytes must be swapped, 0 therwise.
 */

void
bft_file_set_swap_endian(bft_file_t  *const f,
                         const int          swap);

/*
 * Test a file's error or EOF condition.
 *
 * parameters:
 *   f:    <-- bft_file_t descriptor.
 *   line: <-- file line number if available, or 0.
 *
 * returns:
 *   0 if no error, system error code, or -1 if EOF.
 */

int
bft_file_read_check_error(const bft_file_t  *const f,
                          const int                line);

/*
 * Formatted output to a text file (as fprintf()).
 *
 * parameters:
 *   f:      <-- bft_file_t descriptor.
 *   format: <-- format string, as printf() and family.
 *   ... :   <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0'
 *   used to end output strings
 */

int
bft_file_printf(const bft_file_t  *const f,
                const char        *const format,
                ...);

/*
 * Formatted input from a text file (as fgets()).
 *
 * parameters:
 *   s:    --> buffer to which string is to be read.
 *   size: <-- maximum number of characters to be read plus one.
 *   f:    <-- bft_file_t descriptor.
 *   line: <-> file line number if available, or NULL.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 */

char *
bft_file_gets(char              *const s,
              const int                size,
              const bft_file_t  *const f,
              int               *const line);

/*
 * Formatted input from a text file if possible (as fgets()).
 *
 * This function is similar to bft_file_gets(), but failure to read
 * a line du to an end-of-file condition is not considered an error with
 * this variant, which may be used to read text files or sections thereof
 * of unknown length
 *
 * parameters:
 *   s:    --> buffer to which string is to be read.
 *   size: <-- maximum number of characters to be read plus one.
 *   f:    <-- bft_file_t descriptor.
 *   line: <-> file line number if available, or NULL.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 */

char *
bft_file_gets_try(char              *const s,
                  const int                size,
                  const bft_file_t  *const f,
                  int               *const line);

/*
 * Read a binary C or Fortran type record.
 *
 * A Fortran record compatible with most compilers is structured
 * as follows:
 *   - a 4-byte integer indicating the number of bytes in the record.
 *   - the raw data
 *   - a 4-byte integer indicating the number of bytes in the record.
 *
 * A C record contains only the raw data.
 *
 * parameters:
 *   rec:  --> pointer to location receiving data.
 *   size: <-- size of each item of data in bytes.
 *   ni:   <-- number of items to read.
 *   f:    <-- bft_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran
 *   record, if the whole record could not be read, returns 0.
 */

size_t
bft_file_read(void              *const rec,
              const size_t             size,
              const size_t             ni,
              const bft_file_t  *const f);

/*
 * Read a binary C or Fortran type record.
 *
 * This function is similar to bft_file_read(), but failure to read
 * a record due to an end-of-file condition is not considered an error with
 * this variant, which may be used to read records whose presence in the
 * file is unknown.
 *
 * A Fortran record compatible with most compilers is structured
 * as follows:
 *   - a 4-byte integer indicating the number of bytes in the record.
 *   - the raw data
 *   - a 4-byte integer indicating the number of bytes in the record.
 *
 * A C record contains only the raw data.
 *
 * parameters:
 *   rec:  --> pointer to location receiving data.
 *   size: <-- size of each item of data in bytes.
 *   ni:   <-- number of items to read.
 *   f:    <-- bft_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran
 *   record, if the whole record could not be read, returns 0.
 */

size_t
bft_file_read_try(void              *const rec,
                  const size_t             size,
                  const size_t             ni,
                  const bft_file_t  *const f);

/*
 * Write a binary C or Fortran type record.
 *
 * A Fortran record compatible with most compilers is structured
 * as follows:
 *   - a 4-byte integer indicating the number of bytes in the record.
 *   - the raw data
 *   - a 4-byte integer indicating the number of bytes in the record.
 *
 * A C record contains only the raw data.
 *
 * parameters:
 *   rec:  <-- pointer to location containing data.
 *   size: <-- size of each item of data in bytes.
 *   ni:   <-- number of items to write.
 *   f:    <-- bft_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully written.
 */

size_t
bft_file_write(const void        *const rec,
               const size_t             size,
               const size_t             ni,
               const bft_file_t  *const f);

/*
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * parameters:
 *   dest: --> pointer to converted data location.
 *   src:  <-- pointer to source data location.
 *   size: <-- size of each item of data in bytes.
 *   ni:   <-- number of data items.
 */

void
bft_file_swap_endian(void *const          dest,
                     const void    *const src,
                     const size_t         size,
                     const size_t         ni);

/*
 * Create a new directory using default permissions.
 *
 * This function is similar to the POSIX function mkdir(), except that
 * it has no "mode" argument: by default, on a POSIX type system,
 * permissions include read, write, and execute access for the user,
 * group and others, modified by the users umask value (so with a
 * typical configuration, the user will have read, write, and execute
 * pemission, the group and others will only have read and execute
 * permission, but this behavior may be modified).
 *
 * Also, contrary to the usual mkdir(), if the directory already
 * exists (and is truly a directory), this is considered a success
 * and not a failure, and 0 is returned: the aim of this function
 * is to make a directory available, so if it already exists,
 * this is considered acceptable.
 *
 * parameters:
 *   pathname: <-- name of new directory.
 *
 * returns:
 *   0 on success, -1 if an error occured (in which case errno
 *   contains the appropriate error code). If the underlying
 *   system has no mkdir() function or it was not detected
 *   upon BFT configuration, 1 is returned.
 */

int
bft_file_mkdir_default(const char  *const pathname);

/*
 * Check if a file exists and is a regular file.
 *
 * parameters:
 *   name: <-- file name.
 *
 * returns:
 *   1 if file exists and is a regular file, 0 otherwise.
 */

int
bft_file_isreg(const char  *const name);

/*
 * Check if a directory exists.
 *
 * parameters:
 *   name: <-- directory name.
 *
 * returns:
 *   1 if directory exists, 0 otherwise.
 */

int
bft_file_isdir(const char  *const name);

/*
 * List files inside a directory.
 *
 * The array returned must be freed by the caller using BFT_FREE,
 * as well as the individual entries in the array.
 *
 * parameters:
 *   pathname <-- name of directory.
 *
 * returns:
 *   an array of file names in a directory. The last entry is set to NULL.
 *   If no means to list the directory or an error occured, the return
 *    value is simply NULL.
 */

char **
bft_file_listdir(const char  *pathname);

/* Returns the error handler associated with the bft_file_...() functions.
 *
 * returns:
 *   pointer to the error handler function.
 */

bft_error_handler_t *
bft_file_error_handler_get(void);

/*
 * Associates an error handler with the bft_file_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameter:
 *   handler: <-- pointer to the error handler function.
 */

void
bft_file_error_handler_set(bft_error_handler_t *const handler);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFT_FILE_H__ */
