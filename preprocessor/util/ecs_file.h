#ifndef __ECS_FILE_H__
#define __ECS_FILE_H__

/*============================================================================
 * Base file wrapper type and associated functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "ecs_def.h"

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/* ECS file descriptor */

typedef struct _ecs_file_t ecs_file_t;

/* ECS file types */

typedef enum {

  ECS_FILE_TYPE_TEXT,          /* Text file */
  ECS_FILE_TYPE_BINARY,        /* Simple C binary file */
  ECS_FILE_TYPE_FORTRAN_BINARY /* Common Fortran binary file */

} ecs_file_type_t;

/* ECS file modes */

typedef enum {

  ECS_FILE_MODE_READ,   /* Read mode */
  ECS_FILE_MODE_WRITE,  /* Write mode */
  ECS_FILE_MODE_APPEND  /* Append mode */

} ecs_file_mode_t;

/* Offset for ECS file position indicator */

#if defined(SIZEOF_LONG_LONG)
typedef long long ecs_file_off_t;
#else
typedef long ecs_file_off_t;
#endif

/* Possibilities for the third argument of ecs_file_seek() */

typedef enum {

  ECS_FILE_SEEK_SET,   /* Seek from beginning of file */
  ECS_FILE_SEEK_CUR,   /* Seek from current position */
  ECS_FILE_SEEK_END    /* Seek from end of file */

} ecs_file_seek_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Create a `ecs_file_t' file descriptor and open the associated file.
 *
 * The associated file is also opened. By default, data is written
 * or read as big-endian data. This behavior may be modified by
 * ecs_file_set_swap_endian().
 *
 * parameters:
 *   name: <-- file name.
 *   mode: <-- file acces mode: read, write, or append.
 *   type: <-- file type: text, binary, or Fortran binary.
 *
 * returns:
 *   pointer to ecs_file_t file descriptor (NULL in case of failure).
 */

ecs_file_t *
ecs_file_open(const char             *name,
              const ecs_file_mode_t   mode,
              const ecs_file_type_t   type);

/*
 * Destroy a `ecs_file_t' descriptor and close the associated file.
 *
 * The descriptor may only be destroyed if the file was successfully
 * closed. To force destruction of a ecs_file_t descriptor even
 * if the associated file was not closed, use (ecs_file_free_force()).
 *
 * The associated file is only closed if this was not already the case.
 *
 * parameters:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   pointer to ecs_file_t file descriptor (NULL in case of success,
 *   f in case of failure).
 */

ecs_file_t *
ecs_file_free(ecs_file_t  *f);

/*
 * Destroy a `ecs_file_t' descriptor without closing its associated file.
 *
 * parameters:
 *   f:  -> ecs_file_t descriptor.
 *
 * returns:
 *   NULL pointer.
 */

ecs_file_t *
ecs_file_free_descriptor(ecs_file_t *f);

/*
 * Open `ecs_file_t' descriptor's associated file.
 *
 * If the file is already open, this function does nothing.
 *
 * parameters:
 *  f:    <-- ecs_file_t descriptor.
 *  mode: <-- file acces mode: read, write, or append.
 *
 * returns:
 *   0 in case of success, system error code in case of failure
 *   (or Zlib error code in case of Zlib memory allocation problem
 *   for a gzipped file).
 */

int
ecs_file_open_stream(ecs_file_t       *f,
                     ecs_file_mode_t   mode);

/*
 * Close a ecs_file_t file descriptor's associated file.
 *
 * If the file is already closed, this function does nothing.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   0 in case of success, system error code in case of failure
 *   (or Zlib error code in case of a Zlib specific error
 *   for a gzipped file).
 */

int
ecs_file_close_stream(ecs_file_t  *f);

/*
 * Test the end-of-file indicator for a given file.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   0 if the end-of-file has not been reached, or non-zero
 *   (1 or feof() return value) otherwise.
 */

int
ecs_file_eof(const ecs_file_t  *f);

/*
 * Force write of all user-space buffered data for a given file.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   0 upon successful completion, system error code otherwise.
 */

int
ecs_file_flush(ecs_file_t  *f);

/*
 * Obtain the current value of a file's position indicator.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   current value of the file's position indicator, or -1 in case of failure.
 */

ecs_file_off_t
ecs_file_tell(ecs_file_t  *f);

/*
 * Sets the file position indicator to the beginning of the file.
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 */

void
ecs_file_rewind(ecs_file_t  *f);

/*
 * This function may call the libc's fseek() function, or Zlib's gzseek()
 * function. The C 99 standard draft specifies that for a text file, the offset
 * argument to fseek() should be zero or a value returned by an earlier
 * successful call to ftell() (here ecs_file_ftell()) on a stream (here a
 * ecs_file_t structure). Zlib's gzseek() does not support SEEK_END, at least
 * as of version 1.2.1.
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * parameters:
 *   f:      <-- ecs_file_t descriptor.
 *   offset: <-- add to position specified to whence to obtain new position,
 *               measured in characters from the beginning of the file.
 *   whence: <-- beginning if ECS_FILE_SEEK_SET, current if ECS_FILE_SEEK_CUR,
 *               or end-of-file if ECS_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise.
 */

int
ecs_file_seek(ecs_file_t             *f,
              const ecs_file_off_t    offset,
              const ecs_file_seek_t   whence);

/*
 * Return a file's name.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   pointer to file's name.
 */

const char *
ecs_file_get_name(const ecs_file_t  *f);

/*
 * Return a file's type.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   file's type.
 */

ecs_file_type_t
ecs_file_get_type(const ecs_file_t  *f);

/*
 * Change a file's type.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f:    <-> ecs_file_t descriptor.
 *   type: <-- text, binary, or Fortran binary type descriptor.
 */

void
ecs_file_set_type(ecs_file_t             *f,
                  const ecs_file_type_t   type);

/*
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * By default, data is written or read in native format (as regards
 * big-endian or little-endian)..
 *
 * parameter:
 *   f: <-> ecs_file_t descriptor.
 */

void
ecs_file_set_big_endian(ecs_file_t  *f);

/*
 * Return a file's byte-swapping behavior.
 *
 * parameter:
 *   f: <-- ecs_file_t descriptor.
 *
 * returns:
 *   0 if file's endianness is the same as the system's, 1 otherwise.
 */

int
ecs_file_get_swap_endian(const ecs_file_t  *f);

/*
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f:    <-> ecs_file_t descriptor.
 *   swap: <-- 1 if bytes must be swapped, 0 therwise.
 */

void
ecs_file_set_swap_endian(ecs_file_t  *f,
                         const int    swap);

/*
 * Test a file's error or EOF condition.
 *
 * parameters:
 *   f:    <-- ecs_file_t descriptor.
 *   line: <-- file line number if available, or 0.
 *
 * returns:
 *   0 if no error, system error code, or -1 if EOF.
 */

int
ecs_file_read_check_error(const ecs_file_t  *f,
                          const int          line);

/*
 * Formatted output to a text file (as fprintf()).
 *
 * parameters:
 *   f:      <-- ecs_file_t descriptor.
 *   format: <-- format string, as printf() and family.
 *   ... :   <-- variable arguments based on format string.
 *
 * returns:
 *   number of characters printed, not counting the trailing '\0'
 *   used to end output strings
 */

int
ecs_file_printf(const ecs_file_t  *const f,
                const char        *const format,
                ...);

/*
 * Formatted input from a text file (as fgets()).
 *
 * parameters:
 *   s:    --> buffer to which string is to be read.
 *   size: <-- maximum number of characters to be read plus one.
 *   f:    <-- ecs_file_t descriptor.
 *   line: <-> file line number if available, or NULL.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 */

char *
ecs_file_gets(char              *s,
              const int          size,
              const ecs_file_t  *f,
              int               *line);

/*
 * Formatted input from a text file if possible (as fgets()).
 *
 * This function is similar to ecs_file_gets(), but failure to read
 * a line du to an end-of-file condition is not considered an error with
 * this variant, which may be used to read text files or sections thereof
 * of unknown length
 *
 * parameters:
 *   s:    --> buffer to which string is to be read.
 *   size: <-- maximum number of characters to be read plus one.
 *   f:    <-- ecs_file_t descriptor.
 *   line: <-> file line number if available, or NULL.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 */

char *
ecs_file_gets_try(char              *s,
                  const int          size,
                  const ecs_file_t  *f,
                  int               *line);

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
 *   f:    <-- ecs_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran
 *   record, if the whole record could not be read, returns 0.
 */

size_t
ecs_file_read(void              *rec,
              const size_t       size,
              const size_t       ni,
              const ecs_file_t  *f);

/*
 * Read a binary C or Fortran type record.
 *
 * This function is similar to ecs_file_read(), but failure to read
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
 *   f:    <-- ecs_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran
 *   record, if the whole record could not be read, returns 0.
 */

size_t
ecs_file_read_try(void              *rec,
                  const size_t       size,
                  const size_t       ni,
                  const ecs_file_t  *f);

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
 *   f:    <-- ecs_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully written.
 */

size_t
ecs_file_write(const void        *rec,
               const size_t       size,
               const size_t       ni,
               const ecs_file_t  *f);

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
ecs_file_swap_endian(void          *dest,
                     const void    *src,
                     const size_t   size,
                     const size_t   ni);

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
 *   upon ECS configuration, 1 is returned.
 */

int
ecs_file_mkdir_default(const char  *pathname);

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
ecs_file_isreg(const char  *name);

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
ecs_file_isdir(const char  *name);

/*
 * Indicate Zlib version available at run time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and run time in case we use dynamic libraries.
 *
 * returns:
 *   pointer to string indicating Zlib version in use, or NULL
 *   if Zlib support is not available.
 */

const char *
ecs_file_version_zlib(void);

/*
 * Indicate Zlib version available at compilation time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and link time in case we use dynamic libraries.
 *
 * returns:
 *   pointer to string indicating Zlib version at compilation, or NULL
 *   if Zlib support is not available.
 */

const char *
ecs_file_version_build_zlib(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __ECS_FILE_H__ */
