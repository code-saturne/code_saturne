#ifndef __CS_FILE_H__
#define __CS_FILE_H__

/*============================================================================
 * File and directory operations, with parallel file I/O
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* File descriptor */

typedef struct _cs_file_t  cs_file_t;

/* Helper structure for IO serialization */

#if defined(HAVE_MPI)
typedef struct _cs_file_serializer_t cs_file_serializer_t;
#endif

/* File modes */

typedef enum {

  CS_FILE_MODE_READ,   /* Read mode */
  CS_FILE_MODE_WRITE,  /* Write mode */
  CS_FILE_MODE_APPEND  /* Append mode */

} cs_file_mode_t;

/* Possibilities for the third argument of cs_file_seek() */

typedef enum {

  CS_FILE_SEEK_SET,   /* Seek from beginning of file */
  CS_FILE_SEEK_CUR,   /* Seek from current position */
  CS_FILE_SEEK_END    /* Seek from end of file */

} cs_file_seek_t;

/* File access methods */

typedef enum {

  CS_FILE_DEFAULT,
  CS_FILE_STDIO_SERIAL,
  CS_FILE_STDIO_PARALLEL,
  CS_FILE_MPI_INDEPENDENT,
  CS_FILE_MPI_NON_COLLECTIVE,
  CS_FILE_MPI_COLLECTIVE

} cs_file_access_t;

/* MPI-IO file positioning methods */

typedef enum {

  CS_FILE_MPI_EXPLICIT_OFFSETS,
  CS_FILE_MPI_INDIVIDUAL_POINTERS

} cs_file_mpi_positioning_t;

/* Offset for file position indicator (int64_t in C99) */

#if defined(SIZEOF_LONG_LONG)
typedef long long cs_file_off_t;
#else
typedef long cs_file_off_t;
#endif

/*=============================================================================
 * Global variables
 *============================================================================*/

/* names associated with file access methods */

extern const char  *cs_file_access_name[];

/* names associated with MPI-IO positioning */

extern const char  *cs_file_mpi_positioning_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a file descriptor and open the associated file.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by cs_file_set_swap_endian().
 *
 * parameters:
 *   name       <-- file name
 *   mode       <-- file access mode: read, write, or append
 *   method     <-- file access method
 *   hints      <-- associated hints for MPI-IO, or MPI_INFO_NULL
 *   block_comm <-- handle to MPI communicator used for distributed file
 *                  block access (may be a subset of comm if some ranks do
 *                  not directly access distributed data blocks)
 *   comm       <-- handle to main MPI communicator
 *
 * returns:
 *   pointer to cs_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

cs_file_t *
cs_file_open(const char        *name,
             cs_file_mode_t     mode,
             cs_file_access_t   method,
             MPI_Info           hints,
             MPI_Comm           block_comm,
             MPI_Comm           comm);

#else

cs_file_t *
cs_file_open(const char        *name,
             cs_file_mode_t     mode,
             cs_file_access_t   method);

#endif

/*----------------------------------------------------------------------------
 * Create a file descriptor and open the associated file, using the default
 * file communicator and access method.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by cs_file_set_swap_endian().
 *
 * parameters:
 *   name    <-- file name
 *   mode    <-- file access mode: read, write, or append
 *
 * returns:
 *   pointer to cs_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 *----------------------------------------------------------------------------*/

cs_file_t *
cs_file_open_default(const char        *name,
                     cs_file_mode_t     mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a file descriptor and open the associated file, using the
 *        serial IO on the root rank.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by cs_file_set_swap_endian().
 *
 * \param[in]  name   file name
 * \param[in]  mode   file access mode: read, write, or append
 *
 * \return pointer to cs_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

cs_file_t *
cs_file_open_serial(const char      *name,
                    cs_file_mode_t   mode);

/*----------------------------------------------------------------------------
 * Destroy a file descriptor and close the associated file.
 *
 * parameters:
 *   f <-> file descriptor to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_file_t *
cs_file_free(cs_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's name.
 *
 * parameters:
 *   f <-- cs_file_t descriptor
 *
 * returns:
 *   pointer to the file's name.
 *----------------------------------------------------------------------------*/

const char *
cs_file_get_name(const cs_file_t  *f);

/*----------------------------------------------------------------------------
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * parameters:
 *   f <-> cs_file_t descriptor
 *----------------------------------------------------------------------------*/

void
cs_file_set_big_endian(cs_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's byte-swapping behavior.
 *
 * parameters:
 *   f <-- cs_file_t descriptor
 *
 * returns:
 *   0 if file's endianness is the same as the system's, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_file_get_swap_endian(const cs_file_t  *f);

/*----------------------------------------------------------------------------
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f    <-> cs_file_t descriptor
 *   swap <-- 1 if bytes must be swapped, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_file_set_swap_endian(cs_file_t  *f,
                        int         swap);

/*----------------------------------------------------------------------------
 * Read global data from a file, distributing it to all processes
 * associated with that file.
 *
 * parameters:
 *   f    <-- cs_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
cs_file_read_global(cs_file_t  *f,
                    void       *buf,
                    size_t      size,
                    size_t      ni);

/*----------------------------------------------------------------------------
 * Write global data to a file.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The buffers on other ranks are ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   f    <-- cs_file_t descriptor
 *   buf  <-- pointer to location containing data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
cs_file_write_global(cs_file_t   *f,
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
 *   f                <-- cs_file_t descriptor
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
cs_file_read_block(cs_file_t  *f,
                   void       *buf,
                   size_t      size,
                   size_t      stride,
                   cs_gnum_t   global_num_start,
                   cs_gnum_t   global_num_end);

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
 * cs_file_write_block_buffer() instead may be used to avoid an unneeded
 * memory allocation and copy.
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
 *   buf              <-- pointer to location containing data
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
cs_file_write_block(cs_file_t   *f,
                    const void  *buf,
                    size_t       size,
                    size_t       stride,
                    cs_gnum_t    global_num_start,
                    cs_gnum_t    global_num_end);

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
 * This function is intended to be used mainly data that is already a
 * copy of original data (such as data that has been redistributed across
 * processors just for the sake of output), or that is to be deleted after
 * writing, so it may modify the values in its input buffer (notably to
 * convert from little-endian to big-endian of vice-versa if necessary).
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
 *   buf              <-> pointer to location containing data
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
cs_file_write_block_buffer(cs_file_t  *f,
                           void       *buf,
                           size_t      size,
                           size_t      stride,
                           cs_gnum_t   global_num_start,
                           cs_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Update the file pointer according to whence.
 *
 * parameters:
 *   f      <-> cs_file_t descriptor.
 *   offset <-- add to position specified to whence to obtain new position,
 *              measured in characters from the beginning of the file.
 *   whence <-- beginning if CS_FILE_SEEK_SET, current if CS_FILE_SEEK_CUR,
 *               or end-of-file if CS_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise; currently, errors are fatal.
 *----------------------------------------------------------------------------*/

int
cs_file_seek(cs_file_t       *f,
             cs_file_off_t    offset,
             cs_file_seek_t   whence);

/*----------------------------------------------------------------------------
 * Return the position of the file pointer.
 *
 * In parallel, we consider the file pointer to be equal to the highest
 * value of the individual file pointers.
 *
 * parameters:
 *   f <-- cs_file_t descriptor
 *
 * returns:
 *   current position of the file pointer
 *----------------------------------------------------------------------------*/

cs_file_off_t
cs_file_tell(cs_file_t  *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Formatted input from a text file (as fgets()).
 *
 * \param [out]      s     buffer to which string is to be read.
 * \param [in]       size  maximum number of characters to be read plus one.
 * \param [in]       f     ecs_file_t descriptor.
 * \param [in, out]  line  file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */
/*----------------------------------------------------------------------------*/

char *
cs_file_gets(char             *s,
             const int         size,
             const cs_file_t  *f,
             int              *line);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Formatted input from a text file if possible (as fgets()).
 *
 * This function is similar to cs_file_gets(), but failure to read
 * a line due to an end-of-file condition is not considered an error with
 * this variant, which may be used to read text files or sections thereof
 * of unknown length.
 *
 * \param [out]      s     buffer to which string is to be read.
 * \param [in]       size  maximum number of characters to be read plus one.
 * \param [in]       f     cs_file_t descriptor.
 * \param [in, out]  line  file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */
/*----------------------------------------------------------------------------*/

char *
cs_file_gets_try(char              *s,
                 const int          size,
                 const cs_file_t   *f,
                 int               *line);

/*----------------------------------------------------------------------------
 * Dump the metadata of a file structure in human readable form
 *
 * parameters:
 *   f <-- pointer to file
 *----------------------------------------------------------------------------*/

void
cs_file_dump(const cs_file_t  *f);

/*----------------------------------------------------------------------------
 * Free the default options for file access.
 *----------------------------------------------------------------------------*/

void
cs_file_free_defaults(void);

/*----------------------------------------------------------------------------
 * Get the default options for file access.
 *
 * parameters:
 *   mode   <-- file mode for which the default is queried (write and
 *              append use the same method, and are interchangeable here)
 *   access --> default file access method, or NULL
 *   hints  --> MPI-IO hints, or NULL
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_file_get_default_access(cs_file_mode_t     mode,
                           cs_file_access_t  *method,
                           MPI_Info          *hints);

#else

void
cs_file_get_default_access(cs_file_mode_t     mode,
                           cs_file_access_t  *method);

#endif

/*----------------------------------------------------------------------------
 * Set the default options for file access.
 *
 * If the method given contains incompatible values, such as when setting
 * MPI-IO methods when MPI-IO is not available, a "reasonable" default
 * is used instead.
 *
 * parameters:
 *   mode      <-- file mode for which the default is to be set (write and
 *                 append use the same method, and are interchangeable here)
 *   method    <-- default access method to set
 *   hints     <-- MPI-IO hints, or MPI_INFO_NULL
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_file_set_default_access(cs_file_mode_t    mode,
                           cs_file_access_t  method,
                           MPI_Info          hints);

#else

void
cs_file_set_default_access(cs_file_mode_t    mode,
                           cs_file_access_t  method);

#endif

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get default MPI communicator values for file access.
 *
 * A block rank stepping value may be used, allowing the use of a reduced
 * communicator for distributed block reads and writes.
 * If this value is greater than 1, ranks not a multiple of this step must be
 * guaranteed to be empty for block reads and writes with files opened using
 * this default.
 *
 * parameters:
 *   block_rank_step --> MPI rank stepping between non-empty distributed blocks,
 *                       or NULL
 *   block_comm      --> Handle to MPI communicator used for distributed
 *                       file block access, or NULL
 *   comm            --> Handle to main MPI communicator, or NULL
 *----------------------------------------------------------------------------*/

void
cs_file_get_default_comm(int       *block_rank_step,
                         MPI_Comm  *block_comm,
                         MPI_Comm  *comm);

/*----------------------------------------------------------------------------
 * Set default MPI communicator values for file access.
 *
 * A block rank stepping value may be used, allowing the use of a reduced
 * communicator for distributed block reads and writes.
 * If this value is greater than 1, ranks not a multiple of this step must be
 * guaranteed to be empty for block reads and writes with files opened using
 * this default.
 *
 * For each argument, an "out of range" value may be used to avoid modifying
 * the previous default for that argument.
 *
 * parameters:
 *   block_rank_step <-- MPI rank stepping between non-empty blocks for
 *                       file block reads and writes (not set if <= 0)
 *   comm            <-- handle to main MPI communicator
 *                       (not set if MPI_COMM_SELF)
 *----------------------------------------------------------------------------*/

void
cs_file_set_default_comm(int       block_rank_step,
                         MPI_Comm  comm);

/*----------------------------------------------------------------------------
 * Create an MPI communicator for distributed block parallel IO.
 *
 * parameters:
 *   block_rank_step <-- MPI rank stepping between non-empty blocks
 *   comm            <-- Handle to main MPI communicator
 *
 * returns:
 *   communicator associated with IO, MPI_COMM_NULL for ranks not
 *   participating in parallel IO (including ranks participating in IO
 *   where communicator size would be 1)
 *----------------------------------------------------------------------------*/

MPI_Comm
cs_file_block_comm(int       block_rank_step,
                   MPI_Comm  comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Get the positioning method for MPI-IO
 *
 * For details, see cs_file_set_mpi_io_positioning().
 *
 * returns:
 *   positioning method for MPI-IO
 *----------------------------------------------------------------------------*/

cs_file_mpi_positioning_t
cs_file_get_mpi_io_positioning(void);

/*----------------------------------------------------------------------------
 * Set the positioning method for MPI-IO
 *
 * It is not always known whether a performance or robustness difference is
 * to be expected using explicit file offsets or individual file pointers.
 * Perusal of a sampling of ROMIO code would seem to indicate that no
 * difference is to be expected, but this might change with MPI IO variants
 * or file systems, so this advanced setting is made possible.
 *
 * This setting is not available on a per-file basis, though this could be
 * done in the future in the unexpected case of performance results
 * showing this would be useful.
 *
 * parameters:
 *   positioning <-- chosen positioning method for MPI-IO
 *----------------------------------------------------------------------------*/

void
cs_file_set_mpi_io_positioning(cs_file_mpi_positioning_t  positioning);

/*----------------------------------------------------------------------------
 * Print information on default options for file access.
 *----------------------------------------------------------------------------*/

void
cs_file_defaults_info(void);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a cs_file_serializer_t structure.
 *
 * The buf_block_size argument is optional, and may be used when the buffer
 * on rank 0 is larger than (global_num_end - global_num_start)*size*stride
 * bytes. If zero, a block size of (global_num_end - global_num_start) on
 * rank 0 is assumed; a buffer may not be smaller than this, as it must
 * initially contain all data on rank 0's block.
 *
 * parameters:
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   buf_block_size   <-- Local data buffer block size, or 0 for default
 *                        global_num_end - global_num_start
 *                        (only useful on rank 0)
 *   buf              <-- pointer to local block data buffer
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   pointer to new serializer structure
 *----------------------------------------------------------------------------*/

cs_file_serializer_t *
cs_file_serializer_create(size_t      size,
                          size_t      stride,
                          cs_gnum_t   global_num_start,
                          cs_gnum_t   global_num_end,
                          size_t      buf_block_size,
                          void       *buf,
                          MPI_Comm    comm);

/*----------------------------------------------------------------------------
 * Destroy a cs_file_serializer_t structure.
 *
 * parameters:
 *   s <-> pointer to pointer structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
cs_file_serializer_destroy(cs_file_serializer_t  **s);

/*----------------------------------------------------------------------------
 * Advance a cs_file_serializer_t structure.
 *
 * Data from the buffer of the next communicating rank is copied
 * to rank 0 (this is a no-op the first time this function is called,
 * as rank 0 already has its data).
 *
 * On rank 0, the return value may point to the buffer defined when
 * initializing the serializer, or to an aditional buffer if the former is
 * too small to receive data from all ranks.
 *
 * Note also that for ranks > 0, this function always returns NULL,
 * as only one call is needed for those ranks.
 *
 * parameters:
 *   s         <-- pointer to serializer structure
 *   cur_range --> optional start and past-the end global numbers for the
 *                 current block (size: 2), or NULL; only on rank 0
 *
 * returns:
 *   a pointer to the buffer containing new data (first call counts as new),
 *   or NULL if we are finished; always NULL on ranks > 0
 *----------------------------------------------------------------------------*/

void *
cs_file_serializer_advance(cs_file_serializer_t  *s,
                           cs_gnum_t              cur_range[2]);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
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
 *   path: <-- name of new directory.
 *
 * returns:
 *   0 on success, -1 if an error occured (in which case errno
 *   contains the appropriate error code). If the underlying
 *   system has no mkdir() function or it was not detected
 *   upon BFT configuration, 1 is returned.
 *----------------------------------------------------------------------------*/

int
cs_file_mkdir_default(const char  *path);

/*----------------------------------------------------------------------------
 * Check if a file exists and is a regular file.
 *
 * parameters:
 *   path <-- file name.
 *
 * returns:
 *   1 if file exists and is a regular file, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_file_isreg(const char  *path);

/*----------------------------------------------------------------------------
 * Check if a directory exists.
 *
 * parameters:
 *   path <-- directory name.
 *
 * returns:
 *   1 if directory exists, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_file_isdir(const char  *path);

/*----------------------------------------------------------------------------
 * List files inside a directory.
 *
 * The array returned must be freed by the caller using BFT_FREE,
 * as well as the individual entries in the array.
 *
 * parameters:
 *   path <-- name of directory.
 *
 * returns:
 *   an array of file names in a directory. The last entry is set to NULL.
 *   If no means to list the directory or an error occured, the return
 *    value is simply NULL.
 *----------------------------------------------------------------------------*/

char **
cs_file_listdir(const char *path);

/*----------------------------------------------------------------------------
 * Return the size of a file.
 *
 * If the file does not exist, 0 is returned.
 *
 * Note that for some special files, such as files in the Linux /proc
 * directory, this may return 0.
 *
 * parameters
 *   path <-- file path.
 *
 * returns:
 *   size of file.
 *----------------------------------------------------------------------------*/

cs_file_off_t
cs_file_size(const char  *path);

/*----------------------------------------------------------------------------
 * Remove a file if it exists and is a regular file or an empty directory.
 *
 * parameters
 *   path <-- file path.
 *
 * returns:
 *   0 in case of success or if file does not exist, not 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_file_remove(const char  *path);

/*----------------------------------------------------------------------------
 * Check if a file name ends with a specific string (extension)
 *
 * parameters
 *   path <-- file path.
 *   end  <-- string to compare
 *
 * returns:
 *   1 if the path ends with the given string, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_file_endswith(const char  *path,
                 const char  *end);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FILE_H__ */
