/*============================================================================
 * Parallel file I/O
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2007-2011  EDF

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) \
                              && defined(HAVE_UNISTD_H)
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

/*
  Force LARGEFILE_SOURCE if largefiles enabled under 32-bit Linux or Blue Gene
  (otherwise, we may encounter bugs with glibc 2.3 due to fseeko end ftello
  not being correctly defined). Compiling with -D_GNU_SOURCE instead
  of -D_POSIX_C_SOURCE=200112L seems to be another way to solve the problem.
*/

#if (SIZEOF_LONG < 8) && (_FILE_OFFSET_BITS == 64)
# if defined(__linux__) || defined(__blrts__) || defined(__bgp__)
#  if !defined(_POSIX_SOURCE)
#    define _GNU_SOURCE 1
#  endif
#  if !defined(_GNU_SOURCE) && !defined(_LARGEFILE_SOURCE)
#   define _LARGEFILE_SOURCE 1
#  endif
# endif
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_file.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/* FVM file descriptor */

struct _fvm_file_t {

  char              *name;         /* File name */
  fvm_file_mode_t    mode;         /* File mode */
  int                semantics;    /* Preferred file positioning semantics */
  int                rank;         /* MPI rank */
  int                n_ranks;      /* MPI rank */
  _Bool              swap_endian;  /* Swap big-endian and little-endian ? */

  FILE              *sh;           /* Serial file handle */

#if defined(HAVE_MPI)
  MPI_Comm           comm;         /* Associated MPI communicator */
#if defined(HAVE_MPI_IO)
  MPI_File           fh;           /* MPI file handle */
  MPI_Info           info;         /* MPI file info */
  MPI_Offset         offset;       /* MPI file offset */
#endif
#endif

};

#if defined(HAVE_MPI)

/* Helper structure for IO serialization */

struct _fvm_file_serializer_t {

  int          rank_id;        /* Local rank in communicator */
  int          n_ranks;        /* Number of ranks in communicator */

  fvm_gnum_t   range[2];       /* Global start and past-the-end numbers
                                  for local rank */

  size_t       size;           /* datatype size (may include stride) */

  fvm_gnum_t   next_g_num;     /* Next global number */
  int          next_rank_id;   /* Next rank with which we will communicate */

  fvm_lnum_t  *count;          /* Number of elements in each block */

  void        *buf;            /* pointer to external buffer */
  void        *recv_buf;       /* pointer to external buffer if
                                  buf_block_size >= max_block_size,
                                  or to buf otherwise */

  MPI_Comm     comm;           /* Associated MPI communicator */
};

#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default hints */

#if defined(HAVE_MPI_IO)
static fvm_file_hints_t _default_semantics = FVM_FILE_INDIVIDUAL_POINTERS;
#else
static fvm_file_hints_t _default_semantics = 0;
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize an fvm_file_serializer_t structure.
 *
 * The buf_block_size argument is optional, and may be used when the buffer
 * on rank 0 is larger than (global_num_end - global_num_start)*size*stride
 * bytes. If zero, a block size of (global_num_end - global_num_start) on
 * rank 0 is assumed; a buffer may not be smaller than this, as it must
 * initially contain all data on rank 0's block.
 *
 * parameters:
 *   s                <-> pointer to structure that should be initialized
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   buf_block_size   <-- Local data buffer block size, or 0 for default
 *                        global_num_end - global_num_start
 *                        (only useful on rank 0)
 *   buf              <-- pointer to local block data buffer
 *   comm             <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_serializer_init(fvm_file_serializer_t  *s,
                 size_t                  size,
                 fvm_gnum_t              global_num_start,
                 fvm_gnum_t              global_num_end,
                 size_t                  buf_block_size,
                 void                   *buf,
                 MPI_Comm                comm)

{
  fvm_lnum_t l_count = 0;

  s->range[0] = global_num_start;
  s->range[1] = global_num_end;

  s->size = size;

  if (s->range[1] > s->range[0])
    l_count = s->range[1] - s->range[0];

  /* Get local rank and size of the current MPI communicator */

  MPI_Comm_rank(comm, &(s->rank_id));
  MPI_Comm_size(comm, &(s->n_ranks));

  s->next_rank_id = 0;
  s->next_g_num = global_num_start;

  /* Initialize counter */

  if (s->rank_id == 0)
    BFT_MALLOC(s->count, s->n_ranks, fvm_lnum_t);
  else
    s->count = NULL;

  MPI_Gather(&l_count, 1, FVM_MPI_LNUM, s->count, 1, FVM_MPI_LNUM, 0, comm);

  /* Allocate local buffer if necessary, or point to external buffer */

  s->buf = buf;
  s->recv_buf = NULL;

  if (s->rank_id == 0) {
    int i;
    fvm_lnum_t _max_block_size = 0;
    fvm_lnum_t _buf_block_size = FVM_MAX((fvm_lnum_t)buf_block_size, l_count);
    for (i = 0; i < s->n_ranks; i++)
      _max_block_size = FVM_MAX(_max_block_size, s->count[i]);
    if (_max_block_size > _buf_block_size)
      BFT_MALLOC(s->recv_buf, _max_block_size*size, unsigned char);
    else
      s->recv_buf = buf;
  }

  s->comm = comm;
}

/*----------------------------------------------------------------------------
 * Finalize an fvm_file_serializer_t structure.
 *
 * parameters:
 *   s <-- pointer to structure that should be finalized
 *----------------------------------------------------------------------------*/

static void
_serializer_finalize(fvm_file_serializer_t  *s)
{
  s->next_rank_id = 0;
  s->next_g_num = 1;

  if (s->count != NULL)
    BFT_FREE(s->count);

  if (s->recv_buf != s->buf && s->recv_buf != NULL)
    BFT_FREE(s->recv_buf);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * parameters:
 *   dest <-- pointer to converted data location.
 *   src  --> pointer to source data location.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of data items.
 *----------------------------------------------------------------------------*/

static void
_swap_endian(void        *dest,
             const void  *src,
             size_t       size,
             size_t       ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0; i < ni; i++) {

    shift = i * size;

    for (ib = 0; ib < (size / 2); ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

/*----------------------------------------------------------------------------
 * Open a file using standard C IO.
 *
 * parameters:
 *   f    <-- pointer to file handler
 *   mode <-- file acces mode: read, write, or append
 *
 * returns:
 *   0 in case of success, error number in case of failure
 *----------------------------------------------------------------------------*/

static int
_file_open(fvm_file_t       *f,
           fvm_file_mode_t   mode)
{
  int retval = 0;

  assert(f != NULL);

  if (f->sh != NULL)
    return 0;

  /* The file handler exists and the corresponding file is closed */

  f->mode = mode;

  switch (f->mode) {
  case FVM_FILE_MODE_APPEND:
    f->sh = fopen(f->name, "a");
    break;
  case FVM_FILE_MODE_WRITE:
    f->sh = fopen(f->name, "w");
    break;
  default:
    assert(f->mode == FVM_FILE_MODE_READ);
    f->sh = fopen(f->name, "r");
  }

  if (f->sh == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error opening file \"%s\":\n\n"
                "  %s"), f->name, strerror(errno));
    retval = errno;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Close a file using standard C IO.
 *
 * parameters:
 *   f <-> pointer to file handler
 *
 * returns:
 *   0 in case of success, -1 in case of failure
 *----------------------------------------------------------------------------*/

static int
_file_close(fvm_file_t  *f)
{
  int retval = 0;

  if (f->sh != NULL)
    retval = fclose(f->sh);

  if (retval != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error closing file \"%s\":\n\n"
                "  %s"), f->name, strerror(errno));
    retval = errno;
  }
  f->sh = NULL;

  return retval;
}

/*----------------------------------------------------------------------------
 * Read data to a buffer using standard C IO.
 *
 * parameters:
 *   f    <-- fvm_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_read(fvm_file_t  *f,
           void        *buf,
           size_t       size,
           size_t       ni)
{
  size_t retval = 0;

  assert(f->sh != NULL);

  if (ni != 0)
    retval = fread(buf, size, ni, f->sh);

  /* In case of error, determine error type */

  if (retval != ni) {
    int err_num = ferror(f->sh);
    if (err_num != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error reading file \"%s\":\n\n  %s"),
                f->name, strerror(err_num));
    else if (feof(f->sh) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Premature end of file \"%s\""), f->name);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error reading file \"%s\""), f->name);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Write data to a file using standard C IO.
 *
 * parameters:
 *   f    <-- fvm_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_write(fvm_file_t  *f,
            const void  *buf,
            size_t       size,
            size_t       ni)
{
  size_t retval = 0;

  assert(f->sh != NULL);

  if (ni != 0)
    retval = fwrite(buf, size, ni, f->sh);

  /* In case of error, determine error type */

  if (retval != ni) {
    int err_num = ferror(f->sh);
    if (err_num != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error writing file \"%s\":\n\n  %s"),
                f->name, strerror(err_num));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error writing file \"%s\""), f->name);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Sets a file's position indicator using standard C IO.
 *
 * This function may call the libc's fseek() or fseeko() function.
 * The C 99 standard specifies that for a text file, the offset
 * argument to fseek() should be zero or a value returned by an earlier
 * successful call to ftell().
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * parameters:
 *   f      <-> file descriptor.
 *   offset <-- add to position specified to whence to obtain new
 *              position, measured in characters from the beginning of
 *              the file.
 *   whence <-- beginning if FVM_FILE_SEEK_SET, current if
 *              FVM_FILE_SEEK_CUR, or end-of-file if FVM_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise.
 *----------------------------------------------------------------------------*/

static int
_file_seek(fvm_file_t       *f,
           fvm_file_off_t    offset,
           fvm_file_seek_t   whence)
{
  static int _stdio_seek[3] = {SEEK_SET, SEEK_CUR, SEEK_END};

  int _whence = _stdio_seek[whence];
  int retval = 0;

  /* Convert fvm_file_seek to stdio values */

  assert(f != NULL);

  if (f->sh != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support may be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)

    retval = fseeko(f->sh, (off_t)offset, _whence);

    if (retval != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error setting position in file \"%s\":\n\n  %s"),
                f->name, strerror(errno));
# else

    /* Test if offset larger than allowed */

    long _offset = offset;

    if (_offset == offset) {
      retval = fseek(f->sh, (long)offset, _whence);
      if (retval != 0)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error setting position in file \"%s\":\n\n  %s"),
                  f->name, strerror(errno));
    }
    else {
      retval = -1;
      bft_error
        (__FILE__, __LINE__, 0,
         _("Error setting position in file \"%s\":\n\n  %s"),
         f->name,
         _("sizeof(off_t) > sizeof(long) but fseeko() not available"));
    }

# endif /* defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64) */

#else /* SIZEOF_LONG >= 8 */

    /* For 64-bit systems, standard fseek should be enough */

    retval = fseek(f->sh, (long)offset, _whence);
    if (retval != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error setting position in file \"%s\":\n\n  %s"),
                f->name, strerror(errno));

#endif /* SIZEOF_LONG */
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Obtain the current value of a file's position indicator.
 *
 * parameters:
 *   f  <-- file descriptor.
 *
 * returns:
 *   current value of the file's position indicator, or -1 in case of failure.
 *----------------------------------------------------------------------------*/

static fvm_file_off_t
_file_tell(fvm_file_t  *f)
{
  fvm_file_off_t offset = 0;

  assert(f != NULL);

  if (f->sh != NULL) {

    /* For 32-bit systems, large file support may be necessary */

#if (SIZEOF_LONG < 8)

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    offset = ftello(f->sh);
# else
    /*
      Without ftello, ftell will fail above 2 Gigabytes, in which case
      offset == -1 and errno == EOVERFLOW, but should work on smaller
      files. We prefer not to be too strict about fseeko availability, as
      the only 32-bit case without ftello we have encountered is Cygwin
      (for which ftello requires additional non-default libraries), which
      is expected to be used mainly for small cases.
    */
    offset = ftell(f->sh);
# endif

    /* For 64-bit systems, standard ftell should be enough */

#else /* SIZEOF_LONG >= 8 */
    offset = ftell(f->sh);
#endif

  }

  if (offset < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error obtaining position in file \"%s\":\n\n  %s"),
              f->name, strerror(errno));

  return offset;
}

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
 * This version does not use MPI-IO
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_read_block(fvm_file_t  *f,
                 void        *buf,
                 size_t       size,
                 fvm_gnum_t   global_num_start,
                 fvm_gnum_t   global_num_end)
{
  size_t retval = 0;

  if (f->rank == 0)
    retval = _file_read(f,
                        buf,
                        size,
                        (size_t)(global_num_end - global_num_start));

#if defined(HAVE_MPI)

  if (f->comm != MPI_COMM_NULL) {

    MPI_Status status;

    int loc_count = global_num_end - global_num_start;
    int _counts[64];
    int *counts = NULL;

    if (f->rank == 0) {
      if (f->n_ranks < 64)
        counts = _counts;
      else
        BFT_MALLOC(counts, f->n_ranks, int);
    }

    /* Exchange counts */

    MPI_Gather(&loc_count, 1, MPI_INT, counts, 1, MPI_INT, 0, f->comm);

    /* Rank 0 reads data for other ranks from file and distributes it */

    if (f->rank == 0) {

      int dist_rank;
      int _buf_size = global_num_end - global_num_start;
      unsigned char *_buf = NULL;

      /* Allocate exchange buffer */

      for (dist_rank = 1; dist_rank < f->n_ranks; dist_rank++)
        _buf_size = FVM_MAX(_buf_size, counts[dist_rank]);

      BFT_MALLOC(_buf, _buf_size*size, unsigned char);

      /* Loop on distant ranks */

      for (dist_rank = 1; dist_rank < f->n_ranks; dist_rank++) {

        if (counts[dist_rank] == 0)
          continue;

        /* Read data from file */

        counts[dist_rank]
          = (int)_file_read(f, _buf, size, (size_t)counts[dist_rank]);

        /* Send to corresponding rank */

        MPI_Send(_buf, counts[dist_rank]*size, MPI_BYTE, dist_rank,
                 FVM_MPI_TAG, f->comm);

      } /* End of loop on distant ranks */

      BFT_FREE(_buf);

    }

    /* Other ranks receive data from rank 0 */

    else if (loc_count > 0) {

      /* Receive data */

      MPI_Recv(buf, (int)(loc_count*size), MPI_BYTE, 0,
               FVM_MPI_TAG, f->comm, &status);

      MPI_Get_count(&status, MPI_BYTE, &loc_count);
      retval = loc_count / size;

    }

    if (counts != NULL && counts != _counts)
      BFT_FREE(counts);
  }

#endif /* defined(HAVE_MPI) */

  return retval;
}

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
 * This version does not use MPI-IO
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              <-> pointer to location containing data
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully written;
 *----------------------------------------------------------------------------*/

static size_t
_file_write_block(fvm_file_t  *f,
                  void        *buf,
                  size_t       size,
                  fvm_gnum_t   global_num_start,
                  fvm_gnum_t   global_num_end)
{
  size_t retval = 0;

  if (f->n_ranks == 1)
    retval = _file_write(f,
                         buf,
                         size,
                         (size_t)(global_num_end - global_num_start));

#if defined(HAVE_MPI)

  if (f->n_ranks > 1) {

    fvm_file_serializer_t  s;
    fvm_lnum_t *count = NULL;
    void  *write_buf = NULL;

    _serializer_init(&s,
                     size,
                     global_num_start,
                     global_num_end,
                     0,
                     buf,
                     f->comm);

    do {

      int dist_rank = s.next_rank_id;

      write_buf = fvm_file_serializer_advance(&s, NULL);

      if (write_buf != NULL) /* only on rank 0 */
        s.count[dist_rank]
          = (fvm_lnum_t)_file_write(f,
                                    write_buf,
                                    size,
                                    (size_t)(s.count[dist_rank]));

    } while (write_buf != NULL);

    /* Exchange return codes */

    if (s.rank_id == 0)
      count = s.count;
    else
      BFT_MALLOC(count, s.n_ranks, fvm_lnum_t);

    MPI_Scatter(count, 1, MPI_INT, &retval, 1, MPI_INT,  0, f->comm);

    if (s.rank_id != 0)
      BFT_FREE(count);

    _serializer_finalize(&s);
  }

#endif /* defined(HAVE_MPI) */

  return retval;
}

#if defined(HAVE_MPI_IO)

/*----------------------------------------------------------------------------
 * Check for and remove existing file.
 *
 * This is necessary because some MPI-IO implementations seem to append
 * data to an exisiting file even when MPI_MODE_CREATE is specified.
 *
 * parameters:
 *   f <-- fvm_file_t descriptor
 *----------------------------------------------------------------------------*/

static void
_file_clear(fvm_file_t  *f)
{
  int exists = 0;

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) \
                              && defined(HAVE_UNISTD_H)

  struct stat s;

  if (stat(f->name, &s) == 0) {
    if (S_ISREG(s.st_mode) != 0)
      exists = 1;
  }
  if (exists)
    unlink(f->name);

#else

  /* If Posix-type API is not available, revert to basic method */

  FILE *f;

  if ((f = fopen(fic_name, "w")) != NULL)
    fclose(f);

#endif
}

/*----------------------------------------------------------------------------
 * Output MPI error message.
 *
 * This supposes that the default MPI errorhandler is not used
 *
 * parameters:
 *   file_name  <-- file name
 *   error_code <-- associated MPI error code
 *
 * returns:
 *   0 in case of success, system error code in case of failure
 *----------------------------------------------------------------------------*/

static void
_mpi_io_error_message
(
 const char  *file_name,
 int          error_code
)
{
  char buffer[MPI_MAX_ERROR_STRING];
  int  buffer_len;

  MPI_Error_string(error_code, buffer, &buffer_len);

  bft_error(__FILE__, __LINE__, 0,
            _("MPI IO error for file: %s\n"
              "Error type: %s"), file_name, buffer);
}

/*----------------------------------------------------------------------------
 * Open a file using MPI IO.
 *
 * parameters:
 *   f     <-- pointer to file handler
 *   mode  <-- file acces mode: read, write, or append
 *
 * returns:
 *   MPI_SUCCESS in case of success, MPI error code in case of failure
 *----------------------------------------------------------------------------*/

static int
_mpi_file_open(fvm_file_t       *f,
               fvm_file_mode_t   mode)
{
  int amode = MPI_MODE_RDWR;
  MPI_Info  info = MPI_INFO_NULL;
  int retval = 0;

  assert(f != NULL);

  if (f->fh != MPI_FILE_NULL)
    return 0;

  /* Set access mode */

  f->mode = mode;

  if (f->mode == FVM_FILE_MODE_APPEND)
    amode = MPI_MODE_WRONLY | MPI_MODE_APPEND;

  else if (f->mode == FVM_FILE_MODE_WRITE) {
    int rank;
    amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    MPI_Comm_rank(f->comm, &rank);
    if (rank < 1)
      _file_clear(f);
  }

  else if (f->mode == FVM_FILE_MODE_READ)
    amode = MPI_MODE_RDONLY;

  /* Open file */

  retval = MPI_File_open(f->comm, f->name, amode, info, &(f->fh));

  if (retval != MPI_SUCCESS)
    _mpi_io_error_message(f->name, retval);

  if (f->mode == FVM_FILE_MODE_APPEND) {
    retval = MPI_File_get_position(f->fh, &(f->offset));
    if (retval != MPI_SUCCESS)
      _mpi_io_error_message(f->name, retval);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Close a file using MPI IO.
 *
 * parameters:
 *   f <-> pointer to file handler
 *
 * returns:
 *   MPI_SUCCESS in case of success, MPI error code in case of failure
 *----------------------------------------------------------------------------*/

static int
_mpi_file_close(fvm_file_t  *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->fh == MPI_FILE_NULL)
    return 0;

  /* Close file */

  retval = MPI_File_close(&(f->fh));

  if (retval != MPI_SUCCESS)
    _mpi_io_error_message(f->name, retval);

  return retval;
}

/*----------------------------------------------------------------------------
 * Read data to a buffer, distributing a contiguous part of it to each
 * process associated with a file.
 *
 * Each process should receive a block of the data, and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * There are 2 variants, depending on the semantics:
 *   _mpi_file_read_block_eo (using explicit offsets)
 *   _mpi_file_read_block_ip (using individual pointers, setting a file view)
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_mpi_file_read_block_eo(fvm_file_t  *f,
                        void        *buf,
                        size_t       size,
                        fvm_gnum_t   global_num_start,
                        fvm_gnum_t   global_num_end)
{
  fvm_gnum_t global_num_end_last = global_num_end;
  size_t retval = 0;

  MPI_Offset disp;
  MPI_Status status;
  int errcode, count;

  disp = f->offset + ((global_num_start - 1) * size);
  count = (global_num_end - global_num_start)*size;

  errcode = MPI_File_read_at_all(f->fh, disp, buf, count, MPI_BYTE, &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  if (count > 0)
    MPI_Get_count(&status, MPI_BYTE, &count);
  retval = count / size;

  MPI_Bcast(&global_num_end_last, 1, FVM_MPI_GNUM, f->n_ranks-1, f->comm);
  f->offset += ((global_num_end_last - 1) * size);

  return retval;
}

static size_t
_mpi_file_read_block_ip(fvm_file_t  *f,
                        void        *buf,
                        size_t       size,
                        fvm_gnum_t   global_num_start,
                        fvm_gnum_t   global_num_end)
{
  int errcode;
  int lengths[1];
  MPI_Aint disps[1];
  MPI_Status status;
  MPI_Datatype file_type;

  int count = 0;
  char datarep[] = "native";
  fvm_gnum_t global_num_end_last = global_num_end;
  size_t retval = 0;

  lengths[0] = (global_num_end - global_num_start) * size;
  disps[0] = (global_num_start - 1) * size;

  MPI_Type_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
  MPI_Type_commit(&file_type);

  MPI_File_set_view(f->fh, f->offset, MPI_BYTE, file_type, datarep, f->info);

  errcode = MPI_File_read_all(f->fh, buf, (int)(lengths[0]), MPI_BYTE,
                              &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  MPI_Type_free(&file_type);

  if (lengths[0] > 0)
    MPI_Get_count(&status, MPI_BYTE, &count);
  retval = count / size;

  MPI_Bcast(&global_num_end_last, 1, FVM_MPI_GNUM, f->n_ranks-1, f->comm);
  f->offset += ((global_num_end_last - 1) * size);

  return retval;
}

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
 * There are 2 variants, depending on the semantics:
 *   _mpi_file_write_block_eo (using explicit offsets)
 *   _mpi_file_write_block_ip (using individual pointers, setting a file view)
 *
 * parameters:
 *   f                <-- fvm_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_mpi_file_write_block_eo(fvm_file_t  *f,
                         void        *buf,
                         size_t       size,
                         fvm_gnum_t   global_num_start,
                         fvm_gnum_t   global_num_end)
{
  fvm_gnum_t global_num_end_last = global_num_end;
  size_t retval = 0;

  MPI_Offset disp;
  MPI_Status status;
  int errcode, count;

  disp = f->offset + ((global_num_start - 1) * size);
  count = (global_num_end - global_num_start)*size;

  errcode = MPI_File_write_at_all(f->fh, disp, buf, count, MPI_BYTE, &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  if (count > 0)
    MPI_Get_count(&status, MPI_BYTE, &count);
  retval = count / size;

  MPI_Bcast(&global_num_end_last, 1, FVM_MPI_GNUM, f->n_ranks-1, f->comm);
  f->offset += ((global_num_end_last - 1) * size);

  return retval;
}

static size_t
_mpi_file_write_block_ip(fvm_file_t  *f,
                         void        *buf,
                         size_t       size,
                         fvm_gnum_t   global_num_start,
                         fvm_gnum_t   global_num_end)
{
  int lengths[1];
  MPI_Aint disps[1];
  MPI_Status status;
  MPI_Datatype file_type;

  int errcode = MPI_SUCCESS, count = 0;
  char datarep[] = "native";
  fvm_gnum_t global_num_end_last = global_num_end;
  size_t retval = 0;

  lengths[0] = (global_num_end - global_num_start) * size;
  disps[0] = (global_num_start - 1) * size;

  MPI_Type_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
  MPI_Type_commit(&file_type);

  MPI_File_set_view(f->fh, f->offset, MPI_BYTE, file_type, datarep, f->info);

  errcode = MPI_File_write_all(f->fh, buf, (int)(lengths[0]), MPI_BYTE,
                               &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  MPI_Type_free(&file_type);

  if (lengths[0] > 0)
    MPI_Get_count(&status, MPI_BYTE, &count);
  retval = count / size;

  MPI_Bcast(&global_num_end_last, 1, FVM_MPI_GNUM, f->n_ranks-1, f->comm);
  f->offset += ((global_num_end_last - 1) * size);

  return retval;
}

#endif /* defined(HAVE_MPI_IO) */

/*=============================================================================
 * Public function definitions
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
 *   comm  <-- associated MPI communicator
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
              MPI_Comm            comm)

#else

fvm_file_t *
fvm_file_open(const char         *name,
              fvm_file_mode_t     mode,
              fvm_file_hints_t    hints)

#endif
{
  int errcode = 0;
  fvm_file_t * f = NULL;
  fvm_file_hints_t _hints = _default_semantics;

  BFT_MALLOC(f, 1, fvm_file_t);

  f->sh = NULL;

#if defined(HAVE_MPI)
  f->comm = MPI_COMM_NULL;
#if defined(HAVE_MPI_IO)
  f->fh = MPI_FILE_NULL;
  f->info = MPI_INFO_NULL;
  f->offset = 0;
#endif
#endif

  BFT_MALLOC(f->name, strlen(name) + 1, char);
  strcpy(f->name, name);

  f->mode = mode;
  f->semantics = FVM_FILE_NO_MPI_IO;
  f->rank = 0;
  f->n_ranks = 1;

  f->swap_endian = false; /* Use native endianness by default */

  /* Set communicator */

#if defined(HAVE_MPI)
  {
    if (hints != 0)
      _hints = hints;

    if (comm != MPI_COMM_NULL) {
      MPI_Comm_size(comm, &(f->n_ranks));
      if (f->n_ranks > 1) {
        MPI_Comm_dup(comm, &(f->comm));
        MPI_Comm_rank(f->comm, &(f->rank));
      }
      else
        f->comm = MPI_COMM_NULL;
    }
  }
#endif /* defined(HAVE_MPI) */

  /* Use MPI IO ? */

#if defined(HAVE_MPI_IO)
  if (   f->comm != MPI_COMM_NULL
      && !(_hints & FVM_FILE_NO_MPI_IO)) {
    int positioning_mask = (  FVM_FILE_EXPLICIT_OFFSETS
                            | FVM_FILE_INDIVIDUAL_POINTERS);
    if (_hints & positioning_mask)
      f->semantics = _hints & positioning_mask;
    else
      f->semantics = FVM_FILE_INDIVIDUAL_POINTERS;
    f->semantics = f->semantics | (_hints & FVM_FILE_NO_PREDISTRIBUTE);
  }
#endif

  /* Open file. In case of failure, destroy the allocated structure;
     this is only useful with a non-default error handler,
     as the program is terminated by default */

  if ((f->semantics & FVM_FILE_NO_MPI_IO) && f->rank == 0)
    errcode = _file_open(f, f->mode);

#if defined(HAVE_MPI_IO)
  else if (!(f->semantics & FVM_FILE_NO_MPI_IO))
    errcode = _mpi_file_open(f, f->mode);
#endif

  if (errcode != 0)
    f = fvm_file_free(f);

  return f;
}

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
fvm_file_free(fvm_file_t  *f)
{
  fvm_file_t  *_f = f;
  int errcode = 0;

  if (_f->sh != NULL)
    errcode = _file_close(_f);

#if defined(HAVE_MPI_IO)
  else if (_f->fh != MPI_FILE_NULL)
    errcode = _mpi_file_close(_f);
#endif

#if defined(HAVE_MPI)
  if (_f->comm != MPI_COMM_NULL)
    MPI_Comm_free(&(_f->comm));
#endif

  BFT_FREE(_f->name);
  BFT_FREE(_f);

  return NULL;
}

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
fvm_file_get_name(const fvm_file_t  *f)
{
  assert(f != NULL);

  return f->name;
}

/*----------------------------------------------------------------------------
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * parameters:
 *   f <-> fvm_file_t descriptor
 *----------------------------------------------------------------------------*/

void
fvm_file_set_big_endian(fvm_file_t  *f)
{
  unsigned  int_endian;

  /* Check if system is "big-endian" or "little-endian" */

  int_endian = 0;
  *((char *)(&int_endian)) = '\1';

  if (int_endian == 1)
    f->swap_endian = 1;

#if defined(DEBUG) && !defined(NDEBUG)

  else {
    int_endian = 0;
    *((char *) (&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert(int_endian == 1);
  }

#endif
}

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
fvm_file_get_swap_endian(const fvm_file_t  *f)
{
  assert(f != NULL);

  return f->swap_endian;
}

/*----------------------------------------------------------------------------
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f    <-> fvm_file_t descriptor
 *   swap  --> 1 if bytes must be swapped, 0 otherwise
 *----------------------------------------------------------------------------*/

void
fvm_file_set_swap_endian(fvm_file_t  *f,
                         int          swap)
{
  assert(f != NULL);

  f->swap_endian = swap;
}

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
                     size_t       ni)
{
  size_t retval = 0;

  if ((f->semantics & FVM_FILE_NO_MPI_IO)&& f->rank == 0) {
    retval = _file_read(f, buf, size, ni);
  }

#if defined(HAVE_MPI)
  {
    if ((f->semantics & FVM_FILE_NO_MPI_IO) && f->comm != MPI_COMM_NULL) {
      long _retval = retval;
      MPI_Bcast(buf, size*ni, MPI_BYTE, 0, f->comm);
      MPI_Bcast(&_retval, 1, MPI_LONG, 0, f->comm);
      retval = _retval;
    }

#if defined(HAVE_MPI_IO)

    else if (!(f->semantics & FVM_FILE_NO_MPI_IO)) {

      MPI_Status status;
      int errcode = MPI_SUCCESS, count = 0;

      if (f->semantics & FVM_FILE_EXPLICIT_OFFSETS) {
        errcode = MPI_File_read_at_all(f->fh,
                                       f->offset,
                                       buf,
                                       size*ni,
                                       MPI_BYTE,
                                       &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
      }

      else if (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS) {
        MPI_Datatype file_type;
        MPI_Aint disps[1];
        int lengths[1];
        char datarep[] = "native";
        lengths[0] = ni * size;
        disps[0] = 0;
        MPI_Type_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
        MPI_Type_commit(&file_type);
        MPI_File_set_view(f->fh, f->offset, MPI_BYTE, file_type,
                          datarep, f->info);
        errcode = MPI_File_read_all(f->fh,
                                    buf,
                                    size*ni,
                                    MPI_BYTE,
                                    &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
        MPI_Type_free(&file_type);
      }

      if (errcode != MPI_SUCCESS)
        _mpi_io_error_message(f->name, errcode);

      retval = count / size;

      f->offset += count;
    }

#endif /* defined(HAVE_MPI_IO) */
  }
#endif /* defined(HAVE_MPI) */

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf, buf, size, retval);

  return retval;
}

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
                      size_t       ni)
{
  size_t retval = ni;

  unsigned char _copybuf[1024];
  unsigned char *copybuf = _copybuf;
  const void *_buf = buf;

  /* Copy contents to ensure buffer constedness if necessary */

  if (   f->rank == 0
      && (   (f->swap_endian == true && size > 1)
          || !(f->semantics & FVM_FILE_NO_MPI_IO))) {

    if (size*ni > sizeof(_copybuf))
      BFT_MALLOC(copybuf, size*ni, unsigned char);
    memcpy(copybuf, buf, size*ni);

    if (f->swap_endian == true && size > 1)
      _swap_endian(copybuf, copybuf, size, ni);

    _buf = copybuf;
  }

  if ((f->semantics & FVM_FILE_NO_MPI_IO) && f->sh != NULL) {
    retval = _file_write(f,
                         _buf,
                         size,
                         ni);
  }

#if defined(HAVE_MPI_IO)

  if (f->comm != MPI_COMM_NULL && (!(f->semantics & FVM_FILE_NO_MPI_IO))) {

    MPI_Status status;
    int aux[2] = {MPI_SUCCESS, 0}; /* 0: return value; 1: count */

    if (f->semantics & FVM_FILE_EXPLICIT_OFFSETS) {
      if (f->rank == 0) {
        aux[0] = MPI_File_write_at(f->fh,
                                   f->offset,
                                   copybuf,
                                   size*ni,
                                   MPI_BYTE,
                                   &status);
        MPI_Get_count(&status, MPI_BYTE, &(aux[1]));
      }
    }

    else if (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS) {
      MPI_Datatype file_type;
      MPI_Aint disps[1];
      int lengths[1];
      char datarep[] = "native";
      lengths[0] = ni * size;
      disps[0] = 0;
      MPI_Type_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
      MPI_Type_commit(&file_type);
      MPI_File_set_view(f->fh, f->offset, MPI_BYTE, file_type,
                        datarep, f->info);
      if (f->rank == 0) {
        aux[0] = MPI_File_write(f->fh,
                                copybuf,
                                size*ni,
                                MPI_BYTE,
                                &status);
        MPI_Get_count(&status, MPI_BYTE, &(aux[1]));
      }
      MPI_Type_free(&file_type);
    }

    MPI_Bcast(aux, 2, MPI_INT, 0, f->comm);

    if (aux[0] != MPI_SUCCESS)
      _mpi_io_error_message(f->name, aux[0]);

    retval = aux[1] / size;

    f->offset += aux[1];
  }

#endif /* defined(HAVE_MPI_IO) */

  if (copybuf != _copybuf) /* Free allocated memory if necessary */
    BFT_FREE(copybuf);

  return retval;
}

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
 *   the (local) number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvm_file_read_block(fvm_file_t  *f,
                    void        *buf,
                    size_t       size,
                    size_t       stride,
                    fvm_gnum_t   global_num_start,
                    fvm_gnum_t   global_num_end)
{
  size_t retval = 0;

  const fvm_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
  const fvm_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

  if (f->semantics & FVM_FILE_NO_MPI_IO)
    retval = _file_read_block(f,
                              buf,
                              size,
                              _global_num_start,
                              _global_num_end);

#if defined(HAVE_MPI_IO)

  else if (!(f->semantics & FVM_FILE_NO_MPI_IO)) {

    if (f->semantics & FVM_FILE_EXPLICIT_OFFSETS)
      retval = _mpi_file_read_block_eo(f,
                                       buf,
                                       size,
                                       _global_num_start,
                                       _global_num_end);

    else /* if (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS) */
      retval = _mpi_file_read_block_ip(f,
                                       buf,
                                       size,
                                       _global_num_start,
                                       _global_num_end);
  }

#endif /* defined(HAVE_MPI_IO) */

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf, buf, size, retval);

  return retval;
}

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
                     const void  *buf,
                     size_t       size,
                     size_t       stride,
                     fvm_gnum_t   global_num_start,
                     fvm_gnum_t   global_num_end)
{
  size_t retval = 0;

  const size_t bufsize = (global_num_end - global_num_start)*stride*size;

  /* Copy contents to ensure buffer constedness if necessary */

  if (   (f->swap_endian == true && size > 1)
      || (f->n_ranks > 1)
      || !(f->semantics & FVM_FILE_NO_MPI_IO)) {

    unsigned char *copybuf = NULL;

    BFT_MALLOC(copybuf, bufsize, unsigned char);

    memcpy(copybuf, buf, bufsize);

    retval = fvm_file_write_block_buffer(f,
                                         copybuf,
                                         size,
                                         stride,
                                         global_num_start,
                                         global_num_end);

    BFT_FREE(copybuf);
  }

  /* In single-processor case with no byte-swapping, write directly */

  else if (f->sh != NULL) {

    const fvm_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
    const fvm_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

    retval = _file_write(f,
                         buf,
                         size,
                         (_global_num_end - _global_num_start));
  }

  return retval;
}

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
fvm_file_write_block_buffer(fvm_file_t  *f,
                            void        *buf,
                            size_t       size,
                            size_t       stride,
                            fvm_gnum_t   global_num_start,
                            fvm_gnum_t   global_num_end)
{
  size_t retval = 0;

  const fvm_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
  const fvm_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

  /* Swap bytes prior to writing if necessary */

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf,
                 buf,
                 size,
                 (_global_num_end - _global_num_start));

  /* Write to file using chosen method */

  if (f->semantics & FVM_FILE_NO_MPI_IO)
    retval = _file_write_block(f,
                               buf,
                               size,
                               _global_num_start,
                               _global_num_end);

#if defined(HAVE_MPI_IO)

  else if (!(f->semantics & FVM_FILE_NO_MPI_IO)) {

    if (f->semantics & FVM_FILE_EXPLICIT_OFFSETS)
      retval = _mpi_file_write_block_eo(f,
                                        buf,
                                        size,
                                        _global_num_start,
                                        _global_num_end);

    else /* if (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS) */
      retval = _mpi_file_write_block_ip(f,
                                        buf,
                                        size,
                                        _global_num_start,
                                        _global_num_end);
  }

#endif /* defined(HAVE_MPI_IO) */

  return retval;
}

/*----------------------------------------------------------------------------
 * Update the file pointer according to whence.
 *
 * parameters:
 *   f      <-- fvm_file_t descriptor
 *   offset <-- add to position specified to whence to obtain new position,
 *              measured in characters from the beginning of the file
 *   whence <-- beginning if FVM_FILE_SEEK_SET, current if FVM_FILE_SEEK_CUR,
 *              or end-of-file if FVM_FILE_SEEK_END
 *
 * returns:
 *   0 upon success, nonzero otherwise; currently, errors are fatal.
 *----------------------------------------------------------------------------*/

int
fvm_file_seek(fvm_file_t       *f,
              fvm_file_off_t    offset,
              fvm_file_seek_t   whence)
{
  int retval = 0;

  if (f->semantics & FVM_FILE_NO_MPI_IO) {
    if (f->rank == 0)
      retval = _file_seek(f, offset, whence);
  }

#if defined(HAVE_MPI_IO)

  else if (!(f->semantics & FVM_FILE_NO_MPI_IO)) {

    retval = MPI_SUCCESS;

    /* Always update f->offset, regardless of mode */

    switch(whence) {
    case FVM_FILE_SEEK_SET:
      f->offset = offset;
      break;
    case FVM_FILE_SEEK_CUR:
      f->offset += offset;
      break;
    case FVM_FILE_SEEK_END:
      {
        MPI_Offset f_size = 0;
        retval = MPI_File_get_size(f->fh, &f_size);
        f->offset = f_size + offset;
      }
    }

    if (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS)
      retval = MPI_File_seek(f->fh, f->offset, MPI_SEEK_SET);

    if (retval != MPI_SUCCESS)
      _mpi_io_error_message(f->name, retval);
  }

#endif /* defined(HAVE_MPI_IO) */

  return retval;
}

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
fvm_file_tell(fvm_file_t  *f)
{
  fvm_file_off_t retval = 0;

  if (f->semantics & FVM_FILE_NO_MPI_IO) {

    if (f->rank == 0)
      retval = _file_tell(f);

#if defined(HAVE_MPI)
    if (f->comm != MPI_COMM_NULL) {
      int64_t _offset = retval;
      MPI_Bcast(&_offset, 1, fvm_datatype_to_mpi[FVM_INT64], 0, f->comm);
      retval = _offset;
    }
#endif

  }

#if defined(HAVE_MPI_IO)

  else if (!(f->semantics & FVM_FILE_NO_MPI_IO)) {

    /*
      Note that in case of individual file pointers, using
      MPI_File_get_position() and MPI_File_get_byte_offset() should also
      work, but fail after certain collective writes with some processes
      writing zero values (at least on Open MPI 1.2.6), so we prefer to
      keep track of the global offset (which we use to set views anyways).
    */

    retval = f->offset;
  }

#endif /* defined(HAVE_MPI_IO) */

  return retval;
}

/*----------------------------------------------------------------------------
 * Get the default semantics for file access.
 *
 * returns:
 *   current default semantics for file access
 *----------------------------------------------------------------------------*/

fvm_file_hints_t
fvm_file_get_default_semantics(void)
{
  return _default_semantics;
}

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
fvm_file_set_default_semantics(fvm_file_hints_t  hints)
{
  int retval = 0;

  const fvm_file_hints_t mpi_io_hints = (  FVM_FILE_EXPLICIT_OFFSETS
                                         | FVM_FILE_INDIVIDUAL_POINTERS);

#if defined(HAVE_MPI_IO)
  if (   (hints & FVM_FILE_EXPLICIT_OFFSETS)
      && (hints & FVM_FILE_INDIVIDUAL_POINTERS))
    retval = 1;
  else if (   (hints & mpi_io_hints)
           && (hints & FVM_FILE_NO_MPI_IO))
    retval = 1;
#else
  if (hints & mpi_io_hints)
    retval = 1;
#endif

  if (retval == 0)
    _default_semantics = hints;

  return retval;
}

/*----------------------------------------------------------------------------
 * Dump the metadata of a file structure in human readable form
 *
 * parameters:
 *   f <-- pointer to file
 *----------------------------------------------------------------------------*/

void
fvm_file_dump(const fvm_file_t  *f)
{
  const char *mode_name[] = {"FVM_FILE_MODE_READ",
                             "FVM_FILE_MODE_WRITE",
                             "FVM_FILE_MODE_APPEND"};

  if (f == NULL) {
    bft_printf("\n"
               "Null file dump:\n");
    return;
  }

  bft_printf("\n"
             "File name:                \"%s\"\n"
             "Access mode:              %s\n"
             "Semantics:\n"
             "  no_mpi_io:              %d\n"
             "  no_predistribute:       %d\n"
             "  explicit_offsets:       %d\n"
             "  individual_pointers:    %d\n"
             "Rank:                     %d\n"
             "N ranks:                  %d\n"
             "Swap endian:              %d\n"
             "Serial handle:            %p\n",
             f->name, mode_name[f->mode],
             (f->semantics & FVM_FILE_NO_MPI_IO),
             (f->semantics & FVM_FILE_NO_PREDISTRIBUTE) >> 1,
             (f->semantics & FVM_FILE_EXPLICIT_OFFSETS) >> 2,
             (f->semantics & FVM_FILE_INDIVIDUAL_POINTERS) >> 3,
             f->rank, f->n_ranks, (int)(f->swap_endian),
             f->sh);

#if defined(HAVE_MPI)
  bft_printf("Associated communicator:  %llu\n",
             (unsigned long long)(f->comm));
#if defined(HAVE_MPI_IO)
  bft_printf("MPI file handle:          %llu\n"
               "MPI file offset:          %llu\n",
             (unsigned long long)(f->fh),
             (unsigned long long)(f->offset));
#endif
#endif

  bft_printf("\n");
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a fvm_file_serializer_t structure.
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

fvm_file_serializer_t *
fvm_file_serializer_create(size_t        size,
                           size_t        stride,
                           fvm_gnum_t    global_num_start,
                           fvm_gnum_t    global_num_end,
                           size_t        buf_block_size,
                           void         *buf,
                           MPI_Comm      comm)
{
  fvm_file_serializer_t  *s = NULL;

  BFT_MALLOC(s, 1, fvm_file_serializer_t);

  _serializer_init(s,
                   size * stride,
                   global_num_start,
                   global_num_end,
                   buf_block_size,
                   buf,
                   comm);

  return s;
}

/*----------------------------------------------------------------------------
 * Destroy a fvm_file_serializer_t structure.
 *
 * parameters:
 *   s <-- pointer to pointer structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
fvm_file_serializer_destroy(fvm_file_serializer_t  **s)
{
  if (s != NULL) {
    _serializer_finalize(*s);
    BFT_FREE(*s);
  }
}

/*----------------------------------------------------------------------------
 * Advance a fvm_file_serializer_t structure.
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
fvm_file_serializer_advance(fvm_file_serializer_t  *s,
                            fvm_gnum_t              cur_range[2])
{
  MPI_Status status;
  fvm_gnum_t sync_range[2] = {s->next_g_num, 0};

  void *retval = NULL;

  /* Rank 0 receives data */

  if (s->rank_id == 0) {

    int dist_rank = s->next_rank_id;
    int count = 0;

    if (s->next_rank_id >= s->n_ranks)
      return NULL;

    else if (s->next_rank_id != 0) {

      count = s->count[dist_rank];

      /* Forced synchronization */
      sync_range[1] = sync_range[0] + count;
      MPI_Send(&sync_range, 2, FVM_MPI_GNUM, dist_rank, FVM_MPI_TAG, s->comm);

      /* Receive data */
      MPI_Recv(s->recv_buf, (count * s->size), MPI_BYTE, dist_rank,
               FVM_MPI_TAG, s->comm, &status);

      retval = s->recv_buf;
    }

    else { /* First call, rank id 0 */
      count = s->count[dist_rank];
      retval = s->buf;
    }

    /* Update status */

    s->next_rank_id += 1;
    while (s->next_rank_id < s->n_ranks) {
      if (s->count[s->next_rank_id] > 0)
        break;
      else
        s->next_rank_id += 1;
    }

    if (cur_range != NULL) {
      cur_range[0] = s->next_g_num;
      cur_range[1] = cur_range[0] + count;
    }

    s->next_g_num += count;

  }

  /* Other ranks send data */

  else {

    int count = s->range[1] - s->range[0];

    if (count > 0) {

      /* Forced synchronization */
      MPI_Recv(&sync_range, 2, FVM_MPI_GNUM, 0, FVM_MPI_TAG, s->comm, &status);
      count = (sync_range[1] - sync_range[0]);

      if (sync_range[0] != s->range[0] || sync_range[1] != s->range[1])
        bft_error(__FILE__, __LINE__, 0,
                  _("Error serializing data:\n\n"
                    "  requested range: [%llu, %llu[\n"
                    "  local range:     [%llu, %llu["),
                  (unsigned long long)sync_range[0],
                  (unsigned long long)sync_range[1],
                  (unsigned long long)(s->range[0]),
                  (unsigned long long)(s->range[1]));

      /* Send data */
      MPI_Send(s->buf, (count * s->size), MPI_BYTE, 0, FVM_MPI_TAG, s->comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
