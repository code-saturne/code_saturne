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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*
  Force LARGEFILE_SOURCE if largefiles enabled under 32-bit Linux or Blue Gene
  (otherwise, we may encounter bugs with glibc 2.3 due to fseeko end ftello
  not being correctly defined). Compiling with -D_GNU_SOURCE instead
  of -D_POSIX_C_SOURCE=200112L seems to be another way to solve the problem.
*/

#if (SIZEOF_LONG < 8) && (_FILE_OFFSET_BITS == 64)
# if defined(__linux__)
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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H)
# include <sys/stat.h>
# include <sys/types.h>
# if defined(HAVE_UNISTD_H)
#  include <unistd.h>
# endif
#endif /* defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) */

#if defined(HAVE_DIRENT_H)
#include <dirent.h>
#endif

#if defined(WIN32) || defined(_WIN32)
#include <io.h>
#endif

#if defined(HAVE_MPI_IO)
#include <limits.h>
#endif

#if defined(HAVE_ZLIB)
#include <zlib.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_file.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_file.c
        File and directory operations, with parallel IO.

  \typedef cs_file_t
           File descriptor (opaque object)

  \typedef cs_file_off_t
           Offset for file position indicator

  \enum cs_file_mode_t

  \brief File acces modes

  \var CS_FILE_MODE_READ
       Read mode
  \var CS_FILE_MODE_WRITE
       Write mode
  \var CS_FILE_MODE_APPEND
       Append

  \enum cs_file_seek_t

  \brief seek semantics (third argument of \ref cs_file_seek)

  \var CS_FILE_SEEK_SET
       Seek from beginning of file
  \var CS_FILE_SEEK_CUR
       Seek from current position
  \var CS_FILE_SEEK_END
       Seek from end of file

  \enum cs_file_access_t

  \brief Shared file access methods

  \var CS_FILE_STDIO_SERIAL
       Default IO option
  \var CS_FILE_STDIO_SERIAL
       Serial standard C IO (funnelled through rank 0 in parallel)
  \var CS_FILE_STDIO_PARALLEL
       Per-process standard C IO (for reading only)
  \var CS_FILE_MPI_INDEPENDENT
       Non-collective MPI-IO with independent file open and close
       (for reading only)
  \var CS_FILE_MPI_NON_COLLECTIVE
       Non-collective MPI-IO with collective file open and close
  \var CS_FILE_MPI_COLLECTIVE
       Collective MPI-IO

  \enum cs_file_mpi_positioning_t

  \brief MPI-IO positioning methods
  \details It is not always known whether a performance or robustness
          difference is to be expected using explicit file offsets
          or individual file pointers. Perusal of a sampling of ROMIO
          code would seem to indicate that no difference is to be
          expected, but this might change with MPI IO variants
          or file systems, so an advanced setting is made possible.

  \var CS_FILE_MPI_EXPLICIT_OFFSETS
       Use explicit offsets positioning with MPI-IO
  \var CS_FILE_MPI_INDIVIDUAL_POINTERS
       Use individual file pointer positioning with MPI-IO
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* MPI tag for file operations */
#define CS_FILE_MPI_TAG  (int)('C'+'S'+'_'+'F'+'I'+'L'+'E')

/*============================================================================
 * Type definitions
 *============================================================================*/

/* File descriptor */

struct _cs_file_t {

  char              *name;         /* File name */
  cs_file_mode_t     mode;         /* File mode */
  cs_file_access_t   method;       /* File access method */
  int                rank;         /* MPI rank */
  int                n_ranks;      /* MPI rank */
  bool               swap_endian;  /* Swap big-endian and little-endian ? */

  FILE              *sh;           /* Serial file handle */

#if defined(HAVE_ZLIB)
  gzFile             gzh;          /* Zlib (serial) file handle */
#endif

#if defined(HAVE_MPI)
  int                rank_step;    /* Rank step between ranks and io ranks */
  cs_gnum_t         *block_size;   /* Block sizes on IO ranks in case
                                      of rank stepping */
  MPI_Comm           comm;         /* Associated MPI communicator */
  MPI_Comm           io_comm;      /* Associated MPI-IO communicator */
#endif
#if defined(HAVE_MPI_IO)
  MPI_File           fh;           /* MPI file handle */
  MPI_Info           info;         /* MPI file info */
  MPI_Offset         offset;       /* MPI file offset */
#else
  cs_file_off_t      offset;       /* File offset */
#endif

};

/* Associated typedef documentation (for cs_file.h) */

/*!
 * \typedef cs_file_t
 * \brief Pointer to opaque file descriptor
 */

#if defined(HAVE_MPI)

/* Helper structure for IO serialization */

struct _cs_file_serializer_t {

  int          rank_id;        /* Local rank in communicator */
  int          n_ranks;        /* Number of ranks in communicator */

  cs_gnum_t    range[2];       /* Global start and past-the-end numbers
                                  for local rank */

  size_t       size;           /* datatype size (may include stride) */

  cs_gnum_t    next_g_num;     /* Next global number */
  int          next_rank_id;   /* Next rank with which we will communicate */

  cs_lnum_t   *count;          /* Number of elements in each block */

  void        *buf;            /* pointer to external buffer */
  void        *recv_buf;       /* pointer to external buffer if
                                  buf_block_size >= max_block_size,
                                  or to buf otherwise */

  MPI_Comm     comm;           /* Associated MPI communicator */
};

#endif /* defined(HAVE_MPI) */

/* Offset type for zlib */

#if defined(HAVE_ZLIB)

/* Zlib API may be broken when using large file support, as z_off_t
   is based on current off_t, and not on a value fixed at compilation time.
   We redefine prototypes for gzseek() and gztell() ;
   This is ugly, but not as wrong as zlib's logic, and should work with an
   unmodified Zlib (as of Zlib 1.2.11). */

#if defined (SIZEOF_Z_OFF_T)
#  if (SIZEOF_Z_OFF_T == SIZEOF_LONG)
typedef long _cs_z_off_t;
#  elif defined (HAVE_LONG_LONG)
#    if (SIZEOF_Z_OFF_T == SIZEOF_LONG_LONG)
typedef long long _cs_z_off_t;
#    else
#      error "z_off_t returned by zlibCompileFlags() neither long nor long long"
#    endif
#  endif
#else
typedef z_off_t _cs_z_off_t;
#endif

typedef _cs_z_off_t (cs_gzseek_t) (gzFile file,
                                  _cs_z_off_t offset,
                                   int whence);

typedef _cs_z_off_t (cs_gztell_t) (gzFile file);

#endif /* HAVE_ZLIB */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default access */

static cs_file_mpi_positioning_t
  _mpi_io_positioning = CS_FILE_MPI_EXPLICIT_OFFSETS;

static cs_file_access_t _default_access_r = CS_FILE_DEFAULT;
static cs_file_access_t _default_access_w = CS_FILE_DEFAULT;

/* Communicator and hints used for file operations */

#if defined(HAVE_MPI)

static bool     _mpi_defaults_are_set = false;
static int      _mpi_rank_step = 1;
static MPI_Comm _mpi_comm = MPI_COMM_NULL;
static MPI_Comm _mpi_io_comm = MPI_COMM_NULL;
static MPI_Info _mpi_io_hints_r = MPI_INFO_NULL;
static MPI_Info _mpi_io_hints_w = MPI_INFO_NULL;

#endif

#if defined(HAVE_ZLIB)

/* Zlib API broken offset size workaround, continued... */

static cs_gzseek_t  *_cs_gzseek = (cs_gzseek_t *)gzseek;
static cs_gztell_t  *_cs_gztell = (cs_gztell_t *)gztell;

#endif /* HAVE_ZLIB */

/*============================================================================
 * Global variables
 *============================================================================*/

/* names associated with file I/O methods */

const char  *cs_file_access_name[]
  = {N_("default"),
     N_("standard input and output, serial access"),
     N_("standard input and output, parallel access"),
     N_("non-collective MPI-IO, independent file open/close"),
     N_("non-collective MPI-IO, collective file open/close"),
     N_("collective MPI-IO")};

/* names associated with MPI-IO positioning */

#if defined(HAVE_MPI_IO)
const char *cs_file_mpi_positioning_name[] = {N_("explicit offsets"),
                                              N_("individual file pointers")};
#endif

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Evaluate an access method, transforming default to actual value.
 *
 * parameters:
 *   m <-- access method
 *   w <-- true if write access (false for readonly)
 *
 * returns:
 *   actual access method
 *----------------------------------------------------------------------------*/

static cs_file_access_t
_access_method(cs_file_access_t  m,
               bool              w)
{
  cs_file_access_t  _m = m;

  /* Handle default */

  if (_m == CS_FILE_DEFAULT) {

#if defined(HAVE_MPI)
#  if defined(HAVE_MPI_IO)
    _m = CS_FILE_MPI_COLLECTIVE;
#  else
    _m = CS_FILE_STDIO_PARALLEL;
#  endif
#else
    _m = CS_FILE_STDIO_SERIAL;
#endif

  }

  /* Restrict to possible values */

#if defined(HAVE_MPI)
#  if !defined(HAVE_MPI_IO)
  _m = CS_MAX(_m, CS_FILE_STDIO_PARALLEL);
#  endif
  if (cs_glob_mpi_comm == MPI_COMM_NULL)
    _m = CS_FILE_STDIO_SERIAL;
#else
  _m = CS_FILE_STDIO_SERIAL;
#endif

  if (w && _m == CS_FILE_STDIO_PARALLEL)
    _m = CS_FILE_STDIO_SERIAL;

  return _m;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize an cs_file_serializer_t structure.
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
_serializer_init(cs_file_serializer_t  *s,
                 size_t                 size,
                 cs_gnum_t              global_num_start,
                 cs_gnum_t              global_num_end,
                 size_t                 buf_block_size,
                 void                  *buf,
                 MPI_Comm               comm)
{
  cs_lnum_t l_count = 0;

  s->range[0] = global_num_start;
  s->range[1] = global_num_end;

  s->size = size;

  if (s->range[1] > s->range[0])
    l_count = s->range[1] - s->range[0];

  /* Get local rank and size of the current MPI communicator */

  if (comm != MPI_COMM_NULL) {

    MPI_Comm_rank(comm, &(s->rank_id));
    MPI_Comm_size(comm, &(s->n_ranks));

    s->next_rank_id = 0;
    s->next_g_num = global_num_start;

    /* Initialize counter */

    if (s->rank_id == 0)
      BFT_MALLOC(s->count, s->n_ranks, cs_lnum_t);
    else
      s->count = NULL;

    MPI_Gather(&l_count, 1, CS_MPI_LNUM, s->count, 1, CS_MPI_LNUM, 0, comm);

    /* Allocate local buffer if necessary, or point to external buffer */

    s->buf = buf;
    s->recv_buf = NULL;

    if (s->rank_id == 0) {
      int i;
      cs_lnum_t _max_block_size = 0;
      cs_lnum_t _buf_block_size = CS_MAX((cs_lnum_t)buf_block_size, l_count);
      for (i = 0; i < s->n_ranks; i++)
        _max_block_size = CS_MAX(_max_block_size, s->count[i]);
      if (_max_block_size > _buf_block_size)
        BFT_MALLOC(s->recv_buf, _max_block_size*size, unsigned char);
      else
        s->recv_buf = buf;
    }

  }

  else { /* if (comm == MPI_COMM_NULL) */

    s->rank_id = -1;
    s->n_ranks = 0;

    s->next_rank_id = 0;
    s->next_g_num = 0;

    s->count = NULL;

    s->buf = buf;
    s->recv_buf = NULL;

  }

  s->comm = comm;
}

/*----------------------------------------------------------------------------
 * Finalize an cs_file_serializer_t structure.
 *
 * parameters:
 *   s <-- pointer to structure that should be finalized
 *----------------------------------------------------------------------------*/

static void
_serializer_finalize(cs_file_serializer_t  *s)
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
 *
 * returns:
 *   0 in case of success, error number in case of failure
 *----------------------------------------------------------------------------*/

static int
_file_open(cs_file_t  *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->sh != NULL)
    return 0;

  /* Compressed with gzip ? (currently for reading only) */

#if defined(HAVE_ZLIB)

  if (f->gzh != NULL)
    return 0;

  if (f->mode == CS_FILE_MODE_READ) {

    bool gzipped = false;

    size_t l = strlen(f->name);
    if (l > 3 && (strncmp((f->name + l-3), ".gz", 3) == 0))
      gzipped = true;

    if (gzipped) {
      f->gzh = gzopen(f->name, "r");

      if (f->gzh == NULL) {
        const char *err_str
          = (errno == 0) ? zError(Z_MEM_ERROR) : strerror(errno);
        retval = (errno == 0) ? Z_MEM_ERROR : errno;
        bft_error(__FILE__, __LINE__, 0,
                  _("Error opening file \"%s\":\n\n"
                    "  %s"), f->name, err_str);
      }
      return retval;
    }

  }

#endif

  /* The file handler exists and the corresponding file is closed */

  switch (f->mode) {
  case CS_FILE_MODE_APPEND:
    if (f->rank == 0)
      f->sh = fopen(f->name, "ab");
    else
      f->sh = fopen(f->name, "a+b");
    break;
  case CS_FILE_MODE_WRITE:
    if (f->rank == 0)
      f->sh = fopen(f->name, "wb");
    else
      f->sh = fopen(f->name, "a+b");
    break;
  default:
    assert(f->mode == CS_FILE_MODE_READ);
    f->sh = fopen(f->name, "rb");
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
_file_close(cs_file_t  *f)
{
  int retval = 0;

  if (f->sh != NULL)
    retval = fclose(f->sh);

  /* Compressed with gzip ? (currently for reading only) */

#if defined(HAVE_ZLIB)

  else if (f->gzh != NULL) {
    retval = gzclose(f->gzh);
    if (retval != 0) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error closing file \"%s\":\n\n"
                  "  %s"), f->name, gzerror(f->gzh, &retval));
      return retval;
    }
    f->gzh = NULL;
  }

#endif

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
 *   f    <-- cs_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_read(cs_file_t  *f,
           void       *buf,
           size_t      size,
           size_t      ni)
{
  size_t retval = 0;

  if (f->sh != NULL) {

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

#if defined(HAVE_ZLIB)

  else if (f->gzh != NULL) {

    if (ni != 0)
      retval = fread(buf, size, ni, f->sh);

    size_t rec_size = size * ni;

    retval = ((size_t)gzread(f->gzh, buf, rec_size)) / size;

    if (retval != ni) {
      int err_num = 0;
      const char *err_str = gzerror(f->gzh, &err_num);
      if (err_num != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error reading file \"%s\":\n\n  %s"),
                  f->name, err_str);
      else if (gzeof(f->gzh) != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Premature end of file \"%s\""), f->name);
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Error reading file \"%s\""), f->name);
    }

    return retval;

  }

#endif /* defined(HAVE_ZLIB) */

  assert(0);

  return retval;
}

/*----------------------------------------------------------------------------
 * Write data to a file using standard C IO.
 *
 * parameters:
 *   f    <-- cs_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_write(cs_file_t   *f,
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
 *   whence <-- beginning if CS_FILE_SEEK_SET, current if
 *              CS_FILE_SEEK_CUR, or end-of-file if CS_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise.
 *----------------------------------------------------------------------------*/

static int
_file_seek(cs_file_t       *f,
           cs_file_off_t    offset,
           cs_file_seek_t   whence)
{
  static int _stdio_seek[3] = {SEEK_SET, SEEK_CUR, SEEK_END};

  int _whence = _stdio_seek[whence];
  int retval = 0;

  const char err_fmt[] = "Error setting position in file \"%s\":\n\n  %s";

  /* Convert cs_file_seek to stdio values */

  assert(f != NULL);

  if (f->sh != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support may be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)

    retval = fseeko(f->sh, (off_t)offset, _whence);

    if (retval != 0)
      bft_error(__FILE__, __LINE__, errno, _(err_fmt),
                f->name, strerror(errno));
# else

    /* Test if offset larger than allowed */

    long _offset = offset;

    if (_offset == offset) {
      retval = fseek(f->sh, (long)offset, _whence);
      if (retval != 0)
        bft_error(__FILE__, __LINE__, errno, _(err_fmt),
                  f->name, strerror(errno));
    }
    else {
      retval = -1;
      bft_error
        (__FILE__, __LINE__, 0, _(err_fmt),
         f->name,
         _("sizeof(off_t) > sizeof(long) but fseeko() not available"));
    }

# endif /* defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64) */

#else /* SIZEOF_LONG >= 8 */

    /* For 64-bit systems, standard fseek should be enough */

    retval = fseek(f->sh, (long)offset, _whence);
    if (retval != 0)
      bft_error(__FILE__, __LINE__, errno, _(err_fmt),
                f->name, strerror(errno));

#endif /* SIZEOF_LONG */
  }

#if defined(HAVE_ZLIB)

  else if (f->gzh != NULL) {

    retval = _cs_gzseek(f->gzh, (_cs_z_off_t)offset, _whence);

    if (retval != 0) {
      int err_num = 0;
      const char *err_str = gzerror(f->gzh, &err_num);
      if (err_num == 0)
        err_str = "";

      bft_error(__FILE__, __LINE__, 0, _(err_fmt),
                f->name, err_str);
    }
  }

#endif

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

static cs_file_off_t
_file_tell(cs_file_t  *f)
{
  cs_file_off_t offset = 0;

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

#if defined(HAVE_ZLIB)

  else if (f->gzh != NULL) {
    offset = (cs_file_off_t)_cs_gztell(f->gzh);

    if (offset < 0) {
      int err_num = 0;
      const char *err_str = gzerror(f->gzh, &err_num);
      if (err_num == 0)
        err_str = "";

      bft_error(__FILE__, __LINE__, 0,
                _("Error obtaining position in file \"%s\":\n\n  %s"),
                f->name, err_str);
    }
  }

#endif

  return offset;
}

/*----------------------------------------------------------------------------
 * Formatted input from a text file if possible (as fgets()).
 *
 * This function is the base for ecs_file_gets() and ecs_file_gets_try();
 * depending on the allow_eof parameter, failure to read a line due to
 * an end-of-file condition is considered an error or not.
 *
 * parameters:
 *   s:         --> buffer to which string is to be read.
 *   size:      <-- maximum number of characters to be read plus one.
 *   f:         <-- ecs_file_t descriptor.
 *   line:      <-> file line number if available, or NULL.
 *   allow_eof: <-- 1 if EOF is allowed, 0 if considered an error.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 *----------------------------------------------------------------------------*/

static char *
_cs_file_gets(char             *s,
              const int         size,
              const cs_file_t  *f,
              int              *line,
              const int         allow_eof)
{
  char *retval = NULL;

  assert(f != NULL);

  if (f->sh != NULL)
    retval = fgets(s, size, f->sh);

#if defined(HAVE_ZLIB)

  else if (f->gzh != NULL)
    retval = gzgets(f->gzh, s, size);

#endif /* defined(HAVE_ZLIB) */

  else {
    if (cs_glob_n_ranks > 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error: reading from file \"%s\",\n"
                  "       which is not open on rank %d."),
              f->name, cs_glob_rank_id);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error: reading from file \"%s\",\n"
                  "       which is not open."),
                f->name);
  }

  if (retval != NULL) {

    /* Convert Windows type line ending to Unix type line ending if needed */
    int i = strlen(s) - 2;
    if (i > 0) {
      if (s[i] == '\r' && s[i+1] == '\n') {
        s[i] = '\n';
        s[i+1] = '\0';
      }
    }

    if (line != NULL)
      *line += 1;

    return retval;
  }

  /* We should reach this point only in case of a failed read */

  assert(retval == NULL);

  int is_eof = 0;
  if (allow_eof) {
    if (feof(f->sh) != 0)
      is_eof = 1;

#if defined(HAVE_ZLIB)
    else if (gzeof(f->gzh) != 0)
      is_eof = 1;
#endif
  }

  if (allow_eof == 0 || is_eof == 0) {

    const char *err_str = cs_empty_string;

    if (f->sh != NULL) {
      int err_num = ferror(f->sh);
      if (err_num != 0)
        err_str = strerror(err_num);
    }

#if defined(HAVE_ZLIB)

    else if (f->gzh != NULL) {
      int err_num = 0;
      err_str = gzerror(f->gzh, &err_num);
      if (err_num == 0)
        err_str = cs_empty_string;
    }

#endif /* defined(HAVE_ZLIB) */

    if (line != NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Error reading line %d of file \"%s\":\n\n  %s"),
                *line, f->name, err_str);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error reading text file \"%s\":\n\n  %s"),
                f->name, err_str);
  }

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
 * This version does not use MPI-IO
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
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
_file_read_block_s(cs_file_t  *f,
                   void       *buf,
                   size_t      size,
                   cs_gnum_t   global_num_start,
                   cs_gnum_t   global_num_end)
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

    cs_lnum_t loc_count = global_num_end - global_num_start;
    int _counts[64];
    int *counts = NULL;

    MPI_Datatype ent_type = MPI_BYTE;
    size_t _size = size;

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
      cs_lnum_t _buf_size = global_num_end - global_num_start;
      unsigned char *_buf = NULL;

      /* Allocate exchange buffer */

      for (dist_rank = 1; dist_rank < f->n_ranks; dist_rank++)
        _buf_size = CS_MAX(_buf_size, counts[dist_rank]);

      BFT_MALLOC(_buf, _buf_size*size, unsigned char);

      if (_buf_size*size > INT_MAX) {
        MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
        MPI_Type_commit(&ent_type);
        _size = 1;
      }

      /* Loop on distant ranks */

      for (dist_rank = 1; dist_rank < f->n_ranks; dist_rank++) {

        if (counts[dist_rank] == 0)
          continue;

        /* Read data from file */

        counts[dist_rank]
          = (int)_file_read(f, _buf, size, (size_t)counts[dist_rank]);

        /* Send to corresponding rank */

        MPI_Send(_buf, counts[dist_rank]*_size, ent_type, dist_rank,
                 CS_FILE_MPI_TAG, f->comm);

      } /* End of loop on distant ranks */

      BFT_FREE(_buf);

    }

    /* Other ranks receive data from rank 0 */

    else if (loc_count > 0) {

      if (loc_count*size > INT_MAX) {
        MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
        MPI_Type_commit(&ent_type);
        _size = 1;
      }

      /* Receive data */

      MPI_Recv(buf, (int)(loc_count*_size), ent_type, 0,
               CS_FILE_MPI_TAG, f->comm, &status);

      MPI_Get_count(&status, ent_type, &loc_count);
      retval = loc_count / _size;

    }

    if (ent_type != MPI_BYTE)
      MPI_Type_free(&ent_type);

    if (counts != NULL && counts != _counts)
      BFT_FREE(counts);
  }

#endif /* defined(HAVE_MPI) */

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
 * This version does not use MPI-IO
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
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
_file_read_block_p(cs_file_t  *f,
                   void       *buf,
                   size_t      size,
                   cs_gnum_t   global_num_start,
                   cs_gnum_t   global_num_end)
{
  size_t retval = 0;
  cs_gnum_t loc_count = global_num_end - global_num_start;

  if (loc_count > 0) {

    /* Only rank 0 initially opened (to check existence/rights, and
       as all ranks might not participate), so open here if needed */

    cs_file_off_t offset = f->offset + ((global_num_start - 1) * size);

    if (f->sh == NULL)
      _file_open(f);

    if (_file_seek(f, offset, CS_FILE_SEEK_SET) == 0)
      retval = _file_read(f, buf, size, (size_t)loc_count);

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
 * This version does not use MPI-IO
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
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
_file_write_block_s(cs_file_t  *f,
                    void       *buf,
                    size_t      size,
                    cs_gnum_t   global_num_start,
                    cs_gnum_t   global_num_end)
{
  size_t retval = 0;

  if (f->n_ranks == 1)
    retval = _file_write(f,
                         buf,
                         size,
                         (size_t)(global_num_end - global_num_start));

#if defined(HAVE_MPI)

  if (f->n_ranks > 1) {

    cs_file_serializer_t  s;
    cs_lnum_t  local_count;
    cs_lnum_t *count = NULL;
    void  *write_buf = NULL;

    _serializer_init(&s,
                     size,
                     global_num_start,
                     global_num_end,
                     0,
                     buf,
                     f->io_comm);

    do {

      int dist_rank = s.next_rank_id;

      write_buf = cs_file_serializer_advance(&s, NULL);

      if (write_buf != NULL) /* only on rank 0 */
        s.count[dist_rank]
          = (cs_lnum_t)_file_write(f,
                                   write_buf,
                                   size,
                                   (size_t)(s.count[dist_rank]));

    } while (write_buf != NULL);

    /* Exchange return codes */

    if (s.rank_id == 0)
      count = s.count;
    else
      BFT_MALLOC(count, s.n_ranks, cs_lnum_t);

    MPI_Scatter(count, 1, CS_MPI_LNUM,
                &local_count, 1, CS_MPI_LNUM,
                0, f->comm);
    retval = local_count;

    if (s.rank_id != 0)
      BFT_FREE(count);

    _serializer_finalize(&s);
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
 *   f                <-- cs_file_t descriptor
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
_file_write_block_p(cs_file_t  *f,
                    void       *buf,
                    size_t      size,
                    cs_gnum_t   global_num_start,
                    cs_gnum_t   global_num_end)
{
  size_t retval = 0;
  cs_gnum_t loc_count = 0;

  if (global_num_end > global_num_start) {

    loc_count = global_num_end - global_num_start;

    if (f->n_ranks == 1)
      retval = _file_write(f, buf, size, (size_t)loc_count);

#if defined(HAVE_MPI)

    if (f->n_ranks > 1) {

      cs_file_off_t offset = f->offset + ((global_num_start - 1) * size);

      /* Only rank 0 initially opened (to check existence/rights, as
         all ranks might not participate), so open here if needed */

      if (f->sh == NULL)
        _file_open(f);

      if (_file_seek(f, offset, SEEK_SET) == 0)
        retval = _file_write(f, buf, size, (size_t)loc_count);

    }

#endif /* defined(HAVE_MPI) */

  }

  return retval;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Gather blocks sizes across several ranks and allocates matching buffer
 *
 * The caller is responsible for freeing the returned buffer once it is
 * no longer needed.
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-> pointer to global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   pointer to gathered values buffer for gathering rank, NULL for others
 *----------------------------------------------------------------------------*/

static void *
_gather_block_sizes(cs_file_t   *f,
                    size_t       size,
                    cs_gnum_t    global_num_start,
                    cs_gnum_t   *global_num_end)
{
  unsigned char *gather_buf = NULL;

  assert(f != NULL);

  cs_gnum_t _global_num_end = *global_num_end;

  static const int tag = 'f'+'a'+'g'+'g'+'r'+'e'+'g'+'a'+'t'+'e';

  /* Aggregator rank */

  if (f->rank % f->rank_step == 0) {

    f->block_size[0] = _global_num_end - global_num_start;
    size_t block_size = f->block_size[0];

    int rank_end = f->rank + f->rank_step;
    if (rank_end >= f->n_ranks)
      rank_end = f->n_ranks;

    int n_aggr = rank_end - f->rank;

    /* Receive counts */

    for (int i = 1; i < n_aggr; i++) {
      int src_rank = f->rank + i;
      MPI_Status status;
      MPI_Recv(f->block_size + i, 1, CS_MPI_GNUM,
               src_rank, tag, f->comm, &status);
      block_size += f->block_size[i];
    }

    /* Allocate buffer */

    size_t alloc_size = size * (size_t)block_size;
    BFT_MALLOC(gather_buf, alloc_size, unsigned char);

    *global_num_end = global_num_start + block_size;
  }

  /* Sending rank */

  else {

    int dest_rank = f->rank - (f->rank % f->rank_step);
    cs_gnum_t block_size = _global_num_end - global_num_start;
    f->block_size[0] = block_size;

    /* Send counts */

    MPI_Send(&block_size, 1, CS_MPI_GNUM,
             dest_rank, tag, f->comm);

    *global_num_end = global_num_start;  /* For empty message */
  }

  return (void *)gather_buf;
}

/*----------------------------------------------------------------------------
 * Gather blocks across several ranks
 *
 * The caller is responsible for freeing the returned buffer once it is
 * no longer needed.
 *
 * parameters:
 * parameters:
 *   f                <-- cs_file_t descriptor
 *   buf              <-> pointer to location containing data
 *   size             <-- size of each item of data in bytes
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-> pointer to global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   pointer to gathered values buffer for gathering rank, NULL for others
 *----------------------------------------------------------------------------*/

static void *
_gather_blocks(cs_file_t   *f,
               void        *buf,
               size_t       size,
               cs_gnum_t    global_num_start,
               cs_gnum_t   *global_num_end)
{
  assert(f != NULL);

  unsigned char *gather_buf = _gather_block_sizes(f,
                                                  size,
                                                  global_num_start,
                                                  global_num_end);

  MPI_Datatype ent_type = MPI_BYTE;
  size_t _size = size;

  static const int tag = 'f'+'a'+'g'+'g'+'r'+'e'+'g'+'a'+'t'+'e';

  /* Aggregator rank */

  if (f->rank % f->rank_step == 0) {

    int rank_end = f->rank + f->rank_step;
    if (rank_end >= f->n_ranks)
      rank_end = f->n_ranks;

    int n_aggr = rank_end - f->rank;

    /* Precaution for large messages */

    for (int i = 1; i < n_aggr; i++) {
      if (_size > 1 && (size * (size_t)(f->block_size[i])) > INT_MAX) {
        MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
        MPI_Type_commit(&ent_type);
        _size = 1;
      }
    }

    /* Copy local data to gather buffer */

    size_t gather_buf_count = size * (size_t)(f->block_size[0]);
    memcpy(gather_buf, buf, gather_buf_count);

    /* Receive data */

    for (int i = 1; i < n_aggr; i++) {
      int src_rank = f->rank + i;
      MPI_Status status;
      size_t add_buf_count = size * (size_t)(f->block_size[i]);
      int recv_count = f->block_size[i];
      if (recv_count == 0)
        continue;
      if (size * (size_t)recv_count > INT_MAX) {
        MPI_Recv(gather_buf + gather_buf_count, recv_count, ent_type,
                 src_rank, tag, f->comm, &status);
      }
      else {
        recv_count *= size;
        MPI_Recv(gather_buf + gather_buf_count, recv_count, MPI_BYTE,
                 src_rank, tag, f->comm, &status);
      }
      gather_buf_count += add_buf_count;
    }

  }

  /* Sending rank */

  else {

    int dest_rank = f->rank - (f->rank % f->rank_step);
    cs_gnum_t block_size = f->block_size[0];

    size_t message_size = _size * (size_t)block_size;

    /* Precaution for large messages */

    if (message_size > INT_MAX) {
      MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
      MPI_Type_commit(&ent_type);
      _size = 1;
    }

    int send_count = _size * block_size;

    /* Send data */

    if (send_count > 0)
      MPI_Send(buf, send_count, ent_type, dest_rank, tag, f->comm);
  }

  if (ent_type != MPI_BYTE)
    MPI_Type_free(&ent_type);

  return (void *)gather_buf;
}

/*----------------------------------------------------------------------------
 * Gather blocks across several ranks
 *
 * The caller is responsible for freeing the returned buffer once it is
 * no longer needed.
 *
 * parameters:
 * parameters:
 *   f                <-- cs_file_t descriptor
 *   io_buf           <-> pointer to location containing read data
 *   buf              <-> pointer to location containing scattered data
 *   size             <-- size of each item of data in bytes
 *
 * returns:
 *   number of values in block after scatter
 *----------------------------------------------------------------------------*/

static int
_scatter_blocks(cs_file_t   *f,
                void        *io_buf,
                void        *buf,
                size_t       size)
{
  assert(f != NULL);

  MPI_Datatype ent_type = MPI_BYTE;
  size_t _size = size;

  static const int tag = 'f'+'a'+'g'+'g'+'r'+'e'+'g'+'a'+'t'+'e';

  /* Aggregator rank */

  if (f->rank % f->rank_step == 0) {

    int rank_end = f->rank + f->rank_step;
    if (rank_end >= f->n_ranks)
      rank_end = f->n_ranks;

    int n_aggr = rank_end - f->rank;

    /* Precaution for large messages */

    for (int i = 1; i < n_aggr; i++) {
      if (_size > 1 && (size * (size_t)(f->block_size[i])) > INT_MAX) {
        MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
        MPI_Type_commit(&ent_type);
        _size = 1;
      }
    }

    /* Send local data to destination buffer */

    unsigned char *scatter_buf = io_buf;
    size_t scatter_buf_count = size * (size_t)(f->block_size[0]);
    memcpy(buf, io_buf, scatter_buf_count);

    /* Send data */

    for (int i = 1; i < n_aggr; i++) {
      int src_rank = f->rank + i;
      size_t add_buf_count = size * (size_t)(f->block_size[i]);
      int send_count = f->block_size[i];
      if (send_count == 0)
        continue;
      if (size * (size_t)send_count > INT_MAX) {
        MPI_Send(scatter_buf + scatter_buf_count, send_count, ent_type,
                 src_rank, tag, f->comm);
      }
      else {
        send_count *= size;
        MPI_Send(scatter_buf + scatter_buf_count, send_count, MPI_BYTE,
                 src_rank, tag, f->comm);
      }
      scatter_buf_count += add_buf_count;
    }

  }

  /* Receving rank */

  else {

    int dest_rank = f->rank - (f->rank % f->rank_step);
    cs_gnum_t block_size = f->block_size[0];

    size_t message_size = _size * (size_t)block_size;

    /* Precaution for large messages */

    if (message_size > INT_MAX) {
      MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
      MPI_Type_commit(&ent_type);
      _size = 1;
    }

    int recv_count = _size * block_size;
    MPI_Status status;

    /* Receive data */

    if (recv_count > 0)
      MPI_Recv(buf, recv_count, ent_type, dest_rank, tag, f->comm, &status);
  }

  if (ent_type != MPI_BYTE)
    MPI_Type_free(&ent_type);

  return f->block_size[0];
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_MPI_IO)

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
 *   mode  <-- file access mode: read, write, or append
 *
 * returns:
 *   MPI_SUCCESS in case of success, MPI error code in case of failure
 *----------------------------------------------------------------------------*/

static int
_mpi_file_open(cs_file_t       *f,
               cs_file_mode_t   mode)
{
  int amode = MPI_MODE_RDWR;
  int retval = 0;

  assert(f != NULL);

  if (f->fh != MPI_FILE_NULL)
    return 0;

  /* Set access mode */

  f->mode = mode;

  if (f->mode == CS_FILE_MODE_APPEND)
    amode = MPI_MODE_WRONLY | MPI_MODE_APPEND;

  else if (f->mode == CS_FILE_MODE_WRITE) {
    int rank;
    if (f->method == CS_FILE_MPI_INDEPENDENT && f->rank > 0)
      amode = MPI_MODE_WRONLY;
    else
      amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    MPI_Comm_rank(f->comm, &rank);
    if (rank < 1)
      cs_file_remove(f->name);
  }

  else if (f->mode == CS_FILE_MODE_READ)
    amode = MPI_MODE_RDONLY;

  /* Open file (for independent access, only on rank 0 initially) */

  if (f->io_comm != MPI_COMM_NULL) {
    retval = MPI_File_open(f->io_comm, f->name, amode, f->info, &(f->fh));
    if (retval == MPI_SUCCESS)
      retval = MPI_File_get_position(f->fh, &(f->offset));
  }

  if (retval != MPI_SUCCESS)
    _mpi_io_error_message(f->name, retval);

  if (f->mode == CS_FILE_MODE_APPEND)
    f->offset = cs_file_tell(f);

  return retval;
}

/*----------------------------------------------------------------------------
 * Open a file independently of other ranks if required using MPI IO.
 *
 * This function is used in the case of independent file IO, to allow
 * files to be opened only on ranks reading/writing nonempty blocks.
 *
 * parameters:
 *   f     <-- pointer to file handler
 *
 * returns:
 *   MPI_SUCCESS in case of success, MPI error code in case of failure
 *----------------------------------------------------------------------------*/

static int
_mpi_file_ensure_isopen(cs_file_t *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->io_comm != MPI_COMM_NULL && f->fh == MPI_FILE_NULL) {

    int amode = MPI_MODE_RDWR;
    if (f->mode == CS_FILE_MODE_APPEND)
      amode = MPI_MODE_WRONLY | MPI_MODE_APPEND;
    else if (f->mode == CS_FILE_MODE_WRITE)
      amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    else if (f->mode == CS_FILE_MODE_READ)
      amode = MPI_MODE_RDONLY;

    retval = MPI_File_open(MPI_COMM_SELF, f->name, amode, f->info, &(f->fh));
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
_mpi_file_close(cs_file_t  *f)
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
 * There are 3 variants, depending on the semantics:
 *   _mpi_file_read_block_noncoll (non-collective)
 *   _mpi_file_read_block_eo (using explicit offsets)
 *   _mpi_file_read_block_ip (using individual pointers, setting a file view)
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
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
_mpi_file_read_block_noncoll(cs_file_t  *f,
                             void       *buf,
                             size_t      size,
                             cs_gnum_t   global_num_start,
                             cs_gnum_t   global_num_end)
{
  cs_gnum_t gcount = (global_num_end - global_num_start)*size;
  size_t retval = 0;

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > 0) {

    int errcode, count;
    MPI_Status status;

    MPI_Offset disp = f->offset + ((global_num_start - 1) * size);
    MPI_Datatype ent_type = MPI_BYTE;

    if (gcount > INT_MAX) {
      MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
      MPI_Type_commit(&ent_type);
      count = global_num_end - global_num_start;
    }
    else
      count = gcount;

    errcode = _mpi_file_ensure_isopen(f);

    if (errcode == MPI_SUCCESS) {

      if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS)
        errcode = MPI_File_read_at(f->fh, disp, buf, count, ent_type, &status);

      else {
        errcode = MPI_File_seek(f->fh, disp, MPI_SEEK_SET);
        if (errcode == MPI_SUCCESS)
          errcode = MPI_File_read(f->fh, buf, count, ent_type, &status);
      }

    }

    if (errcode != MPI_SUCCESS)
      _mpi_io_error_message(f->name, errcode);

    MPI_Get_count(&status, ent_type, &count);

    if (ent_type != MPI_BYTE) {
      MPI_Type_free(&ent_type);
      retval = count;
    }
    else
      retval = count / size;

  }

  return retval;
}

static size_t
_mpi_file_read_block_eo(cs_file_t  *f,
                        void       *buf,
                        size_t      size,
                        cs_gnum_t   global_num_start,
                        cs_gnum_t   global_num_end)
{
  MPI_Status status;
  int errcode, count;
  cs_gnum_t gcount = (global_num_end - global_num_start)*size;
  MPI_Datatype ent_type = MPI_BYTE;
  MPI_Offset disp = f->offset + ((global_num_start - 1) * size);

  size_t retval = 0;

  assert(gcount == 0 || f->fh != MPI_FILE_NULL);

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > INT_MAX) {
    MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
    MPI_Type_commit(&ent_type);
    count = global_num_end - global_num_start;
  }
  else
    count = gcount;

  errcode = MPI_File_read_at_all(f->fh, disp, buf, count, ent_type, &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  if (count > 0)
    MPI_Get_count(&status, ent_type, &count);

  if (ent_type != MPI_BYTE) {
    MPI_Type_free(&ent_type);
    retval = count;
  }
  else {
    if (count > 0)
      retval = count / size;
    else
      retval = 0;
  }

  return retval;
}

static size_t
_mpi_file_read_block_ip(cs_file_t  *f,
                        void       *buf,
                        size_t      size,
                        cs_gnum_t   global_num_start,
                        cs_gnum_t   global_num_end)
{
  int errcode;
  int lengths[1];
  MPI_Aint disps[1];
  MPI_Status status;
  MPI_Datatype file_type;

  int count = 0;
  char datarep[] = "native";
  MPI_Datatype ent_type = MPI_BYTE;
  cs_gnum_t gcount = (global_num_end - global_num_start) * size;
  cs_gnum_t gdisp = (global_num_start - 1) * size;

  size_t retval = 0;

  assert(gcount == 0 || f->fh != MPI_FILE_NULL);

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > INT_MAX || gdisp > INT_MAX) {
    MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
    MPI_Type_commit(&ent_type);
    lengths[0] = global_num_end - global_num_start;
    disps[0] = global_num_start - 1;
  }
  else {
    lengths[0] = gcount;
    disps[0] = gdisp;
  }

  MPI_Type_create_hindexed(1, lengths, disps, ent_type, &file_type);
  MPI_Type_commit(&file_type);

  MPI_File_set_view(f->fh, f->offset, ent_type, file_type, datarep, f->info);

  errcode = MPI_File_read_all(f->fh, buf, lengths[0], ent_type, &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  MPI_Type_free(&file_type);

  if (lengths[0] > 0)
    MPI_Get_count(&status, ent_type, &count);

  if (ent_type != MPI_BYTE) {
    MPI_Type_free(&ent_type);
    retval = count;
  }
  else
    retval = count / size;

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
 * There are 3 variants, depending on the semantics:
 *   _mpi_file_write_block_noncoll (non-collective)
 *   _mpi_file_write_block_eo (using explicit offsets)
 *   _mpi_file_write_block_ip (using individual pointers, setting a file view)
 *
 * parameters:
 *   f                <-- cs_file_t descriptor
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
_mpi_file_write_block_noncoll(cs_file_t  *f,
                              void       *buf,
                              size_t      size,
                              cs_gnum_t   global_num_start,
                              cs_gnum_t   global_num_end)
{
  cs_gnum_t gcount = (global_num_end - global_num_start)*size;
  size_t retval = 0;

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > 0) {

    int errcode, count;
    MPI_Status status;
    MPI_Offset disp = f->offset + ((global_num_start - 1) * size);
    MPI_Datatype ent_type = MPI_BYTE;

    if (gcount > INT_MAX) {
      MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
      MPI_Type_commit(&ent_type);
      count = global_num_end - global_num_start;
    }
    else
      count = gcount;

    errcode = _mpi_file_ensure_isopen(f);

    if (errcode == MPI_SUCCESS) {

      if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS)
        errcode = MPI_File_write_at(f->fh, disp, buf, count, ent_type, &status);

      else {
        errcode = MPI_File_seek(f->fh, disp, MPI_SEEK_SET);
        if (errcode == MPI_SUCCESS)
          errcode = MPI_File_write(f->fh, buf, count, ent_type, &status);
      }

    }

    if (errcode != MPI_SUCCESS)
      _mpi_io_error_message(f->name, errcode);

    if (count > 0)
      MPI_Get_count(&status, ent_type, &count);

    if (ent_type != MPI_BYTE) {
      MPI_Type_free(&ent_type);
      retval = count;
    }
    else
      retval = count / size;

  }

  return retval;
}

static size_t
_mpi_file_write_block_eo(cs_file_t  *f,
                         void       *buf,
                         size_t      size,
                         cs_gnum_t   global_num_start,
                         cs_gnum_t   global_num_end)
{
  MPI_Status status;
  int errcode, count;

  MPI_Datatype ent_type = MPI_BYTE;
  MPI_Offset disp = f->offset + ((global_num_start - 1) * size);
  cs_gnum_t gcount = (global_num_end - global_num_start)*size;

  size_t retval = 0;

  assert(gcount == 0 || f->fh != MPI_FILE_NULL);

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > INT_MAX) {
    MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
    MPI_Type_commit(&ent_type);
    count = global_num_end - global_num_start;
  }
  else
    count = gcount;

  errcode = MPI_File_write_at_all(f->fh, disp, buf, count, ent_type, &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  if (count > 0)
    MPI_Get_count(&status, ent_type, &count);

  if (ent_type != MPI_BYTE) {
    MPI_Type_free(&ent_type);
    retval = count;
  }
  else
    retval = count / size;

  return retval;
}

static size_t
_mpi_file_write_block_ip(cs_file_t  *f,
                         void       *buf,
                         size_t      size,
                         cs_gnum_t   global_num_start,
                         cs_gnum_t   global_num_end)
{
  int lengths[1];
  MPI_Aint disps[1];
  MPI_Status status;
  MPI_Datatype file_type;

  int errcode = MPI_SUCCESS, count = 0;
  char datarep[] = "native";
  MPI_Datatype ent_type = MPI_BYTE;
  cs_gnum_t gcount = (global_num_end - global_num_start) * size;
  cs_gnum_t gdisp = (global_num_start - 1) * size;

  size_t retval = 0;

  assert(gcount == 0 || f->fh != MPI_FILE_NULL);

  if (f->fh == MPI_FILE_NULL)
    return retval;

  if (gcount > INT_MAX || gdisp > INT_MAX) {
    MPI_Type_contiguous(size, MPI_BYTE, &ent_type);
    MPI_Type_commit(&ent_type);
    lengths[0] = global_num_end - global_num_start;
    disps[0] = global_num_start - 1;
  }
  else {
    lengths[0] = gcount;
    disps[0] = gdisp;
  }

  MPI_Type_create_hindexed(1, lengths, disps, ent_type, &file_type);
  MPI_Type_commit(&file_type);

  MPI_File_set_view(f->fh, f->offset, ent_type, file_type, datarep, f->info);

  errcode = MPI_File_write_all(f->fh, buf, (int)(lengths[0]), ent_type,
                               &status);

  if (errcode != MPI_SUCCESS)
    _mpi_io_error_message(f->name, errcode);

  MPI_Type_free(&file_type);

  if (lengths[0] > 0)
    MPI_Get_count(&status, ent_type, &count);

  if (ent_type != MPI_BYTE) {
    MPI_Type_free(&ent_type);
    retval = count;
  }
  else
    retval = count / size;

  return retval;
}

#endif /* defined(HAVE_MPI_IO) */

/*----------------------------------------------------------------------------
 * Compare strings (qsort function).
 *
 * parameters:
 *   a <-> pointer to first string
 *   b <-> pointer to second string
 *
 * returns:
 *   result of strcmp() on strings
 *----------------------------------------------------------------------------*/

static int
_cs_file_compare_names(const void  *a,
                       const void  *b)
{
  return strcmp(*((const char *const *)a), *((const char *const *)b));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a file descriptor and open the associated file.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by cs_file_set_swap_endian().
 *
 * \param[in]  name        file name
 * \param[in]  mode        file acces mode: read, write, or append
 * \param[in]  method      file access method
 * \param[in]  hints       associated hints for MPI-IO, or MPI_INFO_NULL
 * \param[in]  block_comm  handle to MPI communicator used for distributed file
 *                         block access (may be a subset of comm if some ranks
 *                         do not directly access distributed data blocks)
 * \param[in]  comm        handle to main MPI communicator
 *
 * \return pointer to cs_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

#else

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a file descriptor and open the associated file.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by cs_file_set_swap_endian().
 *
 * \param[in]  name        file name
 * \param[in]  mode        file access mode: read, write, or append
 * \param[in]  method      file access method (currently only C standard-IO
 *                         when built without MPI)
 * \param[in]  hints       associated hints for MPI-IO, or MPI_INFO_NULL
 * \param[in]  block_comm  handle to MPI communicator used for distributed
 *                         file block access (may be a subset of comm if some
 *                         ranks do not directly access distributed data blocks)
 * \param[in]  comm        handle to main MPI communicator
 *
 * \return pointer to cs_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

#endif

#if defined(HAVE_MPI)

cs_file_t *
cs_file_open(const char        *name,
             cs_file_mode_t     mode,
             cs_file_access_t   method,
             MPI_Info           hints,
             MPI_Comm           block_comm,
             MPI_Comm           comm)

#else

cs_file_t *
cs_file_open(const char        *name,
             cs_file_mode_t     mode,
             cs_file_access_t   method)

#endif
{
  int errcode = 0;
  cs_file_t * f = NULL;

  BFT_MALLOC(f, 1, cs_file_t);

  f->sh = NULL;

#if defined(HAVE_ZLIB)
  f->gzh = NULL;
#endif

#if defined(HAVE_MPI)
  f->comm = MPI_COMM_NULL;
  f->io_comm = MPI_COMM_NULL;
#if defined(HAVE_MPI_IO)
  f->fh = MPI_FILE_NULL;
  f->info = hints;
#endif
#endif

  f->offset = 0;

  BFT_MALLOC(f->name, strlen(name) + 1, char);
  strcpy(f->name, name);

  f->mode = mode;
  f->method = method = _access_method(method, (mode != CS_FILE_MODE_READ));

  f->rank = 0;
  f->n_ranks = 1;

  f->swap_endian = false; /* Use native endianness by default */

  /* Set communicator */

#if defined(HAVE_MPI)
  {
    int n_io_ranks = f->n_ranks;

    if (comm != MPI_COMM_NULL) {
      MPI_Comm_size(comm, &(f->n_ranks));
      if (f->n_ranks > 1) {
        f->comm = comm;
        f->io_comm = block_comm;
        MPI_Comm_rank(f->comm, &(f->rank));
        if (f->io_comm != f->comm) {
          int _n_io_ranks = 0;
          if (f->io_comm != MPI_COMM_NULL)
            MPI_Comm_size(f->io_comm, &_n_io_ranks);
          MPI_Allreduce(&_n_io_ranks, &n_io_ranks, 1, MPI_INT, MPI_MAX,
                        f->comm);
        }
      }
      else {
        f->comm = MPI_COMM_NULL;
        f->io_comm = MPI_COMM_NULL;
      }
    }

    if (n_io_ranks < 1)
      n_io_ranks = 1;
    f->rank_step = f->n_ranks / n_io_ranks;
    if (f->n_ranks % n_io_ranks)
      f->rank_step += 1;

    f->block_size = NULL;
    if (f->rank_step > 1) {
      if (f->io_comm != MPI_COMM_NULL)
        BFT_MALLOC(f->block_size, f->rank_step, cs_gnum_t);
      else
        BFT_MALLOC(f->block_size, 1, cs_gnum_t);
    }

    if (f->comm == MPI_COMM_NULL)
      f->method = CS_FILE_STDIO_SERIAL;
  }
#else
  f->method = CS_FILE_STDIO_SERIAL;
#endif

  /* Use MPI IO ? */

#if !defined(HAVE_MPI_IO)
  if (f->method > CS_FILE_STDIO_PARALLEL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error opening file:\n%s\n"
                "MPI-IO is requested, but not available."),
              name);
#endif

  /* Open file. In case of failure, destroy the allocated structure;
     this is only useful with a non-default error handler,
     as the program is terminated by default */

  if (f->method <= CS_FILE_STDIO_PARALLEL && f->rank == 0)
    errcode = _file_open(f);

#if defined(HAVE_MPI_IO)
  if (f->method == CS_FILE_MPI_INDEPENDENT) {
    f->io_comm = MPI_COMM_SELF;
    if (f->rank == 0)
      errcode = _mpi_file_open(f, f->mode);
  }
  else if (f->method > CS_FILE_MPI_INDEPENDENT)
    errcode = _mpi_file_open(f, f->mode);
#endif

  if (errcode != 0)
    f = cs_file_free(f);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a file descriptor and open the associated file, using the
 *        default file communicator and access method.
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
cs_file_open_default(const char      *name,
                     cs_file_mode_t   mode)
{
  cs_file_t *f = NULL;

  if (mode == CS_FILE_MODE_READ) {
#if defined(HAVE_MPI)
    f = cs_file_open(name,
                     mode,
                     _default_access_r,
                     _mpi_io_hints_r,
                     _mpi_io_comm,
                     cs_glob_mpi_comm);
#else
    f = cs_file_open(name,
                     mode,
                     _default_access_r);
#endif
  }
  else {
#if defined(HAVE_MPI)
    f = cs_file_open(name,
                     mode,
                     _default_access_w,
                     _mpi_io_hints_w,
                     _mpi_io_comm,
                     cs_glob_mpi_comm);
#else
    f = cs_file_open(name,
                     mode,
                     _default_access_w);
#endif
  }

  return f;
}

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
                    cs_file_mode_t   mode)
{
  cs_file_t *f = NULL;

#if defined(HAVE_MPI)
  f = cs_file_open(name,
                   mode,
                   CS_FILE_STDIO_SERIAL,
                   MPI_INFO_NULL,
                   MPI_COMM_NULL,
                   cs_glob_mpi_comm);
#else
  f = cs_file_open(name,
                   mode,
                   CS_FILE_STDIO_SERIAL);
#endif

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a file descriptor and close the associated file.
 *
 * \param[in, out]  f  file descriptor to destroy
 */
/*----------------------------------------------------------------------------*/

cs_file_t *
cs_file_free(cs_file_t  *f)
{
  cs_file_t  *_f = f;

  if (_f->sh != NULL)
    _file_close(_f);

#if defined(HAVE_MPI_IO)
  else if (_f->fh != MPI_FILE_NULL)
    _mpi_file_close(_f);
  BFT_FREE(f->block_size);
#endif

  BFT_FREE(_f->name);
  BFT_FREE(_f);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a file's name.
 *
 * \param[in]  f  cs_file_t descriptor
 *
 * \return pointer to the file's name.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_file_get_name(const cs_file_t  *f)
{
  assert(f != NULL);

  return f->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * \param[in, out]  f  cs_file_t descriptor
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_big_endian(cs_file_t  *f)
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a file's byte-swapping behavior.
 *
 * \param[in]  f  cs_file_t descriptor
 *
 * \return 0 if file's endianness is the same as the system's, 1 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_get_swap_endian(const cs_file_t  *f)
{
  assert(f != NULL);

  return f->swap_endian;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a file's byte-swapping behavior.
 *
 * \param[in, out]  f     cs_file_t descriptor
 * \param[in]       swap  1 if bytes must be swapped, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_swap_endian(cs_file_t  *f,
                        int         swap)
{
  assert(f != NULL);

  f->swap_endian = swap;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read global data from a file, distributing it to all processes
 * associated with that file.
 *
 * \param[in]  f     cs_file_t descriptor
 * \param[out] buf   pointer to location receiving data
 * \param[in]  size  size of each item of data in bytes
 * \param[in]  ni    number of items to read
 *
 * \return the number of items (not bytes) sucessfully read;
 *         currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_file_read_global(cs_file_t  *f,
                    void       *buf,
                    size_t      size,
                    size_t      ni)
{
  size_t retval = 0;

  if (f->method <= CS_FILE_STDIO_PARALLEL) {
    if (f->rank == 0) {
      if (_file_seek(f, f->offset, CS_FILE_SEEK_SET) == 0)
        retval = _file_read(f, buf, size, ni);
    }
  }

#if defined(HAVE_MPI_IO)

  else if ((f->method > CS_FILE_STDIO_PARALLEL)) {

    MPI_Status status;
    int errcode = MPI_SUCCESS, count = 0;

    if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS) {
      if (f->rank == 0) {
        errcode = MPI_File_read_at(f->fh,
                                   f->offset,
                                   buf,
                                   size*ni,
                                   MPI_BYTE,
                                   &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
      }
    }

    else {
      MPI_Datatype file_type;
      MPI_Aint disps[1];
      int lengths[1];
      char datarep[] = "native";
      lengths[0] = ni * size;
      disps[0] = 0;
      MPI_Type_create_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
      MPI_Type_commit(&file_type);
      MPI_File_set_view(f->fh, f->offset, MPI_BYTE, file_type,
                        datarep, f->info);
      if (f->rank == 0) {
        errcode = MPI_File_read(f->fh, buf, size*ni, MPI_BYTE, &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
      }
      MPI_Type_free(&file_type);
    }

    if (errcode != MPI_SUCCESS)
      _mpi_io_error_message(f->name, errcode);

    retval = count / size;

  }

#endif /* defined(HAVE_MPI_IO) */

#if defined(HAVE_MPI)
  if (f->comm != MPI_COMM_NULL) {
    long _retval = retval;
    MPI_Bcast(buf, size*ni, MPI_BYTE, 0, f->comm);
    MPI_Bcast(&_retval, 1, MPI_LONG, 0, f->comm);
    retval = _retval;
  }
#endif

  /* Update offset */

  f->offset += (cs_file_off_t)ni * (cs_file_off_t)size;

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf, buf, size, retval);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write global data to a file.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The buffers on other ranks are ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * \param[in]  f     cs_file_t descriptor
 * \param[in]  buf   pointer to location containing data
 * \param[in]  size  size of each item of data in bytes
 * \param[in]  ni    number of items to write
 *
 * \return the number of items (not bytes) sucessfully written;
 *         currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_file_write_global(cs_file_t   *f,
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
          || (f->method > CS_FILE_STDIO_PARALLEL))) {

    if (size*ni > sizeof(_copybuf))
      BFT_MALLOC(copybuf, size*ni, unsigned char);
    memcpy(copybuf, buf, size*ni);

    if (f->swap_endian == true && size > 1)
      _swap_endian(copybuf, copybuf, size, ni);

    _buf = copybuf;
  }

  if (f->rank == 0 && f->sh != NULL && f->method <= CS_FILE_STDIO_PARALLEL) {
    if (f->method == CS_FILE_STDIO_PARALLEL) {
      if (_file_seek(f, f->offset, CS_FILE_SEEK_SET) != 0)
        retval = 0;
    }
    if (retval != 0)
      retval = _file_write(f, _buf, size, ni);
  }

#if defined(HAVE_MPI_IO)

  else if ((f->method > CS_FILE_STDIO_PARALLEL)) {

    MPI_Status status;
    int errcode = MPI_SUCCESS, count = 0;

    if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS) {
      if (f->rank == 0) {
        errcode = MPI_File_write_at(f->fh,
                                    f->offset,
                                    copybuf,
                                    size*ni,
                                    MPI_BYTE,
                                    &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
      }
    }

    else {
      MPI_Datatype file_type;
      MPI_Aint disps[1];
      int lengths[1];
      char datarep[] = "native";
      lengths[0] = ni * size;
      disps[0] = 0;
      MPI_Type_create_hindexed(1, lengths, disps, MPI_BYTE, &file_type);
      MPI_Type_commit(&file_type);
      MPI_File_set_view(f->fh, f->offset, MPI_BYTE,
                        file_type, datarep, f->info);
      if (f->rank == 0) {
        errcode = MPI_File_write(f->fh,
                                 copybuf,
                                 size*ni,
                                 MPI_BYTE,
                                 &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
      }
      MPI_Type_free(&file_type);
    }

    if (errcode != MPI_SUCCESS)
      _mpi_io_error_message(f->name, errcode);

    retval = count / size;

  }

#endif /* defined(HAVE_MPI_IO) */

  if (copybuf != _copybuf) /* Free allocated memory if necessary */
    BFT_FREE(copybuf);

#if defined(HAVE_MPI)
  if (f->comm != MPI_COMM_NULL) {
    long _retval = retval;
    MPI_Bcast(&_retval, 1, MPI_LONG, 0, f->comm);
    retval = _retval;
  }
#endif

  /* Update offset */

  f->offset += (cs_file_off_t)ni * (cs_file_off_t)size;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read data to a buffer, distributing a contiguous part of it to each
 * process associated with a file.
 *
 * Each process should receive a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * \param[in]  f                 cs_file_t descriptor
 * \param[out] buf               pointer to location receiving data
 * \param[in]  size              size of each item of data in bytes
 * \param[in]  stride            number of (interlaced) values per block item
 * \param[in]  global_num_start  global number of first block item
 *                               (1 to n numbering)
 * \param[in]  global_num_end    global number of past-the end block item
 *                               (1 to n numbering)
 *
 * \return the (local) number of items (not bytes) sucessfully read;
 *         currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_file_read_block(cs_file_t  *f,
                   void       *buf,
                   size_t      size,
                   size_t      stride,
                   cs_gnum_t   global_num_start,
                   cs_gnum_t   global_num_end)
{
  size_t retval = 0;

  cs_gnum_t global_num_end_last = global_num_end;

  cs_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
  cs_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

  if (_global_num_end < _global_num_start)
    _global_num_end = _global_num_start;

  void *_buf = buf;

#if defined(HAVE_MPI)
  if (f->rank_step > 1)
    _buf = _gather_block_sizes(f,
                               size,
                               _global_num_start,
                               &_global_num_end);
#endif

  assert(global_num_end >= global_num_start);

  switch(f->method) {

  case CS_FILE_STDIO_SERIAL:
    retval = _file_read_block_s(f,
                                _buf,
                                size,
                                _global_num_start,
                                _global_num_end);
    break;

  case CS_FILE_STDIO_PARALLEL:
    retval = _file_read_block_p(f,
                                _buf,
                                size,
                                _global_num_start,
                                _global_num_end);
    break;

#if defined(HAVE_MPI_IO)

  case CS_FILE_MPI_INDEPENDENT:
  case CS_FILE_MPI_NON_COLLECTIVE:
    retval = _mpi_file_read_block_noncoll(f,
                                          _buf,
                                          size,
                                          _global_num_start,
                                          _global_num_end);
    break;

  case CS_FILE_MPI_COLLECTIVE:

    if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS)
      retval = _mpi_file_read_block_eo(f,
                                       _buf,
                                       size,
                                       _global_num_start,
                                       _global_num_end);
    else
      retval = _mpi_file_read_block_ip(f,
                                       _buf,
                                       size,
                                       _global_num_start,
                                       _global_num_end);
    break;

#endif /* defined(HAVE_MPI_IO) */

  default:
    assert(0);
  }

  /* Update offset */

  assert(f->rank > 0 || global_num_start == 1);

#if defined(HAVE_MPI)
  if (f->n_ranks > 1)
    MPI_Bcast(&global_num_end_last, 1, CS_MPI_GNUM, f->n_ranks-1, f->comm);
#endif

  f->offset += ((global_num_end_last - 1) * size * stride);

#if defined(HAVE_MPI)
  if (f->rank_step > 1) {
    retval = _scatter_blocks(f, _buf, buf, size);
    if (_buf != buf)
      BFT_FREE(_buf);
  }
#endif

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf, buf, size, retval);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write data to a file, each associated process providing a
 * contiguous part of this data.
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
 * \param[in]  f                 cs_file_t descriptor
 * \param[in]  buf               pointer to location containing data
 * \param[in]  size              size of each item of data in bytes
 * \param[in]  stride            number of (interlaced) values per block item
 * \param[in]  global_num_start  global number of first block item
 *                               (1 to n numbering)
 * \param[in]  global_num_end    global number of past-the end block item
 *                               (1 to n numbering)
 *
 * \return the (local) number of items (not bytes) sucessfully written;
 *         currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_file_write_block(cs_file_t   *f,
                    const void  *buf,
                    size_t       size,
                    size_t       stride,
                    cs_gnum_t    global_num_start,
                    cs_gnum_t    global_num_end)
{
  size_t retval = 0;

  const size_t bufsize = (global_num_end - global_num_start)*stride*size;

  /* Copy contents to ensure buffer constedness if necessary */

  bool direct_w = true;

  if (f->swap_endian == true && size > 1)
    direct_w = false;

#if defined(HAVE_MPI)
  if (f->n_ranks > 1) {
    if (f->rank_step > 1 || f->method != CS_FILE_STDIO_PARALLEL)
      direct_w = false;
  }
#endif

  if (direct_w == false) {

    unsigned char *copybuf = NULL;

    BFT_MALLOC(copybuf, bufsize, unsigned char);

    if (copybuf != NULL)
      memcpy(copybuf, buf, bufsize);

    retval = cs_file_write_block_buffer(f,
                                        copybuf,
                                        size,
                                        stride,
                                        global_num_start,
                                        global_num_end);

    BFT_FREE(copybuf);
  }

  /* Using Standard IO with no byte-swapping or serialization, write directly */

  else {

    cs_gnum_t global_num_end_last = global_num_end;

    const cs_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
    const cs_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

    if (_global_num_end > _global_num_start) {

      if (f->sh == NULL)
        _file_open(f);

      retval = _file_write(f,
                           buf,
                           size,
                           (_global_num_end - _global_num_start));

    }

    /* Update offset */

#if defined(HAVE_MPI)
    if (f->n_ranks > 1)
      MPI_Bcast(&global_num_end_last, 1, CS_MPI_GNUM, f->n_ranks-1, f->comm);
#endif

    f->offset += ((global_num_end_last - 1) * size * stride);

  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write data to a file, each associated process providing a
 * contiguous part of this data.
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
 * \param[in]  f                 cs_file_t descriptor
 * \param[in, out]  buf          pointer to location containing data
 * \param[in]  size              size of each item of data in bytes
 * \param[in]  stride            number of (interlaced) values per block item
 * \param[in]  global_num_start  global number of first block item
 *                               (1 to n numbering)
 * \param[in]  global_num_end    global number of past-the end block item
 *                               (1 to n numbering)
 *
 * \return the (local) number of items (not bytes) sucessfully written;
 *         currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_file_write_block_buffer(cs_file_t  *f,
                           void       *buf,
                           size_t      size,
                           size_t      stride,
                           cs_gnum_t   global_num_start,
                           cs_gnum_t   global_num_end)
{
  size_t retval = 0;

  cs_gnum_t global_num_end_last = global_num_end;

  cs_gnum_t _global_num_start = (global_num_start-1)*stride + 1;
  cs_gnum_t _global_num_end = (global_num_end-1)*stride + 1;

  if (_global_num_end < _global_num_start)
    _global_num_end = _global_num_start;

  void *_buf = buf;

  /* Swap bytes prior to writing if necessary */

  if (f->swap_endian == true && size > 1)
    _swap_endian(buf,
                 buf,
                 size,
                 (_global_num_end - _global_num_start));

#if defined(HAVE_MPI)
   if (f->rank_step > 1)
    _buf = _gather_blocks(f, buf, size, _global_num_start,
                          &_global_num_end);
#endif

  /* Write to file using chosen method */

  switch(f->method) {

  case CS_FILE_STDIO_SERIAL:
    retval = _file_write_block_s(f,
                                 _buf,
                                 size,
                                 _global_num_start,
                                 _global_num_end);
    break;

  case CS_FILE_STDIO_PARALLEL:
    retval = _file_write_block_p(f,
                                 _buf,
                                 size,
                                 _global_num_start,
                                 _global_num_end);
    break;

#if defined(HAVE_MPI_IO)

  case CS_FILE_MPI_INDEPENDENT:
  case CS_FILE_MPI_NON_COLLECTIVE:
      retval = _mpi_file_write_block_noncoll(f,
                                             _buf,
                                             size,
                                             _global_num_start,
                                             _global_num_end);
      break;

  case CS_FILE_MPI_COLLECTIVE:
    if (_mpi_io_positioning == CS_FILE_MPI_EXPLICIT_OFFSETS)
      retval = _mpi_file_write_block_eo(f,
                                        _buf,
                                        size,
                                        _global_num_start,
                                        _global_num_end);
    else
      retval = _mpi_file_write_block_ip(f,
                                        _buf,
                                        size,
                                        _global_num_start,
                                        _global_num_end);
    break;

#endif /* defined(HAVE_MPI_IO) */

  default:
    assert(0);
  }

#if defined(HAVE_MPI)
  if (f->rank_step > 1) {
    if (f->rank % f->rank_step == 0) {
      /* Check for inconsistent read sizes */
      int rank_end = f->rank + f->rank_step;
      if (rank_end >= f->n_ranks)
        rank_end = f->n_ranks;
      int n_aggr = rank_end - f->rank;
      cs_gnum_t retval_cmp = 0;
      for (int i = 0; i < n_aggr; i++)
        retval_cmp += f->block_size[i];
      if (retval_cmp != retval)  /* Error in this case */
        f->block_size[0] = retval_cmp;
    }
    retval = f->block_size[0];
    if (_buf != buf)
      BFT_FREE(_buf);
  }

  /* Update offset */

  if (f->n_ranks > 1)
    MPI_Bcast(&global_num_end_last, 1, CS_MPI_GNUM, f->n_ranks-1, f->comm);
#endif

  f->offset += ((global_num_end_last - 1) * size * stride);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the file pointer according to whence.
 *
 * \param[in, out]  f       cs_file_t descriptor
 * \param[in]       offset  add to position specified to whence to obtain
 *                          new position, measured in characters from the
 *                          beginning of the file
 * \param[in]       whence  beginning if CS_FILE_SEEK_SET,
 *                          current if CS_FILE_SEEK_CUR,
 *                          or end-of-file if CS_FILE_SEEK_END
 *
 * \return 0 upon success, nonzero otherwise; currently, errors are fatal.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_seek(cs_file_t       *f,
             cs_file_off_t    offset,
             cs_file_seek_t   whence)
{
  int retval = 0;

  /* Always update f->offset, regardless of mode */

  switch(whence) {

  case CS_FILE_SEEK_SET:

    f->offset = offset;
    break;

  case CS_FILE_SEEK_CUR:

    f->offset += offset;
    break;

  case CS_FILE_SEEK_END:

    if (f->sh != NULL)
      f->offset = cs_file_tell(f) + offset;

#if defined(HAVE_MPI_IO)
    if (f->fh != MPI_FILE_NULL) {
      MPI_Offset f_size = 0;
      retval = MPI_File_get_size(f->fh, &f_size);
      f->offset = f_size + offset;
    }
#endif

#if defined(HAVE_MPI)
  if (f->comm != MPI_COMM_NULL) {
#if defined(MPI_LONG_LONG)
    long long offset_g;
    long long offset_l = f->offset;
    MPI_Datatype  _mpi_datatype_offset = MPI_LONG_LONG;
#else
    long offset_g;
    long offset_l = f->offset;
    MPI_Datatype  _mpi_datatype_offset = MPI_LONG_INT;
#endif
    MPI_Allreduce(&offset_l, &offset_g, 1, _mpi_datatype_offset, MPI_MAX,
                  f->comm);
    f->offset = offset_g;
  }
#endif

  break;
  }

  /* Now update actual file position */

  if (f->sh != NULL)
      retval = _file_seek(f, offset, whence);

#if defined(HAVE_MPI_IO)

  else if (   f->fh != MPI_FILE_NULL
           && _mpi_io_positioning == CS_FILE_MPI_INDIVIDUAL_POINTERS) {

    retval = MPI_File_seek(f->fh, f->offset, MPI_SEEK_SET);

    if (retval != MPI_SUCCESS)
      _mpi_io_error_message(f->name, retval);

  }

#endif /* defined(HAVE_MPI_IO) */

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the position of the file pointer.
 *
 * In parallel, we consider the file pointer to be equal to the highest
 * value of the individual file pointers.
 *
 * \param[in]  f  cs_file_t descriptor
 *
 * \return current position of the file pointer.
 */
/*----------------------------------------------------------------------------*/

cs_file_off_t
cs_file_tell(cs_file_t  *f)
{
  cs_file_off_t retval = f->offset;

  if (f->method == CS_FILE_STDIO_SERIAL && f->rank == 0 && f->sh != NULL)
    retval = _file_tell(f);

#if defined(HAVE_MPI)
  if (f->comm != MPI_COMM_NULL) {
#if defined(MPI_LONG_LONG)
    long long _offset = retval;
    MPI_Datatype  _mpi_datatype_offset = MPI_LONG_LONG;
#else
    long _offset = retval;
    MPI_Datatype  _mpi_datatype_offset = MPI_LONG_INT;
#endif
    MPI_Bcast(&_offset, 1, _mpi_datatype_offset, 0, f->comm);
    retval = _offset;
  }
#endif

  /*
    Note that in case of individual file pointers, using
    MPI_File_get_position() and MPI_File_get_byte_offset() should also
    work, but fail after certain collective writes with some processes
    writing zero values (at least on Open MPI 1.2.6), so we prefer to keep
    track of the global offset (which we need for seeking or views anyways).
  */

  return retval;
}

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
             int              *line)
{
  return _cs_file_gets(s, size, f, line, 0);
}

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
cs_file_gets_try(char             *s,
                 const int         size,
                 const cs_file_t  *f,
                 int              *line)
{
  return _cs_file_gets(s, size, f, line, 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump the metadata of a file structure in human readable form.
 *
 * \param[in]  f  cs_file_t descriptor
 */
/*----------------------------------------------------------------------------*/

void
cs_file_dump(const cs_file_t  *f)
{
  const char *mode_name[] = {"CS_FILE_MODE_READ",
                             "CS_FILE_MODE_WRITE",
                             "CS_FILE_MODE_APPEND"};
  const char *access_name[] = {"CS_FILE_STDIO_SERIAL",
                               "CS_FILE_STDIO_PARALLEL",
                               "CS_FILE_MPI_INDEPENDENT",
                               "CS_FILE_MPI_NON_COLLECTIVE",
                               "CS_FILE_MPI_COLLECTIVE"};

  if (f == NULL) {
    bft_printf("\n"
               "Null file dump:\n");
    return;
  }

#if defined(HAVE_MPI)
  bft_printf("\n"
             "File name:                   \"%s\"\n"
             "Access mode:                 %s\n"
             "Access method:               %s\n"
             "Rank:                        %d\n"
             "N ranks:                     %d\n"
             "rank step:                   %d\n"
             "Swap endian:                 %d\n"
             "Serial handle:               %p\n",
             f->name, mode_name[f->mode], access_name[f->method-1],
             f->rank, f->n_ranks, f->rank_step, (int)(f->swap_endian),
             (const void *)f->sh);
#else
  bft_printf("\n"
             "File name:                   \"%s\"\n"
             "Access mode:                 %s\n"
             "Access method:               %s\n"
             "Rank:                        %d\n"
             "N ranks:                     %d\n"
             "Swap endian:                 %d\n"
             "Serial handle:               %p\n",
             f->name, mode_name[f->mode], access_name[f->method-1],
             f->rank, f->n_ranks, (int)(f->swap_endian), (const void *)f->sh);
#endif

#if defined(HAVE_MPI)
  bft_printf("Associated io communicator:  %llu\n",
             (unsigned long long)(f->io_comm));
  bft_printf("Associated communicator:     %llu\n",
             (unsigned long long)(f->comm));
#if defined(HAVE_MPI_IO)
  bft_printf("MPI file handle:             %llu\n"
             "MPI file offset:             %llu\n",
             (unsigned long long)(f->fh),
             (unsigned long long)(f->offset));
#endif
#endif

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the default options for file access.
 */
/*----------------------------------------------------------------------------*/

void
cs_file_free_defaults(void)
{
  _mpi_io_positioning = CS_FILE_MPI_EXPLICIT_OFFSETS;

  _default_access_r = CS_FILE_DEFAULT;
  _default_access_w = CS_FILE_DEFAULT;

  /* Communicator and hints used for file operations */

#if defined(HAVE_MPI)
  _mpi_defaults_are_set = false;
  _mpi_rank_step = 1;
  _mpi_comm = MPI_COMM_NULL;

  if (_mpi_io_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&_mpi_io_comm);
    _mpi_io_comm = MPI_COMM_NULL;
  }
#endif

#if defined(HAVE_MPI_IO)
#  if MPI_VERSION > 1
  if (_mpi_io_hints_r != MPI_INFO_NULL)
    MPI_Info_free(&_mpi_io_hints_r);
  if (_mpi_io_hints_w != MPI_INFO_NULL)
    MPI_Info_free(&_mpi_io_hints_w);
#  endif /* MPI_VERSION > 1 */
#endif /* defined(HAVE_MPI_IO) */
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the default options for file access.
 *
 * \param[in]    mode    file mode for which the default is queried
 *                       (write and append use the same method, and are
 *                       interchangeable here)
 * \param[out]   method  default file access method, or NULL
 * \param[out]   hints   MPI-IO hints, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_file_get_default_access(cs_file_mode_t     mode,
                           cs_file_access_t  *method,
                           MPI_Info          *hints)
{
  if (mode == CS_FILE_MODE_READ) {
    if (method != NULL)
      *method = _access_method(_default_access_r, false);
    if (hints != NULL)
      *hints = _mpi_io_hints_r;
  }
  else {
    if (method != NULL)
      *method = _access_method(_default_access_w, true);
    if (hints != NULL)
      *hints = _mpi_io_hints_w;
  }
}

#else /* if !defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the default options for file access.
 *
 * \param[in]    mode    file mode for which the default is queried
 *                       (write and append use the same method, and are
 *                       interchangeable here)
 * \param[out]   method  default file access method, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_file_get_default_access(cs_file_mode_t     mode,
                           cs_file_access_t  *method)
{
  if (mode == CS_FILE_MODE_READ) {
    if (method != NULL)
      *method = _access_method(_default_access_r, false);
  }
  else {
    if (method != NULL)
      *method = _access_method(_default_access_w, true);
  }
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the default options for file access.
 *
 * If the method given contains incompatible values, such as when setting
 * MPI-IO methods when MPI-IO is not available, a "reasonable" default
 * is used instead.
 *
 * \param[in]  mode       file mode for which the default is being set
 *                        (write and append use the same method, and are
 *                        interchangeable here)
 * \param[in]  method     default access method to set
 * \param[in]  hints      MPI-IO hints, or MPI_INFO_NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_default_access(cs_file_mode_t    mode,
                           cs_file_access_t  method,
                           MPI_Info          hints)
{
  cs_file_access_t  _method;

  if (mode == CS_FILE_MODE_READ) {
    _method = _access_method(method, false);
    _default_access_r = _method;
  }
  else { /* if (mode == CS_FILE_MODE_WRITE || mode == CS_FILE_MODE_APPEND) */
    _method = _access_method(method, true);
    _default_access_w = _method;
  }

#if defined(HAVE_MPI_IO)
#  if MPI_VERSION > 1

  /* Free previous info objects */

  if (mode == CS_FILE_MODE_READ && _mpi_io_hints_r != MPI_INFO_NULL)
    MPI_Info_free(&_mpi_io_hints_r);
  else if (    (mode == CS_FILE_MODE_WRITE || mode == CS_FILE_MODE_APPEND)
           && _mpi_io_hints_w != MPI_INFO_NULL)
    MPI_Info_free(&_mpi_io_hints_w);

  /* Set info objects */

  if (_method > CS_FILE_STDIO_PARALLEL && hints != MPI_INFO_NULL) {
    if (mode == CS_FILE_MODE_READ)
      MPI_Info_dup(hints, &_mpi_io_hints_r);
    else if (mode == CS_FILE_MODE_WRITE || mode == CS_FILE_MODE_APPEND)
      MPI_Info_dup(hints, &_mpi_io_hints_w);
  }

#  endif /* MPI_VERSION > 1 */
#endif /* defined(HAVE_MPI_IO) */
}

#else /* if !defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the default options for file access.
 *
 * If the method given contains incompatible values, such as when setting
 * MPI-IO methods when MPI-IO is not available, a "reasonable" default
 * is used instead.
 *
 * \param[in]  mode       file mode for which the default is being set
 *                        (write and append use the same method, and are
 *                        interchangeable here)
 * \param[in]  method     default access method to set
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_default_access(cs_file_mode_t    mode,
                           cs_file_access_t  method)
{
  if (mode == CS_FILE_MODE_READ)
    _default_access_r = _access_method(method, false);
  else if (mode == CS_FILE_MODE_WRITE || mode == CS_FILE_MODE_APPEND)
    _default_access_w = _access_method(method, true);
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get default MPI communicator values for file access.
 *
 * A block rank stepping value may be used, allowing the use of a reduced
 * communicator for distributed block reads and writes.
 * If this value is greater than 1, ranks not a multiple of this step must be
 * guaranteed to be empty for block reads and writes with files opened using
 * this default.
 *
 * \param[out]   block_rank_step  MPI rank stepping between non-empty
 *                                distributed blocks, or NULL
 * \param[out]   block_comm       Handle to MPI communicator used for
 *                                distributed file block access, or NULL
 * \param[out]   comm             Handle to main MPI communicator, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_file_get_default_comm(int       *block_rank_step,
                         MPI_Comm  *block_comm,
                         MPI_Comm  *comm)
{
  /* Initialize defauts if not already done */

  if (_mpi_defaults_are_set == false && cs_glob_mpi_comm != MPI_COMM_NULL) {
    cs_file_set_default_comm(0, MPI_COMM_SELF);
    _mpi_defaults_are_set = true;
  }

  /* Return defaults */

  if (block_rank_step != NULL)
    *block_rank_step = _mpi_rank_step;

  if (block_comm != NULL) {
    if (_mpi_comm != MPI_COMM_NULL)
      *block_comm = _mpi_io_comm;
    else
      *block_comm = cs_glob_mpi_comm;
  }

  if (comm != NULL) {
    if (_mpi_comm != MPI_COMM_NULL)
      *comm = _mpi_comm;
    else
      *comm = cs_glob_mpi_comm;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set default MPI communicator values for file access.
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
 * \param[in]  block_rank_step  MPI rank stepping between non-empty blocks for
 *                              file block reads and writes (not set if <= 0)
 * \param[in]  comm             Handle to main MPI communicator
 *                              (not set if MPI_COMM_SELF)
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_default_comm(int       block_rank_step,
                         MPI_Comm  comm)
{
  if (block_rank_step > 0) {
    if (block_rank_step > cs_glob_n_ranks)
      block_rank_step = cs_glob_n_ranks;
    _mpi_rank_step = block_rank_step;
  }

  if (comm != MPI_COMM_SELF)
    _mpi_comm = comm;
  else if (_mpi_defaults_are_set == false)
    _mpi_comm = cs_glob_mpi_comm;

  if (   comm != MPI_COMM_SELF
      || block_rank_step > 0
      || _mpi_defaults_are_set == false) {

    if (_mpi_io_comm != MPI_COMM_NULL) {
      MPI_Comm_free(&_mpi_io_comm);
      _mpi_io_comm = MPI_COMM_NULL;
    }

    if (_mpi_comm != MPI_COMM_NULL) {

      if (_mpi_rank_step < 2) {
        _mpi_rank_step = 1;
        MPI_Comm_dup(_mpi_comm, &_mpi_io_comm);
      }

      else /* Create reduced communicator */
        _mpi_io_comm = cs_file_block_comm(_mpi_rank_step, _mpi_comm);

    }

  }

  _mpi_defaults_are_set = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an MPI communicator for distributed block parallel IO.
 *
 * \param[in]  block_rank_step  MPI rank stepping between non-empty blocks
 * \param[in]  comm             Handle to main MPI communicator
 *
 * \return communicator associated with IO, MPI_COMM_NULL for ranks not
 *         participating in parallel IO (including ranks participating in IO
 *         where communicator size would be 1)
 */
/*----------------------------------------------------------------------------*/

MPI_Comm
cs_file_block_comm(int       block_rank_step,
                   MPI_Comm  comm)
{
  MPI_Comm  new_comm = MPI_COMM_NULL;

  if (comm == MPI_COMM_NULL)
    return new_comm;

  int rank_id, n_ranks;
  MPI_Comm_rank(comm, &rank_id);
  MPI_Comm_size(comm, &n_ranks);
  if (n_ranks < 2) {
    new_comm = MPI_COMM_NULL;
    return new_comm;
  }

  if (block_rank_step > n_ranks)
    block_rank_step = n_ranks;

  if (block_rank_step < 2)
    MPI_Comm_dup(comm, &new_comm);

  /* Create reduced communicator in more general case */

  else {

    int ranges[1][3];
    MPI_Group old_group, new_group;

    MPI_Comm_group(comm, &old_group);

    MPI_Barrier(comm); /* For debugging */

    ranges[0][0] = 0;
    ranges[0][1] = n_ranks - 1;
    ranges[0][2] = block_rank_step;

    MPI_Group_range_incl(old_group, 1, ranges, &new_group);
    MPI_Comm_create(comm, new_group, &new_comm);
    MPI_Group_free(&new_group);

    MPI_Group_free(&old_group);

    if (rank_id % block_rank_step)
      new_comm = MPI_COMM_NULL;

    MPI_Barrier(comm); /* For debugging */

  }

  return new_comm;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the positioning method for MPI-IO
 *
 * For details, see \ref cs_file_set_mpi_io_positioning.
 *
 * \return  positioning method for MPI-IO
 */
/*----------------------------------------------------------------------------*/

cs_file_mpi_positioning_t
cs_file_get_mpi_io_positioning(void)
{
  return _mpi_io_positioning;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the positioning method for MPI-IO
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
 * \param[in]  positioning  chosen positioning method for MPI-IO
 */
/*----------------------------------------------------------------------------*/

void
cs_file_set_mpi_io_positioning(cs_file_mpi_positioning_t  positioning)
{
  _mpi_io_positioning = positioning;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on default options for file access.
 */
/*----------------------------------------------------------------------------*/

void
cs_file_defaults_info(void)
{
#if defined(HAVE_MPI)

  int             log_id;
  cs_log_t logs[] = {CS_LOG_DEFAULT, CS_LOG_PERFORMANCE};

  const char *fmt[4] = {N_("  I/O read method:     %s\n"),
                        N_("  I/O write method:    %s\n"),
                        N_("  I/O read method:     %s (%s)\n"),
                        N_("  I/O write method:    %s (%s)\n")};

  for (log_id = 0; log_id < 2; log_id++)
    cs_log_printf(logs[log_id], "\n");

  for (cs_file_mode_t mode = CS_FILE_MODE_READ;
       mode < CS_FILE_MODE_APPEND;
       mode++) {

    MPI_Info hints;
    cs_file_access_t method;

    cs_file_get_default_access(mode, &method, &hints);

#if defined(HAVE_MPI_IO)
    if (method > CS_FILE_STDIO_PARALLEL) {
      for (log_id = 0; log_id < 2; log_id++)
        cs_log_printf(logs[log_id],
                      _(fmt[mode + 2]),
                      _(cs_file_access_name[method]),
                      _(cs_file_mpi_positioning_name[_mpi_io_positioning]));
    }
#endif
    if (method <= CS_FILE_STDIO_PARALLEL) {
      for (log_id = 0; log_id < 2; log_id++)
        cs_log_printf(logs[log_id],
                      _(fmt[mode]), _(cs_file_access_name[method]));
    }

#if MPI_VERSION > 1

    if (hints != MPI_INFO_NULL) {
      int i, n_keys, flag;
      char *val;
      char key[MPI_MAX_INFO_KEY + 1];
      BFT_MALLOC(val, MPI_MAX_INFO_VAL + 1, char);
      MPI_Info_get_nkeys(hints, &n_keys);
      if (n_keys > 0)
        bft_printf(_("    hints:\n"));
      for (i = 0; i < n_keys; i++) {
        MPI_Info_get_nthkey(hints, i, key);
        MPI_Info_get(hints, key, MPI_MAX_INFO_VAL, val, &flag);
        if (flag) {
          val[MPI_MAX_INFO_VAL] = '\0';
          for (log_id = 0; log_id < 2; log_id++)
            cs_log_printf(logs[log_id],
                          _("      %s: %s\n"), key, val);
        }
      }
      BFT_FREE(val);
    }

#endif /* MPI_VERSION > 1 */

  }

  if (cs_glob_n_ranks > 1) {
    int block_rank_step;
    cs_file_get_default_comm(&block_rank_step, NULL, NULL);
    for (log_id = 0; log_id < 2; log_id++)
      cs_log_printf(logs[log_id],
                    _("  I/O rank step:        %d\n"), block_rank_step);
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

#endif
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a cs_file_serializer_t structure.
 *
 * The buf_block_size argument is optional, and may be used when the buffer
 * on rank 0 is larger than (global_num_end - global_num_start)*size*stride
 * bytes. If zero, a block size of (global_num_end - global_num_start) on
 * rank 0 is assumed; a buffer may not be smaller than this, as it must
 * initially contain all data on rank 0's block.
 *
 * \param[in]  size              size of each item of data in bytes
 * \param[in]  stride            number of (interlaced) values per block item
 * \param[in]  global_num_start  global number of first block item
 *                               (1 to n numbering)
 * \param[in]  global_num_end    global number of past-the end block item
 *                               (1 to n numbering)
 * \param[in]  buf_block_size    Local data buffer block size, or 0 for
 *                               default global_num_end - global_num_start
 *                               (only useful on rank 0)
 * \param[in]  buf               pointer to local block data buffer
 * \param[in]  comm              associated MPI communicator
 *
 * \return pointer to new serializer structure.
 */
/*----------------------------------------------------------------------------*/

cs_file_serializer_t *
cs_file_serializer_create(size_t       size,
                          size_t       stride,
                          cs_gnum_t    global_num_start,
                          cs_gnum_t    global_num_end,
                          size_t       buf_block_size,
                          void        *buf,
                          MPI_Comm     comm)
{
  cs_file_serializer_t  *s = NULL;

  BFT_MALLOC(s, 1, cs_file_serializer_t);

  _serializer_init(s,
                   size * stride,
                   global_num_start,
                   global_num_end,
                   buf_block_size,
                   buf,
                   comm);

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_file_serializer_t structure.
 *
 * \param[in, out]  s  pointer to pointer structure that should be destroyed
 */
/*----------------------------------------------------------------------------*/

void
cs_file_serializer_destroy(cs_file_serializer_t  **s)
{
  if (s != NULL) {
    _serializer_finalize(*s);
    BFT_FREE(*s);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advance a cs_file_serializer_t structure.
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
 * \param[in]   s          pointer to serializer structure
 * \param[out]  cur_range  optional start and past-the end global numbers for
 *                         the current block (size: 2), or NULL; only on rank 0
 *
 * \return a pointer to the buffer containing new data (first call counts as
 *          new), or NULL if we are finished; always NULL on ranks > 0.
 */
/*----------------------------------------------------------------------------*/

void *
cs_file_serializer_advance(cs_file_serializer_t  *s,
                           cs_gnum_t              cur_range[2])
{
  MPI_Status status;
  cs_gnum_t sync_range[2] = {s->next_g_num, 0};

  void *retval = NULL;

  /* Rank 0 receives data */

  if (s->rank_id == 0) {

    int count = 0;

    while (count == 0) {

      int dist_rank = s->next_rank_id;

      count = 0;

      if (s->next_rank_id >= s->n_ranks)
        return NULL;

      else if (s->next_rank_id != 0) {

        count = s->count[dist_rank];

        /* Forced synchronization */
        sync_range[1] = sync_range[0] + count;
        MPI_Send(&sync_range, 2, CS_MPI_GNUM, dist_rank, CS_FILE_MPI_TAG, s->comm);

        /* Receive data */
        MPI_Recv(s->recv_buf, (count * s->size), MPI_BYTE, dist_rank,
                 CS_FILE_MPI_TAG, s->comm, &status);

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

    };

  }

  /* Other ranks send data */

  else {

    int count = s->range[1] - s->range[0];

    if (count > 0) {

      assert(s->rank_id > -1);

      /* Forced synchronization */
      MPI_Recv(&sync_range, 2, CS_MPI_GNUM, 0, CS_FILE_MPI_TAG, s->comm, &status);
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
      MPI_Send(s->buf, (count * s->size), MPI_BYTE, 0, CS_FILE_MPI_TAG, s->comm);

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  MPI_Barrier(comm);
#endif

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new directory using default permissions.
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
 * \param[in]  path  name of new directory.
 *
 * \return 0 on success, -1 if an error occured (in which case errno
 *         contains the appropriate error code). If the underlying
 *         system has no mkdir() function or it was not detected
 *         upon BFT configuration, 1 is returned.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_mkdir_default(const char  *path)
{
  static const char  *str_fail = N_("Failure to create "
                                    "directory \"%s\":\n\n%s");

#if defined(HAVE_MKDIR)

#if defined(WIN32) || defined(_WIN32)

  mkdir(path);
  return 0;

#else

  if (mkdir(path, S_IRWXU|S_IRWXG|S_IRWXO) != 0) {

    if (errno == EEXIST) {

#if defined(HAVE_SYS_STAT_H)

      struct stat buf;

      if (stat(path, &buf) != 0)
        bft_error(__FILE__, __LINE__, 0, _(str_fail),
                  path,
                  _("  A similarly named file or directory exists "
                    "and its status is\n  not available."));
      else if (S_ISDIR(buf.st_mode) != 1)
        bft_error(__FILE__, __LINE__, 0, _(str_fail),
                  path,
                  _("  A similarly named file exists and is "
                    "not a directory."));
      else
        return 0;

#endif

      errno = EEXIST; /* In case modified by stat() */

    }
    else {
      bft_error(__FILE__, __LINE__, errno, _(str_fail),
                path,
                _("  A similarly named file exists and is "
                  "not a directory."));

    }

    return -1;

  } /* End of directory creation failure case */

#endif

  return 0;

#else /* #if defined(HAVE_MKDIR) */

  return 1;

#endif /* #if defined(HAVE_MKDIR) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a file exists and is a regular file.
 *
 * \param[in]  path  file path.
 *
 * \return 1 if file exists and is a regular file, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_isreg(const char  *path)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(path, &s) != 0) {
    if (errno != ENOENT)
      bft_error(__FILE__, __LINE__, errno,
                _("Error querying information for file:\n%s."),
                path);
  }
  else {
    if (S_ISREG(s.st_mode) != 0)
      retval = 1;
  }

#else /* defined(HAVE_SYS_STAT_H) */

  /* If Posix-type API is not available, revert to basic method */

  FILE *f;

  if ((f = fopen(fic_path, "r")) != NULL) {
    retval = 1;
    fclose(f);
  }

#endif /* defined(HAVE_SYS_STAT_H) */

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a directory exists.
 *
 * \param[in]  path  directory path.
 *
 * \return 1 if directory exists, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_isdir(const char  *path)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(path, &s) != 0) {
    if (errno != ENOENT)
      bft_error(__FILE__, __LINE__, errno,
                _("Error querying information for directory:\n%s."),
                path);
  }
  else {
    if (S_ISDIR(s.st_mode) != 0)
      retval = 1;
  }

#else /* defined(HAVE_SYS_STAT_H) */

  /* If Posix-type API is not available,
     consider that directories are not available either */

  retval = 0;

#endif /* defined(HAVE_SYS_STAT_H) */

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief List files inside a directory.
 *
 * The array returned must be freed by the caller using BFT_FREE,
 * as well as the individual entries in the array.
 *
 * \param[in]  path name of directory.
 *
 * \return an array of file names in a directory. The last entry is
 *         set to NULL. If no means to list the directory or an error
 *         occured, the return value is simply NULL.
 */
/*----------------------------------------------------------------------------*/

char **
cs_file_listdir(const char  *path)
{
  char **dirnames = NULL;

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_DIRENT_H)

  struct dirent *ent;
  int n_ent = 0;
  DIR *d = opendir(path);

  if (d == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error opening directory \"%s\":\n\n"
                "  %s"), path, strerror(errno));
    return NULL;
  }

  /* Counting pass */

  while(readdir(d) != NULL)
    n_ent += 1;

  rewinddir(d);

  BFT_MALLOC(dirnames, n_ent + 1, char *);

  n_ent = 0;
  while((ent = readdir(d)) != NULL) {
    BFT_MALLOC(dirnames[n_ent], strlen(ent->d_name) + 1, char);
    strcpy(dirnames[n_ent], ent->d_name);
    n_ent += 1;
  }
  dirnames[n_ent] = NULL;

  closedir(d);

  qsort(dirnames, n_ent, sizeof(char *), &_cs_file_compare_names);

#endif /* defined(HAVE_SYS_TYPES_H) && defined(HAVE_DIRENT_H) */

  return dirnames;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the size of a file.
 *
 * If the file does not exist, 0 is returned.
 *
 * Note also that for some special files, such as files in the Linux /proc
 * directory, this may return 0.
 *
 * \param[in]  path  file path.
 *
 * \return size of file.
 */
/*----------------------------------------------------------------------------*/

cs_file_off_t
cs_file_size(const char  *path)
{
  cs_file_off_t retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(path, &s) != 0) {
    if (errno != ENOENT)
      bft_error(__FILE__, __LINE__, errno,
                _("Error querying information for file:\n%s."),
                path);
  }
  else
    retval = s.st_size;

#else /* defined(HAVE_SYS_STAT_H) */

  /* If Posix-type API is not available, revert to basic method */

  FILE *f;

  if ((f = fopen(fic_path, "r")) != NULL) {

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    if (fseeko(f, 0, SEEK_END) == 0)
      retval = ftello(f);
# else
    if (fseek(f, 0, SEEK_END) == 0)
      retval = ftell(f);
# endif

    fclose(f);
  }

#endif /* defined(HAVE_SYS_STAT_H) */

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove a file if it exists and is a regular file or an empty
 *        directory.
 *
 * \param[in]  path  file path.
 *
 * \return 0 in case of success or if file does not exist,  not 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_remove(const char  *path)
{
  int retval = 0;

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) \
                              && defined(HAVE_UNISTD_H)

  struct stat s;

  if (stat(path, &s) == 0) {
    if (S_ISREG(s.st_mode) != 0) {
      retval = unlink(path);
      if (retval != 0) {
        /* Some error types are accepted */
        if (errno == ENOENT)
          retval = 0;
      }
    }
    else if (S_ISDIR(s.st_mode) != 0) {
      retval = rmdir(path);
      if (retval != 0) {
        /* Some error types are accepted */
        if (   errno == ENOTDIR || errno == EEXIST
            || errno == ENOTEMPTY || errno == EBUSY)
          retval = 0;
      }
    }
  }

#else

  /* If Posix-type API is not available, revert to basic method */

  FILE *f;

  if ((f = fopen(path, "w")) != NULL) {
    fclose(f);
    retval = remove(f);
    /* Some error types are accepted */
    if (errno == ENOENT)
      retval = 0;
  }

#endif

  if (retval != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error removing file \"%s\":\n\n"
                "  %s"), path, strerror(errno));

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if file name/path ends with a specific string
 *
 * The function returns an int: 1 if the file name ends with the
 * given string, 0 otherwise.
 *
 * \param[in]  path  name of file
 * \param[in]  end   end string to test
 *
 * \return  1 if the path ends with the given string, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_file_endswith(const char  *path,
                 const char  *end)
{
  int retval = 0;

  /* If either pointers is NULL, return 0 */
  if (path == NULL || end == NULL)
    retval = 0;

  else {

    const int lpath = strlen(path);
    const int lext  = strlen(end);

    /* If either string is empty, or if the path is shorter than the end
     * string, return 0 */
    if (lpath == 0 || lext == 0)
      retval = 0;

    else if (lext > lpath)
      retval = 0;

    else
      retval = (strcmp(path + (lpath-lext), end) == 0);
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
