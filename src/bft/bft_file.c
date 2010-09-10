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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
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

#include "bft_config_defs.h"

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H)
#include <sys/stat.h>
#include <sys/types.h>
#endif /* defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) */

#if defined(HAVE_DIRENT_H)
#include <dirent.h>
#endif

#if defined(WIN32) || defined(_WIN32)
#include <io.h>
#endif

/*
 * Optional library and BFT headers
 */

#include "bft_error.h"
#include "bft_file.h"
#include "bft_mem.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*
 * BFT file descriptor
 */

struct _bft_file_t {

  FILE              *ptr;         /* File descriptor */

  char              *name;        /* File name */
  bft_file_mode_t    mode;        /* File mode */
  bft_file_type_t    type;        /* Type (text, binary, Fortan binary) */
  int                swp_endian;  /* Swap big-endian and little-endian ? */

};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/* Associated typedef documentation (for bft_file.h) */

/*!
 * \typedef bft_file_t
 * \brief Pointer to opaque file descriptor
 */

/*-----------------------------------------------------------------------------
 * Local macro definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default file I/O error handler.
 *
 * The general error handler used by bft_error() is called (which results in the
 * termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bft_file_error_handler_default(const char  *const file_name,
                                const int          line_num,
                                const int          sys_error_code,
                                const char  *const format,
                                va_list            arg_ptr);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static bft_error_handler_t  *_bft_file_error_handler
                               = (_bft_file_error_handler_default);

/* Message strings that may be used more than once */

static const char * _bft_file_str_b_read_error
                      = N_("Error reading binary file \"%s\":"
                           "\n\n  %s");
static const char * _bft_file_str_b_read_closed_error
                      = N_("Error: reading from closed file \"%s\"");
static const char * _bft_file_str_b_write_error
                      = N_("Error writing binary file \"%s\":"
                           "\n\n  %s");
static const char * _bft_file_str_f_read_error
                      = N_("Error reading Fortran binary file \"%s\":"
                           "\n\n  %s");
static const char * _bft_file_str_f_write_error
                      = N_("Error writing Fortran binary file \"%s\":"
                           "\n\n  %s");
static const char * _bft_file_str_f_rec_too_large
                      = N_("A record is too large to be represented "
                           "in this format (i.e. > 2GB).");

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default file I/O error handler.
 *
 * The general error handler used by bft_error() is called (which results in the
 * termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bft_file_error_handler_default(const char  *const file_name,
                                const int          line_num,
                                const int          sys_error_code,
                                const char  *const format,
                                va_list            arg_ptr)
{
  bft_error_handler_t * general_err_handler;

  general_err_handler = bft_error_handler_get();
  general_err_handler(file_name, line_num, sys_error_code, format, arg_ptr);
}

/*
 * Calls the file I/O error handler (set by bft_file_error_handler_set()
 * or default).
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameters:
 *   file_name:      <-- name of source file from which failed bft_mem_...()
 *                       function was called.
 *   line_num:       <-- line of source file from which failed bft_mem_...()
 *                       function was called.
 *   sys_error_code: <-- error code if error in system or libc call,
 *                       0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   ... :           <-- variable arguments based on format string.
 */

static void
_bft_file_error(const char  *file_name,
                const int    line_num,
                const int    sys_error_code,
                const char  *format,
                ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  _bft_file_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*
 * Return an error message associated with an open file.
 *
 * parameter:
 *   f: <-- file descriptor.
 */

static const char *
_bft_file_error_string(const bft_file_t *const f)
{
  int err_num;

  assert(f != NULL) ;

  if (f->ptr == NULL) return "\0";

  err_num = ferror(f->ptr);

  if (err_num != 0)
    return strerror(err_num);

  else if (feof(f->ptr) != 0)
    return _("Premature end of file.");

  else
    return "\0";
}

/*
 * Formatted input from a text file if possible (as fgets()).
 *
 * This function is the base for bft_file_gets() and bft_file_gets_try();
 * depending on the allow_eof parameter, failure to read a line due to
 * an end-of-file condition is considered an error or not.
 *
 * parameters:
 *   s:         --> buffer to which string is to be read.
 *   size:      <-- maximum number of characters to be read plus one.
 *   f:         <-- bft_file_t descriptor.
 *   line:      <-> file line number if available, or NULL.
 *   allow_eof: <-- 1 if EOF is allowed, 0 if considered an error.
 *
 * returns:
 *   s on success, NULL on error or when end of file occurs and
 *   no characters have been read.
 */

static char *
_bft_file_gets(char              *const s,
               const int                size,
               const bft_file_t  *const f,
               int               *const line,
               const int                allow_eof)
{
  char *retval = NULL;

  assert(f != NULL);

  if (f->ptr != NULL)
    retval = fgets(s, size, f->ptr);

  else
    _bft_file_error(__FILE__, __LINE__, 0,
                    _(_bft_file_str_b_read_closed_error), f->name);

  if (retval == NULL && (allow_eof == 0 || bft_file_eof(f) == 0)) {
    if (line != NULL)
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error reading line %d of text file \"%s\":\n\n  %s"),
                      *line, f->name, _bft_file_error_string(f));
    else
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error reading text file \"%s\":\n\n  %s"),
                      f->name, _bft_file_error_string(f));
  }

  if (retval != NULL && line != NULL)
    *line += 1;

  return retval;
}

/*
 * Reads a "usual" Fortran binary record header or footer.
 *
 * Depending on the allow_eof parameter, failure to read a line due to
 * an end-of-file condition is considered an error or not.
 *
 * parameters:
 *   f:        <-- bft_file_t descriptor.
 *   rec_size: <-- size of corresponding record.
 *
 * returns:
 *   0 on success, 1 if the record format or size does not correspond
 *   to what we expect.
 */

static int
_bft_file_read_fortran_size(const bft_file_t  *const f,
                            const size_t             rec_size,
                            const int                allow_eof)
{
  int32_t  n_bytes_read;
  size_t   retval = 0;

  assert(sizeof(int32_t) == 4);

  if (f->ptr != NULL)
    retval = fread((void *)(&n_bytes_read), sizeof(int32_t), 1, f->ptr);

  else
    _bft_file_error(__FILE__, __LINE__, 0,
                    _(_bft_file_str_b_read_closed_error), f->name);

  if (retval < 1) {
    if (allow_eof == 0 || bft_file_eof(f) == 0)
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_read_error),
                      f->name, _bft_file_error_string(f));
    return 1;
  }

  if (f->swp_endian == 1)
    bft_file_swap_endian(&n_bytes_read, &n_bytes_read, sizeof(int32_t), 1);

  if ((size_t)n_bytes_read != rec_size) {
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error reading Fortran binary file \"%s\":\n\n"
                      "  expected record size: %lu\n"
                      "  read record size:     %lu\n"), f->name,
                    (unsigned long)rec_size, (unsigned long)n_bytes_read);
    return 1;
  }

  return 0;
}

/*
 * Read a binary C or Fortran type record.
 *
 * This function is the base for bft_file_read() and bft_file_read_try();
 * depending on the allow_eof parameter, failure to read a line due to
 * an end-of-file condition is considered an error or not.
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
 *   rec  <-> pointer to location receiving data.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of items to write.
 *   f    <-- bft_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran record,
 *   if the whole record could not be read, returns 0.
 */

static size_t
_bft_file_read(void              *const rec,
               const size_t             size,
               const size_t             ni,
               const bft_file_t  *const f,
               const int                allow_eof)
{
  int32_t  n_bytes;
  size_t   retval;
  size_t   rec_size;

  assert(sizeof(int32_t) == 4);

  /* Check file state */

  assert(f != NULL);
  assert(rec != NULL || ni == 0);
  assert(   f->type == BFT_FILE_TYPE_BINARY
         || f->type == BFT_FILE_TYPE_FORTRAN_BINARY);
  assert(f->mode == BFT_FILE_MODE_READ);

  /* Number of bytes of record to read */

  rec_size = size * ni;

  /* In Fortran binary case, read record header */

  if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY) {

    /* Check that 4 bytes is enough for record */

    n_bytes = (int32_t)rec_size;

    if ((size_t)n_bytes != rec_size) {
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_read_error),
                      f->name, _(_bft_file_str_f_rec_too_large));
      return 0;
    }

    if (_bft_file_read_fortran_size(f, rec_size, allow_eof) != 0)
      return 0;

  }

  /* Read the record proper (C or Fortran) */

  retval = 0;

  if (f->ptr != NULL)
    retval = fread((void *)rec, size, ni, f->ptr);

  else
    _bft_file_error(__FILE__, __LINE__, 0,
                    _(_bft_file_str_b_read_closed_error), f->name);

  if (retval != (size_t)ni) {
    if (allow_eof == 0 || bft_file_eof(f) == 0) {
      if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY) {
        _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_read_error),
                        f->name, _bft_file_error_string(f));
        retval = 0;
      }
      else
        _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_b_read_error),
                        f->name, _bft_file_error_string(f));
    }
    return retval;
  }

  if (f->swp_endian == 1 && size > 1)
    bft_file_swap_endian(rec, rec, size, ni);

  /* In Fortran binary case, read record footer */

  if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY) {

    if (_bft_file_read_fortran_size(f, rec_size, allow_eof) != 0)
      return 0;

  }

  return retval;
}

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
_bft_file_compare_names(const void  *a,
                        const void  *b)
{
  return strcmp(*((const char **)a), *((const char **)b));
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Create a `bft_file_t' file descriptor and open the associated file.
 *
 * By default, data is written or read in native format (as regards
 * big-endian or little-endian). This behavior may be modified by
 * bft_file_set_big_endian() or bft_file_set_swap_endian().
 *
 * \param [in] name file name.
 * \param [in] mode file acces mode: read, write, or append.
 * \param [in] type file type: text, binary, or Fortran binary.
 *
 * \return pointer to bft_file_t file descriptor (NULL in case of failure).
 */

bft_file_t *
bft_file_open(const char             *const name,
              const bft_file_mode_t         mode,
              const bft_file_type_t         type)
{

  bft_file_t * f;

  BFT_MALLOC(f, 1, bft_file_t);

  f->ptr = NULL;

  BFT_MALLOC(f->name, strlen(name) + 1, char);
  strcpy(f->name, name);

  f->type = type;
  f->mode = mode;

  /* Use native endianness by default */

  f->swp_endian = 0;

  /* Open file. In case of failure, destroy the allocated structure;
     this is only useful with a non-default error handler,
     as the program is terminated by default */

  if (bft_file_open_stream(f, mode) != 0)
    f = bft_file_free(f);

  return f;

}

/*!
 * \brief Destroy a `bft_file_t' descriptor and close the associated file.
 *
 * The descriptor may only be destroyed if the file was successfully
 * closed. To force destruction of a bft_file_t descriptor even
 * if the associated file was not closed, use (bft_file_free_force()).
 *
 * The associated file is only closed if this was not already the case.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return pointer to bft_file_t file descriptor (NULL in case of,
 *         success, f in case of failure).
 */

bft_file_t *
bft_file_free(bft_file_t  *f)
{
  if (f != NULL) {

    if (bft_file_close_stream(f) == 0)
      return bft_file_free_descriptor(f);

  }

  return f;
}

/*
 * \brief Destroy a `bft_file_t' descriptor without closing its associated file.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * returns:
 *   NULL pointer.
 */

bft_file_t *
bft_file_free_descriptor(bft_file_t  *f)
{
  if (f != NULL) {
    BFT_FREE(f->name);
    BFT_FREE(f);
  }
  return NULL;
}

/*!
 * \brief Open `bft_file_t' descriptor's associated file.
 *
 * If the file is already open, this function does nothing.
 *
 * \param [in] f    bft_file_t descriptor.
 * \param [in] mode file acces mode: read, write, or append.
 *
 * \return 0 in case of success, system error code in case of failure.
 */

int
bft_file_open_stream(bft_file_t       *const f,
                     bft_file_mode_t   const mode)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL)
    return 0;

  /* The file wrapper exists and the corresponding file is closed */

  f->mode = mode;

  switch (f->type) {

  case BFT_FILE_TYPE_TEXT:

    if (f->mode == BFT_FILE_MODE_APPEND)
      f->ptr = fopen(f->name, "a");

    else if (f->mode == BFT_FILE_MODE_WRITE)
      f->ptr = fopen(f->name, "w");

    else if (f->mode == BFT_FILE_MODE_READ)
      f->ptr = fopen(f->name, "r");

    break;

  default:

    if (f->mode == BFT_FILE_MODE_APPEND)
      f->ptr = fopen(f->name, "ab");

    else if (f->mode == BFT_FILE_MODE_WRITE)
      f->ptr = fopen(f->name, "wb");

    else if (f->mode == BFT_FILE_MODE_READ)
      f->ptr = fopen(f->name, "rb");

    break;

  }

  if (f->ptr == NULL) {
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error opening file \"%s\":\n\n"
                      "  %s"), f->name, strerror(errno));
    retval = errno;
  }

  return retval;
}

/*!
 * \brief Close a bft_file_t file descriptor's associated file.
 *
 * If the file is already closed, this function does nothing.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return 0 in case of success, system error code in case of failure.
 */

int
bft_file_close_stream(bft_file_t  *const f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {
    retval = fclose(f->ptr);
    if (retval != 0) {
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error closing file \"%s\":\n\n"
                        "  %s"), f->name, strerror(errno));
      return errno;
    }
    f->ptr = NULL ;
  }

  return retval;
}

/*!
 * \brief Test the end-of-file indicator for a given file.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return 0 if the end-of-file has not been reached, or non-zero
 *         (1 or feof() return value) otherwise.
 */

int
bft_file_eof(const bft_file_t  *const f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL)
    retval = feof(f->ptr);

  return retval;
}

/*!
 * \brief Force write of all user-space buffered data for a given file.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return 0 upon successful completion, system error code otherwise.
 */

int
bft_file_flush(bft_file_t  *const f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {

    retval = fflush(f->ptr);

    if (retval != 0) {
      retval = errno;
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error flushing file \"%s\":\n\n"
                        "  %s"), f->name, strerror(retval));
    }

  }

  return retval;
}

/*!
 * \brief Obtain the current value of a file's position indicator.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return current value of the file's position indicator, or -1 in case
 *         of failure.
 */

bft_file_off_t
bft_file_tell(bft_file_t  *const f)
{
  bft_file_off_t offset = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support might be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    offset = (bft_file_off_t)ftello(f->ptr);
# else
    /*
      Without ftello, ftell will fail above 2 Gigabytes, in which case
      offset == -1 and errno == EOVERFLOW, but should work on smaller
      files. We prefer not to be too strict about fseeko availability, as
      the only 32-bit case without ftello we have encountered is Cygwin
      (for which ftello requires additional non-default libraries), which
      is expected to be used mainly for small cases.
    */
    offset = (bft_file_off_t)ftell(f->ptr);
#endif

#else /* SIZEOF_LONG >= 8) */

    /* For 64-bit systems, standard ftell should be enough */

    offset = (bft_file_off_t)ftell(f->ptr);

#endif

  }

  if (offset < 0)
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error obtaining position in file \"%s\":\n\n  %s"),
                    f->name, _bft_file_error_string(f));

  return offset;
}

/*!
 * \brief Sets the file position indicator to the beginning of the file.
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * \param [in] f bft_file_t descriptor.
 */

void
bft_file_rewind(bft_file_t  *const f)
{
  assert(f != NULL);

  if (f->ptr != NULL) {

    int retval = 0;

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support might be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    retval = fseeko(f->ptr, 0L, SEEK_SET);
# else
    retval = fseek(f->ptr, 0L, SEEK_SET);
# endif

#else /* SIZEOF_LONG >= 8) */

    /* For 64-bit systems, standard fseek should be enough */

    retval = fseek(f->ptr, 0L, SEEK_SET);

#endif

    if (retval != 0)
      _bft_file_error(__FILE__, __LINE__, errno,
                      _("Error rewinding file \"%s\":\n\n  %s"),
                      f->name, _bft_file_error_string(f));
  }
}

/*!
 * \brief Sets a file's position indicator.
 *
 * This function may call the libc's fseek() function. The C 99 standard draft
 * specifies that for a text file, the offset argument to fseek() should be
 * zero or a value returned by an earlier successful call to ftell()
 * (here bft_file_ftell()) on a stream (here a bft_file_t structure).
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * \param [in] f      bft_file_t descriptor.
 * \param [in] offset add to position specified to whence to obtain new
 *                    position, measured in characters from the beginning of
 *                    the file.
 * \param [in] whence beginning if BFT_FILE_SEEK_SET, current if
 *                    BFT_FILE_SEEK_CUR, or end-of-file if BFT_FILE_SEEK_END.
 *
 * \return 0 upon success, nonzero otherwise.
 */

int
bft_file_seek(bft_file_t             *const f,
              const bft_file_off_t          offset,
              const bft_file_seek_t         whence)
{
  int _whence = BFT_FILE_SEEK_SET;
  int retval = 0;

  assert(f != NULL);

  switch (whence) {
  case BFT_FILE_SEEK_SET:
    _whence = SEEK_SET;
    break;
  case BFT_FILE_SEEK_CUR:
    _whence = SEEK_CUR;
    break;
  case BFT_FILE_SEEK_END:
    _whence = SEEK_END;
    break;
  default:
    _bft_file_error
      (__FILE__, __LINE__, 0,
       _("Invalid offset argument \"%d\" setting position in file\n"
         "\"%s\""),
       (int)whence, f->name);
  }

  if (f->ptr != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support might be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)

    retval = fseeko(f->ptr, (off_t)offset, _whence);

    if (retval != 0)
      _bft_file_error(__FILE__, __LINE__, errno,
                      _("Error setting position in file \"%s\":\n\n  %s"),
                      f->name, _bft_file_error_string(f));

# else

    /* Test if offset larger than allowed */

    long _offset = offset;

    if (_offset == offset) {
      retval = fseek(f->ptr, (long)offset, _whence);
      if (retval != 0)
        _bft_file_error(__FILE__, __LINE__, errno,
                        _("Error setting position in file \"%s\":\n\n  %s"),
                        f->name, _bft_file_error_string(f));
    }
    else {
      retval = -1;
      _bft_file_error
        (__FILE__, __LINE__, errno,
         _("Error setting position in file \"%s\":\n\n  %s"),
         f->name,
         _("sizeof(off_t) > sizeof(long) but fseeko() not available"));
    }

# endif /* defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64) */

#else /* SIZEOF_LONG >= 8) */

    /* For 64-bit systems, standard fseek should be enough */

    retval = fseek(f->ptr, (long)offset, _whence);

    if (retval != 0)
      _bft_file_error(__FILE__, __LINE__, errno,
                      _("Error setting position in file \"%s\":\n\n  %s"),
                      f->name, _bft_file_error_string(f));

#endif /* SIZEOF_LONG */
  }

  return retval;
}

/*!
 * \brief Return a file's name.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return pointer to file's name.
 */

const char *
bft_file_get_name(const bft_file_t  *const f)
{
  assert(f != NULL);

  return f->name;
}

/*!
 * \brief Return a file's type.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return file's type.
 */

bft_file_type_t
bft_file_get_type(const bft_file_t  *const f)
{
  assert(f != NULL);

  return f->type;
}

/*!
 * \brief Change a file's type.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * \param [in, out ] f    bft_file_t descriptor.
 * \param [in]       type text, binary, or Fortran binary type descriptor.
 */

void
bft_file_set_type(bft_file_t             *const f,
                  const bft_file_type_t         type)
{
  assert(f != NULL);

  f->type = type;
}

/*!
 * \brief Ensure that data is read or written in big-endian
 *        (network standard) format.
 *
 * By default, data is written or read in native format (as regards
 * big-endian or little-endian)..
 *
 * \param [in, out] f bft_file_t descriptor.
 */

void
bft_file_set_big_endian(bft_file_t  *const f)
{
  unsigned     int_endian;

  /* Check if system is "big-endian" or "little-endian" */

  int_endian = 0;
  *((char *)(&int_endian)) = '\1';

  if (int_endian == 1)
    f->swp_endian = 1;

#if defined(DEBUG) && !defined(NDEBUG)

  else {
    int_endian = 0;
    *((char *) (&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert(int_endian == 1);
  }

#endif

}

/*!
 * \brief Return a file's byte-swapping behavior.
 *
 * \param [in] f bft_file_t descriptor.
 *
 * \return 0 if file's endianness is the same as the system's, 1 otherwise.
 */

int
bft_file_get_swap_endian(const bft_file_t  *const f)
{
  assert(f != NULL);

  return f->swp_endian;
}

/*!
 * \brief Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * \param [in] f    bft_file_t descriptor.
 * \param [in] swap 1 if bytes must be swapped, 0 otherwise.
 */

void
bft_file_set_swap_endian(bft_file_t  *const f,
                         const int          swap)
{
  assert(f != NULL);

  f->swp_endian = swap;
}

/*!
 * \brief Test a file's error or EOF condition.
 *
 * \param [in]      f    bft_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return 0 if no error, system error code, or -1 if EOF.
 */

int
bft_file_read_check_error(const bft_file_t  *const f,
                          const int                line)
{
  int retval = 0;

  /* Check for possible error */

  if (f->ptr != NULL)
    retval = ferror(f->ptr);

  /* If we have a read error */

  if (retval != 0) {
    if (line > 0)
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error reading line %d of file \"%s\":\n\n"
                        "  %s"), line, f->name, _bft_file_error_string(f));
    else
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error reading file \"%s\":\n\n"
                        "  %s"), f->name, _bft_file_error_string(f));
    return retval;
  }

  /* Check for EOF condition */

  if (f->ptr != NULL)
    retval = feof(f->ptr);

  if (retval != 0) {
    if (line > 0)
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Premature end of file \"%s\" at line %d\n\n"),
                      f->name, line);
    else
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Premature end of file \"%s\"\n\n"), f->name);
    retval = -1;
  }

  return retval;
}

/*!
 * \brief Formatted output to a text file (as fprintf()).
 *
 * \param [in] f      bft_file_t descriptor.
 * \param [in] format format string, as printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int bft_file_printf(const bft_file_t  *const f,
                    const char        *const format,
                    ...)
{
  va_list       arg_ptr;
  int           retval = 0;

  assert(f != NULL) ;
  assert(f->mode == BFT_FILE_MODE_APPEND || f->mode == BFT_FILE_MODE_WRITE);
  assert(f->type == BFT_FILE_TYPE_TEXT);

  if (f->ptr != NULL) {

    va_start(arg_ptr, format);

    retval = vfprintf(f->ptr, format, arg_ptr);

    va_end(arg_ptr);

    if (retval <= 0) {
      retval = errno;
      _bft_file_error(__FILE__, __LINE__, 0,
                      _("Error writing to text file \"%s\":\n\n  %s"),
                      f->name, _bft_file_error_string(f));
    }

  }

  else
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error writing to closed file \"%s\""), f->name);

  return retval;
}

/*!
 * \brief Formatted input from a text file (as fgets()).
 *
 * \param [out]     s    buffer to which string is to be read.
 * \param [in]      size maximum number of characters to be read plus one.
 * \param [in]      f    bft_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */

char *
bft_file_gets(char              *const s,
              const int                size,
              const bft_file_t  *const f,
              int               *const line)
{
  return _bft_file_gets(s, size, f, line, 0);
}

/*!
 * \brief Formatted input from a text file if possible (as fgets()).
 *
 * This function is similar to bft_file_gets(), but failure to read
 * a line due to an end-of-file condition is not considered an error with
 * this variant, which may be used to read text files or sections thereof
 * of unknown length.
 *
 * \param [out]     s    buffer to which string is to be read.
 * \param [in]      size maximum number of characters to be read plus one.
 * \param [in]      f    bft_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */

char *
bft_file_gets_try(char              *const s,
                  const int                size,
                  const bft_file_t  *const f,
                  int               *const line)
{
  return _bft_file_gets(s, size, f, line, 1);
}

/*!
 * \brief Read a binary C or Fortran type record.
 *
 * A Fortran record compatible with most compilers is structured
 * as follows:
 *   - a 4-byte integer indicating the number of bytes in the record.
 *   - the raw data
 *   - a 4-byte integer indicating the number of bytes in the record.
 *
 * A C record contains only the raw data.
 *
 * \param [out] rec  pointer to location receiving data.
 * \param [in]  size size of each item of data in bytes.
 * \param [in]  ni   number of items to read.
 * \param [in]  f    bft_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully read; for a Fortran
 *         record, if the whole record could not be read, returns 0.
 */

size_t
bft_file_read(void              *const rec,
              const size_t             size,
              const size_t             ni,
              const bft_file_t  *const f)
{
  return _bft_file_read(rec, size, ni, f, 0);
}

/*!
 * \brief Read a binary C or Fortran type record.
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
 * \param [out] rec  pointer to location receiving data.
 * \param [in]  size size of each item of data in bytes.
 * \param [in]  ni   number of items to read.
 * \param [in]  f    bft_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully read; for a Fortran
 *         record, if the whole record could not be read, returns 0.
 */

size_t
bft_file_read_try(void              *const rec,
                  const size_t             size,
                  const size_t             ni,
                  const bft_file_t  *const f)
{
  return _bft_file_read(rec, size, ni, f, 1);
}

/*!
 * \brief Write a binary C or Fortran type record.
 *
 * A Fortran record compatible with most compilers is structured
 * as follows:
 *   - a 4-byte integer indicating the number of bytes in the record.
 *   - the raw data
 *   - a 4-byte integer indicating the number of bytes in the record.
 *
 * A C record contains only the raw data.
 *
 * \param [in] rec  pointer to location containing data.
 * \param [in] size size of each item of data in bytes.
 * \param [in] ni   number of items to write.
 * \param [in] f    bft_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully written.
 */

size_t
bft_file_write(const void        *const rec,
               const size_t             size,
               const size_t             ni,
               const bft_file_t  *const f)
{
  int32_t  n_bytes;
  size_t   retval;
  size_t   rec_size;

  assert(sizeof(int32_t) == 4);

  /* Check file state */

  assert(f != NULL);
  assert(rec != NULL || ni == 0);
  assert(   f->type == BFT_FILE_TYPE_BINARY
         || f->type == BFT_FILE_TYPE_FORTRAN_BINARY);
  assert(f->mode == BFT_FILE_MODE_APPEND || f->mode == BFT_FILE_MODE_WRITE);

  if (f->ptr == NULL)
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error writing to closed file \"%s\""), f->name);

  /* Number of bytes of record to write */

  rec_size = size * ni;

  /* In Fortran binary case, write record header */

  if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY) {

    /* Check that 4 bytes is enough for record */

    n_bytes = (int32_t)rec_size;

    if ((size_t)n_bytes != rec_size) {
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_write_error),
                      f->name, _(_bft_file_str_f_rec_too_large));
      return 0;
    }

    if (f->swp_endian == 1)
      bft_file_swap_endian((void *)&n_bytes, (void *)&n_bytes,
                           sizeof(int32_t), 1) ;

    if (fwrite((void *)(&n_bytes), sizeof(int32_t), 1, f->ptr) != 1) {
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_write_error),
                      f->name, _bft_file_error_string(f));
      return 0;
    }

  }

  /* Write the record proper (C or Fortran) */

  if (f->swp_endian == 1 && size > 1) {

    void    *buf;

    BFT_MALLOC(buf, rec_size, unsigned char);

    bft_file_swap_endian(buf, rec, size, ni);

    retval = fwrite((void *)buf, (size_t)size, (size_t)ni, f->ptr);

    BFT_FREE(buf);

  }
  else

    retval = fwrite((const void *)rec, (size_t)size, (size_t)ni, f->ptr);

  if (retval != (size_t)ni) {
    if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY)
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_write_error),
                      f->name, _bft_file_error_string(f));
    else
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_b_write_error),
                      f->name, _bft_file_error_string(f));
    return retval;
  }

  /* In Fortran binary case, write record footer */

  if (f->type == BFT_FILE_TYPE_FORTRAN_BINARY) {

    if (fwrite((void *)(&n_bytes), sizeof(int32_t), 1, f->ptr) != 1) {
      _bft_file_error(__FILE__, __LINE__, 0, _(_bft_file_str_f_write_error),
                      f->name, _bft_file_error_string(f));
      return 0;
    }

  }

  return retval;
}

/*!
 * \brief Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * \param [out] dest pointer to converted data location.
 * \param [in]  src  pointer to source data location.
 * \param [in]  size size of each item of data in bytes.
 * \param [in]  ni   number of data items.
 */

void
bft_file_swap_endian(void *const          dest,
                     const void    *const src,
                     const size_t         size,
                     const size_t         ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0 ; i < ni ; i++) {

    shift = i * size;

    for (ib = 0 ; ib < (size / 2) ; ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

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
 * \param [in] pathname name of new directory.
 *
 * \returns 0 on success, -1 if an error occured (in which case errno
 *          contains the appropriate error code). If the underlying
 *          system has no mkdir() function or it was not detected
 *          upon BFT configuration, 1 is returned.
 */

int
bft_file_mkdir_default(const char  *const pathname)
{
  static const char  *str_fail = N_("Failure to create "
                                    "directory \"%s\":\n\n%s");

#if defined(HAVE_MKDIR)

#if defined(WIN32) || defined(_WIN32)

  mkdir(pathname);
  return 0;

#else

  if (mkdir(pathname, S_IRWXU|S_IRWXG|S_IRWXO) != 0) {

    if (errno == EEXIST) {

#if defined(HAVE_SYS_STAT_H)

      struct stat buf;

      if (stat(pathname, &buf) != 0)
        _bft_file_error(__FILE__, __LINE__, 0, _(str_fail),
                        pathname,
                        _("  A similarly named file or directory exists "
                          "and its status is\n  not available."));
      else if (S_ISDIR(buf.st_mode) != 1)
        _bft_file_error(__FILE__, __LINE__, 0, _(str_fail),
                        pathname,
                        _("  A similarly named file exists and is "
                          "not a directory."));
      else
        return 0;

#endif

      errno = EEXIST; /* In case modified by stat() */

    }
    else {
      _bft_file_error(__FILE__, __LINE__, errno, _(str_fail),
                      pathname,
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

/*!
 * \brief Check if a file exists and is a regular file.
 *
 * \param [in] name file name.
 *
 * \returns 1 if file exists and is a regular file, 0 otherwise.
 */

int
bft_file_isreg(const char  *const name)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(name, &s) != 0) {
    if (errno != ENOENT)
      _bft_file_error(__FILE__, __LINE__, errno,
                      _("Error querying information for file:\n%s."),
                      name) ;
  }
  else {
    if (S_ISREG(s.st_mode) != 0)
      retval = 1;
  }

#else /* defined(HAVE_SYS_STAT_H) */

  /* If Posix-type API is not available, revert to basic method */

  FILE *f;

  if ((f = fopen(fic_name, "r")) != NULL) {
    retval = 1;
    fclose(f) ;
  }

#endif /* defined(HAVE_SYS_STAT_H) */

  return retval;
}


/*!
 * \brief Check if a directory exists.
 *
 * \param [in] name directory name.
 *
 * \returns 1 if directory exists, 0 otherwise.
 */

int
bft_file_isdir(const char  *const name)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(name, &s) != 0) {
    if (errno != ENOENT)
      _bft_file_error(__FILE__, __LINE__, errno,
                      _("Error querying information for directory:\n%s."),
                      name) ;
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

/*!
 * \brief List files inside a directory.
 *
 * The array returned must be freed by the caller using BFT_FREE,
 * as well as the individual entries in the array.
 *
 * \param [in] pathname name of directory.
 *
 * \returns an array of file names in a directory. The last entry is
 *          set to NULL. If no means to list the directory or an error
 *          occured, the return value is simply NULL.
 */

char **
bft_file_listdir(const char  *pathname)
{
  char **dirnames = NULL;

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_DIRENT_H)

  struct dirent *ent;
  int n_ent = 0;
  DIR *d = opendir(pathname);

  if (d == NULL) {
    _bft_file_error(__FILE__, __LINE__, 0,
                    _("Error opening directory \"%s\":\n\n"
                      "  %s"), pathname, strerror(errno));
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

  qsort(dirnames, n_ent, sizeof(char *), &_bft_file_compare_names);

#endif /* defined(HAVE_SYS_TYPES_H) && defined(HAVE_DIRENT_H) */

  return dirnames;
}

/*!
 * \brief Returns the error handler associated with the bft_file_...()
 *        functions.
 *
 * \return pointer to the error handler function.
 */

bft_error_handler_t *
bft_file_error_handler_get(void)
{
  return _bft_file_error_handler;
}

/*!
 * \brief Associates an error handler with the bft_file_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * \param [in] handler pointer to the error handler function.
 */

void
bft_file_error_handler_set(bft_error_handler_t  *const handler)
{
  _bft_file_error_handler = handler;
}

/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through bft_mem_...() (in kB).
 */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
