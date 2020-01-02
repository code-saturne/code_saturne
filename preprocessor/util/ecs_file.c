/*============================================================================
 * Base file wrapper type and associated functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_config.h"

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

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#if defined(HAVE_ZLIB)
#include <zlib.h>
#endif /* defined(_HAVE_ZLIB) */

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H)
#include <sys/stat.h>
#include <sys/types.h>
#endif /* defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) */

#if defined(WIN32) || defined(_WIN32)
#include <io.h>
#endif

/*
 * Optional library and ECS headers
 */

#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"

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
 * ECS file descriptor
 */

struct _ecs_file_t {

  FILE              *ptr;         /* File descriptor */

#if defined(HAVE_ZLIB)
  gzFile             gzptr;       /* Zlib file descriptor */
#endif

  char              *name;        /* File name */
  ecs_file_mode_t    mode;        /* File mode */
  ecs_file_type_t    type;        /* Type (text, binary, Fortan binary) */
  int                swp_endian;  /* Swap big-endian and little-endian ? */

};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/* Associated typedef documentation (for ecs_file.h) */

/*!
 * \typedef ecs_file_t
 * \brief Pointer to opaque file descriptor
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(HAVE_ZLIB)

/* Zlib API may be broken when using large file support, as z_off_t
   is based on current off_t, and not on a value fixed at compilation time.
   We redefine prototypes for gzseek() and gztell() ;
   This is ugly, but should work with an unmodified Zlib (as of Zlib 1.2.1) */

#if defined (SIZEOF_Z_OFF_T)
#  if (SIZEOF_Z_OFF_T == SIZEOF_LONG)
typedef long ecs_z_off_t;
#  elif defined (HAVE_LONG_LONG)
#    if (SIZEOF_Z_OFF_T == SIZEOF_LONG_LONG)
typedef long long ecs_z_off_t;
#    else
#      error "z_off_t returned by zlibCompileFlags() neither long nor long long"
#    endif
#  endif
#else
typedef z_off_t ecs_z_off_t;
#endif

typedef ecs_z_off_t (ecs_gzseek_t) (gzFile file,
                                    ecs_z_off_t offset,
                                    int whence);

typedef ecs_z_off_t (ecs_gztell_t) (gzFile file);

#endif /* HAVE_ZLIB */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local macro definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

/* Message strings that may be used more than once */

static const char * _ecs_file_str_b_read_error
                      = N_("Error reading binary file \"%s\":"
                           "\n\n  %s");
static const char * _ecs_file_str_b_read_closed_error
                      = N_("Error: reading from closed file \"%s\"");
static const char * _ecs_file_str_b_write_error
                      = N_("Error writing binary file \"%s\":"
                           "\n\n  %s");
static const char * _ecs_file_str_f_read_error
                      = N_("Error reading Fortran binary file \"%s\":"
                           "\n\n  %s");
static const char * _ecs_file_str_f_write_error
                      = N_("Error writing Fortran binary file \"%s\":"
                           "\n\n  %s");
static const char * _ecs_file_str_f_rec_too_large
                      = N_("A record is too large to be represented "
                           "in this format (i.e. > 2GB).");

#if defined(HAVE_ZLIB)

/* Zlib API may be broken when using large file support, as z_off_t
   is based on current off_t, and not on a value fixed at compilation time.
   Prototypes for gzseek() and gztell() were redefined above, we
   point to the true functions here.
   This is ugly, but should work with an unmodified Zlib (as of Zlib 1.2.1) */

static ecs_gzseek_t  *_ecs_gzseek = (ecs_gzseek_t *)gzseek;
static ecs_gztell_t  *_ecs_gztell = (ecs_gztell_t *)gztell;

#endif /* HAVE_ZLIB */

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Return an error message associated with an open file.
 *
 * parameter:
 *   f: <-- file descriptor.
 */

static const char *
_ecs_file_error_string(const ecs_file_t  *f)
{
  int err_num;

  assert(f != NULL) ;

#if defined(HAVE_ZLIB)

  if (f->gzptr != NULL) {

    const char *err_str;

    err_str = gzerror(f->gzptr, &err_num);

    if (err_num != 0)
      return err_str;

    else if (gzeof(f->gzptr) != 0)
      return _("Premature end of file.");

    else
      return "\0";

  }

#endif /*defined(HAVE_ZLIB) */

  if (f->ptr == NULL) return "\0";

  err_num = ferror(f->ptr);

  if (err_num != 0)
    return strerror(err_num);

  else if (feof(f->ptr) != 0)
    return _("Premature end of file.");

  else
    return "\0";
}

/*!
 * \brief Formatted output to a text file (as fprintf()).
 *
 * \param [in] f      ecs_file_t descriptor.
 * \param [in] format format string, as printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int ecs_file_printf(const ecs_file_t  *const f,
                    const char        *const format,
                    ...)
{
  va_list       arg_ptr;
  int           retval = 0;

  assert(f != NULL) ;
  assert(f->mode == ECS_FILE_MODE_APPEND || f->mode == ECS_FILE_MODE_WRITE);
  assert(f->type == ECS_FILE_TYPE_TEXT);

  if (f->ptr != NULL) {

    va_start(arg_ptr, format);

    retval = vfprintf(f->ptr, format, arg_ptr);

    va_end(arg_ptr);

    if (retval <= 0) {
      retval = errno;
      ecs_error(__FILE__, __LINE__, 0,
                _("Error writing to text file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));
    }

  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    ecs_error(__FILE__, __LINE__, 0,
              _("ECS library formatted output to gzipped file "
                "not implemented\n\n(file: \"%s\")"), f->name);

#endif /* defined(HAVE_ZLIB) */

  else
    ecs_error(__FILE__, __LINE__, 0,
              _("Error writing to closed file \"%s\""), f->name);

  return retval;
}

/*
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
 */

static char *
_ecs_file_gets(char              *s,
               const int          size,
               const ecs_file_t  *f,
               int               *line,
               const int          allow_eof)
{
  char *retval = NULL;

  assert(f != NULL);

  if (f->ptr != NULL)
    retval = fgets(s, size, f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    retval = gzgets(f->gzptr, s, size);

#endif /* defined(HAVE_ZLIB) */

  else
    ecs_error(__FILE__, __LINE__, 0,
              _(_ecs_file_str_b_read_closed_error), f->name);

  if (retval == NULL && (allow_eof == 0 || ecs_file_eof(f) == 0)) {
    if (line != NULL)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading line %d of text file \"%s\":\n\n  %s"),
                *line, f->name, _ecs_file_error_string(f));
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading text file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));
  }

  if (retval != NULL) {

    int i = strlen(s) - 2;
    if (i > 0) {
      if (s[i] == '\r' && s[i+1] == '\n') {
        s[i] = '\n';
        s[i+1] = '\0';
      }
    }

    if (line != NULL)
      *line += 1;
  }

  return retval;
}

/*
 * Reads a "usual" Fortran binary record header or footer.
 *
 * Depending on the allow_eof parameter, failure to read a line due to
 * an end-of-file condition is considered an error or not.
 *
 * parameters:
 *   f:        <-- ecs_file_t descriptor.
 *   rec_size: <-- size of corresponding record.
 *
 * returns:
 *   0 on success, 1 if the record format or size does not correspond
 *   to what we expect.
 */

static int
_ecs_file_read_fortran_size(const ecs_file_t  *f,
                            const size_t       rec_size,
                            const int          allow_eof)
{
  int32_t  n_bytes_read;
  size_t   retval = 0;

  assert(sizeof(int32_t) == 4);

  if (f->ptr != NULL)
    retval = fread((void *)(&n_bytes_read), sizeof(int32_t), 1, f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL) {
    retval = (size_t)gzread(f->gzptr, (void *)(&n_bytes_read), sizeof(int32_t));
    retval = retval / sizeof(int32_t);
  }

#endif /* defined(HAVE_ZLIB) */

  else
    ecs_error(__FILE__, __LINE__, 0,
              _(_ecs_file_str_b_read_closed_error), f->name);

  if (retval < 1) {
    if (allow_eof == 0 || ecs_file_eof(f) == 0)
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_read_error),
                f->name, _ecs_file_error_string(f));
    return 1;
  }

  if (f->swp_endian == 1)
    ecs_file_swap_endian(&n_bytes_read, &n_bytes_read, sizeof(int32_t), 1);

  if ((size_t)n_bytes_read != rec_size) {
    ecs_error(__FILE__, __LINE__, 0,
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
 * This function is the base for ecs_file_read() and ecs_file_read_try();
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
 *   f    <-- ecs_file_t descriptor.
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; for a Fortran record,
 *   if the whole record could not be read, returns 0.
 */

static size_t
_ecs_file_read(void              *rec,
               const size_t       size,
               const size_t       ni,
               const ecs_file_t  *f,
               const int          allow_eof)
{
  int32_t  n_bytes;
  size_t   retval;
  size_t   rec_size;

  assert(sizeof(int32_t) == 4);

  /* Check file state */

  assert(f != NULL);
  assert(rec != NULL || ni == 0);
  assert(   f->type == ECS_FILE_TYPE_BINARY
         || f->type == ECS_FILE_TYPE_FORTRAN_BINARY);
  assert(f->mode == ECS_FILE_MODE_READ);

  /* Number of bytes of record to read */

  rec_size = size * ni;

  /* In Fortran binary case, read record header */

  if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY) {

    /* Check that 4 bytes is enough for record */

    n_bytes = (int32_t)rec_size;

    if ((size_t)n_bytes != rec_size) {
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_read_error),
                f->name, _(_ecs_file_str_f_rec_too_large));
      return 0;
    }

    if (_ecs_file_read_fortran_size(f, rec_size, allow_eof) != 0)
      return 0;

  }

  /* Read the record proper (C or Fortran) */

  retval = 0;

  if (f->ptr != NULL)

    retval = fread((void *)rec, size, ni, f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)

    retval = ((size_t)gzread(f->gzptr, (void *)(rec), rec_size)) / size;

#endif /* defined(HAVE_ZLIB) */

  else
    ecs_error(__FILE__, __LINE__, 0,
              _(_ecs_file_str_b_read_closed_error), f->name);

  if (retval != (size_t)ni) {
    if (allow_eof == 0 || ecs_file_eof(f) == 0) {
      if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY) {
        ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_read_error),
                  f->name, _ecs_file_error_string(f));
        retval = 0;
      }
      else
        ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_b_read_error),
                  f->name, _ecs_file_error_string(f));
    }
    return retval;
  }

  if (f->swp_endian == 1 && size > 1)
    ecs_file_swap_endian(rec, rec, size, ni);

  /* In Fortran binary case, read record footer */

  if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY) {

    if (_ecs_file_read_fortran_size(f, rec_size, allow_eof) != 0)
      return 0;

  }

  return retval;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Create a `ecs_file_t' file descriptor and open the associated file.
 *
 * By default, data is written or read in native format (as regards
 * big-endian or little-endian). This behavior may be modified by
 * ecs_file_set_big_endian() or ecs_file_set_swap_endian().
 *
 * \param [in] name file name.
 * \param [in] mode file acces mode: read, write, or append.
 * \param [in] type file type: text, binary, or Fortran binary.
 *
 * \return pointer to ecs_file_t file descriptor (NULL in case of failure).
 */

ecs_file_t *
ecs_file_open(const char             *name,
              const ecs_file_mode_t   mode,
              const ecs_file_type_t   type)
{

  ecs_file_t * f;

  ECS_MALLOC(f, 1, ecs_file_t);

  f->ptr = NULL;

#if defined(HAVE_ZLIB)
  f->gzptr = NULL;
#endif

  ECS_MALLOC(f->name, strlen(name) + 1, char);
  strcpy(f->name, name);

  f->type = type;
  f->mode = mode;

  /* Use native endianness by default */

  f->swp_endian = 0;

  /* Open file. In case of failure, destroy the allocated structure;
     this is only useful with a non-default error handler,
     as the program is terminated by default */

  if (ecs_file_open_stream(f, mode) != 0)
    f = ecs_file_free(f);

  return f;

}

/*!
 * \brief Destroy a `ecs_file_t' descriptor and close the associated file.
 *
 * The descriptor may only be destroyed if the file was successfully
 * closed. To force destruction of a ecs_file_t descriptor even
 * if the associated file was not closed, use (ecs_file_free_force()).
 *
 * The associated file is only closed if this was not already the case.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return pointer to ecs_file_t file descriptor (NULL in case of,
 *         success, f in case of failure).
 */

ecs_file_t *
ecs_file_free(ecs_file_t  *f)
{
  if (f != NULL) {

    if (ecs_file_close_stream(f) == 0)
      return ecs_file_free_descriptor(f);

  }

  return f;
}

/*
 * \brief Destroy a `ecs_file_t' descriptor without closing its associated file.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * returns:
 *   NULL pointer.
 */

ecs_file_t *
ecs_file_free_descriptor(ecs_file_t  *f)
{
  if (f != NULL) {
    ECS_FREE(f->name);
    ECS_FREE(f);
  }
  return NULL;
}

/*!
 * \brief Open `ecs_file_t' descriptor's associated file.
 *
 * If the file is already open, this function does nothing.
 *
 * \param [in] f    ecs_file_t descriptor.
 * \param [in] mode file acces mode: read, write, or append.
 *
 * \return 0 in case of success, system error code in case of failure
 *         (or Zlib error code in case of Zlib memory allocation problem
 *         for a gzipped file).
 */

int
ecs_file_open_stream(ecs_file_t       *f,
                     ecs_file_mode_t   mode)
{
  int retval = 0;

#if defined(HAVE_ZLIB)
  int gzipped = 0;
#endif

  assert(f != NULL);

  if (f->ptr != NULL
#if defined(HAVE_ZLIB)
      || f->gzptr != NULL
#endif
      )
    return 0;

  /* The file wrapper exists and the corresponding file is closed */

  f->mode = mode;

  switch (f->type) {

  case ECS_FILE_TYPE_TEXT:

    if (f->mode == ECS_FILE_MODE_APPEND)
      f->ptr = fopen(f->name, "a");

    else if (f->mode == ECS_FILE_MODE_WRITE)
      f->ptr = fopen(f->name, "w");

    else if (f->mode == ECS_FILE_MODE_READ) {

#if defined(HAVE_ZLIB)

      if (strlen(f->name) > 3 &&
          (! strncmp((f->name + strlen(f->name) - 3), ".gz", 3))) {

        f->gzptr = gzopen(f->name, "r");
        gzipped = 1;
        break;
      }

#endif /* defined(HAVE_ZLIB) */

      f->ptr = fopen(f->name, "r");

    }

    break;

  default:

    if (f->mode == ECS_FILE_MODE_APPEND)
      f->ptr = fopen(f->name, "ab");

    else if (f->mode == ECS_FILE_MODE_WRITE)
      f->ptr = fopen(f->name, "wb");

    else if (f->mode == ECS_FILE_MODE_READ) {

#if defined(HAVE_ZLIB)

      if (strlen(f->name) > 3 &&
          (! strncmp((f->name + strlen(f->name) - 3), ".gz", 3))) {

        f->gzptr = gzopen(f->name, "rb");
        gzipped = 1;
        break;
      }

#endif /* defined(HAVE_ZLIB) */

      f->ptr = fopen(f->name, "rb");

    }

    break;

  }

  if (f->ptr == NULL
#if defined(HAVE_ZLIB)
      && f->gzptr == NULL
#endif
      ) {
#if defined(HAVE_ZLIB)
    if (gzipped == 1 && errno == 0) {
      ecs_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n\n"
                  "  %s"), f->name, zError(Z_MEM_ERROR));
      return Z_MEM_ERROR;
    }
#endif
    ecs_error(__FILE__, __LINE__, 0,
              _("Error opening file \"%s\":\n\n"
                "  %s"), f->name, strerror(errno));
    retval = errno;
  }

  return retval;
}

/*!
 * \brief Close a ecs_file_t file descriptor's associated file.
 *
 * If the file is already closed, this function does nothing.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return 0 in case of success, system error code in case of failure
 *         (or Zlib error code in case of a Zlib specific error
 *         for a gzipped file).
 */

int
ecs_file_close_stream(ecs_file_t  *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {
    retval = fclose(f->ptr);
    if (retval != 0) {
      ecs_error(__FILE__, __LINE__, 0,
                _("Error closing file \"%s\":\n\n"
                  "  %s"), f->name, strerror(errno));
      return errno;
    }
    f->ptr = NULL ;
  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL) {
    retval = gzclose(f->gzptr);
    if (retval != 0) {
      ecs_error(__FILE__, __LINE__, 0,
                _("Error closing file \"%s\":\n\n"
                  "  %s"), f->name, gzerror(f->gzptr, &retval));
      return retval;
    }
    f->gzptr = NULL ;
  }

#endif /* defined(HAVE_ZLIB) */

  return retval;
}

/*!
 * \brief Test the end-of-file indicator for a given file.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return 0 if the end-of-file has not been reached, or non-zero
 *         (1 or feof() return value) otherwise.
 */

int
ecs_file_eof(const ecs_file_t  *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL)
    retval = feof(f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    retval = gzeof(f->gzptr);

#endif

  return retval;
}

/*!
 * \brief Force write of all user-space buffered data for a given file.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return 0 upon successful completion, system error code otherwise.
 */

int
ecs_file_flush(ecs_file_t  *f)
{
  int retval = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {

    retval = fflush(f->ptr);

    if (retval != 0) {
      retval = errno;
      ecs_error(__FILE__, __LINE__, 0,
                _("Error flushing file \"%s\":\n\n"
                  "  %s"), f->name, strerror(retval));
    }

  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL) {

    retval = gzflush(f->gzptr, Z_FULL_FLUSH);

    if (retval != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error closing file \"%s\":\n\n"
                  "  %s"), f->name, gzerror(f->gzptr, &retval));

  }

#endif

  return retval;
}

/*!
 * \brief Obtain the current value of a file's position indicator.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return current value of the file's position indicator, or -1 in case
 *         of failure.
 */

ecs_file_off_t
ecs_file_tell(ecs_file_t  *f)
{
  ecs_file_off_t offset = 0;

  assert(f != NULL);

  if (f->ptr != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support might be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    offset = (ecs_file_off_t)ftello(f->ptr);
# else
    /*
      Without ftello, ftell will fail above 2 Gigabytes, in which case
      offset == -1 and errno == EOVERFLOW, but should work on smaller
      files. We prefer not to be too strict about fseeko availability, as
      the only 32-bit case without ftello we have encountered is Cygwin
      (for which ftello requires additional non-default libraries), which
      is expected to be used mainly for small cases.
    */
    offset = (ecs_file_off_t)ftell(f->ptr);
#endif

#else /* SIZEOF_LONG >= 8) */

    /* For 64-bit systems, standard ftell should be enough */

    offset = (ecs_file_off_t)ftell(f->ptr);

#endif

  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    offset = (ecs_file_off_t)_ecs_gztell(f->gzptr);

#endif

  if (offset < 0)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error obtaining position in file \"%s\":\n\n  %s"),
              f->name, _ecs_file_error_string(f));

  return offset;
}

/*!
 * \brief Sets the file position indicator to the beginning of the file.
 *
 * A successful call to this function clears the end-of-file indicator for
 * this file.
 *
 * \param [in] f ecs_file_t descriptor.
 */

void
ecs_file_rewind(ecs_file_t  *f)
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
      ecs_error(__FILE__, __LINE__, errno,
                _("Error rewinding file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));
  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    (void)gzrewind(f->gzptr);

#endif
}

/*!
 * \brief Sets a file's position indicator.
 *
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
 * \param [in] f      ecs_file_t descriptor.
 * \param [in] offset add to position specified to whence to obtain new
 *                    position, measured in characters from the beginning of
 *                    the file.
 * \param [in] whence beginning if ECS_FILE_SEEK_SET, current if
 *                    ECS_FILE_SEEK_CUR, or end-of-file if ECS_FILE_SEEK_END.
 *
 * \return 0 upon success, nonzero otherwise.
 */

int
ecs_file_seek(ecs_file_t             *f,
              const ecs_file_off_t    offset,
              const ecs_file_seek_t   whence)
{
  int _whence = ECS_FILE_SEEK_SET;
  int retval = 0;

  assert(f != NULL);

  switch (whence) {
  case ECS_FILE_SEEK_SET:
    _whence = SEEK_SET;
    break;
  case ECS_FILE_SEEK_CUR:
    _whence = SEEK_CUR;
    break;
  case ECS_FILE_SEEK_END:
    _whence = SEEK_END;
    break;
  default:
    ecs_error
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
      ecs_error(__FILE__, __LINE__, errno,
                _("Error setting position in file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));

# else

    /* Test if offset larger than allowed */

    long _offset = offset;

    if (_offset == offset) {
      retval = fseek(f->ptr, (long)offset, _whence);
      if (retval != 0)
        ecs_error(__FILE__, __LINE__, errno,
                  _("Error setting position in file \"%s\":\n\n  %s"),
                  f->name, _ecs_file_error_string(f));
    }
    else {
      retval = -1;
      ecs_error
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
      ecs_error(__FILE__, __LINE__, errno,
                _("Error setting position in file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));

#endif /* SIZEOF_LONG */
  }

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL) {

    retval = _ecs_gzseek(f->gzptr, (ecs_z_off_t)offset, _whence);

    if (retval != 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error setting position in file \"%s\":\n\n  %s"),
                f->name, _ecs_file_error_string(f));
  }

#endif

  return retval;
}

/*!
 * \brief Return a file's name.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return pointer to file's name.
 */

const char *
ecs_file_get_name(const ecs_file_t  *f)
{
  assert(f != NULL);

  return f->name;
}

/*!
 * \brief Return a file's type.
 *
 * \param [in] f ecs_file_t descriptor.
 *
 * \return file's type.
 */

ecs_file_type_t
ecs_file_get_type(const ecs_file_t  *f)
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
 * \param [in, out ] f    ecs_file_t descriptor.
 * \param [in]       type text, binary, or Fortran binary type descriptor.
 */

void
ecs_file_set_type(ecs_file_t             *f,
                  const ecs_file_type_t   type)
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
 * \param [in, out] f ecs_file_t descriptor.
 */

void
ecs_file_set_big_endian(ecs_file_t  *f)
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
 * \param [in] f ecs_file_t descriptor.
 *
 * \return 0 if file's endianness is the same as the system's, 1 otherwise.
 */

int
ecs_file_get_swap_endian(const ecs_file_t  *f)
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
 * \param [in] f    ecs_file_t descriptor.
 * \param [in] swap 1 if bytes must be swapped, 0 otherwise.
 */

void
ecs_file_set_swap_endian(ecs_file_t  *f,
                         const int    swap)
{
  assert(f != NULL);

  f->swp_endian = swap;
}

/*!
 * \brief Test a file's error or EOF condition.
 *
 * \param [in]      f    ecs_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return 0 if no error, system error code, or -1 if EOF.
 */

int
ecs_file_read_check_error(const ecs_file_t  *f,
                          const int          line)
{
  int retval = 0;

  /* Check for possible error */

  if (f->ptr != NULL)
    retval = ferror(f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    gzerror(f->gzptr, &retval);

#endif /* defined(HAVE_ZLIB) */

  /* If we have a read error */

  if (retval != 0) {
    if (line > 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading line %d of file \"%s\":\n\n"
                  "  %s"), line, f->name, _ecs_file_error_string(f));
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("Error reading file \"%s\":\n\n"
                  "  %s"), f->name, _ecs_file_error_string(f));
    return retval;
  }

  /* Check for EOF condition */

  if (f->ptr != NULL)
    retval = feof(f->ptr);

#if defined(HAVE_ZLIB)

  else if (f->gzptr != NULL)
    retval = gzeof(f->gzptr);

#endif /* defined(HAVE_ZLIB) */

  if (retval != 0) {
    if (line > 0)
      ecs_error(__FILE__, __LINE__, 0,
                _("Premature end of file \"%s\" at line %d\n\n"),
                f->name, line);
    else
      ecs_error(__FILE__, __LINE__, 0,
                _("Premature end of file \"%s\"\n\n"), f->name);
    retval = -1;
  }

  return retval;
}

/*!
 * \brief Formatted input from a text file (as fgets()).
 *
 * \param [out]     s    buffer to which string is to be read.
 * \param [in]      size maximum number of characters to be read plus one.
 * \param [in]      f    ecs_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */

char *
ecs_file_gets(char              *s,
              const int          size,
              const ecs_file_t  *f,
              int               *line)
{
  return _ecs_file_gets(s, size, f, line, 0);
}

/*!
 * \brief Formatted input from a text file if possible (as fgets()).
 *
 * This function is similar to ecs_file_gets(), but failure to read
 * a line due to an end-of-file condition is not considered an error with
 * this variant, which may be used to read text files or sections thereof
 * of unknown length.
 *
 * \param [out]     s    buffer to which string is to be read.
 * \param [in]      size maximum number of characters to be read plus one.
 * \param [in]      f    ecs_file_t descriptor.
 * \param [in, out] line file line number if available, or NULL.
 *
 * \return s on success, NULL on error or when end of file occurs and
 *         no characters have been read.
 */

char *
ecs_file_gets_try(char              *s,
                  const int          size,
                  const ecs_file_t  *f,
                  int               *line)
{
  return _ecs_file_gets(s, size, f, line, 1);
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
 * \param [in]  f    ecs_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully read; for a Fortran
 *         record, if the whole record could not be read, returns 0.
 */

size_t
ecs_file_read(void              *rec,
              const size_t       size,
              const size_t       ni,
              const ecs_file_t  *f)
{
  return _ecs_file_read(rec, size, ni, f, 0);
}

/*!
 * \brief Read a binary C or Fortran type record.
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
 * \param [out] rec  pointer to location receiving data.
 * \param [in]  size size of each item of data in bytes.
 * \param [in]  ni   number of items to read.
 * \param [in]  f    ecs_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully read; for a Fortran
 *         record, if the whole record could not be read, returns 0.
 */

size_t
ecs_file_read_try(void              *rec,
                  const size_t       size,
                  const size_t       ni,
                  const ecs_file_t  *f)
{
  return _ecs_file_read(rec, size, ni, f, 1);
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
 * \param [in] f    ecs_file_t descriptor.
 *
 * \return the number of items (not bytes) sucessfully written.
 */

size_t
ecs_file_write(const void        *rec,
               const size_t       size,
               const size_t       ni,
               const ecs_file_t  *f)
{
  int32_t  n_bytes;
  size_t   retval;
  size_t   rec_size;

  assert(sizeof(int32_t) == 4);

  /* Check file state */

  assert(f != NULL);
  assert(rec != NULL || ni == 0);
  assert(   f->type == ECS_FILE_TYPE_BINARY
         || f->type == ECS_FILE_TYPE_FORTRAN_BINARY);
  assert(f->mode == ECS_FILE_MODE_APPEND || f->mode == ECS_FILE_MODE_WRITE);

  if (f->ptr == NULL)
    ecs_error(__FILE__, __LINE__, 0,
              _("Error writing to closed file \"%s\""), f->name);

  /* Number of bytes of record to write */

  rec_size = size * ni;

  /* In Fortran binary case, write record header */

  if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY) {

    /* Check that 4 bytes is enough for record */

    n_bytes = (int32_t)rec_size;

    if ((size_t)n_bytes != rec_size) {
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_write_error),
                f->name, _(_ecs_file_str_f_rec_too_large));
      return 0;
    }

    if (f->swp_endian == 1)
      ecs_file_swap_endian((void *)&n_bytes, (void *)&n_bytes,
                           sizeof(int32_t), 1) ;

    if (fwrite((void *)(&n_bytes), sizeof(int32_t), 1, f->ptr) != 1) {
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_write_error),
                f->name, _ecs_file_error_string(f));
      return 0;
    }

  }

  /* Write the record proper (C or Fortran) */

  if (f->swp_endian == 1 && size > 1) {

    void    *buf;

    ECS_MALLOC(buf, rec_size, unsigned char);

    ecs_file_swap_endian(buf, rec, size, ni);

    retval = fwrite((void *)buf, (size_t)size, (size_t)ni, f->ptr);

    ECS_FREE(buf);

  }
  else

    retval = fwrite((const void *)rec, (size_t)size, (size_t)ni, f->ptr);

  if (retval != (size_t)ni) {
    if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY)
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_write_error),
                f->name, _ecs_file_error_string(f));
    else
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_b_write_error),
                f->name, _ecs_file_error_string(f));
    return retval;
  }

  /* In Fortran binary case, write record footer */

  if (f->type == ECS_FILE_TYPE_FORTRAN_BINARY) {

    if (fwrite((void *)(&n_bytes), sizeof(int32_t), 1, f->ptr) != 1) {
      ecs_error(__FILE__, __LINE__, 0, _(_ecs_file_str_f_write_error),
                f->name, _ecs_file_error_string(f));
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
ecs_file_swap_endian(void          *dest,
                     const void    *src,
                     const size_t   size,
                     const size_t   ni)
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
 *          upon ECS configuration, 1 is returned.
 */

int
ecs_file_mkdir_default(const char  *pathname)
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
        ecs_error(__FILE__, __LINE__, 0, _(str_fail),
                  pathname,
                  _("  A similarly named file or directory exists "
                    "and its status is\n  not available."));
      else if (S_ISDIR(buf.st_mode) != 1)
        ecs_error(__FILE__, __LINE__, 0, _(str_fail),
                  pathname,
                  _("  A similarly named file exists and is "
                    "not a directory."));
      else
        return 0;

#endif

      errno = EEXIST; /* In case modified by stat() */

    }
    else {
      ecs_error(__FILE__, __LINE__, errno, _(str_fail),
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
ecs_file_isreg(const char  *name)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(name, &s) != 0) {
    if (errno != ENOENT)
      ecs_error(__FILE__, __LINE__, errno,
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
ecs_file_isdir(const char  *name)
{
  int retval = 0;

#if defined(HAVE_SYS_STAT_H)

  struct stat s;

  if (stat(name, &s) != 0) {
    if (errno != ENOENT)
      ecs_error(__FILE__, __LINE__, errno,
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
 * \brief Indicate Zlib version available at run time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and run time in case we use dynamic libraries.
 *
 * \return pointer to string indicating Zlib version in use, or NULL
 *         if Zlib support is not available.
 */

const char *
ecs_file_version_zlib(void)
{
#if defined(HAVE_ZLIB)
  return zlibVersion();
#else
  return NULL;
#endif
}

/*!
 * \brief Indicate Zlib version available at compilation time.
 *
 * It may be useful to compare the Zlib version used at compile
 * and link time in case we use dynamic libraries.
 *
 * \return pointer to string indicating Zlib version at compilation, or NULL
 *         if Zlib support is not available.
 */

const char *
ecs_file_version_build_zlib(void)
{
#if defined(HAVE_ZLIB)
#if defined(ZLIB_VERSION)
  return ZLIB_VERSION;
#else
  return _("unknown");
#endif
#else
  return NULL;
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
