/*============================================================================
 *  Dump of Kernel I/O file for code_saturne
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

#define CS_IGNORE_MPI 1  /* No MPI for this application */

#include "cs_defs.h"

/*
  Force LARGEFILE_SOURCE if large files enabled under 32-bit Linux or Blue Gene
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
#include <locale.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure used to index cs_io_file contents when reading */
/*----------------------------------------------------------*/

typedef struct {

  size_t          size;              /* Current number of entries */
  size_t          max_size;          /* Maximum number of entries */

  /* For each entry, we need 8 values, which we store in h_vals :
   *   0: number of values in section
   *   1: location_id
   *   2: index id
   *   3: number  of values per location
   *   4: index of section name in names array
   *   5: index of type name in types array
   *   6: index of embedded data in data array + 1 if data is
   *      embedded, 0 otherwise
   *   7: associated file id (in case of multiple files)
   */

  long long      *h_vals;            /* Base values associated
                                        with each header */

  long long      *offset;            /* Position of associated data
                                        in file (-1 if embedded) */

  size_t          max_names_size;    /* Maximum size of names array */
  size_t          names_size;        /* Current size of names array */
  char           *names;             /* Array containing section names */

  size_t          max_types_size;    /* Maximum size of type names array */
  size_t          types_size;        /* Current size of type names array */
  char           *types;             /* Array containing section type names */

  size_t          max_data_size;     /* Maximum size of embedded data array */
  size_t          data_size;         /* Current size of data array */

  unsigned char  *data;              /* Array containing embedded data */

} _cs_io_sec_index_t;

/* Main kernel IO state structure */
/*--------------------------------*/

/* File descriptor */

typedef struct {

  const char     *filename;       /* File name */
  FILE           *f;              /* File handle */

  size_t          header_size;    /* Header default size */
  size_t          header_align;   /* Header alignment */
  size_t          body_align;     /* Body alignment */

  size_t          buffer_size;    /* Current size of header buffer */
  unsigned char  *buffer;         /* Header buffer */

  size_t          n_vals;         /* Number of values in section header */
  size_t          location_id;    /* Optional value location id (0 for none) */
  size_t          index_id;       /* Optional index id (0 for none) */
  size_t          n_loc_vals;     /* Optional, number of values per location */
  size_t          type_size;      /* Size of current type */
  const char     *name;           /* Pointer to name field in section header */
  const char     *type_name;      /* Pointer to type field in section header */
  void           *data;           /* Pointer to data in section header */

  long long       offset;         /* Current position in file */
  int             swap_endian;    /* Swap big-endian and little-endian ? */

  _cs_io_sec_index_t  *index;     /* Optional section index (on read) */

} _cs_io_t;

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls _mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define MEM_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) _mem_malloc(_ni, sizeof(_type), \
                             #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls _mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define MEM_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) _mem_realloc(_ptr, _ni, sizeof(_type), \
                              #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#define MEM_FREE(_ptr) \
free(_ptr), _ptr = NULL

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print usage and exit.
 *
 * parameters:
 *   arg_0     <-- name of executable as given by argv[0]
 *   exit_code <-- EXIT_SUCCESS or EXIT_FAILURE
 *----------------------------------------------------------------------------*/

static void
_usage(const char  *arg_0,
       int          exit_code)
{
  printf
    (_("\n"
       "Usage: %s [options] <file_name>\n"
       "   or: %s -d [options] <file_name_1> <file_name_2>\n\n"
       "Dump headers and optionnaly content of a code_saturne Preprocessor,\n"
       "Partitioner, or restart file, or compare two such files.\n\n"
       "Options:\n\n"
       "  -d, --diff        diff mode.\n\n"
       "  -e, --extract     extract mode (extract full section data, with\n"
       "                    no metadata).\n\n"
       "  --f-format <fmt>  define format for floating-point numbers (default:\n"
       "                    \"15.9e\" for floats, \"22.15e\" for doubles).\n"
       "  --location <id>   only output section(s) with given location id.\n\n"
       "  -n <level>        number of first and last elements of each section\n"
       "                    to output (default: print headers only).\n\n"
       "  --section <name>  only consider section matching given criteria.\n\n"
       "  --threshold <val> in diff mode, real value above which a difference is\n"
       "                    considered significant (default: 1e-30).\n\n"
       "  -h, --help        this message.\n\n"),
     arg_0, arg_0);

  exit(exit_code);
}

/*----------------------------------------------------------------------------
 * Abort with error message.
 *
 * parameters:
 *   file_name      <-- name of source file from which function is called.
 *   line_num       <-- line of source file from which function is called.
 *   sys_error_code <-- error code if error in system or libc call,
 *                      0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   ...            <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

static void
_error(const char  *file_name,
       int          line_num,
       int          sys_error_code,
       const char  *format,
       ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  fflush(stdout);

  fprintf(stderr, "\n");

  if (sys_error_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_error_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  va_end(arg_ptr);

  assert(0);

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Allocate memory and check result.
 *
 * This function simply wraps malloc() and calls _error() if it fails.
 *
 * parameters:
 *   ni        <-- number of elements.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 *----------------------------------------------------------------------------*/

static void *
_mem_malloc(size_t       ni,
            size_t       size,
            const char  *var_name,
            const char  *file_name,
            int          line_num)
{
  void  *p_ret;
  size_t  alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  /* Allocate memory and check return */

  p_ret = malloc(alloc_size);

  if (p_ret == NULL)
    _error(file_name, line_num, errno,
           _("Failure to allocate \"%s\" (%lu bytes)"),
           var_name, (unsigned long)alloc_size);

  return p_ret;
}

/*----------------------------------------------------------------------------
 * Allocate memory and check result.
 *
 * This function simply wraps realloc() and calls _error() if it fails.
 *
 * parameters:
 *   ptr       <-- pointer to previous memory location
 *   ni        <-- number of elements.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to reallocated memory.
 *----------------------------------------------------------------------------*/

static void *
_mem_realloc(void        *ptr,
             size_t       ni,
             size_t       size,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
  void  *p_ret;

  size_t realloc_size = ni * size;

  p_ret = realloc(ptr, realloc_size);

  if (size != 0 && p_ret == NULL)
    _error(file_name, line_num, errno,
           _("Failure to reallocate \"%s\" (%lu bytes)"),
           var_name, (unsigned long)realloc_size);

  return p_ret;
}

/*----------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * parameters:
 *   buf  <-> pointer to converted data location.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of data items.
 *----------------------------------------------------------------------------*/

static void
_swap_endian(void        *buf,
             size_t       size,
             size_t       ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)buf;
  const unsigned char  *psrc = (const unsigned char *)buf;

  for (i = 0; i < ni; i++) {

    shift = i * size;

    for (ib = 0; ib < (size / 2); ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }
  }
}

/*----------------------------------------------------------------------------
 * Read data to a buffer using standard C IO.
 *
 * parameters:
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *   inp  <-> input file descriptor
 *
 * returns:
 *   the number of items (not bytes) sucessfully read;
 *----------------------------------------------------------------------------*/

static size_t
_file_read(void        *buf,
           size_t       size,
           size_t       ni,
           _cs_io_t    *inp)
{
  size_t retval = 0;

  assert(inp->f != NULL);

  if (ni != 0)
    retval = fread(buf, size, ni, inp->f);

  inp->offset += size*ni;

  /* In case of error, determine error type */

  if (retval != ni) {
    int err_num = ferror(inp->f);
    if (err_num != 0)
      _error(__FILE__, __LINE__, 0,
             _("Error reading file \"%s\":\n\n  %s"),
             inp->filename, strerror(err_num));
    else if (feof(inp->f) != 0)
      _error(__FILE__, __LINE__, 0,
             _("Premature end of file \"%s\""), inp->filename);
    else
      _error(__FILE__, __LINE__, 0,
             _("Error reading file \"%s\""), inp->filename);
  }

  if (inp->swap_endian)
    _swap_endian(buf, size, ni);

  return retval;
}

/*----------------------------------------------------------------------------
 * Update the file pointer according to whence.
 *
 * parameters:
 *   inp    <-- input file descriptor
 *   offset <-- add to position specified to whence to obtain new position,
 *              measured in characters from the beginning of the file
 *   whence <-- beginning if SEEK_SET, current if SEEK_CUR, or end-of-file
 *              if SEEK_END
 *
 * returns:
 *   0 upon success, nonzero otherwise; currently, errors are fatal.
 *----------------------------------------------------------------------------*/

static int
_file_seek(_cs_io_t   *inp,
           long long   offset,
           int         whence)
{
  int retval = 0;

  assert(inp != NULL);

  if (inp->f != NULL) {

#if (SIZEOF_LONG < 8)

    /* For 32-bit systems, large file support may be necessary */

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)

    retval = fseeko(inp->f, (off_t)offset, whence);

    if (retval != 0)
      _error(__FILE__, __LINE__, errno,
             _("Error setting position in file \"%s\":\n\n  %s"),
             inp->filename, strerror(errno));
# else

    /* Test if offset larger than allowed */

    long _offset = offset;

    if (_offset == offset) {
      retval = fseek(inp->f, (long)offset, whence);
      if (retval != 0)
        _error(__FILE__, __LINE__, errno,
               _("Error setting position in file \"%s\":\n\n  %s"),
               inp->filename, strerror(errno));
    }
    else {
      retval = -1;
      _error
        (__FILE__, __LINE__, 0,
         _("Error setting position in file \"%s\":\n\n  %s"),
         inp->filename,
         _("sizeof(off_t) > sizeof(long) but fseeko() not available"));
    }

# endif /* defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64) */

#else /* SIZEOF_LONG >= 8 */

    /* For 64-bit systems, standard fseek should be enough */

    retval = fseek(inp->f, (long)offset, whence);
    if (retval != 0)
      _error(__FILE__, __LINE__, errno,
             _("Error setting position in file \"%s\":\n\n  %s"),
             inp->filename, strerror(errno));

#endif /* SIZEOF_LONG */
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the position of the file pointer.
 *
 * parameters:
 *   inp <-- input file descriptor
 *
 * returns:
 *   current position of the file pointer
 *----------------------------------------------------------------------------*/

static long long
_file_tell(_cs_io_t  *inp)
{
  long long offset = 0;
  assert(inp != NULL);

  if (inp->f != NULL) {

    /* For 32-bit systems, large file support may be necessary */

#if (SIZEOF_LONG < 8)

# if defined(HAVE_FSEEKO) && (_FILE_OFFSET_BITS == 64)
    offset = ftello(inp->f);
# else
    /*
      Without ftello, ftell will fail above 2 Gigabytes, in which case
      offset == -1 and errno == EOVERFLOW, but should work on smaller
      files. We prefer not to be too strict about fseeko availability, as
      the only 32-bit case without ftello we have encountered is Cygwin
      (for which ftello requires additional non-default libraries), which
      is expected to be used mainly for small cases.
    */
    offset = ftell(inp->f);
# endif

    /* For 64-bit systems, standard ftell should be enough */

#else /* SIZEOF_LONG >= 8 */
    offset = ftell(inp->f);
#endif

  }

  if (offset < 0)
    _error(__FILE__, __LINE__, 0,
           _("Error obtaining position in file \"%s\":\n\n  %s"),
           inp->filename, strerror(errno));

  return offset;
}

/*----------------------------------------------------------------------------
 * Convert a buffer of type uint64_t to size_t
 *
 * parameters:
 *   f   <-- pointer to file
 *   val --> array to which values are read
 *   n   <-- number of values to read
 *----------------------------------------------------------------------------*/

static void
_convert_size(unsigned char  buf[],
              size_t         val[],
              size_t         n)
{
  size_t i;

#if (__STDC_VERSION__ >= 199901L)

  for (i = 0; i < n; i++)
    val[i] = ((uint64_t *)buf)[i];

#else

  if (sizeof(size_t) == 8) {
    for (i = 0; i < n; i++)
      val[i] = ((size_t *)buf)[i];
  }
  else if (sizeof(unsigned long long) == 8) {
    for (i = 0; i < n; i++)
      val[i] = ((unsigned long long *)buf)[i];
  }
  else
    _error(__FILE__, __LINE__, 0,
           "Compilation configuration / porting error:\n"
           "Unable to determine a 64-bit unsigned int type.\n"
           "size_t is %d bits, unsigned long long %d bits",
           sizeof(size_t)*8, sizeof(unsigned long long)*8);

#endif
}

/*----------------------------------------------------------------------------
 * Open input file, testing magic string for type.
 *
 * parameters:
 *   filename <-- file name
 *   mode     <-- 0 for dump, 1 for extract, 2 for diff
 *
 * returns:
 *   File metadata structure
 *----------------------------------------------------------------------------*/

static _cs_io_t
_open_input(const char  *filename,
            int          mode)
{
  unsigned int_endian;
  size_t alignments[3];
  char header_buf[65];

  _cs_io_t inp;

  inp.filename = filename;
  inp.f = NULL;
  inp.header_size = 0;
  inp.header_align = 0;
  inp.body_align = 0;
  inp.buffer_size = 0;
  inp.buffer = NULL;

  inp.n_vals = 0;
  inp.location_id = 0;
  inp.index_id = 0,
  inp.n_loc_vals = 0;
  inp.type_size = 0;
  inp.name = NULL;
  inp.type_name = NULL;
  inp.data = NULL;

  inp.offset = 0;
  inp.swap_endian = 0;

  inp.index = NULL;

  /* Check if system is "big-endian" or "little-endian" */

  int_endian = 0;
  *((char *)(&int_endian)) = '\1';
  if (int_endian == 1)
    inp.swap_endian = 1;

  /* Open file */

  if (mode != 1)
    printf(_("\nOpening input file: \"%s\"\n\n"), filename) ;

  fflush(stdout);

  inp.f = fopen(inp.filename, "rb");

  if (inp.f == NULL)
    _error(__FILE__, __LINE__, 0,
           _("Error opening file \"%s\":\n\n"
             "  %s"), inp.filename, strerror(errno));

  /* Read "magic string" */

  _file_read(header_buf, 1, 64, &inp);

  if (strncmp(header_buf, "Code_Saturne I/O, BE, R0", 64) != 0) {
    header_buf[64] = '\0';
    _error(__FILE__, __LINE__, 0,
           _("File format of \"%s\" is not recognized:\n"
             "First %d bytes: \"%s\"."),
           filename, 64, header_buf);
  }

  _file_read(header_buf, 1, 64, &inp);

  header_buf[64] = '\0';
  if (mode != 1)
    printf(_("  File type: %s\n"), header_buf);

  _file_read(header_buf, 8, 3, &inp);

  _convert_size((unsigned char*)header_buf, alignments, 3);

  inp.header_size = alignments[0];
  inp.header_align = alignments[1];
  inp.body_align = alignments[2];

  if (mode != 1)
    printf(_("\n"
             "  Base header size: %d\n"
             "  Header alignment: %d\n"
             "  Body alignment:   %d\n"),
           (int)(inp.header_size),
           (int)(inp.header_align),
           (int)(inp.body_align));

  inp.buffer_size = inp.header_size;
  MEM_MALLOC(inp.buffer, inp.buffer_size, unsigned char);

  /* Finish */

  return inp;
}

/*----------------------------------------------------------------------------
 * Close input file
 *
 * parameters:
 *   f    <-> pointer to file object
 *   mode <-- 0 for dump, 1 for extract, 2 for diff
 *----------------------------------------------------------------------------*/

static void
_close_input(_cs_io_t  *inp,
             int        mode)
{
  if (inp != NULL) {
    if (inp->f != NULL) {
      int retval = 0;
      if (mode == 0)
        printf(_("\nClosing input: \"%s\"\n\n"), inp->filename);
      retval = fclose(inp->f);
      inp->f = NULL;
      if (retval != 0)
        _error(__FILE__, __LINE__, 0,
               _("Error closing file \"%s\":\n\n"
                 "  %s"), inp->filename, strerror(errno));
    }
    inp->header_size = 0;
    inp->header_align = 0;
    inp->body_align = 0;
    inp->buffer_size = 0;
    inp->filename = NULL;
    MEM_FREE(inp->buffer);
  }
}

/*----------------------------------------------------------------------------
 * Convert an argument to an integer and check its validity
 *
 * parameters:
 *   arg_id  <-- index of argument in argv
 *   argc    <-- number of command line arguments
 *   argv    <-- array of command line arguments
 *
 * returns:
 *   integer value
 *----------------------------------------------------------------------------*/

#if (__STDC_VERSION__ >= 199901L)
static long long
_arg_to_int(int    arg_id,
            int    argc,
            char  *argv[])
{
  char  *start = NULL;
  char  *end = NULL;
  long long  retval = 0;

  if (arg_id < argc) {
    start = argv[arg_id];
    end = start + strlen(start);
    retval = strtoll(start, &end, 0);
    if (end != start + strlen(start))
      _usage(argv[0], EXIT_FAILURE);
  }
  else
    _usage(argv[0], EXIT_FAILURE);

  return retval;
}

#else /* (__STDC_VERSION__ == 1989) */

static long
_arg_to_int(int    arg_id,
            int    argc,
            char  *argv[])
{
  char  *start = NULL;
  char  *end = NULL;
  long  retval = 0;

  if (arg_id < argc) {
    start = argv[arg_id];
    end = start + strlen(start);
    retval = strtol(start, &end, 0);
    if (end != start + strlen(start))
      _usage(argv[0], EXIT_FAILURE);
  }
  else
    _usage(argv[0], EXIT_FAILURE);

  return retval;
}

#endif /* (__STDC_VERSION__) */

/*----------------------------------------------------------------------------
 * Read command line arguments.
 *
 * parameters:
 *   argc             <-- number of command line arguments
 *   argv             <-- array of command line arguments
 *   mode             --> 0 for dump, 1 for extract, 2 for diff
 *   echo             --> echo (verbosity) level
 *   location_id      --> if >= 0 location id filter
 *   sec_name_arg_id  --> if > 0, id of command line argument defining
 *                        section name filter
 *   f_fmt_arg_id     --> if > 0, id of command line argument defining format
 *                        for output of floating-point values
 *   threshold        --> threshold above which 2 floating-point values are
 *                        considered different.
 *   file_name_arg_id --> index of command line arguments defining file names
 *----------------------------------------------------------------------------*/

static void
_read_args(int          argc,
           char       **argv,
           int         *mode,
           size_t      *echo,
           int         *location_id,
           int         *sec_name_arg_id,
           int         *f_fmt_arg_id,
           double      *threshold,
           int          file_name_arg_id[2])
{
  int i = 1;
  int n_files = 0;

  /* Initialize return arguments */

  *echo = 0;
  *mode = 0;
  *location_id = -1;
  *sec_name_arg_id = 0;
  *f_fmt_arg_id = 0;
  *threshold = 1.e-30;
  *file_name_arg_id = 0;

  /* Parse and check command line */

  if (argc < 2)
    _usage(argv[0], EXIT_FAILURE);

  while (i < argc) {

    if (strcmp(argv[i], "--f-format") == 0) {
      i++;
      if (i < argc)
        *f_fmt_arg_id = i;
      else
        _usage(argv[0], EXIT_FAILURE);
    }

    else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--diff") == 0) {
      if (*mode == 0)
        *mode = 2;
      else
        _usage(argv[0], EXIT_SUCCESS);
    }

    else if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--extract") == 0) {
      if (*mode == 0)
        *mode = 1;
      else
        _usage(argv[0], EXIT_SUCCESS);
    }

    else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
      _usage(argv[0], EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      *echo = _arg_to_int(i, argc, argv);
    }

    else if (strcmp(argv[i], "--location") == 0) {
      i++;
      *location_id = _arg_to_int(i, argc, argv);
    }

    else if (strcmp(argv[i], "--section") == 0) {
      i++;
      if (i < argc)
        *sec_name_arg_id = i;
      else
        _usage(argv[0], EXIT_FAILURE);
    }

    else if (strcmp(argv[i], "--threshold") == 0) {
      i++;
      if (i < argc) {
        char *end = argv[i] + strlen(argv[i]);
        *threshold = strtod(argv[i], &end);
        if (end == argv[i])
          _usage(argv[0], EXIT_FAILURE);
      }
      else
        _usage(argv[0], EXIT_FAILURE);
    }

    else {

      if (n_files == 0 || (n_files == 1 && *mode == 2))
        file_name_arg_id[n_files++] = i;
      else
        _usage(argv[0], EXIT_FAILURE);

    }

    i++;
  }

  if ((*mode != 2 && n_files != 1) || (*mode == 2 && n_files != 2))
    _usage(argv[0], EXIT_FAILURE);

  /* At this point, command line seems correct */

  if (*mode != 1)
    printf(_("\n"
             "  .----------------------------.\n"
             "  |   code_saturne file dump   |\n"
             "  `----------------------------'\n"));
}

/*----------------------------------------------------------------------------
 * Echo values depending on type.
 *
 * Type name should already have been checked when this function is called,
 * so no additional check is done here.
 *
 * parameters:
 *   n_values       <-- number of values to echo
 *   n_values_shift <-- shift to second (end) series of values to echo
 *   buffer         <-- pointer to data
 *   type_name      <-- name of data type
 *   f_fmt          <-- if != NULL, format for output of floating-point values
 *----------------------------------------------------------------------------*/

static void
_echo_values(size_t       n_values,
             size_t       n_values_shift,
             const void  *buffer,
             const char   type_name[],
             const char  *f_fmt)
{
  size_t i;

  /* Check type name */

  if (type_name[0] == 'c') {
    const char *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : '%c'\n",
             (unsigned long)(i + n_values_shift),
             _buffer[i]);
  }

#if (__STDC_VERSION__ >= 199901L)

  else if (type_name[0] == 'i' && type_name[1] == '4') {
    const int32_t *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %d\n",
             (unsigned long)(i + n_values_shift),
             (int)(_buffer[i]));
  }

  else if (type_name[0] == 'i' && type_name[1] == '8') {
    const int64_t *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %ld\n",
             (unsigned long)(i + n_values_shift),
             (long)(_buffer[i]));
  }

  else if (type_name[0] == 'u' && type_name[1] == '4') {
    const uint32_t *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %u\n",
             (unsigned long)(i + n_values_shift),
             (unsigned)(_buffer[i]));
  }

  else if (type_name[0] == 'u' && type_name[1] == '8') {
    const uint64_t *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %lu\n",
             (unsigned long)(i + n_values_shift),
             (unsigned long)(_buffer[i]));
  }

#else /* (__STDC_VERSION__ < 199901L) */

  else if (type_name[0] == 'i' && type_name[1] == '4') {
    if (sizeof(int) == 4) {
      const int *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %d\n",
               (unsigned long)(i + n_values_shift),
               _buffer[i]);
    }
    else if (sizeof(short) == 4) {
      const short *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %d\n",
               (unsigned long)(i + n_values_shift),
               (int)(_buffer[i]));
    }
    else
      printf("    int32_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'i' && type_name[1] == '8') {
    if (sizeof(long) == 8) {
      const long *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %ld\n",
               (unsigned long)(i + n_values_shift),
               _buffer[i]);
    }
    else if (sizeof(long long) == 8) {
      const long long *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %ld\n",
               (unsigned long)(i + n_values_shift),
               (long)(_buffer[i]));
    }
    else
      printf("    int64_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'u' && type_name[1] == '4') {
    if (sizeof(unsigned) == 4) {
      const unsigned *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %u\n",
               (unsigned long)(i + n_values_shift),
               _buffer[i]);
    }
    else if (sizeof(unsigned short) == 4) {
      const unsigned short *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %u\n",
               (unsigned long)(i + n_values_shift),
               (unsigned)(_buffer[i]));
    }
    else
      printf("    uint32_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'u' && type_name[1] == '8') {
    if (sizeof(unsigned long) == 8) {
      const unsigned long *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %lu\n",
               (unsigned long)(i + n_values_shift),
               (unsigned long)(_buffer[i]));
    }
    else if (sizeof(unsigned long long) == 8) {
      const unsigned long long *_buffer = buffer;
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %lu\n",
               (unsigned long)(i + n_values_shift),
               (unsigned long)(_buffer[i]));
    }
    else
      printf("    uint64_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

#endif /* (__STDC_VERSION__) */

  else if (type_name[0] == 'r' && type_name[1] == '4') {
    const float *_buffer = buffer;
    if (f_fmt == NULL) {
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %15.9e\n",
               (unsigned long)(i + n_values_shift),
               (double)(_buffer[i]));
    }
    else {
      char format[64] = "    %10lu : %";
      strncat(format, f_fmt, 48);
      strcat(format, "\n");
      for (i = 0; i < n_values; i++)
        printf(format,
               (unsigned long)(i + n_values_shift),
               (double)(_buffer[i]));
    }
  }

  else if (type_name[0] == 'r' && type_name[1] == '8') {
    const double *_buffer = buffer;
    if (f_fmt == NULL) {
      for (i = 0; i < n_values; i++)
        printf("    %10lu : %22.15e\n",
               (unsigned long)(i + n_values_shift),
               _buffer[i]);
    }
    else {
      char format[64] = "    %10lu : %";
      strncat(format, f_fmt, 48);
      strcat(format, "\n");
      for (i = 0; i < n_values; i++)
        printf(format,
               (unsigned long)(i + n_values_shift),
               _buffer[i]);
    }
  }

}

/*----------------------------------------------------------------------------
 * Determine a type's size based on its name.
 *
 * parameters:
 *   type_name  <-- type name
 *
 * returns:
 *   size of data type
 *----------------------------------------------------------------------------*/

static size_t
_type_size_from_name(const char  *type_name)
{
  size_t type_size = 0;

  assert(type_name != NULL);

  /* Check type name and compute size of data */

  if (type_name[0] == 'c' && type_name[1] == ' ')
    type_size = 1;

  else if (   type_name[0] == 'i'
           || type_name[0] == 'u'
           || type_name[0] == 'r') {

    if (type_name[1] == '4')
      type_size = 4;
    else if (type_name[1] == '8')
      type_size = 8;
  }

  if (type_size == 0)
    _error(__FILE__, __LINE__, 0,
           _("Type \"%s\" is not known\n"
             "Known types: \"c \", \"i4\", \"i8\", \"u4\", \"u8\", "
             "\"r4\", \"r8\"."), type_name);
  return type_size;
}

/*----------------------------------------------------------------------------
 * Read section header.
 *
 * parameters:
 *   inp <-> pointer to input object
 *
 * returns:
 *   number of bytes in section body
 *----------------------------------------------------------------------------*/

static size_t
_read_section_header(_cs_io_t  *inp)
{
  size_t body_size = 0;
  size_t header_vals[6];
  unsigned int_endian = 0;

  *((char *)(&int_endian)) = '\1'; /* Determine if we are little-endian */

  assert(inp != NULL);
  assert(inp->f != NULL);
  assert(inp->buffer != NULL);

  /* Position read pointer if necessary */
  /*------------------------------------*/

  {
    long long offset = _file_tell(inp);
    size_t ha = inp->header_align;
    offset += (ha - (offset % ha)) % ha;
    _file_seek(inp, offset, SEEK_SET);
  }

  /* Read header */
  /*-------------*/

  _file_read(inp->buffer, 1, inp->header_size, inp);

  if (int_endian == 1)
    _swap_endian(inp->buffer, 8, 6);

  _convert_size(inp->buffer, header_vals, 6);

  if (header_vals[0] > inp->header_size) {

    if (header_vals[0] > inp->buffer_size) {
      while (header_vals[0] > inp->buffer_size)
        inp->buffer_size *=2;
      MEM_REALLOC(inp->buffer, inp->buffer_size, unsigned char);
    }

    _file_read(inp->buffer + inp->header_size,
               1,
               header_vals[0] - inp->header_size,
               inp);

  }

  /* Set pointers to data fields */

  inp->n_vals = header_vals[1];
  inp->location_id = header_vals[2];
  inp->index_id = header_vals[3];
  inp->n_loc_vals = header_vals[4];
  inp->type_size = 0;
  inp->data = NULL;
  inp->type_name = (char *)(inp->buffer + 48);
  inp->name = (char *)(inp->buffer + 56);

  if (header_vals[1] > 0 && inp->type_name[7] == 'e')
    inp->data = inp->buffer + 56 + header_vals[5];

  inp->type_size = 0;

  if (inp->n_vals > 0) {

    inp->type_size = _type_size_from_name(inp->type_name);

    if (inp->data == NULL)
      body_size = inp->type_size*inp->n_vals;

    else if (int_endian == 1 && inp->type_size > 1)
      _swap_endian(inp->data, inp->type_size, inp->n_vals);
  }

  return body_size;
}

/*----------------------------------------------------------------------------
 * Read section values and print associated info
 *
 * If values are already embedded in the header, no actual reading is done
 *
 * parameters:
 *   inp   <-> pointer to input object
 *   echo  <-- number of values to print
 *   f_fmt <-- if != NULL, format for output of floating-point values
 *----------------------------------------------------------------------------*/

static void
_read_section_values(_cs_io_t    *inp,
                     size_t       echo,
                     const char  *f_fmt)
{
  size_t n_print = 0, n_skip = 0;
  unsigned char  *buffer = NULL;
  const unsigned char  *data = NULL;

  assert(inp->n_vals > 0);

  if (inp->data != NULL)
    printf(_("      Values in header\n"));

  /* Compute number of values to skip */

  if (inp->n_vals > echo*2) {
    n_skip = inp->n_vals - echo*2;
    n_print = echo;
  }
  else {
    n_skip = 0;
    n_print = inp->n_vals;
  }

  /* Position read pointer if non-embedded data is present */

  if (inp->data == NULL) {

    long long offset = _file_tell(inp);
    size_t ba = inp->body_align;
    offset += (ba - (offset % ba)) % ba;
    _file_seek(inp, offset, SEEK_SET);

    /* Allocate buffer */

    if (n_print > 0) {
      MEM_MALLOC(buffer, n_print*inp->type_size, unsigned char);
      _file_read(buffer, inp->type_size, n_print, inp);
      data = buffer;
    }
  }

  else if (n_print > 0)
    data = inp->data;

  /* Print first part of data */

  if (n_print > 0) {

    if (n_skip > 0)
      printf(_("    %lu first and last elements:\n"), (unsigned long)echo);
    else
      printf(_("    elements:\n"));

    _echo_values(n_print,
                 1,
                 data,
                 inp->type_name,
                 f_fmt);

    if (n_skip > 0)
      printf("    ..........   ..........\n");

  }

  /* Move to tail of data and read if necessary */

  if (n_skip > 0) {

    if (inp->data == NULL) {

      long long offset = _file_tell(inp) + n_skip*inp->type_size;
      _file_seek(inp, offset, SEEK_SET);

      if (n_print > 0) {
        _file_read(buffer, inp->type_size, n_print, inp);
        data = buffer;
      }
    }

    else if (n_print > 0)
      data = ((unsigned char *)inp->data) + ((n_print+n_skip)*inp->type_size);

    if (n_print > 0)
      _echo_values(n_print,
                   inp->n_vals - n_print + 1,
                   data,
                   inp->type_name,
                   f_fmt);

  }

  if (buffer != NULL)
    MEM_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Skip section values
 *
 * If values are already embedded in the header, no actual reading is done
 *
 * parameters:
 *   inp   <-> pointer to input object
 *----------------------------------------------------------------------------*/

static void
_skip_section_values(_cs_io_t    *inp)
{
  assert(inp->n_vals > 0);

  /* Position read pointer if non-embedded data is present */

  if (inp->data == NULL) {
    long long offset = _file_tell(inp);
    size_t ba = inp->body_align;
    offset += (ba - (offset % ba)) % ba + inp->n_vals*inp->type_size;
    _file_seek(inp, offset, SEEK_SET);
  }
}

/*----------------------------------------------------------------------------
 * Extract section values and print associated info
 *
 * If values are already embedded in the header, no actual reading is done
 *
 * parameters:
 *   inp   <-> pointer to input object
 *   f_fmt <-- if != NULL, format for output of floating-point values
 *----------------------------------------------------------------------------*/

static void
_extract_section_values(_cs_io_t    *inp,
                        const char  *f_fmt)
{
  size_t i;
  size_t n_vals = inp->n_vals;
  void *data = NULL;
  const char *type_name = inp->type_name;

  /* Position read pointer if non-embedded data is present */

  if (inp->data == NULL) {

    _file_seek(inp, inp->offset, SEEK_SET);

    /* Allocate buffer */

    if (n_vals > 0) {
      MEM_MALLOC(data, n_vals*inp->type_size, unsigned char);
      _file_read(data, inp->type_size, n_vals, inp);
    }
  }

  else if (n_vals > 0)
    data = inp->data;

  /* Print data */

  if (n_vals > 0) {

    /* Check type name */

    if (type_name[0] == 'c') {
      if (inp->location_id == 0) {
        char *_data = NULL;
        MEM_MALLOC(_data, n_vals + 1, char);
        memcpy(_data, data, n_vals);
        for (i = 0; i < n_vals; i++)
          if (_data[i] == '\0')
            _data[i] = '\n';
        _data[n_vals] = '\0';
        printf("%s", _data);
        MEM_FREE(_data);
      }
      else {
        char *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%d\n", (int)(_data[i]));
      }
    }

#if (__STDC_VERSION__ >= 199901L)

    else if (type_name[0] == 'i' && type_name[1] == '4') {
      const int32_t *_data = data;
      for (i = 0; i < n_vals; i++)
        printf("%d\n", (int)(_data[i]));
    }

    else if (type_name[0] == 'i' && type_name[1] == '8') {
      const int64_t *_data = data;
      for (i = 0; i < n_vals; i++)
        printf("%ld\n", (long)(_data[i]));
    }

    else if (type_name[0] == 'u' && type_name[1] == '4') {
      const uint32_t *_data = data;
      for (i = 0; i < n_vals; i++)
        printf("%u\n", (unsigned)(_data[i]));
    }

    else if (type_name[0] == 'u' && type_name[1] == '8') {
      const uint64_t *_data = data;
      for (i = 0; i < n_vals; i++)
        printf("%lu\n", (unsigned long)(_data[i]));
    }

#else /* (__STDC_VERSION__ < 199901L) */

    else if (type_name[0] == 'i' && type_name[1] == '4') {
      if (sizeof(int) == 4) {
        const int *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%d\n",  _data[i]);
      }
      else if (sizeof(short) == 4) {
        const short *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%d\n", (int)(_data[i]));
      }
      else
        printf("    int32_t undefined"
               " (porting error, C99 compiler needed)\n");
    }

    else if (type_name[0] == 'i' && type_name[1] == '8') {
      if (sizeof(long) == 8) {
        const long *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%ld\n", _data[i]);
      }
      else if (sizeof(long long) == 8) {
        const long long *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%ld\n", (long)(_data[i]));
      }
      else
        printf("    int64_t undefined"
               " (porting error, C99 compiler needed)\n");
    }

    else if (type_name[0] == 'u' && type_name[1] == '4') {
      if (sizeof(unsigned) == 4) {
        const unsigned *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%u\n", _data[i]);
      }
      else if (sizeof(unsigned short) == 4) {
        const unsigned short *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%u\n", (unsigned)(_data[i]));
      }
      else
        printf("    uint32_t undefined"
               " (porting error, C99 compiler needed)\n");
    }

    else if (type_name[0] == 'u' && type_name[1] == '8') {
      if (sizeof(unsigned long) == 8) {
        const unsigned long *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%lu\n", (unsigned long)(_data[i]));
      }
      else if (sizeof(unsigned long long) == 8) {
        const unsigned long long *_data = data;
        for (i = 0; i < n_vals; i++)
          printf("%lu\n", (unsigned long)(_data[i]));
      }
      else
        printf("    uint64_t undefined"
               " (porting error, C99 compiler needed)\n");
    }

#endif /* (__STDC_VERSION__) */

    else if (type_name[0] == 'r' && type_name[1] == '4') {
      const float *_data = data;
      if (f_fmt == NULL) {
        for (i = 0; i < n_vals; i++)
          printf("%15.9e\n", (double)(_data[i]));
      }
      else {
        char format[64] = "%";
        strncat(format, f_fmt, 48);
        strcat(format, "\n");
        for (i = 0; i < n_vals; i++)
          printf(format, (double)(_data[i]));
      }
    }

    else if (type_name[0] == 'r' && type_name[1] == '8') {
      const double *_data = data;
      if (f_fmt == NULL) {
        for (i = 0; i < n_vals; i++)
          printf("%22.15e\n", _data[i]);
      }
      else {
        char format[64] = "%";
        strncat(format, f_fmt, 48);
        strcat(format, "\n");
        for (i = 0; i < n_vals; i++)
          printf(format, _data[i]);
      }
    }
  }

  if (inp->data == NULL)
    MEM_FREE(data);
}

/*----------------------------------------------------------------------------
 * Read section.
 *
 * parameters:
 *   inp         <-> pointer to input object
 *   echo        <-- number of values to print
 *   location_id <-- if >= 0 location id filter
 *   sec_name    <-- if != NULL, section name filter
 *   f_fmt       <-- if != NULL, format for output of floating-point values
 *----------------------------------------------------------------------------*/

static void
_read_section(_cs_io_t    *inp,
              int          echo,
              int          location_id,
              const char  *sec_name,
              const char  *f_fmt)
{
  int read_section = 0;

  assert(inp != NULL);
  assert(inp->f != NULL);

  /* Read section header and print basic information */

  _read_section_header(inp);

  if (   (location_id < 0 || (unsigned)location_id == inp->location_id)
      && (sec_name == NULL || !strcmp(sec_name, inp->name)))
    read_section = 1;

  if (read_section) {

    printf(_("\n"
             "  Section:                \"%s\"\n"
             "    Number of values:      %lu\n"),
           inp->name, (unsigned long)(inp->n_vals));

    if (inp->n_vals > 0)
      printf(_("    Type:                 \"%s\"\n"), inp->type_name);

    printf(_("      Location id:         %lu\n"
             "      Index id:            %lu\n"
             "      Values per location: %lu\n"),
           (unsigned long)(inp->location_id),
           (unsigned long)(inp->index_id),
           (unsigned long)(inp->n_loc_vals));
  }

  if (inp->n_vals > 0) {
    if (read_section)
      _read_section_values(inp, echo, f_fmt);
    else
      _skip_section_values(inp);
  }
}

/*----------------------------------------------------------------------------
 * Update an index structure with info from the last header read
 *
 * Also sets the file position for the next read
 *
 * parameters:
 *   inp    <-> input kernel IO structure
 *----------------------------------------------------------------------------*/

static void
_update_index_and_shift(_cs_io_t  *inp)
{
  size_t id = 0;
  size_t new_names_size = 0;
  size_t new_types_size = 0;
  size_t new_data_size = 0;

  _cs_io_sec_index_t  *idx = inp->index;

  if (idx == NULL)
    return;

  /* Reallocate if necessary */

  if (idx->size + 1 == idx->max_size) {
    if (idx->max_size == 0)
      idx->max_size = 32;
    else
      idx->max_size *= 2;
    MEM_REALLOC(idx->h_vals, idx->max_size*8, long long);
    MEM_REALLOC(idx->offset, idx->max_size, long long);
  };

  new_names_size = idx->names_size + strlen(inp->name) + 1;
  new_types_size = idx->types_size + strlen(inp->type_name) + 1;

  if (inp->data != NULL)
    new_data_size = idx->data_size + (inp->n_vals * inp->type_size);

  if (new_names_size > idx->max_names_size) {
    if (idx->max_names_size == 0)
      idx->max_names_size = 128;
    while (new_names_size > idx->max_names_size)
      idx->max_names_size *= 2;
    MEM_REALLOC(idx->names, idx->max_names_size, char);
  }

  if (new_types_size > idx->max_types_size) {
    if (idx->max_types_size == 0)
      idx->max_types_size = 128;
    while (new_types_size > idx->max_types_size)
      idx->max_types_size *= 2;
    MEM_REALLOC(idx->types, idx->max_types_size, char);
  }

  if (new_data_size > idx->max_data_size) {
    if (idx->max_data_size == 0)
      idx->max_data_size = 128;
    while (new_data_size > idx->max_data_size)
      idx->max_data_size *= 2;
    MEM_REALLOC(idx->data, idx->max_data_size, unsigned char);
  }

  /* Set values */

  id = idx->size;

  idx->h_vals[id*8]     = inp->n_vals;
  idx->h_vals[id*8 + 1] = inp->location_id;
  idx->h_vals[id*8 + 2] = inp->index_id;
  idx->h_vals[id*8 + 3] = inp->n_loc_vals;
  idx->h_vals[id*8 + 4] = idx->names_size;
  idx->h_vals[id*8 + 5] = idx->types_size;
  idx->h_vals[id*8 + 6] = 0;

  strcpy(idx->names + idx->names_size, inp->name);
  idx->names[new_names_size - 1] = '\0';
  idx->names_size = new_names_size;

  strcpy(idx->types + idx->types_size, inp->type_name);
  idx->types[new_types_size - 1] = '\0';
  idx->types_size = new_types_size;

  if (inp->data == NULL) {
    long long offset = _file_tell(inp);
    long long data_shift = inp->n_vals * inp->type_size;
    if (inp->body_align > 0) {
      size_t ba = inp->body_align;
      idx->offset[id] = offset + (ba - (offset % ba)) % ba;
    }
    else
      idx->offset[id] = offset;
    _file_seek(inp, idx->offset[id] + data_shift, SEEK_SET);
  }
  else {
    idx->h_vals[id*8 + 6] = idx->data_size + 1;
    memcpy(idx->data + idx->data_size,
           inp->data,
           new_data_size - idx->data_size);
    idx->data_size = new_data_size;
    idx->offset[id] = -1;
  }

  idx->size += 1;
}

/*----------------------------------------------------------------------------
 * Create an index structure to a _cs_io_t structure.
 *
 * parameters:
 *   inp        <-> pointer to cs_io_t structure
 *   end_offset <-- file size
 *----------------------------------------------------------------------------*/

static void
_create_index(_cs_io_t  *inp,
             long long   end_offset)
{
  _cs_io_sec_index_t  *idx = NULL;

  MEM_MALLOC(idx, 1, _cs_io_sec_index_t);

  /* Set structure fields */

  idx->size = 0;
  idx->max_size = 32;

  MEM_MALLOC(idx->h_vals, idx->max_size*8, long long);
  MEM_MALLOC(idx->offset, idx->max_size, long long);

  idx->max_names_size = 256;
  idx->names_size = 0;

  MEM_MALLOC(idx->names, idx->max_names_size, char);

  idx->max_types_size = 64;
  idx->types_size = 0;

  MEM_MALLOC(idx->types, idx->max_names_size, char);

  idx->max_data_size = 256;
  idx->data_size = 0;

  MEM_MALLOC(idx->data, idx->max_data_size, unsigned char);

  /* Add structure */

  inp->index = idx;

  /* Read headers to build index index */

  while (_file_tell(inp) + (long long)(inp->header_size) <= end_offset) {
    _read_section_header(inp);
    _update_index_and_shift(inp);
  }
}

/*----------------------------------------------------------------------------
 * Destroy a cs_io_t structure's optional index structure.
 *
 * parameters:
 *   inp <-> pointer to cs_io_t structure
 *----------------------------------------------------------------------------*/

static void
_destroy_index(_cs_io_t *inp)
{
  _cs_io_sec_index_t *idx = inp->index;

  if (idx == NULL)
    return;

  MEM_FREE(idx->h_vals);
  MEM_FREE(idx->offset);
  MEM_FREE(idx->names);
  MEM_FREE(idx->types);
  MEM_FREE(idx->data);

  MEM_FREE(inp->index);
}

/*----------------------------------------------------------------------------
 * Ready cs_io_t structure to read a given indexed section
 *
 * parameters:
 *   inp        <-> pointer to input object
 *   section_id <-- section id
 *----------------------------------------------------------------------------*/

static void
_set_indexed_section(_cs_io_t  *inp,
                     int        section_id)
{
  const _cs_io_sec_index_t *index = inp->index;
  const long long *h_vals = index->h_vals + section_id*8;

  inp->n_vals = h_vals[0];
  inp->location_id = h_vals[1];
  inp->index_id = h_vals[2];
  inp->n_loc_vals = h_vals[3];
  inp->type_size = 0;
  inp->data = NULL;
  if (h_vals[6] != 0)
    inp->data = index->data + h_vals[6] - 1;
  inp->name = index->names + h_vals[4];
  inp->type_name = index->types + h_vals[5];
  inp->offset = index->offset[section_id];
  inp->type_size = _type_size_from_name(inp->type_name);
}

/*----------------------------------------------------------------------------
 * Extract data from a cs_io_t structure whose index has been built
 *
 * parameters:
 *   inp         <-> pointer to input object
 *   location_id <-- if >= 0 location id filter
 *   sec_name    <-- if != NULL, section name filter
 *   f_fmt       <-- if != NULL, format for output of floating-point values
 *----------------------------------------------------------------------------*/

static void
_find_and_extract_section(_cs_io_t    *inp,
                          int          location_id,
                          const char  *sec_name,
                          const char  *f_fmt)
{
  size_t id;
  int extract_id = -1;
  _cs_io_sec_index_t *index = inp->index;

  /* Find matching sections */

  for (id = 0; id < index->size; id++) {

    int match = 1;
    const long long *h_vals = index->h_vals + id*8;
    const char *_name = index->names + h_vals[4];
    const int _location = h_vals[1];

    if (sec_name != NULL && strcmp(sec_name, _name))
      match = 0;
    if (location_id >= 0 && location_id != _location)
      match = 0;

    if (match == 1) {
      if (extract_id < 0)
        extract_id = id;
      else
        _error(__FILE__, __LINE__, 0,
               _("File \"%s\" contains multiple sections\n"
                 "named \"%s\" with location id %d\n\n"),
               inp->filename, _name, _location);
    }
  }

  /* If section to extract found, output it */

  if (extract_id > -1) {
    _set_indexed_section(inp, extract_id);
    _extract_section_values(inp, f_fmt);
  }
}

/*----------------------------------------------------------------------------
 * Copy data to compare buffer
 *
 * This allows comparing arrays of similar but not identical types.
 *
 * Type name should already have been checked when this function is called,
 * so no additional check is done here.
 *
 * parameters:
 *   dest       --> pointer to comparable (destination) data
 *   buffer     <-- pointer to data
 *   type_named <-- name of data type
 *   n_values   <-- number of values to echo
 *----------------------------------------------------------------------------*/

static void
_copy_to_cmp(void        *dest,
             const void  *buffer,
             const char   type_name[],
             size_t       n_values)
{
  size_t i;

  /* Check type name */

  if (type_name[0] == 'c')
    memcpy(dest, buffer, n_values);

#if (__STDC_VERSION__ >= 199901L)

  else if (type_name[0] == 'i' && type_name[1] == '4') {
    const int32_t *_buffer = buffer;
    long long *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

  else if (type_name[0] == 'i' && type_name[1] == '8') {
    const int64_t *_buffer = buffer;
    long long *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

  else if (type_name[0] == 'u' && type_name[1] == '4') {
    const uint32_t *_buffer = buffer;
    long long *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

  else if (type_name[0] == 'u' && type_name[1] == '8') {
    const uint64_t *_buffer = buffer;
    long long *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

#else /* (__STDC_VERSION__ < 199901L) */

  else if (type_name[0] == 'i' && type_name[1] == '4') {
    if (sizeof(int) == 4) {
      const int *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else if (sizeof(short) == 4) {
      const short *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else
      printf("    int32_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'i' && type_name[1] == '8') {
    if (sizeof(long) == 8) {
      const long *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else if (sizeof(long long) == 8) {
      const long long *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else
      printf("    int64_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'u' && type_name[1] == '4') {
    if (sizeof(unsigned) == 4) {
      const unsigned *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else if (sizeof(unsigned short) == 4) {
      const unsigned short *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else
      printf("    uint32_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

  else if (type_name[0] == 'u' && type_name[1] == '8') {
    if (sizeof(unsigned long) == 8) {
      const unsigned long *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else if (sizeof(unsigned long long) == 8) {
      const unsigned long long *_buffer = buffer;
      long long *_dest = dest;
      for (i = 0; i < n_values; i++)
        _dest[i] = _buffer[i];
    }
    else
      printf("    uint64_t undefined"
             " (porting error, C99 compiler needed)\n");
  }

#endif /* (__STDC_VERSION__) */

  else if (type_name[0] == 'r' && type_name[1] == '4') {
    const float *_buffer = buffer;
    double *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

  else if (type_name[0] == 'r' && type_name[1] == '8') {
    const double *_buffer = buffer;
    double *_dest = dest;
    for (i = 0; i < n_values; i++)
      _dest[i] = _buffer[i];
  }

}

/*----------------------------------------------------------------------------
 * Utility function to print header info for sections with differences.
 *
 * parameters:
 *   inp1  <-- pointer to first input object
 *   inp2  <-- pointer to second input object
 *----------------------------------------------------------------------------*/

static void
_echo_diff_headers(const _cs_io_t  *inp1,
                   const _cs_io_t  *inp2)
{
  if (strcmp(inp1->type_name, inp2->type_name))
    printf(_("  \"%-32s\"; Location: %2lu Size: %llu\n"
             "    Type: %-6s;  |  Type: %-6s; \n"),
           inp1->name,  inp1->location_id,
           (unsigned long long)inp1->n_vals,
           inp1->type_name, inp2->type_name);
  else
    printf(_("  \"%-32s\"; Location: %2lu; Type: %-6s; Size: %llu\n"),
           inp1->name,  inp1->location_id, inp1->type_name,
           (unsigned long long)inp1->n_vals);
}

/*----------------------------------------------------------------------------
 * Compare character data from 2 cs_io_t section buffers in different files
 *
 * parameters:
 *   inp1        <-> pointer to first input object
 *   inp2        <-> pointer to second input object
 *   cmp1        <-- buffer with characters from first file
 *   cmp1        <-- buffer with characters from second file
 *   block_start <-- buffer start id in total array
 *   block_size  <-- buffer size
 *   n_echo      <-- number of values to echo
 *   n_echo_cur  <-> number of values already echoed
 *
 * returns:
 *   number of different values
 *----------------------------------------------------------------------------*/

static size_t
_compare_chars(_cs_io_t    *inp1,
               _cs_io_t    *inp2,
               const char   cmp1[],
               const char   cmp2[],
               size_t       block_start,
               size_t       block_size,
               size_t       n_echo,
               long long   *n_echo_cur)
{
  size_t i;
  size_t n_diffs = 0;

  for (i = 0; i < block_size; i++) {
    if (cmp1[i] != cmp2[i])
      n_diffs++;
  }

  if (n_diffs > 0) {

    if (*n_echo_cur < 0) {
      _echo_diff_headers(inp1, inp2);
      *n_echo_cur = 0;
    }

    for (i = 0; i < block_size && (size_t)(*n_echo_cur) < n_echo; i++) {
      if (cmp1[i] != cmp2[i]) {
        unsigned long long j = block_start + i + 1;
        printf("    %12llu:  %c  | %c\n", j, cmp1[i], cmp2[i]);
        *n_echo_cur += 1;
      }
    }
  }

  return n_diffs;
}

/*----------------------------------------------------------------------------
 * Compare character data from 2 cs_io_t section buffers in different files
 *
 * parameters:
 *   inp1        <-> pointer to first input object
 *   inp2        <-> pointer to second input object
 *   cmp1        <-- buffer with integers from first file
 *   cmp1        <-- buffer with integers from second file
 *   block_start <-- buffer start id in total array
 *   block_size  <-- buffer size
 *   n_echo      <-- number of values to echo
 *   n_echo_cur  <-> number of values already echoed
 *
 * returns:
 *   number of different values
 *----------------------------------------------------------------------------*/

static size_t
_compare_ints(_cs_io_t         *inp1,
              _cs_io_t         *inp2,
              const long long   cmp1[],
              const long long   cmp2[],
              size_t            block_start,
              size_t            block_size,
              size_t            n_echo,
              long long        *n_echo_cur)
{
  size_t i;
  size_t n_diffs = 0;

  for (i = 0; i < block_size; i++) {
    if (cmp1[i] != cmp2[i])
      n_diffs++;
  }

  if (n_diffs > 0) {

    if (*n_echo_cur < 0) {
      _echo_diff_headers(inp1, inp2);
      *n_echo_cur = 0;
    }

    for (i = 0; i < block_size && (size_t)(*n_echo_cur) < n_echo; i++) {
      if (cmp1[i] != cmp2[i]) {
        unsigned long long j = block_start + i + 1;
        printf("    %12llu:  %lld  | %lld\n", j, cmp1[i], cmp2[i]);
        *n_echo_cur += 1;
      }
    }
  }

  return n_diffs;
}

/*----------------------------------------------------------------------------
 * Compare floating-point data from 2 cs_io_t section buffers in different
 * files
 *
 * parameters:
 *   inp1        <-> pointer to first input object
 *   inp2        <-> pointer to second input object
 *   cmp1        <-- buffer with integers from first file
 *   cmp1        <-- buffer with integers from second file
 *   f_fmt       <-- if != NULL, format for output of floating-point values
 *   f_threshold <-- threshold above which 2 floating-point values are
 *                   considered different.
 *   block_start <-- buffer start id in total array
 *   block_size  <-- buffer size
 *   n_echo      <-- number of values to echo
 *   n_echo_cur  <-> number of values already echoed
 *   f_stats     <-> max, total, max relative, total relative difference
 *
 * returns:
 *   number of different values
 *----------------------------------------------------------------------------*/

static size_t
_compare_floats(_cs_io_t         *inp1,
                _cs_io_t         *inp2,
                const double      cmp1[],
                const double      cmp2[],
                const char       *f_fmt,
                double            f_threshold,
                size_t            block_start,
                size_t            block_size,
                size_t            n_echo,
                long long        *n_echo_cur,
                double            f_stats[4])
{
  size_t i;
  size_t n_diffs = 0;

  for (i = 0; i < block_size; i++) {
    double delta = cmp1[i] - cmp2[i];
    if (delta < 0.0)
      delta = -delta;
    if (delta > f_threshold) {
      double delta_r = delta / CS_MAX(CS_ABS(cmp1[i]), CS_ABS(cmp2[i]));
      n_diffs++;
      if (delta > f_stats[0])
        f_stats[0] = delta;
      if (delta_r > f_stats[2])
        f_stats[2] = delta_r;
      f_stats[1] += delta;
      f_stats[3] += delta_r;
    }
  }

  if (n_diffs > 0) {

    char fmt[128] = "    %12llu:  %22.15e  | %22.15e\n";

    if (f_fmt != NULL) {
      strcpy(fmt, "    %12llu:  %");
      strcat(fmt, f_fmt);
      strcat(fmt, "  |  %");
      strcat(fmt, f_fmt);
      strcat(fmt, "\n");
    }

    if (*n_echo_cur < 0) {
      _echo_diff_headers(inp1, inp2);
      *n_echo_cur = 0;
    }

    for (i = 0; i < block_size && (size_t)(*n_echo_cur) < n_echo; i++) {
      double delta = cmp1[i] - cmp2[i];
      if (delta < 0.0)
        delta = -delta;
      if (delta > f_threshold) {
        unsigned long long j = block_start + i + 1;
        printf(fmt, j, cmp1[i], cmp2[i]);
        *n_echo_cur += 1;
      }
    }
  }

  return n_diffs;
}

/*----------------------------------------------------------------------------
 * Compare data from 2 cs_io_t sections in different files
 *
 * This function is expected to be called for sections with identical
 * names and locations.
 *
 * parameters:
 *   inp1        <-> pointer to first input object
 *   inp2        <-> pointer to second input object
 *   id1         <-- id of section in first input
 *   id2         <-- id of section in second input
 *   f_fmt       <-- if != NULL, format for output of floating-point values
 *   f_threshold <-- threshold above which 2 floating-point values are
 *                   considered different.
 *   n_echo      <-- maximum number of differences to output
 *
 * returns:
 *   1 if data differs, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_compare_sections(_cs_io_t    *inp1,
                  _cs_io_t    *inp2,
                  size_t       id1,
                  size_t       id2,
                  const char  *f_fmt,
                  double       f_threshold,
                  size_t       n_echo)
{
  char  compare_type1 = ' ', compare_type2 = ' ';
  int retval = 0;
  const char no_type[] = " ";
  const _cs_io_sec_index_t  *index1 = inp1->index;
  const _cs_io_sec_index_t  *index2 = inp2->index;
  const long long *h_vals1 = index1->h_vals + id1*8;
  const long long *h_vals2 = index2->h_vals + id2*8;
  const char *type1 = no_type;
  const char *type2 = no_type;
  const unsigned long long n_vals1 = h_vals1[0];
  const unsigned long long n_vals2 = h_vals2[0];
  const char *name = index1->names + h_vals1[4];
  const unsigned long location = h_vals1[1];

  /* If both sections have zero elements, they are identical
     (as their names have already been compared) */

  if (n_vals1 == 0 && n_vals2 == 0)
    return 0;

  /* Determine "equivalent" types for comparison; to reduce combinations,
     we will transform floats to doubles, and all integer
     types to 64-bit signed integers */

  if (n_vals1 > 0) {
    type1 = index1->types + h_vals1[5];
    compare_type1 = type1[0];
    if (type1[0] == 'u')
      compare_type1 = 'i';
  }
  if (n_vals2 > 0) {
    type2 = index2->types + h_vals2[5];
    compare_type2 = type2[0];
    if (type2[0] == 'u')
      compare_type2 = 'i';
  }

  /* If number differs or types are incompatible, sections differ */

  if (n_vals1 != n_vals2 || compare_type1 != compare_type2) {
    if (n_vals1 == n_vals2)
      printf(_("  \"%-32s\"; Location: %2lu; Size: %llu\n"
               "    Type: %-6s; |  Type: %-6s\n\n"),
             name,  location, n_vals1, type1, type2);
    else if (!strcmp(type1, type2))
      printf(_("  \"%-32s\"; Location: %2lu; Type: %-6s\n"
               "    Size: %llu  |  Size: %llu\n\n"),
             name,  location, type1, n_vals1, n_vals2);
    else
      printf(_("  \"%-32s\"; Location: %2lu\n"
               "    Type: %-6s; Size: %llu  |  Type: %-6s; Size: %llu\n\n"),
             name,  location, type1, n_vals1, type2, n_vals2);
    return 1;
  }

  /* If sections are comparable, their contents must be compared */

  else {

    unsigned long long n_diffs = 0;
    unsigned long long n_read = 0;
    long long n_echo_cur = -1;
    size_t block_size = n_vals1;
    size_t max_block_size = 2 << 16;
    void *buf1 = NULL, *buf2 = NULL;
    void *cmp1 = NULL, *cmp2 = NULL;
    double f_stats[4] = {0.0, 0.0, 0.0, 0.0};
    const size_t type_size1 = _type_size_from_name(type1);
    const size_t type_size2 = _type_size_from_name(type2);

    _set_indexed_section(inp1, id1);
    _set_indexed_section(inp2, id2);

    if (inp1->data == NULL && inp2->data == NULL && block_size > max_block_size)
      block_size = max_block_size;

    MEM_MALLOC(cmp1, block_size*8, unsigned char);
    MEM_MALLOC(cmp2, block_size*8, unsigned char);

    if (inp1->data == NULL) {
      MEM_MALLOC(buf1, block_size*type_size1, unsigned char);
      _file_seek(inp1, inp1->offset, SEEK_SET);
    }
    else
      buf1 = inp1->data;

    if (inp2->data == NULL) {
      MEM_MALLOC(buf2, block_size*type_size2, unsigned char);
      _file_seek(inp2, inp2->offset, SEEK_SET);
    }
    else
      buf2 = inp2->data;

    for (n_read = 0; n_read < n_vals1; n_read += block_size) {

      if (n_read + block_size > n_vals1)
        block_size = n_vals1 - n_read;

      if (inp1->data == NULL)
        _file_read(buf1, inp1->type_size, block_size, inp1);
      if (inp2->data == NULL)
        _file_read(buf2, inp2->type_size, block_size, inp2);

      _copy_to_cmp(cmp1, buf1, type1, block_size);
      _copy_to_cmp(cmp2, buf2, type2, block_size);

      if (compare_type1 == 'c')
        n_diffs += _compare_chars(inp1, inp2,
                                  cmp1, cmp2,
                                  n_read, block_size,
                                  n_echo, &n_echo_cur);
      else if (compare_type1 == 'i')
        n_diffs += _compare_ints(inp1, inp2,
                                 cmp1, cmp2,
                                 n_read, block_size,
                                 n_echo, &n_echo_cur);
      else if (compare_type1 == 'r')
        n_diffs += _compare_floats(inp1, inp2,
                                   cmp1, cmp2,
                                   f_fmt,
                                   f_threshold,
                                   n_read, block_size,
                                   n_echo, &n_echo_cur,
                                   f_stats);
    }

    if (n_diffs > 0) {
      if (type1[0] == 'r')
        printf(_("    Differences: %llu; Max: %g; Mean: %g; "
                 "Rel Max: %6.2e; Rel Mean: %6.2e\n\n"),
               n_diffs, f_stats[0], f_stats[1]/n_diffs,
               f_stats[2], f_stats[3]/n_diffs);
      else
        printf(_("    Differences: %llu\n\n"), n_diffs);
      retval = 1;
    }

    MEM_FREE(cmp1);
    MEM_FREE(cmp2);

    if (buf1 != inp1->data)
      MEM_FREE(buf1);
    if (buf2 != inp2->data)
      MEM_FREE(buf2);
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Compare data from 2 cs_io_t structures whose indexes have been built
 *
 * parameters:
 *   index <-- pointer to first index
 *   id    <-- id of section in index
 *----------------------------------------------------------------------------*/

static void
_echo_indexed_header(const _cs_io_sec_index_t  *index,
                     const size_t               id)
{
  const long long *h_vals = index->h_vals + id*8;
  const char *name = index->names + h_vals[4];
  const long long n_vals = h_vals[0];

  if (n_vals > 0) {
    const char *type = index->types + h_vals[5];
    const unsigned long location = h_vals[1];
    printf(_("  \"%-32s\"; Type: %-6s; Location: %2lu; Size: %llu\n"),
           name,  type, location, n_vals);
  }
  else
    printf(_("  \"%-32s\"\n"), name);
}

/*----------------------------------------------------------------------------
 * Compare data from 2 cs_io_t structures whose indexes have been built
 *
 * parameters:
 *   inp1        <-> pointer to first input object
 *   inp2        <-> pointer to second input object
 *   location_id <-- if >= 0 location id filter
 *   sec_name    <-- if != NULL, section name filter
 *   f_fmt       <-- if != NULL, format for output of floating-point values
 *   f_threshold <-- threshold above which 2 floating-point values are
 *                   considered different.
 *   n_echo      <-- maximum number of differences to output
 *
 * returns:
 *   0 if contents are identical, 1 if they differ
 *----------------------------------------------------------------------------*/

static int
_compare_files(_cs_io_t    *inp1,
               _cs_io_t    *inp2,
               int          location_id,
               const char  *sec_name,
               const char  *f_fmt,
               double       f_threshold,
               size_t       n_echo)
{
  size_t i, j;
  size_t n_diffs = 0;
  int has_unmatched1 = 0, has_unmatched2 = 0;
  int *compared1 = NULL, *compared2 = NULL;
  _cs_io_sec_index_t *index1 = inp1->index;
  _cs_io_sec_index_t *index2 = inp2->index;

  int retval = 0;

  /* Prepare marker on first file to flag sections with no match in
     second file, and marker on second file to flag sections compared
     during loop on first file. */

  MEM_MALLOC(compared1, index1->size, int);
  for (i = 0; i < index1->size; i++)
    compared1[i] = 0;

  MEM_MALLOC(compared2, index2->size, int);
  for (i = 0; i < index2->size; i++)
    compared2[i] = 0;

  /* Find matching sections */

  for (i = 0; i < index1->size; i++) {

    int match_filter = 1;
    const long long *h_vals1 = index1->h_vals + i*8;
    const char *_name1 = index1->names + h_vals1[4];
    const int _location1 = h_vals1[1];

    if (sec_name != NULL && strcmp(sec_name, _name1))
      match_filter = 0;
    if (location_id >= 0 && location_id != _location1)
      match_filter = 0;

    /* Search for matching section in second file */

    if (match_filter == 1) {

      for (j = 0; j < index2->size; j++) {

        const long long *h_vals2 = index2->h_vals + j*8;
        const char *_name2 = index2->names + h_vals2[4];
        const int _location2 = h_vals2[1];

        /* If matching section is found, compare it */

        if (!strcmp(_name1, _name2) && (_location1 == _location2)) {
          n_diffs += _compare_sections(inp1, inp2, i, j,
                                       f_fmt, f_threshold, n_echo);
          compared1[i] = 1;
          compared2[j] = 1;
        }
      }
    }
    else
      compared1[i] = 1;
  }

  /* List unmatched sections from first file */

  for (i = 0; i < index1->size; i++) {
    if (!compared1[i])
      has_unmatched1 = 1;
  }

  if (has_unmatched1) {

    printf(_("Sections only found in file \"%s\":\n\n"), inp1->filename);

    for (i = 0; i < index1->size; i++) {
      if (!compared1[i])
        _echo_indexed_header(index1, i);
    }

    printf("\n");
  }

  /* List unmatched sections from second file */

  for (i = 0; i < index2->size; i++) {

    const long long *h_vals2 = index2->h_vals + i*8;
    const char *_name2 = index2->names + h_vals2[4];
    const int _location2 = h_vals2[1];

    if (   (sec_name != NULL && strcmp(sec_name, _name2))
        || (location_id >= 0 && location_id != _location2))
      compared2[i] = 1;

    else if (!compared2[i])
      has_unmatched2 = 1;
  }

  if (has_unmatched2) {

    printf(_("Sections only found in file \"%s\":\n\n"), inp2->filename);

    for (i = 0; i < index2->size; i++) {
      if (!compared2[i])
        _echo_indexed_header(index2, i);
    }
  }

  MEM_FREE(compared1);
  MEM_FREE(compared2);

  if (n_diffs > 0 || has_unmatched1 || has_unmatched2)
    retval = 1;

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

int
main (int argc, char *argv[])
{
  int file_name_arg[2] = {0, 0};
  size_t echo = 0;
  int mode = 0, location_id = -1, sec_name_arg_id = 0, f_fmt_arg_id = 0;
  double f_threshold = 1.e-30;

  long long start_offset = 0, end_offset = 0;
  _cs_io_t inp;

  int retval = EXIT_SUCCESS;

  const char *sec_name = NULL, *f_fmt = NULL;

  if (getenv("LANG") != NULL)
    setlocale(LC_ALL,"");
  else
    setlocale(LC_ALL,"C");
  setlocale(LC_NUMERIC,"C");

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);
#endif

  /* Parse command line arguments */

  _read_args(argc,
             argv,
             &mode,
             &echo,
             &location_id,
             &sec_name_arg_id,
             &f_fmt_arg_id,
             &f_threshold,
             file_name_arg);

  if (sec_name_arg_id > 0)
    sec_name = argv[sec_name_arg_id];

  if (f_fmt_arg_id > 0)
    f_fmt = argv[f_fmt_arg_id];

  /* Initialize return arguments */

  inp = _open_input(argv[file_name_arg[0]], mode);

  /* Determine end of file;
     feof() may not work when using seeks,
     so we determine the size of the file first */

  start_offset = _file_tell(&inp);
  _file_seek(&inp, 0, SEEK_END);
  end_offset = _file_tell(&inp);
  _file_seek(&inp, start_offset, SEEK_SET);

  /* Standard dump mode */

  if (mode == 0) {

    /* Read file sections (or portions thereof) */

    while (  _file_tell(&inp) + (long long)(inp.header_size)
           <= end_offset)
      _read_section(&inp, echo, location_id, sec_name, f_fmt);

  }

  /* Extraction mode (build index to check for duplicates and find section) */

  else if (mode == 1) {
    _create_index(&inp, end_offset);
    _find_and_extract_section(&inp, location_id, sec_name, f_fmt);
    _destroy_index(&inp);
  }

  /* Diff mode */

  else if (mode == 2) {

    long long _start_offset = 0, _end_offset = 0;
    _cs_io_t inp2 = _open_input(argv[file_name_arg[1]], mode);

    /* Determine end of file as for first file */
    _start_offset = _file_tell(&inp);
    _file_seek(&inp2, 0, SEEK_END);
    _end_offset = _file_tell(&inp2);
    _file_seek(&inp2, _start_offset, SEEK_SET);

    /* Create indexes for both files */

    _create_index(&inp, end_offset);
    _create_index(&inp2, _end_offset);

    printf("\n");

    retval = _compare_files(&inp, &inp2, location_id, sec_name,
                            f_fmt, f_threshold, echo);

    _destroy_index(&inp2);
    _destroy_index(&inp);

    printf("\n");

    _close_input(&inp2, mode);
  }

  /* Clean-up */

  _close_input(&inp, mode);

  exit(retval);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
