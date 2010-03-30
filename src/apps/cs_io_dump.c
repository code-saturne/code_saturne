/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

/*============================================================================
 *  Dump of Kernel I/O file for Code_Saturne
 *============================================================================*/

/* Detect version of C used (C89 or C99) */

#if !defined(__STDC_VERSION__)
#  define __STDC_VERSION__ 1989
#endif

/* Include configuration file */

#include "cs_config.h"

/*
  Force LARGEFILE_SOURCE if large files enabled under 32-bit Linux or Blue Gene
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
#include <locale.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (__STDC_VERSION__ >= 199901L)
#include <stdint.h>
#endif

/*----------------------------------------------------------------------------
 * Internationalization macros
 *----------------------------------------------------------------------------*/

#if defined(ENABLE_NLS)

#include <libintl.h>
#define _(String) gettext(String)
#define gettext_noop(String) String
#define N_(String) gettext_noop(String)

#else

#define _(String) String
#define N_(String) String
#define textdomain(Domain)
#define bindtextdomain(Package, Directory)

#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
   *   5: index of embedded data in data array + 1 if data is
   *      embedded, 0 otherwise
   *   6: datatype id in file
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
       "  Usage: %s [-n <level>] <file_name>\n\n"
       "  Dump headers and optionnaly content of a Code_Saturne\n"
       "  Preprocessor, Partitioner, or restart file.\n\n"
       "  -n  <level>    number of first and last elements of each section\n"
       "                 to output (default: print headers only).\n\n"
       "  -h             this message.\n\n"),
     arg_0);

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
           _("Compilation configuration / porting error:\n"
             "Unable to determine a 64-bit unsigned int type.\n"
             "size_t is %d bits, unsigned long long %d bits"),
           sizeof(size_t)*8, sizeof(unsigned long long)*8);

#endif
}

/*----------------------------------------------------------------------------
 * Open input file, testing magic string for type.
 *
 * parameters:
 *   filename     <-- file name
 *   header_size  --> header default size
 *   header_align --> header alignment
 *   body_align   --> body alignment
 *
 * returns:
 *   File metadata structure
 *----------------------------------------------------------------------------*/

static _cs_io_t
_open_input(const char *filename)
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

  printf(_("\nOpening input file: \"%s\"\n\n"), filename) ;

  fflush(stdout);

  inp.f = fopen(inp.filename, "r");

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
  printf(_("\n  File type: %s\n"), header_buf);

  _file_read(header_buf, 8, 3, &inp);

  _convert_size((unsigned char*)header_buf, alignments, 3);

  inp.header_size = alignments[0];
  inp.header_align = alignments[1];
  inp.body_align = alignments[2];

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
 *   f <-> pointer to file object
 *----------------------------------------------------------------------------*/

static void
_close_input(_cs_io_t *inp)
{
  if (inp != NULL) {
    if (inp->f != NULL) {
      int retval = 0;
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
 * Read command line arguments.
 *
 * parameters:
 *   argc             <-- number of command line arguments
 *   argv             <-- array of command line arguments
 *   echo             --> echo (verbosity) level
 *   file_name_arg_id --> index of command line arguments defining file name
 *----------------------------------------------------------------------------*/

static void
_read_args(int               argc,
           char            **argv,
           size_t           *echo,
           int              *file_name_arg_id)
{
  int i = 1;

  /* Initialize return arguments */

  *echo = 0;
  *file_name_arg_id = 0;

  /* Parse and check command line */

  if (argc < 2)
    _usage(argv[0], EXIT_FAILURE);

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(argv[0], EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {

      i++;

      if (i >= argc)
        _usage(argv[0], EXIT_FAILURE);

      else {

#if (__STDC_VERSION__ >= 199901L)
        *echo = atoll(argv[i]);
#else
        *echo = atol(argv[i]);
#endif

      }

    }

    else {

      if (*file_name_arg_id == 0)
        *file_name_arg_id = i;
      else
        _usage(argv[0], EXIT_FAILURE);

    }

    i++;
  }

  if (*file_name_arg_id == 0)
    _usage(argv[0], EXIT_FAILURE);

  /* At this point, command line seems correct */

  printf(_("\n"
           "  .----------------------------.\n"
           "  |   Code_Saturne file dump   |\n"
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
 *----------------------------------------------------------------------------*/

static void
_echo_values(size_t       n_values,
             size_t       n_values_shift,
             const void  *buffer,
             const char   type_name[])
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
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %15.9e\n",
             (unsigned long)(i + n_values_shift),
             (double)(_buffer[i]));
  }

  else if (type_name[0] == 'r' && type_name[1] == '8') {
    const double *_buffer = buffer;
    for (i = 0; i < n_values; i++)
      printf("    %10lu : %22.15e\n",
             (unsigned long)(i + n_values_shift),
             _buffer[i]);
  }

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
  int type_name_error = 0;
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
  inp->index_id = header_vals[3],
  inp->n_loc_vals = header_vals[4];
  inp->type_size = 0;
  inp->data = NULL;
  inp->type_name = (char *)(inp->buffer + 48);
  inp->name = (char *)(inp->buffer + 56);

  if (header_vals[1] > 0 && inp->type_name[7] == 'e')
    inp->data = inp->buffer + 56 + header_vals[5];

  inp->type_size = 0;

  if (inp->n_vals > 0) {

    /* Check type name and compute size of data */

    if (inp->type_name[0] == 'c') {
      if (inp->type_name[1] != ' ')
        type_name_error = 1;
      else
        inp->type_size = 1;
    }
    else if (   inp->type_name[0] == 'i'
             || inp->type_name[0] == 'u'
             || inp->type_name[0] == 'r') {

      if (inp->type_name[1] == '4')
        inp->type_size = 4;
      else if (inp->type_name[1] == '8')
        inp->type_size = 8;
      else
        type_name_error = 1;

    }
    else
      type_name_error = 1;

    if (type_name_error)
      _error(__FILE__, __LINE__, 0,
             _("Type \"%s\" is not known\n"
               "Known types: \"c \", \"i4\", \"i8\", \"u4\", \"u8\", "
               "\"r4\", \"r8\"."), inp->type_name);

    else if (inp->data == NULL)
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
 *   inp  <-> pointer to input object
 *   echo <-- number of values to print
 *----------------------------------------------------------------------------*/

static void
_read_section_values(_cs_io_t  *inp,
                     size_t     echo)
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

    _echo_values(n_print, 1, data, inp->type_name);

    if (n_skip > 0)
      printf("    ..........   ..........\n");

  }

  /* Mode to tail of data and read if necessary */

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
      _echo_values(n_print, inp->n_vals - n_print + 1, data, inp->type_name);

  }

  if (buffer != NULL)
    MEM_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Read section.
 *
 * parameters:
 *   inp  <-> pointer to input object
 *   echo <-- number of values to print
 *----------------------------------------------------------------------------*/

static void
_read_section(_cs_io_t  *inp,
              int        echo)
{
  assert(inp != NULL);
  assert(inp->f != NULL);

  /* Read section header and print basic information */

  _read_section_header(inp);

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

  if (inp->n_vals > 0)
    _read_section_values(inp, echo);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

int
main (int argc, char *argv[])
{
  int file_name_arg = 0;
  size_t echo = 0;
  long long start_offset = 0, end_offset = 0;
  _cs_io_t inp;

  if (getenv("LANG") != NULL)
     setlocale(LC_ALL,"");
  else
     setlocale(LC_ALL,"C");
  setlocale(LC_NUMERIC,"C");

#if defined(ENABLE_NLS)
  bindtextdomain (PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);
#endif

  /* Parse command line arguments */

  _read_args(argc, argv, &echo, &file_name_arg);

  inp = _open_input(argv[file_name_arg]);

  /* Determine end of file;
     feof() may not work when using seeks,
     so we determine the size of the file first */

  start_offset = _file_tell(&inp);
  _file_seek(&inp, 0, SEEK_END);
  end_offset = _file_tell(&inp);
  _file_seek(&inp, start_offset, SEEK_SET);

  /* Read file sections (or portions thereof) */

  while (  _file_tell(&inp) + (long long)(inp.header_size)
         <= end_offset)
    _read_section(&inp, echo);

  /* Clean-up */

  _close_input(&inp);

  exit(EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
