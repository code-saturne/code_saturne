/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 *  Low level file I/O utility functions for Preprocessor and restart files
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#undef HAVE_STDINT_H
#if defined(__STDC_VERSION__)
#  if (__STDC_VERSION__ >= 199901L)
#    define HAVE_STDINT_H
#    include <stdint.h>
#  endif
#endif

/*----------------------------------------------------------------------------
 * BFT and FVM library headers
 *----------------------------------------------------------------------------*/

#include <bft_file.h>
#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include <fvm_file.h>

#if defined(HAVE_MPI)
#include <fvm_parall.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_io.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local types and structures
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

  fvm_file_off_t *h_vals;            /* Base values associated
                                        with each header */

  fvm_file_off_t *offset;            /* Position of associated data
                                        in file (-1 if embedded) */

  size_t          max_names_size;    /* Maximum size of names array */
  size_t          names_size;        /* Current size of names array */
  char           *names;             /* Array containing section names */

  size_t          max_data_size;     /* Maximum size of embedded data array */
  size_t          data_size;         /* Current size of data array */
  unsigned char  *data;              /* Array containing embedded data */

  /* An index maintains its own list of files, in case multiple
     files are used to contain the data */

  size_t          n_files;           /* Number of associated files */
  fvm_file_t    **f;                 /* Pointers to associated files */

} cs_io_sec_index_t;

/* Main kernel IO state structure */
/*--------------------------------*/

struct _cs_io_t {

  /* File information */

  fvm_file_t         *f;              /* Pointer to associated file */

  char                contents[64];   /* String describing file contents */

  cs_io_mode_t        mode;           /* File access mode */

  size_t              header_size;    /* Header default size */
  size_t              header_align;   /* Header alignment */
  size_t              body_align;     /* Body alignment */

  cs_io_sec_index_t  *index;          /* Optional section index (on read) */

  /* Current section buffer state */

  size_t              buffer_size;    /* Current size of header buffer */
  unsigned char      *buffer;         /* Header buffer */

  fvm_file_off_t      n_vals;         /* Number of values in section header */
  size_t              location_id;    /* Id of location, or 0 */
  size_t              index_id;       /* Id of index, or 0 */
  size_t              n_loc_vals;     /* Number of values per location */
  size_t              type_size;      /* Size of current type */
  char               *sec_name;       /* Pointer to name in section header */
  char               *type_name;      /* Pointer to type in section header */
  void               *data;           /* Pointer to data in section header
                                         (if embedded; NULL otherwise) */

  /* Other flags */

  long                echo;           /* Data echo level (verbosity) */
};

/*============================================================================
 * Constants and Macros
 *============================================================================*/

#define CS_IO_MPI_TAG     'C'+'S'+'_'+'I'+'O'

/*============================================================================
 * Static global variables
 *============================================================================*/

static char  _type_name_none[] = "  ";
static char  _type_name_char[] = "c ";  /* Character string */
static char  _type_name_i4[] =   "i4";  /* Signed 32 bit integer */
static char  _type_name_i8[] =   "i8";  /* Signed 64 bit integer */
static char  _type_name_u4[] =   "u4";  /* Unsigned 32 bit integer */
static char  _type_name_u8[] =   "u8";  /* Unsigned 64 bit integer */
static char  _type_name_r4[] =   "r4";  /* Single precision real */
static char  _type_name_r8[] =   "r8";  /* Double precsision real */

/* Default hints for files using this API (for MPI-IO) */
#if defined(FVM_HAVE_MPI_IO)
int  cs_glob_io_hints = FVM_FILE_EXPLICIT_OFFSETS;
#else
int  cs_glob_io_hints = 0;
#endif

/* Global pointer on preprocessor data file handle */
cs_io_t  *cs_glob_pp_io = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Default conversion rule from type in file to type in memory.
 *
 * parameters:
 *   type_read <-- type in file
 *
 * returns:
 *   default corresponding type in memory (may need conversion)
 *----------------------------------------------------------------------------*/

static fvm_datatype_t
_type_read_to_elt_type(fvm_datatype_t type_read)
{
  fvm_datatype_t elt_type = FVM_DATATYPE_NULL;

  if (type_read == FVM_INT32 || type_read == FVM_INT64) {
    assert(sizeof(fvm_lnum_t) == 4 || sizeof(fvm_lnum_t) == 8);
    if (sizeof(fvm_lnum_t) == 4)
      elt_type = FVM_INT32;
    else
      elt_type = FVM_INT64;
  }

  else if (type_read == FVM_UINT32 || type_read == FVM_UINT64) {
    assert(sizeof(fvm_gnum_t) == 4 || sizeof(fvm_gnum_t) == 8);
    if (sizeof(fvm_gnum_t) == 4)
      elt_type = FVM_UINT32;
    else
      elt_type = FVM_UINT64;
  }

  else if (type_read == FVM_FLOAT || type_read == FVM_DOUBLE) {
    if (sizeof(cs_real_t) == 4)
      elt_type = FVM_FLOAT;
    else
      elt_type = FVM_DOUBLE;
  }

  else if (type_read == FVM_CHAR)
    elt_type = FVM_CHAR;

  return elt_type;
}

/*----------------------------------------------------------------------------
 * Convert a buffer of type uint64_t to fvm_file_off_t
 *
 * parameters:
 *   buf <-- buffer
 *   val --> array to which values are converted
 *   n   <-- number of values to convert
 *----------------------------------------------------------------------------*/

static void
_convert_to_offset(const unsigned char  buf[],
                   fvm_file_off_t       val[],
                   size_t               n)
{
  size_t i;

#if defined(HAVE_STDINT_H)

  for (i = 0; i < n; i++)
    val[i] = ((const uint64_t *)buf)[i];

#else

  if (sizeof(size_t) == 8) {
    for (i = 0; i < n; i++)
      val[i] = ((const size_t *)buf)[i];
  }
  else if (sizeof(unsigned long long) == 8) {
    for (i = 0; i < n; i++)
      val[i] = ((const unsigned long long *)buf)[i];
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Compilation configuration / porting error:\n"
                "Unable to determine a 64-bit unsigned int type.\n"
                "size_t is %d bits, unsigned long long %d bits"),
              sizeof(size_t)*8, sizeof(unsigned long long)*8);

#endif
}

/*----------------------------------------------------------------------------
 * Convert a buffer of type fvm_file_off_t to uint64_t
 *
 * parameters:
 *   buf --> buffer
 *   val <-- array from which values are converted
 *   n   <-- number of values to convert
 *----------------------------------------------------------------------------*/

static void
_convert_from_offset(unsigned char         buf[],
                     const fvm_file_off_t  val[],
                     size_t                n)
{
  size_t i;

#if defined(HAVE_STDINT_H)

  for (i = 0; i < n; i++)
    ((uint64_t *)buf)[i]=  val[i];

#else

  if (sizeof(size_t) == 8) {
    for (i = 0; i < n; i++)
      ((size_t *)buf)[i] = val[i];
  }
  else if (sizeof(unsigned long long) == 8) {
    for (i = 0; i < n; i++)
      ((unsigned long long *)buf)[i] = val[i];
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Compilation configuration / porting error:\n"
                "Unable to determine a 64-bit unsigned int type.\n"
                "size_t is %d bits, unsigned long long %d bits"),
              sizeof(size_t)*8, sizeof(unsigned long long)*8);

#endif
}

/*----------------------------------------------------------------------------
 * Return an empty kernel IO file structure.
 *
 * parameters:
 *   mode     --> read or write
 *   echo     --> echo on main output (< 0 if none, header if 0,
 *                n first and last elements if n > 0)
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

static cs_io_t *
_cs_io_create(cs_io_mode_t   mode,
              size_t         echo)
{
  cs_io_t  *cs_io = NULL;

  BFT_MALLOC(cs_io, 1, cs_io_t);

  /* Set structure fields */

  cs_io->mode = mode;

  cs_io->f  = NULL;

  memset(cs_io->contents, 0, 64);

  cs_io->header_size = 0;
  cs_io->header_align = 0;
  cs_io->body_align = 0;

  cs_io->index = NULL;

  /* Current section buffer state */

  cs_io->buffer_size = 0;
  cs_io->buffer = NULL;

  cs_io->n_vals = 0;
  cs_io->type_size = 0;
  cs_io->sec_name = NULL;
  cs_io->type_name = NULL;
  cs_io->data = NULL;

  /* Verbosity */

  cs_io->echo = echo;

  return cs_io;
}

/*----------------------------------------------------------------------------
 * Add an empty index structure to a cs_io_t structure.
 *
 * parameters:
 *   inp <-> pointer to cs_io_t structure
 *----------------------------------------------------------------------------*/

static void
_create_index(cs_io_t *inp)
{
  cs_io_sec_index_t  *idx = NULL;

  BFT_MALLOC(idx, 1, cs_io_sec_index_t);

  /* Set structure fields */

  idx->size = 0;
  idx->max_size = 32;

  BFT_MALLOC(idx->h_vals, idx->max_size*8, fvm_file_off_t);
  BFT_MALLOC(idx->offset, idx->max_size, fvm_file_off_t);

  idx->max_names_size = 256;
  idx->names_size = 0;

  BFT_MALLOC(idx->names, idx->max_names_size, char);

  idx->max_data_size = 256;
  idx->data_size = 0;

  BFT_MALLOC(idx->data, idx->max_data_size, unsigned char);

  idx->n_files = 0;
  idx->f = NULL;

  /* Add structure */

  inp->index = idx;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_io_t structure's optional index structure.
 *
 * parameters:
 *   inp <-> pointer to cs_io_t structure
 *----------------------------------------------------------------------------*/

static void
_destroy_index(cs_io_t *inp)
{
  size_t i;
  cs_io_sec_index_t *idx = inp->index;

  if (idx == NULL)
    return;

  BFT_FREE(idx->h_vals);
  BFT_FREE(idx->offset);
  BFT_FREE(idx->names);
  BFT_FREE(idx->data);

  for (i = 0; i < idx->n_files; i++) {
    if (idx->f[i] == inp->f)
      idx->f[i] = NULL;
    else if (idx->f[i] != NULL)
      idx->f[i] = fvm_file_free(idx->f[i]);
  }
  BFT_FREE(idx->f);

  BFT_FREE(inp->index);
}

/*----------------------------------------------------------------------------
 * Update an index structure with info from the last header read
 *
 * Also sets the file position for the next read
 *
 * parameters:
 *   inp    <-> input kernel IO structure
 *   header <-- header structure
 *----------------------------------------------------------------------------*/

static void
_update_index_and_shift(cs_io_t             *inp,
                        cs_io_sec_header_t  *header)
{
  size_t id = 0;
  size_t new_names_size = 0;
  size_t new_data_size = 0;

  cs_io_sec_index_t  *idx = inp->index;

  if (idx == NULL)
    return;

  /* Reallocate if necessary */

  if (idx->size + 1 == idx->max_size) {
    if (idx->max_size == 0)
      idx->max_size = 32;
    else
      idx->max_size *= 2;
    BFT_REALLOC(idx->h_vals, idx->max_size*8, fvm_file_off_t);
    BFT_REALLOC(idx->offset, idx->max_size, fvm_file_off_t);
  };

  new_names_size = idx->names_size + strlen(inp->sec_name) + 1;

  if (inp->data != NULL)
    new_data_size
      = idx->data_size + (  inp->n_vals
                          * fvm_datatype_size[header->type_read]);

  if (new_names_size > idx->max_names_size) {
    if (idx->max_names_size == 0)
      idx->max_names_size = 128;
    while (new_names_size > idx->max_names_size)
      idx->max_names_size *= 2;
    BFT_REALLOC(idx->names, idx->max_names_size, char);
  }

  if (new_data_size > idx->max_data_size) {
    if (idx->max_data_size == 0)
      idx->max_data_size = 128;
    while (new_data_size > idx->max_data_size)
      idx->max_data_size *= 2;
    BFT_REALLOC(idx->data, idx->max_data_size, unsigned char);
  }

  /* Set values */

  id = idx->size;

  idx->h_vals[id*8]     = inp->n_vals;
  idx->h_vals[id*8 + 1] = inp->location_id;
  idx->h_vals[id*8 + 2] = inp->index_id;
  idx->h_vals[id*8 + 3] = inp->n_loc_vals;
  idx->h_vals[id*8 + 4] = idx->names_size;
  idx->h_vals[id*8 + 5] = 0;
  idx->h_vals[id*8 + 6] = header->type_read;
  idx->h_vals[id*8 + 7] = idx->n_files - 1;

  strcpy(idx->names + idx->names_size, inp->sec_name);
  idx->names[new_names_size - 1] = '\0';
  idx->names_size = new_names_size;

  if (inp->data == NULL) {
    fvm_file_off_t offset = fvm_file_tell(inp->f);
    fvm_file_off_t data_shift = inp->n_vals * inp->type_size;
    if (inp->body_align > 0) {
      size_t ba = inp->body_align;
      idx->offset[id] = offset + (ba - (offset % ba)) % ba;
    }
    else
      idx->offset[id] = offset;
    fvm_file_seek(inp->f, idx->offset[id] + data_shift, FVM_FILE_SEEK_SET);
  }
  else {
    idx->h_vals[id*8 + 5] = idx->data_size + 1;
    memcpy(idx->data + idx->data_size,
           inp->data,
           new_data_size - idx->data_size);
    idx->data_size = new_data_size;
    idx->offset[id] = -1;
  }

  idx->size += 1;
}

/*----------------------------------------------------------------------------
 * Open the interface file descriptor and initialize the file by writing
 * or reading a "magic string" used to check the file content type.
 *
 * parameters:
 *   cs_io        <-> kernel IO structure
 *   name         <-- file name
 *   magic_string <-- magic string associated with file content type
 *   hints        <-- file handling method options (0 for default)
 *   comm         <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
static void
_file_open(cs_io_t     *cs_io,
           const char  *name,
           const char  *magic_string,
           int          hints,
           MPI_Comm     comm)
#else
static void
_file_open(cs_io_t     *cs_io,
           const char  *name,
           const char  *magic_string,
           int          hints)
#endif
{
  fvm_file_mode_t f_mode;
  char header_data[128 + 24];
  fvm_file_off_t header_vals[3];

  char  base_header[] = "Code_Saturne I/O, BE, R0";

  /* Prepare file open */

  switch(cs_io->mode) {

  case CS_IO_MODE_READ:
    f_mode = FVM_FILE_MODE_READ;
    break;

  case CS_IO_MODE_WRITE:
    f_mode = FVM_FILE_MODE_WRITE;
    break;

  default:
    assert(   cs_io->mode == CS_IO_MODE_READ
           || cs_io->mode == CS_IO_MODE_WRITE);
    return;
  }


  /* Create interface file descriptor */

#if defined(HAVE_MPI)
  cs_io->f = fvm_file_open(name, f_mode, hints, comm);
#else
  cs_io->f = fvm_file_open(name, f_mode, hints);
#endif

  fvm_file_set_big_endian(cs_io->f);

  /* Write or read a magic string */
  /*------------------------------*/

  if (cs_io->mode == CS_IO_MODE_READ) {

    fvm_file_read_global(cs_io->f, header_data, 1, 128 + 24);

    header_data[63] = '\0';
    header_data[127] = '\0';

    /* If the format does not correspond, we have an error */

    if (strncmp(header_data, base_header, 64) != 0) {

      bft_error(__FILE__, __LINE__, 0,
                _("Error reading file: \"%s\".\n"
                  "File format is not the correct version.\n"
                  "The first 64 bytes expected contain:\n"
                  "\"%s\"\n"
                  "The first 64 bytes read contain:\n"
                  "\"%s\"\n"),
                fvm_file_get_name(cs_io->f), base_header, header_data);

    }

    /* Copy magic string */

    strncpy(cs_io->contents, header_data + 64, 64);
    cs_io->contents[63] = '\0';

    /* If the magic string does not correspond, we have an error */

    if (magic_string != NULL) {
      if (strncmp(cs_io->contents, magic_string, 64) != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error reading file: \"%s\".\n"
                    "The file contents are not of the expected type.\n"
                    "\"%s\" was expected,\n"
                    "\"%s\" was read."),
                  fvm_file_get_name(cs_io->f), magic_string, cs_io->contents);
    }

    /* Now decode the sizes */

    if (fvm_file_get_swap_endian(cs_io->f) == 1)
      bft_file_swap_endian(header_data + 128, header_data + 128, 8, 3);

    _convert_to_offset((unsigned char *)(header_data + 128), header_vals, 3);

    cs_io->header_size = header_vals[0];
    cs_io->header_align = header_vals[1];
    cs_io->body_align = header_vals[2];

  }
  else if (cs_io->mode == CS_IO_MODE_WRITE) {

    size_t n_written = 0;

    memset(header_data, 0, sizeof(header_data));
    strcpy(header_data, base_header);
    strncpy(header_data + 64, magic_string, 64);
    header_data[127] = '\0';

    /* Set default header size and alignments */

    cs_io->header_size = 128;
    cs_io->header_align = 64;
    cs_io->body_align = 64;

    header_vals[0] = cs_io->header_size;
    header_vals[1] = cs_io->header_align;
    header_vals[2] = cs_io->body_align;

    _convert_from_offset((unsigned char *)(header_data + 128), header_vals, 3);

    if (fvm_file_get_swap_endian(cs_io->f) == 1)
      bft_file_swap_endian(header_data + 128, header_data + 128, 8, 3);

    n_written = fvm_file_write_global(cs_io->f, header_data, 1, 128 + 24);

    if (n_written < 128 + 24)
      bft_error(__FILE__, __LINE__, 0,
                _("Error writing the header of file: \"%s\".\n"),
                fvm_file_get_name(cs_io->f));
  }

  cs_io->buffer_size = cs_io->header_size;
  BFT_MALLOC(cs_io->buffer, cs_io->buffer_size, unsigned char);
}

/*----------------------------------------------------------------------------
 * Add fictitious sections for mesh sizes for a legacy restart file.
 *
 * parameters:
 *   inp     <-> kernel IO input structure
 *   sizes   <-- number of cells, interior face, bountery face, and vertices
 *----------------------------------------------------------------------------*/

static void
_file_legacy_add_sizes(cs_io_t     *inp,
                       fvm_lnum_t   sizes[4])
{
  int i;
  cs_io_sec_header_t h;

  char _sec_name[32];
  const char *sec_name[] = {"cells",
                            "interior_faces",
                            "boundary_faces",
                            "vertices"};

  assert(inp->mode == CS_IO_MODE_READ);
  assert(inp->index != NULL);

  /* Common initializations */

  h.n_location_vals = 1;
  h.type_read = FVM_INT32;
  h.elt_type = FVM_INT32;

  /* Add 4 sizes and associated locations as sections with embedded data */

  for (i = 0; i < 4; i++) {

    size_t embedded_val = sizes[i];

    strcpy(_sec_name, sec_name[i]); /* Copy to avoid const warnings */
    inp->sec_name = _sec_name;
    inp->n_vals = 1;
    inp->location_id = i+1;
    inp->index_id = 0;
    inp->n_loc_vals = 1;
    inp->data = &embedded_val;

    h.sec_name = inp->sec_name;
    h.n_vals = inp->n_vals;
    h.location_id = inp->location_id;

    assert(sizeof(size_t) == 4 || sizeof(size_t) == 8);

    if (sizeof(size_t) == 4)
      h.type_read = FVM_UINT32;
    else
      h.type_read = FVM_UINT64;
    h.elt_type = _type_read_to_elt_type(h.type_read);

    /* Todo add to index */

    _update_index_and_shift(inp, &h);

  }
}

/*----------------------------------------------------------------------------
 * Check if a file is a Code_Saturne version 1.1 - 1.3 restart file,
 * and open it if it is.
 *
 * parameters:
 *   inp   <-> kernel IO input structure
 *   name  <-- file name
 *   sizes <-- location sizes
 *   comm  <-- associated MPI communicator
 *
 * returns:
 *   1 if file is a legacy restart file, 0 otherwise.
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
static int
_file_legacy_restart_open(cs_io_t     *inp,
                          const char  *name,
                          fvm_lnum_t   sizes[4],
                          MPI_Comm     comm)
#else
static int
_file_legacy_restart_open(cs_io_t    *inp,
                          const char  *name,
                          fvm_lnum_t   sizes[4])
#endif
{
  char expected_header[] = "Code_Saturne_1.1_bin_reprise\n";
  char header_read[32] = "";
  fvm_lnum_t expected_len = strlen(expected_header);

  fvm_lnum_t n_read = 0;
  int retval = 0;

  assert(inp->mode == CS_IO_MODE_READ);

  /* Create interface file descriptor */

#if defined(HAVE_MPI)
  inp->f = fvm_file_open(name, FVM_FILE_MODE_READ, FVM_FILE_NO_MPI_IO, comm);
#else
  inp->f = fvm_file_open(name, FVM_FILE_MODE_READ, 0);
#endif

  fvm_file_set_big_endian(inp->f);

  /* Read first characters and compare */
  /*-----------------------------------*/

  n_read = fvm_file_read_global(inp->f, header_read, 1, expected_len);

  if (n_read == expected_len) {
    header_read[expected_len] = '\0';
    if (strcmp(header_read, expected_header) == 0)
      retval = 1;
  }

  /* Close file and return if it is not what we are looking for */

  if (retval == 0) {
    inp->f = fvm_file_free(inp->f);
    return retval;
  }

  /* From here on, we know the file is in the legacy restart format */
  /*----------------------------------------------------------------*/

  if (inp->buffer_size == 0) {
    inp->buffer_size = 128;
    BFT_MALLOC(inp->buffer, inp->buffer_size, unsigned char);
    memset(inp->buffer, 0, inp->buffer_size);
  }

  /* Read location sizes */

  n_read = fvm_file_read_global(inp->f, sizes, sizeof(fvm_lnum_t), 4);
  if (n_read < 4) {
    bft_error(__FILE__, __LINE__, 0,
              _("Restart file \"%s\"\n"
                "in format 1.1 is not conforming."),
              fvm_file_get_name(inp->f));

    /* following code will not be reached as long as errors are fatal */
    inp->f = fvm_file_free(inp->f);
    retval = 0;
  }

  /* Update index file section */
  if (inp->index != NULL) {
    BFT_REALLOC(inp->index->f, inp->index->n_files + 1, fvm_file_t *);
    inp->index->f[inp->index->n_files] = inp->f;
    inp->index->n_files += 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a file is a legacy restart file, and open it if it is.
 *
 * parameters:
 *   inp     <-> kernel IO input structure
 *   name    <-- file name
 *   comm    <-- associated MPI communicator
 *
 * returns:
 *   1 if file is a legacy restart file, 0 otherwise.
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
static int
_file_legacy_restart_index(cs_io_t     *inp,
                           const char  *name,
                           MPI_Comm     comm)
#else
static int
_file_legacy_restart_index(cs_io_t     *inp,
                           const char  *name)
#endif
{
  fvm_lnum_t sizes[4] = {0, 0, 0, 0};

  fvm_lnum_t n_read = 0;
  int end_of_file = 0;
  int retval = 0;

  const char incorrect_next_file_msg[]
    = N_("Restart file \"%s\" does not correspond\n"
         "to part %d of the original restart file.");

#if defined(FVM_HAVE_MPI)
  retval = _file_legacy_restart_open(inp, name, sizes, comm);
#else
  retval = _file_legacy_restart_open(inp, name, sizes);
#endif

  if (retval == 0)
    return retval;

  /* From here on, we know the file is in the legacy restart format */

  _file_legacy_add_sizes(inp, sizes);

  /* Now analyze the file */

  inp->header_size = 0;
  inp->header_align = 0;
  inp->body_align = 0;
  inp->data = NULL;

  while (end_of_file == 0) {

    fvm_lnum_t buf[4];
    char *sec_name = NULL;

    /* Read section */

    n_read = fvm_file_read_global(inp->f, buf, sizeof(fvm_lnum_t), 4);

    if (n_read < 4) {
      end_of_file = 1;
      break;
    }

    if (buf[0] + 56 >= (fvm_lnum_t)(inp->buffer_size)) {
      while (buf[0] + 56 >= (fvm_lnum_t)(inp->buffer_size))
        inp->buffer_size *= 2;
      BFT_REALLOC(inp->buffer, inp->buffer_size, unsigned char);
    }
    sec_name = (char *)(inp->buffer + 56);

    n_read = fvm_file_read_global(inp->f, sec_name, 1, buf[0]);
    sec_name[n_read] = '\0';

    if (n_read < buf[0]) {
      end_of_file = 1;
      break;
    }

    /* Now handle section, starting with special records */
    /*---------------------------------------------------*/

    /* If contents continue on another file, switch files */

    if (strcmp(sec_name, "reprise : fic suivant") == 0) {

      size_t ii;
      fvm_lnum_t cmp_sizes[4] = {0, 0, 0, 0};
      size_t _name_len = strlen(name) + strlen("_pxx");
      char *_name = NULL;

      /* Build next restart file name */
      BFT_MALLOC(_name, _name_len + 1, char);
      sprintf(_name, "%s_p%02d", name, (int)(inp->index->n_files + 1));

      /* Open new file */
      end_of_file = 0;
      inp->index->f[inp->index->n_files - 1] = inp->f;
#if defined(FVM_HAVE_MPI)
      retval = _file_legacy_restart_open(inp, _name, cmp_sizes, comm);
#else
      retval = _file_legacy_restart_open(inp, _name, cmp_sizes);
#endif
      if (retval == 0)
        end_of_file = 1;

      for (ii = 0; ii < 4; ii++) {
        if (sizes[ii] != cmp_sizes[ii]) {
          bft_error(__FILE__, __LINE__, 0, _(incorrect_next_file_msg),
                    fvm_file_get_name(inp->f), (int)(inp->index->n_files));
          end_of_file = 1;
        }
      }
      BFT_FREE(_name);
    }

    /* If end of file is indicated */

    else if (strcmp(sec_name, "reprise : fin") == 0)
      end_of_file = 1;

    /* If the beginning of a new file is indicated */

    else if (strcmp(sec_name, "reprise : partie num") == 0) {
      if (buf[0] != (fvm_lnum_t)(inp->index->n_files)) {
        bft_error(__FILE__, __LINE__, 0, _(incorrect_next_file_msg),
                  fvm_file_get_name(inp->f), (int)(inp->index->n_files));
        end_of_file = 1;
      }
      continue;
    }

    /* Standard record */

    else {

      /* Prepare addition to index */

      cs_io_sec_header_t h;

      inp->sec_name = sec_name;
      inp->n_vals = buf[2];
      if (buf[1] > 0)
        inp->n_vals *= buf[2] * sizes[buf[1] - 1];
      inp->location_id = buf[1];
      inp->index_id = 0;
      inp->n_loc_vals = buf[2];

      h.sec_name = inp->sec_name;
      h.n_vals = inp->n_vals;
      h.location_id = inp->location_id;
      h.n_location_vals = inp->n_loc_vals;
      h.type_read = FVM_DATATYPE_NULL;
      if (buf[3] == 0)
        h.type_read = FVM_CHAR;
      else if (buf[3] == 1)
        h.type_read = FVM_INT32;
      else if (buf[3] == 2)
        h.type_read = FVM_DOUBLE;
      h.elt_type = _type_read_to_elt_type(h.type_read);

      inp->type_size = fvm_datatype_size[h.type_read];

      /* Todo add to index */

      _update_index_and_shift(inp, &h);

    }
  }

  return 1;
}

/*----------------------------------------------------------------------------
 * Re-open the interface file descriptor when building an index.
 *
 * parameters:
 *   inp          <-> kernel IO structure
 *   hints        <-- file handling method options (0 for default)
 *   comm         <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
static void
_file_reopen_read(cs_io_t   *inp,
                  int        hints,
                  MPI_Comm   comm)
#else
static void
_file_reopen_read(cs_io_t   *inp,
                  int        hints)
#endif /* FVM_HAVE_MPI */
{
  size_t i;
  char _tmpname[128];
  char *tmpname = _tmpname;

  assert(inp->index != NULL);

  for (i = 0; i < inp->index->n_files; i++) {

    const char *filename = fvm_file_get_name(inp->index->f[i]);

    if (strlen(filename) >= 128)
      BFT_MALLOC(tmpname, strlen(filename) + 1, char);
    strcpy(tmpname, filename);

    inp->index->f[i] = fvm_file_free(inp->index->f[i]);

#if defined(HAVE_MPI)
    inp->index->f[i] = fvm_file_open(tmpname, FVM_FILE_MODE_READ, hints, comm);
#else
    inp->index->f[i] = fvm_file_open(tmpname, FVM_FILE_MODE_READ, hints);
#endif

    fvm_file_set_big_endian(inp->index->f[i]);

    if (tmpname != _tmpname)
      BFT_FREE(tmpname);
  }

  if (inp->index->n_files > 0)
    inp->f = inp->index->f[0];
  else
    inp->f = NULL;
}

/*----------------------------------------------------------------------------
 * Close the interface file.
 *
 * parameters:
 *   cs_io <-> kernel IO structure
 *----------------------------------------------------------------------------*/

static void
_file_close(cs_io_t  *cs_io)
{
  if (cs_io->f != NULL)
    cs_io->f = fvm_file_free(cs_io->f);
}

/*----------------------------------------------------------------------------
 * Echo pending section read or write
 *
 * parameters:
 *   cs_io --> kernel IO structure
 *----------------------------------------------------------------------------*/

static void
_echo_pre(const cs_io_t  *cs_io)
{
  assert(cs_io != NULL);

  switch(cs_io->mode) {

  case CS_IO_MODE_READ:
    bft_printf(_("\n  Section read on \"%s\":\n"),
               fvm_file_get_name(cs_io->f));
    break;

  case CS_IO_MODE_WRITE:
    bft_printf(_("\n  Section written on \"%s\":\n"),
               fvm_file_get_name(cs_io->f));
    break;

  default:
    assert(   cs_io->mode == CS_IO_MODE_READ
           || cs_io->mode == CS_IO_MODE_WRITE);
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Echo a section header
 *
 * parameters:
 *   sec_name  --> section name
 *   n_elts    --> number of elements
 *   elt_type  --> expected (final) element type
 *   type      --> element type in file
 *----------------------------------------------------------------------------*/

static void
_echo_header(const char      *sec_name,
             fvm_gnum_t       n_elts,
             fvm_datatype_t   type)
{
  /* Instructions */

  bft_printf(_("    section name:           \"%s\"\n"
               "    number of elements:     %lu\n"),
             sec_name, (unsigned long)n_elts);

  if (n_elts > 0) {

    char *type_name;

    switch(type) {
    case FVM_CHAR:
      type_name = _type_name_char;
      break;
    case FVM_INT32:
      type_name = _type_name_i4;
      break;
    case FVM_INT64:
      type_name = _type_name_i8;
      break;
    case FVM_UINT32:
      type_name = _type_name_u4;
      break;
    case FVM_UINT64:
      type_name = _type_name_u8;
      break;
    case FVM_FLOAT:
      type_name = _type_name_r4;
      break;
    case FVM_DOUBLE:
      type_name = _type_name_r8;
      break;
    default:
      assert(0);
      type_name = _type_name_none;
    }

    bft_printf(_("    element type name:      \"%s\"\n"), type_name);

  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Partial echo of a section's contents.
 *
 * FVM datatypes must have been converted to the corresponding
 * Code_Saturne compatible datatype before calling this function:
 *   FVM_CHAR                -> char
 *   FVM_INT32 / FVM_INT64   -> fvm_lnum_t / cs_int_t
 *   FVM_UINT32 / FVM_UINT64 -> fvm_gnum_t
 *   FVM_REAL / FVM_FLOAT    -> double / cs_real_t
 *
 * If global_num_start and global_num_end are > 0, the echo shows that
 * a different block is assigned to each processor. Otherwise, the full
 * data is replicated for each processor, and this is also indicated.
 *
 * parameters:
 *   echo             --> echo level
 *   n_elts           --> number of elements
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elt_type         --> element type
 *   elts             --> element buffer
 *----------------------------------------------------------------------------*/

static void
_echo_data(size_t           echo,
           fvm_file_off_t   n_elts,
           fvm_gnum_t       global_num_start,
           fvm_gnum_t       global_num_end,
           fvm_datatype_t   elt_type,
           const void      *elts)
{
  fvm_gnum_t  i;
  fvm_gnum_t  num_shift = 1;
  size_t  _n_elts = n_elts;
  fvm_file_off_t  echo_start = 0;
  fvm_file_off_t  echo_end = 0;
  const char *_loc_glob[] = {N_(" (local)"), ""};
  const char *loc_glob = _loc_glob[1];

  /* Instructions */

  if (n_elts == 0) return;

  if (cs_glob_n_ranks == 1)
    loc_glob = _loc_glob[0];
  else if (global_num_start > 0) {
    loc_glob = _loc_glob[0];
    num_shift = global_num_start;
  }

  if (global_num_start > 0 && global_num_end > 0) {
    assert(global_num_end >= global_num_start);
    _n_elts = global_num_end - global_num_start;
  }

  if (echo * 2 < _n_elts) {
    echo_end = echo;
    bft_printf(_("    %d first and last elements%s:\n"),
               echo, loc_glob);
  }
  else {
    echo_end = _n_elts;
    bft_printf(_("    elements%s:\n"), _(loc_glob));
  }

  /* Note that FVM datatypes will have been converted to
     the corresponding datatype if necessary before
     calling this function, hence the identical treatment
     for different cases. */

  do {

    switch (elt_type) {

    case FVM_INT32:
    case FVM_INT64:
      {
        const fvm_lnum_t *_elts = elts;

        for (i = echo_start ; i < echo_end ; i++)
          bft_printf("    %10lu : %12d\n",
                     (unsigned long)(i + num_shift), *(_elts + i));
      }
      break;

    case FVM_UINT32:
    case FVM_UINT64:
      {
        const fvm_gnum_t *_elts = elts;

        for (i = echo_start ; i < echo_end ; i++)
          bft_printf("    %10lu : %12lu\n",
                     (unsigned long)(i + num_shift),
                     (unsigned long)*(_elts + i));
      }
      break;

    case FVM_FLOAT:
    case FVM_DOUBLE:
      {
        const cs_real_t *_elts = elts;

        for (i = echo_start ; i < echo_end ; i++)
          bft_printf("    %10lu : %12.5e\n",
                     (unsigned long)(i + num_shift), *(_elts + i));
      }
      break;

    case FVM_CHAR:
      {
        const char *_elts = elts;

        for (i = echo_start ; i < echo_end ; i++) {
          if (*(_elts + i) != '\0')
            bft_printf("    %10lu : '%c'\n",
                       (unsigned long)(i + num_shift), *(_elts + i));
          else
            bft_printf("    %10lu : '\\0'\n",
                       (unsigned long)(i + num_shift));
        }
      }
      break;

    default:
      assert(0);

    }

    if (echo_end < n_elts) {
      bft_printf("    ..........   ............\n");
      echo_start = n_elts - echo;
      echo_end = n_elts;
    }
    else {
      assert(echo_end == n_elts);
      echo_end = n_elts + 1;
    }

  } while (echo_end <= n_elts);

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Convert read data.
 *
 * dest_type must have been set to the corresponding
 * Code_Saturne compatible datatype before calling this function and
 * conversion will be done accordingly:
 *   FVM_INT32 / FVM_INT64   -> fvm_lnum_t / cs_int_t
 *   FVM_UINT32 / FVM_UINT64 -> fvm_gnum_t
 *   FVM_REAL / FVM_FLOAT    -> double / cs_real_t
 *
 * parameters:
 *   buffer      --> input buffer
 *   dest        <-- output buffer
 *   n_elts      --> number of elements
 *   buffer_type --> input buffer element type
 *   dest_type   --> output buffer element type
 *----------------------------------------------------------------------------*/

static void
_cs_io_convert_read(void            *buffer,
                    void            *dest,
                    fvm_file_off_t   n_elts,
                    fvm_datatype_t   buffer_type,
                    fvm_datatype_t   dest_type)
{
  fvm_file_off_t ii;
  size_t buffer_type_size = fvm_datatype_size[buffer_type];

  assert(dest_type != buffer_type);

  /* Note that dest will have been set to the corresponding datatype
     before calling this function, hence the identical treatment
     for different cases. */

  switch(dest_type) {

  case FVM_INT32:
  case FVM_INT64:
    {
      fvm_lnum_t *_dest = dest;

      if (   buffer_type == FVM_INT32
          || buffer_type == FVM_INT64) {

        if (sizeof(long) == buffer_type_size) {
          long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
            _dest[ii] = _buffer[ii];
        }
        else if (sizeof(long long) == buffer_type_size) {
          long long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(int) == buffer_type_size) {
          int * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(short) == buffer_type_size) {
          short * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }

      }

      else if (   buffer_type == FVM_UINT32
               || buffer_type == FVM_UINT64) {

        if (sizeof(unsigned long) == buffer_type_size) {
          unsigned long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
            _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned long long) == buffer_type_size) {
          unsigned long long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned int) == buffer_type_size) {
          unsigned int * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned short) == buffer_type_size) {
          unsigned short * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }

      }

      assert(   sizeof(long) == buffer_type_size
             || sizeof(long long) == buffer_type_size
             || sizeof(int) == buffer_type_size
             || sizeof(short) == buffer_type_size);
    }
    break;

  case FVM_UINT32:
  case FVM_UINT64:
    {
      fvm_gnum_t *_dest = dest;

      if (   buffer_type == FVM_INT32
          || buffer_type == FVM_INT64) {

        if (sizeof(long) == buffer_type_size) {
          long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
            _dest[ii] = _buffer[ii];
        }
        else if (sizeof(long long) == buffer_type_size) {
          long long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(int) == buffer_type_size) {
          int * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(short) == buffer_type_size) {
          short * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }

      }

      else if (   buffer_type == FVM_UINT32
               || buffer_type == FVM_UINT64) {

        if (sizeof(unsigned long) == buffer_type_size) {
          unsigned long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
            _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned long long) == buffer_type_size) {
          unsigned long long * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned int) == buffer_type_size) {
          unsigned int * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }
        else if (sizeof(unsigned short) == buffer_type_size) {
          unsigned short * _buffer = buffer;
          for (ii = 0; ii < n_elts; ii++)
          _dest[ii] = _buffer[ii];
        }

      }

      assert(   sizeof(long) == buffer_type_size
             || sizeof(long long) == buffer_type_size
             || sizeof(int) == buffer_type_size
             || sizeof(short) == buffer_type_size);
    }
    break;

  case FVM_FLOAT:
    {
      cs_real_t *_dest = dest;
      double * _buffer = buffer;

      assert(buffer_type == FVM_DOUBLE);

      for (ii = 0; ii < n_elts; ii++)
        _dest[ii] = _buffer[ii];
    }
    break;

  case FVM_DOUBLE:
    {
      cs_real_t *_dest = dest;
      float * _buffer = buffer;

      assert(buffer_type == FVM_FLOAT);

      for (ii = 0; ii < n_elts; ii++)
        _dest[ii] = _buffer[ii];
    }
    break;

  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------
 * Read a section body.
 *
 * If global_num_start and global_num_end are > 0, a different block is
 * assigned to each processor. Otherwise, the full data is replicated
 * for each processor.
 *
 * If location_id > 0 and header->n_location_vals > 1, then
 * global_num_start and global_num_end will be based on location element
 * numbers, so the total number of values read equals
 * (global_num_end - global_num_start) * header->n_location_vals.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elts             <-> pointer to data array, or NULL
 *   inp              --> input kernel IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

static void *
_cs_io_read_body(const cs_io_sec_header_t  *header,
                 fvm_gnum_t                 global_num_start,
                 fvm_gnum_t                 global_num_end,
                 void                      *elts,
                 cs_io_t                   *inp)
{
  size_t  type_size = 0;
  fvm_file_off_t  n_vals = inp->n_vals;
  cs_bool_t  convert_type = false;
  void  *_elts = NULL;
  void  *_buf = NULL;
  size_t  stride = 1;

  assert(inp  != NULL);

  assert(header->n_vals == inp->n_vals);

  assert(global_num_end <= header->n_vals + 1);

  /* Choose global or block mode */

  if (header->n_location_vals > 1)
    stride = header->n_location_vals;

  if (global_num_start > 0 && global_num_end > 0) {
    assert(global_num_end >= global_num_start);
    n_vals = (global_num_end - global_num_start)*stride;
  }

  /* Datatype size given by FVM datatype */

  type_size = fvm_datatype_size[header->type_read];

  /* Assign or allocate */

  _elts = elts;

  if (_elts == NULL && n_vals != 0) {
    if (header->elt_type == FVM_CHAR)
      BFT_MALLOC(_elts, n_vals + 1, char);
    else
      BFT_MALLOC(_elts, n_vals*type_size, char);
  }

  /* Element values */

  if (n_vals != 0 && header->elt_type != header->type_read)
    convert_type = true;

  if (inp->data != NULL)
    _buf = NULL;
  else if (convert_type == true
           && (   fvm_datatype_size[header->type_read]
               != fvm_datatype_size[header->elt_type]))
    BFT_MALLOC(_buf, n_vals*type_size, char);
  else
    _buf = _elts;

  /* Read data from file */

  if (inp->data == NULL) {

    /* Position read pointer if necessary */

    if (inp->body_align > 0) {
      fvm_file_off_t offset = fvm_file_tell(inp->f);
      size_t ba = inp->body_align;
      offset += (ba - (offset % ba)) % ba;
      fvm_file_seek(inp->f, offset, FVM_FILE_SEEK_SET);
    }

    /* Read local or global values */

    if (global_num_start > 0 && global_num_end > 0)
      fvm_file_read_block(inp->f,
                          _buf,
                          type_size,
                          stride,
                          global_num_start,
                          global_num_end);
    else if (n_vals > 0)
      fvm_file_read_global(inp->f,
                           _buf,
                           type_size,
                           n_vals);

  }

  /* If data is embedded in header, simply point to it */

  else {

    if (global_num_start > 0 && global_num_end > 0)
      _buf =   ((unsigned char *)inp->data)
             + (  (global_num_start - 1) * stride
                * fvm_datatype_size[header->type_read]);
    else
      _buf = inp->data;

  }

  /* Convert data if necessary */

  if (convert_type == true) {
    _cs_io_convert_read(_buf,
                        _elts,
                        n_vals,
                        header->type_read,
                        header->elt_type);
    if (   inp->data == NULL
        && _buf != _elts)
      BFT_FREE(_buf);
  }
  else if (inp->data != NULL) {
    memcpy(_elts, _buf, n_vals*fvm_datatype_size[header->type_read]);
    _buf = NULL;
  }

  if (inp->data != NULL)  /* Reset for next read */
    inp->data = NULL;

  /* Add null character at end of string to ensure C-type string */

  if (n_vals != 0 && header->elt_type == FVM_CHAR)
    ((char *)_elts)[header->n_vals] = '\0';

  /* Optional echo */

  if (header->n_vals != 0 && inp->echo > CS_IO_ECHO_HEADERS)
    _echo_data(inp->echo,
               n_vals,
               (global_num_start-1)*stride + 1,
               (global_num_end-1)*stride + 1,
               header->elt_type,
               _elts);

  /* Return read values */

  return _elts;
}

/*----------------------------------------------------------------------------
 * Build an index for a kernel IO file structure in read mode.
 *
 * The magic string may be NULL, if we choose to ignore it.
 *
 * parameters:
 *   inp          <-> empty input kernel IO file structure
 *   name         <-- file name
 *   magic_string <-- magic string associated with file type
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
static void
_cs_io_initialize_with_index(cs_io_t       *inp,
                             const char    *file_name,
                             const char    *magic_string,
                             MPI_Comm       comm)
#else
static void
_cs_io_initialize_with_index(cs_io_t       *inp,
                             const char    *file_name,
                             const char    *magic_string)
#endif /* FVM_HAVE_MPI */
{
  cs_io_sec_header_t  h;
  int  end_reached = 0;

  /* Create interface file descriptor; do not use MPI-IO at this
     stage, as we only read global headers of limited size, and
     a "lighter" method than MPI-IO should be well adapted. */

#if defined(FVM_HAVE_MPI)
  _file_open(inp, file_name, magic_string, FVM_FILE_NO_MPI_IO, comm);
#else
  _file_open(inp, file_name, magic_string, 0);
#endif

  /* Update index file section */

  if (inp->index != NULL) {
    BFT_REALLOC(inp->index->f, inp->index->n_files + 1, fvm_file_t *);
    inp->index->f[inp->index->n_files] = inp->f;
    inp->index->n_files += 1;
  }

  /* Read headers to build index index */

  while (end_reached == 0) {

    end_reached = cs_io_read_header(inp, &h);

    if (end_reached == 0)
      _update_index_and_shift(inp, &h);

  }
}

/*----------------------------------------------------------------------------
 * Write padding with zero bytes if necessary to ensure alignment.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The section data on other ranks is ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   align     <-- required alignment
 *   outp      --> output kernel IO structure
 *----------------------------------------------------------------------------*/

static void
_write_padding(size_t    align,
               cs_io_t  *outp)
{
  if (align > 0) {

    fvm_file_off_t offset = fvm_file_tell(outp->f);
    fvm_file_off_t add_offset = (align - (offset % align)) % align;

    if (add_offset > 0) {

      size_t pad_size = add_offset;
      size_t n_written = 0;

      if (pad_size > outp->buffer_size) {
        while (pad_size > outp->buffer_size)
          outp->buffer_size *=2;
        BFT_REALLOC(outp->buffer, outp->buffer_size, unsigned char);
      }

      memset(outp->buffer, 0, pad_size);

      n_written = fvm_file_write_global(outp->f,
                                        outp->buffer,
                                        1,
                                        pad_size);

      if (pad_size != n_written)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error writing %lu bytes to file \"%s\"."),
                  (unsigned long)pad_size, fvm_file_get_name(outp->f));
    }
  }
}

/*----------------------------------------------------------------------------
 * Write a section header, with possibly embedded data.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The section data on other ranks is ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   section_name     <-- section name
 *   n_vals           <-- total number of values
 *   location_id      <-- id of associated location, or 0
 *   index_id         <-- id of associated index, or 0
 *   n_location_vals  <-- number of values per location
 *   elt_type         <-- element type
 *   elts             <-- pointer to element data, if it may be embedded
 *   outp             --> output kernel IO structure
 *
 * returns:
 *   true if element data is embedded, false otherwise
 *----------------------------------------------------------------------------*/

static cs_bool_t
_write_header(const char      *sec_name,
              fvm_gnum_t       n_vals,
              size_t           location_id,
              size_t           index_id,
              size_t           n_location_vals,
              fvm_datatype_t   elt_type,
              const void      *elts,
              cs_io_t         *outp)
{
  size_t name_size = 0, name_pad_size = 0;
  fvm_file_off_t write_size = 0;
  fvm_file_off_t data_size = n_vals * fvm_datatype_size[elt_type];
  fvm_file_off_t header_vals[6];
  size_t n_written = 0;

  cs_bool_t embed = false;

  assert(outp != NULL);

  if (outp->echo >= CS_IO_ECHO_HEADERS)
    _echo_pre(outp);

  _write_padding(outp->header_align, outp); /* Pad if necessary */

  /* Prepare header data */
  /*---------------------*/

  header_vals[0] = 56;
  header_vals[1] = n_vals;
  header_vals[2] = location_id;
  header_vals[3] = index_id;

  if (n_vals > 0)
    header_vals[4] = n_location_vals;
  else
    header_vals[4] = 0;

  name_size = strlen(sec_name);
  name_pad_size = 8-(name_size%8); /* At least 1 NULL
                                      character with this rule */

  header_vals[5] = name_size + name_pad_size;
  header_vals[0] += (name_size + name_pad_size);

  /* Decide if data is to be embedded */

  if (   n_vals > 0
      && elts != NULL
      && (header_vals[0] + data_size <= (fvm_file_off_t)(outp->header_size))) {
    header_vals[0] += data_size;
    embed = true;
  }

  /* Ensure buffer is big enough for data */

  if (header_vals[0] > (fvm_file_off_t)(outp->buffer_size)) {
    while (header_vals[0] > (fvm_file_off_t)(outp->buffer_size))
      outp->buffer_size *=2;
    BFT_REALLOC(outp->buffer, outp->buffer_size, unsigned char);
  }

  memset(outp->buffer, 0, outp->buffer_size);

  /* Now build buffer */

  _convert_from_offset(outp->buffer, header_vals, 6);

  if (fvm_file_get_swap_endian(outp->f) == 1)
    bft_file_swap_endian(outp->buffer, outp->buffer, 8, 6);

  /* Element type name */

  outp->type_name = (char *)outp->buffer + 48;

  if (n_vals > 0) {
    switch(elt_type) {
    case FVM_FLOAT:
      outp->type_name[0] = 'r';
      if (sizeof(float) == 4)
        outp->type_name[1] = '4';
      else
        outp->type_name[1] = '8';
      assert(sizeof(float) == 4 || sizeof(float) == 8);
      break;
    case FVM_DOUBLE:
      outp->type_name[0] = 'r';
      if (sizeof(double) == 8)
        outp->type_name[1] = '8';
      else {
        outp->type_name[1] = '1';
        outp->type_name[2] = '6';
      }
      assert(sizeof(double) == 8 || sizeof(float) == 16);
      break;
    case FVM_INT32:
      outp->type_name[0] = 'i';
      outp->type_name[1] = '4';
      break;
    case FVM_INT64:
      outp->type_name[0] = 'i';
      outp->type_name[1] = '4';
      break;
    case FVM_UINT32:
      outp->type_name[0] = 'u';
      outp->type_name[1] = '4';
      break;
    case FVM_UINT64:
      outp->type_name[0] = 'u';
      outp->type_name[1] = '4';
      break;
    case FVM_CHAR:
      outp->type_name[0] = 'c';
      outp->type_name[1] = ' ';
      break;
    default:
      break;
    }
  }

  if (embed == true)
    outp->type_name[7] = 'e';

  /* Section name */

  strcpy((char *)(outp->buffer) + 56, sec_name);

  if (embed == true) {

    unsigned char *data =   (unsigned char *)(outp->buffer)
                          + (56 + name_size + name_pad_size);

    memcpy(data, elts, data_size);

    if (   fvm_file_get_swap_endian(outp->f) == 1
        && fvm_datatype_size[elt_type] > 1)
      bft_file_swap_endian(data, data, fvm_datatype_size[elt_type], n_vals);
  }

  /* Now write header data */

  write_size = CS_MAX((fvm_file_off_t)(outp->header_size), header_vals[0]);

  n_written = fvm_file_write_global(outp->f,
                                    outp->buffer,
                                    1,
                                    write_size);

  if (write_size != (fvm_file_off_t)n_written)
    bft_error(__FILE__, __LINE__, 0,
              _("Error writing %lu bytes to file \"%s\"."),
              (unsigned long)write_size, fvm_file_get_name(outp->f));

  if (outp->echo >= CS_IO_ECHO_HEADERS)
    _echo_header(sec_name, n_vals, elt_type);

  return embed;
}

/*----------------------------------------------------------------------------
 * Dump a kernel IO file handle's metadata.
 *
 * parameters:
 *   idx  <-- kernel IO index
 *----------------------------------------------------------------------------*/

static void
_dump_index(const cs_io_sec_index_t  *idx)
{
  size_t ii;

  assert(idx != NULL);

  bft_printf(_(" %lu indexed records:\n"
               "   (name, n_vals, location_id, index_id, n_loc_vals, type, "
               "embed, file_id, offset)\n\n"), (unsigned long)(idx->size));

  for (ii = 0; ii < idx->size; ii++) {

    char embed = 'n';
    fvm_file_off_t *h_vals = idx->h_vals + ii*8;
    const char *name = idx->names + h_vals[4];

    if (h_vals[5] > 0)
      embed = 'y';

    bft_printf(_(" %40s %10lu %2u %2u %2u %6s %c %2u %ld\n"),
               name, (unsigned long)(h_vals[0]),
               (unsigned)(h_vals[1]), (unsigned)(h_vals[2]),
               (unsigned)(h_vals[3]), fvm_datatype_name[h_vals[6]],
               embed, (unsigned)(h_vals[7]),
               (long)(idx->offset[ii]));

  }

  bft_printf(_("\n %u associated file(s):\n"), (unsigned)(idx->n_files));

  for (ii = 0; ii < idx->n_files; ii++)
    bft_printf(_("  \"%s\"\n"), fvm_file_get_name(idx->f[ii]));

  bft_printf("\n");
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a kernel IO file structure.
 *
 * The magic string may be NULL only in read mode, if we choose to ignore it.
 *
 * parameters:
 *   name         <-- file name
 *   magic_string <-- magic string associated with file type
 *   mode         <-- read or write
 *   hints        <-- optional flags for file access method (see fvm_file.h)
 *   echo         <-- echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 int            hints,
                 long           echo,
                 MPI_Comm       comm)
#else
cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 int            hints,
                 long           echo)
#endif /* FVM_HAVE_MPI */
{
  cs_io_t  *cs_io =_cs_io_create(mode, echo);

  /* Info on interface creation */

  if (echo >= CS_IO_ECHO_OPEN_CLOSE) {
    if (mode == CS_IO_MODE_READ)
      bft_printf(_("\n Reading file:        %s\n"), file_name);
    else
      bft_printf(_("\n Writing file:        %s\n"), file_name);
    bft_printf_flush();
  }

  /* Create interface file descriptor */

#if defined(FVM_HAVE_MPI)
  _file_open(cs_io, file_name, magic_string, hints, comm);
#else
  _file_open(cs_io, file_name, magic_string, hints);
#endif

  return cs_io;
}

/*----------------------------------------------------------------------------
 * Initialize a kernel IO file structure in read mode, building an index.
 *
 * The magic string may be NULL, if we choose to ignore it.
 *
 * parameters:
 *   name         <-- file name
 *   magic_string <-- magic string associated with file type
 *   hints        <-- optional flags for file access method (see fvm_file.h)
 *   echo         <-- echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   pointer to kernel IO structure
 *----------------------------------------------------------------------------*/

#if defined(FVM_HAVE_MPI)
cs_io_t *
cs_io_initialize_with_index(const char    *file_name,
                            const char    *magic_string,
                            int            hints,
                            long           echo,
                            MPI_Comm       comm)
#else
cs_io_t *
cs_io_initialize_with_index(const char    *file_name,
                            const char    *magic_string,
                            int            hints,
                            long           echo)
#endif /* FVM_HAVE_MPI */
{
  int retval = 0;

  cs_io_t  *inp =_cs_io_create(CS_IO_MODE_READ, echo);

  /* Info on interface creation */

  if (echo >= CS_IO_ECHO_OPEN_CLOSE) {
    bft_printf(_("\n Reading file:        %s\n"), file_name);
    bft_printf_flush();
  }

  /* Initialize index */

  _create_index(inp);

  /* Test for legacy restart format first */

#if defined(FVM_HAVE_MPI)
  retval = _file_legacy_restart_index(inp, file_name, comm);
#else
  retval = _file_legacy_restart_index(inp, file_name);

#endif

  if (retval == 0) {

#if defined(FVM_HAVE_MPI)
    _cs_io_initialize_with_index(inp, file_name, magic_string, comm);
#else
    _cs_io_initialize_with_index(inp, file_name, magic_string);
#endif

  }

  /* Now reopen all indexed files using hints */

#if defined(FVM_HAVE_MPI)
  _file_reopen_read(inp, hints, comm);
#else
  _file_reopen_read(inp, hints);
#endif

  return inp;
}

/*----------------------------------------------------------------------------
 * Free a kernel IO file structure, closing the associated file.
 *
 * parameters:
 *   cs_io <-> kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_finalize(cs_io_t **cs_io)
{
  cs_io_t *_cs_io = *cs_io;

  if(_cs_io->mode == CS_IO_MODE_WRITE)
    cs_io_write_global("EOF", 0, 0, 0, 0, FVM_DATATYPE_NULL, NULL, _cs_io);

  /* Info on closing of interface file */

  if (_cs_io->echo >= CS_IO_ECHO_OPEN_CLOSE) {
    if (_cs_io->mode == CS_IO_MODE_READ)
      bft_printf(_(" Finished reading:    %s\n"),
                 fvm_file_get_name(_cs_io->f));
    else
      bft_printf(_(" Finished writing:    %s\n"),
                 fvm_file_get_name(_cs_io->f));
    bft_printf_flush();
  }

  if (_cs_io->index != NULL)
    _destroy_index(_cs_io);

  _file_close(_cs_io);

  _cs_io->buffer_size = 0;
  BFT_FREE(_cs_io->buffer);

  BFT_FREE(*cs_io);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a kernel IO structure's name.
 *
 * parameters:
 *   cs_io <-- kernel IO structure
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_name(const cs_io_t  *cs_io)
{
  assert(cs_io != NULL);

  return(fvm_file_get_name(cs_io->f));
}

/*----------------------------------------------------------------------------
 * Return the number of indexed entries in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *
 * returns:
 *   size of index if present, 0 otherwise,
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_index_size(const cs_io_t  *inp)
{
  size_t retval = 0;

  if (inp != NULL && inp->index != NULL)
    retval = inp->index->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the name of an indexed section in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *   id  <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   pointer to section name if id in index range, NULL otherwise
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_indexed_sec_name(const cs_io_t  *inp,
                           size_t          id)
{
  const char *retval = NULL;

  if (inp != NULL && inp->index != NULL) {
    if (id < inp->index->size) {
      size_t name_id = inp->index->h_vals[8*id + 4];
      retval = inp->index->names + name_id;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return header data for an indexed section in a kernel IO structure.
 *
 * parameters:
 *   inp <-- input kernel IO structure
 *   id  <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   section header data (if id not in index range, fields set to zero)
 *----------------------------------------------------------------------------*/

cs_io_sec_header_t
cs_io_get_indexed_sec_header(const cs_io_t  *inp,
                             size_t          id)
{
  cs_io_sec_header_t h;

  h.sec_name = NULL;

  if (inp != NULL && inp->index != NULL) {
    if (id < inp->index->size) {

      size_t name_id = inp->index->h_vals[8*id + 4];

      h.sec_name = inp->index->names + name_id;

      h.n_vals          = inp->index->h_vals[8*id];
      h.location_id     = inp->index->h_vals[8*id + 1];
      h.index_id        = inp->index->h_vals[8*id + 2];
      h.n_location_vals = inp->index->h_vals[8*id + 3];
      h.type_read       = inp->index->h_vals[8*id + 6];
      h.elt_type        = _type_read_to_elt_type(h.type_read);
    }
  }

  if (h.sec_name == NULL) {
    h.n_vals          = 0;
    h.location_id     = 0;
    h.index_id        = 0;
    h.n_location_vals = 0;
    h.type_read       = FVM_DATATYPE_NULL;
    h.elt_type        = h.type_read;
  }

  return h;
}

/*----------------------------------------------------------------------------
 * Return a kernel IO structure's echo (verbosity) level.
 *
 * parameters:
 *   cs_io --> kernel IO structure
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_echo(const cs_io_t  *cs_io)
{
  assert(cs_io != NULL);

  return (size_t)(cs_io->echo);
}

/*----------------------------------------------------------------------------
 * Read a section header.
 *
 * The header values remain valid until the next section header read
 * or the file is closed.
 *
 * parameters:
 *   inp    --> input kernel IO structure
 *   header <-- header structure
 *
 * returns:
 *   0 if a header was read, 1 in case of error or end-of-file
 *----------------------------------------------------------------------------*/

int
cs_io_read_header(cs_io_t             *inp,
                  cs_io_sec_header_t  *header)
{
  int type_name_error = 0;
  fvm_file_off_t body_size = 0;
  fvm_file_off_t header_vals[6];
  size_t n_read = 0;

  assert(inp != NULL);

  if (inp->echo >= CS_IO_ECHO_HEADERS)
    _echo_pre(inp);

  /* Position read pointer if necessary */
  /*------------------------------------*/

  if (inp->header_align > 0) {
    size_t ha = inp->header_align;
    fvm_file_off_t offset = fvm_file_tell(inp->f);
    fvm_file_off_t add_offset = (ha - (offset % ha)) % ha;
    if (add_offset > 0) {
      int errcode = 0;
      errcode = fvm_file_seek(inp->f, add_offset, FVM_FILE_SEEK_CUR);
      if (errcode != 0)
        return 1;
    }
  }

  inp->n_vals = 0;

  /* Read header */
  /*-------------*/

  n_read = fvm_file_read_global(inp->f, inp->buffer, 1, inp->header_size);

  if (n_read < inp->header_size)
    return 1;

  if (fvm_file_get_swap_endian(inp->f) == 1)
    bft_file_swap_endian(inp->buffer, inp->buffer, 8, 6);

  _convert_to_offset(inp->buffer, header_vals, 6);

  if (header_vals[0] > (fvm_file_off_t)(inp->header_size)) {

    size_t n_add = header_vals[0] - inp->header_size;

    if (header_vals[0] > (fvm_file_off_t)(inp->buffer_size)) {
      while (header_vals[0] > (fvm_file_off_t)(inp->buffer_size))
        inp->buffer_size *=2;
      BFT_REALLOC(inp->buffer, inp->buffer_size, unsigned char);
    }

    n_read = fvm_file_read_global(inp->f,
                                  inp->buffer + inp->header_size,
                                  1,
                                  n_add);

    if (n_read < n_add)
      return 1;
  }

  /* Set pointers to data fields */

  inp->n_vals = header_vals[1];
  inp->location_id = header_vals[2];
  inp->index_id = header_vals[3];
  inp->n_loc_vals = header_vals[4];
  inp->type_size = 0;
  inp->data = NULL;
  inp->type_name = (char *)(inp->buffer + 48);
  inp->sec_name = (char *)(inp->buffer + 56);

  if (header_vals[1] > 0 && inp->type_name[7] == 'e')
    inp->data = inp->buffer + 56 + header_vals[5];

  inp->type_size = 0;

  /* Return immediately if we have an end-of file marker */

  if ((inp->n_vals == 0) && (strcmp(inp->sec_name, "EOF") == 0))
    return 1;

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
      bft_error(__FILE__, __LINE__, 0,
                _("Type \"%s\" is not known\n"
                  "Known types: \"c \", \"i4\", \"i8\", \"u4\", \"u8\", "
                  "\"r4\", \"r8\"."), inp->type_name);

    else if (inp->data == NULL)
      body_size = inp->type_size*inp->n_vals;

    else if (fvm_file_get_swap_endian(inp->f) == 1 && inp->type_size > 1)
      bft_file_swap_endian(inp->data,
                           inp->data,
                           inp->type_size,
                           inp->n_vals);
  }

  /* Set externally visible header values */
  /*--------------------------------------*/

  header->sec_name = inp->sec_name;
  header->n_vals = inp->n_vals;
  header->location_id = inp->location_id;
  header->index_id = inp->index_id;
  header->n_location_vals = inp->n_loc_vals;

  /* Initialize data type */
  /*----------------------*/

  assert(sizeof(unsigned long) == 8 || sizeof(unsigned long long) == 8);

  if (header->n_vals != 0) {

    const char *elt_type_name = inp->type_name;

    if (   strcmp(elt_type_name, _type_name_i4) == 0
        || strcmp(elt_type_name, "i ") == 0)
      header->type_read = FVM_INT32;

    else if (strcmp(elt_type_name, _type_name_i8) == 0)
      header->type_read = FVM_INT64;

    else if (strcmp(elt_type_name, _type_name_u4) == 0)
      header->type_read = FVM_UINT32;

    else if (strcmp(elt_type_name, _type_name_u8) == 0)
      header->type_read = FVM_UINT64;

    else if (strcmp(elt_type_name, _type_name_r4) == 0)
      header->type_read = FVM_FLOAT;

    else if (strcmp(elt_type_name, _type_name_r8) == 0)
      header->type_read = FVM_DOUBLE;

    else if (strcmp(elt_type_name, _type_name_char) == 0)
      header->type_read = FVM_CHAR;

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error reading file: \"%s\".\n"
                  "Data type \"%s\" is not recognized."),
                fvm_file_get_name(inp->f), elt_type_name);

    header->elt_type = _type_read_to_elt_type(header->type_read);
  }

  /* Possible echo */

  if (inp->echo >= CS_IO_ECHO_HEADERS)
    _echo_header(header->sec_name,
                 header->n_vals,
                 header->elt_type);

  return 0;
}

/*----------------------------------------------------------------------------
 * Set a kernel IO's state so as to be ready to read an indexed section.
 *
 * The header values and position in the file are set so as to be equivalent
 * to those they would have if the corresponding header had just been read.
 *
 * parameters:
 *   inp    <-> input kernel IO structure
 *   header --> associated header
 *   id     <-- id of section in index (0 to n-1 numbering)
 *
 * returns:
 *   0 in case of success, 1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_io_set_indexed_position(cs_io_t             *inp,
                           cs_io_sec_header_t  *header,
                           size_t               id)
{
  int retval = 0;

  /* Return immediately if out of range */

  if (inp == NULL || inp->index == NULL)
    return 1;
  if (id >= inp->index->size)
    return 1;

  header->sec_name = inp->index->names + inp->index->h_vals[8*id + 4];

  header->n_vals          = inp->index->h_vals[8*id];
  header->location_id     = inp->index->h_vals[8*id + 1];
  header->index_id        = inp->index->h_vals[8*id + 2];
  header->n_location_vals = inp->index->h_vals[8*id + 3];
  header->type_read       = inp->index->h_vals[8*id + 6];
  header->elt_type        = _type_read_to_elt_type(header->type_read);

  inp->n_vals      = header->n_vals;
  inp->location_id = header->location_id;
  inp->index_id    = header->index_id;
  inp->n_loc_vals  = header->n_location_vals;
  inp->type_size   = fvm_datatype_size[header->type_read];

  /* The following values are not taken from the header buffer as
     usual, but are base on the index */

  strcpy((char *)(inp->buffer + 56), header->sec_name);
  inp->sec_name = (char *)(inp->buffer + 56);
  inp->type_name = NULL; /* should not be needed now that datatype is known */

  /* Non-embedded values */

  if (inp->index->h_vals[8*id + 5] == 0) {
    size_t file_id = inp->index->h_vals[8*id + 7];
    fvm_file_off_t offset = inp->index->offset[id];
    inp->f = inp->index->f[file_id];
    retval = fvm_file_seek(inp->f, offset, FVM_FILE_SEEK_SET);
  }

  /* Embedded values */

  else {
    size_t data_id = inp->index->h_vals[8*id + 5] - 1;
    unsigned char *_data = inp->index->data + data_id;
    inp->data = _data;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Set a section's final data type to fvm_lnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   cs_io  --> kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_lnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *cs_io)
{
  assert(header != NULL);

  if (   header->type_read != FVM_INT32
      && header->type_read != FVM_INT64
      && header->type_read != FVM_UINT32
      && header->type_read != FVM_UINT64)
    bft_error(__FILE__, __LINE__, 0,
              _("Error reading file: \"%s\".\n"
                "Type expected for section: "
                "\"%s\" is a signed integer\n"
                "and is not convertible from type read: \"%s\"."),
              fvm_file_get_name(cs_io->f), cs_io->type_name);

  assert(sizeof(fvm_lnum_t) == 4 || sizeof(fvm_lnum_t) == 8);

  if (sizeof(fvm_lnum_t) == 4)
    header->elt_type = FVM_INT32;
  else
    header->elt_type = FVM_INT64;
}

/*----------------------------------------------------------------------------
 * Set a section's final data type to fvm_gnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   cs_io  --> kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_gnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *cs_io)
{
  assert(header != NULL);

  if (   header->type_read != FVM_INT32
      && header->type_read != FVM_INT64
      && header->type_read != FVM_UINT32
      && header->type_read != FVM_UINT64)
    bft_error(__FILE__, __LINE__, 0,
              _("Error reading file: \"%s\".\n"
                "Type expected for section: "
                "\"%s\" is an unsigned integer\n"
                "and is not convertible from type read: \"%s\"."),
              fvm_file_get_name(cs_io->f), cs_io->type_name);

  assert(sizeof(fvm_gnum_t) == 4 || sizeof(fvm_gnum_t) == 8);

  if (sizeof(fvm_gnum_t) == 4)
    header->elt_type = FVM_UINT32;
  else
    header->elt_type = FVM_UINT64;
}

/*----------------------------------------------------------------------------
 * Check that a section's final data type corresponds to cs_real_t.
 *
 * parameters:
 *   header <-- header structure
 *   cs_io  --> kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_assert_cs_real(const cs_io_sec_header_t  *header,
                     const cs_io_t             *cs_io)
{
  assert(header != NULL);

  if (   header->elt_type != FVM_FLOAT
      && header->elt_type != FVM_DOUBLE)
    bft_error(__FILE__, __LINE__, 0,
              _("Error reading file: \"%s\".\n"
                "Type expected for section: \"%s\"\n"
                "is \"r4\" or \"r8\" (real), and not \"%s\"."),
              fvm_file_get_name(cs_io->f), cs_io->type_name);
}

/*----------------------------------------------------------------------------
 * Read a section body and replicate it to all processors.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   elts             <-> pointer to data array, or NULL
 *   cs_io            --> kernel IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_global(const cs_io_sec_header_t  *header,
                  void                      *elts,
                  cs_io_t                   *cs_io)
{
  return _cs_io_read_body(header, 0, 0, elts, cs_io);
}

/*----------------------------------------------------------------------------
 * Read a section body, assigning a different block to each processor.
 *
 * If location_id > 0 and header->n_location_vals > 1, then
 * global_num_start and global_num_end will be based on location element
 * numbers, so the total number of values read equals
 * (global_num_end - global_num_start) * header->n_location_vals.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elts             <-> pointer to data array, or NULL
 *   cs_io            --> kernel IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_block(const cs_io_sec_header_t  *header,
                 fvm_gnum_t                 global_num_start,
                 fvm_gnum_t                 global_num_end,
                 void                      *elts,
                 cs_io_t                   *cs_io)
{
  assert(global_num_start > 0);
  assert(global_num_end >= global_num_start);

  return _cs_io_read_body(header,
                          global_num_start,
                          global_num_end,
                          elts,
                          cs_io);
}

/*----------------------------------------------------------------------------
 * Read a section body, assigning a different block to each processor,
 * when the body corresponds to an index.
 *
 * In serial mode, this function behaves just like cs_io_read_block(),
 * except that it allows only unsigned integer values (fvm_gnum_t).
 *
 * In parallel mode, global_num_end should be set to the past-the-end value
 * of the base data block, the same as for regular data (and not increased
 * by 1 for the last rank, as this will be handled internally).
 * On each rank, the buffer size should be:
 * global_num_end - global_num_start + 1, as the past-the end index
 * for the local block is added automatically.
 *
 * If the array intended to receive the data already exists, we pass an
 * "elt" pointer to this array; this same pointer is then returned.
 * Otherwise, if this pointer is passed as NULL, memory is allocated
 * by this function, and the corresponding pointer is returned. It is
 * the caller's responsibility to free this array.
 *
 * parameters:
 *   header           <-- header structure
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elts             <-> pointer to data array, or NULL
 *   cs_io            --> kernel IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_index_block(cs_io_sec_header_t  *header,
                       fvm_gnum_t           global_num_start,
                       fvm_gnum_t           global_num_end,
                       fvm_gnum_t          *elts,
                       cs_io_t             *cs_io)
{
  fvm_gnum_t _global_num_start = global_num_start;
  fvm_gnum_t _global_num_end = global_num_end;
  fvm_gnum_t *retval = NULL;
  cs_bool_t last_data_rank = false;
  cs_bool_t past_last_data_rank = false;

  assert(global_num_start > 0);
  assert(global_num_end >= global_num_start);

  /* Check type */

  cs_io_set_fvm_gnum(header, cs_io);

  /* Increase _global_num_end by 1 for the last rank */

  if (header->n_vals == global_num_end) {

    _global_num_end += 1;
    last_data_rank = true;

    /* Also shift start values for possibly empty
       blocks past the last rank reading data */

    if (global_num_end <= global_num_start) {
      _global_num_start += 1;
      past_last_data_rank = true;
    }

  }

  retval = _cs_io_read_body(header,
                            _global_num_start,
                            _global_num_end,
                            elts,
                            cs_io);

  /* Exchange past-the-end values */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int needs_safe_algo_loc = 0;
    int needs_safe_algo = 0;
    int rank = cs_glob_rank_id;
    int send_rank = rank - 1;
    int recv_rank = rank + 1;
    MPI_Comm comm = cs_glob_mpi_comm;
    MPI_Status  status;
    fvm_gnum_t  past_last_recv = 0;
    fvm_gnum_t  past_last_send = 0;

    /* Prepare for MPI_Sendrecv */

    if (last_data_rank == true)
      recv_rank = MPI_PROC_NULL;

    if (past_last_data_rank == true) {
      send_rank = MPI_PROC_NULL;
      recv_rank = MPI_PROC_NULL;
    }

    if (rank == 0)
      send_rank = MPI_PROC_NULL;

    /* For empty blocks, past_last_send was initialized to 0,
       and this situation will lead to resorting to a safer
       but slower algorithm */

    if (global_num_end > global_num_start)
      past_last_send = retval[0];

    MPI_Sendrecv(&past_last_send, 1, FVM_MPI_GNUM, send_rank, CS_IO_MPI_TAG,
                 &past_last_recv, 1, FVM_MPI_GNUM, recv_rank, CS_IO_MPI_TAG,
                 comm, &status);

    if (recv_rank != MPI_PROC_NULL && past_last_recv == 0)
      needs_safe_algo_loc = 1;

    /* Check that everything is OK (i.e. no empty intermediate blocks) */

    MPI_Allreduce(&needs_safe_algo_loc, &needs_safe_algo, 1, MPI_INT, MPI_MAX,
                  comm);

    if (needs_safe_algo == 1) {

      int n_ranks = 1;
      int i;
      fvm_gnum_t *past_last = NULL;

      MPI_Comm_size(comm, &n_ranks);

      if (rank == 0)
        BFT_MALLOC(past_last, n_ranks, fvm_gnum_t);

      MPI_Gather(&past_last_send, 1, FVM_MPI_GNUM,
                 &past_last, 1, FVM_MPI_GNUM,
                 0, comm);

      /* Propagate values from higher ranks if necessary */

      for (i = n_ranks - 1; i > 1; i--) {
        if (past_last[i-1] == 0)
          past_last[i-1] = past_last[i];
      }

      MPI_Scatter(&past_last, 1, FVM_MPI_GNUM,
                  &past_last_recv, 1, FVM_MPI_GNUM,
                  0, comm);

      if (rank == 0)
        BFT_FREE(past_last);

    } /* End of condition on safer algorithm */

    if (last_data_rank == false && global_num_end > global_num_start)
      retval[global_num_end - global_num_start] = past_last_recv;

  }

  if (   header->n_vals != 0 && header->n_vals != global_num_end
      && cs_io->echo > CS_IO_ECHO_HEADERS)
    bft_printf(_("    first element for next rank:\n"
                 "    %10lu : %12d\n"),
               (unsigned long)(global_num_end),
               (unsigned long)retval[global_num_end - global_num_start]);

#endif /* defined(HAVE_MPI) */

  return retval;
}

/*----------------------------------------------------------------------------
 * Write a global section.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The section data on other ranks is ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   section_name     <-- section name
 *   n_vals           <-- total number of values
 *   location_id      <-- id of associated location, or 0
 *   index_id         <-- id of associated index, or 0
 *   n_location_vals  <-- number of values per location
 *   elt_type         <-- element type
 *   elts             <-- pointer to element data
 *   outp             <-> output kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_write_global(const char      *sec_name,
                   fvm_gnum_t       n_vals,
                   size_t           location_id,
                   size_t           index_id,
                   size_t           n_location_vals,
                   fvm_datatype_t   elt_type,
                   const void      *elts,
                   cs_io_t         *outp)
{
  cs_bool_t embed = false;

  if (outp->echo >= CS_IO_ECHO_HEADERS)
    _echo_header(sec_name, n_vals, elt_type);

  embed = _write_header(sec_name,
                        n_vals,
                        location_id,
                        index_id,
                        n_location_vals,
                        elt_type,
                        elts,
                        outp);

  if (n_vals > 0 && embed == false) {

    size_t n_written = 0;

    _write_padding(outp->body_align, outp);

    n_written = fvm_file_write_global(outp->f,
                                      elts,
                                      fvm_datatype_size[elt_type],
                                      n_vals);

    if (n_vals != n_written)
      bft_error(__FILE__, __LINE__, 0,
                _("Error writing %lu bytes to file \"%s\"."),
                (unsigned long)n_vals, fvm_file_get_name(outp->f));

  }

  if (n_vals != 0 && outp->echo > CS_IO_ECHO_HEADERS)
    _echo_data(outp->echo, n_vals, 1, n_vals + 1, elt_type, elts);
}

/*----------------------------------------------------------------------------
 * Write a section to file, each associated process providing a contiguous
 * of the section's body.
 *
 * Each process should provide a (possibly empty) block of the body,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * If location_id > 0 and n_location_vals > 1, then global_num_start
 * and global_num_end will be based on location element numbers, so the
 * total number of values read equals
 * (global_num_end - global_num_start) * header->n_location_vals.
 *
 * This function is intended to be used mainly data that is already of
 * copy of original data (such as data that has been redistributed across
 * processors just for the sake of output), or that is to be deleted after
 * writing, so it may modify the values in its input buffer (notably to
 * convert from little-endian to big-endian of vice-versa if necessary).
 *
 * parameters:
 *   section_name     <-- section name
 *   n_g_elts         <-- number of global elements (locations)
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *   location_id      <-- id of associated location, or 0
 *   index_id         <-- id of associated index, or 0
 *   n_location_vals  <-- number of values per location
 *   elt_type         <-- element type
 *                        (1 to n numbering)
 *   elts             <-- pointer to element data
 *   outp             <-> output kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_write_block_buffer(const char      *sec_name,
                         fvm_gnum_t       n_g_elts,
                         fvm_gnum_t       global_num_start,
                         fvm_gnum_t       global_num_end,
                         size_t           location_id,
                         size_t           index_id,
                         size_t           n_location_vals,
                         fvm_datatype_t   elt_type,
                         void            *elts,
                         cs_io_t         *outp)
{
  size_t n_written = 0;
  size_t n_g_vals = n_g_elts;
  size_t n_vals = global_num_end - global_num_start;
  size_t stride = 1;

  if (n_location_vals > 1) {
    stride = n_location_vals;
    n_g_vals *= n_location_vals;
    n_vals *= n_location_vals;
  }

  _write_header(sec_name,
                n_g_vals,
                location_id,
                index_id,
                n_location_vals,
                elt_type,
                NULL,
                outp);

  _write_padding(outp->body_align, outp);

  n_written = fvm_file_write_block_buffer(outp->f,
                                          elts,
                                          fvm_datatype_size[elt_type],
                                          stride,
                                          global_num_start,
                                          global_num_end);

  if (n_vals != n_written)
    bft_error(__FILE__, __LINE__, 0,
              _("Error writing %lu bytes to file \"%s\"."),
              (unsigned long)n_vals, fvm_file_get_name(outp->f));

  if (n_vals != 0 && outp->echo > CS_IO_ECHO_HEADERS)
    _echo_data(outp->echo, n_g_vals,
               (global_num_start-1)*stride + 1,
               (global_num_end -1)*stride + 1,
               elt_type, elts);
}

/*----------------------------------------------------------------------------
 * Print information on default options for file access.
 *----------------------------------------------------------------------------*/

void
cs_io_defaults_info(void)
{
  cs_bool_t  mpi_io = false;
  const char *fmt = N_("  I/O mode:          %s\n");

#if defined(FVM_HAVE_MPI_IO)

  if (cs_glob_n_ranks > 1) {
    if (cs_glob_io_hints & FVM_FILE_EXPLICIT_OFFSETS) {
      bft_printf(_(fmt), _("MPI-IO, explicit offsets"));
      mpi_io = true;
    }
    else if (cs_glob_io_hints & FVM_FILE_INDIVIDUAL_POINTERS) {
      bft_printf(_(fmt), _("MPI-IO, individual file pointers"));
      mpi_io = true;
    }
    if (mpi_io == false || (cs_glob_io_hints & FVM_FILE_NO_MPI_IO))
      bft_printf(_(fmt), _("serial IO\n\n"));
  }
#endif
}

/*----------------------------------------------------------------------------
 * Set the default semantics for file access.
 *
 * Allowed values for mpi_io_mode are:
 *   0: no MPI-IO,
 *   1: MPI-IO with explicit offsets,
 *   2: MPI-IO with individual file pointers
 *
 * Invalid values (for example an MPI-IO mode with no MPI or MPI-IO
 * support) are silently ignored.
 *
 * parameters:
 *   mpi_io_mode <-- mode for default semantics
 *----------------------------------------------------------------------------*/

void
cs_io_set_defaults(int  mpi_io_mode)
{
#if defined(HAVE_MPI)

  if (mpi_io_mode == 0)
    cs_glob_io_hints = FVM_FILE_NO_MPI_IO;
  else if (mpi_io_mode == 1)
    cs_glob_io_hints = FVM_FILE_EXPLICIT_OFFSETS;
  else if (mpi_io_mode == 2)
    cs_glob_io_hints = FVM_FILE_INDIVIDUAL_POINTERS;

#endif

  fvm_file_set_default_semantics(cs_glob_io_hints);
}

/*----------------------------------------------------------------------------
 * Dump a kernel IO file handle's metadata.
 *
 * parameters:
 *   cs_io  <-- kernel IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_dump(const cs_io_t  *cs_io)
{
  assert(cs_io != NULL);

  bft_printf(_("\n\n file contents:\n\n"));

  if (cs_io->f != NULL)
    bft_printf(_("  file: %s\n"), fvm_file_get_name(cs_io->f));

  bft_printf(_("  contents: \"%s\"\n"), cs_io->contents);
  if (cs_io->mode == CS_IO_MODE_READ)
    bft_printf(_("  mode: CS_IO_MODE_READ\n"), cs_io->contents);
  else if (cs_io->mode == CS_IO_MODE_WRITE)
    bft_printf(_("  mode: CS_IO_MODE_WRITE\n"), cs_io->contents);

  bft_printf(_("  default header size: %lu\n"
               "  header alignment:    %lu\n"
               "  body alignment:      %lu\n"
               "  verbosity level:     %ld\n\n"),
             cs_io->header_size, cs_io->header_align, cs_io->body_align,
             cs_io->echo);

  if (cs_io->index != NULL)
    _dump_index(cs_io->index);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
