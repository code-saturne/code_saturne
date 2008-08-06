/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT and FVM library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include <fvm_file.h>

#if defined(_CS_HAVE_MPI)
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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local types and structures
 *============================================================================*/

struct _cs_io_t {

  char            *name;         /* File name */

  fvm_file_t      *f;            /* Pointer to associated file */

  cs_io_mode_t     mode;         /* Communication mode */
  cs_int_t         echo;         /* Data echo level */

};

/*============================================================================
 * Constants and Macros
 *============================================================================*/

#define CS_IO_MPI_TAG     'C'+'S'+'_'+'I'+'O'

/*============================================================================
 * Static global variables
 *============================================================================*/

static char  _cs_io_type_name_char[] = "c ";  /* Character string */
static char  _cs_io_type_name_i4[] =   "i4";  /* Signed 32 bit integer */
static char  _cs_io_type_name_i8[] =   "i8";  /* Signed 64 bit integer */
static char  _cs_io_type_name_u4[] =   "u4";  /* Unsigned 32 bit integer */
static char  _cs_io_type_name_u8[] =   "u8";  /* Unsigned 64 bit integer */
static char  _cs_io_type_name_r4[] =   "r4";  /* Single precision real */
static char  _cs_io_type_name_r8[] =   "r8";  /* Double precsision real */

/* Global pointer on preprocessor data file handle */
cs_io_t  *cs_glob_pp_io = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Open the interface file descriptor and initialize the file by writing
 * or reading a "magic string" used to check the file format.
 *
 * parameters:
 *   pp_io        <-> preprocessor IO structure
 *   name         --> file name
 *   magic_string --> magic string associated with file type
 *----------------------------------------------------------------------------*/

static void
_cs_io_file_open(cs_io_t     *const pp_io,
                 const char  *const name,
                 const char  *const magic_string)
{
  fvm_file_mode_t pp_io_mode;

  int hints = 0;

  /* Prepare file open */

  switch(pp_io->mode) {

  case CS_IO_MODE_READ:
    pp_io_mode = FVM_FILE_MODE_READ;
    break;

  case CS_IO_MODE_WRITE:
    pp_io_mode = FVM_FILE_MODE_WRITE;
    break;

  default:
    assert(   pp_io->mode == CS_IO_MODE_READ
           || pp_io->mode == CS_IO_MODE_WRITE);
  }

  /* Create interface file descriptor */

  pp_io->f = fvm_file_open(name, pp_io_mode, FVM_FILE_TYPE_BINARY, hints);
  fvm_file_set_big_endian(pp_io->f);

  /* Write or read a magic string */
  /*------------------------------*/

  if (pp_io->mode == CS_IO_MODE_READ) {

    char    *magic_str_read;
    size_t   magic_str_size = strlen(magic_string);

    BFT_MALLOC(magic_str_read, magic_str_size + 1, char);

    fvm_file_read_global(magic_str_read, 1, strlen(magic_string), pp_io->f);

    magic_str_read[magic_str_size] = '\0';

    /* If the magic string does not correspond, we have an error */

    if (strcmp(magic_str_read, magic_string) != 0) {

      bft_error(__FILE__, __LINE__, 0,
                _("Erreur à la lecture du fichier de pré traitement : "
                  "\"%s\".\n"
                  "Le format de l'interface n'est pas à la bonne version.\n"
                  "La chaîne magique repère la version du format "
                  "d'interface :\n"
                  "chaîne magique lue      : \"%s\"\n"
                  "chaîne magique actuelle : \"%s\"\n"),
                pp_io->name, magic_str_read, magic_string);

    }

    BFT_FREE(magic_str_read);

  }
  else if (pp_io->mode == CS_IO_MODE_WRITE) {

      bft_error(__FILE__, __LINE__, 0,
                _("Erreur à l'écriture du fichier de pré traitement : "
                  "\"%s\".\n"
                  "Cette fonctionnalité n'est pas encore implémentée."),
                pp_io->name);

  }

}

/*----------------------------------------------------------------------------
 * Close the interface file.
 *
 * parameters:
 *   pp_io <-> preprocessor IO structure
 *----------------------------------------------------------------------------*/

static void
_cs_io_file_close(cs_io_t  *pp_io)
{
  pp_io->f = fvm_file_free(pp_io->f);
}

/*----------------------------------------------------------------------------
 * Echo pending section read or write
 *
 * parameters:
 *   pp_io --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

static void
_cs_io_echo_pre(const cs_io_t  *pp_io)
{
  assert(pp_io != NULL);

  switch(pp_io->mode) {

  case CS_IO_MODE_READ:
    bft_printf(_("\nSection lue sur \"%s\" :\n"), pp_io->name);
    break;

  case CS_IO_MODE_WRITE:
    bft_printf(_("\nSection écrite sur \"%s\" :\n"), pp_io->name);
    break;

  default:
    assert(   pp_io->mode == CS_IO_MODE_READ
           || pp_io->mode == CS_IO_MODE_WRITE);
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Echo a message header
 *
 * parameters:
 *   sec_name  --> section name
 *   n_elts    --> number of elements
 *   elt_type  --> expected (final) element type
 *   type_read --> read element type (convertible to expected type)
 *----------------------------------------------------------------------------*/

static void
_cs_io_echo_header(const char  *sec_name,
                   fvm_gnum_t   n_elts,
                   cs_type_t    type_read)
{
  char sec_name_echo[CS_IO_NAME_LEN + 1];

  /* Instructions */

  strncpy(sec_name_echo, sec_name,  CS_IO_NAME_LEN);
  sec_name_echo[CS_IO_NAME_LEN] = '\0';

  bft_printf(_("    nom de la rubrique    : \"%s\"\n"
               "    nombre d'éléments     : %lu\n"),
             sec_name_echo, (unsigned long)n_elts);

  if (n_elts > 0) {

    char *type_name;

    switch(type_read) {
    case FVM_DATATYPE_NULL:
      type_name = _cs_io_type_name_char;
      break;
    case FVM_INT32:
      type_name = _cs_io_type_name_i4;
      break;
    case FVM_INT64:
      type_name = _cs_io_type_name_i8;
      break;
    case FVM_UINT32:
      type_name = _cs_io_type_name_u4;
      break;
    case FVM_UINT64:
      type_name = _cs_io_type_name_u8;
      break;
    case FVM_FLOAT:
      type_name = _cs_io_type_name_r4;
      break;
    case FVM_DOUBLE:
      type_name = _cs_io_type_name_r8;
      break;
    default:
      assert(0);
    }

    bft_printf(_("    nom du type d'élément : \"%s\"\n"), type_name);

  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Partial echo of a section's contents.
 *
 * FVM datatypes must have been converted to the corresponding
 * Code_Saturne compatible datatype before calling this function:
 *   FVM_DATATYPE_NULL       -> char
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
 *   n_elts           -> number of elements
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *   elt_type         --> element type
 *   elts             --> element buffer
 *----------------------------------------------------------------------------*/

static void
_cs_io_echo_data(size_t       echo,
                 size_t       n_elts,
                 fvm_gnum_t   global_num_start,
                 fvm_gnum_t   global_num_end,
                 cs_type_t    elt_type,
                 const void  *elts)
{
  fvm_gnum_t  i;
  fvm_gnum_t  num_shift = 1;
  size_t  _n_elts = n_elts;
  size_t  echo_start = 0;
  size_t  echo_end = 0;
  const char *_loc_glob[] = {N_(" (locaux)"), ""};
  const char *loc_glob = _loc_glob[1];

  /* Instructions */

  if (n_elts == 0) return;

  if (cs_glob_base_nbr == 1)
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
    bft_printf(_("    %d premiers et derniers éléments%s :\n"),
               echo, loc_glob);
  }
  else {
    echo_end = _n_elts;
    bft_printf(_("    éléments%s :\n"), loc_glob);
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

    case CS_TYPE_char:
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

      assert(   elt_type == CS_TYPE_cs_int_t
             || elt_type == CS_TYPE_cs_real_t
             || elt_type == CS_TYPE_char);

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
                    size_t           n_elts,
                    fvm_datatype_t   buffer_type,
                    fvm_datatype_t   dest_type)
{
  size_t ii;
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
 * Read a message body.
 *
 * If global_num_start and global_num_end are > 0, a different block is
 * assigned to each processor. Otherwise, the full data is replicated
 * for each processor.
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
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

static void *
_cs_io_read_body(const cs_io_sec_header_t  *header,
                 fvm_gnum_t                 global_num_start,
                 fvm_gnum_t                 global_num_end,
                 void                      *elts,
                 cs_io_t                   *pp_io)
{
  size_t      type_size = 0;
  fvm_gnum_t  n_elts = header->n_elts;
  cs_bool_t   convert_type = CS_FALSE;
  void       *_elts = NULL;
  void       *_buf = NULL;

  assert(pp_io  != NULL);

  assert(global_num_end <= header->n_elts + 1);

  /* Choose global or block mode */

  if (global_num_start > 0 && global_num_end > 0) {
    assert(global_num_end >= global_num_start);
    n_elts = global_num_end - global_num_start;
  }

  /* Datatype size given by FVM datatype, except for character type */

  type_size = fvm_datatype_size[header->type_read];
  if (type_size == 0)
    type_size = 1;

  /* Assign or allocate */

  _elts = elts;

  if (_elts == NULL && n_elts != 0) {
    if (header->elt_type == CS_TYPE_char)
      BFT_MALLOC(_elts, n_elts + 1, char);
    else
      BFT_MALLOC(_elts, n_elts*type_size, char);
  }

  /* Element values */

  if (n_elts != 0 && header->elt_type != header->type_read)
    convert_type = CS_TRUE;

  if (convert_type == true
      && (   fvm_datatype_size[header->type_read]
          != fvm_datatype_size[header->elt_type]))
    BFT_MALLOC(_buf, n_elts*type_size, char);
  else
    _buf = _elts;

  if (global_num_start > 0 && global_num_end > 0)
    fvm_file_read_block(_buf,
                        type_size,
                        global_num_start,
                        global_num_end,
                        pp_io->f);
  else
    fvm_file_read_global(_buf,
                         type_size,
                         n_elts,
                         pp_io->f);

  /* Convert data if necessary */

  if (convert_type == true) {
    _cs_io_convert_read(_buf,
                        _elts,
                        n_elts,
                        header->type_read,
                        header->elt_type);
    if (_buf != _elts)
      BFT_FREE(_buf);
  }

  /* Add null character at end of string to ensure C-type string */

  if (n_elts != 0 && header->elt_type == CS_TYPE_char)
    ((char *)_elts)[header->n_elts] = '\0';

  /* Optional echo */

  if (header->n_elts != 0 && pp_io->echo > 0)
    _cs_io_echo_data(pp_io->echo,
                     n_elts,
                     global_num_start,
                     global_num_end,
                     header->elt_type,
                     _elts);

  /* Return read values */

  return _elts;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a preprocessor output file structure.
 *
 * parameters:
 *   name         --> file name
 *   magic_string --> magic string associated with file type
 *   mode         --> read or write
 *   echo         --> echo on main output (< 0 if none, header if 0,
 *                    n first and last elements if n > 0)
 *
 * returns:
 *   pointer to preprocessor IO structure
 *----------------------------------------------------------------------------*/

cs_io_t *
cs_io_initialize(const char    *file_name,
                 const char    *magic_string,
                 cs_io_mode_t   mode,
                 size_t         echo)
{
  cs_io_t  *pp_io = NULL;

  BFT_MALLOC(pp_io, 1, cs_io_t);

  /* Set structure fields */

  BFT_MALLOC(pp_io->name, strlen(file_name) + 1, char);

  strcpy(pp_io->name, file_name);

  pp_io->mode = mode;
  pp_io->echo = echo;

  pp_io->f  = NULL;

  /* Info on interface creation */

  bft_printf(_("\n Lecture du pré traitement :  %s"), file_name);
  bft_printf_flush();

  /* Create interface file descriptor */

  _cs_io_file_open(pp_io, pp_io->name, magic_string);

  return pp_io;
}

/*----------------------------------------------------------------------------
 * Free a preprocessor output file structure, closing the associated file.
 *
 * parameters:
 *   pp_io <-> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_finalize(cs_io_t **pp_io)
{
  cs_io_t *_pp_io = *pp_io;

  /* Info on closing of interface file */

  bft_printf(_("\n Fin de la lecture :  %s\n"), _pp_io->name);
  bft_printf_flush();

  _cs_io_file_close(_pp_io);

  BFT_FREE(_pp_io->name);
  BFT_FREE(*pp_io);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a preprocessor IO structure's name.
 *
 * parameters:
 *   pp_io --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

const char *
cs_io_get_name(const cs_io_t  *pp_io)
{
  assert(pp_io != NULL);

  return(pp_io->name);
}

/*----------------------------------------------------------------------------
 * Return a preprocessor IO structure's echo (verbosity) level.
 *
 * parameters:
 *   pp_io --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

size_t
cs_io_get_echo(const cs_io_t  *pp_io)
{
  assert(pp_io != NULL);

  return (size_t)(pp_io->echo);
}

/*----------------------------------------------------------------------------
 * Read a message header.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_read_header(cs_io_sec_header_t  *header,
                  cs_io_t             *pp_io)
{
  assert(pp_io != NULL);

  header->n_elts = 0;

  if (pp_io->echo >= 0)
    _cs_io_echo_pre(pp_io);


  /* Read from file */
  /*----------------*/

  /* section type name */

  fvm_file_read_global(header->sec_name, 1, CS_IO_NAME_LEN, pp_io->f);

  header->sec_name[CS_IO_NAME_LEN] = '\0';

  /* No true end-of-file detection with current fvm_file API,
     so we guess that we are at the end when a header is empty. */

  if (strlen(header->sec_name) == 0) {
    strcpy(header->sec_name, "EOF");
    return;
  }

  /* number of elements */

  if (sizeof(unsigned long) == 8) {
    unsigned long rec;
    fvm_file_read_global(&rec, 8, 1, pp_io->f);
    header->n_elts = rec;
  }
  else if (sizeof(unsigned long long) == 8) {
    unsigned long long rec;
    fvm_file_read_global(&rec, 8, 1, pp_io->f);
    header->n_elts = rec;
  }

  assert(sizeof(unsigned long) == 8 || sizeof(unsigned long long) == 8);

  if (header->n_elts != 0) {

    char *elt_type_name = header->type_read_name;

    fvm_file_read_global(header->type_read_name, 1, 2, pp_io->f);
    header->type_read_name[2] = '\0';

    if (   strcmp(elt_type_name, _cs_io_type_name_i4) == 0
        || strcmp(elt_type_name, "i ") == 0)
      header->type_read = FVM_INT32;

    else if (strcmp(elt_type_name, _cs_io_type_name_i8) == 0)
      header->type_read = FVM_INT64;

    else if (strcmp(elt_type_name, _cs_io_type_name_u4) == 0)
      header->type_read = FVM_UINT32;

    else if (strcmp(elt_type_name, _cs_io_type_name_u8) == 0)
      header->type_read = FVM_UINT64;

    else if (strcmp(elt_type_name, _cs_io_type_name_r4) == 0)
      header->type_read = FVM_FLOAT;

    else if (strcmp(elt_type_name, _cs_io_type_name_r8) == 0)
      header->type_read = FVM_DOUBLE;

    else if (strcmp(elt_type_name, _cs_io_type_name_char) == 0)
      header->type_read = FVM_DATATYPE_NULL;

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Erreur à la lecture du fichier de pré traitement : "
                  "\"%s\".\n"
                  "Le type de données \"%s\" n'est pas reconnu."),
                pp_io->name, elt_type_name);

    if (header->type_read == FVM_INT32 || header->type_read == FVM_INT64) {
      assert(sizeof(fvm_lnum_t) == 4 || sizeof(fvm_lnum_t) == 8);
      if (sizeof(fvm_lnum_t) == 4)
        header->elt_type = FVM_INT32;
      else
        header->elt_type = FVM_INT64;
    }

    else if (   header->type_read == FVM_UINT32
             || header->type_read == FVM_UINT64) {
      assert(sizeof(fvm_gnum_t) == 4 || sizeof(fvm_gnum_t) == 8);
      if (sizeof(fvm_gnum_t) == 4)
        header->elt_type = FVM_UINT32;
      else
        header->elt_type = FVM_UINT64;
    }

    else if (   header->type_read == FVM_FLOAT
             || header->type_read == FVM_DOUBLE) {
      if (sizeof(cs_real_t) == 4)
        header->elt_type = FVM_FLOAT;
      else
        header->elt_type = FVM_DOUBLE;
    }

    else if (header->type_read == FVM_DATATYPE_NULL)
      header->elt_type = FVM_DATATYPE_NULL;

  }

  else
    header->type_read_name[0] = '\0';

  /* Affichage eventuel */

  if (pp_io->echo >= 0)
    _cs_io_echo_header(header->sec_name,
                       header->n_elts,
                       header->elt_type);
}

/*----------------------------------------------------------------------------
 * Set a message's final data type to fvm_lnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_lnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *pp_io)
{
  assert(header != NULL);

  if (   header->type_read != FVM_INT32
      && header->type_read != FVM_INT64
      && header->type_read != FVM_UINT32
      && header->type_read != FVM_UINT64)
    bft_error(__FILE__, __LINE__, 0,
              _("Erreur à la lecture du fichier de pré traitement : "
                "\"%s\".\n"
                "Le type attendu pour la section : "
                "\"%s\" est un entier signé.\n"
                "et n'est pas convertible à partir du type lu \"%s\"."),
              pp_io->name, header->type_read_name);

  assert(sizeof(fvm_lnum_t) == 4 || sizeof(fvm_lnum_t) == 8);

  if (sizeof(fvm_lnum_t) == 4)
    header->elt_type = FVM_INT32;
  else
    header->elt_type = FVM_INT64;
}

/*----------------------------------------------------------------------------
 * Set a message's final data type to fvm_gnum_t.
 *
 * It the datatype is not compatible, throw an error.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_set_fvm_gnum(cs_io_sec_header_t  *header,
                   const cs_io_t       *pp_io)
{
  assert(header != NULL);

  if (   header->type_read != FVM_INT32
      && header->type_read != FVM_INT64
      && header->type_read != FVM_UINT32
      && header->type_read != FVM_UINT64)
    bft_error(__FILE__, __LINE__, 0,
              _("Erreur à la lecture du fichier de pré traitement : "
                "\"%s\".\n"
                "Le type attendu pour la section : "
                "\"%s\" est un entier non signé.\n"
                "et n'est pas convertible à partir du type lu \"%s\"."),
              pp_io->name, header->type_read_name);

  assert(sizeof(fvm_gnum_t) == 4 || sizeof(fvm_gnum_t) == 8);

  if (sizeof(fvm_gnum_t) == 4)
    header->elt_type = FVM_UINT32;
  else
    header->elt_type = FVM_UINT64;
}

/*----------------------------------------------------------------------------
 * Check that a message's final data type corresponds to cs_real_t.
 *
 * parameters:
 *   header <-- header structure
 *   pp_io  --> preprocessor IO structure
 *----------------------------------------------------------------------------*/

void
cs_io_assert_cs_real(const cs_io_sec_header_t  *header,
                     const cs_io_t             *pp_io)
{
  assert(header != NULL);

  if (   header->elt_type != FVM_FLOAT
      && header->elt_type != FVM_DOUBLE)
    bft_error(__FILE__, __LINE__, 0,
              _("Erreur à la lecture du fichier de pré traitement : "
                "\"%s\".\n"
                "Le type attendu pour la section : "
                "\"%s\".\n"
                "est \"r4\" ou \"r8\" (réel), et non \"%s\"."),
              pp_io->name, header->type_read_name);
}

/*----------------------------------------------------------------------------
 * Read a message body and replicate it to all processors.
 *
 * If global_num_start and global_num_end are > 0, a different block is
 * assigned to each processor. Otherwise, the full data is replicated
 * for each processor.
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
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_global(const cs_io_sec_header_t  *header,
                  void                      *elts,
                  cs_io_t                   *pp_io)
{
  return _cs_io_read_body(header, 0, 0, elts, pp_io);
}

/*----------------------------------------------------------------------------
 * Read a message body, assigning a different block to each processor.
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
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_block(const cs_io_sec_header_t  *header,
                 fvm_gnum_t                 global_num_start,
                 fvm_gnum_t                 global_num_end,
                 void                      *elts,
                 cs_io_t                   *pp_io)
{
  assert(global_num_start > 0);
  assert(global_num_end >= global_num_start);

  return _cs_io_read_body(header,
                          global_num_start,
                          global_num_end,
                          elts,
                          pp_io);
}

/*----------------------------------------------------------------------------
 * Read a message body, assigning a different block to each processor,
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
 *   pp_io            --> preprocessor IO structure
 *
 * returns:
 *   elts if non NULL, or pointer to allocated array otherwise
 *----------------------------------------------------------------------------*/

void *
cs_io_read_index_block(cs_io_sec_header_t  *header,
                       fvm_gnum_t           global_num_start,
                       fvm_gnum_t           global_num_end,
                       fvm_gnum_t          *elts,
                       cs_io_t             *pp_io)
{
  fvm_gnum_t _global_num_start = global_num_end;
  fvm_gnum_t _global_num_end = global_num_end;
  fvm_gnum_t *retval = NULL;
  cs_bool_t last_data_rank = false;
  cs_bool_t past_last_data_rank = false;

  assert(global_num_start > 0);
  assert(global_num_end >= global_num_start);

  /* Check type */

  cs_io_set_fvm_gnum(header, pp_io);

  /* Increase _global_num_end by 1 for the last rank */

  if (header->n_elts == global_num_end) {

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
                            global_num_start,
                            _global_num_end,
                            elts,
                            pp_io);

  /* Exchange past-the-end values */

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_nbr > 1) {

    int needs_safe_algo_loc = 0;
    int needs_safe_algo = 0;
    int rank = cs_glob_base_rang;
    int send_rank = rank - 1;
    int recv_rank = rank + 1;
    MPI_Comm comm = cs_glob_base_mpi_comm;
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

  if (   header->n_elts != 0 && header->n_elts != global_num_end
      && pp_io->echo > 0)
    bft_printf(_("    premier élement rang suivant :\n"
                 "    %10lu : %12d\n"),
               (unsigned long)(global_num_end),
               (unsigned long)retval[global_num_end - global_num_start]);

#endif /* defined(_CS_HAVE_MPI) */

  return retval;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
