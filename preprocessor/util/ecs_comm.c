/*============================================================================
 * Base functions for writing Kernel I/O files.
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

#include "ecs_def.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"

/*----------------------------------------------------------------------------
 *  Headers for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_comm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an output writer
 *----------------------------------------------------------------------------*/

struct _ecs_comm_t {

  char           *name;             /* Writer name */

  ecs_file_t     *f;                /* Associated file structure pointer */

  size_t          header_size;      /* Header default size */
  size_t          header_align;     /* Header alignment */
  size_t          body_align;       /* Body alignment */

  size_t          cur_pos;          /* Current position in file */

};

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write padding zeroes to ensure following alignement if necessary
 *----------------------------------------------------------------------------*/

static void
_write_pad(ecs_comm_t  *comm,
           size_t       alignment)
{
  size_t pad_size = (alignment - (comm->cur_pos % alignment)) % alignment;

  if (pad_size > 0) {

    char padding[128] = "";
    size_t rem_size = pad_size;

    memset(padding, 0, sizeof(padding));

    while (rem_size > 0) {
      size_t write_size = ECS_MIN(rem_size, sizeof(padding));
      ecs_file_write(padding, 1, write_size, comm->f);
      rem_size -= write_size;
    }

    comm->cur_pos += pad_size;
  }
}

/*----------------------------------------------------------------------------
 * Write a record to the interface file
 *----------------------------------------------------------------------------*/

static void
_write_rec(ecs_comm_t  *comm,
           const void  *rec,
           size_t       n_elts,
           ecs_type_t   datatype)
{
  size_t t_size = 0;

  assert(comm != NULL);
  assert(rec  != NULL);

  /* Determine size associated with type */

  switch (datatype) {
  case ECS_TYPE_ecs_int_t:
    t_size = sizeof(ecs_int_t);
    break;
  case ECS_TYPE_ecs_coord_t:
    t_size = sizeof(ecs_coord_t);
    break;
  case ECS_TYPE_ecs_size_t:
    t_size = sizeof(ecs_size_t);
    break;
  case ECS_TYPE_char:
    t_size = sizeof(char);
    break;
  case ECS_TYPE_size_t:
    t_size = 8;
    break;
  default:
    assert(   datatype == ECS_TYPE_ecs_int_t
           || datatype == ECS_TYPE_ecs_coord_t
           || datatype == ECS_TYPE_ecs_size_t
           || datatype == ECS_TYPE_size_t
           || datatype == ECS_TYPE_char);
  }

  /* write record to file */
  /*----------------------*/

  if (datatype != ECS_TYPE_size_t || sizeof(size_t) == 8) {
    ecs_file_write(rec, t_size, n_elts, comm->f);
    comm->cur_pos += n_elts*t_size;
  }

  else { /* if (datatype == ECS_TYPE_size_t && sizeof(size_t) != 8) */

    size_t i;

    assert(sizeof(unsigned long long) == 8);

    for (i = 0; i < n_elts; i++) {
      const unsigned char *rp = rec;
      unsigned long long _rec = *((const size_t *)(rp + i*sizeof(size_t)));
      assert(sizeof(long long) == 8);
      ecs_file_write(&_rec, 8, 1, comm->f);
      comm->cur_pos += 8;
    }

  }
}

/*----------------------------------------------------------------------------
 * Open output file and write base file description headers
 *----------------------------------------------------------------------------*/

static void
_file_open(ecs_comm_t  *comm,
           const char  *file_name)
{
  char header[64] = "";
  size_t sizes[3] = {comm->header_size,
                     comm->header_align,
                     comm->body_align};

  memset(header, 0, 64);

  /* BE stands for big-endian, allowing for future
     native file mode generation, using BE/LE */

  strncpy(header, "Code_Saturne I/O, BE, R0", 63);

  /* Open file */

  comm->f = ecs_file_open(file_name,
                          ECS_FILE_MODE_WRITE,
                          ECS_FILE_TYPE_BINARY);

  ecs_file_set_big_endian(comm->f);

  /* Write header and comment information */

  _write_rec(comm, header, 64, ECS_TYPE_char);

  memset(header, 0, 64);
  strncpy(header, "Face-based mesh definition, R0", 63);

  _write_rec(comm, header, 64, ECS_TYPE_char);

  _write_rec(comm, sizes, 3, ECS_TYPE_size_t);
}

/*----------------------------------------------------------------------------
 * Echo section output data
 *----------------------------------------------------------------------------*/

static void
_echo_header(const char   *name,
             size_t        n_values,
             const char   *datatype_name)
{
  char name_pad[33] = "";
  char type_pad[3] = "";

  if (strlen(name) < 32) {
    memset(name_pad, ' ', sizeof(name_pad));
    name_pad[32 - strlen(name)] = '\0';
  }
  if (strlen(datatype_name) < 2) {
    type_pad[0] = ' '; type_pad[1] = ' ';
    type_pad[2 - strlen(datatype_name)] = '\0';
  }

  if (n_values > 0)
    printf(_("  Wrote: \"%s\"%s; Type: \"%s\"%s; Size: %lu\n"),
           name, name_pad, datatype_name, type_pad,
           (unsigned long)n_values);

  else
    printf(_("  Wrote: \"%s\"\n"),
           name);

  fflush(stdout);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

ecs_comm_t *
ecs_comm_initialize(const char  *file_name)
{
  ecs_comm_t  * comm = NULL;

  /* Create structure */

  ECS_MALLOC(comm, 1, ecs_comm_t);

  ECS_MALLOC(comm->name, strlen(file_name) + 1, char);

  strcpy(comm->name, file_name);

  /* Initialize other fields */

  comm->header_size = 96;
  comm->header_align = 64;
  comm->body_align = 64;
  comm->cur_pos = 0;

  comm->f  = NULL;

  /* Jump a line in log file */

  printf("\n");

  /* Info on interface creation */

  printf(_("  Opening file: %s\n"),
         comm->name);

  printf("\n");

  fflush(stdout);

  /* Create file descriptor */

  _file_open(comm, file_name);

  return comm;
}

/*----------------------------------------------------------------------------
 * Close writer.
 *
 * arguments:
 *   comm <-- pointer to writer structure pointer
 *----------------------------------------------------------------------------*/

void
ecs_comm_finalize(ecs_comm_t **comm)
{
  ecs_comm_t *_comm = *comm;

  printf("\n");

  printf(_("  Closing file: %s\n"),
         _comm->name);

  if (_comm->f != NULL)
    _comm->f = ecs_file_free(_comm->f);

  ECS_FREE(_comm->name);

  ECS_FREE(_comm);

  *comm = _comm;
}

/*----------------------------------------------------------------------------
 * Write a section to the Kernel I/O file.
 *
 * Grid locations and possibly indexes may be assigned to a section by
 * specifying a location id; the first time a given location id appears in
 * the file is considered a declaration. In the same manner, an index id
 * may be specified. Values of zero indicate no location or index base
 * is used. It is up to the calling code to ensure that total number of
 * values, location size, and number of values per location are consistent,
 * as this my be important for code reading the file.
 *
 * arguments:
 *   name              <-- section name
 *   location_id       <-- id of associated location
 *   index_id          <-- id of associated index
 *   n_location_values <-- number of values per location
 *   embed             <-- embed values in header
 *   values            <-- values to write
 *   value_type        <-- type of value to write
 *   comm              <-- Kernel I/O file output structure
 *----------------------------------------------------------------------------*/

void
ecs_comm_write_section(const char  *name,
                       size_t       n_values,
                       size_t       location_id,
                       size_t       index_id,
                       size_t       n_location_values,
                       bool         embed,
                       const void  *values,
                       ecs_type_t   value_type,
                       ecs_comm_t  *comm)
{
  /* Section is defined by:

     section size (in bytes): 8 bytes
     n_elts, location_id, index_id, n_location_values: 4*8 bytes
     name_size (including padding to 8 byte alignment): 8 bytes
     datatype_name: 8 bytes (bit 6 = '\0', bit 7 = embed flag)
     section name: strlen(name) + 1 + padding to 8 bytes
     optional embedded data: n_values*value_type_size
  */

  size_t header_sizes[6] = {56,       /* 6*8 +8 */
                            n_values,
                            location_id,
                            index_id,
                            n_location_values,
                            0};
  size_t name_size = 0;
  size_t name_pad_size = 0;
  size_t data_size = 0;
  char   datatype_name[8];
  char   name_pad[8];

  assert(comm != NULL);
  assert(name != NULL);

  /* Initialization */

  memset(datatype_name, 0, sizeof(datatype_name));
  memset(name_pad, 0, sizeof(name_pad));

  /* Section name */

  name_size = strlen(name);
  name_pad_size = 8-(name_size%8); /* At least 1 NULL
                                      character with this rule */

  header_sizes[0] += (name_size + name_pad_size);
  header_sizes[5] = (name_size + name_pad_size);

  /* Value type name */

  if (n_values > 0) {

    switch(value_type) {

    case ECS_TYPE_ecs_int_t:

      switch(sizeof(ecs_int_t)) {
      case 4:
        strcpy(datatype_name, "i4");
        break;
      case 8:
        strcpy(datatype_name, "i8");
        break;
      default:
        assert(sizeof(ecs_int_t) == 4 || sizeof(ecs_int_t) == 8);
      }
      break;

    case ECS_TYPE_ecs_coord_t:

      switch(sizeof(ecs_coord_t)) {
      case 4:
        strcpy(datatype_name, "r4");
        break;
      case 8:
        strcpy(datatype_name, "r8");
        break;
      default:
        assert(sizeof(ecs_coord_t) == 4 || sizeof(ecs_coord_t) == 8);
      }
      break;

    case ECS_TYPE_ecs_size_t:

      switch(sizeof(ecs_size_t)) {
      case 4:
        strcpy(datatype_name, "u4");
        break;
      case 8:
        strcpy(datatype_name, "u8");
        break;
      default:
        assert(sizeof(ecs_size_t) == 4 || sizeof(ecs_size_t) == 8);
      }
      break;

    case ECS_TYPE_size_t:

      strcpy(datatype_name, "u8");
      break;

    case ECS_TYPE_char:

      strcpy(datatype_name, "c ");
      break;

    default:

      assert(   value_type == ECS_TYPE_ecs_int_t
             || value_type == ECS_TYPE_ecs_coord_t
             || value_type == ECS_TYPE_ecs_size_t
             || value_type == ECS_TYPE_size_t
             || value_type == ECS_TYPE_char);


    } /* End of switch on value_type */

  }

  if (embed == true) {

    datatype_name[7] = 'e';

    if (datatype_name[1] == '4')
      data_size = 4*n_values;
    else if (datatype_name[1] == '8')
      data_size = 8*n_values;
    else
      data_size = n_values;

    header_sizes[0] += data_size;

  }

  /* Only output data if file exists */

  if (comm->f != NULL) {

    /* Align if necessary */

    _write_pad(comm, comm->header_align);

    /* Write sizes */

    _write_rec(comm, header_sizes, 6, ECS_TYPE_size_t);

    /* Type information */

    _write_rec(comm, datatype_name, 8, ECS_TYPE_char);

    /* Section name */

    _write_rec(comm, name, name_size, ECS_TYPE_char);
    _write_rec(comm, name_pad, name_pad_size, ECS_TYPE_char);

    /* Embedded values */

    if (n_values > 0 && embed == true)
      _write_rec(comm, values, n_values, value_type);

    /* Ensure header is at least the size of the basic header block size,
       as reads may read at least this much data */

    if (header_sizes[0] < comm->header_size) {
      char   header_pad[64];
      size_t header_pad_size = comm->header_size - header_sizes[0];
      memset(header_pad, 0, sizeof(header_pad));
      while (header_pad_size > 0) {
        size_t write_size = ECS_MIN(header_pad_size, sizeof(header_pad));
        ecs_file_write(header_pad, 1, write_size, comm->f);
        header_pad_size -= write_size;
        comm->cur_pos += write_size;
      }
    }

    /* If values are not embedded, handle separate value block */

    if (n_values > 0 && embed == false) {

      _write_pad(comm, comm->body_align);
      _write_rec(comm, values, n_values, value_type);

    }

  }

  /* Print message */

  _echo_header(name, n_values, datatype_name);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

