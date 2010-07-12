/*============================================================================
 *  Convert a NUMECA Hex mesh to MED file format.
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2009 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

#include "cs_config.h"

/*----------------------------------------------------------------------------*
 *  C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------*
 *  ECS headers
 *----------------------------------------------------------------------------*/

#include "ecs_def.h"
#include "ecs_file.h"
#include "ecs_mem.h"

/*----------------------------------------------------------------------------*
 *  MED library headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include <med.h>

#ifdef __cplusplus
}
#endif

/*============================================================================
 * Definitions that may not always be provided directly by the system
 *============================================================================*/

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if HAVE_STDDEF_H
# include <stddef.h>
#else
# include <stdlib.h>
#endif

/*
 * Usually stdint.h is included by inttypes.h, but only inttypes.h exists
 * on certain systems, such as Tru64 Unix
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* C99 _Bool type */

#if HAVE_STDBOOL_H
# include <stdbool.h>
#else
# if !HAVE__BOOL
#  ifdef __cplusplus
typedef bool _Bool;
#  else
typedef unsigned char _Bool;
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_tru_false_are_defined 1
#endif

/* int32_t type */

#if !defined(HAVE_INT32_T)
# if (SIZEOF_INT == 4)
typedef int int32_t;
# elif (SIZEOF_SHORT == 4)
typedef short int32_t;
# else
#  error
# endif
#endif

/*============================================================================
 *                            Macro definitions
 *============================================================================*/

/* To read 80 characters per line, adding `\n', `\0', plus safety margin */

#define MAX_LINE_LENGTH  84

#define CROSS_PRODUCT(u, v, w) (\
  u[0] = v[1]*w[2] - w[1]*v[2], \
  u[1] = w[0]*v[2] - v[0]*w[2], \
  u[2] = v[0]*w[1] - w[0]*v[1] )

#define DOT_PRODUCT(v, w) (v[0]*w[0] + v[1]*w[1] + v[2]*w[2])

/*============================================================================
 *                              Private functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Read an integer (format " %d") from a text file with potentially very
 *  long lines. In this case, we read to a buffer, which we slide after
 *  extracting the integer.
 *----------------------------------------------------------------------------*/

static int
_scan_int_from_line(ecs_file_t  *f,                     /* --> File */
                    int         *line_num,              /* <-> Line number */
                    char        *line,                  /* <-> Buffer */
                    int          len)                   /* <-> String length  */
{
  int  p1, p2, res;

  for (p1 = 0;
       (   line[p1] != '0' && line[p1] != '1'
        && line[p1] != '2' && line[p1] != '3'
        && line[p1] != '4' && line[p1] != '5'
        && line[p1] != '6' && line[p1] != '7'
        && line[p1] != '8' && line[p1] != '9');
       p1++) {
    if (line[p1] == '\0') {
      (*line_num)--;
      ecs_file_gets(line, len, f, line_num);
      p1 = 0;
    }
  }

  for (p2 = p1 + 1;
       line[p2] != '\0'  && line[p2] != '\n' && line[p2] != ' ';
       p2++);

  /* If the buffer seems cut, we read the next part */
  if (line[p2] == '\0') {
    (*line_num)--;
    ecs_file_gets(line+p2, len-p2, f, line_num);
  }

  while (line[p2] != '\0' && line[p2] != '\n' && line[p2] != ' ')
    p2++;

  assert (line[p2] != '\0');

  if (sscanf(line+p1, "%d", &res) != 1)
    ecs_error(__FILE__, __LINE__, errno,
              "Error reading line %d of file \"%s\".",
              *line_num, ecs_file_get_name(f));

  /* Slide buffer */

  for (p1 = 0; line[p2] != '\0'; p1++, p2++)
    line[p1] = line[p2];
  line[p1] = '\0';

  return res;
}

/*----------------------------------------------------------------------------*
 * Read connectivity array
 *----------------------------------------------------------------------------*/

static void
_read_elements(ecs_file_t    *hex_file,         /* --> Input file descriptor */
               int           *line_num,         /* <-> Line counter */
               int           *n_vertices,
               int           *n_edges,
               int           *n_faces,
               int           *n_cells,
               double       **vertex_coords,
               int32_t      **edge_vertices,
               int32_t      **face_edges,
               int32_t      **cell_faces,
               int32_t      **face_hierarchy,
               int32_t      **cell_hierarchy)
{
  char     line[MAX_LINE_LENGTH];
  int      n_scan;
  int      elt_count = 0;

  ecs_file_type_t  hex_file_type;

  /* Variables read from file */

  int32_t  _n_read = 0;
  int32_t  *_edge_hierarchy = NULL;

  /* Mesh type */

  hex_file_type = ecs_file_get_type(hex_file);

  /* Vertices */
  /*==========*/

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&_n_read, sizeof(int32_t), 1, hex_file);
    *n_vertices = (int)_n_read;

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH,
                  hex_file, line_num);
    n_scan = sscanf(line, " %d", n_vertices);
    if (n_scan != 1)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));

  }

  ECS_MALLOC((*vertex_coords), *n_vertices * 3, double);

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(*vertex_coords, sizeof(double),
                  *n_vertices * 3, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    static char format_lec[] = "%*d %lg %lg %lg";

    for (elt_count = 0; elt_count < *n_vertices; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH,
                    hex_file, line_num);

      n_scan = sscanf(line, format_lec,
                      &((*vertex_coords)[elt_count * 3    ]),
                      &((*vertex_coords)[elt_count * 3 + 1]),
                      &((*vertex_coords)[elt_count * 3 + 2]));
      if (n_scan != 3)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_vertices; elt_count++)
    printf("vertex  %d: %g %g %g\n", elt_count,
           (*vertex_coords)[elt_count * 3    ],
           (*vertex_coords)[elt_count * 3 + 1],
           (*vertex_coords)[elt_count * 3 + 2]);
#endif

  /* Edges */
  /*=======*/

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&_n_read, sizeof(int32_t), 1, hex_file);
    *n_edges = (int)_n_read;

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    n_scan = sscanf(line, " %d", n_edges);
    if (n_scan != 1)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));

  }

  ECS_MALLOC((*edge_vertices), *n_edges * 2, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(*edge_vertices, sizeof(int32_t),
                  *n_edges * 2, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_edges; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d",
                      &((*edge_vertices)[elt_count * 2    ]),
                      &((*edge_vertices)[elt_count * 2 + 1]));
      if (n_scan != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_edges; elt_count++)
    printf("edge %d: vertices %d and %d\n", elt_count,
           (int)((*edge_vertices)[elt_count * 2]),
           (int)((*edge_vertices)[elt_count * 2 + 1]));
#endif

  /* Edge hierarchy */

  ECS_MALLOC(_edge_hierarchy, *n_edges * 2, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(_edge_hierarchy, sizeof(int32_t),
                  *n_edges * 2, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_edges; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d",
                      &(_edge_hierarchy[elt_count * 2    ]),
                      &(_edge_hierarchy[elt_count * 2 + 1]));
      if (n_scan != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

  ECS_FREE(_edge_hierarchy);

  /* Faces */
  /*=======*/

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&_n_read, sizeof(int32_t), 1, hex_file);
    *n_faces = (int)_n_read;

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    n_scan = sscanf(line, " %d", n_faces);
    if (n_scan != 1)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));

  }

  ECS_MALLOC((*face_edges), *n_faces * 4, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(*face_edges, sizeof(int32_t), *n_faces * 4, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_faces; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH,
                    hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d %d %d",
                      &((*face_edges)[elt_count * 4    ]),
                      &((*face_edges)[elt_count * 4 + 1]),
                      &((*face_edges)[elt_count * 4 + 2]),
                      &((*face_edges)[elt_count * 4 + 3]));
      if (n_scan != 4)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_faces; elt_count++)
    printf("fac  %d : arêtes %d, %d, %d et %d\n", elt_count,
           (*face_edges)[elt_count*4    ], (*face_edges)[elt_count*4 + 1],
           (*face_edges)[elt_count*4 + 2], (*face_edges)[elt_count*4 + 3]);
#endif

  /* Face hierarchy */

  ECS_MALLOC((*face_hierarchy), *n_faces * 2, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read((*face_hierarchy), sizeof(int32_t),
                  *n_faces * 2, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_faces; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d",
                      &((*face_hierarchy)[elt_count * 2    ]),
                      &((*face_hierarchy)[elt_count * 2 + 1]));
      if (n_scan != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_faces; elt_count++)
    printf("face %d: hierarchy %d %d\n", elt_count,
           (*face_hierarchy)[elt_count * 2    ],
           (*face_hierarchy)[elt_count * 2 + 1]);
#endif

  /* Cells */
  /*=======*/

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&_n_read, sizeof(int32_t), 1, hex_file);
    *n_cells = (int)_n_read;

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    n_scan = sscanf(line, " %d", n_cells);
    if (n_scan != 1)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));

  }

  ECS_MALLOC((*cell_faces), *n_cells * 6, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read((*cell_faces), sizeof(int32_t), *n_cells * 6, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_cells; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d %d %d %d %d",
                      &((*cell_faces)[elt_count * 6    ]),
                      &((*cell_faces)[elt_count * 6 + 1]),
                      &((*cell_faces)[elt_count * 6 + 2]),
                      &((*cell_faces)[elt_count * 6 + 3]),
                      &((*cell_faces)[elt_count * 6 + 4]),
                      &((*cell_faces)[elt_count * 6 + 5]));
      if (n_scan != 6)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_cells; elt_count++)
    printf("cell  %d: faces %d, %d, %d, %d, %d and %d\n", elt_count,
           (*cell_faces)[elt_count*6    ], (*cell_faces)[elt_count*6 + 1],
           (*cell_faces)[elt_count*6 + 2], (*cell_faces)[elt_count*6 + 3],
           (*cell_faces)[elt_count*6 + 4], (*cell_faces)[elt_count*6 + 5]);
#endif

  /* Cell hierarchy */

  ECS_MALLOC((*cell_hierarchy), *n_cells * 2, int32_t);

  elt_count = 0;

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read((*cell_hierarchy), sizeof(int32_t), *n_cells * 2, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    for (elt_count = 0; elt_count < *n_cells; elt_count++) {

      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);

      n_scan = sscanf(line, "%*d %d %d",
                      &((*cell_hierarchy)[elt_count * 2    ]),
                      &((*cell_hierarchy)[elt_count * 2 + 1]));
      if (n_scan != 2)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));

    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (elt_count = 0; elt_count < *n_cells; elt_count++)
    printf("cell %d: hierarch %d %d\n", elt_count,
           (*cell_hierarchy)[elt_count * 2    ],
           (*cell_hierarchy)[elt_count * 2 + 1]);
#endif
}

/*----------------------------------------------------------------------------*
 * Read CAD links (references)
 *----------------------------------------------------------------------------*/

static void
_read_references(ecs_file_t    *hex_file,
                 int           *line_num,
                 int32_t       *n_cad_faces,
                 int32_t      **n_cad_face_faces,
                 int32_t     ***cad_face_face)
{
  int   n_scan;
  char  line[MAX_LINE_LENGTH];

  int i, j;

  /* IggHexa variables read */

  int32_t    n_domain_vertices;
  int32_t    n_domain_edges;
  int32_t    n_periodic;

  int32_t   *n_edges;
  int32_t   *n_edgevertex;
  int32_t   *n_facevertex;
  int32_t   *vertex;
  int32_t   *periodic;

  int32_t  **edge;
  int32_t  **edgevertex;
  int32_t  **facevertex;

  ecs_file_type_t  hex_file_type;

  /* Mesh type */

  hex_file_type = ecs_file_get_type(hex_file);

  /* Dimensions */
  /*============*/

  /* Number of CAD corners, curves, and faces */

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&n_domain_vertices, sizeof(int32_t), 1, hex_file);
    ecs_file_read(&n_domain_edges,   sizeof(int32_t), 1, hex_file);
    ecs_file_read(&(*n_cad_faces),   sizeof(int32_t), 1, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    n_scan = sscanf(line, " %d %d %d",
                    &n_domain_vertices, &n_domain_edges, &(*n_cad_faces));
    if (n_scan != 3)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("n_domain_vertices: %d, n_domain_edges: %d, (*n_cad_faces): %d\n",
         n_domain_vertices, n_domain_edges, (*n_cad_faces));
#endif

  ECS_MALLOC((*n_cad_face_faces), (*n_cad_faces), int32_t);
  ECS_MALLOC((*cad_face_face), (*n_cad_faces), int32_t *);

  ECS_MALLOC(n_edges, n_domain_edges, int32_t);
  ECS_MALLOC(n_edgevertex, n_domain_edges, int32_t);
  ECS_MALLOC(n_facevertex, (*n_cad_faces), int32_t);
  ECS_MALLOC(vertex, n_domain_vertices, int32_t);

  ECS_MALLOC(edge, n_domain_edges, int32_t *);
  ECS_MALLOC(edgevertex, n_domain_edges, int32_t *);
  ECS_MALLOC(facevertex, (*n_cad_faces), int32_t *);

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    if (n_domain_edges > 0)
      ecs_file_read(n_edges, sizeof(int32_t), n_domain_edges, hex_file);

    if ((*n_cad_faces) > 0)
      ecs_file_read((*n_cad_face_faces), sizeof(int32_t),
                    (*n_cad_faces), hex_file);

    if (n_domain_edges > 0)
      ecs_file_read(n_edgevertex, sizeof(int32_t), n_domain_edges, hex_file);

    if ((*n_cad_faces) > 0)
      ecs_file_read(n_facevertex, sizeof(int32_t), (*n_cad_faces), hex_file);

    if (n_domain_vertices > 0)
      ecs_file_read(vertex, sizeof(int32_t), n_domain_vertices, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    for (i = 0; i < n_domain_vertices; i++) {
      vertex[i] = _scan_int_from_line(hex_file,
                                      line_num,
                                      line,
                                      MAX_LINE_LENGTH);
    }

  }

  for (i = 0; i < n_domain_edges; i++) {

    /* If the file is binary, n_edges has been read already */
    if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      n_scan = sscanf(line, " %d", &(n_edges[i]));
      if (n_scan != 1)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));
    }

    ECS_MALLOC(edge[i], n_edges[i], int32_t);

    if (hex_file_type == ECS_FILE_TYPE_BINARY) {
      if (n_edges[i] > 0)
        ecs_file_read(edge[i], sizeof(int32_t), n_edges[i], hex_file);

    }
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      for (j = 0; j < n_edges[i]; j++) {
        edge[i][j] = _scan_int_from_line(hex_file,
                                         line_num,
                                         line,
                                         MAX_LINE_LENGTH);
      }
    }

    /* If the file is binary, n_edgevertex has been read already */
    if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      n_scan = sscanf(line, " %d", &(n_edgevertex[i]));
      if (n_scan != 1)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));
    }

    ECS_MALLOC(edgevertex[i], n_edgevertex[i], int32_t);

    if (hex_file_type == ECS_FILE_TYPE_BINARY) {
      if (n_edgevertex[i] > 0)
        ecs_file_read(edgevertex[i], sizeof(int32_t),
                      n_edgevertex[i], hex_file);
    }
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      for (j = 0; j < n_edgevertex[i]; j++) {
        edgevertex[i][j] = _scan_int_from_line(hex_file,
                                               line_num,
                                               line,
                                               MAX_LINE_LENGTH);
      }
    }

  }

  for (i = 0; i < (*n_cad_faces); i++) {

    /* If the file is binary, n_faces has been read already */
    if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      n_scan = sscanf(line, " %d", &((*n_cad_face_faces)[i]));
      if (n_scan != 1)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));
    }

    ECS_MALLOC((*cad_face_face)[i], (*n_cad_face_faces)[i], int32_t);

    if (hex_file_type == ECS_FILE_TYPE_BINARY) {
      if ((*n_cad_face_faces)[i] > 0)
        ecs_file_read((*cad_face_face)[i], sizeof(int32_t),
                      (*n_cad_face_faces)[i], hex_file);
    }
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      for (j = 0; j < (*n_cad_face_faces)[i]; j++) {
        (*cad_face_face)[i][j] = _scan_int_from_line(hex_file,
                                                     line_num,
                                                     line,
                                                     MAX_LINE_LENGTH);
      }
    }

    /* If the file is binary, n_facevertex has been read already */
    if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      n_scan = sscanf(line, " %d", &(n_facevertex[i]));
      if (n_scan != 1)
        ecs_error(__FILE__, __LINE__, errno,
                  "Error reading line %d of file \"%s\".",
                  *line_num, ecs_file_get_name(hex_file));
    }

    ECS_MALLOC(facevertex[i], n_facevertex[i], int32_t);

    if (hex_file_type == ECS_FILE_TYPE_BINARY) {
      if (n_facevertex[i] > 0)
        ecs_file_read(facevertex[i], sizeof(int32_t),
                      n_facevertex[i], hex_file);
    }
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
      for (j = 0; j < n_facevertex[i]; j++) {
        facevertex[i][j] = _scan_int_from_line(hex_file,
                                               line_num,
                                               line,
                                               MAX_LINE_LENGTH);
      }
    }

  }

  /* Periodicity (ignored) */

  if (hex_file_type == ECS_FILE_TYPE_BINARY) {

    ecs_file_read(&n_periodic, sizeof(int32_t), 1, hex_file);

  }
  else if (hex_file_type == ECS_FILE_TYPE_TEXT) {

    ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
    n_scan = sscanf(line, " %d", &n_periodic);
    if (n_scan != 1)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                *line_num, ecs_file_get_name(hex_file));
  }

  ECS_MALLOC(periodic, n_periodic * 2, int32_t);

  if (n_periodic > 0) {

    if (hex_file_type == ECS_FILE_TYPE_BINARY) {
      ecs_file_read(periodic, sizeof(int32_t), n_periodic * 2, hex_file);
    }
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      for (j = 0; j < n_facevertex[i]; j++) {
        ecs_file_gets(line, MAX_LINE_LENGTH, hex_file, line_num);
        n_scan = sscanf(line, "%d %d",
                        &(periodic[i * 2    ]),
                        &(periodic[i * 2 + 1]));
        if (n_scan != 2)
          ecs_error(__FILE__, __LINE__, errno,
                    "Error reading line %d of file \"%s\".",
                    *line_num, ecs_file_get_name(hex_file));
      }
    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)

  printf("n_edges :\n");
  for (i = 0; i < n_domain_edges; i++) {
    printf("%d\n", n_edges[i]);
    for (j = 0; j < n_edges[i]; j++)
      printf("%d", edge[i][j]);
    printf("\n");
  }

  printf("n_edgevertex:\n");
  for (i = 0; i < n_domain_edges; i++) {
    printf("%d\n", n_edgevertex[i]);
    for (j = 0; j < n_edgevertex[i]; j++)
      printf("%d", edgevertex[i][j]);
    printf("\n");
  }

  printf("n_cad_face_faces:\n");
  for (i = 0; i < (*n_cad_faces); i++) {
    printf("%d\n", (*n_cad_face_faces)[i]);
    for (j = 0; j < (*n_cad_face_faces)[i]; j++)
      printf("%d", (*cad_face_face)[i][j]);
    printf("\n");
  }

  printf("n_facevertex :\n");
  for (i = 0; i < (*n_cad_faces); i++) {
    printf("%d\n", n_facevertex[i]);
    for (j = 0; j < n_facevertex[i]; j++)
      printf("%d", facevertex[i][j]);
    printf("\n");
  }

  printf("\n\nvertex :");
  for (i = 0; i < n_domain_vertices; i++)
    printf (" %d", vertex[i]);

  printf("\n\nperiodic : %d\n", n_periodic);
  for (i = 0; i < n_periodic; i++)
    printf ("%d %d\n", periodic[i * 2], periodic[i * 2 + 1]);

#endif

  /* Free unused arrays (we only keep links between mesh faces and CAD faces) */

  ECS_FREE(vertex);

  for (i = 0; i < n_domain_edges; i++)
    ECS_FREE(edge[i]);
  ECS_FREE(edge);
  ECS_FREE(n_edges);

  for (i = 0; i < n_domain_edges; i++)
    ECS_FREE(edgevertex[i]);
  ECS_FREE(edgevertex);
  ECS_FREE(n_edgevertex);

  for (i = 0; i < (*n_cad_faces); i++)
    ECS_FREE(facevertex[i]);
  ECS_FREE(facevertex);
  ECS_FREE(n_facevertex);

  ECS_FREE(periodic);
}

/*----------------------------------------------------------------------------*
 * Assign CAD face references to mesh faces
 *----------------------------------------------------------------------------*/

static int32_t *
_process_references(int         n_faces,
                    int32_t    *face_hierarchy,
                    int32_t    *n_cad_faces,
                    int32_t   **n_cad_face_faces,
                    int32_t  ***cad_face_face)
{
  int   i, j, face_id;
  int   mod_count;

  int32_t  *reference;

  /*
    It seems that a given face may only be associated to one CAD face,
    but in the hierarchical structure, only the face closest to the root
    is given in the associations list; all its sub-faces then have
    the same association by default.
    To be verified: if this is not the case, we should either assign priority
    levels to CAD faces, or use the notion of group (allowing for overlaps).
  */

  /* Initialize array of references and free associations list */

  ECS_MALLOC(reference, n_faces, int32_t);
  for (i = 0; i < n_faces; i++)
    reference[i] = 0;

  for (i = 0; i < *n_cad_faces; i++) {

    for (j = 0; j < (*n_cad_face_faces)[i]; j++) {
      face_id = (*cad_face_face)[i][j];
      if (reference[face_id] == 0)
        reference[face_id] = i + 1;
      else
        ecs_error(__FILE__, __LINE__, 0,
                  "Face %d of IggHexa mesh is attached to at least 2 CAD faces"
                  " (%d and %d).\n"
                  "This case is not currently handled.",
                  face_id, reference[face_id], i);

    }

    ECS_FREE((*cad_face_face)[i]);

  }

  ECS_FREE(*cad_face_face);
  ECS_FREE(*n_cad_face_faces);

  *n_cad_faces = 0;

  /* Recursive propagation to sub-faces of referenced faces */

  do {

    mod_count = 0;

    for (i = 0; i < n_faces; i++) {

      if (reference[i] > 0) {

        for (j = 0; j < 2; j++) {

          face_id = face_hierarchy[i * 2  + j];
          if (face_id > -1 && reference[face_id] == 0) {
            reference[face_id] = reference[i];
            mod_count++;
          }

        }

      }

    }

  } while (mod_count > 0);

  return reference;
}

/*----------------------------------------------------------------------------*
 * Removal of hierarchic cells and associated compaction
 *----------------------------------------------------------------------------*/

static void
_compact_cells(int       *const n_vertices,
               int       *const n_edges,
               int       *const n_faces,
               int       *const n_cells,
               double   **const vertex_coords,
               int32_t  **const edge_vertices,
               int32_t  **const face_edges,
               int32_t  **const cell_faces,
               int32_t  **const face_hierarchy,
               int32_t  **const cell_hierarchy,
               int32_t  **const face_reference)
{
  int   i, j;
  int   elt_count;
  int   mod_count;

  int  *mask;
  int  *renum;

  /* Cells */
  /*=======*/

  /* We only keep cells which do not have sub-cells */

  ECS_MALLOC(mask, *n_cells, int);

  for (i = 0; i < *n_cells; i++) {
    if ((*cell_hierarchy)[i*2] != -1)
      mask[i] = 0;
    else
      mask[i] = 1;
  }

  ECS_FREE(*cell_hierarchy);

  /* Compact cell definitions */

  elt_count = 0;

  for (i = 0; i < *n_cells; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 6; j++)
        (*cell_faces)[elt_count * 6 + j]
          = (*cell_faces)[i * 6 + j];

      elt_count++;

    }

  }

  ECS_FREE(mask);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Number of initial cells (with hierarchy): %d\n"
         "Number of final cells (without hierarchy): %d\n",
         (int)(*n_cells), elt_count);
#endif

  *n_cells = elt_count;

  ECS_REALLOC(*cell_faces, (*n_cells) * 6, int32_t);

  /* Faces */
  /*=======*/

  /* First compaction: keep only faces referenced by compacted cells */

  ECS_MALLOC(mask, *n_faces, int);
  ECS_MALLOC(renum, *n_faces, int);

  for (i = 0; i < *n_faces; i++)
    mask[i] = 0;

  for (i = 0; i < *n_cells; i++) {
    for (j = 0; j < 6; j++)
      mask[(*cell_faces)[i * 6 + j]] = 1;
  }

  do {

    mod_count = 0;

    for (i = 0; i < *n_faces; i++) {

      for (j = 0; j < 2; j++) {

        if (   (mask[i] == 1)
            && ((*face_hierarchy)[i * 2 + j] != -1)
            && (mask[(*face_hierarchy)[i*2 + j]]) == 0) {
          (*face_hierarchy)[i*2 + j] = -1;
          mod_count += 1;
        }

      }

    }

  } while (mod_count != 0);


  elt_count = 0;

  for (i = 0; i < *n_faces; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 4; j++)
        (*face_edges)[elt_count * 4 + j]
          = (*face_edges)[i * 4 + j];

      (*face_reference)[elt_count] = (*face_reference)[i];

      renum[i] = elt_count;

      elt_count++;

    }
    else {

      renum[i] = -1;

    }

  }

  elt_count = 0;

  /* Compaction with hierarchy renumbering */

  for (i = 0; i < *n_faces; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 2; j++) {

        if ((*face_hierarchy)[i * 2 + j] != -1)
          (*face_hierarchy)[elt_count * 2 + j]
            = renum[(*face_hierarchy)[i * 2 + j]];
        else
          (*face_hierarchy)[elt_count * 2 + j] = -1;

      }

      elt_count++;

    }

  }

  /* Renumbering of cell -> faces connectivity */

  for (i = 0; i < *n_cells * 6; i++) {
    (*cell_faces)[i] = renum[(*cell_faces)[i]];
    assert ((*cell_faces)[i] > -1);
  }

  ECS_FREE(renum);
  ECS_FREE(mask);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Number of initial faces (with hierarchy): %d\n"
         "Number of final faces (without hierarchy): %d\n",
         (int)(*n_faces), elt_count);
#endif

  *n_faces = elt_count;

  ECS_REALLOC(*face_edges,     (*n_faces) * 4, int32_t);
  ECS_REALLOC(*face_hierarchy, (*n_faces) * 2, int32_t);
  ECS_REALLOC(*face_reference, (*n_faces),     int32_t);

  /* Edges */
  /*=======*/

  /* First compaction: keep only edges referenced by compacted faces */

  ECS_MALLOC(mask, *n_edges, int);
  ECS_MALLOC(renum, *n_edges, int);

  for (i = 0; i < *n_edges; i++)
    mask[i] = 0;

  for (i = 0; i < *n_faces; i++) {
    for (j = 0; j < 4; j++)
      mask[(*face_edges)[i * 4 + j]] = 1;
  }

  elt_count = 0;

  for (i = 0; i < *n_edges; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 2; j++)
        (*edge_vertices)[elt_count * 2 + j]
          = (*edge_vertices)[i * 2 + j];

      renum[i] = elt_count;

      elt_count++;

    }
    else {

      renum[i] = -1;

    }

  }

  /* Renumbering of face -> edges connectivity */

  for (i = 0; i < *n_faces * 4; i++) {
    (*face_edges)[i] = renum[(*face_edges)[i]];
    assert ((*face_edges)[i] > -1);
  }

  ECS_FREE(renum);
  ECS_FREE(mask);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Number of initial edges (with hierarchy): %d\n"
         "Number of final edges (without hierarchy): %d\n",
         (int)(*n_edges), elt_count);
#endif

  *n_edges = elt_count;

  ECS_REALLOC(*edge_vertices, (*n_edges) * 2, int32_t);

  /* Vertices */
  /*==========*/

  /* First compaction: keep only vertices referenced by compacted edges */

  ECS_MALLOC(mask, *n_vertices, int);
  ECS_MALLOC(renum, *n_vertices, int);

  for (i = 0; i < *n_vertices; i++)
    mask[i] = 0;

  for (i = 0; i < *n_edges; i++) {
    for (j = 0; j < 2; j++)
      mask[(*edge_vertices)[i * 2 + j]] = 1;
  }

  elt_count = 0;

  for (i = 0; i < *n_vertices; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 3; j++)
        (*vertex_coords)[elt_count * 3 + j]
          = (*vertex_coords)[i * 3 + j];

      renum[i] = elt_count;

      elt_count++;

    }
    else {

      renum[i] = -1;

    }
  }

  /* Renumbering of edge -> vertices connectivity */

  for (i = 0; i < *n_edges * 2; i++) {
    (*edge_vertices)[i] = renum[(*edge_vertices)[i]];
    assert ((*edge_vertices)[i] > -1);
  }

  ECS_FREE(renum);
  ECS_FREE(mask);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Number of initial vertices (with hierarchy): %d\n"
         "Number of final vertices (without hierarchy): %d\n",
         (int)(*n_vertices), elt_count);
#endif

  *n_vertices = elt_count;

  ECS_REALLOC(*vertex_coords, (*n_vertices) * 3, double);
}

/*----------------------------------------------------------------------------*
 * Build nodal connectivity
 *----------------------------------------------------------------------------*/

static void
_nodal_connect(int       *n_edges,
               int       *n_faces,
               int        n_cells,
               double    *vertex_coords,
               int32_t  **edge_vertices,
               int32_t  **face_edges,
               int32_t  **fac_val_som,
               int32_t  **cell_faces,
               int32_t  **cell_vertices,
               int32_t  **face_hierarchy,
               int32_t  **face_reference)
{
  int     i, j;
  int     coord_id;
  int     edge_id;
  int     cell_id;
  int     face_id;
  int     l_face_id;
  int     l_vtx_id;
  int     vtx_id;
  int     vtx_id_1;
  int     vtx_id_2;
  int     start_id;
  int     pass;
  int     direct;
  int     elt_count;
  int     mod_count;

  double  cell_ctr_coord[3];
  double  face_ctr_coord[3];
  double  v_1[3], v_2[3], cross_p[3], comp_v[3];

  int    *mask;
  int    *renum;

  /* First processing of faces */
  /*===========================*/

  /* Build nodal connectivity */

  ECS_MALLOC(*fac_val_som, (*n_faces) * 4, int32_t);

  for (face_id = 0; face_id < *n_faces; face_id++) {

    start_id = face_id * 4;

    (*fac_val_som)[start_id    ]
      = (*edge_vertices)[((*face_edges)[start_id]) * 2    ];

    (*fac_val_som)[start_id + 1]
      = (*edge_vertices)[((*face_edges)[start_id]) * 2 + 1];

    /*
      We know vertices 1 and 2, but we do not yet know which order to
      adopt to go to vertex 3.
    */

    vtx_id_1 = (*edge_vertices)[((*face_edges)[start_id + 1]) * 2    ];
    vtx_id_2 = (*edge_vertices)[((*face_edges)[start_id + 1]) * 2 + 1];

    if (   vtx_id_1 == (*fac_val_som)[start_id + 1]
        || vtx_id_2 == (*fac_val_som)[start_id + 1])
      direct = 1;
    else
      direct = -1;

    for (j = 1; j < 3; j++) {

      edge_id = (*face_edges)[start_id + ((4 + (direct * j)) % 4)];

      vtx_id_1 = (*edge_vertices)[edge_id * 2    ];
      vtx_id_2 = (*edge_vertices)[edge_id * 2 + 1];

      if (vtx_id_1 == (*fac_val_som)[start_id + j])
        (*fac_val_som)[start_id + j + 1] = vtx_id_2;
      else if (vtx_id_2 == (*fac_val_som)[start_id + j])
        (*fac_val_som)[start_id + j + 1] = vtx_id_1;
      else
        assert(   vtx_id_1 == (*fac_val_som)[start_id + j]
               || vtx_id_2 == (*fac_val_som)[start_id + j]);

    }

#if 0 && defined(DEBUG) && !defined(NDEBUG)

    printf("Face: %d\n", face_id);
    for (l_vtx_id = 0; l_vtx_id < 4; l_vtx_id++) {
      vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];
      printf("vertex %d: %f %f %f\n", vtx_id,
             vertex_coords[vtx_id * 3    ],
             vertex_coords[vtx_id * 3 + 1],
             vertex_coords[vtx_id * 3 + 2]);
    }

#endif

  }

  /* We do not need the face -> vertices connectivity anymore */

  ECS_FREE(*face_edges);

  /* Edges */
  /*=======*/

  /* We do not need the edges anymore */

  *n_edges = 0;

  ECS_FREE(*edge_vertices);

  /* Cells */
  /*=======*/

  /* Cell -> vertices connectivity */

  ECS_MALLOC(*cell_vertices, n_cells * 8, int32_t);

  for (cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Center of gravity of face centers (close to hexahedron center) */

    cell_ctr_coord[0] = 0.0; cell_ctr_coord[1] = 0.0; cell_ctr_coord[2] = 0.0;

    for (l_face_id = 0; l_face_id < 6; l_face_id++) {

      face_id = (*cell_faces)[cell_id * 6 + l_face_id];

      for (l_vtx_id = 0; l_vtx_id < 4; l_vtx_id++) {

        vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];

        for (coord_id = 0; coord_id < 3; coord_id++)
          cell_ctr_coord[coord_id] += vertex_coords[vtx_id * 3 + coord_id];

      }

    }

    for (coord_id = 0; coord_id < 3; coord_id++)
      cell_ctr_coord[coord_id] /= 24.0;

    /* Center of gravity and normal of first face */

    face_ctr_coord[0] = 0.0; face_ctr_coord[1] = 0.0; face_ctr_coord[2] = 0.0;

    face_id = (*cell_faces)[cell_id * 6];

    for (l_vtx_id = 0; l_vtx_id < 4; l_vtx_id++) {

      vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];

      for (coord_id = 0; coord_id < 3; coord_id++)
        face_ctr_coord[coord_id] += vertex_coords[vtx_id * 3 + coord_id];

    }

    for (coord_id = 0; coord_id < 3; coord_id++)
      face_ctr_coord[coord_id] /= 4.0;

    for (coord_id = 0; coord_id < 3; coord_id++) {
      v_1[coord_id]
        =   vertex_coords[(*fac_val_som)[face_id * 4 + 1] * 3 + coord_id]
          - vertex_coords[(*fac_val_som)[face_id * 4    ] * 3 + coord_id];
      v_2[coord_id]
        =   vertex_coords[(*fac_val_som)[face_id * 4 + 3] * 3 + coord_id]
          - vertex_coords[(*fac_val_som)[face_id * 4    ] * 3 + coord_id];
    }

    CROSS_PRODUCT(cross_p, v_1, v_2);

    for (coord_id = 0; coord_id < 3; coord_id++)
      comp_v[coord_id] = face_ctr_coord[coord_id] - cell_ctr_coord[coord_id];


    /*
      If the normal points outwards, we take the vertices in the opposite
      order (hexahedron "interior" face numbering).

      *     8 x-------x 7
      *      /|      /|
      *     / |     / |
      *  5 x-------x6 |
      *    | 4x----|--x 3
      *    | /     | /
      *    |/      |/
      *  1 x-------x 2

    */

    if (DOT_PRODUCT(comp_v, cross_p) > 0.0) {

      /* Exterior normal */

      for (l_vtx_id = 0; l_vtx_id < 4; l_vtx_id++) {
        vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];
        (*cell_vertices)[cell_id * 8 + 3 - l_vtx_id] = vtx_id;
      }

    }
    else {

      /* Interior normal */

      for (l_vtx_id = 0; l_vtx_id < 4; l_vtx_id++) {
        vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];
        (*cell_vertices)[cell_id * 8 + l_vtx_id] = vtx_id;
      }
    }

    /*
      We now have 4 of 8 vertices; we now search for the face sharing
      vertices 1 and 2, so as to determine vertices 5 and 6, then for
      the face sharing vertices 3 and 4, to determine vertices 7 and 8.
    */

    for (pass = 0; pass < 2; pass++) {

      vtx_id_1 = (*cell_vertices)[cell_id * 8     + (pass * 2)];
      vtx_id_2 = (*cell_vertices)[cell_id * 8 + 1 + (pass * 2)];

      direct = 0;

      for (l_face_id = 1;
           l_face_id < 6 && direct == 0;
           l_face_id++) {

        face_id = (*cell_faces)[cell_id * 6 + l_face_id];

        for (l_vtx_id = 0;
             l_vtx_id < 4;
             l_vtx_id++) {

          vtx_id = (*fac_val_som)[face_id * 4 + l_vtx_id];

          if (vtx_id == vtx_id_1) {

            if ((*fac_val_som)[face_id * 4 + ((l_vtx_id + 1) % 4)]
                == vtx_id_2) {
              direct = 1;
              break;
            }

            else if ((*fac_val_som)[face_id * 4 + ((l_vtx_id - 1 + 4) % 4)]
                     == vtx_id_2) {
              direct = -1;
              break;
            }
          }
        }
      }

      (*cell_vertices)[cell_id * 8 + 4 + (pass * 2)]
        =(*fac_val_som)[face_id * 4
                        + ((l_vtx_id + (4 + (3 * direct))) % 4)];

      (*cell_vertices)[cell_id * 8 + 5 + (pass * 2)]
        =(*fac_val_som)[face_id * 4
                        + ((l_vtx_id + (4 + (2 * direct))) % 4)];

    }

#if 0 && defined(DEBUG) && !defined(NDEBUG)

    printf("Cell: %d\n", cell_id);
    for (l_vtx_id = 0; l_vtx_id < 8; l_vtx_id++) {
      vtx_id = (*cell_vertices)[cell_id * 8 + l_vtx_id];
      printf("vertex %d: %f %f %f\n", vtx_id,
             vertex_coords[vtx_id * 3    ],
             vertex_coords[vtx_id * 3 + 1],
             vertex_coords[vtx_id * 3 + 2]);
    }

#endif

  }

  /* We do not need the cell -> faces connectivity anymore */

  ECS_FREE(*cell_faces);

  /* Second processing of faces */
  /*============================*/

  /*
    We do not need faces anymore which do not intervene in the hierarchy
    and do not bear a reference
  */

  ECS_MALLOC(mask, *n_faces, int);
  ECS_MALLOC(renum, *n_faces, int);

  for (i = 0; i < *n_faces; i++)
    if ((*face_reference)[i] != 0)
      mask[i] = 1;
    else
      mask[i] = 0;

  for (i = 0; i < *n_faces; i++) {
    if ((*face_hierarchy)[i * 2] != -1) {
      mask[i] = 1;
      mask[(*face_hierarchy)[i * 2]] = 1;
    }
  }

  do {

    mod_count = 0;

    for (i = 0; i < *n_faces; i++) {

      if (mask[i] == 1 && (*face_hierarchy)[i * 2 + 1] != -1) {

        if (mask[(*face_hierarchy)[i * 2 + 1]] == 0) {

          mask[(*face_hierarchy)[i * 2 + 1]] = 1;
          mod_count += 1;

        }

      }

    }

  } while (mod_count != 0);

  elt_count = 0;

  for (i = 0; i < *n_faces; i++) {

    if (mask[i] == 1) {

      for (j = 0; j < 4; j++)
        (*fac_val_som)[elt_count * 4 + j]
          = (*fac_val_som)[i * 4 + j];

      for (j = 0; j < 2; j++)
        (*face_hierarchy)[elt_count * 2 + j]
          = (*face_hierarchy)[i * 2 + j];

      (*face_reference)[elt_count] = (*face_reference)[i];

      renum[i] = elt_count;

      elt_count++;

    }
    else {

      renum[i] = -1;

    }
  }

  for (i = 0; i < elt_count * 2; i++) {
    if ((*face_hierarchy)[i] != -1)
      (*face_hierarchy)[i] = renum[(*face_hierarchy)[i]];
  }

  ECS_FREE(mask);
  ECS_FREE(renum);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  printf("Number of initial faces (descending connectivity): %d\n"
         "Number of final faces (nodal connectivity): %d\n",
         *n_faces, elt_count);
#endif

  *n_faces = elt_count;

  ECS_REALLOC(*face_reference, (*n_faces), int32_t);
  ECS_REALLOC(*face_hierarchy, (*n_faces)*2, int32_t);
  ECS_REALLOC(*fac_val_som, (*n_faces)*4, int32_t);

  /*
    Shift connectivities (increment by 1).
    Hierarchies are handled in a specific function.
  */

  for (i = 0; i < n_cells * 8; i++)
    (*cell_vertices)[i] += 1;

  for (i = 0; i < (*n_faces) * 4; i++)
    (*fac_val_som)[i] += 1;
}

/*----------------------------------------------------------------------------*
 * Write mesh strucure
 *----------------------------------------------------------------------------*/

static void
_write_connect(med_idt          f,
               char            *mesh_name,
               int              n_vertices,
               int              n_faces,
               int              n_cells,
               double   **const vertex_coords,
               int      **const face_vertices,
               int      **const cell_vertices)
{
  med_int i;

  char coord_name[3 * MED_TAILLE_PNOM + 1];
  char unit_name[3 * MED_TAILLE_PNOM + 1];

  med_err retval = 0;

  med_int _n_vertices = n_vertices;
  med_int _n_faces = n_faces;
  med_int _n_cells = n_cells;

  med_float  *_vertex_coords = (*vertex_coords);
  med_int    *_face_vertices = NULL;
  med_int    *_cell_vertices = NULL;

  /* Vertices */
  /*----------*/

  assert(sizeof(double) == sizeof(med_float));

  /* Coordinate names and units */

  for (i = 0; i < (med_int)(3 * MED_TAILLE_PNOM); i++)
    coord_name[i] = ' ', unit_name[i] = ' ';

  coord_name[3 * MED_TAILLE_PNOM] =  '\0';
  unit_name[3 * MED_TAILLE_PNOM] = '\0';

  coord_name[0] = 'x';
  coord_name[MED_TAILLE_PNOM] = 'y';
  coord_name[2*MED_TAILLE_PNOM] = 'z';

  /* No vertex families needed for implicit 0 with MED > 2.1,
     so only coordinates are needed */

  retval = MEDcoordEcr(f,
                       mesh_name,
                       3,
                       _vertex_coords,
                       MED_FULL_INTERLACE,
                       _n_vertices,
                       MED_CART,
                       coord_name,
                       unit_name);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDcoordEcr() failed to write coords:\n");

  /* Vertex coordinates are not needed anymore */

  _vertex_coords = NULL;
  ECS_FREE(*vertex_coords);

  /* Faces */
  /*-------*/

  ECS_MALLOC(_face_vertices, _n_faces*4, med_int);

  for (i = 0; i < _n_faces; i++) {
    _face_vertices[i*4] = (*face_vertices)[i*4];
    _face_vertices[i*4 + 1] = (*face_vertices)[i*4 + 1];
    _face_vertices[i*4 + 2] = (*face_vertices)[i*4 + 2];
    _face_vertices[i*4 + 3] = (*face_vertices)[i*4 + 3];
  }

  retval  = MEDconnEcr(f,
                       mesh_name,
                       2, /* entity_dim */
                       _face_vertices,
                       MED_FULL_INTERLACE,
                       _n_faces,
                       MED_MAILLE,
                       MED_QUAD4,
                       MED_NOD);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDconnEcr() failed to write face connectivity:\n");

  ECS_FREE(_face_vertices);

  /* Face definitions are not needed anymore */

  ECS_FREE(*face_vertices);

  /* Cells */
  /*-------*/

  ECS_MALLOC(_cell_vertices, _n_cells*8, med_int);

  for (i = 0; i < _n_cells; i++) {
    _cell_vertices[i*8] = (*cell_vertices)[i*8];
    _cell_vertices[i*8 + 1] = (*cell_vertices)[i*8 + 3];
    _cell_vertices[i*8 + 2] = (*cell_vertices)[i*8 + 2];
    _cell_vertices[i*8 + 3] = (*cell_vertices)[i*8 + 1];
    _cell_vertices[i*8 + 4] = (*cell_vertices)[i*8 + 4];
    _cell_vertices[i*8 + 5] = (*cell_vertices)[i*8 + 7];
    _cell_vertices[i*8 + 6] = (*cell_vertices)[i*8 + 6];
    _cell_vertices[i*8 + 7] = (*cell_vertices)[i*8 + 5];
  }

  retval  = MEDconnEcr(f,
                       mesh_name,
                       3, /* entity_dim */
                       _cell_vertices,
                       MED_FULL_INTERLACE,
                       _n_cells,
                       MED_MAILLE,
                       MED_HEXA8,
                       MED_NOD);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDconnEcr() failed to write cell connectivity:\n");

  ECS_FREE(_cell_vertices);

  /* Cell definitions are not needed anymore */

  ECS_FREE(*cell_vertices);
}

/*----------------------------------------------------------------------------*
 * Write mesh strucure
 *----------------------------------------------------------------------------*/

static void
_write_families(med_idt          f,
                char            *mesh_name,
                int              n_faces,
                int      **const face_reference
)
{
  med_int i;

  char family_name[MED_TAILLE_NOM + 1];
  char att_descr[MED_TAILLE_DESC + 1];

  med_err retval = 0;

  med_int next_family_id = 0, max_reference = 0;
  med_int _n_faces = n_faces;

  med_int    *face_family = NULL;

  size_t     *family_id = NULL;

  /* Initialize MED families data */
  /*------------------------------*/

  for (i = 0; i < MED_TAILLE_DESC; i++)
    att_descr[i] = ' ';
  att_descr[MED_TAILLE_DESC] = '\0';

  sprintf(family_name, "FAMILLE_%d", 0) ;

  retval = MEDfamCr(f,
                    mesh_name,
                    family_name,
                    0,
                    NULL,
                    NULL,
                    NULL,
                    0,
                    NULL,
                    0);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDfamEcr() failed to write family 0:\n");

  /* Build face families */

  max_reference = 0;

  for (i = 0; i < n_faces; i++) {
    if ((*face_reference)[i] > max_reference)
      max_reference = (*face_reference)[i];
  }

  ECS_MALLOC(family_id, max_reference + 1, size_t);

  for (i = 0; i < max_reference + 1; i++)
    family_id[i] = 0;

  for (i = 0; i < _n_faces; i++)
    family_id[(*face_reference)[i]] += 1;

  /* Define families and transform from counter to id */

  next_family_id = -1;

  for (i = 0; i < max_reference + 1; i++) {

    if (family_id[i] > 0) {

      med_int family_val = i;

      family_id[i] = next_family_id;

      sprintf(family_name, "FAMILLE_ELEMENT_%d", (int)next_family_id) ;

      retval = MEDfamCr(f,
                        mesh_name,
                        family_name,
                        next_family_id,
                        &next_family_id, /* not useful */
                        &family_val,
                        att_descr,
                        1,
                        NULL,
                        0);

      if (retval < 0)
        ecs_error(__FILE__, __LINE__, 0,
                  "MEDfamEcr() failed to write family %d:\n",
                  (int)family_val);

      next_family_id--;

    }
  }

  ECS_MALLOC(face_family, _n_faces, med_int);
  for (i = 0; i < _n_faces; i++)
    face_family[i] = family_id[(*face_reference)[i]];

  ECS_FREE(family_id);

  /* Write face family numbers */

  retval  = MEDfamEcr(f,
                      mesh_name,
                      face_family,
                      _n_faces,
                      MED_MAILLE,
                      MED_QUAD4);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDfamEcr() failed to write face families:\n");

  ECS_FREE(face_family);

  /* Face references are not needed anymore */

  ECS_FREE(*face_reference);
}

/*----------------------------------------------------------------------------
 * Transform hierarchy into a face equivalence array
 *----------------------------------------------------------------------------*/

static void
_write_equivalence(med_idt    f,
                   char      *mesh_name,
                   int        n_faces,
                   int32_t  **face_hierarchy)
{
  int     i;
  int     sub_id;
  int     sibling_id;
  int     mod_count;

  med_int   n_equiv;
  med_int  *equiv;

  int      *parent;

  char equiv_name[MED_TAILLE_NOM+1] = "Face Connectivity";
  char equiv_desc[MED_TAILLE_NOM+1] = "parent -> sub face relation";

  _Bool   propagate = false;

  med_err retval = 0;

  /* Initialize parent -> sub face relation */

  ECS_MALLOC(parent, n_faces, int);

  for (i = 0; i < n_faces; i++)
    parent[i] = 0;

  for (i = 0; i < n_faces; i++) {
    sub_id = (*face_hierarchy)[i * 2];
    if (sub_id != -1) {
      parent[sub_id] = i + 1;
      if (propagate == false)
        propagate = true;
    }
  }

  if (propagate == false) {
    ECS_FREE(*face_hierarchy);
    ECS_FREE(parent);
    return;
  }

  /* Recursive propagation from sibling to sibling */

  do {

    for (i = 0, mod_count = 0; i < n_faces; i++) {
      if (parent[i] > 0) {
        sibling_id = (*face_hierarchy)[i * 2 + 1];
        if (sibling_id > -1 && parent[sibling_id] == 0) {
          parent[sibling_id] = parent[i];
          mod_count++;
        }
      }
    }

  } while (mod_count > 0);

  ECS_FREE(*face_hierarchy);

  /* Transform parent relation to equivalence */

  for (i = 0, n_equiv = 0; i < n_faces; i++) {
    if (parent[i] != 0)
      n_equiv++;
  }

  ECS_MALLOC(equiv, n_equiv * 2, med_int);

  for (i = 0, n_equiv = 0; i < n_faces; i++) {
    if (parent[i] != 0) {
      equiv[n_equiv*2] = parent[i];
      equiv[n_equiv*2 + 1] = i + 1;
      n_equiv++;
    }
  }

  ECS_FREE(parent);

  /* Add equivalence to mesh model */

  retval = MEDequivCr(f, mesh_name, equiv_name, equiv_desc);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDequivCr failed for equivalence \"%s\".",
              equiv_name);

  retval = MEDequivEcr(f,
                       mesh_name,
                       equiv_name,
                       equiv,
                       n_equiv,
                       MED_MAILLE,
                       MED_QUAD4);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDequivEcr failed for equivalence \"%s\".",
              equiv_name);

  ECS_FREE(equiv);
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------*
 *  Lecture d'un fichier Hex NUMECA
 *   et affectation des donnees dans la structure de maillage
 *----------------------------------------------------------------------------*/

int
main(int    argc,
     char  *argv[])
{
  int        n_scan;
  char       line[MAX_LINE_LENGTH];
  int        line_num;

  int        n_vertices = 0;
  int        n_edges = 0;
  int        n_faces = 0;
  int        n_cells = 0;
  double    *vertex_coords = NULL;
  int       *edge_vertices = NULL;
  int       *face_edges = NULL;
  int       *face_vertices = NULL;
  int       *cell_faces = NULL;
  int       *cell_vertices = NULL;
  int       *face_hierarchy = NULL;
  int       *cell_hierarchy = NULL;
  int       *face_reference = NULL;

  int        n_cad_faces;
  int       *n_cad_face_faces = NULL;
  int      **cad_face_face = NULL;

  char       mesh_name[] = "Fluid Domain";
  char       mesh_info[MED_TAILLE_DESC + 1] = "Generated by IggHexatoMED";

  ecs_file_type_t  hex_file_type = ECS_FILE_TYPE_BINARY;
  ecs_file_t      *hex_file = NULL;
  med_idt          med_file = -1;
  med_err          retval = 0;

  if  (argc < 3) {
    printf("Usage:\n%s input_file_name output_file_name\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  ecs_mem_init(getenv("IGGHEXA_TO_MED_MEM_LOG"));

  /* Print title */
  /*=============*/

  printf("\n\n"
         "Reading mesh from file in NUMECA Hex format\n"
         "----------------------\n");

  printf("  Mesh file: %s\n\n", argv[1]);

  /* Initialization */
  /*================*/

  line_num = 1;

  /* Open input IggHexa file */
  /*-------------------------*/

  hex_file = ecs_file_open(argv[1],
                           ECS_FILE_MODE_READ,
                           ECS_FILE_TYPE_BINARY);

  ecs_file_set_big_endian(hex_file);

  /* Read header */

  {
    char     *header;
    int32_t   headersize;

    char * version;
    char * fileversion;
    char * tmpdate;

    ecs_file_read(line, sizeof(char), 2, hex_file);

    if (line[0] == '0') {

      ecs_file_free(hex_file);

      hex_file = ecs_file_open(argv[1],
                               ECS_FILE_MODE_READ,
                               ECS_FILE_TYPE_TEXT);
      ecs_file_read(line, sizeof(char), 2, hex_file);

    }

    hex_file_type = ecs_file_get_type(hex_file);

    if (hex_file_type == ECS_FILE_TYPE_BINARY)
      ecs_file_read(&headersize, sizeof(int32_t), 1, hex_file);
    else if (hex_file_type == ECS_FILE_TYPE_TEXT) {
      ecs_file_read(line, sizeof(char), 4, hex_file);
      sscanf(line, "%4d", &headersize);
    }

    ECS_MALLOC(header, headersize + 1, char);
    ECS_MALLOC(version, headersize + 1, char);
    ECS_MALLOC(fileversion, headersize + 1, char);
    ECS_MALLOC(tmpdate, headersize + 1, char);

    if (hex_file_type == ECS_FILE_TYPE_BINARY)
      ecs_file_read(header, sizeof(char), headersize, hex_file);
    else if (hex_file_type == ECS_FILE_TYPE_TEXT)
      ecs_file_gets(header, headersize, hex_file, &line_num);

    n_scan = sscanf(header, "%*s %*s %*s %*s %s %*s %*s %s %s",
                    version, fileversion, tmpdate);
    if (n_scan != 3)
      ecs_error(__FILE__, __LINE__, errno,
                "Error reading line %d of file \"%s\".",
                1, argv[1]);

    printf("  IggHexa %s, format version  %s, %s\n\n",
           version, fileversion, tmpdate);

    ECS_FREE(header);
    ECS_FREE(version);
    ECS_FREE(fileversion);
    ECS_FREE(tmpdate);


    if (hex_file_type == ECS_FILE_TYPE_BINARY)
      ecs_file_read(line, sizeof(char), 1, hex_file);
    else if (hex_file_type == ECS_FILE_TYPE_TEXT)
      ecs_file_gets(line, MAX_LINE_LENGTH,
                    hex_file, &line_num);

    assert (line[0] == '0' || line[0] == '1');

    if (line[0] == '0')
      ecs_error(__FILE__, __LINE__, 0,
                "The IggHexa file %s\ncontains no mesh",
                argv[1]);

  }

  /* Read vertices and connectivity */
  /*--------------------------------*/

  _read_elements(hex_file,
                 &line_num,
                 &n_vertices,
                 &n_edges,
                 &n_faces,
                 &n_cells,
                 &vertex_coords,
                 &edge_vertices,
                 &face_edges,
                 &cell_faces,
                 &face_hierarchy,
                 &cell_hierarchy);

  /* Read flag indicating if we have CAD references */

  if (hex_file_type == ECS_FILE_TYPE_BINARY)
    ecs_file_read(line, sizeof(char), 1, hex_file);
  else if (hex_file_type == ECS_FILE_TYPE_TEXT)
    ecs_file_gets(line, MAX_LINE_LENGTH,
                  hex_file, &line_num);

  assert (line[0] == '0' || line[0] == '1');

  if (line[0] == '0') {
    printf("Warning:\n"
           "--------\n"
           "  The IggHexa file %s\n"
           "  contains no links with CAD.\n"
           "  No face will be referenced",
           argv[1]);
  }

  /* Read CAD references */
  /* ------------------- */

  _read_references(hex_file,
                   &line_num,
                   &n_cad_faces,
                   &n_cad_face_faces,
                   &cad_face_face);

  /* Close input file */
  /*------------------*/

  ecs_file_free(hex_file);

  /* Process then free CAD references */
  /*----------------------------------*/

  face_reference = _process_references(n_faces,
                                       face_hierarchy,
                                       &n_cad_faces,
                                       &n_cad_face_faces,
                                       &cad_face_face);


  /* Prepare mesh */
  /*--------------*/

  _compact_cells(&n_vertices,
                 &n_edges,
                 &n_faces,
                 &n_cells,
                 &vertex_coords,
                 &edge_vertices,
                 &face_edges,
                 &cell_faces,
                 &face_hierarchy,
                 &cell_hierarchy,
                 &face_reference);

  _nodal_connect(&n_edges,
                 &n_faces,
                 n_cells,
                 vertex_coords,
                 &edge_vertices,
                 &face_edges,
                 &face_vertices,
                 &cell_faces,
                 &cell_vertices,
                 &face_hierarchy,
                 &face_reference);

  /* Open output MED file */
  /*----------------------*/

  printf("\n"
         "Writing mesh to file in MED format\n"
         "--------------------\n");

  printf("  Mesh file: %s\n\n", argv[2]);


  med_file = MEDouvrir(argv[2], MED_CREATION) ;

  if (med_file < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "Error opening MED file \"%s\".",
              argv[2]);

  retval = MEDmaaCr(med_file,
                    mesh_name,
                    3,
                    MED_NON_STRUCTURE,
                    mesh_info);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDmaaCr() failed to create a new med_mesh.n");

  retval = MEDdimEspaceCr(med_file,
                          mesh_name,
                          3);

  if (retval < 0)
    ecs_error(__FILE__, __LINE__, 0,
              "MEDdimEspaceCr() failed to set spacial dimension.");

  _write_connect(med_file,
                 mesh_name,
                 n_vertices,
                 n_faces,
                 n_cells,
                 &vertex_coords,
                 &face_vertices,
                 &cell_vertices);

  _write_families(med_file,
                  mesh_name,
                  n_faces,
                  &face_reference);

  _write_equivalence(med_file,
                     mesh_name,
                     n_faces,
                     &face_hierarchy);

  ecs_mem_end();

  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------*/
