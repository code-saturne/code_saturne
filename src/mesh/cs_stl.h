#ifndef __CS_STL_H__
#define __CS_STL_H__

/*============================================================================
 * STL reader and writer.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an STL mesh
 *----------------------------------------------------------------------------*/

typedef struct {

  char           name[20];         /*!< Name identifier of the STL file*/

  char           header[80];       /*!< Header of the STL file */

  cs_lnum_t      n_faces;          /*!< Number of triangles */

  cs_real_3_t   *coords;           /*!< Array of face vertex coordinates
                                   *  coord[n_vertices][3] */

  cs_real_3_t   *coords_ini;       /*!< Array of face vertex coordinates
                                   *  at init coord_ini[n_vertices][3] */

  int            n_seeds;          /*!< Number of prescribed points
                                   *  outside the STL */

  cs_real_t     *seed_coords;      /*!< Coordinates of the reference seed points
                                    *  seed_coords[n_seed][3] */

  bool           is_porous;        /*!< If true the STL is used for porosity
                                      computation. (Default : False) */

  fvm_nodal_t    *ext_mesh;         /*!< Associated external mesh */

} cs_stl_mesh_t ;

/*----------------------------------------------------------------------------
 * Structure containing the number of STL meshes and their associated pointers
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_stl_mesh_t  **mesh_list;  /*!< Array of STL meshes
                                *   size: n_meshes*/

  cs_lnum_t        n_meshes;   /*!< Total number of STL meshes */

  int              writer_id;  /*!< Writer id, if postprocessing needed */

} cs_stl_mesh_info_t ;

/*=============================================================================
 * Static global variables
 *============================================================================*/

extern cs_stl_mesh_info_t  *cs_glob_stl_meshes;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Add a new STL mesh structure.
 *
 * parameters:
 *   path <-- path of the STL mesh
 *
 * returns:
 *   pointer to the new STL mesh structure
 *----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_add(const char  *path);

/*----------------------------------------------------------------------------
 * Return a pointer to a STL mesh based on its name if present.
 *
 * parameters:
 *   name <-- name of the STL mesh
 *
 * returns:
 *   pointer to the STL mesh structure, or NULL
 *----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_get_by_name(const char  *name);

/*----------------------------------------------------------------------------
 * Free all allocated STL mesh structures
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_destroy_all(void);

/*----------------------------------------------------------------------------
 * Read a binary STL file and store its content in a STL mesh structure.
 *
 * parameters:
 *   stl_mesh  <-- pointer to the associated STL mesh structure
 *   path      <-- path to the STL file
 *----------------------------------------------------------------------------*/

void
cs_stl_file_read(cs_stl_mesh_t  *stl_mesh,
                 const char     *path);

/*----------------------------------------------------------------------------
 * Apply a transformation matrix to a STL mesh structure.
 *
 * parameters:
 *   stl_mesh        <-- pointer to the associated STL mesh structure
 *   matrix          <-- transformation matrix
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_transform(cs_stl_mesh_t  *stl_mesh,
                      double          matrix[3][4]);

/*----------------------------------------------------------------------------
 * Apply a transformation matrix to a STL mesh structure, but use
 * the initial coordinates as inputs
 *
 * parameters:
 *   stl_mesh        <-- pointer to the associated STL mesh structure
 *   matrix          <-- transformation matrix
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_transform_from_init(cs_stl_mesh_t  *stl_mesh,
                                double          matrix[3][4]);

/*----------------------------------------------------------------------------
 * Apply a translation to a STL mesh structure.
 *
 * parameters:
 *   stl_mesh        <-- pointer to the associated STL mesh structure
 *   vector          <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_translate(cs_stl_mesh_t  *stl_mesh,
                      double          vector[3]);

/*----------------------------------------------------------------------------
 * Apply a rotation to a STL mesh structure.
 *
 * parameters:
 *   stl_mesh        <-- pointer to the associated STL mesh structure
 *   tehta           <-- rotation angle
 *   axis            <-- rotation axis
 *   center          <-- rotation center
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_rotate(cs_stl_mesh_t  *stl_mesh,
                   double          theta,
                   double          axis[3],
                   double          center[3]);

/*----------------------------------------------------------------------------
 * Apply a scaling to a STL mesh structure.
 *
 * parameters:
 *   stl_mesh        <-- pointer to the associated STL mesh structure
 *   scale           <-- scale factor
 *----------------------------------------------------------------------------*/

void
cs_stl_mesh_scale(cs_stl_mesh_t  *stl_mesh,
                  double          scale);

/*----------------------------------------------------------------------------
 * Set the points outside he STL geometry. Those points will be used as seed
 * to propagate porosity values outside the STL geometry.
 *
 * parameters:
 *   stl_mesh  <-- pointer to the associated STL mesh structure
 *   n_points  <-- number of points
 *   coords    <-- coordinates (x1,y1,z1,...,xn,yn,zn)
 *                 (size : 3*n_point)
 *----------------------------------------------------------------------------*/

void
cs_stl_set_porosity_seed(cs_stl_mesh_t  *stl_mesh,
                         int            n_points,
                         cs_real_t      *coords);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return writer_id used for stl meshes. 0 means unused.
 */
/*----------------------------------------------------------------------------*/

int
cs_stl_post_get_writer_id(void);

/*----------------------------------------------------------------------------
 * Create a new writer that will contains the STL mesh added by the user
 * The writer_id is stored in the cs_glob_stl_meshes structure.
 *
 * case_name        associated case name
 * dir_name         associated directory name
 * fmt_name         associated format name
 * fmt_opts         associated format options string
 * time_dep         FVM_WRITER_FIXED_MESH if mesh definitions
 *                  are fixed, FVM_WRITER_TRANSIENT_COORDS if
 *                  coordinates change, FVM_WRITER_TRANSIENT_CONNECT if
 *                  connectivity changes
 * output_at_start  force output at calculation start if true
 * output_at_end    force output at calculation end if true
 * frequency_n      default output frequency in time-steps, or < 0
 * frequency_t      default output frequency in secon
 *----------------------------------------------------------------------------*/

void
cs_stl_post_init_writer(const char             *case_name,
                        const char             *dir_name,
                        const char             *fmt_name,
                        const char             *fmt_opts,
                        fvm_writer_time_dep_t   time_dep,
                        bool                    output_at_start,
                        bool                    output_at_end,
                        int                     frequency_n,
                        double                  frequency_t);

/*----------------------------------------------------------------------------
 * Associate a STL mesh to the default writer
 *
 * stl_mesh <-- pointer to the associated STL mesh structure
 *----------------------------------------------------------------------------*/

void
cs_stl_post_add_mesh(cs_stl_mesh_t  *stl_mesh);

/*----------------------------------------------------------------------------
 * Write a binary STL file according to a given STL mesh structure.
 *
 * parameters:
 *   stl_mesh  <-- pointer to the associated STL mesh structure
 *   path      <-- path to the STL file
 *----------------------------------------------------------------------------*/

void
cs_stl_file_write(cs_stl_mesh_t  *stl_mesh,
                  const char     *path);

/*----------------------------------------------------------------------------
 * Compute intersection between a STL mesh and the main mesh.
 *
 * parameters:
 *   stl_mesh         <-- pointer to the associated STL mesh structure
 *   n_input          <-- number of cells on which intersection is done
 *   input_index      <-- index of input cells (size: n_input)
 *   n_selected_cells --> number of intersecting cells
 *   selected_cells   --> index of output cells (size: n_output)
 *   tria_in_cell_idx --> start index of triangle intersecting each cell
 *                        (size: n_output)
 *   tria_in_cell_lst --> list of triangles in intersecting cells
 *   max_size         --> maximum size of tria_in_cell_lst array
 *----------------------------------------------------------------------------*/

void
cs_stl_intersection(cs_stl_mesh_t *stl_mesh,
                    cs_lnum_t     n_input,
                    cs_lnum_t     *input_idx,
                    cs_lnum_t     *n_selected_cells,
                    cs_lnum_t     *selected_cells,
                    cs_lnum_t     *tria_in_cell_idx,
                    cs_lnum_t     **tria_in_cell_lst,
                    cs_lnum_t     *max_size);

/*----------------------------------------------------------------------------
 * Refine the mesh following a given STL mesh
 *
 * parameters:
 *   stl_mesh       <-- pointer to the associated STL mesh structure
 *   n_ref          <-- level of refinement
 *   n_add_layer    <-- additional layers between two refinement stage
 *----------------------------------------------------------------------------*/

void
cs_stl_refine(cs_stl_mesh_t *stl_mesh,
              int           n_ref,
              int           n_add_layer);

/*----------------------------------------------------------------------------
 * Compute porosity field according to a given STL mesh
 *
 * parameters:
 *   stl_mesh    <-- pointer to the associated STL mesh structure
 *   n_ref_point <-- number of prescribed points outside the STL
 *   ref_coords  <-- coordinates of the reference points :
 *                   x1,y1,z1,...,xn,yn,zn (size : 3*n_ref_point)
 *   porosity    --> interpolated porosity field
 *   indic       --> indicator of the STL location (NULL if not needed)
 *----------------------------------------------------------------------------*/

void
cs_stl_compute_porosity(cs_stl_mesh_t *stl_mesh,
                        cs_real_t     *porosity,
                        int           *indic);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_STL_H__ */
