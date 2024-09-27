/*============================================================================
 * Cartesian mesh generation
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_block_dist.h"
#include "cs_math.h"

#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_cartesian.h"

BEGIN_C_DECLS

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* parameters for a direction (x, y or z) */
/*----------------------------------------------------------------------------*/

typedef struct {

  /* Law type: Constant, geometric, parabolic or user */
  cs_mesh_cartesian_law_t  law;

  /* Number of cells */
  int                      ncells;

  /* Min and max coordinates */
  cs_real_t                smin;
  cs_real_t                smax;

  /* Progression, used only for geometric or parabolic laws */
  cs_real_t                progression;

  /* Two possibilities :
   *  - If constant law, this is an array of size 1 containing the step
   *  - Else, array of size ncells + 1, containing vertex coordinates.
   */
  cs_real_t               *s;
} _cs_mesh_cartesian_direction_t ;

/*----------------------------------------------------------------------------*/
/* Cartesian mesh parameters structure */
/*----------------------------------------------------------------------------*/

struct _cs_mesh_cartesian_params_t {

  // Name of cartesian block
  char                           *name;

  // Id of the mesh
  int                             id;

  /* Number of direction, set to 3 by default */
  int                             ndir;

  /* Array of size ndir, containing parameters for each direction */
  _cs_mesh_cartesian_direction_t **params;

  /* Index shifting for group id */
  int                             gc_id_shift;

  /* global values */
  cs_gnum_t n_g_cells;
  cs_gnum_t n_g_faces;
  cs_gnum_t n_g_vtx;

  cs_gnum_t n_g_cells_offset;
  cs_gnum_t n_g_faces_offset;
  cs_gnum_t n_g_vtx_offset;

  cs_gnum_t n_cells_on_rank;
  cs_gnum_t n_faces_on_rank;
  cs_gnum_t n_vtx_on_rank;

};

/*============================================================================
 * Private global variables
 *============================================================================*/

/* Flag to tell if a cartesian mesh is to be built */
static int _build_mesh_cartesian = 0;

/* Parameters for structured mesh */
static int _n_structured_meshes = 0;
static cs_mesh_cartesian_params_t **_mesh_params = NULL;

/* Flag to set a maximum number of cartesian blocks */
static int _n_structured_meshes_max = -1;

/*============================================================================
 * Private functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute intersection between two intervals
 *
 * \param[in]  i1 First interval
 * \param[in]  i2 Second interval
 * \param[out] i3 Intersection interval
 */
/*----------------------------------------------------------------------------*/

static void
_intersect_intervals(const cs_gnum_t *i1,
                     const cs_gnum_t *i2,
                     cs_gnum_t *i3)
{
  if (i1[0] > i2[1] || i2[0] > i1[1]) {
    i3[0] = -1;
    i3[1] = -1;
  }
  else {
    i3[0] = (i1[0] < i2[0]) ? i2[0] : i1[0]; // start with max
    i3[1] = (i1[1] < i2[1]) ? i1[1] : i2[1]; // end with min
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a cartesian mesh parameters based on id.
 *
 * \param[in] id  Id of the mesh parameters asked
 *
 * \returns pointer to corresponding mesh parameters.
 *          Raises an error if not found.
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_cartesian_params_t *
_get_structured_mesh_by_id(const int id)
{

  if (id < 0 || id >= _n_structured_meshes)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Out of bound id.\n"));

  return _mesh_params[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get function for structured mesh based on its name.
 *
 * \param[in] name  Name of mesh
 *
 * \returns pointer to corresponding mesh parameters, or NULL if mesh does
 *          not exist
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_cartesian_params_t *
_get_structured_mesh_by_name_try(const char *name)
{
  cs_mesh_cartesian_params_t *retval = NULL;

  if (name != NULL && strlen(name) > 0) {
    for (int i = 0; i < _n_structured_meshes; i++) {
      if (_mesh_params[i]->name != NULL &&
          strcmp(_mesh_params[i]->name, name) == 0) {
        retval = _mesh_params[i];
        break;
      }
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*! \brief Create the mesh parameters structure
 *
 * \param[in] ndir  number of directions
 *
 * \return  pointer to mesh parameters structure
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_cartesian_params_t *
_cs_mesh_cartesian_init(const char *name,
                        const int   ndir)
{
  if (_n_structured_meshes_max > 0 &&
      _n_structured_meshes == _n_structured_meshes_max)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error: Maximum number \"%d\" of cartesian mesh blocks reached.\n"),
       _n_structured_meshes_max);

  cs_mesh_cartesian_params_t *_new_mesh =
    _get_structured_mesh_by_name_try(name);

  if (_new_mesh != NULL)
    bft_error(__FILE__, __LINE__, 0,
              "Error: a mesh with name \"%s\" allready exists.\n",
              name);

  BFT_MALLOC(_new_mesh, 1, cs_mesh_cartesian_params_t);

  _new_mesh->name = NULL;
  if (name != NULL && strlen(name) > 0) {
    size_t _l = strlen(name);
    BFT_MALLOC(_new_mesh->name, _l+1, char);
    strcpy(_new_mesh->name, name);
    _new_mesh->name[_l] = '\0';
  }

  _new_mesh->gc_id_shift = 0;

  /* Global values */
  _new_mesh->n_g_cells = 0;
  _new_mesh->n_g_faces = 0;
  _new_mesh->n_g_vtx   = 0;

  _new_mesh->n_g_cells_offset = 0;
  _new_mesh->n_g_faces_offset = 0;
  _new_mesh->n_g_vtx_offset   = 0;

  _new_mesh->n_cells_on_rank = 0;
  _new_mesh->n_faces_on_rank = 0;
  _new_mesh->n_vtx_on_rank   = 0;

  _new_mesh->ndir = ndir;
  BFT_MALLOC(_new_mesh->params, ndir, _cs_mesh_cartesian_direction_t *);
  for (int i = 0; i < ndir; i++)
    _new_mesh->params[i] = NULL;

  int _id = _n_structured_meshes;
  _new_mesh->id = _id;

  _n_structured_meshes += 1;
  BFT_REALLOC(_mesh_params, _n_structured_meshes ,cs_mesh_cartesian_params_t *);

  _mesh_params[_id] = _new_mesh;

  return _new_mesh;
}

/*----------------------------------------------------------------------------*/
/*! \brief Create parameters for a direction.
 *
 * \param[in] law         1D discreization law : constant, geometric or parabolic
 * \param[in] ncells      Number of cells for this direction
 * \param[in] smin        Min coordinate value for this direction
 * \param[in] smax        Max coordinate value for this direction
 * \param[in] progression Progression value, only used for geometric or
 *                        parabolic laws.
 *
 * \return pointer to direction parameter structure
 */
/*----------------------------------------------------------------------------*/

static _cs_mesh_cartesian_direction_t *
_cs_mesh_cartesian_create_direction(cs_mesh_cartesian_law_t law,
                                    int                     ncells,
                                    cs_real_t               smin,
                                    cs_real_t               smax,
                                    cs_real_t               progression)
{
  _cs_mesh_cartesian_direction_t *dirp = NULL;

  if (smax < smin)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: smax < smin in %s\n"), __func__);

  BFT_MALLOC(dirp, 1, _cs_mesh_cartesian_direction_t);

  dirp->ncells = ncells;
  dirp->smin   = smin;
  dirp->smax   = smax;

  dirp->law = law;

  cs_real_t dir_len = smax-smin;

  if (law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
    dirp->progression = -1.;
    BFT_MALLOC(dirp->s, 1, cs_real_t);

    dirp->s[0] = dir_len / dirp->ncells;
  }
  else if (law == CS_MESH_CARTESIAN_GEOMETRIC_LAW) {
    dirp->progression = progression;
    cs_real_t rho   = dirp->progression;
    cs_real_t rho_n = pow(rho, dirp->ncells);
    cs_real_t dx0   = dir_len * (rho - 1.) / (rho_n - 1.);

    BFT_MALLOC(dirp->s, ncells+1, cs_real_t);

    cs_real_t dx_cur = dx0;
    dirp->s[0] = smin;
    for (int c_id = 0; c_id < ncells; c_id++) {
      dirp->s[c_id+1] = dirp->s[c_id] + dx_cur;
      dx_cur *= rho;
    }

  }
  else if (law == CS_MESH_CARTESIAN_PARABOLIC_LAW) {
    dirp->progression = progression;

    BFT_MALLOC(dirp->s, ncells+1, cs_real_t);
    cs_real_t rho   = dirp->progression;

    cs_real_t dx0 = 0.;

    /* We need to disguish the case of even or odd number of cells */
    int is_even = (ncells % 2 == 0);
    int np = 0;

    if (is_even) {
      np = ncells / 2;
      cs_real_t rho_np = pow(rho, np);
      dx0 = 0.5 * dir_len * (rho - 1.) / (rho_np - 1.);
    } else {
      np = (ncells - 1) / 2;
      cs_real_t rho_np = pow(rho, np);
      cs_real_t rho_np1 = rho_np * rho;
      dx0 = dir_len * (rho - 1.) / (rho_np1 + rho_np - 2.);
    }

    cs_real_t dx_cur = dx0;
    dirp->s[0]      = smin;
    dirp->s[ncells] = smax;

    for (int i = 0; i < np; i++) {
      dirp->s[i+1] = dirp->s[i] + dx_cur;
      dirp->s[ncells-i-1] = dirp->s[ncells-i] - dx_cur;

      dx_cur *= rho;
    }

  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Method not yet implemented for '%s'\n"),
              __func__);
  }

  return dirp;
}

/*----------------------------------------------------------------------------*/
/*! \brief Add a face with x-normal.
 *
 * \param[in]       mp    pointer to cartesian mesh parameters
 * \param[in, out]  mb    pointer to cs_mesh_builder_t structure
 * \param[in]       f_id  id of added face
 * \param[in]       nx    number of cells in x direction
 * \param[in]       ny    number of cells in x direction
 * \param[in]       nz    number of cells in x direction
 * \param[in]       i     i index (x) direction
 * \param[in]       j     j index (x) direction
 * \param[in]       k     k index (x) direction
 */
/*----------------------------------------------------------------------------*/

static void
_add_nx_face(cs_mesh_cartesian_params_t *mp,
             cs_mesh_builder_t          *mb,
             cs_lnum_t                   f_id,
             cs_gnum_t                   nx,
             cs_gnum_t                   ny,
             cs_gnum_t                   nz,
             cs_gnum_t                   i,
             cs_gnum_t                   j,
             cs_gnum_t                   k)
{
  CS_UNUSED(nz);

  /* Global numbering starts at 1! */

  cs_lnum_t i0 = 1 + mp->n_g_vtx_offset;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  /* Face to cell connectivity */
  cs_gnum_t c0 = i + j*nx + k*nx*ny + mp->n_g_cells_offset;

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (i == 0) {
    c_id2 = c0 + 1;
    mb->face_gc_id[f_id] = 2 + mp->gc_id_shift;
  }
  else if (i == nx) {
    c_id1 = c0;
    mb->face_gc_id[f_id] = 3 + mp->gc_id_shift;
  }
  else {
    c_id1 = c0;
    c_id2 = c0 + 1;
  }
  mb->face_cells[2*f_id]     = c_id1;
  mb->face_cells[2*f_id + 1] = c_id2;

  /*  Connectiviy for x-normal faces:
   *
   *  Vtx2        Vtx3
   *  (j,k+1)     (j+1,k+1)
   *
   *   *-----------*       z (k)
   *   |           |       ^
   *   |           |       |
   *   |     *     |       |
   *   |  (i,j,k)  |       .----->y (j)
   *   |           |
   *   *-----------*
   *  Vtx1        Vtx4
   * (j,k)        (j+1,k)
   *
   */
  mb->face_vertices[4*f_id + 3] = i0 + i + j*nxp1     + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 2] = i0 + i + j*nxp1     + (k+1)*nxp1*nyp1;
  mb->face_vertices[4*f_id + 1] = i0 + i + (j+1)*nxp1 + (k+1)*nxp1*nyp1;
  mb->face_vertices[4*f_id + 0] = i0 + i + (j+1)*nxp1 + k*nxp1*nyp1;
}

/*----------------------------------------------------------------------------*/
/*! \brief Add a face with y-normal.
 *
 * \param[in]       mp    pointer to cartesian mesh parameters
 * \param[in, out]  mb    pointer to cs_mesh_builder_t structure
 * \param[in]       f_id  id of added face
 * \param[in]       nx    number of cells in x direction
 * \param[in]       ny    number of cells in x direction
 * \param[in]       nz    number of cells in x direction
 * \param[in]       i     i index (x) direction
 * \param[in]       j     j index (x) direction
 * \param[in]       k     k index (x) direction
 */
/*----------------------------------------------------------------------------*/

static void
_add_ny_face(cs_mesh_cartesian_params_t *mp,
             cs_mesh_builder_t          *mb,
             cs_lnum_t                   f_id,
             cs_gnum_t                   nx,
             cs_gnum_t                   ny,
             cs_gnum_t                   nz,
             cs_gnum_t                   i,
             cs_gnum_t                   j,
             cs_gnum_t                   k)
{
  CS_UNUSED(ny);
  CS_UNUSED(nz);

  /* Global numbering starts at 1!
   * Adding block offset
   */
  const cs_gnum_t i0 = 1 + mp->n_g_vtx_offset;
  const cs_gnum_t c0 = 1 + mp->n_g_cells_offset;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  /* Face to cell connectivity */

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (j == 0) {
    c_id2 = c0 + i + j*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 4 + mp->gc_id_shift;
  } else if (j == ny) {
    c_id1 = c0 + i + (j-1)*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 5 + mp->gc_id_shift;
  } else {
    c_id1 = c0 + i + (j-1)*nx + k*nx*ny;
    c_id2 = c0 + i + j*nx     + k*nx*ny;
  }

  mb->face_cells[2*f_id]     = c_id1;
  mb->face_cells[2*f_id + 1] = c_id2;

  /*  Connectiviy for y-normal faces:
   *
   *  Vtx2        Vtx3
   *  (i+1,k)     (i+1,k+1)
   *
   *   *-----------*       x (i)
   *   |           |       ^
   *   |           |       |
   *   |     *     |       |
   *   |  (i,j,k)  |       .----->z (k)
   *   |           |
   *   *-----------*
   *  Vtx1        Vtx4
   * (i,k)        (i,k+1)
   *
   */
  mb->face_vertices[4*f_id + 3] = i0 + i     + j*nxp1 + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 2] = i0 + (i+1) + j*nxp1 + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 1] = i0 + (i+1) + j*nxp1 + (k+1)*nxp1*nyp1;
  mb->face_vertices[4*f_id + 0] = i0 + i     + j*nxp1 + (k+1)*nxp1*nyp1;
}

/*----------------------------------------------------------------------------*/
/*! \brief Add a face with z-normal.
 *
 * \param[in]       mp    pointer to cartesian mesh parameters
 * \param[in, out]  mb    pointer to cs_mesh_builder_t structure
 * \param[in]       f_id  id of added face
 * \param[in]       nx    number of cells in x direction
 * \param[in]       ny    number of cells in x direction
 * \param[in]       nz    number of cells in x direction
 * \param[in]       i     i index (x) direction
 * \param[in]       j     j index (x) direction
 * \param[in]       k     k index (x) direction
 */
/*----------------------------------------------------------------------------*/

static void
_add_nz_face(cs_mesh_cartesian_params_t *mp,
             cs_mesh_builder_t          *mb,
             cs_lnum_t                   f_id,
             cs_gnum_t                   nx,
             cs_gnum_t                   ny,
             cs_gnum_t                   nz,
             cs_gnum_t                   i,
             cs_gnum_t                   j,
             cs_gnum_t                   k)
{
  /* Global numbering starts at 1! */

  const cs_gnum_t i0 = 1 + mp->n_g_vtx_offset;
  const cs_gnum_t c0 = 1 + mp->n_g_cells_offset;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (k == 0) {
    c_id2 = c0 + i + j*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 6 + mp->gc_id_shift;
  } else if (k == nz) {
    c_id1 = c0 + i + j*nx + (k-1)*nx*ny;
    mb->face_gc_id[f_id] = 7 + mp->gc_id_shift;
  } else {
    c_id1 = c0 + i + j*nx + (k-1)*nx*ny;
    c_id2 = c0 + i + j*nx + k*nx*ny;
  }

  mb->face_cells[2*f_id]     = c_id1;
  mb->face_cells[2*f_id + 1] = c_id2;

  /* Connectiviy for z-normal faces:
   *
   *  Vtx2        Vtx3
   *  (i,j+1)     (i+1,j+1)
   *
   *   *-----------*       y (j)
   *   |           |       ^
   *   |           |       |
   *   |     *     |       |
   *   |  (i,j,k)  |       .----->x (i)
   *   |           |
   *   *-----------*
   *  Vtx1        Vtx4
   * (i,j)        (i+1,j)
   *
   */
  mb->face_vertices[4*f_id + 3] = i0 + i     + j*nxp1     + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 2] = i0 + i     + (j+1)*nxp1 + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 1] = i0 + (i+1) + (j+1)*nxp1 + k*nxp1*nyp1;
  mb->face_vertices[4*f_id + 0] = i0 + (i+1) + j*nxp1     + k*nxp1*nyp1;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of structured meshes to build.
 *
 * \returns number of structured meshes to build.
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_number_of_meshes(void)
{
  return _n_structured_meshes;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to cartesian mesh parameters structure
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \return pointer to cs_mesh_cartesian_params_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_id(const int id)
{
  cs_mesh_cartesian_params_t *retval = _get_structured_mesh_by_id(id);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get function for structured mesh based on its name.
 *
 * \param[in] name  Name of mesh
 *
 * \returns pointer to corresponding mesh parameters, or NULL if mesh
 *          does not exist.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_name_try(const char *name)
{
  cs_mesh_cartesian_params_t *retval = _get_structured_mesh_by_name_try(name);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get function for structured mesh based on its name.
 *
 * \param[in] name  Name of mesh
 *
 * \returns pointer to corresponding mesh parameters. Raises error if mesh does
 * not exist.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_name(const char *name)
{

  if (name == NULL || strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error: Empty name string.\n");

  cs_mesh_cartesian_params_t *retval = _get_structured_mesh_by_name_try(name);

  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "Error: cartesian mesh \"%s\" does not exist.\n",
              name);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*! \brief Create cartesian mesh structure
 *
 * \param[in] name  Name of mesh to create
 *
 * \returns pointer to newly created mesh parameters
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_create(const char *name)
{
  cs_mesh_cartesian_params_t *retval = _cs_mesh_cartesian_init(name, 3);

  _build_mesh_cartesian = 1;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh with a constant step in all
 *         directions
 *
 * \param[in] name    Name of mesh to create
 * \param[in] ncells  Array of size 3 containing number of cells in each
 *                    direction
 * \param[in] xyz     Array of size 6 containing min values, followed by
 *                    max values for the three directions.
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_define_simple(const char *name,
                                int         ncells[3],
                                cs_real_t   xyz[6])
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_create(name);

  assert(mp->ndir == 3);
  for (int idim = 0; idim < 3; idim++)
    mp->params[idim] =
      _cs_mesh_cartesian_create_direction(CS_MESH_CARTESIAN_CONSTANT_LAW,
                                          ncells[idim],
                                          xyz[idim],
                                          xyz[idim+3],
                                          -1.);

  return mp->id;
}

/*----------------------------------------------------------------------------*/
/*! \brief Define directions parameters based on a user input
 *
 * \param[in] mp         Pointer to mesh parameters
 * \param[in] idir       Direction index. 0->X, 1->Y, 2->Z
 * \param[in] ncells     Number of cells for the direction
 * \param[in] vtx_coord  Array of size ncells+1 containing 1D coordinate values
 *                       for vertices on the given direction
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_user(cs_mesh_cartesian_params_t *mp,
                                  int                         idir,
                                  int                         ncells,
                                  cs_real_t                   vtx_coord[])
{
  assert(mp != NULL);

  if (mp->params[idir] != NULL)
    bft_error(__FILE__, __LINE__, 0,
              "Error: %d-th component was allready defined for this mesh.\n",
              idir);

  _cs_mesh_cartesian_direction_t *dirp = NULL;
  BFT_MALLOC(dirp, 1, _cs_mesh_cartesian_direction_t);

  dirp->ncells = ncells;
  dirp->law    = CS_MESH_CARTESIAN_USER_LAW;

  BFT_MALLOC(dirp->s, ncells + 1, cs_real_t);
  for (int i = 0; i < ncells+1; i++)
    dirp->s[i] = vtx_coord[i];

  dirp->smin = vtx_coord[0];
  dirp->smax = vtx_coord[ncells];

  dirp->progression = -1.;

  mp->params[idir] = dirp;

}

/*----------------------------------------------------------------------------*/
/*! \brief Define direction parameters based on a piecewise definition. Each
 *         part follows a geometric (or uniform) sequence. To get the uniform
 *         sequence, set the amplification factor to 1 in the wanted part.
 *
 *         A direction is split in several parts. Each part contains a number
 *         of cells, its starting and ending position (stored in a compact way)
 *         inside part_coords, the amplification factor (f) between the first
 *         and last cell size of each part. Notice that if f = 1, this leads to
 *         a uniform refinement. If f > 1, (resp f < 1) this leads to a growing
 *         (resp. decreasing) geometric progression of the cell size when
 *         moving along the direction of increasing coordinates.
 *
 * \param[in] mp            Pointer to mesh parameters
 * \param[in] idir          Direction index. 0->X, 1->Y, 2->Z
 * \param[in] n_parts       Number of parts to define the direction
 * \param[in] part_coords   Position delimiting each part (size = n_parts + 1)
 * \param[in] n_part_cells  Number of cells in each part (size = n_parts)
 * \param[in] amp_factors   Amplification factor in each part (size = n_parts)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_geom_by_part(cs_mesh_cartesian_params_t *mp,
                                          int                         idir,
                                          int                         n_parts,
                                          const cs_real_t  part_coords[],
                                          const cs_lnum_t  n_part_cells[],
                                          const cs_real_t  amp_factors[])
{
  if (n_parts == 0)
    return;

  assert(part_coords != NULL && n_part_cells != NULL && amp_factors != NULL);

  /* Compute the cumulated number of cells along this direction */

  cs_lnum_t  n_tot_cells = 0;
  for (int i = 0; i < n_parts; i++)
    n_tot_cells += n_part_cells[i];

  /* There are n_cells + 1 coordinates to define */

  cs_real_t  *vtx_coord = NULL;
  BFT_MALLOC(vtx_coord, n_tot_cells + 1, cs_real_t);

  vtx_coord[0] = part_coords[0];

  cs_lnum_t  shift = 0;

  for (int i = 0; i < n_parts; i++) {

    const cs_lnum_t  _n_cells = n_part_cells[i];
    const cs_real_t  part_length = part_coords[i+1] - part_coords[i];

    cs_real_t  *_coord = vtx_coord + shift;

    if (fabs(amp_factors[i] - 1.0) < 1e-6) {

      /* Simple case: uniform refinement for this part */

      const cs_real_t  dx = part_length/_n_cells;

      for (cs_lnum_t ix = 1; ix < _n_cells; ix++)
        _coord[ix] = part_coords[i] + ix*dx;

    }
    else { /* geometric progression (or sequence) */

      const cs_real_t  common_ratio = pow(amp_factors[i], 1./(_n_cells-1));
      const cs_real_t  l0 = part_length *
        (1 - common_ratio) / (1 - pow(common_ratio, _n_cells));

      cs_real_t  coef = l0;
      for (cs_lnum_t ix = 1; ix < _n_cells; ix++) {
        _coord[ix] = _coord[ix-1] + coef;
        coef *= common_ratio;
      }

    }

    /* Ending coordinates for this part */

    _coord[_n_cells] = part_coords[i+1];

    /* Update the shift value */

    shift += _n_cells;

  } /* Loop on parts */

  /* Finally, one relies on the user-defined API to build the direction */

  cs_mesh_cartesian_define_dir_user(mp, idir, n_tot_cells, vtx_coord);

  BFT_FREE(vtx_coord);
}

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh based on a CSV file.
 *         CSV file needs to contain :
 *         (1) First line which is empty or contains a header
 *         (2) Second line containing number of vertices per direction:
 *             format is 'nx;ny;nz'
 *         (3) Third line is empty or contains a header
 *         (4) Fourth line and onwards contains vertices coordinates for each
 *             direction. Format is "X1[i];X2[i];X3[i]" for index i.
 *             If current vertex index is beyond max value for a given
 *             direction, an empty value is expected.
 *             For example, if for index 'j' the first direction is empty,
 *             format is : ';X2[j];X3[j]'
 *
 * \param[in] name           Name of new mesh
 * \param[in] csv_file_name  name of CSV file containing mesh information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_from_csv(const char  *name,
                                  const char  *csv_file_name)
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_create(name);

  const int _ndim = 3;

  /* Read CSV file */
  FILE *f = fopen(csv_file_name, "r");

  char line[128];

  int ln     = 0;
  int vtx_id = 0;

  cs_real_t *s[3] = {NULL, NULL, NULL};
  int nc[3] = {0,0,0};

  /* Read the file lines one by one */
  while (fgets(line, 128, f))
  {
    if (ln == 0 || ln == 2) {
      /* First and third lines contain header or are empty */
      ln += 1;
      continue;

    } else if (ln == 1) {
      /* Second line contains values : <nx>;<ny>;<nz> */
      sscanf(line, "%d;%d;%d", &nc[0], &nc[1], &nc[2]);

      for (int i = 0; i < _ndim; i++)
        BFT_MALLOC(s[i], nc[i], cs_real_t);

      ln += 1;
      continue;

    }
    else {
      /* Fourth line and beyond contain values for vertices coordinates */

      char *n = NULL;
      char *c = line;

      int idim = 0;
      while (true) {
        n = strchr(c, ';');
        if (n != NULL) {
          size_t l_c = strlen(c);
          size_t l_n = strlen(n);

          if (l_c > l_n) {
            char tmp[40];
            memcpy(tmp, c, l_c - l_n);
            tmp[l_c-l_n] = '\0';

            s[idim][vtx_id] = atof(tmp);
          }

          c = n + 1;
        }
        else {
          if (strlen(c) > 1 && strcmp(c, "\n") && strcmp(c, "\r\n"))
            s[idim][vtx_id] = atof(c);

          break;
        }
        idim += 1;
      }
      vtx_id += 1;
    }
  }

  for (int i = 0; i < _ndim; i++)
    cs_mesh_cartesian_define_dir_user(mp, i, nc[i]-1, s[i]);

  for (int i = 0; i < _ndim; i++)
    BFT_FREE(s[i]);

  fclose(f);
}

/*----------------------------------------------------------------------------*/
/*! \brief Define parameters for a given direction.
 *
 * \param[in] mp           Pointer to mesh parameters
 * \param[in] idim         Geometrical direction: 0->X, 1->Y, 2->Z
 * \param[in] law          1D discretization law: constant, geometric or
 *                         parabolic
 * \param[in] ncells       Number of cells for this direction
 * \param[in] smin         Min coordinate value for this direction
 * \param[in] smax         Max coordinate value for this direction
 * \param[in] progression  Progression value, only used for geometric or
 *                         parabolic laws.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_params(cs_mesh_cartesian_params_t  *mp,
                                    int                          idim,
                                    cs_mesh_cartesian_law_t      law,
                                    int                          ncells,
                                    cs_real_t                    smin,
                                    cs_real_t                    smax,
                                    cs_real_t                    progression)
{
  assert(mp != NULL);

  cs_mesh_cartesian_law_t _law = law;
  cs_real_t _p   = progression;

  /* Sanity check that min != max */
  if (cs_math_fabs(smin-smax) < 1.e-6) {
    const char *dirnames[3] = {"X", "Y", "Z"};

    bft_error(__FILE__, __LINE__, 0,
              _("Error: min and max values for direction '%s' are equal in"
                " cartesian mesh definition.\n"),
              dirnames[idim]);
  }

  /* Sanity check for progression value */
  if (law == CS_MESH_CARTESIAN_GEOMETRIC_LAW ||
      law == CS_MESH_CARTESIAN_PARABOLIC_LAW) {
    if (cs_math_fabs(progression - 1.) < 1.e-6) {
      bft_printf("Warning: \n");
      if (law == CS_MESH_CARTESIAN_GEOMETRIC_LAW)
        bft_printf("A geometric law was defined ");
      else
        bft_printf("A parabolic law was defined ");
      bft_printf("for direction #%d using a unitary progression (p=%f).\n",
                 idim+1, progression);

      bft_printf("A constant step law is set for this direction.\n");

      _law = CS_MESH_CARTESIAN_CONSTANT_LAW;
      _p   = -1.;
    }
  }

  if (mp->params[idim] != NULL) {
    bft_printf("Warning: You are modifying parameters for direction \"%d\""
               "which was allready defined.\n",
               idim);
    bft_printf_flush();
    BFT_FREE(mp->params[idim]);
  }

  assert(idim < mp->ndir);
  mp->params[idim] = _cs_mesh_cartesian_create_direction(_law,
                                                         ncells,
                                                         smin,
                                                         smax,
                                                         _p);
}

/*----------------------------------------------------------------------------*/
/*! \brief Indicate if a cartesian mesh is to be built.
 *
 * \return 1 if mesh needs to be built, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_need_build(void)
{
  int retval = _build_mesh_cartesian;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get name of structured mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns Name of the mesh
 */
/*----------------------------------------------------------------------------*/

const char *
cs_mesh_cartesian_get_name(int  id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get group class id shift of cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns shift value
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_gc_id_shift(int  id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->gc_id_shift;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set group class id shift of cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] shift Value of shift index
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_set_gc_id_shift(int  id,
                                  int  shift)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  mp->gc_id_shift = shift;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of cells of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of cells
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_cells(int  id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->n_g_cells;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of faces of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of faces
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_faces(int  id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->n_g_faces;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of vertices of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of vertices
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_vtx(int  id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->n_g_vtx;
}

/*----------------------------------------------------------------------------*/
/*! \brief Get number of cells in a given direction.
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] idim  Index of direction: 0->X, 1->Y, 2->Z
 *
 * \return Number of cells in corresponding direction (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_ncells(int  id,
                             int  idim)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  return mp->params[idim]->ncells;
}

/*----------------------------------------------------------------------------*/
/*! \brief Build unstructured connectivity needed for partitionning.
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] m     pointer to cs_mesh_t structure
 * \param[in] mb    pointer to cs_mesh_builder_t structure
 * \param[in] echo  verbosity flag
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_block_connectivity(int                 id,
                                     cs_mesh_t          *m,
                                     cs_mesh_builder_t  *mb,
                                     long                echo)
{
  CS_UNUSED(echo);

  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);

  const cs_gnum_t nx = cs_mesh_cartesian_get_ncells(id, 0);
  const cs_gnum_t ny = cs_mesh_cartesian_get_ncells(id, 1);
  const cs_gnum_t nz = cs_mesh_cartesian_get_ncells(id, 2);

  const cs_gnum_t nxp1 = nx + 1;
  const cs_gnum_t nyp1 = ny + 1;
  const cs_gnum_t nzp1 = nz + 1;

  /* Compute global values and distribution */
  const cs_gnum_t n_g_cells = mp->n_g_cells;;
  const cs_gnum_t n_g_vtx   = mp->n_g_vtx;
  const cs_gnum_t n_g_faces = mp->n_g_faces;

  cs_lnum_t n_cells = (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]);
  cs_lnum_t n_faces = (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]);
  cs_lnum_t n_vertices = (  mb->vertex_bi.gnum_range[1]
                          - mb->vertex_bi.gnum_range[0]);

  // Compute the intervals of cells/faces/vertices of this mesh which are no
  // this rank
  const cs_gnum_t mp_g_c_range[2] = {mp->n_g_cells_offset + 1,
                                     mp->n_g_cells_offset + 1 + n_g_cells};
  cs_gnum_t       _rank_c_range[2] = { 0, 0 };
  _intersect_intervals(mb->cell_bi.gnum_range, mp_g_c_range, _rank_c_range);

  const cs_gnum_t mp_g_f_range[2] = {mp->n_g_faces_offset + 1,
                                     mp->n_g_faces_offset + 1 + n_g_faces};
  cs_gnum_t       _rank_f_range[2] = { 0, 0 };
  _intersect_intervals(mb->face_bi.gnum_range, mp_g_f_range, _rank_f_range);

  const cs_gnum_t mp_g_v_range[2] = {mp->n_g_vtx_offset + 1,
                                     mp->n_g_vtx_offset + 1 + n_g_vtx};
  cs_gnum_t       _rank_v_range[2] = { 0, 0 };
  _intersect_intervals(mb->vertex_bi.gnum_range, mp_g_v_range, _rank_v_range);

  // Number of cells on this rank
  mp->n_cells_on_rank = _rank_c_range[1] - _rank_c_range[0];
  mp->n_faces_on_rank = _rank_f_range[1] - _rank_f_range[0];
  mp->n_vtx_on_rank   = _rank_v_range[1] - _rank_v_range[0];

  /* Compute local offset on rank */
  cs_gnum_t _rank_c_offset = 0;
  cs_gnum_t _rank_f_offset = 0;
  cs_gnum_t _rank_v_offset = 0;

  if (id > 0) {
    for (int i = 0; i < id; i++) {
      _rank_c_offset += _mesh_params[i]->n_cells_on_rank;
      _rank_f_offset += _mesh_params[i]->n_faces_on_rank;
      _rank_v_offset += _mesh_params[i]->n_vtx_on_rank;
    }
  }

  /* --------- */
  /* Group ids */
  /* --------- */

  if (mb->cell_gc_id == NULL)
    BFT_MALLOC(mb->cell_gc_id, n_cells, int);

  for (cs_gnum_t i = 0; i < mp->n_cells_on_rank; i++)
    mb->cell_gc_id[i + _rank_c_offset] = mp->gc_id_shift + 1;

  if (mb->face_gc_id == NULL)
    BFT_MALLOC(mb->face_gc_id, n_faces, int);

  // Default face group is 8
  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_gc_id[i + _rank_f_offset] = 1;  /* default family */

  /* number of vertices per face array */
  if (mb->face_vertices_idx == NULL) {
    BFT_MALLOC(mb->face_vertices_idx, n_faces + 1, cs_lnum_t);
    /* First value is always 0 */
    mb->face_vertices_idx[0] = 0;
  }

  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_vertices_idx[_rank_f_offset + i + 1] =
      mb->face_vertices_idx[_rank_f_offset + i] + 4;

  /* Face to cell connectivity using global numbering */
  if (mb->face_cells == NULL)
    BFT_MALLOC(mb->face_cells, 2*n_faces, cs_gnum_t);
  if (mb->face_vertices == NULL)
    BFT_MALLOC(mb->face_vertices, 4*n_faces, cs_gnum_t);

  /* Global numbering starts at 1! */

  cs_lnum_t f_id = _rank_f_offset;
  cs_gnum_t g_f_num = 1 + mp->n_g_faces_offset;

  /* We should find a better way of filtering what is built on the
     current rank, but currently ignore everything which is out of range */
  cs_gnum_t g_f_num_min = _rank_f_range[0];
  cs_gnum_t g_f_num_max = _rank_f_range[1];

  /* X normal faces */

  for (cs_gnum_t k = 0; k < nz && g_f_num < g_f_num_max; k++) {
    for (cs_gnum_t j = 0; j < ny && g_f_num < g_f_num_max; j++) {
      for (cs_gnum_t i = 0; i < nxp1 && g_f_num < g_f_num_max; i++) {

        if (g_f_num >= g_f_num_min) {
          _add_nx_face(mp, mb, f_id, nx, ny, nz, i, j, k);
          f_id += 1;
        }
        g_f_num += 1;

      }
    }
  }

  /* Y normal faces */
  for (cs_gnum_t k = 0; k < nz && g_f_num < g_f_num_max; k++) {
    for (cs_gnum_t j = 0; j < nyp1 && g_f_num < g_f_num_max; j++) {
      for (cs_gnum_t i = 0; i < nx && g_f_num < g_f_num_max; i++) {

        if (g_f_num >= g_f_num_min) {
          _add_ny_face(mp, mb, f_id, nx, ny, nz, i, j, k);
          f_id += 1;
        }
        g_f_num += 1;

      }
    }
  }

  /* Z normal faces */
  for (cs_gnum_t k = 0; k < nzp1 && g_f_num < g_f_num_max; k++) {
    for (cs_gnum_t j = 0; j < ny && g_f_num < g_f_num_max; j++) {
      for (cs_gnum_t i = 0; i < nx && g_f_num < g_f_num_max; i++) {

        if (g_f_num >= g_f_num_min) {
          _add_nz_face(mp, mb, f_id, nx, ny, nz, i, j, k);
          f_id += 1;
        }
        g_f_num += 1;

      }
    }
  }

  BFT_REALLOC(mb->vertex_coords, 3*(_rank_v_offset + n_vertices), cs_real_t);

  /* We should find a better way of filtering what is built on the
     current rank, but currently ignore everything which is out of range */
  cs_gnum_t g_v_num_min = _rank_v_range[0];
  cs_gnum_t g_v_num_max = _rank_v_range[1];

  /* Vertex coords */
  cs_lnum_t v_id = 0;

  for (cs_gnum_t k = 0; k < nzp1; k++) {
    for (cs_gnum_t j = 0; j < nyp1; j++) {
      for (cs_gnum_t i = 0; i < nxp1; i++) {
        cs_gnum_t g_v_num = 1 + mp->n_g_vtx_offset + i + j*nxp1 + k*nxp1*nyp1;

        if (g_v_num >= g_v_num_min && g_v_num < g_v_num_max) {

          /* X coord */
          cs_gnum_t ijk[3] = { i, j, k };
          for (cs_lnum_t idim = 0; idim < 3; idim++) {
            /* Constant step: xyz[idim] = xyzmin[idim] + ijk*dx[idim] */
            if (mp->params[idim]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
              mb->vertex_coords[3*(_rank_v_offset+v_id) + idim]
                = mp->params[idim]->smin + ijk[idim] * mp->params[idim]->s[0];
            }
            /* Non constant step: We allready stored the vertices in dx,
             * since dx[j+1] - dx[j] == dx of cell j */
            else {
              mb->vertex_coords[3*(_rank_v_offset+v_id) + idim]
                = mp->params[idim]->s[ijk[idim]];
            }
          }
          v_id++;
        }

      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute all global values for meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_finalize_definition(void)
{
  for (int i = 0; i < _n_structured_meshes; i++) {
    cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(i);
    cs_real_t nxyz[3] = {0.};
    for (int j = 0; j < 3; j++)
      nxyz[j] = mp->params[j]->ncells;

    mp->n_g_cells = nxyz[0] * nxyz[1] * nxyz[2];
    mp->n_g_faces = (nxyz[0] + 1) * nxyz[1] * nxyz[2]
                  + (nxyz[1] + 1) * nxyz[2] * nxyz[0]
                  + (nxyz[2] + 1) * nxyz[0] * nxyz[1];
    mp->n_g_vtx   = (nxyz[0] + 1) * (nxyz[1] + 1) * (nxyz[2] + 1);

    /* If multiple blocks, compute offset */
    if (i > 0) {
      cs_mesh_cartesian_params_t *mp_m1 = _get_structured_mesh_by_id(i-1);
      mp->n_g_cells_offset = mp_m1->n_g_cells_offset + mp_m1->n_g_cells;
      mp->n_g_faces_offset = mp_m1->n_g_faces_offset + mp_m1->n_g_faces;
      mp->n_g_vtx_offset   = mp_m1->n_g_vtx_offset   + mp_m1->n_g_vtx;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Destroy cartesian mesh parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_params_destroy(void)
{
  if (_mesh_params == NULL)
    return;

  for (int i = 0; i < _n_structured_meshes; i++) {
    for (int j = 0; j < _mesh_params[i]->ndir; j++) {
      BFT_FREE(_mesh_params[i]->params[j]->s);
      BFT_FREE(_mesh_params[i]->params[j]);
    }
    BFT_FREE(_mesh_params[i]->params);

    BFT_FREE(_mesh_params[i]);
  }
  BFT_FREE(_mesh_params);
  _n_structured_meshes = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set maximum number of cartesian blocks (by default is set to None)
 *
 * \param[in] n_blocks  maximum number of cartesian blocks which can be created
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_set_max_number_of_blocks(int  n_blocks)
{
  if (n_blocks < _n_structured_meshes)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Max number of cartesian mesh blocks was set to \"%d\""
                " using \"%s\" while \"%d\" allready exist.\n"),
              n_blocks, __func__, _n_structured_meshes);

  _n_structured_meshes_max = n_blocks;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
