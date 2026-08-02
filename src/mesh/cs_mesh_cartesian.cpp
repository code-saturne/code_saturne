/*============================================================================
 * Cartesian mesh generation
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_block_dist.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_builder.h"
#include "mesh/cs_mesh_cartesian.h"

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
} _cs_mesh_cartesian_direction_t;

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

  /* O-grid cylinder parameters */
  int       ogrid_cylinder_mode;
  cs_real_t ogrid_r_core;
  cs_real_t ogrid_r_outer;
  int       ogrid_nr;
  cs_real_t ogrid_r_prog;

};

/*============================================================================
 * Private global variables
 *============================================================================*/

/* Flag to tell if a cartesian mesh is to be built */
static int _build_mesh_cartesian = 0;

/* Parameters for structured mesh */
static int _n_structured_meshes = 0;
static cs_mesh_cartesian_params_t **_mesh_params = nullptr;

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
    i3[0] = 0;
    i3[1] = 0;
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
 * \returns pointer to corresponding mesh parameters, or nullptr if mesh does
 *          not exist
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_cartesian_params_t *
_get_structured_mesh_by_name_try(const char *name)
{
  cs_mesh_cartesian_params_t *retval = nullptr;

  if (name != nullptr && strlen(name) > 0) {
    for (int i = 0; i < _n_structured_meshes; i++) {
      if (_mesh_params[i]->name != nullptr &&
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

  if (_new_mesh != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: a mesh with name \"%s\" allready exists.\n",
              name);

  CS_MALLOC(_new_mesh, 1, cs_mesh_cartesian_params_t);

  _new_mesh->name = nullptr;
  if (name != nullptr && strlen(name) > 0) {
    size_t _l = strlen(name);
    CS_MALLOC(_new_mesh->name, _l+1, char);
    strcpy(_new_mesh->name, name);
    _new_mesh->name[_l] = '\0';
  }

  _new_mesh->gc_id_shift = 0;
  _new_mesh->ogrid_cylinder_mode = 0;
  _new_mesh->ogrid_r_core = -1.0;
  _new_mesh->ogrid_r_outer = -1.0;
  _new_mesh->ogrid_nr = 0;
  _new_mesh->ogrid_r_prog = 1.0;

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
  CS_MALLOC(_new_mesh->params, ndir, _cs_mesh_cartesian_direction_t *);
  for (int i = 0; i < ndir; i++)
    _new_mesh->params[i] = nullptr;

  int _id = _n_structured_meshes;
  _new_mesh->id = _id;

  _n_structured_meshes += 1;
  CS_REALLOC(_mesh_params, _n_structured_meshes ,cs_mesh_cartesian_params_t *);

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
  _cs_mesh_cartesian_direction_t *dirp = nullptr;

  if (smax < smin)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: smax < smin in %s\n"), __func__);

  CS_MALLOC(dirp, 1, _cs_mesh_cartesian_direction_t);

  dirp->ncells = ncells;
  dirp->smin   = smin;
  dirp->smax   = smax;

  dirp->law = law;

  cs_real_t dir_len = smax-smin;

  if (law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
    dirp->progression = -1.;
    CS_MALLOC(dirp->s, 1, cs_real_t);

    if (dirp->ncells > 0)
      dirp->s[0] = dir_len / dirp->ncells;
    else
      dirp->s[0] = smin;
  }
  else if (law == CS_MESH_CARTESIAN_GEOMETRIC_LAW) {
    dirp->progression = progression;
    cs_real_t rho   = dirp->progression;
    cs_real_t rho_n = pow(rho, dirp->ncells);
    cs_real_t dx0   = dir_len * (rho - 1.) / (rho_n - 1.);

    CS_MALLOC(dirp->s, ncells+1, cs_real_t);

    cs_real_t dx_cur = dx0;
    dirp->s[0] = smin;
    for (int c_id = 0; c_id < ncells; c_id++) {
      dirp->s[c_id+1] = dirp->s[c_id] + dx_cur;
      dx_cur *= rho;
    }

  }
  else if (law == CS_MESH_CARTESIAN_PARABOLIC_LAW) {
    dirp->progression = progression;

    CS_MALLOC(dirp->s, ncells+1, cs_real_t);
    cs_real_t rho   = dirp->progression;

    cs_real_t dx0 = 0.;

    /* We need to disguish the case of even or odd number of cells */
    int is_even = (ncells % 2 == 0);
    int np = 0;

    if (is_even) {
      np = ncells / 2;
      cs_real_t rho_np = pow(rho, np);
      dx0 = 0.5 * dir_len * (rho - 1.) / (rho_np - 1.);
    }
    else {
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
  }
  else if (j == ny) {
    c_id1 = c0 + i + (j-1)*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 5 + mp->gc_id_shift;
  }
  else {
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
  }
  else if (k == nz) {
    c_id1 = c0 + i + j*nx + (k-1)*nx*ny;
    mb->face_gc_id[f_id] = 7 + mp->gc_id_shift;
  }
  else {
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
 * \returns pointer to corresponding mesh parameters, or nullptr if mesh
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

  if (name == nullptr || strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              "Error: Empty name string.\n");

  cs_mesh_cartesian_params_t *retval = _get_structured_mesh_by_name_try(name);

  if (retval == nullptr)
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
  assert(mp != nullptr);

  if (mp->params[idir] != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: %d-th component was allready defined for this mesh.\n",
              idir);

  _cs_mesh_cartesian_direction_t *dirp = nullptr;
  CS_MALLOC(dirp, 1, _cs_mesh_cartesian_direction_t);

  dirp->ncells = ncells;
  dirp->law    = CS_MESH_CARTESIAN_USER_LAW;

  CS_MALLOC(dirp->s, ncells + 1, cs_real_t);
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

  assert(part_coords != nullptr && n_part_cells != nullptr && amp_factors != nullptr);

  /* Compute the cumulated number of cells along this direction */

  cs_lnum_t  n_tot_cells = 0;
  for (int i = 0; i < n_parts; i++)
    n_tot_cells += n_part_cells[i];

  /* There are n_cells + 1 coordinates to define */

  cs_real_t  *vtx_coord = nullptr;
  CS_MALLOC(vtx_coord, n_tot_cells + 1, cs_real_t);

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

  CS_FREE(vtx_coord);
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

  cs_real_t *s[3] = {nullptr, nullptr, nullptr};
  int nc[3] = {0,0,0};

  /* Read the file lines one by one */
  while (fgets(line, 128, f))
  {
    if (ln == 0 || ln == 2) {
      /* First and third lines contain header or are empty */
      ln += 1;
      continue;
    }
    else if (ln == 1) {
      /* Second line contains values : <nx>;<ny>;<nz> */
      sscanf(line, "%d;%d;%d", &nc[0], &nc[1], &nc[2]);

      for (int i = 0; i < _ndim; i++)
        CS_MALLOC(s[i], nc[i], cs_real_t);

      ln += 1;
      continue;
    }
    else {
      /* Fourth line and beyond contain values for vertices coordinates */

      char *n = nullptr;
      char *c = line;

      int idim = 0;
      while (true) {
        n = strchr(c, ';');
        if (n != nullptr) {
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
    CS_FREE(s[i]);

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
  assert(mp != nullptr);

  cs_mesh_cartesian_law_t _law = law;
  cs_real_t _p   = progression;

  /* Sanity check that min != max */
  if (cs::abs(smin-smax) < 1.e-6) {
    const char *dirnames[3] = {"X", "Y", "Z"};

    bft_error(__FILE__, __LINE__, 0,
              _("Error: min and max values for direction '%s' are equal in"
                " cartesian mesh definition.\n"),
              dirnames[idim]);
  }

  /* Sanity check for progression value */
  if (law == CS_MESH_CARTESIAN_GEOMETRIC_LAW ||
      law == CS_MESH_CARTESIAN_PARABOLIC_LAW) {
    if (cs::abs(progression - 1.) < 1.e-6) {
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

  if (mp->params[idim] != nullptr) {
    bft_printf("Warning: You are modifying parameters for direction \"%d\""
               "which was allready defined.\n",
               idim);
    bft_printf_flush();
    CS_FREE(mp->params[idim]->s);
    CS_FREE(mp->params[idim]);
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
 * \brief Set parameters for O-grid cylinder mode
 *
 * \param[in] mp       Pointer to mesh parameters
 * \param[in] enable   Enable/disable flag
 * \param[in] r_core   Core radius
 * \param[in] r_outer  Outer cylinder radius
 * \param[in] nr       Number of radial layers
 * \param[in] r_prog   Radial geometric progression factor
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_set_ogrid_cylinder(cs_mesh_cartesian_params_t *mp,
                                     int                         enable,
                                     cs_real_t                   r_outer,
                                     int                         nr,
                                     cs_real_t                   r_prog)
{
  assert(mp != nullptr);

  mp->ogrid_cylinder_mode = enable;
  mp->ogrid_r_core = -1.0;  /* Computed automatically from core bounds */
  mp->ogrid_r_outer = r_outer;
  mp->ogrid_nr = nr;
  mp->ogrid_r_prog = r_prog;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get ogrid_cylinder_mode parameter
 *
 * \param[in] id  Id of the cartesian mesh
 * \return 1 if O-grid cylinder mode is enabled, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_ogrid_cylinder_mode(int id)
{
  cs_mesh_cartesian_params_t *mp = _get_structured_mesh_by_id(id);
  return mp->ogrid_cylinder_mode;
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

/*----------------------------------------------------------------------------*/
/*! \brief Helper to map 2D core vertex indices to global 2D index.
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_ogrid_g_c(cs_gnum_t nx, cs_gnum_t ny, cs_gnum_t i, cs_gnum_t j)
{
  return i + (nx + 1) * j;
}

/*----------------------------------------------------------------------------*/
/*! \brief Helper to map 2D crown vertex indices to global 2D index.
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_ogrid_g_o(cs_gnum_t nx, cs_gnum_t ny, cs_gnum_t nr, cs_gnum_t s, cs_gnum_t r)
{
  cs_gnum_t n_gamma = 2 * nx + 2 * ny;
  cs_gnum_t n_center = (nx + 1) * (ny + 1);

  if (r == 0) {
    cs_gnum_t i_val = 0, j_val = 0;
    if (s < nx) {
      i_val = s;
      j_val = 0;
    }
    else if (s < nx + ny) {
      i_val = nx;
      j_val = s - nx;
    }
    else if (s < 2 * nx + ny) {
      i_val = 2 * nx + ny - s;
      j_val = ny;
    }
    else {
      i_val = 0;
      j_val = 2 * nx + 2 * ny - s;
    }
    return _ogrid_g_c(nx, ny, i_val, j_val);
  }
  else {
    return n_center + (r - 1) * n_gamma + s;
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Helper to map 2D core cell indices to global 2D cell index.
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_ogrid_c_c(cs_gnum_t                   nx,
           [[maybe_unused]] cs_gnum_t  ny,
           cs_gnum_t                   i,
           cs_gnum_t                   j)
{
  return i + nx * j;
}

/*----------------------------------------------------------------------------*/
/*! \brief Helper to map 2D crown cell indices to global 2D cell index.
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_ogrid_c_o(cs_gnum_t nx, cs_gnum_t ny, cs_gnum_t s, cs_gnum_t r)
{
  cs_gnum_t n_gamma = 2 * nx + 2 * ny;
  return nx * ny + r * n_gamma + s;
}

/*----------------------------------------------------------------------------*/
/*! \brief Build unstructured conformed 5-block O-grid cylinder connectivity.
 */
/*----------------------------------------------------------------------------*/

static void
_cs_mesh_cartesian_block_connectivity_ogrid(cs_mesh_cartesian_params_t *mp,
                                            cs_mesh_t                  *m,
                                            cs_mesh_builder_t          *mb)
{
  CS_UNUSED(m);

  cs_gnum_t nx = mp->params[0]->ncells;
  cs_gnum_t ny = mp->params[1]->ncells;
  cs_gnum_t nz = mp->params[2]->ncells;
  cs_gnum_t nr = mp->ogrid_nr;
  cs_gnum_t n_gamma = 2 * nx + 2 * ny;
  cs_gnum_t n_center = (nx + 1) * (ny + 1);

  cs_gnum_t n_g_cells = mp->n_g_cells;
  cs_gnum_t n_g_vtx   = mp->n_g_vtx;
  cs_gnum_t n_g_faces = mp->n_g_faces;

  cs_lnum_t n_cells = (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]);
  cs_lnum_t n_faces = (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]);
  cs_lnum_t n_vertices = (mb->vertex_bi.gnum_range[1]
                          - mb->vertex_bi.gnum_range[0]);

  const cs_gnum_t mp_g_c_range[2] = {mp->n_g_cells_offset + 1,
                                     mp->n_g_cells_offset + 1 + n_g_cells};
  cs_gnum_t _rank_c_range[2] = {0, 0};
  _intersect_intervals(mb->cell_bi.gnum_range, mp_g_c_range, _rank_c_range);

  const cs_gnum_t mp_g_f_range[2] = {mp->n_g_faces_offset + 1,
                                     mp->n_g_faces_offset + 1 + n_g_faces};
  cs_gnum_t _rank_f_range[2] = {0, 0};
  _intersect_intervals(mb->face_bi.gnum_range, mp_g_f_range, _rank_f_range);

  const cs_gnum_t mp_g_v_range[2] = {mp->n_g_vtx_offset + 1,
                                     mp->n_g_vtx_offset + 1 + n_g_vtx};
  cs_gnum_t _rank_v_range[2] = {0, 0};
  _intersect_intervals(mb->vertex_bi.gnum_range, mp_g_v_range, _rank_v_range);

  mp->n_cells_on_rank = _rank_c_range[1] - _rank_c_range[0];
  mp->n_faces_on_rank = _rank_f_range[1] - _rank_f_range[0];
  mp->n_vtx_on_rank   = _rank_v_range[1] - _rank_v_range[0];

  cs_gnum_t _rank_c_offset = 0;
  cs_gnum_t _rank_f_offset = 0;
  cs_gnum_t _rank_v_offset = 0;

  if (mp->id > 0) {
    for (int i = 0; i < mp->id; i++) {
      _rank_c_offset += _mesh_params[i]->n_cells_on_rank;
      _rank_f_offset += _mesh_params[i]->n_faces_on_rank;
      _rank_v_offset += _mesh_params[i]->n_vtx_on_rank;
    }
  }

  if (mb->cell_gc_id == nullptr)
    CS_MALLOC(mb->cell_gc_id, n_cells, int);

  for (cs_gnum_t i = 0; i < mp->n_cells_on_rank; i++)
    mb->cell_gc_id[i + _rank_c_offset] = mp->gc_id_shift + 1;

  if (mb->face_gc_id == nullptr)
    CS_MALLOC(mb->face_gc_id, n_faces, int);

  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_gc_id[i + _rank_f_offset] = 1;

  if (mb->face_vertices_idx == nullptr) {
    CS_MALLOC(mb->face_vertices_idx, n_faces + 1, cs_lnum_t);
    mb->face_vertices_idx[0] = 0;
  }

  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_vertices_idx[_rank_f_offset + i + 1] =
      mb->face_vertices_idx[_rank_f_offset + i] + 4;

  if (mb->face_cells == nullptr)
    CS_MALLOC(mb->face_cells, 2*n_faces, cs_gnum_t);
  if (mb->face_vertices == nullptr)
    CS_MALLOC(mb->face_vertices, 4*n_faces, cs_gnum_t);

  cs_real_t r_core = mp->ogrid_r_core;
  cs_real_t r_outer = mp->ogrid_r_outer;
  cs_real_t q = mp->ogrid_r_prog;

  CS_REALLOC(mb->vertex_coords, 3*(_rank_v_offset + n_vertices), cs_real_t);

  cs_gnum_t g_v_num_min = _rank_v_range[0];
  cs_gnum_t g_v_num_max = _rank_v_range[1];

  cs_lnum_t v_id = 0;

  for (cs_gnum_t k = 0; k <= nz; k++) {
    cs_real_t z_orig = 0.0;
    if (mp->params[2]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
      z_orig = mp->params[2]->smin + k * mp->params[2]->s[0];
    } else {
      z_orig = mp->params[2]->s[k];
    }

    cs_gnum_t n_vtx_2d = n_center + nr * n_gamma;
    for (cs_gnum_t g_2d = 0; g_2d < n_vtx_2d; g_2d++) {
      cs_gnum_t g_v_num = 1 + mp->n_g_vtx_offset + g_2d + k * n_vtx_2d;

      if (g_v_num >= g_v_num_min && g_v_num < g_v_num_max) {
        cs_real_t x_phys = 0.0;
        cs_real_t y_phys = 0.0;

        if (g_2d < n_center) {
          cs_gnum_t i = g_2d % (nx + 1);
          cs_gnum_t j = g_2d / (nx + 1);

          if (mp->params[0]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            x_phys = mp->params[0]->smin + i * mp->params[0]->s[0];
          } else {
            x_phys = mp->params[0]->s[i];
          }

          if (mp->params[1]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            y_phys = mp->params[1]->smin + j * mp->params[1]->s[0];
          } else {
            y_phys = mp->params[1]->s[j];
          }
        } else {
          cs_gnum_t rem = g_2d - n_center;
          cs_gnum_t r = 1 + rem / n_gamma;
          cs_gnum_t s = rem % n_gamma;

          cs_real_t eta_r = 0.0;
          if (q > 1.0 + 1e-9 || q < 1.0 - 1e-9) {
            eta_r = (pow(q, r) - 1.0) / (pow(q, nr) - 1.0);
          } else {
            eta_r = (cs_real_t)r / (cs_real_t)nr;
          }

          cs_real_t x_int = 0.0, y_int = 0.0;
          cs_gnum_t i_val = 0, j_val = 0;
          if (s < nx) {
            i_val = s;
            j_val = 0;
          } else if (s < nx + ny) {
            i_val = nx;
            j_val = s - nx;
          } else if (s < 2 * nx + ny) {
            i_val = 2 * nx + ny - s;
            j_val = ny;
          } else {
            i_val = 0;
            j_val = 2 * nx + 2 * ny - s;
          }

          if (mp->params[0]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            x_int = mp->params[0]->smin + i_val * mp->params[0]->s[0];
          } else {
            x_int = mp->params[0]->s[i_val];
          }

          if (mp->params[1]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            y_int = mp->params[1]->smin + j_val * mp->params[1]->s[0];
          } else {
            y_int = mp->params[1]->s[j_val];
          }

          cs_real_t theta = 0.0;

          if (s < nx) {
            theta = -3.0 * cs_math_pi / 4.0 + (cs_real_t)s / (cs_real_t)nx * cs_math_pi / 2.0;
          }
          else if (s < nx + ny) {
            cs_gnum_t p = s - nx;
            theta = -cs_math_pi / 4.0 + (cs_real_t)p / (cs_real_t)ny * cs_math_pi / 2.0;
          }
          else if (s < 2 * nx + ny) {
            cs_gnum_t p = s - (nx + ny);
            theta = cs_math_pi / 4.0 + (cs_real_t)p / (cs_real_t)nx * cs_math_pi / 2.0;
          }
          else {
            cs_gnum_t p = s - (2 * nx + ny);
            theta = 3.0 * cs_math_pi / 4.0 + (cs_real_t)p / (cs_real_t)ny * cs_math_pi / 2.0;
          }

          cs_real_t x_ext = r_outer * cos(theta);
          cs_real_t y_ext = r_outer * sin(theta);

          x_phys = (1.0 - eta_r) * x_int + eta_r * x_ext;
          y_phys = (1.0 - eta_r) * y_int + eta_r * y_ext;
        }

        mb->vertex_coords[3 * (_rank_v_offset + v_id) + 0] = x_phys;
        mb->vertex_coords[3 * (_rank_v_offset + v_id) + 1] = y_phys;
        mb->vertex_coords[3 * (_rank_v_offset + v_id) + 2] = z_orig;
        v_id++;
      }
    }
  }

  cs_gnum_t g_f_num_min = _rank_f_range[0];
  cs_gnum_t g_f_num_max = _rank_f_range[1];
  cs_lnum_t f_id = _rank_f_offset;
  cs_gnum_t g_f_num = 1 + mp->n_g_faces_offset;

  cs_gnum_t i0 = 1 + mp->n_g_vtx_offset;
  cs_gnum_t c0 = 1 + mp->n_g_cells_offset;
  cs_gnum_t n_vtx_2d = n_center + nr * n_gamma;
  cs_gnum_t n_cells_2d = nx * ny + nr * n_gamma;

  #define G_V(g2d, l) (i0 + (g2d) + (l) * n_vtx_2d)
  #define G_C(g2d, l) (c0 + (g2d) + (l) * n_cells_2d)

  for (cs_gnum_t k = 0; k < nz; k++) {
    /* Lateral/Vertical faces: X-normal core faces */
    for (cs_gnum_t j = 0; j < ny; j++) {
      for (cs_gnum_t i = 0; i <= nx; i++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          mb->face_vertices[4 * f_id + 0] = G_V(_ogrid_g_c(nx, ny, i, j+1), k);
          mb->face_vertices[4 * f_id + 1] = G_V(_ogrid_g_c(nx, ny, i, j+1), k+1);
          mb->face_vertices[4 * f_id + 2] = G_V(_ogrid_g_c(nx, ny, i, j), k+1);
          mb->face_vertices[4 * f_id + 3] = G_V(_ogrid_g_c(nx, ny, i, j), k);

          cs_gnum_t c_id1 = 0;
          cs_gnum_t c_id2 = 0;

          if (i == 0) {
            cs_gnum_t s_val = 2 * nx + 2 * ny - 1 - j;
            c_id1 = G_C(_ogrid_c_o(nx, ny, s_val, 0), k);
            c_id2 = G_C(_ogrid_c_c(nx, ny, 0, j), k);
          }
          else if (i == nx) {
            cs_gnum_t s_val = nx + j;
            c_id1 = G_C(_ogrid_c_c(nx, ny, nx-1, j), k);
            c_id2 = G_C(_ogrid_c_o(nx, ny, s_val, 0), k);
          }
          else {
            c_id1 = G_C(_ogrid_c_c(nx, ny, i-1, j), k);
            c_id2 = G_C(_ogrid_c_c(nx, ny, i, j), k);
          }

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;

          if (i == 0 || i == nx)
            mb->face_gc_id[f_id] = 1 + mp->gc_id_shift;
          else
            mb->face_gc_id[f_id] = 1;

          f_id++;
        }
        g_f_num++;
      }
    }

    /* Lateral/Vertical faces: Y-normal core faces */
    for (cs_gnum_t j = 0; j <= ny; j++) {
      for (cs_gnum_t i = 0; i < nx; i++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          mb->face_vertices[4 * f_id + 0] = G_V(_ogrid_g_c(nx, ny, i, j), k+1);
          mb->face_vertices[4 * f_id + 1] = G_V(_ogrid_g_c(nx, ny, i+1, j), k+1);
          mb->face_vertices[4 * f_id + 2] = G_V(_ogrid_g_c(nx, ny, i+1, j), k);
          mb->face_vertices[4 * f_id + 3] = G_V(_ogrid_g_c(nx, ny, i, j), k);

          cs_gnum_t c_id1 = 0;
          cs_gnum_t c_id2 = 0;

          if (j == 0) {
            cs_gnum_t s_val = i;
            c_id1 = G_C(_ogrid_c_o(nx, ny, s_val, 0), k);
            c_id2 = G_C(_ogrid_c_c(nx, ny, i, 0), k);
          }
          else if (j == ny) {
            cs_gnum_t s_val = 2 * nx + ny - 1 - i;
            c_id1 = G_C(_ogrid_c_c(nx, ny, i, ny-1), k);
            c_id2 = G_C(_ogrid_c_o(nx, ny, s_val, 0), k);
          }
          else {
            c_id1 = G_C(_ogrid_c_c(nx, ny, i, j-1), k);
            c_id2 = G_C(_ogrid_c_c(nx, ny, i, j), k);
          }

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;

          if (j == 0 || j == ny)
            mb->face_gc_id[f_id] = 1 + mp->gc_id_shift;
          else
            mb->face_gc_id[f_id] = 1;

          f_id++;
        }
        g_f_num++;
      }
    }

    /* Lateral/Vertical faces: Radial-normal/circumferential outer faces */
    for (cs_gnum_t r = 1; r <= nr; r++) {
      for (cs_gnum_t s = 0; s < n_gamma; s++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          cs_gnum_t sp = (s + 1) % n_gamma;

          mb->face_vertices[4 * f_id + 0] =
            G_V(_ogrid_g_o(nx, ny, nr, s, r), k);
          mb->face_vertices[4 * f_id + 1] =
            G_V(_ogrid_g_o(nx, ny, nr, sp, r), k);
          mb->face_vertices[4 * f_id + 2] =
            G_V(_ogrid_g_o(nx, ny, nr, sp, r), k+1);
          mb->face_vertices[4 * f_id + 3] =
            G_V(_ogrid_g_o(nx, ny, nr, s, r), k+1);

          cs_gnum_t c_id1 = 0;
          cs_gnum_t c_id2 = 0;

          if (r == nr) {
            c_id1 = G_C(_ogrid_c_o(nx, ny, s, nr-1), k);
            mb->face_gc_id[f_id] = 2 + mp->gc_id_shift;
          } else {
            c_id1 = G_C(_ogrid_c_o(nx, ny, s, r-1), k);
            c_id2 = G_C(_ogrid_c_o(nx, ny, s, r), k);
            mb->face_gc_id[f_id] = 1;
          }

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;

          f_id++;
        }
        g_f_num++;
      }
    }

    /* Lateral/Vertical faces: Radial/transverse outer faces */
    for (cs_gnum_t r = 0; r < nr; r++) {
      for (cs_gnum_t s = 0; s < n_gamma; s++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          mb->face_vertices[4 * f_id + 0] = G_V(_ogrid_g_o(nx, ny, nr, s, r+1), k);
          mb->face_vertices[4 * f_id + 1] = G_V(_ogrid_g_o(nx, ny, nr, s, r+1), k+1);
          mb->face_vertices[4 * f_id + 2] = G_V(_ogrid_g_o(nx, ny, nr, s, r), k+1);
          mb->face_vertices[4 * f_id + 3] = G_V(_ogrid_g_o(nx, ny, nr, s, r), k);

          cs_gnum_t sm = (s - 1 + n_gamma) % n_gamma;
          cs_gnum_t c_id1 = G_C(_ogrid_c_o(nx, ny, sm, r), k);
          cs_gnum_t c_id2 = G_C(_ogrid_c_o(nx, ny, s, r), k);

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;
          mb->face_gc_id[f_id] = 1;

          f_id++;
        }
        g_f_num++;
      }
    }
  }

  /* Z-normal/Horizontal faces loop */
  for (cs_gnum_t k = 0; k <= nz; k++) {
    /* Horizontal faces of Core block */
    for (cs_gnum_t j = 0; j < ny; j++) {
      for (cs_gnum_t i = 0; i < nx; i++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          mb->face_vertices[4 * f_id + 0] = G_V(_ogrid_g_c(nx, ny, i+1, j), k);
          mb->face_vertices[4 * f_id + 1] = G_V(_ogrid_g_c(nx, ny, i+1, j+1), k);
          mb->face_vertices[4 * f_id + 2] = G_V(_ogrid_g_c(nx, ny, i, j+1), k);
          mb->face_vertices[4 * f_id + 3] = G_V(_ogrid_g_c(nx, ny, i, j), k);

          cs_gnum_t c_id1 = 0;
          cs_gnum_t c_id2 = 0;

          if (k == 0) {
            c_id2 = G_C(_ogrid_c_c(nx, ny, i, j), 0);
            mb->face_gc_id[f_id] = 4 + mp->gc_id_shift;
          } else if (k == nz) {
            c_id1 = G_C(_ogrid_c_c(nx, ny, i, j), nz-1);
            mb->face_gc_id[f_id] = 5 + mp->gc_id_shift;
          } else {
            c_id1 = G_C(_ogrid_c_c(nx, ny, i, j), k-1);
            c_id2 = G_C(_ogrid_c_c(nx, ny, i, j), k);
            mb->face_gc_id[f_id] = 1;
          }

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;

          f_id++;
        }
        g_f_num++;
      }
    }

    /* Horizontal faces of Outer blocks */
    for (cs_gnum_t r = 0; r < nr; r++) {
      for (cs_gnum_t s = 0; s < n_gamma; s++) {
        if (g_f_num >= g_f_num_min && g_f_num < g_f_num_max) {
          cs_gnum_t sp = (s + 1) % n_gamma;

          mb->face_vertices[4 * f_id + 0] =
            G_V(_ogrid_g_o(nx, ny, nr, s, r), k);
          mb->face_vertices[4 * f_id + 1] =
            G_V(_ogrid_g_o(nx, ny, nr, s, r+1), k);
          mb->face_vertices[4 * f_id + 2] =
            G_V(_ogrid_g_o(nx, ny, nr, sp, r+1), k);
          mb->face_vertices[4 * f_id + 3] =
            G_V(_ogrid_g_o(nx, ny, nr, sp, r), k);

          cs_gnum_t c_id1 = 0;
          cs_gnum_t c_id2 = 0;

          if (k == 0) {
            c_id2 = G_C(_ogrid_c_o(nx, ny, s, r), 0);
            mb->face_gc_id[f_id] = 4 + mp->gc_id_shift;
          } else if (k == nz) {
            c_id1 = G_C(_ogrid_c_o(nx, ny, s, r), nz-1);
            mb->face_gc_id[f_id] = 5 + mp->gc_id_shift;
          } else {
            c_id1 = G_C(_ogrid_c_o(nx, ny, s, r), k-1);
            c_id2 = G_C(_ogrid_c_o(nx, ny, s, r), k);
            mb->face_gc_id[f_id] = 1;
          }

          mb->face_cells[2 * f_id] = c_id1;
          mb->face_cells[2 * f_id + 1] = c_id2;

          f_id++;
        }
        g_f_num++;
      }
    }
  }

  #undef G_V
  #undef G_C
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

  if (mp->ogrid_cylinder_mode) {
    _cs_mesh_cartesian_block_connectivity_ogrid(mp, m, mb);
    return;
  }

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

  if (mb->cell_gc_id == nullptr)
    CS_MALLOC(mb->cell_gc_id, n_cells, int);

  for (cs_gnum_t i = 0; i < mp->n_cells_on_rank; i++)
    mb->cell_gc_id[i + _rank_c_offset] = mp->gc_id_shift + 1;

  if (mb->face_gc_id == nullptr)
    CS_MALLOC(mb->face_gc_id, n_faces, int);

  // Default face group is 8
  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_gc_id[i + _rank_f_offset] = 1;  /* default family */

  /* number of vertices per face array */
  if (mb->face_vertices_idx == nullptr) {
    CS_MALLOC(mb->face_vertices_idx, n_faces + 1, cs_lnum_t);
    /* First value is always 0 */
    mb->face_vertices_idx[0] = 0;
  }

  for (cs_gnum_t i = 0; i < mp->n_faces_on_rank; i++)
    mb->face_vertices_idx[_rank_f_offset + i + 1] =
      mb->face_vertices_idx[_rank_f_offset + i] + 4;

  /* Face to cell connectivity using global numbering */
  if (mb->face_cells == nullptr)
    CS_MALLOC(mb->face_cells, 2*n_faces, cs_gnum_t);
  if (mb->face_vertices == nullptr)
    CS_MALLOC(mb->face_vertices, 4*n_faces, cs_gnum_t);

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

  CS_REALLOC(mb->vertex_coords, 3*(_rank_v_offset + n_vertices), cs_real_t);

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
          cs_real_t x_orig = 0.0;
          cs_real_t y_orig = 0.0;
          cs_real_t z_orig = 0.0;

          if (mp->params[0]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            x_orig = mp->params[0]->smin + i * mp->params[0]->s[0];
          } else {
            x_orig = mp->params[0]->s[i];
          }

          if (mp->params[1]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            y_orig = mp->params[1]->smin + j * mp->params[1]->s[0];
          } else {
            y_orig = mp->params[1]->s[j];
          }

          if (mp->params[2]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
            z_orig = mp->params[2]->smin + k * mp->params[2]->s[0];
          } else {
            z_orig = mp->params[2]->s[k];
          }

          if (mp->ogrid_cylinder_mode) {
            cs_real_t s_min_x = mp->params[0]->smin;
            cs_real_t s_max_x = mp->params[0]->smax;
            cs_real_t s_min_y = mp->params[1]->smin;
            cs_real_t s_max_y = mp->params[1]->smax;

            cs_real_t x_mid = 0.5 * (s_min_x + s_max_x);
            cs_real_t y_mid = 0.5 * (s_min_y + s_max_y);
            cs_real_t Lx = 0.5 * (s_max_x - s_min_x);
            cs_real_t Ly = 0.5 * (s_max_y - s_min_y);

            cs_real_t x_norm = (x_orig - x_mid) / Lx;
            cs_real_t y_norm = (y_orig - y_mid) / Ly;

            cs_real_t d = fabs(x_norm) > fabs(y_norm) ?
                          fabs(x_norm) : fabs(y_norm);
            cs_real_t a = 0.5;
            if (mp->ogrid_r_core > 0.0 && mp->ogrid_r_outer > 0.0) {
              a = mp->ogrid_r_core / mp->ogrid_r_outer;
            }

            cs_real_t x_mapped = 0.0;
            cs_real_t y_mapped = 0.0;

            if (d < 1e-12) {
              x_mapped = 0.0;
              y_mapped = 0.0;
            } else if (d <= a) {
              x_mapped = x_norm;
              y_mapped = y_norm;
            } else {
              cs_real_t t = (d - a) / (1.0 - a);
              cs_real_t u = 0.0;

              /* Top sector */
              if (y_norm == d) {
                u = x_norm / y_norm;
                x_mapped = (1.0 - t) * a * u + t * sin(u * cs_math_pi / 4.0);
                y_mapped = (1.0 - t) * a + t * cos(u * cs_math_pi / 4.0);
              }
              /* Bottom sector */
              else if (-y_norm == d) {
                u = -x_norm / y_norm;
                x_mapped = (1.0 - t) * a * u + t * sin(u * cs_math_pi / 4.0);
                y_mapped = -((1.0 - t) * a + t * cos(u * cs_math_pi / 4.0));
              }
              /* Right sector */
              else if (x_norm == d) {
                u = y_norm / x_norm;
                x_mapped = (1.0 - t) * a + t * cos(u * cs_math_pi / 4.0);
                y_mapped = (1.0 - t) * a * u + t * sin(u * cs_math_pi / 4.0);
              }
              /* Left sector */
              else {
                u = -y_norm / x_norm;
                x_mapped = -((1.0 - t) * a + t * cos(u * cs_math_pi / 4.0));
                y_mapped = (1.0 - t) * a * u + t * sin(u * cs_math_pi / 4.0);
              }
            }

            x_orig = x_mid + x_mapped * Lx;
            y_orig = y_mid + y_mapped * Ly;
          }

          mb->vertex_coords[3 * (_rank_v_offset + v_id) + 0] = x_orig;
          mb->vertex_coords[3 * (_rank_v_offset + v_id) + 1] = y_orig;
          mb->vertex_coords[3 * (_rank_v_offset + v_id) + 2] = z_orig;

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

    if (mp->ogrid_cylinder_mode) {
      cs_gnum_t nx = nxyz[0];
      cs_gnum_t ny = nxyz[1];
      cs_gnum_t nz = nxyz[2];

      cs_real_t r_core = mp->ogrid_r_core;
      if (r_core <= 0.0) {
        r_core = 0.5 * (mp->params[0]->smax - mp->params[0]->smin);
        mp->ogrid_r_core = r_core;
      }

      cs_real_t r_outer = mp->ogrid_r_outer;
      cs_real_t L_r = r_outer - r_core;

      cs_real_t dx = nx > 0 ?
                     (mp->params[0]->smax - mp->params[0]->smin) / nx : 1.0;
      cs_real_t dy = ny > 0 ?
                     (mp->params[1]->smax - mp->params[1]->smin) / ny : 1.0;
      cs_real_t h_core = dx < dy ? dx : dy;
      if (h_core < 1e-12) {
        h_core = 1.0;
      }

      cs_gnum_t nr = mp->ogrid_nr;
      cs_real_t q = mp->ogrid_r_prog;

      /* Compute automatically if nr is 0 */
      if (nr == 0) {
        if (cs::abs(q - 1.0) > 1e-9) {
          /* If q < 1.0, the infinite geometric series sum is h_core / (1 - q).
             If this sum is smaller than L_r, there is no real solution.
             To prevent this and keep the starting cell size exactly equal
             to h_core, we automatically adjust q to a valid decreasing
             progression factor. */
          if (q < 1.0 && q <= 1.0 - h_core / L_r) {
            q = 1.0 - 0.8 * h_core / L_r;
            mp->ogrid_r_prog = q;
          }

          cs_real_t arg = 1.0 + L_r * (q - 1.0) / h_core;
          if (arg > 1e-9) {
            nr = (cs_gnum_t)round(log(arg) / log(q));
          } else {
            nr = 1;
          }
        } else {
          nr = (cs_gnum_t)round(L_r / h_core);
        }
        if (nr < 1) {
          nr = 1;
        }
        mp->ogrid_nr = nr;
      }

      cs_gnum_t n_gamma = 2 * nx + 2 * ny;
      cs_gnum_t n_center = (nx + 1) * (ny + 1);

      mp->n_g_cells = nz * (nx * ny + nr * n_gamma);
      mp->n_g_vtx = (nz + 1) * (n_center + nr * n_gamma);

      cs_gnum_t e_2d = 2 * nx * ny + 2 * nr * n_gamma + n_gamma / 2;
      mp->n_g_faces = nz * e_2d + (nz + 1) * (nx * ny + nr * n_gamma);
    }

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
  if (_mesh_params == nullptr)
    return;

  for (int i = 0; i < _n_structured_meshes; i++) {
    for (int j = 0; j < _mesh_params[i]->ndir; j++) {
      CS_FREE(_mesh_params[i]->params[j]->s);
      CS_FREE(_mesh_params[i]->params[j]);
    }
    CS_FREE(_mesh_params[i]->params);

    CS_FREE(_mesh_params[i]);
  }
  CS_FREE(_mesh_params);
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
