/*============================================================================
 * Cartesian mesh generation
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

  /* Number of direction, set to 3 by default */
  int                             ndir;

  /* Array of size ndir, containing parameters for each direction */
  _cs_mesh_cartesian_direction_t **params;

};

/*============================================================================
 * Private global variables
 *============================================================================*/

static int _build_mesh_cartesian = 0;

static cs_mesh_cartesian_params_t *_mesh_params = NULL;

/*============================================================================
 * Private functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Create the mesh parameters structure
 *
 * \param[in] ndir  number of directions
 *
 * \return  pointer to mesh parameters structure.
 */
/*----------------------------------------------------------------------------*/

static cs_mesh_cartesian_params_t *
_cs_mesh_cartesian_init(int ndir)
{
  if (_mesh_params != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: cartesian mesh parameters were allready defined!\n"));

  BFT_MALLOC(_mesh_params, 1, cs_mesh_cartesian_params_t);

  _mesh_params->ndir = ndir;
  BFT_MALLOC(_mesh_params->params, ndir, _cs_mesh_cartesian_direction_t *);

  return _mesh_params;
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
_add_nx_face(cs_mesh_builder_t  *mb,
             cs_lnum_t           f_id,
             cs_gnum_t           nx,
             cs_gnum_t           ny,
             cs_gnum_t           nz,
             cs_gnum_t           i,
             cs_gnum_t           j,
             cs_gnum_t           k)
{
  CS_UNUSED(nz);

  /* Global numbering starts at 1! */

  cs_lnum_t i0 = 1;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  /* Face to cell connectivity */
  cs_gnum_t c0 = i + j*nx + k*nx*ny;

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (i == 0) {
    c_id2 = c0 + 1;
    mb->face_gc_id[f_id] = 1;
  }
  else if (i == nx) {
    c_id1 = c0;
    mb->face_gc_id[f_id] = 2;
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
_add_ny_face(cs_mesh_builder_t  *mb,
             cs_lnum_t           f_id,
             cs_gnum_t           nx,
             cs_gnum_t           ny,
             cs_gnum_t           nz,
             cs_gnum_t           i,
             cs_gnum_t           j,
             cs_gnum_t           k)
{
  CS_UNUSED(ny);
  CS_UNUSED(nz);

  /* Global numbering starts at 1! */

  cs_lnum_t i0 = 1;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  /* Face to cell connectivity */

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (j == 0) {
    c_id2 = i0 + i + j*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 3;
  } else if (j == ny) {
    c_id1 = i0 + i + (j-1)*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 4;
  } else {
    c_id1 = i0 + i + (j-1)*nx + k*nx*ny;
    c_id2 = i0 + i + j*nx     + k*nx*ny;
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
_add_nz_face(cs_mesh_builder_t  *mb,
             cs_lnum_t           f_id,
             cs_gnum_t           nx,
             cs_gnum_t           ny,
             cs_gnum_t           nz,
             cs_gnum_t           i,
             cs_gnum_t           j,
             cs_gnum_t           k)
{
  /* Global numbering starts at 1! */

  cs_lnum_t i0 = 1;

  cs_gnum_t nxp1 = nx+1, nyp1 = ny+1;

  cs_gnum_t c_id1 = 0;
  cs_gnum_t c_id2 = 0;

  if (k == 0) {
    c_id2 = i0 + i + j*nx + k*nx*ny;
    mb->face_gc_id[f_id] = 5;
  } else if (k == nz) {
    c_id1 = i0 + i + j*nx + (k-1)*nx*ny;
    mb->face_gc_id[f_id] = 6;
  } else {
    c_id1 = i0 + i + j*nx + (k-1)*nx*ny;
    c_id2 = i0 + i + j*nx + k*nx*ny;
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
 * \brief Return pointer to cartesian mesh parameters structure
 *
 * \return pointer to cs_mesh_cartesian_params_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_get_params(void)
{
  return _mesh_params;
}

/*----------------------------------------------------------------------------*/
/*! \brief Create cartesian mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_create(void)
{
  _cs_mesh_cartesian_init(3);

  _build_mesh_cartesian = 1;
}

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh with a constant step in all
 *         directions
 *
 * \param[in] ncells  Array of size 3 containing number of cells in each
 *                    direction
 * \param[in] xyz     Array of size 6 containing min values, followed by
 *                    max values for the three directions.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_simple(int        ncells[3],
                                cs_real_t  xyz[6])
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_get_params();

  if (mp == NULL) {
    cs_mesh_cartesian_create();
    mp = cs_mesh_cartesian_get_params();
  }

  for (int idim = 0; idim < 3; idim++)
    mp->params[idim] =
      _cs_mesh_cartesian_create_direction(CS_MESH_CARTESIAN_CONSTANT_LAW,
                                          ncells[idim],
                                          xyz[idim],
                                          xyz[idim+3],
                                          -1.);
}

/*----------------------------------------------------------------------------*/
/*! \brief Define directions parameters based on a user input
 *
 * \param[in] idir       Direction index. 0->X, 1->Y, 2->Z
 * \param[in] ncells     Number of cells for the direction
 * \param[in] vtx_coord  Array of size ncells+1 containing 1D coordinate values
 *                       for vertices on the given direction
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_user(int       idir,
                                  int       ncells,
                                  cs_real_t vtx_coord[])
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_get_params();

  if (mp == NULL) {
    cs_mesh_cartesian_create();
    mp = cs_mesh_cartesian_get_params();
  }

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
 * \param[in] idir          Direction index. 0->X, 1->Y, 2->Z
 * \param[in] n_parts       Number of parts to define the direction
 * \param[in] part_coords   Position delimiting each part (size = n_parts + 1)
 * \param[in] n_part_cells  Number of cells in each part (size = n_parts)
 * \param[in] amp_factors   Amplification factor in each part (size = n_parts)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_geom_by_part(int                idir,
                                          int                n_parts,
                                          const cs_real_t    part_coords[],
                                          const cs_lnum_t    n_part_cells[],
                                          const cs_real_t    amp_factors[])
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

  cs_mesh_cartesian_define_dir_user(idir, n_tot_cells, vtx_coord);

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
 * \param[in] csv_file_name  name of CSV file containing mesh information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_from_csv(const char *csv_file_name)
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_get_params();

  const int _ndim = 3;

  if (mp == NULL) {
    cs_mesh_cartesian_create();
    mp = cs_mesh_cartesian_get_params();
  }

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
    cs_mesh_cartesian_define_dir_user(i, nc[i]-1, s[i]);

  for (int i = 0; i < _ndim; i++)
    BFT_FREE(s[i]);

  fclose(f);
}

/*----------------------------------------------------------------------------*/
/*! \brief Define parameters for a given direction.
 *
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
cs_mesh_cartesian_define_dir_params(int                     idim,
                                    cs_mesh_cartesian_law_t law,
                                    int                     ncells,
                                    cs_real_t               smin,
                                    cs_real_t               smax,
                                    cs_real_t               progression)
{
  cs_mesh_cartesian_params_t *mp = cs_mesh_cartesian_get_params();

  if (mp == NULL) {
    cs_mesh_cartesian_create();
    mp = cs_mesh_cartesian_get_params();
  }

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
/*! \brief Get number of cells in a given direction.
 *
 * \param[in] idim  Index of direction: 0->X, 1->Y, 2->Z
 *
 * \return Number of cells in corresponding direction (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_ncells(int idim)
{
  return _mesh_params->params[idim]->ncells;
}

/*----------------------------------------------------------------------------*/
/*! \brief Build unstructured connectivity needed for partitionning.
 *
 * \param[in] m     pointer to cs_mesh_t structure
 * \param[in] mb    pointer to cs_mesh_builder_t structure
 * \param[in] echo  verbosity flag
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_connectivity(cs_mesh_t          *m,
                               cs_mesh_builder_t  *mb,
                               long                echo)
{
  CS_UNUSED(echo);

  cs_mesh_cartesian_params_t *mp = _mesh_params;

  /* Number of cells per direction */
  cs_gnum_t nx = mp->params[0]->ncells;
  cs_gnum_t ny = mp->params[1]->ncells;
  cs_gnum_t nz = mp->params[2]->ncells;

  /* Number of vertices per direction */
  cs_gnum_t nxp1 = nx+1;
  cs_gnum_t nyp1 = ny+1;
  cs_gnum_t nzp1 = nz+1;

  /* Compute global values and distribution */
  cs_gnum_t n_g_cells = nx * ny * nz;
  cs_gnum_t n_g_vtx = nxp1 * nyp1 * nzp1;
  cs_gnum_t n_g_faces = 3*n_g_cells + nx*ny + nx*nz + ny*nz;

  mb->n_g_faces = n_g_faces;
  mb->n_g_face_connect_size = n_g_faces * 2;

  m->n_g_cells = n_g_cells;
  m->n_g_vertices = n_g_vtx;

  cs_mesh_builder_define_block_dist(mb,
                                    cs_glob_rank_id,
                                    cs_glob_n_ranks,
                                    mb->min_rank_step,
                                    0,
                                    m->n_g_cells,
                                    mb->n_g_faces,
                                    m->n_g_vertices);

  cs_lnum_t n_cells = (mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0]);
  cs_lnum_t n_faces = (mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0]);
  cs_lnum_t n_vertices = (  mb->vertex_bi.gnum_range[1]
                          - mb->vertex_bi.gnum_range[0]);

  /* Group ids */
  BFT_REALLOC(mb->cell_gc_id, n_cells, int);
  for (cs_lnum_t i = 0; i < n_cells; i++)
    mb->cell_gc_id[i] = 7;

  BFT_REALLOC(mb->face_gc_id, n_faces, int);
  for (cs_lnum_t i = 0; i < n_faces; i++)
    mb->face_gc_id[i] = 7;

  /* number of vertices per face array */
  BFT_REALLOC(mb->face_vertices_idx, n_faces+1, cs_lnum_t);

  mb->face_vertices_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_faces; i++)
    mb->face_vertices_idx[i+1] = mb->face_vertices_idx[i] + 4;

  /* Face to cell connectivity using global numbering */
  BFT_REALLOC(mb->face_cells, 2*n_faces, cs_gnum_t);
  BFT_REALLOC(mb->face_vertices, 4*n_faces, cs_gnum_t);

  /* Global numbering starts at 1! */

  cs_lnum_t f_id = 0;
  cs_gnum_t g_f_num = 1;

  /* We should find a better way of filtering what is built on the
     current rank, but currently ignore everything which is out of range */
  cs_gnum_t g_f_num_min = mb->face_bi.gnum_range[0];
  cs_gnum_t g_f_num_max = mb->face_bi.gnum_range[1];

  /* X normal faces */

  for (cs_gnum_t k = 0; k < nz && g_f_num < g_f_num_max; k++) {
    for (cs_gnum_t j = 0; j < ny && g_f_num < g_f_num_max; j++) {
      for (cs_gnum_t i = 0; i < nxp1 && g_f_num < g_f_num_max; i++) {

        if (g_f_num >= g_f_num_min) {
          _add_nx_face(mb, f_id, nx, ny, nz, i, j, k);
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
          _add_ny_face(mb, f_id, nx, ny, nz, i, j, k);
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
          _add_nz_face(mb, f_id, nx, ny, nz, i, j, k);
          f_id += 1;
        }
        g_f_num += 1;

      }
    }
  }

  BFT_REALLOC(mb->vertex_coords, n_vertices*3, cs_real_t);

  /* We should find a better way of filtering what is built on the
     current rank, but currently ignore everything which is out of range */
  cs_gnum_t g_v_num_min = mb->vertex_bi.gnum_range[0];
  cs_gnum_t g_v_num_max = mb->vertex_bi.gnum_range[1];

  /* Vertex coords */
  cs_lnum_t v_id = 0;

  for (cs_gnum_t k = 0; k < nzp1; k++) {
    for (cs_gnum_t j = 0; j < nyp1; j++) {
      for (cs_gnum_t i = 0; i < nxp1; i++) {
        cs_gnum_t g_v_num = 1 + i + j*nxp1 + k*nxp1*nyp1;

        if (g_v_num >= g_v_num_min && g_v_num < g_v_num_max) {

          /* X coord */
          cs_lnum_t ijk[3] = {i,j,k};
          for (cs_lnum_t idim = 0; idim < 3; idim++) {
            /* Constant step: xyz[idim] = xyzmin[idim] + ijk*dx[idim] */
            if (mp->params[idim]->law == CS_MESH_CARTESIAN_CONSTANT_LAW) {
              mb->vertex_coords[3*v_id + idim]
                = mp->params[idim]->smin + ijk[idim] * mp->params[idim]->s[0];
            }
            /* Non constant step: We allready stored the vertices in dx,
             * since dx[j+1] - dx[j] == dx of cell j */
            else {
              mb->vertex_coords[3*v_id + idim] = mp->params[idim]->s[ijk[idim]];
            }
          }

          v_id++;
        }

      }
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

  for (int i = 0; i < _mesh_params->ndir; i++) {
    BFT_FREE(_mesh_params->params[i]->s);
    BFT_FREE(_mesh_params->params[i]);
  }
  BFT_FREE(_mesh_params->params);

  BFT_FREE(_mesh_params);
  _mesh_params = NULL;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
