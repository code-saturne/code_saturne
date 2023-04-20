/*============================================================================
 * Postprocessing utilities based on MEDCoupling functions
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

#include <string>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_defs.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

#include "cs_medcoupling_utils.h"
#include "cs_medcoupling_mesh.hxx"

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
#include <MEDCoupling_version.h>
#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingNormalizedUnstructuredMesh.txx>
#include "Interpolation2D3D.hxx"
#include <MEDFileMesh.hxx>

using namespace MEDCoupling;
#endif

#include "cs_medcoupling_postprocess.h"

/* -------------------------------------------------------------------------- */
/* Internal structures */
/* -------------------------------------------------------------------------- */

/* Enum type to define integral/mean computations mode */
typedef enum {
  CS_MEDCPL_INT_WEIGHT_NONE,
  CS_MEDCPL_INT_WEIGHT_SCALAR,
  CS_MEDCPL_INT_WEIGHT_VECTOR,
  CS_MEDCPL_INT_WEIGHT_SCALAR_VECTOR,

  CS_MEDCPL_INT_WEIGHT_N_TYPES
} cs_medcoupling_int_weight_t;

/* Slice structure */
struct _medcoupling_slice_t {
  char *name;                  /* Name of the intersection */

  cs_medcoupling_mesh_t *csm;  /* MED representation of the local mesh */

  cs_real_t normal[3];         /* Normal vector of the surface */
  cs_real_t origin[3];         /* Origin point of the surface */

  cs_lnum_t  n_elts;           /* Number of elements considered by the selection criteria */
  cs_lnum_t *elt_ids;          /* Ids of the selected elements */
  cs_real_t *surface;          /* Intersected surface for eacu element */
  cs_real_t  total_surface;    /* Total intersected surface */
};

/*============================================================================
 * Private global variables
 *============================================================================*/

static int _n_slices = 0; /* Number of defined intersections */
static cs_medcoupling_slice_t **_slices = NULL;

/*============================================================================
 * Private functions
 *============================================================================*/

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
/* -------------------------------------------------------------------------- */
/*!
 * \brief Get a slice by name. Returns NULL if not found.
 *
 * \param[in] name  name of the searched intersection
 *
 * \return pointer to slice structure. NULL if not found.
 */
/* -------------------------------------------------------------------------- */

static inline cs_medcoupling_slice_t *
_get_slice_try(const char *name)
{
  cs_medcoupling_slice_t *retval = NULL;

  if (_n_slices > 0 && name != NULL) {
    for (int i = 0; i < _n_slices; i++) {
      if (strcmp(name, _slices[i]->name) == 0) {
        retval = _slices[i];
        break;
      }
    }
  }

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Allocate a new pointer for a slice.
 *
 * \return A newly allocated and initialised pointer.
 */
/* -------------------------------------------------------------------------- */

static inline cs_medcoupling_slice_t *
_allocate_new_slice()
{
  cs_medcoupling_slice_t *_si = NULL;
  BFT_MALLOC(_si, 1, cs_medcoupling_slice_t);

  _si->name = NULL;
  _si->csm  = NULL;

  for (int i = 0; i < 3; i++) {
    _si->normal[i] = 0.;
    _si->origin[i] = 0.;
  }

  _si->n_elts = 0;
  _si->elt_ids = NULL;

  _si->surface = NULL;
  _si->total_surface = 0.;

  return _si;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Add a new slice.
 *
 * \param[in] name               Name of the new intersection
 * \param[in] selection_criteria Selection criteria for cells to intersect
 * \param[in] origin             Origin point of the surface
 * \param[in] normal             Normal vector of the surface
 *
 * \return pointer to newly created slice
 */
/* -------------------------------------------------------------------------- */

static inline cs_medcoupling_slice_t *
_add_slice(const char      *name,
           const char      *selection_criteria,
           const cs_real_t  origin[],
           const cs_real_t  normal[])
{

  cs_medcoupling_slice_t *_si =
    _get_slice_try(name);

  if (_si != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: a slice with name \"%s\" allready exists.\n"),
              name);

  _si = _allocate_new_slice();

  BFT_MALLOC(_si->name, strlen(name) + 1, char);
  strcpy(_si->name, name);

  cs_math_3_normalize(normal, _si->normal);

  memcpy(_si->origin, origin, sizeof(cs_real_3_t));

  _si->csm = cs_medcoupling_mesh_from_base(cs_glob_mesh,
                                           name,
                                           selection_criteria,
                                           3,
                                           0);

  if (_n_slices == 0)
    BFT_MALLOC(_slices, 1, cs_medcoupling_slice_t *);
  else
    BFT_REALLOC(_slices,
                _n_slices + 1,
                cs_medcoupling_slice_t *);

  _slices[_n_slices] = _si;

  _n_slices += 1;

  return _si;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute intersected surface
 *
 * \param[in] si Slice
 * \param[in] m  Pointer to MEDCouplingUMesh representing the slice
 */
/* -------------------------------------------------------------------------- */

static inline void
_compute_slice(cs_medcoupling_slice_t *si,
               MEDCouplingUMesh       *m)
{

  /* We can only intersect the elements which are inside the selection criteria */
  si->n_elts = si->csm->n_elts;
  BFT_MALLOC(si->elt_ids, si->n_elts, cs_lnum_t);
  for (cs_lnum_t e_id = 0; e_id < si->n_elts; e_id++)
    si->elt_ids[e_id] = si->csm->new_to_old[e_id];

  BFT_MALLOC(si->surface,  si->n_elts, cs_real_t);
  cs_array_real_fill_zero(si->n_elts, si->surface);

  /* Get global arrays */
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_with_ghost = cs_glob_mesh->n_cells_with_ghosts;

  cs_lnum_t *_flag = NULL;
  BFT_MALLOC(_flag, n_cells_with_ghost, cs_lnum_t);
  memset(_flag, 0, n_cells_with_ghost * sizeof(cs_lnum_t));

  if (si->n_elts > 0) {
    MEDCouplingNormalizedUnstructuredMesh<3,3> tm_wrapper(si->csm->med_mesh);
    MEDCouplingNormalizedUnstructuredMesh<3,3> sm_wrapper(m);

    /* Compute the intersection matrix between source and target meshes */
    std::vector<std::map<mcIdType, double> > mat;
    INTERP_KERNEL::Interpolation2D3D interpolator;
    interpolator.setPrecision(1e-12);
    interpolator.interpolateMeshes(sm_wrapper, tm_wrapper, mat, "P0P0");

    /* Loop on the different elements of the target mesh.
     * For each element, we sum all intersected volumes to retrieve the total
     * intersected volume per cell.
     * The iterator map contains two elements:
     * -> first  : which is the index of the intersected cell in source mesh
     * -> second : which the intersection volume
     */
    const cs_lnum_t n_elts = si->n_elts;
    const cs_lnum_t *connec = si->csm->new_to_old;
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = connec[e_id];
      for (std::map<mcIdType, double>::iterator it = mat[e_id].begin();
           it != mat[e_id].end();
           ++it) {
        si->surface[e_id] += it->second;
        _flag[c_id] = 1; // Flag cell as intersected
      }
    }
  }

  // ----------------------------------------------------------------------
  // Sanity check for when the plane intersects with a face and counts both
  // cells (when only one should).
  // Hence we apply a "0.5" coefficient on the cell
  // ----------------------------------------------------------------------

  if (cs_glob_mesh->halo != NULL)
    cs_halo_sync_num(cs_glob_mesh->halo, CS_HALO_STANDARD, _flag);

  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;
  const cs_lnum_2_t *i_face_cells = cs_glob_mesh->i_face_cells;

  const cs_real_3_t *face_normal =
    (cs_real_3_t *)cs_glob_mesh_quantities->i_face_normal;
  const cs_real_3_t *face_cog =
    (cs_real_3_t *)cs_glob_mesh_quantities->i_face_cog;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];
    if (_flag[c_id1] == 1 && _flag[c_id2] == 1) {
      // If both cells are intersected double check :
      // - If face center is inside the plane
      // - If face normal and plane normal are the same
      cs_real_t d_face_plane = cs_math_3_distance_dot_product(si->origin,
                                                              face_cog[f_id],
                                                              si->normal);
      cs_real_t nxn[3] = {0.};
      cs_math_3_cross_product(si->normal, face_normal[f_id], nxn);
      cs_real_t norm_nxn = cs_math_3_norm(nxn);

      if (cs_math_fabs(d_face_plane) < 1.e-8 && norm_nxn < 1.e-8) {
        if (c_id1 < n_cells)
          si->surface[c_id1] *= 0.5;

        if (c_id2 < n_cells)
          si->surface[c_id2] *= 0.5;
      }
    }
  }

  BFT_FREE(_flag);

  // ----------------------------------------------------------------------
  // Compute total surface
  // ----------------------------------------------------------------------

  si->total_surface = 0.;
  for (cs_lnum_t e_id = 0; e_id < si->n_elts; e_id++)
    si->total_surface += si->surface[e_id];

  cs_parall_sum(1, CS_REAL_TYPE, &(si->total_surface));

  return;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute the local integral contribution. Parallel sum needs to be done
 *        in another call.
 *
 * \param[in] si        Slice pointer
 * \param[in] scalar    Scalar array (size on n_cells)
 * \param[in] weight_s  Scalar weight array (size n_cells)
 * \param[in] weight_v  Vector weight array (size n_cells) in cs_real_3_t *
 * \param[in] int_l     Local integral value
 * \param[in] w_l       Local integrated weight value (if normalization is used)
 *
 */
/* -------------------------------------------------------------------------- */

template <const cs_medcoupling_int_weight_t iw_type>
void
_compute_scalar_integral_l(cs_medcoupling_slice_t *si,
                           const cs_real_t        *scalar,
                           const cs_real_t        *weight_s,
                           const cs_real_3_t      *weight_v,
                           cs_real_t              *int_l,
                           cs_real_t              *w_l)
{
  assert(si != NULL);
  assert(scalar != NULL);

  cs_real_t _int_l = 0.;
  cs_real_t _w_l   = 0.;

  /* --------------------------- */
  /* Avoid warnings for template */
  /* --------------------------- */
  switch (iw_type) {
  case CS_MEDCPL_INT_WEIGHT_NONE:
    {
      CS_NO_WARN_IF_UNUSED(weight_s);
      CS_NO_WARN_IF_UNUSED(weight_v);
      break;
    }
  case CS_MEDCPL_INT_WEIGHT_SCALAR:
    {
      CS_NO_WARN_IF_UNUSED(weight_v);
      break;
    }
  case CS_MEDCPL_INT_WEIGHT_VECTOR:
    {
      CS_NO_WARN_IF_UNUSED(weight_s);
      break;
    }
  default:
    break; // We use all arguments in this case.
  }

  if (si->n_elts > 0) {
    const cs_lnum_t  n_elts  = si->n_elts;
    const cs_lnum_t *elt_ids = si->elt_ids;
    const cs_real_t *surf    = si->surface;

    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = elt_ids[e_id];

      if (iw_type == CS_MEDCPL_INT_WEIGHT_NONE) {
        _int_l += surf[e_id] * scalar[c_id];
      }
      else if (iw_type == CS_MEDCPL_INT_WEIGHT_SCALAR) {
        _int_l += surf[e_id] * weight_s[c_id] * scalar[c_id];
        _w_l += surf[e_id] * weight_s[c_id];
      }
      else if (iw_type == CS_MEDCPL_INT_WEIGHT_VECTOR) {
        cs_real_t dotp = cs_math_3_dot_product(weight_v[c_id], si->normal);
        _int_l += surf[e_id] * dotp * scalar[c_id];
        _w_l += surf[e_id] * dotp;
      }
      else if (iw_type == CS_MEDCPL_INT_WEIGHT_SCALAR_VECTOR) {
        cs_real_t dotp = cs_math_3_dot_product(weight_v[c_id], si->normal);
        _int_l += surf[e_id] * weight_s[c_id] * dotp * scalar[c_id];
        _w_l += surf[e_id] * weight_s[c_id] * dotp;
      }
    }
  }

  *int_l = _int_l;

  if (iw_type > CS_MEDCPL_INT_WEIGHT_NONE)
    *w_l   = _w_l;

  return;
}

#endif

BEGIN_C_DECLS

/*============================================================================
 * Public functions
 *============================================================================*/

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get pointer to a slice based on id
 *
 * \param[in] id index of slice
 *
 * \return pointer to slice. Raises an error if index is out of
 * bounds.
 */
/* -------------------------------------------------------------------------- */

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_id(const int id)
{
  if (id < 0 || id >= _n_slices)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Requested id \"%d\" does not exist.\n"), id);


  return _slices[id];
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get pointer to slice based on name. Returns NULL if
 * not found.
 *
 * \param[in] name  Name of the slice structure
 *
 * \return pointer to slice, NULL if not found.
 */
/* -------------------------------------------------------------------------- */

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name_try(const char *name)
{
  if (name == NULL || strcmp(name, "") == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: An empty name was provided.\n"));

  cs_medcoupling_slice_t *retval = NULL;

  for (int i = 0; i < _n_slices; i++) {
    if (strcmp(name, _slices[i]->name) == 0) {
      retval = _slices[i];
      break;
    }
  }

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get pointer to slice based on name, raises an error
 * if not found.
 *
 * \param[in] name  Name of the slice structure
 *
 * \return pointer to slice, raises error if not found.
 */
/* -------------------------------------------------------------------------- */

cs_medcoupling_slice_t *
cs_medcoupling_slice_by_name(const char *name)
{
  cs_medcoupling_slice_t *retval =
    cs_medcoupling_slice_by_name_try(name);

  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: intersection with name \"%s\" was not found.\n"),
              name);

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Add a slice based on a plane.
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] length1             Length along the first axis of the plane
 * \param[in] length2             Length along the second axis of the plane
 *
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_postprocess_add_plane_slice(const char *name,
                                           const char *selection_criteria,
                                           const cs_real_t origin[],
                                           const cs_real_t normal[],
                                           const cs_real_t length1,
                                           const cs_real_t length2)
{
#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(selection_criteria);
  CS_NO_WARN_IF_UNUSED(origin);
  CS_NO_WARN_IF_UNUSED(normal);
  CS_NO_WARN_IF_UNUSED(length1);
  CS_NO_WARN_IF_UNUSED(length2);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  MEDCouplingUMesh *umesh = cs_medcoupling_create_plane_mesh(origin,
                                                             normal,
                                                             length1,
                                                             length2);

  cs_medcoupling_slice_t *si =
    _add_slice(name, selection_criteria, origin, normal);

  _compute_slice(si, umesh);

  umesh->decrRef();
#endif
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Add a slice based on a disc
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] radius              Radius of the disc
 * \param[in] n_sectors           Number of sectors for discretization. If negative
 *                                default value (36) is used.
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_postprocess_add_disc_slice(const char *name,
                                          const char *selection_criteria,
                                          const cs_real_t origin[],
                                          const cs_real_t normal[],
                                          const cs_real_t radius,
                                          const int       n_sectors)
{
#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(selection_criteria);
  CS_NO_WARN_IF_UNUSED(origin);
  CS_NO_WARN_IF_UNUSED(normal);
  CS_NO_WARN_IF_UNUSED(radius);
  CS_NO_WARN_IF_UNUSED(n_sectors);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  MEDCouplingUMesh *umesh = cs_medcoupling_create_disc_mesh(origin,
                                                            normal,
                                                            radius,
                                                            n_sectors);

  cs_medcoupling_slice_t *si =
    _add_slice(name, selection_criteria, origin, normal);

  _compute_slice(si, umesh);

  umesh->decrRef();
#endif
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Add a slice based on an annulus
 *
 * \param[in] name                Name of the slice
 * \param[in] selection_criteria  Selection criteria for cells to intersect
 * \param[in] origin              Coordinates of origin point of slice
 * \param[in] normal              Normal vector of the slice
 * \param[in] radius1             Inner radius of the annulus (hole)
 * \param[in] radius2             Outer radius of the annulus
 * \param[in] n_sectors           Number of sectors for discretization. If negative
 *                                default value (36) is used.
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_postprocess_add_annulus_slice(const char *name,
                                             const char *selection_criteria,
                                             const cs_real_t origin[],
                                             const cs_real_t normal[],
                                             const cs_real_t radius1,
                                             const cs_real_t radius2,
                                             const int       n_sectors)
{
#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(selection_criteria);
  CS_NO_WARN_IF_UNUSED(origin);
  CS_NO_WARN_IF_UNUSED(normal);
  CS_NO_WARN_IF_UNUSED(radius1);
  CS_NO_WARN_IF_UNUSED(radius2);
  CS_NO_WARN_IF_UNUSED(n_sectors);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  MEDCouplingUMesh *umesh = cs_medcoupling_create_annulus_mesh(origin,
                                                               normal,
                                                               radius1,
                                                               radius2,
                                                               n_sectors);

  cs_medcoupling_slice_t *si =
    _add_slice(name, selection_criteria, origin, normal);

  _compute_slice(si, umesh);

  umesh->decrRef();
#endif
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get number cells that may be intersected by the slice.
 *
 * \param[in] name  Name of the slice
 *
 * \return Number of elements
 */
/* -------------------------------------------------------------------------- */

cs_lnum_t
cs_medcoupling_slice_get_n_elts(const char *name)
{
  cs_medcoupling_slice_t *si =
    cs_medcoupling_slice_by_name(name);

  return si->n_elts;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get list of ids of the elements which may be intersected.
 *
 * \param[in] name  Name of the slice
 *
 * \return Pointer to list of ids (cs_lnum_t *). Do not deallocate!
 */
/* -------------------------------------------------------------------------- */

cs_lnum_t *
cs_medcoupling_slice_get_elt_ids(const char *name)
{
  cs_medcoupling_slice_t *si =
    cs_medcoupling_slice_by_name(name);

  return si->elt_ids;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get list of intersection surfaces for each cell intersected.
 *
 * \param[in] name  Name of the slice
 *
 * \return Pointer to list of intersection surfaces (cs_real_t *)
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
cs_medcoupling_slice_get_surfaces(const char *name)
{
  cs_medcoupling_slice_t *si =
    cs_medcoupling_slice_by_name(name);

  return si->surface;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Get total intersection surface between a slice and volume mesh
 *
 * \param[in] name  Name of the slice
 *
 * \return Value of total intersection surface
 */
/* -------------------------------------------------------------------------- */

cs_real_t
cs_medcoupling_slice_get_total_surface(const char *name)
{
  cs_medcoupling_slice_t *si =
    cs_medcoupling_slice_by_name(name);

  return si->total_surface;
}


/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute integral of a scalar over a slice.
 *
 * \param[in] name    Name of the slice
 * \param[in] scalar  Array of scalar values (size n_cells)
 *
 * \return Global integrated value. A cs_parall_sum is used.
 */
/* -------------------------------------------------------------------------- */

cs_real_t
cs_medcoupling_slice_scalar_integral(const char *name,
                                     cs_real_t  *scalar)
{
  cs_real_t retval = 0.;
#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(scalar);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  cs_medcoupling_slice_t * si =
    cs_medcoupling_slice_by_name(name);

  _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_NONE>
    (si, scalar, NULL, NULL, &retval, NULL);

  cs_parall_sum(1, CS_REAL_TYPE, &retval);
#endif

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute mean value of a scalar over a slice.
 *
 * \param[in] name    Name of the slice
 * \param[in] scalar  Array of scalar values (size n_cells)
 *
 * \return Global integrated value. A cs_parall_sum is used.
 */
/* -------------------------------------------------------------------------- */

cs_real_t
cs_medcoupling_slice_scalar_mean(const char *name,
                                 cs_real_t  *scalar)
{
  cs_real_t _s   = cs_medcoupling_slice_get_total_surface(name);
  cs_real_t _int = cs_medcoupling_slice_scalar_integral(name, scalar);

  cs_real_t retval = _int / _s;

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute integral of a scalar over a slice using a scalar and/or vectorial
 *        weights. If NULL is provided for both weights, the non-weighted
 *        function is called.
 *
 * \param[in] name      Name of the slice
 * \param[in] scalar    Array of scalar values (size n_cells)
 * \param[in] weight_s  Scalar weight array (size n_cells)
 * \param[in] weight_v  Vectorial weight array (size n_cells), cs_real_3_t pointer
 *
 * \return Computed integral value over entire slice (parallel)
 */
/* -------------------------------------------------------------------------- */

cs_real_t
cs_medcoupling_slice_scalar_integral_weighted(const char  *name,
                                              cs_real_t   *scalar,
                                              cs_real_t   *weight_s,
                                              cs_real_3_t *weight_v)
{
  cs_real_t retval = 0.;

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(scalar);
  CS_NO_WARN_IF_UNUSED(weight_s);
  CS_NO_WARN_IF_UNUSED(weight_v);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  /* If no weight use simpler function */
  if (weight_v == NULL && weight_s == NULL) {
    retval = cs_medcoupling_slice_scalar_integral(name, scalar);
  }
  else {
    cs_real_t _int_l;
    cs_real_t _weight_l;

    cs_medcoupling_slice_t * si =
      cs_medcoupling_slice_by_name(name);

    if (weight_s != NULL && weight_v == NULL) {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_SCALAR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }
    else if (weight_s == NULL && weight_v != NULL) {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_VECTOR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }
    else {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_SCALAR_VECTOR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }

    cs_parall_sum(1, CS_REAL_TYPE, &_int_l);

    retval = _int_l;
  }
#endif

  return retval;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Compute mean of a scalar over a slice using a scalar and/or vectorial
 *        weights. If NULL is provided for both weights, the non-weighted
 *        function is called.
 *
 * \param[in] name      Name of the slice
 * \param[in] scalar    Array of scalar values (size n_cells)
 * \param[in] weight_s  Scalar weight array (size n_cells)
 * \param[in] weight_v  Vectorial weight array (size n_cells), cs_real_3_t pointer
 *
 * \return Computed mean value over entire slice (parallel)
 */
/* -------------------------------------------------------------------------- */

cs_real_t
cs_medcoupling_slice_scalar_mean_weighted(const char  *name,
                                          cs_real_t   *scalar,
                                          cs_real_t   *weight_s,
                                          cs_real_3_t *weight_v)
{
  cs_real_t retval = 0.;

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(scalar);
  CS_NO_WARN_IF_UNUSED(weight_s);
  CS_NO_WARN_IF_UNUSED(weight_v);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without MEDCoupling support"));
#else
  /* If no weight use simpler function */
  if (weight_v == NULL && weight_s == NULL) {
    retval = cs_medcoupling_slice_scalar_mean(name, scalar);
  }
  else {
    cs_real_t _int_l;
    cs_real_t _weight_l;

    cs_medcoupling_slice_t * si =
      cs_medcoupling_slice_by_name(name);

    if (weight_s != NULL && weight_v == NULL) {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_SCALAR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }
    else if (weight_s == NULL && weight_v != NULL) {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_VECTOR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }
    else {
      _compute_scalar_integral_l<CS_MEDCPL_INT_WEIGHT_SCALAR_VECTOR>
        (si, scalar, weight_s, weight_v, &_int_l, &_weight_l);
    }

    cs_real_t work[2] = {_int_l, _weight_l};
    cs_parall_sum(2, CS_REAL_TYPE, work);

    if (cs_math_fabs(work[1]) < 1.e-12)
      retval = work[0];
    else
      retval = work[0] / work[1];

  }
#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all slices
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_slice_destroy_all()
{

  for (int i = 0; i < _n_slices; i++) {
    cs_medcoupling_slice_t *_s = _slices[i];

    BFT_FREE(_s->name);
    BFT_FREE(_s->elt_ids);
    BFT_FREE(_s->surface);

    BFT_FREE(_s);
  }

  BFT_FREE(_slices);
  _n_slices = 0;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
