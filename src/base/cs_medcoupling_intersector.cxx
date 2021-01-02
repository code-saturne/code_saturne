/*============================================================================
 * Interpolation using MEDCoupling Intersector.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_file.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_medcoupling_utils.hxx"
#include "cs_medcoupling_intersector.h"

/*----------------------------------------------------------------------------
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)

#include <MEDCoupling_version.h>

#include <MEDFileMesh.hxx>
#include <MEDCouplingUMesh.hxx>

#include <MEDFileField1TS.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldFloat.hxx>
#include <MEDCouplingFieldDouble.hxx>
#include <MEDFileFieldMultiTS.hxx>

#include <MEDCouplingRemapper.hxx>

#include <MEDLoader.hxx>

#include <MEDCouplingNormalizedUnstructuredMesh.txx>
#include "Interpolation3D.hxx"

using namespace MEDCoupling;
#endif

/*----------------------------------------------------------------------------
 *  Intersector structure
 *----------------------------------------------------------------------------*/

struct _cs_medcoupling_intersector_t {

  char                           *name;
  char                           *medfile_path;
  char                           *interp_method;

  cs_medcoupling_mesh_t          *local_mesh;

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)
  MEDCouplingUMesh               *source_mesh;
#else
  void                           *source_mesh;
#endif

  int                             matrix_needs_update;
  cs_real_t                      *vol_intersect;

};

/*============================================================================
 * Private global variables
 *============================================================================*/

static int                             _n_intersects = 0;
static cs_medcoupling_intersector_t  **_intersects   = NULL;

#if defined(HAVE_MEDCOUPLING) && defined(HAVE_MEDCOUPLING_LOADER)

/*----------------------------------------------------------------------------*/
/*!
 * \brief create a cs_medcoupling_intersector_t object
 *
 * \return pointer to the newly created object
 */
/*----------------------------------------------------------------------------*/

static cs_medcoupling_intersector_t *
_create_intersector(void)
{

  cs_medcoupling_intersector_t *mi = NULL;
  BFT_MALLOC(mi, 1, cs_medcoupling_intersector_t);

  mi->name                = NULL;
  mi->medfile_path        = NULL;
  mi->interp_method       = NULL;
  mi->local_mesh          = NULL;
  mi->source_mesh         = NULL;
  mi->vol_intersect       = NULL;
  mi->matrix_needs_update = 1;

  return mi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a cs_medcoupling_intersector with given parameters
 *
 * \param[in] mi               pointer to the cs_medcoupling_intersector_t struct
 * \param[in] name             name of the intersector
 * \param[in] medfile_path     path to the MED file
 * \param[in] interp_method    interpolation method (P0P0, P1P0, ..)
 * \param[in] select_criteria  selection criteria
 *
 */
/*----------------------------------------------------------------------------*/

void
_allocate_intersector(cs_medcoupling_intersector_t *mi,
                      const char                   *name,
                      const char                   *medfile_path,
                      const char                   *interp_method,
                      const char                   *select_criteria)
{

  BFT_MALLOC(mi->name, strlen(name)+1, char);
  strcpy(mi->name, name);

  BFT_MALLOC(mi->medfile_path, strlen(medfile_path)+1, char);
  strcpy(mi->medfile_path, medfile_path);

  BFT_MALLOC(mi->interp_method, strlen(interp_method)+1, char);
  strcpy(mi->interp_method, interp_method);

  mi->local_mesh = cs_medcoupling_mesh_create(name, select_criteria, 3);
  cs_medcoupling_mesh_copy_from_base(cs_glob_mesh, mi->local_mesh, 1);

  mi->matrix_needs_update = 1;

  MEDCoupling::MEDFileUMesh *mesh = MEDCoupling::MEDFileUMesh::New(medfile_path);
  mi->source_mesh = mesh->getMeshAtLevel(0);

  BFT_MALLOC(mi->vol_intersect, cs_glob_mesh->n_cells, cs_real_t);
  for (cs_lnum_t e_id = 0; e_id < cs_glob_mesh->n_cells; e_id++)
    mi->vol_intersect[e_id] = 0.;

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief destroy a given intersector
 *
 * \param[in] mi  pointer to the cs_medcoupling_intersector_t struct
 */
/*----------------------------------------------------------------------------*/

void
_destroy_intersector(cs_medcoupling_intersector_t *mi)
{

  BFT_FREE(mi->name);
  BFT_FREE(mi->medfile_path);
  BFT_FREE(mi->interp_method);
  BFT_FREE(mi->source_mesh);
  BFT_FREE(mi->vol_intersect);
  cs_medcoupling_mesh_destroy(mi->local_mesh);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute intersection matrix and update the intersection array
 *
 * \param[in] mi  pointer to the cs_medcoupling_intersector_t struct
 */
/*----------------------------------------------------------------------------*/

void
_compute_intersection_volumes(cs_medcoupling_intersector_t *mi)
{

  /* If local mesh is empty, nothing to do... */
  cs_lnum_t n_elts = mi->local_mesh->n_elts;
  if (n_elts > 0 && mi->matrix_needs_update) {

    /* initialize the pointer */
    for (cs_lnum_t c_id = 0; c_id < cs_glob_mesh->n_cells; c_id++)
      mi->vol_intersect[c_id] = 0.;

    /* Matrix for the target mesh */
    MEDCouplingNormalizedUnstructuredMesh<3,3>
      tMesh_wrapper(mi->local_mesh->med_mesh);

    /* Matrix for the source mesh, based on the bbox of the target mesh */
    const cs_real_t *bbox = mi->local_mesh->bbox;

    const DataArrayIdType *subcells
      = mi->source_mesh->getCellsInBoundingBox(bbox, 1.05);

    MEDCouplingNormalizedUnstructuredMesh<3,3>
      sMesh_wrapper(mi->source_mesh->buildPartOfMySelf(subcells->begin(),
                                                       subcells->end(),
                                                       true));

    /* Compute the intersection matrix between source and target meshes */
    std::vector<std::map<mcIdType, double> > mat;
    INTERP_KERNEL::Interpolation3D interpolator;

    interpolator.interpolateMeshes(sMesh_wrapper,
                                   tMesh_wrapper,
                                   mat,
                                   mi->interp_method);

    /* Loop on the different elements of the target mesh.
     * For each element, we sum all intersected volumes to retrieve the total
     * intersected volume per cell.
     * The iterator map contains two elements:
     * -> first  : which is the index of the intersected cell in source mesh
     * -> second : which the intersection volume
     */
    const cs_lnum_t *connec = mi->local_mesh->new_to_old;
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = connec[e_id];
      for (std::map<mcIdType, double>::iterator it = mat[e_id].begin();
           it != mat[e_id].end();
           ++it)
        mi->vol_intersect[c_id] += it->second;
    }

  }


  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief dump a medcoupling mesh
 *
 * \param[in] m         MEDCouplingUMesh to dump
 * \param[in] prefix    folder where the file is to be written
 * \param[in] filename  name of the file to write
 */
/*----------------------------------------------------------------------------*/

void
_dump_medcoupling_mesh(MEDCouplingUMesh *m,
                       const char       *prefix,
                       const char       *filename)
{

#if defined(WIN32) || defined(_WIN32)
  static const char _dir_separator = '\\';
#else
  static const char _dir_separator = '/';
#endif

  const char _medfiles[] = "medfiles";
  const char _ext[]      = ".med";

  const char *subdir = prefix;

  if (cs_glob_rank_id < 1) {

    /* Sanity check to ensure the subdirectory is defined */
    if (subdir != NULL) {
      if (strlen(subdir) == 0)
        subdir = NULL;
    }
    if (subdir == NULL)
      subdir = _medfiles;

    /* Creat the subdirectory */
    cs_file_mkdir_default(subdir);

    size_t lsdir = strlen(subdir);
    size_t lname = strlen(filename);

    size_t lext = 0;
    if (cs_file_endswith(filename, _ext) == 0)
      lext = strlen(_ext);
    char *fname = NULL;
    BFT_MALLOC(fname, lsdir + lname + lext + 2, char);

    strcpy(fname, subdir);
    fname[lsdir] = _dir_separator;
    fname[lsdir+1] = '\0';
    strcat(fname, filename);
    if (lext != 0)
      strcat(fname, _ext);
    fname[lsdir+lname+lext+1] = '\0';

    WriteUMesh(fname, m, true);
  }

  return;
}

#endif

/* ========================================================================== */

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a MEDCoupling intersector.
 *
 * \param[in] name             name of the intersector
 * \param[in] medfile_path     path to the MED file
 * \param[in] interp_method    interpolation method (P0P0, P1P0, ..)
 * \param[in] select_criteria  selection criteria
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_add(const char  *name,
                               const char  *medfile_path,
                               const char  *interp_method,
                               const char  *select_criteria)
{

#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(name);
  CS_NO_WARN_IF_UNUSED(medfile_path);
  CS_NO_WARN_IF_UNUSED(interp_method);
  CS_NO_WARN_IF_UNUSED(select_criteria);

  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  if (_n_intersects == 0)
    BFT_MALLOC(_intersects, _n_intersects + 1, cs_medcoupling_intersector_t *);
  else
    BFT_REALLOC(_intersects, _n_intersects + 1, cs_medcoupling_intersector_t *);


  cs_medcoupling_intersector_t *mi = _create_intersector();
  _allocate_intersector(mi,
                        name,
                        medfile_path,
                        interp_method,
                        select_criteria);

  _intersects[_n_intersects] = mi;

  _n_intersects++;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a given MEDCoupling intersector.
 *
 * \param[in]  mi  pointer to the cs_medcoupling_intersector_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_destroy(cs_medcoupling_intersector_t  *mi)
{

#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(mi);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  _destroy_intersector(mi);

  BFT_FREE(mi);
#endif
  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all allocated intersectors.
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_destroy_all(void)
{

#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  for (int i = 0; i < _n_intersects; i++)
    cs_medcoupling_intersector_destroy(_intersects[i]);

  BFT_FREE(_intersects);
#endif

  return;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a MEDCoupling intersector using its id.
 *
 * \param[in] id  id of the intersector
 *
 * \return pointer to the cs_medcoupling_intersector_t or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_intersector_t *
cs_medcoupling_intersector_by_id(int id)
{

  cs_medcoupling_intersector_t *mi = NULL;

  if (id > -1 && id < _n_intersects)
    mi = _intersects[id];

  return mi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an intersector by name
 *
 * \param[in] name  name of the intersector
 *
 * \return pointer to the cs_medcoupling_intersector_t or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_medcoupling_intersector_t *
cs_medcoupling_intersector_by_name(const char  *name)
{
  cs_medcoupling_intersector_t *mi = NULL;

  if (_n_intersects > 0) {
    for (int i = 0; i < _n_intersects; i++) {
      if (strcmp(name, _intersects[i]->name) == 0) {
        mi = _intersects[i];
        break;
      }
    }
  }

  return mi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the intersection volumes between the source mesh and
 * code mesh
 *
 * \param[in] mi            pointer to the cs_medcoupling_intersector_t struct
 *
 * \return a pointer to the array containing the intersected volume of each cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_medcoupling_intersect_volumes(cs_medcoupling_intersector_t  *mi)
{
  cs_real_t *retval = NULL;

#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(mi);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  /* Compute intersection */
  _compute_intersection_volumes(mi);

  /* Reset intersector matrix status */
  mi->matrix_needs_update = 0;

  /* Return intersected volumes */
  retval = mi->vol_intersect;
#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] mi           pointer to the cs_medcoupling_intersector_t struct
 * \param[in] translation  translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_translate(cs_medcoupling_intersector_t  *mi,
                                     cs_real_t             translation[3])
{
#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(mi);
  CS_NO_WARN_IF_UNUSED(translation);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  mi->source_mesh->translate(translation);
  mi->matrix_needs_update = 1;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief rotate the mesh
 *
 * \param[in] mi         pointer to the cs_medcoupling_intersector_t struct
 * \param[in] invariant  Invariant point
 * \param[in] axis       Rotation axis
 * \param[in] angle      angle (in radians)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_rotate(cs_medcoupling_intersector_t  *mi,
                                  cs_real_t                      invariant[3],
                                  cs_real_t                      axis[3],
                                  cs_real_t                      angle)
{
#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(mi);
  CS_NO_WARN_IF_UNUSED(invariant);
  CS_NO_WARN_IF_UNUSED(axis);
  CS_NO_WARN_IF_UNUSED(angle);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  mi->source_mesh->rotate(invariant, axis, angle);
  mi->matrix_needs_update = 1;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] mi      pointer to the cs_medcoupling_intersector_t struct
 * \param[in] prefix  subdir prefix
 */
/*----------------------------------------------------------------------------*/

void
cs_medcoupling_intersector_dump_mesh(cs_medcoupling_intersector_t  *mi,
                                     const char                    *prefix)
{

#if !defined(HAVE_MEDCOUPLING) || !defined(HAVE_MEDCOUPLING_LOADER)
  CS_NO_WARN_IF_UNUSED(mi);
  CS_NO_WARN_IF_UNUSED(prefix);
  bft_error(__FILE__, __LINE__, 0,
            _("Error: This function cannot be called without "
              "MEDCoupling support.\n"));
#else
  _dump_medcoupling_mesh(mi->source_mesh, prefix, mi->name);
#endif

  return;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS

