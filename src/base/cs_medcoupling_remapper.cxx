/*============================================================================
 * Interpolation using MEDCoupling Remapper.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

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
 * MEDCOUPLING library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MEDCOUPLING_LOADER)

#include <MEDCoupling_version.h>

#include <MEDFileMesh.hxx>
#include <MEDCouplingUMesh.hxx>

#include <MEDFileField1TS.hxx>
#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldFloat.hxx>
#include <MEDCouplingFieldDouble.hxx>

#include <MEDCouplingRemapper.hxx>

#include <MEDLoader.hxx>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

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
#include "cs_medcoupling_remapper.hxx"

using namespace MEDCoupling;

/*----------------------------------------------------------------------------
 * ParaMEDMED field structure
 *----------------------------------------------------------------------------*/

struct _cs_medcoupling_remapper_t {

  char                     *name;
  char                     *medfile_path;
  char                    **field_names;

  char                     *interp_method;

  int                       n_fields;

  cs_medcoupling_mesh_t    *target_mesh;

  MEDCouplingUMesh         *bbox_source_mesh;

  MEDCouplingFieldDouble  **source_fields;

  MEDCouplingRemapper      *remapper;     /* MEDCoupling remapper */

};

/*============================================================================
 * Private global variables
 *============================================================================*/

static int                          _n_remappers = 0;
static cs_medcoupling_remapper_t  **_remapper = NULL;

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Read a MEDCoupling field from a MEDFile and convert it to
 *          MEDCouplingFieldDouble
 *
 * \param[in] medfile_path  path to the med file
 * \param[in] field_name    name of the field to load
 * \param[in] iteration     associated time iteration
 * \param[in] order         associated time order
 *
 * \return  pointer to the new MEDCouplingFieldDouble struct
 *
 */
/* -------------------------------------------------------------------------- */

static MEDCouplingFieldDouble *
_cs_medcoupling_read_field_real(const char *medfile_path,
                                const char *field_name,
                                int         iteration,
                                int         order)
{

  MCAuto<MEDFileAnyTypeField1TS> f(MEDFileAnyTypeField1TS::New(medfile_path,field_name,iteration,order));
  MCAuto<MEDFileMesh> mesh(MEDFileMesh::New(medfile_path,f->getMeshName()));

  /* Case 1: Field is a allready a double */
  {
    MCAuto<MEDFileField1TS> f1(MEDCoupling::DynamicCast<MEDFileAnyTypeField1TS,MEDFileField1TS>(f));
    if(f1.isNotNull()) {
      MEDCouplingFieldDouble *dble_field(f1->field(mesh));
      return dble_field;
    }
  }

  /* Case 2: Field is a float, and we convert it to double */
  {
    MCAuto<MEDFileFloatField1TS> f1(MEDCoupling::DynamicCast<MEDFileAnyTypeField1TS,MEDFileFloatField1TS>(f));
    if(f1.isNotNull()) {
      MEDCouplingFieldFloat *float_field(f1->field(mesh));
      MEDCouplingFieldDouble *dble_field = float_field->convertToDblField();
      return dble_field;
    }
  }

}

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Creates a new cs_medcoupling_remapper_t struct
 *
 * \param[in] name             name of the new remapper
 * \param[in] elt_dim          element dimension
 * \param[in] select_criteria  selection criteria for the elements
 * \param[in] medfile_path     path to the med file
 * \param[in] n_fields         number of fields to load
 * \param[in] field_names      names of the fields to load
 * \param[in] iteration        time iteration to load
 * \param[in] order            iteration order to load
 *
 * \return  pointer to the new cs_medcoupling_remapper_t struct
 *
 */
/* -------------------------------------------------------------------------- */

static cs_medcoupling_remapper_t *
_cs_paramedmem_create_remapper(const char   *name,
                               int           elt_dim,
                               const char   *select_criteria,
                               const char   *medfile_path,
                               int           n_fields,
                               const char  **field_names,
                               int           iteration,
                               int           order)
{
  cs_medcoupling_remapper_t *r = NULL;
  BFT_MALLOC(r, 1, cs_medcoupling_remapper_t);

  r->n_fields = n_fields;

  BFT_MALLOC(r->name, strlen(name)+1, char);
  strcpy(r->name, name);

  // Store fields and medfile info in case updates are needed

  BFT_MALLOC(r->medfile_path, strlen(medfile_path)+1, char);
  strcpy(r->medfile_path, medfile_path);

  BFT_MALLOC(r->field_names, n_fields, char *);
  for (int i = 0; i < n_fields; i++) {
    BFT_MALLOC(r->field_names[i], strlen(field_names[i])+1, char);
    strcpy(r->field_names[i], field_names[i]);
  }

  // New MEDCoupling UMesh linked to Code_Saturne mesh
  cs_medcoupling_mesh_t *new_mesh = cs_medcoupling_mesh_create(name,
                                                               select_criteria,
                                                               elt_dim);

  cs_mesh_t *parent_mesh = cs_glob_mesh;
  cs_medcoupling_mesh_copy_from_base(parent_mesh, new_mesh, 1);
  r->target_mesh = new_mesh;

  // MEDCoupling remapper (sequential interpolation)

  r->remapper = new MEDCouplingRemapper;
  r->remapper->setPrecision(1.0e-12);
  r->remapper->setIntersectionType(INTERP_KERNEL::Triangulation);

  // Read the fields from the medfile

  BFT_MALLOC(r->source_fields, n_fields, MEDCouplingFieldDouble *);
  for (int ii = 0; ii < n_fields; ii++) {
    r->source_fields[ii] = _cs_medcoupling_read_field_real(medfile_path,
                                                           field_names[ii],
                                                           iteration,
                                                           order);

  }

  // Set the interpolation type (P0P0 or P1P0) based on source_fields type
  BFT_MALLOC(r->interp_method, 5, char);
  if (r->source_fields[0]->getTypeOfField() == MEDCoupling::ON_CELLS) {
    r->interp_method = "P0P0";
  } else if (r->source_fields[0]->getTypeOfField() == MEDCoupling::ON_NODES) {
    r->interp_method = "P1P0";
  }

  // REduced file mesh: to improve the interpolation performance,
  //                    we use a reduced mesh based only on the cells which
  //                    are intersected by the local mesh bounding box

  r->bbox_source_mesh
    = dynamic_cast<MEDCouplingUMesh *>(r->source_fields[0]->getMesh());

  return r;
}

/*----------------------------------------------------------------------------
 * Add a new remapper.
 *
 * parameters:
 *   name            <-- new remapper name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------- */
/*!
 * \brief   Add a new remapper to the list
 *
 * \param[in] name             name of the new remapper
 * \param[in] elt_dim          element dimension
 * \param[in] select_criteria  selection criteria for the elements
 * \param[in] medfile_path     path to the med file
 * \param[in] n_fields         number of fields to load
 * \param[in] field_names      names of the fields to load
 * \param[in] iteration        time iteration to load
 * \param[in] order            iteration order to load
 *
 */
/* -------------------------------------------------------------------------- */

static void
_cs_paramedmem_add_remapper(const char   *name,
                            int           elt_dim,
                            const char   *select_criteria,
                            const char   *medfile_path,
                            int           n_fields,
                            const char  **field_names,
                            int           iteration,
                            int           order)
{
  // Allocate or reallocate if needed

  if (_remapper == NULL)
    BFT_MALLOC(_remapper, 1, cs_medcoupling_remapper_t *);
  else
    BFT_REALLOC(_remapper, _n_remappers+1, cs_medcoupling_remapper_t *);

  // Initialize new remapper, and update number of remappers

  _remapper[_n_remappers] = _cs_paramedmem_create_remapper(name,
                                                           elt_dim,
                                                           select_criteria,
                                                           medfile_path,
                                                           n_fields,
                                                           field_names,
                                                           iteration,
                                                           order);

  _n_remappers++;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Interpolate values for a given field without using the reduced bbox
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] field_id     id of the field to interpolate (in the list given before)
 * \param[in] default_val  value to apply for elements not intersected by
 *                         source mesh
 *
 * \return  pointer to cs_real_t array containing new values
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
_cs_medcoupling_remapper_copy_values_no_bbox(cs_medcoupling_remapper_t *r,
                                             int                        field_id,
                                             double                     default_val)
{

  cs_lnum_t n_elts = r->target_mesh->n_elts;

  cs_real_t *new_vals;
  BFT_MALLOC(new_vals, n_elts, cs_real_t);

  for (int i = 0; i < n_elts; i++) {
    new_vals[i] = default_val;
  }

  if (n_elts > 0) {

    MEDCouplingFieldDouble *source_field = r->source_fields[field_id];
    source_field->setNature(IntensiveMaximum);

    MEDCouplingFieldDouble *target_field
      = r->remapper->transferField(source_field, default_val);


    const double *val_ptr = target_field->getArray()->getConstPointer();

    for (int i = 0; i < n_elts; i++) {
      new_vals[i] = val_ptr[i];
    }
  }

  return new_vals;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Interpolate values for a given field using the reduced bbox
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] field_id     id of the field to interpolate (in the list given before)
 * \param[in] default_val  value to apply for elements not intersected by
 *                         source mesh
 *
 * \return  pointer to cs_real_t array containing new values
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
_cs_medcoupling_remapper_copy_values_with_bbox(cs_medcoupling_remapper_t *r,
                                               int                        field_id,
                                               double                     default_val)
{

  cs_lnum_t n_elts = r->target_mesh->n_elts;
  cs_lnum_t n_elts_loc = cs_glob_mesh->n_cells;

  cs_real_t *new_vals;
  BFT_MALLOC(new_vals, n_elts_loc, cs_real_t);
  for (int i = 0; i < n_elts_loc; i++) {
    new_vals[i] = default_val;
  }

  if (n_elts > 0) {
    // List of subcells intersecting the local mesh bounding box
    const cs_real_t *rbbox = r->target_mesh->bbox;

    const DataArrayInt *subcells
      = r->bbox_source_mesh->getCellsInBoundingBox(rbbox,
                                                     1.1);
      // Construct the subfields based on the subcells list
    MEDCouplingFieldDouble *source_field
      = r->source_fields[field_id]->buildSubPart(subcells);

    // Set the nature of the field
    source_field->setNature(IntensiveMaximum); /* TODO options */

    // Interpolate the new values
    MEDCouplingFieldDouble *target_field
      = r->remapper->transferField(source_field, default_val);

    // Generate the output array
    const double *val_ptr = target_field->getArray()->getConstPointer();
    int npts = target_field->getNumberOfValues();

    const cs_lnum_t *r_elt_list = r->target_mesh->elt_list;

    if (r_elt_list != NULL) {
      const cs_lnum_t *r_new_connec = r->target_mesh->new_to_old;

      for (int i = 0; i < npts; i++) {
        int e_id = r_new_connec[i];
        new_vals[e_id] = val_ptr[i];
      }
    } else {
      for (int i = 0; i < npts; i++) {
        new_vals[i] = val_ptr[i];
      }
    }

  }

  return new_vals;

}

/* -------------------------------------------------------------------------- */
/*!
 * \brief update the interpolation matrix without using the reduced bbox
 *
 * \param[in] r  pointer to the cs_medcoupling_remapper_t struct
 *
 */
/* -------------------------------------------------------------------------- */

void
_cs_medcoupling_remapper_setup_no_bbox(cs_medcoupling_remapper_t *r)
{

  cs_lnum_t n_elts = r->target_mesh->n_elts;

  if (n_elts > 0) {

    MEDCouplingFieldDouble *source_field = r->source_fields[0];

    source_field->setNature(IntensiveMaximum);

    r->remapper->setPrecision(1.e-12);
    r->remapper->setIntersectionType(INTERP_KERNEL::Triangulation);

    r->remapper->prepare(source_field->getMesh(),
                         r->target_mesh->med_mesh,
                         r->interp_method);
  }

}

/* -------------------------------------------------------------------------- */
/*!
 * \brief update the interpolation matrix using the reduced bbox
 *
 * \param[in] r  pointer to the cs_medcoupling_remapper_t struct
 *
 */
/* -------------------------------------------------------------------------- */

void
_cs_medcoupling_remapper_setup_with_bbox(cs_medcoupling_remapper_t  *r)
{
  cs_lnum_t n_elts = r->target_mesh->n_elts;

  if (n_elts > 0) {
    // List of subcells intersecting the local mesh bounding box
    const cs_real_t *rbbox = r->target_mesh->bbox;

    const DataArrayInt *subcells
      = r->bbox_source_mesh->getCellsInBoundingBox(rbbox,
                                                     1.1);

      // Construction of a subfield and the submesh associated with it.
    MEDCouplingFieldDouble *source_field
      = r->source_fields[0]->buildSubPart(subcells);

    source_field->setNature(IntensiveMaximum); /* TODO options */

    // Update the remapper structure and interpolation matrix
    // TODO allow settings for precision and interpolation type
    r->remapper->setPrecision(1.e-12);
    r->remapper->setIntersectionType(INTERP_KERNEL::Triangulation);

    r->remapper->prepare(source_field->getMesh(),
                         r->target_mesh->med_mesh,
                         r->interp_method);
  }
}
/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/* -------------------------------------------------------------------------- */
/*!
 * \brief get a remapper by its id
 *
 * \param[in] r_id  id of the remapper
 *
 * \return  pointer to cs_medcoupling_remapper_t struct
 */
/* -------------------------------------------------------------------------- */

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_id(int  r_id)
{
  if (r_id < _n_remappers) {
    cs_medcoupling_remapper_t *r = _remapper[r_id];
    return r;
  } else {
    return NULL;
  }
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief get a remapper by its name
 *
 * \param[in] name  name of the remapper
 *
 * \return  pointer to cs_medcoupling_remapper_t struct
 */
/* -------------------------------------------------------------------------- */

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_name_try(const char  *name)
{
  if (_n_remappers > 0) {
    for (int r_id = 0; r_id < _n_remappers; r_id++) {
      const char *r_name = _remapper[r_id]->name;
      if (strcmp(r_name, name) == 0) {
        return _remapper[r_id];

      }
    }
  }

  return NULL;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief initialize a remapper based on a set of given arguments
 *
 * \param[in] name             name of the new remapper
 * \param[in] elt_dim          element dimension
 * \param[in] select_criteria  selection criteria for the elements
 * \param[in] medfile_path     path to the med file
 * \param[in] n_fields         number of fields to load
 * \param[in] field_names      names of the fields to load
 * \param[in] iteration        time iteration to load
 * \param[in] order            iteration order to load
 *
 * \return  id of the new remapper
 */
/* -------------------------------------------------------------------------- */

int
cs_medcoupling_remapper_initialize(const char   *name,
                                   int           elt_dim,
                                   const char   *select_criteria,
                                   const char   *medfile_path,
                                   int           n_fields,
                                   const char  **field_names,
                                   int           iteration,
                                   int           order)
{
  _cs_paramedmem_add_remapper(name,
                              elt_dim,
                              select_criteria,
                              medfile_path,
                              n_fields,
                              field_names,
                              iteration,
                              order);

  int r_id = _n_remappers - 1;

  return r_id;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief set and load a given time iteration from the MED file
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] iteration    time iteration to load
 * \param[in] order        iteration order to load
 *
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_remapper_set_iteration(cs_medcoupling_remapper_t  *r,
                                      int                         iteration,
                                      int                         order)
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i] = _cs_medcoupling_read_field_real(r->medfile_path,
                                                          r->field_names[i],
                                                          iteration,
                                                          order);
  }
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief update the interpolation matrix of the remapper
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 *
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_remapper_setup(cs_medcoupling_remapper_t  *r)
{
  cs_lnum_t n_elts = r->target_mesh->n_elts;

  if (n_elts > 0) {
    // List of subcells intersecting the local mesh bounding box
    const cs_real_t *rbbox = r->target_mesh->bbox;

    if (rbbox == NULL) {
      _cs_medcoupling_remapper_setup_no_bbox(r);

    } else {
      _cs_medcoupling_remapper_setup_with_bbox(r);
    }

  }
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Interpolate values for a given field
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] field_id     id of the field to interpolate (in the list given before)
 * \param[in] default_val  value to apply for elements not intersected by
 *                         source mesh
 *
 * \return  pointer to cs_real_t array containing the new values
 *
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
cs_medcoupling_remapper_copy_values(cs_medcoupling_remapper_t  *r,
                                    int                         field_id,
                                    double                      default_val)
{

  cs_real_t *new_vals = NULL;

  if (r->target_mesh->elt_dim == 2) {
    new_vals
      = _cs_medcoupling_remapper_copy_values_no_bbox(r, field_id, default_val);
  } else if (r->target_mesh->elt_dim == 3) {
    new_vals
      = _cs_medcoupling_remapper_copy_values_with_bbox(r, field_id, default_val);
  }

  return new_vals;
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief translate the mesh using a given vector
 *
 * \param[in] r            pointer to the cs_medcoupling_remapper_t struct
 * \param[in] translation  translation vector
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_remapper_translate(cs_medcoupling_remapper_t  *r,
                                  cs_real_t                   translation[3])
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i]->getMesh()->translate(translation);
  }
}

/* -------------------------------------------------------------------------- */
/*!
 * \brief Rotate the mesh using a center point, axis and angle
 *
 * \param[in] r          pointer to the cs_medcoupling_remapper_t struct
 * \param[in] invariant  coordinates of the invariant point
 * \param[in] axis       rotation axis vector
 * \param[in] angle      rotation angle in radians
 *
 */
/* -------------------------------------------------------------------------- */

void
cs_medcoupling_remapper_rotate(cs_medcoupling_remapper_t  *r,
                               cs_real_t                   invariant[3],
                               cs_real_t                   axis[3],
                               cs_real_t                   angle)
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i]->getMesh()->rotate(invariant, axis, angle);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif // HAVE_MEDCOUPLING_LOADER
