/*============================================================================
 * Interpolation using MEDCoupling Remapper.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Read a MEDCoupling field (float or double) from a MEDFile and convert it to
 * MEDCouplingFieldDouble * object
 *
 * parameters:
 *   medfile_path   <-- path to med file
 *   field_name     <-- field name
 *   iteration      <-- associated iteration
 *   order          <-- associated iteration order
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Create a new cs_medcoupling_remapper_t * object.
 *
 * parameters:
 *   name            <-- new object name
 *   elt_dim         <-- element dimension
 *   select_criteria <-- selection criteria
 *   medfile_path    <-- path of associated MED file
 *   n_fields        <-- number of fields
 *   field_names     <-- associated field names
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *
 * return:
 *   new remapper object
 *----------------------------------------------------------------------------*/

static cs_medcoupling_remapper_t *
_cs_paramedmem_create_remapper(const char   *name,
                               int           elt_dim,
                               const char   *select_criteria,
                               const char   *medfile_path,
                               int           n_fields,
                               const char  **field_names,
                               int           iteration,
                               int           iteration_order)
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

  cs_medcoupling_mesh_copy_from_base(parent_mesh, new_mesh);

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
                                                           iteration_order);

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

static void
_cs_paramedmem_add_remapper(const char   *name,
                            int           elt_dim,
                            const char   *select_criteria,
                            const char   *medfile_path,
                            int           n_fields,
                            const char  **field_names,
                            int           iteration,
                            int           iteration_order)
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
                                                           iteration_order);

  _n_remappers++;
}

/*----------------------------------------------------------------------------
 * Copy interpolated values to a new array when no bbox is available.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   field_id        <-- id of given field
 *   r               <-- pointer to remapper object
 *   default_val     <-- default value
 *
 * return:
 *   pointer to allocated values array
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Copy interpolated values to a new array when a bbox is available.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   field_id        <-- id of given field
 *   r               <-- pointer to remapper object
 *   default_val     <-- default value
 *
 * return:
 *   pointer to allocated values array
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Internal function: Creating the interpolation matrix when no bbox is available.
 *
 * This step is separated from the interpolation step since it only needs
 * to be done once per mesh, while interpolation can be done for several
 * fields.
 *
 * parameters:
 *   r               <-- remapper object
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Internal function: Creating the interpolation matrix when a bbox is available.
 *
 * This step is separated from the interpolation step since it only needs
 * to be done once per mesh, while interpolation can be done for several
 * fields.
 *
 * parameters:
 *   r               <-- remapper object
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Return remapper associated with a given id
 *
 * parameters:
 *   id <-- remapper id
 *
 * return:
 *   pointer to remapper
 *----------------------------------------------------------------------------*/

cs_medcoupling_remapper_t *
cs_medcoupling_remapper_by_id(int  r_id)
{
  cs_medcoupling_remapper_t *r = _remapper[r_id];

  return r;
}

/*----------------------------------------------------------------------------
 * Return remapper associated with a given name
 *
 * parameters:
 *   name <-- remapper name
 *
 * return:
 *   pointer to remapper, or NULL
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Create or update update the list of remappers in the case where
 * several remappers may be needed.
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
 *
 * return:
 *   id of the newly added remapper within the list
 *----------------------------------------------------------------------------*/

int
cs_medcoupling_remapper_initialize(const char   *name,
                                   int           elt_dim,
                                   const char   *select_criteria,
                                   const char   *medfile_path,
                                   int           n_fields,
                                   const char  **field_names,
                                   int           iteration,
                                   int           iteration_order)
{
  _cs_paramedmem_add_remapper(name,
                              elt_dim,
                              select_criteria,
                              medfile_path,
                              n_fields,
                              field_names,
                              iteration,
                              iteration_order);

  int r_id = _n_remappers - 1;

  return r_id;
}

/*----------------------------------------------------------------------------
 * Update field values (if several time steps are available in the MED file).
 *
 * parameters:
 *   r               <-- remapper object
 *   iteration       <-- associated iteration
 *   iteration_order <-- associated iteration order
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_set_iteration(cs_medcoupling_remapper_t  *r,
                                      int                         iteration,
                                      int                         iteration_order)
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i] = _cs_medcoupling_read_field_real(r->medfile_path,
                                                          r->field_names[i],
                                                          iteration,
                                                          iteration_order);
  }
}

/*----------------------------------------------------------------------------
 * Create the interpolation matrix.
 *
 * This step is separated from the interpolation step since it only needs
 * to be done once per mesh, while interpolation can be done for several
 * fields.
 *
 * parameters:
 *   r               <-- remapper object
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Copy interpolated values to a new array.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   field_id        <-- id of given field
 *   r               <-- pointer to remapper object
 *   default_val     <-- default value
 *
 * return:
 *   pointer to allocated values array
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Translate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r           <-- pointer to remapper object
 *   translation <-- translation vector
 *----------------------------------------------------------------------------*/

void
cs_medcoupling_remapper_translate(cs_medcoupling_remapper_t  *r,
                                  cs_real_t                   translation[3])
{
  for (int i = 0; i < r->n_fields; i++) {
    r->source_fields[i]->getMesh()->translate(translation);
  }
}

/*----------------------------------------------------------------------------
 * Rotate the mapped source mesh.
 *
 * Caution: cs_medcoupling_remapper_prepare() must to be called after this
 * function in order to update the interpolation matrix.
 *
 * parameters:
 *   r         <-- pointer to remapper object
 *   invariant <-- coordinates of invariant point
 *   axis      <-- rotation axis vector
 *   angle     <-- rotation angle in radians
 *----------------------------------------------------------------------------*/

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
