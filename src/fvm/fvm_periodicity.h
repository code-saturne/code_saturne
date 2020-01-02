#ifndef __FVM_PERIODICITY_H__
#define __FVM_PERIODICITY_H__

/*============================================================================
 * Main structure for handling of periodicities
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Periodicity types
 *----------------------------------------------------------------------------*/

typedef enum {

  FVM_PERIODICITY_NULL,
  FVM_PERIODICITY_TRANSLATION,
  FVM_PERIODICITY_ROTATION,
  FVM_PERIODICITY_MIXED

} fvm_periodicity_type_t;

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef struct _fvm_periodicity_t fvm_periodicity_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of periodicity type */

extern const char  *fvm_periodicity_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty periodicity definition structure.
 *
 * To create a composed periodicity, use fvm_periodicity_compose().
 *
 * parameters:
 *   equiv_tolerance <-- relative tolerance for identification of
 *                       possibly equivalent directions
 *
 * returns:
 *   pointer to the created fvm_periodicity_t structure.
 *---------------------------------------------------------------------------*/

fvm_periodicity_t *
fvm_periodicity_create(double  equiv_tolerance);

/*----------------------------------------------------------------------------
 * Destroy an fvm_periodicity_t structure.
 *
 * parameters:
 *   this_periodicity  <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *---------------------------------------------------------------------------*/

fvm_periodicity_t *
fvm_periodicity_destroy(fvm_periodicity_t  *this_periodicity);

/*----------------------------------------------------------------------------
 * Return number of transformations associated with peridocity.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *
 * returns:
 *   number of periodicity transformations
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_get_n_transforms(const fvm_periodicity_t  *this_periodicity);

/*----------------------------------------------------------------------------
 * Return number of periodicity combination levels.
 *
 * Combination level is 0 for single periodicity, and 1 or 2 for
 * combinations of periodicities (2 max in 3d).
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *
 * returns:
 *   number of periodicity combination levels (1 to 3)
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_get_n_levels(const fvm_periodicity_t  *this_periodicity);

/*----------------------------------------------------------------------------
 * Return index of periodicity transform combination levels.
 *
 * Combination level is 0 for single periodicity, and 1 or 2 for
 * combinations of periodicities (2 max in 3d).
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   tr_level_index   --> start index of first transform of each level
 *                        (size: 3 + 1)
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_get_tr_level_idx(const fvm_periodicity_t  *this_periodicity,
                                 int                       tr_level_index[4]);

/*----------------------------------------------------------------------------
 * Add a general periodicity direction.
 *
 * For each direction defined, two transformations are defined: direct
 * and reverse. The id of the reverse transformation is equal to the
 * id of the direct transformation + 1.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   external_num     <-- external number (1 to n) associated with direction
 *   type             <-- transformation type (translation or rotation)
 *   matrix           <-- transformation matrix (3x4 matrix, 3 first rows
 *                        of a homogeneous coordinates transformation
 *                        matrix, with last row = [0 0 0 1])
 *
 * returns:
 *   id of the associated direct transformation.
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_add_by_matrix(fvm_periodicity_t       *this_periodicity,
                              int                      external_num,
                              fvm_periodicity_type_t   type,
                              double                   matrix[3][4]);

/*----------------------------------------------------------------------------
 * Add a translation-type periodicity direction.
 *
 * For each direction defined, two transformations are defined: direct
 * and reverse. The id of the reverse transformation is equal to the
 * id of the direct transformation + 1.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   external_num     <-- external number (1 to n) associated with direction
 *   translation      <-- components of translation vector (3)
 *
 * returns:
 *   id of the associated direct transformation.
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_add_translation(fvm_periodicity_t  *this_periodicity,
                                int                 external_num,
                                const double        translation[3]);

/*----------------------------------------------------------------------------
 * Add a rotation type periodicity direction.
 *
 * For each direction defined, two transformations are defined: direct
 * and reverse. The id of the reverse transformation is equal to the
 * id of the direct transformation + 1.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   external_num     <-- external number (1 to n) associated with direction
 *   angle            <-- rotation angle in degrees
 *   axis             <-- components of rotation axis direction vector (3)
 *   invariant_point  <-- components of invariant point (3)
 *
 * returns:
 *   id of the associated direct transformation.
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_add_rotation(fvm_periodicity_t  *this_periodicity,
                             int                 external_num,
                             double              angle,
                             const double        axis[3],
                             const double        invariant_point[3]);

/*----------------------------------------------------------------------------
 * Return the periodicity transform id associated with an external number.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   external_num     <-- external number (1 to n) associated with direction
 *   direction        <-- 1 for direct transformation, -1 for reverse
 *
 * returns:
 *  transformation id (0 to n-1), or -1 if not found.
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_get_transform_id(const fvm_periodicity_t  *this_periodicity,
                                 int                       external_num,
                                 int                       direction);

/*----------------------------------------------------------------------------
 * Return a periodicity transformation's type.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *
 * returns:
 *   periodicity transform's type
 *---------------------------------------------------------------------------*/

fvm_periodicity_type_t
fvm_periodicity_get_type(const fvm_periodicity_t  *this_periodicity,
                         int                       tr_id);

/*----------------------------------------------------------------------------
 * Return the periodicity transform reverse's id.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *
 * returns:
 *  reverse transformation id (0 to n-1), or -1 if not found.
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_get_reverse_id(const fvm_periodicity_t  *this_periodicity,
                               int                       tr_id);

/*----------------------------------------------------------------------------
 * Return a periodicity transformation's parents.
 *
 * A standard transformation has combination level 0, and no parents.
 * Transformations built from combination of 2 parents have combination
 * level 1, and those built from the combination of 3 parents have
 * combination level 2.
 *
 * Level 2 transformations are built from a combination of a level 0
 * transformation with a level 1 transformation (itself a combination of
 * 2 level 0 transformations).
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *   parent_ids       --> parent ids, or -1 (size: 2)
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_get_parent_ids(const fvm_periodicity_t  *this_periodicity,
                               int                       tr_id,
                               int                       parent_ids[2]);

/*----------------------------------------------------------------------------
 * Return a periodicity transformation's component ids.
 *
 * A standard transformation has combination level 0, and no parents.
 * Transformations built from combination of 2 parents have combination
 * level 1, and those built from the combination of 3 parents have
 * combination level 2.
 *
 * Level 2 transformations are built from a combination of a level 0
 * transformation with a level 1 transformation (itself a combination of
 * 2 level 0 transformations). Component ids allow direct access to the 3
 * corresponding level 0 combinations.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *   component_ids    --> component ids, or -1 (size: 3)
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_get_components(const fvm_periodicity_t  *this_periodicity,
                               int                       tr_id,
                               int                       component_ids[3]);

/*----------------------------------------------------------------------------
 * Return a periodicity transformation's first equivalent transform id.
 *
 * If multiple transformations are equivalent, each will point to
 * the previous equivalent transformation, defining a reversed linked list.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *
 * returns:
 *  id of first equivalent transformation
 *---------------------------------------------------------------------------*/

int
fvm_periodicity_get_equiv_id(const fvm_periodicity_t  *this_periodicity,
                             int                       tr_id);

/*----------------------------------------------------------------------------
 * Return a periodicity transformation's matrix.
 *
 * parameters:
 *   this_periodicity <-- pointer to periodicity structure
 *   tr_id            <-- id of transformation we are interested in
 *   matrix           --> coefficients of transformation matrix
 *
 * returns:
 *   periodicity transform's matrix
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_get_matrix(const fvm_periodicity_t  *this_periodicity,
                           int                       tr_id,
                           double                    matrix[3][4]);

/*----------------------------------------------------------------------------
 * Complete periodicity information with combined transformations.
 *
 * This function should only be called once, after all base periodicity
 * transforms have been defined. It returns immediately if combined
 * transforms are already defined.
 *
 * parameters:
 *   this_periodicity <-> pointer to the periodicity structure
 *   abort_on_error   <-- 0: non-commuting combinations are discarded
 *                        1: abort in presence of non-commuting combinations
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_combine(fvm_periodicity_t  *this_periodicity,
                        int                 abort_on_error);

/*----------------------------------------------------------------------------
 * Dump fvm_periodicity_t structure
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_dump(const fvm_periodicity_t  *this_periodicity);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_PERIODICITY_H__ */
