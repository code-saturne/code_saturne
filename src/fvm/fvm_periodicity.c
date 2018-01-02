/*============================================================================
 * Main structure for handling of periodicity
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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_periodicity.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local structure defining a periodic transformation.
 *----------------------------------------------------------------------------*/

/*
 * Notes:
 *
 * In 3 D, we may have at most a composition of 3 translations, 2 rotations,
 * or 1 translation + 1 rotation (as geometrically equivalent transformations
 * may not be composed).
 * Note that in the translation + rotation case, the order of operations is
 * important. We represent these affine transformations using homogeneous
 * transformations. This allows us to represent all combinations. For example,
 * for a translation A, followed by a rotation B around an axis containing an
 * invariant point which may be obtained by translating the origin by C,
 * the transform of a point x would be defined as: y = (C.B.inv(C))Ax.
 *
 * For vector variables, contrary to coordinates, it is essential to
 * apply only the rotational part of a transformation, hence only the
 * first lines and columns of the matrix (the part corresponding to
 * rotation) will be used.
 */

typedef struct {

  fvm_periodicity_type_t  type; /* Transformation type */

  int     external_num;   /* Given periodicity number (1 to n), with sign
                             indicating direction; 0 for transformations
                             defined by composition */

  int     reverse_id;     /* Id of reverse transformation */
  int     parent_ids[2];  /* Id's of possible parent transformations
                             (which may in turn have parents), or -1 */
  int     equiv_id;       /* Id of first equivalent transformation
                             (own id if no equivalent transformation) */

  double  m[3][4];        /* Matrix of associated transformation
                             (3x4 matrix, 3 first rows of a homogeneous
                             coordinates transformation matrix,
                             with last row = [0 0 0 1]) */

} _transform_t;

/*----------------------------------------------------------------------------
 * Main structure defining all periodicity transformations
 *----------------------------------------------------------------------------*/

struct _fvm_periodicity_t {

  int             n_transforms;    /* Number of transformations */
  _transform_t  **transform;       /* List of transformations */

  int             n_levels;        /* Number of periodicity combination
                                      levels; 0 to 2) */
  int             tr_level_idx[4]; /* Start index of transformations of each
                                      combination level */

  double          equiv_tolerance; /* Relative tolerance for identification
                                      of possibly equivalent directions */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of periodicity type */

const char  *fvm_periodicity_type_name[] = {N_("null"),
                                            N_("translation"),
                                            N_("rotation"),
                                            N_("mixed")};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Combine transformation matrixes.
 *
 * parameters:
 *   a <-- first transformation matrix
 *   b <-- second transformation matrix
 *   c --> combined transformation matrix
 *---------------------------------------------------------------------------*/

static void
_combine_tr_matrixes(const double a[3][4],
                     const double b[3][4],
                     double       c[3][4])
{
  c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
  c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
  c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
  c[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3];

  c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
  c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
  c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
  c[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3];

  c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
  c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
  c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
  c[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3];
}

/*----------------------------------------------------------------------------
 * Check if a combination of transformations is possible.
 *
 * If a transformation direction would appear more than once in a combination,
 * that combination is ignored. If a combination is not commutative (implying
 * at least one rotation component non-orthogonal to a translation or
 * another rotation), the function also returns immediatedly with a specific
 * value.
 *
 * Note that the second transform argument may itself represent a
 * combination of two previous transforms (this may not occur with the
 * first transform due to the way combinations are built in the calling loop).
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   tr_id_0          <-- id of first transformation
 *   tr_id_1          <-- id of second transformation
 *   abort_on_error   <-- abort in case of non-commutative combination ?
 *
 * returns: 0 if combination possible, 1 if combination ignored,
 *          2 if not commutative.
 *---------------------------------------------------------------------------*/

static int
_combine_transforms_check(const fvm_periodicity_t  *this_periodicity,
                          int                       tr_id_0,
                          int                       tr_id_1,
                          _Bool                     abort_on_error)
{
  int i;
  int eq[3];
  int rev[3] = {-1, -1, -1};
  const _transform_t  *tr0 = this_periodicity->transform[tr_id_0];
  const _transform_t  *tr1 = this_periodicity->transform[tr_id_1];

  assert(tr0->parent_ids[0] < 0 && tr0->parent_ids[1] < 0);

  eq[0] = tr_id_0; eq[1] = tr_id_1; eq[2] = -1;

  /* Check that all components are independant */

  if (tr1->parent_ids[1] > -1) {
    assert(tr1->parent_ids[0] > -1);
    eq[1] = tr1->parent_ids[0];
    eq[2] = tr1->parent_ids[1];
  }

  for (i = 0; i < 3; i++) {
    if (eq[i] > -1) {
      eq[i] = (this_periodicity->transform[eq[i]])->equiv_id;
      rev[i] = (this_periodicity->transform[eq[i]])->reverse_id;
    }
  }

  if (   eq[0] == eq[1] || eq[0] == rev[1]
      || rev[0] == eq[1] || rev[0] == rev[1]
      || eq[0] == eq[2] || eq[0] == rev[2]
      || rev[0] == eq[2] || rev[0] == rev[2])
    return 1;

  /* Check for commutativity if necessary */

  if (   tr0->type != FVM_PERIODICITY_TRANSLATION
      || tr1->type != FVM_PERIODICITY_TRANSLATION) {

    int j;
    double m0[3][4], m1[3][4];

    _combine_tr_matrixes(tr0->m, tr1->m, m0);
    _combine_tr_matrixes(tr1->m, tr0->m, m1);

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 4; j++) {

        if (  CS_ABS(m0[i][j] - m1[i][j])
            > this_periodicity->equiv_tolerance) {

          if (abort_on_error) {
            int ext0 = CS_ABS(tr0->external_num);
            int ext1 = CS_ABS(tr1->external_num);
            if (ext1 != 0) {
              bft_error(__FILE__, __LINE__, 0,
                        _("Periodicity transforms %d and %d\n"
                          "(based on directions %d and %d)\n"
                          "are not commutative and may not be combined\n"),
                        tr_id_0, tr_id_1, ext0, ext1);
            }
            else {
              int ext2;
              int j_0 = tr1->parent_ids[0];
              int j_1 = tr1->parent_ids[1];
              ext1 = (this_periodicity->transform[j_0])->external_num;
              ext2 = (this_periodicity->transform[j_1])->external_num;
              ext1 = CS_ABS(ext1);
              ext2 = CS_ABS(ext2);
              bft_error(__FILE__, __LINE__, 0,
                        _("Periodicity transforms %d and %d\n"
                          "(based on directions %d, %d %d)\n"
                          "are not commutative and may not be combined\n"),
                        tr_id_0, tr_id_1, ext0, ext1, ext2);
            }
          } /* end of abort_on_error case */

          return 2;

        }

      }
    }

  } /* End of commutativity check */

  return 0;
}

/*----------------------------------------------------------------------------
 * Sort component transform ids.
 *
 * Transformations built from combination of 2 parents have 2 component ids,
 * and those built from the combination of 3 parents (1 parent plus
 * a combination of 2 parents) have 3 component ids.
 *
 * parameters:
 *   comp_ids  <-> component ids, or -1 (size: 3)
 *
 * returns:
 *   periodicity transform's combination level, or -1 in case of error
 *---------------------------------------------------------------------------*/

static inline void
_sort_component_ids(int  comp_ids[3])
{
  if (comp_ids[1] > -1 && comp_ids[0] > comp_ids[1]) {
    int tmp_id = comp_ids[0];
    comp_ids[0] = comp_ids[1];
    comp_ids[1] = tmp_id;
  }
  if (comp_ids[2] > -1 && comp_ids[2] < comp_ids[1]) {
    int tmp_id = comp_ids[2];
    comp_ids[2] = comp_ids[1];
    if (tmp_id < comp_ids[0]) {
      comp_ids[1] = comp_ids[0];
      comp_ids[0] = tmp_id;
    }
    else
      comp_ids[1] = tmp_id;
  }
}

/*----------------------------------------------------------------------------
 * Determine component equivalent ids for the combination of transformations.
 *
 * Transformations built from combination of 2 parents have 2 component ids,
 * and those built from the combination of 3 parents (1 parent plus
 * a combination of 2 parents) have 3 component ids.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   tr_id_0          <-- id of first transformation
 *   tr_id_1          <-- id of second transformation
 *   comp_ids         --> component ids, or -1 (size: 3)
 *---------------------------------------------------------------------------*/

static void
_component_equiv_ids(fvm_periodicity_t  *this_periodicity,
                     int                 tr_id_0,
                     int                 tr_id_1,
                     int                 comp_id[3])
{
  int i;
  _transform_t  *tr;

  assert(this_periodicity != NULL);

  /* Initialization */

  tr  = this_periodicity->transform[tr_id_1];
  comp_id[0] = tr_id_0;
  if (tr->parent_ids[0] > -1) {
    comp_id[1] = tr->parent_ids[0];
    comp_id[2] = tr->parent_ids[1];
  }
  else {
    comp_id[1] = tr_id_1;
    comp_id[2] = -1;
  }

  /* Replace component ids by their equivalents */

  for (i = 0; i < 3 && comp_id[i] > -1; i++)
    comp_id[i] = (this_periodicity->transform[comp_id[i]])->equiv_id;

  /* Order component ids */

  _sort_component_ids(comp_id);
}

/*----------------------------------------------------------------------------
 * Combine transormations if possible.
 *
 * If a transformation direction would appear more than once in a combination,
 * that combination is ignored. If a combination is not commutative (implying
 * at least one rotation component non-orthogonal to a translation or
 * another rotation), the function also returns immediatedly with a specific
 * value.
 *
 * Note that the reverse id's of combined transformations will be determined
 * later, once all combinations are built.
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *   tr_id_0          <-- id of first transformation
 *   tr_id_1          <-- id of second transformation
 *   tr               --> pointer to the transformation built from
 *                        combination of transforms tr_id_0 and tr_id_1
 *---------------------------------------------------------------------------*/

static void
_combine_transforms(fvm_periodicity_t  *this_periodicity,
                    int                 tr_id_0,
                    int                 tr_id_1,
                    int                 tr_id)
{
  int i, eq_0, eq_1;
  int lv_1, level;

  const _transform_t  *tr0 = this_periodicity->transform[tr_id_0];
  const _transform_t  *tr1 = this_periodicity->transform[tr_id_1];
  _transform_t  *tr  = this_periodicity->transform[tr_id];

  assert(tr_id_0 < this_periodicity->tr_level_idx[1]);

  for (lv_1 = 0;
       lv_1 < 3 && tr_id_1 > this_periodicity->tr_level_idx[lv_1+1];
       lv_1++);
  assert(lv_1 < 3);

  level = lv_1 + 1; /* lv_0 + lv_1 + 1, with lv_0 = 0 */

  assert(this_periodicity != NULL);

  /* Initialize transformation */

  if (tr0->type == tr1->type)
    tr->type = tr0->type;
  else
    tr->type = FVM_PERIODICITY_MIXED;

  tr->external_num = 0;

  tr->reverse_id = -1; /* Will be re-computed later */
  tr->parent_ids[0] = tr_id_0;
  tr->parent_ids[1] = tr_id_1;

  /* Determine previous equivalent transformation id if applicable */

  tr->equiv_id = tr_id;

  eq_0 = tr0->equiv_id;
  eq_1 = tr1->equiv_id;

  if (eq_0 != tr_id_0 || eq_1 != tr_id_1) {

    int comp_id_ref[3];

    _component_equiv_ids(this_periodicity, tr_id_0, tr_id_1, comp_id_ref);

    for (i = this_periodicity->tr_level_idx[level];
         i < tr_id;
         i++) {

      int comp_id_cmp[3];
      _transform_t  *tr_cmp  = this_periodicity->transform[i];

      _component_equiv_ids(this_periodicity,
                           tr_cmp->parent_ids[0],
                           tr_cmp->parent_ids[1],
                           comp_id_cmp);

      if (   comp_id_cmp[0] == comp_id_ref[0]
          && comp_id_cmp[1] == comp_id_ref[1]
          && comp_id_cmp[2] == comp_id_ref[2]) {

        tr->equiv_id = i;
        break;

      }
    }
  }

  /* Compute combined transformation matrix */

  _combine_tr_matrixes(tr0->m, tr1->m, tr->m);

  /* Update periodicity level info */

  if (this_periodicity->n_levels < level + 1)
    this_periodicity->n_levels = level + 1;

  for (i = level + 1; i < 4; i++)
    this_periodicity->tr_level_idx[i] = tr_id + 1;

}

/*----------------------------------------------------------------------------
 * Dump fvm_periodicity_t structure
 *
 *   this_periodicity <-- pointer to the periodicity structure
 *---------------------------------------------------------------------------*/

static void
_transform_dump(_transform_t  *this_transform,
                int            id)
{
  _transform_t  *tr = this_transform;

  bft_printf("\n"
             "  Transform:           %d\n"
             "  Type:                %s\n"
             "  External_num         %d\n"
             "  Reverse id           %d\n"
             "  Parent ids           %d %d\n"
             "  First equivalent id  %d\n",
             id,
             fvm_periodicity_type_name[tr->type],
             tr->external_num,
             tr->reverse_id,
             tr->parent_ids[0], tr->parent_ids[1],
             tr->equiv_id);

  bft_printf("  Matrix:              %12.5g %12.5g %12.5g %12.5g\n"
             "                       %12.5g %12.5g %12.5g %12.5g\n"
             "                       %12.5g %12.5g %12.5g %12.5g\n",
             tr->m[0][0], tr->m[0][1], tr->m[0][2], tr->m[0][3],
             tr->m[1][0], tr->m[1][1], tr->m[1][2], tr->m[1][3],
             tr->m[2][0], tr->m[2][1], tr->m[2][2], tr->m[2][3]);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
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
fvm_periodicity_create(double  equiv_tolerance)
{
  int i;
  fvm_periodicity_t  *period = NULL;

  BFT_MALLOC(period, 1, fvm_periodicity_t);

  period->n_transforms = 0;
  period->transform = NULL;

  period->n_levels = 1;

  for (i = 0; i < 4; i++)
    period->tr_level_idx[i] = 0;

  period->equiv_tolerance = equiv_tolerance;

  return period;
}

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
fvm_periodicity_destroy(fvm_periodicity_t  *this_periodicity)
{
  int  i;

  if (this_periodicity == NULL)
    return NULL;

  for (i = 0; i < this_periodicity->n_transforms; i++)
    BFT_FREE(this_periodicity->transform[i]);

  BFT_FREE(this_periodicity->transform);

  BFT_FREE(this_periodicity);

  return NULL;
}

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
fvm_periodicity_get_n_transforms(const fvm_periodicity_t  *this_periodicity)
{
  int retval = 0;

  if (this_periodicity != NULL)
    retval = this_periodicity->n_transforms;

  return retval;
}

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
fvm_periodicity_get_n_levels(const fvm_periodicity_t  *this_periodicity)
{
  int retval = 0;

  if (this_periodicity != NULL)
    retval = this_periodicity->n_levels;

  return retval;
}

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
                                 int                       tr_level_index[4])
{
  int i;

  for (i = 0; i < 4; i++)
    tr_level_index[i] = 0;

  if (this_periodicity != NULL) {
    for (i = 0; i < 4; i++)
      tr_level_index[i] = this_periodicity->tr_level_idx[i];
  }
}

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
                              double                   matrix[3][4])
{
  int  sub_transf, i, j, k;

  _transform_t  *transform = NULL;

  if (this_periodicity == NULL)
    return -1;

  BFT_REALLOC(this_periodicity->transform,
              this_periodicity->n_transforms + 2,
              _transform_t *);

  for (sub_transf = 0; sub_transf < 2; sub_transf++) {

    /* Create new transformation */

    BFT_MALLOC(transform, 1, _transform_t);

    this_periodicity->transform[this_periodicity->n_transforms] = transform;

    transform->type = type;

    if (sub_transf == 0) {
      transform->external_num = external_num;
      transform->reverse_id = this_periodicity->n_transforms + 1;
    }
    else {
      transform->external_num = - external_num;
      transform->reverse_id = this_periodicity->n_transforms - 1;
    }

    this_periodicity->n_transforms += 1;

    for (i = 1; i < 4; i++)
      this_periodicity->tr_level_idx[i] = this_periodicity->n_transforms;

    for (i = 0; i < 2; i++)
      transform->parent_ids[i] = -1;

    /* Compute transformation matrix */

    if (sub_transf == 0) {
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 4; j++) {
          transform->m[i][j] = matrix[i][j];
        }
      }
    }
    else {
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
          transform->m[i][j] = matrix[j][i]; /* transpose rotation part */
        }
      }
      for (i = 0; i < 3; i++) {
        transform->m[i][3] = 0.;
        for (j = 0; j < 3; j++)
          transform->m[i][3] -= matrix[j][i]*matrix[j][3];
      }
    }

    /* Define equivalent periodicities if necessary */

    transform->equiv_id = this_periodicity->n_transforms - 1;

    for (i = 0; i < this_periodicity->n_transforms - 1; i++) {
      const _transform_t *tr = this_periodicity->transform[i];
      _Bool is_equiv = true;
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 4; k++) {
          if (  CS_ABS(transform->m[j][k] - tr->m[j][k])
              > this_periodicity->equiv_tolerance)
            is_equiv = false;
        }
      }
      if (is_equiv == true) {
        transform->equiv_id = i;
        break;
      }
    }
  }

  /* Return id of first (direct) added transformation */

  return this_periodicity->n_transforms - 2;
}

/*----------------------------------------------------------------------------
 * Add a translation type periodicity direction.
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
                                const double        translation[3])
{
  double  matrix[3][4] = {{1., 0., 0., 0.},
                          {0., 1., 0., 0.},
                          {0., 0., 1., 0.}};

  matrix[0][3] = translation[0];
  matrix[1][3] = translation[1];
  matrix[2][3] = translation[2];

  return fvm_periodicity_add_by_matrix(this_periodicity,
                                       external_num,
                                       FVM_PERIODICITY_TRANSLATION,
                                       matrix);
}

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
                             const double        invariant_point[3])
{
  int  i, j;
  double norm;
  double direction[3];
  double rot[3][3];
  double matrix[3][4];

  const double pi = 4 * atan(1);
  const double theta = pi * angle/180.;
  const double cost = cos(theta);
  const double sint = sin(theta);
  const double onemcost = (1.0 - cost);

  /* Compute the rotation matrix, using formula:
   *  R = (1-cos(theta))axis.transp(axis) + cos(theta)I + sin(theta)V
   *
   *           [ 0            -direction(3)  direction(2)]
   *  with V = [ direction(3)       0       -direction(1)]
   *           [-direction(2)  direction(1)       0      ]
   */

  norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

  direction[0] = axis[0] / norm;
  direction[1] = axis[1] / norm;
  direction[2] = axis[2] / norm;

  /* first row of rotation maxtrix */
  rot[0][0] = onemcost*direction[0]*direction[0] + cost;
  rot[0][1] = onemcost*direction[0]*direction[1] - sint*direction[2];
  rot[0][2] = onemcost*direction[0]*direction[2] + sint*direction[1];

  /* second row of rotation maxtrix */
  rot[1][0] = onemcost*direction[1]*direction[0] + sint*direction[2];
  rot[1][1] = onemcost*direction[1]*direction[1] + cost;
  rot[1][2] = onemcost*direction[1]*direction[2] - sint*direction[0];

  /* third row of rotation maxtrix */
  rot[2][0] = onemcost*direction[2]*direction[0] - sint*direction[1];
  rot[2][1] = onemcost*direction[2]*direction[1] + sint*direction[0];
  rot[2][2] = onemcost*direction[2]*direction[2] + cost;

  /* Now compute full rotation matrix in homogeneous coordinates,
   * accounting for invariant point of coordiantes t[], with the formula:
   *
   *     [1 0 0 t[0]] [r[0][0] r[0][1] r[0][3] 0] [1 0 0 -t[0]]
   * M = [0 1 0 t[1]].[r[1][0] r[1][1] r[1][3] 0].[0 1 0 -t[1]]
   *     [0 0 1 t[2]] [r[2][0] r[2][1] r[2][3] 0] [0 0 1 -t[2]]
   *     [0 0 0 1   ] [0       0       0       1] [0 0 0  1]
   */

  for (i = 0; i < 3; i++) {       /* rotation part of matrix */
    for (j = 0; j < 3; j++) {
      matrix[i][j] = rot[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    matrix[i][3] = invariant_point[i];
    for (j = 0; j < 3; j++)
      matrix[i][3] -= rot[i][j]*invariant_point[j];
  }

  /* Clip "zero" values */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 4; j++) {
      if (CS_ABS(matrix[i][j]) < 1.e-16)
        matrix[i][j] = 0.;
    }
  }

  /* Now add periodicity */

  return fvm_periodicity_add_by_matrix(this_periodicity,
                                       external_num,
                                       FVM_PERIODICITY_ROTATION,
                                       matrix);
}

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
                                 int                       direction)
{
  int  i;

  assert(direction == 1 || direction == -1);

  if (this_periodicity != NULL) {

    for (i = 0; i < this_periodicity->n_transforms; i++) {

      _transform_t  *transform = this_periodicity->transform[i];

      if (transform->external_num == external_num * direction) {
        if (   (direction > 0 && transform->reverse_id > i)
            || (direction < 0 && transform->reverse_id < i))
          return i;
      }

    }

  }

  return -1;
}

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
                         int                       tr_id)
{
  if (this_periodicity == NULL)
    return FVM_PERIODICITY_NULL;

  if (tr_id < 0 || tr_id >= this_periodicity->n_transforms)
    return FVM_PERIODICITY_NULL;

  return (this_periodicity->transform[tr_id])->type;
}

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
                               int                       tr_id)
{
  int  retval = -1;

  if (this_periodicity != NULL) {

    if (tr_id > -1 && tr_id < this_periodicity->n_transforms)
      retval = this_periodicity->transform[tr_id]->reverse_id;

  }

  return retval;
}

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
                               int                       parent_ids[2])
{
  _transform_t  *transform;

  /* Initialization */

  if (parent_ids == NULL)
    return;

  parent_ids[0] = -1;
  parent_ids[1] = -1;

  if (this_periodicity == NULL)
    return;

  if (tr_id < 0 || tr_id >= this_periodicity->n_transforms)
    return;

  /* Return parent ids */

  transform = this_periodicity->transform[tr_id];

  parent_ids[0] = transform->parent_ids[0];
  parent_ids[1] = transform->parent_ids[1];
}

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
                               int                       component_ids[3])
{
  int tr_id_1;
  _transform_t  *transform;

  /* Initialization */

  if (component_ids == NULL)
    return;

  if (this_periodicity == NULL)
    return;

  if (tr_id < 0 || tr_id >= this_periodicity->n_transforms)
    return;

  /* Component ids */

  transform = this_periodicity->transform[tr_id];

  tr_id_1 = transform->parent_ids[1];

  if (tr_id_1 > -1) { /* 2 or more components */

    component_ids[0] = transform->parent_ids[0];
    if (tr_id_1 < this_periodicity->tr_level_idx[1]) {
      component_ids[1] = transform->parent_ids[1];
      component_ids[2] = -1;
    }
    else {
      _transform_t  *tr_1;
      assert(tr_id_1 < this_periodicity->tr_level_idx[2]);
      tr_1 = this_periodicity->transform[tr_id_1];
      component_ids[1] = tr_1->parent_ids[0];
      component_ids[2] = tr_1->parent_ids[1];
    }

    _sort_component_ids(component_ids);
  }

  else { /* Only 1 component */

    component_ids[0] = tr_id;
    component_ids[1] = -1;
    component_ids[2] = -1;

  }
}

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
                             int                       tr_id)
{
  int  retval = -1;

  if (this_periodicity != NULL) {

    if (tr_id > -1 && tr_id < this_periodicity->n_transforms)
      retval = this_periodicity->transform[tr_id]->equiv_id;

  }

  return retval;
}

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
                           double                    matrix[3][4])
{
  int i, j;
  _transform_t  *transform = NULL;

  if (this_periodicity != NULL)
    if (tr_id >= 0 && tr_id < this_periodicity->n_transforms) {

      transform = this_periodicity->transform[tr_id];

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 4; j++)
          matrix[i][j] = transform->m[i][j];
      }
    }

  if (transform == NULL) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 4; j++)
        matrix[i][j] = 0.0;
    }
  }
}

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
                        int                 abort_on_error)
{
  int i, j, tr_count;
  int j_start, j_end;
  int level;
  int n_level_0, n_level_1;
  int n_transforms_max;
  int retval;

  if (this_periodicity == NULL)
    return;

  /* Check that only level 0 periodicities are defined */

  if (  this_periodicity->tr_level_idx[2]
      > this_periodicity->tr_level_idx[1])
    return;

  n_level_0 = this_periodicity->n_transforms;
  n_level_1 = 0; /* temporary value */

  /* Build combinations */

  for (level = 1; level < 3; level++) {

    tr_count = this_periodicity->n_transforms;

    /* Count possible combinations for this level */

    if (level == 1)
      n_transforms_max = this_periodicity->n_transforms + (n_level_0*n_level_0);
    else /* if (level == 2) */
      n_transforms_max = this_periodicity->n_transforms + (n_level_0*n_level_1);

    BFT_REALLOC(this_periodicity->transform,
                n_transforms_max,
                _transform_t *);

    for (i = 0; i < n_level_0; i++) {

      _transform_t  *tr0 = this_periodicity->transform[i];

      if (level == 1) {
        j_start = i+1;
        j_end = n_level_0;
      }
      else { /* if level == 2 */
        j_start = n_level_0;
        j_end = n_level_0 + n_level_1;
      }

      for (j = j_start; j < j_end; j++) {

        _transform_t  *tr1 = this_periodicity->transform[j];

        /* Avoid combining transforms from the same periodicity direction */
        if (tr0->reverse_id == j || tr1->reverse_id == i)
          continue;

        /* Ensure that for any combination i, {j,k}, i < j < k
           to avoid multiple equivalent definitions */

        if (tr1->parent_ids[0] > -1 && tr1->parent_ids[0] < i)
          continue;

        /* Combine independent transformations */

        retval = _combine_transforms_check(this_periodicity,
                                           i,
                                           j,
                                           (_Bool)abort_on_error);

        if (retval == 0) {

          BFT_MALLOC(this_periodicity->transform[tr_count], 1, _transform_t);

          _combine_transforms(this_periodicity,
                              i,
                              j,
                              tr_count);

          /* Increment number of transforms */

          tr_count += 1;

        }

      }

    }

    /* Determine which new transformations are the reverses of others */

    for (i = this_periodicity->n_transforms; i < tr_count; i++) {

      _transform_t *tr = this_periodicity->transform[i];

      int rev_0 = (this_periodicity->transform[tr->parent_ids[0]])->reverse_id;
      int rev_1 = (this_periodicity->transform[tr->parent_ids[1]])->reverse_id;

      for (j = i; j < tr_count; j++) {

        _transform_t *tr_cmp = this_periodicity->transform[j];

        if (   (   tr_cmp->parent_ids[0] == rev_0
                && tr_cmp->parent_ids[1] == rev_1)
            || (   tr_cmp->parent_ids[0] == rev_1
                && tr_cmp->parent_ids[1] == rev_0)) {
          tr->reverse_id = j;
          tr_cmp->reverse_id = i;
        }
      }

    }

    /* Update number of periodic transforms */

    this_periodicity->n_transforms = tr_count;

    if (level == 1)
      n_level_1 = tr_count - n_level_0;
  }

  BFT_REALLOC(this_periodicity->transform,
              this_periodicity->n_transforms,
              _transform_t *);
}

/*----------------------------------------------------------------------------
 * Dump fvm_periodicity_t structure
 *
 * parameters:
 *   this_periodicity <-- pointer to the periodicity structure
 *---------------------------------------------------------------------------*/

void
fvm_periodicity_dump(const fvm_periodicity_t  *this_periodicity)
{
  int  i;
  int  level = 0;

  bft_printf("\n"
             "Periodicity:          %p\n", (const void *)this_periodicity);

  if (this_periodicity == NULL) {
    bft_printf("\n");
    return;
  }

  bft_printf("Number of transforms  %d\n"
             "Number of levels  %d\n"
             "Levels index      %d %d %d %d\n"
             "Equivalence tolerance %12.5g\n",
             this_periodicity->n_transforms,
             this_periodicity->n_levels,
             this_periodicity->tr_level_idx[0],
             this_periodicity->tr_level_idx[1],
             this_periodicity->tr_level_idx[2],
             this_periodicity->tr_level_idx[3],
             this_periodicity->equiv_tolerance);

  for (i = 0; i < this_periodicity->n_transforms; i++) {

    if (i == this_periodicity->tr_level_idx[level]) {
      bft_printf("\n  Combination level %d\n", level);
      level += 1;
    }

    _transform_dump(this_periodicity->transform[i], i);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
