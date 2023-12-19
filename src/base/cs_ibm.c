/*============================================================================
 * Define immersed boundaries based on user inputs (experimental).
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_cell_to_vertex.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_medcoupling_intersector.h"
#include "cs_meg_prototypes.h"
#include "cs_mesh_adjacencies.h"
#include "cs_parameters.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_stl.h"
#include "cs_turbomachinery.h"
#include "cs_velocity_pressure.h"
#include "cs_vertex_to_cell.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ibm.h"

/*----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

BEGIN_C_DECLS

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_ibm.c
        Define immersed boundaries based on user inputs (experimental).
        Cloud of points are dealt with porosity from scan files.
*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_porosity_ibm_opt_t _porosity_ibm_opt = {
  .porosity_mode = 0
};

cs_porosity_ibm_opt_t *cs_glob_porosity_ibm_opt
= &(_porosity_ibm_opt);

/* Pointer to cs_ibm_t structure for the main initialization */
cs_ibm_t  *cs_ibm = NULL;

/* Names of algorithms */
const char *_ibm_algo_names[] = {"CS_IBM_ALGO_NONE",
                                 "CS_IBM_ALGO_CUT_CELLS",
                                 "CS_IBM_ALGO_MEDCOUPLING",
                                 "CS_IBM_ALGO_STL"};

const char *_ibm_obj_property_names[] = {"density",
                                         "mass",
                                         "inertia matrix",
                                         "cp",
                                         "lambda",
                                         "stiffness",
                                         "damping",
                                         "Young module",
                                         "Inertia momentum",
                                         "cross section",
                                         "rayleigh_coeff_a",
                                         "rayleigh_coeff_b"};

const char *_ibm_obj_init_vals_names[] = {"Equilibrium Center of gravity",
                                          "Center of gravity",
                                          "Angle",
                                          "Velocity",
                                          "Acceleration",
                                          "Angular velocity",
                                          "Fluid force"};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_immersed_boundaries(void);

void
cs_f_porosity_ibm_get_pointer(int **ibm_porosity_mode);

/*----------------------------------------------------------------------------
 * Wrapper immersed boundary function, intended for use by Fortran wrapper only.
 *----------------------------------------------------------------------------*/

void
cs_f_immersed_boundaries(void)
{
  cs_immersed_boundaries(cs_glob_mesh, cs_glob_mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Get pointer
 *----------------------------------------------------------------------------*/

void
cs_f_porosity_ibm_get_pointer(int **ibm_porosity_mode)
{
  *ibm_porosity_mode
    = &(_porosity_ibm_opt.porosity_mode);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the face porosity depending on the neighbouring cells
 *          porosities based on geometry.
 *
 * \param[in]  alphai       value at neighbouring cell i
 * \param[in]  alphaj       value at neighbouring cell j
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_geom_face_fraction(cs_real_t alphai,
                    cs_real_t alphaj)
{
  cs_real_t alpi = cs_math_fmin(alphai, alphaj);
  cs_real_t alpj = cs_math_fmax(alphai, alphaj);

  cs_real_t alpij = 0.5 * (alpi + alpj);

  if (alpj < 1.e-10)
    return alpij;

  if (alpi > 1. - 1.e-10)
    return alpij;

  /* Four cases are possible : faces are cut each one on Y, one on X and one
   * on Y, one on X and so the other one, or Y and X. The minimax of
   * initialization of alpi alpj guarantees the positivity of the slope. */

  /* Case YY */
  cs_real_t x1 = 0.;
  cs_real_t x2 = 0.;

  /* Case YX */
  cs_real_t bb = 4. * alpj - 6.;
  cs_real_t cc = 1. + 4. * alpi * (1. - alpj);
  cs_real_t delta = bb * bb - 4. * cc;

  if (delta >= 0.)
    x2 = 0.5 * (-bb - sqrt(delta));
  else
    x2 = -1.;

  if (x2 >= alpi && x2 <= 2. * alpi && x2 <= 2. * alpj - 1.)
    return x2;

  x1 = 0.5 * (alpi + alpj);
  if (x1 <= 2. * alpi && x1 >= 2. * alpj - 1.)
    return x1;

  /* Case XX */
  cs_real_t aa = 1. - (alpi + alpj);
  bb = 2. * alpi;
  cc = -alpi;
  delta = bb * bb - 4. * aa * cc;

  if (delta >= 0. && cs_math_fabs(aa) > 0.) {
    cs_real_t sqrt_delta = sqrt(delta);
    cs_real_t den = 0.5 / aa;
    x1 = (-bb + sqrt_delta) * den;
    x2 = (-bb - sqrt_delta) * den;
  } else {
    x1 = -1.;
    x2 = -1.;
  }

  if (x1 >= 2. * alpi && x1 <= 2. * alpj - 1.)
    return x1;
  else if (x2 >= 2. * alpi && x2 <= 2 * alpj - 1.)
    return x2;

  /* Case XY */
  bb = 4. * alpi;
  cc = -4. * alpi * alpj;
  delta = bb * bb - 4. * cc;

  if (delta >= 0.)
    x1 = 0.5 * (-bb + sqrt(delta));
  else
    x1 = -1.;

  if (x1 >= 2. * alpi && x1 >= 2. * alpj - 1. && x1 <= alpj)
    return x1;

  cs_real_t eps = 1.e-10;
  if (alpi < 1.-eps && alpj  > 1.-eps)
    return alpj;
  if (alpi < eps && alpj  > eps)
    return alpi;

  return alpij;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by dichotomy the length of the immersed part of a line
 *          between two points (i.e. the distance between the point in the
 *          solid and the immersed boundary) based on the cut-cell method.
 *
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_imm_lgth_cutcell(cs_real_3_t x1,
                  cs_real_3_t x2,
                  cs_real_t   t,
                  int         num_object)
{
  cs_real_t length = 0.;

  int ipenal1 = cs_ibm_object_compute_cut_porosity(1, x1, t, num_object);
  int ipenal2 = cs_ibm_object_compute_cut_porosity(2, x2, t, num_object);

  if (ipenal1 + ipenal2 == 2)
    return length;
  else if (ipenal1 + ipenal2 == 0)
    return cs_math_3_distance(x1, x2);
  else {
    cs_real_3_t xx1, xx2, x3, xout;
    int sens = 0;

    for (int idim = 0; idim < 3; idim++) {
      xx1[idim] = x1[idim];
      xx2[idim] = x2[idim];
      xout[idim] = xx1[idim];
    }

    if (ipenal1 == 1) {
      sens = 1;
      for (int idim = 0; idim < 3; idim++)
        xout[idim] = xx2[idim];
    }

    for (int isou = 0; isou < 10; isou++) {
      for (int idim = 0; idim < 3; idim++)
        x3[idim] = 0.5 * (xx1[idim] + xx2[idim]);

      int ipenal3 = cs_ibm_object_compute_cut_porosity(3, x3, t, num_object);

      if (ipenal3 == sens)
        for (int idim = 0; idim < 3; idim++)
          xx1[idim] = x3[idim];
      else
        for (int idim = 0; idim < 3; idim++)
          xx2[idim] = x3[idim];
    }

    return cs_math_3_distance(xx1, xout);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by dichotomy the length of the immersed part of a line
 *          between two points (i.e. the distance between the point in the
 *          solid and the immersed boundary) based on the input porosities.
 *
 * \param[in]  x1           point 1
 * \param[in]  por1         porosity at 1
 * \param[in]  x2           point 2
 * \param[in]  por2         porosity at 2
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_imm_lgth_poro(cs_real_3_t x1,
               cs_real_t   por1,
               cs_real_3_t x2,
               cs_real_t   por2)
{
  /* Both considered as solids, returns 0 */
  if (por1 < 0.5 && por2 < 0.5)
    return 0.;

  /* Both considered as fluids, returns total distance */
  else if (por1 > 0.5 && por2 > 0.5)
    return cs_math_3_distance(x1, x2);

  /* Fluid-solid cases, return weighted distance  */
  else if (por1 < 0.5 && por2 >= 0.5)
    return (por2 - 0.5) / (por2 - por1) * cs_math_3_distance(x1, x2);

  else if (por2 < 0.5 && por1 >= 0.5)
    return (por1 - 0.5) / (por1 - por2) * cs_math_3_distance(x1, x2);

  /* Medium case */
  else
    return 0.5 * cs_math_3_distance(x1, x2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume of a tetrahedron described by its vertices
 *          (x1,x2,x3,x4).
 *
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_tetra_vol(cs_real_3_t x1,
           cs_real_3_t x2,
           cs_real_3_t x3,
           cs_real_3_t x4)
{
  /* Volume of the tetrahedron x1, x2, x3, x4 */
  cs_real_t tetra_vol = cs_math_fabs(( (x1[0]-x4[0])*(x2[1]-x4[1])
                                     - (x1[1]-x4[1])*(x2[0]-x4[0]) )
                                     * (x3[2]-x4[2])
                                 + (   (x1[1]-x4[1])*(x2[2]-x4[2])
                                     - (x1[2]-x4[2])*(x2[1]-x4[1]) )
                                     * (x3[0]-x4[0])
                                 + (   (x1[2]-x4[2])*(x2[0]-x4[0])
                                     - (x1[0]-x4[0])*(x2[2]-x4[2]) )
                                     * (x3[1]-x4[1]) )
                                 * cs_math_1ov6;

  return tetra_vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume of a pyramid with a quadrangle base
 *          (x1,x2,x3,x4) and apex (x5).
 *
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 * \param[in]  x5           point 5 (apex)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_pyram_vol(cs_real_3_t x1,
           cs_real_3_t x2,
           cs_real_3_t x3,
           cs_real_3_t x4,
           cs_real_3_t x5)
{
  /* Volume of the pyramid of base x1, x2, x3, x4 and apex x5 */
  cs_real_3_t xc;

  for (int i = 0; i < 3; i++)
    xc[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);

  cs_real_t vol12 = _tetra_vol(x1, x2, xc, x5);
  cs_real_t vol23 = _tetra_vol(x2, x3, xc, x5);
  cs_real_t vol34 = _tetra_vol(x3, x4, xc, x5);
  cs_real_t vol41 = _tetra_vol(x4, x1, xc, x5);

  cs_real_t pyram_vol = vol12 + vol23 + vol34 + vol41;

  return pyram_vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume of a prism with a quadrangle base
 *          (x3,x4,x5,x6) and edge (x1 x2).
 *
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 * \param[in]  x5           point 5
 * \param[in]  x6           point 6
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_prism_vol(cs_real_3_t x1,
           cs_real_3_t x2,
           cs_real_3_t x3,
           cs_real_3_t x4,
           cs_real_3_t x5,
           cs_real_3_t x6)
{
  /* Volume of the prism of base x3, x4, x5, x6 and edge x1 x2 */
  /* The two triangles are x1 x3 x6 and x2 x4 x5 */

  cs_real_3_t xc;

  for (int i = 0; i < 3; i++)
    xc[i] = (x1[i]+x2[i]+x3[i]+x4[i]+x5[i]+x6[i]) * cs_math_1ov6;

  cs_real_t vol136c  = _tetra_vol(x1, x3, x6, xc);
  cs_real_t vol245c  = _tetra_vol(x2, x4, x5, xc);
  cs_real_t vol1256c = _pyram_vol(x1, x2, x5, x6, xc);
  cs_real_t vol1243c = _pyram_vol(x1, x2, x4, x3, xc);
  cs_real_t vol3456c = _pyram_vol(x3, x4, x5, x6, xc);

  cs_real_t prism_vol = vol136c + vol245c + vol1256c +  vol1243c + vol3456c;

  return prism_vol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the volume and center of gravity of a tetrahedron described
 *          by its vertices (x1,x2,x3,x4) but truncated due to solid parts
 *          identified with porosities, through a recursive approach.
 *
 * \param[out] vol          volume
 * \param[out] cog          center of gravity
 * \param[in]  x1           point 1
 * \param[in]  por1         porosity at point 1
 * \param[in]  x2           point 2
 * \param[in]  por2         porosity at point 2
 * \param[in]  x3           point 3
 * \param[in]  por3         porosity at point 3
 * \param[in]  x4           point 4
 * \param[in]  por4         porosity at point 4
 * \param[in]  icut         number of sub-cut for cells in cut-cells algorithm
 */
/*----------------------------------------------------------------------------*/

static void
_tetra_vol_poro(cs_real_t   *vol,
                 cs_real_3_t  cog,
                 cs_real_3_t  x1,
                 cs_real_t    por1,
                 cs_real_3_t  x2,
                 cs_real_t    por2,
                 cs_real_3_t  x3,
                 cs_real_t    por3,
                 cs_real_3_t  x4,
                 cs_real_t    por4,
                 int          icut)
{
  /* Mean porosity */
  cs_real_t porc = (por1 + por2 + por3 + por4) * 0.25;

  /* Check if some vertices are considered as solid */
  int ipenal1 = 0;
  if (por1 < 0.5)
    ipenal1 = 1;

  int ipenal2 = 0;
  if (por2 < 0.5)
    ipenal2 = 1;

  int ipenal3 = 0;
  if (por3 < 0.5)
    ipenal3 = 1;

  int ipenal4 = 0;
  if (por4 < 0.5)
    ipenal4 = 1;

  int ipenalc = 0;
  if (porc < 0.5)
    ipenalc = 1;

  int cpt4 = ipenal1 + ipenal2 + ipenal3 + ipenal4;
  int cpt5 = cpt4 + ipenalc;

  *vol = 0.;

  /* Check if the volume is considered as completely solid, null volume */
  if (cpt5 == 5)
    return;
  /* Check if no solid vertex, complete volume */
  else if (cpt5 == 0 && icut > 0) {

    *vol = _tetra_vol(x1, x2, x3, x4);

    for (int i = 0; i < 3; i++)
      cog[i] += 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]) * *vol;

    return;

  /* Check if at least one solid vertex, subdivision of the edges and
   * distribution of porosities before recursive call */
  } else if (icut > 0) {
    icut--;

    cs_real_3_t x12, x13, x14, x23, x24, x34;
    cs_real_t por12, por13, por14, por23, por24, por34;
    for (int i = 0; i < 3; i++) {
      x12[i] = 0.5 * (x1[i] + x2[i]);
      x13[i] = 0.5 * (x1[i] + x3[i]);
      x14[i] = 0.5 * (x1[i] + x4[i]);
      x23[i] = 0.5 * (x2[i] + x3[i]);
      x24[i] = 0.5 * (x2[i] + x4[i]);
      x34[i] = 0.5 * (x3[i] + x4[i]);
    }
    por12 = 0.5 * (por1 + por2);
    por13 = 0.5 * (por1 + por3);
    por14 = 0.5 * (por1 + por4);
    por23 = 0.5 * (por2 + por3);
    por24 = 0.5 * (por2 + por4);
    por34 = 0.5 * (por3 + por4);

    cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
    _tetra_vol_poro(&vol1, cog, x1,  por1,  x12, por12, x13, por13,
                    x14, por14, icut);
    _tetra_vol_poro(&vol2, cog, x2,  por2,  x23, por23, x12, por12,
                    x24, por24, icut);
    _tetra_vol_poro(&vol3, cog, x3,  por3,  x13, por13, x23, por23,
                    x34, por34, icut);
    _tetra_vol_poro(&vol4, cog, x4,  por4,  x14, por14, x24, por24,
                    x34, por34, icut);
    _tetra_vol_poro(&vol5, cog, x12, por12, x13, por13, x34, por34,
                    x23, por23, icut);
    _tetra_vol_poro(&vol6, cog, x12, por12, x13, por13, x34, por34,
                    x14, por14, icut);
    _tetra_vol_poro(&vol7, cog, x34, por34, x24, por24, x12, por12,
                    x23, por23, icut);
    _tetra_vol_poro(&vol8, cog, x34, por34, x24, por24, x12, por12,
                    x14, por14, icut);

    *vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
    return;

  } else {
    if (cpt4 == 0) {
      /* No penalized point -> full tetrahedron */

      *vol = _tetra_vol(x1, x2, x3, x4);

      for (int i = 0; i < 3; i++)
        cog[i] += 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]) * *vol;

    } else if (cpt4 == 1) {
      /* One penalized point tetrahedron - tetra from the penalized point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      *vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 1) {
        cs_real_t l12 = cs_math_3_norm(d12);
        cs_real_t l13 = cs_math_3_norm(d13);
        cs_real_t l14 = cs_math_3_norm(d14);
        cs_real_t l12a = l12 - _imm_lgth_poro(x1, por1, x2, por2);
        cs_real_t l13a = l13 - _imm_lgth_poro(x1, por1, x3, por3);
        cs_real_t l14a = l14 - _imm_lgth_poro(x1, por1, x4, por4);

        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd13 = l13a / l13;
        cs_real_t lbd14 = l14a / l14;

        cs_real_t volp = *vol * lbd12 * lbd13 * lbd14;

        cs_real_3_t x12, x13, x14;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x12[i] = x1[i] + lbd12 * d12[i];
          x13[i] = x1[i] + lbd13 * d13[i];
          x14[i] = x1[i] + lbd14 * d14[i];
          cogp[i] = 0.25 * (x1[i] + x12[i] + x13[i] + x14[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * *vol - cogp[i] * volp;
        }

        *vol *= (1. - lbd12 * lbd13 * lbd14);

        return;

      } else if (ipenal2 == 1) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t l21 = cs_math_3_norm(d21);
        cs_real_t l23 = cs_math_3_norm(d23);
        cs_real_t l24 = cs_math_3_norm(d24);
        cs_real_t l21a = l21 - _imm_lgth_poro(x2, por2, x1, por1);
        cs_real_t l23a = l23 - _imm_lgth_poro(x2, por2, x3, por3);
        cs_real_t l24a = l24 - _imm_lgth_poro(x2, por2, x4, por4);

        cs_real_t lbd21 = l21a / l21;
        cs_real_t lbd23 = l23a / l23;
        cs_real_t lbd24 = l24a / l24;

        cs_real_t volp = *vol * lbd21 * lbd23 * lbd24;

        cs_real_3_t x21, x23, x24;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x21[i] = x2[i] + lbd21 * d21[i];
          x23[i] = x2[i] + lbd23 * d23[i];
          x24[i] = x2[i] + lbd24 * d24[i];
          cogp[i] = 0.25 * (x2[i] + x21[i] + x23[i] + x24[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * *vol - cogp[i] * volp;
        }

        *vol *= (1. - lbd21 * lbd23 * lbd24);

        return;

      } else if (ipenal3 == 1) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t l31 = cs_math_3_norm(d31);
        cs_real_t l32 = cs_math_3_norm(d32);
        cs_real_t l34 = cs_math_3_norm(d34);
        cs_real_t l31a = l31 - _imm_lgth_poro(x3, por3, x1, por1);
        cs_real_t l32a = l32 - _imm_lgth_poro(x3, por3, x2, por2);
        cs_real_t l34a = l34 - _imm_lgth_poro(x3, por3, x4, por4);

        cs_real_t lbd31 = l31a / l31;
        cs_real_t lbd32 = l32a / l32;
        cs_real_t lbd34 = l34a / l34;

        cs_real_t volp = *vol * lbd31 * lbd32 * lbd34;

        cs_real_3_t x31, x32, x34;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x31[i] = x3[i] + lbd31 * d31[i];
          x32[i] = x3[i] + lbd32 * d32[i];
          x34[i] = x3[i] + lbd34 * d34[i];
          cogp[i] = 0.25 * (x3[i] + x31[i] + x32[i] + x34[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * *vol - cogp[i] * volp;
        }

        *vol *= (1. - lbd31 * lbd32 * lbd34);

        return;

      } else if (ipenal4 == 1) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t l41 = cs_math_3_norm(d41);
        cs_real_t l42 = cs_math_3_norm(d42);
        cs_real_t l43 = cs_math_3_norm(d43);
        cs_real_t l41a = l41 - _imm_lgth_poro(x4, por4, x1, por1);
        cs_real_t l42a = l42 - _imm_lgth_poro(x4, por4, x2, por2);
        cs_real_t l43a = l43 - _imm_lgth_poro(x4, por4, x3, por3);

        cs_real_t lbd41 = l41a / l41;
        cs_real_t lbd42 = l42a / l42;
        cs_real_t lbd43 = l43a / l43;

        cs_real_t volp = *vol * lbd41 * lbd42 * lbd43;

        cs_real_3_t x41, x42, x43;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x41[i] = x4[i] + lbd41 * d41[i];
          x42[i] = x4[i] + lbd42 * d42[i];
          x43[i] = x4[i] + lbd43 * d43[i];
          cogp[i] = 0.25 * (x4[i] + x41[i] + x42[i] + x43[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * *vol - cogp[i] * volp;
        }

        *vol *= (1. - lbd41 * lbd42 * lbd43);

        return;

      } else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_poro\n");

    } else if (cpt4 == 3) {
      /* Three penalized points : tetrahedron from the non penalized
       * fourth point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      *vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 0) {
        cs_real_t l12 = cs_math_3_norm(d12);
        cs_real_t l13 = cs_math_3_norm(d13);
        cs_real_t l14 = cs_math_3_norm(d14);
        cs_real_t l12a = _imm_lgth_poro(x1, por1, x2, por2 );
        cs_real_t l13a = _imm_lgth_poro(x1, por1, x3, por3 );
        cs_real_t l14a = _imm_lgth_poro(x1, por1, x4, por4 );

        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd13 = l13a / l13;
        cs_real_t lbd14 = l14a / l14;

        *vol *= lbd12 * lbd13 * lbd14;

        cs_real_3_t x12, x13, x14;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x12[i] = x1[i] + lbd12 * d12[i];
          x13[i] = x1[i] + lbd13 * d13[i];
          x14[i] = x1[i] + lbd14 * d14[i];
          cogp[i] = 0.25 * (x1[i] + x12[i] + x13[i] + x14[i]);
          cog[i] += cogp[i] * *vol;
        }

        return;

      } else if (ipenal2 == 0) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t l21 = cs_math_3_norm(d21);
        cs_real_t l23 = cs_math_3_norm(d23);
        cs_real_t l24 = cs_math_3_norm(d24);
        cs_real_t l21a = _imm_lgth_poro(x2, por2, x1, por1);
        cs_real_t l23a = _imm_lgth_poro(x2, por2, x3, por3);
        cs_real_t l24a = _imm_lgth_poro(x2, por2, x4, por4);

        cs_real_t lbd21 = l21a / l21;
        cs_real_t lbd23 = l23a / l23;
        cs_real_t lbd24 = l24a / l24;

        *vol *= lbd21 * lbd23 * lbd24;

        cs_real_3_t x21, x23, x24;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x21[i] = x2[i] + lbd21 * d21[i];
          x23[i] = x2[i] + lbd23 * d23[i];
          x24[i] = x2[i] + lbd24 * d24[i];
          cogp[i] = 0.25 * (x2[i] + x21[i] + x23[i] + x24[i]);
          cog[i] += cogp[i] * *vol;
        }

        return;

      } else if (ipenal3 == 0) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t l31 = cs_math_3_norm(d31);
        cs_real_t l32 = cs_math_3_norm(d32);
        cs_real_t l34 = cs_math_3_norm(d34);
        cs_real_t l31a = _imm_lgth_poro(x3, por3, x1, por1);
        cs_real_t l32a = _imm_lgth_poro(x3, por3, x2, por2);
        cs_real_t l34a = _imm_lgth_poro(x3, por3, x4, por4);

        cs_real_t lbd31 = l31a / l31;
        cs_real_t lbd32 = l32a / l32;
        cs_real_t lbd34 = l34a / l34;

        *vol *= lbd31 * lbd32 * lbd34;

        cs_real_3_t x31, x32, x34;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x31[i] = x3[i] + lbd31 * d31[i];
          x32[i] = x3[i] + lbd32 * d32[i];
          x34[i] = x3[i] + lbd34 * d34[i];
          cogp[i] = 0.25 * (x3[i] + x31[i] + x32[i] + x34[i]);
          cog[i] += cogp[i] * *vol;
        }

        return;

      } else if (ipenal4 == 0) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t l41 = cs_math_3_norm(d41);
        cs_real_t l42 = cs_math_3_norm(d42);
        cs_real_t l43 = cs_math_3_norm(d43);
        cs_real_t l41a = _imm_lgth_poro(x4, por4, x1, por1);
        cs_real_t l42a = _imm_lgth_poro(x4, por4, x2, por2);
        cs_real_t l43a = _imm_lgth_poro(x4, por4, x3, por3);

        cs_real_t lbd41 = l41a / l41;
        cs_real_t lbd42 = l42a / l42;
        cs_real_t lbd43 = l43a / l43;

        *vol *= lbd41 * lbd42 * lbd43;

        cs_real_3_t x41, x42, x43;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x41[i] = x4[i] + lbd41 * d41[i];
          x42[i] = x4[i] + lbd42 * d42[i];
          x43[i] = x4[i] + lbd43 * d43[i];
          cogp[i] = 0.25 * (x4[i] + x41[i] + x42[i] + x43[i]);
          cog[i] += cogp[i] * *vol;
        }

        return;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    } else if (cpt4 == 2) {
      /* Two penalized points : more complex */
      /* One puts in y1 and y2 the two non penalized points and in y3 and y4
       * the penalized points */
      cs_real_3_t y1, y2, y3, y4;
      cs_real_t porr1 = por1;
      cs_real_t porr2 = por2;
      cs_real_t porr3 = por3;
      cs_real_t porr4 = por4;

      if (ipenal1 == 0 && ipenal2 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x2[i];
          y3[i] = x3[i];
          y4[i] = x4[i];
        }

        porr1 = por1;
        porr2 = por2;
        porr3 = por3;
        porr4 = por4;

      } else if (ipenal1 == 0 && ipenal3 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x3[i];
          y3[i] = x2[i];
          y4[i] = x4[i];
        }

        porr1 = por1;
        porr2 = por3;
        porr3 = por2;
        porr4 = por4;

      } else if (ipenal1 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x4[i];
          y3[i] = x2[i];
          y4[i] = x3[i];
        }

        porr1 = por1;
        porr2 = por4;
        porr3 = por2;
        porr4 = por3;

      } else if (ipenal2 == 0 && ipenal3 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x3[i];
          y3[i] = x1[i];
          y4[i] = x4[i];
        }

        porr1 = por2;
        porr2 = por3;
        porr3 = por1;
        porr4 = por4;

      } else if (ipenal2 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x3[i];
        }

        porr1 = por2;
        porr2 = por4;
        porr3 = por1;
        porr4 = por3;

      } else if (ipenal3 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x3[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x2[i];
        }

        porr1 = por3;
        porr2 = por4;
        porr3 = por1;
        porr4 = por2;

      } else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

      cs_real_3_t d13, d14, d23, d24;
      for (int i = 0; i < 3; i++) {
        d13[i] = y3[i] - y1[i];
        d14[i] = y4[i] - y1[i];
        d23[i] = y3[i] - y2[i];
        d24[i] = y4[i] - y2[i];
      }

      cs_real_t l13 = cs_math_3_norm(d13);
      cs_real_t l13a = _imm_lgth_poro(y1, porr1, y3, porr3);
      cs_real_t lbd13 = l13a / l13;

      cs_real_t l14 = cs_math_3_norm(d14);
      cs_real_t l14a = _imm_lgth_poro(y1, porr1, y4, porr4);
      cs_real_t lbd14 = l14a / l14;

      cs_real_t l23 = cs_math_3_norm(d23);
      cs_real_t l23a = _imm_lgth_poro(y2, porr2, y3, porr3);
      cs_real_t lbd23 = l23a / l23;

      cs_real_t l24 = cs_math_3_norm(d24);
      cs_real_t l24a = _imm_lgth_poro(y2, porr2, y4, porr4);
      cs_real_t lbd24 = l24a / l24;

      cs_real_3_t y13, y14, y23, y24;
      for (int i = 0; i < 3; i++) {
        y13[i] = y1[i] + lbd13 * d13[i];
        y14[i] = y1[i] + lbd14 * d14[i];
        y23[i] = y2[i] + lbd23 * d23[i];
        y24[i] = y2[i] + lbd24 * d24[i];
      }


      *vol = _prism_vol(y1, y2, y13, y23, y24, y14);

      cs_real_3_t cogp;
      for (int i = 0; i < 3; i++) {
        cogp[i] = (y1[i] + y2[i] + y13[i] + y23[i] + y24[i] + y14[i])
                  * cs_math_1ov6;
        cog[i] += cogp[i] * *vol;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by dichotomy the area of the immersed part of a triangle
 *          (x1,x2,x3) based on the input porosities at its vertices through
 *          a recursive approach.
 *
 * \param[in]  x1           point 1
 * \param[in]  por1         porosity at 1
 * \param[in]  x2           point 2
 * \param[in]  por2         porosity at 2
 * \param[in]  x3           point 3
 * \param[in]  por3         porosity at 3
 * \param[in]  icut         number of sub-cut for cells in cut-cells algorithm
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_tri_surf_trunc(cs_real_3_t x1,
                cs_real_t por1,
                cs_real_3_t x2,
                cs_real_t por2,
                cs_real_3_t x3,
                cs_real_t por3,
                int icut)
{
  /* Mean porosity */
  cs_real_t por4 = (por1 + por2 + por3) * cs_math_1ov3;

  /* Check if some vertices are considered as solid */
  int ipenal1 = 0;
  if (por1 < 0.5)
    ipenal1 = 1;

  int ipenal2 = 0;
  if (por2 < 0.5)
    ipenal2 = 1;

  int ipenal3 = 0;
  if (por3 < 0.5)
    ipenal3 = 1;

  int ipenal4 = 0;
  if (por4 < 0.5)
    ipenal4 = 1;

  int cpt3 = ipenal1 + ipenal2 + ipenal3;
  int cpt4 = cpt3 + ipenal4;

  cs_real_t surf = 0.;

  /* Check if the area is considered as completely solid, null are */
  if (cpt4 == 4)
    return surf;
  /* Check if no solid vertex, complete area */
  else if (cpt4 == 0) {

    cs_real_3_t dx12, dx13;
    cs_real_3_t pvec;
    for (int i = 0; i < 3; i++) {
      dx12[i] = x2[i] - x1[i];
      dx13[i] = x3[i] - x1[i];
    }

    pvec[0] = dx12[1]*dx13[2] - dx12[2]*dx13[1];
    pvec[1] = dx12[2]*dx13[0] - dx12[0]*dx13[2];
    pvec[2] = dx12[0]*dx13[1] - dx12[1]*dx13[0];

    surf = 0.5 * cs_math_3_norm(pvec);
    return surf;

  /* Check if at least one solid vertex, subdivision of the edges and
   * distribution of porosities before recursive call */
  } else if (icut > 0) {
    cs_real_3_t x12, x23, x13;
    icut--;

    for (int i = 0; i < 3; i++) {
      x12[i] = 0.5 * (x1[i] + x2[i]);
      x13[i] = 0.5 * (x1[i] + x3[i]);
      x23[i] = 0.5 * (x2[i] + x3[i]);
    }
    cs_real_t por12 = 0.5*(por1+por2);
    cs_real_t por13 = 0.5*(por1+por3);
    cs_real_t por23 = 0.5*(por2+por3);

    surf += _tri_surf_trunc(x1,  por1,  x12, por12, x13, por13, icut);
    surf += _tri_surf_trunc(x2,  por2,  x12, por12, x23, por23, icut);
    surf += _tri_surf_trunc(x3,  por3,  x13, por13, x23, por23, icut);
    surf += _tri_surf_trunc(x12, por12, x13, por13, x23, por23, icut);

    return surf;

  } else {
    if (cpt3 != 1 && cpt3 != 2)
      return surf;

    if (cpt3 == 1) {
      if (ipenal1 == 1) {

        cs_real_3_t dx12, dx13;
        for (int i = 0; i < 3; i++) {
          dx12[i] = x2[i] - x1[i];
          dx13[i] = x3[i] - x1[i];
        }

        cs_real_3_t pvec;
        pvec[0] = dx12[1]*dx13[2] - dx12[2]*dx13[1];
        pvec[1] = dx12[2]*dx13[0] - dx12[0]*dx13[2];
        pvec[2] = dx12[0]*dx13[1] - dx12[1]*dx13[0];

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t l12 = cs_math_3_norm(dx12);
        cs_real_t l13 = cs_math_3_norm(dx13);

        cs_real_t l12a = l12 - _imm_lgth_poro(x1, por1, x2, por2);
        cs_real_t l13a = l13 - _imm_lgth_poro(x1, por1, x3, por3);
        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd13 = l13a / l13;

        surf = surft * (1. - lbd12 * lbd13);

        return surf;

      } else if (ipenal2 == 1) {

        cs_real_3_t dx12, dx23;
        for (int i = 0; i < 3; i++) {
          dx12[i] = x2[i] - x1[i];
          dx23[i] = x3[i] - x2[i];
        }

        cs_real_3_t pvec;
        pvec[0] = dx12[1]*dx23[2] - dx12[2]*dx23[1];
        pvec[1] = dx12[2]*dx23[0] - dx12[0]*dx23[2];
        pvec[2] = dx12[0]*dx23[1] - dx12[1]*dx23[0];

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t l12 = cs_math_3_norm(dx12);
        cs_real_t l23 = cs_math_3_norm(dx23);

        cs_real_t l12a = l12 - _imm_lgth_poro(x1, por1, x2, por2);
        cs_real_t l23a = l23 - _imm_lgth_poro(x2, por2, x3, por3);
        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd23 = l23a / l23;

        surf = surft * (1. - lbd12 * lbd23);

        return surf;

      } else {

        cs_real_3_t dx13, dx23;
        for (int i = 0; i < 3; i++) {
          dx13[i] = x3[i] - x1[i];
          dx23[i] = x3[i] - x2[i];
        }

        cs_real_3_t pvec;
        pvec[0] = dx13[1]*dx23[2] - dx13[2]*dx23[1];
        pvec[1] = dx13[2]*dx23[0] - dx13[0]*dx23[2];
        pvec[2] = dx13[0]*dx23[1] - dx13[1]*dx23[0];

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t l13 = cs_math_3_norm(dx13);
        cs_real_t l23 = cs_math_3_norm(dx23);

        cs_real_t l13a = l13 - _imm_lgth_poro(x1, por1, x3, por3);
        cs_real_t l23a = l23 - _imm_lgth_poro(x2, por2, x3, por3);
        cs_real_t lbd13 = l13a / l13;
        cs_real_t lbd23 = l23a / l23;

        surf = surft * (1. - lbd13 * lbd23);

        return surf;
      }
    } else if (cpt3 == 2) {
      if (ipenal1 == 0) {
        cs_real_t l12 = _imm_lgth_poro(x1, por1, x2, por2);
        cs_real_t l13 = _imm_lgth_poro(x1, por1, x3, por3);

        cs_real_3_t dx12;
        for (int i = 0; i < 3; i++)
          dx12[i] = x2[i] - x1[i];

        cs_real_t l12t = cs_math_3_norm(dx12);
        cs_real_t lamb12 = l12 / l12t;
        lamb12 = cs_math_fmax(lamb12, 0.);
        lamb12 = cs_math_fmin(lamb12, 1.);

        for (int i = 0; i < 3; i++)
          dx12[i] *= lamb12;

        cs_real_3_t dx13;
        for (int i = 0; i < 3; i++)
          dx13[i] = x3[i] - x1[i];

        cs_real_t l13t = cs_math_3_norm(dx13);
        cs_real_t lamb13 = l13 / l13t;
        lamb13 = cs_math_fmax(lamb13, 0.);
        lamb13 = cs_math_fmin(lamb13, 1.);

        for (int i = 0; i < 3; i++)
          dx13[i] *= lamb13;

        cs_real_3_t pvec;
        pvec[0] = dx12[1]*dx13[2] - dx12[2]*dx13[1];
        pvec[1] = dx12[2]*dx13[0] - dx12[0]*dx13[2];
        pvec[2] = dx12[0]*dx13[1] - dx12[1]*dx13[0];

        surf = 0.5 * cs_math_3_norm(pvec);

        return surf;

      } else if (ipenal2 == 0) {
        cs_real_t l12 = _imm_lgth_poro(x1, por1, x2, por2);
        cs_real_t l23 = _imm_lgth_poro(x2, por2, x3, por3);

        cs_real_3_t dx12;
        for (int i =  0; i < 3; i++)
          dx12[i] = x2[i] - x1[i];

        cs_real_t l12t = cs_math_3_norm(dx12);
        cs_real_t lamb12 = l12 / l12t;
        lamb12 = cs_math_fmax(lamb12, 0.);
        lamb12 = cs_math_fmin(lamb12, 1.);

        for (int i =  0; i < 3; i++)
          dx12[i] *= lamb12;

        cs_real_3_t dx23;
        for (int i =  0; i < 3; i++)
          dx23[i] = x3[i] - x2[i];

        cs_real_t l23t = cs_math_3_norm(dx23);
        cs_real_t lamb23 = l23 / l23t;
        lamb23 = cs_math_fmax(lamb23, 0.);
        lamb23 = cs_math_fmin(lamb23, 1.);
        for (int i =  0; i < 3; i++)
          dx23[i] *= lamb23;

        cs_real_3_t pvec;
        pvec[0] = dx12[1]*dx23[2] - dx12[2]*dx23[1];
        pvec[1] = dx12[2]*dx23[0] - dx12[0]*dx23[2];
        pvec[2] = dx12[0]*dx23[1] - dx12[1]*dx23[0];

        surf = 0.5 * cs_math_3_norm(pvec);

        return surf;

      } else {
        cs_real_t l13 = _imm_lgth_poro(x1, por1, x3, por3);
        cs_real_t l23 = _imm_lgth_poro(x2, por2, x3, por3);

        cs_real_3_t dx13;
        for (int i = 0; i < 3; i++)
          dx13[i] = x3[i] - x1[i];

        cs_real_t l13t = cs_math_3_norm(dx13);
        cs_real_t lamb13 = l13 / l13t;
        lamb13 = cs_math_fmax(lamb13, 0.);
        lamb13 = cs_math_fmin(lamb13, 1.);

        for (int i = 0; i < 3; i++)
          dx13[i] *= lamb13;

        cs_real_3_t dx23;
        for (int i = 0; i < 3; i++)
          dx23[i] = x3[i] - x2[i];

        cs_real_t l23t = cs_math_3_norm(dx23);
        cs_real_t lamb23 = l23 / l23t;
        lamb23 = cs_math_fmax(lamb23, 0.);
        lamb23 = cs_math_fmin(lamb23, 1.);

        for (int i = 0; i < 3; i++)
          dx23[i] *=lamb23;

        cs_real_3_t pvec;
        pvec[0] = dx13[1]*dx23[2] - dx13[2]*dx23[1];
        pvec[1] = dx13[2]*dx23[0] - dx13[0]*dx23[2];
        pvec[2] = dx13[0]*dx23[1] - dx13[1]*dx23[0];

        surf = 0.5 * cs_math_3_norm(pvec);

        return surf;
      }
    }
  }
  return surf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute (using bisection) the volume of the immersed part of a
 *        tetrahedron described by its vertices (x1,x2,x3,x4) based on the
 *        cut-cell method and recursive approach.
 *
 * \param[out] vol          volume
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 * \param[in]  t            time value for the current time step
 * \param[in]  icut         number of bisections
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

static void
_tetra_vol_cutcell(cs_real_t    *vol,
                   cs_real_3_t   x1,
                   cs_real_3_t   x2,
                   cs_real_3_t   x3,
                   cs_real_3_t   x4,
                   cs_real_t     t,
                   int           icut,
                   int           num_object)
{
  /*
   * One subdivision level is enough
   * (triangles based on cell and edge centers). */

  cs_real_3_t cog1234;
  for (int i = 0; i < 3; i++)
    cog1234[i] = (x1[i] + x2[i] + x3[i] + x4[i]) * 0.25;

  /* Check if some vertices are considered as solid */
  int ipenal1 = cs_ibm_object_compute_cut_porosity(1, x1, t, num_object);

  int ipenal2 = cs_ibm_object_compute_cut_porosity(1, x2, t, num_object);

  int ipenal3 = cs_ibm_object_compute_cut_porosity(1, x3, t, num_object);

  int ipenal4 = cs_ibm_object_compute_cut_porosity(1, x4, t, num_object);

  int ipenalc = cs_ibm_object_compute_cut_porosity(1, cog1234, t, num_object);

  int cpt4 = ipenal1 + ipenal2 + ipenal3 + ipenal4;
  int cpt5 = cpt4 + ipenalc;

  *vol = 0.;

  /* Check if the volume is considered as completely solid, null volume */
  if (cpt5 == 5)
    return;
  /* Check if no solid vertex, complete volume */
  else if (cpt5 == 0 && icut > 0) {

    *vol = _tetra_vol(x1, x2, x3, x4);

    return;

  /* Check if at least one solid vertex, subdivision of the edges and
   * distribution of porosities before recursive call */
  }
  else if (icut > 0) {
    icut--;

    cs_real_3_t x12, x13, x14, x23, x24, x34;
    for (int i = 0; i < 3; i++) {
      x12[i] = 0.5 * (x1[i] + x2[i]);
      x13[i] = 0.5 * (x1[i] + x3[i]);
      x14[i] = 0.5 * (x1[i] + x4[i]);
      x23[i] = 0.5 * (x2[i] + x3[i]);
      x24[i] = 0.5 * (x2[i] + x4[i]);
      x34[i] = 0.5 * (x3[i] + x4[i]);
    }

    cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
    _tetra_vol_cutcell(&vol1, x1,  x12, x13, x14, t, icut, num_object);
    _tetra_vol_cutcell(&vol2, x2,  x23, x12, x24, t, icut, num_object);
    _tetra_vol_cutcell(&vol3, x3,  x13, x23, x34, t, icut, num_object);
    _tetra_vol_cutcell(&vol4, x4,  x14, x24, x34, t, icut, num_object);
    _tetra_vol_cutcell(&vol5, x12, x13, x34, x23, t, icut, num_object);
    _tetra_vol_cutcell(&vol6, x12, x13, x34, x14, t, icut, num_object);
    _tetra_vol_cutcell(&vol7, x34, x24, x12, x23, t, icut, num_object);
    _tetra_vol_cutcell(&vol8, x34, x24, x12, x14, t, icut, num_object);

    *vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
    return;

  } else {
    if (cpt4 == 0) {
      /* No penalized point -> full tetrahedron */

      *vol = _tetra_vol(x1, x2, x3, x4);

    } else if (cpt4 == 1) {
      /* One penalized point tetrahedron - tetra from the penalized point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      *vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 1) {
        cs_real_t l12 = cs_math_3_norm(d12);
        cs_real_t l13 = cs_math_3_norm(d13);
        cs_real_t l14 = cs_math_3_norm(d14);
        cs_real_t l12a = l12 - _imm_lgth_cutcell(x1, x2, t, num_object);
        cs_real_t l13a = l13 - _imm_lgth_cutcell(x1, x3, t, num_object);
        cs_real_t l14a = l14 - _imm_lgth_cutcell(x1, x4, t, num_object);

        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd13 = l13a / l13;
        cs_real_t lbd14 = l14a / l14;

        *vol *= (1. - lbd12 * lbd13 * lbd14);

        return;

      } else if (ipenal2 == 1) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t l21 = cs_math_3_norm(d21);
        cs_real_t l23 = cs_math_3_norm(d23);
        cs_real_t l24 = cs_math_3_norm(d24);
        cs_real_t l21a = l21 - _imm_lgth_cutcell(x2, x1, t, num_object);
        cs_real_t l23a = l23 - _imm_lgth_cutcell(x2, x3, t, num_object);
        cs_real_t l24a = l24 - _imm_lgth_cutcell(x2, x4, t, num_object);

        cs_real_t lbd21 = l21a / l21;
        cs_real_t lbd23 = l23a / l23;
        cs_real_t lbd24 = l24a / l24;

        *vol *= (1. - lbd21 * lbd23 * lbd24);

        return;

      } else if (ipenal3 == 1) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t l31 = cs_math_3_norm(d31);
        cs_real_t l32 = cs_math_3_norm(d32);
        cs_real_t l34 = cs_math_3_norm(d34);
        cs_real_t l31a = l31 - _imm_lgth_cutcell(x3, x1, t, num_object);
        cs_real_t l32a = l32 - _imm_lgth_cutcell(x3, x2, t, num_object);
        cs_real_t l34a = l34 - _imm_lgth_cutcell(x3, x4, t, num_object);

        cs_real_t lbd31 = l31a / l31;
        cs_real_t lbd32 = l32a / l32;
        cs_real_t lbd34 = l34a / l34;

        *vol *= (1. - lbd31 * lbd32 * lbd34);

        return;

      } else if (ipenal4 == 1) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t l41 = cs_math_3_norm(d41);
        cs_real_t l42 = cs_math_3_norm(d42);
        cs_real_t l43 = cs_math_3_norm(d43);
        cs_real_t l41a = l41 - _imm_lgth_cutcell(x4, x1, t, num_object);
        cs_real_t l42a = l42 - _imm_lgth_cutcell(x4, x2, t, num_object);
        cs_real_t l43a = l43 - _imm_lgth_cutcell(x4, x3, t, num_object);

        cs_real_t lbd41 = l41a / l41;
        cs_real_t lbd42 = l42a / l42;
        cs_real_t lbd43 = l43a / l43;

        *vol *= (1. - lbd41 * lbd42 * lbd43);

        return;

      } else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    } else if (cpt4 == 3) {
      /* Three penalized points : tetrahedron from the non penalized
       * fourth point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      *vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 0) {
        cs_real_t l12 = cs_math_3_norm(d12);
        cs_real_t l13 = cs_math_3_norm(d13);
        cs_real_t l14 = cs_math_3_norm(d14);
        cs_real_t l12a = _imm_lgth_cutcell(x1, x2, t, num_object);
        cs_real_t l13a = _imm_lgth_cutcell(x1, x3, t, num_object);
        cs_real_t l14a = _imm_lgth_cutcell(x1, x4, t, num_object);

        cs_real_t lbd12 = l12a / l12;
        cs_real_t lbd13 = l13a / l13;
        cs_real_t lbd14 = l14a / l14;

        *vol *= lbd12 * lbd13 * lbd14;

        return;

      } else if (ipenal2 == 0) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t l21 = cs_math_3_norm(d21);
        cs_real_t l23 = cs_math_3_norm(d23);
        cs_real_t l24 = cs_math_3_norm(d24);
        cs_real_t l21a = _imm_lgth_cutcell(x2, x1, t, num_object);
        cs_real_t l23a = _imm_lgth_cutcell(x2, x3, t, num_object);
        cs_real_t l24a = _imm_lgth_cutcell(x2, x4, t, num_object);

        cs_real_t lbd21 = l21a / l21;
        cs_real_t lbd23 = l23a / l23;
        cs_real_t lbd24 = l24a / l24;

        *vol *= lbd21 * lbd23 * lbd24;

        return;

      } else if (ipenal3 == 0) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t l31 = cs_math_3_norm(d31);
        cs_real_t l32 = cs_math_3_norm(d32);
        cs_real_t l34 = cs_math_3_norm(d34);
        cs_real_t l31a = _imm_lgth_cutcell(x3, x1, t, num_object);
        cs_real_t l32a = _imm_lgth_cutcell(x3, x2, t, num_object);
        cs_real_t l34a = _imm_lgth_cutcell(x3, x4, t, num_object);

        cs_real_t lbd31 = l31a / l31;
        cs_real_t lbd32 = l32a / l32;
        cs_real_t lbd34 = l34a / l34;

        *vol *= lbd31 * lbd32 * lbd34;

        return;

      } else if (ipenal4 == 0) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t l41 = cs_math_3_norm(d41);
        cs_real_t l42 = cs_math_3_norm(d42);
        cs_real_t l43 = cs_math_3_norm(d43);
        cs_real_t l41a = _imm_lgth_cutcell(x4, x1, t, num_object);
        cs_real_t l42a = _imm_lgth_cutcell(x4, x2, t, num_object);
        cs_real_t l43a = _imm_lgth_cutcell(x4, x3, t, num_object);

        cs_real_t lbd41 = l41a / l41;
        cs_real_t lbd42 = l42a / l42;
        cs_real_t lbd43 = l43a / l43;

        *vol *= lbd41 * lbd42 * lbd43;

        return;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    } else if (cpt4 == 2) {
      /* Two penalized points : more complex */
      /* One puts in y1 et y2 the two non penalized points and in y3 and y4
       * the penalized points */
      cs_real_3_t y1, y2, y3, y4;
      if (ipenal1 == 0 && ipenal2 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x2[i];
          y3[i] = x3[i];
          y4[i] = x4[i];
        }
      else if (ipenal1 == 0 && ipenal3 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x3[i];
          y3[i] = x2[i];
          y4[i] = x4[i];
        }
      else if (ipenal1 == 0 && ipenal4 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x4[i];
          y3[i] = x2[i];
          y4[i] = x3[i];
        }
      else if (ipenal2 == 0 && ipenal3 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x3[i];
          y3[i] = x1[i];
          y4[i] = x4[i];
        }
      else if (ipenal2 == 0 && ipenal4 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x3[i];
        }
      else if (ipenal3 == 0 && ipenal4 == 0)
        for (int i = 0; i < 3; i++) {
          y1[i] = x3[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x2[i];
        }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

      cs_real_3_t d13, d14, d23, d24;
      for (int i = 0; i < 3; i++) {
        d13[i] = y3[i] - y1[i];
        d14[i] = y4[i] - y1[i];
        d23[i] = y3[i] - y2[i];
        d24[i] = y4[i] - y2[i];
      }

      cs_real_t l13 = cs_math_3_norm(d13);
      cs_real_t l13a = _imm_lgth_cutcell(y1, y3, t, num_object);
      cs_real_t lbd13 = l13a / l13;

      cs_real_t l14 = cs_math_3_norm(d14);
      cs_real_t l14a = _imm_lgth_cutcell(y1, y4, t, num_object);
      cs_real_t lbd14 = l14a / l14;

      cs_real_t l23 = cs_math_3_norm(d23);
      cs_real_t l23a = _imm_lgth_cutcell(y2, y3, t, num_object);
      cs_real_t lbd23 = l23a / l23;

      cs_real_t l24 = cs_math_3_norm(d24);
      cs_real_t l24a = _imm_lgth_cutcell(y2, y4, t, num_object);
      cs_real_t lbd24 = l24a / l24;

      cs_real_3_t y13, y14, y23, y24;
      for (int i = 0; i < 3; i++) {
        y13[i] = y1[i] + lbd13 * d13[i];
        y14[i] = y1[i] + lbd14 * d14[i];
        y23[i] = y2[i] + lbd23 * d23[i];
        y24[i] = y2[i] + lbd24 * d24[i];
      }

      *vol = _prism_vol(y1, y2, y13, y23, y24, y14);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Call to user function to know whether we are in the solid or not
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

static int
_penal_glob(const cs_lnum_t   c_id,
            const cs_real_3_t xyz,
            const cs_real_t   t,
            const int         num_object)
{
  /* Call to user function to know whether we are in the solid or not */
  int ipenal = cs_ibm_object_compute_cut_porosity(c_id, xyz, t, num_object);

  return ipenal;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell porosity using the cut cell method
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_cut_porosity(const cs_mesh_t *mesh,
                           const cs_mesh_quantities_t *mesh_quantities,
                           int  *comp_cell)
{

  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  cs_real_t *cell_vol  = mesh_quantities->cell_vol;

  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)mesh_quantities->cell_cen;
  const cs_real_3_t *i_face_cog
    = (const cs_real_3_t *)mesh_quantities->i_face_cog;
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)mesh_quantities->b_face_cog;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  cs_real_t t_cur = cs_glob_time_step->t_cur;
  int icut = cs_ibm->nb_cut_cells;
  cs_real_t voltot = 0;

  int *nbvtx, *nbvtx_in, *cog_in, *node_in;
  BFT_MALLOC(nbvtx, n_cells_ext, int);
  BFT_MALLOC(nbvtx_in, n_cells_ext, int);
  BFT_MALLOC(cog_in, n_cells_ext, int);
  BFT_MALLOC(node_in, mesh->n_vertices, int);
  memset(nbvtx, 0, n_cells_ext*sizeof(int));
  memset(nbvtx_in, 0, n_cells_ext*sizeof(int));
  memset(cog_in, 0, n_cells_ext*sizeof(int));

  for (int v_id = 0; v_id < mesh->n_vertices; v_id++)
    node_in[v_id] = -1;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (comp_cell[c_id0] + comp_cell[c_id1] > 0) {

      int num_objecti = -1;
      int num_objectj = -1;
      int num_objectij = -1;

      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
          vtx_id < i_face_vtx_idx[f_id + 1]; vtx_id++) {
        cs_lnum_t v_id = i_face_vtx[vtx_id];

        cs_real_3_t xn;
        for (int idim = 0; idim < 3; idim++)
          xn[idim] = vtx_crd[v_id][idim];

        nbvtx[c_id0]++;
        nbvtx[c_id1]++;

        if (node_in[v_id] < 0)
          node_in[v_id] = _penal_glob(c_id0, xn, t_cur, num_objectij);
        int in = node_in[v_id];

        nbvtx_in[c_id0] += in;
        nbvtx_in[c_id1] += in;
      }

      cs_real_3_t xif, xjf;
      for (int idim = 0; idim < 3; idim++) {
        xif[idim] = 0.5 * (i_face_cog[f_id][idim] + cell_cen[c_id0][idim]);
        xjf[idim] = 0.5 * (i_face_cog[f_id][idim] + cell_cen[c_id1][idim]);
      }

      if (cog_in[c_id0] == 0)
        cog_in[c_id0] += _penal_glob(c_id0, xif, t_cur, num_objecti);
      if (cog_in[c_id1] == 0)
        cog_in[c_id1] += _penal_glob(c_id1, xjf, t_cur, num_objectj);
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    if (comp_cell[c_id] > 0) {

      int num_objecti = -1;

      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
          vtx_id < b_face_vtx_idx[f_id + 1]; vtx_id++) {
        cs_lnum_t v_id = b_face_vtx[vtx_id];

        cs_real_3_t xn;
        for (int idim = 0; idim < 3; idim++)
          xn[idim] = vtx_crd[v_id][idim];

        nbvtx[c_id]++;

        if (node_in[v_id] < 0)
          node_in[v_id] = _penal_glob(c_id, xn, t_cur, num_objecti);

        int in = node_in[v_id];
        nbvtx_in[c_id] += in;
      }

      cs_real_3_t xif;
      for (int idim = 0; idim < 3; idim++)
        xif[idim] = 0.5 * (b_face_cog[f_id][idim] + cell_cen[c_id][idim]);

      if (cog_in[c_id] == 0)
        cog_in[c_id] += _penal_glob(c_id, xif, t_cur, num_objecti);
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (cog_in[c_id] == 0) {
      int num_object = -1;

      cog_in[c_id] += _penal_glob(c_id, cell_cen[c_id], t_cur, num_object);
    }

  cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, nbvtx);
  cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, nbvtx_in);
  cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, cog_in);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    if (comp_cell[c_id] > 0)
      CS_F_(poro)->val[c_id] = 0.;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (comp_cell[c_id] > 0)
      if (nbvtx_in[c_id] == 0 && cog_in[c_id] == 0)
        CS_F_(poro)->val[c_id] = 1.;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    bool computi = true;
    if (comp_cell[c_id0] == 0)
      computi = false;

    if ((nbvtx[c_id0] == nbvtx_in[c_id0] && cog_in[c_id0] > 0)
         || (nbvtx_in[c_id0] == 0 && cog_in[c_id0] == 0))
      computi = false;

    bool computj = true;
    if (comp_cell[c_id1] == 0)
      computj = false;

    if ((nbvtx[c_id1] == nbvtx_in[c_id1] && cog_in[c_id1] > 0)
         || (nbvtx_in[c_id1] == 0 && cog_in[c_id1] == 0))
      computj = false;

    if (computi || computj) {
      cs_real_3_t xf,xi,xj;
      for (int i = 0; i < 3; i++) {
        xf[i] = i_face_cog[f_id][i];
        xi[i] = cell_cen[c_id0][i];
        xj[i] = cell_cen[c_id1][i];
      }

      int num_objecti = -1;
      int num_objectj = -1;

      /* For each face, one builds the tetrahedra base on an edge,
       * the cog of the face and the cog of the cell (i or j) */
      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
           vtx_id < i_face_vtx_idx[f_id + 1] - 1; vtx_id++) {
        cs_lnum_t vtx_id1 = i_face_vtx[vtx_id];
        cs_lnum_t vtx_id2 = i_face_vtx[vtx_id + 1];

        cs_real_3_t x1, x2;
        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
        }

        if (computi) {
          _tetra_vol_cutcell(&voltot, x1, x2, xf, xi, t_cur,
                            icut, num_objecti);
          CS_F_(poro)->val[c_id0] += voltot;
        }

        if (computj) {
          _tetra_vol_cutcell(&voltot, x1, x2, xf, xj, t_cur,
                            icut, num_objectj);
          CS_F_(poro)->val[c_id1] += voltot;
        }
      }

      cs_lnum_t vtx_id1
        = i_face_vtx[i_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2
        = i_face_vtx[i_face_vtx_idx[f_id]];

      cs_real_3_t x1, x2;
      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
      }

      if (computi) {
        _tetra_vol_cutcell(&voltot, x1, x2, xf, xi, t_cur, icut, num_objecti);
        CS_F_(poro)->val[c_id0] += voltot;
      }

      if (computj) {
        _tetra_vol_cutcell(&voltot, x1, x2, xf, xj, t_cur, icut, num_objectj);
        CS_F_(poro)->val[c_id1] += voltot;
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    bool computi = true;
    if (comp_cell[c_id] == 0)
      computi = false;

    if ((     nbvtx[c_id] == nbvtx_in[c_id] && cog_in[c_id] > 0)
         || ( nbvtx_in[c_id] == 0 && cog_in[c_id] == 0))
      computi = false;

    if (computi) {
      cs_real_3_t xf, xi;
      for (int i = 0; i < 3; i++) {
        xf[i] = b_face_cog[f_id][i];
        xi[i] = cell_cen[c_id][i];
      }

      int num_objecti = -1;

      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
           vtx_id < b_face_vtx_idx[f_id + 1] - 1; vtx_id++) {
        cs_lnum_t vtx_id1 = b_face_vtx[vtx_id];
        cs_lnum_t vtx_id2 = b_face_vtx[vtx_id + 1];

        cs_real_3_t x1, x2;
        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
        }

        _tetra_vol_cutcell(&voltot, x1, x2, xf, xi, t_cur, icut, num_objecti);
        CS_F_(poro)->val[c_id] += voltot;
      }

      cs_lnum_t vtx_id1
        = b_face_vtx[b_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2
        = b_face_vtx[b_face_vtx_idx[f_id]];

      cs_real_3_t x1, x2;
      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
      }

      _tetra_vol_cutcell(&voltot, x1, x2, xf, xi, t_cur, icut, num_objecti);
      CS_F_(poro)->val[c_id] += voltot;
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (comp_cell[c_id] > 0)
      if (nbvtx_in[c_id] != 0 || cog_in[c_id] != 0) {
        CS_F_(poro)->val[c_id] /= cell_vol[c_id];
        CS_F_(poro)->val[c_id] = cs_math_fmax(CS_F_(poro)->val[c_id], 0.);
        CS_F_(poro)->val[c_id] = cs_math_fmin(CS_F_(poro)->val[c_id], 1.);

        if (CS_F_(poro)->val[c_id] < 1.e-5)
          CS_F_(poro)->val[c_id] = 0.;
      }
  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(poro)->val);

  BFT_FREE(cog_in);
  BFT_FREE(nbvtx);
  BFT_FREE(nbvtx_in);
  BFT_FREE(node_in);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell cog (and cell porosity from porosity at vertices)
 *         Cut method into sub-tetras
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   por_vtx            vertex porosity
 * \param[in]   por_init           initialization porosity
 * \param[in]   comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_cog(const cs_mesh_t            *mesh,
                  const cs_mesh_quantities_t *mesh_quantities,
                  cs_real_t                  *v_poro,
                  cs_real_t                  *por_init,
                  int                        *comp_cell)
{
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  cs_real_t *cell_vol  = mesh_quantities->cell_vol;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mesh_quantities->cell_cen;
  cs_real_3_t *cell_f_cen = (cs_real_3_t *)mesh_quantities->cell_f_cen;
  cs_real_3_t *cell_s_cen = (cs_real_3_t *)mesh_quantities->cell_s_cen;
  const cs_real_3_t *i_face_cog
    = (const cs_real_3_t *)mesh_quantities->i_face_cog;
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)mesh_quantities->b_face_cog;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  int nt_cur = cs_glob_time_step->nt_cur;
  int nt_prev = cs_glob_time_step->nt_prev;
  bool comp_all_cog = true;
  if (!cs_restart_present() && nt_cur > 1)
    comp_all_cog = false;
  if (cs_restart_present() && nt_cur != nt_prev)
    comp_all_cog = false;

  int icut = cs_ibm->nb_cut_cells;

  cs_real_t *c_poro;
  BFT_MALLOC(c_poro, n_cells_ext, cs_real_t);

  cs_vertex_to_cell(CS_VERTEX_TO_CELL_SHEPARD, 0, 1, NULL,
                    v_poro, c_poro);

  cs_real_t voltot = 0.;

  cs_real_t *porbis;
  BFT_MALLOC(porbis, n_cells_ext, cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    if (comp_cell[c_id] > 0 || comp_all_cog) {
      porbis[c_id] = 0.;
      for (int i = 0; i < 3; i++) {
        cell_f_cen[c_id][i] = 0.;
        cell_s_cen[c_id][i] = 0.;
      }
    }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    bool computi = true;
    bool computj = true;

    if (!comp_all_cog) {
      if (comp_cell[c_id0] == 0)
        computi = false;
      if (comp_cell[c_id1] == 0)
        computj = false;
    }

    if (computi || computj) {
      cs_real_3_t xf,xi,xj;
      for (int i = 0; i < 3; i++) {
        xf[i] = i_face_cog[f_id][i];
        xi[i] = cell_cen[c_id0][i];
        xj[i] = cell_cen[c_id1][i];
      }

      cs_real_t pori = c_poro[c_id0];
      cs_real_t porj = c_poro[c_id1];

      cs_real_t porf = 0.;
      cs_real_t cpt = 0.;
      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
           vtx_id < i_face_vtx_idx[f_id + 1]; vtx_id++) {
         cs_lnum_t inod = i_face_vtx[vtx_id];
         porf += v_poro[inod];
         cpt += 1;
      }
      porf /= cpt;

      /* For each face, one builds the tetrahedra base on an edge,
       * the cog of the face and the cog of the cell (i or j) */
      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
           vtx_id < i_face_vtx_idx[f_id + 1] - 1; vtx_id++) {
        cs_lnum_t vtx_id1 = i_face_vtx[vtx_id];
        cs_lnum_t vtx_id2 = i_face_vtx[vtx_id + 1];

        cs_real_3_t x1, x2;
        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
        }

        cs_real_t por1 = v_poro[vtx_id1];
        cs_real_t por2 = v_poro[vtx_id2];

        if (computi) {
          _tetra_vol_poro(&voltot, cell_f_cen[c_id0], x1, por1,
                          x2, por2, xf, porf, xi, pori, icut);
          porbis[c_id0] += voltot;
          _tetra_vol_poro(&voltot, cell_s_cen[c_id0], x1, 1.-por1,
                          x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);
        }

        if (computj) {
          _tetra_vol_poro(&voltot, cell_f_cen[c_id1], x1, por1,
                          x2, por2, xf, porf, xj, porj, icut);
          porbis[c_id1] += voltot;
          _tetra_vol_poro(&voltot, cell_s_cen[c_id1], x1, 1.-por1,
                          x2, 1.-por2, xf, 1.-porf, xj, 1.-porj, icut);
        }
      }

      cs_lnum_t vtx_id1 = i_face_vtx[i_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2 = i_face_vtx[i_face_vtx_idx[f_id]];

      cs_real_3_t x1, x2;
      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
      }

      cs_real_t por1 = v_poro[vtx_id1];
      cs_real_t por2 = v_poro[vtx_id2];

      if (computi) {
        _tetra_vol_poro(&voltot, cell_f_cen[c_id0], x1, por1,
                        x2, por2, xf, porf, xi, pori, icut);
        porbis[c_id0] += voltot;
        _tetra_vol_poro(&voltot, cell_s_cen[c_id0], x1, 1.-por1,
                        x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);
      }

      if (computj) {
        _tetra_vol_poro(&voltot, cell_f_cen[c_id1], x1, por1,
                        x2, por2, xf, porf, xj, porj, icut);
        porbis[c_id1] += voltot;
        _tetra_vol_poro(&voltot, cell_s_cen[c_id1], x1, 1.-por1,
                        x2, 1.-por2, xf, 1.-porf, xj, 1.-porj, icut);
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    bool computi = true;

    if (!comp_all_cog)
      if (comp_cell[c_id] == 0)
        computi = false;

    if (computi) {
      cs_real_3_t xf, xi;
      for (int i = 0; i < 3; i++) {
        xf[i] = b_face_cog[f_id][i];
        xi[i] = cell_cen[c_id][i];
      }

      cs_real_t pori = c_poro[c_id];

      cs_real_t porf = 0.;
      cs_real_t cpt = 0.;
      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
           vtx_id < b_face_vtx_idx[f_id + 1]; vtx_id++) {
         cs_lnum_t inod = b_face_vtx[vtx_id];
         porf += v_poro[inod];
         cpt += 1;
      }
      porf /= cpt;

      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
           vtx_id < b_face_vtx_idx[f_id + 1] - 1; vtx_id++) {
        cs_lnum_t vtx_id1 = b_face_vtx[vtx_id];
        cs_lnum_t vtx_id2 = b_face_vtx[vtx_id + 1];

        cs_real_3_t x1, x2;
        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
        }

        cs_real_t por1 = v_poro[vtx_id1];
        cs_real_t por2 = v_poro[vtx_id2];

        _tetra_vol_poro(&voltot, cell_f_cen[c_id], x1, por1,
                        x2, por2, xf, porf, xi, pori, icut);
        porbis[c_id] += voltot;
        _tetra_vol_poro(&voltot, cell_s_cen[c_id], x1, 1.-por1,
                        x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);
      }

      cs_lnum_t vtx_id1 = b_face_vtx[b_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2 = b_face_vtx[b_face_vtx_idx[f_id]];

      cs_real_3_t x1, x2;
      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
      }

      cs_real_t por1 = v_poro[vtx_id1];
      cs_real_t por2 = v_poro[vtx_id2];

      _tetra_vol_poro(&voltot, cell_f_cen[c_id], x1, por1,
                      x2, por2, xf, porf, xi, pori, icut);
      porbis[c_id] += voltot;
      _tetra_vol_poro(&voltot, cell_s_cen[c_id], x1, 1.-por1,
                      x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);

    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (comp_cell[c_id] > 0 || comp_all_cog) {
      for (int i = 0; i < 3; i++)
        cell_f_cen[c_id][i] /= cs_math_fmax(porbis[c_id],
                                            cs_math_epzero * cell_vol[c_id]);
      for (int i = 0; i < 3; i++)
        cell_s_cen[c_id][i] /= cs_math_fmax(cell_vol[c_id]-porbis[c_id],
                                            cs_math_epzero * cell_vol[c_id]);

      porbis[c_id] /= cell_vol[c_id];
      porbis[c_id] = cs_math_fmax( porbis[c_id], 0.);
      porbis[c_id] = cs_math_fmin( porbis[c_id], 1.);

      if (porbis[c_id] < 1.e-5) {
        porbis[c_id] = 0.;
        for (int i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = cell_cen[c_id][i];
      }

      if (porbis[c_id] > 1. - 1.e-5) {
        for (int i = 0; i < 3; i++)
          cell_s_cen[c_id][i] = cell_cen[c_id][i];
      }

      if (por_init[c_id] < 1.e-5)
        for (int i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = cell_cen[c_id][i];
    }

  /* Recompute cognew for cells with almost null porosities */
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_f_cen, 3);
  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, porbis);

  if (mesh->n_init_perio > 0)
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_STANDARD,
                              (cs_real_t *)cell_f_cen);

  cs_lnum_t size_weight = n_cells_ext;
  if (mesh->n_vertices > size_weight)
    size_weight = mesh->n_vertices;

  cs_real_t *weight;
  BFT_MALLOC(weight, size_weight, cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < mesh->n_vertices; c_id++)
    weight[c_id] = 1.;

  for (cs_lnum_t c_id = 0.; c_id < n_cells_ext; c_id++)
    if (comp_cell[c_id] > 0 || comp_all_cog) {
      weight[c_id] = 0.;

      if (porbis[c_id] < 1.e-5)
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id][idir] = 0.;
    }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (comp_cell[c_id0] > 0 || comp_all_cog)
      if (porbis[c_id0] < 1.e-5 && porbis[c_id1] >= 1.e-5) {
        weight[c_id0] += porbis[c_id1];
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id0][idir] += porbis[c_id1] * i_face_cog[f_id][idir];
      }

    if (comp_cell[c_id1] > 0 || comp_all_cog)
      if (porbis[c_id1] < 1.e-5 && porbis[c_id0] >= 1.e-5) {
        weight[c_id1] += porbis[c_id0];
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id1][idir] += porbis[c_id0] * i_face_cog[f_id][idir];
      }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (comp_cell[c_id] > 0 || comp_all_cog)
      if (porbis[c_id] < 1.e-5) {
        if (weight[c_id] > 0.01)
          for (int idir = 0; idir < 3; idir++)
            cell_f_cen[c_id][idir] /= weight[c_id];
        else
          for (int idir = 0; idir < 3; idir++)
            cell_f_cen[c_id][idir] = cell_cen[c_id][idir];
      }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (comp_cell[c_id] > 0 || comp_all_cog)
      if (por_init[c_id] < 1.e-5)
        for (int i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = cell_cen[c_id][i];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cs_ibm->prob_dim == CS_IBM_2D_X) {
      cell_f_cen[c_id][0] = cell_cen[c_id][0];
      cell_s_cen[c_id][0] = cell_cen[c_id][0];
    } else if (cs_ibm->prob_dim == CS_IBM_2D_Y) {
      cell_f_cen[c_id][1] = cell_cen[c_id][1];
      cell_s_cen[c_id][1] = cell_cen[c_id][1];
    } else if (cs_ibm->prob_dim == CS_IBM_2D_Z) {
      cell_f_cen[c_id][2] = cell_cen[c_id][2];
      cell_s_cen[c_id][2] = cell_cen[c_id][2];
    }

    if (CS_F_(poro)->val[c_id] > 1.-1.e-5)
      for (int i = 0; i < 3; i++)
        cell_s_cen[c_id][i] = cell_cen[c_id][i];
  }

  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_f_cen, 3);
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_s_cen, 3);

  if (mesh->n_init_perio > 0) {
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_STANDARD,
                              (cs_real_t *)cell_f_cen);
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_STANDARD,
                              (cs_real_t *)cell_s_cen);
  }

  /* -> porosity = porbis + volume conservation */
  if (cs_ibm->porosity_from_nodes) {
    cs_real_t volpor = 0.;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0)
        volpor += cell_vol[c_id] * CS_F_(poro)->val[c_id];

    cs_parall_sum(1, CS_REAL_TYPE, &volpor);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0)
        CS_F_(poro)->val[c_id] = porbis[c_id];

    if (cs_ibm->ensure_isovol)
      for (int iter = 1; iter < 4; iter++) {
        cs_real_t aa = 0.;
        cs_real_t bb = 0.;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          if (comp_cell[c_id] > 0) {
            cs_real_t porloc =  CS_F_(poro)->val[c_id];
            cs_real_t pump = cs_math_fabs(porloc * (1. - porloc));
            aa += pump * cell_vol[c_id];
            bb += porloc * cell_vol[c_id];
          }

        cs_parall_sum(1, CS_REAL_TYPE, &aa);
        cs_parall_sum(1, CS_REAL_TYPE, &bb);

        cs_real_t beta = (volpor - bb) / cs_math_fmax(aa, 1.e-20);

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          if (comp_cell[c_id] > 0) {
            cs_real_t porloc = CS_F_(poro)->val[c_id];
            cs_real_t pump = cs_math_fabs(porloc * (1. - porloc));
            CS_F_(poro)->val[c_id] += beta * pump;
            CS_F_(poro)->val[c_id] = cs_math_fmax(CS_F_(poro)->val[c_id], 0.);
            CS_F_(poro)->val[c_id] = cs_math_fmin(CS_F_(poro)->val[c_id], 1.);
          }
      }

    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(poro)->val);
  }

  BFT_FREE(porbis);
  BFT_FREE(c_poro);
  BFT_FREE(weight);
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute boundary faces porosity
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   v_poro             vertex porosity
 * \param[out]  bfpro_poro         boundary face porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_b_fac_porosity(const cs_mesh_t            *mesh,
                        const cs_mesh_quantities_t *mesh_quantities,
                        cs_real_t                  *c_poro,
                        cs_real_t                  *v_poro,
                        cs_real_t                  *bfpro_poro)
{
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *b_face_surf  = mesh_quantities->b_face_surf;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mesh_quantities->b_face_normal;

  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  int icut = cs_ibm->nb_cut_cells;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    bfpro_poro[f_id] = c_poro[c_id];

    cs_real_3_t cog;
    for (int i = 0; i < 3; i++)
      cog[i] = 0.;

    cs_real_t porc = 0.;

    int cpt = 0;
    int ptin = 0;

    for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
         vtx_id < b_face_vtx_idx[f_id + 1]; vtx_id++) {
      cs_lnum_t v_id = b_face_vtx[vtx_id];

      cpt++;
      if (v_poro[v_id] < 0.5)
        ptin++;

      for (int i = 0; i < 3; i++)
        cog[i] += vtx_crd[v_id][i];

      porc += v_poro[v_id];
    }

    if (ptin == cpt)
      bfpro_poro[f_id] = 0.;

    else if (ptin == 0)
      bfpro_poro[f_id] = 1.;

    else {
      for (int i = 0; i < 3; i++)
        cog[i] /= (cs_real_t)(cpt);
      porc /= (cs_real_t)(cpt);

      cs_real_t surftot = b_face_surf[f_id];
      cs_real_t surfpor = 0.;

      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
          vtx_id < b_face_vtx_idx[f_id + 1] - 1; vtx_id++) {
        cs_lnum_t vtx_id1 = b_face_vtx[vtx_id];
        cs_lnum_t vtx_id2 = b_face_vtx[vtx_id + 1];

        cs_real_3_t x1, x2;
        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
        }
        cs_real_t por1 = v_poro[vtx_id1];
        cs_real_t por2 = v_poro[vtx_id2];

        surfpor += _tri_surf_trunc(x1, por1, x2, por2, cog, porc, icut);
      }

      cs_lnum_t id1_ = b_face_vtx_idx[f_id];
      cs_lnum_t id2_ = b_face_vtx_idx[f_id + 1] - 1;
      cs_lnum_t vtx_id1 = b_face_vtx[id1_];
      cs_lnum_t vtx_id2 = b_face_vtx[id2_];

      cs_real_3_t x1, x2;
      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
      }
      cs_real_t por1 = v_poro[vtx_id1];
      cs_real_t por2 = v_poro[vtx_id2];

      surfpor += _tri_surf_trunc(x1, por1, x2, por2, cog, porc, icut);

      cs_real_t porb = surfpor / surftot;
      bfpro_poro[f_id] = cs_math_fmin(porb, c_poro[c_id]);
    }

    /* Porosity treatment in the solid */
    if (cs_ibm->solid_porosity[c_id] > 1.e-5) {
      cs_real_t porb =  bfpro_poro[f_id];
      cs_real_t porbis = cs_math_fmin(cs_ibm->solid_porosity[c_id],
                                      c_poro[c_id]);
      bfpro_poro[f_id] = cs_math_fmax(porb, porbis);
    }

    /* 2D cases treatment */
    if (cs_ibm->prob_dim != CS_IBM_3D) {
      cs_real_3_t rn;
      for (int i = 0; i < 3; i++)
        rn[i] = b_face_normal[f_id][i]/b_face_surf[f_id];

      if (cs_ibm->prob_dim ==  CS_IBM_2D_X)
        if (rn[0]*rn[0] > rn[1]*rn[1] + rn[2]*rn[2])
          bfpro_poro[f_id] = c_poro[c_id];

      if (cs_ibm->prob_dim ==  CS_IBM_2D_Y)
        if (rn[1]*rn[1] > rn[0]*rn[0] + rn[2]*rn[2])
          bfpro_poro[f_id] = c_poro[c_id];

      if (cs_ibm->prob_dim ==  CS_IBM_2D_Z)
        if (rn[2]*rn[2] > rn[0]*rn[0] + rn[1]*rn[1])
          bfpro_poro[f_id] = c_poro[c_id];
    }
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute internal faces porosity
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   v_poro             vertex porosity
 * \param[out]  ifpro_poro         internal face porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_i_fac_porosity(const cs_mesh_t            *mesh,
                        const cs_mesh_quantities_t *mesh_quantities,
                        cs_real_t                  *c_poro,
                        cs_real_t                  *v_poro,
                        cs_real_t                  *ifpro_poro)
{
  CS_NO_WARN_IF_UNUSED(mesh_quantities);
  CS_NO_WARN_IF_UNUSED(v_poro);

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t pori = 0.5 * (c_poro[c_id0] + CS_F_(poro)->val_pre[c_id0]);
    cs_real_t porj = 0.5 * (c_poro[c_id1] + CS_F_(poro)->val_pre[c_id1]);

    cs_real_t porij = _geom_face_fraction(pori, porj);
    cs_real_t porimin = cs_math_fmin(c_poro[c_id0],
                                     CS_F_(poro)->val_pre[c_id0]);
    cs_real_t porjmin = cs_math_fmin(c_poro[c_id1],
                                     CS_F_(poro)->val_pre[c_id1]);

    if (porij < 1.e-5 || porimin < 1.e-5 || porjmin < 1.e-5)
      porij = 0.;

    cs_real_t porijmin = cs_math_fmin(pori, porj);
    porij = cs_math_fmin(porij, 50.*porijmin);

    ifpro_poro[f_id] = porij;
  }


  /* If for a cell, only one face porosity is positive,
   * one cancels face porosity */
  cs_lnum_t *cpt;
  BFT_MALLOC(cpt, n_cells_ext, cs_lnum_t);

  for (int ii = 0; ii < 4; ii++) {
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
      cpt[c_id] = 0;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      if (ifpro_poro[f_id] > 1.e-5) {
        cpt[c_id0]++;
        cpt[c_id1]++;
      }
    }

    cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, cpt);

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      if (cpt[c_id0] <= 1 || cpt[c_id1] <= 1)
        ifpro_poro[f_id] = 0;
    }
  }

  BFT_FREE(cpt);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Iso-volume porosity check
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_iso_vol_porosity(const cs_mesh_t            *mesh,
                          const cs_mesh_quantities_t *mesh_quantities,
                          cs_real_t                  *c_poro,
                          int                        *comp_cell)
{
  cs_lnum_t n_cells    = mesh->n_cells;
  cs_real_t *cell_vol  = mesh_quantities->cell_vol;

  static int ipass = 0;
  if (   cs_restart_present()
      && cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev)
    ipass = 1;
  ipass ++;

  if (ipass == 1) {
    cs_real_t volpor = 0.;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0)
        volpor += cell_vol[c_id] * c_poro[c_id];

    cs_parall_sum(1, CS_DOUBLE, &volpor);
    cs_ibm->isovol = volpor;

  } else if (ipass > 1) {
    cs_real_t aa = 0.;
    cs_real_t bb = 0.;
    cs_real_t volpor = cs_ibm->isovol;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0) {
        cs_real_t porloc = c_poro[c_id];
        cs_real_t pump = cs_math_fabs(porloc * (1. - porloc));
        aa += pump * cell_vol[c_id];
        bb += porloc * cell_vol[c_id];
      }

    cs_parall_sum(1, CS_DOUBLE, &aa);
    cs_parall_sum(1, CS_DOUBLE, &bb);

    cs_real_t beta = (volpor - bb) / cs_math_fmax(aa, 1.e-20);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0) {
        cs_real_t porloc = c_poro[c_id];
        cs_real_t pump = cs_math_fabs(porloc * (1. - porloc));
        c_poro[c_id] += beta * pump;
        c_poro[c_id] = cs_math_fmax(c_poro[c_id], 0.);
        c_poro[c_id] = cs_math_fmin(c_poro[c_id], 1.);
      }

    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_poro);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account internal solid porosity
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_solid_porosity(const cs_mesh_t            *mesh,
                        const cs_mesh_quantities_t *mesh_quantities,
                        cs_real_t                  *c_poro,
                        int                        *comp_cell)
{
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mesh_quantities->cell_cen;
  cs_real_3_t *cell_f_cen = (cs_real_3_t *)mesh_quantities->cell_f_cen;

  memset(cs_ibm->solid_porosity, 0, n_cells_ext*sizeof(cs_real_t));

  cs_real_t t_cur = cs_glob_time_step->t_cur;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    int num_object = -1;

    cs_user_ibm_solid_por(c_id,
                          cell_cen[c_id],
                          t_cur,
                          num_object);
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (comp_cell[c_id] > 0) {
      cs_real_t porsol = cs_ibm->solid_porosity[c_id];
      cs_real_t por1 =  c_poro[c_id];
      cs_real_t por2 = por1 + porsol * (1. - por1);
      c_poro[c_id] = por2;

      for (int i = 0; i < 3; i++) {
        cs_real_t xyz1 = cell_f_cen[c_id][i];
        cell_f_cen[c_id][i] = xyz1 + porsol * (cell_cen[c_id][i] - xyz1);
      }
    }
  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_poro);
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_f_cen, 3);

  if (mesh->n_init_perio > 0)
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_STANDARD,
                              (cs_real_t *)cell_f_cen);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute list of cells for porosity computing
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[out]  comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_list_porosity(const cs_mesh_t            *mesh,
                            const cs_mesh_quantities_t *mesh_quantities,
                            cs_lnum_t                  *comp_cell)
{
  CS_NO_WARN_IF_UNUSED(mesh_quantities);

  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    comp_cell[c_id] = 1;

  int nt_cur = cs_glob_time_step->nt_cur;
  int nt_prev = cs_glob_time_step->nt_prev;
  if (   (nt_cur > 1 && nt_cur > nt_prev)
      || cs_ibm->porosity_user_source_term_modification) {
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
      comp_cell[c_id] = 0;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
           vtx_id < i_face_vtx_idx[f_id + 1]; vtx_id++) {
        cs_lnum_t v_id = i_face_vtx[vtx_id];

        int iok = 0;
        for (int idim = 0; idim < 3; idim++) {
          if (   vtx_crd[v_id][idim] >= cs_ibm->xyzmin_moving_porosity[idim]
              && vtx_crd[v_id][idim] <= cs_ibm->xyzmax_moving_porosity[idim])
            iok++;
        }

        if (iok == 3) {
          comp_cell[c_id0] = 1;
          comp_cell[c_id1] = 1;
        }
      }
    }

    cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, comp_cell);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute porosity solid surface
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   i_poro             internal face porosity
 * \param[in]   b_poro             boundary face porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_solid_surface_vector(const cs_mesh_t            *mesh,
                              const cs_mesh_quantities_t *mesh_quantities,
                              cs_real_t                  *i_poro,
                              cs_real_t                  *b_poro)
{
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mesh_quantities->i_face_normal;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mesh_quantities->b_face_normal;
  cs_real_t *c_w_face_surf
    = (cs_real_t *)mesh_quantities->c_w_face_surf;
  cs_real_3_t *c_w_face_normal
    = (cs_real_3_t *)mesh_quantities->c_w_face_normal;

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int i = 0; i < 3; i++)
      c_w_face_normal[c_id][i] = 0;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t alpij = i_poro[f_id] - 1.;

    for (int idim = 0; idim < 3; idim++) {
      c_w_face_normal[c_id0][idim] += alpij  * i_face_normal[f_id][idim];
      c_w_face_normal[c_id1][idim] -= alpij  * i_face_normal[f_id][idim];
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t alpij = b_poro[f_id] - 1.;

    for (int idim = 0; idim < 3; idim++)
      c_w_face_normal[c_id][idim] += alpij  * b_face_normal[f_id][idim];
  }

  cs_real_3_t sx;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (int i = 0; i < 3; i++)
      sx[i] = c_w_face_normal[c_id][i];
    c_w_face_surf[c_id] = cs_math_3_norm(sx);
  }


  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)c_w_face_normal, 3);
  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_w_face_surf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute porosity solid surface cog (dist and ponderation coeff)
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   v_poro             vertex porosity
 * \param[in]   i_poro             internal face porosity
 * \param[in]   b_poro             boundary face porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_solid_surface_cog(const cs_mesh_t            *mesh,
                           const cs_mesh_quantities_t *mesh_quantities,
                           cs_real_t                  *c_poro,
                           cs_real_t                  *v_poro,
                           cs_real_t                  *i_poro,
                           cs_real_t                  *b_poro)
{
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)mesh->i_face_cells;

  cs_real_t *cell_vol = mesh_quantities->cell_vol;
  cs_real_t *dist = mesh_quantities->i_dist;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2f = ma->cell_i_faces;

  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)mesh_quantities->cell_cen;
  const cs_real_3_t *cell_f_cen
    = (const cs_real_3_t *)mesh_quantities->cell_f_cen;
  cs_real_3_t *cell_s_cen = (cs_real_3_t *)mesh_quantities->cell_s_cen;
  const cs_real_3_t *i_face_cog
    = (const cs_real_3_t *)mesh_quantities->i_face_cog;
  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mesh_quantities->i_face_normal;
  const cs_real_t *i_face_surf  = mesh_quantities->i_face_surf;
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)mesh_quantities->b_face_cog;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mesh_quantities->b_face_normal;
  const cs_real_t *b_face_surf  = mesh_quantities->b_face_surf;
  const cs_real_t *c_w_face_surf
    = (const cs_real_t *)mesh_quantities->c_w_face_surf;
  const cs_real_3_t *c_w_face_normal
    = (const cs_real_3_t *)mesh_quantities->c_w_face_normal;
  cs_real_3_t *c_w_face_cog
    = (cs_real_3_t *)mesh_quantities->c_w_face_cog;
  cs_real_t *c_w_dist_inv = (cs_real_t *)mesh_quantities->c_w_dist_inv;
  cs_real_t *i_f_weight = mesh_quantities->i_f_weight;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t num = 0.;
    cs_real_t den = 0.;
    cs_real_3_t xip, xjp, ipf, ipjp;
    for (int idim = 0; idim < 3; idim++) {
      xip[idim] = cell_f_cen[c_id0][idim];
      xjp[idim] = cell_f_cen[c_id1][idim];
      ipf[idim] = i_face_cog[f_id][idim] - xip[idim];
      ipjp[idim] = xjp[idim] - xip[idim];
      num += ipf[idim] * i_face_normal[f_id][idim];
      den += ipjp[idim] * i_face_normal[f_id][idim];
    }

    cs_real_t lambda = 0.5;
    if (cs_math_fabs(den) > 1.e-20)
      lambda = num / den;

    lambda = cs_math_fmin(cs_math_fmax(lambda, 0.001), 0.999);
    cs_real_t weight_loc = 1. - lambda;

    cs_real_t i_face_normal_unit[3];
    for (int idim = 0; idim < 3; idim++)
      i_face_normal_unit[idim] = i_face_normal[f_id][idim]/i_face_surf[f_id];

    cs_real_t dist_ipjp
      = cs_math_fabs(cs_math_3_distance_dot_product(xip, xjp,
                                                    i_face_normal_unit));
    dist_ipjp = cs_math_fmax(dist_ipjp, 0.5 * dist[f_id]);

    // Security
    if (cs_math_fmin(CS_F_(poro)->val[c_id0],CS_F_(poro)->val[c_id1]) > 0.9999)
      dist_ipjp = dist[f_id];

    weight_loc = cs_math_fmax(weight_loc, 0.2);
    i_f_weight[f_id] = cs_math_fmin(weight_loc, 0.8);
  }

  // TODO: Check

  cs_real_33_t *cut = NULL;
  BFT_MALLOC(cut, n_cells_ext, cs_real_33_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        cut[c_id][k][l] = 0.;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (int i = 0; i < 3; i++)
      cut[c_id][i][i] -= cell_vol[c_id] * c_poro[c_id];

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t alpij = i_poro[f_id];

    cs_real_t num = 0.;
    cs_real_t den = 0.;
    cs_real_3_t xip, xjp, ipf, ipjp, cogfac2;
    for (int idim = 0; idim < 3; idim++) {
      xip[idim] = cell_f_cen[c_id0][idim];
      xjp[idim] = cell_f_cen[c_id1][idim];
      ipf[idim] = i_face_cog[f_id][idim] - xip[idim];
      ipjp[idim] = xjp[idim] - xip[idim];
      num += ipf[idim] * i_face_normal[f_id][idim];
      den += ipjp[idim] * i_face_normal[f_id][idim];
    }

    cs_real_t lambda = 0.5;
    if (cs_math_fabs(den) > 1.e-20)
      lambda = num / den;
    lambda = cs_math_fmin(cs_math_fmax(lambda, 0.001), 0.999);
    cs_real_t weight_loc = 1. - lambda;

    if (i_poro[f_id] < 0.99)
      for (int idim = 0; idim < 3; idim++)
        cogfac2[idim] = weight_loc * xip[idim] + (1. - weight_loc) * xjp[idim];
    else
      for (int idim = 0; idim < 3; idim++)
        cogfac2[idim] = i_face_cog[f_id][idim];

    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++) {
        cut[c_id0][k][l] += alpij * cogfac2[k] * i_face_normal[f_id][l];
        cut[c_id1][k][l] -= alpij * cogfac2[k] * i_face_normal[f_id][l];
      }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t alpij = b_poro[f_id];

    cs_real_t num = 0.;
    cs_real_t den = 0.;
    cs_real_3_t xip, ipf, b_face_cog_tp;
    for (int idim = 0; idim < 3; idim++) {
      xip[idim] = cell_f_cen[c_id][idim];
      ipf[idim] = b_face_cog[f_id][idim] - xip[idim];
      num += ipf[idim] * b_face_normal[f_id][idim];
      den += cs_math_pow2(b_face_normal[f_id][idim]);
    }

    cs_real_t lambda = num / den;
    for (int idim = 0; idim < 3; idim++)
      b_face_cog_tp[idim] = xip[idim] + lambda * b_face_normal[f_id][idim];

    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        cut[c_id][k][l] += alpij * b_face_cog_tp[k] * b_face_normal[f_id][l];
  }

  cs_real_t *denom2;
  BFT_MALLOC(denom2, n_cells_ext, cs_real_t);

  cs_real_3_t *coord_node_out;
  BFT_MALLOC(coord_node_out, n_cells_ext, cs_real_3_t);
  cs_real_3_t *coord_node_min;
  BFT_MALLOC(coord_node_min, n_cells_ext, cs_real_3_t);
  cs_real_3_t *coord_node_max;
  BFT_MALLOC(coord_node_max, n_cells_ext, cs_real_3_t);

  /* Compute positions out, min and max */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int idim = 0; idim < 3; idim++) {
      coord_node_out[c_id][idim] = 0.;
      coord_node_min[c_id][idim] = 1.e20;
      coord_node_max[c_id][idim] =-1.e20;
    }

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    denom2[c_id] = 0.;

  cs_real_t pornmin = 0.5;
  cs_real_t pornmax = 0.52;

  /* Loop on the internal faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
         vtx_id < i_face_vtx_idx[f_id + 1]; vtx_id++) {
      cs_lnum_t v_id = i_face_vtx[vtx_id];

      if (v_poro[v_id] < pornmax) {
        cs_real_t pds = (pornmax - v_poro[v_id]) / (pornmax - pornmin);
        pds = cs_math_fmin(pds, 1.);

        for (int idim = 0; idim < 3; idim++)
          coord_node_out[c_id0][idim] += pds * (1. - v_poro[v_id])
                                             * vtx_crd[v_id][idim];

        for (int idim = 0; idim < 3; idim++)
          coord_node_out[c_id1][idim] += pds * (1. - v_poro[v_id])
                                             * vtx_crd[v_id][idim];

        denom2[c_id0] += pds * (1. - v_poro[v_id]);
        denom2[c_id1] += pds * (1. - v_poro[v_id]);
      }

      for (int idim = 0; idim < 3; idim++)
        coord_node_min[c_id0][idim] = cs_math_fmin(coord_node_min[c_id0][idim],
                                                   vtx_crd[v_id][idim]);
      for (int idim = 0; idim < 3; idim++)
        coord_node_min[c_id1][idim] = cs_math_fmin(coord_node_min[c_id1][idim],
                                                   vtx_crd[v_id][idim]);

      for (int idim = 0; idim < 3; idim++)
        coord_node_max[c_id0][idim] = cs_math_fmax(coord_node_max[c_id0][idim],
                                                   vtx_crd[v_id][idim]);
      for (int idim = 0; idim < 3; idim++)
        coord_node_max[c_id1][idim] = cs_math_fmax(coord_node_max[c_id1][idim],
                                                   vtx_crd[v_id][idim]);
    }
  } /* End of loop on internal faces */

  /* Loop on the boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
         vtx_id < b_face_vtx_idx[f_id + 1]; vtx_id++) {
      cs_lnum_t v_id = b_face_vtx[vtx_id];

      if (v_poro[v_id] < pornmax) {
        cs_real_t pds = (pornmax - v_poro[v_id]) / (pornmax - pornmin);
        pds = cs_math_fmin(pds, 1.);

        for (int idim = 0; idim < 3; idim++)
          coord_node_out[c_id][idim] += pds * (1. - v_poro[v_id])
                                            * vtx_crd[v_id][idim];

        denom2[c_id] += pds * (1. - v_poro[v_id]);
      }

      for (int idim = 0; idim < 3; idim++)
        coord_node_min[c_id][idim] = cs_math_fmin(coord_node_min[c_id][idim],
                                                  vtx_crd[v_id][idim]);

      for (int idim = 0; idim < 3; idim++)
        coord_node_max[c_id][idim] = cs_math_fmax(coord_node_max[c_id][idim],
                                                  vtx_crd[v_id][idim]);
    }
  } /* End of loop on boundary faces */

  /* Finalization */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (int idim = 0; idim < 3; idim++)
      coord_node_out[c_id][idim] /= cs_math_fmax(denom2[c_id], 1.e-20);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (int i = 0; i < 3; i++)
      c_w_face_cog[c_id][i] = cell_cen[c_id][i];

    if (c_poro[c_id] < 1.e-5)
      for (int i = 0; i < 3; i++)
        c_w_face_cog[c_id][i] = cell_f_cen[c_id][i];

    cs_real_3_t nx;
    for (int i = 0; i < 3; i++)
      nx[i] = c_w_face_normal[c_id][i];

    cs_real_t n2 = cs_math_3_square_norm(nx);

    if (n2 > 1.e-20) {
      cs_real_3_t xp;
      for (int i = 0; i < 3; i++)
        xp[i] = cs_math_3_dot_product(cut[c_id][i], nx) / n2;

      int iok = 1;

      const cs_lnum_t s_id = c2c_idx[c_id];
      const cs_lnum_t e_id = c2c_idx[c_id+1];

      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t f_id = c2f[i];
        cs_real_t psca = cs_math_3_distance_dot_product(xp, i_face_cog[f_id],
                                                        i_face_normal[f_id]);
        psca *= i_poro[f_id];

        if (psca < -1.e-10)
          iok = 0;
      }

      int iok2 = 1;
      if (c_poro[c_id] > 0.75) {
        cs_real_t psca2
          = (xp[0]-cell_cen[c_id][0])*(cell_f_cen[c_id][0]-cell_cen[c_id][0])
          + (xp[1]-cell_cen[c_id][1])*(cell_f_cen[c_id][1]-cell_cen[c_id][1])
          + (xp[2]-cell_cen[c_id][2])*(cell_f_cen[c_id][2]-cell_cen[c_id][2]);
        cs_real_t psca3
          = cs_math_3_distance_dot_product(cell_f_cen[c_id], xp, nx);

        if (psca2 > 0. || psca3 > 0.)
          iok2 = 0;
      }

      if ((iok == 0 || iok2 == 0) && denom2[c_id] >= 1.e-20)
        for (int i = 0; i < 3; i++)
          xp[i] = cell_cen[c_id][i]
                + 2. * (c_poro[c_id] - 0.5)
                     * (coord_node_out[c_id][i] - cell_cen[c_id][i]);

      for (int i = 0; i < 3; i++) {
        xp[i] = cs_math_fmax(xp[i], coord_node_min[c_id][i]);
        xp[i] = cs_math_fmin(xp[i], coord_node_max[c_id][i]);
      }

      for (int i = 0; i < 3; i++)
        c_w_face_cog[c_id][i] = xp[i];
    }

    cs_real_3_t xf_xc, xs_xc;
    for (int i = 0; i < 3; i++) {
      xf_xc[i] = cell_f_cen[c_id][i] - c_w_face_cog[c_id][i];
      xs_xc[i] = cell_s_cen[c_id][i] - c_w_face_cog[c_id][i];
    }

    cs_real_t psca3 = cs_math_3_dot_product(xf_xc, xs_xc);

    if (psca3 > 1.e-20)
      for (int i = 0; i < 3; i++)
        cell_s_cen[c_id][i] = c_w_face_cog[c_id][i];
  }

  /* Readjustment of cog_cut */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)   {
    cs_real_t poro = CS_F_(poro)->val[c_id];

    if (poro > 1.e-6 && poro < 1. - 1.e-6) {
      cs_real_3_t nxyz;
      for (int i = 0; i < 3; i++)
        nxyz[i] = c_w_face_normal[c_id][i];

      cs_real_t nn = cs_math_fmax(cs_math_3_norm(nxyz), 1.e-20);

      for (int i = 0; i < 3; i++)
        nxyz[i] /= nn;

      cs_real_3_t xyzc, xyzp;
      for (int i = 0; i < 3; i++) {
        xyzc[i] = cell_f_cen[c_id][i];
        xyzp[i] = c_w_face_cog[c_id][i];
      }

      cs_real_t psca = cs_math_3_distance_dot_product(xyzp, xyzc, nxyz);
      if (psca < 0.) {
        for (int i = 0; i < 3; i++)
          xyzp[i] = cell_cen[c_id][i] + 0.9 * (xyzc[i] - cell_cen[c_id][i]);

        for (int i = 0; i < 3; i++)
          c_w_face_cog[c_id][i] = xyzp[i];
      }
    }
  }

  BFT_FREE(cut);
  BFT_FREE(coord_node_out);
  BFT_FREE(coord_node_min);
  BFT_FREE(coord_node_max);
  BFT_FREE(denom2);

  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)c_w_face_cog, 3);
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_s_cen, 3);

  // Compute the "characteristic size" vector of each cell
  cs_real_3_t *cell_length3D, *weight3D;
  BFT_MALLOC(cell_length3D, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(weight3D, n_cells_ext, cs_real_3_t);
  memset(cell_length3D, 0, n_cells_ext*sizeof(cs_real_3_t));
  memset(weight3D, 0, n_cells_ext*sizeof(cs_real_3_t));

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t i_face_normal_unit_cp = 0.;
    for (int idim = 0; idim < 3; idim++) {
      i_face_normal_unit_cp = i_face_normal[f_id][idim]/i_face_surf[f_id];
      cell_length3D[c_id0][idim]
        += cs_math_fabs((cell_cen[c_id1][idim] - cell_cen[c_id0][idim])
                        * i_face_normal_unit_cp);
      cell_length3D[c_id1][idim]
        += cs_math_fabs((cell_cen[c_id1][idim] - cell_cen[c_id0][idim])
                        * i_face_normal_unit_cp);

      weight3D[c_id0][idim] += cs_math_fabs(i_face_normal_unit_cp);
      weight3D[c_id1][idim] += cs_math_fabs(i_face_normal_unit_cp);
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    cs_real_t b_face_normal_unit_cp = 0.;
    for (int idim = 0; idim < 3; idim++) {
      b_face_normal_unit_cp = b_face_normal[f_id][idim]/b_face_surf[f_id];
      cell_length3D[c_id][idim]
        += cs_math_fabs(2. * (b_face_cog[f_id][idim] - cell_cen[c_id][idim])
                        * b_face_normal_unit_cp);
      weight3D[c_id][idim] += cs_math_fabs(b_face_normal_unit_cp);
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (int i = 0; i < 3; i++)
      cell_length3D[c_id][i] /= weight3D[c_id][i];

  /* c_w_dist_inv : inverse distance from new cell_cog to immersed wall cog */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_3_t npx, xip, xpp;
    for (int idim = 0; idim < 3; idim++) {
      xip[idim] = cell_f_cen[c_id][idim];
      xpp[idim] = c_w_face_cog[c_id][idim];
      npx[idim] = c_w_face_normal[c_id][idim];
    }

    cs_real_t nn = c_w_face_surf[c_id];
    nn = cs_math_fmax(nn, 1.e-20);

    for (int idim = 0; idim < 3; idim++)
      npx[idim] /= nn;

    cs_real_t dwall
      = cs_math_fabs(cs_math_3_distance_dot_product(xpp, xip, npx));

    cs_real_3_t dipxyz;
    for (int i = 0; i < 3; i++)
      dipxyz[i] = xpp[i] - xip[i];

    cs_real_t dipn = cs_math_3_norm(dipxyz);

    if (dipn < 1.e-10) {
      for (int i = 0; i < 3; i++)
        dipxyz[i] = 1.e-10;

      dipn = cs_math_3_norm(dipxyz);
    }

    cs_real_t dcell_geom =  cs_math_fabs(cell_length3D[c_id][0] * dipxyz[0])
                          + cs_math_fabs(cell_length3D[c_id][1] * dipxyz[1])
                          + cs_math_fabs(cell_length3D[c_id][2] * dipxyz[2]);
    dcell_geom /= dipn;

    c_w_dist_inv[c_id]
      = 1./cs_math_fmax(dwall, 1.e-3 * pow(cell_vol[c_id], cs_math_1ov3));
    c_w_dist_inv[c_id]
      = cs_math_fmin(c_w_dist_inv[c_id], 1./(0.1 * dcell_geom));

  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_w_dist_inv);

  BFT_FREE(cell_length3D);
  BFT_FREE(weight3D);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new immersed boundary object with a given name and method.
 *
 * \param[in] name    name of the new object
 * \param[in] method  method used to compute the object porosity
 *
 * \return pointer to new object structure
 */
/*----------------------------------------------------------------------------*/

static cs_ibm_object_t *
_create_ibm_object(const char          *name,
                   cs_ibm_algo_type_t   method)
{

  cs_ibm_object_t *new_obj = NULL;

  BFT_MALLOC(new_obj, 1, cs_ibm_object_t);

  /* Set name */
  if (name == NULL || strcmp(name, "") == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty name provided for IBM object creation.\n"));

  new_obj->name = NULL;
  BFT_MALLOC(new_obj->name, strlen(name) + 1, char);
  strcpy(new_obj->name, name);

  /* Method */
  new_obj->method = method;

  /* Pointer to medcoupling or stl mesh structures */
  new_obj->cutcell_func = NULL;
  new_obj->stl = NULL;
  new_obj->mi  = NULL;

  for (int i = 0; i < CS_N_IBM_OBJ_PROP_TYPES; i++)
    new_obj->property_defs[i] = NULL;

  for (int i = 0; i < CS_N_IBM_OBJ_INIT_TYPES; i++)
    new_obj->init_vals_defs[i] = NULL;

  return new_obj;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Adds a new immersed boundary object with a given name and method.
 *
 * \param[in] name    name of the new object
 * \param[in] method  method used to compute the object porosity
 *
 * \return id of the new object (int)
 */
/*----------------------------------------------------------------------------*/

static int
_add_ibm_object(const char          *name,
                cs_ibm_algo_type_t   method)
{

  /* Check that the chosen algorithm is the correct one */
  if (cs_ibm->algo_choice == CS_IBM_ALGO_NONE)
    cs_ibm->algo_choice = method;
  if (cs_ibm->algo_choice != method)
    bft_error(__FILE__, __LINE__, 0,
              _("Current approach requires all objects be defined using the "
                "same method.\n You tried to define an object using the \"%s\" "
                "algorithm while the global algorithm is \"%s\".\n"),
              _ibm_algo_names[method],
              _ibm_algo_names[cs_ibm->algo_choice]);

  /* If object allready exists, exit the function */
  cs_ibm_object_t *obj = cs_ibm_object_by_name_try(name);

  if (obj != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating object: object \"%s\" already exists.\n"),
              name);

  int new_obj_id = cs_ibm->n_objects;

  if (new_obj_id == 0)
    BFT_MALLOC(cs_ibm->objects, new_obj_id + 1, cs_ibm_object_t *);
  else
    BFT_REALLOC(cs_ibm->objects, new_obj_id + 1, cs_ibm_object_t *);

  cs_ibm->objects[new_obj_id] = _create_ibm_object(name, method);

  cs_ibm->n_objects += 1;

  return new_obj_id;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_ibm_object_t structure.
 *
 * \param[in] obj  pointer to structure to free
 */
/*----------------------------------------------------------------------------*/

static void
_free_ibm_object(cs_ibm_object_t *obj)
{

  BFT_FREE(obj->name);

  if (obj->cutcell_func != NULL)
    obj->cutcell_func = NULL;

  BFT_FREE(obj);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new property definition for an object.
 *
 * \param[in] obj       pointer to object
 * \param[in] ppty_id   property id (si enum for list)
 * \param[in] n_vals    number of values (array size)
 * \param[in] vals      array of values
 */
/*----------------------------------------------------------------------------*/

static void
_ibm_object_define_property_def(cs_ibm_object_t               *obj,
                                cs_ibm_object_property_type_t  ppty_id,
                                int                            n_vals,
                                cs_real_t                     *vals)
{

  assert(ppty_id < CS_N_IBM_OBJ_PROP_TYPES);
  cs_xdef_t *def = obj->property_defs[ppty_id];

  if (def != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Property \"%s\" was already set for object \"%s\".\n"),
              _ibm_obj_property_names[ppty_id],
              obj->name);

  def = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                              n_vals,                /* Array size */
                              -1,                    /* Zone id */
                                                    /* -1 since no zone */
                              CS_FLAG_STATE_UNIFORM, /* Uniform value */
                              0,                     /* No meta data */
                              vals);                 /* Pointer to values */

  obj->property_defs[ppty_id] = def;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set initial value xdef for an object
 *
 * \param[in]  obj      Pointer to object
 * \param[in]  p_id     Property id (enum)
 * \param[in]  n_vals   Number of values (array size)
 * \param[in]  vals     Values array (size n_vals)
 */
/*----------------------------------------------------------------------------*/

static void
_ibm_object_define_initial_val_def(cs_ibm_object_t             *obj,
                                   cs_ibm_object_init_param_t   p_id,
                                   int                          n_vals,
                                   cs_real_t                   *vals)
{
  assert(p_id >= 0 && p_id < CS_N_IBM_OBJ_INIT_TYPES);

  cs_xdef_t *def = obj->init_vals_defs[p_id];

  if (def != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Initial value of \"%s\" was already set for object \"%s\".\n"),
              _ibm_obj_init_vals_names[p_id],
              obj->name);

  def = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                              n_vals,                /* Array size */
                              -1,                    /* Zone id */
                                                     /* -1 since no zone */
                              CS_FLAG_STATE_UNIFORM, /* Uniform value */
                              0,                     /* No meta data */
                              vals);                 /* Pointer to values */

  obj->init_vals_defs[p_id] = def;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get xdef value at object element
 *
 * \param[in] def     Pointer to cs_xdef_t definition
 * \param[in] elt_id  Element id (int)
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_get_xdef_val_at_object_elt(const cs_xdef_t *def,
                            const int        elt_id)
{
  assert(elt_id < def->dim);

  cs_real_t retval = -1.e30;

  const cs_real_t *cval = (cs_real_t *)def->context;

  if (def->dim == 1) {
    retval = cval[0];
  }
  else {
    retval = cval[elt_id];
  }

  return retval;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its id.
 *
 * \param[in] obj_id  id of the object
 *
 * \return pointer to object structure.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_id(int obj_id)
{
  assert(obj_id > -1 && obj_id < cs_ibm->n_objects);

  cs_ibm_object_t *obj = cs_ibm->objects[obj_id];

  return obj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Try to get an object based on its name. Returns NULL if not found
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure, NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name_try(const char *name)
{
  assert(name != NULL);

  cs_ibm_object_t *obj = NULL;

  for (int i = 0; i < cs_ibm->n_objects; i++) {
    if (strcmp(name, cs_ibm->objects[i]->name) == 0) {
      obj = cs_ibm->objects[i];
      break;
    }
  }

  return obj;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its name. Error if not found
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure, NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name(const char *name)
{
  assert(name != NULL);

  cs_ibm_object_t *obj = NULL;

  for (int i = 0; i < cs_ibm->n_objects; i++) {
    if (strcmp(name, cs_ibm->objects[i]->name) == 0) {
      obj = cs_ibm->objects[i];
      break;
    }
  }

  if (obj == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Object \"%s\" does not exist.\n"),
              name);
  return obj;
}

/*----------------------------------------------------------------------------
 * Create an empty cs_ibm structure
 *
 * returns:
 *   pointer to created cs_ibm structure
 *----------------------------------------------------------------------------*/

cs_ibm_t *
cs_ibm_create(void)
{
  cs_ibm_t  *ibm = NULL;

  BFT_MALLOC(ibm, 1, cs_ibm_t);

  ibm->n_objects = 0;
  ibm->objects = NULL;

  ibm->prob_dim       = CS_IBM_3D;
  ibm->algo_choice    = CS_IBM_ALGO_CUT_CELLS;
  ibm->wall_condition = CS_IBM_WALL_LAW_WALL_CONDITION;
  ibm->nb_cut_cells   = 1;
  ibm->nb_cut_faces   = 1;
  ibm->porosity_user_source_term_modification = false;
  ibm->solid_porosity = NULL;
  ibm->isovol         = 0.;
  ibm->ensure_isovol  = false;
  ibm->porosity_from_nodes = false;

  for (int i = 0; i < 3; i++) {
    ibm->xyzmin_moving_porosity[i] = -1.e20;
    ibm->xyzmax_moving_porosity[i] = 1.e20;
  }

  return (ibm);
}

/*----------------------------------------------------------------------------
 * Destroy a cs_ibm structure
 *
 * ibm <-- pointer to a cs_ibm_t structure
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

void
cs_ibm_finalize(void)
{
  if (cs_ibm != NULL) {

  BFT_FREE(cs_ibm->solid_porosity);

  /* Clean objects */
  for (int iobj = 0; iobj < cs_ibm->n_objects; iobj++) {
    _free_ibm_object(cs_ibm->objects[iobj]);
  }
  BFT_FREE(cs_ibm->objects);

  BFT_FREE(cs_ibm);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Define immersed boundaries in time and space
 *          (solid(s) interior part).
 *
 *  This function is called several times during each time step.
 *
 *  Ipenal: 1 means only solid and 0 only fluid.
 *
 *  Warning, porosity values have to be 0 or 1.
 *
 * The solid(s) velocities and porosity are prescirbed within the user function
 * (cs_user_ibm).
 *
 *
 * \param[in]  mesh                 pointer to associated mesh structure
 * \param[in]  mesh_quantities      pointer to associated mesh quantities
 */
/*----------------------------------------------------------------------------*/

void cs_immersed_boundaries(const cs_mesh_t *mesh,
                            const cs_mesh_quantities_t *mesh_quantities)
{
  cs_lnum_t n_cells     = mesh->n_cells;
  cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  cs_lnum_t n_i_faces   = mesh->n_i_faces;
  cs_lnum_t n_b_faces   = mesh->n_b_faces;

  cs_real_3_t *restrict i_f_face_normal =
     (cs_real_3_t *restrict)mesh_quantities->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal =
     (cs_real_3_t *restrict)mesh_quantities->b_f_face_normal;
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)mesh->b_face_cells;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)mesh_quantities->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mesh_quantities->b_face_normal;
  cs_real_t *restrict i_f_face_surf
    = (cs_real_t *restrict)mesh_quantities->i_f_face_surf;
  cs_real_t *restrict b_f_face_surf
    = (cs_real_t *restrict)mesh_quantities->b_f_face_surf;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_real_3_t *cell_f_cen
    = (const cs_real_3_t *)mesh_quantities->cell_f_cen;
  cs_real_t *cell_vol  = mesh_quantities->cell_vol;
  cs_real_t *cell_f_vol  = mesh_quantities->cell_f_vol;

  static int ipass = 0;
  ipass++;

  if (ipass == 1) {
    /* Initialization for time/space immersed boundaries */
    if (cs_glob_porosity_ibm_opt->porosity_mode > 0)
      cs_ibm = cs_ibm_create();
  }

  /* Value pointers for fields */
  cs_real_t *ifpro_poro = cs_field_by_name("i_face_porosity")->val;
  cs_real_t *bfpro_poro = cs_field_by_name("b_face_porosity")->val;

  /* Structure members allocation */
  if (cs_ibm->solid_porosity == NULL)
    BFT_MALLOC(cs_ibm->solid_porosity, n_cells_ext, cs_real_t);

  /* First call to user function to determine the problem dimension
   * and chosen algo for porosity */

  /* Initialize with undefined values to force user choice */
  if (ipass == 1) {
    cs_ibm_user_parameters();
    cs_ibm_init_writer();
  }

  int nt_cur = cs_glob_time_step->nt_cur;
  int nt_prev = cs_glob_time_step->nt_prev;
  if (   cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT
      || nt_cur == 0 || (cs_restart_present() && nt_cur == nt_prev)
      || (_porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID
          && cs_ibm->porosity_user_source_term_modification) )  {

    if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_TRANSIENT) {
      BFT_REALLOC(cs_ibm->solid_porosity, n_cells_ext, cs_real_t);
    }

    cs_real_3_t *gradp;

    BFT_MALLOC(gradp, n_cells_ext, cs_real_3_t);

    int hyd_p_flag = cs_glob_velocity_pressure_param->iphydr;
    cs_real_3_t *f_ext = (hyd_p_flag == 1) ?
      (cs_real_3_t *)cs_field_by_name_try("volume_forces")->val:NULL;

    bool use_previous_t = false;
    int inc = 1;
    cs_field_gradient_potential(CS_F_(p),
                                use_previous_t,
                                inc,
                                hyd_p_flag,
                                f_ext,
                                gradp);

    /* Initializer porosity at the previous value */
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
      CS_F_(poro)->val_pre[c_id] = CS_F_(poro)->val[c_id];

    /* List of cells for which one has to recompute porosity -> comp_cell */
    int *comp_cell;

    BFT_MALLOC(comp_cell, n_cells_ext, int);
    _compute_cell_list_porosity(mesh, mesh_quantities, comp_cell);

    /* Compute cell porosity */
    if (   nt_cur <=1
       || !cs_ibm->porosity_user_source_term_modification
       || _porosity_ibm_opt.porosity_mode != CS_IBM_FIXED_SOLID) {

      /* Compute cell porosity with a cut-cell method
       * from the cs_user_ibm function */
      if (cs_ibm->algo_choice == CS_IBM_ALGO_CUT_CELLS) {

        _compute_cell_cut_porosity(mesh, mesh_quantities, comp_cell);

      /* Compute cell porosity from a file */
      /* Object logic in this section */
      } else {

        /* User imposed rotations/translations */
        cs_user_ibm_object_transformations(cs_glob_time_step->t_cur);

        /* Some initialization */
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          CS_F_(poro)->val[c_id] = 1.;

        /* Local declarations */
        cs_real_t *obj_vol_f_tot = NULL;
        BFT_MALLOC(obj_vol_f_tot, n_cells_ext, cs_real_t);
        memset(obj_vol_f_tot, 0, n_cells_ext*sizeof(cs_real_t));

        for (int o_id = 0; o_id < cs_ibm->n_objects; o_id++) {
          cs_ibm_object_t *ibm_obj = cs_ibm_object_by_id(o_id);

          /* Increment total solid volume fraction using the current object.
           * Update object indicator array if needed. */
          cs_ibm_object_compute_intersect_vol(ibm_obj,
                                              mesh,
                                              cell_vol,
                                              obj_vol_f_tot,
                                              NULL);

        }

        /* Once total volume of solids in each cell is computed,
         * substract the result from the porosity. */
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          CS_F_(poro)->val[c_id] -= obj_vol_f_tot[c_id];

        /* Clip sync porosity values */
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          CS_F_(poro)->val[c_id] =
            cs_math_fmin(cs_math_fmax(CS_F_(poro)->val[c_id], 0.), 1.);

        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(poro)->val);

        BFT_FREE(obj_vol_f_tot);
      }
    }

    /* Possible modification of porosity by the user */
    if (cs_ibm->porosity_user_source_term_modification &&
        (nt_cur >= 2 || _porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID))
      cs_user_ibm_modify(mesh, mesh_quantities);

    /* One guarantees the same volume at each iteration
     * for cells whose porosity did not change */
    if (nt_cur > 1 && nt_cur > nt_prev)
      if (cs_ibm->ensure_isovol)
        if (!cs_ibm->porosity_user_source_term_modification)
          _compute_iso_vol_porosity(mesh, mesh_quantities,
                                    CS_F_(poro)->val, comp_cell);

    /* Porosity projection at vertices */
    cs_real_t *v_poro;
    BFT_MALLOC(v_poro, mesh->n_vertices, cs_real_t);

    int ncycle = 1;

    /* Save cell cog before modification */
    cs_real_3_t *cog_save;
    BFT_MALLOC(cog_save, n_cells, cs_real_3_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      for (int idir = 0; idir < 3; idir++)
        cog_save[c_id][idir] = cell_f_cen[c_id][idir];

    for (int ii = 0; ii < ncycle; ii++) {

      /* Update val_pre at the first time step */
      if (   cs_glob_time_step->nt_cur <= 1
          && _porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID)
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          CS_F_(poro)->val_pre[c_id] = CS_F_(poro)->val[c_id];

      cs_cell_to_vertex(CS_CELL_TO_VERTEX_SHEPARD, 0, 1, false, NULL,
                        CS_F_(poro)->val, NULL, v_poro);

      /* Compute cell centers of gravity */
      _compute_cell_cog(mesh, mesh_quantities,
                        v_poro, CS_F_(poro)->val, comp_cell);

      /* One guarantees the same volume at each iteration
       * for cells whose porosity did not change */
      if (nt_cur > 1 && nt_cur > nt_prev)
        if (cs_ibm->ensure_isovol)
          if (!cs_ibm->porosity_user_source_term_modification)
            _compute_iso_vol_porosity(mesh, mesh_quantities,
                                      CS_F_(poro)->val, comp_cell);

      /* Take into account inner porosity of solid */
      _compute_solid_porosity(mesh, mesh_quantities,
                              CS_F_(poro)->val, comp_cell);

      cs_cell_to_vertex(CS_CELL_TO_VERTEX_SHEPARD, 0, 1, false, NULL,
                        CS_F_(poro)->val, NULL, v_poro);

      /* Boundary face porosity */
      _compute_b_fac_porosity(mesh, mesh_quantities,
                              CS_F_(poro)->val, v_poro, bfpro_poro);

      /* Internal face porosity */
      _compute_i_fac_porosity(mesh, mesh_quantities,
                              CS_F_(poro)->val, v_poro, ifpro_poro);

      /* Compute solid surface vector */
      _compute_solid_surface_vector(mesh, mesh_quantities,
                                    ifpro_poro, bfpro_poro);

      /* Compute cog solid face + weight + dist + dist_wall */
      _compute_solid_surface_cog(mesh, mesh_quantities,
                                 CS_F_(poro)->val, v_poro,
                                ifpro_poro, bfpro_poro);
    }

    BFT_FREE(comp_cell);
    BFT_FREE(v_poro);

    /* Pressure update after moving cell cog */
    cs_real_3_t dcog;
    if (ipass >= 3) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (int idir = 0; idir < 3; idir++)
          dcog[idir] = cell_f_cen[c_id][idir] - cog_save[c_id][idir];

        CS_F_(p)->val[c_id] += cs_math_3_dot_product(dcog, gradp[c_id]);
        CS_F_(p)->val_pre[c_id] = CS_F_(p)->val[c_id];
      }

      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(p)->val);
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(p)->val_pre);
    }

    BFT_FREE(cog_save);
    BFT_FREE(gradp);
  }

  /* Porosity */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    cs_real_t porosity = CS_F_(poro)->val[c_id];
    if (porosity < cs_math_epzero) {
      CS_F_(poro)->val[c_id] = 0.;
      mesh_quantities->c_disable_flag[c_id] = 1;
    }
    cell_f_vol[c_id] = CS_F_(poro)->val[c_id] * cell_vol[c_id];
  }

  /* Set interior face values */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t c_id0 = i_face_cells[face_id][0];
    cs_lnum_t c_id1 = i_face_cells[face_id][1];

    cs_real_t face_porosity = CS_MIN(CS_F_(poro)->val[c_id0],
                                     CS_F_(poro)->val[c_id1]);

    for (cs_lnum_t i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

    mesh_quantities->i_f_face_surf[face_id]
      = cs_math_3_norm(i_f_face_normal[face_id]);

    if (mesh_quantities->i_f_face_factor != NULL) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mesh_quantities->i_f_face_factor[face_id][0]
      //    = cpro_porosi[c_id0] / face_porosity;
      //  mesh_quantities->i_f_face_factor[face_id][1]
      //    = cpro_porosi[c_id1] / face_porosity;
      //}
      //else {
        mesh_quantities->i_f_face_factor[face_id][0] = 1.;
        mesh_quantities->i_f_face_factor[face_id][1] = 1.;
      //}
    }

  }

  /* Set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    cs_lnum_t c_id = b_face_cells[face_id];

    cs_real_t face_porosity = CS_F_(poro)->val[c_id];

    for (cs_lnum_t i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

    mesh_quantities->b_f_face_surf[face_id]
      = cs_math_3_norm(b_f_face_normal[face_id]);

    if (mesh_quantities->b_f_face_factor != NULL) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mesh_quantities->b_f_face_factor[face_id]
      //    = cpro_porosi[c_id] / face_porosity;
      //}
      //else {
        mesh_quantities->b_f_face_factor[face_id] = 1.;
      //}
    }
  }

  /* Set interior face values */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t c_id0 = i_face_cells[face_id][0];
    cs_lnum_t c_id1 = i_face_cells[face_id][1];

    if (   mesh_quantities->c_disable_flag[c_id0] == 1
        || mesh_quantities->c_disable_flag[c_id1] == 1) {
      i_f_face_normal[face_id][0] = 0.;
      i_f_face_normal[face_id][1] = 0.;
      i_f_face_normal[face_id][2] = 0.;
      i_f_face_surf[face_id] = 0.;
    }
  }

  /* Set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    cs_lnum_t c_id0 = b_face_cells[face_id];

    if (mesh_quantities->c_disable_flag[c_id0] == 1) {
      b_f_face_normal[face_id][0] = 0.;
      b_f_face_normal[face_id][1] = 0.;
      b_f_face_normal[face_id][2] = 0.;
      b_f_face_surf[face_id] = 0.;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a point is solid or fluid based on the cut-cell method.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

int
cs_ibm_object_compute_cut_porosity(const cs_lnum_t    c_id,
                                   const cs_real_3_t  xyz,
                                   const cs_real_t    t,
                                   const int          num_object)
{

  int retval = 0;

  for (int i = 0; i < cs_ibm->n_objects; i++) {
    cs_ibm_object_t *obj = cs_ibm->objects[i];
    if (obj->method == CS_IBM_ALGO_CUT_CELLS) {
      retval = obj->cutcell_func(c_id, xyz, t, num_object);

      /* Stop at first intersection */
      if (retval > 0)
        break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an object from a file using STL or MED formats
 *
 * \param[in] name          name of the object
 * \param[in] method        Porosity computation method
 * \param[in] file_name     file name
 * \param[in] solve_fsi     Is the object used in the FSI resolution ?
 *                          (currently ignored)
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_object_from_file(const char          *name,
                            cs_ibm_algo_type_t   method,
                            const char          *file_name,
                            bool                 solve_fsi)
{
  CS_UNUSED(solve_fsi);

  int obj_id = _add_ibm_object(name,
                               method);

  cs_ibm_object_t *obj = cs_ibm_object_by_id(obj_id);

  /* STL/MED objects are only rigid for the moment! */
  /*obj->solve_fsi = solve_fsi;
  // TODO create a different function?
  if (solve_fsi) {
    BFT_MALLOC(obj->fsi_index, 1, int);
    obj->fsi_index[0] = cs_fsi_object->number;
    cs_fsi_object->number += 1;
  }*/

  /* STL */
  if (method == CS_IBM_ALGO_STL) {
    obj->stl = cs_stl_mesh_add(name);
    cs_stl_file_read(obj->stl, file_name);
    obj->stl->is_porous = true;
  } else if (method == CS_IBM_ALGO_MEDCOUPLING) {
    const char *sel_crit = "all[]";
    const char *intersect_method   = "P0P0";

    cs_medcoupling_intersector_add_vol(name,
                                       file_name,
                                       intersect_method,
                                       sel_crit);
    obj->mi = cs_medcoupling_intersector_by_name(name);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an object from a function used in the cutcell algorithm
 *
 * \param[in] name          name of the object
 * \param[in] cutcell_func  pointer to the cutcell function of the object
 * \param[in] solve_fsi     Is the object used in the FSI resolution ?
 *                          (currently ignored)
 * \param[in] n_nodes       Number of nodes if the object is deformable
 *                          (currently ignored)
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_object_from_func(const char        *name,
                            cs_cutcell_func_t *cutcell_func,
                            bool               solve_fsi,
                            int                n_nodes)
{
  CS_UNUSED(solve_fsi);
  CS_UNUSED(n_nodes);

  int obj_id = _add_ibm_object(name,
                               CS_IBM_ALGO_CUT_CELLS);

  cs_ibm_object_t *obj = cs_ibm_object_by_id(obj_id);

  // TODO create a different function?
  /*obj->solve_fsi = solve_fsi;
  if (solve_fsi) {
    int id0 = cs_fsi_object->number;
    if (n_nodes > 1) {
      obj->is_deformable = true;
      BFT_MALLOC(obj->fsi_index, n_nodes - 1, int);
      for (int i = 0; i < n_nodes - 1; i++)
        obj->fsi_index[i] = id0 + i;

      cs_fsi_object->number  += n_nodes - 1;

      obj->deform_id = cs_fsi_object->n_solid;
      cs_fsi_object->n_solid += 1;
      obj->n_nodes = n_nodes;

    } else {
      obj->is_deformable = false;
      BFT_MALLOC(obj->fsi_index, 1, int);
      obj->fsi_index[0] = id0;

      cs_fsi_object->number += 1;

    }
  }*/

  obj->cutcell_func = cutcell_func;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define exterior points for an stl object.
 *
 * \param[in] name          name of the object
 * \param[in] n_pts         number of points
 * \param[in] pts_coords    coordinates of the points
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_stl_define_ext_points(const char      *name,
                             const int        n_pts,
                             cs_real_t       *pts_coords)
{

  cs_ibm_object_t *obj = cs_ibm_object_by_name(name);

  if (obj->method != CS_IBM_ALGO_STL)
    bft_error(__FILE__, __LINE__, 0,
              _("You can't define exterior points to a non stl object.\n"));

  cs_stl_set_porosity_seed(obj->stl,
                           n_pts,
                           pts_coords);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotate an object based on the STL or MED algorithms
 *
 * \param[in] name          name of the object
 * \param[in] angle         angle of rotation
 * \param[in] axis          axis of rotation
 * \param[in] center        center of rotation
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_rotate(const char *name,
                     cs_real_t   angle,
                     cs_real_t   axis[3],
                     cs_real_t   center[3])
{
  cs_ibm_object_t *obj = cs_ibm_object_by_name(name);

  switch(obj->method) {

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      cs_medcoupling_intersector_rotate(obj->mi,
                                        center,
                                        axis,
                                        angle);
    }
    break;

  case CS_IBM_ALGO_STL:
    {
      cs_stl_mesh_rotate(obj->stl,
                         angle,
                         axis,
                         center);
    }
    break;

  default:
    {
      bft_error(__FILE__, __LINE__, 0,
                "Object %s definition method is neither STL nor MEDCoupling\n",
                name);
    }
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new constant property definition for an object.
 *
 * \param[in] obj       pointer to object
 * \param[in] ppty_id   property id (si enum for list)
 * \param[in] val       property constant value
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_set_property_const(cs_ibm_object_t               *obj,
                                 cs_ibm_object_property_type_t  ppty_id,
                                 cs_real_t                      val)
{
  _ibm_object_define_property_def(obj, ppty_id, 1, &val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Translate an object based on the STL or MED algorithms
 *
 * \param[in] name          name of the object
 * \param[in] vector        translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_translate(const char *name,
                        cs_real_t   vector[3])
{
  cs_ibm_object_t *obj = cs_ibm_object_by_name(name);

  switch(obj->method) {

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      cs_medcoupling_intersector_translate(obj->mi, vector);
    }
    break;

  case CS_IBM_ALGO_STL:
    {
      cs_stl_mesh_translate(obj->stl, vector);
    }
    break;

  default:
    {
      bft_error(__FILE__, __LINE__, 0,
                "Object %s was not defined using MEDCoupling or STL\n", name);
    }
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scale an object based on a factor
 *
 * \param[in] name          name of the object
 * \param[in] factor        scaling factor
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_scale(const char *name,
                    cs_real_t   factor)
{
  cs_ibm_object_t *obj = cs_ibm_object_by_name(name);

  switch(obj->method) {

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      cs_medcoupling_intersector_scale_auto(obj->mi, factor);
    }
    break;

  case CS_IBM_ALGO_STL:
    {
      cs_stl_mesh_scale(obj->stl, factor);
    }
    break;

  default:
    {
      bft_error(__FILE__, __LINE__, 0,
                "Object %s was not defined using MEDCoupling or STL\n", name);
    }
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply user parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_user_parameters(void)
{
  cs_ibm->algo_choice = CS_IBM_ALGO_NONE;

  // TODO: Add gui call
  cs_user_ibm_parameters();

  // TODO: Add gui call
  cs_user_ibm_define_objects();
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Init writers for STL or MED objects.
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_init_writer(void)
{
  /* Get the time control from the default writer (writer_id = -1) */
  cs_time_control_t *tc = cs_post_get_time_control(-1);

  bool      output_at_start = tc->at_start;
  bool      output_at_end   = tc->at_end;
  int       interval_n      = tc->interval_nt;
  cs_real_t interval_t      = tc->interval_t;

  switch(cs_ibm->algo_choice) {

  case CS_IBM_ALGO_STL:
    {
      cs_stl_post_init_writer("STL_OBJECTS",
                              "postprocessing",
                              "Ensight Gold",
                              "",
                              FVM_WRITER_TRANSIENT_COORDS,
                              output_at_start,
                              output_at_end,
                              interval_n,
                              interval_t);
    }
    break;

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      cs_mi_post_init_writer("MED_OBJECTS",
                             "postprocessing",
                             "Ensight Gold",
                             "",
                             FVM_WRITER_TRANSIENT_COORDS,
                             output_at_start,
                             output_at_end,
                             interval_n,
                             interval_t);
    }

  /* Do nothing */
  default:
    break;

  }

  /* Attach meshes to writers */
  for (int i = 0; i < cs_ibm->n_objects; i++) {
    /* STL Objects */
    if (cs_ibm->objects[i]->method == CS_IBM_ALGO_STL)
      cs_stl_post_add_mesh(cs_ibm->objects[i]->stl);

    /* MED Objects */
    else if (cs_ibm->objects[i]->method == CS_IBM_ALGO_MEDCOUPLING)
      cs_mi_post_add_mesh(cs_ibm->objects[i]->mi);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transform an object from its initial state using a transformation
 *         matrix.
 *
 * \param[in] obj     pointer to object structure
 * \param[in] matrix  transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_transform_from_init(cs_ibm_object_t *obj,
                                  cs_real_34_t     matrix)
{
  switch(obj->method) {

  case CS_IBM_ALGO_STL:
    {
      cs_stl_mesh_transform_from_init(obj->stl, matrix);
    }
    break;

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      cs_medcoupling_intersector_transform_from_init(obj->mi, matrix);
    }
    break;

  default:
    break;
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the volume fraction of an object over all cells.
 *
 * \param[in]  obj            pointer to object structure
 * \param[in]  m              pointer to mesh structure
 * \param[in]  cell_vol       pointer to cell volume array
 * \param[out] obj_frac_tot   array containing the total vol fraction of solids
 * \param[in]  indic          indicator array (currently ignored)
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_compute_intersect_vol(cs_ibm_object_t            *obj,
                                    const cs_mesh_t            *m,
                                    const cs_real_t            *cell_vol,
                                    cs_real_t                  *obj_frac_tot,
                                    int                        *indic)
{
  CS_UNUSED(indic);

  cs_real_t *wfrac = NULL;
  cs_lnum_t *windic = NULL;
  BFT_MALLOC(wfrac , m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(windic, m->n_cells_with_ghosts, cs_lnum_t);

  switch(obj->method) {

  case CS_IBM_ALGO_STL:
    {
      CS_NO_WARN_IF_UNUSED(cell_vol);

      /* Compute porosity from updated STL*/
      cs_stl_compute_porosity(obj->stl,
                              wfrac,
                              windic);

      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
        cs_real_t obj_frac = 1. - wfrac[c_id];

        obj_frac_tot[c_id] += obj_frac;
      }

    }
    break;

  case CS_IBM_ALGO_MEDCOUPLING:
    {
      /* Calling the intersection routine for the "o_id"_th object.
       * the function takes as input the intersector (mi).
       * Result is absolute volume, which needs to be divided
       * by the cell volume to obtain the object volume fraction of
       * the cell.
       */
      cs_real_t *obj_volume = cs_medcoupling_intersect_volumes(obj->mi);

      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
        cs_real_t obj_frac = obj_volume[c_id] / cell_vol[c_id];

        obj_frac_tot[c_id] += obj_frac;
      }
      cs_real_t obj_frac_tot_clip = 0.;
      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
        obj_frac_tot_clip = obj_frac_tot[c_id];
        if (obj_frac_tot_clip > 0.999999) {
          obj_frac_tot_clip = 1.;
        }
        obj_frac_tot[c_id] = obj_frac_tot_clip;
      }

    }
    break;

  default:
    break;

  }

  BFT_FREE(wfrac);
  BFT_FREE(windic);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define space immersed boundaries on a set of zones defined by the user
 *         in the GUI.
 *
 * \param[in]  mesh_quantities  pointer to associated mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void cs_ibm_volumic_zone(const cs_mesh_quantities_t *mesh_quantities)
{
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mesh_quantities->cell_cen;
  int n_v_zones = cs_volume_zone_n_zones();

  cs_tree_node_t *tn_p = cs_tree_get_node(cs_glob_tree,
                                   "thermophysical_models/porosities/porosity");

  /* Loop on all zones */
  for (int i = 0; i < n_v_zones; i++) {
    const cs_zone_t *z = cs_volume_zone_by_id(i);

    /* If zone is defined as Porous zone we use the formulae defined by the
     * user in the GUI. */
    if (z->type & CS_VOLUME_ZONE_POROSITY) {

      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", i);
      cs_tree_node_t *tn_zp = cs_tree_node_get_sibling_with_tag(tn_p,
                                                                "zone_id",
                                                                z_id_str);

      const char *formula = cs_tree_node_get_child_value_str(tn_zp, "formula");

      if (formula != NULL) {
        cs_field_t *f = CS_F_(poro);
        cs_meg_volume_function(z->name,
                               z->n_elts,
                               z->elt_ids,
                               cell_cen,
                               f->name,
                               &(f->val));
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
