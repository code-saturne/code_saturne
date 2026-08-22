/*============================================================================
 * Define immersed boundaries based on user inputs (experimental).
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

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_writer.h"

#include "base/cs_array.h"
#include "alge/cs_cell_to_vertex.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_mem.h"
#include "base/cs_medcoupling_intersector.h"
#include "meg/cs_meg_prototypes.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_post.h"
#include "base/cs_prototypes.h"
#include "mesh/cs_stl.h"
#include "base/cs_turbomachinery.h"
#include "base/cs_velocity_pressure.h"
#include "alge/cs_vertex_to_cell.h"
#include "base/cs_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_ibm.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_ibm.cpp
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
cs_ibm_t  *cs_ibm = nullptr;

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
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the face porosity depending on the neighbouring cells
 *          porosities based on geometry.
 *
 * \param[in]  alphai       value at neighbouring cell i
 * \param[in]  alphaj       value at neighbouring cell j
 * \param[in]  sizei        distance face - center of gravity cell i
 * \param[in]  sizej        distance face - center of gravity cell j
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_geom_face_fraction(cs_real_t alphai,
                    cs_real_t alphaj,
                    cs_real_t sizei,
                    cs_real_t sizej)
{
  cs_real_t alpi = alphai;
  cs_real_t alpj = alphaj;
  cs_real_t li = sizei;
  cs_real_t lj = sizej;

  if (alphaj < alphai) {
    alpi = alphaj;
    alpj = alphai;
    li = sizej;
    lj = sizei;
  }

  cs_real_t li_ov_lj = li/lj;
  cs_real_t lj_ov_li = lj/li;

  cs_real_t a, b, x1, x2, y1, y2, aa, bb, cc, delta, b1, b2;
  cs_real_t rap;

  /* Small porosity treatment */
  cs_real_t eps = 1.e-10;
  if (alpi < 1.-eps && alpj  > 1.-eps)
    return alpj;
  if (alpi < eps && alpj  > eps)
    return alpi;

  /* Initialize the solution */
  rap = li_ov_lj;
  cs_real_t alpij = (alpi + rap*alpj) / (1. + rap);

  if (alpj < 1.e-10)
    return alpij;

  if (alpi > 1. - 1.e-10)
    return alpij;

  if (alpi < 1.e-3 && alpj > 1.-1.e-3) return alpij;

  /* Four crossing cases are possible : faces are cut each one on Y, one on X
   * one on Y, one on X the other one too, or Y and X. The minimax of
   * initialization of the alpi alpj guarantees the positivity of the slope. */

  /* Case YY */
  b = alpij;
  a = (b - alpi) * 2./li;
  y1 = b - a*li;
  y2 = b + a*lj;

  if (y1 >= 0. && y2 <= 1.)
    return b;

  /* Case YX */
  rap = lj_ov_li;
  aa = 1.;
  bb = -2. - 4.*(1. - alpj)*rap;
  cc = 1. + 4 * alpi*(1. - alpj)*rap;
  delta = bb*bb - 4.*aa*cc;

  if (delta >= 0.) {
    b = (-bb - sqrt(delta)) / (2.*aa);
    if (b >= alpi && b <= alpj) {
      a = (b - alpi)*2./li;
      if (a > 0.) {
        y1 = b - a*li;
        x2 = (1. - b)/a;
        if (y1 >= 0. && x2 <= lj)
          return b;
      }
    }
  }

  /* Case XX */
  aa = (1. - alpj)*lj - alpi*li;
  bb = 2.*alpi*li;
  cc = -alpi*li;
  delta = bb*bb - 4.*aa*cc;

  if (delta >= 0.) {
    if (cs::abs(aa) < 1.e-20) {
      return 0.5;
    }
    else {
      b1 = (-bb - sqrt(delta)) / (2.*aa);
      b2 = (-bb + sqrt(delta)) / (2.*aa);
      if (b1 >= alpi && b1 <= alpj) {
        a = b1*b1/(2.*alpi*li);
        if (a > 0) {
          x1 = -b1/a;
          x2 = (1. - b1)/a;
          if (x1 >= -li && x2 <= lj)
            return b1;
        }
      }
      if (b2 >= alpi && b2 <= alpj) {
        a = b2*b2/(2.*alpi*li);
        if (a > 0.) {
          x1 = -b2/a;
          x2 = (1. - b2)/a;
          if (x1 >= -li && x2 <= lj)
            return b2;
        }
      }
    }
  }

  /* Case XY */
  rap = li_ov_lj;
  aa = 1.;
  bb = 4.*alpi*rap;
  cc = -4.*alpi*alpj*rap;
  delta = bb*bb - 4.*aa*cc;
  if (delta >= 0.) {
    b = (-bb + sqrt(delta)) / (2.*aa);
    if (b >= alpi && b <= alpj) {
      a = (alpj - b)*2./lj;
      if (a > 0.) {
        y2 = b + a*lj;
        x1 = -b/a;
        if (x1 >= -li && y2 <= 1.)
          return b;
      }
    }
  }

  return alpij;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by dichotomy the length of the immersed part of a line
 *          between two points (i.e. the distance between the point in the
 *          solid and the immersed boundary) based on the cut-cell method.
 *
 * \param[in]  ipenal1      penalization
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_imm_lgth_cutcell(int         ipenal1,
                  cs_real_3_t x1,
                  cs_real_3_t x2,
                  cs_real_t   t,
                  int         num_object)
{
  cs_real_3_t xx1, xx2, x3;
  for (int idim = 0; idim < 3; idim++) {
    xx1[idim] = x1[idim];
    xx2[idim] = x2[idim];
  }

  cs_real_t l1 = 0;
  cs_real_t l2 = 1.;

  /* Method of dichotomy (be careful not to lower the value 20, pressure
   * oscillations otherwise) */
  constexpr int nsou = 20;
  for (int isou = 0; isou < nsou; isou++) {
    for (int idim = 0; idim < 3; idim++)
      x3[idim] = 0.5*(xx1[idim] + xx2[idim]);

    int ipenal3 = cs_ibm_object_compute_cut_porosity(3, x3, t, num_object);

    if (ipenal3 == ipenal1) {
      for (int idim = 0; idim < 3; idim++)
        xx1[idim] = x3[idim];
      l1 = 0.5*(l1 + l2);
    }
    else {
      for (int idim = 0; idim < 3; idim++)
        xx2[idim] = x3[idim];
      l2 = 0.5*(l1 + l2);
    }
  }

  if (ipenal1 == 1)
    return 1. - l1;
  else
    return l1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute by dichotomy the length of the immersed part of a line
 *          between two points (i.e. the distance between the point in the
 *          solid and the immersed boundary) based on the input porosities.
 *
 * \param[in]  por1         porosity at 1
 * \param[in]  por2         porosity at 2
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_imm_lgth_poro(cs_real_t   por1,
               cs_real_t   por2)
{
  /* Both considered as solids, returns 0 */
  if (por1 < 0.5 && por2 < 0.5)
    return 0.;

  /* Both considered as fluids, returns total distance */
  else if (por1 > 0.5 && por2 > 0.5)
    return 1.;

  /* Fluid-solid cases, return weighted distance  */
  else if (por1 < 0.5 && por2 >= 0.5)
    return (por2 - 0.5) / (por2 - por1);

  else if (por2 < 0.5 && por1 >= 0.5)
    return (por1 - 0.5) / (por1 - por2);

  /* Medium case */
  else
    return 0.5;
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
  cs_real_3_t dx14, dx24, dx34;
  for (int i = 0; i < 3; i++) {
    dx14[i] = x1[i] - x4[i];
    dx24[i] = x2[i] - x4[i];
    dx34[i] = x3[i] - x4[i];
  }

  cs_real_t tetra_vol
  = cs::abs(cs_math_3_triple_product(dx14, dx24, dx34)) * cs_math_1ov6;

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

static cs_real_t
_pyram_vol_cog(cs_real_3_t x1,
               cs_real_3_t x2,
               cs_real_3_t x3,
               cs_real_3_t x4,
               cs_real_3_t x5,
               cs_real_3_t cogvol)
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

  for (int i = 0; i < 3; i++) {
    cogvol[i] = 0.25 * (x1[i] + x2[i] + xc[i] + x5[i]) * vol12
              + 0.25 * (x2[i] + x3[i] + xc[i] + x5[i]) * vol23
              + 0.25 * (x3[i] + x4[i] + xc[i] + x5[i]) * vol34
              + 0.25 * (x4[i] + x1[i] + xc[i] + x5[i]) * vol41;
  }

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

  cs_real_t xc[3];

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
 * \brief  Compute the volume of a prism with a quadrangle base
 *          (x3,x4,x5,x6) and edge (x1 x2), and its center of gravity.
 *
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 * \param[in]  x5           point 5
 * \param[in]  x6           point 6
 * \param[in]  cogvol       center of gravity
 */
/*----------------------------------------------------------------------------*/


static cs_real_t
_prism_vol_cog
(
  cs_real_3_t x1,    /*!<[in] */
  cs_real_3_t x2,    /*!<[in] */
  cs_real_3_t x3,    /*!<[in] */
  cs_real_3_t x4,    /*!<[in] */
  cs_real_3_t x5,    /*!<[in] */
  cs_real_3_t x6,    /*!<[in] */
  cs_real_3_t cogvol /*!<[out] */
)
{
  /* Volume of the prism of base x3, x4, x5, x6 and edge x1 x2 */
  /* The two triangles are x1 x3 x6 and x2 x4 x5 */
  cs_real_3_t xc;

  for (int i = 0; i < 3; i++)
    xc[i] = (x1[i]+x2[i]+x3[i]+x4[i]+x5[i]+x6[i]) * cs_math_1ov6;

  cs_real_3_t cogvolpyram;

  for (int i = 0; i < 3; i++)
    cogvol[i] = 0.;

  cs_real_t vol136c  = _tetra_vol(x1, x3, x6, xc);
  for (int i = 0; i < 3; i++)
    cogvol[i] += 0.25 * (x1[i] + x3[i] + x6[i] + xc[i]) * vol136c;

  cs_real_t vol245c  = _tetra_vol(x2, x4, x5, xc);
  for (int i = 0; i < 3; i++)
    cogvol[i] += 0.25 * (x2[i] + x4[i] + x5[i] + xc[i]) * vol245c;

  cs_real_t vol1256c = _pyram_vol_cog(x1, x2, x5, x6, xc, cogvolpyram);
  for (int i = 0; i < 3; i++)
    cogvol[i] += cogvolpyram[i];

  cs_real_t vol1243c = _pyram_vol_cog(x1, x2, x4, x3, xc, cogvolpyram);
  for (int i = 0; i < 3; i++)
    cogvol[i] += cogvolpyram[i];

  cs_real_t vol3456c = _pyram_vol_cog(x3, x4, x5, x6, xc, cogvolpyram);
  for (int i = 0; i < 3; i++)
    cogvol[i] += cogvolpyram[i];

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
_tetra_vol_poro(cs_lnum_t    numcel,
                cs_real_t   &vol,
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

  vol = 0.;

  /* Check if the volume is considered as completely solid, null volume */
  if (cpt5 == 5)
    return;

  /* Check if no solid vertex, complete volume */
  else if (cpt5 == 0 && icut > 0) {

    vol = _tetra_vol(x1, x2, x3, x4);

    for (int i = 0; i < 3; i++)
      cog[i] += 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]) * vol;

    return;

  /* Check if at least one solid vertex, subdivision of the edges and
   * distribution of porosities before recursive call */
  }
  else if (icut > 0) {
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
    _tetra_vol_poro(numcel, vol1, cog, x1,  por1,  x12, por12, x13, por13,
                    x14, por14, icut);
    _tetra_vol_poro(numcel, vol2, cog, x2,  por2,  x23, por23, x12, por12,
                    x24, por24, icut);
    _tetra_vol_poro(numcel, vol3, cog, x3,  por3,  x13, por13, x23, por23,
                    x34, por34, icut);
    _tetra_vol_poro(numcel, vol4, cog, x4,  por4,  x14, por14, x24, por24,
                    x34, por34, icut);
    _tetra_vol_poro(numcel, vol5, cog, x12, por12, x13, por13, x34, por34,
                    x23, por23, icut);
    _tetra_vol_poro(numcel, vol6, cog, x12, por12, x13, por13, x34, por34,
                    x14, por14, icut);
    _tetra_vol_poro(numcel, vol7, cog, x34, por34, x24, por24, x12, por12,
                    x23, por23, icut);
    _tetra_vol_poro(numcel, vol8, cog, x34, por34, x24, por24, x12, por12,
                    x14, por14, icut);

    vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
    return;

  }
  else {
    if (cpt4 == 0) {

      /* No penalized point -> full tetrahedron */
      vol = _tetra_vol(x1, x2, x3, x4);

      for (int i = 0; i < 3; i++)
        cog[i] += 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]) * vol;

    }
    else if (cpt4 == 1) {

      /* One penalized point tetrahedron - tetra from the penalized point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 1) {
        cs_real_t lbd12 = 1. - _imm_lgth_poro(por2, por1);
        cs_real_t lbd13 = 1. - _imm_lgth_poro(por3, por1);
        cs_real_t lbd14 = 1. - _imm_lgth_poro(por4, por1);

        cs_real_t volp = vol * lbd12 * lbd13 * lbd14;

        cs_real_3_t x12, x13, x14;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x12[i] = x1[i] + lbd12 * d12[i];
          x13[i] = x1[i] + lbd13 * d13[i];
          x14[i] = x1[i] + lbd14 * d14[i];
          cogp[i] = 0.25 * (x1[i] + x12[i] + x13[i] + x14[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * vol - cogp[i] * volp;
        }

        vol *= (1. - lbd12 * lbd13 * lbd14);

        return;

      }
      else if (ipenal2 == 1) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t lbd21 = 1. - _imm_lgth_poro(por1, por2);
        cs_real_t lbd23 = 1. - _imm_lgth_poro(por3, por2);
        cs_real_t lbd24 = 1. - _imm_lgth_poro(por4, por2);

        cs_real_t volp = vol * lbd21 * lbd23 * lbd24;

        cs_real_3_t x21, x23, x24;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x21[i] = x2[i] + lbd21 * d21[i];
          x23[i] = x2[i] + lbd23 * d23[i];
          x24[i] = x2[i] + lbd24 * d24[i];
          cogp[i] = 0.25 * (x2[i] + x21[i] + x23[i] + x24[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * vol - cogp[i] * volp;
        }

        vol *= (1. - lbd21 * lbd23 * lbd24);

        return;

      }
      else if (ipenal3 == 1) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t lbd31 = 1. - _imm_lgth_poro(por1, por3);
        cs_real_t lbd32 = 1. - _imm_lgth_poro(por2, por3);
        cs_real_t lbd34 = 1. - _imm_lgth_poro(por4, por3);

        cs_real_t volp = vol * lbd31 * lbd32 * lbd34;

        cs_real_3_t x31, x32, x34;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x31[i] = x3[i] + lbd31 * d31[i];
          x32[i] = x3[i] + lbd32 * d32[i];
          x34[i] = x3[i] + lbd34 * d34[i];
          cogp[i] = 0.25 * (x3[i] + x31[i] + x32[i] + x34[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * vol - cogp[i] * volp;
        }

        vol *= (1. - lbd31 * lbd32 * lbd34);

        return;

      }
      else if (ipenal4 == 1) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t lbd41 = 1. - _imm_lgth_poro(por1, por4);
        cs_real_t lbd42 = 1. - _imm_lgth_poro(por2, por4);
        cs_real_t lbd43 = 1. - _imm_lgth_poro(por3, por4);

        cs_real_t volp = vol * lbd41 * lbd42 * lbd43;

        cs_real_3_t x41, x42, x43;
        cs_real_3_t cogp, cogtot;
        for (int i = 0; i < 3; i++) {
          x41[i] = x4[i] + lbd41 * d41[i];
          x42[i] = x4[i] + lbd42 * d42[i];
          x43[i] = x4[i] + lbd43 * d43[i];
          cogp[i] = 0.25 * (x4[i] + x41[i] + x42[i] + x43[i]);
          cogtot[i] = 0.25 * (x1[i] + x2[i] + x3[i] + x4[i]);
          cog[i] += cogtot[i] * vol - cogp[i] * volp;
        }

        vol *= (1. - lbd41 * lbd42 * lbd43);

        return;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_poro\n");

    }
    else if (cpt4 == 3) {

      /* Three penalized points : tetrahedron from the non penalized
       * fourth point */
      cs_real_3_t d12, d13, d14;
      for (int i = 0; i < 3; i++) {
        d12[i] = x2[i] - x1[i];
        d13[i] = x3[i] - x1[i];
        d14[i] = x4[i] - x1[i];
      }

      vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 0) {
        cs_real_t lbd12 = _imm_lgth_poro(por1, por2 );
        cs_real_t lbd13 = _imm_lgth_poro(por1, por3 );
        cs_real_t lbd14 = _imm_lgth_poro(por1, por4 );

        vol *= lbd12 * lbd13 * lbd14;

        cs_real_3_t x12, x13, x14;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x12[i] = x1[i] + lbd12 * d12[i];
          x13[i] = x1[i] + lbd13 * d13[i];
          x14[i] = x1[i] + lbd14 * d14[i];
          cogp[i] = 0.25 * (x1[i] + x12[i] + x13[i] + x14[i]);
          cog[i] += cogp[i] * vol;
        }

        return;

      }
      else if (ipenal2 == 0) {
        cs_real_3_t d21, d23, d24;
        for (int i = 0; i < 3; i++) {
          d21[i] = x1[i] - x2[i];
          d23[i] = x3[i] - x2[i];
          d24[i] = x4[i] - x2[i];
        }

        cs_real_t lbd21 = _imm_lgth_poro(por2, por1);
        cs_real_t lbd23 = _imm_lgth_poro(por2, por3);
        cs_real_t lbd24 = _imm_lgth_poro(por2, por4);

        vol *= lbd21 * lbd23 * lbd24;

        cs_real_3_t x21, x23, x24;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x21[i] = x2[i] + lbd21 * d21[i];
          x23[i] = x2[i] + lbd23 * d23[i];
          x24[i] = x2[i] + lbd24 * d24[i];
          cogp[i] = 0.25 * (x2[i] + x21[i] + x23[i] + x24[i]);
          cog[i] += cogp[i] * vol;
        }

        return;

      }
      else if (ipenal3 == 0) {
        cs_real_3_t d31, d32, d34;
        for (int i = 0; i < 3; i++) {
          d31[i] = x1[i] - x3[i];
          d32[i] = x2[i] - x3[i];
          d34[i] = x4[i] - x3[i];
        }

        cs_real_t lbd31 = _imm_lgth_poro(por3, por1);
        cs_real_t lbd32 = _imm_lgth_poro(por3, por2);
        cs_real_t lbd34 = _imm_lgth_poro(por3, por4);

        vol *= lbd31 * lbd32 * lbd34;

        cs_real_3_t x31, x32, x34;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x31[i] = x3[i] + lbd31 * d31[i];
          x32[i] = x3[i] + lbd32 * d32[i];
          x34[i] = x3[i] + lbd34 * d34[i];
          cogp[i] = 0.25 * (x3[i] + x31[i] + x32[i] + x34[i]);
          cog[i] += cogp[i] * vol;
        }

        return;

      }
      else if (ipenal4 == 0) {
        cs_real_3_t d41, d42, d43;
        for (int i = 0; i < 3; i++) {
          d41[i] = x1[i] - x4[i];
          d42[i] = x2[i] - x4[i];
          d43[i] = x3[i] - x4[i];
        }

        cs_real_t lbd41 = _imm_lgth_poro(por4, por1);
        cs_real_t lbd42 = _imm_lgth_poro(por4, por2);
        cs_real_t lbd43 = _imm_lgth_poro(por4, por3);

        vol *= lbd41 * lbd42 * lbd43;

        cs_real_3_t x41, x42, x43;
        cs_real_3_t cogp;
        for (int i = 0; i < 3; i++) {
          x41[i] = x4[i] + lbd41 * d41[i];
          x42[i] = x4[i] + lbd42 * d42[i];
          x43[i] = x4[i] + lbd43 * d43[i];
          cogp[i] = 0.25 * (x4[i] + x41[i] + x42[i] + x43[i]);
          cog[i] += cogp[i] * vol;
        }

        return;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    }
    else if (cpt4 == 2) {

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

      }
      else if (ipenal1 == 0 && ipenal3 == 0) {
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

      }
      else if (ipenal1 == 0 && ipenal4 == 0) {
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

      }
      else if (ipenal2 == 0 && ipenal3 == 0) {
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

      }
      else if (ipenal2 == 0 && ipenal4 == 0) {
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

      }
      else if (ipenal3 == 0 && ipenal4 == 0) {
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

      cs_real_t lbd13 = _imm_lgth_poro(porr1, porr3);
      cs_real_t lbd14 = _imm_lgth_poro(porr1, porr4);
      cs_real_t lbd23 = _imm_lgth_poro(porr2, porr3);
      cs_real_t lbd24 = _imm_lgth_poro(porr2, porr4);

      cs_real_3_t y13, y14, y23, y24;
      for (int i = 0; i < 3; i++) {
        y13[i] = y1[i] + lbd13 * d13[i];
        y14[i] = y1[i] + lbd14 * d14[i];
        y23[i] = y2[i] + lbd23 * d23[i];
        y24[i] = y2[i] + lbd24 * d24[i];
      }

      cs_real_3_t cogp;
      vol = _prism_vol_cog(y1, y2, y13, y23, y24, y14, cogp);

      for (int i = 0; i < 3; i++)
        cog[i] += cogp[i];
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
  constexpr cs_real_t c_1ov3 = 1./3.;
  cs_real_t por4 = (por1 + por2 + por3) * c_1ov3;

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
    for (int i = 0; i < 3; i++) {
      dx12[i] = x2[i] - x1[i];
      dx13[i] = x3[i] - x1[i];
    }

    cs_real_3_t pvec;
    cs_math_3_cross_product(dx12, dx13, pvec);

    surf = 0.5 * cs_math_3_norm(pvec);
    return surf;

  /* Check if at least one solid vertex, subdivision of the edges and
   * distribution of porosities before recursive call */
  }
  else if (icut > 0) {
    cs_real_3_t x12, x23, x13;
    icut--;

    for (int i = 0; i < 3; i++) {
      x12[i] = 0.5*(x1[i] + x2[i]);
      x13[i] = 0.5*(x1[i] + x3[i]);
      x23[i] = 0.5*(x2[i] + x3[i]);
    }

    cs_real_t por12 = 0.5*(por1 + por2);
    cs_real_t por13 = 0.5*(por1 + por3);
    cs_real_t por23 = 0.5*(por2 + por3);

    surf += _tri_surf_trunc(x1,  por1,  x12, por12, x13, por13, icut);
    surf += _tri_surf_trunc(x2,  por2,  x12, por12, x23, por23, icut);
    surf += _tri_surf_trunc(x3,  por3,  x13, por13, x23, por23, icut);
    surf += _tri_surf_trunc(x12, por12, x13, por13, x23, por23, icut);

    return surf;

  }
  else {
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
        cs_math_3_cross_product(dx12, dx13, pvec);

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t lbd12 = _imm_lgth_poro(por2, por1);
        cs_real_t lbd13 = _imm_lgth_poro(por3, por1);

        surf = surft * (1. - lbd12 * lbd13);

        return surf;

      }
      else if (ipenal2 == 1) {

        cs_real_3_t dx12, dx23;
        for (int i = 0; i < 3; i++) {
          dx12[i] = x2[i] - x1[i];
          dx23[i] = x3[i] - x2[i];
        }

        cs_real_3_t pvec;
        cs_math_3_cross_product(dx12, dx23, pvec);

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t lbd12 = _imm_lgth_poro(por2, por1);
        cs_real_t lbd23 = _imm_lgth_poro(por3, por2);

        surf = surft * (1. - lbd12 * lbd23);

        return surf;

      }
      else {

        cs_real_3_t dx13, dx23;
        for (int i = 0; i < 3; i++) {
          dx13[i] = x3[i] - x1[i];
          dx23[i] = x3[i] - x2[i];
        }

        cs_real_3_t pvec;
        cs_math_3_cross_product(dx13, dx23, pvec);

        cs_real_t surft = 0.5 * cs_math_3_norm(pvec);

        cs_real_t lbd13 = _imm_lgth_poro(por3, por1);
        cs_real_t lbd23 = _imm_lgth_poro(por3, por2);

        surf = surft * (1. - lbd13 * lbd23);

        return surf;
      }
    }
    else if (cpt3 == 2) {
      if (ipenal1 == 0) {
        cs_real_t lamb12 = _imm_lgth_poro(por1, por2);
        cs_real_t lamb13 = _imm_lgth_poro(por1, por3);

        cs_real_3_t dx12, dx13;
        for (int i = 0; i < 3; i++) {
          dx12[i] = lamb12 * (x2[i] - x1[i]);
          dx13[i] = lamb13 * (x3[i] - x1[i]);
        }

        cs_real_3_t pvec;
        cs_math_3_cross_product(dx12, dx13, pvec);

        surf = 0.5 * cs_math_3_norm(pvec);

        return surf;

      }
      else if (ipenal2 == 0) {
        cs_real_t lamb12 = _imm_lgth_poro(por1, por2);
        cs_real_t lamb23 = _imm_lgth_poro(por2, por3);

        cs_real_3_t dx12, dx23;
        for (int i =  0; i < 3; i++) {
          dx12[i] = lamb12 * (x2[i] - x1[i]);
          dx23[i] = lamb23 * (x3[i] - x2[i]);
        }

        cs_real_3_t pvec;
        cs_math_3_cross_product(dx12, dx23, pvec);

        surf = 0.5 * cs_math_3_norm(pvec);

        return surf;

      }
      else {
        cs_real_t lamb13 = _imm_lgth_poro(por1, por3);
        cs_real_t lamb23 = _imm_lgth_poro(por2, por3);

        cs_real_3_t dx13, dx23;
        for (int i = 0; i < 3; i++) {
          dx13[i] = lamb13 * (x3[i] - x1[i]);
          dx23[i] = lamb23 * (x3[i] - x2[i]);
        }

        cs_real_3_t pvec;
        cs_math_3_cross_product(dx13, dx23, pvec);

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
 * \param[in]  numcel       volume
 * \param[out] vol          volume
 * \param[in]  x1           point 1
 * \param[in]  x2           point 2
 * \param[in]  x3           point 3
 * \param[in]  x4           point 4
 * \param[in]  t            time value for the current time step
 * \param[in]  icut         number of bisections
 * \param[in]  nraf         number of refinement
 * \param[in]  lambda       number of refinement
 * \param[in]  lambda2      number of refinement
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

static void
_tetra_vol_cutcell(cs_lnum_t     numcel,
                   cs_real_t    &vol,
                   cs_real_3_t   x1,
                   cs_real_3_t   x2,
                   cs_real_3_t   x3,
                   cs_real_3_t   x4,
                   cs_real_t     t,
                   int           icut,
                   int           nraf,
                   cs_span_2d<cs_real_t>  lambda,
                   cs_span_2d<cs_real_t>  lambda2,
                   int           num_object)
{
  vol = 0.;

  /* If sub cut activated */
  if (icut > 0) {

    cs_real_3_t xnd[4];
    for (int i = 0; i < 3; i++) {
      xnd[0][i] = x1[i];
      xnd[1][i] = x2[i];
      xnd[2][i] = x3[i];
      xnd[3][i] = x4[i];
    }

    int in = 0;
    int inmax = 0;

    /* Loop on the vertices of the tetrahedron */
    for (int pt = 0; pt < 4; pt++) {
      in += cs_ibm_object_compute_cut_porosity(1, xnd[pt], t, num_object);
      inmax++;

      if (in > 0 && in < inmax) {
        icut--;
        cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
        cs_real_3_t xca[6];

        for (int i = 0; i < 3; i++) {
          xca[0][i] = 0.5*(x1[i] + x2[i]);
          xca[1][i] = 0.5*(x1[i] + x3[i]);
          xca[2][i] = 0.5*(x1[i] + x4[i]);
          xca[3][i] = 0.5*(x2[i] + x3[i]);
          xca[4][i] = 0.5*(x2[i] + x4[i]);
          xca[5][i] = 0.5*(x3[i] + x4[i]);
        }

        _tetra_vol_cutcell(numcel, vol1, x1,  xca[0], xca[1], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol2, x2,  xca[3], xca[0], xca[4], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol3, x3,  xca[1], xca[3], xca[5], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol4, x4,  xca[2], xca[4], xca[5], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol5, xca[0], xca[1], xca[5], xca[3], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol6, xca[0], xca[1], xca[5], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol7, xca[5], xca[4], xca[0], xca[3], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol8, xca[5], xca[4], xca[0], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);

        vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
        return;
      }
    }

    /* Center of gravity of the tetrahedron */
    cs_real_3_t cog1234;
    for (int i = 0; i < 3; i++)
      cog1234[i] = (x1[i] + x2[i] + x3[i] + x4[i])*0.25;

    in += cs_ibm_object_compute_cut_porosity(1, cog1234, t, num_object);
    inmax++;

    if (in > 0 && in < inmax) {
      icut--;
      cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
      cs_real_3_t xca[6];

      for (int i = 0; i < 3; i++) {
        xca[0][i] = 0.5*(x1[i] + x2[i]);
        xca[1][i] = 0.5*(x1[i] + x3[i]);
        xca[2][i] = 0.5*(x1[i] + x4[i]);
        xca[3][i] = 0.5*(x2[i] + x3[i]);
        xca[4][i] = 0.5*(x2[i] + x4[i]);
        xca[5][i] = 0.5*(x3[i] + x4[i]);
      }

      _tetra_vol_cutcell(numcel, vol1, x1,  xca[0], xca[1], xca[2], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol2, x2,  xca[3], xca[0], xca[4], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol3, x3,  xca[1], xca[3], xca[5], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol4, x4,  xca[2], xca[4], xca[5], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol5, xca[0], xca[1], xca[5], xca[3], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol6, xca[0], xca[1], xca[5], xca[2], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol7, xca[5], xca[4], xca[0], xca[3], t,
                         icut, nraf, lambda, lambda2, num_object);
      _tetra_vol_cutcell(numcel, vol8, xca[5], xca[4], xca[0], xca[2], t,
                         icut, nraf, lambda, lambda2, num_object);

      vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
      return;
    }

    /* Loop on the edges */
    if (true) {

      cs_real_3_t xdeb[6], xfin[6];
      cs_real_t xa[3];

      for (int i = 0; i < 3; i++) {
        xdeb[0][i] = x1[i];
        xfin[0][i] = x2[i];
        xdeb[1][i] = x1[i];
        xfin[1][i] = x3[i];
        xdeb[2][i] = x1[i];
        xfin[2][i] = x4[i];
        xdeb[3][i] = x2[i];
        xfin[3][i] = x3[i];
        xdeb[4][i] = x2[i];
        xfin[4][i] = x4[i];
        xdeb[5][i] = x3[i];
        xfin[5][i] = x4[i];
      }

      for (int edge = 0; edge < 6; edge++) {
        for (int iraf = 1; iraf < cs_ibm->nb_cut_edges; iraf++) {
          for (int i = 0; i < 3; i++)
            xa[i] = lambda(iraf, 0) * xdeb[edge][i]
                  + lambda(iraf, 1) * xfin[edge][i];

          in += cs_ibm_object_compute_cut_porosity(1, xa, t, num_object);
          inmax++;

          if (in > 0 && in < inmax) {
            icut--;
            cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
            cs_real_3_t xca[6];

            for (int i = 0; i < 3; i++) {
              xca[0][i] = 0.5*(x1[i] + x2[i]);
              xca[1][i] = 0.5*(x1[i] + x3[i]);
              xca[2][i] = 0.5*(x1[i] + x4[i]);
              xca[3][i] = 0.5*(x2[i] + x3[i]);
              xca[4][i] = 0.5*(x2[i] + x4[i]);
              xca[5][i] = 0.5*(x3[i] + x4[i]);
            }
            _tetra_vol_cutcell(numcel, vol1, x1,  xca[0], xca[1], xca[2], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol2, x2,  xca[3], xca[0], xca[4], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol3, x3,  xca[1], xca[3], xca[5], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol4, x4,  xca[2], xca[4], xca[5], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol5, xca[0], xca[1], xca[5], xca[3], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol6, xca[0], xca[1], xca[5], xca[2], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol7, xca[5], xca[4], xca[0], xca[3], t,
                               icut, nraf, lambda, lambda2, num_object);
            _tetra_vol_cutcell(numcel, vol8, xca[5], xca[4], xca[0], xca[2], t,
                               icut, nraf, lambda, lambda2, num_object);

            vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
            return;
          }
        }
      }
    }

    if (in == 0)
      vol = _tetra_vol(x1, x2, x3, x4);
    return;

  }
  else {

    int ipenal1 = cs_ibm_object_compute_cut_porosity(1, x1, t, num_object);
    int ipenal2 = cs_ibm_object_compute_cut_porosity(1, x2, t, num_object);
    int ipenal3 = cs_ibm_object_compute_cut_porosity(1, x3, t, num_object);
    int ipenal4 = cs_ibm_object_compute_cut_porosity(1, x4, t, num_object);

    // Last test before volume computation
    if (icut == 0) {

      cs_real_3_t x;
      int pb = 0;

      /* pb 12 */
      if (ipenal1 == ipenal2) {
        for (int k = 1; k < nraf; k++) {

          for (int i = 0; i < 3; i++)
            x[i] = lambda2(k, 0) * x1[i] + lambda2(k, 1) * x2[i];

          if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal1) {
            pb++;
            break;
          }
        }
      }

      /* pb 13 */
      if (pb == 0) {
        if (ipenal1 == ipenal3) {
          for (int k = 1; k < nraf; k++) {

            for (int i = 0; i < 3; i++)
              x[i] = lambda2(k, 0) * x1[i] + lambda2(k, 1) * x3[i];

            if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal1) {
              pb++;
              break;
            }
          }
        }

        /* pb 14 */
        if (ipenal1 == ipenal4) {
          for (int k = 1; k < nraf; k++) {

            for (int i = 0; i < 3; i++)
              x[i] = lambda2(k, 0) * x1[i] + lambda2(k, 1) * x4[i];

            if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal1) {
              pb++;
              break;
            }
          }
        }

        /* pb 23 */
        if (ipenal2 == ipenal3) {
          for (int k = 1; k < nraf; k++) {

            for (int i = 0; i < 3; i++)
              x[i] = lambda2(k, 0) * x2[i] + lambda2(k, 1) * x3[i];

            if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal2) {
              pb++;
              break;
            }
          }
        }

        /* pb 24 */
        if (ipenal2 == ipenal4) {
          for (int k = 1; k < nraf; k++) {

            for (int i = 0; i < 3; i++)
              x[i] = lambda2(k, 0) * x2[i] + lambda2(k, 1) * x4[i];

            if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal2) {
              pb++;
              break;
            }
          }
        }

        /* pb 34 */
        if (ipenal3 == ipenal4) {
          for (int k = 1; k < nraf; k++) {

            for (int i = 0; i < 3; i++)
              x[i] = lambda2(k, 0) * x3[i] + lambda2(k, 1) * x4[i];

            if (cs_ibm_object_compute_cut_porosity(1, x, t, num_object) != ipenal3) {
              pb++;
              break;
            }
          }
        }
      }

      if (pb > 0) {
        icut--;
        cs_real_t vol1, vol2, vol3, vol4, vol5, vol6, vol7, vol8;
        cs_real_3_t xca[6];

        for (int i = 0; i < 3; i++) {
          xca[0][i] = 0.5*(x1[i] + x2[i]);
          xca[1][i] = 0.5*(x1[i] + x3[i]);
          xca[2][i] = 0.5*(x1[i] + x4[i]);
          xca[3][i] = 0.5*(x2[i] + x3[i]);
          xca[4][i] = 0.5*(x2[i] + x4[i]);
          xca[5][i] = 0.5*(x3[i] + x4[i]);
        }
        _tetra_vol_cutcell(numcel, vol1, x1,  xca[0], xca[1], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol2, x2,  xca[3], xca[0], xca[4], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol3, x3,  xca[1], xca[3], xca[5], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol4, x4,  xca[2], xca[4], xca[5], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol5, xca[0], xca[1], xca[5], xca[3], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol6, xca[0], xca[1], xca[5], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol7, xca[5], xca[4], xca[0], xca[3], t,
                           icut, nraf, lambda, lambda2, num_object);
        _tetra_vol_cutcell(numcel, vol8, xca[5], xca[4], xca[0], xca[2], t,
                           icut, nraf, lambda, lambda2, num_object);

        vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8;
        return;
      }
    }

    int cpt4 = ipenal1 + ipenal2 + ipenal3 + ipenal4;

    if (cpt4 == 0) {
      /* No penalized point -> full tetrahedron */
      vol = _tetra_vol(x1, x2, x3, x4);

    } else if (cpt4 == 4) {
      vol = 0;

    } else if (cpt4 == 1) {
      /* One penalized point tetrahedron - tetra from the penalized point */
      vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 1) {
        cs_real_t lbd12 = 1.-_imm_lgth_cutcell(ipenal1, x1, x2, t, num_object);
        cs_real_t lbd13 = 1.-_imm_lgth_cutcell(ipenal1, x1, x3, t, num_object);
        cs_real_t lbd14 = 1.-_imm_lgth_cutcell(ipenal1, x1, x4, t, num_object);

        vol *= (1. - lbd12 * lbd13 * lbd14);

        return;

      }
      else if (ipenal2 == 1) {
        cs_real_t lbd21 = 1.-_imm_lgth_cutcell(ipenal2, x2, x1, t, num_object);
        cs_real_t lbd23 = 1.-_imm_lgth_cutcell(ipenal2, x2, x3, t, num_object);
        cs_real_t lbd24 = 1.-_imm_lgth_cutcell(ipenal2, x2, x4, t, num_object);

        vol *= (1. - lbd21 * lbd23 * lbd24);

        return;

      }
      else if (ipenal3 == 1) {
        cs_real_t lbd31 = 1.-_imm_lgth_cutcell(ipenal3, x3, x1, t, num_object);
        cs_real_t lbd32 = 1.-_imm_lgth_cutcell(ipenal3, x3, x2, t, num_object);
        cs_real_t lbd34 = 1.-_imm_lgth_cutcell(ipenal3, x3, x4, t, num_object);

        vol *= (1. - lbd31 * lbd32 * lbd34);

        return;

      }
      else if (ipenal4 == 1) {
        cs_real_t lbd41 = 1.-_imm_lgth_cutcell(ipenal4, x4, x1, t, num_object);
        cs_real_t lbd42 = 1.-_imm_lgth_cutcell(ipenal4, x4, x2, t, num_object);
        cs_real_t lbd43 = 1.-_imm_lgth_cutcell(ipenal4, x4, x3, t, num_object);

        vol *= (1. - lbd41 * lbd42 * lbd43);

        return;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    }
    else if (cpt4 == 3) {

      /* Three penalized points : tetrahedron from the non penalized
       * fourth point */
      vol = _tetra_vol(x1, x2, x3, x4);

      if (ipenal1 == 0) {
        cs_real_t lbd12 = _imm_lgth_cutcell(ipenal1, x1, x2, t, num_object);
        cs_real_t lbd13 = _imm_lgth_cutcell(ipenal1, x1, x3, t, num_object);
        cs_real_t lbd14 = _imm_lgth_cutcell(ipenal1, x1, x4, t, num_object);

        vol *= lbd12 * lbd13 * lbd14;

        return;

      }
      else if (ipenal2 == 0) {
        cs_real_t lbd21 = _imm_lgth_cutcell(ipenal2, x2, x1, t, num_object);
        cs_real_t lbd23 = _imm_lgth_cutcell(ipenal2, x2, x3, t, num_object);
        cs_real_t lbd24 = _imm_lgth_cutcell(ipenal2, x2, x4, t, num_object);

        vol *= lbd21 * lbd23 * lbd24;

        return;

      }
      else if (ipenal3 == 0) {
        cs_real_t lbd31 = _imm_lgth_cutcell(ipenal3, x3, x1, t, num_object);
        cs_real_t lbd32 = _imm_lgth_cutcell(ipenal3, x3, x2, t, num_object);
        cs_real_t lbd34 = _imm_lgth_cutcell(ipenal3, x3, x4, t, num_object);

        vol *= lbd31 * lbd32 * lbd34;

        return;

      }
      else if (ipenal4 == 0) {
        cs_real_t lbd41 = _imm_lgth_cutcell(ipenal4, x4, x1, t, num_object);
        cs_real_t lbd42 = _imm_lgth_cutcell(ipenal4, x4, x2, t, num_object);
        cs_real_t lbd43 = _imm_lgth_cutcell(ipenal4, x4, x3, t, num_object);

        vol *= lbd41 * lbd42 * lbd43;

        return;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Error in function _tetra_vol_cutcell\n");

    }
    else if (cpt4 == 2) {

      /* Two penalized points : more complex */
      /* One puts in y1 et y2 the two non penalized points and in y3 and y4
       * the penalized points */
      cs_real_3_t y1, y2, y3, y4;
      int ipenal1_init = 0;
      int ipenal2_init = 0;

      if (ipenal1 == 0 && ipenal2 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x2[i];
          y3[i] = x3[i];
          y4[i] = x4[i];
        }
      } else if (ipenal1 == 0 && ipenal3 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x3[i];
          y3[i] = x2[i];
          y4[i] = x4[i];
        }
      } else if (ipenal1 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x1[i];
          y2[i] = x4[i];
          y3[i] = x2[i];
          y4[i] = x3[i];
        }
      } else if (ipenal2 == 0 && ipenal3 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x3[i];
          y3[i] = x1[i];
          y4[i] = x4[i];
        }
      } else if (ipenal2 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x2[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x3[i];
        }
      } else if (ipenal3 == 0 && ipenal4 == 0) {
        for (int i = 0; i < 3; i++) {
          y1[i] = x3[i];
          y2[i] = x4[i];
          y3[i] = x1[i];
          y4[i] = x2[i];
        }
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

      cs_real_t lbd13 = _imm_lgth_cutcell(ipenal1_init, y1, y3, t, num_object);
      cs_real_t lbd14 = _imm_lgth_cutcell(ipenal1_init, y1, y4, t, num_object);
      cs_real_t lbd23 = _imm_lgth_cutcell(ipenal2_init, y2, y3, t, num_object);
      cs_real_t lbd24 = _imm_lgth_cutcell(ipenal2_init, y2, y4, t, num_object);

      cs_real_3_t y13, y14, y23, y24;
      for (int i = 0; i < 3; i++) {
        y13[i] = y1[i] + lbd13 * d13[i];
        y14[i] = y1[i] + lbd14 * d14[i];
        y23[i] = y2[i] + lbd23 * d23[i];
        y24[i] = y2[i] + lbd24 * d24[i];
      }

      vol = _prism_vol(y1, y2, y13, y23, y24, y14);
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
 * \param[in]   comp_cell          list of cells to recompute porosity
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_cut_porosity(const cs_mesh_t *mesh,
                           const int  *comp_cell)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_t *cell_vol  = mq_g->cell_vol;

  const cs_real_3_t *cell_cen = mq_g->cell_cen;
  const cs_real_3_t *i_face_cog = mq_g->i_face_cog;
  const cs_real_3_t *b_face_cog = mq_g->b_face_cog;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  cs_real_t t_cur = cs_glob_time_step->t_cur;
  int icut = cs_ibm->nb_cut_cells;
  int cut_edges = cs_ibm->nb_cut_edges*4;

  cs_real_t den_cut_edges = 1. / (cs_real_t)(cut_edges);

  cs_array<cs_real_t> lambda(cut_edges + 1), lambda0(cut_edges + 1);

  for (int iraf = 0; iraf <= cut_edges; iraf++)
    lambda0[iraf] = iraf * den_cut_edges;

  lambda[0] = lambda0[0];
  lambda[1] = lambda0[cut_edges];
  lambda[2] = lambda0[cut_edges/2];

  for (int k = 3; k < cut_edges/2+2; k++)
    lambda[k] = lambda0[k-2];

  for (int k = cut_edges/2+2; k <= cut_edges; k++)
    lambda[k] = lambda0[k-1];

  cs_real_t voltot = 0;

  /* Array poro_compute :
     -> 0 completely fluid cell
     -> 1 completely solid cell
     -> 2 mixed cell (porosity to compute) */

  cs_array<int> poro_compute(n_cells_ext);
  poro_compute.set_to_val(-1);

  /* Cell centers */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    int num_object = -1;
    poro_compute[c_id] = _penal_glob(c_id, cell_cen[c_id], t_cur, num_object);
  }

  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, poro_compute);

  cs_real_t *poro_val = CS_F_(poro)->val;

  cs_lnum_t v_id1, v_id2;
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (comp_cell[c_id0] + comp_cell[c_id1] > 0) {

      int num_objecti = -1;
      int num_objectj = -1;
      int num_objectij = -1;

      /* Go through the edges */
      for (cs_lnum_t vtx_id = i_face_vtx_idx[f_id];
          vtx_id < i_face_vtx_idx[f_id + 1]; vtx_id++) {
        v_id1 = i_face_vtx[vtx_id];

        if (vtx_id < mesh->i_face_vtx_idx[f_id + 1] - 1)
          v_id2 = mesh->i_face_vtx_lst[vtx_id + 1];
        else
          v_id2 = mesh->i_face_vtx_lst[mesh->i_face_vtx_idx[f_id]];

        cs_real_3_t x1, x2, x12;
        for (int i = 0; i < 3; i++) {
          x1[i] = mesh->vtx_coord[3*v_id1 + i];
          x2[i] = mesh->vtx_coord[3*v_id2 + i];
        }

        for (int iraf = 0; iraf < cut_edges; iraf++) {
          for (int i = 0; i < 3; i++)
            x12[i] = lambda[iraf] * x1[i] + (1.-lambda[iraf]) * x2[i];

          int in12 = _penal_glob(c_id0, x12, t_cur, num_objectij);

          if (poro_compute[c_id0] < 2 && in12 != poro_compute[c_id0])
            poro_compute[c_id0] = 2;

          if (poro_compute[c_id1] < 2 && in12 != poro_compute[c_id1])
            poro_compute[c_id1] = 2;

          if (poro_compute[c_id0] > 1 && poro_compute[c_id1] > 1)
            break;
        }

        if (poro_compute[c_id0] > 1 && poro_compute[c_id1] > 1)
          break;
      }

      /* Face center */
      cs_real_3_t xif, xjf;
      for (int idim = 0; idim < 3; idim++) {
        xif[idim] = 0.5 * (i_face_cog[f_id][idim] + cell_cen[c_id0][idim]);
        xjf[idim] = 0.5 * (i_face_cog[f_id][idim] + cell_cen[c_id1][idim]);
      }

      if (poro_compute[c_id0] < 2) {
        int in12 = _penal_glob(c_id0, xif, t_cur, num_objecti);
        if (in12 != poro_compute[c_id0])
          poro_compute[c_id0] = 2;
      }

      if (poro_compute[c_id1] < 2) {
        int in12 = _penal_glob(c_id1, xjf, t_cur, num_objectj);
        if (in12 != poro_compute[c_id1])
          poro_compute[c_id1] = 2;
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    if (comp_cell[c_id] > 0  && poro_compute[c_id] < 2) {

      int num_objecti = -1;

      /* Go through the edges */
      for (cs_lnum_t vtx_id = b_face_vtx_idx[f_id];
          vtx_id < b_face_vtx_idx[f_id + 1]; vtx_id++) {
        v_id1 = b_face_vtx[vtx_id];

        if (vtx_id < mesh->b_face_vtx_idx[f_id + 1] - 1)
          v_id2 = mesh->b_face_vtx_lst[vtx_id + 1];
        else
          v_id2 = mesh->b_face_vtx_lst[mesh->b_face_vtx_idx[f_id]];

        cs_real_3_t x1, x2, x12;
        for (int i = 0; i < 3; i++) {
          x1[i] = mesh->vtx_coord[3*v_id1 + i];
          x2[i] = mesh->vtx_coord[3*v_id2 + i];
        }

        for (int iraf = 0; iraf < cut_edges; iraf++) {
          for (int i = 0; i < 3; i++)
            x12[i] = lambda[iraf] * x1[i] + (1.-lambda[iraf]) * x2[i];

          if (poro_compute[c_id] < 2) {
            int in12 = _penal_glob(c_id, x12, t_cur, num_objecti);
            if (in12 != poro_compute[c_id])
              poro_compute[c_id] = 2;
          }

          if (poro_compute[c_id] > 1)
            break;
        }

        if (poro_compute[c_id] > 1)
          break;
      }

      /* Face center */
      if (poro_compute[c_id] < 2) {

        cs_real_3_t xif;
        for (int idim = 0; idim < 3; idim++)
          xif[idim] = 0.5 * (b_face_cog[f_id][idim] + cell_cen[c_id][idim]);

        if (poro_compute[c_id] < 2) {
          int in12 = _penal_glob(c_id, xif, t_cur, num_objecti);
          if (in12 != poro_compute[c_id])
            poro_compute[c_id] = 2;
        }
      }
    }
  }

  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, poro_compute);

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    if (comp_cell[c_id] > 0)
      poro_val[c_id] = 0.;

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    if (comp_cell[c_id] > 0)
      if (poro_compute[c_id] == 0)
        poro_val[c_id] = 1.;

  cs_real_3_t x1, x2, x12;

  int nraf = 5;
  cs_array_2d<cs_real_t> lbd_tetporvol(cs_ibm->nb_cut_edges, 2);
  cs_array_2d<cs_real_t> lbd_tetporvol2(nraf, 2);

  cs_real_t den_cut_edge = 1. / (cs_real_t)(cs_ibm->nb_cut_edges);
  for (int k = 1; k < cs_ibm->nb_cut_edges; k++) {
    lbd_tetporvol(k, 0) = k*den_cut_edge;
    lbd_tetporvol(k, 1) = 1. - lbd_tetporvol(k, 0);
  }

  cs_real_t den_raf = 1. / (cs_real_t)(nraf);
  for (int k = 1; k < nraf; k++) {
    lbd_tetporvol2(k, 0) = k*den_raf;
    lbd_tetporvol2(k, 1) = 1. - lbd_tetporvol2(k, 0);
  }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (poro_compute[c_id0] > 1 || poro_compute[c_id1] > 1) {
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

        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
          x12[i] = 0.5*(x1[i]+x2[i]);
        }

        if (poro_compute[c_id0] > 1) {
          _tetra_vol_cutcell(c_id0, voltot, x1, x12, xf, xi, t_cur,
                             icut, nraf, lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
          poro_val[c_id0] += voltot;

          _tetra_vol_cutcell(c_id0, voltot, x12, x2, xf, xi, t_cur,
                             icut, nraf, lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
          poro_val[c_id0] += voltot;
        }

        if (poro_compute[c_id1] > 1) {
          _tetra_vol_cutcell(c_id1, voltot, x1, x12, xf, xj, t_cur,
                             icut, nraf, lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objectj);
          poro_val[c_id1] += voltot;

          _tetra_vol_cutcell(c_id1, voltot, x12, x2, xf, xj, t_cur,
                             icut, nraf, lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objectj);
          poro_val[c_id1] += voltot;
        }
      }

      cs_lnum_t vtx_id1
        = i_face_vtx[i_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2
        = i_face_vtx[i_face_vtx_idx[f_id]];

      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
        x12[i] = 0.5*(x1[i]+x2[i]);
      }

      if (poro_compute[c_id0] > 1) {
        _tetra_vol_cutcell(c_id0, voltot, x1, x12, xf, xi, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
        poro_val[c_id0] += voltot;

        _tetra_vol_cutcell(c_id0, voltot, x12, x2, xf, xi, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
        poro_val[c_id0] += voltot;
      }

      if (poro_compute[c_id1] > 1) {
        _tetra_vol_cutcell(c_id1, voltot, x1, x12, xf, xj, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objectj);
        poro_val[c_id1] += voltot;

        _tetra_vol_cutcell(c_id1, voltot, x12, x2, xf, xj, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objectj);
        poro_val[c_id1] += voltot;
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    if (poro_compute[c_id] > 1) {
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

        for (int i = 0; i < 3; i++) {
          x1[i] = vtx_crd[vtx_id1][i];
          x2[i] = vtx_crd[vtx_id2][i];
          x12[i] = 0.5*(x1[i] + x2[i]);
        }

        _tetra_vol_cutcell(c_id, voltot, x1, x12, xf, xi, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
        poro_val[c_id] += voltot;

        _tetra_vol_cutcell(c_id, voltot, x12, x2, xf, xi, t_cur, icut, nraf,
                           lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
        poro_val[c_id] += voltot;
      }

      cs_lnum_t vtx_id1
        = b_face_vtx[b_face_vtx_idx[f_id + 1] - 1];
      cs_lnum_t vtx_id2
        = b_face_vtx[b_face_vtx_idx[f_id]];

      for (int i = 0; i < 3; i++) {
        x1[i] = vtx_crd[vtx_id1][i];
        x2[i] = vtx_crd[vtx_id2][i];
        x12[i] = 0.5*(x1[i] + x2[i]);
      }

      _tetra_vol_cutcell(c_id, voltot, x1, x12, xf, xi, t_cur, icut, nraf,
                         lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
      poro_val[c_id] += voltot;

      _tetra_vol_cutcell(c_id, voltot, x12, x2, xf, xi, t_cur, icut, nraf,
                         lbd_tetporvol.view(), lbd_tetporvol2.view(), num_objecti);
      poro_val[c_id] += voltot;
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (comp_cell[c_id] > 0) {
      if (poro_compute[c_id] > 1) {
        poro_val[c_id] /= cell_vol[c_id];
        poro_val[c_id] = cs::clamp(poro_val[c_id], 0., 1.);

        if (poro_val[c_id] < 1.e-5)
          poro_val[c_id] = 0.;
      }
    }
    else {
      poro_val[c_id] = cs::clamp(poro_val[c_id], 0., 1.);
    }
  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, poro_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell cog (and cell porosity from porosity at vertices)
 *         Cut method into sub-tetras
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   v_poro             vertex porosity
 * \param[in]   por_init           initialization porosity
 * \param[in]   comp_cell          list of cells to recompute porosity
 * \param[in]   icycle
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_cog(const cs_mesh_t      *mesh,
                  cs_mesh_quantities_t *mesh_quantities,
                  const cs_real_t      *v_poro,
                  const cs_real_t      *por_init,
                  const int            *comp_cell,
                  const int             icycle)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_t *cell_vol = mq_g->cell_vol;

  const cs_real_3_t *cell_cen = mq_g->cell_cen;
  cs_real_3_t *cell_f_cen = mesh_quantities->cell_cen;
  cs_real_3_t *cell_s_cen
    = (cs_real_3_t *)(cs_field_by_name("cell_s_cen")->val);
  const cs_real_3_t *i_face_cog = mq_g->i_face_cog;
  const cs_real_3_t *b_face_cog = mq_g->b_face_cog;

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
  if (icycle > 0)
    comp_all_cog = false;

  int icut = cs_ibm->nb_cut_cells;

  cs_array<cs_real_t> c_poro(n_cells_ext);

  cs_vertex_to_cell<1>(CS_VERTEX_TO_CELL_SHEPARD, 0, nullptr,
                       v_poro, c_poro);

  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, c_poro);

  cs_real_t voltot = 0.;

  cs_array<cs_real_t> porbis(n_cells_ext);
  porbis.zero();

  cs_array<bool> compute_cell(n_cells_ext);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    compute_cell[c_id] = true;

  if (!comp_all_cog)
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
      if (comp_cell[c_id] < 1)
        compute_cell[c_id] = false;

  cs_real_t *poro_val = CS_F_(poro)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    if (compute_cell[c_id]) {
      for (int i = 0; i < 3; i++) {
        cell_f_cen[c_id][i] = 0.;
        cell_s_cen[c_id][i] = 0.;
      }
    }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    bool computi = compute_cell[c_id0];
    bool computj = compute_cell[c_id1];

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
          _tetra_vol_poro(c_id0, voltot, cell_f_cen[c_id0], x1, por1,
                          x2, por2, xf, porf, xi, pori, icut);
          porbis[c_id0] += voltot;
          _tetra_vol_poro(c_id0, voltot, cell_s_cen[c_id0], x1, 1.-por1,
                          x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);
        }

        if (computj) {
          _tetra_vol_poro(c_id1, voltot, cell_f_cen[c_id1], x1, por1,
                          x2, por2, xf, porf, xj, porj, icut);
          porbis[c_id1] += voltot;
          _tetra_vol_poro(c_id1, voltot, cell_s_cen[c_id1], x1, 1.-por1,
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
        _tetra_vol_poro(c_id0, voltot, cell_f_cen[c_id0], x1, por1,
                        x2, por2, xf, porf, xi, pori, icut);
        porbis[c_id0] += voltot;
        _tetra_vol_poro(c_id0, voltot, cell_s_cen[c_id0], x1, 1.-por1,
                        x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);
      }

      if (computj) {
        _tetra_vol_poro(c_id1, voltot, cell_f_cen[c_id1], x1, por1,
                        x2, por2, xf, porf, xj, porj, icut);
        porbis[c_id1] += voltot;
        _tetra_vol_poro(c_id1, voltot, cell_s_cen[c_id1], x1, 1.-por1,
                        x2, 1.-por2, xf, 1.-porf, xj, 1.-porj, icut);
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    bool computi = compute_cell[c_id];

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

        _tetra_vol_poro(c_id, voltot, cell_f_cen[c_id], x1, por1,
                        x2, por2, xf, porf, xi, pori, icut);
        porbis[c_id] += voltot;
        _tetra_vol_poro(c_id, voltot, cell_s_cen[c_id], x1, 1.-por1,
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

      _tetra_vol_poro(c_id, voltot, cell_f_cen[c_id], x1, por1,
                      x2, por2, xf, porf, xi, pori, icut);
      porbis[c_id] += voltot;
      _tetra_vol_poro(c_id, voltot, cell_s_cen[c_id], x1, 1.-por1,
                      x2, 1.-por2, xf, 1.-porf, xi, 1.-pori, icut);

    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (compute_cell[c_id]) {
      for (int i = 0; i < 3; i++)
        cell_f_cen[c_id][i] /= cs::max(porbis[c_id],
                                       cs_math_epzero * cell_vol[c_id]);
      for (int i = 0; i < 3; i++)
        cell_s_cen[c_id][i] /= cs::max(cell_vol[c_id]-porbis[c_id],
                                       cs_math_epzero * cell_vol[c_id]);

      porbis[c_id] /= cell_vol[c_id];
      porbis[c_id] = cs::clamp(porbis[c_id], 0., 1.);

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

  /* Recompute cell_f_cen for cells with almost null porosities */
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                           (cs_real_t *)cell_f_cen, 3);
  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, porbis);

  if (mesh->n_init_perio > 0)
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                              (cs_real_t *)cell_f_cen);

  cs_lnum_t size_weight = n_cells_ext;
  if (mesh->n_vertices > size_weight)
    size_weight = mesh->n_vertices;

  cs_array<cs_real_t> weight(size_weight);
  weight.set_to_val(1.);

  for (cs_lnum_t c_id = 0.; c_id < n_cells_ext; c_id++)
    if (compute_cell[c_id]) {
      weight[c_id] = 0.;

      if (porbis[c_id] < 1.e-5)
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id][idir] = 0.;
    }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (compute_cell[c_id0])
      if (porbis[c_id0] < 1.e-5 && porbis[c_id1] >= 1.e-5) {
        weight[c_id0] += porbis[c_id1];
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id0][idir] += porbis[c_id1] * i_face_cog[f_id][idir];
      }

    if (compute_cell[c_id1])
      if (porbis[c_id1] < 1.e-5 && porbis[c_id0] >= 1.e-5) {
        weight[c_id1] += porbis[c_id0];
        for (int idir = 0; idir < 3; idir++)
          cell_f_cen[c_id1][idir] += porbis[c_id0] * i_face_cog[f_id][idir];
      }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (compute_cell[c_id])
      if (porbis[c_id] < 1.e-5) {
        if (weight[c_id] > 0.01)
          for (int idir = 0; idir < 3; idir++)
            cell_f_cen[c_id][idir] /= weight[c_id];
        else
          for (int idir = 0; idir < 3; idir++)
            cell_f_cen[c_id][idir] = cell_cen[c_id][idir];
      }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (compute_cell[c_id])
      if (por_init[c_id] < 1.e-5)
        for (int i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = cell_cen[c_id][i];

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cs_ibm->prob_dim == CS_IBM_2D_X) {
      cell_f_cen[c_id][0] = cell_cen[c_id][0];
      cell_s_cen[c_id][0] = cell_cen[c_id][0];
    }
    else if (cs_ibm->prob_dim == CS_IBM_2D_Y) {
      cell_f_cen[c_id][1] = cell_cen[c_id][1];
      cell_s_cen[c_id][1] = cell_cen[c_id][1];
    }
    else if (cs_ibm->prob_dim == CS_IBM_2D_Z) {
      cell_f_cen[c_id][2] = cell_cen[c_id][2];
      cell_s_cen[c_id][2] = cell_cen[c_id][2];
    }

    if (poro_val[c_id] > 1.-1.e-5) {
      for (int i = 0; i < 3; i++) {
        cell_s_cen[c_id][i] = cell_cen[c_id][i];
        cell_f_cen[c_id][i] = cell_cen[c_id][i];
      }
    }
  }

  cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                           (cs_real_t *)cell_f_cen, 3);
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                           (cs_real_t *)cell_s_cen, 3);

  if (mesh->n_init_perio > 0) {
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                              (cs_real_t *)cell_f_cen);
    cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                              (cs_real_t *)cell_s_cen);
  }

  /* -> porosity = porbis + volume conservation */
  if (cs_ibm->porosity_from_nodes) {
    cs_real_t volpor = 0.;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0)
        volpor += cell_vol[c_id] * poro_val[c_id];

    cs::parall::sum(volpor);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0)
        poro_val[c_id] = porbis[c_id];

    if (cs_ibm->ensure_isovol)
      for (int iter = 1; iter < 4; iter++) {
        cs_real_t aa = 0.;
        cs_real_t bb = 0.;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          if (comp_cell[c_id] > 0) {
            cs_real_t porloc =  poro_val[c_id];
            cs_real_t pump = cs::abs(porloc * (1. - porloc));
            aa += pump * cell_vol[c_id];
            bb += porloc * cell_vol[c_id];
          }

        cs::parall::sum(aa, bb);

        cs_real_t beta = (volpor - bb) / cs::max(aa, 1.e-20);

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          if (comp_cell[c_id] > 0) {
            cs_real_t porloc = poro_val[c_id];
            cs_real_t pump = cs::abs(porloc * (1. - porloc));
            poro_val[c_id] += beta * pump;
            poro_val[c_id] = cs::clamp(poro_val[c_id], 0., 1.);
          }
      }

    cs_halo_sync(mesh->halo, CS_HALO_STANDARD, poro_val);
  }
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
_compute_b_fac_porosity(const cs_mesh_t      *mesh,
                        cs_mesh_quantities_t *mesh_quantities,
                        const cs_real_t      *c_poro,
                        const cs_real_t      *v_poro,
                        cs_real_t            *bfpro_poro)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;

  const cs_real_t *b_face_surf  = mq_g->b_face_surf;
  const cs_real_3_t *b_face_cog = mq_g->b_face_cog;
  cs_real_3_t *b_f_face_cog = mesh_quantities->b_face_cog;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mq_g->b_face_normal;
  const cs_nreal_3_t *b_face_u_normal = mq_g->b_face_u_normal;
  const cs_real_3_t *cell_f_cen = mesh_quantities->cell_cen;
  const cs_real_t *b_dist = mq_g->b_dist;
  cs_real_t *b_f_dist = mesh_quantities->b_dist;
  cs_real_3_t *b_f_face_normal = (cs_real_3_t *)mesh_quantities->b_face_normal;
  cs_nreal_3_t *b_f_face_u_normal = mesh_quantities->b_face_u_normal;
  cs_real_t *b_f_face_surf = mesh_quantities->b_face_surf;
  int *c_disable_flag = mesh_quantities->c_disable_flag;

  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  int icut = cs_ibm->nb_cut_cells;

  const cs_real_t *poro_val_pre = CS_F_(poro)->val_pre;

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    bfpro_poro[f_id] = 0.5*(c_poro[c_id] + poro_val_pre[c_id]);

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

    if (ptin == 0)
      bfpro_poro[f_id] = 1.;

    else if (ptin == cpt)
      bfpro_poro[f_id] = 0.;

    else if (ptin < cpt) {
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
      bfpro_poro[f_id] = cs::min(porb, c_poro[c_id]);
    }

    /* 2D cases treatment */
    if (cs_ibm->prob_dim != CS_IBM_3D) {
      cs_real_3_t rn;
      for (int i = 0; i < 3; i++)
        rn[i] = b_face_u_normal[f_id][i];

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

  /* Only for cut-cell (cf. internal face treatment) */
  if (cs_ibm->algo_choice == CS_IBM_ALGO_CUT_CELLS) {
    cs_lnum_t node[100];
    cs_real_3_t xn[100];
    cs_real_3_t xp[100];
    int penal[100];

    cs_real_3_t cog;
    cs_real_3_t cogf;
    cs_real_t surf;
    cs_real_t t_cur = cs_glob_time_step->t_cur;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t c_id = b_face_cells[f_id];

      int num_object = -1;
      cs_real_t rnx = cs::abs(b_face_u_normal[f_id][0]);
      cs_real_t rny = cs::abs(b_face_u_normal[f_id][1]);
      cs_real_t rnz = cs::abs(b_face_u_normal[f_id][2]);

      /* Cell I totally fluid -> cog = b_face_cog */
      if (c_poro[c_id] >= 0.9999) {
        bfpro_poro[f_id] = 1.;
        for (int k = 0; k < 3; k++)
          b_f_face_cog[f_id][k] = b_face_cog[f_id][k];

      }  else if (cs_ibm->prob_dim == CS_IBM_2D_X && rnx > rny + rnz) {
        bfpro_poro[f_id] = c_poro[c_id];
        b_f_face_cog[f_id][0] = b_face_cog[f_id][0];
        b_f_face_cog[f_id][1] = cell_f_cen[c_id][1];
        b_f_face_cog[f_id][2] = cell_f_cen[c_id][2];


      }  else if (cs_ibm->prob_dim == CS_IBM_2D_Y && rny > rnx + rnz) {
        bfpro_poro[f_id] = c_poro[c_id];
        b_f_face_cog[f_id][0] = cell_f_cen[c_id][0];
        b_f_face_cog[f_id][1] = b_face_cog[f_id][1];
        b_f_face_cog[f_id][2] = cell_f_cen[c_id][2];

      }  else if (cs_ibm->prob_dim == CS_IBM_2D_Z && rnz > rnx + rny) {
        bfpro_poro[f_id] = c_poro[c_id];
        b_f_face_cog[f_id][0] = cell_f_cen[c_id][0];
        b_f_face_cog[f_id][1] = cell_f_cen[c_id][1];
        b_f_face_cog[f_id][2] = b_face_cog[f_id][2];

        /* For safety */
      } else if (c_poro[c_id] < 1.e-10) {
        bfpro_poro[f_id] = 0.;
        for (int k = 0; k < 3; k++)
          b_f_face_cog[f_id][k] = b_face_cog[f_id][k];

      /* Otherwise loop on the edges */
      }  else  {

        int nd = -1;
        int sumt = mesh->b_face_vtx_idx[f_id + 1] - mesh->b_face_vtx_idx[f_id];
        int sump = 0;

        for (cs_lnum_t v_id = mesh->b_face_vtx_idx[f_id]; v_id < mesh->b_face_vtx_idx[f_id + 1]; v_id++) {
          nd++;
          if (nd >= 100)
            bft_error(__FILE__, __LINE__, 0, "Arrays declared above should be resized above 100\n");

          node[nd] = mesh->b_face_vtx_lst[v_id];
          for (int k = 0; k < 3; k++)
            xn[nd][k] = mesh->vtx_coord[3*node[nd] + k];

          int in = _penal_glob(c_id, xn[nd], t_cur, num_object);
          penal[nd] = in;

          if (in > 0)
            sump++;
        }

        /* If all points are solid */
        if (sump == sumt) {
          cs_real_t lambda = cs_math_3_distance_dot_product(cell_f_cen[c_id], b_face_cog[f_id], b_face_u_normal[f_id]);

          bfpro_poro[f_id] = 0;
          for (int k = 0; k < 3; k++)
            b_f_face_cog[f_id][k] = cell_f_cen[c_id][k]
                                  + lambda * b_face_u_normal[f_id][k];

        /* If all points are fluid */
        } else if (sump == 0) {
          bfpro_poro[f_id] = 1;
          for (int k = 0; k < 3; k++)
            b_f_face_cog[f_id][k] = b_face_cog[f_id][k];

        }
        else {
        /* Compute cog and bfpro_poro */

          int cpt = -1;
          /* Loop on the edges */
          for (int ar = 0; ar < sumt; ar++) {
            int pt1 = ar;
            int pt2 = ar + 1;
            if (ar == sumt - 1)
              pt2 = 0;

            /* First point is fluid, second point is solid */
            if (penal[pt1] == 0 && penal[pt2] == 1) {
              cpt++;
              cs_real_t lambda = _imm_lgth_cutcell(penal[pt1], xn[pt1], xn[pt2], t_cur, num_object);

              for (int k = 0; k < 3; k++)
                xp[cpt][k] = lambda * xn[pt2][k] + (1. - lambda) * (xn[pt1][k]);
            }

            /* First point is solid, second point is fluid */
            if (penal[pt1] == 1 && penal[pt2] == 0) {
              cpt++;
              cs_real_t lambda = _imm_lgth_cutcell(penal[pt1], xn[pt1], xn[pt2], t_cur, num_object);

              for (int k = 0; k < 3; k++)
                xp[cpt][k] = lambda * xn[pt1][k] + (1. - lambda) * (xn[pt2][k]);
            }

            /* Second point is fluid */
            if (penal[pt2] == 0) {
              cpt++;
              for (int k = 0; k < 3; k++)
                xp[cpt][k] = xn[pt2][k];
            }

            if (cpt >= 100 || pt2 >= 100)
              bft_error(__FILE__, __LINE__, 0, "Arrays declared above should be resized above 100\n");
          }

          cpt++;
          for (int k = 0; k < 3; k++)
            xp[cpt][k] = xp[0][k];

          /* Center of gravity */
          for (int k = 0; k < 3; k++)
            cog[k] = 0;

          for (int ipt = 0; ipt < cpt; ipt++)
            for (int k = 0; k < 3; k++)
              cog[k] += xp[ipt][k];

          for (int k = 0; k < 3; k++)
            cog[k] /= (cs_real_t)cpt;

          /* Loop on the triangles */
          for (int k = 0; k < 3; k++)
            cogf[k] = 0.;

          surf = 0.;

          /* Increment surface and cogfac */
          for (int tr = 0; tr < cpt; tr++) {
            cs_real_3_t cogloc;
            cs_real_3_t vectloc;

            for (int k = 0; k < 3; k++)
              cogloc[k] = (xp[tr][k] + xp[tr+1][k] + cog[k])*cs_math_1ov3;

            vectloc[0] = (xp[tr][1]-cog[1]) * (xp[tr+1][2]-cog[2])
                       - (xp[tr][2]-cog[2]) * (xp[tr+1][1]-cog[1]);
            vectloc[1] = (xp[tr][2]-cog[2]) * (xp[tr+1][0]-cog[0])
                       - (xp[tr][0]-cog[0]) * (xp[tr+1][2]-cog[2]);
            vectloc[2] = (xp[tr][0]-cog[0]) * (xp[tr+1][1]-cog[1])
                       - (xp[tr][1]-cog[1]) * (xp[tr+1][0]-cog[0]);

            cs_real_t surfloc = cs_math_3_norm(vectloc);
            surf += surfloc;

            for (int k = 0; k < 3; k++)
              cogf[k] += surfloc * cogloc[k];
          }

          /* Finalize cogfac */
          if (surf >= 2. * 1.e-6 * b_face_surf[f_id]) {
            for (int k = 0; k < 3; k++)
              cogf[k] /= surf;

            /* Final storage */
            for (int k = 0; k < 3; k++)
              b_f_face_cog[f_id][k] = cogf[k];

            cs_real_t poros_loc = 0.5*surf / b_face_surf[f_id];

            bfpro_poro[f_id] = poros_loc;
          }
          else {
            /* Standard initialization of  */
            for (int k = 0; k < 3; k++)
              b_f_face_cog[f_id][k] = b_face_cog[f_id][k];
            bfpro_poro[f_id] = 0.;
          }
        }
      }
    }
  }
  else {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t c_id = b_face_cells[f_id];

      cs_real_t rnx = cs::abs(b_face_u_normal[f_id][0]);
      cs_real_t rny = cs::abs(b_face_u_normal[f_id][1]);
      cs_real_t rnz = cs::abs(b_face_u_normal[f_id][2]);

      /* Cell I totally fluid -> cog = b_face_cog */
      if (c_poro[c_id] >= 0.9999) {
        for (int k = 0; k < 3; k++)
          b_f_face_cog[f_id][k] = b_face_cog[f_id][k];

      }  else if (cs_ibm->prob_dim == CS_IBM_2D_X && rnx > rny + rnz) {
        b_f_face_cog[f_id][0] = b_face_cog[f_id][0];
        b_f_face_cog[f_id][1] = cell_f_cen[c_id][1];
        b_f_face_cog[f_id][2] = cell_f_cen[c_id][2];

      }  else if (cs_ibm->prob_dim == CS_IBM_2D_Y && rny > rnx + rnz) {
        b_f_face_cog[f_id][0] = cell_f_cen[c_id][0];
        b_f_face_cog[f_id][1] = b_face_cog[f_id][1];
        b_f_face_cog[f_id][2] = cell_f_cen[c_id][2];

      }  else if (cs_ibm->prob_dim == CS_IBM_2D_Z && rnz > rnx + rny) {
        b_f_face_cog[f_id][0] = cell_f_cen[c_id][0];
        b_f_face_cog[f_id][1] = cell_f_cen[c_id][1];
        b_f_face_cog[f_id][2] = b_face_cog[f_id][2];

      /* Ad hoc interpolation otherwise */
      }
      else {
        cs_real_t num = 0.;
        cs_real_t den = 0.;
        cs_real_3_t xip, ipf;
        for (int k = 0; k < 3; k++) {
          xip[k] = cell_f_cen[c_id][k];
          ipf[k] = b_face_cog[f_id][k] - xip[k];
          num += ipf[k] * b_face_normal[f_id][k];
          den += cs_math_pow2(b_face_normal[f_id][k]);
        }

        if (den > 1.e-20) {
          cs_real_t lambda = num / den;
          for (int k = 0; k < 3; k++)
            b_f_face_cog[f_id][k] = xip[k] + lambda * b_face_normal[f_id][k];

        }
        else {
          for (int k = 0; k < 3; k++)
            b_f_face_cog[f_id][k] = b_face_cog[f_id][k];
        }
      }
    }
  }

  /* Compute b_f_dist and face normal and surface*/
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    //if (comp_fac[f_id]) {

      cs_real_t dist_ipf
        = cs::abs(cs_math_3_distance_dot_product(cell_f_cen[c_id],
                                                 b_f_face_cog[f_id],
                                                 b_face_u_normal[f_id]));

      dist_ipf = cs::max(dist_ipf, 0.5 * b_dist[f_id]);

      b_f_dist[f_id] = dist_ipf;
    //}

    if (c_disable_flag[c_id] == 1 || bfpro_poro[f_id] < cs_ibm->min_poro) {
      b_f_face_normal[f_id][0] = 0.;
      b_f_face_normal[f_id][1] = 0.;
      b_f_face_normal[f_id][2] = 0.;
      b_f_face_u_normal[f_id][0] = 0.;
      b_f_face_u_normal[f_id][1] = 0.;
      b_f_face_u_normal[f_id][2] = 0.;
      b_f_face_surf[f_id] = 0.;
    }
    else {
      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_normal[f_id][i] = bfpro_poro[f_id] * b_face_normal[f_id][i];

      b_f_face_surf[f_id] = cs_math_3_norm(b_f_face_normal[f_id]);
    }

    if (mesh_quantities->b_f_face_factor != nullptr) {
      //FIXME
      //if (face_porosity > cs_math_epzero) {
      //  mesh_quantities->b_f_face_factor[f_id] = c_poro[c_id] / bfpro_poro[f_id];
      //}
      //else {
      mesh_quantities->b_f_face_factor[f_id] = 1.;
      //}
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
_compute_i_fac_porosity(const cs_mesh_t      *mesh,
                        cs_mesh_quantities_t *mesh_quantities,
                        const cs_real_t      *c_poro,
                        cs_real_t            *ifpro_poro)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_3_t *i_face_cog = mq_g->i_face_cog;
  const cs_real_3_t *cell_cen = mq_g->cell_cen;
  const cs_real_3_t *i_face_normal = (cs_real_3_t *)mq_g->i_face_normal;
  cs_real_3_t *i_f_face_normal = (cs_real_3_t *)mesh_quantities->i_face_normal;
  cs_nreal_3_t *i_f_face_u_normal = mesh_quantities->i_face_u_normal;
  cs_real_t *i_f_face_surf = mesh_quantities->i_face_surf;
  int *c_disable_flag = mesh_quantities->c_disable_flag;

  const cs_real_t *poro_val_pre = CS_F_(poro)->val_pre;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t pori = 0.5 * (c_poro[c_id0] + poro_val_pre[c_id0]);
    cs_real_t porj = 0.5 * (c_poro[c_id1] + poro_val_pre[c_id1]);

    cs_real_t sizei = cs_math_3_distance(i_face_cog[f_id], cell_cen[c_id0]);
    cs_real_t sizej = cs_math_3_distance(i_face_cog[f_id], cell_cen[c_id1]);

    cs_real_t porij = _geom_face_fraction(pori, porj, sizei, sizej);

    cs_real_t porimin = cs::max(c_poro[c_id0],
                                poro_val_pre[c_id0]);
    cs_real_t porjmin = cs::max(c_poro[c_id1],
                                poro_val_pre[c_id1]);

    if (porij < 1.e-6 || porimin < 1.e-6 || porjmin < 1.e-6)
      porij = 0.;

    ifpro_poro[f_id] = porij;
  }

  /* If only one face porosity is positive for a cell, one cancels face
   * porosity */
  cs_array<cs_lnum_t> cpt(n_cells_ext);

  for (int ii = 0; ii < 4; ii++) {
    cpt.zero();

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      if (ifpro_poro[f_id] > 1.e-5) {
        cpt[c_id0]++;
        cpt[c_id1]++;
      }
    }

    cs_halo_sync(mesh->halo, CS_HALO_STANDARD, cpt);

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      if (cpt[c_id0] <= 1 || cpt[c_id1] <= 1)
        ifpro_poro[f_id] = 0;
    }
  }

  cpt.zero();

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (ifpro_poro[f_id] > 1.e-5) {
      cpt[c_id0]++;
      cpt[c_id1]++;
    }
  }

  cs_real_t *cvara_poro = CS_F_(poro)->val_pre;
  cs_real_t *cvar_poro  = CS_F_(poro)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (cpt[c_id] <= 1)
      cvar_poro[c_id] = 0.;

  cs_halo_sync(mesh->halo, CS_HALO_STANDARD, cvar_poro);

  if (cs_glob_time_step->nt_cur <= 1)
    cs_array_real_copy(n_cells_ext, cvar_poro, cvara_poro);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (   c_disable_flag[c_id0] == 1
        || c_disable_flag[c_id1] == 1 || ifpro_poro[f_id] < cs_ibm->min_poro) {
      i_f_face_normal[f_id][0] = 0.;
      i_f_face_normal[f_id][1] = 0.;
      i_f_face_normal[f_id][2] = 0.;
      i_f_face_u_normal[f_id][0] = 0.;
      i_f_face_u_normal[f_id][1] = 0.;
      i_f_face_u_normal[f_id][2] = 0.;
      i_f_face_surf[f_id] = 0.;
    }
    else {
      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[f_id][i] = ifpro_poro[f_id] * i_face_normal[f_id][i];

      i_f_face_surf[f_id] = cs_math_3_norm(i_f_face_normal[f_id]);
    }

    if (mesh_quantities->i_f_face_factor != nullptr) {
      //FIXME
      //if (porij > cs_math_epzero) {
      //  mesh_quantities->i_f_face_factor[f_id][0] = c_poro[c_id0] / porij;
      //  mesh_quantities->i_f_face_factor[f_id][1] = c_poro[c_id1] / porij;
      //}
      //else {
        mesh_quantities->i_f_face_factor[f_id][0] = 1.;
        mesh_quantities->i_f_face_factor[f_id][1] = 1.;
      //}
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute internal faces center of gravity
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 * \param[in]   c_poro             cell porosity
 * \param[in]   v_poro             vertex porosity
 */
/*----------------------------------------------------------------------------*/
static void
_compute_i_fac_cog(const cs_mesh_t      *mesh,
                   cs_mesh_quantities_t *mesh_quantities,
                   const cs_real_t      *c_poro,
                   const cs_real_t      *v_poro,
                   bool                 *comp_fac)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_t *i_face_surf = mq_g->i_face_surf;
  const cs_real_t *i_dist   = mq_g->i_dist;
  cs_real_t *i_f_dist = mesh_quantities->i_dist;
  const cs_real_3_t *i_face_cog = mq_g->i_face_cog;
  cs_real_3_t *i_f_face_cog = mesh_quantities->i_face_cog;
  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mq_g->i_face_normal;
  const cs_nreal_3_t *i_face_u_normal = mq_g->i_face_u_normal;
  cs_real_t *i_f_weight = mesh_quantities->weight;
  const cs_real_3_t *cell_f_cen = mesh_quantities->cell_cen;

  /* Compute weight and inv_dist */
  constexpr cs_real_t weight_min = 0.001;
  constexpr cs_real_t weight_max = 0.999;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (comp_fac[f_id]) {
      cs_real_t num = 0.;
      cs_real_t den = 0.;
      cs_real_3_t xip, xjp, ipf, ipjp;
      for (int k = 0; k < 3; k++) {
        xip[k] = cell_f_cen[c_id0][k];
        xjp[k] = cell_f_cen[c_id1][k];

        ipf[k] = i_face_cog[f_id][k] - xip[k];
        ipjp[k] = xjp[k] - xip[k];

        num += ipf[k] * i_face_normal[f_id][k];
        den += ipjp[k] * i_face_normal[f_id][k];
      }

      cs_real_t lambda = 0.5;
      if (cs::abs(den) > 1.e-20)
        lambda = num / den;

      lambda = cs::clamp(lambda, weight_min, weight_max);
      i_f_weight[f_id] = 1. - lambda;

      cs_real_t dist_ipjp = cs::abs(cs_math_3_distance_dot_product(xip, xjp, i_face_u_normal[f_id]));
      dist_ipjp = cs::max(dist_ipjp, 0.5 * i_dist[f_id]);

      i_f_dist[f_id] = dist_ipjp;
    }
  }

  cs_lnum_t node[100];
  cs_real_t porn[100];
  cs_real_3_t xn[100];
  cs_real_3_t xp[100];
  cs_real_3_t cog;
  cs_real_3_t cogf;
  cs_real_t surf;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    if (comp_fac[f_id]) {
      /* Cell I or J totally fluid -> cog = i_face_cog */
      if (c_poro[c_id0] > 0.99999 || c_poro[c_id1] > 0.99999) {
        for (int k = 0; k < 3; k++)
          i_f_face_cog[f_id][k] = i_face_cog[f_id][k];

      /* Loop on the edges */
      }  else {
        int nd = -1;
        int sumt = mesh->i_face_vtx_idx[f_id + 1] - mesh->i_face_vtx_idx[f_id];
        int sump = 0;

        for (cs_lnum_t v_id = mesh->i_face_vtx_idx[f_id]; v_id < mesh->i_face_vtx_idx[f_id + 1]; v_id++) {
          nd++;
          if (nd >= 100)
            bft_error(__FILE__, __LINE__, 0, "Arrays declared above should be resized above 100\n");

          node[nd] = mesh->i_face_vtx_lst[v_id];
          porn[nd] = v_poro[node[nd]];
          if (v_poro[node[nd]] >= 0.5)
            sump++;
        }
        /* If all points are fluid */
        if (sump == sumt) {
          for (int k = 0; k < 3; k++)
            i_f_face_cog[f_id][k] = i_face_cog[f_id][k];

        /* If all points are solid */
        } else if (sump == 0) {
          cs_real_t pnd = i_f_weight[f_id];
          for (int k = 0; k < 3; k++) {
            cs_real_t xip = cell_f_cen[c_id0][k];
            cs_real_t xjp = cell_f_cen[c_id1][k];

            cs_real_t xf = pnd * xip + (1. - pnd) * xjp;
            i_f_face_cog[f_id][k] = xf;
          }
        /* Compute the new internal face cog */
        }
        else {
          /* Store node coordinates */
          for (int nod = 0; nod < sumt; nod++)
            for (int k = 0; k < 3; k++)
              xn[nod][k] = mesh->vtx_coord[3*node[nod] + k];

          /* Initialize node counter */
          int cpt = -1;

          /* Loop on the edges */
          for (int ar = 0; ar < sumt; ar++) {

            int pt1 = ar;
            int pt2 = ar + 1;
            if (ar == sumt - 1)
              pt2 = 0;

            /* First point fluid, second point solid */
            if (porn[pt1] >= 0.5 && porn[pt2] < 0.5) {
              cpt++;
              cs_real_t lambda = (0.5 - porn[pt2]) / (porn[pt1] - porn[pt2]);

              for (int k = 0; k < 3; k++)
                xp[cpt][k] = lambda * xn[pt1][k] + (1. - lambda) * (xn[pt2][k]);
            }

            /* First point solid, second point fluid */
            if (porn[pt1] < 0.5 && porn[pt2] >= 0.5) {
              cpt++;
              cs_real_t lambda = (0.5 - porn[pt1]) / (porn[pt2] - porn[pt1]);
              for (int k = 0; k < 3; k++)
                xp[cpt][k] = lambda * xn[pt2][k] + (1. - lambda) * (xn[pt1][k]);
            }

            /* Second point fluid */
            if (porn[pt2] >= 0.5) {
              cpt++;
              for (int k = 0; k < 3; k++)
                xp[cpt][k] = xn[pt2][k];
            }

            if (cpt >= 100 || pt2 >= 100)
              bft_error(__FILE__, __LINE__, 0, "Arrays declared above should be resized above 100\n");
          }

          cpt++;
          for (int k = 0; k < 3; k++)
            xp[cpt][k] = xp[0][k];

          /* Center of gravity */
           for (int k = 0; k < 3; k++)
            cog[k] = 0;
          for (int ipt = 0; ipt < cpt; ipt++)
            for (int k = 0; k < 3; k++)
              cog[k] += xp[ipt][k];

          for (int k = 0; k < 3; k++)
            cog[k] /= (cs_real_t)cpt;

           /* Loop on triangles */
          for (int k = 0; k < 3; k++)
            cogf[k] = 0.;

          surf = 0.;

          /* Increment surface et i_face_cog */
          for (int tr = 0; tr < cpt; tr++) {
            cs_real_3_t cogloc;
            cs_real_3_t vectloc;

            for (int k = 0; k < 3; k++)
              cogloc[k] = (xp[tr][k] + xp[tr+1][k] + cog[k])*cs_math_1ov3;

            vectloc[0] = (xp[tr][1]-cog[1]) * (xp[tr+1][2]-cog[2])
                       - (xp[tr][2]-cog[2]) * (xp[tr+1][1]-cog[1]);
            vectloc[1] = (xp[tr][2]-cog[2]) * (xp[tr+1][0]-cog[0])
                       - (xp[tr][0]-cog[0]) * (xp[tr+1][2]-cog[2]);
            vectloc[2] = (xp[tr][0]-cog[0]) * (xp[tr+1][1]-cog[1])
                       - (xp[tr][1]-cog[1]) * (xp[tr+1][0]-cog[0]);

            cs_real_t surfloc = cs_math_3_norm(vectloc);
            surf += surfloc;

            for (int k = 0; k < 3; k++)
              cogf[k] += surfloc * cogloc[k];
          }

          /* Finalize i_face_cog */
          if (surf >= 2 * 1.e-6 * i_face_surf[f_id]) {
            for (int k = 0; k < 3; k++)
              cogf[k] /= surf;

            /* Final storage */
            for (int k = 0; k < 3; k++)
              i_f_face_cog[f_id][k] = cogf[k];

          }
          else {
            cs_real_t pnd = i_f_weight[f_id];
            for (int k = 0; k < 3; k++) {
              cs_real_t xip = cell_f_cen[c_id0][k];
              cs_real_t xjp = cell_f_cen[c_id1][k];

              cs_real_t xf = pnd * xip + (1. - pnd) * xjp;
              i_f_face_cog[f_id][k] = xf;
            }
          }
        }
      }

      if (cs::abs(i_f_dist[f_id]) > 1e-20) {
        /* Distance between the face center of gravity
           and the neighbor cell center
           and dot-product with the normal */
        cs_real_t dist2f = cs_math_3_distance_dot_product(i_f_face_cog[f_id],
                                                          cell_f_cen[c_id1],
                                                          i_face_u_normal[f_id]);
        i_f_weight[f_id] = dist2f / i_f_dist[f_id];
      }
      else {
        i_f_weight[f_id] = 0.5;
      }
    }
  }
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
_compute_iso_vol_porosity(const cs_mesh_t      *mesh,
                          cs_real_t            *c_poro,
                          const int            *comp_cell)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells    = mesh->n_cells;

  const cs_real_t *cell_vol  = mq_g->cell_vol;

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

    cs::parall::sum(volpor);
    cs_ibm->isovol = volpor;

  }
  else if (ipass > 1) {
    cs_real_t aa = 0.;
    cs_real_t bb = 0.;
    cs_real_t volpor = cs_ibm->isovol;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0) {
        cs_real_t porloc = c_poro[c_id];
        cs_real_t pump = cs::abs(porloc * (1. - porloc));
        aa += pump * cell_vol[c_id];
        bb += porloc * cell_vol[c_id];
      }

    cs::parall::sum(aa, bb);

    cs_real_t beta = (volpor - bb) / cs::max(aa, 1.e-20);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      if (comp_cell[c_id] > 0) {
        cs_real_t porloc = c_poro[c_id];
        cs_real_t pump = cs::abs(porloc * (1. - porloc));
        c_poro[c_id] += beta * pump;
        c_poro[c_id] = cs::clamp(c_poro[c_id], 0., 1.);
      }

    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_poro);
  }
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
_compute_cell_list_porosity(const cs_mesh_t *mesh,
                            cs_lnum_t       *comp_cell)
{
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  cs_array_lnum_fill_zero(n_cells_ext, comp_cell);

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

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    for (cs_lnum_t ind = mesh->b_face_vtx_idx[f_id]; ind < mesh->b_face_vtx_idx[f_id + 1]; ind++) {
      cs_lnum_t v_id = mesh->b_face_vtx_lst[ind];

      int iok = 0;
      for (int idim = 0; idim < 3; idim++) {
        if (   vtx_crd[v_id][idim] >= cs_ibm->xyzmin_moving_porosity[idim]
            && vtx_crd[v_id][idim] <= cs_ibm->xyzmax_moving_porosity[idim])
          iok++;
      }

      if (iok == 3)
        comp_cell[c_id] = 1;
    }
  }

  cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, comp_cell);

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
_compute_solid_surface_vector(const cs_mesh_t      *mesh,
                              cs_mesh_quantities_t *mesh_quantities,
                              const cs_real_t      *i_poro,
                              const cs_real_t      *b_poro)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;

  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mq_g->i_face_normal;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mq_g->b_face_normal;
  cs_real_t *c_w_face_surf = mesh_quantities->c_w_face_surf;
  cs_real_3_t *c_w_face_normal
    = (cs_real_3_t *)mesh_quantities->c_w_face_normal;

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int i = 0; i < 3; i++)
      c_w_face_normal[c_id][i] = 0;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t alpij = i_poro[f_id];
    for (int idim = 0; idim < 3; idim++) {
      c_w_face_normal[c_id0][idim] += alpij * i_face_normal[f_id][idim];
      c_w_face_normal[c_id1][idim] -= alpij * i_face_normal[f_id][idim];
    }

  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    cs_real_t alpij = b_poro[f_id];
    for (int idim = 0; idim < 3; idim++)
      c_w_face_normal[c_id][idim] += alpij * b_face_normal[f_id][idim];

  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    c_w_face_surf[c_id] = cs_math_3_norm(c_w_face_normal[c_id]);
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
_compute_solid_surface_cog(const cs_mesh_t      *mesh,
                           cs_mesh_quantities_t *mesh_quantities,
                           const cs_real_t      *c_poro,
                           const cs_real_t      *v_poro,
                           const cs_real_t      *i_poro,
                           const cs_real_t      *b_poro)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_t *cell_vol = mq_g->cell_vol;
  const cs_real_3_t *cell_cen = mq_g->cell_cen;
  const cs_real_3_t *cell_f_cen = mesh_quantities->cell_cen;
  cs_real_3_t *cell_s_cen
    = (cs_real_3_t *)(cs_field_by_name("cell_s_cen")->val);
  const cs_real_3_t *i_face_cog = mq_g->i_face_cog;
  const cs_real_3_t *i_f_face_cog = mesh_quantities->i_face_cog;
  const cs_real_3_t *i_face_normal
    = (const cs_real_3_t *)mq_g->i_face_normal;
  const cs_nreal_3_t *i_face_u_normal = mq_g->i_face_u_normal;
  const cs_real_3_t *b_face_cog = mq_g->b_face_cog;
  const cs_real_3_t *b_f_face_cog = mesh_quantities->b_face_cog;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)mq_g->b_face_normal;
  const cs_nreal_3_t *b_face_u_normal = mq_g->b_face_u_normal;
  const cs_real_t *c_w_face_surf = mesh_quantities->c_w_face_surf;
  const cs_real_3_t *c_w_face_normal
    = (const cs_real_3_t *)mesh_quantities->c_w_face_normal;
  cs_real_3_t *c_w_face_cog = (cs_real_3_t *)mesh_quantities->c_w_face_cog;
  cs_real_t *c_w_dist_inv = mesh_quantities->c_w_dist_inv;

  const cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t *i_face_vtx = mesh->i_face_vtx_lst;
  const cs_lnum_t *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_3_t *vtx_crd = (const cs_real_3_t *)mesh->vtx_coord;

  cs_real_t *poro_val = CS_F_(poro)->val;
  const cs_real_t *poro_val_pre = CS_F_(poro)->val_pre;

  static bool ipass = true;
  static int nbface_max = 24;
  int nbfm7 = 7 * nbface_max;

  while (ipass) {
    ipass = false;

    if (cs_ibm->nfc == nullptr)
      CS_MALLOC(cs_ibm->nfc, n_cells_ext, cs_lnum_t);

    if (cs_ibm->xandnfc == nullptr)
      CS_MALLOC(cs_ibm->xandnfc, nbfm7 * n_cells_ext, cs_real_t);
    else
      CS_REALLOC(cs_ibm->xandnfc, nbfm7 * n_cells_ext, cs_real_t);

    cs_array_lnum_fill_zero(n_cells_ext, cs_ibm->nfc);

    cs_gnum_t error = 0;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id0 = i_face_cells[f_id][0];
      cs_lnum_t c_id1 = i_face_cells[f_id][1];

      for (int idim = 0; idim < 3; idim++) {
        cs_ibm->xandnfc[nbfm7*c_id0 + 7*cs_ibm->nfc[c_id0] + idim]     =  i_face_cog[f_id][idim];
        cs_ibm->xandnfc[nbfm7*c_id0 + 7*cs_ibm->nfc[c_id0] + idim + 3] =  i_face_u_normal[f_id][idim];
        cs_ibm->xandnfc[nbfm7*c_id1 + 7*cs_ibm->nfc[c_id1] + idim]     =  i_face_cog[f_id][idim];
        cs_ibm->xandnfc[nbfm7*c_id1 + 7*cs_ibm->nfc[c_id1] + idim + 3] = -i_face_u_normal[f_id][idim];
      }

      cs_ibm->xandnfc[nbfm7*c_id0 + 7*cs_ibm->nfc[c_id0] + 6] = i_poro[f_id];
      cs_ibm->xandnfc[nbfm7*c_id1 + 7*cs_ibm->nfc[c_id1] + 6] = i_poro[f_id];

      if (cs_ibm->nfc[c_id0] < nbface_max)
        cs_ibm->nfc[c_id0]++;
      else
        error = 1;

      if (cs_ibm->nfc[c_id1] < nbface_max)
        cs_ibm->nfc[c_id1]++;
      else
        error = 1;
    }

    cs::parall::sum(error);
    /* Increase nbface_max if there is memory issue */
    if (error > 0) {
      ipass = true;
      nbface_max += 10;
      nbfm7 = 7 * nbface_max;
    }
  }

  cs_array_3d<cs_real_t> cut(n_cells_ext, 3, 3);
  cut.zero();

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int i = 0; i < 3; i++)
      cut(c_id, i, i) -= cell_vol[c_id] * c_poro[c_id];

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    cs_real_t alpij = i_poro[f_id];

    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
        cut(c_id0, k, l) += alpij * i_f_face_cog[f_id][k] * i_face_normal[f_id][l];
        cut(c_id1, k, l) -= alpij * i_f_face_cog[f_id][k] * i_face_normal[f_id][l];
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];
    cs_real_t alpij = b_poro[f_id];

    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        cut(c_id, k, l) += alpij * b_f_face_cog[f_id][k] * b_face_normal[f_id][l];
  }

  cs_array<cs_real_t> denom2(n_cells_ext);

  cs_array_2d<cs_real_t> coord_node_out(n_cells_ext, 3);
  cs_array_2d<cs_real_t> coord_node_min(n_cells_ext, 3);
  cs_array_2d<cs_real_t> coord_node_max(n_cells_ext, 3);

  /* Compute positions out, min and max */
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    denom2[c_id] = 0.;
    for (cs_lnum_t i = 0; i < 3; i++) {
      coord_node_out(c_id, i) = 0.;
      coord_node_min(c_id, i) = 1.e20;
      coord_node_max(c_id, i) = -1.e20;
    }
  }

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
        pds = cs::min(pds, 1.);

        for (int idim = 0; idim < 3; idim++)
          coord_node_out(c_id0, idim) += pds * (1. - v_poro[v_id])
                                             * vtx_crd[v_id][idim];

        for (int idim = 0; idim < 3; idim++)
          coord_node_out(c_id1, idim) += pds * (1. - v_poro[v_id])
                                             * vtx_crd[v_id][idim];

        denom2[c_id0] += pds * (1. - v_poro[v_id]);
        denom2[c_id1] += pds * (1. - v_poro[v_id]);
      }

      for (int idim = 0; idim < 3; idim++) {
        coord_node_min(c_id0, idim) = cs::min(coord_node_min(c_id0, idim),
                                              vtx_crd[v_id][idim]);
      }
      for (int idim = 0; idim < 3; idim++) {
        coord_node_min(c_id1, idim) = cs::min(coord_node_min(c_id1, idim),
                                              vtx_crd[v_id][idim]);

      }

      for (int idim = 0; idim < 3; idim++)
        coord_node_max(c_id0, idim) = cs::max(coord_node_max(c_id0, idim),
                                              vtx_crd[v_id][idim]);
      for (int idim = 0; idim < 3; idim++)
        coord_node_max(c_id1, idim) = cs::max(coord_node_max(c_id1, idim),
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
        pds = cs::min(pds, 1.);

        for (int idim = 0; idim < 3; idim++)
          coord_node_out(c_id, idim) += pds * (1. - v_poro[v_id])
                                            * vtx_crd[v_id][idim];

        denom2[c_id] += pds * (1. - v_poro[v_id]);
      }

      for (int idim = 0; idim < 3; idim++) {
        coord_node_min(c_id, idim) = cs::min(coord_node_min(c_id, idim),
                                             vtx_crd[v_id][idim]);
      }

      for (int idim = 0; idim < 3; idim++)
        coord_node_max(c_id, idim) = cs::max(coord_node_max(c_id, idim),
                                             vtx_crd[v_id][idim]);
    }
  } /* End of loop on boundary faces */

  /* Finalization */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    for (int idim = 0; idim < 3; idim++)
      coord_node_out(c_id, idim) /= cs::max(denom2[c_id], 1.e-20);

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
      for (int i = 0; i < 3; i++) {
        xp[i] = cs_math_3_dot_product(cut.sub_array(c_id, i), nx) / n2;
      }

      int iok = 1;
      int nfc_loc = cs_ibm->nfc[c_id];

      cs_real_t *xfc_loc = cs_ibm->xandnfc;

      for (int i = 0; i < nfc_loc; i++) {
        cs_real_t psca = (xfc_loc[nbfm7*c_id + 7*i    ] - xp[0]) * xfc_loc[nbfm7*c_id + 7*i + 3]
                       + (xfc_loc[nbfm7*c_id + 7*i + 1] - xp[1]) * xfc_loc[nbfm7*c_id + 7*i + 4]
                       + (xfc_loc[nbfm7*c_id + 7*i + 2] - xp[2]) * xfc_loc[nbfm7*c_id + 7*i + 5];
        psca *= xfc_loc[nbfm7*c_id + 7*i + 6];

        if (psca < -1.e-10)
          iok = 0;
      }

      int iok2 = 1;
      if (c_poro[c_id] > 0.75) {
        cs_real_t psca3
          = cs_math_3_distance_dot_product(cell_f_cen[c_id], xp, nx);

        if (psca3 > 0.)
          iok2 = 0;
      }

      if ((iok == 0 || iok2 == 0) && denom2[c_id] >= 1.e-20)
        for (int i = 0; i < 3; i++) {

          xp[i] = cell_cen[c_id][i]
                + 2. * (c_poro[c_id] - 0.5)
                     * (coord_node_out(c_id, i) - cell_cen[c_id][i]);
        }

      for (int i = 0; i < 3; i++) {
        xp[i] = cs::clamp(xp[i], coord_node_min(c_id, i), coord_node_max(c_id, i));
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
  cs_real_t smin = 1.e-20;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)   {
    cs_real_t poro = 0.5*(poro_val_pre[c_id] + poro_val[c_id]);

    if (poro > 1.e-6 && poro < 1. - 1.e-6) {
      cs_real_3_t nxyz;
      for (int i = 0; i < 3; i++)
        nxyz[i] = c_w_face_normal[c_id][i];

      cs_real_t nn = cs::max(cs_math_3_norm(nxyz), smin);

      for (int i = 0; i < 3; i++)
        nxyz[i] /= nn;

      cs_real_t psca = cs_math_3_distance_dot_product(c_w_face_cog[c_id],
                                                      cell_f_cen[c_id], nxyz);
      if (psca < 0.) {
        for (int i = 0; i < 3; i++)
          c_w_face_cog[c_id][i] = cell_cen[c_id][i]
                                + 0.9 * (cell_f_cen[c_id][i] - cell_cen[c_id][i]);
      }
    }
  }

  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)c_w_face_cog, 3);
  cs_halo_sync_var_strided(mesh->halo, CS_HALO_STANDARD,
                           (cs_real_t *)cell_s_cen, 3);

  // Compute the "characteristic size" vector of each cell
  cs_array_2d<cs_real_t> cell_length_3d(n_cells_ext, 3);
  cs_array_2d<cs_real_t> weight_3d(n_cells_ext, 3);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      cell_length_3d(c_id, i) = 0.;
      weight_3d(c_id, i) = 0.;
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id0 = i_face_cells[f_id][0];
    cs_lnum_t c_id1 = i_face_cells[f_id][1];

    for (int idim = 0; idim < 3; idim++) {
      cell_length_3d(c_id0, idim)
        += cs::abs((cell_cen[c_id1][idim] - cell_cen[c_id0][idim])
                   * i_face_u_normal[f_id][idim]);
      cell_length_3d(c_id1, idim)
        += cs::abs((cell_cen[c_id1][idim] - cell_cen[c_id0][idim])
                   * i_face_u_normal[f_id][idim]);

      weight_3d(c_id0, idim) += cs::abs(i_face_u_normal[f_id][idim]);
      weight_3d(c_id1, idim) += cs::abs(i_face_u_normal[f_id][idim]);
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    for (int idim = 0; idim < 3; idim++) {
      cell_length_3d(c_id, idim)
        += cs::abs(2. * (b_face_cog[f_id][idim] - cell_cen[c_id][idim])
                   * b_face_u_normal[f_id][idim]);
      weight_3d(c_id, idim) += cs::abs(b_face_u_normal[f_id][idim]);
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
    for (int i = 0; i < 3; i++)
      cell_length_3d(c_id, i) /= weight_3d(c_id, i);

  /* c_w_dist_inv : inverse distance from cell_f_cen to c_w_face_cog */
  cs_real_t loc_min = 1.e-20;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_3_t npx, xip, xpp;
    for (int idim = 0; idim < 3; idim++) {
      xip[idim] = cell_f_cen[c_id][idim];
      xpp[idim] = c_w_face_cog[c_id][idim];
      npx[idim] = c_w_face_normal[c_id][idim];
    }

    cs_real_t nn = c_w_face_surf[c_id];
    nn = cs::max(nn, loc_min);

    for (int idim = 0; idim < 3; idim++)
      npx[idim] /= nn;

    cs_real_t dwall
      = cs::abs(cs_math_3_distance_dot_product(xpp, xip, npx));

    cs_real_t poro = 0.5*(poro_val[c_id] + poro_val_pre[c_id]);

    for (int idim = 0; idim < 3; idim++)
      npx[idim] = cs::max(cs::abs(npx[idim]), loc_min);

    cs_real_t lambdax = cell_length_3d(c_id, 0) / npx[0];
    cs_real_t lambday = cell_length_3d(c_id, 1) / npx[1];
    cs_real_t lambdaz = cell_length_3d(c_id, 2) / npx[2];

    cs_real_t lambda = cs::min(lambdax, lambday);
    lambda = cs::min(lambda, lambdaz);

    cs_real_t dcell_geom_max = cs_math_3_norm(cell_length_3d.sub_array(c_id));

    cs_real_t dcell_geom = cs::min(lambda, dcell_geom_max) * 0.5
                         * cs::max(poro, 1.e-3);

    c_w_dist_inv[c_id]
      = 1./cs::max(dwall, 0.5 * dcell_geom);
    c_w_dist_inv[c_id]
      = cs::max(c_w_dist_inv[c_id], 1./(3. * dcell_geom));
  }

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, c_w_dist_inv);

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
  cs_ibm_object_t *new_obj = nullptr;

  CS_MALLOC(new_obj, 1, cs_ibm_object_t);

  /* Set name */
  if (name == nullptr || strcmp(name, "") == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty name provided for IBM object creation.\n"));

  new_obj->name = nullptr;
  CS_MALLOC(new_obj->name, strlen(name) + 1, char);
  strcpy(new_obj->name, name);

  /* Method */
  new_obj->method = method;

  /* Pointer to medcoupling or stl mesh structures */
  new_obj->cutcell_func = nullptr;
  new_obj->stl = nullptr;
  new_obj->mi  = nullptr;

  for (int i = 0; i < CS_N_IBM_OBJ_PROP_TYPES; i++)
    new_obj->property_defs[i] = nullptr;

  for (int i = 0; i < CS_N_IBM_OBJ_INIT_TYPES; i++)
    new_obj->init_vals_defs[i] = nullptr;

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

  if (obj != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating object: object \"%s\" already exists.\n"),
              name);

  int new_obj_id = cs_ibm->n_objects;

  if (new_obj_id == 0)
    CS_MALLOC(cs_ibm->objects, new_obj_id + 1, cs_ibm_object_t *);
  else
    CS_REALLOC(cs_ibm->objects, new_obj_id + 1, cs_ibm_object_t *);

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
  CS_FREE(obj->name);

  if (obj->cutcell_func != nullptr)
    obj->cutcell_func = nullptr;

  CS_FREE(obj);
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

  if (def != nullptr)
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

  if (def != nullptr)
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
 * \brief Try to get an object based on its name. Returns null if not found
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure, null if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name_try(const char *name)
{
  assert(name != nullptr);

  cs_ibm_object_t *obj = nullptr;

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
 * \return pointer to object structure, null if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name(const char *name)
{
  assert(name != nullptr);

  cs_ibm_object_t *obj = nullptr;

  for (int i = 0; i < cs_ibm->n_objects; i++) {
    if (strcmp(name, cs_ibm->objects[i]->name) == 0) {
      obj = cs_ibm->objects[i];
      break;
    }
  }

  if (obj == nullptr)
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
  cs_ibm_t  *ibm = nullptr;

  CS_MALLOC(ibm, 1, cs_ibm_t);

  ibm->n_objects = 0;
  ibm->objects = nullptr;

  ibm->prob_dim       = CS_IBM_3D;
  ibm->algo_choice    = CS_IBM_ALGO_CUT_CELLS;
  ibm->wall_condition = CS_IBM_WALL_LAW_WALL_CONDITION;
  ibm->min_poro       = 0.001;
  ibm->nb_cut_cells   = 2;
  ibm->nb_cut_edges   = 5;
  ibm->porosity_user_source_term_modification = false;
  ibm->nfc            = nullptr;
  ibm->xandnfc        = nullptr;
  ibm->isovol         = 0.;
  ibm->ensure_isovol  = false;
  ibm->remove_bad_cells = true;
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
 *   nullptr
 *----------------------------------------------------------------------------*/

void
cs_ibm_finalize(void)
{
  if (cs_ibm != nullptr) {

    /* Clean objects */
    for (int iobj = 0; iobj < cs_ibm->n_objects; iobj++)
      _free_ibm_object(cs_ibm->objects[iobj]);

    CS_FREE(cs_ibm->objects);

    CS_FREE(cs_ibm);
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

void cs_immersed_boundaries(cs_mesh_t      *mesh,
                            cs_mesh_quantities_t *mesh_quantities)
{
  cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells     = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = mesh->n_b_faces;

  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;

  const cs_real_3_t *cell_f_cen = mesh_quantities->cell_cen;
  const cs_real_t *cell_vol = mq_g->cell_vol;
  cs_real_t *cell_f_vol = mesh_quantities->cell_vol;
  int *c_disable_flag = mesh_quantities->c_disable_flag;

  cs_real_t timdeb = cs_timer_wtime();

  static int ipass = 0;
  ipass++;

  if (ipass == 1) {

    /* Initialization for time/space immersed boundaries */
    if (cs_glob_porosity_ibm_opt->porosity_mode > 0)
      cs_ibm = cs_ibm_create();
  }

  cs_array<int> error_cell(n_cells_ext);
  error_cell.zero();
  cs_array<int> final_error_cell(n_cells);
  final_error_cell.zero();

  /* Value pointers for fields */
  cs_real_t *ifpro_poro = cs_field("i_face_porosity")->val;
  cs_real_t *bfpro_poro = cs_field("b_face_porosity")->val;

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
      || nt_cur == 0
      || (cs_restart_present() && nt_cur == nt_prev)
      || (   _porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID
          && cs_ibm->porosity_user_source_term_modification)) {

    cs_real_3_t *gradp;

    CS_MALLOC(gradp, n_cells_ext, cs_real_3_t);

    int hyd_p_flag = cs_glob_velocity_pressure_param->iphydr;
    cs_real_3_t *f_ext = (hyd_p_flag == 1) ?
      (cs_real_3_t *)cs_field_by_name_try("volume_forces")->val : nullptr;

    bool use_previous_t = false;
    int inc = 1;

    if (nt_cur == nt_prev) {

      /* FIXME: gradp computation below needs cell_f_vol wich is computed
         after porosity. If possible, gradp could be computed after cell_f_vol.
         The following bad alternative solution is to copy cell_vol
         in cell_f_vol for initialization ...*/

      cs_array_real_copy(n_cells,
                         (const cs_real_t *)cell_vol,
                         cell_f_vol);
    }

    if (ipass >= 3) {
      cs_field_gradient_potential(CS_F_(p),
                                  use_previous_t,
                                  inc,
                                  hyd_p_flag,
                                  f_ext,
                                  gradp);
    }

    /* Initialize porosity at the previous value */
    cs_real_t *poro_val_pre = CS_F_(poro)->val_pre;
    cs_real_t *poro_val = CS_F_(poro)->val;

    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
      poro_val_pre[c_id] = poro_val[c_id];

    /* List of cells for which one has to recompute porosity -> comp_cell */
    int *comp_cell;

    CS_MALLOC(comp_cell, n_cells_ext, int);
    _compute_cell_list_porosity(mesh, comp_cell);

    /* Compute cell porosity */
    if (   nt_cur <=1
       || !cs_ibm->porosity_user_source_term_modification
       || _porosity_ibm_opt.porosity_mode != CS_IBM_FIXED_SOLID) {

      /* Compute cell porosity with a cut-cell method
       * from the cs_user_ibm function */
      if (cs_ibm->algo_choice == CS_IBM_ALGO_CUT_CELLS) {

        _compute_cell_cut_porosity(mesh, comp_cell);

        for(int c_id = 0; c_id < n_cells_ext; c_id++)
          if (poro_val[c_id] < 1.e-4)
            poro_val[c_id] = 0.;

      /* Compute cell porosity from a file */
      /* Object logic in this section */
      }
      else {

        /* User imposed rotations/translations */
        cs_user_ibm_object_transformations(cs_glob_time_step->t_cur);

        /* Some initialization */
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          poro_val[c_id] = 1.;

        /* Local declarations */
        cs_real_t *obj_vol_f_tot = nullptr;
        CS_MALLOC(obj_vol_f_tot, n_cells_ext, cs_real_t);
        cs_array_real_fill_zero(n_cells_ext, obj_vol_f_tot);

        for (int o_id = 0; o_id < cs_ibm->n_objects; o_id++) {
          cs_ibm_object_t *ibm_obj = cs_ibm_object_by_id(o_id);

          /* Increment total solid volume fraction using the current object.
           * Update object indicator array if needed. */
          cs_ibm_object_compute_intersect_vol(ibm_obj,
                                              mesh,
                                              cell_vol,
                                              obj_vol_f_tot,
                                              nullptr);

        }

        /* Once total volume of solids in each cell is computed,
         * substract the result from the porosity. */
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          poro_val[c_id] -= obj_vol_f_tot[c_id];

        /* Clip sync porosity values */
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          poro_val[c_id] =
            cs::min(cs::max(poro_val[c_id], 0.), 1.);

        cs_halo_sync(mesh->halo, CS_HALO_STANDARD, poro_val);

        CS_FREE(obj_vol_f_tot);
      }
      /* Additional safety */
      if (_porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID)
        cs_array_real_copy(n_cells_ext, poro_val, CS_F_(poro)->val_pre);
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
          _compute_iso_vol_porosity(mesh, poro_val, comp_cell);

    /* Disable cell and update fluid volume */
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      if (poro_val[c_id] < cs_ibm->min_poro) {
        poro_val[c_id] = 0.;
        c_disable_flag[c_id] = 1;
      }
      cell_f_vol[c_id] = poro_val[c_id] * cell_vol[c_id];
    }

    /* Porosity projection at vertices */
    cs_real_t *v_poro;
    CS_MALLOC(v_poro, mesh->n_vertices, cs_real_t);

    /* Automatic sub-cycles activated to remove bad cells */
    int ncycle = 1;
    bool remove_bad_cells = false;
    if (   _porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID
        && cs_ibm->remove_bad_cells) {
      remove_bad_cells = true;
      ncycle = 1000000;
    }

    /* Save cell cog before modification */
    cs_real_3_t *cog_save;
    CS_MALLOC(cog_save, n_cells, cs_real_3_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      for (int idir = 0; idir < 3; idir++)
        cog_save[c_id][idir] = cell_f_cen[c_id][idir];

    bool leave_cycle = false;
    for (int icycle = 0; icycle < ncycle; icycle++) {

      /* Update val_pre at the first time step */
      if (   cs_glob_time_step->nt_cur <= 1
          && _porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID)
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          poro_val_pre[c_id] = poro_val[c_id];

      /* Projection at vertices */
      cs_real_t *c_poro;
      CS_MALLOC(c_poro, n_cells_ext, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        c_poro[c_id] = 0.5*(poro_val_pre[c_id] + poro_val[c_id]);

      cs_cell_to_vertex<1>(CS_CELL_TO_VERTEX_SHEPARD, 0, false, nullptr,
                           c_poro, nullptr, v_poro);

      /* Compute cell centers of gravity */
      _compute_cell_cog(mesh, mesh_quantities,
                        v_poro, poro_val, comp_cell, icycle);

      /* One guarantees the same volume at each iteration
       * for cells whose porosity did not change */
      if (nt_cur > 1 && nt_cur > nt_prev)
        if (cs_ibm->ensure_isovol)
          if (!cs_ibm->porosity_user_source_term_modification)
            _compute_iso_vol_porosity(mesh, poro_val, comp_cell);

      /* Boundary face porosity */
      _compute_b_fac_porosity(mesh, mesh_quantities,
                              poro_val, v_poro, bfpro_poro);

      cs_array<bool> comp_fac(n_i_faces);
      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id0 = i_face_cells[f_id][0];
        cs_lnum_t c_id1 = i_face_cells[f_id][1];

        comp_fac[f_id] = false;
        if (comp_cell[c_id0] + comp_cell[c_id1] > 0)
          comp_fac[f_id] = true;
      }

      /* Internal face porosity */
      _compute_i_fac_porosity(mesh, mesh_quantities, c_poro, ifpro_poro);

      /* Compute solid surface vector */
      _compute_solid_surface_vector(mesh, mesh_quantities,
                                    ifpro_poro, bfpro_poro);

      /* Compute cog internal faces + weight + dist */
      _compute_i_fac_cog(mesh, mesh_quantities, c_poro, v_poro, comp_fac);

      /* Compute cog solid face + dist_wall */
      _compute_solid_surface_cog(mesh, mesh_quantities,
                                 poro_val, v_poro,
                                 ifpro_poro, bfpro_poro);

      cs_mesh_quantities_compute_face_vectors
        (mesh->dim,
         n_i_faces,
         n_b_faces,
         i_face_cells,
         b_face_cells,
         mesh_quantities->b_face_u_normal,
         (const cs_real_t *)(mesh_quantities->i_face_cog),
         (const cs_real_t *)(mesh_quantities->b_face_cog),
         (const cs_real_t *)(mesh_quantities->cell_cen),
         mesh_quantities->weight,
         mesh_quantities->b_dist,
         mesh_quantities->diipb,
         (cs_real_t *)(mesh_quantities->dofij));

      cs_real_3_t *dofij = mesh_quantities->dofij;

      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

        const cs_lnum_t c_id1 = i_face_cells[f_id][0];
        const cs_lnum_t c_id2 = i_face_cells[f_id][1];

        if (   c_disable_flag[c_id1] == 1
            || c_disable_flag[c_id2] == 1
            || ifpro_poro[f_id] < cs_ibm->min_poro) {
          dofij[f_id][0] = 0.;
          dofij[f_id][1] = 0.;
          dofij[f_id][2] = 0.;
        }

      }

      cs_mesh_quantities_compute_face_sup_vectors
        (n_cells,
         n_i_faces,
         i_face_cells,
         mesh_quantities->i_face_u_normal,
         (const cs_real_3_t *)(mesh_quantities->i_face_normal),
         (const cs_real_3_t *)(mesh_quantities->i_face_cog),
         (const cs_real_3_t *)(mesh_quantities->cell_cen),
         mesh_quantities->cell_vol,
         mesh_quantities->i_dist,
         mesh_quantities->diipf,
         mesh_quantities->djjpf);

      cs_lnum_t nb_cell_errors = 0;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        final_error_cell[c_id] += error_cell[c_id];

      /* Update comp_cell for the next cycle */
      cs_array<int> comp_prev(n_cells_ext);
      comp_prev.zero();

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (error_cell[c_id] < 1)
          comp_cell[c_id] = 0;
        else {
          comp_cell[c_id] = 1;
          comp_prev[c_id] = 1;
        }
      }

      cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, comp_prev);

      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id0 = i_face_cells[f_id][0];
        cs_lnum_t c_id1 = i_face_cells[f_id][1];

        if (comp_prev[c_id0] > 0)
          comp_cell[c_id1] = 1;
        if (comp_prev[c_id1] > 0)
          comp_cell[c_id0] = 1;
      }

      cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, comp_cell);

      CS_FREE(c_poro);

      if (leave_cycle)
        break;

      if (remove_bad_cells && nb_cell_errors == 0 && !leave_cycle )
        leave_cycle = true;
    }

    if (cs_log_default_is_active()) {

      cs_lnum_t total_nb_cell_errors = 0;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        if (final_error_cell[c_id] > 0)
          total_nb_cell_errors++;

      cs::parall::sum(total_nb_cell_errors);

      bft_printf("\n   ** Number of removed bad cells:  %d\n", total_nb_cell_errors);
      bft_printf("      ----------------------------\n");
    }

    CS_FREE(comp_cell);
    CS_FREE(v_poro);

    /* Pressure update after moving cell cog */
    cs_real_3_t dcog;
    if (ipass >= 3) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (int idir = 0; idir < 3; idir++)
          dcog[idir] = cell_f_cen[c_id][idir] - cog_save[c_id][idir];

        CS_F_(p)->val[c_id] += cs_math_3_dot_product(dcog, gradp[c_id]);
      }

      cs_field_synchronize(CS_F_(p), CS_HALO_STANDARD);
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, CS_F_(p)->val_pre);
    }

    CS_FREE(cog_save);
    CS_FREE(gradp);

  }

  /* Immersed boundary cells number */
  cs_lnum_t n_ib_cells = 0;

  cs_real_t *c_poro = CS_F_(poro)->val;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (c_poro[c_id] < 1. && c_poro[c_id] > 0.) {
      n_ib_cells++;
    }
  }

  /* Immersed boundary cells number */
  cs_lnum_t n_ib_cells_verif = 0;

  /* Connectivity ib_cell to cells */
  cs_lnum_t *ibcell_cells; // disable_cells included
  CS_MALLOC(ibcell_cells, n_ib_cells, cs_lnum_t);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (c_poro[c_id] < 1. && c_poro[c_id] > 0.) {
      ibcell_cells[n_ib_cells_verif] = c_id;
      n_ib_cells_verif++;
    }
  }

  // True connectivity skipping the disable cells
  cs_lnum_t *ibcell_cells_filt;
  CS_MALLOC(ibcell_cells_filt, n_ib_cells, cs_lnum_t);

  cs_lnum_t n_ib_cells_filt = 0;
  cs_lnum_t *ib_cells_filt;
  CS_MALLOC(ib_cells_filt, n_ib_cells, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_ib_cells; i++) {
    const cs_lnum_t c_id = ibcell_cells[i];
    if (c_disable_flag[c_id] == 0) {
      ibcell_cells_filt[n_ib_cells_filt] = c_id;
      ib_cells_filt[n_ib_cells_filt] = i;
      n_ib_cells_filt++;
    }
  }

  cs_real_3_t *c_w_face_normal
    = (cs_real_3_t *)mesh_quantities->c_w_face_normal;
  // Correction of orientation of c_w_face_normal to point towards solid region
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      for (cs_lnum_t i = 0; i < 3; i++)
        c_w_face_normal[c_id][i] = - c_w_face_normal[c_id][i];
  }

  cs_compute_corr_grad_lin(mesh, mesh_quantities);

  cs_porous_model_convert_cell_to_boundary(n_ib_cells_filt, ibcell_cells_filt);

  /* Free memory */
  CS_FREE(ibcell_cells);

  CS_FREE(ibcell_cells_filt);
  CS_FREE(ib_cells_filt);

  if (cs_log_default_is_active()) {
    bft_printf("\n   ** Immersed boundary calculation:  %5.2f s\n", cs_timer_wtime() - timdeb);
    bft_printf("      ------------------------------\n\n");
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

  /* Sanity check */
  if (file_name == nullptr || strlen(file_name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Object '%s' was added with an empty file name...\n",
              __func__, name);

  int obj_id = _add_ibm_object(name, method);

  cs_ibm_object_t *obj = cs_ibm_object_by_id(obj_id);

  /* STL/MED objects are only rigid for the moment! */
  /*obj->solve_fsi = solve_fsi;
  // TODO create a different function?
  if (solve_fsi) {
    CS_MALLOC(obj->fsi_index, 1, int);
    obj->fsi_index[0] = cs_fsi_object->number;
    cs_fsi_object->number += 1;
  }*/

  /* STL */
  if (method == CS_IBM_ALGO_STL) {
    obj->stl = cs_stl_mesh_add(name);
    cs_stl_file_read(obj->stl, file_name);
    obj->stl->is_porous = true;
  }
  else if (method == CS_IBM_ALGO_MEDCOUPLING) {
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

  /* sanity check */
  if (cutcell_func == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Object '%s' is defined a nullptr instead of a valid function.\n",
              __func__, name);

  int obj_id = _add_ibm_object(name,
                               CS_IBM_ALGO_CUT_CELLS);

  cs_ibm_object_t *obj = cs_ibm_object_by_id(obj_id);

  // TODO create a different function?
  /*obj->solve_fsi = solve_fsi;
  if (solve_fsi) {
    int id0 = cs_fsi_object->number;
    if (n_nodes > 1) {
      obj->is_deformable = true;
      CS_MALLOC(obj->fsi_index, n_nodes - 1, int);
      for (int i = 0; i < n_nodes - 1; i++)
        obj->fsi_index[i] = id0 + i;

      cs_fsi_object->number  += n_nodes - 1;

      obj->deform_id = cs_fsi_object->n_solid;
      cs_fsi_object->n_solid += 1;
      obj->n_nodes = n_nodes;

    }
    else {
      obj->is_deformable = false;
      CS_MALLOC(obj->fsi_index, 1, int);
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

  /* Default options */
  bool      output_at_start = tc->at_start;
  bool      output_at_end   = tc->at_end;
  int       interval_n      = tc->interval_nt;
  cs_real_t interval_t      = tc->interval_t;
  fvm_writer_time_dep_t time_dep = FVM_WRITER_TRANSIENT_COORDS;

  /* If non-moving objects, we force a single output only */
  if (_porosity_ibm_opt.porosity_mode == CS_IBM_FIXED_SOLID) {
    output_at_start = true;
    output_at_end   = false;
    interval_n      = -1;
    interval_t      = -1.0;
    time_dep = FVM_WRITER_FIXED_MESH;
  }

  switch(cs_ibm->algo_choice) {

  case CS_IBM_ALGO_STL:
    {
      cs_stl_post_init_writer("STL_OBJECTS",
                              "postprocessing",
                              "Ensight Gold",
                              "",
                              time_dep,
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
                             time_dep,
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

  cs_array<cs_real_t> wfrac(m->n_cells_with_ghosts);
  cs_array<cs_lnum_t> windic(m->n_cells_with_ghosts);

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

  return;
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
  const cs_real_3_t *cell_cen = mesh_quantities->cell_cen;
  int n_v_zones = cs_volume_zone_n_zones();

  cs_tree_node_t *tn_p
    = cs_tree_get_node(cs_glob_tree,
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

      if (formula != nullptr) {
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
