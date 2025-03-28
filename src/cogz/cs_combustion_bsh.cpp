/*============================================================================
 * Burke Schumann combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"
#include "cogz/cs_combustion_bsh.h"

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_bsh.cpp
        Burke Schumann combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! Burke Schumann combustion model constants */

#define N_Z 80
#define N_XR 5
#define N_ZVAR 10
#define N_VAR_BSH 7

/*! Burke Schumann combustion model variables */

cs_real_t        cs_bsh_coeff_therm[7][2][5];
static cs_real_t bsh_lib[N_XR][N_Z][N_VAR_BSH];
static cs_real_t turb_bsh_lib[CS_BSH_NVAR_TURB][N_XR][N_ZVAR][N_Z];

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute temperature.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_temperature
(
  cs_real_t h[N_Z][N_XR],
  cs_real_t t[N_Z][N_XR],
  cs_real_t yspece[N_Z][CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS]
)
{
  cs_lnum_t iter, ie, k, icp;
  cs_real_t som1, som2, dvar, var1, var;

  cs_combustion_gas_model_t *cm    = cs_glob_combustion_gas_model;
  const cs_lnum_t            ngase = cm->n_gas_el_comp;

  cs_real_t deltat = 1.0e-7;

  for (int ixr = 0; ixr < N_XR; ixr++) {
    for (int iz = 0; iz < N_Z; iz++) {

      var  = t[iz][ixr];
      iter = 0;

      do {

        iter++;
        var1 = var;
        som1 = 0.0;
        som2 = 0.0;
        icp  = 1;

        if (var1 > 1000.0)
          icp = 0;

        for (ie = 0; ie < ngase; ie++) {
          som1 += cm->coeff_therm[5][icp][ie] * yspece[iz][ie];
          for (k = 0; k < 5; k++) {
            som1 += cm->coeff_therm[k][icp][ie] / (k + 1) * pow(var1, k + 1)
                    * yspece[iz][ie];
            som2 += cm->coeff_therm[k][icp][ie] * pow(var1, k) * yspece[iz][ie];
          }
        }

        var  = var1 - (som1 - h[iz][ixr]) / som2;
        dvar = fabs(var - var1);

      } while (dvar > deltat && iter < 10000);

      t[iz][ixr] = var;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_density
(
  cs_real_t rho[N_Z][N_XR],
  cs_real_t t[N_Z][N_XR],
  cs_real_t yspece[N_Z][CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS]
)
{
  cs_lnum_t iz, ixr, ie;
  cs_real_t som1, p0;

  cs_combustion_gas_model_t *cm    = cs_glob_combustion_gas_model;
  const cs_lnum_t            ngase = cm->n_gas_el_comp;

  p0 = 101325.0;

  for (iz = 0; iz < N_Z; iz++) {
    for (ixr = 0; ixr < N_XR; ixr++) {
      som1 = 0.0;
      for (ie = 0; ie < ngase; ie++) {
        som1 += yspece[iz][ie] / cm->wmole[ie];
      }
      rho[iz][ixr] = p0 * (1.0 / (8.31433 * t[iz][ixr] * som1));
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief 1D interpolation.
 */
/*----------------------------------------------------------------------------*/

static void
_interp1d(const cs_real_t *xdata,
          const cs_real_t *ydata,
          const cs_real_t *xval,
          cs_real_t       *yval,
          cs_lnum_t        xdata_size,
          cs_lnum_t        xval_size)
{
  cs_lnum_t inputindex, dataindex;
  cs_real_t weight;

  for (inputindex = 0; inputindex < xval_size; ++inputindex) {
    dataindex = 0;
    while (xdata[dataindex] < xval[inputindex] && dataindex < xdata_size - 1) {
      dataindex++;
    }

    weight = (xval[inputindex] - xdata[dataindex - 1])
             / (xdata[dataindex] - xdata[dataindex - 1]);

    yval[inputindex]
      = (1.0 - weight) * ydata[dataindex - 1] + weight * ydata[dataindex];

    if (xdata[0] > xval[inputindex]) {
      yval[inputindex] = ydata[0];
    }

    if (xdata[xdata_size - 1] < xval[inputindex]) {
      yval[inputindex] = ydata[xdata_size - 1];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief 1D interpolation for the scalar.
 */
/*----------------------------------------------------------------------------*/

static void
_interp1d_scalar(const cs_real_t *xdata,
                 const cs_real_t *ydata,
                 cs_real_t        xval,
                 cs_real_t       *yval,
                 cs_lnum_t        xdata_size)
{
  cs_lnum_t dataindex = 0;
  cs_real_t weight;

  while (xdata[dataindex] < xval && dataindex < xdata_size - 1) {
    dataindex++;
  }

  weight
    = (xval - xdata[dataindex - 1]) / (xdata[dataindex] - xdata[dataindex - 1]);
  *yval = (1.0 - weight) * ydata[dataindex - 1] + weight * ydata[dataindex];

  if (xdata[0] >= xval) {
    *yval = ydata[0];
  }

  if (xdata[xdata_size - 1] <= xval) {
    *yval = ydata[xdata_size - 1];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Partition.
 */
/*----------------------------------------------------------------------------*/

static void
_partition_zz(cs_lnum_t m0, cs_lnum_t n0, cs_real_t eps, cs_real_t *zz)
{
  cs_lnum_t nsubrange, ii, jj;
  cs_real_t dz;

  nsubrange = (cs_lnum_t)log10(0.1 / eps);

  zz[0] = eps;

  for (ii = 1; ii <= nsubrange; ii++) {
    dz = (pow(10.0, ii) * eps - pow(10.0, ii - 1) * eps) / n0;
    for (jj = 1; jj <= n0; jj++) {
      zz[(ii - 1) * n0 + jj] = zz[(ii - 1) * n0 + jj - 1] + dz;
    }
  }

  dz = (0.9 - 0.1) / m0;
  for (ii = nsubrange * n0; ii < nsubrange * n0 + m0; ii++) {
    zz[ii + 1] = zz[ii] + dz;
  }

  for (ii = 1; ii <= nsubrange; ii++) {
    dz = (pow(10.0, nsubrange - ii + 1) * eps - pow(10.0, nsubrange - ii) * eps)
         / n0;
    for (jj = 1; jj <= n0; jj++) {
      zz[n0 * nsubrange + m0 + (ii - 1) * n0 + jj]
        = zz[m0 + n0 * nsubrange + (ii - 1) * n0 + jj - 1] + dz;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the phi/beta integral.
 */
/*----------------------------------------------------------------------------*/

static void
_integral_phi_beta(cs_real_t *z,
                   cs_real_t *phi,
                   cs_real_t  z_mean,
                   cs_real_t  z_var,
                   cs_real_t *phi_integral,
                   cs_lnum_t  size_z)
{
  cs_real_t  sum1, sum2, sum3, denom, beta_pdf;
  cs_real_t  zz_max, aa, bb, eps, aa_max, bb_max, dz;
  cs_real_t *zz, *phi_;
  cs_lnum_t  m0, n0, dim_zz, ii;

  if (fabs(z_var) <= 1e-7) {
    _interp1d_scalar(z, phi, z_mean, phi_integral, size_z);
    return;
  }

  aa_max = 500.0;
  bb_max = 500.0;

  aa     = z_mean * (z_mean * (1.0 - z_mean) / z_var - 1.0);
  bb     = (1.0 - z_mean) * (z_mean * (1.0 - z_mean) / z_var - 1.0);
  zz_max = 1.0 / (1.0 + (bb - 1.0) / (aa - 1.0));

  if (aa > aa_max) {
    aa = aa_max;
    bb = (aa - 1.0 - zz_max * (aa - 2.0)) / zz_max;
  }

  if (bb > bb_max) {
    bb = bb_max;
    aa = (1.0 + zz_max * (bb - 2.0)) / (1.0 - zz_max);
  }

  m0  = 500;
  n0  = 50;
  eps = 1.0e-6;

  dim_zz = m0 + 2 * n0 * (cs_lnum_t)log10(0.1 / eps) + 1;
  CS_MALLOC(zz, dim_zz, cs_real_t);
  for (ii = 0; ii < dim_zz; ii++) {
    zz[ii] = 0.0;
  }

  CS_MALLOC(phi_, dim_zz, cs_real_t);
  for (ii = 0; ii < dim_zz; ii++) {
    phi_[ii] = 0.0;
  }

  _partition_zz(m0, n0, eps, zz);

  sum1 = pow(eps, aa) / aa;
  sum3 = pow(eps, bb) / bb;
  sum2 = 0.0;

  for (ii = 1; ii < dim_zz; ii++) {
    dz       = zz[ii] - zz[ii - 1];
    beta_pdf = pow(zz[ii - 1], aa - 1.0) * pow(1.0 - zz[ii - 1], bb - 1.0);
    sum2 += 0.5 * dz * beta_pdf;
    beta_pdf = pow(zz[ii], aa - 1.0) * pow(1.0 - zz[ii], bb - 1.0);
    sum2 += 0.5 * dz * beta_pdf;
  }

  denom = sum1 + sum2 + sum3;

  _interp1d(z, phi, zz, phi_, size_z, dim_zz);

  sum1 = phi_[0] * pow(eps, aa) / aa / denom;
  sum3 = phi_[dim_zz - 1] * pow(eps, bb) / bb / denom;
  sum2 = 0.0;

  for (ii = 1; ii < dim_zz; ii++) {
    dz = zz[ii] - zz[ii - 1];
    beta_pdf
      = pow(zz[ii - 1], aa - 1.0) * pow(1.0 - zz[ii - 1], bb - 1.0) / denom;
    sum2 += 0.5 * dz * phi_[ii - 1] * beta_pdf;
    beta_pdf = pow(zz[ii], aa - 1.0) * pow(1.0 - zz[ii], bb - 1.0) / denom;
    sum2 += 0.5 * dz * phi_[ii] * beta_pdf;
  }

  *phi_integral = sum1 + sum2 + sum3;

  CS_FREE(zz);
  CS_FREE(phi_);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the fluid properties from the Burke-Schumann combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_burke_schumann_properties(cs_real_t z_m_0,
                                     cs_real_t zvar_0,
                                     cs_real_t xr_m_0,
                                     cs_real_t phi_t[CS_BSH_NVAR_TURB])
{
  cs_real_t  xdata[N_Z];  // largest of NZ, N_ZVAR, N_XR;
  assert(N_Z >= N_ZVAR && N_Z >= N_XR);

  cs_real_t  phi_4[CS_BSH_NVAR_TURB][N_XR][N_ZVAR];
  cs_real_t  phi_3[CS_BSH_NVAR_TURB][N_XR];
  cs_real_t  weight_z_m, weight_z_var, weight_xr_m;
  cs_lnum_t  i, j, k, dataindex;

  // Start with z_m
  for (i = 0; i < N_Z; i++) {
    xdata[i] = turb_bsh_lib[0][0][0][i];
  }
  dataindex = 0;

  while (xdata[dataindex] < z_m_0 && dataindex < N_Z-1) {
    dataindex++;
  }

  if (dataindex > 0) {
    weight_z_m = (z_m_0 - xdata[dataindex - 1])
                 / (xdata[dataindex] - xdata[dataindex - 1]);

    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        for (k = 0; k < N_ZVAR; k++) {
          phi_4[i][j][k]
            = (1.0 - weight_z_m) * turb_bsh_lib[i][j][k][dataindex - 1]
              + weight_z_m * turb_bsh_lib[i][j][k][dataindex];
        }
      }
    }
  }

  if (xdata[0] >= z_m_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        for (k = 0; k < N_ZVAR; k++) {
          phi_4[i][j][k] = turb_bsh_lib[i][j][k][0];
        }
      }
    }
  }

  if (xdata[N_Z - 1] <= z_m_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        for (k = 0; k < N_ZVAR; k++) {
          phi_4[i][j][k] = turb_bsh_lib[i][j][k][N_Z - 1];
        }
      }
    }
  }

  // Then interpolate over z_var
  for (i = 0; i < N_ZVAR; i++) {
    xdata[i] = phi_4[1][0][i];
  }
  dataindex = 0;

  while (xdata[dataindex] < zvar_0 && dataindex < N_ZVAR) {
    dataindex++;
  }

  if (dataindex > 0) {
    weight_z_var = (zvar_0 - xdata[dataindex - 1])
                   / (xdata[dataindex] - xdata[dataindex - 1]);
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        phi_3[i][j] = (1.0 - weight_z_var) * phi_4[i][j][dataindex - 1]
                      + weight_z_var * phi_4[i][j][dataindex];
      }
    }
  }

  if (xdata[0] >= zvar_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        phi_3[i][j] = phi_4[i][j][0];
      }
    }
  }

  if (xdata[N_ZVAR - 1] <= zvar_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      for (j = 0; j < N_XR; j++) {
        phi_3[i][j] = phi_4[i][j][N_ZVAR - 1];
      }
    }
  }

  // Then interpolate over XR_m
  for (i = 0; i < N_XR; i++) {
    xdata[i] = phi_3[2][i];
  }
  dataindex = 0;

  while (xdata[dataindex] < xr_m_0 && dataindex < N_XR-1) {
    dataindex++;
  }

  if (dataindex > 0) {
    weight_xr_m = (xr_m_0 - xdata[dataindex - 1])
                  / (xdata[dataindex] - xdata[dataindex - 1]);
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      phi_t[i] = (1.0 - weight_xr_m) * phi_3[i][dataindex - 1]
                 + weight_xr_m * phi_3[i][dataindex];
    }
  }

  if (xdata[0] >= xr_m_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      phi_t[i] = phi_3[i][0];
    }
  }

  if (xdata[N_XR - 1] <= xr_m_0) {
    for (i = 0; i < CS_BSH_NVAR_TURB; i++) {
      phi_t[i] = phi_3[i][N_XR - 1];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the parameters needed for the Burke-Schumann combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_burke_schumann(void)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  // Variables of the Burke Schumann model
  const cs_lnum_t N_Z_left  = 30;
  const cs_lnum_t N_Z_right = N_Z - N_Z_left;
  cs_real_t       z[N_Z], dz_left, dz_right;

  // Space discretization of the defect enthalpy
  const int ngasem = CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS;
  const cs_real_t xr1[N_XR] = {0.0, 0.2, 0.4, 0.6, 0.8};
  const cs_real_t q = 0.5;

  cs_real_t hunburnt, had;
  cs_real_t yspecg[N_Z][CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  cs_real_t yspece[N_Z][CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
  cs_real_t yspec[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
  cs_real_t ytot;
  cs_real_t xr[N_Z][N_XR], h[N_Z][N_XR], t[N_Z][N_XR], rho[N_Z][N_XR];

  // Space discretization of the mixture fraction
  z[0]            = 0.0;
  z[N_Z_left - 1] = cm->fs[0];

  dz_left  = pow(cm->fs[0], q) / (cs_real_t)(N_Z_left - 1);
  dz_right = pow(1.0 - cm->fs[0], q) / (cs_real_t)N_Z_right;

  for (int ii = 1; ii < N_Z_left; ii++) {
    z[N_Z_left - ii - 1] = z[N_Z_left - 1] - pow(ii * dz_left, 1.0 / q);
  }

  for (int ii = 1; ii <= N_Z_right; ii++) {
    z[N_Z_left + ii - 1] = z[N_Z_left - 1] + pow(ii * dz_right, 1.0 / q);
  }

  // Calculation of the mass fraction of species
  for (int ii = 0; ii < N_Z; ii++) {
    if (z[ii] <= cm->fs[0]) {
      // Oxidiser side
      yspecg[ii][0] = 0.0;
      yspecg[ii][1] = 1.0 - z[ii] / cm->fs[0];
      yspecg[ii][2] = 1.0 - yspecg[ii][0] - yspecg[ii][1];
    }
    else {
      // Fuel side
      yspecg[ii][0] = (z[ii] - cm->fs[0]) / (1.0 - cm->fs[0]);
      yspecg[ii][1] = 0.0;
      yspecg[ii][2] = 1.0 - yspecg[ii][0] - yspecg[ii][1];
    }
  }

  for (int ii = 0; ii < N_Z; ii++) {
    for (int nn = 0; nn < cm->n_gas_el_comp - 1; nn++) {
      yspece[ii][nn] = 0.0;
      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        yspece[ii][nn] += cm->coefeg[igg][nn] * yspecg[ii][igg];
      }
    }
  }

  for (int ii = 0; ii < N_Z; ii++) {
    ytot = 0.0;
    for (int nn = 0; nn < cm->n_gas_el_comp - 1; nn++) {
      ytot += yspece[ii][nn];
    }
    yspece[ii][cm->n_gas_el_comp - 1] = 1.0 - ytot;
  }

  // Calculation of the enthalpy

  for (int ii = 0; ii < ngasem; ii++)
    yspec[ii] = 0;

  yspec[0]   = yspece[N_Z - 1][0];
  cm->hinfue = cs_compute_burke_schumann_enthalpy(cm->tinfue, yspec);

  yspec[0]   = 0.0;
  yspec[1]   = yspece[0][1];
  yspec[4]   = yspece[0][cm->n_gas_el_comp - 1];
  cm->hinoxy = cs_compute_burke_schumann_enthalpy(cm->tinoxy, yspec);

  for (int ii = 0; ii < N_Z; ii++) {
    had = z[ii] * cm->hinfue + (1.0 - z[ii]) * cm->hinoxy;

    for (int nn = 0; nn < cm->n_gas_el_comp; nn++) {
      yspec[nn] = yspece[ii][nn];
    }

    hunburnt = cs_compute_burke_schumann_enthalpy(cm->tinoxy, yspec);

    for (int mm = 0; mm < N_XR; mm++) {
      h[ii][mm] = xr1[mm] * (hunburnt - had) + had;
    }
    for (int mm = 0; mm < N_XR; mm++) {
      xr[ii][mm] = -(h[ii][mm] - had);
    }
  }

  // Calculation of temperature
  for (int ii = 0; ii < N_Z; ii++) {
    for (int mm = 0; mm < N_XR; mm++) {
      t[ii][mm] = cm->tinoxy;
    }
  }
  _compute_temperature(h, t, yspece);

  // Calculation of density
  _compute_density(rho, t, yspece);

  for (int mm = 0; mm < N_XR; mm++) {
    for (int ii = 0; ii < N_Z; ii++) {
      bsh_lib[mm][ii][0] = z[ii];
      bsh_lib[mm][ii][1] = xr[ii][mm];
      bsh_lib[mm][ii][2] = rho[ii][mm];
      bsh_lib[mm][ii][3] = t[ii][mm];
      bsh_lib[mm][ii][4] = yspecg[ii][0];
      bsh_lib[mm][ii][5] = yspecg[ii][1];
      bsh_lib[mm][ii][6] = yspecg[ii][2];
    }
  }

  cs_lnum_t idx, k;
  cs_real_t zvarmin, zvarmax, dzvar;
  cs_real_t zvar[N_Z * N_ZVAR], phi[N_Z], z_m[N_Z];

  // Initialize z array from bsh_lib
  for (int ii = 0; ii < N_Z; ii++) {
    z[ii] = bsh_lib[0][ii][0];
  }

  // Discretization on z_m
  for (int ii = 0; ii < N_Z; ii++) {
    z_m[ii] = z[ii];
  }

  z_m[0]       = z[0] + 1.0e-6;
  z_m[N_Z - 1] = z[N_Z - 1] - 1.0e-6;

  // Discretization on zvar
  for (int ii = 0; ii < N_Z; ii++) {

    zvarmin = 1.0e-7;
    zvarmax = z_m[ii] * (1.0 - z_m[ii]) * 0.999999;
    dzvar   = sqrt(zvarmax - zvarmin) / (cs_real_t)(N_ZVAR - 1);

    zvar[ii * N_ZVAR] = zvarmin;
    for (int jj = 1; jj < N_ZVAR; jj++) {
      zvar[ii * N_ZVAR + jj]
        = zvar[ii * N_ZVAR] + (cs_real_t)(jj * dzvar * jj * dzvar);
    }
  }

  for (int ii = 0; ii < CS_BSH_NVAR_TURB; ii++) {
    for (int jj = 0; jj < N_XR; jj++) {
      for (int kk = 0; kk < N_ZVAR; kk++) {
        for (int ll = 0; ll < N_Z; ll++) {
          turb_bsh_lib[ii][jj][kk][ll] = 0.0;
        }
      }
    }
  }

  // Laminar base
  for (int ii = 0; ii < N_Z; ii++) {      // Loop on the mixture fraction
    for (int jj = 0; jj < N_ZVAR; jj++) { // Loop on the variance
      for (int nn = 0; nn < N_XR; nn++) { // Loop on the defect enthalpy

        turb_bsh_lib[0][nn][jj][ii] = z_m[ii];
        turb_bsh_lib[1][nn][jj][ii] = zvar[ii * N_ZVAR + jj];

        for (k = 1; k < 7; k++) {
          for (idx = 0; idx < N_Z; idx++) {
            phi[idx] = bsh_lib[nn][idx][k];
          }
          if (k == 2) {
            for (idx = 0; idx < N_Z; idx++) {
              phi[idx] = 1.0 / phi[idx];
            }
          }
          _integral_phi_beta(z,
                             phi,
                             z_m[ii],
                             zvar[ii * N_ZVAR + jj],
                             &turb_bsh_lib[k + 1][nn][jj][ii],
                             N_Z);
        }

        turb_bsh_lib[3][nn][jj][ii] = 1.0 / turb_bsh_lib[3][nn][jj][ii];

        // T^2 calculation for the T variance
        for (idx = 0; idx < N_Z; idx++) {
          phi[idx] = bsh_lib[nn][idx][3] * bsh_lib[nn][idx][3];
        }
        _integral_phi_beta(z,
                           phi,
                           z_m[ii],
                           zvar[ii * N_ZVAR + jj],
                           &turb_bsh_lib[CS_BSH_NVAR_TURB - 2][nn][jj][ii],
                           N_Z);

        // T^4 calculation for the radiative emission
        for (idx = 0; idx < N_Z; idx++) {
          phi[idx] = pow(bsh_lib[nn][idx][3], 4.0);
        }
        _integral_phi_beta(z,
                           phi,
                           z_m[ii],
                           zvar[ii * N_ZVAR + jj],
                           &turb_bsh_lib[CS_BSH_NVAR_TURB - 1][nn][jj][ii],
                           N_Z);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute enthalpy using the Burke-Schumann model.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_compute_burke_schumann_enthalpy
(
 cs_real_t t,
 cs_real_t yspec[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS]
)
{
  cs_combustion_gas_model_t *cm    = cs_glob_combustion_gas_model;
  const int                  ngase = cm->n_gas_el_comp;

  // Determination of the set of coefficients used
  int icp = 1;
  if (t > 1000.0)
    icp = 0;

  cs_real_t h = 0.0;

  for (int ne = 0; ne < ngase; ne++) {
    cs_real_t he = cm->coeff_therm[5][icp][ne];
    for (int i = 0; i < 5; i++) {
      he += cm->coeff_therm[i][icp][ne] / (i + 1) * pow(t, i + 1);
    }
    h += he * yspec[ne];
  }
  return h;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
