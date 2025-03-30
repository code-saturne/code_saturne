/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "base/cs_time_step.h"
#include "base/cs_timer.h"

#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "rayt/cs_rad_transfer_rcfsk.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/* \file cs_rad_transfer_rcfsk.c.f90 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local static variables
 *============================================================================*/

static int ipass_rcfsk = 0;

static cs_real_t *kgdatabase;
static cs_real_t *agdatabase;

/*=============================================================================
 * Local const variables
 *============================================================================*/

const int nq_rcfsk    = 10; /* nq_rcfsk = nxsgg */
const int nt_rcfsk    = 28;
const int nconc_rcfsk = 9;
const int nfvs_rcfsk  = 6;

const cs_real_t t_kg_rcfsk[28]
  = {300.0,  400.0,  500.0,  600.0,  700.0,  800.0,  900.0,
     1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0,
     1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0,
     2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0};

const cs_real_t x_kg_rcfsk[9]
  = {1.e-5, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0};
const cs_real_t fvs_kg_rcfsk[6]
 = {1.e-15, 0.5e-6, 1.e-6, 1.5e-6, 2.0e-6, 1.e-5};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert a line of values in a file into an array of cs_real_t and
 *        returns the number of values read
 *
 * \param[in]  radfile     pointer to file for reading values
 * \param[in]  values      array for values storage
 * \param[in]  nvalues     number of values read
 */
/*----------------------------------------------------------------------------*/

static inline void
_line_to_array_rcfsk(FILE       *radfile,
                     cs_real_t   values[],
                     int        *nvalues)
{
  char line[256];
  fgets(line, 256, radfile);
  int index = 0;

  /* now read all info along the line */
  while (strlen(line) > 1) {
    char temp[256];
    /* store next value in string format */
    sscanf(line, "%s", temp);
    /* read and convert to double precision */
    sscanf(temp, "%lf", &(values[index]));
    /* increase index in array */
    index++;

    int l = strlen(temp);
    int i = 0;
    while (line[i] == ' ')
      i++;
    snprintf(temp, 256, "%s", &line[l+i]);
    strcpy(line, temp);
  }

  *nvalues = index;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief number of segments
 *
 * \param[in]     t             Local Temperature
 * \param[in]     xco2          CO2 volume fraction
 * \param[in]     xh2o          H2O volume fraction
 * \param[in]     fvs
 * \param[inout]  itx[4][4]     itx[0]: TPlanck; iTx[1]: Tloc;
 *                              iTx[2]: xCO2; iTx[3] = xH2O
 */
/*----------------------------------------------------------------------------*/

inline static void
_gridposnbsg1_rcfsk(cs_real_t t,
                    cs_real_t xco2,
                    cs_real_t xh2o,
                    cs_real_t fvs,
                    int       itx[4][2])
{

  int ita[2]   = { 0, 0 };
  int ico2a[2] = { 0, 0 };
  int ih2oa[2] = { 0, 0 };
  int ifvsa[2] = { 0, 0 };

  /* Soot volume fraction interpolation: determine xfvs */

  int ifvs  = 0;
  int jfvs  = nfvs_rcfsk - 1;
  int ipfvs = (ifvs + jfvs) / 2;

  while (jfvs - ifvs > 1) {
    if (fvs < fvs_kg_rcfsk[ipfvs])
      jfvs = ipfvs;
    else
      ifvs = ipfvs;
    ipfvs = (ifvs + jfvs) / 2;
  }

  if (ifvs < 0)
    ifvs = 0;
  else if (ifvs > nfvs_rcfsk - 2)
    ifvs = nfvs_rcfsk - 2;

  ifvsa[0] = ifvs;
  ifvsa[1] = ifvs + 1;

  /* Mole fraction interpolation: determine xCO2 */

  int i  = 0;
  int j  = nconc_rcfsk - 1;
  int ip = (i + j) / 2;

  while (j - i > 1) {
    if (xco2 < x_kg_rcfsk[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  if (i < 0)
    i = 0;
  else if (i > nconc_rcfsk - 2)
    i = nconc_rcfsk - 2;

  ico2a[0] = i;
  ico2a[1] = i + 1;

  /* Mole fraction interpolation: determine xH2O */

  i  = 0;
  j  = nconc_rcfsk - 1;
  ip = (i + j) / 2;

  while (j - i > 1) {
    if (xh2o < x_kg_rcfsk[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  if (i < 0)
    i = 0;
  else if (i > nconc_rcfsk - 2)
    i = nconc_rcfsk - 2;

  ih2oa[0] = i;
  ih2oa[1] = i + 1;

  /* Temperature interpolation */

  i  = 0;
  j  = nt_rcfsk - 1;
  ip = (i + j) / 2;

  while (j - i > 1) {
    if (t < t_kg_rcfsk[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  if (i < 0)
    i = 0;
  else if (i > nt_rcfsk - 2)
    i = nt_rcfsk - 2;

  ita[0] = i;
  ita[1] = i + 1;

  /* attribution itx */
  for (int k = 0; k < 2; k++) {
    itx[0][k] = ita[k];
    itx[1][k] = ico2a[k];
    itx[2][k] = ih2oa[k];
    itx[3][k] = ifvsa[k];
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief interpolation4d
 *
 * \param[in]     val_cal
 * \param[in]     t             Local temperature
 * \param[in]     xco2          Reference CO2 volume fraction
 * \param[in]     xh2o          Reference H2O volume fraction
 * \param[in]     fvs           Reference soot volume fraction
 * \param[out]    phi
 */
/*----------------------------------------------------------------------------*/

inline static void
_interpolation4d_rcfsk(int       val_cal,
                       cs_real_t t,
                       cs_real_t xco2,
                       cs_real_t xh2o,
                       cs_real_t fvs,
                       cs_real_t phi[])
{
  int itx[4][2];
  int nix = 2;
  int nit = 2;
  int nis = 2;

  cs_real_t *karray, *kint1, *kint2, *kint3;

  CS_MALLOC(karray, 2 * 2 * 2 * 2 * nq_rcfsk, cs_real_t);
  CS_MALLOC(kint1, 2 * 2 * 2 * nq_rcfsk, cs_real_t);
  CS_MALLOC(kint2, 2 * 2 * nq_rcfsk, cs_real_t);
  CS_MALLOC(kint3, 2 * nq_rcfsk, cs_real_t);

  /* number of interpolation points along t & x */
  _gridposnbsg1_rcfsk(t, xco2, xh2o, fvs, itx);

  /* Attribute interpolation point indexes along T & x */

  int ita[2]   = { 0, 0 };
  int ico2a[2] = { 0, 0 };
  int ih2oa[2] = { 0, 0 };
  int ifvsa[2] = { 0, 0 };
  int i, is, ih2o, ico2, it, ig;

  for (i = 0; i < nit; i++) {
    ita[i] = itx[0][i];
  }

  for (i = 0; i < nix; i++) {
    ico2a[i] = itx[1][i];
    ih2oa[i] = itx[2][i];
  }

  for (i = 0; i < nis; i++) {
    ifvsa[i] = itx[3][i];
  }

  for (is = 0; is < nis; is++) {
    for (ih2o = 0; ih2o < nix; ih2o++) {
      for (ico2 = 0; ico2 < nix; ico2++) {
        for (it = 0; it < nit; it++) {
          for (ig = 0; ig < nq_rcfsk; ig++) {
            if (val_cal == 1) {
              karray[is + ih2o * 2 + ico2 * 2 * 2 + it * 2 * 2 * 2
                     + ig * 2 * 2 * 2 * 2]
                = agdatabase[ifvsa[is] + ih2oa[ih2o] * nfvs_rcfsk
                             + ico2a[ico2] * nfvs_rcfsk * nconc_rcfsk
                             + ita[it] * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                             + ig * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk * nt_rcfsk];
            }
            else if (val_cal == 2) {
              karray[is + ih2o * 2 + ico2 * 2 * 2 + it * 2 * 2 * 2
                     + ig * 2 * 2 * 2 * 2]
                = kgdatabase[ifvsa[is] + ih2oa[ih2o] * nfvs_rcfsk
                             + ico2a[ico2] * nfvs_rcfsk * nconc_rcfsk
                             + ita[it] * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                             + ig * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk * nt_rcfsk];
            }
          }
        }
      }
    }
  }

  /* interpolation on SOOT */

  const cs_real_t wxfvs
    = (fvs - fvs_kg_rcfsk[ifvsa[0]]) / (fvs_kg_rcfsk[ifvsa[1]] - fvs_kg_rcfsk[ifvsa[0]]);

  for (ih2o = 0; ih2o < nix; ih2o++) {
    for (ico2 = 0; ico2 < nix; ico2++) {
      for (it = 0; it < nit; it++) {
        for (ig = 0; ig < nq_rcfsk; ig++) {
          kint1[ih2o + ico2 * 2 + it * 2 * 2 + ig * 2 * 2 * 2]
            = wxfvs
                * karray[1 + ih2o * 2 + ico2 * 2 * 2 + it * 2 * 2 * 2
                         + ig * 2 * 2 * 2 * 2]
              + (1.0 - wxfvs)
                  * karray[0 + ih2o * 2 + ico2 * 2 * 2 + it * 2 * 2 * 2
                           + ig * 2 * 2 * 2 * 2];
        }
      }
    }
  }

  /* Interpolation on XH2O */

  const cs_real_t wxh2o
    = (xh2o - x_kg_rcfsk[ih2oa[0]]) / (x_kg_rcfsk[ih2oa[1]] - x_kg_rcfsk[ih2oa[0]]);

  for (ico2 = 0; ico2 < nix; ico2++) {
    for (it = 0; it < nit; it++) {
      for (ig = 0; ig < nq_rcfsk; ig++) {
        kint2[ico2 + it * 2 + ig * 2 * 2]
          = wxh2o * kint1[1 + ico2 * 2 + it * 2 * 2 + ig * 2 * 2 * 2]
            + (1.0 - wxh2o) * kint1[0 + ico2 * 2 + it * 2 * 2 + ig * 2 * 2 * 2];
      }
    }
  }

  /* Interpolation on XCO2 */

  const cs_real_t wxco2
    = (xco2 - x_kg_rcfsk[ico2a[0]]) / (x_kg_rcfsk[ico2a[1]] - x_kg_rcfsk[ico2a[0]]);

  for (it = 0; it < nit; it++) {
    for (ig = 0; ig < nq_rcfsk; ig++) {
      kint3[it + ig * 2] = wxco2 * kint2[1 + it * 2 + ig * 2 * 2]
                           + (1.0 - wxco2) * kint2[0 + it * 2 + ig * 2 * 2];
    }
  }

  /* Interpolation on T */

  const cs_real_t wt = (t - t_kg_rcfsk[ita[0]]) / (t_kg_rcfsk[ita[1]] - t_kg_rcfsk[ita[0]]);

  for (ig = 0; ig < nq_rcfsk; ig++) {
    phi[ig] = wt * kint3[1 + ig * 2] + (1.0 - wt) * kint3[0 + ig * 2];
  }

  /* Free memory */
  CS_FREE(karray);
  CS_FREE(kint1);
  CS_FREE(kint2);
  CS_FREE(kint3);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Determine the radiation coefficients of the RCFSK model
 *         as well as the corresponding weights. This is applicable
 *         for the case of gas combustion only.
 *
 * \param[in]     pco2        CO2 volume fraction
 * \param[in]     ph2o        H2O volume fraction
 * \param[in]     fvsloc      soot volume fraction
 * \param[in]     teloc       gas temperature
 * \param[out]    kloc        radiation coefficient of the i different gases
 * \param[out]    aloc        weights of the i different gases in cells
 * \param[out]    alocb       weights of the i different gases at boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_rcfsk(const cs_real_t  *restrict pco2,
                      const cs_real_t  *restrict ph2o,
                      const cs_real_t  *restrict fvsloc,
                      const cs_real_t  *restrict teloc,
                      cs_real_t        *restrict kloc,
                      cs_real_t        *restrict aloc,
                      cs_real_t        *restrict alocb)
{
  ipass_rcfsk++;
  if (ipass_rcfsk == 1) { /* Read parameters files */

    FILE       *radfile     = NULL;
    const char *pathdatadir = cs_base_get_pkgdatadir();
    char        filepath[256];

    CS_MALLOC(kgdatabase,
              nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk * nt_rcfsk * nq_rcfsk,
              cs_real_t);
    CS_MALLOC(agdatabase,
              nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk * nt_rcfsk * nq_rcfsk,
              cs_real_t);

    /* Read k-distributions kgdatabase (cm-1) */
    {
      snprintf(filepath, 256, "%s/data/thch/dp_radiat_MFS_RCFSK", pathdatadir);
      radfile = fopen(filepath, "r");

      for (int csoot = 0; csoot < nfvs_rcfsk; csoot++) {
        for (int ch2o = 0; ch2o < nconc_rcfsk; ch2o++) {
          for (int cco2 = 0; cco2 < nconc_rcfsk; cco2++) {
            for (int it = 0; it < nt_rcfsk; it++) {
              for (int ig = 0; ig < nq_rcfsk; ig++) {
                cs_real_t temp[2] = { 0. };
                int       nvalues;
                _line_to_array_rcfsk(radfile, temp, &nvalues);
                assert(nvalues == 2);
                kgdatabase[csoot + ch2o * nfvs_rcfsk
                           + cco2 * nfvs_rcfsk * nconc_rcfsk
                           + it * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                           + ig * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                               * nt_rcfsk]
                  = temp[0];
                agdatabase[csoot + ch2o * nfvs_rcfsk
                           + cco2 * nfvs_rcfsk * nconc_rcfsk
                           + it * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                           + ig * nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk
                               * nt_rcfsk]
                  = temp[1];
              }
            }
          }
        }
      }
      fclose(radfile);
    }

    const int ntot
      = nfvs_rcfsk * nconc_rcfsk * nconc_rcfsk * nt_rcfsk * nq_rcfsk;

    for (int it = 0; it < ntot; it++)
      kgdatabase[it] *= 100.;

  } /* End read parameters file */

  /* Determination of the local absorption coefficient and local function a */

  cs_real_t *kgfsk, *agfsk;
  CS_MALLOC(kgfsk, cs_glob_rad_transfer_params->nwsgg, cs_real_t);
  CS_MALLOC(agfsk, cs_glob_rad_transfer_params->nwsgg, cs_real_t);

  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {

    for (int i = 0; i < nq_rcfsk; i++) {
      kgfsk[i] = 0.;
    }

    const cs_real_t tloc    = cs::max(teloc[iel], 300.0);
    const cs_real_t xh2oloc = cs::max(ph2o[iel], 1e-5);
    const cs_real_t xco2loc = cs::max(pco2[iel], 1e-5);
    const cs_real_t sootloc = cs::min(cs::max(fvsloc[iel], 1e-15), 1e-5);

    _interpolation4d_rcfsk(2, tloc, xco2loc, xh2oloc, sootloc, kgfsk);
    _interpolation4d_rcfsk(1, tloc, xco2loc, xh2oloc, sootloc, agfsk);

    for (int i = 0; i < cs_glob_rad_transfer_params->nwsgg; i++) {
      kloc[i * cs_glob_mesh->n_cells + iel] = kgfsk[i];
      aloc[i * cs_glob_mesh->n_cells + iel] = agfsk[i];
    }
  }

  CS_FREE(kgfsk);
  CS_FREE(agfsk);

  /* Determination of the function a for boundary faces  */

  cs_field_t *f_bound_t = cs_field_by_name_try("boundary_temperature");
  cs_real_t  *tpfsck    = f_bound_t->val;
  cs_real_t  *agb;
  cs_lnum_t  *ifabor = cs_glob_mesh->b_face_cells;

  CS_MALLOC(agb, cs_glob_rad_transfer_params->nwsgg, cs_real_t);

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {

    const cs_lnum_t cell_id = ifabor[ifac];

    const cs_real_t tbord   = cs::max(tpfsck[ifac], 300.0);
    const cs_real_t xh2oloc = cs::max(ph2o[cell_id], 1e-5);
    const cs_real_t xco2loc = cs::max(pco2[cell_id], 1e-5);
    const cs_real_t sootloc = cs::min(cs::max(fvsloc[cell_id], 1e-15), 1.e-5);

    _interpolation4d_rcfsk(1, tbord, xco2loc, xh2oloc, sootloc, agb);

    for (int i = 0; i < cs_glob_rad_transfer_params->nwsgg; i++) {
      alocb[i * cs_glob_mesh->n_b_faces + ifac] = agb[i];
    }
  }

  CS_FREE(agb);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
