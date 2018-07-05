/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_thermal_model.h"
#include "cs_physical_constants.h"

#include "cs_prototypes.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_adf_models.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/* \file cs_rad_transfer_adf_models.c.f90 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local static variables
 *============================================================================*/

/* Common statics */
static cs_int_t ntsto = 0;
static cs_int_t ipass = 0;
static cs_real_t *tsto = NULL;
static cs_real_t *asto = NULL;
static cs_real_t *ksto2 = NULL;

/* ADF08 model specifics */
static cs_int_t nysto = 0;
static cs_real_t *ysto = NULL;

/* ADF50 model specifics */
static cs_int_t nxh2osto = 0;
static cs_real_t *xh2osto = NULL;
static cs_real_t *ksto1 = NULL;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert a line of values in a file into an array of cs_real_t and
 *         returns the number of values read
 *
 * Parameters
 * \param[in]  radfile     pointer to file for reading values
 * \param[in]  values      array for values storage
 * \param[in]  nvalues     number of values read
 */
/*----------------------------------------------------------------------------*/

static inline void
_line_to_array(FILE      *radfile,
               cs_real_t  values[],
               cs_int_t  *nvalues)
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
 * \brief Determine the radiation coefficients of the ADF 08 model
 *        as well as the corresponding weights.
 *
 * \param[in]     pco2       CO2 volume fraction
 * \param[in]     ph2o       H2O volume fraction
 * \param[in]     teloc      gas temperature
 * \param[out]    kloc       radiation coefficient of the i different gases
 * \param[out]    aloc       weights of the i different gases in cells
 * \param[out]    alocb      weights of the i different gases at boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_adf08(const cs_real_t  pco2[],
                      const cs_real_t  ph2o[],
                      const cs_real_t  teloc[],
                      cs_real_t        kloc[],
                      cs_real_t        aloc[],
                      cs_real_t        alocb[])
{
  cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;
  cs_lnum_t ncells = cs_glob_mesh->n_cells;

  int nwsgg = cs_glob_rad_transfer_params->nwsgg;
  cs_real_t  tkelvi = 273.15;

  cs_real_t  tref, xh2oref, rt, rx;

  /* Memory allocation and initialization */

  cs_real_t *y, *tpaadf;
  cs_field_t *f_b_temp = cs_field_by_name_try("boundary_temperature");
  BFT_MALLOC(y, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);

  if (cs_glob_thermal_model->itpscl == 2) {

    BFT_MALLOC(tpaadf, nfabor, cs_real_t );
    for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++)
      tpaadf[ifac] = f_b_temp->val[ifac] + tkelvi;

  }
  else
    tpaadf = f_b_temp->val;

  /* Absorption coefficient of gas mix (m-1)
     --------------------------------------- */

  ipass++;

  /* The ADF data base is read only once during the very first iteration */

  if (ipass == 1) {

    const char *pathdatadir = cs_base_get_pkgdatadir();
    char filepath[256];
    snprintf(filepath, 256, "%s/data/thch/dp_radiat_ADF8", pathdatadir);
    FILE *radfile = fopen(filepath, "r");

    char line[256];
    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    fscanf(radfile, "%d", &ntsto);

    /* Number of tabulated gas phase temperatures    */
    BFT_MALLOC(tsto, ntsto, cs_real_t);

    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    /* Tabulated gas phase temperatures    */
    int itsto = 0;
    while (itsto < ntsto - 1) {
      cs_real_t temp[20] = {0.};
      int nvalues;
      _line_to_array(radfile, temp, &nvalues);
      for (int i = 0; i < nvalues; i++) {
        tsto[itsto] = temp[i];
        itsto++;
      }
    }

    fgets(line, 256, radfile);

    /* Number of tabulated ratios ph2o/pco2     */
    fscanf(radfile, "%d", &nysto);

    BFT_MALLOC(ysto, nysto, cs_real_t);

    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    /* Tabulated ratios ph2o/pco2     */
    int iysto = 0;
    while (iysto < nysto - 1) {
      cs_real_t temp[20] = {0.};
      int nvalues;
      _line_to_array(radfile, temp, &nvalues);
      for (int i = 0; i < nvalues; i++) {
        ysto[iysto] = temp[i];
        iysto++;
      }
    }

    fgets(line, 256, radfile);

    /* Reference temperature and h2o volume fraction */
    fscanf(radfile, "%lf %lf", &tref, &xh2oref);
    fgets(line, 256, radfile);

    BFT_MALLOC(asto,  nwsgg * nysto * ntsto, cs_real_t);
    BFT_MALLOC(ksto2, nwsgg * nysto * ntsto, cs_real_t);

    /*    READING THE PARAMETERS */
    fgets(line, 256, radfile);
    for (int i = 0; i < nwsgg; i++) {

      fgets(line, 256, radfile);

      for (int j = 0; j < ntsto; j++) {
        cs_real_t *temp;
        BFT_MALLOC(temp, 2*nysto, cs_real_t);
        int nvalues;
        _line_to_array(radfile, temp, &nvalues);
        assert(nvalues == nysto * 2);
        for (int k = 0; k < nysto; k++) {
          ksto2[i + k * nwsgg + j * nysto * nwsgg] = temp[k];
          /* ksto2(i,k,j): Radiation coefficient of the i-th grey gas, */
          /*               the k-th ratio of ph2o/pco2, and */
          /*               the j-th tabulated temperature   */
          asto[i + k * nwsgg + j * nysto * nwsgg] = temp[k + nysto];
          /* asto2(i,k,j): Weight of the i-th grey gas,     */
          /*               the k-th ratio of ph2o/pco2, and */
          /*               the j-th tabulated temperature   */
        }
        BFT_FREE(temp);
      }
    }
  }

  int it, ix, l;
  for (cs_lnum_t iel = 0; iel < ncells; iel++) {
    if (pco2[iel] > 0.0)
      y[iel] = ph2o[iel] / pco2[iel];
    else
      y[iel] = ysto[nysto - 1];

    /* Interpolation temperature */
    if (teloc[iel] <= tsto[0]) {
      rt = 0.0;
      it = 0;
    }
    else if (teloc[iel] >= tsto[ntsto - 1]) {
      rt = 1.0;
      it = ntsto - 2;
    }
    else {
      l = 0;
      while (teloc[iel] > tsto[l])
        l++;
      it = l - 1;
      rt = (teloc[iel] - tsto[it]) / (tsto[it + 1] - tsto[it]);
    }

    /* Interpolation H2O-molefraction */
    if (y[iel] <= ysto[0]) {
      rx = 0.0;
      ix = 0;
    }
    else if (y[iel] >= ysto[nysto - 1]) {
      rx = 1.0;
      ix = nysto - 2;
    }
    else {
      l = 0;
      while (y[iel] > ysto[l])
        l++;
      ix = l - 1;
      rx = (y[iel] - ysto[ix]) / (ysto[ix + 1] - ysto[ix]);
    }

    /* Absortion Coefficient */

    for (cs_int_t i = 0; i < nwsgg; i++) {
      cs_real_t kmloc =    (1.0 - rt) * (1.0 - rx)
                         * ksto2[i + ix * nwsgg + it * nysto * nwsgg]
                       +   (1.0 - rt) * (rx)
                         * ksto2[i + (ix + 1) * nwsgg + it * nysto * nwsgg]
                       +   rt * (1.0 - rx)
                         * ksto2[i + ix * nysto + (it + 1) * nysto * nwsgg]
                       + rt * rx
                         * ksto2[i + (ix + 1) * nwsgg + (it + 1) * nysto * nwsgg];

      kloc[iel + i * ncells] =  pco2[iel] * kmloc * 100.0
                              * (cs_glob_fluid_properties->p0 / 100000.0);

      /* Local radiation coefficient of the i-th grey gas    */
      aloc[iel + i * ncells]
        =   (1.0 - rt) * (1.0 - rx) * asto[i + ix * nwsgg + it * nysto * nwsgg]
          + (1.0 - rt) * (rx) * asto[i + (ix + 1) * nwsgg + it * nysto * nwsgg]
          + rt * (1.0 - rx) * asto[i + ix * nysto + (it + 1) * nysto * nwsgg]
          + rt * rx * asto[i + (ix + 1) * nwsgg + (it + 1) * nysto * nwsgg];
      /* Local weight of the i-th grey gas   */
    }

  }

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {
    cs_lnum_t iel = cs_glob_mesh->b_face_cells[ifac];
    if (pco2[iel] > 0.0)
      y[iel] = ph2o[iel] / pco2[iel];
    else
      y[iel] = ysto[nysto - 1];

    /* Interpolation temperature */
    if (tpaadf[ifac] <= tsto[0]) {
      rt = 0.0;
      it = 0;
    }
    else if (tpaadf[ifac] >= tsto[ntsto - 1]) {
      rt = 1.0;
      it = ntsto - 2;
    }
    else {
      l = 0;
      while (tpaadf[ifac] > tsto[l])
        l++;
      it = l - 1;
      rt = (tpaadf[ifac] - tsto[it]) / (tsto[it + 1] - tsto[it]);
    }

    /* Interpolation H2O-molefraction */
    if (y[iel] <= ysto[0]) {
      rx = 0.0;
      ix = 0;
    }
    else if (y[iel] >= ysto[nysto - 1]) {
      rx = 1.0;
      ix = nysto - 2;
    }
    else {
      l = 0;
      while (y[iel] > ysto[l])
        l++;
      ix = l - 1;
      rx = (y[iel] - ysto[ix]) / (ysto[ix + 1] - ysto[ix]);
    }

    /* Absortion Coefficient     */

    for (int i = 0; i < nwsgg; i++)
      alocb[ifac + i * nfabor]
        =   (1.0 - rt) * (1.0 - rx) * asto[i + ix * nwsgg + it * nysto * nwsgg]
          + (1.0 - rt) * rx * asto[i + (ix + 1) * nwsgg + it * nysto * nwsgg]
          + rt * (1.0 - rx) * asto[i + ix * nwsgg + (it + 1) * nysto * nwsgg]
          + rt * rx * asto[i + (ix + 1) * nwsgg + (it + 1) * nysto * nwsgg];

    /* Local weight of the i-th grey gas   */

  }

  BFT_FREE(tpaadf);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine the radiation coefficients of the ADF 50 model
 *        as well as the corresponding weights.
 *
 * \param[in]     pco2       CO2 volume fraction
 * \param[in]     ph2o       H2O volume fraction
 * \param[in]     teloc      gas temperature
 * \param[out]    kloc       radiation coefficient of the i different gases
 * \param[out]    aloc       weights of the i different gases in cells
 * \param[out]    alocb      weights of the i different gases at boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_adf50(const cs_real_t  pco2[],
                      const cs_real_t  ph2o[],
                      const cs_real_t  teloc[],
                      cs_real_t        kloc[],
                      cs_real_t        aloc[],
                      cs_real_t        alocb[])
{
  cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;
  cs_lnum_t ncells = cs_glob_mesh->n_cells;

  int nwsgg = cs_glob_rad_transfer_params->nwsgg;
  cs_real_t  tkelvi = 273.15;

  cs_real_t  tref, xh2oref, rt, rx;

  /* Memory allocation and initialization */

  cs_real_t *tpaadf;
  cs_field_t *f_b_temp = cs_field_by_name_try("boundary_temperature");

  if (cs_glob_thermal_model->itpscl == 2) {

    BFT_MALLOC(tpaadf, nfabor, cs_real_t );
    for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++)
      tpaadf[ifac] = f_b_temp->val[ifac] + tkelvi;

  }
  else
    tpaadf = f_b_temp->val;

  /* Absorption coefficient of gas mix (m-1)
     --------------------------------------- */

  ipass++;

  /* The ADF data base is read only once during the very first iteration */

  if (ipass == 1) {

    const char *pathdatadir = cs_base_get_pkgdatadir();
    char filepath[256];
    snprintf(filepath, 256, "%s/data/thch/dp_radiat_ADF50", pathdatadir);
    FILE *radfile = fopen(filepath, "r");

    char line[256];
    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    fscanf(radfile, "%d", &ntsto);

    /* Number of tabulated gas phase temperatures    */
    BFT_MALLOC(tsto, ntsto, cs_real_t);

    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    /* Tabulated gas phase temperatures    */
    int itsto = 0;
    while (itsto < ntsto - 1) {
      cs_real_t temp[20] = {0.};
      int nvalues;
      _line_to_array(radfile, temp, &nvalues);
      for (int i = 0; i < nvalues; i++) {
        tsto[itsto] = temp[i];
        itsto++;
      }
    }

    fgets(line, 256, radfile);

    /* Number of tabulated h2o volume fractions  */
    fscanf(radfile, "%d", &nxh2osto);

    BFT_MALLOC(xh2osto, nxh2osto, cs_real_t);

    fgets(line, 256, radfile);
    fgets(line, 256, radfile);
    /* Tabulated ratios ph2o/pco2     */
    int ixh2osto = 0;
    while (ixh2osto < nxh2osto - 1) {
      cs_real_t temp[20] = {0.};
      int nvalues;
      _line_to_array(radfile, temp, &nvalues);
      for (int i = 0; i < nvalues; i++) {
        xh2osto[ixh2osto] = temp[i];
        ixh2osto++;
      }
    }

    fgets(line, 256, radfile);

    /* Reference temperature and h2o volume fraction */
    fscanf(radfile, "%lf %lf", &tref, &xh2oref);
    fgets(line, 256, radfile);

    BFT_MALLOC(asto,  nwsgg * nxh2osto * ntsto, cs_real_t);
    BFT_MALLOC(ksto1, nwsgg * nxh2osto, cs_real_t);
    BFT_MALLOC(ksto2, nwsgg * nxh2osto * ntsto, cs_real_t);

    /*    READING THE PARAMETERS */
    fgets(line, 256, radfile);
    for (int i = 0; i < nwsgg; i++) {

      fgets(line, 256, radfile);

      for (int j = 0; j < ntsto; j++) {
        cs_real_t *temp;
        BFT_MALLOC(temp, 2*nxh2osto, cs_real_t);
        int nvalues;
        _line_to_array(radfile, temp, &nvalues);
        assert(nvalues == nxh2osto * 2);

        ksto1[i + j * nwsgg] = temp[0];
        /*ksto1(i,j): Radiation coefficient CO2 of the i-th grey gas,
         *            and the j-th tabulated temperature.*/
        for (int k = 0; k < nxh2osto; k++) {
          ksto2[i + k * nwsgg + j * nxh2osto * nwsgg] = temp[k + 1];
          /* ksto2(i,k,j): Radiation coefficient of h2o of the i-th grey gas,
           *               the k-th h2o volume fraction, and
           *               the j-th tabulated temperature.*/
          asto[i + k * nwsgg + j * nxh2osto * nwsgg] = temp[k + nxh2osto + 1];
          /* asto2(i,k,j): Weight of the i-th grey gas,
           *               the k-th h2o volume fraction, and
           *               the j-th tabulated temperature   */
        }
        BFT_FREE(temp);
      }
    }
  }

  int it, ix, l;
  for (cs_lnum_t iel = 0; iel < ncells; iel++) {

    /* Interpolation temperature */
    if (teloc[iel] <= tsto[0]) {
      rt = 0.0;
      it = 0;
    }
    else if (teloc[iel] >= tsto[ntsto - 1]) {
      rt = 1.0;
      it = ntsto - 2;
    }
    else {
      l = 0;
      while (teloc[iel] > tsto[l])
        l++;
      it = l - 1;
      rt = (teloc[iel] - tsto[it]) / (tsto[it + 1] - tsto[it]);
    }

    /* Interpolation H2O-molefraction */
    if (ph2o[iel] <= xh2osto[0]) {
      rx = 0.0;
      ix = 0;
    }
    else if (ph2o[iel] >= xh2osto[nxh2osto - 1]) {
      rx = 1.0;
      ix = nxh2osto - 2;
    }
    else {
      l = 0;
      while (ph2o[iel] > xh2osto[l])
        l++;
      ix = l - 1;
      rx = (ph2o[iel] - xh2osto[ix]) / (xh2osto[ix + 1] - xh2osto[ix]);
    }

    /* Absortion Coefficient */
    for (cs_int_t i = 0; i < nwsgg; i++) {
      cs_real_t kco2loc =  ksto1[i + it * nwsgg]
                         + rt * (  ksto1[i + (it+1) * nwsgg]
                                 - ksto1[i + it * nwsgg]);
      cs_real_t kh2oloc =   (1.0 - rt) * (1.0 - rx)
                          * ksto2[i + ix * nwsgg + it * nxh2osto * nwsgg]
                         +  (1.0 - rt) * (rx)
                          * ksto2[i + (ix + 1) * nwsgg + it * nxh2osto * nwsgg]
                         +  rt * (1.0 - rx)
                          * ksto2[i + ix * nxh2osto + (it + 1) * nxh2osto * nwsgg]
                         +  rt * rx
                          * ksto2[i + (ix + 1) * nwsgg + (it + 1) * nxh2osto * nwsgg];

      kloc[iel + i * ncells] =  (pco2[iel] * kco2loc + ph2o[iel] * kh2oloc) * 100.0
                              * (cs_glob_fluid_properties->p0 / 100000.0);

      /* Local radiation coefficient of the i-th grey gas */
      aloc[iel + i * ncells] =   (1.0 - rt) * (1.0 - rx)
                               * asto[i + ix * nwsgg + it * nxh2osto * nwsgg]
                              +  (1.0 - rt) * (rx)
                               * asto[i + (ix + 1) * nwsgg + it * nxh2osto * nwsgg]
                              +  rt * (1.0 - rx)
                               * asto[i + ix * nxh2osto + (it + 1) * nxh2osto * nwsgg]
                              +  rt * rx
                               * asto[i + (ix + 1) * nwsgg + (it + 1) * nxh2osto * nwsgg];
                              /* Local weight of the i-th grey gas */
    }

  }

  for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++) {
    cs_lnum_t iel = cs_glob_mesh->b_face_cells[ifac];

    /* Interpolation temperature */
    if (tpaadf[ifac] <= tsto[0]) {
      rt = 0.0;
      it = 0;
    }
    else if (tpaadf[ifac] >= tsto[ntsto - 1]) {
      rt = 1.0;
      it = ntsto - 2;
    }
    else {
      l = 0;
      while (tpaadf[ifac] > tsto[l])
        l++;
      it = l - 1;
      rt = (tpaadf[ifac] - tsto[it]) / (tsto[it + 1] - tsto[it]);
    }

    /* Interpolation H2O-molefraction */
    if (ph2o[iel] <= xh2osto[0]) {
      rx = 0.0;
      ix = 0;
    }
    else if (ph2o[iel] >= xh2osto[nxh2osto - 1]) {
      rx = 1.0;
      ix = nxh2osto - 2;
    }
    else {
      l = 0;
      while (ph2o[iel] > xh2osto[l])
        l++;
      ix = l - 1;
      rx = (ph2o[iel] - xh2osto[ix]) / (xh2osto[ix + 1] - xh2osto[ix]);
    }

    /* Absortion Coefficient     */

    for (int i = 0; i < nwsgg; i++)
      alocb[ifac + i * nfabor] =   (1.0 - rt) * (1.0 - rx)
                                  * asto[i + ix * nwsgg + it * nxh2osto * nwsgg]
                                 +  (1.0 - rt) * rx
                                  * asto[i + (ix + 1) * nwsgg + it * nxh2osto * nwsgg]
                                 +  rt * (1.0 - rx)
                                  * asto[i + ix * nwsgg + (it + 1) * nxh2osto * nwsgg]
                                 +  rt * rx
                                  * asto[i + (ix + 1) * nwsgg + (it + 1) * nxh2osto * nwsgg];
                                 /* Local weight of the i-th grey gas   */

  }

  BFT_FREE(tpaadf);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
