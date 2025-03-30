/*============================================================================
 * Combustion model parameters.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"
#include "base/cs_base.h"
#include "base/cs_file.h"
#include "base/cs_dispatch.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "comb/cs_coal.h"
#include "cogz/cs_combustion_bsh.h"
#include "cogz/cs_combustion_gas.h"

#include "pprt/cs_combustion_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_model.cpp
        Combustion  model selection parameters.
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Enthalpy and Cp based on the JANAF band.
 *
 * \param[in]   ncoel        number of elementary constituents
 * \param[in]   ngazem       number of elementary constituents
 * \param[in]   npo          number of interpolation points
 * \param[in]   nomcoel      names of elementary constituants
 * \param[out]  ehcoel       enthalpy for each elementary species
 *                           (for point i and species j, ehcoel[i*ngazem + j])
 * \param[out]  cpcoel       cp for each elementary species
 *                           (for point i and species j, cpcoel[i*ngazem + j])
 * \param[out]  coeff_therm  coefficients for the Burke-Scumann model,
 *                           or null otherwise
 * \param[in]   wmolce       molar mass of each species
 * \param[in]   th           temperature in K
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_enthalpy_and_cp_from_janaf(int           ncoel,
                                         int           ngazem,
                                         int           npo,
                                         const char    nomcoel[][13],
                                         double        ehcoel[],
                                         double        cpcoel[],
                                         double        coeff_therm[][2][5],
                                         const double  wmolce[],
                                         const double  th[])
{
#  undef MAX_ELEMENTARY_COMPONENTS
#  define MAX_ELEMENTARY_COMPONENTS 20
  if (ngazem > MAX_ELEMENTARY_COMPONENTS)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Maximum number of elmentary components handled is %d,\n"
                "but %d requested."),
              __func__, ngazem, MAX_ELEMENTARY_COMPONENTS);

  char nomesp[13] = "";

  int  icoeff[MAX_ELEMENTARY_COMPONENTS];
  char buf[128];

  double  wcoeff[7][2];
  double  coeff[7][2][MAX_ELEMENTARY_COMPONENTS];

  /* Read JANAF thermochemical data
     ------------------------------ */

  const char *datadir = cs_base_get_pkgdatadir();
  const char sub_path[] = "/data/thch/JANAF";

  char *pathdatadir;
  BFT_MALLOC(pathdatadir, strlen(datadir) + strlen(sub_path) + 1, char);
  sprintf(pathdatadir, "%s%s", datadir, sub_path);

#if defined(HAVE_MPI)

  cs_file_t *impjnf = cs_file_open(pathdatadir,
                                   CS_FILE_MODE_READ,
                                   CS_FILE_STDIO_SERIAL,
                                   MPI_INFO_NULL,
                                   MPI_COMM_NULL,
                                   MPI_COMM_NULL);

#else

  cs_file_t *impjnf = cs_file_open(pathdatadir,
                                   CS_FILE_MODE_READ,
                                   CS_FILE_DEFAULT);

#endif

  BFT_FREE(pathdatadir);
  int line = 0;

  /* Initialization */

  for (int ne = 0; ne < ngazem; ne++) {
   icoeff[ne]= 0;
   for (int inicff = 0; inicff < 2; inicff++) {
     for (int injcff = 0; injcff < 7; injcff++)
       coeff[injcff][inicff][ne] = 0;
   }
  }

  for (int ne = 0; ne < ncoel; ne++) {
    for (int nt = 0; nt < npo; nt++) {
      cpcoel[nt*ngazem + ne] = 0;
      ehcoel[nt*ngazem + ne] = 0;
    }
  }

  char *s = cs_file_gets(buf, 127, impjnf, &line);  /* dummy read for 1st line */

  /* Read temperature data */

  double  tlim[3];
  s = cs_file_gets(buf, 127, impjnf, &line);
  sscanf(s, "%lf %lf %lf", tlim, tlim+1, tlim+2);

  /* Read chemical species with partial storage */

  while (true) {

    s = cs_file_gets(buf, 127, impjnf, &line);

    char *p = strtok(s, " \t");
    strncpy(nomesp, p, 12); nomesp[12] = '\0';

    if (strcmp(nomesp, "END") == 0) break;

    s = cs_file_gets(buf, 127, impjnf, &line);
    sscanf(s, "%lf %lf %lf %lf %lf",
           &wcoeff[0][0], &wcoeff[1][0], &wcoeff[2][0],
           &wcoeff[3][0], &wcoeff[4][0]);

    s = cs_file_gets(buf, 127, impjnf, &line);
    sscanf(s, "%lf %lf %lf %lf %lf",
           &wcoeff[5][0], &wcoeff[6][0],
           &wcoeff[0][1], &wcoeff[1][1], &wcoeff[2][1]);

    s = cs_file_gets(buf, 127, impjnf, &line);
    sscanf(s, "%lf %lf %lf %lf",
           &wcoeff[3][1], &wcoeff[4][1],
           &wcoeff[5][1], &wcoeff[6][1]);

    /* We store the coefficients only if the considered species
       is part of the example */

    for (int ne = 0; ne < ncoel; ne++) {
      if (strcmp(nomcoel[ne], nomesp) == 0) {
        icoeff[ne] = 1;
        for (int inicff = 0; inicff < 2; inicff++) {
          for (int injcff = 0; injcff < 7; injcff++)
            coeff[injcff][inicff][ne] = wcoeff[injcff][inicff];
        }
      }
    }

  }  /* end while */

  /* Finish reading if all data has beeen stored */

  impjnf = cs_file_free(impjnf);

  /* Test and possible stop */

  int iok = 0;
  for (int ne = 0; ne < ncoel; ne++) {
    if (icoeff[ne] == 0) {
      iok += 1;
      bft_printf(_("@@ Error: species \'%s\' not in JANAF\n"), nomcoel[ne]);
    }
  }

  if (iok != 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: combustion data input:\n"
                " The parametric file references %d species\n"
                " not in the data/thch/JANAF file\n."
                "\n"
                " Check the setup parameters."),
              __func__, iok);

  /* Compute enthalpies and Cp
     ------------------------- */

  /* perfect gas constants in J/mol/K */

  for (int nt = 0; nt < npo; nt++) {

    /* Determination of the set of coefficients used */
    int ind;
    if (th[nt] > tlim[1])
      ind = 0;
    else
      ind = 1;

    for (int ne = 0; ne < ncoel; ne++) {
      ehcoel[nt*ngazem + ne]  = coeff[5][ind][ne] + coeff[0][ind][ne] * th[nt];
      cpcoel[nt*ngazem + ne]  =                     coeff[0][ind][ne];
      double cth = th[nt];
      double ctc = 1.;

      /* In the JANAF table, coefficients are adimensional (CP/R,H/R) */
      for (int nc = 1; nc < 5; nc++) {
        cth = cth * th[nt];
        ctc = ctc * th[nt];
        ehcoel[nt*ngazem + ne] += coeff[nc][ind][ne] * cth / (double)(nc+1);
        cpcoel[nt*ngazem + ne] += coeff[nc][ind][ne] * ctc;
      }

      /* Compute CP and H for each species */
      ehcoel[nt*ngazem + ne] *= cs_physical_constants_r / wmolce[ne];
      cpcoel[nt*ngazem + ne] *= cs_physical_constants_r / wmolce[ne];
    }

  }

  /* Compute coeff_therm for Burk-Schumann model */

  if (coeff_therm != nullptr) {
    for (int ne = 0; ne < ncoel; ne++) {
      for (int inicff = 0; inicff < 2; inicff++) {
        for (int injcff = 0; injcff < 7; injcff++) {
          coeff_therm[injcff][inicff][ne] = cs_physical_constants_r
                                            * coeff[injcff][inicff][ne]
                                            / wmolce[ne];
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rectangle-Dirac pdf parameters.
 *
 * From P. Plion & A. Escaich
 *
 * \param[in]       n_cells       number of cells
 * \param[out]      indpdf        indicator for pdf integration or mean value
 * \param[out]      tpdf          indicator for pdf shape:
 *                               - 0: Dirac at mean value
 *                               - 1: rectangle
 *                               - 2: Dirac's peak at \f$ f_{min} \f$
 *                               - 3: Dirac's peak at \f$ f_{max} \f$
 *                               - 4: rectangle and 2 Dirac's pics
 * \param[in]       fm            mean mixture fraction at cell centers
 * \param[in, out]  fp2m          mean mixture fraction variance at cell centers
 * \param[in, out]  fmini         mixture fraction low boundary
 * \param[in]       fmaxi         mixture fraction high boundary
 * \param[out]      dirmin        Dirac's peak value at \f$ f_{min} \f$
 * \param[out]      dirmax        Dirac's peak value at \f$ f_{max} \f$
 * \param[out]      fdeb          abscissa of rectangle low boundary
 * \param[out]      ffin          abscissa of rectangle high boundary
 * \param[out]      hrec          rectangle height
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_dirac_pdf(cs_lnum_t         n_cells,
                        int               indpdf[],
                        cs_real_t         tpdf[],
                        const cs_real_t   fm[],
                        cs_real_t         fp2m[],
                        const cs_real_t   fmini[],
                        const cs_real_t   fmaxi[],
                        cs_real_t         dirmin[],
                        cs_real_t         dirmax[],
                        cs_real_t         fdeb[],
                        cs_real_t         ffin[],
                        cs_real_t         hrec[])
{
  /* Initialization
   * -------------- */

  bool log_active = cs_log_default_is_active();

  // Parameter relative to variance
  cs_real_t t1 = 1.e-08;
  // Parameter relative to mean
  cs_real_t t2 = 5.e-07;

  cs_real_t epzero = cs_math_epzero;

  cs_host_context ctx;

  /* Preliminary computations
   * ------------------------ */

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

    tpdf  [c_id] = 0.0;
    dirmin[c_id] = 0.0;
    dirmax[c_id] = 0.0;
    fdeb  [c_id] = 0.0;
    ffin  [c_id] = 0.0;
    hrec  [c_id] = 0.0;

    // Change parameters T1 and T2 to acccount for fact that
    // FMINI < FM < FMAXI
    cs_real_t delta_t = fmaxi[c_id]-fmini[c_id];
    cs_real_t t1mod = t1 * (delta_t * delta_t);
    cs_real_t t2mod = t2 * delta_t;

    if (   (fp2m[c_id] > t1mod)
        && (fm[c_id] >= (fmini[c_id] + t2mod))
        && (fm[c_id] <= (fmaxi[c_id] - t2mod)))
      indpdf[c_id] = 1;
    else
      indpdf[c_id] = 0;

  });

  // Clip variance

  cs_real_t fp2mmin1 = HUGE_VALF, fp2mmax1 = -HUGE_VALF;
  cs_gnum_t nfp2 = 0;

  if (log_active) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      fp2mmin1 = cs::min(fp2mmin1, fp2m[c_id]);
      fp2mmax1 = cs::max(fp2mmax1, fp2m[c_id]);
    }
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t fp2max = (fmaxi[c_id]-fm[c_id]) * (fm[c_id]-fmini[c_id]);
    if (fp2m[c_id] > fp2max+1.e-20) {
      fp2m[c_id] = fp2max;
      nfp2 += 1;
    }
  }

  if (log_active) {
    cs_parall_counter(&nfp2, 1);

    cs_real_t fp2mmin2 = HUGE_VALF, fp2mmax2 = -HUGE_VALF;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      fp2mmin2 = cs::min(fp2mmin2, fp2m[c_id]);
      fp2mmax2 = cs::max(fp2mmax2, fp2m[c_id]);
    }
    cs_real_t rval[4] = {-fp2mmin1, fp2mmax1, -fp2mmin2, fp2mmax2};
    cs_parall_max(4, CS_REAL_TYPE, rval);
    fp2mmin1 = -rval[0]; fp2mmax1 = rval[1];
    fp2mmin2 = -rval[2]; fp2mmax2 = rval[3];

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("  pppdfr: variance clipping points: %llu\n"),
       (unsigned long long)nfp2);

    if (nfp2 > 0) {
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("     Variance before clipping min and max: %g %g\n"
           "     Variance after  clipping min and max: %g %g\n"),
         fp2mmin1, fp2mmax1, fp2mmin2, fp2mmax2);
    }
  }

  /* Compute parameters of probability density function
   * -------------------------------------------------- */

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {

    if (indpdf[c_id] == 1) {
      cs_real_t f_mid = (fmini[c_id] + fmaxi[c_id])*0.5;
      cs_real_t fm_m_fmini = fm[c_id] - fmini[c_id];

      if (   (   (fm[c_id] <= f_mid)
              && (fp2m[c_id] <= cs_math_pow2(fm_m_fmini)/3.))
          || (   (fm[c_id] > f_mid)
              && (fp2m[c_id] <= cs_math_pow2(fmaxi[c_id] -fm[c_id])/3.))) {
        // Rectangle only

        tpdf[c_id] = 1.0;

        hrec[c_id]   = sqrt(3.0*fp2m[c_id]);
        dirmin[c_id] = 0.0;
        dirmax[c_id] = 0.0;
        fdeb[c_id] = fm[c_id] - hrec[c_id];
        ffin[c_id] = fm[c_id] + hrec[c_id];
      }
      else if (   (fm[c_id] <= f_mid)
               && (fp2m[c_id] <= (  fm_m_fmini
                                  * (2.0*fmaxi[c_id] + fmini[c_id] - 3.0*fm[c_id])
                                  / 3.0))) {
        // Rectangle and Dirac at FMINI

        tpdf[c_id] = 2.0;

        fdeb[c_id]   = fmini[c_id];
        dirmax[c_id] = 0.0;
        ffin[c_id]   =   fmini[c_id]
                       + 1.5*(cs_math_pow2(fm_m_fmini) + fp2m[c_id])
                            /(fm_m_fmini);
        dirmin[c_id] =   (3.0*fp2m[c_id] - cs_math_pow2(fm_m_fmini))
                       / (3.*(cs_math_pow2(fm_m_fmini) + fp2m[c_id]));
      }

      else if (   (fm[c_id]  > f_mid)
               && (fp2m[c_id] <= (  (fmaxi[c_id] - fm[c_id])
                                  * (3.0*fm[c_id]-fmaxi[c_id]-2.0*fmini[c_id])
                                  / 3.0))) {

        // Rectangle and Dirac at FMAXI
        // (correct: HI/81/02/03/A has an error p 12).

        tpdf[c_id] = 3.0;

        ffin[c_id]   = fmaxi[c_id];
        dirmin[c_id] = 0.;
        fdeb[c_id]   =   fmini[c_id]
                       + 3.0*(  (cs_math_pow2(fm_m_fmini) + fp2m[c_id])
                              +  cs_math_pow2(fmaxi[c_id] - fmini[c_id])
                          - 4.0 * fm_m_fmini * (fmaxi[c_id] - fmini[c_id]))
                         / (2.0*(fm[c_id] - fmaxi[c_id]));
        dirmax[c_id] =   (3.0*fp2m[c_id] - cs_math_pow2(fm[c_id] - fmaxi[c_id]))
                       / (3.0*(cs_math_pow2(fm[c_id] - fmaxi[c_id]) +fp2m[c_id]));
      }
      else {
        // Rectangle and 2 Diracs

        tpdf  [c_id] = 4.0;

        fdeb[c_id]   = fmini[c_id];
        ffin[c_id]   = fmaxi[c_id];
        dirmax[c_id] =   3.0*(cs_math_pow2(fm_m_fmini) +fp2m[c_id])
                       / cs_math_pow2(fmaxi[c_id] - fmini[c_id])
                       -2.0 * (fm_m_fmini)
                            / (fmaxi[c_id] - fmini[c_id]);
        dirmin[c_id] =  dirmax[c_id] + 1.0
                       - 2.0*(fm[c_id]-fmini[c_id])/(fmaxi[c_id]-fmini[c_id]);
      }

      if (fabs(ffin[c_id] - fdeb[c_id]) > epzero) {
        hrec[c_id] = (1.0-dirmin[c_id]-dirmax[c_id]) / (ffin[c_id]-fdeb[c_id]);
      }
      else {
        cs_real_t t3 = sqrt(3.*t1*cs_math_pow2(fmaxi[c_id]-fmini[c_id]));
        fdeb[c_id] = fmin(fmaxi[c_id], fmax(fmini[c_id], fm[c_id] - t3));
        ffin[c_id] = fmin(fmaxi[c_id], fmax(fmini[c_id], fm[c_id] + t3));
        if (fabs(ffin[c_id] - fdeb[c_id]) > epzero)
          hrec[c_id] = (1.0-dirmin[c_id]-dirmax[c_id]) / (ffin[c_id] - fdeb[c_id]);
        else
          hrec[c_id] = 0.0;
      }
    }
    else  {
      tpdf[c_id] = 0.;

      dirmin[c_id] = 0.;
      dirmax[c_id] = 0.;
      fdeb[c_id]   = 0.;
      ffin[c_id]   = 0.;
      hrec[c_id]   = 0.;
    }

  });

  // Check: if Hrec <= 0 we pass without the PDF

  cs_gnum_t nbspdf = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (hrec[c_id] <= epzero && indpdf[c_id] == 1) {
      indpdf[c_id] = 0;
      nbspdf += 1;
    }
  }

  if (log_active) {
    cs_parall_counter(&nbspdf, 1);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("  pppdfr: switch off PDF %llu\n\n"),
       (unsigned long long)nbspdf);

    /* Logging
     * ------- */

    cs_gnum_t n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0, n6 = n_cells;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (indpdf[c_id] == 1) {
        n1 += 1;
        if (dirmin[c_id] > epzero && dirmax[c_id] < epzero)
          n2 += 1;
        else if (dirmin[c_id] < epzero && dirmax[c_id] > epzero)
          n3 += 1;
        else if (dirmin[c_id] > epzero && dirmax[c_id] > epzero)
          n4 += 1;
        else if (dirmin[c_id] < epzero && dirmax[c_id] < epzero)
          n5 += 1;
      }
    }

    cs_parall_sum_scalars(n1, n2, n3, n4, n5, n6);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("Rectangle PDF - Dirac peaks\n"
         "Mean, variance of transported tracer\n"
         "Number of turbulent points (using the PDFs)   = %lu\n"
         "Number of computation points                  = %lu\n"),
       (unsigned long)n1, (unsigned long)n6);

    cs_log_printf
      (CS_LOG_DEFAULT,
       _(" Nb points with rectangle PDF without Dirac              = %lu\n"
         " - - - - - - - - - -- - - - and Dirac in FMINI           = %lu\n"
         " - - - - - - - - - -- - - - - - - - - -  FMAXI           = %lu\n"
         " - - - - - - - - - - - - - - - Diracs in FMINI and FMAXI = %lu\n"),
       (unsigned long)n5, (unsigned long)n2,
       (unsigned long)n3, (unsigned long)n4);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
