/*============================================================================
 * Combustion model parameters.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"
#include "cs_base.h"
#include "cs_file.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal.h"
#include "cs_combustion_bsh.h"
#include "cs_combustion_gas.h"

#include "cs_combustion_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_model.c
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
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_coal_model_map(void);

void
cs_f_ppthch_get_pointers(int     **ngaze,
                         int     **ngazg,
                         int     **nato,
                         int     **nrgaz,
                         int     **iic,
                         int     **iico2,
                         int     **iio2,
                         int     **npo,
                         double  **wmole,
                         double  **wmolg,
                         double  **wmolat,
                         double  **xco2,
                         double  **xh2o,
                         double  **fs,
                         double  **th,
                         double  **cpgazg);

void
cs_f_coincl_get_pointers(int     **isoot,
                         int     **ngazfl,
                         int     **nki,
                         int     **nxr,
                         int     **nzm,
                         int     **nzvar,
                         int     **nlibvar,
                         int     **ikimid,
                         int     **mode_fp2m,
                         bool    **use_janaf,
                         double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot,
                         double  **lsp_fuel,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy,
                         double  **fmin_lwc,
                         double  **fmax_lwc,
                         double  **hmin_lwc,
                         double  **hmax_lwc);

void
cs_f_ppcpfu_get_pointers(double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2);

void
cs_f_combustion_model_get_pointers(double  **srrom);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical model flags.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ngaze  --> pointer to number of elementary species
 *   ngazg  --> pointer to number of global species
 *   nato   --> pointer to number of atomic species
 *   iic    --> pointer to rank of C in gas composition
 *   wmole  --> pointer to molar mass of elementary gas components
 *   wmolg  --> pointer to molar mass of global species
 *   wmolat --> pointer to molar mass of atomic species
 *   xco2   --> pointer to molar coefficient of co2
 *   xh2o   --> pointer to molar coefficient of h2o
 *   fs     --> pointer to mixing rate at the stoichiometry
 *----------------------------------------------------------------------------*/

void
cs_f_ppthch_get_pointers(int     **ngaze,
                         int     **ngazg,
                         int     **nato,
                         int     **nrgaz,
                         int     **iic,
                         int     **iico2,
                         int     **iio2,
                         int     **npo,
                         double  **wmole,
                         double  **wmolg,
                         double  **wmolat,
                         double  **xco2,
                         double  **xh2o,
                         double  **fs,
                         double  **th,
                         double  **cpgazg)
{
  *npo = NULL;
  *wmolg  = NULL;
  *th = NULL;
  *cpgazg = NULL;

  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *ngaze  = &(cm->n_gas_el_comp);
    *ngazg  = &(cm->n_gas_species);
    *nato   = &(cm->n_atomic_species);
    *nrgaz  = &(cm->n_reactions);
    *npo    = &(cm->n_tab_points);
    *iic    = &(cm->iic);
    *iio2   = &(cm->iio2);
    *iico2  = &(cm->iico2);
    *wmole  = cm->wmole;
    *wmolg  = cm->wmolg;
    *wmolat = cm->wmolat;
    *xco2   = &(cm->xco2);
    *xh2o   = &(cm->xh2o);
    *fs     = cm->fs;
    *th     = cm->th;
    *cpgazg = (double *)cm->cpgazg;

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *ngaze  = &(cm->n_gas_el_comp);
    *ngazg  = &(cm->n_gas_species);
    *nato   = &(cm->n_atomic_species);
    *nrgaz  = &(cm->n_reactions);

    *wmole  = cm->wmole;
    *wmolat = cm->wmolat;
    *xco2   = &(cm->xco2);
    *xh2o   = &(cm->xh2o);

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical model (coincl).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   coefeg --> pointer to conversion coefficients
 *   compog --> pointer to conversion coefficients
 *----------------------------------------------------------------------------*/

void
cs_f_coincl_get_pointers(int     **isoot,
                         int     **ngazfl,
                         int     **nki,
                         int     **nxr,
                         int     **nzm,
                         int     **nzvar,
                         int     **nlibvar,
                         int     **ikimid,
                         int     **mode_fp2m,
                         bool    **use_janaf,
                         double  **coefeg,
                         double  **compog,
                         double  **xsoot,
                         double  **rosoot,
                         double  **lsp_fuel,
                         double  **hinfue,
                         double  **hinoxy,
                         double  **pcigas,
                         double  **tinfue,
                         double  **tinoxy,
                         double  **fmin_lwc,
                         double  **fmax_lwc,
                         double  **hmin_lwc,
                         double  **hmax_lwc)
{
  *isoot  = NULL;
  *ngazfl = NULL;
  *nki = NULL;
  *nxr = NULL;
  *nzm = NULL;
  *nzvar = NULL;
  *nlibvar = NULL;
  *ikimid = NULL;
  *mode_fp2m = NULL;
  *use_janaf = NULL;
  *coefeg = NULL;
  *compog = NULL;
  *xsoot  = NULL;
  *rosoot = NULL;
  *lsp_fuel = NULL;
  *hinfue = NULL;
  *tinfue = NULL;
  *hinoxy = NULL;
  *tinoxy = NULL;
  *pcigas = NULL;
  *fmin_lwc = NULL;
  *fmax_lwc = NULL;
  *hmin_lwc = NULL;
  *hmax_lwc = NULL;

  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *isoot  = &(cm->isoot);
    *ngazfl = &(cm->ngazfl);
    *nki = &(cm->nki);
    *nxr = &(cm->nxr);
    *nzm = &(cm->nzm);
    *nzvar = &(cm->nzvar);
    *nlibvar = &(cm->nlibvar);
    *ikimid = &(cm->ikimid);
    *mode_fp2m = &(cm->mode_fp2m);
    *use_janaf = &(cm->use_janaf);
    *coefeg = &(cm->coefeg[0][0]);
    *compog = &(cm->compog[0][0]);
    *xsoot  = &(cm->xsoot);
    *rosoot = &(cm->rosoot);
    *lsp_fuel = &(cm->lsp_fuel);
    *hinfue = &(cm->hinfue);
    *tinfue = &(cm->tinfue);
    *hinoxy = &(cm->hinoxy);
    *tinoxy = &(cm->tinoxy);
    *pcigas = &(cm->pcigas);
    *fmin_lwc = &(cm->fmin_lwc);
    *fmax_lwc = &(cm->fmax_lwc);
    *hmin_lwc = &(cm->hmin_lwc);
    *hmax_lwc = &(cm->hmax_lwc);

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *pcigas = &(cm->pcigas);

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to members of combustion model (ppcpfu).
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ieqco2 --> pointer to cm->ieqco2
 *   oxyo2  --> pointer to cm->oxyo2
 *   oxyn2  --> pointer to cm->oxyn2
 *   oxyh2o --> pointer to cm->oxyh2o
 *   oxyco2 --> pointer to cm->oxyco2
 *----------------------------------------------------------------------------*/

void
cs_f_ppcpfu_get_pointers(double  **oxyo2,
                         double  **oxyn2,
                         double  **oxyh2o,
                         double  **oxyco2)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *oxyo2 = NULL;
    *oxyn2 =  cm->oxyn2;
    *oxyh2o = cm->oxyh2o;
    *oxyco2 = cm->oxyco2;

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *oxyo2 =  cm->oxyo2;
    *oxyn2 =  cm->oxyn2;
    *oxyh2o = cm->oxyh2o;
    *oxyco2 = cm->oxyco2;

  }
}

/*----------------------------------------------------------------------------
 * Get pointers to generic physical model pointers (pincl)
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_combustion_model_get_pointers(double  **srrom)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    *srrom = &(cm->srrom);

  }
  else if (cs_glob_coal_model != NULL) {

    cs_coal_model_t  *cm = cs_glob_coal_model;

    *srrom = &(cm->srrom);

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Enthalpy and Cp based on the JANAF band.
 *
 * \param[in]   ncoel    number of elementary constituents
 * \param[in]   ngazem   number of elementary constituents
 * \param[in]   npo      number of interpolation points
 * \param[in]   nomcoel  names of elementary constituants
 * \param[out]  ehcoel   enthalpy for each elementary species
 *                       (for point i and species j, ehcoel[i*ngazem + j])
 * \param[out]  cpcoel   cp for each elementary species
 *                       (for point i and species j, cpcoel[i*ngazem + j])
 * \param[in]   wmolce   molar mass of each species
 * \param[in]   th       temperature in K
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_enthalpy_and_cp_from_janaf(int           ncoel,
                                         int           ngazem,
                                         int           npo,
                                         const char    nomcoel[][13],
                                         double        ehcoel[],
                                         double        cpcoel[],
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

  for (int ne = 0; ne < ncoel; ne++) {
    for (int inicff = 0; inicff < 2; inicff++) {
      for (int injcff = 0; injcff < 7; injcff++) {
        coeff_therm[injcff][inicff][ne] = 0.0;
      }
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
          for (int injcff = 0; injcff < 7; injcff++) {
            coeff[injcff][inicff][ne] = wcoeff[injcff][inicff];
            coeff_therm[injcff][inicff][ne] = wcoeff[injcff][inicff];
          }
        }
      }
    }

  }  /* end while */

  /* Finish reading if all data has beeen stored */

  impjnf = cs_file_free(impjnf);

  for (int ne = 0; ne < ncoel; ne++) {
    for (int inicff = 0; inicff < 2; inicff++) {
      for (int injcff = 0; injcff < 7; injcff++) {
        coeff_therm[injcff][inicff][ne] = cs_physical_constants_r
                                          * coeff_therm[injcff][inicff][ne]
                                          / wmolce[ne];
      }
    }
  }

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
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
