/*============================================================================
 * Coal combustion model:  setup data from input.
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

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_assert.h"
#include "cs_coal.h"
#include "cs_combustion_model.h"
#include "cs_field.h"
#include "cs_file.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal_read_data.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_read_data.c

  \brief Coal combustion model: setup data from input.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finds the solution of the linear system A*XX=BB , dim(A)=ndim
 *        using the Gauss elimination method with partial pivoting
 *        A,B - modified data ; N - data ; X - result
 *
 * \param[in]   ndim   system dimension
 * \param[in]   aa     system matrix
 * \param[in]   bb     right-hand side
 * \param[in]   xx     solution
 */
/*----------------------------------------------------------------------------*/

static int
_coal_solve_matrix(int     ndim,
                   double  aa[],
                   double  bb[],
                   double  xx[])
{
  int retcode = 0;

  double  pp, ss;

  const double epsil = 1.e-10;

  /* initialization of the solution and test to ensure the elements
     of the matrix diagonal are non zero */

  for (int ii = 0; ii < ndim; ii++) {
    int iw = ii;
    double ww = cs_math_fabs(aa[ii*ndim + ii]);
    for (int jj = ii; jj < ndim; jj++) {
      if (cs_math_fabs(aa[ii*ndim + jj]) > ww) {
        iw = jj;
        ww = cs_math_fabs(aa[ii*ndim + jj]);
      }
    }
    if (ww <= epsil)  {
      retcode = 1;
      break;
    }

    for (int jj = ii; jj < ndim; jj++) {
      pp = aa[jj*ndim + ii];
      aa[jj*ndim + ii] = aa[jj*ndim + iw];
      aa[jj*ndim + iw] = pp;
    }

    pp = bb[ii];
    bb[ii] = bb[iw];
    bb[iw] = pp;

    for (int jj = ii+1; jj < ndim; jj++) {
      pp = aa[ii*ndim + jj] / aa[ii*ndim + ii];
      for (int kk = ii+1; kk < ndim; kk++) {
        aa[kk*ndim + jj] -= pp * aa[kk*ndim + ii];
      }
      bb[jj] -= pp * bb[ii];
    }
  }

  if (retcode != 1) {
    if (cs_math_fabs(aa[ndim*ndim - 1]) < epsil) {
      retcode = 1;
    }
    else {
      xx[ndim - 1] = bb[ndim - 1]/aa[ndim*ndim - 1];
      for (int ii = ndim-2; ii >= 0; ii--) {
        pp = 1. / aa[ii*ndim + ii];
        ss = 0.;
        for (int jj = ii+1; jj < ndim; jj++) {
          ss += aa[jj*ndim + ii] * xx[jj];
        }
        xx[ii] = pp*(bb[ii]-ss);
      }
    }
  }

  return retcode;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup coal model data based on in input.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_read_data(void)
{
  cs_coal_model_t  *cm = cs_glob_coal_model;

  int     atcoel[5][CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

  double  wmolce[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];

  double  ychxle[CS_COMBUSTION_MAX_COALS];
  double  ychxlo[CS_COMBUSTION_MAX_COALS];
  double  ycoch1[CS_COMBUSTION_MAX_COALS];
  double  ycoch2[CS_COMBUSTION_MAX_COALS];
  double  nnch[CS_COMBUSTION_MAX_COALS];
  double  nnckle[CS_COMBUSTION_MAX_COALS];
  double  nhckle[CS_COMBUSTION_MAX_COALS];
  double  ncckle[CS_COMBUSTION_MAX_COALS];
  double  nncklo[CS_COMBUSTION_MAX_COALS];
  double  nhcklo[CS_COMBUSTION_MAX_COALS];
  double  nccklo[CS_COMBUSTION_MAX_COALS];
  double  wchx1c[CS_COMBUSTION_MAX_COALS];
  double  wchx2c[CS_COMBUSTION_MAX_COALS];

  /* NOx model
     --------- */

  /* Terminology:
     MF: mass fraction
     VM: volatile matter */

  // N mol of H, O, N and S in coal
  double  nalpha[CS_COMBUSTION_MAX_COALS];
  double  nbeta[CS_COMBUSTION_MAX_COALS];
  double  ntheta[CS_COMBUSTION_MAX_COALS];
  double  nomega[CS_COMBUSTION_MAX_COALS];

  // N mol of H, O, N and S inchar1 (devolatilization of light VM)
  double  ngama1[CS_COMBUSTION_MAX_COALS];
  double  ndelt1[CS_COMBUSTION_MAX_COALS];
  double  nkapp1[CS_COMBUSTION_MAX_COALS];
  double  nzeta1[CS_COMBUSTION_MAX_COALS];

  // N mol of H, O, N and S in char2 (devolatiliaation of heavy VM)
  double  ngama2[CS_COMBUSTION_MAX_COALS];
  double  ndelt2[CS_COMBUSTION_MAX_COALS];
  double  nkapp2[CS_COMBUSTION_MAX_COALS];
  double  nzeta2[CS_COMBUSTION_MAX_COALS];

  // COMPOSITION of VM
  // H/C ratio of light VM, MF of light VM, H/C ratio of heavy VM, MF of heavy VM
  double  nchx1[CS_COMBUSTION_MAX_COALS];
  // double  ny1ch[CS_COMBUSTION_MAX_COALS];
  double  nchx2[CS_COMBUSTION_MAX_COALS];
  // double  ny2ch[CS_COMBUSTION_MAX_COALS];

  // N mol of CHx1, CO, H2O, H2S, HCN, NH3 of light MV
  double  noxa1[CS_COMBUSTION_MAX_COALS];
  double  noxb1[CS_COMBUSTION_MAX_COALS];
  double  noxc1[CS_COMBUSTION_MAX_COALS];
  double  noxd1[CS_COMBUSTION_MAX_COALS];
  double  noxe1[CS_COMBUSTION_MAX_COALS];
  double  noxf1[CS_COMBUSTION_MAX_COALS];

  // N mol of CHx1, CO, H2O, H2S, HCN, NH3 of heavy VM
  double  noxa2[CS_COMBUSTION_MAX_COALS];
  double  noxb2[CS_COMBUSTION_MAX_COALS];
  double  noxc2[CS_COMBUSTION_MAX_COALS];
  double  noxd2[CS_COMBUSTION_MAX_COALS];
  double  noxe2[CS_COMBUSTION_MAX_COALS];
  double  noxf2[CS_COMBUSTION_MAX_COALS];

  // Composition of products of heterogeneous reaction
  // N mol of O2, CO, HCN and CO (char1 combustion)
  // double  noxh1[CS_COMBUSTION_MAX_COALS];
  double  noxi1[CS_COMBUSTION_MAX_COALS];
  double  noxj1[CS_COMBUSTION_MAX_COALS];
  double  noxk1[CS_COMBUSTION_MAX_COALS];

  // N mol of O2, CO, HCN and CO (char2 combustion)
  // double  noxh2[CS_COMBUSTION_MAX_COALS];
  double  noxi2[CS_COMBUSTION_MAX_COALS];
  double  noxj2[CS_COMBUSTION_MAX_COALS];
  double  noxk2[CS_COMBUSTION_MAX_COALS];

  /* composition of reactive coal: alpha(c) = hch(c)/cch(c)
   *                               beta (c) = och(c)/cch(c)
   *                               theta(c) = sch(c)/cch(c)
   *                               omega (c) = nch(c)/cch(c) */
  double alpha[CS_COMBUSTION_MAX_COALS];
  double beta[CS_COMBUSTION_MAX_COALS];
  double theta[CS_COMBUSTION_MAX_COALS];
  double omega[CS_COMBUSTION_MAX_COALS];

  /* composition of coke: gamma(c) = hck(c)/cck(c)
   *                      delta(c) = ock(c)/cck(c)
   *                      kappa(c) = sck(c)/cck(c)
   *                      zeta(c) = nck(c)/cck(c) */
  double gamma[CS_COMBUSTION_MAX_COALS];
  double delta[CS_COMBUSTION_MAX_COALS];
  double kappa[CS_COMBUSTION_MAX_COALS];
  double zeta[CS_COMBUSTION_MAX_COALS];

  const int  ngazem = CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS;
  const int  npot = CS_COMBUSTION_COAL_MAX_TABULATION_POINTS;
  double *ehcoel, *cpcoel;
  BFT_MALLOC(ehcoel, ngazem * npot, double);
  BFT_MALLOC(cpcoel, ngazem * npot, double);

  const int n_coals = cm->n_coals;

  /* "Read" specific data
   * ==================== */

  const int ncoel = 13;

  if (ncoel > CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS) {
    bft_error
      (__FILE__, __LINE__, 0,
       _(" %s: pulverized coal model:\n\n"
         " The current number of species is %d\n"
         " but must be less than or equal to %d\n"
         "\n"
         " Check the setup parameters."),
       __func__, CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS, ncoel);
  }

  /* Number of enthalpy-temperature tabulation points */

  cm->n_tab_points = 8;
  cs_assert(cm->n_tab_points <= CS_COMBUSTION_COAL_MAX_TABULATION_POINTS);

  /* Build names of elementary constituents */

  const char *el_comp_names[] = {"CH4", "C2H4", "CO", "H2S", "H2", "HCN", "NH3",
                                 "O2", "CO2", "H2O", "SO2", "N2", "C(S)"};

  char  nomcoel[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS][13];

  for (int ice = 0; ice < ncoel; ice++)
    strncpy(nomcoel[ice], el_comp_names[ice], 13);
  for (int ice = ncoel; ice < ngazem; ice++)
    nomcoel[ice][0] = '\0';

  /* Min and max temperatures */
  double  tmin = 300., tmax = 2400.;

  /* Number of atomic species  (C, H, O, S et N) */

  const int nato = 5;
  cs_assert(nato <= CS_COMBUSTION_COAL_MAX_ATOMIC_SPECIES);

  /* Molar masse of atomic species */

  cm->wmolat[0] = 0.012;
  cm->wmolat[1] = 0.001;
  cm->wmolat[2] = 0.016;
  cm->wmolat[3] = 0.014;
  cm->wmolat[4] = 0.032;

  /* Composition of elemantary constituants based on elementary species */

  // CH4
  atcoel[0][0] = 1;
  atcoel[1][0] = 4;
  atcoel[2][0] = 0;
  atcoel[3][0] = 0;
  atcoel[4][0] = 0;

  // C2H4
  atcoel[0][1] = 2;
  atcoel[1][1] = 4;
  atcoel[2][1] = 0;
  atcoel[3][1] = 0;
  atcoel[4][1] = 0;

  // CO
  atcoel[0][2] = 1;
  atcoel[1][2] = 0;
  atcoel[2][2] = 1;
  atcoel[3][2] = 0;
  atcoel[4][2] = 0;

  // H2S
  atcoel[0][3] = 0;
  atcoel[1][3] = 2;
  atcoel[2][3] = 0;
  atcoel[3][3] = 0;
  atcoel[4][3] = 1;

  // H2
  atcoel[0][4] = 0;
  atcoel[1][4] = 2;
  atcoel[2][4] = 0;
  atcoel[3][4] = 0;
  atcoel[4][4] = 0;

  // HCN
  atcoel[0][5] = 1;
  atcoel[1][5] = 1;
  atcoel[2][5] = 0;
  atcoel[3][5] = 1;
  atcoel[4][5] = 0;

  // NH3
  atcoel[0][6] = 0;
  atcoel[1][6] = 3;
  atcoel[2][6] = 0;
  atcoel[3][6] = 1;
  atcoel[4][6] = 0;

  // O2
  atcoel[0][7] = 0;
  atcoel[1][7] = 0;
  atcoel[2][7] = 2;
  atcoel[3][7] = 0;
  atcoel[4][7] = 0;

  // CO2
  atcoel[0][8] = 1;
  atcoel[1][8] = 0;
  atcoel[2][8] = 2;
  atcoel[3][8] = 0;
  atcoel[4][8] = 0;

  // H2O
  atcoel[0][9] = 0;
  atcoel[1][9] = 2;
  atcoel[2][9] = 1;
  atcoel[3][9] = 0;
  atcoel[4][9] = 0;

  // SO2
  atcoel[0][10] = 0;
  atcoel[1][10] = 0;
  atcoel[2][10] = 2;
  atcoel[3][10] = 0;
  atcoel[4][10] = 1;

  // N2
  atcoel[0][11] = 0;
  atcoel[1][11] = 0;
  atcoel[2][11] = 0;
  atcoel[3][11] = 2;
  atcoel[4][11] = 0;

  // C(S)
  atcoel[0][12] = 1;
  atcoel[1][12] = 0;
  atcoel[2][12] = 0;
  atcoel[3][12] = 0;
  atcoel[4][12] = 0;

  /* Compute molar masses of elementary constituants */

  for (int ice = 0; ice < ncoel; ice++) {
    wmolce[ice] = 0;
    for (int iat = 0; iat < nato; iat++)
      wmolce[ice] += atcoel[iat][ice] * cm->wmolat[iat];
  }

  /* Compute additional data
     ======================= */

  /* Temperature discretization */

  for (int it = 0; it < cm->n_tab_points; it++) {
    cm->th[it] = it * (tmax-tmin) / (cm->n_tab_points - 1) + tmin;
  }

  /* Compute enthalpies for the different usual species */

  cs_combustion_enthalpy_and_cp_from_janaf(ncoel,
                                           ncoel,
                                           cm->n_tab_points,
                                           nomcoel,
                                           ehcoel,
                                           cpcoel,
                                           wmolce,
                                           cm->th);

  BFT_FREE(cpcoel);

  /* Compute enthalpy - temperature tabulation for gas mix */

  /* Number of gaseous constituants
     Attention: we also count CH4 and the CH2 monomere
     we add: H2S, H2 et SO2 (we go from  5 to 8). */

  cm->n_gas_el_comp = 2 + 10 + 2*n_coals;
  cm->n_gas_species = 2 + 10;

  /* Define indexes in wmole and ehgaze.
     Remark: these indexes will also be used for the IYMI pointers
     array (property indices).

     For now, these indexes are 1-based, but can be converted to
     0-based when the Fortran-C conversion is finished. */

  int s_id  = 0;
  int ichx1 = s_id++;
  int ichx2 = s_id++;
  int ico   = s_id++;
  int ih2s  = s_id++;
  int ihy   = s_id++;
  int ihcn  = s_id++;
  int inh3  = s_id++;
  int io2   = s_id++;
  int ico2  = s_id++;
  int ih2o  = s_id++;
  int iso2  = s_id++;
  int in2   = s_id++;

  cm->ichx1 = ichx1 + 1;
  cm->ichx2 = ichx2 + 1;
  cm->ico   = ico + 1;
  cm->ih2s  = ih2s + 1;
  cm->ihy   = ihy + 1;
  cm->ihcn  = ihcn + 1;
  cm->inh3  = inh3 + 1;
  cm->io2   = io2 + 1;
  cm->ico2  = ico2 + 1;
  cm->ih2o  = ih2o + 1;
  cm->iso2  = iso2 + 1;
  cm->in2   = in2 + 1;

  for (int icha = 0; icha < n_coals; icha++) {
    cm->ichx1c[icha] = s_id + icha + 1;
    cm->ichx2c[icha] = s_id + n_coals + icha + 1;
  }

  /* Fill ehgaze and wmole from ehcoel and wmolce */

  /* Attention:
     We take by default CH4 for CHX1m
     and the monomere CH2 (plays the role of C2H4) for CHX2m */

  const int ntp = cm->n_tab_points;

  for (int it = 0; it < ntp; it++) {
    cm->ehgaze[it][ichx1] = ehcoel[ncoel*it];
    cm->ehgaze[it][ichx2] = ehcoel[ncoel*it + 1];
    cm->ehgaze[it][ico]   = ehcoel[ncoel*it + 2];
    cm->ehgaze[it][ih2s]  = ehcoel[ncoel*it + 3];
    cm->ehgaze[it][ihy]   = ehcoel[ncoel*it + 4];
    cm->ehgaze[it][ihcn]  = ehcoel[ncoel*it + 5];
    cm->ehgaze[it][inh3]  = ehcoel[ncoel*it + 6];
    cm->ehgaze[it][io2]   = ehcoel[ncoel*it + 7];
    cm->ehgaze[it][ico2]  = ehcoel[ncoel*it + 8];
    cm->ehgaze[it][ih2o]  = ehcoel[ncoel*it + 9];
    cm->ehgaze[it][iso2]  = ehcoel[ncoel*it + 10];
    cm->ehgaze[it][in2]   = ehcoel[ncoel*it + 11];
  }
  cm->wmole[ichx1] = wmolce[0];

  int iatc = cs_coal_atom_id_c;
  int iath = cs_coal_atom_id_h;
  int iato = cs_coal_atom_id_o;
  int iatn = cs_coal_atom_id_n;
  int iats = cs_coal_atom_id_s;

  /* Attention: we must use the molar mass of the CH2 monomere
     and not that of C2H4 */
  cm->wmole[ichx2] = cm->wmolat[iatc] + cm->wmolat[iath]*2.;
  cm->wmole[ico]   = wmolce[2];
  cm->wmole[ih2s]  = wmolce[3];
  cm->wmole[ihy]   = wmolce[4];
  cm->wmole[ihcn]  = wmolce[5];
  cm->wmole[inh3]  = wmolce[6];
  cm->wmole[io2]   = wmolce[7];
  cm->wmole[ico2]  = wmolce[8];
  cm->wmole[ih2o]  = wmolce[9];
  cm->wmole[iso2]  = wmolce[10];
  cm->wmole[in2]   = wmolce[11];

  /* Build enthalpy - temperature tabulation for the solid
     reactive coal, coke, and ashes */

  /* Numbur of solid components */

  cm->nsolid = 4*n_coals;

  /* Definition of ich, ick et iash pointers (1-based for Fortran) */

  s_id = 1;
  for (int icha = 0; icha < n_coals; icha++) {
    cm->ich[icha] = s_id++;
    cm->ick[icha] = s_id++;
    cm->iash[icha] = s_id++;
    cm->iwat[icha] = s_id++;
  }

  /* Composition of reactive coal
     under the forme CH(alpha)O(beta)N(theta)S(omega)
     and compute molar mass */

  for (int icha = 0; icha < n_coals; icha++) {

    alpha[icha] =    (cm->hch[icha] / cm->wmolat[iath])
                   / (cm->cch[icha] / cm->wmolat[iatc]);
    beta[icha]  =    (cm->och[icha] / cm->wmolat[iato])
                   / (cm->cch[icha] / cm->wmolat[iatc]);

    /* Currently, nitrous species are not considered in the combustion
       of the gas phase. So we must ensure that no nitrous specie is
       released in the gas phase. */

    theta[icha]  =   ((cm->nch[icha]-cm->nch[icha]) / cm->wmolat[iatn])
                   / (cm->cch[icha] / cm-> wmolat[iatc]);

    omega[icha] =   (cm->sch[icha] / cm->wmolat[iats])
                  / (cm->cch[icha] / cm->wmolat[iatc]);

    /* Correction of correlations theta - sulfur and omega - nitrogen */

    int ich_icha = cm->ich[icha]-1;
    cm->wmols[ich_icha] =   cm->wmolat[iatc]
                          + alpha[icha] * cm->wmolat[iath]
                          + beta[icha]  * cm->wmolat[iato]
                          + theta[icha]  * cm->wmolat[iatn]
                          + omega[icha] * cm->wmolat[iats];
  }

  /* Transformation distribution coefs of nitrogen in HCN etNH3 */

  if (cm->ieqnox) {
    for (int icha = 0; icha < n_coals; icha++) {
      if (cm->nch[icha] > 0) {
        double som1 = cm->crepn1[icha][0]+cm->crepn1[icha][1];
        double som2 = cm->crepn2[icha][0]+cm->crepn2[icha][1];
        if (som1 < 0. || som2 < 0.) {
          bft_error
            (__FILE__, __LINE__, 0,
             _(" %s: pulverized coal model:\n\n"
               "Nitrogen was chosen in the coal composition but the\n"
               "distribution between HCN and NH3 is zero for coal %d.\n"
               "\n"
               " Check the setup parameters."), __func__, icha+1);
        }
        cm->crepn1[icha][0] /= som1;
        cm->crepn1[icha][1] /= som1;
        cm->crepn2[icha][0] /= som2;
        cm->crepn2[icha][1] /= som2;
      }
    }
  }

  /* Transformation by Schaff's formula of the PCI over dry on PCI over pure
     when ipci = 1 */

  /* Compute PCI (pure) */

  const double fcor = 1.2;
  const cs_real_t xcal2j = 4.1855e0;

  for (int icha = 0; icha < n_coals; icha++) {

    double pcisec = 0., pcipur, pcibrut, xwatpc, pcspur;

    xwatpc = cm->xwatch[icha] * 100.;

    /* PCI(dry) ---> PCI(pure) */
    if (cm->ipci[icha] == 1) {
      pcisec = cm->pcich[icha] / 1000. / xcal2j;
      // Compute PCI over pure
      // DS remark: no usage of xashsec instead of xashpc
      // double xashpc = cm->xashch[icha] * 100.;
      pcipur =  ( pcisec+5.95*((1.-fcor)*cm->xashsec[icha])) *100.
               / (100. -fcor*cm->xashsec[icha]);
      pcipur = pcipur * 1000. * xcal2j; // kcal/kg to J/kg
      cm->pcich[icha] = pcipur;
    }

    /* PCI(raw) ---> PCI(pure) */
    else if (cm->ipci[icha] == 2) {
      pcibrut = cm->pcich[icha] / 1000. /xcal2j;
      pcipur =   ( pcibrut + 5.95*( xwatpc+(1.-fcor)
                                   *cm->xashsec[icha])) *100.
               /  (100. - xwatpc -fcor*cm->xashsec[icha]);
      pcipur = pcipur * 1000. * xcal2j;  // kcal/kg to J/kg
      cm->pcich[icha] = pcipur;
    }

    else if (cm->ipci[icha] >= 3) {

      double pcsbrut, pcssec;

      if (cm->ipci[icha] == 3) {   /* PCS(pure) ---> PCI(pure) */
        pcspur = cm->pcich[icha] / 1000.;  // J/kg to KJ/kg
        pcssec = pcspur * (100.-cm->xashsec[icha]) / 100.;
        pcisec = pcssec -226. * cm->hch[icha];
      }

      else if (cm->ipci[icha] == 4) {   /* PCS(dry) ---> PCI(pure) */
        pcisec = cm->pcich[icha] / 1000.;  // J/kg to KJ/kg
      }
      else if (cm->ipci[icha] == 5) {   /* PCS(raw)---> PCI(pure) */
        pcsbrut = cm->pcich[icha] / 1000.; // J/kg to KJ/kg
        pcssec = pcsbrut*100. / (100.-xwatpc);
        pcisec = pcssec -226. * cm->hch[icha];
      }
      else if (cm->ipci[icha] == 6) {   /* IGT correlation */
        // Compute PCS(dry) in KJ/kg
        pcssec =    340.94 * cm->cch[icha]
                 + 1322.98 * cm->hch[icha]
                 + 68.3844 * cm->sch[icha]
                 - 119.86  * (cm->och[icha] + cm->nch[icha])
                 - 15.305  * cm->xashsec[icha];
        pcisec = pcssec -226. * cm->hch[icha];
      }

      // Compute PCI(dry) from PCS(dry)

      pcisec = pcisec / xcal2j;  // KJ/kg to kcal/kg

      // Compute PCI over pure
      pcipur =   ( pcisec+5.95*((1.-fcor)*cm->xashsec[icha])) *100.
               / (100. - fcor*cm->xashsec[icha]);
      pcipur = pcipur*1000.*xcal2j;  // kcal/kg to J/kg

      cm->pcich[icha] = pcipur;

    }
  }

  // Compute: H02 for CH, CK, ASH and WT

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;
    double hco20 = cm->ehgaze[0][ico2] * cm->wmole[ico2];
    double hh2o0 = cm->ehgaze[0][ih2o] * cm->wmole[ih2o];
    double ho20  = cm->ehgaze[0][io2] * cm->wmole[io2];
    double hso20 = cm->ehgaze[0][iso2] * cm->wmole[iso2];

    // Correction of correlation theta - sulfur

    cm->h02ch[icha] = cm->pcich[icha]
                      + (  hco20 + alpha[icha]/2.*hh2o0
                                 + omega[icha]   *hso20
                         - ( 1. + alpha[icha]/4.
                                - beta[icha]/2.
                               + omega[icha]) *ho20) / cm->wmols[ich_icha];

    cm->eh0sol[ich_icha] = cm->h02ch[icha];
  }

  // For now H02CK an H02ASH are zero

  for (int icha = 0; icha < n_coals; icha++) {
    int ick_icha = cm->ick[icha]-1;
    int iash_icha = cm->iash[icha]-1;
    cm->eh0sol[ick_icha] = 0.;
    cm->eh0sol[iash_icha] = 0.;
  }

  /* Compute H02 for liquid water: based on JANAF
     H0 = -15871802.2
     CP = 4181.35 */

  for (int icha = 0; icha < n_coals; icha++) {
    int iwat_icha = cm->iwat[icha]-1;
    cm->eh0sol[iwat_icha] = -15871802.2;
  }

  /* Build enthalpie/temperature tabulation for coal
     For now, same as that for gas */

  // Number of tabulations: npoc

  cm->npoc = cm->n_tab_points;

  for (int i = 0; i < cm->npoc; i++)
    cm->thc[i] = cm->th[i];

  /* Compute ehsoli for reactive coal
     If        cp2ch > 0 : HCH = H02CH + CP2CH(T2-TREFTH)
     Otherwise           : we consider PCI constant whatever T2 */

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;

    if (cm->cp2ch[icha] > cs_coal_epsilon) {
      for (int it = 0; it < cm->npoc; it++) {
        cm->ehsoli[it][ich_icha] =   cm->eh0sol[ich_icha]
                                   + cm->cp2ch[icha]*(cm->thc[it]-cs_coal_trefth);
      }
    }
    else {
      for (int it = 0; it < cm->npoc; it++) {
        double hco20 = cm->ehgaze[it][ico2] * cm->wmole[ico2];
        double hh2o0 = cm->ehgaze[it][ih2o] * cm->wmole[ih2o];
        double ho20  = cm->ehgaze[it][io2]  * cm->wmole[io2];
        double hso20 = cm->ehgaze[it][iso2] * cm->wmole[iso2];

        // Correction of correlation theta - sulfur
        cm->ehsoli[it][ich_icha]
          =   cm->pcich[icha]
            + (  hco20 + alpha[icha]/2.*hh2o0
                       + omega[icha] *hso20
               - (1.   + alpha[icha]/4.
                       - beta[icha]/2.
                       + omega[icha]) *ho20) / cm->wmols[ich_icha];
      }
    }

    gamma[icha] = 0;
    delta[icha] = 0;
    kappa[icha] = 0;
    zeta[icha] = 0;
  }

  /* Computation relative to coke

     By default
     Coke = solid carbon
     gamma = 0, delta = 0, kappa = 0, zeta = 0
     If    cp2ch > 0 : HCK = H02CH + CP2CH(T2-TREFTH)
     Else:           : HCK = Enthalpu of pure carbon */

  for (int icha = 0; icha < n_coals; icha++) {
    int ick_icha = cm->ick[icha]-1;
    cm->wmols[ick_icha] = wmolce[ncoel - 1];
    if (cm->cp2ch[icha] > cs_coal_epsilon) {
      for (int it = 0; it < cm->npoc; it++) {
        cm->ehsoli[it][ick_icha]
          =   cm->eh0sol[ick_icha]
            + cm->cp2ch[icha] * (cm->thc[it] - cs_coal_trefth);
      }
    }
    else {
      for (int it = 0; it < cm->npoc; it++) {
        cm->ehsoli[it][ick_icha] = ehcoel[ncoel*it + ncoel-1];
      }
    }
  }

  BFT_FREE(ehcoel);

  /* Coke = CH(gamma)O(delta)N(kappa)S(zeta) */

  for (int icha = 0; icha < n_coals; icha++) {
    int ick_icha = cm->ick[icha]-1;
    if (cm->pcick[icha] > cs_coal_epsilon) {
      gamma[icha] = cm->hck[icha] / cm->cck[icha];
      delta[icha] = cm->ock[icha] / cm->cck[icha];

      // Correction of correlation kappa - sulfur.
      kappa[icha] = cm->nck[icha] / cm->cck[icha];

     // Correction of correlation zeta - nitrogen.
      zeta [icha] = cm->sck[icha] / cm->cck[icha];

     // Correction of correlations kappa - sulfur and zeta - nitrogen.
     cm->wmols[ick_icha] =   cm->wmolat[iatc]
                           + gamma[icha] * cm->wmolat[iath]
                           + delta[icha] * cm->wmolat[iato]
                           + kappa[icha] * cm->wmolat[iatn]
                           + zeta[icha]  * cm->wmolat[iats];

     // Consider PCI constant whatever T.
     for (int it = 0; it < cm->npoc; it++) {
       double hco20 = cm->ehgaze[it][ico2] * cm->wmole[ico2];
       double hh2o0 = cm->ehgaze[it][ih2o] * cm->wmole[ih2o];
       double hso20 = cm->ehgaze[it][iso2] * cm->wmole[iso2];
       double ho20  = cm->ehgaze[it][io2]  * cm->wmole[io2];

       // Correction of correlation kappa - sulfur.
       cm->ehsoli[it][ick_icha]
         =   cm->pcick[icha]
           + (hco20 + gamma[icha]/2.*hh2o0
                    + zeta[icha]    *hso20
           - (1.    + gamma[icha]/4.
                    - delta[icha]/2.
                    + zeta[icha]) *ho20) / cm->wmols[ick_icha];
     }
    }
  }

  // Computation relative to ashes

  for (int icha = 0; icha < n_coals; icha++) {
    int iash_icha = cm->iash[icha]-1;
    if (cm->cp2ch[icha] > cs_coal_epsilon) {
      for (int it = 0; it < cm->npoc; it++) {
        cm->ehsoli[it][iash_icha]
          =   cm->eh0sol[iash_icha]
            + cm->cp2ch[icha] * (cm->thc[it] - cs_coal_trefth);
      }
    }
    else {
      for (int it = 0; it < cm->npoc; it++) {
        cm->ehsoli[it][iash_icha]
          =   cm->h0ashc[icha]
            + cm->cpashc[icha] * (cm->thc[it] - cs_coal_trefth);
      }
    }
    cm->wmols[iash_icha] = 0;
  }

  // Computation relative to water

  for (int icha = 0; icha < n_coals; icha++) {
    int iwat_icha = cm->iwat[icha]-1;
    const double cp2wat = 4181.35;

    for (int it = 0; it < cm->npoc; it++) {
      cm->ehsoli[it][iwat_icha]
        =   cm->eh0sol[iwat_icha]
          + cp2wat * (cm->thc[it] - cs_coal_trefth);
    }
  }

  // Computations relative to volatile matter

  /* We check that the same option is used for Y1 and Y2 for each coal */

  for (int icha = 0; icha < n_coals; icha++) {
    if (cm->iy1ch[icha] != cm->iy2ch[icha]) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "A different option is used for Y1 and Y2 for coal %d\n"
           "  iy1ch = %d, iy2ch = %d\n"
           "\n"
           " Check the setup parameters."), __func__,
         icha+1, cm->iy1ch[icha], cm->iy2ch[icha]);
    }
    if (   cm->iy1ch[icha] < 0 || cm->iy1ch[icha] > 2
        || cm->iy2ch[icha] < 0 || cm->iy2ch[icha] > 2) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "An unavailable option is used for Y1 or Y2 for coal %d\n"
           "  iy1ch = %d, iy2ch = %d\n"
           "\n"
           " Check the setup parameters."), __func__,
         icha+1, cm->iy1ch[icha], cm->iy2ch[icha]);
    }
  }

  // Light volatile matter : [CH(CHX1); CO]
  //   CH(alpha)O(beta)N(omega)S(theta) --> A1 CH(CHX1) + B1 CO + D1 H2S
  //                            + E1 HCN + F1 NH3
  //                            + (1-A1-B1-E1) CH(gamma)O(delta)N(zeta)S(kappa)
  // If IY1CH = 0 CHX1 fixed , Y1CH computed
  // If IY1CH = 1 Y1CH fixed , CHX1 computed
  // If IY1CH = 2 Y1CH fixed , CHX1 fixed, we add water
  //   CH(alpha)O(beta)N(omega)S(theta) --> A1 CH(CHX1) + B1 CO + C1 H2O
  //                            + D1 H2S + E1 HCN + F1 NH3
  //                            + (1-A1-B1-E1) CH(gamma)O(delta)N(zeta)S(kappa)

  for (int icha = 0; icha < n_coals; icha++) {

    double  sm[8], solu[8], mat[8][8];

    mat[0][0] = -gamma[icha];
    mat[1][0] = -gamma[icha];
    mat[2][0] = 2.;
    mat[3][0] = 2.;
    mat[4][0] = 1.-gamma[icha];
    mat[5][0] = 3.;
    mat[6][0] = 1.;
    mat[7][0] = 0.;

    mat[0][1] = -delta[icha];
    mat[1][1] = 1.-delta[icha];
    mat[2][1] = 1.;
    mat[3][1] = 0.;
    mat[4][1] = -delta[icha];
    mat[5][1] = 0.;
    mat[6][1] = 0.;
    mat[7][1] = 0.;

    mat[0][2] = -kappa[icha];
    mat[1][2] = -kappa[icha];
    mat[2][2] = 0.;
    mat[3][2] = 0.;
    mat[4][2] = 1.-kappa[icha];
    mat[5][2] = 1.;
    mat[6][2] = 0.;
    mat[7][2] = 0.;

    mat[0][3] = 0.;
    mat[1][3] = 0.;
    mat[2][3] = 0.;
    mat[3][3] = 0.;
    mat[4][3] = cm->crepn1[icha][1];
    mat[5][3] =-cm->crepn1[icha][0];
    mat[6][3] = 0.;
    mat[7][3] = 0.;

    mat[0][4] = -zeta[icha];
    mat[1][4] = -zeta[icha];
    mat[2][4] = 0.;
    mat[3][4] = 1.;
    mat[4][4] = -zeta[icha];
    mat[5][4] = 0.;
    mat[6][4] = 0.;
    mat[7][4] = 0.;

    mat[0][5] = cm->wmolat[iatc];
    mat[1][5] = cm->wmole[ico];
    mat[2][5] = cm->wmole[ih2o];
    mat[3][5] = cm->wmole[ih2s];
    mat[4][5] = cm->wmole[ihcn];
    mat[5][5] = cm->wmole[inh3];
    mat[6][5] = cm->wmolat[iath];

    // Correction of correlations theta - sulfur and omega - nitrogen.
    mat[7][5] = -(cm->wmolat[iatc] + alpha[icha]*cm->wmolat[iath]
                                   + beta[icha] *cm->wmolat[iato]
                                   + theta[icha] *cm->wmolat[iatn]
                                   + omega[icha]*cm->wmolat[iats]);

    // Right-hand side
    sm[0] = alpha[icha] - gamma[icha];
    sm[1] = beta[icha]  - delta[icha];
    sm[2] = theta[icha]  - kappa[icha];
    sm[3] = 0.;
    sm[4] = omega[icha] - zeta[icha];
    sm[5] = 0.;

    // Complete the  matrix and RHS base on value of IY1CH

    if (cm->iy1ch[icha] == 0) {
      cm->chx1[icha] = 4.;

      mat[0][6] = -cm->chx1[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = 0.;
      sm[7] = 0.;
    }
    else if (cm->iy1ch[icha] == 1)  {
      mat[0][6] = 0.;
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 0.;
      mat[7][6] = 1.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = cm->y1ch[icha];
      sm[7] = 0.;
    }
    else if (cm->iy1ch[icha] == 2) {
      cm->chx1[icha] = 4.;

      mat[0][6] = -cm->chx1[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 0.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 1.;

      sm[6] = 0.;
      sm[7] = cm->y1ch[icha];
    }

    int ierr = _coal_solve_matrix(8, (double *)mat, sm, solu);

    if (ierr == 0) {
      cm->a1[icha] = solu[0];
      cm->b1[icha] = solu[1];
      cm->c1[icha] = solu[2];
      cm->d1[icha] = solu[3];
      cm->e1[icha] = solu[4];
      cm->f1[icha] = solu[5];
      if (cm->iy1ch[icha] == 0)
        cm->y1ch[icha] = solu[7];
      else if (cm->iy1ch[icha] == 1)
        cm->chx1[icha] = solu[6]/solu[0];
    }
    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "Non invertible matrix for volatile matter computation for coal %d.\n"
           "\n"
           " Check the setup parameters."), __func__, icha+1);
    }

  }

  // Heavy volatile matter: [CH(CHX2); CO]
  //   CH(alpha)O(beta)N(omega)S(theta) --> A2 CH(CHX1) + B2 CO + D2 H2S
  //                            + E2 HCN + F2 NH3
  //                            + (1-A2-B2-E2) CH(gamma)O(delta)N(zeta)S(kappa)
  //     Si IY2CH = 0 CHX2 fixe, Y2CH calcule
  //     Si IY2CH = 1 Y2CH fixe, CHX2 calcule
  //     Si IY2CH = 2 Y2CH fixe, CHX2 fixe, on ajoute de l'eau
  //   CH(alpha)O(beta)N(omega)S(theta) --> A2 CH(CHX1) + B2 CO + C2 H2O
  //                            + D2 H2S + E2 HCN + F2 NH3
  //                            + (1-A2-B2-E2) CH(gamma)O(delta)N(zeta)S(kappa)

  for (int icha = 0; icha < n_coals; icha++) {

    double  sm[8], solu[8], mat[8][8];

    mat[0][0] = -gamma[icha];
    mat[1][0] = -gamma[icha];
    mat[2][0] = 2.;
    mat[3][0] = 2.;
    mat[4][0] = 1.-gamma[icha];
    mat[5][0] = 3.;
    mat[6][0] = 1.;
    mat[7][0] = 0.;

    mat[0][1] = -delta[icha];
    mat[1][1] = 1.-delta[icha];
    mat[2][1] = 1.;
    mat[3][1] = 0.;
    mat[4][1] = -delta[icha];
    mat[5][1] = 0.;
    mat[6][1] = 0.;
    mat[7][1] = 0.;

    mat[0][2] = -kappa[icha];
    mat[1][2] = -kappa[icha];
    mat[2][2] = 0.;
    mat[3][2] = 0.;
    mat[4][2] = 1.-kappa[icha];
    mat[5][2] = 1.;
    mat[6][2] = 0.;
    mat[7][2] = 0.;

    mat[0][3] = 0.;
    mat[1][3] = 0.;
    mat[2][3] = 0.;
    mat[3][3] = 0.;
    mat[4][3] =  cm->crepn2[icha][1];
    mat[5][3] = -cm->crepn2[icha][0];
    mat[6][3] = 0.;
    mat[7][3] = 0.;

    mat[0][4] = -zeta[icha];
    mat[1][4] = -zeta[icha];
    mat[2][4] = 0.;
    mat[3][4] = 1.;
    mat[4][4] = -zeta[icha];
    mat[5][4] = 0.;
    mat[6][4] = 0.;
    mat[7][4] = 0.;

    mat[0][5] = cm->wmolat[iatc];
    mat[1][5] = cm->wmole[ico];
    mat[2][5] = cm->wmole[ih2o];
    mat[3][5] = cm->wmole[ih2s];
    mat[4][5] = cm->wmole[ihcn];
    mat[5][5] = cm->wmole[inh3];
    mat[6][5] = cm->wmolat[iath];

    // Correction of correlations theta - sulfur and omega - nitrogen.
    mat[7][5] = -(cm->wmolat[iatc] + alpha[icha]*cm->wmolat[iath]
                                   + beta[icha] *cm->wmolat[iato]
                                   + theta[icha]*cm->wmolat[iatn]
                                   + omega[icha]*cm->wmolat[iats]);

    // RHS

    sm[0] = alpha[icha] - gamma[icha];
    sm[1] = beta [icha] - delta[icha];
    sm[2] = theta[icha] - kappa[icha];
    sm[3] = 0.;
    sm[4] = omega[icha] - zeta[icha];
    sm[5] = 0.;

    // Complete matrix and RHS based on iy2ch

    if (cm->iy2ch[icha] == 0) {
      cm->chx2[icha] = 2.;

      mat[0][6] = -cm->chx2[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = 0.;
      sm[7] = 0.;
    }
    else if (cm->iy2ch[icha] == 1) {
      mat[0][6] = 0.;
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 0.;
      mat[7][6] = 1.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = cm->y2ch[icha];
      sm[7] = 0.;
    }
    else if (cm->iy2ch[icha] == 2) {
      cm->chx2[icha] = 2.;

      mat[0][6] = -cm->chx2[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 0.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 1.;

      sm[6] = 0.;
      sm[7] = cm->y2ch[icha];
    }

    int ierr = _coal_solve_matrix(8, (double *)mat, sm, solu);

    if (ierr == 0) {
      cm->a2[icha] = solu[0];
      cm->b2[icha] = solu[1];
      cm->c2[icha] = solu[2];
      cm->d2[icha] = solu[3];
      cm->e2[icha] = solu[4];
      cm->f2[icha] = solu[5];
      if (cm->iy2ch[icha] == 0)
        cm->y2ch[icha] = solu[7];
      else if (cm->iy2ch[icha] == 1)
        cm->chx2[icha] = solu[6]/solu[0];
    }
    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "Non invertible matrix for volatile matter computation for coal %d.\n"
           "\n"
           " Check the setup parameters."), __func__, icha+1);
    }

  }

  // Compute ehgaze and wmole for species CH(CHX1) and CH(CHX2)

  // CH(CHX1) species

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;
    int ick_icha = cm->ick[icha]-1;
    int ichx1c_icha = cm->ichx1c[icha]-1;
    cm->wmole[ichx1c_icha]
      = cm->wmolat[iatc] + cm->chx1[icha]*cm->wmolat[iath];
    if (cm->iy1ch[icha] == 0 || cm->iy1ch[icha] == 2) {
      for (int it = 0; it < cm->n_tab_points; it++) {
        cm->ehgaze[it][ichx1c_icha] = cm->ehgaze[it][ichx1];
      }
    }
    else {
      // We assume D(HDEV,1,ICHA) = 0
      for (int it = 0; it < cm->n_tab_points; it++) {
        double den1 = cm->a1[icha]*cm->wmole[ichx1c_icha];
        cm->ehgaze[it][ichx1c_icha]
          =  (     (                         cm->ehsoli[it][ich_icha]
                     - (1.-cm->y1ch[icha]) * cm->ehsoli[it][ick_icha])
                 * cm->wmols[ich_icha]
              - cm->b1[icha] * cm->wmole[ico] * cm->ehgaze[it][ico])
           / den1;
      }
    }
  }

  // CH(CHX2) species

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;
    int ick_icha = cm->ick[icha]-1;
    int ichx2c_icha = cm->ichx2c[icha]-1;
    cm->wmole[ichx2c_icha]
      = cm->wmolat[iatc] + cm->chx2[icha]*cm->wmolat[iath];
    if (cm->iy2ch[icha] == 0 || cm->iy2ch[icha] == 2) {
      for (int it = 0; it < cm->n_tab_points; it++) {
        cm->ehgaze[it][ichx2c_icha] = cm->ehgaze[it][ichx2];
      }
    }
    else {
      // We assume D(HDEV,2,ICHA) = 0
      for (int it = 0; it < cm->n_tab_points; it++) {
        double den2 = cm->a2[icha]*cm->wmole[ichx2c_icha];
        cm->ehgaze[it][ichx2c_icha]
          =  (     (                         cm->ehsoli[it][ich_icha]
                     - (1.-cm->y2ch[icha]) * cm->ehsoli[it][ick_icha])
                 * cm->wmols[ich_icha]
              - cm->b2[icha] * cm->wmole[ico] * cm->ehgaze[it][ico])
           / den2;
      }
    }
  }

  // Computation fo the different classes

  for (int icla = 0; icla < cm->nclacp; icla++) {
    int ichcor_icla = cm->ichcor[icla] - 1;

    cm->dia2mn[icla] = 0;                        // min diameter (m)
    cm->rho20[icla]  = cm->rho0ch[ichcor_icla];  // density (kg/m3)
    cm->rho2mn[icla] = cm->rho20[icla] * cm->xashch[ichcor_icla];

    // initial particle mass
    cm->xmp0[icla]   =    cm->rho20[icla]
                       * cs_math_pi * cs_math_pow3(cm->diam20[icla])/6.;

    // mass of ashes of particule
    cm->xmash[icla] = cm->xmp0[icla] * cm->xashch[ichcor_icla];
  }

  /* New NOx model (Marcus Charwath and Sandro Dal Secco)
     ---------------------------------------------------- */

  // Composition of reactive coal under the form
  // CH(nalpha)O(nbeta)N(ntheta)S(nomega)

  for (int icha = 0; icha < n_coals; icha++) {
    nalpha[icha] =   (cm->hch[icha] / cm->wmolat[iath])
                   / (cm->cch[icha] / cm->wmolat[iatc]);
    nbeta[icha]  =   (cm->och[icha] / cm->wmolat[iato])
                   / (cm->cch[icha] / cm->wmolat[iatc]);
    ntheta[icha] =   (cm->nch[icha] / cm->wmolat[iatn])
                   / (cm->cch[icha] / cm->wmolat[iatc]);
    nomega[icha] =   (cm->sch[icha] / cm->wmolat[iats])
                   / (cm->cch[icha] / cm->wmolat[iatc]);
  }

  // Composition of char under the form
  // CH(noxgamma)O(noxdelta)N(noxkappa)S(noxzeta)

  for (int icha = 0; icha < n_coals; icha++) {

    // Mass fraction of nitrogen of coal (only one computation at start)
    nnch[icha] =   (ntheta[icha] * cm->wmolat[iatn])
                 / (               cm->wmolat[iatc]
                    + nalpha[icha]*cm->wmolat[iath]
                    + nbeta[icha] *cm->wmolat[iato]
                    + ntheta[icha]*cm->wmolat[iatn]
                    + nomega[icha]*cm->wmolat[iats]);

    /* Test on the char composition. If all coefficients are zero,
       they need to be computed automatically. */

    if ((  fabs(cm->cck[icha]) + fabs(cm->hck[icha]) +  fabs(cm->ock[icha])
         + fabs(cm->nck[icha]) + fabs(cm->sck[icha])) <= 0.) {

      // Mass fraction of nitrogen in  char i (function of Yi)
      nnckle[icha] = cm->repnle[icha] * nnch[icha] / (1.-cm->y1ch[icha]);
      nncklo[icha] = cm->repnlo[icha] * nnch[icha] / (1.-cm->y2ch[icha]);

      // In case of release of HCN in heterogeneaou combustion, the mass
      // fraction of hydrogen in char i is (function of Yi):
      nhckle[icha] =   cm->repnck[icha] * nnckle[icha]
                     * cm->wmolat[iath] / cm->wmolat[iatn];
      nhcklo[icha] =   cm->repnck[icha] * nncklo[icha]
                     * cm->wmolat[iath] / cm->wmolat[iatn];

      // Mass fraction of carbon of char (function of Yi)
      ncckle[icha]  = 1. - nnckle[icha] - nhckle[icha];
      nccklo[icha]  = 1. - nncklo[icha] - nhcklo[icha];

      //  Coefficients (function of Yi)
      nkapp1[icha] =   (cm->wmolat[iatc] - cm->wmolat[iatc]*ncckle[icha])
                     / (ncckle[icha]*(  cm->repnck[icha]*cm->wmolat[iath]
                                      +                  cm->wmolat[iatn]));
      nkapp2[icha] =   (cm->wmolat[iatc] - cm->wmolat[iatc]*nccklo[icha])
                    / (nccklo[icha]*(  cm->repnck[icha]*cm->wmolat[iath]
                                     +                  cm->wmolat[iatn]));

      ngama1[icha] = cm->repnck[icha] * nkapp1[icha];
      ngama2[icha] = cm->repnck[icha] * nkapp2[icha];

      ndelt1[icha] = 0.;
      ndelt2[icha] = 0.;

      nzeta1 [icha] = 0.;
      nzeta2 [icha] = 0.;

    }

    // Composition defined by user
    else {

      // Ensure the carbon coefficient is not zero.
      if (fabs(cm->cck[icha]) <= 0.) {
        bft_error
          (__FILE__, __LINE__, 0,
           _(" %s: pulverized coal model:\n\n"
             "Wrong elementary coke balance for coal %d.\n"
             "\n"
             " Check the setup parameters."), __func__, icha+1);
      }

      ngama1[icha] =   (cm->hck[icha] / cm->wmolat[iath])
                     / (cm->cck[icha] / cm->wmolat[iatc]);
      ngama2[icha] = ngama1[icha];

      ndelt1[icha] =   (cm->ock[icha] / cm->wmolat[iato])
                     / (cm->cck[icha] / cm->wmolat[iatc]);
      ndelt2[icha] = ndelt1[icha];

      nkapp1[icha] =   (cm->nck[icha] / cm->wmolat[iatn])
                     / (cm->cck[icha] / cm->wmolat[iatc]);
      nkapp2[icha] = nkapp1[icha];

      nzeta1[icha] =   (cm->sck[icha] / cm->wmolat[iats])
                     / (cm->cck[icha] / cm->wmolat[iatc]);
      nzeta2[icha] = nzeta1[icha];

    }

    // Coefficients of first pyrolysis reaction.

    double  sm[8], solu[8], mat[8][8];

    mat[0][0] = -ngama1[icha];
    mat[1][0] = -ngama1[icha];
    mat[2][0] = 2.;
    mat[3][0] = 2.;
    mat[4][0] = 1.-ngama1[icha];
    mat[5][0] = 3.;
    mat[6][0] = 1.;
    mat[7][0] = 0.;

    mat[0][1] = -ndelt1[icha];
    mat[1][1] = 1.-ndelt1[icha];
    mat[2][1] = 1.;
    mat[3][1] = 0.;
    mat[4][1] = -ndelt1[icha];
    mat[5][1] = 0.;
    mat[6][1] = 0.;
    mat[7][1] = 0.;

    mat[0][2] = -nkapp1[icha];
    mat[1][2] = -nkapp1[icha];
    mat[2][2] = 0.;
    mat[3][2] = 0.;
    mat[4][2] = 1.-nkapp1[icha];
    mat[5][2] = 1.;
    mat[6][2] = 0.;
    mat[7][2] = 0.;

    mat[0][3] = 0.;
    mat[1][3] = 0.;
    mat[2][3] = 0.;
    mat[3][3] = 0.;
    mat[4][3] = cm->crepn1[icha][1];
    mat[5][3] = -cm->crepn1[icha][0];
    mat[6][3] = 0.;
    mat[7][3] = 0.;

    mat[0][4] = -nzeta1[icha];
    mat[1][4] = -nzeta1[icha];
    mat[2][4] = 0.;
    mat[3][4] = 1.;
    mat[4][4] = -nzeta1[icha];
    mat[5][4] = 0.;
    mat[6][4] = 0.;
    mat[7][4] = 0.;

    mat[0][5] = cm->wmolat[iatc];
    mat[1][5] = cm->wmole[ico];
    mat[2][5] = cm->wmole[ih2o];
    mat[3][5] = cm->wmole[ih2s];
    mat[4][5] = cm->wmole[ihcn];
    mat[5][5] = cm->wmole[inh3];
    mat[6][5] = cm->wmolat[iath];
    mat[7][5] = -(cm->wmolat[iatc] + nalpha[icha]*cm->wmolat[iath]
                                   + nbeta[icha] *cm->wmolat[iato]
                                   + ntheta[icha]*cm->wmolat[iatn]
                                   + nomega[icha]*cm->wmolat[iats]);

    sm[0] = nalpha[icha]-ngama1[icha];
    sm[1] = nbeta[icha] -ndelt1[icha];
    sm[2] = ntheta[icha]-nkapp1[icha];
    sm[3] = 0.;
    sm[4] = nomega[icha]-nzeta1[icha];
    sm[5] = 0.;

    // Complete matrix and RHS basd on value of iy1ch

    if (cm->iy1ch[icha] == 0) {
      nchx1[icha] = 4.;

      mat[0][6] = -nchx1[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = 0.;
      sm[7] = 0.;
    }
    else if (cm->iy1ch[icha] == 1) {
      mat[0][6] = 0.;
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 0.;
      mat[7][6] = 1.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = cm->y1ch[icha];
      sm[7] = 0.;
    }
    else if (cm->iy1ch[icha] == 2) {
      nchx1[icha] = 4.;

      mat[0][6] = -nchx1[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 0.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 1.;

      sm[6] = 0.;
      sm[7] = cm->y1ch[icha];
    }

    int ierr = _coal_solve_matrix(8, (double *)mat, sm, solu);

    if (ierr == 0) {
      noxa1[icha] = solu[0];
      noxb1[icha] = solu[1];
      noxc1[icha] = solu[2];
      noxd1[icha] = solu[3];
      noxe1[icha] = solu[4];
      noxf1[icha] = solu[5];
      if (cm->iy1ch[icha] == 0) {
        // ny1ch[icha] = solu[7];
      }
      else if (cm->iy1ch[icha] == 1)
        nchx1[icha] = solu[6]/solu[0];
    }
    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "Non invertible matrix for volatile matter computation for coal %d.\n"
           "\n"
           " Check the setup parameters."), __func__, icha+1);
    }

    // Molar mass of CHX1
    wchx1c[icha] = cm->wmolat[iatc] + nchx1[icha]*cm->wmolat[iath];

    // Mass fractions of HCN, NH3, CHx1 in light volatile matter.
    cm->yhcnle[icha] =    (noxe1[icha]*cm->wmole[ihcn])
                       /  (  noxa1[icha]*(              cm->wmolat[iatc]
                                          + nchx1[icha]*cm->wmolat[iath])
                           + noxb1[icha]*cm->wmole[ico]
                           + noxc1[icha]*cm->wmole[ih2o]
                           + noxd1[icha]*cm->wmole[ih2s]
                           + noxe1[icha]*cm->wmole[ihcn]
                           + noxf1[icha]*cm->wmole[inh3]);

    cm->ynh3le[icha] =   (noxf1[icha]*cm->wmole[inh3])
                       / (  noxa1[icha]*(              cm->wmolat[iatc]
                                         + nchx1[icha]*cm->wmolat[iath])
                          + noxb1[icha]*cm->wmole[ico]
                          + noxc1[icha]*cm->wmole[ih2o]
                          + noxd1[icha]*cm->wmole[ih2s]
                          + noxe1[icha]*cm->wmole[ihcn]
                          + noxf1[icha]*cm->wmole[inh3]);

    ychxle[icha] =  (noxa1[icha]*(              cm->wmolat[iatc]
                                  + nchx1[icha]*cm->wmolat[iath]))
                  / (  noxa1[icha]*(              cm->wmolat[iatc]
                                    + nchx1[icha]*cm->wmolat[iath])
                     + noxb1[icha]*cm->wmole[ico]
                     + noxc1[icha]*cm->wmole[ih2o]
                     + noxd1[icha]*cm->wmole[ih2s]
                     + noxe1[icha]*cm->wmole[ihcn]
                     + noxf1[icha]*cm->wmole[inh3]);

    // Coefficients of heterogeneous reaction (Char1)
    noxj1[icha] = ngama1[icha];
    noxi1[icha] = 1. - noxj1[icha];
    noxk1[icha] = nkapp1[icha] - noxj1[icha];
    // noxh1[icha] = (noxi1[icha] + noxk1[icha])/2.;

    // Mass fractions massiques of HCN, NO and CO per kg Char1.
    ycoch1[icha] =   (noxi1[icha]*cm->wmole[ico])
                   / (  noxi1[icha]*cm->wmole[ico]
                      + noxj1[icha]*cm->wmole[ihcn]
                      + noxk1[icha]*3.e-2);
    cm->yhcnc1[icha] =   (noxj1[icha]*cm->wmole[ihcn])
                       / (  noxi1[icha]*cm->wmole[ico]
                          + noxj1[icha]*cm->wmole[ihcn]
                          + noxk1[icha]*3.e-2);
    cm->ynoch1[icha] =   (noxk1[icha]*3.e-2)
                       / (  noxi1[icha]*cm->wmole[ico]
                          + noxj1[icha]*cm->wmole[ihcn]
                          + noxk1[icha]*3.e-2);

    // Coefficients of second pyrolysis reaction.

    mat[0][0] = -ngama2[icha];
    mat[1][0] = -ngama2[icha];
    mat[2][0] = 2.;
    mat[3][0] = 2.;
    mat[4][0] = 1.-ngama2[icha];
    mat[5][0] = 3.;
    mat[6][0] = 1.;
    mat[7][0] = 0.;

    mat[0][1] = -ndelt2[icha];
    mat[1][1] = 1.-ndelt2[icha];
    mat[2][1] = 1.;
    mat[3][1] = 0.;
    mat[4][1] = -ndelt2[icha];
    mat[5][1] = 0.;
    mat[6][1] = 0.;
    mat[7][1] = 0.;

    mat[0][2] = -nkapp2[icha];
    mat[1][2] = -nkapp2[icha];
    mat[2][2] = 0.;
    mat[3][2] = 0.;
    mat[4][2] = 1.-nkapp2[icha];
    mat[5][2] = 1.;
    mat[6][2] = 0.;
    mat[7][2] = 0.;

    mat[0][3] = 0.;
    mat[1][3] = 0.;
    mat[2][3] = 0.;
    mat[3][3] = 0.;
    mat[4][3] =  cm->crepn2[icha][1];
    mat[5][3] = -cm->crepn2[icha][0];
    mat[6][3] = 0.;
    mat[7][3] = 0.;

    mat[0][4] = -nzeta2[icha];
    mat[1][4] = -nzeta2[icha];
    mat[2][4] = 0.;
    mat[3][4] = 1.;
    mat[4][4] = -nzeta2[icha];
    mat[5][4] = 0.;
    mat[6][4] = 0.;
    mat[7][4] = 0.;

    mat[0][5] = cm->wmolat[iatc];
    mat[1][5] = cm->wmole[ico];
    mat[2][5] = cm->wmole[ih2o];
    mat[3][5] = cm->wmole[ih2s];
    mat[4][5] = cm->wmole[ihcn];
    mat[5][5] = cm->wmole[inh3];
    mat[6][5] = cm->wmolat[iath];
    mat[7][5] = -(               cm->wmolat[iatc]
                  + nalpha[icha]*cm->wmolat[iath]
                  + nbeta[icha] *cm->wmolat[iato]
                  + ntheta[icha]*cm->wmolat[iatn]
                  + nomega[icha]*cm->wmolat[iats]);

    sm[0] = nalpha[icha]-ngama2[icha];
    sm[1] = nbeta[icha] -ndelt2[icha];
    sm[2] = ntheta[icha]-nkapp2[icha];
    sm[3] = 0.;
    sm[4] = nomega[icha]-nzeta2[icha];
    sm[5] = 0.;

    // Complete matrix and RHQ based on value of iy2ch.

    if (cm->iy2ch[icha] == 0) {
      nchx2[icha] = 2.;

      mat[0][6] = -nchx2[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;
      sm[6] = 0.;
      sm[7] = 0.;
    }
    else if (cm->iy2ch[icha] == 1) {
      mat[0][6] = 0.;
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 0.;
      mat[7][6] = 1.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 1.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 0.;

      sm[6] = cm->y2ch[icha];
      sm[7] = 0.;
    }
    else if (cm->iy2ch[icha] == 2) {
      nchx2[icha] = 2.;

      mat[0][6] = -nchx2[icha];
      mat[1][6] = 0.;
      mat[2][6] = 0.;
      mat[3][6] = 0.;
      mat[4][6] = 0.;
      mat[5][6] = 0.;
      mat[6][6] = 1.;
      mat[7][6] = 0.;

      mat[0][7] = 0.;
      mat[1][7] = 0.;
      mat[2][7] = 0.;
      mat[3][7] = 0.;
      mat[4][7] = 0.;
      mat[5][7] = 0.;
      mat[6][7] = 0.;
      mat[7][7] = 1.;

      sm[6] = 0.;
      sm[7] = cm->y2ch[icha];
    }

    ierr = _coal_solve_matrix(8, (double *)mat, sm, solu);

    if (ierr == 0) {
      noxa2[icha] = solu[0];
      noxb2[icha] = solu[1];
      noxc2[icha] = solu[2];
      noxd2[icha] = solu[3];
      noxe2[icha] = solu[4];
      noxf2[icha] = solu[5];
      if (cm->iy2ch[icha] == 0) {
        // ny2ch[icha] = solu[7];
      }
      else if (cm->iy2ch[icha] == 1)
        nchx2[icha] = solu[6]/solu[0];
    }
    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "Non invertible matrix for volatile matter computation for coal %d.\n"
           "\n"
           " Check the setup parameters."), __func__, icha+1);
    }

    // Molar mass of CHX2
    wchx2c[icha] = cm->wmolat[iatc] + nchx2[icha]*cm->wmolat[iath];

    // Mass fractions of HCN,NH3,CHx2 in heavy volatile matter.
    cm->yhcnlo[icha] =   (noxe2[icha]*cm->wmole[ihcn])
                       / (  noxa2[icha]*(              cm->wmolat[iatc]
                                         + nchx2[icha]*cm->wmolat[iath])
                          + noxb2[icha]*cm->wmole[ico]
                          + noxc2[icha]*cm->wmole[ih2o]
                          + noxd2[icha]*cm->wmole[ih2s]
                          + noxe2[icha]*cm->wmole[ihcn]
                          + noxf2[icha]*cm->wmole[inh3]);

    cm->ynh3lo[icha] =   (noxf2[icha]*cm->wmole[inh3])
                       / (  noxa2[icha]*(              cm->wmolat[iatc]
                                         + nchx2[icha]*cm->wmolat[iath])
                          + noxb2[icha]*cm->wmole[ico]
                          + noxc2[icha]*cm->wmole[ih2o]
                          + noxd2[icha]*cm->wmole[ih2s]
                          + noxe2[icha]*cm->wmole[ihcn]
                          + noxf2[icha]*cm->wmole[inh3]);

    ychxlo[icha] =   (noxa2[icha]*(              cm->wmolat[iatc]
                                   + nchx2[icha]*cm->wmolat[iath]))
                   /
                     (  noxa2[icha]*(              cm->wmolat[iatc]
                                     + nchx2[icha]*cm->wmolat[iath])
                      + noxb2[icha]*cm->wmole[ico]
                      + noxc2[icha]*cm->wmole[ih2o]
                      + noxd2[icha]*cm->wmole[ih2s]
                      + noxe2[icha]*cm->wmole[ihcn]
                      + noxf2[icha]*cm->wmole[inh3]);

    // Coefficients of products of heterogeneous reaction (Char2)
    noxj2[icha] = ngama2[icha];
    noxi2[icha] = 1. - noxj2[icha];
    noxk2[icha] = nkapp2[icha] - noxj2[icha];
    // noxh2[icha] = (noxi2[icha] + noxk2[icha])/2.;

    // Mass fractions of HCN, NO and CO per kg Char2.
    ycoch2[icha] =   (  noxi2[icha]*cm->wmole[ico])
                   / (  noxi2[icha]*cm->wmole[ico]
                      + noxj2[icha]*cm->wmole[ihcn]
                      + noxk2[icha]*3.e-2);
    cm->yhcnc2[icha] =   (noxj2[icha]*cm->wmole[ihcn])
                       / (  noxi2[icha]*cm->wmole[ico]
                          + noxj2[icha]*cm->wmole[ihcn]
                          + noxk2[icha]*3.e-2);
    cm->ynoch2[icha] =   (noxk2[icha]*3.e-2)
                       / (  noxi2[icha]*cm->wmole[ico]
                          + noxj2[icha]*cm->wmole[ihcn]
                          + noxk2[icha]*3.e-2);
  }

  // Molar masses of light and heavy volatile matter

  double ychx1t = 0., ychx2t = 0.;

  for (int icha = 0; icha < n_coals; icha++) {
    ychx1t = ychx1t + ychxle[icha];
    ychx2t = ychx2t + ychxlo[icha];
  }

  for (int icha = 0; icha < n_coals; icha++) {
    cm->wmchx1 += (ychxle[icha]/ychx1t * wchx1c[icha]);
    cm->wmchx2 += (ychxlo[icha]/ychx2t * wchx2c[icha]);
  }

  // Discretization of temperature

  cm->teno[0] =  300.;
  cm->teno[1] =  500.;
  cm->teno[2] = 1000.;
  cm->teno[3] = 1500.;
  cm->teno[4] = 2000.;
  cm->teno[5] = 2500.;
  cm->teno[6] = 3000.;

  // Reburning model
  // Kinetic Constant tabulated by temperature

  double  kf1[7]
    = {1.58e+12, 1.85e+12, 1.68e+12, 1.45e+12, 1.26e+12, 1.13e+12, 1.02e+12};

  double kf2[7]
    = {4.10e+13, 4.10e+13, 4.10e+13, 4.10e+13, 4.10e+13, 4.10e+13, 4.10e+13};

  double kf3[7]
    = {1.90e+13, 1.90e+13, 1.90e+13, 1.90e+13, 1.90e+13, 1.90e+13, 1.90e+13};

  double kf4[7]
    = {8.62e+04, 2.84e+08, 2.04e+11, 2.43e+12, 9.61e+12, 2.38e+13, 4.60e+13};

  double kr4[7]
    = {2.05e+04, 3.36e+07, 1.18e+10, 1.20e+11, 4.88e+11, 1.33e+12, 2.92e+12};

  double kf5[7]
    = {5.80e+07, 4.98e+09, 2.31e+11, 1.10e+12, 2.74e+12, 5.14e+12, 8.26e+12};

  double kr5[7]
    = {1.79e+01, 4.76e+05, 1.73e+09, 3.78e+10, 2.11e+11, 6.57e+11, 1.50e+12};

  double kf6[7]
    = {1.86e+14, 1.91e+14, 1.90e+14, 1.81e+14, 1.72e+14, 1.66e+14, 1.63e+14};

  double kr6[7]
    = {5.86e+11, 4.72e+12, 2.26e+13, 3.80e+13, 4.94e+13, 5.78e+13, 6.42e+13};

  double kf7[7]
    = {1.65e+14, 1.65e+14, 1.65e+14, 1.65e+14, 1.65e+14, 1.65e+14, 1.65e+14};

  double kr7[7]
    = {3.23e-03, 2.40e+04, 3.47e+09, 1.89e+11, 1.45e+12, 5.03e+12, 1.17e+13};

  // Constant similar to the one in FLUENT theory
  const double  pflue = 1.e-4;

  for (int ii = 0; ii < 7; ii++) {

    // Chi2 is computed as the kinetics in the FLUENT manual
    cm->chi2[ii] =  (4.52e5 * pow(cm->teno[ii], 1.6)
                            * exp(-80815./cs_physical_constants_r/cm->teno[ii]))
                  / (1.02e5 * pow(cm->teno[ii], 1.6)
                            * exp(-13802./cs_physical_constants_r/cm->teno[ii]));

    // JJ indicates the H/C ratio of fuel (3=CH4;2=CH3,etc.)

    for (int jj = 0; jj < 4; jj++) {

      if (jj == 1) {
        cm->ka[ii][jj] = kf1[ii] * (kr6[ii] / kf6[ii]) * 1.e-6*pflue;
        cm->kb[ii][jj] = kf2[ii] * 1.e-6*pflue;
        cm->kc[ii][jj] = kf3[ii] * (kf7[ii] / kr7[ii]) * 1.e-6*pflue;
      }
      else if (jj == 2) {
        cm->ka[ii][jj] = kf1[ii] * 1.e-6 * pflue;
        cm->kb[ii][jj] = kf2[ii] * kf6[ii] / kr6[ii]  * 1.e-6*pflue;
        cm->kc[ii][jj] = kf3[ii] * (kf6[ii]*kf7[ii]) / (kr6[ii]*kr7[ii])
                                 * 1.e-6*pflue;
      }
      else if (jj == 3) {
        cm->ka[ii][jj] = kf1[ii] * kf5[ii] / kr5[ii]  * 1.e-6*pflue;
        cm->kb[ii][jj] = kf2[ii] * (kf5[ii]*kf6[ii]) / (kr5[ii]*kr6[ii])
                                 * 1.e-6*pflue;
        cm->kc[ii][jj] = kf3[ii] *   (kf5[ii]*kf6[ii]*kf7[ii])
                                   / (kr5[ii]*kr6[ii]* kr7[ii]) * 1.e-6*pflue;
      }
      else if (jj == 4) {
        cm->ka[ii][jj] = kf1[ii] * (kf4[ii]*kf5[ii]) / (kr4[ii]*kr5[ii])
                                 * 1.e-6*pflue;
        cm->kb[ii][jj] = kf2[ii] *   (kf4[ii]*kf5[ii]*kf6[ii])
                                   / (kr4[ii]*kr5[ii]*kr6[ii]) * 1.e-6*pflue;
        cm->kc[ii][jj] = kf3[ii] *   (kf4[ii]*kf5[ii]*kf6[ii]*kf7[ii])
                                   / (kr4[ii]*kr5[ii]*kr6[ii]*kr7[ii])
                                 * 1.e-6*pflue;
      }

    }

  }

  /* Logging
     ======= */

  cs_log_separator(CS_LOG_DEFAULT);

  cs_log_printf(CS_LOG_DEFAULT,
                _("** Coals data summary **\n"
                  "------------------------\n\n"));

  const char coal_number_fmt[]
    = N_("----------------------\n"
         " Coal number: %d\n"
         "----------------------\n");

  const char normed_coal_fmt[]
    = N_("\n"
         "  Normalized coal composition\n"
         "  ---------------------------\n");

  const char normed_char_fmt[]
    = N_("\n"
         "  Normalized char composition\n"
         "  ---------------------------\n");

  const char normed_char_fmt_i[]
    = N_("\n"
         "  Normalized char composition (reaction %d)\n"
         "  ---------------------------------------\n");

  const char devol_model_cst_fmt[]
    = N_("\n"
         "  Devolatization model constants\n"
         "  ------------------------------\n");

  const char vol_matters_compo_fmt[]
    = N_("\n"
         "  Volatile materials composition\n"
         "  ------------------------------\n");

  const char normed_comp_fmt[]
    = N_("      C: %12.4e   H: %12.4e   O: %12.4e  N: %12.4e   S: %12.4e\n");

  const char vm_comp_fmt[]
    = N_("             CHx            CO            H2O             H2S"
         "            HCN            NH3\n"
         "            -------------------------------------------------"
         "-------------------------------------\n"
         "      VM1   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e\n"
         "      VM2   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e   %12.4e\n\n");

  for (int icha = 0; icha < n_coals; icha++) {
    // Coal number
    cs_log_printf(CS_LOG_DEFAULT, _(coal_number_fmt), icha+1);

    // Normed composition of coal
    cs_log_printf(CS_LOG_DEFAULT, _(normed_coal_fmt));
    double molcch = cm->cch[icha] / cm->cch[icha];
    cs_log_printf(CS_LOG_DEFAULT, _(normed_comp_fmt),
                  molcch, alpha[icha], beta[icha], theta[icha], omega[icha]);

    // Normed composition of char
    cs_log_printf(CS_LOG_DEFAULT, _(normed_char_fmt));
    double molcck = 1;
    if (cm->cck[icha] > 0.)
      molcck = cm->cck[icha] / cm->cck[icha];
    cs_log_printf(CS_LOG_DEFAULT, _(normed_comp_fmt),
                  molcck, gamma[icha], delta[icha], kappa[icha], zeta[icha]);

    // Devolatilization model options
    cs_log_printf(CS_LOG_DEFAULT, _(devol_model_cst_fmt));
    if (cm->iy1ch[icha] == 0)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 0: Y1 and Y2 computed\n"));
    else if (cm->iy1ch[icha]  == 1)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 1: Y1 and Y2 fixed\n"));
    else if  (cm->iy1ch[icha] == 2)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 2: CHX1 and CHX2 fixed\n"));

    // Y1 and Y2
    cs_log_printf(CS_LOG_DEFAULT, _("      Y1_ch = %12.4e   Y2_ch = %12.4e\n"),
                  cm->y1ch[icha], cm->y2ch[icha]);

    // Volatile materials Composition
    cs_log_printf(CS_LOG_DEFAULT, _(vol_matters_compo_fmt));
    cs_log_printf(CS_LOG_DEFAULT, _("      CHX1 = %12.4e   CHX2 = %12.4e\n"),
                  cm->chx1[icha], cm->chx2[icha]);
    if (cm->chx1[icha] <= 0. || cm->chx2[icha] <= 0) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           " Some values of CHX1 and CHX2 are negative or zero.\n"
           "\n"
           " Check the setup parameters."), __func__);
    }

    // VM composition
    cs_log_printf(CS_LOG_DEFAULT, _(vm_comp_fmt),
                  cm->a1[icha], cm->b1[icha], cm->c1[icha],
                  cm->d1[icha], cm->e1[icha], cm->f1[icha],
                  cm->a2[icha], cm->b2[icha], cm->c2[icha],
                  cm->d2[icha], cm->e2[icha], cm->f2[icha]);
  }

  // NOx model

  cs_log_printf(CS_LOG_DEFAULT,
                _("** NOx model **\n"
                  "---------------\n\n"));

  cs_log_printf(CS_LOG_DEFAULT,
                _("  ieqnox : %3d (NOx model activation)\n"),
                cm->ieqnox);
  if (cm->ieqnox)
    cs_log_printf(CS_LOG_DEFAULT,
                  _("  imdnox : %3d (NOx formation Features)\n"
                    "  irb    : %3d (Reburning model)\n"),
                  cm->imdnox, cm->irb);

  cs_log_printf(CS_LOG_DEFAULT, "\n");

  // Number of Coal and option Y
  for (int icha = 0; icha < n_coals; icha++) {
    // Coal number
    cs_log_printf(CS_LOG_DEFAULT, _(coal_number_fmt), icha+1);

    // Normed composition of coal
    cs_log_printf(CS_LOG_DEFAULT, _(normed_coal_fmt));
    double molcch = cm->cch[icha] / cm->cch[icha];
    cs_log_printf(CS_LOG_DEFAULT, _(normed_comp_fmt),
                  molcch, nalpha[icha], nbeta[icha],
                  ntheta[icha], nomega[icha]);

    // Composition of Char Reaction 1
    cs_log_printf(CS_LOG_DEFAULT, _(normed_char_fmt_i), 1);
    double molcck = 1.;
    if (cm->cck[icha] > 0.)
      molcck = cm->cck[icha] / cm->cck[icha];
    cs_log_printf(CS_LOG_DEFAULT, _(normed_comp_fmt),
                  molcck, ngama1[icha], ndelt1[icha],
                  nkapp1[icha], nzeta1[icha]);
    // Composition of Char Reaction 2
    cs_log_printf(CS_LOG_DEFAULT, _(normed_char_fmt_i), 2);
    cs_log_printf(CS_LOG_DEFAULT, _(normed_comp_fmt),
                  molcck, ngama2[icha], ndelt2[icha],
                  nkapp2[icha], nzeta2[icha]);

    // Devolatilization model options (remark: already output in previous loop)
    cs_log_printf(CS_LOG_DEFAULT, _(devol_model_cst_fmt));
    if (cm->iy1ch[icha] == 0)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 0: Y1 and Y2 computed\n"));
    else if (cm->iy1ch[icha]  == 1)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 1: Y1 and Y2 fixed\n"));
    else if  (cm->iy1ch[icha] == 2)
      cs_log_printf(CS_LOG_DEFAULT, _("      Option 2: CHX1 and CHX2 fixed\n"));

    // Y1 and Y2 (remark: already output in previous loop)
    cs_log_printf(CS_LOG_DEFAULT, _("      Y1_ch = %12.4e   Y2_ch = %12.4e\n"),
                  cm->y1ch[icha], cm->y2ch[icha]);

    // Volatile materials Composition
    cs_log_printf(CS_LOG_DEFAULT, _(vol_matters_compo_fmt));
    cs_log_printf(CS_LOG_DEFAULT, _("      CHX1 = %12.4e   CHX2 = %12.4e\n"),
                  nchx1[icha], nchx2[icha]);
    if (nchx1[icha] <= 0. || nchx2[icha] <= 0) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           " Some values of CHX1 and CHX2 are negative or zero.\n"
           "\n"
           " Check the setup parameters."), __func__);
    }

    // VM composition
    cs_log_printf(CS_LOG_DEFAULT, _(vm_comp_fmt),
                  noxa1[icha], noxb1[icha], noxc1[icha],
                  noxd1[icha], noxe1[icha], noxf1[icha],
                  noxa2[icha], noxb2[icha], noxc2[icha],
                  noxd2[icha], noxe2[icha], noxf2[icha]);

    // Products of heterogeneos reaction

    cs_log_printf(CS_LOG_DEFAULT,
                  _("  Heterogeneous reaction products\n"
                    "  -------------------------------\n"
                    "           CO             HCN            NO\n"
                    "      R1:  %12.4e   %12.4e   %12.4e\n"
                    "      R2:  %12.4e   %12.4e   %12.4e\n"),
                  ycoch1[icha], cm->yhcnc1[icha], cm->ynoch1[icha],
                  ycoch2[icha], cm->yhcnc2[icha], cm->ynoch2[icha]);
  }

  // Compute de DHdev

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "** Delta Hdev **\n"
                  "----------------\n\n"));

  cs_log_printf(CS_LOG_DEFAULT,
                _("      coal   light          heavy          PCI          \n"
                  "      --------------------------------------------------\n"));

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;
    int ick_icha = cm->ick[icha]-1;
    int ichx1c_icha = cm->ichx1c[icha]-1;
    int ichx2c_icha = cm->ichx2c[icha]-1;

    // Coal enthalpy

    double hchar = cm->ehsoli[0][ich_icha];
    double hcoke = cm->ehsoli[0][ick_icha];

    //  Enthalpy of volatile matter

    double ehvol1
      =   (  cm->a1[icha]*cm->ehgaze[0][ichx1c_icha]*cm->wmole[ichx1c_icha]
           + cm->b1[icha]*cm->ehgaze[0][ico]*cm->wmole[ico]
           + cm->c1[icha]*cm->ehgaze[0][ih2o]*cm->wmole[ih2o]
           + cm->d1[icha]*cm->ehgaze[0][ih2s]*cm->wmole[ih2s])
        / (  cm->a1[icha]*cm->wmole[ichx1c_icha]
           + cm->b1[icha]*cm->wmole[ico]
           + cm->c1[icha]*cm->wmole[ih2o]
           + cm->d1[icha]*cm->wmole[ih2s]);

    double ehvol2
      =   (  cm->a2[icha]*cm->ehgaze[0][ichx2c_icha]*cm->wmole[ichx2c_icha]
           + cm->b2[icha]*cm->ehgaze[0][ico]*cm->wmole[ico]
           + cm->c2[icha]*cm->ehgaze[0][ih2o]*cm->wmole[ih2o]
           + cm->d2[icha]*cm->ehgaze[0][ih2s]*cm->wmole[ih2s])
        / (  cm->a2[icha]*cm->wmole[ichx2c_icha]
           + cm->b2[icha]*cm->wmole[ico]
           + cm->c2[icha]*cm->wmole[ih2o]
           + cm->d2[icha]*cm->wmole[ih2s]);

    double dhvol1 = hchar -     cm->y1ch[icha] *ehvol1
                          - (1.-cm->y1ch[icha])*hcoke;
    double dhvol2 = hchar -     cm->y2ch[icha] *ehvol2
                          - (1.-cm->y2ch[icha])*hcoke;

    cs_log_printf(CS_LOG_DEFAULT,
                  _("      %2d   %12.4e   %12.4e   %12.4e\n"),
                  icha+1, dhvol1, dhvol2, cm->pcich[icha]);
  }

  // enthalpy/temperature law

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  "** Enthlapy/Temperature law **\n"
                  "------------------------------\n\n"));

  for (int icha = 0; icha < n_coals; icha++) {
    int ich_icha = cm->ich[icha]-1;
    int ick_icha = cm->ick[icha]-1;
    int iash_icha = cm->iash[icha]-1;
    int iwat_icha = cm->iwat[icha]-1;

    cs_log_printf(CS_LOG_DEFAULT, _(coal_number_fmt), icha+1);
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("      Temp          h_coal        h_coke        h_ash        h_wat\n"
         "      ------------------------------------------------------------------\n"));
    for (int ii = 0; ii < cm->npoc; ii++) {
      cs_log_printf(CS_LOG_DEFAULT,
                    _("      %12.4e  %12.4e  %12.4e  %12.4e  %12.4e\n"),
                    cm->thc[ii],
                    cm->ehsoli[ii][ich_icha],
                    cm->ehsoli[ii][ick_icha],
                    cm->ehsoli[ii][iash_icha],
                    cm->ehsoli[ii][iwat_icha]);
    }
  }

  // Compute the AiFj terms: number of moles of i per kg of j at start

  // double wmco  = cm->wmole[ico];
  double wmo2  = cm->wmole[io2];
  double wmco2 = cm->wmole[ico2];
  double wmh2o = cm->wmole[ih2o];
  double wmn2  = cm->wmole[in2];
  double wmc   = cm->wmolat[iatc];

  double dmf3 = (  cm->oxyo2[0] *wmo2  + cm->oxyn2[0] *wmn2
                 + cm->oxyh2o[0]*wmh2o + cm->oxyco2[0]*wmco2);
  if (dmf3 <= 0.) {
    bft_error
      (__FILE__, __LINE__, 0,
       _(" %s: pulverized coal model:\n\n"
         "The composition of oxydant 1 is incorrect\n"
         "  O2  : %g\n"
         "  N2  : %g\n"
         "  H2O : %g\n"
         "  CO2 : %g\n"
         "\n"
         " Check the setup parameters."), __func__,
       cm->oxyo2[0], cm->oxyn2[0], cm->oxyh2o[0], cm->oxyco2[0]);
  }

  cm->af3[io2]  = cm->oxyo2[0]  / dmf3;
  cm->af3[in2]  = cm->oxyn2[0]  / dmf3;
  cm->af3[ih2o] = cm->oxyh2o[0] / dmf3;
  cm->af3[ico2] = cm->oxyco2[0] / dmf3;

  if (cm->noxyd >= 2) {
    double dmf4  = (  cm->oxyo2[1]*wmo2   + cm->oxyn2[1]*wmn2
                    + cm->oxyh2o[1]*wmh2o + cm->oxyco2[1]*wmco2);
    if (dmf4 <= 0.) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "The composition of oxydant 2 is incorrect\n"
           "  O2  : %g\n"
           "  N2  : %g\n"
           "  H2O : %g\n"
           "  CO2 : %g\n"
           "\n"
           " Check the setup parameters."), __func__,
         cm->oxyo2[1], cm->oxyn2[1], cm->oxyh2o[1], cm->oxyco2[1]);
    }

    cm->af4[io2]  = cm->oxyo2[1]  / dmf4;
    cm->af4[in2]  = cm->oxyn2[1]  / dmf4;
    cm->af4[ih2o] = cm->oxyh2o[1] / dmf4;
    cm->af4[ico2] = cm->oxyco2[1] / dmf4;
  }

  if (cm->noxyd == 3) {
    double dmf5  = (  cm->oxyo2[2]*wmo2   + cm->oxyn2[2]*wmn2
                    + cm->oxyh2o[2]*wmh2o + cm->oxyco2[2]*wmco2);
    if (dmf5 <= 0.) {
      bft_error
        (__FILE__, __LINE__, 0,
         _(" %s: pulverized coal model:\n\n"
           "The composition of oxydant 3 is incorrect\n"
           "  O2  : %g\n"
           "  N2  : %g\n"
           "  H2O : %g\n"
           "  CO2 : %g\n"
           "\n"
           " Check the setup parameters."), __func__,
         cm->oxyo2[2], cm->oxyn2[2], cm->oxyh2o[2], cm->oxyco2[2]);

      cm->af5[io2]  = cm->oxyo2[2]  / dmf5;
      cm->af5[in2]  = cm->oxyn2[2]  / dmf5;
      cm->af5[ih2o] = cm->oxyh2o[2] / dmf5;
      cm->af5[ico2] = cm->oxyco2[2] / dmf5;
    }
  }

  // vapor
  cm->af6[ih2o]  = 1. / wmh2o;

  // coke per o2
  cm->af7[ico]   =  1.0 / wmc;
  cm->af7[io2]   = -0.5 / wmc;

  // coke per co2
  cm->af8[ico]   =  2. / wmc;
  cm->af8[ico2]  = -1. /wmc;

  // coke per h2o
  cm->af9[ih2o]  = -1. / wmc;
  cm->af9[ico]   =  1. / wmc;
  cm->af9[ihy]   =  1. / wmc;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
