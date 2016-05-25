#ifndef __CS_RAD_TRANSFER_H__
#define __CS_RAD_TRANSFER_H__

/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

typedef struct {

  /*! Radiative transfer model:
    - 0: none
    - 1: DOM
    - 2: P1 */
  int      iirayo;

  /*! Phase which radiates (bulk by default, but may be coal class or fuel
      droplets phase) */
  int      nrphas;

  /*! Verbosity level for the calculation of the wall temperatures */
  int      iimpar;

  /*! Verbosity level in the listing concerning the solution of
    the radiative transfer equation:
    - 0: no display
    - 1: standard
    - 2: complete */
  int      iimlum;

  /*! When gas or coal combustion is activated, \ref imodak indicates whether
    the absorption coefficient shall be calculated ``automatically'' (=1)
    or read from the data file (=0) */
  int      imodak;

  /*! ADF model:
    - 0 no ADF model
    - 1 ADF model with 8 intervals of wave length
    - 2 ADF model with 50 intervals of wave length */
  int      imoadf;

  /*! P1 model transparency warnings counter */
  int iwrp1t;

  /*! FSCK model:
    - 0 no FSCK model
    - 1 FSCK model activated */
  int  imfsck;

  /*! For the P-1 model, percentage of cells for which we allow the optical
      thickness to exceed unity, thish this should be avoided.
      (more precisely, where \f$ KL \f$ is lower than 1, where \f$ K \f$ is
      the absorption coefficient of the medium and \f$ L \f$ is a
      characteristic length of the domain). */
  double  xnp1mx;

  /*! Indicates the method used to calculate the radiative source term:
    - 0: semi-analytic calculation (compulsory with transparent media)
    - 1: conservative calculation
    - 2: semi-analytic calculation corrected in order to be globally
         conservative */
  int  idiver;

  /*! Index of the quadrature and number of directions for a single octant.

      Sn quadrature (n(n+2) directions)
      - 1: S4 (24 directions)
      - 2: S6 (48 directions)
      - 3: S8 (80 directions)

      Tn quadrature (8n^2 directions)
      - 4: T2 (32 directions)
      - 5: T4 (128 directions)
      - 6: Tn (8*ndirec^2 directions) */
  int  i_quadrature;

  //> Parameter assiociated to the Tn
  int      ndirec;

  //> For the Tn quadrature, \ref ndirec squared
  int      ndirs;

  //--> directions of angular values of the quadrature sx, sy, sz
  //    and weight of the solid angle associated

  cs_real_3_t  *sxyz;
  cs_real_t    *angsol;

  /*! Indicates whether the radiation variables should be initialized */
  int  restart;

  //> Period of the radiation module.
  //> The radiation module is called every \ref nfreqr time steps (more precisely,
  //> every time \ref optcal::ntcabs "ntcabs" is a multiple of \ref nfreqr).
  //> Also, in order to have proper initialization of the variables, whatever
  //> the value of \ref nfreqr, the radiation module is called at
  //> the first time step of a calculation (restart or not).
  //> Useful if and only if the radiation module is activated}
  int      nfreqr;

  //> Spectral radiation models (ADF and FSCK)
  //> Number of ETRs to solve
  int      nwsgg;
  //> Weights of the Gaussian quadrature
  cs_real_t *wq;

  //--> Informations sur les zones frontieres

  // NBZRDM Nombre max. de  zones frontieres
  // NOZRDM Numero max. des zones frontieres

  int      nbzrdm;

  int      nozrdm;

  // NZFRAD Nombre de zones de bord (sur le proc courant)
  // ILZRAY Liste des numeros de zone de bord (du proc courant)
  // NOZARM Numero de zone de bord atteint max
  //   exemple zones 1 4 2 : NZFRAD=3,NOZARM=4

  int      nozarm;
  int      nzfrad;
  int     *ilzrad;

  //--> Types de condition pour les temperatures de paroi :
  //       ITPIMP Profil de temperature imposee
  //       IPGRNO Parois grises ou noires
  //       IPREFL Parois reflechissante
  //       IFGRNO Flux de conduction impose dans la paroi
  //                   ET paroi non reflechissante (EPS non nul)
  //       IFREFL Flux de conduction impose dans la paroi
  //                   ET paroi reflechissante     (EPS = 0)
  //       ITPT1D Resolution de l'equation de la chaleur (module tp1d)

  int      itpimp, ipgrno, iprefl, ifgrno, ifrefl, itpt1d;

} cs_rad_transfer_params_t;

extern cs_rad_transfer_params_t *cs_glob_rad_transfer_params;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize radiative transfer module.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_H__ */
