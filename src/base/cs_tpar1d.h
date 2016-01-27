#ifndef __CS_TPAR1D_H__
#define __CS_TPAR1D_H__

/*============================================================================
 * Modelling the thermal wall with 1D approach
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create the 1D mesh for each face and initialize the temperature
 *
 * Fortran interface:
 *
 * SUBROUTINE  MAIT1D
 * ******************
 *
 * INTEGER          NFPT1D         : <-  : number of coupled faces
 * INTEGER          NPPT1D(NFPT1D) : <-  : number of mesh points for each face
 * DOUBLE PRECISION EPPT1D(NFPT1D) : <-  : wall thickness for each face
 * DOUBLE PRECISION RGPT1D(NFPT1D) : <-  : mesh geometric ratio for each face
 * DOUBLE PRECISION TPPT1D(NFPT1D) : <-  : temperature initizalition value
 *----------------------------------------------------------------------------*/

void CS_PROCF (mait1d,MAIT1D)
(
 cs_int_t   *nf,
 cs_int_t    n[],
 cs_real_t   e[],
 cs_real_t   r[],
 cs_real_t   tp[]
);

/*----------------------------------------------------------------------------
 * Solve the 1D equation for a given face
 *
 * Fortran interface:
 *
 * SUBROUTINE  TPAR1D
 * ******************
 *
 * INTEGER          II     : <-  : face number
 * INTEGER          ICLT1D : <-  : type of exterior boundary condition
 * DOUBLE PRECISION TBORD  : <-  : fluid temperature at the boundary
 * DOUBLE PRECISION HBORD  : <-  : exchange coefficient for the fluid
 * DOUBLE PRECISION QINC   : <-  : incident radiative flux at the boundary
 *                         :     : at the boundary
 * DOUBLE PRECISION EPS    : <-  : emissivity (epsilon)
 * DOUBLE PRECISION TET1D  : <-  : temperature on the exterior boundary
 *                         :     : (Dirichlet boundary condition)
 * DOUBLE PRECISION HET1D  : <-  : exchange coefficient on the exterior wall
 * DOUBLE PRECISION FET1D  : <-  : flux on the exterior wall
 *                         :     : (Neumann boundary condition)
 * DOUBLE PRECISION LAMT1D : <-  : conductivity (lambda)
 * DOUBLE PRECISION RCPT1D : <-  : rho*Cp product
 * DOUBLE PRECISION DTPT1D : <-> : time-step for the solid resolution
 * DOUBLE PRECISION TPPT1D : <-> : physical temperature at the fluid/solid
 *                         :     : interface
 *----------------------------------------------------------------------------*/

void CS_PROCF (tpar1d,TPAR1D)
(
 cs_int_t *ii,
 cs_int_t *icdcle,
 cs_real_t *tf,
 cs_real_t *hf,
 cs_real_t *qinc,
 cs_real_t *eps,
 cs_real_t *te,
 cs_real_t *he,
 cs_real_t *fe,
 cs_real_t *lb,
 cs_real_t *rocp,
 cs_real_t *dtf,
 cs_real_t *tp
);

/*----------------------------------------------------------------------------
 * Read the restart file of the 1D-wall thermal module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECT1D
 * *****************
 *
 * CHARACTER        NOMSUI : <-- : Name of the restart file
 * INTEGER          LNGNOM : <-- : Name length
 * INTEGER          NFPT1D : <-- : Number of coupled faces
 * INTEGER          NFPT1T : <-- : Total number of coupled faces
 * INTEGER          NMXT1D : <-- : Max number of points on the 1D meshes
 * INTEGER          NFABOR : <-- : Number of boundary faces
 * INTEGER          NPPT1D : <-- : Number of points of each face 1D-mesh
 * INTEGER          IFPT1D : <-- : Indirection array for 1D-module faces
 * DOUBLE PRECISION EPPT1D : <-- : Wall thickness of each face
 * DOUBLE PRECISION RGPT1D : <-- : Geometric reason associated to faces
 * DOUBLE PRECISION TPPT1D : --> : Wall temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (lect1d, LECT1D)
(
 const char       *const nomsui,
 const cs_int_t   *const lngnom,
 const cs_int_t   *const nfpt1d,
 const cs_int_t   *const nfpt1t,
 const cs_int_t   *const nmxt1d,
 const cs_int_t   *const nfabor,
 const cs_int_t   *const nppt1d,
 const cs_int_t   *const ifpt1d,
 const cs_real_t  *const eppt1d,
 const cs_real_t  *const rgpt1d,
       cs_real_t  *const tppt1d
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Write the restart file of the 1D-wall thermal module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECT1D
 * *****************
 *
 * CHARACTER        NOMSUI : <-- : Name of the restart file
 * INTEGER          LNGNOM : <-- : Name length
 * INTEGER          NFPT1D : <-- : Number of coupled faces
 * INTEGER          NMXT1D : <-- : Max number of points on the 1D meshes
 * INTEGER          NFABOR : <-- : Number of boundary faces
 * DOUBLE PRECISION TPPT1D : --> : Wall temperature
 * INTEGER          IFPT1D : <-- : Indirection array for 1D-module faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrt1d, ECRT1D)
(
 const char       *const nomsui,
 const cs_int_t   *const lngnom,
 const cs_int_t   *const nfpt1d,
 const cs_int_t   *const nmxt1d,
 const cs_int_t   *const nfabor,
 const cs_real_t  *const tppt1d,
 const cs_int_t   *const ifpt1d
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Free allocated memory
 *----------------------------------------------------------------------------*/

void CS_PROCF (lbrt1d, LBRT1D)(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TPAR1D_H__ */

