/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_LAGRANG_TRACKING_H__
#define __CS_LAGRANG_TRACKING_H__

/*============================================================================
 * Utilitarian functions for the diphasic lagrangian module
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate cs_lagr_particle_set_t structure and initialize useful buffers.
 *
 * parameters:
 *  n_particles_max     -->  local max. number of particles
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagbeg, LAGBEG)(const cs_int_t   *n_particles_max,
                          const cs_int_t   *iphyla,
                          const cs_int_t   *nvls,
                          const cs_int_t   *nbclst);

/*----------------------------------------------------------------------------
 * Get variables and parameters associated to each particles and keep it in
 * a new structure
 *
 * parameters:
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (prtget, PRTGET)(const cs_int_t   *nbpmax,  /* n_particles max. */
                          const cs_int_t   *nbpart,  /* number of current particles */
                          const cs_real_t  *dnbpar,  /* particle total weight */
                          cs_int_t          liste[],
                          cs_int_t         *nbvis,
                          const cs_real_t   ettp[],
                          const cs_real_t   ettpa[],
                          const cs_int_t    itepa[],
                          const cs_real_t   tepa[],
                          const cs_int_t    ibord[],
                          const cs_int_t    indep[],
                          const cs_int_t   *jisor,
                          const cs_int_t   *jrpoi,
                          const cs_int_t   *jrtsp,
                          const cs_int_t   *jdp,
                          const cs_int_t   *jmp,
                          const cs_int_t   *jxp,
                          const cs_int_t   *jyp,
                          const cs_int_t   *jzp,
                          const cs_int_t   *jup,
                          const cs_int_t   *jvp,
                          const cs_int_t   *jwp,
                          const cs_int_t   *juf,
                          const cs_int_t   *jvf,
                          const cs_int_t   *jwf,
                          const cs_int_t   *jtaux,
                          const cs_int_t   *jryplu,
                          const cs_int_t   *jdfac,
                          const cs_int_t   *jimark,
                          cs_int_t         *idepst
);

/*----------------------------------------------------------------------------
 * Put variables and parameters associated to each particles into FORTRAN
 * arrays.
 *
 * parameters:
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (prtput, PRTPUT)(const cs_int_t   *nbpmax,  /* n_particles max. */
                          cs_int_t         *nbpart,  /* number of current particles */
                          cs_real_t        *dnbpar,  /* particle total weight */
                          cs_int_t         *nbpout,  /* number of outgoing particles */
                          cs_real_t        *dnbpou,  /* outgoing particle total weight */
                          cs_int_t         *nbperr,  /* number of failed particles */
                          cs_real_t        *dnbper,  /* failed particles total weight */
                          cs_int_t         *nbpdep,  /* number of depositing particles during the timestep*/
                          cs_real_t        *dnbdep,  /* depositing particles total weight during the timestep */
                          cs_int_t          liste[],
                          cs_int_t         *nbvis,
                          cs_real_t         ettp[],
                          cs_real_t         ettpa[],
                          cs_int_t          itepa[],
                          cs_real_t         tepa[],
                          cs_int_t          ibord[],
                          const cs_int_t   *jisor,
                          const cs_int_t   *jrpoi,
                          const cs_int_t   *jrtsp,
                          const cs_int_t   *jdp,
                          const cs_int_t   *jmp,
                          const cs_int_t   *jxp,
                          const cs_int_t   *jyp,
                          const cs_int_t   *jzp,
                          const cs_int_t   *jup,
                          const cs_int_t   *jvp,
                          const cs_int_t   *jwp,
                          const cs_int_t   *juf,
                          const cs_int_t   *jvf,
                          const cs_int_t   *jwf,
                          const cs_int_t   *jtaux,
                          const cs_int_t   *jryplu,
                          const cs_int_t   *jdfac,
                          const cs_int_t   *jimark,
                          cs_int_t         *idepst
);

/*----------------------------------------------------------------------------
 * Get variables and parameters associated to each particles and keep it in
 * a new structure
 *
 * parameters:
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (getbdy, GETBDY)(const cs_int_t    *nflagm,
                          const cs_int_t    *nfrlag,
                          const cs_int_t    *injcon,
                          const cs_int_t     ilflag[],
                          const cs_int_t     iusncl[],
                          const cs_int_t     iusclb[],
                          const cs_int_t     iusmoy[],
                          const cs_real_t    deblag[],
                          const cs_int_t     ifrlag[]);

/*----------------------------------------------------------------------------
 * Displacement of particles.
 *
 * parameters:
 *  p_n_particles     <->  pointer to the number of particles
 *  scheme_order    -->  current order of the scheme used for Lagragian
 *
 *
 *
 *----------------------------------------------------------------------------*/

void
CS_PROCF (dplprt, DPLPRT)(cs_int_t        *p_n_particles,
                          cs_real_t       *p_parts_weight,
                          cs_int_t        *p_scheme_order,
                          cs_real_t        boundary_stat[],
                          const cs_int_t  *iensi3,
                          const cs_int_t  *nvisbr,
                          const cs_int_t  *inbr,
                          const cs_int_t  *inbrbd,
                          const cs_int_t  *iflm,
                          const cs_int_t  *iflmbd,
                          const cs_int_t  *iang,
                          const cs_int_t  *iangbd,
                          const cs_int_t  *ivit,
                          const cs_int_t  *ivitbd,
                          const cs_int_t  *nusbor,
                          cs_int_t         iusb[],
                          cs_real_t        visc_length[],
                          cs_real_t        dlgeo[],
                          cs_real_t        rtp[],
                          const cs_int_t  *iu,
                          const cs_int_t  *iv,
                          const cs_int_t  *iw,
                          cs_int_t        *idepst);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Delete cs_lagr_particle_set_t structure and delete other useful buffers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_destroy(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGRANG_TRACKING_H__ */
