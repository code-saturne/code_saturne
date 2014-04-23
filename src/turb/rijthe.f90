!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine rijthe &
!================

 ( nscal  ,                                                       &
   ivar   ,                                                       &
   rtpa   ,                                                       &
   gradro , smbr   )

!===============================================================================
! FONCTION :
! ----------

! TERMES DE GRAVITE
!   POUR Rij et EPSILON
! VAR  = R11 R22 R33 R12 R13 R23 EP

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! i  ! <-- ! variable number                                !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at previous time step)                       !
! gradro(ncelet,3) ! tr ! <-- ! tableau de travail pour grad rom               !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use dimens, only: nvar
use cstphy
use mesh

!===============================================================================

implicit none

! Arguments

integer          nscal
integer          ivar

double precision rtpa(ncelet,nflown:nvar)
double precision gradro(ncelet,3)
double precision smbr(ncelet)

! Local variables

integer          iel

double precision uns3, const, kseps, csttmp
double precision prdtur, r1t, r2t, r3t
double precision g11, g22, g33, g12, g13, g23, gkks3
double precision g11p, g22p, g33p
double precision phit11, phit22, phit33, phit12, phit13, phit23
double precision aa, bb


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! EBRSM
if (iturb.eq.32) then
  csttmp = cebmr6
else
  csttmp = crij3
endif

if(iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

const = -1.5d0*cmu/prdtur
uns3  = 1.d0/3.d0

!===============================================================================
! 2. TERMES POUR RIJ :
!      ROM*VOLUME*dRij/dt =
!                     ... + (Gij - CRIJ3*(Gij-Delta ij Gkk/3))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!===============================================================================


if     (ivar.eq.ir11) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit11 = -csttmp*(g11-gkks3)

    smbr(iel) = smbr(iel) + (g11+phit11)*volume(iel)

  enddo

elseif (ivar.eq.ir22) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit22 = -csttmp*(g22-gkks3)

    smbr(iel) = smbr(iel) + (g22+phit22)*volume(iel)

  enddo

elseif (ivar.eq.ir33) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit33 = -csttmp*(g33-gkks3)

    smbr(iel) = smbr(iel) + (g33+phit33)*volume(iel)

  enddo

elseif (ivar.eq.ir12) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g12 = const*kseps*     (r1t*gy+r2t*gx)

    phit12 = -csttmp* g12

    smbr(iel) = smbr(iel) + (g12+phit12)*volume(iel)

  enddo

elseif (ivar.eq.ir13) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g13 = const*kseps*     (r1t*gz+r3t*gx)

    phit13 = -csttmp* g13

    smbr(iel) = smbr(iel) + (g13+phit13)*volume(iel)

  enddo

elseif (ivar.eq.ir23) then

  do iel = 1, ncel

    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    kseps = (rtpa(iel,ir11)+rtpa(iel,ir22)+rtpa(iel,ir33))  &
           /(2.d0*rtpa(iel,iep))

    g23 = const*kseps*(r2t*gz+r3t*gy)

    phit23 = -csttmp* g23

    smbr(iel) = smbr(iel) + (g23+phit23)*volume(iel)

  enddo

!===============================================================================
! 3. TERMES POUR EPSILON :
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*(EPS/K)*MAX(0,(Gkk/2))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!            On simplifie (EPS/K) en notant
!                GijP = -(1.5 CMU/PRDTUR)         (Rit Gj + Rjt Gi)
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*        MAX(0,(GkkP/2))*VOLUME
!===============================================================================


elseif (ivar.eq.iep ) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11)*gradro(iel,1)                            &
        + rtpa(iel,ir12)*gradro(iel,2)                            &
        + rtpa(iel,ir13)*gradro(iel,3)
    r2t = rtpa(iel,ir12)*gradro(iel,1)                            &
        + rtpa(iel,ir22)*gradro(iel,2)                            &
        + rtpa(iel,ir23)*gradro(iel,3)
    r3t = rtpa(iel,ir13)*gradro(iel,1)                            &
        + rtpa(iel,ir23)*gradro(iel,2)                            &
        + rtpa(iel,ir33)*gradro(iel,3)

    g11p = const*2.d0*(r1t*gx)
    g22p = const*2.d0*(r2t*gy)
    g33p = const*2.d0*(r3t*gz)

    aa = 0.d0
    bb = 0.5d0*(g11p+g22p+g33p)
    smbr(iel) = smbr(iel) + ce1*max(aa,bb)*volume(iel)

  enddo

endif

return

end subroutine
