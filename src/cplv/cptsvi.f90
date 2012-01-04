!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine cptsvi &
!================

 ( ncelet , ncel   , numtra ,                                     &
   rtp    , propce , volume ,                                     &
   smbrs , rovsdt ,                                               &
   xf1m   , xf2m   ,                                              &
   w1     )

!===============================================================================
! FONCTION :
! --------
! ROUTINE PHYSIQUE PARTICULIERE : FLAMME CP
!  TERMES DE TRANSFERT DE MASSE ENTRE LA PHASE CONTINUE
!  ET LA PHASE DISPERSEE RELATIF A LA VARIANCE DES DIFFERENTS
!  TRACEURS (BILANS EXPLICITE ET IMPLICITE)

!  ATTENTION CETTE SUBROUTINE N'EST OPERATIONNELLE QUE
!              POUR NUMTRA = 3 ou 4

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! numtra           ! e  ! <-- ! numero du traceur concerne (1,2,3,4)           !
! rtp         ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant precedent)                !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! xf1m, xf2m       ! tr ! <-- ! somme de f1,f2,f3,f4 sur l'ensemble            !
!                  !    !     ! sur l'ensemble des charbons et                 !
!                  !    !     ! des classes                                    !
! w1(ncelet)       ! tr ! --- ! tableau de travail                             !
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
use cstphy
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , numtra

double precision rtp(ncelet,*), propce(ncelet,*)
double precision volume(ncelet)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision xf1m(ncelet)  , xf2m(ncelet)
double precision w1(ncelet)

! Local variables

integer          iel    , icla   , icha
integer          ipcrom , ixchcl , ixckcl , ixnpcl
integer          ipcgd1 , ipcgd2 , ipcght , ipcdia

double precision gamdev1 , gamdev2 , gmdevt , gamhet  , d2s3
double precision fdev(4) , fhet(4) , fsd(4)  , fsh(4) , ftrac(4)
double precision diamdv , diamht , x2 , xmp , diamad

!===============================================================================

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --> Calcul de X1 = 1. - SOMME(X2) dans W1

do iel = 1, ncel
  w1(iel) = 1.d0
enddo

do icla = 1, nclacp
  ixchcl = isca(ixch(icla))
  ixckcl = isca(ixck(icla))
  ixnpcl = isca(inp(icla ))
  do iel = 1, ncel
    x2 = rtp(iel,ixchcl)+rtp(iel,ixckcl)                          &
       + rtp(iel,ixnpcl)*xmash(icla)
    w1(iel) = w1(iel) - x2
  enddo
enddo

ipcrom = ipproc(irom)

! Calcul de F1 et F2 sur l'ensemble des charbons

do iel = 1, ncel
  xf1m(iel) = 0.d0
  xf2m(iel) = 0.d0
enddo
do icha = 1, ncharb
  do iel = 1, ncel
    xf1m(iel) = xf1m(iel) + rtp(iel,isca(if1m(icha)))
    xf2m(iel) = xf2m(iel) + rtp(iel,isca(if2m(icha)))
  enddo
enddo

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. CALCULS DES BILANS EXPLICITE ET IMPLICITE
!===============================================================================

do icla = 1, nclacp

  ixchcl = isca(ixch(icla))
  ixckcl = isca(ixck(icla))
  ixnpcl = isca(inp(icla ))
  ipcgd1 = ipproc(igmdv1(icla))
  ipcgd2 = ipproc(igmdv2(icla))
  ipcght = ipproc(igmhet(icla))
  ipcdia = ipproc(idiam2(icla))

  do iel = 1, ncel

! --> Calculs preliminaires

    ftrac(1) = xf1m(iel)           / w1(iel)
    ftrac(2) = xf2m(iel)           / w1(iel)
    ftrac(3) = rtp(iel,isca(if3m)) / w1(iel)
    ftrac(4) = 1.d0-ftrac(1)-ftrac(2)-ftrac(3)

    xmp = xmp0(icla)*rtp(iel,ixnpcl)
    x2  = rtp(iel,ixchcl) + rtp(iel,ixckcl) +                     &
          rtp(iel,ixnpcl)*xmash(icla)

    if ( xmp .gt. epsicp .and.                                    &
         x2  .gt. epsicp        ) then

! --> Prise en compte des TS interfaciaux relatifs
!     au phenomene de devolatilisation

      if (rtp(iel,ixchcl).gt.epsicp) then

!     On prend comme diametre le diametre initial de la particule
!     la devolatilisation ayant lieu a diametre constant

        diamdv = diam20(icla)
        fsd(1) = 1.d0 - (1.d0-ftrac(1))                           &
             * exp( ( rtp(iel,ixchcl)*propce(iel,ipcgd1)  )       &
                   /( 2.d0*pi*2.77d-4*diamdv                      &
                      *rtp(iel,ixnpcl)*propce(iel,ipcrom) ) )

        fsd(2) = 1.d0 - (1.d0-ftrac(2))                           &
             * exp( ( rtp(iel,ixchcl)*propce(iel,ipcgd2)  )       &
                   /( 2.d0*pi*2.77d-4*diamdv                      &
                      *rtp(iel,ixnpcl)*propce(iel,ipcrom) ) )

        fsd(3) = ftrac(3) * (1.d0 -   fsd(1) -   fsd(2))          &
                        / (1.d0 - ftrac(1) - ftrac(2))
        fsd(4) = ftrac(4) * (1.d0 -   fsd(1) -   fsd(2))          &
                        / (1.d0 - ftrac(1) - ftrac(2))

        gamdev1 = - propce(iel,ipcrom) * rtp(iel,ixchcl)          &
                  * propce(iel,ipcgd1)
        gamdev2 = - propce(iel,ipcrom) * rtp(iel,ixchcl)          &
                  * propce(iel,ipcgd2)
        gmdevt = gamdev1+gamdev2

        fdev(1) = gamdev1 / gmdevt
        fdev(2) = gamdev2 / gmdevt
        fdev(3) = zero
        fdev(4) = zero

        if (  (fsd(numtra)-ftrac(numtra))                         &
           *(2.d0*fdev(numtra)-fsd(numtra)-ftrac(numtra))         &
          .gt. epsicp ) then
          smbrs(iel) = smbrs(iel) + gmdevt                        &
                    *volume(iel)*(fsd(numtra)-ftrac(numtra))      &
                    *(2.d0*fdev(numtra)-fsd(numtra)-ftrac(numtra))
!                ROVSDT(IEL) = ROVSDT(IEL) + ZERO
        endif

      endif

! --> Prise en compte des TS interfaciaux relatifs
!     au phenomene de combustion heterogene

!     On prend comme diametre le diametre du coke

      diamht = propce(iel,ipcdia)
      diamad = propce(iel,ipcdia) / diam20(icla)
      fsh(3) = 1.d0
      if ( diamad.gt.epsicp ) then
        fsh(3) = 1.d0 - (1.d0-ftrac(3))                           &
               * exp( ( (rtp(iel,ixckcl)**d2s3)                   &
                        *propce(iel,ipcght)       )               &
                     /( 2.d0*pi*2.77d-4*diamht                    &
                        *rtp(iel,ixnpcl)*propce(iel,ipcrom) ) )
      endif

      fsh(1) = ftrac(1) * (1.d0-fsh(3))/(1.d0-ftrac(3))
      fsh(2) = ftrac(2) * (1.d0-fsh(3))/(1.d0-ftrac(3))
      fsh(4) = ftrac(4) * (1.d0-fsh(3))/(1.d0-ftrac(3))

      gamhet = - propce(iel,ipcrom) * (rtp(iel,ixckcl)**d2s3)     &
                 * propce(iel,ipcght)

      fhet(1) = zero
      fhet(2) = zero
      fhet(3) = 1.d0
      fhet(4) = zero

      if (  (fsh(numtra)-ftrac(numtra))                           &
           *(2.d0*fhet(numtra)-fsh(numtra)-ftrac(numtra))         &
         .gt. epsicp ) then
        smbrs(iel) = smbrs(iel) + gamhet                          &
                    *volume(iel)*(fsh(numtra)-ftrac(numtra))      &
                    *(2.d0*fhet(numtra)-fsh(numtra)-ftrac(numtra))
!              ROVSDT(IEL) = ROVSDT(IEL) + ZERO
      endif

    endif

  enddo

enddo


!----
! FIN
!----

return

end subroutine
