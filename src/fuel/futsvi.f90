!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine futsvi &
!================

 ( ncelet , ncel   , numtra ,                                     &
   rtp    , propce , volume ,                                     &
   smbrs  , rovsdt ,                                              &
   w1     )

!===============================================================================
! FONCTION :
! --------
! ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
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
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant precedent)                !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
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
use fuincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , numtra

double precision rtp(ncelet,*), propce(ncelet,*)
double precision volume(ncelet)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision w1(ncelet)

! Local variables

integer          iel    , icla
integer          ipcrom
integer          ixckcl , ixnpcl , ipcght , ipcdia

double precision gamvap , gamhet  , d2s3
double precision fvap(4), fsd(4)  , fsh(4) , ftrac(4)
double precision fhet(4)
double precision diamht , diamad

integer          ipcte1 , ipcte2
double precision t2mt1

!===============================================================================

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! ---- W1 = X1

do iel = 1, ncel
  w1(iel) = 1.d0
enddo
do icla = 1, nclafu
  do iel = 1, ncel
    w1(iel) = w1(iel) - rtp(iel,isca(iyfol(icla)))
  enddo
enddo

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. CALCULS DES BILANS EXPLICITE ET IMPLICITE
!===============================================================================


ipcrom = ipproc(irom)
ixckcl = isca(ifhtf)
ipcte1 = ipproc(itemp1)

do icla = 1, nclafu

  ixnpcl = isca(ing(icla))
  ipcght = ipproc(igmhtf(icla))
  ipcdia = ipproc(idiam3(icla))
  ipcte2 = ipproc(itemp3(icla))

  do iel = 1, ncel

    if ( rtp(iel,ixnpcl) .gt. epsifl ) then

! --> Calculs preliminaires

      ftrac(1) = rtp(iel,isca(ifvap)) / w1(iel)
      ftrac(3) = rtp(iel,isca(ifhtf)) / w1(iel)
      ftrac(4) = 1.d0-ftrac(1)-ftrac(3)

! --> Prise en compte des TS interfaciaux relatifs
!     au phenomene d'evaporation

!           GAMVAP = PROPCE(IEL,IPPROC(IGMEVA(ICLA)))*PROPCE(IEL,IPCROM)

      t2mt1  = propce(iel,ipcte2)-propce(iel,ipcte1)
      gamvap = -propce(iel,ipproc(igmeva(icla)))*t2mt1            &
               *propce(iel,ipcrom)

      fsd(1) = 0.d0
      fsd(3) = 0.d0
      fsd(4) = 0.d0

      fvap(1) = 1.d0
      fvap(3) = 0.d0
      fvap(4) = 0.d0

      if ( (fsd(numtra)-ftrac(numtra))                            &
          *(2.d0*fvap(numtra)-fsd(numtra)-ftrac(numtra))          &
              .gt. epsifl ) then
        smbrs(iel)  = smbrs(iel) + gamvap                         &
              *volume(iel)*(fsd(numtra)-ftrac(numtra))            &
              *(2.d0*fvap(numtra)-fsd(numtra)-ftrac(numtra))
      endif

! --> Prise en compte des TS interfaciaux relatifs
!     au phenomene de combustion heterogene

      gamhet = -propce(iel,ipcrom)*propce(iel,ipcght)

      diamht = ( ( rtp(iel,isca(iyfol(icla)))                     &
                 /(rtp(iel,isca(ing(icla)))*rho0fl)               &
                  -pi*(dinikf(icla)**3)*xinkf/6.d0  )             &
                   *6.d0/(pi*(1.d0-xinkf)) )
      diamht = abs(diamht)**(1.d0/3.d0)

      diamad = propce(iel,ipcdia) / dinikf(icla)

      fsh(3) = 1.d0
      if ( diamad.gt.epsifl ) then
        fsh(3) = 1.d0 - (1.d0-ftrac(3))                           &
                * exp( propce(iel,ipcght)                         &
                      /( 2.d0*pi*2.77d-4*diamht                   &
                        *rtp(iel,ixnpcl)*propce(iel,ipcrom) ) )
      endif

      fsh(1) = ftrac(1) * (1.d0-fsh(3))/(1.d0-ftrac(3))
      fsh(4) = ftrac(4) * (1.d0-fsh(3))/(1.d0-ftrac(3))

      fhet(1) = zero
      fhet(3) = 1.d0
      fhet(4) = zero

      if (  (fsh(numtra)-ftrac(numtra))                           &
           *(2.d0*fhet(numtra)-fsh(numtra)-ftrac(numtra))         &
             .gt. epsifl ) then
        smbrs(iel)  = smbrs(iel) + gamhet                         &
              *volume(iel)*(fsh(numtra)-ftrac(numtra))            &
              *(2.d0*fhet(numtra)-fsh(numtra)-ftrac(numtra))
      endif

    endif

  enddo

enddo

!----
! FIN
!----

return

end subroutine
