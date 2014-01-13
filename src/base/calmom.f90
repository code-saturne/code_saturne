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

subroutine calmom &
!================

 ( ncel   , ncelet ,        &
   rtp    , dt     , propce )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA MOYENNE TEMPORELLE DE CORRELATIONS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,*   )    !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use optcal
use numvar

!===============================================================================

implicit none

! Arguments

integer          ncel   , ncelet
double precision rtp(ncelet,*) , dt(ncelet), propce(ncelet,*)

! Local variables

integer          ii, iel, imom, icmom, idtcm
integer          jmom, icumut

double precision, allocatable, dimension(:) :: travm

!===============================================================================

!     PHASE DE CUMUL DES TERMES
!     =========================

! Allocate a temporary array
allocate(travm(ncelet))

! Boucle et test sur les moments a calculer
do imom = 1, nbmomt
  if (ntcabs.ge.ntdmom(imom).and.ttcabs.ge.ttdmom(imom)) then

!   Position dans PROPCE du tableau de cumul des moments
    icmom = ipproc(icmome(imom))

!   Tableau de travail pour recueillir la correlation par cellule
    do iel = 1, ncel
      travm(iel) = 1.d0
    enddo

!   Calcul de la correlation (on suppose qu'il y a au moins une variable)
    do ii = 1, idgmom(imom)

!     Si la variable est dans RTP
      if(idfmom(ii,imom).gt.0) then
        do iel = 1, ncel
          travm(iel) = travm(iel) * rtp(iel,idfmom(ii,imom))
        enddo
!     Si la variable est dans PROPCE
      elseif(idfmom(ii,imom).lt.0) then
        do iel = 1, ncel
          travm(iel) = travm(iel) * propce(iel,ipproc(-idfmom(ii,imom)))
        enddo
      endif

    enddo

!   Incrementation du tableau de cumul des correlations
    do iel = 1, ncel
      propce(iel,icmom) = propce(iel,icmom) + dt(iel)*travm(iel)
    enddo

!   Incrementation du tableau de cumul du temps (ncel ou constante)

!     Si plusieurs moyennes partagent le meme cumul temporel,
!     il ne faut incrementer que pour une seule (la premiere)

!     On decide si on doit incrementer
    icumut = 1
    if (imom.gt.1) then
      do jmom = 1, imom-1
        if(idtmom(jmom).eq.idtmom(imom)) then
          icumut = 0
        endif
      enddo
    endif

!     On incremente
    if(icumut.eq.1) then
      if(idtmom(imom).gt.0) then
        idtcm = ipproc(icdtmo(idtmom(imom)))
        do iel = 1, ncel
          propce(iel,idtcm) = propce(iel,idtcm) + dt(iel)
        enddo
      elseif(idtmom(imom).lt.0) then
        idtcm = -idtmom(imom)
        dtcmom(idtcm) = dtcmom(idtcm) + dt(1)
      endif
    endif

  endif
enddo

! Free memory
deallocate(travm)

return

end subroutine
