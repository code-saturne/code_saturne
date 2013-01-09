!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine catsma &
!================

 ( ncelet , ncel   , ncesmp , iterns , isnexp ,                   &
   thetv  ,                                                       &
   icetsm , itpsmp ,                                              &
   volume , pvara  , smcelp , gamma  ,                            &
   tsexp  , tsimp  , gapinj )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES TERMES SOURCES IMPLICITE ET EXPLICITE
!  VENANT DES SOURCES DE MASSE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iterns           ! e  ! <-- ! numero d'iteration sur navsto                  !
! isnexp           ! e  ! <-- ! indicateur pour extrapolation des              !
!                  !    !     !  termes sources de la phase traitee            !
! thetv            ! r  ! <-- ! theta sch. pour la variable                    !
!                  !    !     !    thetv  v(n+1) + (1-thetv) v(n)              !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itpsmp(ncesmp    ! te ! <-- ! type de source de masse pour la                !
!                  !    !     !  variable traitee (cf. ustsma)                 !
! volume(ncel)     ! tr ! <-- ! volume des cellules                            !
! pvara(ncel)      ! tr ! <-- ! valeur de la variable en debut de pas          !
!                  !    !     !  de temps                                      !
! smcelp(ncesmp    ! tr ! <-- ! valeur de la variable associee a la            !
!                  !    !     !  source de masse                               !
! gamma(ncesmp)    ! tr ! <-- ! valeur du flux de masse                        !
! tsexp(ncel)      ! tr ! <-- ! part du terme source explicite                 !
!                  !    !     !  lineaire par rapport a la variable            !
! tsimp(ncel)      ! tr ! <-- ! valeur associee a tsexp qui ira dans           !
!                  !    !     !  la matrice                                    !
! gapinj(ncel)     ! tr ! --> ! part du terme source explicite                 !
!                  !    !     !  independant de la variable                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

implicit none

! Variables

integer          ncelet, ncel  , ncesmp, iterns, isnexp
integer          icetsm(ncesmp), itpsmp(ncesmp)
double precision thetv
double precision volume(ncelet)
double precision pvara (ncelet)
double precision smcelp(ncesmp), gamma (ncesmp)
double precision tsexp (ncelet), tsimp (ncelet), gapinj(ncelet)

! Local variables

integer ii, iel

!===============================================================================

! Explication des tests GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1 :
!     Si on enleve de la matiere ou si on entre a la valeur de la
!       cellule alors l'equation de IVAR n'est pas modifiee
!     Sinon, on ajoute le terme source GAMMA*(f_i-f^(n+1))

!     Dans TSIMP, on ajoute le terme qui ira sur la diagonale,
!       soit Gamma (*theta en cas d'ordre 2)
!     Dans TSEXP on ajoute le terme correspondant du second membre
!       cad Gamma * Pvar (avec Pvar)
!     Dans GAPINJ on place le terme Gamma Pinj qui ira au second membre

!     La distinction entre TSEXP et W1 (qui vont finalement tous les
!       deux au second membre) sert pour l'ordre 2 en temps.

if(iterns.eq.1) then
  do iel = 1, ncel
    gapinj(iel) = 0.d0
  enddo
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      tsexp(iel) = tsexp(iel) - volume(iel)*gamma(ii) * pvara(iel)
      gapinj(iel) = volume(iel)*gamma(ii) * smcelp(ii)
    endif
  enddo
endif

!     Sur la diagonale
if(isnexp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      tsimp(iel) = tsimp(iel) + volume(iel)*gamma(ii)*thetv
    endif
  enddo
else
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      tsimp(iel) = tsimp(iel) + volume(iel)*gamma(ii)
    endif
  enddo
endif

return
end subroutine
