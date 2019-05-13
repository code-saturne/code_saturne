!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file catsmt.f90
!> \brief Compute explicit and implicit source terms coming from mass source.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     iterns        iteration number on Navier-Stoke
!> \param[in]     isnexp        sources terms of treated phasis extrapolation
!>                              indicator
!> \param[in]     icetsm        source mass cells pointer
!> \param[in]     itpsmp        mass source type for the working variable
!>                              (see \ref cs_user_mass_source_terms)
!> \param[in]     volume        cells volume
!> \param[in]     pvara         variable value at time step beginning
!> \param[in]     smcelp        value of the variable associated with mass source
!> \param[in]     gamma         flow mass value
!> \param[in,out] tsexp         explicit source term part linear in the variable
!> \param[in,out] tsimp         associated value with \c tsexp
!>                              to be stored in the matrix
!> \param[out]    gapinj        explicit source term part independant
!>                              of the variable
!______________________________________________________________________________

subroutine catsmt &
 ( ncelet , ncel   , ncesmp , iterns , isnexp ,                   &
   icetsm , itpsmp ,                                              &
   volume , pvara  , smcelp , gamma  ,                            &
   tsexp  , tsimp  , gapinj )


!===============================================================================
! Module files
!===============================================================================

implicit none

! Variables

integer          ncelet, ncel  , ncesmp, iterns, isnexp
integer          icetsm(ncesmp), itpsmp(ncesmp)
double precision volume(ncelet)
double precision pvara(6,ncelet)
double precision smcelp(ncesmp,6), gamma (ncesmp)
double precision tsexp(6,ncelet), tsimp(6,6,ncelet), gapinj(6,ncelet)

! Local variables

integer ii, iel, isou

!===============================================================================

! Explication des tests GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1 :
!     Si on enleve de la matiere ou si on entre a la valeur de la
!       cellule alors l'equation de IVAR n'est pas modifiee
!     Sinon, on ajoute le terme source GAMMA*(f_i-f^(n+1))

!     Dans TSIMP, on ajoute le terme qui ira sur la diagonale,
!       soit Gamma
!     Dans TSEXP on ajoute le terme correspondant du second membre
!       cad Gamma * Pvar (avec Pvar)
!     Dans GAPINJ on place le terme Gamma Pinj qui ira au second membre

!     La distinction entre TSEXP et W1 (qui vont finalement tous les
!       deux au second membre) sert pour l'ordre 2 en temps.

if(iterns.eq.1) then
  do iel = 1, ncel
    do isou = 1, 6
      gapinj(isou,iel) = 0.d0
    enddo
  enddo
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 6
        tsexp(isou,iel) = tsexp(isou,iel)-volume(iel)*gamma(ii)*pvara(isou,iel)
        gapinj(isou,iel) = volume(iel)*gamma(ii) * smcelp(ii,isou)
      enddo
    endif
  enddo
endif

!     Sur la diagonale
if(isnexp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 6
        tsimp(isou,isou,iel) = tsimp(isou,isou,iel)+volume(iel)*gamma(ii)
      enddo
    endif
  enddo
else
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 6
        tsimp(isou,isou,iel) = tsimp(isou,isou,iel)+volume(iel)*gamma(ii)
      enddo
    endif
  enddo
endif

return
end subroutine
