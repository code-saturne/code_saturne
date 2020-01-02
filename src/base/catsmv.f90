!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file catsmv.f90
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
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     isnexp        sources terms of treated phase extrapolation
!>                              indicator
!> \param[in]     icetsm        source mass cells pointer
!> \param[in]     itpsmp        mass source type for the working variable
!>                              (cf. \ref cs_user_mass_source_terms)
!> \param[in]     cell_f_vol    cells fluid volume
!> \param[in]     vela          variable value at time step beginning
!> \param[in]     smcelv        value of the variable associated to the mass
!>                              source; NOT INTERLEAVED
!> \param[in]     gamma         mass flow value
!> \param[in,out] tsexpv        explicit source term part linear in the variable
!> \param[in,out] tsimpv        associated value with \c tsexp
!>                              to be stored in the matrix
!> \param[out]    gavinj        explicit source term part independant
!>                              of the variable
!______________________________________________________________________________

subroutine catsmv &
 ( ncelet , ncel   , ncesmp , iterns , isnexp ,                   &
   icetsm , itpsmp ,                                              &
   cell_f_vol      , vela   , smcelv , gamma  ,                   &
   tsexpv , tsimpv , gavinj )

!===============================================================================
! Module files
!===============================================================================

implicit none

! Variables

integer          ncelet, ncel  , ncesmp, iterns, isnexp
integer          icetsm(ncesmp), itpsmp(ncesmp)
double precision cell_f_vol(ncelet)
double precision vela  (3,ncelet)
double precision smcelv(ncesmp,3), gamma (ncesmp)
double precision tsexpv(3,ncelet), tsimpv(3,3,ncelet), gavinj(3,ncelet)

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
    do isou = 1, 3
      gavinj(isou,iel) = 0.d0
    enddo
  enddo
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 3
        tsexpv(isou,iel) = tsexpv(isou,iel)-cell_f_vol(iel)*gamma(ii)*vela(isou,iel)
        gavinj(isou,iel) = cell_f_vol(iel)*gamma(ii) * smcelv(ii,isou)
      enddo
    endif
  enddo
endif

!     Sur la diagonale
if(isnexp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 3
        tsimpv(isou,isou,iel) = tsimpv(isou,isou,iel)+cell_f_vol(iel)*gamma(ii)
      enddo
    endif
  enddo
else
  do ii = 1, ncesmp
    iel = icetsm(ii)
    if (gamma(ii).gt.0.d0 .and. itpsmp(ii).eq.1) then
      do isou = 1, 3
        tsimpv(isou,isou,iel) = tsimpv(isou,isou,iel)+cell_f_vol(iel)*gamma(ii)
      enddo
    endif
  enddo
endif

return
end subroutine
