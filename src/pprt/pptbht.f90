!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine pptbht &
 ( ncoel  ,                                                       &
   nomcel , ehcoel , cpcoel , wmolce )

!===============================================================================
! FONCTION :
! --------

!         PHYSIQUES PARTICULIERES

!           CALCUL DE L'ENTHALPIE ET DU CP
!                    A PARTIR DE LA BANDE DE JANAF


! Arguments
!_______________.____._____.________________________________________________.
!    nom        !type!mode !                   role                         !
!_______________!____!_____!________________________________________________!
! ncoel         ! e  ! <-- ! nombre de const. elem.                         !
! nomcel(ngazem)! a  ! <-- ! nom des constituants elementaires              !
! ehcoel        ! tr !  <- ! enthalpie pour chaque constituant              !
! (ngazem,npot) !    !     !                elementaire                     !
! cpcoel        ! tr !  <- ! cp pour chaque constituant                     !
! (ngazem,npot) !    !     !                elementaire                     !
! wmolce(ngazem)! tr !  <- ! masse molaire de chaque constituant            !
!_______________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use ppincl

!===============================================================================

implicit none

! Arguments

integer          ncoel

character(len=12) :: nomcel(ngazem)

double precision ehcoel(ngazem,npot) , cpcoel(ngazem,npot)
double precision wmolce (ngazem)

! Local variables

integer  ii, jj, kk
character(len=13, kind=c_char) :: c_name_sp
character(len=13*ngazem, kind=c_char) :: c_nomcel

integer(c_int) :: c_ncoel, c_ngazem, c_npo

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_combustion_enthalpy_and_cp_from_janaf &
    (c_ncoel, c_ngazem, c_npo, c_nomcel, ehcoel, cpcoel, wmolce, th)  &
    bind(C, name='cs_combustion_enthalpy_and_cp_from_janaf')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: c_ncoel, c_ngazem, c_npo
    character(kind=c_char, len=1), dimension(*), intent(in) :: c_nomcel
    real(kind=c_double), dimension(*), intent(in) :: ehcoel, cpcoel, wmolce, th
  end subroutine cs_combustion_enthalpy_and_cp_from_janaf

end interface

!===============================================================================

! Transform name to C array of names

do ii = 1, 13*ngazem
  c_nomcel(ii:ii) = c_null_char
enddo

do ii = 1, ncoel
  c_name_sp = trim(nomcel(ii)) // c_null_char
  do jj = 1, len(nomcel(ii))
    kk = (ii-1)*13 + jj
    c_nomcel(kk:kk) = c_name_sp(jj:jj)
  enddo
enddo

! Other variables

c_ncoel = ncoel
c_ngazem = ngazem
c_npo = npo

call cs_combustion_enthalpy_and_cp_from_janaf &
  (c_ncoel, c_ngazem, c_npo, c_nomcel, ehcoel, cpcoel, wmolce, th)

end subroutine
