!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine cou1di &
!================

 ( nfabor ,                                                       &
   isvtb  , icodcl ,                                              &
   rcodcl )

!===============================================================================

! FONCTION :
! ---------

! LECTURE DE DONNEES RELATIVES A UN COUPLAGE PAROI 1D

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! isvtb            ! e  ! <-- ! numero du scalaire couple                      !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar)   !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar)   !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use dimens, only: nvar
use entsor
use pointe
use field
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nfabor
integer          isvtb  , icodcl(nfabor,nvar)
double precision rcodcl(nfabor,nvar,3)

! Local variables


integer          ii , ivar
integer          ifac
integer          icldef
double precision, dimension(:), allocatable :: h_b
double precision, dimension(:), pointer :: b_temp

integer, dimension(:), pointer :: ifpt1d
double precision, dimension(:), pointer :: tppt1d

!===============================================================================

interface

  subroutine b_t_to_h(nlst, lstfac, t_b, h_b)

    use mesh, only: nfabor
    implicit none

    integer :: nlst
    integer, dimension(nlst) :: lstfac
    double precision, dimension(nfabor), intent(in) :: t_b
    double precision, dimension(nfabor), intent(out), target :: h_b

  end subroutine b_t_to_h

 end interface

!===============================================================================

! Get the 1D wall thermal module arrays
call cs_1d_wall_thermal_get_faces(ifpt1d)
call cs_1d_wall_thermal_get_temp(tppt1d)

! Update boundary temperature field for radiative transfer

if (iirayo.ge.1.and.nfpt1d.gt.0) then

  call field_get_val_s(itempb, b_temp)

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
      if (itpscl.eq.2) then
        b_temp(ifac) = tppt1d(ii) - tkelvi
      else
        b_temp(ifac) = tppt1d(ii)
      endif
    endif
  enddo

endif

! Sans specification, une face couplee est une face de type paroi FIXME, ca pose
! probleme si on a une couleur pour plusieurs types de face, paroi + autre
! Ces conditions peuvent etre ecrasees sur la face est couplee au flux radiatif dans
! cs_user_radiative_transfer_bcs.

icldef = 5

ivar = isca(isvtb)

do ii = 1, nfpt1d

   ifac = ifpt1d(ii)

   if ((icodcl(ifac,ivar).ne.1 .and.                           &
        icodcl(ifac,ivar).ne.5 .and.                           &
        icodcl(ifac,ivar).ne.6).and.                           &
       itypfb(ifac).eq.iparoi      ) icodcl(ifac,ivar) = icldef

   rcodcl(ifac,ivar,1) = tppt1d(ii)
   rcodcl(ifac,ivar,2) = rinfin
   rcodcl(ifac,ivar,3) = 0.d0

enddo

! Conversion eventuelle temperature -> enthalpie

if (isvtb.eq.iscalt .and. itherm.eq.2) then

  allocate(h_b(nfabor))

  do ii = 1, nfabor
    h_b(ii) = 0
  enddo

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    h_b(ifac) = tppt1d(ii)
  enddo

  call b_t_to_h(nfpt1d, ifpt1d, h_b, h_b)

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    rcodcl(ifac,ivar,1) = h_b(ifac)
  enddo

  deallocate(h_b)

endif

end subroutine


