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

subroutine cou1do &
 ( cvcst ,  hbord  , tbord  )

!===============================================================================
! FONCTION :
! ---------

! ECRITURE DE DONNEES RELATIVES A UN COUPLAGE AVEC SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! hbord(nfabor)    ! ra ! <-> ! coefficients d'echange aux bords               !
! tbord(nfabor)    ! ra ! <-> ! temperatures aux bords                         !
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
use entsor
use optcal
use cstphy
use cstnum
use field
use parall
use period
use pointe
use mesh
use cs_cf_bindings
use cs_c_bindings
use radiat

!===============================================================================

implicit none

! Arguments

double precision cvcst
double precision hbord(nfabor),tbord(nfabor)

!     VARIABLES LOCALES

integer          iappel
integer          ifac, iel , ii

double precision energ, cvt

double precision, dimension(:), allocatable :: wa
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cpro_cp, cpro_cv, cpro_rho

integer, dimension(:), pointer :: ifpt1d

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine b_h_to_t(h_b, t_b)

    use mesh, only: nfabor
    implicit none

    double precision, dimension(nfabor), intent(in) :: h_b
    double precision, dimension(nfabor), intent(out), target :: t_b

  end subroutine b_h_to_t

 end interface

!===============================================================================

! Get the 1D wall thermal module arrays

call cs_1d_wall_thermal_get_faces(ifpt1d)

! Conversion to temperature for enthalpy or energy
! (check for surface couplings to make sure it is needed)

! In enthalpy formulation, transform to temperatures for SYRTHES
!  To conserve flux Phi = (lambda/d     ) Delta T
!                or Phi = (lambda/(d Cp)) Delta H
!  recall      hbord = lambda/d.

if (itherm.eq.2) then

  if (icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
  endif

  ! Temperature near boundary faces
  allocate(wa(nfabor))
  call b_h_to_t(tbord, wa)

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    iel  = ifabor(ifac)
    tbord(ifac) = wa(ifac)
  enddo

else if (itherm.eq.3) then

  call field_get_val_v(ivarfl(iu), vel)
  call field_get_val_s(icrom, cpro_rho)
  if (icv.ge.0) then
    call field_get_val_s(icv, cpro_cv)
  endif

  ! Epsilon sup for perfect gas at cells
  allocate(wa(ncelet))
  call cs_cf_thermo_eps_sup(cpro_rho, wa, ncel)

  do ii = 1, nfpt1d
    ifac  = ifpt1d(ii)
    iel   = ifabor(ifac)
    energ = tbord(ifac)
    cvt   = energ                                               &
                  -(0.5d0*(  vel(1,iel)**2                      &
                           + vel(2,iel)**2                      &
                           + vel(3,iel)**2)                     &
                    + wa(iel) )
    if (icv.ge.0) then
      tbord(ifac) = cvt/cpro_cv(iel)
    else
      tbord(ifac) = cvt/cvcst
    endif
  enddo

endif

!     Mise a jour des conditions aux limites externes du module 1D
iappel = 3

call cs_user_1d_wall_thermal(iappel, isuit1)

iappel = 3
call cs_1d_wall_thermal_check(iappel, isuit1)

! coupling with radiative module
if (iirayo.ge.1) then

  do ii = 1, nfpt1d

    ifac = ifpt1d(ii)

    ! FIXME pour gerer les faces qui ne sont pas des parois ou beps n'est pas
    ! renseigne, par exemple si une meme couleur est utilisee pour designer
    ! plusieurs faces, paroi + autre
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then

      call cs_1d_wall_thermal_solve(ii-1, tbord(ifac), hbord(ifac))

    endif

  enddo

! without coupling with radiative module
else

  do ii = 1, nfpt1d

    ifac = ifpt1d(ii)

    call cs_1d_wall_thermal_solve(ii-1, tbord(ifac), hbord(ifac))

  enddo

endif

if (itherm .gt. 1) deallocate(wa)

return
end subroutine
