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

subroutine coupbo &
!================

 ( itherm ,                                                       &
   cvcst  ,                                                       &
   hbord  , tbord  )

!===============================================================================
! Purpose:
! --------

! Send data relative to a SYRTHES coupling

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itherm           ! i  ! <-- ! 1: tparoi is a temperature ; 2: enthalpy       !
!                  !    !     ! 3: total energy (compressible)                 !
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! hbord(nfabor)    ! ra ! <-- ! coefficients d'echange aux bords               !
! tbord(nfabor)    ! ra ! <-- ! temperatures aux bords                         !
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
use cstphy
use mesh
use field
use optcal, only: iporos
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          itherm

double precision cvcst

double precision hbord(nfabor),tbord(nfabor)

! Local variables

integer          nbccou, inbcou, inbcoo, nbfcou, ifac, iloc, iel
integer          mode, flag
integer          iepsel, ifinwa
integer          issurf
double precision energ, cvt

integer, dimension(:), allocatable :: lfcou
double precision, dimension(:), allocatable :: tfluid, hparoi, wa
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cpro_cp, cpro_cv, cpro_rho
double precision, dimension(:), pointer :: porosi

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

! Map field arrays

!===============================================================================

! Initialize variables to avoid compiler warnings

iepsel = 0
ifinwa = 0

!===============================================================================
! SYRTHES coupling: compute fluid temperature and exchange coefficient
!===============================================================================

! Get number of coupling cases

call nbcsyr(nbccou)

flag = 0
do inbcou = 1, nbccou
  inbcoo = inbcou
  call tsursy(inbcoo, issurf)
  if (issurf.eq.1) then
    flag = 1
    exit
  endif
enddo

if (flag.eq.0) return

! Memory management to build arrays (large enough for all cases)

allocate(lfcou(nfabor))
allocate(tfluid(nfabor))
allocate(hparoi(nfabor))

! Prepare conversion to temperature for enthalpy or energy
! (check for surface couplings to make sure it is needed,
! exit earlier otherwise)

if (itherm.eq.2) then

  if (icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
  endif
  ! Temperature near boundary faces
  allocate(wa(nfabor))
  call b_h_to_t(tbord, wa)

else if (itherm.eq.3) then

  call field_get_val_v(ivarfl(iu), vel)
  if (icv.ge.0) then
    call field_get_val_s(icv, cpro_cv)
  endif
  call field_get_val_s(icrom, cpro_rho)
  ! Epsilon sup for perfect gas at cells
  allocate(wa(ncelet))
  call cs_cf_thermo_eps_sup(cpro_rho, wa, ncel)

endif

!---> Loop on couplings

do inbcou = 1, nbccou

  inbcoo = inbcou

  ! Test if this coupling is a surface coupling
  ! This is a surface coupling if issurf = 1

  call tsursy(inbcoo, issurf)

  if (issurf.eq.1) then

    mode = 0 ! Surface coupling

    ! Number of boundary faces per coupling case

    call nbesyr(inbcoo, mode, nbfcou)

    ! Loop on coupled faces to compute coefficients

    inbcoo = inbcou
    call leltsy(inbcoo, mode, lfcou)

    if (itherm.eq.1) then

      do iloc = 1, nbfcou

        ifac = lfcou(iloc)

        ! Saved fluid temperatures and exchange coefficients
        tfluid(iloc) = tbord(ifac)
        hparoi(iloc) = hbord(ifac)

      enddo

    ! In enthalpy formulation, transform to temperatures for SYRTHES
    !  To conserve flux Phi = (lambda/d     ) Delta T
    !                or Phi = (lambda/(d Cp)) Delta H
    !  recall      hbord = lambda/d.
    !  Conservation is not guaranteed, so we add a warning.

    else if (itherm.eq.2) then

      do iloc = 1, nbfcou

        ifac = lfcou(iloc)
        iel  = ifabor(ifac)
        tfluid(iloc) = wa(ifac)
        hparoi(iloc) = hbord(ifac)
      enddo

    else if (itherm.eq.3) then

      ! In energy formulation, transform to temperatures for SYRTHES
      !  To conserve flux Phi = (lambda/d     ) Delta T
      !                or Phi = (lambda/(d Cp)) Delta H
      !  Recall      hbord = lambda/ d
      !  Note that Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
      !  and  that Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
      !    (the difference is thus Cv Delta T)

      ! Modify temperature and exchange coefficient

      ! Compute e - CvT

      do iloc = 1, nbfcou
        ifac  = lfcou(iloc)
        iel   = ifabor(ifac)
        energ = tbord(ifac)
        cvt   = energ                                               &
                      -(0.5d0*(  vel(1,iel)**2                      &
                               + vel(2,iel)**2                      &
                               + vel(3,iel)**2)                     &
                        + wa(iel) )
        if (icv.ge.0) then
          tfluid(iloc) = cvt/cpro_cv(iel)
        else
          tfluid(iloc) = cvt/cvcst
        endif
        hparoi(iloc) = hbord(ifac)
      enddo

    endif

    ! Fluxes are multiplied by porosity if present.
    ! Here as the flux is expressed as h.(Tw-Tf), the exchange coefficient
    ! is multipled by the porosity.

    if (iporos.ge.1) then
      call field_get_val_s(ipori, porosi)
      do iloc = 1, nbfcou
        ifac = lfcou(iloc)
        iel = ifabor(ifac)
        hparoi(iloc) = hparoi(iloc)*porosi(iel)
      enddo
    endif

    ! Send fluid temperature and exchange coefficient

    inbcoo = inbcou
    call varsyo(inbcoo, mode, lfcou, tfluid, hparoi)

  endif ! This coupling is a surface coupling

enddo ! Loop on all syrthes couplings

! Free memory

if (itherm .gt. 1) deallocate(wa)

deallocate(hparoi)
deallocate(tfluid)
deallocate(lfcou)

!===============================================================================
! End of boundary couplings
!===============================================================================

return

end subroutine
