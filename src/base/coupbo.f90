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

subroutine coupbo &
!================

 ( ncv    , ientha ,                                              &
   cvcst  , cv     ,                                              &
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
! ncv              ! i  ! <-- ! dimension de cv (ncelet ou 1)                  !
! ientha           ! i  ! <-- ! 1 si tparoi est une enthalpie                  !
!                  ! i  ! <-- ! 2 si tparoi est une energie                    !
!                  !    !     !    (compressible)                              !
! cvcst            ! r  ! <-- ! chaleur specifique si constante                !
! cv(ncv)          ! ra ! <-- ! chaleur specifique si variable                 !
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

!===============================================================================

implicit none

! Arguments

integer          ncv    , ientha

double precision cvcst

double precision cv(ncv)
double precision hbord(nfabor),tbord(nfabor)

! Local variables

integer          nbccou, inbcou, inbcoo, nbfcou, ifac, iloc, iel
integer          mode, flag
integer          iepsel, iepsfa, igamag, ixmasm, ifinwa
integer          issurf
double precision enthal, temper, energ, cvt

integer, dimension(:), allocatable :: lfcou
double precision, dimension(:), allocatable :: tfluid, hparoi, wa
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cpro_cp

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
if(icp.gt.0 .and. ientha.eq.1) then
  call field_get_val_s(iprpfl(icp), cpro_cp)
endif

!===============================================================================

! Initialize variables to avoid compiler warnings

iepsel = 0
iepsfa = 0
igamag = 0
ixmasm = 0
ifinwa = 0

!===============================================================================
! SYRTHES coupling: compute fluid temperature and exchange coefficient
!===============================================================================

! Get number of coupling cases

call nbcsyr(nbccou)
!==========

!---> Loop on couplings

do inbcou = 1, nbccou

  inbcoo = inbcou

  ! Test if this coupling is a surface coupling
  ! This is a surface coupling if issurf = 1

  call tsursy(inbcoo, issurf)
  !==========

  if (issurf.eq.1) then

    mode = 0 ! Surface coupling

    ! Number of boundary faces per coupling case

    call nbesyr(inbcoo, mode, nbfcou)
    !==========

    ! Memory management to build arrays
    allocate(lfcou(nbfcou))
    allocate(tfluid(nbfcou))
    allocate(hparoi(nbfcou))

    ! Compressible: coupling with energy

    if (ientha .eq. 2) then
      iepsel = 1
      iepsfa = iepsel + ncelet
      igamag = iepsfa + nfabor
      ixmasm = igamag + ncelet
      ifinwa = ixmasm + ncelet
      allocate(wa(ifinwa))
    endif

    ! Loop on coupled faces to compute coefficients

    inbcoo = inbcou
    call leltsy(inbcoo, mode, lfcou)
    !==========

    do iloc = 1, nbfcou

      ifac = lfcou(iloc)

      ! Saved fluid temperatures and exchange coefficients
      tfluid(iloc) = tbord(ifac)
      hparoi(iloc) = hbord(ifac)

    enddo

    ! In enthalpy formulation, transform to temperatures for SYRTHES
    !  To conserve flux Phi = (lambda/d     ) Delta T
    !                or Phi = (lambda/(d Cp)) Delta H
    !  we multiply hbord = lambda/(d Cp) by Cp in the adjacent cell.
    !  Conservation is not guaranteed, so we add a warning.

    if (ientha.eq.1) then

      write(nfecra,1000)
      flag = 1
      do iloc = 1, nbfcou
        ifac = lfcou(iloc)
        iel  = ifabor(ifac)
        enthal = tfluid(iloc)
        call usthht(flag, enthal, temper)
        !==========
        tfluid(iloc) = temper
        if (icp.gt.0) then
          hparoi(iloc) = hparoi(iloc)*cpro_cp(iel)
        else
          hparoi(iloc) = hparoi(iloc)*cp0
        endif
      enddo

    else if (ientha.eq.2) then

      ! In energy formulation, transform to temperatures for SYRTHES
      !  To conserve flux Phi = (lambda/d     ) Delta T
      !                or Phi = (lambda/(d Cp)) Delta H
      !  we multiply hbord = lambda/(d Cp) by Cv in the adjacent cell.
      !  Note that Ei = Cv Ti + 1/2 Ui*Ui + Epsilon_sup_i
      !  and  that Ep = Cv Tp + 1/2 Ui*Ui + Epsilon_sup_i
      !    (the difference is thus Cv Delta T)

      ! Modify temperature and exchange coefficient

      ! Compute e - CvT

      ! At cell centers
      call cf_thermo_eps_sup(wa(iepsel), ncel)
      !=====================

      ! At boundary faces centers
      call cf_thermo_eps_sup(wa(iepsfa), nfabor)
      !=====================

      do iloc = 1, nbfcou
        ifac  = lfcou(iloc)
        iel   = ifabor(ifac)
        energ = tfluid(iloc)
        cvt   = energ                                               &
                      -(0.5d0*(  vel(1,iel)**2                      &
                               + vel(2,iel)**2                      &
                               + vel(3,iel)**2)                     &
                        + wa(iepsel+iel-1)           )
        if (ncv.eq.ncelet) then
          tfluid(iloc) = cvt/cv(iel)
          hparoi(iloc) = hparoi(iloc)*cv(iel)
        else
          tfluid(iloc) = cvt/cvcst
          hparoi(iloc) = hparoi(iloc)*cvcst
        endif
      enddo

    endif

    ! Send fluid temperature and exchange coefficient

    inbcoo = inbcou
    call varsyo(inbcoo, mode, lfcou, tfluid, hparoi)
    !==========

    ! Free memory
    if (ientha .eq. 2) deallocate(wa)
    deallocate(hparoi)
    deallocate(tfluid)
    deallocate(lfcou)

  endif ! This coupling is a surface coupling

enddo ! Loop on all syrthes couplings

!===============================================================================
! End of boundary couplings
!===============================================================================

return

! Formats

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE SYRTHES AVEC CALCUL EN ENTHALPIE   ',/,&
'@    =========                                               ',/,&
'@      OPTION NON VALIDEE - CONTACTER L''EQUIPE DE DVPT      ',/,&
'@                                                            ')

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ WARNING: SYRTHES COUPLING WITH ENTHALPY CALCULATION     ',/,&
'@    ========                                                ',/,&
'@      OPTION NOT VALIDATED - CONTACT THE SUPPORT            ',/,&
'@                                                            ')

#endif

end subroutine
