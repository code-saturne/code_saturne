!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine coupbi &
!================

 ( nfabor , nvar   , nscal  , nphas  ,                            &
   icodcl ,                                                       &
   rcodcl )

!===============================================================================
! Purpose:
! --------

! Read data relative to a SYRTHES coupling

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! icodcl           ! te ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> flux density                          !
!                  !    !     ! = 4   -> slip and u.n=0 (velocity)             !
!                  !    !     ! = 5   -> friction and u.n=0 (velocity)         !
!                  !    !     ! = 6   -> rugosity and u.n=0 (velocity)         !
!                  !    !     ! = 9   -> free inlet/outlet (velocity)          !
! rcodcl           ! tr ! --> ! boundary condition values                      !
!  (nfabor, nvar)  !    !     ! rcodcl(1) = dirichlet value                    !
!                  !    !     ! rcodcl(2) = exchange coefficient value         !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value (negative       !
!                  !    !     !  if gain) in W/m2 or rugosity height (m)       !
!                  !    !     !  if icodcl=6                                   !
!                  !    !     ! for velocities (vistl+visct)*gradu             !
!                  !    !     ! for pressure              dt*gradp             !
!                  !    !     ! for scalars                                    !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
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
use entsor
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer          nfabor, nvar, nscal, nphas
integer          icodcl(nfabor,nvar)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ll, nbccou, inbcou, inbcoo, nbfcou
integer          ifac, iloc, iscal , iphas
integer          icldef
integer          mode
double precision temper, enthal

integer, dimension(:), allocatable :: lfcou
double precision, dimension(:), allocatable :: thpar

!===============================================================================

!===============================================================================
! SYRTHES coupling: get wall temperature
!===============================================================================

! Get number of coupling cases

call nbcsyr (nbccou)
!==========

!---> Loop on couplings: get "tparoi" array for each coupling and apply
!                        matching boundary condition.

do inbcou = 1, nbccou

  ! Number of boundary faces per coupling case
  inbcoo = inbcou
  call nbfsyr (inbcoo, nbfcou)
  !==========

  ! Memory management to receive array
  allocate(lfcou(nbfcou))
  allocate(thpar(nbfcou))

  ! Read wall temperature and interpolate if necessary.
  call varsyi (inbcou, thpar(1))
  !==========

  ! Prescribe wall temperature
  inbcoo = inbcou
  call lfasyr(inbcoo, lfcou(1))
  !==========

  ! Default condition
  icldef = 5

  do iscal = 1, nscal

    if (icpsyr(iscal).eq.1) then

      ! For scalars coupled with SYRTHES, prescribe a Dirichlet
      ! condition at coupled faces.
      ! For the time being, pass here only once, as only one scalar is
      ! coupled with SYRTHES.
      ! For the compressible module, solve in energy, but save the
      ! temperature separately, for BC's to be clearer.

      ll = isca(iscal)
      if (ippmod(icompf).ge.0) then
        iphas = 1
        if (iscal.eq.ienerg(iphas)) then
          ll = isca(itempk(iphas))
        else
          write(nfecra,9000)ienerg(iphas),iscal
          call csexit (1)
        endif
      endif


      do iloc = 1, nbfcou

        ifac = lfcou(iloc)

        if ((icodcl(ifac,ll) .ne. 1) .and.                        &
            (icodcl(ifac,ll) .ne. 5) .and.                        &
            (icodcl(ifac,ll) .ne. 6)) icodcl(ifac,ll) = icldef

        rcodcl(ifac,ll,1) = thpar(iloc)
        rcodcl(ifac,ll,2) = rinfin
        rcodcl(ifac,ll,3) = 0.d0

      enddo

      ! Possible temperature -> enthalpy conversion

      if (iscsth(iscal).eq.2) then

        do iloc = 1, nbfcou

          ifac = lfcou(iloc)

          temper = rcodcl(ifac,ll,1)
          mode   = -1
          call usthht(mode,enthal,temper)
          !==========
          rcodcl(ifac,ll,1) = enthal

        enddo

      endif

    endif

  enddo

  deallocate(thpar)
  deallocate(lfcou)

enddo

!===============================================================================
! End of boundary couplings
!===============================================================================

return

! Formats

#if defined(_CS_LANG_FR)

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU COUPLAGE SYRTHES              ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le module compressible, seul le scalaire ', i10       ,/,&
'@    peut etre couple a SYRTHES. Ici, on cherche a coupler   ',/,&
'@    le scalaire ', i10                                       ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN SYRTHES COUPLING                      ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  With the compressible module, only the scalar ', i10       ,/,&
'@    may be coupled with SYRTHES. Here, one tries to couple  ',/,&
'@    with the scalar ', i10                                   ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine
