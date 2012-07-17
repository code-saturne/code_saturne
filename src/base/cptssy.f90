!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine cptssy &
!================

 ( nvar   , nscal  ,                                              &
   iscal  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! --------

! Compute the source term (implicit and/or explicit part) for a volume
! coupling with SYRTHES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iscal            ! i  ! <-- ! index number of the current scalar             !
! dt(ncelet)       ! ra ! <-- ! time-step                                      !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! crvexp           ! ra ! --> ! explicit part of the source term               !
! crvimp           ! ra ! --> ! implicit part of the source term               !
!__________________!____!_____!________________________________________________!
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
use cstnum
use mesh
use optcal

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

integer          nbccou, inbcou, inbcoo, ncecpl, iloc, iel
integer          mode  , isvol , ipccp , ivart
double precision tsexp, tsimp, cp, cecoef

integer, dimension(:), allocatable :: lcecpl
double precision, dimension(:), allocatable :: tfluid, ctbimp, ctbexp

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Get number of coupling cases and do something if there is
! at least one SYRTHES coupling

call nbcsyr(nbccou)
!==========

if (nbccou.lt.1) then
  return
endif

! Define source term only for the temperature scalar

if (iscalt.ne.iscal) then
  return
endif

if (iscalt.eq.iscal) then
  if (iscsth(iscalt).ne.-1) then
    write(nfecra, 1000)
  endif
endif

!===============================================================================
! 2. Compute implicit and explicit source terms
!===============================================================================

!---> Loop on couplings

do inbcou = 1, nbccou

  inbcoo = inbcou

  ! Test if this coupling is a surface coupling
  ! This is a surface coupling if isvol = 1

  call tvolsy(inbcoo, isvol)
  !==========

  if (isvol.eq.1) then

    mode = 1 ! Volume coupling
    ivart = isca(iscalt)

    ! Number of cells per coupling case

    call nbesyr(inbcoo, mode, ncecpl)
    !==========

    ! Memory management to build arrays
    allocate(lcecpl(ncecpl))
    allocate(tfluid(ncecpl))
    allocate(ctbimp(ncecpl))
    allocate(ctbexp(ncecpl))

    ! Get list of cells implied in this coupling

    call leltsy(inbcoo, mode, lcecpl)
    !==========

    ! Loop on coupled cells to initialize arrays

    do iloc = 1, ncecpl

      iel = lcecpl(iloc)
      tfluid(iloc) = rtpa(iel, ivart)
      ctbimp(iloc) = 0.0d0
      ctbexp(iloc) = 0.0d0

    enddo

    ! Compute implicit and explicit contribution to source terms

    call ctbvsy(inbcoo, tfluid, ctbimp, ctbexp)

    ! Loop on coupled cells to compute crvexp and crvimp

    do iloc = 1, ncecpl

      iel = lcecpl(iloc)

      ! Get Cp constant or variable

      if (icp.gt.0) then
        ipccp = ipproc(icp)
        cp = propce(iel, ipccp)
      else
        cp = cp0
      endif

      !FIXME
      cecoef = volume(iel)/cp
      tsexp = (ctbexp(iloc) - ctbimp(iloc)*tfluid(iloc))* cecoef
      tsimp = ctbimp(iloc) * cecoef

      crvexp(iel) = crvexp(iel) + tsexp
      crvimp(iel) = crvimp(iel) + tsimp

    enddo

    ! Free memory
    deallocate(tfluid)
    deallocate(ctbimp)
    deallocate(ctbexp)
    deallocate(lcecpl)

  endif ! This coupling is a surface coupling

enddo ! Loop on all syrthes couplings

!===============================================================================
! End of term source computation
!===============================================================================

return

! Formats

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE VOLUMIQUE SYRTHES :                ',/,&
'@      LA TEMPERATURE N''EST PAS CONFIGUREE EN °C.           ',/,&
'@    =========                                               ',/,&
'@    Le calcul continue.                                     ',/,&
'@                                                            ')

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@ @@ WARNING: SYRTHES VOLUME COUPLING:                       ',/,&
'@      THE TEMPERATURE IS NOT CONFIGURED IN °C.              ',/,&
'@    ========                                                ',/,&
'@    The calculation continue.                               ',/,&
'@                                                            ')

#endif

end subroutine
