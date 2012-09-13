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

!===============================================================================
! Function:
! ---------

!> \file vistnv.f90
!>
!> \brief This function computes the equivalent tensor viscosity at faces for
!> a 3x3 symetric tensor.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     imvisf        method to compute the viscosity at faces:
!>                               - 0: arithmetic
!>                               - 1: harmonic
!> \param[in]     w1            cell viscosity symmetric tensor
!> \param[in]     iwarnp        verbosity
!> \param[out]    weighf        inner face weight between cells i and j
!>                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
!>                                       {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
!>                              and
!>                              \f$ \frac{\vect{FJ} \cdot \tens{K}_\cellj}
!>                                       {\norm{\tens{K}_\cellj \cdot \vect{S}}^2} \f$
!> \param[out]    weighb        boundary face weight
!>                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
!>                                       {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
!> \param[out]    viscf         inner face tensor viscosity
!>                               (times surface divided distance)
!> \param[out]    viscb         inner face viscosity
!>                               (surface, must be consistent with flux BCs)
!_______________________________________________________________________________

subroutine vistnv &
!================
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use optcal, only: iporos
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          imvisf, iwarnp

double precision w1(6,ncelet)
double precision viscf(3,3,nfac), viscb(nfabor)

! Local variables

integer          ifac, iel, ii, jj, isou, jsou
double precision visci(3,3), viscj(3,3)
double precision distbf
double precision poroi, poroj, pnd

!===============================================================================

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call syntis(w1)
endif

! Without porosity
if (iporos.eq.0) then

  ! Arithmetic mean
  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      visci(1,1) = w1(1,ii)
      visci(2,2) = w1(2,ii)
      visci(3,3) = w1(3,ii)
      visci(1,2) = w1(4,ii)
      visci(2,1) = w1(4,ii)
      visci(2,3) = w1(5,ii)
      visci(3,2) = w1(5,ii)
      visci(1,3) = w1(6,ii)
      visci(3,1) = w1(6,ii)

      viscj(1,1) = w1(1,jj)
      viscj(2,2) = w1(2,jj)
      viscj(3,3) = w1(3,jj)
      viscj(1,2) = w1(4,jj)
      viscj(2,1) = w1(4,jj)
      viscj(2,3) = w1(5,jj)
      viscj(3,2) = w1(5,jj)
      viscj(1,3) = w1(6,jj)
      viscj(3,1) = w1(6,jj)

      do isou = 1, 3
        do jsou = 1, 3
          viscf(isou,jsou,ifac) = 0.5d0*(visci(isou,jsou)+viscj(isou,jsou)) &
                                * surfan(ifac)/dist(ifac)
        enddo
      enddo

    enddo

    do ifac = 1, nfabor
      ii = ifabor(ifac)
      viscb(ifac) = surfbn(ifac)
    enddo

  ! Harmonic mean
  else
!TODO
    call csexit(1)
  endif

! With porosity
else

  ! Arithmetic mean
  if (imvisf.eq.0) then

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      poroi = porosi(ii)
      poroj = porosi(jj)

      visci(1,1) = w1(1,ii)*poroi
      visci(2,2) = w1(2,ii)*poroi
      visci(3,3) = w1(3,ii)*poroi
      visci(1,2) = w1(4,ii)*poroi
      visci(2,1) = w1(4,ii)*poroi
      visci(2,3) = w1(5,ii)*poroi
      visci(3,2) = w1(5,ii)*poroi
      visci(1,3) = w1(6,ii)*poroi
      visci(3,1) = w1(6,ii)*poroi

      viscj(1,1) = w1(1,jj)*poroj
      viscj(2,2) = w1(2,jj)*poroj
      viscj(3,3) = w1(3,jj)*poroj
      viscj(1,2) = w1(4,jj)*poroj
      viscj(2,1) = w1(4,jj)*poroj
      viscj(2,3) = w1(5,jj)*poroj
      viscj(3,2) = w1(5,jj)*poroj
      viscj(1,3) = w1(6,jj)*poroj
      viscj(3,1) = w1(6,jj)*poroj

      do isou = 1, 3
        do jsou = 1, 3
          viscf(isou,jsou,ifac) = 0.5d0*(visci(isou,jsou)+viscj(isou,jsou)) &
                                * surfan(ifac)/dist(ifac)
        enddo
      enddo

    enddo

    do ifac = 1, nfabor
      ii = ifabor(ifac)
      viscb(ifac) = surfbn(ifac)*porosi(ii)
    enddo

  ! Harmonic mean
  else
!TODO
    call csexit(1)
  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
