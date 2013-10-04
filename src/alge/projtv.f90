!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file projtv.f90
!>
!> \brief This function projects the external source termes to the faces
!> in coherence with itrmav.f90 for the improved hydrostatic pressure
!> algorithm (iphydr=1).
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     init           indicator
!>                               - 1 initialize the mass flux to 0
!>                               - 0 otherwise
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     frcxt         body force creating the hydrostatic pressure
!> \param[in]     cofbfp        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
!> \param[in]     weighf        internal face weight between cells i j in case
!>                               of tensor diffusion
!> \param[in,out] flumas        mass flux at interior faces
!> \param[in,out] flumab        mass flux at boundary faces
!_______________________________________________________________________________

subroutine projtv &
 ( init   , nswrgp , ircflp ,                                     &
   nfecra ,                                                       &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
   viscf  , viscb  ,                                              &
   viscel ,                                                       &
   weighf ,                                                       &
   flumas , flumab )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use mesh

!===============================================================================

implicit none

! Arguments

integer          init
integer          nswrgp , ircflp
integer          nfecra

double precision frcxt(3, ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision viscel(6,ncelet)
double precision weighf(2,nfac)
double precision cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, i
double precision distbf,surfn
double precision diippf(3), djjppf(3)
double precision visci(3,3), viscj(3,3)
double precision fikdvi, fjkdvi

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

if (init.eq.1) then
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo

elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit(1)
endif

!===============================================================================
! 2. Update mass flux without reconstruction technics
!===============================================================================

if (nswrgp.le.1) then

  ! ---> Contribution from interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas(ifac) =  flumas(ifac)                                     &
         + viscf(ifac)*(                                             &
           (cdgfac(1,ifac)-xyzcen(1,ii))*frcxt(1, ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*frcxt(2, ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*frcxt(3, ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*frcxt(1, jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*frcxt(2, jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*frcxt(3, jj) )

  enddo

  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn             &
         *cofbfp(ifac)*(frcxt(1, ii)*surfbo(1,ifac)                  &
         +frcxt(2, ii)*surfbo(2,ifac)+frcxt(3, ii)*surfbo(3,ifac) )

  enddo

else
!===============================================================================
! 3. Update mass flux WITH reconstruction technics
!===============================================================================

  ! ---> Contribution from interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    ! Recompute II" and JJ"
    !----------------------

    visci(1,1) = viscel(1,ii)
    visci(2,2) = viscel(2,ii)
    visci(3,3) = viscel(3,ii)
    visci(1,2) = viscel(4,ii)
    visci(2,1) = viscel(4,ii)
    visci(2,3) = viscel(5,ii)
    visci(3,2) = viscel(5,ii)
    visci(1,3) = viscel(6,ii)
    visci(3,1) = viscel(6,ii)

    ! IF.Ki.S / ||Ki.S||^2
    fikdvi = weighf(1,ifac)

    ! II" = IF + FI"
    do i = 1, 3
      diippf(i) = cdgfac(i,ifac)-xyzcen(i,ii)          &
                - fikdvi*( visci(i,1)*surfac(1,ifac)   &
                         + visci(i,2)*surfac(2,ifac)   &
                         + visci(i,3)*surfac(3,ifac) )
    enddo

    viscj(1,1) = viscel(1,jj)
    viscj(2,2) = viscel(2,jj)
    viscj(3,3) = viscel(3,jj)
    viscj(1,2) = viscel(4,jj)
    viscj(2,1) = viscel(4,jj)
    viscj(2,3) = viscel(5,jj)
    viscj(3,2) = viscel(5,jj)
    viscj(1,3) = viscel(6,jj)
    viscj(3,1) = viscel(6,jj)

    ! FJ.Kj.S / ||Kj.S||^2
    fjkdvi = weighf(2,ifac)

    ! JJ" = JF + FJ"
    do i = 1, 3
      djjppf(i) = cdgfac(i,ifac)-xyzcen(i,jj)          &
                + fjkdvi*( viscj(i,1)*surfac(1,ifac)   &
                         + viscj(i,2)*surfac(2,ifac)   &
                         + viscj(i,3)*surfac(3,ifac) )
    enddo

    flumas(ifac) = flumas(ifac)                                               &
                 + viscf(ifac)*(                                              &
                                 frcxt(1, ii)*(cdgfac(1,ifac)-xyzcen(1,ii))   &
                               + frcxt(2, ii)*(cdgfac(2,ifac)-xyzcen(2,ii))   &
                               + frcxt(3, ii)*(cdgfac(3,ifac)-xyzcen(3,ii))   &
                               - frcxt(1, jj)*(cdgfac(1,ifac)-xyzcen(1,jj))   &
                               - frcxt(2, jj)*(cdgfac(2,ifac)-xyzcen(2,jj))   &
                               - frcxt(3, jj)*(cdgfac(3,ifac)-xyzcen(3,jj))   &
                               )                                              &
                 + viscf(ifac)*ircflp*(                                       &
                                      - frcxt(1, ii)*diippf(1)                &
                                      - frcxt(2, ii)*diippf(2)                &
                                      - frcxt(3, ii)*diippf(3)                &
                                      + frcxt(1, jj)*djjppf(1)                &
                                      + frcxt(2, jj)*djjppf(2)                &
                                      + frcxt(3, jj)*djjppf(3)                &
                                      )

  enddo

  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    surfn = surfbn(ifac)
    distbf = distb(ifac)

    ! FIXME: wrong if dirichlet and viscel is really a tensor
    flumab(ifac) = flumab(ifac)                                                   &
                 + viscb(ifac)*distbf/surfn*cofbfp(ifac)*(                        &
                                      frcxt(1, ii)*surfbo(1,ifac)                 &
                                    + frcxt(2, ii)*surfbo(2,ifac)                 &
                                    + frcxt(3, ii)*surfbo(3,ifac) )

  enddo
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format('PROJTV appele avec INIT =',i10)

#else

 1000 format('PROJTV called with INIT =',i10)

#endif

!----
! End
!----

return

end subroutine
