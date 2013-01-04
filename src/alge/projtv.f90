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
!> \param[in]     inc           indicator
!>                               - 0 when solving an increment
!>                               - 1 otherwise
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     iwarnp        verbosity
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     fextx,y,z     body force creating the hydrostatic pressure
!> \param[in]     pvar          solved variable (pressure)
!> \param[in]     cofaf         boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbf         boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
!> \param[in]     weighf        internal face weight between cells i j in case
!>                               of tensor diffusion
!> \param[in]     weighb        boundary face weight for cells i in case
!>                               of tensor diffusion
!> \param[in,out] flumas        mass flux at interior faces
!> \param[in,out] flumab        mass flux at boundary faces
!_______________________________________________________________________________

subroutine projtv &
!================

 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , nswrgu , imligu ,                   &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu ,                                              &
   fextx  , fexty  , fextz  ,                                     &
   cofbfp ,                                                       &
   flumas , flumab , viscf  , viscb  ,                            &
   visel  )

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

integer          nvar   , nscal
integer          init   , inc    , imrgra
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu


double precision pnd
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision visel(3,ncelet)
double precision cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, iii
double precision dijpfx,dijpfy,dijpfz
double precision diipx,diipy,diipz
double precision djjpx,djjpy,djjpz
double precision distbf,surfn

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

if( nswrgu.le.1 ) then

  ! ---> Contribution from interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas(ifac) =  flumas(ifac)                                  &
         + viscf(ifac)*(                                          &
           (cdgfac(1,ifac)-xyzcen(1,ii))*fextx(ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*fexty(ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*fextz(ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*fextx(jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*fexty(jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*fextz(jj) )

  enddo


  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn          &
         *cofbfp(ifac)*(fextx(ii)*surfbo(1,ifac)                  &
         +fexty(ii)*surfbo(2,ifac)+fextz(ii)*surfbo(3,ifac) )

  enddo


else
!===============================================================================
! 3. Update mass flux WITH reconstruction technics
!===============================================================================

  ! ---> Contribution from interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pnd = pond(ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    surfn = surfan(ifac)

!     calcul de II' et JJ'
    diipx = cdgfac(1,ifac)-xyzcen(1,ii)-(1.d0-pnd)*dijpfx
    diipy = cdgfac(2,ifac)-xyzcen(2,ii)-(1.d0-pnd)*dijpfy
    diipz = cdgfac(3,ifac)-xyzcen(3,ii)-(1.d0-pnd)*dijpfz
    djjpx = cdgfac(1,ifac)-xyzcen(1,jj)+pnd*dijpfx
    djjpy = cdgfac(2,ifac)-xyzcen(2,jj)+pnd*dijpfy
    djjpz = cdgfac(3,ifac)-xyzcen(3,jj)+pnd*dijpfz

    flumas(ifac) =  flumas(ifac)                                  &
         + viscf(ifac)*(                                          &
           (cdgfac(1,ifac)-xyzcen(1,ii))*fextx(ii)                &
          +(cdgfac(2,ifac)-xyzcen(2,ii))*fexty(ii)                &
          +(cdgfac(3,ifac)-xyzcen(3,ii))*fextz(ii)                &
          -(cdgfac(1,ifac)-xyzcen(1,jj))*fextx(jj)                &
          -(cdgfac(2,ifac)-xyzcen(2,jj))*fexty(jj)                &
          -(cdgfac(3,ifac)-xyzcen(3,jj))*fextz(jj) )              &
         +surfn/dist(ifac)*0.5d0*(                                &
       (djjpx-diipx)*(visel(1,ii)*fextx(ii)+visel(1,jj)*fextx(jj))  &
      +(djjpy-diipy)*(visel(2,ii)*fexty(ii)+visel(2,jj)*fexty(jj))  &
      +(djjpz-diipz)*(visel(3,ii)*fextz(ii)+visel(3,jj)*fextz(jj)))

  enddo


  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    surfn = surfbn(ifac)
    distbf = distb(ifac)

    flumab(ifac) = flumab(ifac)+viscb(ifac)*distbf/surfn          &
         *cofbfp(ifac)*(fextx(ii)*surfbo(1,ifac)                  &
         +fexty(ii)*surfbo(2,ifac)+fextz(ii)*surfbo(3,ifac) )

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
