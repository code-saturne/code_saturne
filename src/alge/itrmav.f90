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

!> \file itrmav.f90
!>
!> \brief This function adds the explicit part of the pressure gradient
!> to the mass flux for a tensorial diffusion for the pressure
!> field \f$ P \f$.
!>
!> More precisely, the mass flux side \f$ \dot{m}_\fij \f$ is updated as
!> follows:
!> \f[
!> \dot{m}_\fij = \dot{m}_\fij -
!>              \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
!> \f]
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
!> \param[in]     iccocg        indicator
!>                               - 1 re-compute cocg matrix (for iterativ gradients)
!>                               - 0 otherwise
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     iphydr        indicator
!>                               - 1 hydrostatic pressure taken into account
!>                               - 0 otherwise
!> \param[in]     iwarnp        verbosity
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
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

subroutine itrmav &
!================

 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydp , iwarnp , nfecra ,                                     &
   epsrgp , climgp , extrap ,                                     &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp , cofafp , cofbfp , viscf  , viscb  , &
   viscel ,                                                       &
   weighf , weighb ,                                              &
   flumas , flumab )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          init   , inc    , imrgra , iccocg
integer          nswrgp , imligp
integer          ircflp
integer          iwarnp , iphydp , nfecra
double precision epsrgp , climgp , extrap


double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision viscel(6,ncelet)
double precision weighf(2,nfac), weighb(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)

! Local variables

integer          ifac, ii, jj, i
double precision pfac
double precision pi, pj
double precision diippf(3), djjppf(3), pipp, pjpp
double precision visci(3,3), viscj(3,3)
double precision fikdvi, fjkdvi

double precision rvoid(1)

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

if (init.ge.1) then
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo
elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

! Handle parallelism and periodicity

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(pvar)
endif

!===============================================================================
! 2. Update mass flux without reconstruction technics
!===============================================================================

if (nswrgp.le.1) then

  ! ---> Contribution from interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas(ifac) = flumas(ifac) + viscf(ifac)*(pvar(ii) - pvar(jj))

  enddo

  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    pfac = inc*cofafp(ifac) + cofbfp(ifac)*pvar(ii)

    flumab(ifac) = flumab(ifac) + viscb(ifac)*pfac

  enddo

endif

!===============================================================================
! 3. Update mass flux WITH reconstruction technics
!===============================================================================

if (nswrgp.gt.1) then

  ! Allocate a work array for the gradient calculation
  allocate(grad(ncelet,3))

  call grdpot &
  !==========
 ( ipr    , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp ,                                     &
   grad   )

  ! ---> Periodicity and parallelism treatment of symmetric tensors

  if (irangp.ge.0.or.iperio.eq.1) then
    call syntis(viscel)
  endif

  ! Mass flow through interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    pi = pvar(ii)
    pj = pvar(jj)

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

    ! p in I" and J"
    pipp = pi + ircflp*( grad(ii,1)*diippf(1)   &
                       + grad(ii,2)*diippf(2)   &
                       + grad(ii,3)*diippf(3))
    pjpp = pj + ircflp*( grad(jj,1)*djjppf(1)   &
                       + grad(jj,2)*djjppf(2)   &
                       + grad(jj,3)*djjppf(3))

    flumas(ifac) = flumas(ifac) + viscf(ifac)*(pipp - pjpp)

  enddo

  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    pi = pvar(ii)

    ! Recompute II"
    !--------------

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
    fikdvi = weighb(ifac)

    ! II" = IF + FI"
    do i = 1, 3
      diippf(i) = cdgfbo(i,ifac) - xyzcen(i,ii)        &
                - fikdvi*( visci(i,1)*surfbo(1,ifac)   &
                         + visci(i,2)*surfbo(2,ifac)   &
                         + visci(i,3)*surfbo(3,ifac) )
    enddo

    pipp = pi                               &
         + ircflp*( grad(ii,1)*diippf(1)    &
                  + grad(ii,2)*diippf(2)    &
                  + grad(ii,3)*diippf(3))


    pfac = inc*cofafp(ifac) + cofbfp(ifac)*pipp

    flumab(ifac) = flumab(ifac) + viscb(ifac)*pfac

  enddo

  ! Free memory
  deallocate(grad)

endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format('ITRMAS appele avec INIT = ',i10)

#else

 1000 format('ITRMAS called with INIT = ',i10)

#endif

!----
! End
!----

return

end subroutine
