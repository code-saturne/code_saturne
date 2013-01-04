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

 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iphydp , iwarnp , nfecra ,                                     &
   epsrgp , climgp , extrap ,                                     &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp , cofafp , cofbfp , viscf  , viscb  , &
   visel  ,                                                       &
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

integer          nvar   , nscal
integer          init   , inc    , imrgra , iccocg
integer          nswrgp , imligp
integer          iwarnp , iphydp , nfecra
double precision epsrgp , climgp , extrap


double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision visel(3,ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)

! Local variables

integer          ifac, ii, jj, iij, iii
double precision pfac,pip
double precision dpxf  , dpyf  , dpzf
double precision dijpfx, dijpfy, dijpfz
double precision diipbx, diipby, diipbz
double precision dijx  , dijy  , dijz

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
    call synvin(visel)
  endif

  ! Mass flow through interior faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dpxf = 0.5d0*(visel(1,ii)*grad(ii,1) + visel(1,jj)*grad(jj,1))
    dpyf = 0.5d0*(visel(2,ii)*grad(ii,2) + visel(2,jj)*grad(jj,2))
    dpzf = 0.5d0*(visel(3,ii)*grad(ii,3) + visel(3,jj)*grad(jj,3))

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

!---> DIJ = IJ - (IJ.N) N
    dijx = (xyzcen(1,jj)-xyzcen(1,ii))-dijpfx
    dijy = (xyzcen(2,jj)-xyzcen(2,ii))-dijpfy
    dijz = (xyzcen(3,jj)-xyzcen(3,ii))-dijpfz

    flumas(ifac) = flumas(ifac)                                   &
     + viscf(ifac)*( pvar(ii) -pvar(jj) )                         &
     + ( dpxf * dijx                                              &
     +   dpyf * dijy                                              &
     +   dpzf * dijz )*surfan(ifac)/dist(ifac)

  enddo

  ! ---> Contribution from boundary faces

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

    pip = pvar(ii) + grad(ii,1)*diipbx + grad(ii,2)*diipby + grad(ii,3)*diipbz
    pfac = inc*cofafp(ifac) + cofbfp(ifac)*pip

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
