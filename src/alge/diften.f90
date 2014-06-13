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

!> \file diften.f90
!>
!> \brief This function adds the explicit part of the diffusion
!> terms with a symmetric tensor diffusivity for a transport equation of a
!> scalar field \f$ \varia \f$.
!>
!> More precisely, the right hand side \f$ Rhs \f$ is updated as
!> follows:
!> \f[
!> Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
!>      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Warning:
!> - \f$ Rhs \f$ has already been initialized before calling diften!
!> - mind the sign minus
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idtvar        indicator of the temporal scheme
!> \param[in]     ivar          index of the current variable
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
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
!> \param[in]     ipp           index of the variable for post-processing
!> \param[in]     iwarnp        verbosity
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvar          solved variable (current time step)
!> \param[in]     pvara         solved variable (previous time step)
!> \param[in]     coefap        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbp        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafp        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfp        boundary condition array for the diffusion
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
!> \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine diften &
 ( idtvar , ivar   , nswrgp , imligp , ircflp ,          &
   inc    , imrgra , iccocg , ipp    , iwarnp , epsrgp , &
   climgp , extrap , relaxp , thetap ,                   &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp , &
   viscf  , viscb  , viscel ,                            &
   weighf , weighb ,                                     &
   smbrp  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use parall
use period
use cplsat
use optcal, only: iporos
use mesh

!===============================================================================

implicit none

! Arguments

integer          idtvar
integer          ivar   , nswrgp , imligp
integer          ircflp
integer          inc    , imrgra , iccocg
integer          iwarnp , ipp

double precision epsrgp , climgp, extrap, relaxp , thetap

double precision pvar  (ncelet), pvara(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision viscf (nfac)  , viscb (nfabor)
double precision, target :: viscel(6,ncelet)
double precision weighf(2,nfac), weighb(nfabor)
double precision smbrp (ncelet)

! Local variables

integer          ifac,ii,jj,infac,iel, ig, it,i
integer          isou
double precision pfacd,flux,fluxi,fluxj
double precision pi, pj, pia, pja
double precision pir,pjr,pippr,pjppr
double precision diippf(3), djjppf(3), pipp, pjpp
double precision visci(3,3), viscj(3,3)
double precision fikdvi, fjkdvi

double precision, pointer, dimension(:,:) :: viscce
double precision, dimension(:,:), allocatable, target :: w2
double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

viscce => null()

! Allocate work arrays
allocate(grad(ncelet,3))

! Without porosity
if (iporos.eq.0) then
  viscce => viscel(:,:)

! With porosity
else if (iporos.eq.1) then
  allocate(w2(6, ncelet))
  do iel = 1, ncel
    do isou = 1, 6
      w2(isou, iel) = porosi(iel)*viscel(isou, iel)
    enddo
  enddo
  viscce => w2(:,:)

! With tensorial porosity
else if (iporos.eq.2) then
  allocate(w2(6, ncelet))
  do iel = 1, ncel
    call symmetric_matrix_product(w2(1, iel), porosf(1, iel), viscel(1, iel))
  enddo
  viscce => w2(:,:)
endif

! ---> Periodicity and parallelism treatment of symmetric tensors

if (irangp.ge.0.or.iperio.eq.1) then
  call syntis(viscce)
endif

!===============================================================================
! 2. Compute the diffusive part with reconstruction technics
!===============================================================================

! ======================================================================
! ---> Compute the gradient of the current variable if needed
! ======================================================================

if (ircflp.eq.1) then

  call grdcel &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   pvar   , coefap , coefbp ,                                     &
   grad   )

else
  !$omp parallel do
  do iel = 1, ncelet
    grad(iel,1) = 0.d0
    grad(iel,2) = 0.d0
    grad(iel,3) = 0.d0
  enddo
endif

! ======================================================================
! ---> Contribution from interior faces
! ======================================================================

infac = 0

if (ncelet.gt.ncel) then
  !$omp parallel do if(ncelet - ncel > thr_n_min)
  do iel = ncel+1, ncelet
    smbrp(iel) = 0.d0
  enddo
endif

! Steady
if (idtvar.lt.0) then

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, i, visci, viscj,    &
    !$omp                     pipp, pjpp, pippr, pjppr,         &
    !$omp                     fluxi, fluxj, fikdvi,             &
    !$omp                     pi, pj, pir, pjr, pia, pja)       &
    !$omp             reduction(+:infac)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        ! in parallel, face will be counted by one and only one rank
        if (ii.le.ncel) then
          infac = infac+1
        endif

        pi = pvar(ii)
        pj = pvar(jj)
        pia = pvara(ii)
        pja = pvara(jj)

        ! Recompute II" and JJ"
        !----------------------

        visci(1,1) = viscce(1,ii)
        visci(2,2) = viscce(2,ii)
        visci(3,3) = viscce(3,ii)
        visci(1,2) = viscce(4,ii)
        visci(2,1) = viscce(4,ii)
        visci(2,3) = viscce(5,ii)
        visci(3,2) = viscce(5,ii)
        visci(1,3) = viscce(6,ii)
        visci(3,1) = viscce(6,ii)

        ! IF.Ki.S / ||Ki.S||^2
        fikdvi = weighf(1,ifac)

        ! II" = IF + FI"
        do i = 1, 3
          diippf(i) = cdgfac(i,ifac)-xyzcen(i,ii)          &
                    - fikdvi*( visci(i,1)*surfac(1,ifac)   &
                             + visci(i,2)*surfac(2,ifac)   &
                             + visci(i,3)*surfac(3,ifac) )
        enddo

        viscj(1,1) = viscce(1,jj)
        viscj(2,2) = viscce(2,jj)
        viscj(3,3) = viscce(3,jj)
        viscj(1,2) = viscce(4,jj)
        viscj(2,1) = viscce(4,jj)
        viscj(2,3) = viscce(5,jj)
        viscj(3,2) = viscce(5,jj)
        viscj(1,3) = viscce(6,jj)
        viscj(3,1) = viscce(6,jj)

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

        pir = pi/relaxp - (1.d0-relaxp)/relaxp * pia
        pjr = pj/relaxp - (1.d0-relaxp)/relaxp * pja

        ! pr in I" and J"
        pippr = pir + ircflp*( grad(ii,1)*diippf(1)   &
                             + grad(ii,2)*diippf(2)   &
                             + grad(ii,3)*diippf(3))
        pjppr = pjr + ircflp*( grad(jj,1)*djjppf(1)   &
                             + grad(jj,2)*djjppf(2)   &
                             + grad(jj,3)*djjppf(3))


        fluxi = viscf(ifac)*(pippr - pjpp)
        fluxj = viscf(ifac)*(pipp - pjppr)

        smbrp(ii) = smbrp(ii) - fluxi
        smbrp(jj) = smbrp(jj) + fluxj

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, i, visci, viscj, fikdvi, fjkdvi,  &
    !$omp                     pipp, pjpp, diippf, djjppf,                     &
    !$omp                     flux, pi, pj)                                   &
    !$omp             reduction(+:infac)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        ! in parallel, face will be counted by one and only one rank
        if (ii.le.ncel) then
          infac = infac+1
        endif

        pi = pvar(ii)
        pj = pvar(jj)

        ! Recompute II" and JJ"
        !----------------------

        visci(1,1) = viscce(1,ii)
        visci(2,2) = viscce(2,ii)
        visci(3,3) = viscce(3,ii)
        visci(1,2) = viscce(4,ii)
        visci(2,1) = viscce(4,ii)
        visci(2,3) = viscce(5,ii)
        visci(3,2) = viscce(5,ii)
        visci(1,3) = viscce(6,ii)
        visci(3,1) = viscce(6,ii)

        ! IF.Ki.S / ||Ki.S||^2
        fikdvi = weighf(1,ifac)

        ! II" = IF + FI"
        do i = 1, 3
          diippf(i) = cdgfac(i,ifac)-xyzcen(i,ii)          &
                    - fikdvi*( visci(i,1)*surfac(1,ifac)   &
                             + visci(i,2)*surfac(2,ifac)   &
                             + visci(i,3)*surfac(3,ifac) )
        enddo

        viscj(1,1) = viscce(1,jj)
        viscj(2,2) = viscce(2,jj)
        viscj(3,3) = viscce(3,jj)
        viscj(1,2) = viscce(4,jj)
        viscj(2,1) = viscce(4,jj)
        viscj(2,3) = viscce(5,jj)
        viscj(3,2) = viscce(5,jj)
        viscj(1,3) = viscce(6,jj)
        viscj(3,1) = viscce(6,jj)

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

        flux = viscf(ifac)*(pipp -pjpp)

        smbrp(ii) = smbrp(ii) - thetap*flux
        smbrp(jj) = smbrp(jj) + thetap*flux

      enddo
    enddo
  enddo
endif

! ======================================================================
! ---> Contribution from boundary faces
! ======================================================================

! Steady
if (idtvar.lt.0) then

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, i, visci, fikdvi,                        &
    !$omp                     pir, pippr, pfacd, flux, pi, pia)                  &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)
        pia = pvara(ii)

        pir = pi/relaxp - (1.d0-relaxp)/relaxp*pia

        ! Recompute II"
        !--------------

        visci(1,1) = viscce(1,ii)
        visci(2,2) = viscce(2,ii)
        visci(3,3) = viscce(3,ii)
        visci(1,2) = viscce(4,ii)
        visci(2,1) = viscce(4,ii)
        visci(2,3) = viscce(5,ii)
        visci(3,2) = viscce(5,ii)
        visci(1,3) = viscce(6,ii)
        visci(3,1) = viscce(6,ii)

        ! IF.Ki.S / ||Ki.S||^2
        fikdvi = weighb(ifac)

        ! II" = IF + FI"
        do i = 1, 3
          diippf(i) = cdgfbo(i,ifac) - xyzcen(i,ii)        &
                    - fikdvi*( visci(i,1)*surfbo(1,ifac)   &
                             + visci(i,2)*surfbo(2,ifac)   &
                             + visci(i,3)*surfbo(3,ifac) )
        enddo

        pippr = pir                             &
              + ircflp*( grad(ii,1)*diippf(1)   &
                       + grad(ii,2)*diippf(2)   &
                       + grad(ii,3)*diippf(3))

        pfacd = inc*cofafp(ifac) +cofbfp(ifac)*pippr

        flux = viscb(ifac)*pfacd
        smbrp(ii) = smbrp(ii) - flux

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, visci, fikdvi,                          &
    !$omp                     diippf, pipp, pfacd, flux, pi)                    &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)

        ! Recompute II"
        !--------------

        visci(1,1) = viscce(1,ii)
        visci(2,2) = viscce(2,ii)
        visci(3,3) = viscce(3,ii)
        visci(1,2) = viscce(4,ii)
        visci(2,1) = viscce(4,ii)
        visci(2,3) = viscce(5,ii)
        visci(3,2) = viscce(5,ii)
        visci(1,3) = viscce(6,ii)
        visci(3,1) = viscce(6,ii)

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


        pfacd = inc*cofafp(ifac) + cofbfp(ifac)*pipp

        flux = viscb(ifac)*pfacd
        smbrp(ii) = smbrp(ii) - thetap * flux

      enddo
    enddo
  enddo

endif

! Free memory
deallocate(grad)
if (allocated(w2)) deallocate(w2)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
