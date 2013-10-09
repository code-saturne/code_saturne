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

!> \file diftnv.f90
!>
!> \brief This function adds the explicit part of the diffusion
!> terms with a symmetric tensor diffusivity for a transport equation of a
!> vector field \f$ \vect{\varia} \f$.
!>
!> More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
!> follows:
!> \f[
!> \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
!>      - \tens{\mu}_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Warning:
!> - \f$ \vect{Rhs} \f$ has already been initialized before calling diftnv!
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
!> \param[in]     ippu          index of the variable for post-processing
!> \param[in]     iwarnp        verbosity
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvar          solved variable (current time step)
!> \param[in]     pvara         solved variable (previous time step)
!> \param[in]     coefav        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     viscf         \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine diftnv &
 ( idtvar , ivar   , nswrgp , imligp , ircflp ,          &
   inc    , imrgra ,                                     &
   ippu   , iwarnp ,                                     &
   epsrgp ,                                              &
   climgp , relaxp , thetap ,                            &
   pvar   , pvara  ,                                     &
   coefav , coefbv , cofafv , cofbfv ,                   &
   viscf  , viscb  ,                                     &
   rhs    )

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
use mesh

!===============================================================================

implicit none

! Arguments

integer          idtvar
integer          ivar   , nswrgp , imligp
integer          ircflp
integer          inc    , imrgra
integer          iwarnp , ippu

double precision epsrgp , climgp, relaxp , thetap

double precision pvar  (3  ,ncelet)
double precision pvara (3  ,ncelet)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision cofbfv(3,3,nfabor)
double precision viscf (3,3,nfac)  , viscb (nfabor)
double precision rhs(3,ncelet)

! Local variables

integer          ifac,ii,jj,infac,iel, ig, it,isou, jsou
logical          ilved
double precision pfacd,flux,fluxi,fluxj, pnd
double precision pi, pj, pia, pja
double precision pir,pipr(3),pjpr(3)
double precision dpvf(3)
double precision dijpfv(3)
double precision diipfv(3), djjpfv(3), pip(3), pjp(3)
double precision diipbv(3)

double precision, allocatable, dimension(:,:,:) :: gradv

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate work arrays
allocate(gradv(3,3,ncelet))

! Initialize variables to avoid compiler warnings

pi  = 0.d0
pj  = 0.d0
pia = 0.d0
pja = 0.d0

!===============================================================================
! 2. Compute the diffusive part with reconstruction technics
!===============================================================================

! ======================================================================
! ---> Compute the gradient of the current variable if needed
! ======================================================================

if (ircflp.eq.1) then

  ilved = .true.

  call grdvec &
  !==========
( ivar   , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , epsrgp , climgp ,                                     &
  ilved  ,                                                       &
  pvar   , coefav , coefbv ,                                     &
  gradv )

else
  !$omp parallel do private(isou, jsou)
  do iel = 1, ncelet
    do isou =1, 3
      do jsou = 1, 3
        gradv(jsou,isou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

! ======================================================================
! ---> Contribution from interior faces
! ======================================================================

infac = 0

if (ncelet.gt.ncel) then
  !$omp parallel do private(isou) if(ncelet -ncel > thr_n_min)
  do iel = ncel+1, ncelet
    do isou = 1, 3
      rhs(isou,iel) = 0.d0
    enddo
  enddo
endif

! Steady
if (idtvar.lt.0) then

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, isou, jsou, pnd, dijpfv,     &
    !$omp                     diipfv, djjpfv, dpvf, pi, pj,              &
    !$omp                     pia, pja, pip, pjp, pipr, pjpr,            &
    !$omp                     fluxi, fluxj)                  &
    !$omp             reduction(+:infac)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        ! in parallel, face will be counted by one and only one rank
        if (ii.le.ncel) then
          infac = infac+1
        endif

        do jsou = 1, 3
          dijpfv(jsou) = dijpf(jsou,ifac)
        enddo

        pnd = pond(ifac)

        ! Recompute II' and JJ' at this level
        do jsou = 1, 3
          diipfv(jsou) =   cdgfac(jsou,ifac) - (xyzcen(jsou,ii)           &
                         + (1.d0-pnd) * dijpfv(jsou))
          djjpfv(jsou) =   cdgfac(jsou,ifac) -  xyzcen(jsou,jj)           &
                         + pnd  * dijpfv(jsou)
        enddo

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          do jsou = 1, 3
            dpvf(jsou) = 0.5d0*(gradv(jsou,isou,ii) + gradv(jsou,isou,jj))
          enddo

          pi  = pvar (isou,ii)
          pj  = pvar (isou,jj)

          pia = pvara(isou,ii)
          pja = pvara(isou,jj)

          ! reconstruction only if IRCFLP = 1
          pip(isou) = pi + ircflp*(dpvf(1)*diipfv(1)        &
                                  +dpvf(2)*diipfv(2)        &
                                  +dpvf(3)*diipfv(3))
          pjp(isou) = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                                  +dpvf(2)*djjpfv(2)        &
                                  +dpvf(3)*djjpfv(3))

          pipr(isou) = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
                     + ircflp*(dpvf(1)*diipfv(1)                 &
                              +dpvf(2)*diipfv(2)                 &
                              +dpvf(3)*diipfv(3))
          pjpr(isou) = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
                     + ircflp*(dpvf(1)*djjpfv(1)                 &
                              +dpvf(2)*djjpfv(2)                 &
                              +dpvf(3)*djjpfv(3))

        enddo

        do isou = 1, 3

          fluxi = viscf(isou,1,ifac)*(pipr(1) - pjp(1))  &
                + viscf(isou,2,ifac)*(pipr(2) - pjp(2))  &
                + viscf(isou,3,ifac)*(pipr(3) - pjp(3))
          fluxj = viscf(isou,1,ifac)*(pip(1) - pjpr(1))  &
                + viscf(isou,2,ifac)*(pip(2) - pjpr(2))  &
                + viscf(isou,3,ifac)*(pip(3) - pjpr(3))

          rhs(isou,ii) = rhs(isou,ii) - fluxi
          rhs(isou,jj) = rhs(isou,jj) + fluxj

        enddo

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, isou, jsou, pnd, dijpfv,     &
    !$omp                     diipfv, djjpfv, dpvf, pi, pj,              &
    !$omp                     pip, pjp, flux)                            &
    !$omp             reduction(+:infac)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        ! in parallel, face will be counted by one and only one rank
        if (ii.le.ncel) then
          infac = infac+1
        endif

        do jsou = 1, 3
          dijpfv(jsou) = dijpf(jsou,ifac)
        enddo

        pnd = pond(ifac)

        ! Recompute II' and JJ' at this level
        do jsou = 1, 3
          diipfv(jsou) =   cdgfac(jsou,ifac) - (xyzcen(jsou,ii)             &
                         + (1.d0-pnd) * dijpfv(jsou))
          djjpfv(jsou) =   cdgfac(jsou,ifac) -  xyzcen(jsou,jj)             &
                         +  pnd * dijpfv(jsou)
        enddo

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          do jsou = 1, 3
            dpvf(jsou) = 0.5d0*(gradv(jsou,isou,ii) + gradv(jsou,isou,jj))
          enddo

          pi = pvar(isou,ii)
          pj = pvar(isou,jj)

          pip(isou) = pi + ircflp*(dpvf(1)*diipfv(1)        &
                                  +dpvf(2)*diipfv(2)        &
                                  +dpvf(3)*diipfv(3))
          pjp(isou) = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                                  +dpvf(2)*djjpfv(2)        &
                                  +dpvf(3)*djjpfv(3))

        enddo

        do isou = 1, 3

          flux = viscf(isou,1,ifac)*(pip(1) - pjp(1))  &
               + viscf(isou,2,ifac)*(pip(2) - pjp(2))  &
               + viscf(isou,3,ifac)*(pip(3) - pjp(3))

          rhs(isou,ii) = rhs(isou,ii) - thetap*flux
          rhs(isou,jj) = rhs(isou,jj) + thetap*flux

        enddo

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
    !$omp parallel do private(ifac, ii, isou, jsou, diipbv,               &
    !$omp                     pfacd, pir, pipr, flux)                     &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        do jsou = 1, 3
          diipbv(jsou) = diipb(jsou,ifac)
        enddo

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          pfacd = inc*cofafv(isou,ifac)

          !coefu and cofuf are matrices
          do jsou = 1, 3
            pir  = pvar(jsou,ii)/relaxp - (1.d0-relaxp)/relaxp*pvara(jsou,ii)

            pipr(jsou) = pir +ircflp*( gradv(1,jsou,ii)*diipbv(1)         &
                                     + gradv(2,jsou,ii)*diipbv(2)         &
                                     + gradv(3,jsou,ii)*diipbv(3))
            pfacd = pfacd + cofbfv(isou,jsou,ifac)*pipr(jsou)
          enddo

          flux = viscb(ifac)*pfacd
          rhs(isou,ii) = rhs(isou,ii) - flux

        enddo ! isou

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, isou, jsou, diipbv,    &
    !$omp                     pfacd, pir, flux, pi)            &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        do jsou = 1, 3
          diipbv(jsou) = diipb(jsou,ifac)
        enddo

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          pfacd = inc*cofafv(isou,ifac)

          !coefu and cofuf are matrices
          do jsou = 1, 3
            pir = pvar(jsou,ii) + ircflp*( gradv(1,jsou,ii)*diipbv(1)        &
                                         + gradv(2,jsou,ii)*diipbv(2)        &
                                         + gradv(3,jsou,ii)*diipbv(3))
            pfacd = pfacd + cofbfv(isou,jsou,ifac)*pir
          enddo

          flux = viscb(ifac)*pfacd
          rhs(isou,ii) = rhs(isou,ii) - thetap * flux

        enddo ! isou

      enddo
    enddo
  enddo

endif ! idtvar

! Free memory
deallocate(gradv)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

#else

#endif

!----
! End
!----

return

end subroutine
