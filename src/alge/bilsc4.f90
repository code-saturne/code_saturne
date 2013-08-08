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

!> \file bilsc4.f90
!>
!> \brief This function adds the explicit part of the convection/diffusion
!> terms of a transport equation of a vector field \f$ \vect{\varia} \f$.
!>
!> More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
!> follows:
!> \f[
!> \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
!>        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
!>      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Remark:
!> if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
!> + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
!> the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
!>
!> Warning:
!> - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
!> - mind the sign minus
!>
!> Options:
!> - blencp = 0: upwind scheme for the advection
!> - blencp = 1: no upwind scheme except in the slope test
!> - ischcp = 0: second order
!> - ischcp = 1: centred
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idtvar        indicator of the temporal scheme
!> \param[in]     ivar          index of the current variable
!> \param[in]     iconvp        indicator
!>                               - 1 convection,
!>                               - 0 sinon
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 sinon
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     ischcp        indicator
!>                               - 1 centred
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
!> \param[in]     inc           indicator
!>                               - 0 when solving an increment
!>                               - 1 otherwise
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     ivisep        indicator to take \f$ \divv
!>                               \left(\mu \gradt \transpose{\vect{a}} \right)
!>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
!>                               - 1 take into account,
!>                               - 0 otherwise
!> \param[in]     ippu          index of the variable for post-processing
!> \param[in]     ippv          index of the variable for post-processing
!> \param[in]     ippw          index of the variable for post-processing
!> \param[in]     iwarnp        verbosity
!> \param[in]     blencp        fraction of upwinding
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
!> \param[in]     pvar          vitesse resolue (instant courant)
!> \param[in]     pvara         vitesse resolue (instant precedent)
!> \param[in]     coefav        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     secvif        secondary viscosity at interior faces
!> \param[in]     secvib        secondary viscosity at boundary faces
!> \param[in,out] smbr          right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine bilsc4 &
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   ippu   , ippv   , ippw   , iwarnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscf  , viscb  , secvif , secvib ,          &
   smbr   )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
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
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , ivisep
integer          iwarnp , ippu   , ippv   , ippw

double precision blencp , epsrgp , climgp, extrap, relaxp , thetap
double precision pvar  (3  ,ncelet)
double precision pvara (3  ,ncelet)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision cofbfv(3,3,nfabor)
double precision flumas(nfac)  , flumab(nfabor)
double precision viscf (nfac)  , viscb (nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision smbr(3,ncelet)


! Local variables

character*80     chaine
character*8      cnom
integer          ifac,ii,jj,infac,iel,iupwin, ig, it
integer          iiu,iiv,iiw
integer          iitytu
integer          iir11,iir22,iir33
integer          iir12,iir13,iir23
integer          isou, jsou, ityp
logical          ilved
double precision pfac,pfacd,flui,fluj,flux,fluxi,fluxj
double precision vfac(3)
double precision difv(3), djfv(3)
double precision pi , pj, pia, pja
double precision pif,pjf,pip,pjp,pir,pjr,pipr,pjpr
double precision pifr,pjfr,pifri,pifrj,pjfri,pjfrj
double precision testi,testj,testij
double precision dpvf(3)
double precision dcc, ddi, ddj, tesqck
double precision dijpfv(3)
double precision diipfv(3)
double precision djjpfv(3)
double precision diipbv(3)
double precision pnd, distf, srfan
double precision unsvol, visco, grdtrv, tgrdfl, secvis

double precision, dimension(:,:,:), allocatable :: gradv, gradva
double precision, dimension(:), allocatable :: bndcel


!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(gradv(3,3,ncelet))
allocate(gradva(3,3,ncelet))

! Initialize variables to avoid compiler warnings

pif = 0.d0
pjf = 0.d0
pifri = 0.d0
pifrj = 0.d0
pjfri = 0.d0
pjfrj = 0.d0

pi  = 0.d0
pj  = 0.d0
pia = 0.d0
pja = 0.d0

! Memoire

chaine = nomvar(ippu)
cnom   = chaine(1:8)

if (iwarnp.ge.2) then
  if (ischcp.eq.1) then
    write(nfecra,1000) cnom, '  CENTERED ', (1.d0-blencp)*100.d0
  else
    write(nfecra,1000) cnom, ' 2ND ORDER ', (1.d0-blencp)*100.d0
  endif
endif

iupwin = 0
if (blencp.eq.0.d0) iupwin = 1

!===============================================================================
! 2.  CALCUL DU BILAN AVEC TECHNIQUE DE RECONSTRUCTION
!===============================================================================

! ======================================================================
! ---> CALCUL DU GRADIENT DE VITESSE
! ======================================================================
!    DUDX sert a la fois pour la reconstruction des flux et pour le test
!    de pente. On doit donc le calculer :
!        - quand on a de la diffusion et qu'on reconstruit les flux
!        - quand on a de la convection SOLU
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on reconstruit les flux
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on n'a pas shunte le test de pente

if ((idiffp.ne.0 .and. ircflp.eq.1) .or. ivisep.eq.1 .or.         &
    (iconvp.ne.0 .and. iupwin.eq.0 .and.                          &
    (ischcp.eq.0 .or.  ircflp.eq.1 .or. isstpp.eq.0))) then


  ilved = .true.

  call grdvec &
  !==========
( ivar   , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
  ilved  ,                                                       &
  pvar   , coefav , coefbv ,                                     &
  gradv )

else
  !$omp parallel do private(isou, jsou)
  do iel = 1, ncelet
    do isou =1, 3
      do jsou = 1, 3
        gradv(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

! ======================================================================
! ---> Compute uncentred gradient gradva for the slope test
! ======================================================================

!$omp parallel do private(isou, jsou)
do iel = 1, ncelet
  do jsou = 1, 3
    do isou =1, 3
      gradva(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

if (iconvp.gt.0.and.iupwin.eq.0.and.isstpp.eq.0) then

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, isou, jsou, difv, djfv, pif, pjf, &
    !$omp                     pfac, vfac)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        do jsou = 1, 3
          difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
          djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
        enddo
        !-----------------
        ! X-Y-Z component, p=u, v, w
        do isou = 1, 3
          pif = pvar(isou,ii)
          pjf = pvar(isou,jj)
          do jsou = 1, 3
            pif = pif + gradv(isou,jsou,ii)*difv(jsou)
            pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
          enddo

          pfac = pjf
          if (flumas(ifac).gt.0.d0) pfac = pif

          ! U gradient
          do jsou = 1, 3
            vfac(jsou) = pfac*surfac(jsou,ifac)

            gradva(isou,jsou,ii) = gradva(isou,jsou,ii) + vfac(jsou)
            gradva(isou,jsou,jj) = gradva(isou,jsou,jj) - vfac(jsou)
          enddo
        enddo

      enddo
    enddo
  enddo

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, isou, jsou, diipbv, pfac)   &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        do jsou = 1, 3
          diipbv(jsou) = diipb(jsou,ifac)
        enddo
        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1,3
          pfac = inc*coefav(isou,ifac)
          !coefu is a matrix
          do jsou =  1, 3
            pfac = pfac + coefbv(isou,jsou,ifac)*(pvar(jsou,ii)    &
                        + gradv(jsou,1,ii)*diipbv(1)               &
                        + gradv(jsou,2,ii)*diipbv(2)               &
                        + gradv(jsou,3,ii)*diipbv(3))
          enddo

          do jsou = 1, 3
            gradva(isou,jsou,ii) = gradva(isou,jsou,ii) +pfac*surfbo(jsou,ifac)
          enddo
        enddo

      enddo
    enddo
  enddo

  !$omp parallel do private(isou, jsou, unsvol)
  do iel = 1, ncel
    unsvol = 1.d0/volume(iel)
    do isou = 1, 3
      do jsou = 1, 3
        gradva(isou,jsou,iel) = gradva(isou,jsou,iel)*unsvol
      enddo
    enddo
  enddo

  ! Handle parallelism and periodicity

  if (irangp.ge.0.or.iperio.eq.1) then
    call syntin (gradva)
    !==========
  endif

endif

! ======================================================================
! ---> Contribution from interior faces
! ======================================================================

infac = 0

if (ncelet.gt.ncel) then
  !$omp parallel do private(isou) if(ncelet -ncel > thr_n_min)
  do iel = ncel+1, ncelet
    do isou = 1, 3
      smbr(isou,iel) = 0.d0
    enddo
  enddo
endif

! --> Pure upwind flux
! =====================

if (iupwin.eq.1) then

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, pnd, dijpfv,     &
      !$omp                     diipfv, djjpfv, flui, fluj, dpvf, pi, pj,  &
      !$omp                     pia, pja, pip, pjp, pipr, pjpr,            &
      !$omp                     pifr, pjfr, fluxi, fluxj)                  &
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

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            ! reconstruction only if IRCFLP = 1
            pi  = pvar (isou,ii)
            pj  = pvar (isou,jj)

            pia = pvara(isou,ii)
            pja = pvara(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
                 + ircflp*(dpvf(1)*diipfv(1)                 &
                          +dpvf(2)*diipfv(2)                 &
                          +dpvf(3)*diipfv(3))
            pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
                 + ircflp*(dpvf(1)*djjpfv(1)                 &
                          +dpvf(2)*djjpfv(2)                 &
                          +dpvf(3)*djjpfv(3))

            pifr = pi /relaxp - (1.d0-relaxp)/relaxp * pia
            pjfr = pj /relaxp - (1.d0-relaxp)/relaxp * pja

            fluxi = iconvp*(flui*pifr + fluj*pj - flumas(ifac)*pi)    &
                  + idiffp*viscf(ifac)*(pipr -pjp)
            fluxj = iconvp*(flui*pi + fluj*pjfr - flumas(ifac)*pj)    &
                  + idiffp*viscf(ifac)*(pip -pjpr)

            smbr(isou,ii) = smbr(isou,ii) - fluxi
            smbr(isou,jj) = smbr(isou,jj) + fluxj

          enddo

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, pnd, dijpfv,     &
      !$omp                     diipfv, djjpfv, flui, fluj, dpvf, pi, pj,  &
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

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            pi = pvar(isou,ii)
            pj = pvar(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            flux =   iconvp*(flui*pi +fluj*pj)          &
                   + idiffp*viscf(ifac)*(pip -pjp)

            smbr(isou,ii) = smbr(isou,ii) - thetap*(flux - iconvp*flumas(ifac)*pi)
            smbr(isou,jj) = smbr(isou,jj) + thetap*(flux - iconvp*flumas(ifac)*pj)

          enddo

        enddo
      enddo
    enddo

  endif


! --> Flux with no slope test
! ============================

elseif (isstpp.eq.1) then

  if (ischcp.lt.0 .or. ischcp.gt.1) then
    write(nfecra,9000) ischcp
    call csexit(1)
  endif

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, dijpfv, pnd,          &
      !$omp                     diipfv, djjpfv, flui, fluj, difv, djfv, dpvf,   &
      !$omp                     pi, pj, pia, pja, pip, pjp, pipr, pjpr,         &
      !$omp                     pir, pjr, pifri, pjfri, pifrj, pjfrj,           &
      !$omp                     fluxi, fluxj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          do jsou = 1, 3
            dijpfv(jsou) = dijpf(jsou,ifac)
          enddo

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level
          do jsou = 1, 3
            diipfv(jsou) =   cdgfac(jsou,ifac) - (xyzcen(jsou,ii)               &
                           + (1.d0-pnd) * dijpfv(jsou))
            djjpfv(jsou) =   cdgfac(jsou,ifac) -  xyzcen(jsou,jj)               &
                           + pnd  * dijpfv(jsou)
          enddo

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          ! For second order, define IF
          if (ischcp.eq.0) then
            do jsou = 1, 3
              difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
              djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
            enddo
          endif

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            pi = pvar (isou,ii)
            pj = pvar (isou,jj)

            pia = pvara(isou,ii)
            pja = pvara(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
                 + ircflp*(dpvf(1)*diipfv(1)                 &
                          +dpvf(2)*diipfv(2)                 &
                          +dpvf(3)*diipfv(3))
            pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
                 + ircflp*(dpvf(1)*djjpfv(1)                 &
                          +dpvf(2)*djjpfv(2)                 &
                          +dpvf(3)*djjpfv(3))

            pir = pi /relaxp - (1.d0 - relaxp)/relaxp* pia
            pjr = pj /relaxp - (1.d0 - relaxp)/relaxp* pja

            ! Centered
            ! --------

            if (ischcp.eq.1) then

              pifri = pnd*pipr + (1.d0-pnd)*pjp
              pjfri = pifri
              pifrj = pnd*pip  + (1.d0-pnd)*pjpr
              pjfrj = pifrj

            ! Second order
            ! ------------

            else ! if (ischcp.eq.0) then
              ! dif* is already defined

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pifri = pir + difv(1)*gradv(isou,1,ii)      &
                          + difv(2)*gradv(isou,2,ii)      &
                          + difv(3)*gradv(isou,3,ii)
              pifrj = pi  + difv(1)*gradv(isou,1,ii)      &
                          + difv(2)*gradv(isou,2,ii)      &
                          + difv(3)*gradv(isou,3,ii)

              pjfrj = pjr + djfv(1)*gradv(isou,1,jj)      &
                          + djfv(2)*gradv(isou,2,jj)      &
                          + djfv(3)*gradv(isou,3,jj)
              pjfri = pj  + djfv(1)*gradv(isou,1,jj)      &
                          + djfv(2)*gradv(isou,2,jj)      &
                          + djfv(3)*gradv(isou,3,jj)

            endif

            ! Blending
            ! --------

            pifri = blencp*pifri+(1.d0-blencp)*pir
            pifrj = blencp*pifrj+(1.d0-blencp)*pif
            pjfri = blencp*pjfri+(1.d0-blencp)*pjf
            pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr

            ! Flux
            ! ----

            fluxi =   iconvp*(flui*pifri + fluj*pjfri - flumas(ifac)*pi)    &
                    + idiffp*viscf(ifac)*(pipr -pjp)
            fluxj =   iconvp*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj)    &
                    + idiffp*viscf(ifac)*(pip -pjpr)

            ! Assembly
            ! --------

            smbr(isou,ii) = smbr(isou,ii) - fluxi
            smbr(isou,jj) = smbr(isou,jj) + fluxj

          enddo ! isou

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, dijpfv, pnd,          &
      !$omp                     diipfv, djjpfv, flui, fluj, difv, djfv, dpvf,   &
      !$omp                     pi, pj, pip, pjp, pif, pjf, flux)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          do jsou = 1, 3
            dijpfv(jsou) = dijpf(jsou,ifac)
          enddo

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level
          do jsou = 1, 3
            diipfv(jsou) =   cdgfac(jsou,ifac) - (xyzcen(jsou,ii)               &
                           + (1.d0-pnd) * dijpfv(jsou))
            djjpfv(jsou) =   cdgfac(jsou,ifac) -  xyzcen(jsou,jj)               &
                           + pnd  * dijpfv(jsou)
          enddo

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          ! For second order, define IF
          if (ischcp.eq.0) then
            do jsou = 1, 3
              difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
              djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
            enddo
          endif

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            pi = pvar(isou,ii)
            pj = pvar(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            ! Centered
            ! --------

            if (ischcp.eq.1) then

              pif = pnd*pip +(1.d0-pnd)*pjp
              pjf = pif


            ! Second order
            ! ------------

            else ! if (ischcp.eq.0) then
              ! dif* is already defined

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pif = pi
              pjf = pj
              do jsou = 1, 3
                pif = pif + gradv(isou,jsou,ii)*difv(jsou)
                pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
              enddo

            endif

            ! Blending
            ! --------

            pif = blencp*pif+(1.d0-blencp)*pi
            pjf = blencp*pjf+(1.d0-blencp)*pj

            ! Flux
            ! ----

            flux =   iconvp*(flui*pif +fluj*pjf)          &
                   + idiffp*viscf(ifac)*(pip -pjp)

            ! Assembly
            ! --------

            smbr(isou,ii) = smbr(isou,ii) - thetap*(flux - iconvp*flumas(ifac)*pi)
            smbr(isou,jj) = smbr(isou,jj) + thetap*(flux - iconvp*flumas(ifac)*pj)

          enddo ! isou

        enddo
      enddo
    enddo

  endif

! --> Flux with slope test
! =========================

else

  if (ischcp.lt.0 .or. ischcp.gt.1) then
    write(nfecra,9000) ischcp
    call csexit(1)
  endif

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, dijpfv, pnd,     &
      !$omp                     distf, srfan, diipfv, djjpfv, flui, fluj,  &
      !$omp                     difv, djfv, dpvf, pi, pj, pia, pja,        &
      !$omp                     pip, pjp, pipr, pjpr, pir, pjr, testij,    &
      !$omp                     testi, testj, dcc, ddi, ddj, tesqck,       &
      !$omp                     pifri, pifrj, pjfri, pjfrj, fluxi, fluxj)  &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          do jsou = 1, 3
            dijpfv(jsou) = dijpf(jsou,ifac)
          enddo

          pnd   = pond(ifac)
          distf = dist(ifac)
          srfan = surfan(ifac)

          ! Recompute II' and JJ' at this level
          do jsou = 1, 3
            diipfv(jsou) =    cdgfac(jsou,ifac) - (xyzcen(jsou,ii)         &
                           + (1.d0-pnd) * dijpfv(jsou))
            djjpfv(jsou) =    cdgfac(jsou,ifac) -  xyzcen(jsou,jj)         &
                           + pnd * dijpfv(jsou)
          enddo

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          ! For second order, define IF
          if (ischcp.eq.0) then
            do jsou = 1, 3
              difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
              djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
            enddo
          endif

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            pi  = pvar (isou,ii)
            pj  = pvar (isou,jj)

            pia = pvara(isou,ii)
            pja = pvara(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            pipr = pi /relaxp - (1.d0-relaxp)/relaxp * pia   &
                 + ircflp*(dpvf(1)*diipfv(1)                 &
                          +dpvf(2)*diipfv(2)                 &
                          +dpvf(3)*diipfv(3))
            pjpr = pj /relaxp - (1.d0-relaxp)/relaxp * pja   &
                 + ircflp*(dpvf(1)*djjpfv(1)                 &
                          +dpvf(2)*djjpfv(2)                 &
                          +dpvf(3)*djjpfv(3))

            pir = pi /relaxp - (1.d0 - relaxp)/relaxp*pia
            pjr = pj /relaxp - (1.d0 - relaxp)/relaxp*pja

            ! Slope test
            ! ----------

            testi = gradva(isou,1,ii)*surfac(1,ifac)         &
                  + gradva(isou,2,ii)*surfac(2,ifac)         &
                  + gradva(isou,3,ii)*surfac(3,ifac)
            testj = gradva(isou,1,jj)*surfac(1,ifac)         &
                  + gradva(isou,2,jj)*surfac(2,ifac)         &
                  + gradva(isou,3,jj)*surfac(3,ifac)
            testij= gradva(isou,1,ii)*gradva(isou,1,jj)      &
                  + gradva(isou,2,ii)*gradva(isou,2,jj)      &
                  + gradva(isou,3,ii)*gradva(isou,3,jj)

            if (flumas(ifac).gt.0.d0) then
              dcc = gradv(isou,1,ii)*surfac(1,ifac)    &
                  + gradv(isou,2,ii)*surfac(2,ifac)    &
                  + gradv(isou,3,ii)*surfac(3,ifac)
              ddi = testi
              ddj = (pj - pi)/distf *srfan
            else
              dcc = gradv(isou,1,jj)*surfac(1,ifac)    &
                  + gradv(isou,2,jj)*surfac(2,ifac)    &
                  + gradv(isou,3,jj)*surfac(3,ifac)
              ddi = (pj - pi)/distf *srfan
              ddj = testj
            endif
            tesqck = dcc**2 -(ddi-ddj)**2

            ! Upwind
            ! ------

            if (tesqck.le.0.d0 .or. testij.le.0.d0) then

              pifri = pir
              pifrj = pi
              pjfri = pj
              pjfrj = pjr
              ! in parallel, face will be counted by one and only one rank
              if (ii.le.ncel) then
                infac = infac+1
              endif

            else

              ! Centered
              ! --------

              if (ischcp.eq.1) then

                pifri = pnd*pipr +(1.d0-pnd)*pjp
                pjfri = pifri
                pifrj = pnd*pip  +(1.d0-pnd)*pjpr
                pjfrj = pifrj

              ! Second order
              ! ------------

              else ! if (ischcp.eq.0) then
                ! difv already defined

                ! leave reconstruction of PIF and PJF even if IRCFLP=0
                ! otherwise, it is the same as using upwind
                pifri = pir + difv(1)*gradv(isou,1,ii)      &
                            + difv(2)*gradv(isou,2,ii)      &
                            + difv(3)*gradv(isou,3,ii)
                pifrj = pi  + difv(1)*gradv(isou,1,ii)      &
                            + difv(2)*gradv(isou,2,ii)      &
                            + difv(3)*gradv(isou,3,ii)

                pjfrj = pjr + djfv(1)*gradv(isou,1,jj)      &
                            + djfv(2)*gradv(isou,2,jj)      &
                            + djfv(3)*gradv(isou,3,jj)
                pjfri = pj  + djfv(1)*gradv(isou,1,jj)      &
                            + djfv(2)*gradv(isou,2,jj)      &
                            + djfv(3)*gradv(isou,3,jj)

              endif

            endif

            ! Blending
            ! --------

            pifri = blencp*pifri+(1.d0-blencp)*pir
            pifrj = blencp*pifrj+(1.d0-blencp)*pi
            pjfri = blencp*pjfri+(1.d0-blencp)*pj
            pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr

            ! Flux
            ! ----

            fluxi =    iconvp*(flui*pifri + fluj*pjfri - flumas(ifac)*pi)    &
                    + idiffp*viscf(ifac)*(pipr -pjp)
            fluxj =    iconvp*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj)    &
                    + idiffp*viscf(ifac)*(pip -pjpr)

            ! Assembly
            ! --------

            smbr(isou,ii) = smbr(isou,ii) - fluxi
            smbr(isou,jj) = smbr(isou,jj) + fluxj

          enddo ! isou

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, isou, jsou, dijpfv, pnd,     &
      !$omp                     distf, srfan, diipfv, djjpfv, flui, fluj,  &
      !$omp                     difv, djfv, dpvf, pi, pj, pip, pjp,        &
      !$omp                     testi, testj, testij, dcc, ddi, ddj,       &
      !$omp                     tesqck, pif, pjf, flux)                    &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          do jsou = 1, 3
            dijpfv(jsou) = dijpf(jsou,ifac)
          enddo

          pnd   = pond(ifac)
          distf = dist(ifac)
          srfan = surfan(ifac)

          ! Recompute II' and JJ' at this level
          do jsou = 1, 3
            diipfv(jsou) =    cdgfac(jsou,ifac) - (xyzcen(jsou,ii)         &
                           + (1.d0-pnd) * dijpfv(jsou))
            djjpfv(jsou) =    cdgfac(jsou,ifac) -  xyzcen(jsou,jj)         &
                           + pnd * dijpfv(jsou)
          enddo

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          ! For second order, define IF
          if (ischcp.eq.0) then
            do jsou = 1, 3
              difv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,ii)
              djfv(jsou) = cdgfac(jsou,ifac) - xyzcen(jsou,jj)
            enddo
          endif

          !-----------------
          ! X-Y-Z components, p=u, v, w
          do isou = 1, 3

            do jsou = 1, 3
              dpvf(jsou) = 0.5d0*(gradv(isou,jsou,ii) + gradv(isou,jsou,jj))
            enddo

            pi = pvar(isou,ii)
            pj = pvar(isou,jj)

            pip = pi + ircflp*(dpvf(1)*diipfv(1)        &
                              +dpvf(2)*diipfv(2)        &
                              +dpvf(3)*diipfv(3))
            pjp = pj + ircflp*(dpvf(1)*djjpfv(1)        &
                              +dpvf(2)*djjpfv(2)        &
                              +dpvf(3)*djjpfv(3))

            ! Slope test
            ! ----------

            testi = gradva(isou,1,ii)*surfac(1,ifac)    &
                  + gradva(isou,2,ii)*surfac(2,ifac)    &
                  + gradva(isou,3,ii)*surfac(3,ifac)
            testj = gradva(isou,1,jj)*surfac(1,ifac)    &
                  + gradva(isou,2,jj)*surfac(2,ifac)    &
                  + gradva(isou,3,jj)*surfac(3,ifac)
            testij = gradva(isou,1,ii)*gradva(isou,1,jj) &
                   + gradva(isou,2,ii)*gradva(isou,2,jj) &
                   + gradva(isou,3,ii)*gradva(isou,3,jj)

            if (flumas(ifac).gt.0.d0) then
              dcc = gradv(isou,1,ii)*surfac(1,ifac)     &
                  + gradv(isou,2,ii)*surfac(2,ifac)     &
                  + gradv(isou,3,ii)*surfac(3,ifac)
              ddi = testi
              ddj = (pj - pi)/distf *srfan
            else
              dcc = gradv(isou,1,jj)*surfac(1,ifac)    &
                  + gradv(isou,2,jj)*surfac(2,ifac)    &
                  + gradv(isou,3,jj)*surfac(3,ifac)
              ddi = (pj - pi)/distf *srfan
              ddj = testj
            endif
            tesqck = dcc**2 -(ddi-ddj)**2

            ! Upwind
            ! ------

            if (tesqck.le.0.d0 .or. testij.le.0.d0) then

              pif = pi
              pjf = pj
              ! in parallel, face will be counted by one and only one rank
              if (ii.le.ncel) then
                infac = infac+1
              endif

            else

              ! Centered
              ! --------

              if (ischcp.eq.1) then

                pif = pnd*pip +(1.d0-pnd)*pjp
                pjf = pif

              ! Second order
              ! ------------

              else ! if (ischcp.eq.0) then
                ! dif* is already defined

                pif = pi
                pjf = pj
                do jsou = 1, 3
                  pif = pif + gradv(isou,jsou,ii)*difv(jsou)
                  pjf = pjf + gradv(isou,jsou,jj)*djfv(jsou)
                enddo

                ! leave reconstruction of PIF and PJF even if IRCFLP=0
                ! otherwise, it is the same as using upwind

              endif

            endif

            ! Blending
            ! --------

            pif = blencp*pif+(1.d0-blencp)*pi
            pjf = blencp*pjf+(1.d0-blencp)*pj

            ! Flux
            ! ----

            flux =   iconvp*(flui*pif +fluj*pjf)       &
                   + idiffp*viscf(ifac)*(pip -pjp)

            ! Assembly
            ! --------

            smbr(isou,ii) = smbr(isou,ii) - thetap*(flux - iconvp*flumas(ifac)*pi)
            smbr(isou,jj) = smbr(isou,jj) + thetap*(flux - iconvp*flumas(ifac)*pj)

          enddo ! isou

        enddo
      enddo
    enddo

  endif ! idtvar

endif ! iupwin


if (iwarnp.ge.2) then
  if (irangp.ge.0) call parcpt(infac)
  write(nfecra,1100) cnom, infac, nfacgb
endif

! ======================================================================
! ---> Contribution from boundary faces
! ======================================================================

! Steady
if (idtvar.lt.0) then

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, isou, jsou, diipbv, flui, fluj,         &
    !$omp                     pfac, pfacd, pir, pipr, flux, pi, pia)            &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        do jsou = 1, 3
          diipbv(jsou) = diipb(jsou,ifac)
        enddo

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          pfac  = inc*coefav(isou,ifac)
          pfacd = inc*cofafv(isou,ifac)

          !coefu and cofuf are matrices
          do jsou = 1, 3
            pir  = pvar(jsou,ii)/relaxp - (1.d0-relaxp)/relaxp*pvara(jsou,ii)

            pipr = pir +ircflp*( gradv(jsou,1,ii)*diipbv(1)         &
                               + gradv(jsou,2,ii)*diipbv(2)         &
                               + gradv(jsou,3,ii)*diipbv(3))
            pfac  = pfac  + coefbv(isou,jsou,ifac)*pipr
            pfacd = pfacd + cofbfv(isou,jsou,ifac)*pipr
          enddo

          pi  = pvar(isou,ii)
          pia = pvara(isou,ii)

          pir  = pi/relaxp - (1.d0-relaxp)/relaxp*pia
          pipr = pir +ircflp*( gradv(isou,1,ii)*diipbv(1)           &
                             + gradv(isou,2,ii)*diipbv(2)           &
                             + gradv(isou,3,ii)*diipbv(3))

          flux = iconvp*(flui*pir + fluj*pfac - flumab(ifac)*pi)     &
               + idiffp*viscb(ifac)*pfacd
          smbr(isou,ii) = smbr(isou,ii) - flux

        enddo ! isou

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, isou, jsou, diipbv, flui, fluj,         &
    !$omp                     pfac, pfacd, pip, flux, pi)                       &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        do jsou = 1, 3
          diipbv(jsou) = diipb(jsou,ifac)
        enddo

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        !-----------------
        ! X-Y-Z components, p=u, v, w
        do isou = 1, 3

          pfac  = inc*coefav(isou,ifac)
          pfacd = inc*cofafv(isou,ifac)

          !coefu and cofuf are matrices
          do jsou = 1, 3
            pip = pvar(jsou,ii) + ircflp*( gradv(jsou,1,ii)*diipbv(1)        &
                                         + gradv(jsou,2,ii)*diipbv(2)        &
                                         + gradv(jsou,3,ii)*diipbv(3))
            pfac  = pfac  + coefbv(isou,jsou,ifac)*pip
            pfacd = pfacd + cofbfv(isou,jsou,ifac)*pip
          enddo

          pi = pvar(isou,ii)

          flux = iconvp*((flui-flumab(ifac))*pi + fluj*pfac)                 &
               + idiffp*viscb(ifac)*pfacd
          smbr(isou,ii) = smbr(isou,ii) - thetap * flux

        enddo ! isou

      enddo
    enddo
  enddo

endif ! idtvar

!===============================================================================
! 3.  Computation of the transpose grad(vel) term and grad(-2/3 div(vel))
!===============================================================================

if (ivisep.eq.1) then

  ! We do not know what condition to put in the inlets and the outlets, so we
  ! assume that there is an equilibrium. Moreover, cells containing a coupled
  ! are removed.

  ! Allocate a temporary array
  allocate(bndcel(ncelet))

  !$omp parallel do
  do iel = 1, ncelet
    bndcel(iel) = 1.d0
  enddo

  !$omp parallel do private(ityp) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    ityp = itypfb(ifac)
    if (    ityp.eq.isolib                          &
       .or. ityp.eq.ientre                          &
       .or.(ifaccp.eq.1.and.ityp.eq.icscpl)) bndcel(ifabor(ifac)) = 0.d0
  enddo

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(bndcel)
    !==========
  endif

  ! ---> Interior faces

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, isou, jsou, pnd, secvis,     &
    !$omp                     visco, grdtrv, tgrdfl, flux)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        pnd = pond(ifac)
        secvis = secvif(ifac)
        visco = viscf(ifac)

        grdtrv =        pnd*(gradv(1,1,ii)+gradv(2,2,ii)+gradv(3,3,ii))   &
               + (1.d0-pnd)*(gradv(1,1,jj)+gradv(2,2,jj)+gradv(3,3,jj))

        ! We need to compute trans_grad(u).IJ which is equal to IJ.grad(u)

        do isou = 1, 3

          tgrdfl = dijpf(1,ifac) * (        pnd*gradv(1,isou,ii)         &
                                   + (1.d0-pnd)*gradv(1,isou,jj))        &
                 + dijpf(2,ifac) * (        pnd*gradv(2,isou,ii)         &
                                   + (1.d0-pnd)*gradv(2,isou,jj))        &
                 + dijpf(3,ifac) * (        pnd*gradv(3,isou,ii)         &
                                   + (1.d0-pnd)*gradv(3,isou,jj))

          flux = visco*tgrdfl + secvis*grdtrv*surfac(isou,ifac)

          smbr(isou,ii) = smbr(isou,ii) + flux*bndcel(ii)
          smbr(isou,jj) = smbr(isou,jj) - flux*bndcel(jj)

        enddo

      enddo
    enddo
  enddo

  ! ---> Boundary FACES
  !      the whole flux term of the stress tensor is already taken into account
  !      (so, no corresponding term in forbr)
  !TODO in theory we should take the normal component into account (the
  !tangential one is modeled by the wall law)

  !Free memory
  deallocate(bndcel)

endif

! Free memory
deallocate(gradva)
deallocate(gradv)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(1X,A8,' : CONVECTION EN ',A11,                       &
                               ' BLENDING A ',F4.0,' % D''UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES UPWIND SUR ',                &
                               I10,' FACES INTERNES ')
 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS bilsc4                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE bilsc4 POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(1X,A8,' : CONVECTION IN ',A11,                       &
                            ' BLENDING WITH ',F4.0,' % OF UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES WITH UPWIND ON ',            &
                               I10,' INTERIOR FACES ')
 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN bilsc4                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF bilsc4 FOR ',A8 ,' WITH ISCHCP = ',I10          ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return

end subroutine
