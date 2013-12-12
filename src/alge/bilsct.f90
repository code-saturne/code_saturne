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

!> \file bilsct.f90
!>
!> \brief This function adds the explicit part of the convection/diffusion
!> terms of a transport equation of a scalar field \f$ \varia \f$ such as the
!> temperature.
!>
!> More precisely, the right hand side \f$ Rhs \f$ is updated as
!> follows:
!> \f[
!> Rhs = Rhs + \sum_{\fij \in \Facei{\celli}}      \left(
!>        C_p\dot{m}_\ij \varia_\fij
!>      - \lambda_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
!> \f]
!>
!> Warning:
!> \f$ Rhs \f$ has already been initialized before calling bilsct!
!>
!> Options for the convective scheme:
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
!> \param[in]     iccocg        indicator
!>                               - 1 re-compute cocg matrix (for iterativ gradients)
!>                               - 0 otherwise
!> \param[in]     ipp           index of the variable for post-processing
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
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
!> \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$
!_______________________________________________________________________________

subroutine bilsct &
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscf  , viscb  , xcpp   ,                   &
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
use mesh
use field
use numvar, only: ivarfl

!===============================================================================

implicit none

! Arguments

integer          idtvar
integer          ivar   , iconvp , idiffp , nswrgp , imligp
integer          ircflp , ischcp , isstpp
integer          inc    , imrgra , iccocg
integer          iwarnp , ipp

double precision blencp , epsrgp , climgp, extrap, relaxp , thetap

double precision pvar (ncelet), pvara(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf (nfac), viscb (nfabor)
double precision smbrp(ncelet)
double precision xcpp(ncelet)

! Local variables

character*80     chaine
character*8      cnom
integer          ifac,ii,jj,infac,iel,iupwin, ig, it
double precision pfac,pfacd,flui,fluj,flux,fluxi,fluxj
double precision difx,dify,difz,djfx,djfy,djfz
double precision pi, pj, pia, pja
double precision pif,pjf,pip,pjp,pir,pjr,pipr,pjpr
double precision pifri,pifrj,pjfri,pjfrj
double precision testi,testj,testij
double precision dpxf,dpyf,dpzf
double precision dcc, ddi, ddj, tesqck
double precision dijpfx, dijpfy, dijpfz
double precision diipfx, diipfy, diipfz
double precision djjpfx, djjpfy, djjpfz
double precision diipbx, diipby, diipbz
double precision pnd, distf, srfan
double precision pfac1, pfac2, pfac3, unsvol

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpdxa, dpdya, dpdza

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate work arrays
allocate(grad(ncelet,3))
allocate(dpdxa(ncelet), dpdya(ncelet), dpdza(ncelet))

! Initialize variables to avoid compiler warnings

pif = 0.d0
pjf = 0.d0
pifri = 0.d0
pifrj = 0.d0
pjfri = 0.d0
pjfrj = 0.d0

! Memoire

if (ivar.gt.0) then
  call field_get_name(ivarfl(ivar), chaine)
else
  chaine = nomvar(ipp)
endif
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
! ---> CALCUL DU GRADIENT DE P
! ======================================================================
!    GRAD sert a la fois pour la reconstruction des flux et pour le test
!    de pente. On doit donc le calculer :
!        - quand on a de la diffusion et qu'on reconstruit les flux
!        - quand on a de la convection SOLU
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on reconstruit les flux
!        - quand on a de la convection, qu'on n'est pas en upwind pur
!          et qu'on n'a pas shunte le test de pente

if ((idiffp.ne.0 .and. ircflp.eq.1) .or.                          &
    (iconvp.ne.0 .and. iupwin.eq.0 .and.                          &
       (ischcp.eq.0 .or. ircflp.eq.1 .or. isstpp.eq.0))) then

  call grdcel                                                     &
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
! ---> CALCUL DU GRADIENT DECENTRE DPDXA, DPDYA, DPDZA POUR TST DE PENTE
! ======================================================================

!$omp parallel do
do iel = 1, ncelet
  dpdxa(iel) = 0.d0
  dpdya(iel) = 0.d0
  dpdza(iel) = 0.d0
enddo

if (iconvp.gt.0.and.iupwin.eq.0.and.isstpp.eq.0) then

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, difx, dify, difz, djfx, djfy, djfz, &
    !$omp                     pif, pjf, pfac, pfac1, pfac2, pfac3)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        difx = cdgfac(1,ifac) - xyzcen(1,ii)
        dify = cdgfac(2,ifac) - xyzcen(2,ii)
        difz = cdgfac(3,ifac) - xyzcen(3,ii)
        djfx = cdgfac(1,ifac) - xyzcen(1,jj)
        djfy = cdgfac(2,ifac) - xyzcen(2,jj)
        djfz = cdgfac(3,ifac) - xyzcen(3,jj)

        pif = pvar(ii) + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
        pjf = pvar(jj) + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

        pfac = pjf
        if (flumas(ifac).gt.0.d0) pfac = pif

        pfac1 = pfac*surfac(1,ifac)
        pfac2 = pfac*surfac(2,ifac)
        pfac3 = pfac*surfac(3,ifac)

        dpdxa(ii) = dpdxa(ii) + pfac1
        dpdya(ii) = dpdya(ii) + pfac2
        dpdza(ii) = dpdza(ii) + pfac3

        dpdxa(jj) = dpdxa(jj) - pfac1
        dpdya(jj) = dpdya(jj) - pfac2
        dpdza(jj) = dpdza(jj) - pfac3

      enddo
    enddo
  enddo

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, pfac) &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)
        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)
        pfac =   inc*coefap(ifac)                                           &
               + coefbp(ifac) * (  pvar(ii)          + diipbx*grad(ii,1)    &
                                 + diipby*grad(ii,2) + diipbz*grad(ii,3))
        dpdxa(ii) = dpdxa(ii) + pfac*surfbo(1,ifac)
        dpdya(ii) = dpdya(ii) + pfac*surfbo(2,ifac)
        dpdza(ii) = dpdza(ii) + pfac*surfbo(3,ifac)

      enddo
    enddo
  enddo

  !$omp parallel do private(unsvol)
  do iel = 1, ncel
    unsvol = 1.d0/volume(iel)
    dpdxa(iel) = dpdxa(iel)*unsvol
    dpdya(iel) = dpdya(iel)*unsvol
    dpdza(iel) = dpdza(iel)*unsvol
  enddo

  ! Synchronization for parallelism or periodicity

  if (irangp.ge.0 .or. iperio.eq.1) then
    call synvec(dpdxa, dpdya, dpdza)
    !==========
  endif

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

! --> Pure upwind flux
! =====================

if (iupwin.eq.1) then

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, pipr, pjpr,         &
      !$omp                     flui, fluj, pif, pjf, fluxi, fluxj,             &
      !$omp                     pi, pj, pir, pjr, pia, pja)                     &
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

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          ! reconstruction only if IRCFLP = 1
          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0-relaxp)/relaxp * pia
          pjr = pj/relaxp - (1.d0-relaxp)/relaxp * pja

          pipr = pir                                                   &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                   &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          pif  = pi
          pjf  = pj

          fluxi = iconvp*xcpp(ii)*(flui*pir + fluj*pjf - flumas(ifac)*pi)       &
                + idiffp*viscf(ifac)*(pipr - pjp)
          fluxj = iconvp*xcpp(jj)*(flui*pif + fluj*pjr - flumas(ifac)*pj)       &
                + idiffp*viscf(ifac)*(pip - pjpr)

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, flui, fluj,         &
      !$omp                     pif, pjf, fluxi, fluxj, pi, pj)                 &
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

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          pif = pi
          pjf = pj

          fluxi = iconvp*xcpp(ii)*(flui*pif +fluj*pjf -flumas(ifac)*pi) &
                + idiffp*viscf(ifac)*(pip -pjp)
          fluxj = iconvp*xcpp(jj)*(flui*pif +fluj*pjf -flumas(ifac)*pj) &
                + idiffp*viscf(ifac)*(pip -pjp)

          smbrp(ii) = smbrp(ii) - thetap * fluxi
          smbrp(jj) = smbrp(jj) + thetap * fluxj

        enddo
      enddo
    enddo

  endif


! --> Flux with no slope test
! ============================

else if (isstpp.eq.1) then

  if (ischcp.lt.0 .or. ischcp.gt.1) then
    write(nfecra,9000) ischcp
    call csexit(1)
  endif

  ! Steady
  if (idtvar.lt.0) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, pipr, pjpr, flui,   &
      !$omp                     fluj, pir, pjr, pifri, pjfri, pifrj, pjfrj,     &
      !$omp                     difx, dify, difz, djfx, djfy, djfz,             &
      !$omp                     fluxi, fluxj, pi, pj, pia, pja)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          pi = pvar(ii)
          pj = pvar(jj)
          pia = pvara(ii)
          pja = pvara(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0 - relaxp)/relaxp*pia
          pjr = pj/relaxp - (1.d0 - relaxp)/relaxp*pja

          pipr = pir                                                    &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                    &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))


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

            difx = cdgfac(1,ifac) - xyzcen(1,ii)
            dify = cdgfac(2,ifac) - xyzcen(2,ii)
            difz = cdgfac(3,ifac) - xyzcen(3,ii)
            djfx = cdgfac(1,ifac) - xyzcen(1,jj)
            djfy = cdgfac(2,ifac) - xyzcen(2,jj)
            djfz = cdgfac(3,ifac) - xyzcen(3,jj)

            ! leave reconstruction of PIF and PJF even if IRCFLP=0
            ! otherwise, it is the same as using upwind
            pifri = pir + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pifrj = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pjfrj = pjr + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)
            pjfri = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

          endif

          ! Blending
          ! --------

          pifri = blencp*pifri+(1.d0-blencp)*pir
          pifrj = blencp*pifrj+(1.d0-blencp)*pi
          pjfri = blencp*pjfri+(1.d0-blencp)*pj
          pjfrj = blencp*pjfrj+(1.d0-blencp)*pjr

          ! Flux
          ! ----

          fluxi =   iconvp*xcpp(ii)*(flui*pifri + fluj*pjfri - flumas(ifac)*pi) &
                  + idiffp*viscf(ifac)*(pipr -pjp)
          fluxj =   iconvp*xcpp(jj)*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj) &
                  + idiffp*viscf(ifac)*(pip -pjpr)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz, &
      !$omp                     dpxf, dpyf, dpzf, pip, pjp, flui, fluj, pif,    &
      !$omp                     pjf, difx, dify, difz, djfx, djfy, djfz, fluxi, &
      !$omp                     fluxj, pi, pj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd = pond(ifac)

          pi = pvar(ii)
          pj = pvar(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) + abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) - abs(flumas(ifac)))

          ! Centered
          ! --------

          if (ischcp.eq.1) then

            pif = pnd*pip +(1.d0-pnd)*pjp
            pjf = pif

          ! Second order
          ! ------------

          else ! if (ischcp.eq.0) then

            difx = cdgfac(1,ifac) - xyzcen(1,ii)
            dify = cdgfac(2,ifac) - xyzcen(2,ii)
            difz = cdgfac(3,ifac) - xyzcen(3,ii)
            djfx = cdgfac(1,ifac) - xyzcen(1,jj)
            djfy = cdgfac(2,ifac) - xyzcen(2,jj)
            djfz = cdgfac(3,ifac) - xyzcen(3,jj)

            ! leave reconstruction of PIF and PJF even if IRCFLP=0
            ! otherwise, it is the same as using upwind
            pif = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
            pjf = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

          endif

          ! Blending
          ! --------

          pif = blencp*pif+(1.d0-blencp)*pi
          pjf = blencp*pjf+(1.d0-blencp)*pj

          ! Flux
          ! ----

          fluxi = iconvp*xcpp(ii)*(flui*pif +fluj*pjf - flumas(ifac)*pi) &
                + idiffp*viscf(ifac)*(pip -pjp)
          fluxj = iconvp*xcpp(jj)*(flui*pif +fluj*pjf - flumas(ifac)*pj) &
                + idiffp*viscf(ifac)*(pip -pjp)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - thetap * fluxi
          smbrp(jj) = smbrp(jj) + thetap * fluxj

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
      !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,      &
      !$omp                     distf, srfan, diipfx, diipfy, diipfz, djjpfx,   &
      !$omp                     djjpfy, djjpfz, dpxf, dpyf, dpzf, pip, pjp,     &
      !$omp                     pipr, pjpr, flui, fluj, pir, pjr, testi, testj, &
      !$omp                     testij, dcc, ddi, ddj, tesqck, pifri, pjfri,    &
      !$omp                     pifrj, pjfrj, difx, dify, difz, djfx, djfy,     &
      !$omp                     djfz, fluxi, fluxj, pi, pj, pia, pja)           &
      !$omp             reduction(+:infac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd    = pond(ifac)
          distf  = dist(ifac)
          srfan  = surfan(ifac)

          pi = pvar(ii)
          pj = pvar(jj)
          pia = pvara(ii)
          pja = pvara(jj)

          ! Recompute II' and JJ' at this level
          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          pir = pi/relaxp - (1.d0 - relaxp)/relaxp*pia
          pjr = pj/relaxp - (1.d0 - relaxp)/relaxp*pja

          pipr = pir                                                  &
               + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjpr = pjr                                                  &
               + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))


          ! Slope test
          ! ----------

          testi =   dpdxa(ii)*surfac(1,ifac) + dpdya(ii)*surfac(2,ifac)  &
                  + dpdza(ii)*surfac(3,ifac)
          testj =   dpdxa(jj)*surfac(1,ifac) + dpdya(jj)*surfac(2,ifac)  &
                  + dpdza(jj)*surfac(3,ifac)
          testij=   dpdxa(ii)*dpdxa(jj)      + dpdya(ii)*dpdya(jj)       &
                  + dpdza(ii)*dpdza(jj)

          if (flumas(ifac).gt.0.d0) then
            dcc =   grad(ii,1)*surfac(1,ifac) +grad(ii,2)*surfac(2,ifac)    &
                  + grad(ii,3)*surfac(3,ifac)
            ddi = testi
            ddj = (pj-pi)/distf *srfan
          else
            dcc =   grad(jj,1)*surfac(1,ifac) +grad(jj,2)*surfac(2,ifac)    &
                  + grad(jj,3)*surfac(3,ifac)
            ddi = (pj-pi)/distf *srfan
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

              difx = cdgfac(1,ifac) - xyzcen(1,ii)
              dify = cdgfac(2,ifac) - xyzcen(2,ii)
              difz = cdgfac(3,ifac) - xyzcen(3,ii)
              djfx = cdgfac(1,ifac) - xyzcen(1,jj)
              djfy = cdgfac(2,ifac) - xyzcen(2,jj)
              djfz = cdgfac(3,ifac) - xyzcen(3,jj)

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pifri = pir + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pifrj = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pjfrj = pjr + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)
              pjfri = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

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

          fluxi = iconvp*xcpp(ii)*(flui*pifri + fluj*pjfri - flumas(ifac)*pi) &
                + idiffp*viscf(ifac)*(pipr -pjp)
          fluxj = iconvp*xcpp(jj)*(flui*pifrj + fluj*pjfrj - flumas(ifac)*pj) &
                + idiffp*viscf(ifac)*(pip -pjpr)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - fluxi
          smbrp(jj) = smbrp(jj) + fluxj

        enddo
      enddo
    enddo

  ! Unsteady
  else

    do ig = 1, ngrpi
       !$omp parallel do private(ifac, ii, jj, dijpfx, dijpfy, dijpfz, pnd,     &
       !$omp                     distf, srfan, diipfx, diipfy, diipfz, djjpfx,  &
       !$omp                     djjpfy, djjpfz, dpxf, dpyf, dpzf, pip, pjp,    &
       !$omp                     flui, fluj, testi, testj, testij, dcc, ddi,    &
       !$omp                     ddj, tesqck, pif, pjf, difx, dify, difz,       &
       !$omp                     djfx, djfy, djfz, fluxi, fluxj, pi, pj)        &
       !$omp             reduction(+:infac)
       do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)

          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)

          dijpfx = dijpf(1,ifac)
          dijpfy = dijpf(2,ifac)
          dijpfz = dijpf(3,ifac)

          pnd    = pond(ifac)
          distf  = dist(ifac)
          srfan  = surfan(ifac)

          pi = pvar(ii)
          pj = pvar(jj)

          ! Recompute II' and JJ' at this level

          diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
          diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
          diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
          djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) + pnd * dijpfx
          djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) + pnd * dijpfy
          djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) + pnd * dijpfz

          dpxf = 0.5d0*(grad(ii,1) + grad(jj,1))
          dpyf = 0.5d0*(grad(ii,2) + grad(jj,2))
          dpzf = 0.5d0*(grad(ii,3) + grad(jj,3))

          pip = pi + ircflp*(dpxf*diipfx+dpyf*diipfy+dpzf*diipfz)
          pjp = pj + ircflp*(dpxf*djjpfx+dpyf*djjpfy+dpzf*djjpfz)

          flui = 0.5d0*(flumas(ifac) +abs(flumas(ifac)))
          fluj = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))

          ! Slope test
          ! ----------

          testi =   dpdxa(ii)*surfac(1,ifac) + dpdya(ii)*surfac(2,ifac)  &
                  + dpdza(ii)*surfac(3,ifac)
          testj =   dpdxa(jj)*surfac(1,ifac) + dpdya(jj)*surfac(2,ifac)  &
                  + dpdza(jj)*surfac(3,ifac)
          testij =   dpdxa(ii)*dpdxa(jj)    + dpdya(ii)*dpdya(jj)        &
                   + dpdza(ii)*dpdza(jj)

          if (flumas(ifac).gt.0.d0) then
            dcc =   grad(ii,1)*surfac(1,ifac) + grad(ii,2)*surfac(2,ifac)    &
                  + grad(ii,3)*surfac(3,ifac)
            ddi = testi
            ddj = (pj-pi)/distf *srfan
          else
            dcc =   grad(jj,1)*surfac(1,ifac) + grad(jj,2)*surfac(2,ifac)    &
                  + grad(jj,3)*surfac(3,ifac)
            ddi = (pj-pi)/distf *srfan
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

              difx = cdgfac(1,ifac) - xyzcen(1,ii)
              dify = cdgfac(2,ifac) - xyzcen(2,ii)
              difz = cdgfac(3,ifac) - xyzcen(3,ii)
              djfx = cdgfac(1,ifac) - xyzcen(1,jj)
              djfy = cdgfac(2,ifac) - xyzcen(2,jj)
              djfz = cdgfac(3,ifac) - xyzcen(3,jj)

              ! leave reconstruction of PIF and PJF even if IRCFLP=0
              ! otherwise, it is the same as using upwind
              pif = pi + difx*grad(ii,1)+dify*grad(ii,2)+difz*grad(ii,3)
              pjf = pj + djfx*grad(jj,1)+djfy*grad(jj,2)+djfz*grad(jj,3)

            endif

          endif

          ! Blending
          ! --------

          pif = blencp*pif+(1.d0-blencp)*pi
          pjf = blencp*pjf+(1.d0-blencp)*pj

          ! Flux
          ! ----

          fluxi = iconvp*xcpp(ii)*(flui*pif + fluj*pjf - flumas(ifac)*pi)  &
                + idiffp*viscf(ifac)*(pip-pjp)
          fluxj = iconvp*xcpp(jj)*(flui*pif + fluj*pjf - flumas(ifac)*pj)  &
                + idiffp*viscf(ifac)*(pip-pjp)

          ! Assembly
          ! --------

          smbrp(ii) = smbrp(ii) - thetap *fluxi
          smbrp(jj) = smbrp(jj) + thetap *fluxj

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
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, flui, fluj,     &
    !$omp                     pir, pipr, pfac, pfacd, flux, pi, pia)            &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)
        pia = pvara(ii)

        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        pir  = pi/relaxp - (1.d0-relaxp)/relaxp*pia
        pipr = pir                                                            &
             + ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

        pfac  = inc*coefap(ifac) +coefbp(ifac)*pipr
        pfacd = inc*cofafp(ifac) +cofbfp(ifac)*pipr

        flux = iconvp*xcpp(ii)*(flui*pir + fluj*pfac - flumab(ifac)*pi )      &
             + idiffp*viscb(ifac)*pfacd
        smbrp(ii) = smbrp(ii) - flux

      enddo
    enddo
  enddo

! Unsteady
else

  do ig = 1, ngrpb
    !$omp parallel do private(ifac, ii, diipbx, diipby, diipbz, flui, fluj,     &
    !$omp                     pip, pfac, pfacd, flux, pi)                       &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

        ii = ifabor(ifac)

        pi = pvar(ii)

        diipbx = diipb(1,ifac)
        diipby = diipb(2,ifac)
        diipbz = diipb(3,ifac)

        ! Remove decentering for coupled faces
        if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
          flui = 0.0d0
          fluj = flumab(ifac)
        else
          flui = 0.5d0*(flumab(ifac) +abs(flumab(ifac)))
          fluj = 0.5d0*(flumab(ifac) -abs(flumab(ifac)))
        endif

        pip = pi                                                       &
            + ircflp*(grad(ii,1)*diipbx+grad(ii,2)*diipby+grad(ii,3)*diipbz)

        pfac  = inc*coefap(ifac) + coefbp(ifac)*pip
        pfacd = inc*cofafp(ifac) + cofbfp(ifac)*pip

        flux = iconvp*xcpp(ii)*((flui - flumab(ifac))*pi + fluj*pfac)         &
             + idiffp*viscb(ifac)*pfacd
        smbrp(ii) = smbrp(ii) - thetap * flux

      enddo
    enddo
  enddo

endif

! Free memory
deallocate(grad)
deallocate(dpdxa, dpdya, dpdza)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(1X,A8,' : CONVECTION EN ',A11,                             &
                               ' BLENDING A ',F4.0,' % D''UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES UPWIND SUR ',                      &
                               I10,' FACES INTERNES ')
 9000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS bilsct                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE bilsct POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(1X,A8,' : CONVECTION IN ',A11,                             &
                            ' BLENDING WITH ',F4.0,' % OF UPWIND')
 1100 format(1X,A8,' : ',I10,' FACES WITH UPWIND ON ',                  &
                               I10,' INTERIOR FACES ')
 9000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN bilsct                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF bilsct FOR ',A8 ,' WITH ISCHCP = ',I10          ,/,&
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
