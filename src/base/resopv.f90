!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file resopv.f90
!>
!> \brief This subroutine performs the pressure correction step of the Navier
!> Stokes equations for incompressible or slightly compressible flows for
!> the coupled velocity components solver.
!>
!> This function solves the following Poisson equation on the pressure:
!> \f[
!>     D \left( \Delta t, \delta p \right) =
!> \divs \left( \rho \vect{\widetilde{u}}\right)
!>     - \Gamma^n
!>     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
!> \f]
!> The mass flux is then updated as follows:
!> \f[
!>  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
!>                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
!> \f]
!>
!> Remarks:
!> - an iterative process is used to solve the Poisson equation.
!> - if the coefficient arak is set to 1, the the Rhie & Chow filter is
!>   activated.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     isostd        indicator of standard outlet and index
!>                               of the reference outlet face
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     coefa_dp      boundary conditions for the pressure increment
!> \param[in]     coefb_dp      boundary conditions for the pressure increment
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     dfrcxt        variation of the external forces
!> \param[in]                    making the hydrostatic pressure
!> \param[in]     tpucou        non scalar time step in case of
!>                               velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                               the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     drtp          tableau de travail pour increment
!> \param[in]     tslagr        coupling term for the Lagrangian module
!> \param[in]     trava         tableau de travail pour couplage
!_______________________________________________________________________________

subroutine resopv &
 ( nvar   , ncesmp ,                                              &
   icetsm , isostd ,                                              &
   dt     , rtp    , rtpa   , vel    ,                            &
   propce ,                                                       &
   coefav , coefbv , coefa_dp        , coefb_dp ,                 &
   smacel ,                                                       &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  ,                                              &
   drtp   , tslagr ,                                              &
   trava  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use cstphy
use cstnum
use optcal
use pointe, only: itypfb, b_head_loss
use albase
use parall
use period
use mltgrd
use lagpar
use lagran
use cplsat
use mesh
use field
use field_operator
use cs_f_interfaces

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          ncesmp

integer          icetsm(ncesmp)
integer          isostd(nfabor+1)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision smacel(ncesmp,nvar)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision tpucou(6, ncelet), trav(3,ncelet)
double precision viscf(nfac), viscb(ndimfb)
double precision drtp(ncelet)
double precision tslagr(ncelet,*)
double precision trava(ndim,ncelet)
double precision coefav(3  ,ndimfb)
double precision coefbv(3,3,ndimfb)
double precision vel   (3  ,ncelet)
double precision coefa_dp(ndimfb)
double precision coefb_dp(ndimfb)

! Local variables

character*80     chaine
integer          lchain
integer          iccocg, inc   , iprev, init  , isym  , ipol  , isqrt
integer          ii, iel   , ifac  , ifac0 , iel0
integer          ireslp, nswmpr
integer          isweep, niterf, icycle
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          iflmas, iflmab
integer          ipp
integer          idiffp, iconvp, ndircp
integer          nitmap, imgrp , ncymap, nitmgp
integer          iinvpe, indhyd
integer          itypfl
integer          iesdep
integer          nagmax, npstmg
integer          isou  , ibsize, iesize
integer          imucpp, idftnp, iswdyp
integer          iescap, ircflp, ischcp, isstpp, ivar, ncymxp, nitmfp
integer          nswrsp
integer          insqrt

integer          icvflb, f_id0
integer          ivoid(1)

double precision residu, phydr0
double precision ardtsr, arsr  , unsara, thetap
double precision dtsrom, unsvom, romro0
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1, relaxp
double precision hint, qimp, qimpv(3), epsrsp, blencp
double precision ressol, rnorm2
double precision nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk, alph, beta
double precision visci(3,3), fikis, viscis, distfi
double precision cfl, kpdc, rho, pimp, bpmasf

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam, xam
double precision, allocatable, dimension(:) :: res, divu, presa
double precision, dimension(:,:), allocatable :: gradp
double precision, allocatable, dimension(:) :: coefaf_dp, coefbf_dp
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp
double precision, allocatable, dimension(:) :: rhs, rovsdt
double precision, allocatable, dimension(:) :: velflx, velflb, dpvar
double precision, allocatable, dimension(:,:) :: coefar, cofafr
double precision, allocatable, dimension(:,:,:) :: coefbr, cofbfr
double precision, allocatable, dimension(:) :: adxk, adxkm1, dpvarm1, rhs0
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: frchy, dfrchy
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:), pointer :: coefaf_p, coefbf_p
double precision, allocatable, dimension(:) :: iflux, bflux
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, croma

!===============================================================================

!===============================================================================
! 1. Initialisations
!===============================================================================

! Initializations to avoid compiler warnings
rnorm2 = 0.d0

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac))
allocate(res(ncelet), presa(ncelet), divu(ncelet))
allocate(rhs(ncelet), rovsdt(ncelet))
allocate(iflux(nfac), bflux(ndimfb))
iswdyp = iswdyn(ipr)
if (iswdyp.ge.1) allocate(adxk(ncelet), adxkm1(ncelet),   &
                          dpvarm1(ncelet), rhs0(ncelet))
if (icalhy.eq.1) allocate(frchy(ndim,ncelet), dfrchy(ndim,ncelet))

! Diffusive flux Boundary conditions for delta P
allocate(coefaf_dp(ndimfb), coefbf_dp(ndimfb))

! --- Writing
call field_get_name(ivarfl(ipr), chaine)
lchain = 16
ipp    = ipprtp(ipr)

f_id0 = -1

! --- Boundary conditions

call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)
call field_get_coefaf_s(ivarfl(ipr), coefaf_p)
call field_get_coefbf_s(ivarfl(ipr), coefbf_p)

! --- Physical quantities
call field_get_val_s(icrom, crom)
if (icalhy.eq.1.or.idilat.gt.1) then
  call field_get_val_prev_s(icrom, croma)
endif
call field_get_val_s(ibrom, brom)

call field_get_key_int(ivarfl(ipr), kimasf, iflmas)
call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! --- Solving options
isym  = 1
if( iconv (ipr).gt.0 ) then
  isym  = 2
endif

! Matrix block size
ibsize = 1
iesize = 1

if (iresol(ipr).eq.-1) then
  ireslp = 0
  ipol   = 0
  if( iconv(ipr).gt.0 ) then
    ireslp = 1
    ipol   = 0
  endif
else
  ireslp = mod(iresol(ipr)+10000,1000)
  ipol   = (iresol(ipr)-ireslp)/1000
endif

isqrt = 1

!===============================================================================
! 2. Norm residual
!===============================================================================

if(irnpnw.ne.1) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      unsvom = -1.d0/volume(iel)
      trav(1,iel) = trav(1,iel)*unsvom + frcxt(1 ,iel) + dfrcxt(1 ,iel)
      trav(2,iel) = trav(2,iel)*unsvom + frcxt(2 ,iel) + dfrcxt(2 ,iel)
      trav(3,iel) = trav(3,iel)*unsvom + frcxt(3 ,iel) + dfrcxt(3 ,iel)
    enddo
  else
    if(isno2t.gt.0) then
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = crom(iel)-ro0
        trav(1,iel) = (trav(1,iel)+trava(1,iel))*unsvom + romro0*gx
        trav(2,iel) = (trav(2,iel)+trava(2,iel))*unsvom + romro0*gy
        trav(3,iel) = (trav(3,iel)+trava(3,iel))*unsvom + romro0*gz
      enddo
    else
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = crom(iel)-ro0
        trav(1,iel) = trav(1,iel)*unsvom + romro0*gx
        trav(2,iel) = trav(2,iel)*unsvom + romro0*gy
        trav(3,iel) = trav(3,iel)*unsvom + romro0*gz
      enddo
    endif
  endif
  do iel = 1, ncel
    dtsrom = dt(iel)/crom(iel)
    do isou = 1, 3
      trav(isou,iel) = vel(isou,iel) +dtsrom*trav(isou,iel)
    enddo
  enddo

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(trav)
    !==========
  endif


! ON NE RECONSTRUIT PAS POUR GAGNER DU TEMPS
!   EPSRGR N'EST DONC PAS UTILISE

  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = 1
  imligp = imligr(iu )
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  itypfl = 1
  if (idilat.eq.4) itypfl = 0

  call inimav &
  !==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   iflux  , bflux  )

  init = 1
  call divmas(init, iflux, bflux, res)

  ! --- Weakly compressible algorithm: semi analytic scheme
  if (idilat.eq.4) then
    do iel = 1, ncel
      res(iel) = res(iel)*crom(iel)
    enddo
  endif

  if (ncesmp.gt.0) then
    do ii = 1, ncesmp
      iel = icetsm(ii)
      res(iel) = res(iel) -volume(iel)*smacel(ii,ipr)
    enddo
  endif

! ---> LAGRANGIEN : COUPLAGE RETOUR

  if (iilagr.eq.2 .and. ltsmas.eq.1) then

    do iel = 1, ncel
      res(iel) = res(iel) -tslagr(iel,itsmas)
    enddo

  endif

  call prodsc(ncel,isqrt,res,res,rnormp)

  if(iwarni(ipr).ge.2) then
    write(nfecra,1300)chaine(1:16) ,rnormp
  endif
  dervar(ipp) = rnormp
  nbivar(ipp) = 0

else

  if(iwarni(ipr).ge.2) then
    write(nfecra,1300)chaine(1:16) ,rnormp
  endif
  dervar(ipp) = rnormp
  nbivar(ipp) = 0

endif

!===============================================================================
! 3. Compute a approximated pressure increment if needed
!    that is when there is buoyancy terms (gravity and variable density)
!    with a free outlet.
!===============================================================================

! Standard initialization
do ifac = 1, nfac
  iflux(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  coefa_dp(ifac) = 0.d0
  coefaf_dp(ifac) = 0.d0
  coefb_dp(ifac) = coefb_p(ifac)
  coefbf_dp(ifac) = coefbf_p(ifac)
  bflux(ifac) = 0.d0
enddo

! Compute a pseudo hydrostatic pressure increment stored temporarly
! in rtp(., ipr) with Homogeneous Neumann BCs everywhere
if (iphydr.eq.1.and.icalhy.eq.1) then

  ifcsor = isostd(nfabor+1)
  if (irangp.ge.0) then
    call parcmx (ifcsor)
  endif

  ! This computation is needed only if there are outlet faces
  if (ifcsor.le.0) then
    indhyd = 0
  else

    ! Work arrays for BCs
    allocate(coefap(ndimfb), cofafp(ndimfb))
    allocate(coefbp(ndimfb), cofbfp(ndimfb))

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      if (idften(ipr).eq.1) then
        hint = dt(iel)/distb(ifac)
      ! Symmetric tensor diffusivity
      elseif (idften(ipr).eq.6) then

        visci(1,1) = tpucou(1,iel)
        visci(2,2) = tpucou(2,iel)
        visci(3,3) = tpucou(3,iel)
        visci(1,2) = tpucou(4,iel)
        visci(2,1) = tpucou(4,iel)
        visci(2,3) = tpucou(5,iel)
        visci(3,2) = tpucou(5,iel)
        visci(1,3) = tpucou(6,iel)
        visci(3,1) = tpucou(6,iel)

        ! ||Ki.S||^2
        viscis = ( visci(1,1)*surfbo(1,ifac)       &
                 + visci(1,2)*surfbo(2,ifac)       &
                 + visci(1,3)*surfbo(3,ifac))**2   &
               + ( visci(2,1)*surfbo(1,ifac)       &
                 + visci(2,2)*surfbo(2,ifac)       &
                 + visci(2,3)*surfbo(3,ifac))**2   &
               + ( visci(3,1)*surfbo(1,ifac)       &
                 + visci(3,2)*surfbo(2,ifac)       &
                 + visci(3,3)*surfbo(3,ifac))**2

        ! IF.Ki.S
        fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                )*surfbo(1,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                )*surfbo(2,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                )*surfbo(3,ifac)

        distfi = distb(ifac)

        ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
        ! NB: eps =1.d-1 must be consistent with vitens.f90
        fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

        hint = viscis/surfbn(ifac)/fikis

      endif

      ! LOCAL Neumann Boundary Conditions on the hydrostatic pressure
      !--------------------------------------------------------------

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )

    enddo


    ! External forces containing bouyancy force ONLY
    do iel = 1, ncel
      dronm1 = (croma(iel)-ro0)
      drom   = (crom(iel)-ro0)
      frchy(1,iel)  = dronm1*gx
      frchy(2,iel)  = dronm1*gy
      frchy(3,iel)  = dronm1*gz
      dfrchy(1,iel) = drom  *gx - frchy(1,iel)
      dfrchy(2,iel) = drom  *gy - frchy(2,iel)
      dfrchy(3,iel) = drom  *gz - frchy(3,iel)
    enddo

    ! Parallelism and periodicity treatment
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(frchy)
      call synvin(dfrchy)
    endif

    call calhyd &
    !==========
    ( indhyd ,                                &
      frchy  , dfrchy ,                       &
      rtp(1,ipr)   , iflux , bflux ,          &
      coefap , coefbp ,                       &
      cofafp , cofbfp ,                       &
      viscf  , viscb  ,                       &
      dam    , xam    ,                       &
      drtp   , rhs    )

    ! Free memory
    deallocate(coefap, cofafp)
    deallocate(coefbp, cofbfp)

  endif

else

  indhyd = 0

endif

! Compute the BCs for the pressure increment
! (first we set the BCs of a standard pressure increment,
!  that are (A = 0, B_dp = B_p) for the gradient BCs
!  Then the A_dp is set thank to the pre-computed hydrostatic pressure
!  so that the pressure increment will be 0 on the reference outlet face.
if (iphydr.eq.1.or.iifren.eq.1) then


  if (indhyd.eq.1) then
    ifac0 = isostd(nfabor+1)
    if (ifac0.le.0) then
      phydr0 = 0.d0
    else
      iel0 = ifabor(ifac0)
      phydr0 = rtp(iel0,ipr)                                 &
           +(cdgfbo(1,ifac0)-xyzcen(1,iel0))*dfrcxt(1 ,iel0) &
           +(cdgfbo(2,ifac0)-xyzcen(2,iel0))*dfrcxt(2 ,iel0) &
           +(cdgfbo(3,ifac0)-xyzcen(3,iel0))*dfrcxt(3 ,iel0)
    endif

    if (irangp.ge.0) then
      call parsom (phydr0)
    endif
  endif

  ! If hydrostatic pressure increment or free entrance Inlet
  if (indhyd.eq.1.or.iifren.eq.1) then

    do ifac = 1, nfabor
      if (isostd(ifac).eq.1) then
        iel=ifabor(ifac)

        if (indhyd.eq.1) then
          coefa_dp(ifac) =  rtp(iel, ipr)                                 &
                       + (cdgfbo(1,ifac)-xyzcen(1,iel))*dfrcxt(1 ,iel)  &
                       + (cdgfbo(2,ifac)-xyzcen(2,iel))*dfrcxt(2 ,iel)  &
                       + (cdgfbo(3,ifac)-xyzcen(3,iel))*dfrcxt(3 ,iel)  &
                       -  phydr0
        endif

        ! Diffusive flux BCs
        if (idften(ipr).eq.1) then
          hint = dt(iel)/distb(ifac)

        ! Symmetric tensor diffusivity
        elseif (idften(ipr).eq.6) then

          visci(1,1) = tpucou(1,iel)
          visci(2,2) = tpucou(2,iel)
          visci(3,3) = tpucou(3,iel)
          visci(1,2) = tpucou(4,iel)
          visci(2,1) = tpucou(4,iel)
          visci(2,3) = tpucou(5,iel)
          visci(3,2) = tpucou(5,iel)
          visci(1,3) = tpucou(6,iel)
          visci(3,1) = tpucou(6,iel)

          ! ||Ki.S||^2
          viscis = ( visci(1,1)*surfbo(1,ifac)       &
                   + visci(1,2)*surfbo(2,ifac)       &
                   + visci(1,3)*surfbo(3,ifac))**2   &
                 + ( visci(2,1)*surfbo(1,ifac)       &
                   + visci(2,2)*surfbo(2,ifac)       &
                   + visci(2,3)*surfbo(3,ifac))**2   &
                 + ( visci(3,1)*surfbo(1,ifac)       &
                   + visci(3,2)*surfbo(2,ifac)       &
                   + visci(3,3)*surfbo(3,ifac))**2

          ! IF.Ki.S
          fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                  )*surfbo(1,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                  )*surfbo(2,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                  )*surfbo(3,ifac)

          distfi = distb(ifac)

          ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
          ! NB: eps =1.d-1 must be consistent with vitens.f90
          fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

          hint = viscis/surfbn(ifac)/fikis

        endif

        ! Free entrance boundary face (Bernoulli condition to link the pressure
        ! increment and the predicted velocity)
        if (itypfb(ifac).eq.ifrent) then

          ! Boundary mass flux of the predicted velocity
          bpmasf = vel(1, iel)*surfbo(1, ifac)    &
                 + vel(2, iel)*surfbo(2, ifac)    &
                 + vel(3, iel)*surfbo(3, ifac)

          ! Ingoing mass Flux, Bernoulli relation ship is used
          if (bpmasf.le.0.d0) then

            ! Head loss of the fluid outside the domain, between infinity and
            ! the entrance
            kpdc = b_head_loss(ifac)
            rho = brom(ifac)
            cfl = -(bmasfl(ifac)/surfbn(ifac)*dt(iel))    &
                / (2.d0*rho*distb(ifac))*(1.d0 + kpdc)

            pimp = - rtp(iel, ipr)                                              &
                 - 0.5d0*(1.d0 + kpdc)*bmasfl(ifac)*bpmasf/surfbn(ifac)**2

            call set_convective_outlet_scalar &
                 !==================
               ( coefa_dp(ifac), coefaf_dp(ifac),             &
                 coefb_dp(ifac), coefbf_dp(ifac),             &
                 pimp        , cfl         , hint )

          else
            coefaf_dp(ifac) = - hint*coefa_dp(ifac)
          endif

        else
          coefaf_dp(ifac) = - hint*coefa_dp(ifac)
        endif

      endif
    enddo
  endif

endif

!===============================================================================
! 4. Building of the linear system to solve
!===============================================================================

! ---> TERME INSTATIONNAIRE

do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! ---> Face diffusivity
if (idiff(ipr).ge.1) then

  ! Scalar diffusivity
  if (idften(ipr).eq.1) then
    call viscfa &
    !==========
 ( imvisf ,                                                       &
   dt     ,                                                       &
   viscf  , viscb  )

  ! Tensor diffusivity
  else if (idften(ipr).eq.6) then

    ! Allocate temporary arrays
    allocate(weighf(2,nfac))
    allocate(weighb(ndimfb))

    iwarnp = iwarni(ipr)

    call vitens &
    !==========
   ( tpucou , iwarnp ,             &
     weighf , weighb ,             &
     viscf  , viscb  )

  endif
else
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo
endif

iconvp = iconv (ipr)
idiffp = idiff (ipr)
ndircp = ndircl(ipr)

thetap = 1.d0
imucpp = 0

call matrix &
!==========
 ( iconvp , idiffp , ndircp , isym   ,                            &
   thetap , imucpp ,                                              &
   coefb_dp , coefbf_dp     , rovsdt ,                            &
   imasfl , bmasfl , viscf  , viscb  ,                            &
   rvoid  , dam    , xam    )

! Strengthen the diagonal
if (idilat.eq.3) then
  do iel = 1, ncel
    dam(iel) = dam(iel) + epsdp*volume(iel)/dt(iel)
  enddo
endif

! Free memory
deallocate(iflux, bflux)

!===============================================================================
! 5. Mass flux initialization
!===============================================================================

! --- Flux de masse predit et premiere composante Rhie et Chow

! Allocate a work array for the gradient calculation
allocate(gradp(3,ncelet))

iccocg = 1
iprev  = 1
inc    = 1

call field_gradient_potential(ivarfl(ipr), iprev, imrgra, inc,    &
                              iccocg, iphydr,                     &
                              frcxt, gradp)

do iel = 1, ncelet
  do isou = 1, 3
    trav(isou,iel) = gradp(isou,iel)
  enddo
enddo

if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(1,iel) = trav(1,iel) - frcxt(1 ,iel)
    trav(2,iel) = trav(2,iel) - frcxt(2 ,iel)
    trav(3,iel) = trav(3,iel) - frcxt(3 ,iel)
  enddo
endif

! --- Weakly compressible algorithm: semi analytic scheme
!     The RHS contains rho div(u*) and not div(rho u*)
!     so this term will be add afterwards
if (idilat.eq.4) then
  if (idften(ipr).eq.1) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/crom(iel))
      do isou = 1, 3
        trav(isou,iel) = ardtsr*trav(isou,iel)
      enddo
    enddo
  else if (idften(ipr).eq.6) then
    do iel = 1, ncel
      arsr  = arak/crom(iel)

      trav(1,iel) = arsr*(                                 &
                           tpucou(1,iel)*(trav(1,iel))     &
                         + tpucou(4,iel)*(trav(2,iel))     &
                         + tpucou(6,iel)*(trav(3,iel))     &
                         )
      trav(2,iel) = arsr*(                                 &
                           tpucou(4,iel)*(trav(1,iel))     &
                         + tpucou(2,iel)*(trav(2,iel))     &
                         + tpucou(5,iel)*(trav(3,iel))     &
                         )
      trav(3,iel) = arsr*(                                 &
                           tpucou(6,iel)*(trav(1,iel))     &
                         + tpucou(5,iel)*(trav(2,iel))     &
                         + tpucou(3,iel)*(trav(3,iel))     &
                         )

    enddo
  endif

! Standard algorithm
else
  if (idften(ipr).eq.1) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/crom(iel))
      do isou = 1, 3
        trav(isou,iel) = vel(isou,iel) + ardtsr*trav(isou,iel)
      enddo
    enddo
  else if (idften(ipr).eq.6) then
    do iel = 1, ncel
      arsr  = arak/crom(iel)

      trav(1,iel) = vel(1,iel)                             &
                  + arsr*(                                 &
                           tpucou(1,iel)*(trav(1,iel))     &
                         + tpucou(4,iel)*(trav(2,iel))     &
                         + tpucou(6,iel)*(trav(3,iel))     &
                         )
      trav(2,iel) = vel(2,iel)                             &
                  + arsr*(                                 &
                           tpucou(4,iel)*(trav(1,iel))     &
                         + tpucou(2,iel)*(trav(2,iel))     &
                         + tpucou(5,iel)*(trav(3,iel))     &
                         )
      trav(3,iel) = vel(3,iel)                             &
                  + arsr*(                                 &
                           tpucou(6,iel)*(trav(1,iel))     &
                         + tpucou(5,iel)*(trav(2,iel))     &
                         + tpucou(3,iel)*(trav(3,iel))     &
                         )

    enddo
  endif
endif

! ---> Traitement du parallelisme et de la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(trav)
endif

init   = 1
inc    = 1
! BCs will be taken into account after in idilat=4
if (idilat.eq.4) inc = 0
iflmb0 = 1
if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
nswrgp = nswrgr(iu )
imligp = imligr(iu )
iwarnp = iwarni(ipr)
epsrgp = epsrgr(iu )
climgp = climgr(iu )
itypfl = 1

call inimav &
!==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   imasfl , bmasfl )


! --- Projection aux faces des forces exterieures

if (iphydr.eq.1) then
  init   = 0
  inc    = 0
  iccocg = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  ircflp = ircflu(ipr)

  ! Scalar diffusivity
  if (idften(ipr).eq.1) then
    call projts &
    !==========
 ( init   , nswrgp ,                                              &
   dfrcxt ,                                                       &
   coefbf_p ,                                                     &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     )

  ! Tensor diffusivity
  else if (idften(ipr).eq.6) then

    call projtv &
    !==========
  ( init   , nswrgp , ircflp ,                                     &
    dfrcxt ,                                                       &
    coefbf_p ,                                                     &
    viscf  , viscb  ,                                              &
    tpucou ,                                                       &
    weighf ,                                                       &
    imasfl , bmasfl )

  endif
endif

init   = 0
inc    = 1
iccocg = 1

if (arak.gt.0.d0) then

! --- Prise en compte de Arak : la viscosite face est multipliee
!       Le pas de temps aussi. On retablit plus bas.
  do ifac = 1, nfac
    viscf(ifac) = arak*viscf(ifac)
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = arak*viscb(ifac)
  enddo

  ! On annule la viscosite facette pour les faces couplees pour ne pas modifier
  ! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
  ! de pression et le filtre sont annules.
  if (nbrcpl.gt.0) then
    do ifac = 1, nfabor
      if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
        viscb(ifac) = 0.d0
      endif
    enddo
  endif

  ! Scalar diffusivity
  !-------------------
  if (idften(ipr).eq.1) then
    do iel = 1, ncel
      dt(iel) = arak*dt(iel)
    enddo

    nswrgp = nswrgr(ipr )
    imligp = imligr(ipr )
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr )
    climgp = climgr(ipr )
    extrap = extrag(ipr )

    call itrmas &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   frcxt  ,                                                       &
   rtpa(1,ipr)  ,                                                 &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                      &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   imasfl , bmasfl )

    ! Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      iccocg = 1
      nswrgp = nswrgr(ipr)
      imligp = imligr(ipr)
      iwarnp = iwarni(ipr)
      epsrgp = epsrgr(ipr)
      climgp = climgr(ipr)

      ! A 0 boundary coefficient coefbf_dp is passed to projts
      ! to cancel boundary terms
      allocate(cofbfp(ndimfb))
      do ifac = 1,nfabor
        cofbfp(ifac) = 0.d0
      enddo

      call projts &
      !==========
 ( init   , nswrgp ,                                              &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     )

      deallocate(cofbfp)

    endif

    ! --- Correction du pas de temps
    unsara = 1.d0/arak
    do iel = 1, ncel
      dt(iel) = dt(iel)*unsara
    enddo

  ! Tensor diffusivity
  !-------------------
  else if (idften(ipr).eq.6) then

    do iel = 1, ncel
      tpucou(1,iel) = arak*tpucou(1,iel)
      tpucou(2,iel) = arak*tpucou(2,iel)
      tpucou(3,iel) = arak*tpucou(3,iel)
      tpucou(4,iel) = arak*tpucou(4,iel)
      tpucou(5,iel) = arak*tpucou(5,iel)
      tpucou(6,iel) = arak*tpucou(6,iel)
    enddo

    nswrgp = nswrgr(ipr)
    imligp = imligr(ipr)
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr)
    climgp = climgr(ipr)
    extrap = extrag(ipr)
    ircflp = ircflu(ipr)

    call itrmav &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   frcxt  ,                                                       &
   rtpa(1,ipr)  ,                                                 &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                      &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   weighf , weighb ,                                              &
   imasfl , bmasfl )

   ! Projection du terme source pour oter la partie hydrostat de la pression
   if (iphydr.eq.1) then
     init   = 0
     inc    = 0
     nswrgp = nswrgr(ipr)
     imligp = imligr(ipr)
     iwarnp = iwarni(ipr)
     epsrgp = epsrgr(ipr)
     climgp = climgr(ipr)

     ! A 0 boundary coefficient coefbf_dp is passed to projtv
     ! to cancel boundary terms
     allocate(cofbfp(ndimfb))
     do ifac = 1, nfabor
       cofbfp(ifac) = 0.d0
     enddo

     call projtv &
     !==========
   ( init   , nswrgp , ircflp ,                                     &
     frcxt  ,                                                       &
     cofbfp ,                                                       &
     viscf  , viscb  ,                                              &
     tpucou ,                                                       &
     weighf ,                                                       &
     imasfl, bmasfl )

     deallocate(cofbfp)

   endif

    ! --- Correction du pas de temps
    unsara = 1.d0/arak
    do iel = 1, ncel
      tpucou(1,iel) = unsara*tpucou(1,iel)
      tpucou(2,iel) = unsara*tpucou(2,iel)
      tpucou(3,iel) = unsara*tpucou(3,iel)
      tpucou(4,iel) = unsara*tpucou(4,iel)
      tpucou(5,iel) = unsara*tpucou(5,iel)
      tpucou(6,iel) = unsara*tpucou(6,iel)
    enddo

  endif

  ! --- Correction de la viscosite aux faces
  do ifac = 1, nfac
    viscf(ifac) = viscf(ifac)*unsara
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = viscb(ifac)*unsara
  enddo

  ! If Saturne/Saturne coupling, re-set the boundary face viscosity to
  ! the non-zero value
  if (nbrcpl.gt.0) then
    if (idiff(ipr).ge.1) then
      if (idften(ipr).eq.1) then
        call viscfa &
        !==========
     ( imvisf ,                                                       &
       dt     ,                                                       &
       viscf  , viscb  )

      ! Tensor diffusivity
      else if (idften(ipr).eq.6) then

        iwarnp = iwarni(ipr)

        call vitens &
        !==========
       ( tpucou , iwarnp ,             &
         weighf , weighb ,             &
         viscf  , viscb  )

      endif
    else
      do ifac = 1, nfac
        viscf(ifac) = 0.d0
      enddo
      do ifac = 1, nfabor
        viscb(ifac) = 0.d0
      enddo
    endif

  endif

endif

!===============================================================================
! 6. Preparation of the Algebraic Multigrid
!===============================================================================

if (imgr(ipr).gt.0) then

  ! --- Building of the mesh hierarchy

  iwarnp = iwarni(ipr)
  nagmax = nagmx0(ipr)
  npstmg = ncpmgr(ipr)

  call clmlga &
  !==========
 ( chaine(1:16) ,    lchain ,                                     &
   isym   , ibsize , iesize , nagmax , npstmg , iwarnp ,          &
   ngrmax , ncegrm ,                                              &
   rlxp1  ,                                                       &
   dam    , xam    )

endif

!===============================================================================
! 7. Solving (Loop over the non-orthogonalities)
!===============================================================================

! --- Number of sweeps
nswmpr = nswrsm(ipr)

! --- Variables are set to 0
!       rtp(.,IPR) is the increment of the pressure
!       drtp       is the increment of the increment between sweeps
!       divu       is the initial divergence of the predicted mass flux
do iel = 1, ncel
  rtp(iel,ipr) = 0.d0
  drtp(iel) = 0.d0
  presa(iel) = 0.d0
enddo

relaxp = relaxv(ipr)

! --- Initial divergence
init = 1

call divmas(init, imasfl , bmasfl , divu)

! --- Weakly compressible algorithm: semi analytic scheme
!     1. The RHS contains rho div(u*) and not div(rho u*)
!     2. Add dilatation source term to rhs
!     3. The mass flux is completed by rho u* . S

if (idilat.eq.4) then

  allocate(velflx(nfac), velflb(ndimfb))

  ! 1. The RHS contains rho div(u*) and not div(rho u*)
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)

  itypfl = 0

  call inimav &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   velflx , velflb )

  call divmas(init, velflx , velflb , res)

  do iel = 1, ncel
    divu(iel) = divu(iel) + res(iel)*crom(iel)
  enddo

  ! 2. Add the dilatation source term D(rho)/Dt
  do iel = 1, ncel
    divu(iel) = divu(iel) + propce(iel,ipproc(iustdy(itsrho)))
  enddo

  ! 3. The mass flux is completed by rho u* . S
  init   = 0
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)

  itypfl = 1

  call inimav &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   imasfl , bmasfl )

endif

! --- Termes sources de masse
if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    divu(iel) = divu(iel) -volume(iel)*smacel(ii,ipr)
  enddo
endif

! --- Source term associated to the mass aggregation
if (idilat.eq.2.or.idilat.eq.3) then
  do iel = 1, ncel
    drom = crom(iel) - croma(iel)
    divu(iel) = divu(iel) + drom*volume(iel)/dt(iel)
  enddo
endif

! ---> Termes sources Lagrangien
if (iilagr.eq.2 .and. ltsmas.eq.1) then
  do iel = 1, ncel
    divu(iel) = divu(iel) -tslagr(iel,itsmas)
  enddo
endif

! --- Initial right hand side
do iel = 1, ncel
  rhs(iel) = - divu(iel)
enddo

! --- Add eps*pressure*volume/dt in the right hand side
!     to strengthen the diagonal for the low-Mach algo.
if (idilat.eq.3) then
  do iel = 1, ncel
    rhs(iel) = rhs(iel) - epsdp*volume(iel)/dt(iel)*rtp(iel,ipr)
  enddo
endif

! --- Right hand side residual
call prodsc(ncel,isqrt,rhs,rhs,residu)

rnsmbr(ipp) = residu

isweep = 1

! Writing
if (iwarni(ipr).ge.2) then
  write(nfecra,1400)chaine(1:16),isweep,residu, relaxp
endif

! Dynamic relaxation initialization
!----------------------------------
if (iswdyp.ge.1) then

  do iel = 1, ncelet
    adxkm1(iel) = 0.d0
    adxk(iel) = 0.d0
  enddo

  ! ||A.dx^0||^2 = 0
  nadxk = 0.d0

  rnorm2 = rnormp**2

  iccocg = 1
  init = 1
  inc  = 0
  if (iphydr.eq.1.or.iifren.eq.1) inc = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  extrap = extrag(ipr)
  ircflp = ircflu(ipr)

  if (idften(ipr).eq.1) then

    call itrgrp &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   rhs0   )

  else if (idften(ipr).eq.6) then

    call itrgrv &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   weighf , weighb ,                                              &
   rhs0   )

  endif

endif

! Reconstruction loop (beginning)
!--------------------------------

do while (isweep.le.nswmpr.and.residu.gt.epsrsm(ipr)*rnormp)

  ! Solving on the increment drtp
  !------------------------------
  if (iswdyp.eq.0) then
    do iel = 1, ncel
      drtp(iel) = 0.d0
    enddo
  else
    do iel = 1, ncel
      dpvarm1(iel) = drtp(iel)
      drtp(iel) = 0.d0
    enddo
  endif

  nitmap = nitmax(ipr)
  imgrp  = imgr  (ipr)
  ncymap = ncymax(ipr)
  nitmgp = nitmgf(ipr)
  iwarnp = iwarni(ipr)
  epsilp = epsilo(ipr)

  ! The pressure is a scalar => no problem for the periodicity of rotation
  ! (iinvpe=1)
  iinvpe = 1

  ! Solver reisudal
  ressol = residu

  call invers &
  !==========
 ( chaine(1:16)    , isym   , ibsize , iesize ,                   &
   ipol   , ireslp , nitmap , imgrp  ,                            &
   ncymap , nitmgp ,                                              &
   iwarnp , niterf , icycle , iinvpe ,                            &
   epsilp , rnormp , ressol ,                                     &
   dam    , xam    , rhs    , drtp   )

  ! Dynamic relaxation of the system
  !---------------------------------
  if (iswdyp.ge.1) then

    ! Computation of the variable ralaxation coefficient

    !$omp parallel do
    do iel = 1, ncelet
      adxkm1(iel) = adxk(iel)
      adxk(iel) = - rhs0(iel)
    enddo

    iccocg = 1
    init = 0
    inc  = 0
    if (iphydr.eq.1.or.iifren.eq.1) inc = 1
    nswrgp = nswrgr(ipr)
    imligp = imligr(ipr)
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr)
    climgp = climgr(ipr)
    extrap = extrag(ipr)

    if (idften(ipr).eq.1) then

      call itrgrp &
      !==========
   ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
     iwarnp ,                                                       &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt ,                                                       &
     drtp   ,                                                       &
     coefa_dp  , coefb_dp  ,                                        &
     coefaf_dp , coefbf_dp ,                                        &
     viscf  , viscb  ,                                              &
     dt     , dt     , dt     ,                                     &
     adxk   )

    else if (idften(ipr).eq.6) then

      call itrgrv &
      !==========
   ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , iwarnp ,                                              &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt ,                                                       &
     drtp   ,                                                       &
     coefa_dp  , coefb_dp  ,                                        &
     coefaf_dp , coefbf_dp ,                                        &
     viscf  , viscb  ,                                              &
     tpucou ,                                                       &
     weighf , weighb ,                                              &
     adxk   )

    endif

    do iel = 1, ncel
      adxk(iel) = - adxk(iel)
    enddo

    insqrt = 0

    ! ||E.dx^(k-1)-E.0||^2
    nadxkm1 = nadxk

    ! ||E.dx^k-E.0||^2
    call prodsc(ncel , insqrt , adxk , adxk , nadxk)

    ! < E.dx^k-E.0; r^k >
    call prodsc(ncel , insqrt , rhs , adxk , paxkrk)

    ! Relaxation with respect to dx^k and dx^(k-1)
    if (iswdyp.ge.2) then

      ! < E.dx^(k-1)-E.0; r^k >
      call prodsc(ncel , insqrt , rhs , adxkm1 , paxm1rk)

      ! < E.dx^(k-1)-E.0; E.dx^k -E.0 >
      call prodsc(ncel , insqrt , adxk, adxkm1 , paxm1ax)

      if (nadxkm1.gt.1.d-30*rnorm2.and.                    &
         (nadxk*nadxkm1-paxm1ax**2).gt.1.d-30*rnorm2) then
        beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-paxm1ax**2)
      else
        beta = 0.d0
      endif

    else
      beta = 0.d0
      paxm1ax = 1.d0
      paxm1rk = 0.d0
      paxm1ax = 0.d0
    endif

    ! The first sweep is not relaxed
    if (isweep.eq.1) then
      alph = 1.d0
      beta = 0.d0
    elseif (isweep.eq.2) then
      beta = 0.d0
      alph = -paxkrk/max(nadxk, 1.d-30*rnorm2)
    else
      alph = -(paxkrk + beta*paxm1ax)/max(nadxk, 1.d-30*rnorm2)
    endif

    ! Writing
    if (iwarnp.ge.2) then
      write(nfecra,1200) chaine(1:16), isweep, alph, beta, &
                         paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax
    endif

  endif

  ! Update the increment of pressure
  !---------------------------------

  if (iswdyp.eq.0) then
    if (idtvar.ge.0.and.isweep.le.nswmpr.and.residu.gt.epsrsm(ipr)*rnormp) then
      do iel = 1, ncel
        presa(iel) = rtp(iel,ipr)
        rtp(iel,ipr) = presa(iel) + relaxv(ipr)*drtp(iel)
      enddo
    ! If it is the last sweep, update with the total increment
    else
      do iel = 1, ncel
        presa(iel) = rtp(iel,ipr)
        rtp(iel,ipr) = presa(iel) + drtp(iel)
      enddo
    endif
  elseif (iswdyp.eq.1) then
     do iel = 1, ncel
      presa(iel) = rtp(iel,ipr)
      rtp(iel,ipr) = presa(iel) + alph*drtp(iel)
    enddo
  elseif (iswdyp.ge.2) then
    do iel = 1, ncel
      presa(iel) = rtp(iel,ipr)
      rtp(iel,ipr) = presa(iel) + alph*drtp(iel) + beta*dpvarm1(iel)
    enddo
  endif

  ! --- Update the right hand side and update the residual
  !      rhs^{k+1} = - div(rho u^n) - D(dt, delta delta p^{k+1})
  !-------------------------------------------------------------

  iccocg = 1
  init = 1
  inc  = 0
  if (iphydr.eq.1.or.iifren.eq.1) inc = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  extrap = extrag(ipr)

  if (idften(ipr).eq.1) then

    call itrgrp &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   rtp(1,ipr)      ,                                              &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   rhs    )

  else if (idften(ipr).eq.6) then

    call itrgrv &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   rtp(1,ipr)      ,                                              &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   weighf , weighb ,                                              &
   rhs    )

  endif

  do iel = 1, ncel
    rhs(iel) = - divu(iel) - rhs(iel)
  enddo

  ! --- Add eps*pressure*volume/dt in the right hand side
  !     to strengthen the diagonal for the low-Mach algo.
  if (idilat.eq.3) then
    do iel = 1, ncel
      rhs(iel) = rhs(iel) - epsdp*volume(iel)/dt(iel)*rtp(iel,ipr)
    enddo
  endif

  ! --- Convergence test
  call prodsc(ncel,isqrt,rhs,rhs,residu)

  ! Writing
  if (iwarni(ipr).ge.2) then
    write(nfecra,1400)chaine(1:16),isweep,residu, relaxp
  endif

  ! Writing
  nbivar(ipp) = nbivar(ipp) + niterf
  if (abs(rnormp).gt.0.d0) then
    resvar(ipp) = residu/rnormp
  else
    resvar(ipp) = 0.d0
  endif

  ! Writing
  if (iwarnp.ge.3) then
    write(nfecra,1500) chaine(1:16), isweep, residu, rnormp, niterf
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

! Writing
if(iwarni(ipr).ge.2) then
  if(isweep.gt.nswmpr) then
     write(nfecra,1600) chaine(1:16),nswmpr
  endif
endif

! --- Compute the indicator, taken the volume into account (L2 norm)
!     or not
if(iescal(iesder).gt.0) then
  iesdep = ipproc(iestim(iesder))
  do iel = 1, ncel
    propce(iel,iesdep) = abs(rhs(iel))/volume(iel)
  enddo
  if(iescal(iesder).eq.2) then
    do iel = 1, ncel
      propce(iel,iesdep) = propce(iel,iesdep)*sqrt(volume(iel))
    enddo
  endif
endif

! Update the mass flux
!---------------------

! On annule la viscosite facette pour les faces couplees pour ne pas modifier
! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
! de pression et le filtre sont annules.
if (nbrcpl.ge.0) then
  do ifac = 1, nfabor
    if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
      viscb(ifac) = 0.d0
    endif
  enddo
endif

iccocg = 1
init = 0
inc  = 0
! In case of hydrostatic pressure, inc is set to 1 to take explicit
! boundary conditions on the pressure (coefa)
if (iphydr.eq.1.or.iifren.eq.1) inc = 1
nswrgp = nswrgr(ipr)
imligp = imligr(ipr)
iwarnp = iwarni(ipr)
epsrgp = epsrgr(ipr)
climgp = climgr(ipr)
extrap = extrag(ipr)

if (idften(ipr).eq.1) then
  call itrmas &
  !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   presa  ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0

  call itrmas &
  !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   imasfl , bmasfl )

else if (idften(ipr).eq.6) then

  call itrmav &
  !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   presa  ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   weighf , weighb ,                                              &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0
  ircflp = 0

  call itrmav &
  !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa_dp  , coefb_dp  ,                                        &
   coefaf_dp , coefbf_dp ,                                        &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   weighf , weighb ,                                              &
   imasfl , bmasfl )

endif

!===============================================================================
! 8. Suppression of the mesh hierarchy
!===============================================================================

if (imgr(ipr).gt.0) then
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

!===============================================================================
! 9. Weakly compressible algorithm: semi analytic scheme
!    2nd step solving a convection diffusion equation
!===============================================================================

if (idilat.eq.4) then

  deallocate(gradp)
  allocate(gradp(ncelet,3))

  ! Allocate temporary arrays
  allocate(dpvar(ncelet))
  allocate(coefar(3,ndimfb), cofafr(3,ndimfb))
  allocate(coefbr(3,3,ndimfb), cofbfr(3,3,ndimfb))

  ! --- Convective flux: dt/rho grad(rho)
  inc = 1
  iccocg = 1
  ivar   = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  ! Dirichlet Boundary Condition on rho
  !------------------------------------

  do ifac = 1, nfabor
    coefa_dp(ifac) = brom(ifac)
    coefb_dp(ifac) = 0.d0
  enddo

  call grdcel &
  !==========
  (ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   crom   ,                                                       &
   coefa_dp        , coefb_dp        ,                            &
   gradp  )

  ! --- dt/rho * grad rho
  do iel = 1, ncel
    do isou = 1, 3
      trav(isou,iel) = gradp(iel,isou) * dt(iel) / crom(iel)
    enddo
  enddo

  ! --- (dt/rho * grad rho) . S

  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  itypfl = 0

  ! --- Viscosity
  call viscfa (imvisf, dt, viscf, viscb)

  ! --- Boundary Conditions for the convective flux
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! Neumann Boundary Conditions
    !----------------------------

    qimpv(1) = 0.d0
    qimpv(2) = 0.d0
    qimpv(3) = 0.d0

    if (idften(ipr).eq.1) then
      hint = dt(iel)/distb(ifac)

    ! Symmetric tensor diffusivity
    elseif (idften(ipr).eq.6) then

      visci(1,1) = tpucou(1,iel)
      visci(2,2) = tpucou(2,iel)
      visci(3,3) = tpucou(3,iel)
      visci(1,2) = tpucou(4,iel)
      visci(2,1) = tpucou(4,iel)
      visci(2,3) = tpucou(5,iel)
      visci(3,2) = tpucou(5,iel)
      visci(1,3) = tpucou(6,iel)
      visci(3,1) = tpucou(6,iel)

      ! ||Ki.S||^2
      viscis = ( visci(1,1)*surfbo(1,ifac)       &
               + visci(1,2)*surfbo(2,ifac)       &
               + visci(1,3)*surfbo(3,ifac))**2   &
             + ( visci(2,1)*surfbo(1,ifac)       &
               + visci(2,2)*surfbo(2,ifac)       &
               + visci(2,3)*surfbo(3,ifac))**2   &
             + ( visci(3,1)*surfbo(1,ifac)       &
               + visci(3,2)*surfbo(2,ifac)       &
               + visci(3,3)*surfbo(3,ifac))**2

      ! IF.Ki.S
      fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
              )*surfbo(1,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
              )*surfbo(2,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
              )*surfbo(3,ifac)

      distfi = distb(ifac)

      ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
      ! NB: eps =1.d-1 must be consistent with vitens.f90
      fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

      hint = viscis/surfbn(ifac)/fikis

    endif

    call set_neumann_vector &
         !=================
       ( coefar(1,ifac)  , cofafr(1,ifac)  ,             &
         coefbr(1,1,ifac), cofbfr(1,1,ifac),             &
         qimpv           , hint             )

  enddo

  call inimav &
  !==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   trav   ,                                                       &
   coefar , coefbr ,                                              &
   velflx , velflb )

  ! --- Boundary condition for the pressure increment
  do ifac = 1, nfabor
   coefa_dp(ifac) = 0.d0
   coefaf_dp(ifac) = 0.d0
  enddo

  ! --- Convective source term
  do iel = 1, ncel
    rhs(iel) = 0.d0
  enddo

  ivar   = ipr
  iconvp = 1
  idiffp = 0
  nswrsp = 1
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  inc    = 1
  iccocg = 1
  ipp    = ipprtp(ivar)
  iwarnp = iwarni(ivar)
  imucpp = 0
  idftnp = idften(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = 1.d0
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,ipr)      , rtp(1,ipr)      ,                            &
   coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
   velflx , velflb , viscf  , viscb  , rvoid  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rhs   )

  ! --- Initialization of the variable to solve
  do iel = 1, ncel
    rovsdt(iel) = 340.d0/dt(iel) * volume(iel)
    drtp(iel)   = 0.d0
    dpvar(iel)  = 0.d0
    rhs(iel)    = - rhs(iel)
  enddo

  ! --- Solve the convection diffusion equation

  idiffp = 1
  ireslp = 1
  ipol   = 0
  ! To reinforce the diagonal
  ndircp = 0
  nitmap = nitmax(ivar)
  nswrsp = nswrsm(ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  iescap = 0
  imucpp = 0
  idftnp = idften(ivar)
  iswdyp = iswdyn(ivar)
  imgrp  = 0
  ncymxp = ncymax(ivar)
  nitmfp = nitmgf(ivar)
  ipp    = ipprtp(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsilp = epsilo(ivar)
  epsrsp = epsrsm(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)
  ! all boundary convective flux with upwind
  icvflb = 0

  ! --- Solve the convection diffusion equation

  call codits &
  !==========
   ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
     imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
     ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
     imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
     blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
     relaxp , thetap ,                                              &
     drtp   , drtp   ,                                              &
     coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
     velflx , velflb ,                                              &
     viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
     weighf , weighb ,                                              &
     icvflb , ivoid  ,                                              &
     rovsdt , rhs    , drtp   , dpvar  ,                            &
     rvoid  , rvoid  )

  ! --- Update the increment of Pressure

  do iel = 1, ncel
    rtp(iel,ipr) = rtp(iel,ipr) + drtp(iel)
    ! Remove the last increment
    drtp(iel) = drtp(iel) - dpvar(iel)
  enddo

  ! --- Update the Mass flux

  init   = 0
  inc    = 1
  iccocg = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  extrap = extrag(ipr)

  if (idften(ipr).eq.1) then
    call itrmas &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   imasfl , bmasfl )

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrgp = 0
    inc = 0

    call itrmas &
    !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   dpvar  ,                                                       &
   coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   imasfl , bmasfl )

  else if (idften(ipr).eq.6) then

    call itrmav &
    !==========
   ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , iwarnp ,                                              &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt ,                                                       &
     drtp   ,                                                       &
     coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
     viscf  , viscb  ,                                              &
     tpucou ,                                                       &
     weighf , weighb ,                                              &
     imasfl , bmasfl )

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrgp = 0
    inc = 0

    call itrmav &
    !==========
   ( init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , iwarnp ,                                              &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt ,                                                       &
     dpvar  ,                                                       &
     coefa_dp , coefb_p       , coefaf_dp       , coefbf_p ,        &
     viscf  , viscb  ,                                              &
     tpucou ,                                                       &
     weighf , weighb ,                                              &
     imasfl , bmasfl )

  endif

  ! Free memory
  deallocate(dpvar)
  deallocate(coefar, coefbr)
  deallocate(cofafr, cofbfr)
  deallocate(velflx, velflb)

endif

!===============================================================================
! 10. Update the pressure field
!===============================================================================

if (idtvar.lt.0) then
  do iel = 1, ncel
    rtp(iel,ipr) = rtpa(iel,ipr) + relaxv(ipr)*rtp(iel,ipr)
  enddo
else
  do iel = 1, ncel
    rtp(iel,ipr) = rtpa(iel,ipr) + rtp(iel,ipr)
  enddo
endif

! Free memory
deallocate(dam, xam)
deallocate(res, divu, presa)
deallocate(gradp)
deallocate(coefaf_dp, coefbf_dp)
deallocate(rhs, rovsdt)
if (allocated(weighf)) deallocate(weighf, weighb)
if (iswdyp.ge.1) deallocate(adxk, adxkm1, dpvarm1, rhs0)
if (icalhy.eq.1) deallocate(frchy, dfrchy)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : RESIDU DE NORMALISATION =', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6,  &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION'                       ,/,&
'@    ========='                                                    ,/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint'                ,/,&
'@' )

#else

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : NORMED RESIDUALS = ', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6, &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ WARNING: ', A16,' PRESSURE STEP'                              ,/,&
'@    ========'                                                     ,/,&
'@  Maximum number of iterations ',I10   ,' reached'                ,/,&
'@'                                                              )

#endif

!----
! End
!----

return

end subroutine resopv
