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

!> \file navstv.f90
!>
!> \brief Solving of NS equations for incompressible or slightly compressible
!> flows for one time step. Both convection-diffusion and continuity steps are
!> performed.  The velocity components are solved together in once.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        index of the iteration on Navier-Stokes
!> \param[in]     icvrge        indicator of convergence
!> \param[in]     itrale        number of the current ALE iteration
!> \param[in]     isostd        indicator of standar outlet
!>                               +index of the reference face
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     frcxt         external force generating the hydrostatic
!>                              pressure
!> \param[in]     prhyd         hydrostatic pressure predicted at cell centers
!> \param[in]     trava         work array for pressure velocity coupling
!> \param[in]     ximpa         work array for pressure velocity coupling
!> \param[in]     uvwk          work array for pressure velocity coupling
!_______________________________________________________________________________


subroutine navstv &
 ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
   isostd ,                                                       &
   dt     , rtp    , rtpa   , propce ,                            &
   frcxt  , prhyd  ,                                              &
   trava  , ximpa  , uvwk   )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb, nproce
use numvar
use entsor
use cstphy
use cstnum
use optcal
use pointe
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use lagran, only: iilagr
use lagdim, only: ntersl
use turbomachinery
use ptrglo
use field
use cavitation

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iterns , icvrge , itrale

integer          isostd(nfabor+1)

double precision, pointer, dimension(:)   :: dt
double precision, pointer, dimension(:,:) :: rtp, rtpa, propce
double precision, pointer, dimension(:,:) :: frcxt
double precision, pointer, dimension(:) :: prhyd
double precision, pointer, dimension(:,:) :: trava, uvwk
double precision, pointer, dimension(:,:,:) :: ximpa

! Local variables

integer          iccocg, inc, iel, iel1, iel2, ifac, imax, imaxt
integer          ii    , inod, itypfl
integer          isou, ivar, iitsm
integer          init
integer          iflmas, iflmab
integer          iflmb0
integer          nswrgp, imligp, iwarnp
integer          nbrval, iappel, iescop
integer          ndircp, icpt
integer          numcpl
integer          iflvoi, iflvob
double precision rnorm , rnormt, rnorma, rnormi, vitnor
double precision dtsrom, unsrom, surf  , rhom, rovolsdt
double precision epsrgp, climgp, extrap, xyzmax(3)
double precision thetap, xdu, xdv, xdw
double precision xxp0 , xyp0 , xzp0
double precision rhofac, dtfac, ddepx , ddepy, ddepz
double precision xnrdis, xnrtmp
double precision vitbox, vitboy, vitboz
double precision t1, t2, t3, t4
double precision visclc, visctc
double precision distbf, srfbnf, hint
double precision rnx, rny, rnz, rcodcx, rcodcy, rcodcz, rcodcn

double precision, allocatable, dimension(:,:,:), target :: viscf
double precision, allocatable, dimension(:), target :: viscb
double precision, allocatable, dimension(:,:,:), target :: wvisfi
double precision, allocatable, dimension(:), target :: wvisbi
double precision, allocatable, dimension(:) :: drtp
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:) :: esflum, esflub
double precision, allocatable, dimension(:) :: intflx, bouflx
double precision, allocatable, dimension(:) :: secvif, secvib

double precision, dimension(:,:), allocatable :: gradp
double precision, dimension(:,:), allocatable :: mshvel
double precision, dimension(:), allocatable :: coefa_dp, coefb_dp
double precision, dimension(:), allocatable :: xinvro

double precision, dimension(:,:), pointer :: grdphd
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:,:,:), pointer :: viscfi
double precision, dimension(:), pointer :: viscbi
double precision, dimension(:,:), pointer :: dttens
double precision, dimension(:,:), pointer :: dfrcxt

double precision, dimension(:,:), pointer :: coefau, cofafu, claale
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu, clbale

double precision, dimension(:), pointer :: coefa_p
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, croma, viscl, visct
double precision, dimension(:), pointer :: ivoifl, bvoifl
double precision, dimension(:), pointer :: coavoi, cobvoi
double precision, dimension(:,:), pointer :: trav

double precision, dimension(:), pointer :: cvar_pr, cvara_pr
double precision, dimension(:), pointer :: cpro_prtot

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

! Allocate temporary arrays for the velocity-pressure resolution
if (idften(iu).eq.1) then
  allocate(viscf(1, 1, nfac), viscb(ndimfb))
else if (idften(iu).eq.6) then
  allocate(viscf(3, 3, nfac), viscb(ndimfb))
endif

allocate(trav(3,ncelet))

! Allocate other arrays, depending on user options

! Array for delta p gradient boundary conditions
allocate(coefa_dp(ndimfb), coefb_dp(ndimfb))

allocate(dfrcxt(3,ncelet))
if (iphydr.eq.2) then
  allocate(grdphd(ncelet,ndim))
else
  grdphd => rvoid2
endif
if (idften(iu).eq.1) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(1,1,nfac), wvisbi(ndimfb))
    viscfi => wvisfi(:,:,1:nfac)
    viscbi => wvisbi(1:ndimfb)
  else
    viscfi => viscf(:,:,1:nfac)
    viscbi => viscb(1:ndimfb)
  endif
else if(idften(iu).eq.6) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(3,3,nfac), wvisbi(ndimfb))
    viscfi => wvisfi(1:3,1:3,1:nfac)
    viscbi => wvisbi(1:ndimfb)
  else
    viscfi => viscf(1:3,1:3,1:nfac)
    viscbi => viscb(1:ndimfb)
  endif
endif

if (ivisse.eq.1) then
  allocate(secvif(nfac),secvib(ndimfb))
endif

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_val_prev_v(ivarfl(iu), vela)

! Map some specific field arrays
if (idtten.ge.0) then
  call field_get_val_v(idtten, dttens)
else
  dttens => rvoid2
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_prev_s(ivarfl(ipr), cvara_pr)

! Allocate work arrays
allocate(w1(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
if (irnpnw.eq.1) allocate(w10(ncelet))

if (iwarni(iu).ge.1) then
  write(nfecra,1000)
endif

! Initialize variables to avoid compiler warnings

ivar = 0
iflmas = 0
imax = 0

! Memory

if (nterup.gt.1) then

  !$omp parallel do private(isou)
  do iel = 1,ncelet
    do isou = 1, 3
    !     La boucle sur NCELET est une securite au cas
    !       ou on utiliserait UVWK par erreur a ITERNS = 1
      uvwk(isou,iel) = vel(isou,iel)
    enddo
  enddo

  ! Calcul de la norme L2 de la vitesse
  if (iterns.eq.1) then
    xnrtmp = 0.d0
    !$omp parallel do reduction(+:xnrtmp)
    do iel = 1, ncel
      xnrtmp = xnrtmp +(vela(1,iel)**2        &
                      + vela(2,iel)**2        &
                      + vela(3,iel)**2)       &
                      * volume(iel)
    enddo
    xnrmu0 = xnrtmp
    if (irangp.ge.0) then
      call parsom (xnrmu0)
      !==========
    endif
    ! En cas de couplage entre deux instances de Code_Saturne, on calcule
    ! la norme totale de la vitesse
    ! Necessaire pour que l'une des instances ne stoppe pas plus tot que les autres
    ! (il faudrait quand meme verifier les options numeriques, ...)
    do numcpl = 1, nbrcpl
      call tbrcpl ( numcpl, 1, 1, xnrmu0, xnrdis )
      !==========
      xnrmu0 = xnrmu0 + xnrdis
    enddo
    xnrmu0 = sqrt(xnrmu0)
  endif

  ! On assure la periodicite ou le parallelisme de UVWK et la pression
  ! (cette derniere vaut la pression a l'iteration precedente)
  if (iterns.gt.1) then
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(uvwk(1,1))
      !==========
      call synsca(cvara_pr)
      !==========
    endif
  endif

endif

! --- Physical quantities
call field_get_val_s(iprpfl(iviscl), viscl)
call field_get_val_s(iprpfl(ivisct), visct)

! Initialize timers
t1 = 0.d0
t2 = 0.d0
t3 = 0.d0
t4 = 0.d0

!===============================================================================
! 1. Prediction of the mass flux in case of Low Mach compressible algorithm
!===============================================================================

if ((idilat.eq.2.or.idilat.eq.3).and. &
    (ntcabs.gt.1.or.isuite.gt.0)) then

  call predfl(nvar, ncetsm, icetsm, dt, smacel)
  !==========

endif

!===============================================================================
! 2. Hydrostatic pressure prediction in case of Low Mach compressible algorithm
!===============================================================================

if (iphydr.eq.2) then

  call prehyd(prhyd, grdphd)
  !==========

endif

!===============================================================================
! 3. Pressure resolution and computation of mass flux for compressible flow
!===============================================================================

if ( ippmod(icompf).ge.0 ) then

  if(iwarni(iu).ge.1) then
    write(nfecra,1080)
  endif

  call cfmspr &
  !==========
  ( nvar   , nscal  , iterns ,                                     &
    ncepdc , ncetsm , icepdc , icetsm , itypsm ,                   &
    dt     , rtp    , rtpa   , propce , vela   ,                   &
    ckupdc , smacel )

endif

!===============================================================================
! 4. Compute liquid-vapour mass transfer term for cavitating flows
!===============================================================================

if (icavit.ge.0) then

  call cavitation_compute_source_term (rtpa(:,ipr), rtpa(:,ivoidf))
  !==================================

endif

!===============================================================================
! 5. Velocity prediction step
!===============================================================================

iappel = 1

! Id of the mass flux
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)

! Pointers to the mass fluxes
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Pointers to properties
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

! Pointers to BC coefficients
call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

call predvv &
!==========
( iappel ,                                                       &
  nvar   , nscal  , iterns ,                                     &
  ncepdc , ncetsm ,                                              &
  icepdc , icetsm , itypsm ,                                     &
  dt     , vel    , vela   ,                                     &
  propce ,                                                       &
  imasfl , bmasfl ,                                              &
  tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
  ckupdc , smacel , frcxt  , grdphd ,                            &
  trava  , ximpa  , uvwk   , dfrcxt , dttens ,  trav  ,          &
  viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
  w1     , w7     , w8     , w9     , w10    )

! --- Sortie si pas de pression continuite
!       on met a jour les flux de masse, et on sort

if (iprco.le.0) then

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)

  call inimav                                                     &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom,                                                    &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   imasfl , bmasfl )

  ! In the ALE framework, we add the mesh velocity
  if (iale.eq.1) then

    allocate(mshvel(3,ncelet))

    do iel = 1, ncelet
      mshvel(1,iel) = rtp(iel,iuma)
      mshvel(2,iel) = rtp(iel,ivma)
      mshvel(3,iel) = rtp(iel,iwma)
    enddo

    ! One temporary array needed for internal faces, in case some internal vertices
    !  are moved directly by the user
    allocate(intflx(nfac), bouflx(ndimfb))

    itypfl = 1
    init   = 1
    inc    = 1
    iflmb0 = 1
    nswrgp = nswrgr(iuma)
    imligp = imligr(iuma)
    iwarnp = iwarni(iuma)
    epsrgp = epsrgr(iuma)
    climgp = climgr(iuma)

    call field_get_coefa_v(ivarfl(iuma), claale)
    call field_get_coefb_v(ivarfl(iuma), clbale)

    call inimav &
    !==========
  ( ivarfl(iuma)    , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
    iwarnp ,                                                       &
    epsrgp , climgp ,                                              &
    crom, brom,                                                    &
    mshvel ,                                                       &
    claale , clbale ,                                              &
    intflx , bouflx )

    ! Here we need of the opposite of the mesh velocity.
    !$omp parallel do if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      bmasfl(ifac) = bmasfl(ifac) - bouflx(ifac)
    enddo

    !$omp parallel do private(ddepx, ddepy, ddepz, icpt, ii, inod, &
    !$omp                     iel1, iel2, dtfac, rhofac)
    do ifac = 1, nfac
      ddepx = 0.d0
      ddepy = 0.d0
      ddepz = 0.d0
      icpt  = 0
      do ii = ipnfac(ifac),ipnfac(ifac+1)-1
        inod = nodfac(ii)
        icpt = icpt + 1
        ddepx = ddepx + disala(1,inod) + xyzno0(1,inod)-xyznod(1,inod)
        ddepy = ddepy + disala(2,inod) + xyzno0(2,inod)-xyznod(2,inod)
        ddepz = ddepz + disala(3,inod) + xyzno0(3,inod)-xyznod(3,inod)
      enddo
      ! Compute the mass flux using the nodes displacement
      if (iflxmw.eq.0) then
        ! For inner vertices, the mass flux due to the mesh displacement is
        !  recomputed from the nodes displacement
        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)
        dtfac = 0.5d0*(dt(iel1) + dt(iel2))
        rhofac = 0.5d0*(crom(iel1) + crom(iel2))
        imasfl(ifac) = imasfl(ifac) - rhofac*(                    &
              ddepx*surfac(1,ifac)                                &
             +ddepy*surfac(2,ifac)                                &
             +ddepz*surfac(3,ifac) )/dtfac/icpt
      else
        imasfl(ifac) = imasfl(ifac) - intflx(ifac)
      endif
    enddo

    ! Free memory
    deallocate(intflx, bouflx)
    deallocate(mshvel)

  endif

  ! Ajout de la vitesse du solide dans le flux convectif,
  ! si le maillage est mobile (solide rigide)
  ! En turbomachine, on connait exactement la vitesse de maillage a ajouter
  if (imobil.eq.1 .or. iturbo.eq.1 .or. iturbo.eq.2) then
    !$omp parallel do private(iel1, iel2, dtfac, rhofac, vitbox, vitboy, vitboz)
    do ifac = 1, nfac
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      if (irotce(iel1).ne.0 .or. irotce(iel2).ne.0) then
        dtfac  = 0.5d0*(dt(iel1) + dt(iel2))
        rhofac = 0.5d0*(crom(iel1) + crom(iel2))
        vitbox = rotax(2)*cdgfac(3,ifac) - rotax(3)*cdgfac(2,ifac)
        vitboy = rotax(3)*cdgfac(1,ifac) - rotax(1)*cdgfac(3,ifac)
        vitboz = rotax(1)*cdgfac(2,ifac) - rotax(2)*cdgfac(1,ifac)
        imasfl(ifac) = imasfl(ifac) - rhofac*(  vitbox*surfac(1,ifac)  &
                                              + vitboy*surfac(2,ifac)  &
                                              + vitboz*surfac(3,ifac))
      endif
    enddo
    !$omp parallel do private(iel, dtfac, rhofac, vitbox, vitboy, vitboz) &
    !$omp          if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      if (irotce(iel).ne.0) then
        dtfac  = dt(iel)
        rhofac = brom(ifac)
        vitbox = rotax(2)*cdgfbo(3,ifac) - rotax(3)*cdgfbo(2,ifac)
        vitboy = rotax(3)*cdgfbo(1,ifac) - rotax(1)*cdgfbo(3,ifac)
        vitboz = rotax(1)*cdgfbo(2,ifac) - rotax(2)*cdgfbo(1,ifac)
        bmasfl(ifac) = bmasfl(ifac) - rhofac*(  vitbox*surfbo(1,ifac)  &
                                              + vitboy*surfbo(2,ifac)  &
                                              + vitboz*surfbo(3,ifac))
      endif
    enddo
  endif

  ! Free memory
  !--------------
  deallocate(coefa_dp, coefb_dp)

  return

endif

!===============================================================================
! 6. Update mesh for unsteady turbomachinery computations
!===============================================================================

if (iturbo.eq.2 .and. iterns.eq.1) then

  ! Update mesh

  call turbomachinery_update_mesh (ttcmob, rs_ell(1))
  !==============================

  call dmtmps(t1)
  !==========

  do ifac = 1, nfabor
    isympa(ifac) = 1
  enddo

  ! Scratch and resize temporary internal faces arrays

  deallocate(viscf)
  if (idften(iu).eq.1) then
    allocate(viscf(1, 1, nfac))
  else if (idften(iu).eq.6) then
    allocate(viscf(3, 3, nfac))
  endif

  if (allocated(wvisfi)) then
    deallocate(viscfi)

    if (idften(iu).eq.1) then
      if (itytur.eq.3.and.irijnu.eq.1) then
        allocate(wvisfi(1,1,nfac))
        viscfi => wvisfi(:,:,1:nfac)
      else
        viscfi => viscf(:,:,1:nfac)
      endif
    else if(idften(iu).eq.6) then
      if (itytur.eq.3.and.irijnu.eq.1) then
        allocate(wvisfi(3,3,nfac))
        viscfi => wvisfi(1:3,1:3,1:nfac)
      else
        viscfi => viscf(1:3,1:3,1:nfac)
      endif
    endif

  endif

  if (allocated(secvif)) then
    deallocate(secvif)
    allocate(secvif(nfac))
  endif

  ! Scratch, resize and initialize main internal faces properties array

  call turbomachinery_reinit_i_face_fields

  ! Update local pointers on "internal faces" fields

  call field_get_val_s(iflmas, imasfl)

  if (irangp.ge.0 .or. iperio.eq.1) then

    ! Scratch and resize work arrays

    deallocate(w1, w7, w8, w9)
    allocate(w1(ncelet), w7(ncelet), w8(ncelet), w9(ncelet))
    if (allocated(w10)) then
      deallocate(w10)
      allocate(w10(ncelet))
    endif

    ! Resize auxiliary arrays (pointe module)

    call resize_aux_arrays
    !=====================

    ! Resize main real array

    call resize_main_real_array ( dt , rtp , rtpa , propce )
    !==========================

    ! Update turbomachinery module

    call turbomachinery_update
    !=========================

    ! Update field mappings ("owner" fields handled by turbomachinery_update)

    call fldtri(nproce, dt, rtpa, rtp, propce)
    !==========

    ! Resize other arrays related to the velocity-pressure resolution

    call resize_vec_real_array(trav)
    call resize_vec_real_array(dfrcxt)

    ! Resize other arrays, depending on user options

    if (iilagr.gt.0) &
      call resize_n_sca_real_arrays(ntersl, tslagr)

    if (iphydr.eq.1) then
      call resize_vec_real_array(frcxt)
    elseif (iphydr.eq.2) then
      call resize_sca_real_array(prhyd)
      call resize_vec_real_array_ni(grdphd)
    endif

    if (nterup.gt.1) then
      call resize_vec_real_array(trava)
      call resize_vec_real_array(uvwk)
      call resize_tens_real_array(ximpa)
    endif

    ! Update local pointers on "cells" fields

    call field_get_val_s(icrom, crom)

    call field_get_val_s(iprpfl(iviscl), viscl)
    call field_get_val_s(iprpfl(ivisct), visct)

    call field_get_val_v(ivarfl(iu), vel)
    call field_get_val_prev_v(ivarfl(iu), vela)

    call field_get_val_s(ivarfl(ipr), cvar_pr)
    call field_get_val_prev_s(ivarfl(ipr), cvara_pr)

    if (idtten.ge.0) call field_get_val_v(idtten, dttens)

  endif

  ! Update the Dirichlet wall boundary conditions for velocity (based on the
  ! solid body rotation on the new mesh).
  ! Note that the velocity BC update is made only if the user has not specified
  ! any specific Dirichlet condition for velocity.

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    if (coftur(ifac).lt.rinfin*0.5d0) then

      ! --- Physical Propreties
      visclc = viscl(iel)
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)

      ! Unit normal
      rnx = surfbo(1,ifac)/srfbnf
      rny = surfbo(2,ifac)/srfbnf
      rnz = surfbo(3,ifac)/srfbnf

      if (itytur.eq.3) then
        hint =   visclc         /distbf
      else
        hint = ( visclc+visctc )/distbf
      endif

      rcodcx = rotax(2)*cdgfbo(3,ifac) - rotax(3)*cdgfbo(2,ifac)
      rcodcy = rotax(3)*cdgfbo(1,ifac) - rotax(1)*cdgfbo(3,ifac)
      rcodcz = rotax(1)*cdgfbo(2,ifac) - rotax(2)*cdgfbo(1,ifac)

      ! Gradient boundary conditions (Dirichlet)
      !-----------------------------
      rcodcn = rcodcx*rnx+rcodcy*rny+rcodcz*rnz

      coefau(1,ifac) = (1.d0-coftur(ifac))*(rcodcx - rcodcn*rnx) + rcodcn*rnx
      coefau(2,ifac) = (1.d0-coftur(ifac))*(rcodcy - rcodcn*rny) + rcodcn*rny
      coefau(3,ifac) = (1.d0-coftur(ifac))*(rcodcz - rcodcn*rnz) + rcodcn*rnz

      ! Flux boundary conditions (Dirichlet)
      !-------------------------

      cofafu(1,ifac) = -hfltur(ifac)*(rcodcx - rcodcn*rnx) - hint*rcodcn*rnx
      cofafu(2,ifac) = -hfltur(ifac)*(rcodcy - rcodcn*rny) - hint*rcodcn*rny
      cofafu(3,ifac) = -hfltur(ifac)*(rcodcz - rcodcn*rnz) - hint*rcodcn*rnz

    endif

  enddo

  call dmtmps(t2)
  !==========

  rs_ell(2) = t2-t1

endif

!===============================================================================
! 7. Pressure correction step
!===============================================================================

if (iwarni(iu).ge.1) then
  write(nfecra,1200)
endif

! Allocate temporary arrays for the pressure resolution
allocate(drtp(ncelet))

if (ippmod(icompf).lt.0) then

  call resopv &
  !==========
( nvar   , ncetsm ,                                              &
  icetsm , isostd ,                                              &
  dt     , vel    ,                                              &
  propce ,                                                       &
  coefau , coefbu , coefa_dp        , coefb_dp ,                 &
  smacel ,                                                       &
  frcxt  , dfrcxt , dttens , trav   ,                            &
  viscf  , viscb  ,                                              &
  drtp   , tslagr ,                                              &
  trava  )

endif

!===============================================================================
! 8. Mesh velocity solving (ALE)
!===============================================================================

if (iale.eq.1) then

  if (itrale.gt.nalinf) then
    call alelav(rtp, rtpa, propce)
    !==========
  endif

endif

!===============================================================================
! 9. Update of the fluid velocity field
!===============================================================================

! Mass flux initialization for cavitating flows
if (icavit.ge.0) then
  do ifac = 1, nfac
    imasfl(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    bmasfl(ifac) = 0.d0
  enddo
endif

if (ippmod(icompf).lt.0) then

  ! irevmc = 0: Only the standard method is available for the coupled
  !              version of navstv.

  if (irevmc.eq.0) then

    ! The predicted velocity is corrected by the cell gradient of the
    ! pressure increment.

    ! GRADIENT DE L'INCREMENT TOTAL DE PRESSION

    if (idtvar.lt.0) then
      !$omp parallel do
      do iel = 1, ncel
        drtp(iel) = (cvar_pr(iel) -cvara_pr(iel)) / relaxv(ipr)
      enddo
    else
      !$omp parallel do
      do iel = 1, ncel
        drtp(iel) = cvar_pr(iel) -cvara_pr(iel)
      enddo
    endif

    ! --->    TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(drtp)
      !==========
    endif

    iccocg = 1
    inc = 0
    if (iphydr.eq.1.or.iifren.eq.1) inc = 1
    nswrgp = nswrgr(ipr)
    imligp = imligr(ipr)
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr)
    climgp = climgr(ipr)
    extrap = extrag(ipr)

    !Allocation
    allocate(gradp(ncelet,3))

    if (icavit.lt.0) then
      call grdpot &
      !==========
      ( ipr    , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
        iwarnp , epsrgp , climgp , extrap ,                            &
        dfrcxt ,                                                       &
        drtp   , coefa_dp        , coefb_dp        ,                   &
        gradp  )
    else
      allocate(xinvro(ncelet))
      do iel = 1, ncel
        xinvro(iel) = 1.d0/crom(iel)
      enddo

      call grdpre &
      !==========
      ( ipr    , imrgra , inc    , iccocg , nswrgp , imligp ,          &
        iwarnp , epsrgp , climgp , extrap ,                            &
        drtp   , xinvro , coefa_dp , coefb_dp ,                        &
        gradp  )

      deallocate(xinvro)
    endif

    thetap = thetav(ipr)
    !$omp parallel do private(isou)
    do iel = 1, ncelet
      do isou = 1, 3
        trav(isou,iel) = gradp(iel,isou)
      enddo
    enddo

    !Free memory
    deallocate(gradp)

    ! Update the velocity field
    !--------------------------
    thetap = thetav(ipr)

    ! Specific handling of hydrostatic pressure
    !------------------------------------------
    if (iphydr.eq.1) then

      ! Scalar diffusion for the pressure
      if (idften(ipr).eq.1) then
        !$omp parallel do private(dtsrom, isou)
        do iel = 1, ncel
          dtsrom = thetap*dt(iel)/crom(iel)
          do isou = 1, 3
            vel(isou,iel) = vel(isou,iel)                            &
                 + dtsrom*(dfrcxt(isou, iel)-trav(isou,iel))
          enddo
        enddo

      ! Tensorial diffusion for the pressure
      else if (idften(ipr).eq.6) then
        !$omp parallel do private(unsrom)
        do iel = 1, ncel
          unsrom = thetap/crom(iel)

            vel(1, iel) = vel(1, iel)                                             &
                 + unsrom*(                                                &
                   dttens(1,iel)*(dfrcxt(1, iel)-trav(1,iel))     &
                 + dttens(4,iel)*(dfrcxt(2, iel)-trav(2,iel))     &
                 + dttens(6,iel)*(dfrcxt(3, iel)-trav(3,iel))     &
                 )
            vel(2, iel) = vel(2, iel)                                             &
                 + unsrom*(                                                &
                   dttens(4,iel)*(dfrcxt(1, iel)-trav(1,iel))     &
                 + dttens(2,iel)*(dfrcxt(2, iel)-trav(2,iel))     &
                 + dttens(5,iel)*(dfrcxt(3, iel)-trav(3,iel))     &
                 )
            vel(3, iel) = vel(3, iel)                                             &
                 + unsrom*(                                                &
                   dttens(6,iel)*(dfrcxt(1 ,iel)-trav(1,iel))     &
                 + dttens(5,iel)*(dfrcxt(2 ,iel)-trav(2,iel))     &
                 + dttens(3,iel)*(dfrcxt(3 ,iel)-trav(3,iel))     &
                 )
        enddo
      endif

      ! Update external forces for the computation of the gradients
      !$omp parallel do
      do iel=1,ncel
        frcxt(1 ,iel) = frcxt(1 ,iel) + dfrcxt(1 ,iel)
        frcxt(2 ,iel) = frcxt(2 ,iel) + dfrcxt(2 ,iel)
        frcxt(3 ,iel) = frcxt(3 ,iel) + dfrcxt(3 ,iel)
      enddo
      if (irangp.ge.0.or.iperio.eq.1) then
        call synvin(frcxt)
        !==========
      endif
      ! Update of the Dirichlet boundary conditions on the
      ! pressure for the outlet
      call field_get_coefa_s(ivarfl(ipr), coefa_p)
      !$omp parallel do if(nfabor > thr_n_min)
      do ifac = 1, nfabor
        if (isostd(ifac).eq.1) then
          coefa_p(ifac) = coefa_p(ifac) + coefa_dp(ifac)
        endif
      enddo


      ! Standard handling of hydrostatic pressure
      !------------------------------------------
    else

      ! Scalar diffusion for the pressure
      if (idften(ipr).eq.1) then

      !$omp parallel do private(dtsrom, isou)
      do iel = 1, ncel
        dtsrom = thetap*dt(iel)/crom(iel)
        do isou = 1, 3
          vel(isou,iel) = vel(isou,iel) - dtsrom*trav(isou,iel)
        enddo
       enddo

      ! Tensorial diffusion for the pressure
      else if (idften(ipr).eq.6) then

      !$omp parallel do private(unsrom)
      do iel = 1, ncel
        unsrom = thetap/crom(iel)

          vel(1, iel) = vel(1, iel)                              &
                      - unsrom*(                                 &
                                 dttens(1,iel)*(trav(1,iel))     &
                               + dttens(4,iel)*(trav(2,iel))     &
                               + dttens(6,iel)*(trav(3,iel))     &
                               )
          vel(2, iel) = vel(2, iel)                              &
                      - unsrom*(                                 &
                                 dttens(4,iel)*(trav(1,iel))     &
                               + dttens(2,iel)*(trav(2,iel))     &
                               + dttens(5,iel)*(trav(3,iel))     &
                               )
          vel(3, iel) = vel(3, iel)                              &
                      - unsrom*(                                 &
                                 dttens(6,iel)*(trav(1,iel))     &
                               + dttens(5,iel)*(trav(2,iel))     &
                               + dttens(3,iel)*(trav(3,iel))     &
                               )
        enddo

      endif
    endif
  endif

endif

! In the ALE framework, we add the mesh velocity
if (iale.eq.1) then

  allocate(mshvel(3,ncelet))

  !$omp parallel do
  do iel = 1, ncelet
    mshvel(1,iel) = rtp(iel,iuma)
    mshvel(2,iel) = rtp(iel,ivma)
    mshvel(3,iel) = rtp(iel,iwma)
  enddo

  ! One temporary array needed for internal faces, in case some internal vertices
  !  are moved directly by the user
  allocate(intflx(nfac), bouflx(ndimfb))

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  nswrgp = nswrgr(iuma)
  imligp = imligr(iuma)
  iwarnp = iwarni(iuma)
  epsrgp = epsrgr(iuma)
  climgp = climgr(iuma)

  call field_get_coefa_v(ivarfl(iuma), claale)
  call field_get_coefb_v(ivarfl(iuma), clbale)

  call inimav &
  !==========
( ivarfl(iuma)    , itypfl ,                                     &
  iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
  iwarnp ,                                                       &
  epsrgp , climgp ,                                              &
  crom, brom,                                                    &
  mshvel ,                                                       &
  claale , clbale ,                                              &
  intflx , bouflx )

  ! Here we need of the opposite of the mesh velocity.
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    bmasfl(ifac) = bmasfl(ifac) - bouflx(ifac)
  enddo

  !$omp parallel do private(ddepx, ddepy, ddepz, icpt, ii, inod, &
  !$omp                     iel1, iel2, dtfac, rhofac)
  do ifac = 1, nfac
    ddepx = 0.d0
    ddepy = 0.d0
    ddepz = 0.d0
    icpt  = 0
    do ii = ipnfac(ifac),ipnfac(ifac+1)-1
      inod = nodfac(ii)
      icpt = icpt + 1
      ddepx = ddepx + disala(1,inod) + xyzno0(1,inod)-xyznod(1,inod)
      ddepy = ddepy + disala(2,inod) + xyzno0(2,inod)-xyznod(2,inod)
      ddepz = ddepz + disala(3,inod) + xyzno0(3,inod)-xyznod(3,inod)
    enddo
    ! Compute the mass flux using the nodes displacement
    if (iflxmw.eq.0) then
      ! For inner vertices, the mass flux due to the mesh displacement is
      !  recomputed from the nodes displacement
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      dtfac = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(crom(iel1) + crom(iel2))
      imasfl(ifac) = imasfl(ifac) - rhofac*(      &
            ddepx*surfac(1,ifac)                                &
           +ddepy*surfac(2,ifac)                                &
           +ddepz*surfac(3,ifac) )/dtfac/icpt
    else
      imasfl(ifac) = imasfl(ifac) - intflx(ifac)
    endif
  enddo

  ! Free memory
  deallocate(intflx, bouflx)
  deallocate(mshvel)

endif

!FIXME for me we should do that before predvv
! Ajout de la vitesse du solide dans le flux convectif,
! si le maillage est mobile (solide rigide)
! En turbomachine, on connait exactement la vitesse de maillage a ajouter
if (imobil.eq.1 .or. iturbo.eq.1 .or. iturbo.eq.2) then

  if (iturbo.eq.1.or.iturbo.eq.2) call dmtmps(t3)
                                  !==========

  !$omp parallel do private(iel1, iel2, dtfac, rhofac, vitbox, vitboy, vitboz)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    if (irotce(iel1).ne.0 .or. irotce(iel2).ne.0) then
      dtfac  = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(crom(iel1) + crom(iel2))
      vitbox = rotax(2)*cdgfac(3,ifac) - rotax(3)*cdgfac(2,ifac)
      vitboy = rotax(3)*cdgfac(1,ifac) - rotax(1)*cdgfac(3,ifac)
      vitboz = rotax(1)*cdgfac(2,ifac) - rotax(2)*cdgfac(1,ifac)
      imasfl(ifac) = imasfl(ifac) - rhofac*(  vitbox*surfac(1,ifac)  &
                                            + vitboy*surfac(2,ifac)  &
                                            + vitboz*surfac(3,ifac))
    endif
  enddo
  !$omp parallel do private(iel, dtfac, rhofac, vitbox, vitboy, vitboz) &
  !$omp             if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    if (irotce(iel).ne.0) then
      dtfac  = dt(iel)
      rhofac = brom(ifac)
      vitbox = rotax(2)*cdgfbo(3,ifac) - rotax(3)*cdgfbo(2,ifac)
      vitboy = rotax(3)*cdgfbo(1,ifac) - rotax(1)*cdgfbo(3,ifac)
      vitboz = rotax(1)*cdgfbo(2,ifac) - rotax(2)*cdgfbo(1,ifac)
      bmasfl(ifac) = bmasfl(ifac) - rhofac*(  vitbox*surfbo(1,ifac)  &
                                            + vitboy*surfbo(2,ifac)  &
                                            + vitboz*surfbo(3,ifac))
    endif
  enddo

  if (iturbo.eq.1.or.iturbo.eq.2) call dmtmps(t4)
                                  !==========

  rs_ell(2) = rs_ell(2) + t4-t3

endif

!===============================================================================
! 10. Cavitation: void fraction solving and update the mixture density/viscosity
!      and mass flux (resopv solved the convective flux of void fraction, divU)
!===============================================================================

if (icavit.ge.0) then

  ! Void fraction solving

  call resvoi(dt, rtp, rtpa, propce)
  !==========

  ! Halo synchronization

  call synsca(rtp(:,ivoidf))
  !==========

  ! Get the void fraction boundary conditions

  call field_get_coefa_s(ivarfl(ivoidf), coavoi)
  call field_get_coefb_s(ivarfl(ivoidf), cobvoi)

  ! Get the convective flux of the void fraction

  call field_get_key_int(ivarfl(ivoidf), kimasf, iflvoi)
  call field_get_key_int(ivarfl(ivoidf), kbmasf, iflvob)
  call field_get_val_s(iflvoi, ivoifl)
  call field_get_val_s(iflvob, bvoifl)

  ! Update mixture density/viscosity and mass flux

  call cavitation_update_phys_prop &
  !===============================
 ( rtp(:,ivoidf), coavoi, cobvoi, ivoifl, bvoifl, &
   crom, brom, propce(:,ipproc(iviscl)), imasfl, bmasfl )

  ! Verbosity

  if (mod(ntcabs,ntlist).eq.0.and.iterns.eq.nterup) then

    call field_get_val_prev_s(icrom, croma)

    call cavitation_print_mass_budget &
    !================================
   ( crom, croma, dt, imasfl, bmasfl )

  endif

endif

!===============================================================================
! 11. Compute error estimators for correction step and the global algo
!===============================================================================

if (iescal(iescor).gt.0.or.iescal(iestot).gt.0) then

  ! Allocate temporary arrays
  allocate(esflum(nfac), esflub(nfabor))

  ! ---> ECHANGE DES VITESSES ET PRESSION EN PERIODICITE ET PARALLELISME

  !    Pour les estimateurs IESCOR et IESTOT, la vitesse doit etre echangee.

  !    Pour l'estimateur IESTOT, la pression doit etre echangee aussi.

  !    Cela ne remplace pas l'echange du debut de pas de temps
  !     a cause de cs_user_extra_operations qui vient plus tard et des calculs suite)


  ! --- Vitesse

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(vel)
    !==========
  endif

  !  -- Pression

  if (iescal(iestot).gt.0) then

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(cvar_pr)
      !==========
    endif

  endif

  ! ---> CALCUL DU FLUX DE MASSE DEDUIT DE LA VITESSE REACTUALISEE

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)

  call inimav                                                     &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom,                                                    &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   esflum , esflub )


  ! ---> CALCUL DE L'ESTIMATEUR CORRECTION : DIVERGENCE DE ROM * U (N + 1)
  !                                          - GAMMA

  if (iescal(iescor).gt.0) then
    init = 1
    call divmas(init, esflum, esflub, w1)

    if (ncetsm.gt.0) then
      !$omp parallel do private(iel) if(ncetsm > thr_n_min)
      do iitsm = 1, ncetsm
        iel = icetsm(iitsm)
        w1(iel) = w1(iel)-volume(iel)*smacel(iitsm,ipr)
      enddo
    endif

    if (iescal(iescor).eq.2) then
      iescop = ipproc(iestim(iescor))
      !$omp parallel do
      do iel = 1, ncel
        propce(iel,iescop) =  abs(w1(iel))
      enddo
    elseif (iescal(iescor).eq.1) then
      iescop = ipproc(iestim(iescor))
      !$omp parallel do
      do iel = 1, ncel
        propce(iel,iescop) =  abs(w1(iel)) / volume(iel)
      enddo
    endif
  endif


  ! ---> CALCUL DE L'ESTIMATEUR TOTAL

  if (iescal(iestot).gt.0) then

    !   INITIALISATION DE TRAV AVEC LE TERME INSTATIONNAIRE

    !$omp parallel do private(rovolsdt, isou)
    do iel = 1, ncel
      rovolsdt = crom(iel)*volume(iel)/dt(iel)
      do isou = 1, 3
        trav(isou,iel) = rovolsdt *                               &
                 ( vela(isou,iel)- vel(isou,iel) )
      enddo
    enddo

    !   APPEL A PREDVV AVEC RTP ET RTP AU LIEU DE RTP ET RTPA
    !                  AVEC LE FLUX DE MASSE RECALCULE
    iappel = 2
    call predvv &
    !==========
 ( iappel ,                                                       &
   nvar   , nscal  , iterns , ncepdc , ncetsm ,                   &
   icepdc , icetsm , itypsm ,                                     &
   dt     , vel    , vel    ,                                     &
   propce ,                                                       &
   esflum , esflub ,                                              &
   tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
   ckupdc , smacel , frcxt  , grdphd ,                            &
   trava  , ximpa  , uvwk   , dfrcxt , dttens , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
   w1     , w7     , w8     , w9     , w10    )

  endif

  deallocate(esflum, esflub)
endif

!===============================================================================
! 12. Loop on the velocity/Pressure coupling (PISO)
!===============================================================================

if (nterup.gt.1) then
! TEST DE CONVERGENCE DE L'ALGORITHME ITERATIF
! On initialise ICVRGE a 1 et on le met a 0 si on n'a pas convergee

  icvrge = 1

  xnrtmp = 0.d0
  !$omp parallel do reduction(+:xnrtmp) private(xdu, xdv, xdw)
  do iel = 1,ncel
    xdu = vel(1,iel) - uvwk(1,iel)
    xdv = vel(2,iel) - uvwk(2,iel)
    xdw = vel(3,iel) - uvwk(3,iel)
    xnrtmp = xnrmu +(xdu**2 + xdv**2 + xdw**2) * volume(iel)
  enddo
  xnrmu = xnrtmp
  ! --->    TRAITEMENT DU PARALLELISME

  if (irangp.ge.0) call parsom (xnrmu)
                   !==========
  ! -- >    TRAITEMENT DU COUPLAGE ENTRE DEUX INSTANCES DE CODE_SATURNE
  do numcpl = 1, nbrcpl
    call tbrcpl ( numcpl, 1, 1, xnrmu, xnrdis )
    !==========
    xnrmu = xnrmu + xnrdis
  enddo
  xnrmu = sqrt(xnrmu)

  ! Indicateur de convergence du point fixe
  if (xnrmu.ge.epsup*xnrmu0) icvrge = 0

endif

! ---> RECALAGE DE LA PRESSION SUR UNE PRESSION A MOYENNE NULLE
!  On recale si on n'a pas de Dirichlet. Or le nombre de Dirichlets
!  calcule dans typecl.F est NDIRCL si IDIRCL=1 et NDIRCL-1 si IDIRCL=0
!  (ISTAT vaut toujours 0 pour la pression)

if (idircl(ipr).eq.1) then
  ndircp = ndircl(ipr)
else
  ndircp = ndircl(ipr)-1
endif
if (ndircp.le.0) then
  call prmoy0 &
  !==========
( ncelet , ncel   , volume , cvar_pr )
endif

! Calcul de la pression totale IPRTOT : (definie comme propriete )
! En compressible, la pression resolue est deja la pression totale

if (ippmod(icompf).lt.0) then
  call field_get_val_s(iprpfl(iprtot), cpro_prtot)
  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)
  do iel=1,ncel
    cpro_prtot(iel)= cvar_pr(iel)               &
         + ro0*( gx*(xyzcen(1,iel)-xxp0)               &
         + gy*(xyzcen(2,iel)-xyp0)                     &
         + gz*(xyzcen(3,iel)-xzp0) )                   &
         + p0 - pred0
  enddo
endif

!===============================================================================
! 13. Printing
!===============================================================================

if (iwarni(iu).ge.1) then

  write(nfecra,2000)

  rnorm = -1.d0
  do iel = 1, ncel
    rnorm  = max(rnorm,abs(cvar_pr(iel)))
  enddo
  if (irangp.ge.0) call parmax (rnorm)
                   !==========
  write(nfecra,2100)rnorm

  rnorm = -1.d0
  imax = 1
  !$omp parallel private(vitnor, rnormt, imaxt)
  rnormt = -1.d0
  !$omp do
  do iel = 1, ncel
    vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
    if (vitnor.ge.rnormt) then
      rnormt = vitnor
      imaxt  = iel
    endif
  enddo
  !$omp critical
  if (rnormt .gt. rnorm) then
    rnormt = rnorm
    imax = imaxt
  endif
  !$omp end critical
  !$omp end parallel

  xyzmax(1) = xyzcen(1,imax)
  xyzmax(2) = xyzcen(2,imax)
  xyzmax(3) = xyzcen(3,imax)

  if (irangp.ge.0) then
    nbrval = 3
    call parmxl (nbrval, rnorm, xyzmax)
    !==========
  endif

  write(nfecra,2200) rnorm,xyzmax(1),xyzmax(2),xyzmax(3)

  ! Pour la periodicite et le parallelisme, rom est echange dans phyvar

  rnorma = -grand
  rnormi =  grand
  !$omp parallel do reduction(max: rnorma) reduction(min: rnormi)         &
  !$omp             private(iel1, iel2, surf, rhom, rnorm)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    surf = surfan(ifac)
    rhom = (crom(iel1)+crom(iel2))*0.5d0
    rnorm = imasfl(ifac)/(surf*rhom)
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    !==========
    call parmin (rnormi)
    !==========
  endif
  write(nfecra,2300)rnorma, rnormi

  rnorma = -grand
  rnormi =  grand
  do ifac = 1, nfabor
    rnorm = bmasfl(ifac)/(surfbn(ifac)*brom(ifac))
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    !==========
    call parmin (rnormi)
    !==========
  endif
  write(nfecra,2400)rnorma, rnormi

  rnorm = 0.d0
  !$omp parallel do reduction(+: rnorm) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    rnorm = rnorm + bmasfl(ifac)
  enddo

  if (irangp.ge.0) call parsom (rnorm)
                   !==========

  write(nfecra,2500)rnorm

  write(nfecra,2001)

  if (nterup.gt.1) then
    if (icvrge.eq.0) then
      write(nfecra,2600) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
      if (iterns.eq.nterup) then
        write(nfecra,2603)
        write(nfecra,2001)
      endif
    else
      write(nfecra,2602) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
    endif
  endif

endif

if (iturbo.eq.2) then
  if (mod(ntcabs,ntlist).eq.0.and.iterns.eq.nterup) then
     write(nfecra,3000) rs_ell(1), rs_ell(1) + rs_ell(2)
     rs_ell(1) = 0.d0
     rs_ell(2) = 0.d0
  endif
endif

! Free memory
deallocate(viscf, viscb)
deallocate(drtp)
deallocate(trav)
deallocate(dfrcxt)
deallocate(w1)
deallocate(w7, w8, w9)
if (allocated(w10)) deallocate(w10)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
if (allocated(secvif)) deallocate(secvif, secvib)
if (iphydr.eq.2) deallocate(grdphd)

! Free memory
!--------------
deallocate(coefa_dp, coefb_dp)

!--------
! Formats
!--------
#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VITESSE                             ',/,&
'      --------------------------                             ',/)
 1080 format(/,                                                   &
'   ** RESOLUTION DE L''EQUATION DE MASSE                     ',/,&
'      ----------------------------------                     ',/)
 1200 format(/,                                                   &
'   ** RESOLUTION POUR LA PRESSION CONTINUITE                 ',/,&
'      --------------------------------------                 ',/)
 2000 format(/,' APRES PRESSION CONTINUITE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Pression max.',E12.4   ,' (max. de la valeur absolue)       ',/)
 2200 format(                                                           &
' Vitesse  max.',E12.4   ,' en',3E11.3                         ,/)
 2300 format(                                                           &
' Vitesse  en face interne max.',E12.4   ,' ; min.',E12.4        )
 2400 format(                                                           &
' Vitesse  en face de bord max.',E12.4   ,' ; min.',E12.4        )
 2500 format(                                                           &
' Bilan de masse   au bord   ',E14.6                             )
 2600 format(                                                           &
' Informations Point fixe a l''iteration :',I10                ,/)
 2601 format('norme = ',E12.4,' norme 0 = ',E12.4,' toler  = ',E12.4 ,/)
 2602 format(                                                           &
' Convergence du point fixe a l''iteration ',I10               ,/)
 2603 format(                                                           &
' Non convergence du couplage vitesse pression par point fixe  ' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)
 3000 format(/,                                                     &
'   ** INFORMATION SUR LE TRAITEMENT ROTOR/STATOR INSTATIONNAIRE',/,&
'      ---------------------------------------------------------',/,&
' Temps dedie a la mise a jour du maillage (s) :',F12.4,          /,&
' Temps total                              (s) :',F12.4,          /)

#else

 1000 format(/,                                                   &
'   ** SOLVING VELOCITY'                                       ,/,&
'      ----------------'                                       ,/)
 1080 format(/,                                                   &
'   ** SOLVING MASS BALANCE EQUATION                          ',/,&
'      -----------------------------                          ',/)
 1200 format(/,                                                   &
'   ** SOLVING CONTINUITY PRESSURE'                            ,/,&
'      ---------------------------'                            ,/)
 2000 format(/,' AFTER CONTINUITY PRESSURE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Max. pressure',E12.4   ,' (max. absolute value)'             ,/)
 2200 format(                                                           &
' Max. velocity',E12.4   ,' en',3E11.3                         ,/)
 2300 format(                                                           &
' Max. velocity at interior face',E12.4   ,' ; min.',E12.4       )
 2400 format(                                                           &
' Max. velocity at boundary face',E12.4   ,' ; min.',E12.4       )
 2500 format(                                                           &
' Mass balance  at boundary  ',E14.6                             )
 2600 format(                                                           &
' Fixed point informations at iteration:',I10                  ,/)
 2601 format('norm = ',E12.4,' norm 0 = ',E12.4,' toler  = ',E12.4   ,/)
 2602 format(                                                           &
' Fixed point convergence at iteration ',I10                   ,/)
 2603 format(                                                           &
' Non convergence of fixed point for velocity pressure coupling' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)
 3000 format(/,                                             &
'   ** INFORMATION ON UNSTEADY ROTOR/STATOR TREATMENT',/,&
'      ----------------------------------------------',/,&
' Time dedicated to mesh update (s):',F12.4,           /,&
' Global time                   (s):',F12.4,           /)

#endif

!----
! End
!----

return

end subroutine
