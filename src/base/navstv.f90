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
!> \param[in]     tpucou        velocity-pressure coupling (interleaved)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     coefa, coefb  boundary conditions
!> \param[in]     frcxt         external force generating the hydrostatic
!>                              pressure
!> \param[in]     prhyd         hydrostatic pressure predicted at cell centers
!> \param[in]     tslagr        terme de couplage retour du
!>                                  lagrangien
!> \param[in]     trava         work array for pressure velocity coupling
!> \param[in]     ximpa         work array for pressure velocity coupling
!> \param[in]     uvwk          work array for pressure velocity coupling
!_______________________________________________________________________________


subroutine navstv &
 ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
   isostd ,                                                       &
   dt     , tpucou , rtp    , rtpa   , propce , propfa , propfb , &
   tslagr , coefa  , coefb  , frcxt  , prhyd  ,                   &
   trava  , ximpa  , uvwk   )

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
use pointe
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iterns , icvrge , itrale

integer          isostd(nfabor+1)

double precision dt(ncelet), tpucou(ncelet,3), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision frcxt(3,ncelet)
double precision prhyd(ncelet)
double precision trava(ndim,ncelet),ximpa(ndim,ndim,ncelet)
double precision uvwk(ndim,ncelet)

! Local variables

integer          iccocg, inc, iel, iel1, iel2, ifac, imax, imaxt
integer          ii    , inod, itypfl
integer          isou, ivar, iitsm
integer          iclipr, iclipf
integer          icliup, iclivp, icliwp, init
integer          iflmas, iflmab, ipcrom, ipbrom
integer          iflms1, iflmb1, iflmb0
integer          nswrgp, imligp, iwarnp
integer          nbrval, iappel, iescop
integer          ndircp, icpt
integer          numcpl
double precision rnorm , rnormt, rnorma, rnormi, vitnor
double precision dtsrom, unsrom, surf  , rhom, rovolsdt
double precision epsrgp, climgp, extrap, xyzmax(3)
double precision thetap, xdu, xdv, xdw
double precision xxp0 , xyp0 , xzp0
double precision rhofac, dtfac, ddepx , ddepy, ddepz
double precision xnrdis
double precision vitbox, vitboy, vitboz

double precision rvoid(1)

double precision, allocatable, dimension(:,:,:), target :: viscf
double precision, allocatable, dimension(:), target :: viscb
double precision, allocatable, dimension(:,:,:), target :: wvisfi
double precision, allocatable, dimension(:), target :: wvisbi
double precision, allocatable, dimension(:) :: drtp, smbr
double precision, allocatable, dimension(:,:) :: trav
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:,:) :: dfrcxt
double precision, allocatable, dimension(:,:) :: grdphd
double precision, allocatable, dimension(:) :: esflum, esflub
double precision, allocatable, dimension(:) :: intflx, bouflx
double precision, allocatable, dimension(:) :: secvif, secvib

double precision, pointer, dimension(:,:,:) :: viscfi => null()
double precision, pointer, dimension(:) :: viscbi => null()

double precision, dimension(:,:), allocatable :: gradp
double precision, dimension(:,:), allocatable :: vel
double precision, dimension(:,:), allocatable :: vela
double precision, dimension(:,:), allocatable :: mshvel
double precision, dimension(:), allocatable :: coefap

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

! Allocate temporary arrays for the velocity-pressure resolution
if (idften(iu).eq.1) then
  allocate(viscf(1, 1, nfac), viscb(nfabor))
else if (idften(iu).eq.6) then
  allocate(viscf(3, 3, nfac), viscb(nfabor))
endif

allocate(drtp(ncelet))
allocate(trav(3,ncelet))
allocate(vela(3,ncelet))
allocate(vel(3,ncelet))

! Allocate other arrays, depending on user options

! Array for delta p gradient boundary conditions
allocate(coefap(nfabor))

allocate(dfrcxt(3,ncelet))
if (iphydr.eq.2) allocate(grdphd(ncelet,ndim))
if (iescal(iestot).gt.0) allocate(esflum(nfac), esflub(nfabor))
if (idften(iu).eq.1) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(1,1,nfac), wvisbi(nfabor))
    viscfi => wvisfi(:,:,1:nfac)
    viscbi => wvisbi(1:nfabor)
  else
    viscfi => viscf(:,:,1:nfac)
    viscbi => viscb(1:nfabor)
  endif
else if(idften(iu).eq.6) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(3,3,nfac), wvisbi(nfabor))
    viscfi => wvisfi(1:3,1:3,1:nfac)
    viscbi => wvisbi(1:nfabor)
  else
    viscfi => viscf(1:3,1:3,1:nfac)
    viscbi => viscb(1:nfabor)
  endif
endif

if (ivisse.eq.1) then
  allocate(secvif(nfac),secvib(nfabor))
endif

! Allocate work for the ALE module
if (iale.eq.1) allocate(mshvel(3,ncelet))

! Allocate work arrays
allocate(w1(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
if (irnpnw.eq.1) allocate(w10(ncelet))

! Interleaved value of vel and vela
!$omp parallel do
do iel = 1, ncelet
  vel (1,iel) = rtp (iel,iu)
  vel (2,iel) = rtp (iel,iv)
  vel (3,iel) = rtp (iel,iw)
  vela(1,iel) = rtpa(iel,iu)
  vela(2,iel) = rtpa(iel,iv)
  vela(3,iel) = rtpa(iel,iw)
enddo

if (iwarni(iu).ge.1) then
  write(nfecra,1000)
endif

! Initialize variables to avoid compiler warnings

ivar = 0
iflmas = 0
ipcrom = 0
imax = 0

! Memoire

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
    xnrmu0 = 0.d0
    !$omp parallel do reduction(+:xnrmu0)
    do iel = 1, ncel
      xnrmu0 = xnrmu0 +(vela(1,iel)**2        &
                      + vela(2,iel)**2        &
                      + vela(3,iel)**2)       &
                      * volume(iel)
    enddo
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
      call synsca(rtpa(1,ipr))
      !==========
    endif
  endif

endif

!===============================================================================
! 1. Prediction of the mass flux in case of Low Mach compressible algorithm
!===============================================================================

if ((idilat.eq.2.or.idilat.eq.3).and. &
    (ntcabs.gt.1.or.isuite.gt.0)) then

  call predfl &
  !==========
  ( nvar   , nscal  , ncetsm ,                            &
    icetsm ,                                              &
    dt     , rtp    , rtpa   ,                            &
    propce , propfa , propfb ,                            &
    smacel )

endif

!===============================================================================
! 2. Hydrostatic pressure prediction in case of Low Mach compressible algorithm
!===============================================================================

if (iphydr.eq.2) then

  call prehyd &
  !==========
  ( nvar   , nscal  ,                                     &
    dt     , rtp    , rtpa   ,                            &
    propce , propfa , propfb ,                            &
    prhyd  , grdphd )

endif

!===============================================================================
! 3.  ETAPE DE PREDICTION DES VITESSES
!===============================================================================

iappel = 1
iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))

call predvv &
!==========
( iappel ,                                                       &
  nvar   , nscal  , iterns ,                                     &
  ncepdc , ncetsm ,                                              &
  icepdc , icetsm , itypsm ,                                     &
  dt     , rtp    , rtpa   , vel    , vela   ,                   &
  propce , propfa , propfb ,                                     &
  propfa(1,iflmas), propfb(1,iflmab),                            &
  tslagr , coefa  , coefb  , coefau , coefbu , cofafu , cofbfu , &
  ckupdc , smacel , frcxt  , grdphd ,                            &
  trava  , ximpa  , uvwk   , dfrcxt , dttens ,  trav  ,          &
  viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
  w1     , w7     , w8     , w9     , w10    )

! --- Sortie si pas de pression continuite
!       on met a jour les flux de masse, et on sort

if (iprco.le.0) then

  icliup = iclrtp(iu,icoef)
  iclivp = iclrtp(iv,icoef)
  icliwp = iclrtp(iw,icoef)

  iflmas = ipprof(ifluma(iu))
  iflmab = ipprob(ifluma(iu))
  ipcrom = ipproc(irom)
  ipbrom = ipprob(irom)

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
  extrap = extrag(iu)

  call inimav                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iu     , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   propfa(1,iflmas), propfb(1,iflmab)  )

  ! In the ALE framework, we add the mesh velocity
  if (iale.eq.1) then

    do iel = 1, ncelet
      mshvel(1,iel) = rtp(iel,iuma)
      mshvel(2,iel) = rtp(iel,ivma)
      mshvel(3,iel) = rtp(iel,iwma)
    enddo

    ! One temporary array needed for internal faces, in case some internal vertices
    !  are moved directly by the user
    allocate(intflx(nfac), bouflx(nfabor))

    iflmas = ipprof(ifluma(iu))
    iflmab = ipprob(ifluma(iu))
    ipcrom = ipproc(irom)
    ipbrom = ipprob(irom)

    itypfl = 1
    init   = 1
    inc    = 1
    iflmb0 = 1
    nswrgp = nswrgr(iuma)
    imligp = imligr(iuma)
    iwarnp = iwarni(iuma)
    epsrgp = epsrgr(iuma)
    climgp = climgr(iuma)
    extrap = extrag(iuma)

    call inimav &
    !==========
  ( nvar   , nscal  ,                                              &
    iu     , itypfl ,                                              &
    iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
    iwarnp , nfecra ,                                              &
    epsrgp , climgp , extrap ,                                     &
    propce(1,ipcrom), propfb(1,ipbrom),                            &
    mshvel ,                                                       &
    claale , clbale ,                                              &
    intflx , bouflx )

    ! Here we need of the opposite of the mesh velocity.
    !$omp parallel do if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      propfb(ifac,iflmab) = propfb(ifac,iflmab) - bouflx(ifac)
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
        rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
        propfa(ifac,iflmas) = propfa(ifac,iflmas) - rhofac*(      &
              ddepx*surfac(1,ifac)                                &
             +ddepy*surfac(2,ifac)                                &
             +ddepz*surfac(3,ifac) )/dtfac/icpt
      else
        propfa(ifac,iflmas) = propfa(ifac,iflmas) - intflx(ifac)
      endif
    enddo

    ! Free memory
    deallocate(intflx, bouflx)

  endif

  ! Ajout de la vitesse du solide dans le flux convectif,
  ! si le maillage est mobile (solide rigide)
  ! En turbomachine, on connaît exactement la vitesse de maillage à ajouter
  if (imobil.eq.1) then

    iflmas = ipprof(ifluma(iu))
    iflmab = ipprob(ifluma(iu))
    ipcrom = ipproc(irom)
    ipbrom = ipprob(irom)

    !$omp parallel do private(iel1, iel2, dtfac, rhofac, vitbox, vitboy, vitboz)
    do ifac = 1, nfac
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      dtfac  = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
      vitbox = omegay*cdgfac(3,ifac) - omegaz*cdgfac(2,ifac)
      vitboy = omegaz*cdgfac(1,ifac) - omegax*cdgfac(3,ifac)
      vitboz = omegax*cdgfac(2,ifac) - omegay*cdgfac(1,ifac)
      propfa(ifac,iflmas) = propfa(ifac,iflmas) - rhofac*(        &
           vitbox*surfac(1,ifac) + vitboy*surfac(2,ifac) + vitboz*surfac(3,ifac) )
    enddo
    !$omp parallel do private(iel, dtfac, rhofac, vitbox, vitboy, vitboz) &
    !$omp          if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      dtfac  = dt(iel)
      rhofac = propfb(ifac,ipbrom)
      vitbox = omegay*cdgfbo(3,ifac) - omegaz*cdgfbo(2,ifac)
      vitboy = omegaz*cdgfbo(1,ifac) - omegax*cdgfbo(3,ifac)
      vitboz = omegax*cdgfbo(2,ifac) - omegay*cdgfbo(1,ifac)
      propfb(ifac,iflmab) = propfb(ifac,iflmab) - rhofac*(        &
           vitbox*surfbo(1,ifac) + vitboy*surfbo(2,ifac) + vitboz*surfbo(3,ifac) )
    enddo

  endif

  ! Interleaved values of vel and vela

  !$omp parallel do
  do iel = 1, ncelet
    rtp (iel,iu) = vel (1,iel)
    rtp (iel,iv) = vel (2,iel)
    rtp (iel,iw) = vel (3,iel)
    rtpa(iel,iu) = vela(1,iel)
    rtpa(iel,iv) = vela(2,iel)
    rtpa(iel,iw) = vela(3,iel)
    ! Store the diagonal part of dttens for postprocessing purpose
    if (ipucou.eq.1 .or. ncpdct.gt.0) then
      tpucou(iel,1) = dttens(1,iel)
      tpucou(iel,2) = dttens(2,iel)
      tpucou(iel,3) = dttens(3,iel)
    endif
  enddo

  ! Free memory
  !--------------
  deallocate(vel)
  deallocate(vela)
  deallocate(coefap)

  return

endif

!===============================================================================
! 4.  ETAPE DE PRESSION/CONTINUITE ( VITESSE/PRESSION )
!===============================================================================

if (iwarni(iu).ge.1) then
  write(nfecra,1200)
endif

call resopv &
!==========
 ( nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   icepdc , icetsm , itypsm , isostd ,                            &
   dt     , rtp    , rtpa   , vel    , vela   ,                   &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  , coefau , coefbu , coefap ,                   &
   ckupdc , smacel ,                                              &
   frcxt  , dfrcxt , dttens , trav   ,                            &
   viscf  , viscb  , viscfi , viscbi ,                            &
   drtp   , tslagr ,                                              &
   trava  )

!===============================================================================
! 5.  RESOLUTION DE LA VITESSE DE MAILLAGE EN ALE
!===============================================================================

if (iale.eq.1) then

  if (itrale.gt.nalinf) then

    call alelav &
    !==========
  ( nvar   , nscal  ,                                              &
    dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
    coefa  , coefb  )

  endif

endif

!===============================================================================
! 6.  REACTUALISATION DU CHAMP DE VITESSE
!===============================================================================

iclipr = iclrtp(ipr,icoef)
iclipf = iclrtp(ipr,icoeff)
icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))
ipcrom = ipproc(irom  )
ipbrom = ipprob(irom  )

! irevmc = 0: Only the standard method is available for the coupled
!              version of navstv.

if (irevmc.eq.0) then

  ! The predicted velocity is corrected by the cell gradient of the
  ! pressure increment.

  ! GRADIENT DE L'INCREMENT TOTAL DE PRESSION

  if (idtvar.lt.0) then
    !$omp parallel do
    do iel = 1, ncel
      drtp(iel) = (rtp(iel,ipr) -rtpa(iel,ipr)) / relaxv(ipr)
    enddo
  else
    !$omp parallel do
    do iel = 1, ncel
      drtp(iel) = rtp(iel,ipr) -rtpa(iel,ipr)
    enddo
  endif

  ! --->    TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(drtp)
    !==========
  endif

  iccocg = 1
  inc = 0
  if (iphydr.eq.1) inc = 1
  nswrgp = nswrgr(ipr)
  imligp = imligr(ipr)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(ipr)
  climgp = climgr(ipr)
  extrap = extrag(ipr)

  !Allocation
  allocate(gradp(ncelet,3))

  call grdpot &
  !==========
 ( ipr , imrgra , inc    , iccocg , nswrgp , imligp , iphydr ,    &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   dfrcxt ,                                                       &
   drtp   , coefap        , coefb(1,iclipr)  ,                    &
   gradp  )

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
        dtsrom = thetap*dt(iel)/propce(iel,ipcrom)
        do isou = 1, 3
          vel(isou,iel) = vel(isou,iel)                            &
                        + dtsrom*(dfrcxt(isou, iel)-trav(isou,iel))
        enddo
      enddo

    ! Tensorial diffusion for the pressure
    else if (idften(ipr).eq.6) then
      !$omp parallel do private(unsrom)
      do iel = 1, ncel
        unsrom = thetap/propce(iel,ipcrom)

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
    ! Update of the Direchlet boundary conditions on the pressure for the outlet
    iclipr = iclrtp(ipr,icoef)
    iclipf = iclrtp(ipr,icoeff)
    !$omp parallel do if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      if (isostd(ifac).eq.1) then
        coefa(ifac,iclipr) = coefa(ifac,iclipr) + coefap(ifac)
      endif
    enddo


  ! Standard handling of hydrostatic pressure
  !------------------------------------------
  else

    ! Scalar diffusion for the pressure
    if (idften(ipr).eq.1) then

      !$omp parallel do private(dtsrom, isou)
      do iel = 1, ncel
        dtsrom = thetap*dt(iel)/propce(iel,ipcrom)
        do isou = 1, 3
          vel(isou,iel) = vel(isou,iel) - dtsrom*trav(isou,iel)
        enddo
      enddo

    ! Tensorial diffusion for the pressure
    else if (idften(ipr).eq.6) then

      !$omp parallel do private(unsrom)
      do iel = 1, ncel
        unsrom = thetap/propce(iel,ipcrom)

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

! In the ALE framework, we add the mesh velocity
if (iale.eq.1) then

  !$omp parallel do
  do iel = 1, ncelet
    mshvel(1,iel) = rtp(iel,iuma)
    mshvel(2,iel) = rtp(iel,ivma)
    mshvel(3,iel) = rtp(iel,iwma)
  enddo

  ! One temporary array needed for internal faces, in case some internal vertices
  !  are moved directly by the user
  allocate(intflx(nfac), bouflx(nfabor))

  iflmas = ipprof(ifluma(iu))
  iflmab = ipprob(ifluma(iu))
  ipcrom = ipproc(irom)
  ipbrom = ipprob(irom)

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  nswrgp = nswrgr(iuma)
  imligp = imligr(iuma)
  iwarnp = iwarni(iuma)
  epsrgp = epsrgr(iuma)
  climgp = climgr(iuma)
  extrap = extrag(iuma)

  call inimav &
  !==========
( nvar   , nscal  ,                                              &
  iuma   , itypfl ,                                              &
  iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  propce(1,ipcrom), propfb(1,ipbrom),                            &
  mshvel ,                                                       &
  claale , clbale ,                                              &
  intflx , bouflx )

  ! Here we need of the opposite of the mesh velocity.
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    propfb(ifac,iflmab) = propfb(ifac,iflmab) - bouflx(ifac)
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
      rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
      propfa(ifac,iflmas) = propfa(ifac,iflmas) - rhofac*(      &
            ddepx*surfac(1,ifac)                                &
           +ddepy*surfac(2,ifac)                                &
           +ddepz*surfac(3,ifac) )/dtfac/icpt
    else
      propfa(ifac,iflmas) = propfa(ifac,iflmas) - intflx(ifac)
    endif
  enddo

  ! Free memory
  deallocate(intflx, bouflx)

endif

!FIXME for me we should do that before predvv
! Ajout de la vitesse du solide dans le flux convectif,
! si le maillage est mobile (solide rigide)
! En turbomachine, on connaît exactement la vitesse de maillage à ajouter

if (imobil.eq.1) then

  iflmas = ipprof(ifluma(iu))
  iflmab = ipprob(ifluma(iu))
  ipcrom = ipproc(irom)
  ipbrom = ipprob(irom)

  !$omp parallel do private(iel1, iel2, dtfac, rhofac, vitbox, vitboy, vitboz)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    dtfac  = 0.5d0*(dt(iel1) + dt(iel2))
    rhofac = 0.5d0*(propce(iel1,ipcrom) + propce(iel2,ipcrom))
    vitbox = omegay*cdgfac(3,ifac) - omegaz*cdgfac(2,ifac)
    vitboy = omegaz*cdgfac(1,ifac) - omegax*cdgfac(3,ifac)
    vitboz = omegax*cdgfac(2,ifac) - omegay*cdgfac(1,ifac)
    propfa(ifac,iflmas) = propfa(ifac,iflmas) - rhofac*(        &
         vitbox*surfac(1,ifac) + vitboy*surfac(2,ifac) + vitboz*surfac(3,ifac) )
  enddo
  !$omp parallel do private(iel, dtfac, rhofac, vitbox, vitboy, vitboz) &
  !$omp             if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    dtfac  = dt(iel)
    rhofac = propfb(ifac,ipbrom)
    vitbox = omegay*cdgfbo(3,ifac) - omegaz*cdgfbo(2,ifac)
    vitboy = omegaz*cdgfbo(1,ifac) - omegax*cdgfbo(3,ifac)
    vitboz = omegax*cdgfbo(2,ifac) - omegay*cdgfbo(1,ifac)
    propfb(ifac,iflmab) = propfb(ifac,iflmab) - rhofac*(        &
         vitbox*surfbo(1,ifac) + vitboy*surfbo(2,ifac) + vitboz*surfbo(3,ifac) )
  enddo
endif

!===============================================================================
! 7.  CALCUL D'UN ESTIMATEUR D'ERREUR DE L'ETAPE DE CORRECTION ET TOTAL
!===============================================================================

if (iescal(iescor).gt.0.or.iescal(iestot).gt.0) then

  ! ---> REPERAGE DES VARIABLES

  icliup = iclrtp(iu,icoef)
  iclivp = iclrtp(iv,icoef)
  icliwp = iclrtp(iw,icoef)

  ipcrom = ipproc(irom)
  ipbrom = ipprob(irom)


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
      call synsca(rtp(1,ipr))
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
  extrap = extrag(iu)

  call inimav                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iu     , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   esflum , esflub )


  ! ---> CALCUL DE L'ESTIMATEUR CORRECTION : DIVERGENCE DE ROM * U (N + 1)
  !                                          - GAMMA

  if (iescal(iescor).gt.0) then
    init = 1
    call divmas(ncelet, ncel, nfac, nfabor, init, nfecra,         &
    !==========
                ifacel, ifabor, esflum, esflub, w1)

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
      rovolsdt = propce(iel,ipcrom)*volume(iel)/dt(iel)
      do isou = 1, 3
        trav(isou,iel) = rovolsdt *                               &
                 ( vela(isou,iel)- vel(isou,iel) )
      enddo
    enddo

    !   APPEL A PREDUV AVEC RTP ET RTP AU LIEU DE RTP ET RTPA
    !                  AVEC LE FLUX DE MASSE RECALCULE
    iappel = 2
    call predvv &
    !==========
 ( iappel ,                                                       &
   nvar   , nscal  , iterns ,                                     &
   ncepdc , ncetsm ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtp    , vel    , vel    ,                   &
   propce , propfa , propfb ,                                     &
   esflum , esflub ,                                              &
   tslagr , coefa  , coefb  , coefau , coefbu , cofafu , cofbfu , &
   ckupdc , smacel , frcxt  , grdphd ,                            &
   trava  , ximpa  , uvwk   , dfrcxt , dttens , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
   w1     , w7     , w8     , w9     , w10    )

  endif

endif

!===============================================================================
! 8.  TRAITEMENT DU POINT FIXE SUR LE SYSTEME VITESSE/PRESSION
!===============================================================================

if (nterup.gt.1) then
! TEST DE CONVERGENCE DE L'ALGORITHME ITERATIF
! On initialise ICVRGE a 1 et on le met a 0 si on n'a pas convergee

  icvrge = 1

  xnrmu = 0.d0
  !$omp parallel do reduction(+:xnrmu0) private(xdu, xdv, xdw)
  do iel = 1,ncel
    xdu = vel(1,iel) - uvwk(1,iel)
    xdv = vel(2,iel) - uvwk(2,iel)
    xdw = vel(3,iel) - uvwk(3,iel)
    xnrmu = xnrmu +(xdu**2 + xdv**2 + xdw**2)     &
                                * volume(iel)
  enddo
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
( ncelet , ncel   , nfac   , nfabor ,                         &
  volume , rtp(1,ipr) )
endif

! Calcul de la pression totale IPRTOT : (definie comme propriete )
! En compressible, la pression resolue est deja la pression totale

if (ippmod(icompf).lt.0) then
  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)
  do iel=1,ncel
    propce(iel,ipproc(iprtot))= rtp(iel,ipr)           &
         + ro0*( gx*(xyzcen(1,iel)-xxp0)               &
         + gy*(xyzcen(2,iel)-xyp0)                     &
         + gz*(xyzcen(3,iel)-xzp0) )                   &
         + p0 - pred0
  enddo
endif

!===============================================================================
! 9.  IMPRESSIONS
!===============================================================================

iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))
ipcrom = ipproc(irom)
ipbrom = ipprob(irom)

if (iwarni(iu).ge.1) then

  write(nfecra,2000)

  rnorm = -1.d0
  do iel = 1, ncel
    rnorm  = max(rnorm,abs(rtp(iel,ipr)))
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
    rhom = (propce(iel1,ipcrom)+propce(iel2,ipcrom))*0.5d0
    rnorm = propfa(ifac,iflmas)/(surf*rhom)
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
    rnorm = propfb(ifac,iflmab)/(surfbn(ifac)*propfb(ifac,ipbrom))
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
    rnorm = rnorm + propfb(ifac,iflmab)
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

! Free memory
deallocate(viscf, viscb)
deallocate(drtp)
deallocate(trav)
deallocate(dfrcxt)
if (allocated(grdphd)) deallocate(grdphd)
if (allocated(esflum)) deallocate(esflum, esflub)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
deallocate(w1)
deallocate(w7, w8, w9)
if (allocated(w10)) deallocate(w10)
if (allocated(mshvel)) deallocate(mshvel)
if (allocated(secvif)) deallocate(secvif, secvib)

! Interleaved values of vel and vela

!$omp parallel do
do iel = 1, ncelet
  rtp (iel,iu) = vel (1,iel)
  rtp (iel,iv) = vel (2,iel)
  rtp (iel,iw) = vel (3,iel)
  rtpa(iel,iu) = vela(1,iel)
  rtpa(iel,iv) = vela(2,iel)
  rtpa(iel,iw) = vela(3,iel)

  ! Store the diagonal part of dttens for postprocessing purpose
  if (ipucou.eq.1 .or. ncpdct.gt.0) then
    tpucou(iel,1) = dttens(1,iel)
    tpucou(iel,2) = dttens(2,iel)
    tpucou(iel,3) = dttens(3,iel)
  endif
enddo

! Free memory
!--------------
deallocate(vel)
deallocate(vela)
deallocate(coefap)

!--------
! Formats
!--------
#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VITESSE                             ',/,&
'      --------------------------                             ',/)
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

#else

 1000 format(/,                                                   &
'   ** SOLVING VELOCITY'                                       ,/,&
'      ----------------'                                       ,/)
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

#endif

!----
! End
!----

return

end subroutine
