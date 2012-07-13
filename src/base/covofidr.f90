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

subroutine covofidr &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  , itspdv ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt , rhosav)

!===============================================================================
! FONCTION :
! ----------

! Solving the advection/diffusion equation (with source terms) for a scalar
! quantity with a drift velocity over a time step

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use pointe, only: porosi
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use lagpar
use lagran
use radiat
use mesh
use parall

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal  , itspdv

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet), rhosav(ncelet)


! Local variables

character*80     chaine
integer          ivar , ivaraux
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt, iii, iiun, ibcl
integer          ivarsc
integer          iiscav, iicp
integer          iclvar, iclvaf
integer          ipcrom, ipcroa, ipcrho, ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp   , iptsca, ipcvso
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp, ipcvis
integer          iclvit, iclvib
integer          iflvit, iflvib

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp

double precision rvoid(1)

double precision normsurf, aa, nut

integer iflmb0, ii, imaspe, w12, w13, w23

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w7, w8, w9

double precision, allocatable, dimension(:) :: energy, dissip, taufpt

double precision, allocatable, dimension(:) :: grad, dbrow, viscdyn, tempf

double precision, allocatable, dimension(:,:) :: gradu, gradv, gradw, gradc

double precision, allocatable, dimension(:,:) :: drift

double precision, allocatable, dimension(:,:) :: coefa1, coefb1

double precision alpha0, beta0, gamma0, knu, lg, cuning, beta1
double precision ckol, grdpx, grdpy, grdpz, grdsn

integer          icliup, iclivp, icliwp, idimte, itenso, jj
integer          mode



!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet))
allocate(w2(ncelet))
allocate(w3(ncelet))
allocate(w7(ncelet))
allocate(w8(ncelet))
allocate(w9(ncelet))
allocate(drift(ncel,3))
allocate(coefa1(nfabor,3))
allocate(coefb1(nfabor,3))
allocate(energy(ncel))
allocate(dissip(ncel))
allocate(taufpt(ncel))
allocate(gradc(ncelet,3))

allocate(tempf(ncel))
allocate(viscdyn(ncel))
allocate(dbrow(ncel))

if (itytur.ne.3) then
   allocate(gradu(ncelet,3))
   allocate(gradv(ncelet,3))
   allocate(gradw(ncelet,3))
endif


! test sur la viscosité (constante ou non) 
if (ivivar.gt.0) then
   ipcvis = ipproc(iviscl)
else
   ipcvis = 0
endif

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
ippvar = ipprtp(ivar)


! --- Numero des conditions aux limites
iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

! --- Numero des grandeurs physiques

ipcrom = ipproc(irom)
if (iroma .gt. 0) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif

ipcvst = ipproc(ivisct)
iflmas = ipprof(ifluma(ivar))
iflmab = ipprob(ifluma(ivar))

if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

icliup = iclrtp(iu,icoef)
iclivp = iclrtp(iv,icoef)
icliwp = iclrtp(iw,icoef)

! --- Numero du terme source dans PROPCE si extrapolation
if(isso2t(iscal).gt.0) then
  iptsca = ipproc(itssca(iscal))
else
  iptsca = 0
endif

!     S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = thetav(ivar )

chaine = nomvar(ippvar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation


do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrs(iel) = 0.d0
enddo


do iel = 1, ncel
   if (iscalt.gt.0) then
      if (iscsth(iscalt).eq.-1) then
         tempf(iel) = rtpa(iel,isca(iscalt))
      else if ( iscsth(iscalt).eq.1 ) then
         tempf(iel) = rtpa(iel,isca(iscalt))
      else if ( iscsth(iscalt).eq.2 ) then
         mode = 1
         call usthht(mode,rtpa(iel,isca(iscalt)),tempf)
         !==========
         tempf(iel) = tempf(iel)+tkelvi
      endif
   else
      tempf(iel) = t0
   endif

   if(ivivar.gt.0) then
      viscdyn(iel) = propce(iel,ipproc(iviscl))
   else
      viscdyn(iel) = viscl0
   endif

enddo

do iel = 1, ncel
   !        Terme source utilisateur
   smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtpa(iel,ivar)
   !        Diagonale
   rovsdt(iel) = max(-rovsdt(iel),zero)
enddo


! --> Non stationnary term and mass aggregation term
do iel = 1, ncel
   rovsdt(iel) = rovsdt(iel)                                        &
        + istat(ivar)*propce(iel,ipcrom)*volume(iel)/dt(iel)
enddo
 

!-->------------------------------------------
!--> Implementation of the Zaichik model of 
!--> of aerosol transport
!

!-------------------------------------------------------------------------- 
! Cuningham correction for sub-micron particles
!-------------------------------------------------------------------------- 
! Knudsen number
! lg : Mean free path of air molecules

lg  = 7.0d-8
knu = 2 * lg/(diapart(iscal))

alpha0 = 1.257d0
beta0  = 0.4d0
gamma0 = 1.1d0

cuning = 1.d0 + knu *(alpha0 + beta0*exp(-gamma0/knu))

if (diapart(iscal).le.1.d-6) then
   do iel = 1, ncel
      taupae(iscal,iel) = taupae(iscal,iel) * cuning
   enddo
endif


! --> Step 1 : All zeros

do ifac = 1,nfabor

 propfb(ifac,iflmab) = 0.d0 
 viscb(ifac) = 0.d0

enddo

do ifac = 1, nfac

   viscf(ifac) = 0.d0
   propfa(ifac,iflmas) = 0.0

enddo

!-----------------------------------------------------
! Diffusion
!-----------------------------------------------------

if (idrift(iscal).gt.0) then

! --> Step 1 : Diffusion: Turbulence effect

! 
! Retrieve of the turbulent kinetic energy with
! respect to the chosen turbulence model

   
   if (itytur.eq.2 .or. iturb.eq.50) then
      do iel = 1,ncel
         energy(iel) = rtp(iel,ik)
         dissip(iel) = rtp(iel,iep)
      enddo
   else if (itytur.eq.3) then
      do iel = 1,ncel
         energy(iel) = 0.5d0*( rtp(iel,ir11)    &
              +rtp(iel,ir22)                    &
              +rtp(iel,ir33))
         dissip(iel) = rtp(iel,iep)
      enddo
   else if (iturb.eq.60) then
      do iel = 1,ncel
         energy(iel) = rtp(iel,ik)
         dissip(iel) = cmu * energy(iel) * rtp(iel,iomg)
      enddo
   else
      do iel = 1,ncel
         energy(iel) = 0.d0
         dissip(iel) = 0.d0
      enddo
   endif
! 
! Calculation of taufp, interaction time
! between vortices and particle
!
! ckol : Kolmogorov constant

   ckol   = 2.1d0

   beta1  = 0.5d0 + 3 * ckol/4.d0

   if (itytur.eq.2 .or. iturb.eq.50                    &
        .or. iturb.eq.60) then
      do iel = 1,ncel
         taufpt(iel) =                                 &
              (3/2.d0)*(cmu/sigmas(iscal))*energy(iel)/dissip(iel)
      enddo
   else if (itytur.eq.3) then
      do iel = 1,ncel
         taufpt(iel) = (energy(iel)/dissip(iel))/beta1
      enddo
   else
      do iel = 1,ncel
         taufpt(iel) = 0.d0
      enddo
   endif


!
! Retrieve the velocity gradient when a 
! second-order turbulence model (Rij) is
! not used

   if (itytur.ne.3) then
      inc    = 1
      iccocg = 1
      nswrgp = nswrgr(iu)
      imligp = imligr(iu)
      iwarnp = iwarni(iu)
      epsrgp = epsrgr(iu)
      climgp = climgr(iu)
      extrap = extrag(iu)
      !
      call grdcel &
           !==========
           ( iu  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
           iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
           rtpa(1,iu)    , coefa(1,icliup) , coefb(1,icliup) ,         &
           gradu )

      inc    = 1
      iccocg = 1
      nswrgp = nswrgr(iv)
      imligp = imligr(iv)
      iwarnp = iwarni(iv)
      epsrgp = epsrgr(iv)
      climgp = climgr(iv)
      extrap = extrag(iv)

      call grdcel &
           !==========
           ( iv  , imrgra , inc    , iccocg , nswrgp , imligp ,       &
           iwarnp , nfecra , epsrgp , climgp , extrap ,               &
           rtpa(1,iv)   , coefa(1,iclivp) , coefb(1,iclivp) ,         &
           gradv )

      inc    = 1
      iccocg = 1
      nswrgp = nswrgr(iw)
      imligp = imligr(iw)
      iwarnp = iwarni(iw)
      epsrgp = epsrgr(iw)
      climgp = climgr(iw)
      extrap = extrag(iw)

      call grdcel &
           !==========
           ( iw  , imrgra , inc    , iccocg , nswrgp , imligp ,         &
           iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
           rtpa(1,iw)   , coefa(1,icliwp) , coefb(1,icliwp) ,          &
           gradw)
   endif

!
! ---> Calculation of grad(C)
!
   ivaraux = ivar
   ivar = 2

   inc    = 1
   iccocg = 1
   nswrgp = nswrgr(ivar)
   imligp = imligr(ivar)
   iwarnp = iwarni(ivar)
   epsrgp = epsrgr(ivar)
   climgp = climgr(ivar)
   extrap = extrag(ivar)

   ivar = ivaraux

   call grdcel &
   !==========
        ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
        iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
        rtpa(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,         &
        gradc)

! ---> Calculation of de (Extradiag).grad(C)
!
! iiso = 1 : taking into account the extradiagonal 
!            terms in the Zaichik diffusion tensor
!

   if (iiso .eq. 1) then

      if (itytur.eq.3) then  ! Rij model

         do iel = 1, ncel

            if (dissip(iel).gt.0.d0) then
               nut = cmu * energy(iel)**2 / dissip(iel)
            else
               nut = 0.d0
            endif

            w12 = rhosav(iel) * taufpt(iel) * rtp(iel,ir12)
            w13 = rhosav(iel) * taufpt(iel) * rtp(iel,ir13)
            w23 = rhosav(iel) * taufpt(iel) * rtp(iel,ir23)

            w7(iel) = gradc(iel,2) * w12 + gradc(iel,3) * w13
            w8(iel) = gradc(iel,1) * w12 + gradc(iel,3) * w23
            w9(iel) = gradc(iel,1) * w13 + gradc(iel,2) * w23

         enddo

      else

         do iel = 1, ncel

            if (dissip(iel).gt.0.d0) then
               nut = cmu * energy(iel)**2 / dissip(iel)
            else
               nut = 0.d0
            endif

            w12 = rhosav(iel) * taufpt(iel) *                       &
                 (- nut *(gradu(iel,2) + gradv(iel,1)))

            w13 = rhosav(iel) * taufpt(iel) *                       &
                 (- nut *(gradu(iel,3) + gradw(iel,1)))

            w23 = rhosav(iel) * taufpt(iel) *                       &
                 (- nut *(gradv(iel,3) + gradw(iel,2)))

            w7(iel) = gradc(iel,2) * w12 + gradc(iel,3) * w13
            w8(iel) = gradc(iel,1) * w12 + gradc(iel,3) * w23
            w9(iel) = gradc(iel,1) * w13 + gradc(iel,2) * w23
            
         enddo
      endif

   else

      do iel = 1, ncel

         w7(iel) = 0.d0
         w8(iel) = 0.d0
         w9(iel) = 0.d0

      enddo

   endif


   ! ---> Assemble { E.grad(C) } .S at the (internal and boundary) faces

   call vectds                                                      &
  !!==========
  ( w7     , w8     , w9 ,                                    &
    viscf  , viscb )

   init = 1

   call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
        ifacel,ifabor,viscf,viscb,w1)

   !    If extrapolation of the source terms

   if(isto2t.gt.0) then
      do iel = 1, ncel
         propce(iel,iptsca) = propce(iel,iptsca) + w1(iel)
      enddo
      !     Sinon
   else
      do iel = 1, ncel
         smbrs(iel) = smbrs(iel) + w1(iel)
      enddo
   endif

   ! ---> Calculation of de (diag).grad(C)
   ! -->  Brownian motion effect
   !
   if (ibrow .eq. 0) then

      do iel = 1, ncel
         dbrow(iel) = 0.d0
      enddo

   else

      do iel = 1, ncel

         dbrow(iel) = kboltz * tempf(iel) * cuning /              &
              (3.d0 * pi * diapart(iscal) * viscdyn(iel))


      enddo
   endif

   if (itytur.eq.3) then

      do iel = 1, ncel
         w1(iel) = rhosav(iel) * (taufpt(iel) * rtp(iel,ir11) + dbrow(iel))
         w2(iel) = rhosav(iel) * (taufpt(iel) * rtp(iel,ir22) + dbrow(iel))
         w3(iel) = rhosav(iel) * (taufpt(iel) * rtp(iel,ir33) + dbrow(iel))
      enddo

   else

      do iel = 1, ncel

         if (dissip(iel).gt.0.d0) then
            nut = cmu * energy(iel)**2 / dissip(iel)
         else
            nut = 0.d0
         endif

         w1(iel) = rhosav(iel)*taufpt(iel)               &
              * ((2/3.d0) * energy(iel) - nut * 2 *      & 
              gradu(iel,1) + dbrow(iel))

         w2(iel) = rhosav(iel)*taufpt(iel) *              &
              ((2/3.d0) * energy(iel) - nut * 2 *         &
              gradv(iel,2) + dbrow(iel))

         w3(iel) = rhosav(iel)*taufpt(iel)               &
              * ((2/3.d0) * energy(iel) - nut * 2 *      &
              gradw(iel,3) + dbrow(iel))

         if (w1(iel).le.0.0d0) then                                       

            w1(iel) = rhosav(iel) *                       &
                 (taufpt(iel) * (2/3.d0)*energy(iel) + dbrow(iel))

         endif

         if (w2(iel).le.0.0d0) then   
                                     
            w2(iel) = rhosav(iel) *                       &
                 (taufpt(iel) * (2/3.d0)*energy(iel) + dbrow(iel))

         endif

         if (w3(iel).le.0.0d0) then                        
                
            w3(iel) = rhosav(iel)*(taufpt(iel) *          & 
                 (2/3.d0)*energy(iel) + dbrow(iel))
         endif

      enddo
   endif


   ! --->  Management of the parallelism

   if(irangp.ge.0) then
!!      call synvec (gradc)
!!$      !==========
!!$      call parcom (w1)
!!$      !==========
!!$      call parcom (w2)
!!$      !==========
!!$      call parcom (w3)
!!$      !==========
   endif

   ! -->   Management of the periodicity (to be verified)

!!$   if(iperio.eq.1) then
!!$
!!$      idimte = 1
!!$      itenso = 0
!!$
!!$      call percom                                                   &
!!$                                !==========
!!$ !          ( idimte , itenso ,                                           &
!!$ !          w4   , w4   , w4,                                           &
!!$ !          w5   , w5   , w5,                                           &
!!$ !          w6   , w6   , w6)
!!$
!!$      call percom                                                   &
!!$                                !==========
!!$           ( idimte , itenso ,                                           &
!!$           w1   , w1   , w1  ,                                         &
!!$           w2   , w2   , w2  ,                                         &
!!$           w3   , w3   , w3  )
!!$   endif

   do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      grdpx = 0.5d0 * (gradc(ii,1) + gradc(jj,1))
      grdpy = 0.5d0 * (gradc(ii,2) + gradc(ii,2))
      grdpz = 0.5d0 * (gradc(ii,3) + gradc(ii,3))

      grdsn = grdpx*surfac(1,ifac)+grdpy*surfac(2,ifac)             &
           + grdpz*surfac(3,ifac)

      grdpx = grdpx-grdsn*surfac(1,ifac)/surfan(ifac)
      grdpy = grdpy-grdsn*surfac(2,ifac)/surfan(ifac)
      grdpz = grdpz-grdsn*surfac(3,ifac)/surfan(ifac)

      viscf(ifac)= 0.5d0 * (                                          &
           (w1(ii) + w1(jj)) * grdpx * surfac(1,ifac)                &
        +  (w2(ii) + w2(jj)) * grdpy * surfac(2,ifac)                &
        +  (w3(ii) + w3(jj)) * grdpz * surfac(3,ifac))

   enddo


   do ifac = 1, nfabor
      viscb(ifac) = 0.d0
   enddo

   init = 1

   call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
        ifacel,ifabor,viscf,viscb,w7)

   !     Si on extrapole les termes sources
   if(isto2t.gt.0) then
      do iel = 1, ncel
         propce(iel,iptsca) = propce(iel,iptsca) + w7(iel)
      enddo
      !     Sinon
   else
      do iel = 1, ncel
         smbrs(iel) = smbrs(iel) + w7(iel)
      enddo
   endif

   ! ---> Orthotropic viscosity for the implicit part

   call visort &
                                !==========
        ( imvisf ,                 &
        w1     , w2     , w3  ,   &
        viscf  , viscb  )

else

! No diffusitivity
!
do ifac = 1, nfac
   viscf(ifac) = 0.d0
enddo

do ifac = 1, nfabor
   viscb(ifac) = 0.d0
enddo
!
!!$
!!$   if(ipcvsl.eq.0)then
!!$      do iel = 1, ncel
!!$         w1(iel) = visls0(iscal)                                     &
!!$          + idifft(ivar)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
!!$      enddo
!!$   else
!!$      do iel = 1, ncel
!!$         w1(iel) = propce(iel,ipcvsl)                                &
!!$           + idifft(ivar)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
!!$      enddo
!!$   endif
!!$
!!$   call viscfa &
!!$                                !==========
!!$        ( imvisf ,        &
!!$        w1     ,        &
!!$        viscf  , viscb  )

endif



!-----------------------------------------------------
! Convection
!-----------------------------------------------------

do ii = 1,3
   do iel = 1, ncel
      drift(iel,ii) = 0.d0
   enddo
enddo


! --> Step 1: Transport by the fluid velocity

do iel = 1, ncel
   drift(iel,1) = drift(iel,1) + rtp(iel,iu)
   drift(iel,2) = drift(iel,2) + rtp(iel,iv)  
   drift(iel,3) = drift(iel,3) + rtp(iel,iw)  
enddo


! --> Step 2: Sedimentation effect : V_drift = tau_p * g 

do iel = 1, ncel
   drift(iel,1) = drift(iel,1) + gx * taupae(iscal,iel)
   drift(iel,2) = drift(iel,2) + gy * taupae(iscal,iel)
   drift(iel,3) = drift(iel,3) + gz * taupae(iscal,iel)
enddo


!
! --> Step 3: Particle-trajectory deviation term
!
if (idrift(iscal).eq.2 .and. itstde .eq. 1) then

   do iel = 1, ncelet
      w1(iel) = 0.d0
      w2(iel) = 0.d0
      w3(iel) = 0.d0
      
   enddo

! 3 appels à bilsc2 : ivar remplacé par iu, puis iv puis iw

  thetap = thetav(ivar)

! 1er appel : le 2nd membre est W1

  do iel = 1, ncel
     w1(iel) =  - propce(iel,ipcrom) * volume(iel)                             &
          *(rtp(iel,iu)-rtpa(iel,iu)) / dt(iel)
  enddo

  iconvp = 1
  idiffp = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  ircflp = ircflu(iu)
  ischcp = ischcv(iu)
  isstpp = isstpc(iu)
  
  inc    = 1
  iccocg = 1
  ipp    = ipprtp(iu)
  iwarnp = iwarni(iu)
  blencp = blencv(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)
  iclvit = iclrtp(iu,icoef)
  iclvib = iclrtp(iu,icoeff)
  iflvit = ipprof(ifluma(iu))
  iflvib = ipprob(ifluma(iu))
  relaxp = relaxv(iu)
  
call bilsc2                                                       &
!==========
(  nvar   , nscal  ,                                              &
   idtvar , iu     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,iu),rtpa(1,iu),coefa(1,iclvit) , coefb(1,iclvit) ,       &
                     coefa(1,iclvib) , coefb(1,iclvib) ,          &
   propfa(1,iflvit), propfb(1,iflvib), viscf  , viscb  ,          &
   w1     )


! 2e appel : le 2nd membre est W2
!
  do iel = 1, ncel
     w2(iel) = -propce(iel,ipcrom)*volume(iel)                             &
          *(rtp(iel,iv)-rtpa(iel,iv))/dt(iel)
  enddo

  iconvp = 1
  idiffp = 0
  nswrgp = nswrgr(iv)
  imligp = imligr(iv)
  ircflp = ircflu(iv)
  ischcp = ischcv(iv)
  isstpp = isstpc(iv)

  inc    = 1
  iccocg = 1
  ipp    = ipprtp(iv)
  iwarnp = iwarni(iv)
  blencp = blencv(iv)
  epsrgp = epsrgr(iv)
  climgp = climgr(iv)
  extrap = extrag(iv)
  iclvit = iclrtp(iv,icoef)
  iclvib = iclrtp(iv,icoeff)
  iflvit = ipprof(ifluma(iv ))
  iflvib = ipprob(ifluma(iv ))
  relaxp = relaxv(iv)


call bilsc2 &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , iv     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,iv),rtpa(1,iv),coefa(1,iclvit) , coefb(1,iclvit) ,       &
                     coefa(1,iclvib) , coefb(1,iclvib) ,          &
   propfa(1,iflvit), propfb(1,iflvib), viscf  , viscb  ,          &
   w2     )

! 3e appel : le 2nd membre est W3

  do iel = 1, ncel
         w3(iel) = -propce(iel,ipcrom)*volume(iel)                             &
                 *(rtp(iel,iw)-rtpa(iel,iw))/dt(iel)
  enddo

  iconvp = 1
  idiffp = 0
  nswrgp = nswrgr(iw)
  imligp = imligr(iw)
  ircflp = ircflu(iw)
  ischcp = ischcv(iw)
  isstpp = isstpc(iw)

  inc    = 1
  iccocg = 1
  ipp    = ipprtp(iw)
  iwarnp = iwarni(iw)
  blencp = blencv(iw)
  epsrgp = epsrgr(iw)
  climgp = climgr(iw)
  extrap = extrag(iw)
  iclvit = iclrtp(iw,icoef)
  iclvib = iclrtp(iw,icoeff)
  iflvit = ipprof(ifluma(iw ))
  iflvib = ipprob(ifluma(iw ))
  relaxp = relaxv(iw)

  call bilsc2                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , iw     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,iw),rtpa(1,iw),coefa(1,iclvit) , coefb(1,iclvit) ,       &
                     coefa(1,iclvib) , coefb(1,iclvib) ,          &
   propfa(1,iflvit), propfb(1,iflvib), viscf  , viscb  ,          &
   w3     )
  

  write(*,*) "w1(1) = ",w1(1)
  write(*,*) "w2(1) = ",w2(1)
  write(*,*) "w3(1) = ",w3(1)


  do iel = 1, ncel
     drift(iel,1) = drift(iel,1)                                 &
        - (gradu(iel,1) +  gradu(iel,2) +  gradu(iel,3))         &
        + taupae(iscal,iel)*w1(iel)/(propce(iel,ipcrom)*volume(iel))     

     drift(iel,2) =  drift(iel,2)                                 &
        - (gradv(iel,1) +  gradv(iel,2) +  gradv(iel,3))         &
        + taupae(iscal,iel)*w2(iel)/(propce(iel,ipcrom)*volume(iel))   


     drift(iel,3) = drift(iel,3)                                 &
        - (gradw(iel,1) +  gradw(iel,2) +  gradw(iel,3))         &
        + taupae(iscal,iel)*w3(iel)/(propce(iel,ipcrom)*volume(iel))    

  enddo

endif












!===============================================================================
! 3. RESOLUTION
!===============================================================================



do ifac = 1, nfabor

  coefa1(ifac,1) = 0.d0
  coefa1(ifac,2) = 0.d0
  coefa1(ifac,3) = 0.d0

  coefb1(ifac,1) = 0.d0
  coefb1(ifac,2) = 0.d0
  coefb1(ifac,3) = 0.d0      

enddo 

! On ne le met pas a zero en paroi 
iflmb0 = 0
! On l'initialise a 0
init   = 1
! On prend en compte les Dirichlet
inc    = 1
! On recalcule les gradients complets
iccocg = 1
! Calcul de flux std (pas de divrij)
imaspe = 1

ivaraux = ivar
ivar = 2

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)

ivar = ivaraux

call inimas                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   ivar   , ivar   , ivar   , imaspe ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipproc(irom))  ,  propfb(1, ipprob(irom))  ,          &
   drift(1,1)   ,   drift(1,2)   ,   drift(1,3)   ,               &
   coefa1(1,1) , coefa1(1,2) ,coefa1(1,3)  ,                      &
   coefb1(1,1) , coefb1(1,2) ,coefb1(1,3)  ,                      &
   propfa(1,iflmas) ,  propfb(1,iflmab) )




call codits                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , rtp(1,ivar)     ,                            &
   rvoid  )

!===============================================================================
! 4. IMPRESSIONS ET CLIPPINGS
!===============================================================================

if(ivarsc.gt.0) then
  iii = ivarsc
else
! Valeur bidon
  iii = 1
endif

call clpsca                                                       &
!==========
 ( ncelet , ncel   , nvar   , nscal  , iscal  ,                   &
   propce , rtp(1,iii)      , rtp    )


! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

if (iwarni(ivar).ge.2) then
  if(nswrsm(ivar).gt.1) then
    ibcl = 1
  else
    ibcl = 0
  endif
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                       &
            - istat(ivar)*(rhosav(iel)/dt(iel))*volume(iel)&
                *(rtp(iel,ivar)-rtpa(iel,ivar))*ibcl
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

! Free memory
deallocate(w1)
deallocate(w2)
deallocate(w3)
deallocate(w7)
deallocate(w8)
deallocate(w9)
deallocate(drift)
deallocate(coefa1)
deallocate(coefb1)
deallocate(taufpt)
deallocate(energy)
deallocate(dissip)
deallocate(gradc)
deallocate(tempf)
deallocate(viscdyn)
deallocate(dbrow)

if (itytur.ne.3) then
   deallocate(gradu)
   deallocate(gradv)
   deallocate(gradw)
endif


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/,                                                   &
'   ** SOLVING VARIABLE ',A8                                   ,/,&
'      ----------------'                                       ,/)
 1200 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

#endif

!----
! FIN
!----

return

end subroutine
