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

!> \file resopv.f90
!>
!> \brief This subroutine perform the pressure correction step of the Navier
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
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icetsm        index of cells with mass source terms
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] propfa        physical properties at interior face centers
!> \param[in,out] propfb        physical properties at boundary face centers
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)


!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     isostd        indicator of standard outlet and index
!>                               of the reference outlet face
!> \param[in]     idtsca        indicator of non scalar time step
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] propfa        physical properties at interior face centers
!> \param[in,out] propfb        physical properties at boundary face centers
!> \param[in]     coefa, coefb  boundary conditions
!>
!> \param[in]     ckupdc        work array for the head loss
!>
!> \param[in]     smacel        values of variables associated to the mass
!>                               source (for the pressure, smacel
!>                               is the mass flux)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     dfrcxt        variation of the external forces
!> \param[in]                    making the hydrostatic pressure
!> \param[in]     tpucou        non scalar time step in case of
!>                               velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                               the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     viscfi        idem viscf pour increments
!> \param[in]     viscbi        idem viscb pour increments
!> \param[in]     drtp          tableau de travail pour increment
!> \param[in]     tslagr        coupling term for teh Lagrangian module
!> \param[in]     frchy         tableau de travail
!> \param[in]     dfrchy        tableau de travail
!> \param[in]     trava         tableau de travail pour couplage
!_______________________________________________________________________________

subroutine resopv &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm , isostd , idtsca ,                   &
   dt     , rtp    , rtpa   , vel    , vela   ,                   &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  , coefav , coefbv ,                            &
   ckupdc , smacel ,                                              &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  , viscfi , viscbi ,                            &
   drtp   , tslagr ,                                              &
   frchy  , dfrchy , trava  )

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
use pointe, only: itypfb
use albase
use parall
use period
use mltgrd
use lagpar
use lagran
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          isostd(nfabor+1)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(ncelet,3), dfrcxt(ncelet,3)
double precision tpucou(ndim,ncelet), trav(3,ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision drtp(ncelet)
double precision tslagr(ncelet,*)
double precision frchy(ncelet,ndim), dfrchy(ncelet,ndim)
double precision trava(ndim,ncelet)
double precision coefav(3  ,ndimfb)
double precision coefbv(3,3,ndimfb)
double precision vel   (3  ,ncelet)
double precision vela  (3  ,ncelet)

! Local variables

character*80     chaine
integer          lchain
integer          iccocg, inc   , init  , isym  , ipol  , isqrt
integer          ii, iel   , ifac  , ifac0 , iel0
integer          ireslp, nswmpr
integer          isweep, niterf, icycle
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          iclipf
integer                  iclipr, icliup, iclivp, icliwp
integer          ipcrom, ipcroa, ipbrom, iflmas, iflmab
integer          ipp
integer          idiffp, iconvp, ndircp
integer          nitmap, imgrp , ncymap, nitmgp
integer          iinvpe, indhyd
integer          iesdep
integer          idtsca
integer          nagmax, npstmg
integer          isou  , ibsize
double precision residu, resold, phydr0
double precision ardtsr, arsr  , unsara, thetap
double precision dtsrom, unsvom, romro0
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1, tcrite, relaxp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: res, divu, presa
double precision, dimension(:,:), allocatable :: gradp
double precision, allocatable, dimension(:) :: rhs, rovsdt

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac,2))
allocate(res(ncelet), presa(ncelet), divu(ncelet))
allocate(rhs(ncelet), rovsdt(ncelet))

! --- Writting
ipp    = ipprtp(ipr)

! --- Boundary conditions
!     (ICLRTP(IPR,ICOEFF) pointe vers ICLRTP(IPR,ICOEF) si IPHYDR=0)
iclipr = iclrtp(ipr,icoef)
iclipf = iclrtp(ipr,icoeff)
icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

! --- Physical quantities
ipcrom = ipproc(irom)
if (icalhy.eq.1.or.idilat.gt.1) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif
ipbrom = ipprob(irom  )
iflmas = ipprof(ifluma(ipr))
iflmab = ipprob(ifluma(ipr))

! --- Solving options
isym  = 1
if( iconv (ipr).gt.0 ) then
  isym  = 2
endif

if (iresol(ipr).eq.-1) then
  ireslp = 0
  ipol   = 0
  if( iconv(ipr).gt.0 ) then
    ireslp = 1
    ipol   = 0
  endif
else
  ireslp = mod(iresol(ipr),1000)
  ipol   = (iresol(ipr)-ireslp)/1000
endif

isqrt = 1

!===============================================================================
! 2.  RESIDU DE NORMALISATION
!===============================================================================

if(irnpnw.ne.1) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      unsvom = -1.d0/volume(iel)
      trav(1,iel) = trav(1,iel)*unsvom + frcxt(iel,1) + dfrcxt(iel,1)
      trav(2,iel) = trav(2,iel)*unsvom + frcxt(iel,2) + dfrcxt(iel,2)
      trav(3,iel) = trav(3,iel)*unsvom + frcxt(iel,3) + dfrcxt(iel,3)
    enddo
  else
    if(isno2t.gt.0) then
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = propce(iel,ipcrom)-ro0
        trav(1,iel) = (trav(1,iel)+trava(1,iel))*unsvom + romro0*gx
        trav(2,iel) = (trav(2,iel)+trava(2,iel))*unsvom + romro0*gy
        trav(3,iel) = (trav(3,iel)+trava(3,iel))*unsvom + romro0*gz
      enddo
    else
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = propce(iel,ipcrom)-ro0
        trav(1,iel) = trav(1,iel)*unsvom + romro0*gx
        trav(2,iel) = trav(2,iel)*unsvom + romro0*gy
        trav(3,iel) = trav(3,iel)*unsvom + romro0*gz
      enddo
    endif
  endif
  do iel = 1, ncel
    dtsrom = dt(iel)/propce(iel,ipcrom)
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
  nswrgp  = 1
  imligp = imligr(iu )
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  extrap = extrag(iu )

  call inimav                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iu     ,                                                       &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) )

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
       ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),res)

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
    chaine = nomvar(ipp)
    write(nfecra,1300)chaine(1:16) ,rnormp
  endif
  dervar(ipp) = rnormp
  nbivar(ipp) = 0

else

  if(iwarni(ipr).ge.2) then
    chaine = nomvar(ipp)
    write(nfecra,1300)chaine(1:16) ,rnormp
  endif
  dervar(ipp) = rnormp
  nbivar(ipp) = 0

endif

!===============================================================================
! 3.  CALCUL DE L'INCREMENT DE PRESSION HYDROSTATIQUE (SI NECESSAIRE)
!===============================================================================

if (iphydr.eq.1) then
!     L'INCREMENT EST STOCKE PROVISOIREMENT DANS RTP(.,IPRIPH)
!     on resout une equation de Poisson avec des conditions de
!     flux nul partout
!     Ce n'est utile que si on a des faces de sortie
  ifcsor = isostd(nfabor+1)
  if (irangp.ge.0) then
    call parcmx (ifcsor)
  endif

  if (ifcsor.le.0) then
    indhyd = 0
  else
    do ifac=1,nfabor
      coefa(ifac,iclipf) = 0.d0
      coefb(ifac,iclipf) = 1.d0
    enddo

    if (icalhy.eq.1) then


!     Il serait necessaire de communiquer pour periodicite et parallelisme
!      sur le vecteur dfrchy(IEL,1) dfrchy(IEL,2) dfrchy(IEL,3)
!     On peut economiser la communication tant que dfrchy ne depend que de
!      RHO et RHO n-1 qui ont ete communiques auparavant.
!     Exceptionnellement, on fait donc le calcul sur ncelet.
      do iel = 1, ncelet
        dronm1 = (propce(iel,ipcroa)-ro0)
        drom   = (propce(iel,ipcrom)-ro0)
        frchy(iel,1)  = dronm1*gx
        frchy(iel,2)  = dronm1*gy
        frchy(iel,3)  = dronm1*gz
        dfrchy(iel,1) = drom  *gx - frchy(iel,1)
        dfrchy(iel,2) = drom  *gy - frchy(iel,2)
        dfrchy(iel,3) = drom  *gz - frchy(iel,3)
      enddo

      call calhyd &
      !==========
 ( nvar   , nscal  ,                                              &
   indhyd ,                                                       &
   frchy (1,1) , frchy (1,2) , frchy (1,3) ,                      &
   dfrchy(1,1) , dfrchy(1,2) , dfrchy(1,3) ,                      &
   rtp(1,ipr)   , propfa(1,iflmas), propfb(1,iflmab),             &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , rhs    )
    else
      indhyd = 0
    endif

  endif
endif


!===============================================================================
! 4. Building of the linear system to solve
!===============================================================================

! ---> TERME INSTATIONNAIRE

do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! ---> "VITESSE" DE DIFFUSION FACETTE

if( idiff(ipr).ge. 1 ) then
  if (idtsca.eq.0) then
    call viscfa                                                   &
    !==========
 ( imvisf ,                                                       &
   dt     ,                                                       &
   viscf  , viscb  )
  else
    !interleaved version of visort
    call viortv                                                   &
    !==========
 ( imvisf ,                                                       &
   tpucou ,                                                       &
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

call matrix &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp ,                                     &
   isym   , nfecra ,                                              &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefb(1,iclipr) , rovsdt ,                                     &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   dam    , xam    )

! Strengthen the diagonal
if (idilat.eq.3) then
  do iel = 1, ncel
    dam(iel) = dam(iel) + epsdp*volume(iel)/dt(iel)
  enddo
endif

!===============================================================================
! 5. Mass flux initialization
!===============================================================================

! --- Flux de masse predit et premiere composante Rhie et Chow

! On annule la viscosite facette pour les faces couplees pour ne pas modifier
! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
! de pression et le filtre sont annules.
if (nbrcpl.ge.1) then
  do ifac = 1, nfabor
    if (ifaccp.eq.1.and.itypfb(ifac).eq.icscpl) then
      viscb(ifac) = 0.d0
    endif
  enddo
endif

! Allocate a work array for the gradient calculation
allocate(gradp(ncelet,3))

iccocg = 1
inc    = 1
nswrgp = nswrgr(ipr)
imligp = imligr(ipr)
iwarnp = iwarni(ipr)
epsrgp = epsrgr(ipr)
climgp = climgr(ipr)
extrap = extrag(ipr)

call grdpot &
!==========
 ( ipr , imrgra , inc    , iccocg , nswrgp , imligp , iphydr ,    &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr)  ,            &
   gradp  )

do iel = 1, ncelet
  do isou = 1, 3
    trav(isou,iel) = gradp(iel,isou)
  enddo
enddo

! Free memory
deallocate(gradp)

if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(1,iel) = trav(1,iel) - frcxt(iel,1)
    trav(2,iel) = trav(2,iel) - frcxt(iel,2)
    trav(3,iel) = trav(3,iel) - frcxt(iel,3)
  enddo
endif

if (idtsca.eq.0) then
  do iel = 1, ncel
    ardtsr  = arak*(dt(iel)/propce(iel,ipcrom))
    do isou = 1, 3
      trav(isou,iel) = vel(isou,iel) + ardtsr*trav(isou,iel)
    enddo
  enddo
else
  do iel=1,ncel
    arsr  = arak/propce(iel,ipcrom)
    do isou = 1, 3
      trav(isou,iel) = vel(isou,iel) + arsr*tpucou(isou,iel)*trav(isou,iel)
    enddo
  enddo
endif

! ---> Traitement du parallelisme et de la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(trav)
endif

init   = 1
inc    = 1
iflmb0 = 1
if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
nswrgp = nswrgr(iu )
imligp = imligr(iu )
iwarnp = iwarni(ipr)
epsrgp = epsrgr(iu )
climgp = climgr(iu )
extrap = extrag(iu )

call inimav                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   iu     ,                                                       &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) )


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

  if (idtsca.eq.0) then
    call projts                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   coefb(1,iclipr) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     )
  else
    call projtv                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , nswrgp , imligp ,                   &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   coefb(1,iclipr) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou )
  endif
endif

init   = 0
inc    = 1
iccocg = 1

if(arak.gt.0.d0) then

! --- Prise en compte de Arak : la viscosite face est multipliee
!       Le pas de temps aussi. On retablit plus bas.
  do ifac = 1, nfac
    viscf(ifac) = arak*viscf(ifac)
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = arak*viscb(ifac)
  enddo

  if (idtsca.eq.0) then
    do iel = 1, ncel
      dt(iel) = arak*dt(iel)
    enddo

    nswrgp = nswrgr(ipr )
    imligp = imligr(ipr )
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr )
    climgp = climgr(ipr )
    extrap = extrag(ipr )

    call itrmas                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab))

!     Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      iccocg = 1
      nswrgp = nswrgr(ipr)
      imligp = imligr(ipr)
      iwarnp = iwarni(ipr)
      epsrgp = epsrgr(ipr)
      climgp = climgr(ipr)
!     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
      do ifac = 1,nfabor
        coefb(ifac,iclipf) = 1.d0
      enddo

      call projts                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     )

    endif

! --- Correction du pas de temps
    unsara = 1.d0/arak
    do iel = 1, ncel
      dt(iel) = dt(iel)*unsara
    enddo

  else

    do iel = 1, ncel
      tpucou(1,iel) = arak*tpucou(1,iel)
      tpucou(2,iel) = arak*tpucou(2,iel)
      tpucou(3,iel) = arak*tpucou(3,iel)
    enddo

    nswrgp = nswrgr(ipr )
    imligp = imligr(ipr )
    iwarnp = iwarni(ipr )
    epsrgp = epsrgr(ipr )
    climgp = climgr(ipr )
    extrap = extrag(ipr )

    call itrmav                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   propfa(1,iflmas), propfb(1,iflmab) )

!     Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      nswrgp = nswrgr(ipr)
      imligp = imligr(ipr)
      iwarnp = iwarni(ipr)
      epsrgp = epsrgr(ipr)
      climgp = climgr(ipr)
!     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
      do ifac = 1,nfabor
        coefb(ifac,iclipf) = 1.d0
      enddo

      call projtv                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , nswrgp , imligp ,                   &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou )

    endif

! --- Correction du pas de temps
    unsara = 1.d0/arak
    do iel = 1, ncel
      tpucou(1,iel) = unsara*tpucou(1,iel)
      tpucou(2,iel) = unsara*tpucou(2,iel)
      tpucou(3,iel) = unsara*tpucou(3,iel)
    enddo

  endif

! --- Correction de la viscosite aux faces
  do ifac = 1, nfac
    viscf(ifac) = viscf(ifac)*unsara
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = viscb(ifac)*unsara
  enddo

endif

!     Calcul des CL pour l'increment de pression
!     On commence par affecter les CL classiques
!     (COEFA=0 et COEFB=COEFB(P), puis on change
!     les CL de sortie en mettant COEFA a l'increment
!     de pression hydrostatique, decale pour valoir 0
!     sur la face de reference
if (iphydr.eq.1) then
  do ifac=1,nfabor
    coefa(ifac,iclipf) = 0.d0
  enddo
  if (indhyd.eq.1) then
    ifac0 = isostd(nfabor+1)
    if (ifac0.le.0) then
      phydr0 = 0.d0
    else
      iel0 = ifabor(ifac0)
      phydr0 = rtp(iel0,ipr)                                &
           +(cdgfbo(1,ifac0)-xyzcen(1,iel0))*dfrcxt(iel0,1) &
           +(cdgfbo(2,ifac0)-xyzcen(2,iel0))*dfrcxt(iel0,2) &
           +(cdgfbo(3,ifac0)-xyzcen(3,iel0))*dfrcxt(iel0,3)
    endif
    if (irangp.ge.0) then
      call parsom (phydr0)
    endif
    do ifac=1,nfabor
      if (isostd(ifac).eq.1) then
        iel=ifabor(ifac)
        coefa(ifac,iclipf) = rtp(iel,ipr)                   &
             +(cdgfbo(1,ifac)-xyzcen(1,iel))*dfrcxt(iel,1)  &
             +(cdgfbo(2,ifac)-xyzcen(2,iel))*dfrcxt(iel,2)  &
             +(cdgfbo(3,ifac)-xyzcen(3,iel))*dfrcxt(iel,3)  &
             - phydr0
      endif
    enddo
  endif
endif


!===============================================================================
! 6. Preparation of the Algebraic Multigrid
!===============================================================================

if (imgr(ipr).gt.0) then

  ! --- Building of the mesh hierarchy

  chaine = nomvar(ipp)
  iwarnp = iwarni(ipr)
  nagmax = nagmx0(ipr)
  npstmg = ncpmgr(ipr)
  lchain = 16

  call clmlga                                                     &
  !==========
 ( chaine(1:16) ,   lchain ,                                      &
   ncelet , ncel   , nfac   ,                                     &
   isym   , nagmax , npstmg , iwarnp ,                            &
   ngrmax , ncegrm ,                                              &
   rlxp1  ,                                                       &
   dam    , xam    )

endif

!===============================================================================
! 7. Solving (Loop over the non-orthogonalities)
!===============================================================================

! --- Numbre of sweeps
nswmpr = nswrsm(ipr)

! --- Variables are set to 0
!       rtp(.,IPR) is the increment of the pressure
!       drtp       is the increment of the increment between sweeps
!       divu       is the initial divergence of the predicted mass flux
do iel = 1,ncel
  rtp(iel,ipr) = 0.d0
  drtp(iel) = 0.d0
  presa(iel) = 0.d0
enddo

relaxp = relaxv(ipr)

! --- Divergence initiale
init = 1
call divmas                                                       &
  (ncelet,ncel,nfac,nfabor,init,nfecra,                           &
   ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),divu)

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
    drom = propce(iel,ipcrom) - propce(iel,ipcroa)
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
  chaine = nomvar(ipp)
  if (rnormp.gt.0.d0) then
    write(nfecra,1440)chaine(1:16),isweep,residu/rnormp, relaxp
  else
    write(nfecra,1440)chaine(1:16),isweep,residu, relaxp
  endif
endif

! Dynamic relaxation criterion
! (Test to modify if needed: must be scticter than
! the test in the conjugate gradient)
if (swpdyn.eq.1) then
  tcrite = 100.d0*epsilo(ipr)*rnormp
else
  tcrite = 10.d0*epsrsm(ipr)*rnormp
endif

! Reconstruction loop (beginning)
!--------------------------------

do while (isweep.le.nswmpr.and.residu.gt.tcrite)

  ! --- Solving on the increment drtp
  do iel = 1, ncel
    drtp(iel) = 0.d0
  enddo

  chaine = nomvar(ipp)
  nitmap = nitmax(ipr)
  imgrp  = imgr  (ipr)
  ncymap = ncymax(ipr)
  nitmgp = nitmgf(ipr)
  iwarnp = iwarni(ipr)
  epsilp = epsilo(ipr)

  ! The pressure is a scalar => no problem for the periodicity of rotation
  ! (iinvpe=1)
  iinvpe = 1
  ibsize = 1

  call invers &
  !==========
 ( chaine(1:16)    , isym   , ibsize ,                            &
   ipol   , ireslp , nitmap , imgrp  ,                            &
   ncymap , nitmgp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnormp , residu ,                                     &
   dam    , xam    , rhs    , drtp   )

  ! Writing
  nbivar(ipp) = nbivar(ipp) + niterf
  if (abs(rnormp).gt.0.d0) then
    resvar(ipp) = residu/rnormp
  else
    resvar(ipp) = 0.d0
  endif

  ! Update the increment of pressure
  if (idtvar.ge.0.and.isweep.le.nswmpr.and.residu.gt.tcrite) then
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

  isweep = isweep + 1

  ! --- Update the right hand side if needed:
  !      rhs^{k+1} = - div(rho u^n) - D(dt, delta delta p^{k+1})

  if (isweep.le.nswmpr) then
    iccocg = 1
    init = 1
    inc  = 0
    if (iphydr.eq.1) inc = 1
    nswrgp = nswrgr(ipr)
    imligp = imligr(ipr)
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr)
    climgp = climgr(ipr)
    extrap = extrag(ipr)

    if (idtsca.eq.0) then

      call itrgrp &
      !==========
   ( nvar   , nscal  ,                                              &
     init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
     iwarnp , nfecra ,                                              &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
     rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
     viscf  , viscb  ,                                              &
     dt     , dt     , dt     ,                                     &
     rhs    )

    else
      !interleaved tpucou array
      call itrgrv &
      !==========
   ( nvar   , nscal  ,                                              &
     init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
     iwarnp , nfecra ,                                              &
     epsrgp , climgp , extrap ,                                     &
     dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
     rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
     viscf  , viscb  ,                                              &
     tpucou ,                                                       &
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

    ! Dynamic relaxation criterion
    if (swpdyn.eq.1) then
      if (isweep.gt.2) then

        if ((residu + 0.001d0*residu).gt.resold) then
          relaxv(ipr) = max(0.8d0*relaxp, 0.1d0)
        endif

      endif
      resold = residu
    endif


    ! Writing
    if (iwarni(ipr).ge.2) then
      chaine = nomvar(ipp)
      if (rnormp.gt.0.d0) then
        write(nfecra,1440)chaine(1:16),isweep,residu/rnormp, relaxp
      else
        write(nfecra,1440)chaine(1:16),isweep,residu, relaxp
      endif
    endif

  endif

enddo
! --- Reconstruction loop (end)

! Writing
if(iwarni(ipr).ge.2) then
  if(isweep.gt.nswmpr) then
     chaine = nomvar(ipp)
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

iccocg = 1
init = 0
inc  = 0
! In case of hydrostatic pressure, inc is set to 1 to take explicit
! boundary conditions on the pressure (coefa)
if (iphydr.eq.1) inc = 1
nswrgp = nswrgr(ipr)
imligp = imligr(ipr)
iwarnp = iwarni(ipr)
epsrgp = epsrgr(ipr)
climgp = climgr(ipr)
extrap = extrag(ipr)

if (idtsca.eq.0) then
  call itrmas &
  !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   presa  , coefa(1,iclipf) , coefb(1,iclipr) ,                   &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab))

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0

  call itrmas &
  !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp   , coefa(1,iclipr) , coefb(1,iclipr) ,                   &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab))

else
  ! tpucou array is interleaved
  call itrmav &
  !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   presa  , coefa(1,iclipf) , coefb(1,iclipr) ,                   &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   propfa(1,iflmas), propfb(1,iflmab))

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0

  call itrmav &
  !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp   , coefa(1,iclipr) , coefb(1,iclipr) ,                   &
   viscf  , viscb  ,                                              &
   tpucou ,                                                       &
   propfa(1,iflmas), propfb(1,iflmab))

endif

! Update the pressure
if (idtvar.lt.0) then
  do iel = 1, ncel
    rtp(iel,ipr) = rtpa(iel,ipr) + relaxv(ipr)*rtp(iel,ipr)
  enddo
else
  do iel = 1, ncel
    rtp(iel,ipr) = rtpa(iel,ipr) + rtp(iel,ipr)
  enddo
endif

!===============================================================================
! 8. Suppression of the mesh hierarchy
!===============================================================================

if (imgr(ipr).gt.0) then
  chaine = nomvar(ipp)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

! Free memory
deallocate(dam, xam)
deallocate(res, divu, presa)
deallocate(rhs, rovsdt)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1300 format(1X,A16,' : RESIDU DE NORMALISATION =', E14.6)
 1440 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6,&
             ', RELAXP = ',E14.6)
 1600 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION                  ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

#else

 1300 format(1X,A16,' : NORMED RESIDUALS = ', E14.6)
 1440 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6,&
             ', RELAXP = ',E14.6)
 1600 format(                                                     &
'@'                                                            ,/,&
'@ @@ WARNING: ', A16,' PRESSURE STEP '                        ,/,&
'@    ========'                                                ,/,&
'@  Maximum number of iterations ',I10   ,' reached'           ,/,&
'@'                                                              )

#endif

!----
! End
!----

return

end subroutine
