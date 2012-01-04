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

subroutine resolp &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm , isostd , idtsca ,                   &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  , viscfi , viscbi ,                            &
   drtp   , smbr   , rovsdt , tslagr ,                            &
   frchy  , dfrchy , trava  )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
! SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! isostd           ! te ! <-- ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
! idtsca           ! e  ! <-- ! indicateur de pas de temps non scalai          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! frcxt(ncelet,3)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
!dfrcxt(ncelet,3)  ! tr ! <-- ! variation de force exterieure                  !
!                  !    !     !  generant lapression hydrostatique             !
! tpucou           ! tr ! <-- ! couplage vitesse pression                      !
! (ncelel,ndim)    !    !     !                                                !
! trav(ncelet,3    ! tr ! <-- ! smb pour normalisation de residu               !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! viscfi(nfac)     ! tr ! --- ! idem viscf pour increments                     !
! viscbi(nfabor    ! tr ! --- ! idem viscb pour increments                     !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!  (ncelet,*)      !    !     !   lagrangien                                   !
! frchy(ncelet     ! tr ! --- ! tableau de travail                             !
!  ndim  )         !    !     !                                                !
! dfrchy(ncelet    ! tr ! --- ! tableau de travail                             !
!  ndim  )         !    !     !                                                !
! trava            ! tr ! <-- ! tableau de travail pour couplage               !
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
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)
double precision frchy(ncelet,ndim), dfrchy(ncelet,ndim)
double precision trava(ncelet,ndim)

! Local variables

character*80     chaine
integer          lchain
integer          iccocg, inc   , init  , isym  , ipol  , isqrt
integer          ii, iel   , ifac  , ifac0 , iel0
integer          ireslp, nswrp , nswmpr
integer          isweep, niterf, icycle
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          iclipf
integer                  iclipr, icliup, iclivp, icliwp
integer          ipcrom, ipcroa, ipbrom, iflmas, iflmab
integer          ipp
integer          idiffp, iconvp, ndircp
integer          nitmap, imgrp , ncymap, nitmgp
integer          iinvpe, imaspe, indhyd
integer          iesdep
integer          idtsca
integer          iagmax, nagmax, npstmg
integer          ibsize

double precision residu, phydr0
double precision ardtsr, arsr  , arakph, unsara, thetap
double precision dtsrom, unsvom, romro0
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: w1, w7
double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac,2))
allocate(w1(ncelet), w7(ncelet))

! --- Memoire

! --- Impressions
ipp    = ipprtp(ipr)

! --- Conditions aux limites
!     (ICLRTP(IPR,ICOEFF) pointe vers ICLRTP(IPR,ICOEF) si IPHYDR=0)
iclipr = iclrtp(ipr,icoef)
iclipf = iclrtp(ipr,icoeff)
icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

! --- Grandeurs physiques
ipcrom = ipproc(irom  )
if(icalhy.eq.1) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif
ipbrom = ipprob(irom  )
iflmas = ipprof(ifluma(ipr))
iflmab = ipprob(ifluma(ipr))

! --- Options de resolution
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

arakph = arak

isqrt = 1

!===============================================================================
! 2.  RESIDU DE NORMALISATION
!===============================================================================

if(irnpnw.ne.1) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      unsvom = -1.d0/volume(iel)
      trav(iel,1) = trav(iel,1)*unsvom + frcxt(iel,1) + dfrcxt(iel,1)
      trav(iel,2) = trav(iel,2)*unsvom + frcxt(iel,2) + dfrcxt(iel,2)
      trav(iel,3) = trav(iel,3)*unsvom + frcxt(iel,3) + dfrcxt(iel,3)
    enddo
  else
    if(isno2t.gt.0) then
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = propce(iel,ipcrom)-ro0
        trav(iel,1) = (trav(iel,1)+trava(iel,1))*unsvom + romro0*gx
        trav(iel,2) = (trav(iel,2)+trava(iel,2))*unsvom + romro0*gy
        trav(iel,3) = (trav(iel,3)+trava(iel,3))*unsvom + romro0*gz
      enddo
    else
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = propce(iel,ipcrom)-ro0
        trav(iel,1) = trav(iel,1)*unsvom + romro0*gx
        trav(iel,2) = trav(iel,2)*unsvom + romro0*gy
        trav(iel,3) = trav(iel,3)*unsvom + romro0*gz
      enddo
    endif
  endif
  do iel = 1, ncel
    dtsrom = dt(iel)/propce(iel,ipcrom)
    trav(iel,1) = rtp(iel,iu) +dtsrom*trav(iel,1)
    trav(iel,2) = rtp(iel,iv) +dtsrom*trav(iel,2)
    trav(iel,3) = rtp(iel,iw) +dtsrom*trav(iel,3)
  enddo

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(trav(1,1), trav(1,2), trav(1,3))
    !==========
  endif


! ON NE RECONSTRUIT PAS POUR GAGNER DU TEMPS
!   EPSRGR N'EST DONC PAS UTILISE

  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrp  = 1
  imligp = imligr(iu )
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  extrap = extrag(iu )

  imaspe = 1

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu  , iv  , iw  , imaspe ,                                     &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   propfa(1,iflmas), propfb(1,iflmab) )

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
       ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

  if (ncesmp.gt.0) then
    do ii = 1, ncesmp
      iel = icetsm(ii)
      w1(iel) = w1(iel) -volume(iel)*smacel(ii,ipr)
    enddo
  endif

! ---> LAGRANGIEN : COUPLAGE RETOUR

  if (iilagr.eq.2 .and. ltsmas.eq.1) then

    do iel = 1, ncel
      w1(iel) = w1(iel) -tslagr(iel,itsmas)
    enddo

  endif

  call prodsc(ncelet,ncel,isqrt,w1,w1,rnormp)

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
!      sur le vecteur DFRCHY(IEL,1) DFRCHY(IEL,2) DFRCHY(IEL,3)
!     On peut economiser la communication tant que DFRCHY ne depend que de
!      RHO et RHO n-1 qui ont ete communiques auparavant.
!     Exceptionnellement, on fait donc le calcul sur NCELET.
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
   drtp   , smbr   )
    else
      indhyd = 0
    endif

  endif
endif


!===============================================================================
! 4.  PREPARATION DE LA MATRICE DU SYSTEME A RESOUDRE
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
    call visort                                                   &
    !==========
 ( imvisf ,                                                       &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
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
call matrix                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp ,                                     &
   isym   , nfecra ,                                              &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefb(1,iclipr) , rovsdt ,                                     &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   dam    , xam    )

!===============================================================================
! 5.  INITIALISATION DU FLUX DE MASSE
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
allocate(grad(ncelet,3))

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
   grad   )


if (iphydr.eq.1) then
  do iel = 1, ncel
    grad(iel,1) = grad(iel,1) - frcxt(iel,1)
    grad(iel,2) = grad(iel,2) - frcxt(iel,2)
    grad(iel,3) = grad(iel,3) - frcxt(iel,3)
  enddo
endif

if (idtsca.eq.0) then
  do iel = 1, ncel
    ardtsr  = arakph*(dt(iel)/propce(iel,ipcrom))
    grad(iel,1) = rtp(iel,iu) + ardtsr*grad(iel,1)
    grad(iel,2) = rtp(iel,iv) + ardtsr*grad(iel,2)
    grad(iel,3) = rtp(iel,iw) + ardtsr*grad(iel,3)
  enddo
else
  do iel=1,ncel
    arsr  = arakph/propce(iel,ipcrom)
    grad(iel,1) = rtp(iel,iu)+arsr*tpucou(iel,1)*grad(iel,1)
    grad(iel,2) = rtp(iel,iv)+arsr*tpucou(iel,2)*grad(iel,2)
    grad(iel,3) = rtp(iel,iw)+arsr*tpucou(iel,3)*grad(iel,3)
  enddo
endif

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(grad(1,1), grad(1,2), grad(1,3))
  !==========
endif


init   = 1
inc    = 1
iccocg = 1
iflmb0 = 1
if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
nswrgp = nswrgr(iu )
imligp = imligr(iu )
iwarnp = iwarni(ipr)
epsrgp = epsrgr(iu )
climgp = climgr(iu )
extrap = extrag(iu )

imaspe = 1

call inimas &
!==========
 ( nvar   , nscal  ,                                              &
   iu  , iv  , iw  , imaspe ,                                     &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   grad(1,1) , grad(1,2) , grad(1,3) ,                            &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   propfa(1,iflmas), propfb(1,iflmab) )

! Free memory
deallocate(grad)


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
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     )
  endif
endif

init   = 0
inc    = 1
iccocg = 1

if(arakph.gt.0.d0) then

! --- Prise en compte de Arak : la viscosite face est multipliee
!       Le pas de temps aussi. On retablit plus bas.
  do ifac = 1, nfac
    viscf(ifac) = arakph*viscf(ifac)
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = arakph*viscb(ifac)
  enddo

  if (idtsca.eq.0) then
    do iel = 1, ncel
      dt(iel) = arakph*dt(iel)
    enddo

    nswrgp = nswrgr(ipr )
    imligp = imligr(ipr )
    iwarnp = iwarni(ipr)
    epsrgp = epsrgr(ipr )
    climgp = climgr(ipr )
    extrap = extrag(ipr )
    call itrmas &
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
    unsara = 1.d0/arakph
    do iel = 1, ncel
      dt(iel) = dt(iel)*unsara
    enddo

  else

    do iel = 1, ncel
      tpucou(iel,1) = arakph*tpucou(iel,1)
      tpucou(iel,2) = arakph*tpucou(iel,2)
      tpucou(iel,3) = arakph*tpucou(iel,3)
    enddo

    nswrgp = nswrgr(ipr )
    imligp = imligr(ipr )
    iwarnp = iwarni(ipr )
    epsrgp = epsrgr(ipr )
    climgp = climgr(ipr )
    extrap = extrag(ipr )
    call itrmas &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
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
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     )

    endif

! --- Correction du pas de temps
    unsara = 1.d0/arakph
    do iel = 1, ncel
      tpucou(iel,1) = unsara*tpucou(iel,1)
      tpucou(iel,2) = unsara*tpucou(iel,2)
      tpucou(iel,3) = unsara*tpucou(iel,3)
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
! 6.  PREPARATION DU MULTIGRILLE ALGEBRIQUE
!===============================================================================

if (imgr(ipr).gt.0) then

! --- Creation de la hierarchie de maillages

  chaine = nomvar(ipp)
  iwarnp = iwarni(ipr)
  iagmax = iagmx0(ipr)
  nagmax = nagmx0(ipr)
  npstmg = ncpmgr(ipr)
  lchain = 16

  call clmlga                                                     &
  !==========
 ( chaine(1:16) ,   lchain ,                                      &
   ncelet , ncel   , nfac   ,                                     &
   isym   , iagmax , nagmax , npstmg , iwarnp ,                   &
   ngrmax , ncegrm ,                                              &
   rlxp1  ,                                                       &
   dam    , xam    )

endif

!===============================================================================
! 7.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
!===============================================================================

! --- Nombre de sweeps
nswmpr = nswrsm(ipr)

! --- Mise a zero des variables
!       RTP(.,IPR) sera l'increment de pression cumule
!       DRTP       sera l'increment d'increment a chaque sweep
!       W7         sera la divergence du flux de masse predit
do iel = 1,ncel
  rtp(iel,ipr) = 0.d0
  drtp(iel) = 0.d0
  smbr(iel) = 0.d0
enddo

! --- Divergence initiale
init = 1
call divmas                                                       &
  (ncelet,ncel,nfac,nfabor,init,nfecra,                           &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w7)

! --- Termes sources de masse
if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    w7(iel) = w7(iel) -volume(iel)*smacel(ii,ipr)
  enddo
endif

! ---> Termes sources Lagrangien
if (iilagr.eq.2 .and. ltsmas.eq.1) then
  do iel = 1, ncel
    w7(iel) = w7(iel) -tslagr(iel,itsmas)
  enddo
endif

! --- Boucle de reconstruction : debut
do 100 isweep = 1, nswmpr

! --- Mise a jour du second membre
!     (signe "-" a cause de celui qui est implicitement dans la matrice)
  do iel = 1, ncel
    smbr(iel) = - w7(iel) - smbr(iel)
  enddo

! --- Test de convergence du calcul

  call prodsc(ncelet,ncel,isqrt,smbr,smbr,residu)
  if (iwarni(ipr).ge.2) then
     chaine = nomvar(ipp)
     write(nfecra,1400)chaine(1:16),isweep,residu
  endif
  if (isweep.eq.1) rnsmbr(ipp) = residu

!  Test a modifier eventuellement
! (il faut qu'il soit plus strict que celui de gradco)
  if( residu .le. 10.d0*epsrsm(ipr)*rnormp ) then
! --- Si convergence, calcul de l'indicateur
!                     mise a jour du flux de masse et sortie


! --- Calcul d'indicateur, avec prise en compte
!       du volume (norme L2) ou non

    if(iescal(iesder).gt.0) then
      iesdep = ipproc(iestim(iesder))
      do iel = 1, ncel
        propce(iel,iesdep) = abs(smbr(iel))/volume(iel)
      enddo
      if(iescal(iesder).eq.2) then
        do iel = 1, ncel
          propce(iel,iesdep) =                                    &
            propce(iel,iesdep)*sqrt(volume(iel))
        enddo
      endif
    endif



    iccocg = 1
    init = 0
    inc  = 0
    if (iphydr.eq.1) inc = 1
! --- en cas de prise en compte de Phydro, on met INC=1 pour prendre en
!     compte les CL de COEFA(.,ICLIPF)
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
   rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab))

    else

      call itrmas &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

    endif

    goto 101

  endif

! --- Resolution implicite sur l'increment d'increment DRTP
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
  ibsize = 1

! ---> TRAITEMENT PERIODICITE
!     (La pression est un scalaire,
!                 pas de pb pour la rotation: IINVPE=1)
 iinvpe = 1

  call invers                                                     &
  !==========
 ( chaine(1:16)    , isym   , ibsize ,                            &
   ipol   , ireslp , nitmap , imgrp  ,                            &
   ncymap , nitmgp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnormp   , residu ,                                   &
   dam    , xam    , smbr   , drtp   )

  nbivar(ipp) = nbivar(ipp) + niterf
  if(abs(rnormp).gt.0.d0) then
    resvar(ipp) = residu/rnormp
  else
    resvar(ipp) = 0.d0
  endif

  if( isweep.eq.nswmpr ) then

! --- Si dernier sweep :
!       Calcul d'estimateur
!       Incrementation du flux de masse
!         avec reconstruction a partir de (dP)^(NSWMPR-1)
!       Puis on rajoute la correction en (d(dP))^(NSWMPR)
!         sans reconstruction pour assurer la divergence nulle


! --- Calcul d'indicateur, avec prise en compte
!       du volume (norme L2) ou non

    if(iescal(iesder).gt.0) then
      iesdep = ipproc(iestim(iesder))
      do iel = 1, ncel
        propce(iel,iesdep) = abs(smbr(iel))/volume(iel)
      enddo
      if(iescal(iesder).eq.2) then
        do iel = 1, ncel
          propce(iel,iesdep) =                                    &
            propce(iel,iesdep)*sqrt(volume(iel))
        enddo
      endif
    endif

! --- Incrementation du flux de masse et correction

    iccocg = 1
    init = 0
    inc  = 0
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
   rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

!     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
!     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
!     est donc sans importance
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp            , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

    else

      call itrmas &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

!     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
!     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
!     est donc sans importance
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp            , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

    endif

!     Mise a jour de l'increment de pression
    do iel = 1, ncel
      rtp(iel,ipr) = rtp(iel,ipr) + drtp(iel)
    enddo

  else

! --- Si ce n'est pas le dernier sweep
!       Mise a jour de l'increment de pression et calcul direct de la
!       partie en gradient d'increment de pression du second membre
!       (avec reconstruction)

    if (idtvar.ge.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = rtp(iel,ipr) + relaxv(ipr)*drtp(iel)
      enddo
    else
      do iel = 1, ncel
        rtp(iel,ipr) = rtp(iel,ipr) + drtp(iel)
      enddo
    endif

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
   dt          , dt          , dt          ,                      &
   smbr   )

    else

      call itrgrp &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipr)   , coefa(1,iclipf) , coefb(1,iclipr) ,             &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   smbr   )

    endif

  endif

 100  continue
! --- Boucle de reconstruction : fin

if(iwarni(ipr).ge.2) then
   chaine = nomvar(ipp)
   write( nfecra,1600)chaine(1:16),nswmpr
endif

 101  continue

!     REACTUALISATION DE LA PRESSION

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
! 7.  SUPPRESSION DE LA HIERARCHIE DE MAILLAGES
!===============================================================================

if (imgr(ipr).gt.0) then
  chaine = nomvar(ipp)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

! Free memory
deallocate(dam, xam)
deallocate(w1, w7)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1300 format(1X,A16,' : RESIDU DE NORMALISATION =', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
 1600 format(                                                     &
'@                                                            ',/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION                  ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

#else

 1300 format(1X,A16,' : NORMED RESIDUALS = ', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6)
 1600 format(                                                     &
'@'                                                            ,/,&
'@ @@ WARNING: ', A16,' PRESSURE STEP '                        ,/,&
'@    ========'                                                ,/,&
'@  Maximum number of iterations ',I10   ,' reached'           ,/,&
'@'                                                              )

#endif

!----
! FIN
!----

return

end subroutine
