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

subroutine resolp &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm , isostd , idtsca ,                   &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , coefap ,                                     &
   ckupdc , smacel ,                                              &
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
! frcxt(3,ncelet)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
!dfrcxt(3,ncelet)  ! tr ! <-- ! variation de force exterieure                  !
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
use pointe, only: itypfb, dttens
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
double precision coefap(nfabor)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)
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
integer          iinvpe, imaspe, indhyd, itypfl
integer          iesdep
integer          idtsca
integer          nagmax, npstmg
integer          ibsize, iesize
integer          iescap, ircflp, ischcp, isstpp, ivar, ncymxp, nitmfp
integer          nswrsp
integer          imucpp, idftnp, iswdyp

double precision residu, phydr0
double precision ardtsr, arsr  , unsara, thetap
double precision dtsrom, unsvom, romro0
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1
double precision hint, qimp, epsrsp, blencp, relaxp

double precision rvoid(1)

double precision, allocatable, dimension(:)   :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:)   :: w1, w7
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:)   :: cofafp, coefbp, cofbfp
double precision, allocatable, dimension(:)   :: velflx, velflb, dpvar
double precision, allocatable, dimension(:)   :: coefav, cofafv, coefbv, cofbfv
double precision, allocatable, dimension(:,:) :: frchy, dfrchy

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac,2))
allocate(w1(ncelet), w7(ncelet))
if (icalhy.eq.1) allocate(frchy(ndim,ncelet), dfrchy(ndim,ncelet))

! Boundary conditions for delta P
allocate(cofafp(nfabor), coefbp(nfabor), cofbfp(nfabor))

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
if (icalhy.eq.1.or.idilat.gt.1) then
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
! 2.  RESIDU DE NORMALISATION
!===============================================================================

if(irnpnw.ne.1) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      unsvom = -1.d0/volume(iel)
      trav(iel,1) = trav(iel,1)*unsvom + frcxt(1 ,iel) + dfrcxt(1 ,iel)
      trav(iel,2) = trav(iel,2)*unsvom + frcxt(2 ,iel) + dfrcxt(2 ,iel)
      trav(iel,3) = trav(iel,3)*unsvom + frcxt(3 ,iel) + dfrcxt(3 ,iel)
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
  itypfl = 1
  if (idilat.eq.4) itypfl = 0

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu  , iv  , iw  , imaspe , itypfl ,                            &
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

  ! --- Weakly compressible algorithm: semi analytic scheme
  if (idilat.eq.4) then
    do iel = 1, ncel
      w1(iel) = w1(iel)*propce(iel,ipcrom)
    enddo
  endif

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

  call prodsc(ncel,isqrt,w1,w1,rnormp)

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

do ifac = 1, nfabor
  coefap(ifac) = 0.d0
  cofafp(ifac) = 0.d0
enddo

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

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      ! Neumann Boundary Conditions
      !----------------------------

      if (idften(ipr).eq.1) then
        hint = dt(iel)/distb(ifac)
      else if (idften(ipr).eq.3) then
        hint = ( dttens(1, iel)*surfbo(1,ifac)**2              &
               + dttens(2, iel)*surfbo(2,ifac)**2              &
               + dttens(3, iel)*surfbo(3,ifac)**2              &
               ) / (surfbn(ifac)**2 * distb(ifac))
      endif

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )

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
        frchy(1,iel)  = dronm1*gx
        frchy(2,iel)  = dronm1*gy
        frchy(3,iel)  = dronm1*gz
        dfrchy(1,iel) = drom  *gx - frchy(1,iel)
        dfrchy(2,iel) = drom  *gy - frchy(2,iel)
        dfrchy(3,iel) = drom  *gz - frchy(3,iel)
      enddo

      call calhyd &
      !==========
 ( nvar   , nscal  ,                                              &
   indhyd ,                                                       &
   frchy  , dfrchy ,                                              &
   rtp(1,ipr)   , propfa(1,iflmas), propfb(1,iflmab),             &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
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
imucpp = 0

call matrix &
!==========
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp ,                                     &
   isym   , nfecra ,                                              &
   thetap , imucpp ,                                              &
   ifacel , ifabor ,                                              &
   coefb(1,iclipr) , coefb(1,iclipf) , rovsdt ,                   &
   propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  ,          &
   rvoid  , dam    , xam    )

! Strengthen the diagonal
if (idilat.eq.3) then
  do iel = 1, ncel
    dam(iel) = dam(iel) + epsdp*volume(iel)/dt(iel)
  enddo
endif

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
   frcxt  ,                                                       &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr)  ,            &
   grad   )


if (iphydr.eq.1) then
  do iel = 1, ncel
    grad(iel,1) = grad(iel,1) - frcxt(1 ,iel)
    grad(iel,2) = grad(iel,2) - frcxt(2 ,iel)
    grad(iel,3) = grad(iel,3) - frcxt(3 ,iel)
  enddo
endif

! --- Weakly compressible algorithm: semi analytic scheme
!     The RHS contains rho div(u*) and not div(rho u*)
!     so this term will be add afterwards
if (idilat.eq.4) then
  if (idtsca.eq.0) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/propce(iel,ipcrom))
      grad(iel,1) = ardtsr*grad(iel,1)
      grad(iel,2) = ardtsr*grad(iel,2)
      grad(iel,3) = ardtsr*grad(iel,3)
    enddo
  else
    do iel=1,ncel
      arsr  = arak/propce(iel,ipcrom)
      grad(iel,1) =arsr*tpucou(iel,1)*grad(iel,1)
      grad(iel,2) =arsr*tpucou(iel,2)*grad(iel,2)
      grad(iel,3) =arsr*tpucou(iel,3)*grad(iel,3)
    enddo
  endif

! Standard case
else
  if (idtsca.eq.0) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/propce(iel,ipcrom))
      grad(iel,1) = rtp(iel,iu) + ardtsr*grad(iel,1)
      grad(iel,2) = rtp(iel,iv) + ardtsr*grad(iel,2)
      grad(iel,3) = rtp(iel,iw) + ardtsr*grad(iel,3)
    enddo
  else
    do iel=1,ncel
      arsr  = arak/propce(iel,ipcrom)
      grad(iel,1) = rtp(iel,iu)+arsr*tpucou(iel,1)*grad(iel,1)
      grad(iel,2) = rtp(iel,iv)+arsr*tpucou(iel,2)*grad(iel,2)
      grad(iel,3) = rtp(iel,iw)+arsr*tpucou(iel,3)*grad(iel,3)
    enddo
  endif
endif

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(grad(1,1), grad(1,2), grad(1,3))
  !==========
endif

init   = 1
inc    = 1
if (idilat.eq.4) inc = 0
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
itypfl = 1

call inimas &
!==========
 ( nvar   , nscal  ,                                              &
   iu     , iv     , iw     , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   grad(1,1)       , grad(1,2)       , grad(1,3)       ,          &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
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
    call projts &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   dfrcxt ,                                                       &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     )
  else
    call projts &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   dfrcxt ,                                                       &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     )
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
    call itrmas &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   frcxt  ,                                                       &
   rtpa(1,ipr)  ,                                                 &
   coefa(1,iclipr) , coefb(1,iclipr) ,                            &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
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
      ! A 0 boundary coefficient cofbfp is passed to projts
      ! to cancel boundary terms
      do ifac = 1,nfabor
        cofbfp(ifac) = 0.d0
      enddo

      call projts &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
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
      tpucou(iel,1) = arak*tpucou(iel,1)
      tpucou(iel,2) = arak*tpucou(iel,2)
      tpucou(iel,3) = arak*tpucou(iel,3)
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
   frcxt  ,                                                       &
   rtpa(1,ipr)  ,                                                 &
   coefa(1,iclipr) , coefb(1,iclipr) ,                            &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
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
      ! A 0 boundary coefficient cofbfp is passed to projts
      ! to cancel boundary terms
      do ifac = 1,nfabor
        cofbfp(ifac) = 0.d0
      enddo

      call projts &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     )

    endif

! --- Correction du pas de temps
    unsara = 1.d0/arak
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
    coefap(ifac) = 0.d0
    cofafp(ifac) = 0.d0
  enddo

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
    do ifac=1,nfabor
      if (isostd(ifac).eq.1) then
        iel=ifabor(ifac)
        coefap(ifac) = rtp(iel,ipr)                          &
             +(cdgfbo(1,ifac)-xyzcen(1,iel))*dfrcxt(1 ,iel)  &
             +(cdgfbo(2,ifac)-xyzcen(2,iel))*dfrcxt(2 ,iel)  &
             +(cdgfbo(3,ifac)-xyzcen(3,iel))*dfrcxt(3 ,iel)  &
             - phydr0
        if (idften(ipr).eq.1) then
          hint = dt(iel)/distb(ifac)
        else if (idften(ipr).eq.3) then
          hint = ( dttens(1, iel)*surfbo(1,ifac)**2              &
                 + dttens(2, iel)*surfbo(2,ifac)**2              &
                 + dttens(3, iel)*surfbo(3,ifac)**2              &
                 ) / (surfbn(ifac)**2 * distb(ifac))
        endif
        cofafp(ifac) = - hint*coefap(ifac)
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
  nagmax = nagmx0(ipr)
  npstmg = ncpmgr(ipr)
  lchain = 16

  call clmlga &
  !==========
 ( chaine(1:16) ,    lchain ,                                     &
   isym   , ibsize , iesize , nagmax , npstmg , iwarnp ,          &
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

! --- Initial divergence
init = 1

call divmas &
!==========
 ( ncelet , ncel   , nfac  , nfabor , init   , nfecra ,          &
   ifacel , ifabor ,                                             &
   propfa(1,iflmas), propfb(1,iflmab)        , w7 )

! --- Weakly compressible algorithm: semi analytic scheme
!     1. The RHS contains rho div(u*) and not div(rho u*)
!     2. Add dilatation source term to rhs
!     3. The mass flux is completed by rho u* . S
if (idilat.eq.4) then

  allocate(velflx(nfac), velflb(nfabor))

  ! 1. The RHS contains rho div(u*) and not div(rho u*)
  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  imaspe = 1
  itypfl = 0

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu     , iv     , iw     , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iu)       , rtp(1,iv)       , rtp(1,iw)       ,          &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
   velflx , velflb )

  call divmas &
  !==========
 ( ncelet , ncel   , nfac  , nfabor , init   , nfecra ,           &
   ifacel , ifabor ,                                              &
   velflx , velflb , w1 )

  do iel = 1, ncel
    w7(iel) = w7(iel) + w1(iel)*propce(iel,ipcrom)
  enddo

  ! 2. Add dilatation source term
  do iel = 1, ncel
    w7(iel) = w7(iel) + propce(iel,ipproc(iustdy(itsrho)))
  enddo

  ! 3. The mass flux is completed by rho u* . S
  init   = 0
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu     , iv     , iw     , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iu)       , rtp(1,iv)       , rtp(1,iw)       ,          &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
   propfa(1,iflmas), propfb(1,iflmab) )

endif

! --- Termes sources de masse
if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    w7(iel) = w7(iel) - volume(iel)*smacel(ii,ipr)
  enddo
endif

! --- Source term associated to the mass aggregation
if ((idilat.eq.2.or.idilat.eq.3).and. &
    (ntcabs.gt.1).and.(isuite.gt.0)) then
  do iel = 1, ncel
    drom = propce(iel,ipcrom) - propce(iel,ipcroa)
    w7(iel) = w7(iel) + drom*volume(iel)/dt(iel)
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

! --- Rajout de eps*pressure*volume/dt dans le second membre
!     pour le renforcement de la diagonale de l'algorithme
!     bas Mach.

  if (idilat.eq.3) then
    do iel = 1, ncel
      smbr(iel) = smbr(iel) - epsdp*volume(iel)/dt(iel)*rtp(iel,ipr)
    enddo
  endif

! --- Test de convergence du calcul

  call prodsc(ncel,isqrt,smbr,smbr,residu)
  if (iwarni(ipr).ge.2) then
     chaine = nomvar(ipp)
     write(nfecra,1400)chaine(1:16),isweep,residu
  endif
  if (isweep.eq.1) rnsmbr(ipp) = residu

  if (residu .le. epsrsm(ipr)*rnormp ) then
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
   dfrcxt ,                                                       &
   rtp(1,ipr)   ,                                                 &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
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
   dfrcxt ,                                                       &
   rtp(1,ipr)   ,                                                 &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
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
  iesize = 1

! ---> TRAITEMENT PERIODICITE
!     (La pression est un scalaire,
!                 pas de pb pour la rotation: IINVPE=1)
  iinvpe = 1

  call invers &
  !==========
 ( chaine(1:16)    , isym   , ibsize , iesize ,                   &
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
   dfrcxt ,                                                       &
   rtp(1,ipr)   ,                                                 &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa(1,iclipr) , coefb(1,iclipr) ,                            &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
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
   dfrcxt ,                                                       &
   rtp(1,ipr)   ,                                                 &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab))

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefa(1,iclipr) , coefb(1,iclipr) ,                            &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
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
   dfrcxt ,                                                       &
   rtp(1,ipr)      ,                                              &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
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
   dfrcxt ,                                                       &
   rtp(1,ipr)      ,                                              &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
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

!===============================================================================
! 7.  SUPPRESSION DE LA HIERARCHIE DE MAILLAGES
!===============================================================================

if (imgr(ipr).gt.0) then
  chaine = nomvar(ipp)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

!===============================================================================
! 8. Weakly compressible algorithm: semi analytic scheme
!    2nd step solving a convection diffusion equation
!===============================================================================

if (idilat.eq.4) then

  ! Allocate temporary arrays
  allocate(dpvar(ncelet))
  allocate(coefav(nfabor), coefbv(nfabor))
  allocate(cofafv(nfabor), cofbfv(nfabor))

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
    coefap(ifac) = propfb(ifac,ipbrom)
    coefbp(ifac) = 0.d0
  enddo

  call grdcel &
  !==========
  (ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), coefap , coefbp , grad   )

  ! --- dt/rho * grad rho
  do ii = 1, 3
    do iel = 1, ncel
      grad(iel,ii) = grad(iel,ii) * dt(iel) / propce(iel,ipcrom)
    enddo
  enddo

  ! --- (dt/rho * grad rho) . S

  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  extrap = extrag(iu)

  imaspe = 1
  itypfl = 0

  ! Boundary Conditions for the convective flux
  do ifac = 1, nfabor

     iel = ifabor(ifac)

     qimp = 0.d0
     hint = dt(iel)/distb(ifac)

     call set_neumann_scalar &
          !=================
        ( coefav(ifac), cofafv(ifac),             &
          coefbv(ifac), cofbfv(ifac),             &
          qimp        , hint        )

  enddo

  call inimas &
  !==========
   ( nvar   , nscal  ,                                              &
     ipcrom , ipcrom , ipcrom , imaspe , itypfl ,                   &
     iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
     iwarnp , nfecra ,                                              &
     epsrgp , climgp , extrap ,                                     &
     propce(1,ipcrom), propfb(1,ipbrom),                            &
     grad(1,1)       , grad(1,2)   , grad(1,3)   ,                  &
     coefav          , coefav      , coefav      ,                  &
     coefbv          , coefbv      , coefbv      ,                  &
     velflx , velflb )

  ! --- Viscosity
  call viscfa (imvisf, dt, viscf, viscb)

  ! --- Boundary condition for the pressure increment
  do ifac = 1, nfabor
   coefap(ifac) = 0.d0
   cofafp(ifac) = 0.d0
  enddo

  ! --- Convective source term
  do iel = 1, ncel
    smbr(iel) = 0.d0
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

  call bilsca &
  !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp , imucpp , idftnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,ipr)      , rtp(1,ipr)      ,                            &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   velflx , velflb , viscf  , viscb  , rvoid  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   smbr   )

  ! --- Initialization of the variable to solve
  do iel = 1, ncel
    rovsdt(iel) = 340.d0/dt(iel) * volume(iel)
    drtp(iel)   = 0.d0
    dpvar(iel)  = 0.d0
    smbr(iel)   = - smbr(iel)
  enddo

  ! --- Solve the convection diffusion equation

  ivar   = ipr
  iconvp = 1
  idiffp = 1
  ireslp = 1
  ipol   = 0
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

  ! --- Solve the convection diffusion equation

  call codits &
  !==========
   ( nvar   , nscal  ,                                              &
     idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
     imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
     ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
     imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
     blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
     relaxp , thetap ,                                              &
     drtp   , drtp   ,                                              &
     coefap , coefb(1,iclipr) ,                                     &
     cofafp , coefb(1,iclipf) ,                                     &
     velflx , velflb ,                                              &
     viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
     rvoid  , rvoid  ,                                              &
     rovsdt , smbr   , drtp   , dpvar  ,                            &
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

  if (idtsca.eq.0) then
    call itrmas &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab))

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrp = 0
    inc = 0

    call itrmas &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   dpvar  ,                                                       &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
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
   dfrcxt ,                                                       &
   drtp   ,                                                       &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
   propfa(1,iflmas), propfb(1,iflmab))

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrp = 0
    inc = 0

    call itrmas &
    !==========
 ( nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   dfrcxt ,                                                       &
   dpvar  ,                                                       &
   coefap , coefb(1,iclipr) ,                                     &
   cofafp , coefb(1,iclipf) ,                                     &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
   propfa(1,iflmas), propfb(1,iflmab))

  endif

  ! Free memory
  deallocate(dpvar)
  deallocate(velflx, velflb)
  deallocate(coefav, coefbv)
  deallocate(cofafv, cofbfv)

endif

!===============================================================================
! 9. Update the pressure field
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
deallocate(grad)
deallocate(dam, xam)
deallocate(w1, w7)
if (icalhy.eq.1) deallocate(frchy, dfrchy)
deallocate(cofafp, coefbp, cofbfp)

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
