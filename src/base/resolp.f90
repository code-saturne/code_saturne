!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine resolp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm , isostd , idtsca ,                   &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  , viscfi , viscbi ,                            &
   dam    , xam    ,                                              &
   drtp   , smbr   , rovsdt , tslagr ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , frchy  , dfrchy , coefu  , trava  , &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!  (ncelet,*)      !    !     !   lagrangien                                   !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! frchy(ncelet     ! tr ! --- ! tableau de travail                             !
!  ndim  )         !    !     !                                                !
! dfrchy(ncelet    ! tr ! --- ! tableau de travail                             !
!  ndim  )         !    !     !                                                !
! coefu(nfab,3)    ! tr ! --- ! tableau de travail                             !
! trava            ! tr ! <-- ! tableau de travail pour couplage               !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use pointe
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

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          isostd(nfabor+1)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(ncelet,3), dfrcxt(ncelet,3)
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision frchy(ncelet,ndim), dfrchy(ncelet,ndim)
double precision coefu(nfabor,3), trava(ncelet,ndim)
double precision ra(*)

! Local variables

character*80     chaine
integer          lchain
integer          idebia, idebra
integer          iccocg, inc   , init  , isym  , ipol  , isqrt
integer          ii, iel   , ifac  , ifac0 , iel0
integer          ireslp, nswrp , nswmpr
integer          isweep, niterf, icycle
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          ipriph, iuiph , iviph , iwiph ,iclipf
integer                  iclipr, icliup, iclivp, icliwp
integer          ipcrom, ipcroa, ipbrom, iflmas, iflmab
integer          ipp
integer          iismph
integer          idiffp, iconvp, ndircp
integer          nitmap, imgrp , ncymap, nitmgp
integer          iinvpe, imaspe, indhyd
integer          iesdep
integer          idtsca
integer          iagmax, nagmax, npstmg
double precision residu, phydr0
double precision ardtsr, arsr  , arakph, unsara, thetap
double precision dtsrom, unsvom, romro0, ro0iph
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1

double precision rvoid(1)

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Memoire
idebia = idbia0
idebra = idbra0

! --- Impressions
ipp    = ipprtp(ipr)

! --- Variables
ipriph = ipr
iuiph  = iu
iviph  = iv
iwiph  = iw

! --- Conditions aux limites
!     (ICLRTP(IPRIPH,ICOEFF) pointe vers ICLRTP(IPRIPH,ICOEF) si IPHYDR=0)
iclipr = iclrtp(ipriph,icoef)
iclipf = iclrtp(ipriph,icoeff)
icliup = iclrtp(iuiph ,icoef)
iclivp = iclrtp(iviph ,icoef)
icliwp = iclrtp(iwiph ,icoef)

iismph = iisymp

! --- Grandeurs physiques
ipcrom = ipproc(irom  )
if(icalhy.eq.1) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif
ipbrom = ipprob(irom  )
iflmas = ipprof(ifluma(ipriph))
iflmab = ipprob(ifluma(ipriph))

ro0iph = ro0

! --- Options de resolution
isym  = 1
if( iconv (ipriph).gt.0 ) then
  isym  = 2
endif

if (iresol(ipriph).eq.-1) then
  ireslp = 0
  ipol   = 0
  if( iconv(ipriph).gt.0 ) then
    ireslp = 1
    ipol   = 0
  endif
else
  ireslp = mod(iresol(ipriph),1000)
  ipol   = (iresol(ipriph)-ireslp)/1000
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
        romro0 = propce(iel,ipcrom)-ro0iph
        trav(iel,1) = (trav(iel,1)+trava(iel,1))*unsvom + romro0*gx
        trav(iel,2) = (trav(iel,2)+trava(iel,2))*unsvom + romro0*gy
        trav(iel,3) = (trav(iel,3)+trava(iel,3))*unsvom + romro0*gz
      enddo
    else
      do iel = 1, ncel
        unsvom = -1.d0/volume(iel)
        romro0 = propce(iel,ipcrom)-ro0iph
        trav(iel,1) = trav(iel,1)*unsvom + romro0*gx
        trav(iel,2) = trav(iel,2)*unsvom + romro0*gy
        trav(iel,3) = trav(iel,3)*unsvom + romro0*gz
      enddo
    endif
  endif
  do iel = 1, ncel
    dtsrom = dt(iel)/propce(iel,ipcrom)
    trav(iel,1) = rtp(iel,iuiph) +dtsrom*trav(iel,1)
    trav(iel,2) = rtp(iel,iviph) +dtsrom*trav(iel,2)
    trav(iel,3) = rtp(iel,iwiph) +dtsrom*trav(iel,3)
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
  imligp = imligr(iuiph )
  iwarnp = iwarni(ipriph)
  epsrgp = epsrgr(iuiph )
  climgp = climgr(iuiph )
  extrap = extrag(iuiph )

  imaspe = 1

  call inimas                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   iuiph  , iviph  , iwiph  , imaspe ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia(iismph) ,                                                   &
   ia     ,                                                       &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
       ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

  if (ncesmp.gt.0) then
    do ii = 1, ncesmp
      iel = icetsm(ii)
      w1(iel) = w1(iel) -volume(iel)*smacel(ii,ipriph)
    enddo
  endif

! ---> LAGRANGIEN : COUPLAGE RETOUR

  if (iilagr.eq.2 .and. ltsmas.eq.1) then

    do iel = 1, ncel
      w1(iel) = w1(iel) -tslagr(iel,itsmas)
    enddo

  endif

  call prodsc(ncelet,ncel,isqrt,w1,w1,rnormp)

  if(iwarni(ipriph).ge.2) then
    chaine = nomvar(ipp)
    write(nfecra,1300)chaine(1:16) ,rnormp
  endif
  dervar(ipp) = rnormp
  nbivar(ipp) = 0

else

  if(iwarni(ipriph).ge.2) then
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
        dronm1 = (propce(iel,ipcroa)-ro0iph)
        drom   = (propce(iel,ipcrom)-ro0iph)
        frchy(iel,1)  = dronm1*gx
        frchy(iel,2)  = dronm1*gy
        frchy(iel,3)  = dronm1*gz
        dfrchy(iel,1) = drom  *gx - frchy(iel,1)
        dfrchy(iel,2) = drom  *gy - frchy(iel,2)
        dfrchy(iel,3) = drom  *gz - frchy(iel,3)
      enddo

      call calhyd                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   indhyd ,                                                       &
   ia     ,                                                       &
   frchy (1,1) , frchy (1,2) , frchy (1,3) ,                      &
   dfrchy(1,1) , dfrchy(1,2) , dfrchy(1,3) ,                      &
   rtp(1,ipriph)   , propfa(1,iflmas), propfb(1,iflmab),          &
   coefa(1,iclipf) , coefb(1,iclipf) ,                            &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbr   ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , rovsdt ,                            &
   ra     )
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

if( idiff(ipriph).ge. 1 ) then
  if (idtsca.eq.0) then
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   dt     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )
  else
    call visort                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   viscf  , viscb  ,                                              &
   ra     )
  endif
else
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo
endif

iconvp = iconv (ipriph)
idiffp = idiff (ipriph)
ndircp = ndircl(ipriph)

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
    if (ifaccp.eq.1.and.ia(iitypf-1+ifac).eq.icscpl) then
      viscb(ifac) = 0.d0
    endif
  enddo
endif

iccocg = 1
inc    = 1
nswrgp = nswrgr(ipriph)
imligp = imligr(ipriph)
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(ipriph)
climgp = climgr(ipriph)
extrap = extrag(ipriph)

call grdpot                                                       &
!==========
 ( ipriph , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rvoid  ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr)  ,         &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
!        ---------   ---------   ---------
   ra     )


if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(iel,1) = trav(iel,1) - frcxt(iel,1)
    trav(iel,2) = trav(iel,2) - frcxt(iel,2)
    trav(iel,3) = trav(iel,3) - frcxt(iel,3)
  enddo
endif

if (idtsca.eq.0) then
  do iel = 1, ncel
    ardtsr  = arakph*(dt(iel)/propce(iel,ipcrom))
    trav(iel,1) = rtp(iel,iuiph) + ardtsr*trav(iel,1)
    trav(iel,2) = rtp(iel,iviph) + ardtsr*trav(iel,2)
    trav(iel,3) = rtp(iel,iwiph) + ardtsr*trav(iel,3)
  enddo
else
  do iel=1,ncel
    arsr  = arakph/propce(iel,ipcrom)
    trav(iel,1) = rtp(iel,iuiph)+arsr*tpucou(iel,1)*trav(iel,1)
    trav(iel,2) = rtp(iel,iviph)+arsr*tpucou(iel,2)*trav(iel,2)
    trav(iel,3) = rtp(iel,iwiph)+arsr*tpucou(iel,3)*trav(iel,3)
  enddo
endif

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(trav(1,1), trav(1,2), trav(1,3))
  !==========
endif


init   = 1
inc    = 1
iccocg = 1
iflmb0 = 1
if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
nswrgp = nswrgr(iuiph )
imligp = imligr(iuiph )
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(iuiph )
climgp = climgr(iuiph )
extrap = extrag(iuiph )

imaspe = 1

call inimas                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   iuiph  , iviph  , iwiph  , imaspe ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia(iismph) ,                                                   &
   ia     ,                                                       &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )


! --- Projection aux faces des forces exterieures

if (iphydr.eq.1) then
  init   = 0
  inc    = 0
  iccocg = 1
  nswrgp = nswrgr(ipriph)
  imligp = imligr(ipriph)
  iwarnp = iwarni(ipriph)
  epsrgp = epsrgr(ipriph)
  climgp = climgr(ipriph)

  if (idtsca.eq.0) then
    call projts                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   coefb(1,iclipr) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   ra     )
  else
    call projts                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   coefb(1,iclipr) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
   ra     )
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

    nswrgp = nswrgr(ipriph )
    imligp = imligr(ipriph )
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph )
    climgp = climgr(ipriph )
    extrap = extrag(ipriph )
    call itrmas                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      iccocg = 1
      nswrgp = nswrgr(ipriph)
      imligp = imligr(ipriph)
      iwarnp = iwarni(ipriph)
      epsrgp = epsrgr(ipriph)
      climgp = climgr(ipriph)
!     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
      do ifac = 1,nfabor
        coefb(ifac,iclipf) = 1.d0
      enddo

      call projts                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   ia     ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   ra     )

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

    nswrgp = nswrgr(ipriph )
    imligp = imligr(ipriph )
    iwarnp = iwarni(ipriph )
    epsrgp = epsrgr(ipriph )
    climgp = climgr(ipriph )
    extrap = extrag(ipriph )
    call itrmas                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      iccocg = 1
      nswrgp = nswrgr(ipriph)
      imligp = imligr(ipriph)
      iwarnp = iwarni(ipriph)
      epsrgp = epsrgr(ipriph)
      climgp = climgr(ipriph)
!     on passe avec un pseudo coefB=1, pour avoir 0 aux faces de bord
      do ifac = 1,nfabor
        coefb(ifac,iclipf) = 1.d0
      enddo

      call projts                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   ia     ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   coefb(1,iclipf) ,                                              &
   propfa(1,iflmas), propfb(1,iflmab) ,                           &
   viscf  , viscb  ,                                              &
   tpucou(1,1)     , tpucou(1,2)     , tpucou(1,3)     ,          &
   ra     )

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
      phydr0 = rtp(iel0,ipriph)                             &
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
        coefa(ifac,iclipf) = rtp(iel,ipriph)                &
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

if (imgr(ipriph).gt.0) then

! --- Creation de la hierarchie de maillages

  chaine = nomvar(ipp)
  iwarnp = iwarni(ipriph)
  iagmax = iagmx0(ipriph)
  nagmax = nagmx0(ipriph)
  npstmg = ncpmgr(ipriph)
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
nswmpr = nswrsm(ipriph)

! --- Mise a zero des variables
!       RTP(.,IPR) sera l'increment de pression cumule
!       DRTP       sera l'increment d'increment a chaque sweep
!       W7         sera la divergence du flux de masse predit
do iel = 1,ncel
  rtp(iel,ipriph) = 0.d0
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
    w7(iel) = w7(iel) -volume(iel)*smacel(ii,ipriph)
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
  if (iwarni(ipriph).ge.2) then
     chaine = nomvar(ipp)
     write(nfecra,1400)chaine(1:16),isweep,residu
  endif
  if (isweep.eq.1) rnsmbr(ipp) = residu

!  Test a modifier eventuellement
! (il faut qu'il soit plus strict que celui de gradco)
  if( residu .le. 10.d0*epsrsm(ipriph)*rnormp ) then
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
    nswrgp = nswrgr(ipriph)
    imligp = imligr(ipriph)
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph)
    climgp = climgr(ipriph)
    extrap = extrag(ipriph)
    if (idtsca.eq.0) then
      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt     , dt     , dt     ,                                     &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    else

      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    endif

    goto 101

  endif

! --- Resolution implicite sur l'increment d'increment DRTP
  do iel = 1, ncel
    drtp(iel) = 0.d0
  enddo

  chaine = nomvar(ipp)
  nitmap = nitmax(ipriph)
  imgrp  = imgr  (ipriph)
  ncymap = ncymax(ipriph)
  nitmgp = nitmgf(ipriph)
  iwarnp = iwarni(ipriph)
  epsilp = epsilo(ipriph)

! ---> TRAITEMENT PERIODICITE
!     (La pression est un scalaire,
!                 pas de pb pour la rotation: IINVPE=1)
 iinvpe = 1

  call invers                                                     &
  !==========
 ( chaine(1:16)    , isym   , ipol   , ireslp , nitmap , imgrp  , &
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
    nswrgp = nswrgr(ipriph)
    imligp = imligr(ipriph)
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph)
    climgp = climgr(ipriph)
    extrap = extrag(ipriph)

    if (idtsca.eq.0) then
      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
!     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
!     est donc sans importance
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp            , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    else

      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     pour le dernier increment, on ne reconstruit pas, on n'appelle donc
!     pas GRDCEL. La valeur des DFRCXT (qui doit normalement etre nul)
!     est donc sans importance
    iccocg = 0
    nswrp = 0
    inc = 0

      call itrmas                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrp  , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   drtp            , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   propfa(1,iflmas), propfb(1,iflmab),                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    endif

!     Mise a jour de l'increment de pression
    do iel = 1, ncel
      rtp(iel,ipriph) = rtp(iel,ipriph) + drtp(iel)
    enddo

  else

! --- Si ce n'est pas le dernier sweep
!       Mise a jour de l'increment de pression et calcul direct de la
!       partie en gradient d'increment de pression du second membre
!       (avec reconstruction)

    if (idtvar.ge.0) then
      do iel = 1, ncel
        rtp(iel,ipriph) = rtp(iel,ipriph)                         &
                        + relaxv(ipriph)*drtp(iel)
      enddo
    else
      do iel = 1, ncel
        rtp(iel,ipriph) = rtp(iel,ipriph)                         &
                        + drtp(iel)
      enddo
    endif

    iccocg = 1
    init = 1
    inc  = 0
    if (iphydr.eq.1) inc = 1
    nswrgp = nswrgr(ipriph)
    imligp = imligr(ipriph)
    iwarnp = iwarni(ipriph)
    epsrgp = epsrgr(ipriph)
    climgp = climgr(ipriph)
    extrap = extrag(ipriph)

    if (idtsca.eq.0) then

      call itrgrp                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   dt          , dt          , dt          ,                      &
   smbr   ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    else

      call itrgrp                                                 &
      !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   dfrcxt(1,1),dfrcxt(1,2),dfrcxt(1,3),                           &
   rtp(1,ipriph)   , coefa(1,iclipf) , coefb(1,iclipr) ,          &
   viscf  , viscb  ,                                              &
   tpucou(1,1) , tpucou(1,2) , tpucou(1,3) ,                      &
   smbr   ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    endif

  endif

 100  continue
! --- Boucle de reconstruction : fin

if(iwarni(ipriph).ge.2) then
   chaine = nomvar(ipp)
   write( nfecra,1600)chaine(1:16),nswmpr
endif

 101  continue

!     REACTUALISATION DE LA PRESSION

if (idtvar.lt.0) then
  do iel = 1, ncel
    rtp(iel,ipriph) = rtpa(iel,ipriph)                            &
         + relaxv(ipriph)*rtp(iel,ipriph)
  enddo
else
  do iel = 1, ncel
    rtp(iel,ipriph) = rtpa(iel,ipriph) + rtp(iel,ipriph)
  enddo
endif

!===============================================================================
! 7.  SUPPRESSION DE LA HIERARCHIE DE MAILLAGES
!===============================================================================

if (imgr(ipriph).gt.0) then
  chaine = nomvar(ipp)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
  !==========
endif

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
