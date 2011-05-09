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

subroutine cfqdmv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   flumas , flumab ,                                              &
   coefa  , coefb  , ckupdc , smacel , frcxt  , dfrcxt ,          &
   tpucou , trav   , viscf  , viscb  , viscfi , viscbi ,          &
   dam    , xam    ,                                              &
   drtp   , smbr   , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
! SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)

! AU PREMIER APPEL,  ON EFFECTUE LA PREDICITION DES VITESSES
!               ET  ON CALCULE UN ESTIMATEUR SUR LA VITESSE PREDITE

! AU DEUXIEME APPEL, ON CALCULE UN ESTIMATEUR GLOBAL
!               POUR NAVIER-STOKES :
!   ON UTILISE TRAV, SMBR ET LES TABLEAUX DE TRAVAIL
!   ON APPELLE BILSC2 AU LIEU DE CODITS
!   ON REMPLIT LE PROPCE ESTIMATEUR IESTOT
!   CE DEUXIEME APPEL INTERVIENT DANS NAVSTO APRES RESOLP
!   LORS DE CE DEUXIEME APPEL
!     RTPA ET RTP SONT UN UNIQUE TABLEAU (= RTP)
!     LE FLUX DE MASSE EST LE FLUX DE MASSE DEDUIT DE LA VITESSE
!      AU CENTRE CONTENUE DANS RTP

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iphas            ! i  ! <-- ! phase number                                   !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! flumas           ! tr ! <-- ! flux de masse aux faces internes               !
!  (nfac  )        !    !     !   (distinction iappel=1 ou 2)                  !
! flumab           ! tr ! <-- ! flux de masse aux faces de bord                !
!  (nfabor  )      !    !     !    (distinction iappel=1 ou 2)                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
!dfrcxt(ncelet,    ! tr ! <-- ! variation de force exterieure                  !
!   3,nphas)       !    !     !  generant lapression hydrostatique             !
! tpucou           ! tr ! --> ! couplage vitesse pression                      !
! (ncelel,ndim)    !    !     !                                                !
! trav(ncelet,3    ! tr ! --> ! smb qui servira pour normalisation             !
!                  !    !     !  dans resolp                                   !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! viscfi(nfac)     ! tr ! --- ! idem viscf pour increments                     !
! viscbi(nfabor    ! tr ! --- ! idem viscb pour increments                     !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! coefu(nfab,3)    ! tr ! --- ! tableau de travail                             !
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
use pointe
use numvar
use entsor
use cstphy
use cstnum
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision flumas(nfac), flumab(nfabor)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(ncelet,3,nphas), dfrcxt(ncelet,3,nphas)
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision coefu(nfabor,3)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          iel   , ielpdc, ifac  , ivar  , isou  , iii
integer          iccocg, inc   , init  , iphydp, ii
integer          ireslp, nswrgp, imligp, iwarnp, ipp
integer          ipriph, ikiph , iuiph , iviph , iwiph
integer          iclik , iclvar, iclvaf, iclipr
integer          ipcrom, ipcvis, ipcvst
integer          iconvp, idiffp, ndircp, nitmap, nswrsp
integer          ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          idiaex, iterns
integer          iifru
integer          maxelt, ils
double precision rnorm , vitnor
double precision romvom, rtprom
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision vit1  , vit2  , vit3  , thetap, pfac, pfac1
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , pbord , diipbx, diipby, diipbz, pip, xkb
!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Initialize variables to avoid compiler warnings

iclik = 0

! Memoire

idebia = idbia0
idebra = idbra0

ipriph = ipr
iuiph  = iu
iviph  = iv
iwiph  = iw
if(itytur.eq.2 .or. iturb.eq.50                     &
     .or. iturb.eq.60) then
  ikiph  = ik
endif

if(itytur.eq.2 .or. iturb.eq.50                     &
     .or. iturb.eq.60) then
  iclik  = iclrtp(ikiph ,icoef)
endif

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)

!     Indicateur flux de bord Rusanov
if(iifbru.gt.0) then
  iifru = iifbru+(iphas-1)*nfabor
else
  iifru = 1
endif

!===============================================================================
! 2.  GRADIENT DE PRESSION ET GRAVITE
!===============================================================================

! ---> PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE :

if (iphydr.eq.1) then

!     on doit pouvoir adapter l'option iphydr au compressible,
!       mais noter plusieurs points
!       - on dispose de la masse volumique au pas de temps courant et au
!         pas de temps précédent (au pas de temps précédent dans propce
!         en particulier)
!       - la correction de pression est ici généree par la resolution de
!         l'energie (la pression change sans que rho ne change)
!       - si l'objectif se limite a adapter le calcul de grad p pour
!         qu'il compense le rho0 g, noter quand meme que l'on ne resout
!         pas en rho-rho0 et que P est tjrs cohérent avec rho (par la
!         thermo)

  do iel = 1, ncel

! variation de force (utilise dans resolp)

    if(igrdpp.gt.0) then
      rtprom = rtp(iel,isca(irho))
    else
      rtprom = rtpa(iel,isca(irho))
    endif

    dfrcxt(iel,1,iphas) = rtprom*gx - frcxt(iel,1,iphas)
    dfrcxt(iel,2,iphas) = rtprom*gy - frcxt(iel,2,iphas)
    dfrcxt(iel,3,iphas) = rtprom*gz - frcxt(iel,3,iphas)
  enddo
!     Ajout eventuel des pertes de charges
  if (ncepdp.gt.0) then
    do ielpdc = 1, ncepdp
      iel=icepdc(ielpdc)
      vit1   = rtp(iel,iuiph)
      vit2   = rtp(iel,iviph)
      vit3   = rtp(iel,iwiph)
      cpdc11 = ckupdc(ielpdc,1)
      cpdc22 = ckupdc(ielpdc,2)
      cpdc33 = ckupdc(ielpdc,3)
      cpdc12 = ckupdc(ielpdc,4)
      cpdc13 = ckupdc(ielpdc,5)
      cpdc23 = ckupdc(ielpdc,6)
      dfrcxt(iel,1,iphas) = dfrcxt(iel,1,iphas)                   &
 -rtp(iel,isca(irho))*(cpdc11*vit1+cpdc12*vit2+cpdc13*vit3)
      dfrcxt(iel,2,iphas) = dfrcxt(iel,2,iphas)                   &
 -rtp(iel,isca(irho))*(cpdc12*vit1+cpdc22*vit2+cpdc23*vit3)
      dfrcxt(iel,3,iphas) = dfrcxt(iel,3,iphas)                   &
 -rtp(iel,isca(irho))*(cpdc13*vit1+cpdc23*vit2+cpdc33*vit3)
    enddo
  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(dfrcxt(1,1,iphas), dfrcxt(1,2,iphas), dfrcxt(1,3,iphas))
    !==========
  endif

endif

!       Fin du test sur IPHYDR


! ---> PRISE EN COMPTE DU GRADIENT DE PRESSION

iccocg = 1
inc    = 1
nswrgp = nswrgr(ipriph)
imligp = imligr(ipriph)
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(ipriph)
climgp = climgr(ipriph)
extrap = extrag(ipriph)

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ipriph , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   frcxt(1,1,iphas), frcxt(1,2,iphas), frcxt(1,3,iphas),          &
   rtp(1,ipriph)   , coefa(1,iclrtp(ipriph,icoef))  ,             &
                     coefb(1,iclrtp(ipriph,icoef))  ,             &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )


if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(iel,1) = ( frcxt(iel,1,iphas) - w1(iel) )*volume(iel)
    trav(iel,2) = ( frcxt(iel,2,iphas) - w2(iel) )*volume(iel)
    trav(iel,3) = ( frcxt(iel,3,iphas) - w3(iel) )*volume(iel)
  enddo
else
  do iel = 1, ncel

    if(igrdpp.gt.0) then
      rtprom = rtp(iel,isca(irho))
    else
      rtprom = rtpa(iel,isca(irho))
    endif

    trav(iel,1) = ( rtprom*gx - w1(iel) )*volume(iel)
    trav(iel,2) = ( rtprom*gy - w2(iel) )*volume(iel)
    trav(iel,3) = ( rtprom*gz - w3(iel) )*volume(iel)
  enddo
endif

!    Calcul des efforts aux parois (partie 2/5), si demande
!    La pression a la face est calculee comme dans gradrc/gradmc
if (ineedf.eq.1) then
  iclipr = iclrtp(ipriph,icoef)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)
    pip =  rtpa(iel,ipriph)                                       &
         +diipbx*w1(iel) +diipby*w2(iel)                          &
         +diipbz*w3(iel)
    pfac = coefa(ifac,iclipr) +coefb(ifac,iclipr)*pip
    pfac1= rtpa(iel,ipriph)                                       &
         +(cdgfbo(1,ifac)-xyzcen(1,iel))*w1(iel)                  &
         +(cdgfbo(2,ifac)-xyzcen(2,iel))*w2(iel)                  &
         +(cdgfbo(3,ifac)-xyzcen(3,iel))*w3(iel)
    pfac = coefb(ifac,iclipr)*(extrag(ipriph)*pfac1               &
         +(1.d0-extrag(ipriph))*pfac)                             &
         +(1.d0-coefb(ifac,iclipr))*pfac
    do isou = 1, 3
      ra(iforbr+(ifac-1)*ndim + isou-1) =                         &
           ra(iforbr+(ifac-1)*ndim + isou-1)                      &
           + pfac*surfbo(isou,ifac)
    enddo
  enddo
endif


!     Elimination du flux au bord associé au gradient de pression :
!       il est pris en compte par les conditions aux limites dans
!       le flux de Rusanov

if(iifbru.gt.0) then

  do ifac = 1, nfabor

    if(ia(iifru+ifac-1).eq.1) then

      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      pip = rtp(iel,ipriph)                                       &
           +(w1(iel)*diipbx+w2(iel)*diipby+w3(iel)*diipbz)

      pbord = coefa(ifac,iclrtp(ipriph,icoef))                    &
           + coefb(ifac,iclrtp(ipriph,icoef))*pip

      trav(iel,1) = trav(iel,1) + pbord*surfbo(1,ifac)
      trav(iel,2) = trav(iel,2) + pbord*surfbo(2,ifac)
      trav(iel,3) = trav(iel,3) + pbord*surfbo(3,ifac)

    endif

  enddo

endif

!     Flux de C    .L. associé à Rusanov (PROPFB contient la contribution
!       de - div(rho u u) - grad P si on est passé dans cfrusb
!       ou 0 sinon).
!     Pour ne pas ajouter le flux div(rho u u) deux fois, on a remplace
!       codits et bilsc2 par cfcdts et cfbsc2 qui ne different des
!       precedents que par les indicateurs qui permettent de
!       ne pas prendre en compte le flux convectif aux faces de bord
!       pour lesquelles on est passe dans cfrusb

do ifac = 1, nfabor
  iel = ifabor(ifac)
  trav(iel,1) =  trav(iel,1) - propfb(ifac,ipprob(ifbrhu))
  trav(iel,2) =  trav(iel,2) - propfb(ifac,ipprob(ifbrhv))
  trav(iel,3) =  trav(iel,3) - propfb(ifac,ipprob(ifbrhw))
enddo


! ---> 2/3 RHO * GRADIENT DE K SI k epsilon
!      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
!           CA COMPLIQUERAIT LA GESTION DES CL ...

if( (itytur.eq.2 .or. iturb.eq.50                   &
     .or. iturb.eq.60) .and.igrhok.eq.1) then
  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ikiph)
  imligp = imligr(ikiph)
  epsrgp = epsrgr(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

  iwarnp = iwarni(iuiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   ikiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w6     , w6     , w6     ,                                     &
   rtp(1,ikiph)    , coefa(1,iclik)  , coefb(1,iclik)  ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )

  d2s3 = 2.d0/3.d0
  do iel = 1, ncel
    romvom = -rtp(iel,isca(irho))*volume(iel)*d2s3
    trav(iel,1) = trav(iel,1) + w1(iel) * romvom
    trav(iel,2) = trav(iel,2) + w2(iel) * romvom
    trav(iel,3) = trav(iel,3) + w3(iel) * romvom
  enddo

!    Calcul des efforts aux parois (partie 3/5), si demande
  if (ineedf.eq.1) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      xkb = rtpa(iel,ikiph) + diipbx*w1(iel)                      &
           + diipby*w2(iel) + diipbz*w3(iel)
      xkb = coefa(ifac,iclik)+coefb(ifac,iclik)*xkb
      xkb = d2s3*propce(iel,ipcrom)*xkb
      do isou = 1, 3
        ra(iforbr+(ifac-1)*ndim + isou-1) =                       &
             ra(iforbr+(ifac-1)*ndim + isou-1)                    &
             + xkb*surfbo(isou,ifac)
      enddo
    enddo
  endif

endif



! ---> TERMES DE GRADIENT TRANSPOSE

if (ivisse.eq.1) then

  call vissec                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   trav   ,                                                       &
!        ------
   viscf  , viscb  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

endif


! ---> TERMES DE PERTES DE CHARGE
!     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT

if((ncepdp.gt.0).and.(iphydr.eq.0)) then

  idiaex = 1
  call tsepdc                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ncepdp ,                                                       &
   iphas  , idiaex ,                                              &
   icepdc ,                                                       &
   ia     ,                                                       &
   rtp    , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

endif


! ---> - DIVERGENCE DE RIJ

if(itytur.eq.3 ) then

  do isou = 1, 3

    if(isou.eq.1) ivar = iuiph
    if(isou.eq.2) ivar = iviph
    if(isou.eq.3) ivar = iwiph

    call divrij                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   isou   , ivar   , iphas  ,                                     &
   ia     ,                                                       &
   rtp    , propce , propfa , propfb ,                            &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )

    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,              &
                                   ifacel,ifabor,viscf,viscb,w1)

    do iel = 1, ncel
       trav(iel,isou) = trav(iel,isou) -  w1(iel)
    enddo

  enddo

endif



! ---> "VITESSE" DE DIFFUSION FACETTE
!      SI ON FAIT AUTRE CHOSE QUE DU K EPS, IL FAUDRA LA METTRE
!        DANS LA BOUCLE

if( idiff(iuiph).ge. 1 ) then

! --- Si la vitesse doit etre diffusee, on calcule la viscosite
!       pour le second membre (selon Rij ou non)

  if (itytur.eq.3) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
    enddo
  endif

  call viscfa                                                     &
  !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

!     Quand on n'est pas en Rij ou que irijnu = 0, les tableaux
!       VISCFI, VISCBI se trouvent remplis par la meme occasion
!       (ils sont confondus avec VISCF, VISCB)
!     En Rij avec irijnu = 1, on calcule la viscosite increment
!       de la matrice dans VISCFI, VISCBI

  if(itytur.eq.3 .and. irijnu.eq.1) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
    enddo
    call viscfa                                                   &
    !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w1     ,                                                       &
   viscfi , viscbi ,                                              &
   ra     )
  endif

else

! --- Si la vitesse n'a pas de diffusion, on annule la viscosite
!      (matrice et second membre sont dans le meme tableau,
!       sauf en Rij avec IRIJNU = 1)

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  if(itytur.eq.3.and.irijnu.eq.1) then
    do ifac = 1, nfac
      viscfi(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscbi(ifac) = 0.d0
    enddo
  endif

endif


! 2.2  RESOLUTION IMPLICITE NON COUPLEE DES 3 COMPO. DE VITESSES
! ==============================================================

! ---> BOUCLE SUR LES DIRECTIONS DE L'ESPACE (U, V, W)


! Remarque : On suppose que le couplage vitesse pression
!  n'est valable que pour une seule phase.

do isou = 1, 3

  if(isou.eq.1) then
    ivar = iuiph
  endif
  if(isou.eq.2) then
    ivar = iviph
  endif
  if(isou.eq.3) then
    ivar = iwiph
  endif
  ipp  = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef)
  iclvaf = iclrtp(ivar,icoeff)


! ---> TERMES SOURCES UTILISATEURS

  do iel = 1, ncel
    smbr  (iel) = 0.d0
    drtp  (iel) = 0.d0
  enddo

  maxelt = max(ncelet, nfac, nfabor)
  ils    = idebia
  ifinia = ils + maxelt
  call iasize('cfqdmv',ifinia)

  call ustsns                                                     &
  !==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   ivar   , iphas  ,                                              &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbr   , drtp   ,                                              &
!        ------   ------
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )


  do iel = 1, ncel
    rovsdt(iel) = max(-drtp(iel),zero)
    smbr  (iel) = smbr(iel) + drtp(iel) * rtp(iel,ivar)
  enddo


! ---> TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                           ifacel,ifabor,flumas,flumab,w1)



! ---> AJOUT DANS LE TERME SOURCE ET DANS LE TERME INSTATIONNAIRE

  do iel = 1, ncel
    smbr(iel) = smbr  (iel) +                                     &
         trav(iel,isou)+iconv(ivar)*w1(iel)*rtpa(iel,ivar)
  enddo

  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel) +                                   &
         istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)     &
         -iconv(ivar)*w1(iel)
  enddo


! ---> PERTES DE CHARGE

  if (ncepdp.gt.0) then
    do ielpdc = 1, ncepdp
      iel = icepdc(ielpdc)
      rovsdt(iel) = rovsdt(iel) +                                 &
           propce(iel,ipcrom)*volume(iel)*ckupdc(ielpdc,isou)
    enddo
  endif


! --->  TERMES DE SOURCE DE MASSE

  if (ncesmp.gt.0) then
    iterns = 1
    call catsma ( ncelet , ncel , ncesmp , iterns ,               &
                  isno2t, thetav(ivar),                    &
                  icetsm , itypsm(1,ivar)  ,                      &
                  volume , rtp(1,ivar) , smacel(1,ivar) ,         &
                  smacel(1,ipr) , smbr , rovsdt , w1)
  endif



! ---> PARAMETRES POUR LA RESOLUTION DU SYSTEME

  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  ireslp = iresol(ivar)
  ndircp = ndircl(ivar)
  nitmap = nitmax(ivar)
!MO        IMRGRA
  nswrsp = nswrsm(ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  imgrp  = imgr  (ivar)
  ncymxp = ncymax(ivar)
  nitmfp = nitmgf(ivar)
!MO        IPP
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsilp = epsilo(ivar)
  epsrsp = epsrsm(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  thetap = thetav(ivar)
  iescap = 0


! ---> FIN DE LA CONSTRUCTION ET DE LA RESOLUTION DU SYSTEME

  call cfcdts                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ivar   , iconvp , idiffp , ireslp , ndircp ,  nitmap ,         &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , iifbru ,                            &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   ia(iifru) ,                                                    &
   ia     ,                                                       &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     flumas , flumab ,                            &
   viscfi , viscbi , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )


!     PAS DE COUPLAGE INSTATIONNAIRE EN COMPRESSIBLE

enddo

! --->  FIN DE LA BOUCLE SUR U, V, W,


! ---> IMPRESSION DE NORME

if (iwarni(iuiph).ge.2) then
  rnorm = -1.d0
  do iel = 1, ncel
    vitnor =                                                      &
     sqrt(rtp(iel,iuiph)**2+rtp(iel,iviph)**2+rtp(iel,iwiph)**2)
    rnorm = max(rnorm,vitnor)
  enddo
  if (irangp.ge.0) call parmax (rnorm)
                             !==========
  write(nfecra,1100) rnorm
endif


!--------
! FORMATS
!--------

 1100 format(/,                                                   &
 1X,'Vitesse maximale apres qdm ',E12.4)
!----
! FIN
!----

return
end subroutine
