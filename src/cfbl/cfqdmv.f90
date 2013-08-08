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

subroutine cfqdmv &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   flumas , flumab ,                                              &
   coefa  , coefb  , ckupdc , smacel , frcxt  ,                   &
   tpucou )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
! frcxt(3,ncelet)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! tpucou           ! tr ! --> ! couplage vitesse pression                      !
! (ncelel,ndim)    !    !     !                                                !
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
use pointe, only: forbr
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

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision flumas(nfac), flumab(nfabor)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(3, ncelet)
double precision tpucou(ncelet,ndim)

! Local variables

integer          iel   , ielpdc, ifac  , ivar  , isou
integer          iccocg, inc   , init  , ii
integer          ireslp, nswrgp, imligp, iwarnp, ipp
integer          iclik , iclvar, iclvaf, iclipr
integer          ipcrom, ipcvis, ipcvst
integer          iconvp, idiffp, ndircp, nitmap, nswrsp
integer          ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          idiaex, iterns
double precision rnorm , vitnor
double precision romvom, rtprom
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision vit1  , vit2  , vit3  , thetap, pfac, pfac1
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , pbord , diipbx, diipby, diipbz, pip, xkb

double precision rvoid(1)

double precision, allocatable, dimension(:), target :: viscf, viscb
double precision, allocatable, dimension(:), target :: wvisfi, wvisbi
double precision, allocatable, dimension(:) :: drtp, smbr, rovsdt
double precision, allocatable, dimension(:,:) :: trav, grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: dfrcxt

double precision, pointer, dimension(:) :: viscfi => null(), viscbi => null()

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate temporary arrays for the velocity-pressure resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(drtp(ncelet), smbr(ncelet), rovsdt(ncelet))
allocate(trav(ncelet,3))
allocate(dfrcxt(3,ncelet))

! Allocate other arrays, depending on user options
if (itytur.eq.3.and.irijnu.eq.1) then
  allocate(wvisfi(nfac), wvisbi(nfabor))
  viscfi => wvisfi(1:nfac)
  viscbi => wvisbi(1:nfabor)
else
  viscfi => viscf(1:nfac)
  viscbi => viscb(1:nfabor)
endif

! Allocate work arrays
allocate(w1(ncelet))

! Initialize variables to avoid compiler warnings

iclik = 0

! Memoire


if(itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) then
  iclik = iclrtp(ik ,icoef)
endif

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)

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

    dfrcxt(1, iel) = rtprom*gx - frcxt(1, iel)
    dfrcxt(2, iel) = rtprom*gy - frcxt(2, iel)
    dfrcxt(3, iel) = rtprom*gz - frcxt(3, iel)
  enddo
!     Ajout eventuel des pertes de charges
  if (ncepdp.gt.0) then
    do ielpdc = 1, ncepdp
      iel=icepdc(ielpdc)
      vit1   = rtp(iel,iu)
      vit2   = rtp(iel,iv)
      vit3   = rtp(iel,iw)
      cpdc11 = ckupdc(ielpdc,1)
      cpdc22 = ckupdc(ielpdc,2)
      cpdc33 = ckupdc(ielpdc,3)
      cpdc12 = ckupdc(ielpdc,4)
      cpdc23 = ckupdc(ielpdc,5)
      cpdc13 = ckupdc(ielpdc,6)
      dfrcxt(1 ,iel) = dfrcxt(1 ,iel)                   &
 -rtp(iel,isca(irho))*(cpdc11*vit1+cpdc12*vit2+cpdc13*vit3)
      dfrcxt(2 ,iel) = dfrcxt(2 ,iel)                   &
 -rtp(iel,isca(irho))*(cpdc12*vit1+cpdc22*vit2+cpdc23*vit3)
      dfrcxt(3 ,iel) = dfrcxt(3 ,iel)                   &
 -rtp(iel,isca(irho))*(cpdc13*vit1+cpdc23*vit2+cpdc33*vit3)
    enddo
  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(dfrcxt)
  endif

endif

!       Fin du test sur IPHYDR


! ---> PRISE EN COMPTE DU GRADIENT DE PRESSION

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
   rtp(1,ipr)   , coefa(1,iclrtp(ipr,icoef))  ,                   &
                  coefb(1,iclrtp(ipr,icoef))  ,                   &
   grad   )


if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(iel,1) = ( frcxt(1,iel) - grad(iel,1) )*volume(iel)
    trav(iel,2) = ( frcxt(2,iel) - grad(iel,2) )*volume(iel)
    trav(iel,3) = ( frcxt(3,iel) - grad(iel,3) )*volume(iel)
  enddo
else
  do iel = 1, ncel

    if(igrdpp.gt.0) then
      rtprom = rtp(iel,isca(irho))
    else
      rtprom = rtpa(iel,isca(irho))
    endif

    trav(iel,1) = ( rtprom*gx - grad(iel,1) )*volume(iel)
    trav(iel,2) = ( rtprom*gy - grad(iel,2) )*volume(iel)
    trav(iel,3) = ( rtprom*gz - grad(iel,3) )*volume(iel)
  enddo
endif

!    Calcul des efforts aux parois (partie 2/5), si demande
!    La pression a la face est calculee comme dans gradrc/gradmc
if (ineedf.eq.1) then
  iclipr = iclrtp(ipr,icoef)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)
    pip =  rtpa(iel,ipr) &
        +diipbx*grad(iel,1) +diipby*grad(iel,2) +diipbz*grad(iel,3)
    pfac = coefa(ifac,iclipr) +coefb(ifac,iclipr)*pip
    pfac1= rtpa(iel,ipr)                                          &
         +(cdgfbo(1,ifac)-xyzcen(1,iel))*grad(iel,1)              &
         +(cdgfbo(2,ifac)-xyzcen(2,iel))*grad(iel,2)              &
         +(cdgfbo(3,ifac)-xyzcen(3,iel))*grad(iel,3)
    pfac = coefb(ifac,iclipr)*(extrag(ipr)*pfac1                  &
         +(1.d0-extrag(ipr))*pfac)                                &
         +(1.d0-coefb(ifac,iclipr))*pfac
    do isou = 1, 3
      forbr(isou,ifac) = forbr(isou,ifac) + pfac*surfbo(isou,ifac)
    enddo
  enddo
endif


!     Elimination du flux au bord associé au gradient de pression :
!       il est pris en compte par les conditions aux limites dans
!       le flux de Rusanov

do ifac = 1, nfabor

  if(ifbrus(ifac).eq.1) then

    iel = ifabor(ifac)

    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)

    pip = rtp(iel,ipr) &
        +diipbx*grad(iel,1) +diipby*grad(iel,2) +diipbz*grad(iel,3)

    pbord = coefa(ifac,iclrtp(ipr,icoef))                    &
         + coefb(ifac,iclrtp(ipr,icoef))*pip

    trav(iel,1) = trav(iel,1) + pbord*surfbo(1,ifac)
    trav(iel,2) = trav(iel,2) + pbord*surfbo(2,ifac)
    trav(iel,3) = trav(iel,3) + pbord*surfbo(3,ifac)

  endif

enddo

! Free memory
deallocate(grad)

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

if( (itytur.eq.2 .or. iturb.eq.50 .or. iturb.eq.60) .and. igrhok.eq.1) then

  ! Allocate a work array for the gradient calculation
  allocate(grad(ncelet,3))

  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ik)
  imligp = imligr(ik)
  epsrgp = epsrgr(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  iwarnp = iwarni(iu)

  call grdcel                                                     &
  !==========
 ( ik  , imrgra , inc    , iccocg , nswrgp , imligp ,             &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtp(1,ik)    , coefa(1,iclik)  , coefb(1,iclik)  ,             &
   grad   )

  d2s3 = 2.d0/3.d0
  do iel = 1, ncel
    romvom = -rtp(iel,isca(irho))*volume(iel)*d2s3
    trav(iel,1) = trav(iel,1) + grad(iel,1) * romvom
    trav(iel,2) = trav(iel,2) + grad(iel,2) * romvom
    trav(iel,3) = trav(iel,3) + grad(iel,3) * romvom
  enddo

!    Calcul des efforts aux parois (partie 3/5), si demande
  if (ineedf.eq.1) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      xkb = rtpa(iel,ik) &
          + diipbx*grad(iel,1) + diipby*grad(iel,2) + diipbz*grad(iel,3)
      xkb = coefa(ifac,iclik)+coefb(ifac,iclik)*xkb
      xkb = d2s3*propce(iel,ipcrom)*xkb
      do isou = 1, 3
        forbr(isou,ifac) = forbr(isou,ifac) + xkb*surfbo(isou,ifac)
      enddo
    enddo
  endif

  ! Free memory
  deallocate(grad)

endif



! ---> TERMES DE GRADIENT TRANSPOSE

if (ivisse.eq.1) then

  call vissec                                                     &
  !==========
 ( rtpa   , propce , propfb ,                                     &
   coefa  , coefb  ,                                              &
   trav   ,                                                       &
   viscf  , viscb  )

endif


! ---> TERMES DE PERTES DE CHARGE
!     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT

if((ncepdp.gt.0).and.(iphydr.eq.0)) then

  idiaex = 1
  call tsepdc                                                     &
  !==========
 ( ncepdp ,                                                       &
   idiaex ,                                                       &
   icepdc ,                                                       &
   rtp    , propce , ckupdc , trav   )

endif


! ---> - DIVERGENCE DE RIJ

if(itytur.eq.3 ) then

  do isou = 1, 3

    if(isou.eq.1) ivar = iu
    if(isou.eq.2) ivar = iv
    if(isou.eq.3) ivar = iw

    call divrij                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   isou   , ivar   ,                                              &
   rtp    , propce , propfb ,                                     &
   coefa  , coefb  ,                                              &
   viscf  , viscb  )

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

if( idiff(iu).ge. 1 ) then

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
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

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
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscfi , viscbi )
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
    ivar = iu
  endif
  if(isou.eq.2) then
    ivar = iv
  endif
  if(isou.eq.3) then
    ivar = iw
  endif
  ipp  = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef)
  iclvaf = iclrtp(ivar,icoeff)


! ---> TERMES SOURCES UTILISATEURS

  do iel = 1, ncel
    smbr  (iel) = 0.d0
    drtp  (iel) = 0.d0
  enddo

  call ustsns                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , propce , propfb ,                            &
   ckupdc , smacel , smbr   , drtp   )

  do iel = 1, ncel
    rovsdt(iel) = max(-drtp(iel),zero)
    smbr  (iel) = smbr(iel) + drtp(iel) * rtp(iel,ivar)
  enddo

! ---> AJOUT DANS LE TERME SOURCE ET DANS LE TERME INSTATIONNAIRE

  do iel = 1, ncel
    smbr(iel) = smbr(iel) + trav(iel,isou)
  enddo

  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel) +                                   &
         istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)
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
 ( nvar   , nscal  ,                                              &
   ivar   , iconvp , idiffp , ireslp , ndircp ,  nitmap ,         &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap , thetap , &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     flumas , flumab ,                            &
   viscfi , viscbi , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   rvoid  )


!     PAS DE COUPLAGE INSTATIONNAIRE EN COMPRESSIBLE

enddo

! --->  FIN DE LA BOUCLE SUR U, V, W,


! ---> IMPRESSION DE NORME

if (iwarni(iu).ge.2) then
  rnorm = -1.d0
  do iel = 1, ncel
    vitnor =                                                      &
     sqrt(rtp(iel,iu)**2+rtp(iel,iv)**2+rtp(iel,iw)**2)
    rnorm = max(rnorm,vitnor)
  enddo
  if (irangp.ge.0) call parmax (rnorm)
                             !==========
  write(nfecra,1100) rnorm
endif

! Free memory
deallocate(viscf, viscb)
deallocate(drtp, smbr, rovsdt)
deallocate(trav)
deallocate(dfrcxt)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
deallocate(w1)

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
