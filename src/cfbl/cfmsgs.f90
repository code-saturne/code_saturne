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

subroutine cfmsgs &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   flumas , flumab , flabgs , flbbgs ,                            &
   trflms , trflmb )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DU "FLUX DE MASSE" AUX FACES
!   POUR LA RESOLUTION DE LA MASSE VOLUMIQUE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! flumas(nfac)     ! tr ! --> ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! --> ! flux de masse aux faces de bord                !
! trflms(nfac)     ! tr ! --- ! tableau de travail                             !
! trflmb(nfabor    ! tr ! --- ! tableau de travail                             !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision flumas(nfac), flumab(nfabor)
double precision flabgs(nfac), flbbgs(nfabor)
double precision trflms(nfac), trflmb(nfabor)

! Local variables

integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, ii, jj
integer          ipp
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp
integer          ircflp, ischcp, isstpp

integer          iirom , iiromb
integer          ivar0 , imvis1, iccfth, imodif, isou
integer          imaspe, iflmb0, itypfl
integer          icliup, iclivp, icliwp, iclvar
integer          itsqdm, iiun  , iextts

double precision epsrgp, climgp, extrap, blencp
double precision flui  , fluj  , pfac  , thetv

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10, w11, w12
double precision, allocatable, dimension(:,:) :: coefu
double precision, allocatable, dimension(:,:) :: coefuf
!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(coefu(nfabor,3))
allocate(coefuf(nfabor,3))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet), w11(ncelet), w12(ncelet))


! --- Numero des variables de calcul
!     Masse volumique
ivar   = isca(iscal)

!     Masse volumique dans PROPCE
iirom  = ipproc(irom  )
iiromb = ipprob(irom  )

! ---> Initialisation du flux de masse

do ifac = 1, nfac
  flumas(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  flumab(ifac) = 0.d0
enddo
do ifac = 1, nfac
  flabgs(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  flbbgs(ifac) = 0.d0
enddo

!===============================================================================
! 2. FLUX DE MASSE AUX FACES
!===============================================================================

!     2.1 TERMES SOURCES DE L'EQUATION DE QDM
!     =======================================

!     FX -> W5 , FY -> W6 , FZ -> W7

!     Apres premiers tests, il semble (tube à choc en double détente)
!       que ce ne soit pas une bonne idée de prendre en compte dans
!       l'équation de la masse tous les termes de l'équation de la
!       quantité de mouvement (en particulier le terme convectif, mais
!       les termes de diffusion, en gradient transposé, source de masse
!       et utilisateur n'ont pas été vraiment testés, dans la mesure ou
!       ils étaient nuls pour tous les cas unitaires considérés).
!       Bien entendu, s'agissant de tests préliminaires, il est possible
!       que le comportement insatisfaisant observé soit le fruit d'une
!       simple erreur de codage (mais on ne l'a pas trouvée...).
!     On propose donc par sécurité de ne pas prendre en compte les termes
!       source de l'équation de la quantité de mouvement, hormis le terme
!       de gravité (car il contrebalance le gradient de pression et son
!       effet est bien visible lorsqu'on est dans une situation
!       d'équilibre).
!     Cependant, pour l'avenir, on conserve ici le codage
!       préliminaire qui a permis d'effectuer les tests (version 1.1.0.h).
!       On l'encapsule dans un test qui le désactive obligatoirement
!       (ainsi, pas de question utilisateur intempestive, et toujours
!       la possibilité de poursuivre les tests si l'on a le temps  ou si
!       on rencontre des pbs...).
!     Noter que, avec ces termes, le système devient difficile à analyser
!       sur le papier (alors que sans ces termes, on est dans la
!       configuration Euler + gravité).

! --- Initialisation
do iel = 1, ncel
  w10(iel) = 0.d0
enddo
do iel = 1, ncel
  w11(iel) = 0.d0
enddo
do iel = 1, ncel
  w12(iel) = 0.d0
enddo


!     Test sur les termes source de qdm
itsqdm = 0
if(itsqdm.ne.0) then


! --- Terme source utilisateur

!     Suivant X
  call ustsns                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iu  ,                                                          &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfb ,                            &
   ckupdc , smacel , w10    , w9     )

!     Suivant Y
  call ustsns                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iv  ,                                                          &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfb ,                            &
   ckupdc , smacel , w11    , w9     )

!     Suivant Z
  call ustsns                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iw  ,                                                          &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfb ,                            &
   ckupdc , smacel , w12    , w9     )


! --- Terme de convection de quantite de mouvement
  if(iconv(iu).ge.1) then

    icliup = iclrtp(iu ,icoef)
    iclivp = iclrtp(iv ,icoef)
    icliwp = iclrtp(iw ,icoef)

    init   = 1
    inc    = 1
    iccocg = 1
    iflmb0 = 1
    nswrgp = nswrgr(iu)
    imligp = imligr(iu)
    iwarnp = iwarni(iu)
    epsrgp = epsrgr(iu)
    climgp = climgr(iu)
    extrap = extrag(iu)

    imaspe = 1
    itypfl = 1

!     Calcul du flux de masse
    call inimas                                                   &
    !==========
 ( iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,iirom) , propfb(1,iiromb),                            &
   rtpa (1,iu)  , rtpa (1,iv)  , rtpa (1,iw)  ,                   &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
   flumas , flumab )

!     Calcul du terme convecte suivant les 3 directions
!       sans reconstruction
    do isou = 1, 3
      if(isou.eq.1) ivar0  = iu
      if(isou.eq.1) iclvar = icliup
      if(isou.eq.2) ivar0  = iv
      if(isou.eq.2) iclvar = iclivp
      if(isou.eq.3) ivar0  = iw
      if(isou.eq.3) iclvar = icliwp

      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
        fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
        trflms(ifac) = -(flui*rtpa(ii,ivar0)+fluj*rtpa(jj,ivar0))
      enddo

      do ifac = 1, nfabor
        ii = ifabor(ifac)
        flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
        fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
        pfac = coefa(ifac,iclvar)                                 &
             + coefb(ifac,iclvar)*rtpa(ii,ivar0)
        trflmb(ifac) = - ( flui*rtpa(ii,ivar0) + fluj*pfac )
      enddo

      init = 0
      if(isou.eq.1) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w10)
      elseif(isou.eq.2) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w11)
      elseif(isou.eq.3) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w12)
      endif

    enddo

  endif


! --- Terme de viscosite

  if( idiff(iu).ge.1 ) then

    do iel = 1, ncelet
      w8(iel) = 1.d0
      w9(iel) = 0.d0
    enddo

    call cfdivs                                                   &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfb ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w10    , w8     , w9     , w9     )
!        ------

    call cfdivs                                                   &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfb ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w11    , w9     , w8     , w9     )
!        ------

    call cfdivs                                                   &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfb ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w12    , w9     , w9     , w8     )
!        ------

  endif


! --- Terme source de masse
!     On met tout en explicite pour l'instant ... a modifier lors des
!      tests sur ce terme (cf. remarque en debut de fichier)
  iiun   = 1
  iextts = 0
  thetv  = 0.d0
  do isou = 1, 3
    if(isou.eq.1) then
      ivar0  = iu
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w10   , w1     , w2 )
      do iel = 1, ncel
        w10(iel) = w10(iel) + w2(iel)
      enddo

    elseif(isou.eq.2) then
      ivar0  = iv
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w11   , w1     , w2 )
      do iel = 1, ncel
        w11(iel) = w11(iel) + w2(iel)
      enddo

    elseif(isou.eq.3) then
      ivar0  = iw
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w12   , w1     , w2 )
      do iel = 1, ncel
        w12(iel) = w12(iel) + w2(iel)
      enddo

    endif

  enddo

endif
!     Fin du Test sur les termes source de qdm


! --- Terme de forces volumiques (gravite)
do iel = 1, ncel
  w5(iel) = gx + w10(iel)/rtpa(iel,ivar)
  w6(iel) = gy + w11(iel)/rtpa(iel,ivar)
  w7(iel) = gz + w12(iel)/rtpa(iel,ivar)
enddo


! 2.2 CALCUL DU TERME DE DIFFUSION D'ENTROPIE
! ===========================================

! --- Calcul de l'entropie au centre des cellules et affectation a W1
iccfth = 6
imodif = 0
call cfther                                                       &
!==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   w1     , w8     , w9     , w10    )

! --- Communication de l'entropie
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w1)
  !==========
endif

! --- Calcul de dt*Beta/Rho au centre des cellules et affectation a W2
iccfth = 162
imodif = 0
call cfther                                                       &
!==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   w2     , w8     , w9     , w10    )

! --- Pour la condition au bord sur l'entropie
!     COEFA=COEFA(.,ITEMPK) et COEFB=COEFB(.,ITEMPK)


do iel = 1, ncel
  w2(iel) =  dt(iel) * w2(iel)
enddo

! --- Calcul de dt*Beta aux faces et affectation a TRFLMS et TRFLMB
imvis1 = 1

call viscfa                                                       &
!==========
 ( imvis1 ,                                                       &
   w2     ,                                                       &
   trflms , trflmb )

! --- Calcul du flux de diffusion

!     Conditions aux limites de flux nul pour l'entropie
!       on utilise COEFU comme tableau de travail)
do ifac = 1, nfabor
  coefu(ifac,1) = 0.d0
  coefu(ifac,2) = 1.d0
  coefuf(ifac,1) = 0.d0
  coefuf(ifac,2) = 0.d0
enddo


inc =1
iccocg = 1
ivar0  = 0
iconvp = 0
idiffp = 1
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
ipp    = ipprtp(ivar)
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

call cfbsc3                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   ivar0  , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   w1     , coefu(1,1)      , coefu(1,2)      ,                   &
            coefuf(1,1)     , coefuf(1,2)     ,                   &
   trflms , trflmb , trflms , trflmb ,                            &
   flabgs , flbbgs )


! 2.3 CALCUL DU FLUX DE MASSE AUX FACES
! =====================================

! --- Calcul des "vitesses" de convection au centre des cellules

do iel = 1, ncel
  w10(iel) = rtpa(iel,iu) + dt(iel)*w5(iel)
  w11(iel) = rtpa(iel,iv) + dt(iel)*w6(iel)
  w12(iel) = rtpa(iel,iw) + dt(iel)*w7(iel)
enddo

! --- Calcul du flux par appel a INIMAS

!     Pour éviter d'appliquer une condition limite inadaptée, on impose
!       un simple flux nul. Noter que cela ne sert qu'a reconstruire des
!       gradients. La valeur au bord importe peu, dans la mesure ou le
!       flux de bord est ecrasé ensuite.

!     Pour éviter d'avoir recours a d'autres tableaux, on simule un flux
!       nul avec :
!       TRFLMB = 1 pour COEFB = 1 et
!       INC    = 0 pour COEFA = 0

!     On prend ROM = W1 = 1 et ROMB = TRFLMB = 1
!       (et ROMB sert aussi pour COEFB=1)
do iel = 1, ncel
  w1(iel) = 1.d0
enddo
do ifac = 1, nfabor
  trflmb(ifac) = 1.d0
enddo

icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

init   = 1
inc    = 0
!              ^ Comme indiqué ci-dessus, pour un flux nul
iccocg = 1
ivar0  = 0
iflmb0 = 1
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

imaspe = 1
itypfl = 1

call inimas                                                       &
!==========
 ( ivar0  , ivar0  , ivar0  , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   w1     , trflmb ,                                              &
   w10    , w11    , w12    ,                                     &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   trflmb          , trflmb          , trflmb          ,          &
   flumas , flumab )

! Free memory
deallocate(coefu)
deallocate(coefuf)
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)
deallocate(w10, w11, w12)

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine
