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

subroutine cfmsgs &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   flumas , flumab , flabgs , flbbgs ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   trflms , trflmb , coefu  , xam    ,                            &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! flumas(nfac)     ! tr ! --> ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! --> ! flux de masse aux faces de bord                !
! w1..12(ncelet    ! tr ! --- ! tableau de travail                             !
! trflms(nfac)     ! tr ! --- ! tableau de travail                             !
! trflmb(nfabor    ! tr ! --- ! tableau de travail                             !
! coefu(nfabo,3    ! tr ! --- ! tableau de travail cl de la qdm                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision flumas(nfac), flumab(nfabor)
double precision flabgs(nfac), flbbgs(nfabor)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)
double precision w4(ncelet) , w5(ncelet) , w6(ncelet)
double precision w7(ncelet) , w8(ncelet) , w9(ncelet)
double precision w10(ncelet), w11(ncelet), w12(ncelet)
double precision trflms(nfac), trflmb(nfabor)
double precision coefu(nfabor,3), xam(nfac,2)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          ivar  , iphas
integer          ifac  , iel
integer          init  , inc   , iccocg, ii, jj
integer          ipp
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp
integer          ircflp, ischcp, isstpp

integer          iirom , iiromb
integer          ivar0 , imvis1, iccfth, imodif, isou
integer          imaspe, iflmb0, iismph
integer          icliup, iclivp, icliwp, iclvar
integer          iuiph , iviph , iwiph
integer          itsqdm, iiun  , iextts
integer          maxelt, ils

double precision epsrgp, climgp, extrap, blencp
double precision flui  , fluj  , pfac  , thetv

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero de phase associee au scalaire traite
iphas  = 1

! --- Numero des variables de calcul
!     Masse volumique
ivar   = isca(iscal)
!     Vitesses
iuiph  = iu
iviph  = iv
iwiph  = iw

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

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
call iasize('cfmsgs',ifinia)

!     Suivant X
  call ustsns                                                     &
  !==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iuiph  , iphas  ,                                              &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w10    , w9     ,                                              &
!        ------   ------
   w8     , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Suivant Y
  call ustsns                                                     &
  !==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iviph  , iphas  ,                                              &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w11    , w9     ,                                              &
!        ------   ------
   w8     , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!     Suivant Z
  call ustsns                                                     &
  !==========
 ( ifinia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iwiph  , iphas  ,                                              &
   maxelt , ia(ils),                                              &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w12    , w9     ,                                              &
!        ------   ------
   w8     , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )


! --- Terme de convection de quantite de mouvement
  if(iconv(iuiph).ge.1) then

    icliup = iclrtp(iuiph ,icoef)
    iclivp = iclrtp(iviph ,icoef)
    icliwp = iclrtp(iwiph ,icoef)

    init   = 1
    inc    = 1
    iccocg = 1
    iflmb0 = 1
    iismph = iisymp+nfabor*(iphas-1)
    nswrgp = nswrgr(iuiph)
    imligp = imligr(iuiph)
    iwarnp = iwarni(iuiph)
    epsrgp = epsrgr(iuiph)
    climgp = climgr(iuiph)
    extrap = extrag(iuiph)

    imaspe = 1

!     Calcul du flux de masse
    call inimas                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   iuiph  , iviph  , iwiph  , imaspe , iphas  ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia(iismph) ,                                                   &
   ia     ,                                                       &
   propce(1,iirom) , propfb(1,iiromb),                            &
   rtpa (1,iuiph)  , rtpa (1,iviph)  , rtpa (1,iwiph)  ,          &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
   flumas , flumab ,                                              &
!        ------   ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )

!     Calcul du terme convecte suivant les 3 directions
!       sans reconstruction
    do isou = 1, 3
      if(isou.eq.1) ivar0  = iuiph
      if(isou.eq.1) iclvar = icliup
      if(isou.eq.2) ivar0  = iviph
      if(isou.eq.2) iclvar = iclivp
      if(isou.eq.3) ivar0  = iwiph
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

  if( idiff(iuiph).ge.1 ) then

    do iel = 1, ncelet
      w8(iel) = 1.d0
      w9(iel) = 0.d0
    enddo

    call cfdivs                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w10    , w8     , w9     , w9     ,                            &
!        ------
   w7     ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    call cfdivs                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w11    , w9     , w8     , w9     ,                            &
!        ------
   w7     ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

    call cfdivs                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w12    , w9     , w9     , w8     ,                            &
!        ------
   w7     ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

  endif


! --- Terme source de masse
!     On met tout en explicite pour l'instant ... a modifier lors des
!      tests sur ce terme (cf. remarque en debut de fichier)
  iiun   = 1
  iextts = 0
  thetv  = 0.d0
  do isou = 1, 3
    if(isou.eq.1) then
      ivar0  = iuiph
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr),                         &
   w10   , w1     , w2 )
      do iel = 1, ncel
        w10(iel) = w10(iel) + w2(iel)
      enddo

    elseif(isou.eq.2) then
      ivar0  = iviph
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr),                         &
   w11   , w1     , w2 )
      do iel = 1, ncel
        w11(iel) = w11(iel) + w2(iel)
      enddo

    elseif(isou.eq.3) then
      ivar0  = iwiph
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr),                         &
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
call uscfth                                                       &
!==========
 ( nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w8     , w9     , w10    )

! --- Communication de l'entropie
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w1)
  !==========
endif

! --- Calcul de dt*Beta/Rho au centre des cellules et affectation a W2
iccfth = 162
imodif = 0
call uscfth                                                       &
!==========
 ( nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
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
 ( idebia , idebra ,                                              &
   imvis1 ,                                                       &
   ia     ,                                                       &
   w2     ,                                                       &
   trflms , trflmb ,                                              &
!        ------   ------
   ra     )

! --- Calcul du flux de diffusion

!     Conditions aux limites de flux nul pour l'entropie
!       on utilise COEFU comme tableau de travail)
do ifac = 1, nfabor
  coefu(ifac,1) = 0.d0
  coefu(ifac,2) = 1.d0
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
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ivar0  , iconvp , idiffp , nswrgp , imligp , ircflp ,          &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap ,                            &
   ia     ,                                                       &
   w1     , coefu(1,1)      , coefu(1,2)      ,                   &
            coefu(1,1)      , coefu(1,2)      ,                   &
   trflms , trflmb , trflms , trflmb ,                            &
   flabgs , flbbgs ,                                              &
!        ------   ------
   w2     , w3     , w4     , w8     , w9     , w10    ,          &
   ra     )


! 2.3 CALCUL DU FLUX DE MASSE AUX FACES
! =====================================

! --- Calcul des "vitesses" de convection au centre des cellules

do iel = 1, ncel
  w10(iel) = rtpa(iel,iuiph) + dt(iel)*w5(iel)
  w11(iel) = rtpa(iel,iviph) + dt(iel)*w6(iel)
  w12(iel) = rtpa(iel,iwiph) + dt(iel)*w7(iel)
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

icliup = iclrtp(iuiph ,icoef)
iclivp = iclrtp(iviph ,icoef)
icliwp = iclrtp(iwiph ,icoef)

init   = 1
inc    = 0
!              ^ Comme indiqué ci-dessus, pour un flux nul
iccocg = 1
ivar0  = 0
iflmb0 = 1
iismph = iisymp+nfabor*(iphas-1)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

imaspe = 1

call inimas                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ivar0  , ivar0  , ivar0  , imaspe , iphas  ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia(iismph) ,                                                   &
   ia     ,                                                       &
   w1     , trflmb ,                                              &
   w10    , w11    , w12    ,                                     &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   trflmb          , trflmb          , trflmb          ,          &
   flumas , flumab ,                                              &
!        ------   ------
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   ra     )

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine
