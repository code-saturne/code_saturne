!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine lagune &
!================

 ( idbia0 , idbra0 ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  , indep  , ibord  ,          &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor , &
   taup   , tlag   , piil   , bx     , vagaus , tsuf   , tsup   , &
   tsvar  , tempct , tsfext , cpgd1  , cpgd2  , cpght  ,          &
   gradpr , gradvf , croule , brgaus , terbru ,                   &
   w1     , w2     , w3     , auxl   , auxl2  ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   Sous-programme principal du module de modelisation Lagrangienne
!   des ecoulements diphasiques a inclusions dispersees.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! icocel           ! te ! --> ! connectivite cellules -> faces                 !
!   (lndnod)       !    !     !    face de bord si numero negatif              !
! itycel           ! te ! --> ! connectivite cellules -> faces                 !
!   (ncelet+1)     !    !     !    pointeur du tableau icocel                  !
! ifrlag           ! te ! --> ! numero de zone de la face de bord              !
!   (nfabor)       !    !     !  pour le module lagrangien                     !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! indep            ! te ! --> ! pour chaque particule :                        !
!   (nbpmax)       !    !     !   numero de la cellule de depart               !
! ibord            ! te ! --> ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant et prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! --> ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! --> ! moyennes statistiques                          !
!(ncelet,nvlsta    !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! tr ! --> ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! tr ! --> ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
! taup(nbpmax)     ! tr ! --> ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! --> ! temps caracteristique fluide                   !
! piil(nbpmax,3    ! tr ! --> ! terme dans l'integration des eds up            !
! bx(nbpmax,3,2    ! tr ! --> ! caracteristiques de la turbulence              !
! vagaus           ! tr ! --> ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! tsup(nbpmax,3    ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse des particules                    !
! tsuf(nbpmax,3    ! tr ! --> ! prediction 1er sous-pas pour                   !
!                  !    !     !   la vitesse du fluide vu                      !
! tsvar            ! tr ! --> ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable courante, utilise pour la           !
! tempct           ! tr ! --> ! temps caracteristique thermique                !
! (nbpmax,2)       !    !     !                                                !
! tsfext(nbpmax    ! tr ! --> ! forces externes                                !
! cpgd1,cpgd2,     ! tr ! --> ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpmax    !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
! gradpr(ncel,3    ! tr ! --> ! gradient de pression                           !
! gradvf(ncel,9    ! tr ! --> ! gradient de vitesse fluide                     !
! croule           ! tr ! --> ! fonction d'importance pour roulette            !
!   (ncelet)       !    !     !   russe                                        !
! w1..w3(ncelet    ! tr ! --- ! tableaux de travail                            !
! auxl(nbpmax,3    ! tr ! --- ! tableau de travail                             !
! auxl2            ! tr ! --- ! tableau de travail                             !
!    (nbpmax,7)    !    !     !                                                !
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
use optcal
use entsor
use cstphy
use cstnum
use pointe
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          lndnod
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          icocel(lndnod) , itycel(ncelet+1)
integer          ifrlag(nfabor) , itepa(nbpmax,nivep)
integer          indep(nbpmax) , ibord(nbpmax)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision taup(nbpmax) , tlag(nbpmax,3) , piil(nbpmax,3)
double precision vagaus(nbpmax,*) , bx(nbpmax,3,2)
double precision tsuf(nbpmax,3) , tsup(nbpmax,3)
double precision tsvar(nbpmax,nvp1)
double precision tempct(nbpmax,2) , tsfext(nbpmax)
double precision cpgd1(nbpmax) , cpgd2(nbpmax) , cpght(nbpmax)
double precision brgaus(nbpmax,*) , terbru(nbpmax)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision croule(ncelet)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision auxl(nbpmax,3) , auxl2(nbpmax,7)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifinia, ifinra

integer          ip     , npt    , iok
integer          nfin   , npars  , iel    , ivf
integer          npar1  , npar2
integer          iforce , iitslg
integer          modntl

double precision dnpars

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 0.  GESTION MEMOIRE ET COMPTEUR DE PASSAGE
!===============================================================================

ipass = ipass + 1

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

iplar = iplar + 1
iplas = iplas + 1

nbpnew = 0
npcsup = 0
npclon = 0
npkill = 0
npencr = 0
nbpout = 0
nbperr = 0

dnbpnw = 0.d0
dnpcsu = 0.d0
dnpclo = 0.d0
dnpkil = 0.d0
dnpenc = 0.d0
dnbpou = 0.d0
dnbper = 0.d0

!-->Sur Champ fige Lagrangien : RTPA = RTP
!   Rem : cette boucle pourrait etre faite au 1er passage
!         mais la presence de usproj incite a la prudence...

if (iilagr.eq.3) then
  do ivf = 1,nvar
    do iel = 1,ncel
      rtpa(iel,ivf) = rtp(iel,ivf)
    enddo
  enddo
endif

!-->au premier passage relatif :

if (iplar.eq.1) then

!       Connectivite cellules -> faces

  call lagdeb                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   icocel , itycel ,                                              &
   ia     ,                                                       &
   ra     )

endif


!===============================================================================
! 2.  MISE A JOUR DES NOUVELLES PARTICULES ENTREES DANS LE DOMAINE
!===============================================================================

! Au premier pas de temps on initalise les particules avec RTP et
! non RTPA car RTPA = initialisation

if ( ntcabs.eq.1 ) then

  call lagent                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel     , icocel      ,                                     &
   ia(iitypf) , ia(iitrif)  , ifrlag , itepa  ,                   &
   ia     ,                                                       &
   ra(isrfbn) , dt     , rtp    , propce , propfa , propfb ,      &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , auxl   , w1     , w2     , w3     , &
   ra     )

else

  call lagent                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel     , icocel      ,                                     &
   ia(iitypf) , ia(iitrif)  , ifrlag , itepa  ,                   &
   ia     ,                                                       &
   ra(isrfbn) , dt     , rtpa   , propce , propfa , propfb ,      &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , auxl   , w1     , w2     , w3     , &
   ra     )
endif

!===============================================================================
! 2.1 CALCUL DE LA FONCTION D'IMPORTANCE POUR LA ROULETTE RUSSE
!===============================================================================

if (iroule.ge.1) then

  call uslaru                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iitypf) , ia(iitrif)  , itepa ,                             &
   ia     ,                                                       &
   ra(isrfbn)  , dt     , rtpa   , propce , propfa , propfb ,     &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , croule , auxl ,                     &
   ra(idipar) , ra(iyppar) ,                                      &
   w1     , w2     , w3     ,                                     &
   ra     )

  iok = 0
  do iel = 1,ncel
    if (croule(iel).le.0.d0) iok = iok + 1
  enddo
  if (iok.ne.0) then
    write(nfecra,9001)
    call csexit (1)
    !==========
  endif

endif

!===============================================================================
! 3.  GESTION DU TEMPS QUI PASSE...
!===============================================================================

!-->Gestion du pas de temps Lagrangien

dtp = dtref

!-->Incrementation du TEMPS COURANT LAGRANGIEN

ttclag = ttclag + dtp

!-->Test pour savoir si le domaine contient des particules

if (nbpart.eq.0) goto 20

!-->On enregistre l'element de depart de la particule

do ip = 1,nbpart
  indep(ip) = itepa(ip,jisor)
enddo

!===============================================================================
! 4.  GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
!===============================================================================

! Au premier pas de temps on calcul les gradient avec RTP et
! non RTPA car RTPA = initialisation (gradients nuls)

if ( ntcabs.eq.1 ) then

  call laggra                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   rtp    , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

else

  call laggra                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   rtpa   , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

endif

!===============================================================================
! 4.  Initialisation des variables aleatoires gaussiennes
!===============================================================================

!---> CALCUL DES TIRAGES ALEATOIRES
!     remarque : NORMALEN est dans le fichier ZUFALL.F
!     ^^^^^^^^

if (idistu.eq.1) then
  do ivf = 1,nvgaus
    call normalen(nbpart, vagaus(1,ivf))
  enddo
else
  do ivf = 1,nvgaus
    do ip = 1,nbpmax
      vagaus(ip,ivf) = 0.d0
    enddo
  enddo
endif

!---> CALCUL DES TIRAGES ALEATOIRES POUR LE MVT BROWNIEN

if ( lamvbr .eq. 1 ) then

  do ivf = 1,nbrgau
    call normalen(nbpart, brgaus(1,ivf))
  enddo

endif

!===============================================================================
! 5. PROGRESSION DES PARTICULES
!===============================================================================

 10   continue

nor = mod(nor,nordre)
nor = nor + 1

!---> Recopie des resultats de l'etape precedente :

if (nor.eq.1) then

  do ivf = 1,nvp
    do ip = 1,nbpart
      ettpa(ip,ivf) = ettp(ip,ivf)
    enddo
  enddo

endif

!-----> CALCUL GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE
!       EN N+1 (avec RTP)

if (nor.eq.2 .and. iilagr.ne.3) then

  call laggra                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   ia     ,                                                       &
   rtp    , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

endif

!-----> CALCUL DES CARACTERISTIQUES DES PARTICULES

if (nor.eq.1) then

!      sous pas de temps n (avec RTPA)

  call lagcar                                                     &
  !==========
   ( idebia , idebra ,                                            &
     nvar   , nscal  , nphas  ,                                   &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     itepa  ,                                                     &
     ia     ,                                                     &
     dt     , rtpa   , propce , propfa , propfb ,                 &
     ettp   , ettpa  , tepa   , taup   , tlag   ,                 &
     piil   , bx     , tempct , statis ,                          &
     gradpr , gradvf , w1     , w2     , auxl(1,1)  ,             &
     ra     )

else

!     sous pas de temps n+1 (avec RTP)

  call lagcar                                                     &
  !==========
   ( idebia , idebra ,                                            &
     nvar   , nscal  , nphas  ,                                   &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     itepa  ,                                                     &
     ia     ,                                                     &
     dt     , rtp    , propce , propfa , propfb ,                 &
     ettp   , ettpa  , tepa   , taup   , tlag   ,                 &
     piil   , bx     , tempct , statis ,                          &
     gradpr , gradvf , w1     , w2     , auxl(1,1) ,              &
     ra     )

endif


!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     POSITION, VITESSE FLUIDE, VITESSE PARTICULE

call lagesp                                                       &
!==========
   ( idebia , idebra ,                                            &
     nvar   , nscal  , nphas  ,                                   &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     itepa  , ibord  ,                                            &
     ia     ,                                                     &
     dt     , rtpa   , rtp    , propce , propfa , propfb ,        &
     ettp   , ettpa  , tepa   ,                                   &
     statis , stativ , taup   , tlag   , piil   ,                 &
     tsuf   , tsup   , bx     , tsfext ,                          &
     vagaus , gradpr , gradvf , brgaus , terbru ,                 &
     auxl(1,1) , auxl2 ,                                          &
     ra     )

!---> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     LIEES AUX PHYSIQUES PARTICULIERES PARTICULAIRES

if ( iphyla.eq.1 .or. iphyla.eq.2 ) then

  if ( nor.eq.1 ) then
    call lagphy                                                   &
    !==========
    ( idebia , idebra ,                                           &
      nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      ntersl , nvlsta , nvisbr ,                                  &
      itepa  , ibord  ,                                           &
      ia     ,                                                    &
      dt     , rtpa   , propce , propfa , propfb ,                &
      ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,       &
      tsvar  , auxl   , cpgd1  , cpgd2  , cpght  ,                &
      w1     , w2     , w3     ,                                  &
      ra     )
  else
    call lagphy                                                   &
    !==========
    ( idebia , idebra ,                                           &
      nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      ntersl , nvlsta , nvisbr ,                                  &
      itepa  , ibord  ,                                           &
      ia     ,                                                    &
      dt     , rtp    , propce , propfa , propfb ,                &
      ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,       &
      tsvar  , auxl   , cpgd1  , cpgd2  , cpght  ,                &
      w1     , w2     , w3     ,                                  &
      ra     )
  endif

endif

!===============================================================================
! 6.  Couplage Retour - Calcul des termes sources
!===============================================================================

if (iilagr.eq.2 .and. nor.eq.nordre) then

  ifinia = idebia
  iitslg = idebra
  ifinra = iitslg + ntersl*nbpmax
  call rasize('lagune',ifinra)
  !==========

  call lagcou                                                     &
  !==========
   ( ifinia , ifinra ,                                            &
     nvar   , nscal  , nphas  ,                                   &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     itepa  , indep  , ibord  ,                                   &
     ia     ,                                                     &
     rtp    , propce ,                                            &
     ettp   , ettpa  , tepa   , taup   ,                          &
     tempct , tsfext , tslagr ,                                   &
     cpgd1  , cpgd2  , cpght  ,                                   &
     ra(iitslg)      , w1     , w2   ,                            &
     auxl(1,1) , auxl(1,2)   , auxl(1,3) ,                        &
     ra     )

endif

!===============================================================================
! 7.  Reperage des particules - Traitement des conditions aux limites
!     pour la position des particules
!===============================================================================

if (nor.eq.1) then

  call lagcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iitypf)      , ia(iitrif)      ,                            &
   icocel , itycel , ifrlag , itepa  , ibord  , indep  ,          &
   ia     ,                                                       &
   ra(isrfbn)      ,                                              &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , parbor , auxl   ,                   &
   ra     )

  if (ierr.eq.1) then
    call lagerr
    !==========
    goto 20
  endif

endif

!===============================================================================
! 9.  ELIMINATION DES PARTICULES QUI SONT SORTIES DU DOMAINE
!===============================================================================

!     ATTENTION : NBPOUT contient les particules sorties de facon
!                 normal + les particules sorties en erreur de reperage.

if (nor.eq.nordre) then

  call lageli                                                     &
  !==========
 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   npars  ,                                                       &
   itepa  ,                                                       &
   ia     ,                                                       &
   dnpars ,                                                       &
   ettp   , ettpa  , tepa   ,                                     &
   ra     )

  nbpout = npars
  dnbpou = dnpars

endif

!===============================================================================
! 10.  TEMPS DE SEJOUR
!===============================================================================

if (nor.eq.nordre) then

  do npt = 1,nbpart
    if ( itepa(npt,jisor).gt.0 ) then
      tepa(npt,jrtsp) = tepa(npt,jrtsp) + dtp
    endif
  enddo

endif

!===============================================================================
! 11.  CALCUL STATISTIQUES
!===============================================================================

if (nor.eq.nordre .and. istala.eq.1 .and. iplas.ge.idstnt) then

  call lagsta                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , tepa   , statis , stativ ,                            &
   w1     ,                                                       &
   ra     )

endif

!===============================================================================
! 12.  Equation de Poisson
!===============================================================================

if (nor.eq.nordre .and. ilapoi.eq.1) then

  call lagpoi                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , statis ,                                     &
   w1     , w2     , w3     ,                                     &
   ra     )

endif

!===============================================================================
! 13.  Methode de reduction de variances : Clonage/Fusion des particules
!===============================================================================

if ( nor.eq.nordre .and. iroule.ge.1 ) then

  call lagrus                                                     &
  !==========
   ( idebia , idebra ,                                            &
     ncelet , ncel   ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     itepa  , indep  ,                                            &
     ia     ,                                                     &
     ettp   , ettpa  , tepa   , croule ,                          &
     ra     )

  if (npclon.gt.0) then

    npar1 = nbpart - npclon + 1
    npar2 = nbpart

    call lagipn                                                   &
    !==========
    ( idebia , idebra ,                                           &
      ncelet , ncel   ,                                           &
      nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      npar1  , npar2  ,                                           &
      itepa  ,                                                    &
      ia     ,                                                    &
      rtp    ,                                                    &
      ettp   , tepa   , vagaus ,                                  &
      w1     , w2     , w3     ,                                  &
      ra     )

  endif

endif

!===============================================================================
! 14. UN AUTRE TOUR ?
!===============================================================================

if (nordre.eq.2 .and. nor.eq.1) goto 10

!===============================================================================
! 15. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
!     EN FIN D'ITERATION LAGRANGIENNE
!===============================================================================

call uslast                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   statis , stativ ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! 16. Visualisations
!===============================================================================

 20   continue

nfin = 0

!-->Stockage des trajectoires au format Ensight Gold

if (iensi1.eq.1) then

  iforce = 0

  call enslag                                                     &
  !==========
   ( idebia , idebra ,                                            &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nfin   , iforce ,                                            &
     itepa  ,                                                     &
     ettp   , tepa   , ra )
endif

if (iensi2.eq.1) then
  call enswaf                                                     &
  !==========
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     nfin   ,                                                     &
     itepa  ,                                                     &
     ettp   , tepa , auxl   )
endif

!===============================================================================
! 17. NOMBRE DE PARITICULES PERDUES (SUITES COMPRISES)
!===============================================================================

nbpert = nbpert + nbperr

!===============================================================================
! 18. ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
!        - nombre de particules dans le domaine
!        - nombre de particules entrantes
!        - nombre de particules sorties
!        - ...

!===============================================================================

if (ipass.eq.1) then
   modntl = 0
elseif(ntlal.gt.0) then
   modntl = mod(ntcabs,ntlal)
elseif(ntlal.eq.-1.and.ntcabs.eq.ntmabs) then
   modntl = 0
else
   modntl = 1
endif

if (modntl.eq.0) then
   call lagaff                                                    &
   !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , statis , &
   w1     , w2     , w3     ,                                     &
   ra     )

endif

!===============================================================================

!--------
! FORMATS
!--------

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LA TECHNIQUE DE CLONAGE/FUSION DES PARTICULES           ',/,&
'@      EST ENCLENCHEE AVEC UNE FONCTION D''IMPORTANCE        ',/,&
'@      COMPORTANT DES VALEURS NEGATIVES OU NULLES            ',/,&
'@      (LAGUNE).                                             ',/,&
'@                                                            ',/,&
'@    LES ELEMENTS DU TABLEAU CROULE DOIVENT STRICTEMENT      ',/,&
'@      POSITIFS.                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les valeurs de CROULE dans la subroutine USLARU. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
