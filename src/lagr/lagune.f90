!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dlgeo  ,                                                       &
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
! dlgeo            ! tr ! --> ! tableau contenant les donnees geometriques     !
! (nfabor,ngeol)   !    !     ! pour le sous-modele de depot                   !
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
use parall
use period
use pointe
use lagpar
use lagran
use mesh
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          lndnod
integer          nvar   , nscal
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
double precision dlgeo(nfabor,ngeol)
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
integer          modntl , iromf

double precision dnpars

integer          ifac , nn , ifab , ifap , kfap
integer          n10,n20,n30,n50,n100,nmax
integer          ius

double precision distp , d1 , px,py,pz, lvisq, visccf, romf
double precision tvisq, ustar, ustarmoy

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
nbpdep = 0

dnbpnw = 0.d0
dnpcsu = 0.d0
dnpclo = 0.d0
dnpkil = 0.d0
dnpenc = 0.d0
dnbpou = 0.d0
dnbper = 0.d0
dnbdep = 0.d0

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

!
! --> if the deposition model is activated
!

  if (idepst.ge.1) then

     ustarmoy = 0.d0
     ius = 0

    ! boundary faces data

     call laggeo                                                  &
     !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   ia     , dlgeo  , ra     )

    ! the mesh elements yplus checking

     n10  = 0
     n20  = 0
     n30  = 0
     n50  = 0
     n100 = 0
     nmax = 0

     do ifac=1, nfabor

       if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

         distp = 0.d0
         iel = ifabor(ifac)

      ! the density pointer according to the flow location

         if ( ippmod(icp3pl).ge.0 .or. ippmod(icfuel).ge.0 ) then
           iromf = ipproc(irom1)
         else
           iromf = ipproc(irom)
         endif

         romf = propce(iel,iromf)
         visccf = propce(iel,ipproc(iviscl)) / romf

         do kfap = itycel(iel), itycel(iel+1)-1

           ifap = icocel(kfap)

           if (ifap.gt.0) then

             do nn = ipnfac(ifap), ipnfac(ifap+1)-1

               px = xyznod(1,nodfac(nn))
               py = xyznod(2,nodfac(nn))
               pz = xyznod(3,nodfac(nn))
               d1 = abs( px*dlgeo(ifac,1)+py*dlgeo(ifac,2)            &
                    +pz*dlgeo(ifac,3)+   dlgeo(ifac,4) )              &
                    /sqrt( dlgeo(ifac,1)*dlgeo(ifac,1)                &
                    +dlgeo(ifac,2)*dlgeo(ifac,2)                      &
                    +dlgeo(ifac,3)*dlgeo(ifac,3) )

               if ( d1 .gt. distp ) then
                 distp = d1
               endif

             enddo

           else

             ifab = -ifap

             do nn = ipnfbr(ifab), ipnfbr(ifab+1)-1

               px = xyznod(1,nodfbr(nn))
               py = xyznod(2,nodfbr(nn))
               pz = xyznod(3,nodfbr(nn))

               d1 = abs( px*dlgeo(ifac,1)+py*dlgeo(ifac,2)           &
                        +pz*dlgeo(ifac,3)+ dlgeo(ifac,4))            &
                  /sqrt( dlgeo(ifac,1)*dlgeo(ifac,1)                 &
                       + dlgeo(ifac,2)*dlgeo(ifac,2)                 &
                       + dlgeo(ifac,3)*dlgeo(ifac,3))

               if ( d1.gt.distp) then
                 distp = d1
               endif

             enddo

           endif

         enddo

         ustar = uetbor(ifac)

         if (ustar.gt.0.d0) then

           ustarmoy = ustarmoy + ustar
           ius = ius + 1

           lvisq = visccf / ustar


           distp = distp/lvisq

           if ( distp .le. 10.d0 ) then
             n10 = n10+1
           else if ( distp .le. 20.d0 ) then
             n20 = n20+1
           else if ( distp .le. 30.d0 ) then
             n30 = n30+1
           else if ( distp .le. 50.d0 ) then
             n50 = n50+1
           else if ( distp .le. 100.d0 ) then
             n100 = n100+1
           else
             nmax = nmax +1
           endif

         endif

       endif

     enddo

     ustarmoy = ustarmoy / ius

! the mesh edge yplus and average friction velocity display

     write(nfecra,4100) nfabor,n10,n20,n30,n50,n100,nmax,ustarmoy
!
  endif

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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel , icocel ,                                              &
   itypfb , itrifb , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , auxl   , w1     , w2     , w3     , &
   ra     )

else

  call lagent                                                     &
  !==========
 ( idebia , idebra ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itycel , icocel ,                                              &
   itypfb , itrifb , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
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
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa ,                                      &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , croule , auxl ,                     &
   dispar , yplpar ,                                              &
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
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   rtp    , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   ra     )

else

  call laggra                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   rtpa   , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
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
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   rtp    , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   ra     )

endif

!-----> CALCUL DES CARACTERISTIQUES DES PARTICULES

if (nor.eq.1) then

!      sous pas de temps n (avec RTPA)

  call lagcar                                                     &
  !==========
   ( idebia , idebra ,                                            &
     nvar   , nscal  ,                                            &
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
     nvar   , nscal  ,                                            &
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
     nvar   , nscal  , lndnod ,                                   &
     nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ntersl , nvlsta , nvisbr ,                                   &
     icocel , itycel , ifrlag,                                    &
     itepa  , ibord  ,                                            &
     ia     ,                                                     &
     dlgeo  ,                                                     &
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
     nvar   , nscal  ,                                            &
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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb ,                                              &
   icocel , itycel , ifrlag , itepa  , ibord  , indep  ,          &
   ia     ,                                                       &
   dlgeo  ,                                                       &
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
   nvar   , nscal  ,                                              &
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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , statis ,                                     &
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
    ( ncelet , ncel   ,                                           &
      nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
      npar1  , npar2  ,                                           &
      itepa  ,                                                    &
      ia     ,                                                    &
      rtp    ,                                                    &
      ettp   , tepa   , vagaus ,                                  &
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
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   statis , stativ ,                                              &
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
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , statis , &
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

 4100 format(                                                     &
'                                                               '/,&
'   ** LAGRANGIAN MODULE:  '                                     /,&
'   ** Check of the mesh for the deposition submodel  '         ,/,&
'      ---------------------------------------------  '         ,/,&
'                                                               '/,&
' Number of boundary faces                        ',I10         ,/,&
' Number of boundary faces with 0  < y^+ < 10     ',I10         ,/,&
' Number of boundary faces with 10 < y^+ < 20     ',I10         ,/,&
' Number of boundary faces with 20 < y^+ < 30     ',I10         ,/,&
' Number of boundary faces with 30 < y^+ < 50     ',I10         ,/,&
' Number of boundary faces with 50 < y^+ < 100    ',I10         ,/,&
' Number of boundary faces with y^+ > 100         ',I10         ,/,&
'                                                               '/,&
'   ** Mean friction velocity  (ustar) =  ',F7.3                ,/,&
'---------------------------------------------------------------  ',/)

!----
! FIN
!----

end subroutine
