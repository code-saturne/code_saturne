!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine lagich &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , volume ,          &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   cpgd1  , cpgd2  , cpght  ,                                     &
   skp1   , skp2   , skglob ,                                     &
   gamhet , deltah , tempf  ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LE CHARBON

!        - Temperature              (JHP)
!        - Masse de charbon reactif (JMCH)
!        - Masse de coke            (JMCK)

!     ET CALCUL DU DIAMETRE DU COEUR RETRECISSANT (JRDCK)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! cpgd1,cpgd2,     ! tr ! --> ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpmax    !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
! sk1,sk2,         ! tr ! --- ! tableaux de travail                            !
! skglob(nbpmax    !    !     !                                                !
! gamhet(nbpmax    ! tr ! --- ! tableau de travail                             !
! deltah(nbpmax    ! tr ! --- ! tableau de travail                             !
! tempf(ncelet)    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "cstphy.h"
include "cstnum.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "cpincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp , nvp1 , nvep , nivep
integer          ntersl , nvlsta , nvisbr
integer          itepa(nbpmax,nivep) , ibord(nbpmax)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision skp1(nbpmax) , skp2(nbpmax) , skglob(nbpmax)
double precision gamhet(nbpmax) , deltah(nbpmax) , tempf(ncelet)
double precision cpgd1(nbpmax), cpgd2(nbpmax), cpght(nbpmax)
double precision ra(*)

! VARIABLES LOCALES

integer          npt , iel , icha , mode , ige
integer          iromf , iphas
double precision aux1 , aux2 , aux3 , aux4 , aux5 , aux6
double precision ter1 , ter2 , ter3 , diamp2, d2
double precision tpk , tfk , skc , skdd , se , po2
double precision ho2tf , hctp , hcotp , den , sherw
double precision coef , mp0 , d6spi , dpis6 , d1s3 , d2s3
double precision gamdv1(ncharm2) , gamdv2(ncharm2)
double precision f1mc(ncharm2) , f2mc(ncharm2)
double precision coefe(ngazem)

double precision precis
parameter ( precis = 1.d-15 )

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

iphas = ilphas

d6spi = 6.d0 / pi
dpis6 = pi / 6.d0
d1s3  = 1.d0 / 3.d0
d2s3  = 2.d0 / 3.d0

! --- Si couplage retour thermique :

if ( ltsthe.eq.1 .and. nor.eq.1 ) then
  do npt = 1,nbpart
    cpgd1(npt) = 0.d0
    cpgd2(npt) = 0.d0
    cpght(npt) = 0.d0
  enddo
endif
if ( ltsthe.eq.1 ) then
  coef = 1.d0 / dble(nordre)
endif

!===============================================================================
! 2. Pointeur sur la masse volumique en fonction de l'ecoulement
!===============================================================================

if ( ippmod(icp3pl).ge.0 ) then
  iromf = ipproc(irom1)
else
  iromf = ipproc(irom(iphas))
endif

!===============================================================================
! 3. Temperature moyenne Fluide en Kelvin dans le cas ou la phase
!    porteuse est une flamme de charbon pulverise
!===============================================================================

if ( ippmod(icp3pl).ge.0 .or. ippmod(icpl3c).ge.0 ) then

   do iel = 1,ncel
     tempf(iel) = propce(iel,ipproc(itemp1))
   enddo

else
  write(nfecra,1000) iphyla, ippmod(icpl3c), ippmod(icp3pl)
  call csexit (1)
  !==========
endif

!===============================================================================
! 4. Calcul des constantes de vitesses SPK1 et SPK2 du transfert
!    de masse par devolatilisation avec des lois d'Arrhenius
!===============================================================================

!     RR --> Constante des gaz parfaits en J/mol/K

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    tpk = ettp(npt,jhp) + tkelvi
    aux1 = 1.d0 / (rr*tpk)
    skp1(npt) = a1ch(icha) * exp( -e1ch(icha) * aux1)
    skp2(npt) = a2ch(icha) * exp( -e2ch(icha) * aux1)
  endif
enddo

!===============================================================================
! 5. Calcul de la masse volumique du coke
!    On suppose pour le calcul de la masse volumique du coke que
!    la devolatilisation a lieu a volume constant
!===============================================================================

! --- Initialisation

do icha = 1,ncharm
  gamdv1(icha) = 0.d0
  gamdv2(icha) = 0.d0
  rhock(icha) = rho0ch(icha)
enddo

! --- Calcul de l'integrale de GAMDV1 et GAMDV2 pour chaque charbon

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    aux1 = skp1(npt) * y1ch(icha) * ettp(npt,jmch)
    aux2 = skp2(npt) * y2ch(icha) * ettp(npt,jmch)

    gamdv1(icha) = gamdv1(icha) + aux1
    gamdv2(icha) = gamdv2(icha) + aux2

! --- Couplage retour thermique

    if ( ltsthe.eq.1 ) then
      cpgd1(npt) = cpgd1(npt) + coef*aux1
      cpgd2(npt) = cpgd2(npt) + coef*aux2
    endif

  endif
enddo

! --- Calcul de la masse volumique moyenne du coke

do icha = 1,ncharb
  den = y2ch(icha)*gamdv1(icha) + y1ch(icha)*gamdv2(icha)
  if ( den.gt.precis ) then
    rhock(icha) = rho0ch(icha)                                    &
      *( y2ch(icha)*gamdv1(icha)+y1ch(icha)*gamdv2(icha)          &
        -y1ch(icha)*y2ch(icha)*(gamdv1(icha)+gamdv2(icha)) )      &
       / den
  endif
enddo

!===============================================================================
! 6. Calcul du diametre du coeur retrecissant
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)

    tepa(npt,jrdck) =                                             &
             ( (d6spi / ( 1.d0-xashch(icha)) )                    &
              *( ettp(npt,jmch)/rho0ch(icha)                      &
              +ettp(npt,jmck)/rhock(icha) ) )**d1s3
  endif
enddo

!===============================================================================
! 7. Calcul de la constante globale de reaction
!===============================================================================

! ---  Hypothese Sherwood = 2.

sherw = 2.d0

do npt = 1,nbpart

  if (itepa(npt,jisor).gt.0) then

    icha = itepa(npt,jinch)

    tpk = ettp(npt,jhp) + tkelvi

! --- Coefficient de cinetique chimique de formation de CO
!       en (kg.m-2.s-1.atm(-n))

    skc = ahetch(icha)                                            &
      * exp(-ehetch(icha)*4185.d0 / (rr*tpk) )

! --- Coefficient de diffusion en  (Kg/m2/s/atm) et constante
!     globale de reaction

    if ( tepa(npt,jrdck).gt.epsicp ) then
      skdd = sherw * 2.53d-7 * (tpk**0.75d0) / tepa(npt,jrdck)
      skglob(npt) = (skc*skdd) / (skc+skdd)
    else
      skglob(npt) = skc
    endif

  endif
enddo

!===============================================================================
! 8. Calcul de la GAMMAhet , GAMMACH et 0.5(MO2/MC)*(HO2(Tp)-HO2(TF))
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then

    icha = itepa(npt,jinch)
    iel  = itepa(npt,jisor)

    tpk = ettp(npt,jhp) + tkelvi
    tfk = ettp(npt,jtf) + tkelvi

! --- Calcul de la pression partielle en oxygene (atm)
!                                                 ---
!       PO2 = RHO1*RR*T*YO2/MO2

    po2 = propce(iel,iromf) * rr * tempf(iel)                     &
        * propce(iel,ipproc(iym1(io2))) / wmole(io2) / prefth

! --- Calcul de (Surface efficace)/(Mck**2/3) : SE

    se =  ( pi*(1.d0-xashch(icha)) )**d1s3                        &
        * ( 6.d0/rhock(icha)       )**d2s3

! --- Calcul de la GamHET/(Mck**2/3)

    gamhet(npt) = se * po2 * skglob(npt)

! --- Pas de combustion heterogene si Mch/Mp >= 1.D-3

    if ( ettpa(npt,jmch).ge.(1.d-3*ettpa(npt,jmp)) ) then
      gamhet(npt) = 0.d0
    endif

! --- Couplage retour thermique

    if ( ltsthe.eq.1 ) then
      cpght(npt) = cpght(npt)                                     &
                 + coef * gamhet(npt)                             &
                 * ( ettp(npt,jmck)**d2s3 )
    endif

! --- Calcul de Hc(Tp)-Mco/Mc Hco2(Tp)+0.5Mo2/Mc Ho2(Tf)

!        Calcul de Hcoke(TP)

    hctp = h02ch(icha) + ettp(npt,jcp)*(tpk-trefth)

!        Calcul de MCO/MC HCO(TP)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(ico) = wmole(ico) / wmolat(iatc)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    mode      = -1
    call cpthp1 ( mode , hcotp , coefe  , f1mc , f2mc ,  tpk )
    !==========

!        Calcul de MO2/MC/2. HO2(TF)

    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = wmole(io2) / wmolat(iatc) / 2.d0

    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    mode      = -1
    call cpthp1 ( mode  , ho2tf , coefe , f1mc  , f2mc , tfk )
    !==========

    deltah(npt) = hcotp - ho2tf - hctp

  endif

enddo

!===============================================================================
! 9. Integration Masse de Charbon reactif
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      aux1 = exp(-(skp1(npt)+skp2(npt))*dtp)
      ettp(npt,jmch) = ettpa(npt,jmch)*aux1

! Clipping
      if ( ettp(npt,jmch).lt.precis ) then
        ettp(npt,jmch) = 0.d0
      endif

    endif
  enddo

else if (nor.eq.2) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0 .and. ibord(npt).eq.0) then

      aux1 = exp(-(skp1(npt)+skp2(npt))*dtp)
      ettp(npt,jmch) = 0.5d0 * ( ettp(npt,jmch)                   &
                         + ettpa(npt,jmch)*aux1 )

! Clipping
      if ( ettp(npt,jmch).lt.precis ) then
        ettp(npt,jmch) = 0.d0
      endif

    endif
  enddo
endif

!===============================================================================
! 10. Integration Masse de Coke
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      icha = itepa(npt,jinch)

      aux1 = skp1(npt) * (1.d0-y1ch(icha)) * ettpa(npt,jmch)      &
           + skp2(npt) * (1.d0-y2ch(icha)) * ettpa(npt,jmch)

      if ( ettpa(npt,jmck).gt.precis ) then

        aux2 = ettpa(npt,jmck)**d1s3
        aux3 = aux2 + d2s3 *gamhet(npt) *dtp

        ter1 = (aux1*aux2-gamhet(npt)*ettpa(npt,jmck)) *dtp/aux3

        tsvar(npt,jmck) = 0.5d0 * ter1

        ettp(npt,jmck) = ettpa(npt,jmck) + ter1

      else
        ter1 = aux1 * dtp
        tsvar(npt,jmck) = 0.5d0 * ter1
        ettp(npt,jmck) = ettpa(npt,jmck) + ter1
      endif
!  Clipping
      if ( ettp(npt,jmck).lt.0.d0 ) then
        ettp(npt,jmck) = 0.d0
      endif

    endif
  enddo

else if (nor.eq.2) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0 .and. ibord(npt).eq.0) then

      icha = itepa(npt,jinch)

      aux1 = skp1(npt) * (1.d0-y1ch(icha)) * ettp(npt,jmch)       &
           + skp2(npt) * (1.d0-y2ch(icha)) * ettp(npt,jmch)

      if ( ettpa(npt,jmck).gt.precis ) then

        aux2 = ettpa(npt,jmck)**d1s3
        aux3 = aux2 + d2s3 *gamhet(npt) *dtp

        ter1 = ( aux1 *aux2 -gamhet(npt) *ettpa(npt,jmck))        &
             * dtp / aux3

        ettp(npt,jmck) = ettpa(npt,jmck)                          &
                            +tsvar(npt,jmck)+0.5d0*ter1

      else
        ter1 = aux1*dtp
        ettp(npt,jmck) = ettpa(npt,jmck)                          &
                           + tsvar(npt,jmck) + 0.5d0*ter1
      endif
! Clipping
      if ( ettp(npt,jmck).lt.0.d0 ) then
        ettp(npt,jmck) = 0.d0
      endif

    endif
  enddo
endif

!===============================================================================
! 11. Integration de la temperature des grains de charbon
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      icha = itepa(npt,jinch)
      iel   = itepa(npt,jisor)

      d2 = ettp(npt,jdp)*ettp(npt,jdp)
      diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)       &
              +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)

      aux1 = tempct(npt,1)*diamp2/d2

!    Combustion heterogene
      aux2 = ( -gamhet(npt) *(ettp(npt,jmck)**d2s3) *deltah(npt) )

!    Rayonnement

      aux3 = pi * diamp2                                          &
           * ( propce(iel,ipproc(ilumi))                          &
           - 4.d0*stephn*((ettp(npt,jhp)+tkelvi)**4) )

      aux4 = ettpa(npt,jtf) + aux1*(aux2+aux3)

      aux5 = dtp/aux1
      aux6 = exp(-aux5)

      ter1 = ettpa(npt,jhp) * aux6
      ter2 = aux4 * (1.d0-aux6)
      ter3 = aux4 * ( -aux6+(1.d0-aux6) / aux5 )

      tsvar(npt,jhp) = 0.5d0 * ter1 + ter3
      ettp(npt,jhp) = ter1 + ter2

    endif
  enddo

else if (nor.eq.2) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0 .and. ibord(npt).eq.0) then


      icha  = itepa(npt,jinch)
      iel   = itepa(npt,jisor)

      d2 = ettp(npt,jdp)*ettp(npt,jdp)
      diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)       &
             + (1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)

      aux1 = tempct(npt,1)*diamp2/d2

!    Combustion heterogene
      aux2 = -gamhet(npt) *(ettp(npt,jmck)**d2s3) *deltah(npt)

!    Rayonnement

      aux3 = pi * diamp2                                          &
           * ( propce(iel,ipproc(ilumi))                          &
           - 4.d0*stephn*((ettp(npt,jhp)+tkelvi)**4) )

      aux4 = ettp(npt,jtf) + aux1*(aux2+aux3)

      aux5 = dtp / aux1
      aux6 = exp(-aux5)

      ter1 = ettpa(npt,jhp) * aux6
      ter2 = aux4 * ( 1.d0-((1.d0-aux6)/aux5) )

      ettp(npt,jhp) = tsvar(npt,jhp) + 0.5d0*ter1 + ter2

    endif
  enddo
endif

!===============================================================================
! 12. Mise a jour du diametre du coeur retrecissant
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)

    tepa(npt,jrdck) = ( (d6spi/(1.d0-xashch(icha)) )              &
                    *   ( ettp(npt,jmch)/rho0ch(icha)             &
                        + ettp(npt,jmck)/rhock(icha) ) )**d1s3
  endif
enddo

!===============================================================================
! 13. Calcul du diametre des grains de charbon
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    ettp(npt,jdp) = (xashch(icha)*(tepa(npt,jrd0p)**2)            &
                         + (1.d0-xashch(icha))                    &
                          *(tepa(npt,jrdck)**2) )**0.5d0
  endif
enddo

!===============================================================================
! 14. Calcul de la masse des grains de charbon
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    mp0  = dpis6 * (tepa(npt,jrd0p)**3) * rho0ch(icha)
    ettp(npt,jmp) = ettp(npt,jmch) + ettp(npt,jmck)               &
                      + xashch(icha)*mp0
  endif
enddo

!===============================================================================

!=======
! FORMAT
!=======

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON        ',/,&
'@      EST ACTIVE (LAGICH), ALORS QU''AUCUNE PHYSIQUE        ',/,&
'@      PARTICULIERE SUR LA COMBUSTION DU CHABON PULVERISE    ',/,&
'@      N''EST PAS ENCLENCHE (USPPMO).                        ',/,&
'@                                                            ',/,&
'@       IPHYLA = ', I10                                       ,/,&
'@       IPPMOD(ICPL3C) = ', I10                               ,/,&
'@       IPPMOD(ICP3PL) = ', I10                               ,/,&
'@                                                            ',/,&
'@  Le transport Lagrangien de particule de charbon doit      ',/,&
'@   etre couple avec la combustion d''une flamme de charbon  ',/,&
'@   pulverise en phase continue.                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,&
'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end
