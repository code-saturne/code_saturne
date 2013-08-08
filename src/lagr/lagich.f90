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

subroutine lagich &
!================

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  ,                                              &
   dt     , rtp    , propce ,                                     &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   cpgd1  , cpgd2  , cpght  ,                                     &
   skp1   , skp2   , skglob )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LE CHARBON

!        - Temperature              (JHP)
!        - Masse d eau              (JMWAT)
!        - Masse de charbon reactif (JMCH)
!        - Masse de coke            (JMCK)

!     ET CALCUL DU DIAMETRE DU COEUR RETRECISSANT (JRDCK)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp , nvp1 , nvep , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep) , ibord(nbpmax)

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision cpgd1(nbpmax), cpgd2(nbpmax), cpght(nbpmax)
double precision skp1(nbpmax) , skp2(nbpmax) , skglob(nbpmax)


! Local variables

integer          npt , iel , icha , mode , ige
integer          iromf
double precision aux1 , aux2 , aux3 , aux4 , aux5 , aux6
double precision ter1 , ter2 , ter3 , diamp2, dd2
double precision lv, tebl, tlimit, tmini, tsat, fwatsat
double precision tpk , tfk , skc , skdd , se , po2
double precision ho2tf , hctp , hcotp , den , sherw
double precision coef , mp0 , d6spi , dpis6 , d1s3 , d2s3, mv
double precision gamdv1(ncharm2) , gamdv2(ncharm2)
double precision f1mc(ncharm2) , f2mc(ncharm2)
double precision coefe(ngazem)

double precision, allocatable, dimension(:) :: tempf
double precision, allocatable, dimension(:) :: fwat, gamhet, deltah

double precision precis
parameter ( precis = 1.d-15 )

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(tempf(ncelet))
allocate(fwat(nbpmax) , gamhet(nbpmax) , deltah(nbpmax))

! Initialize variables to avoid compiler warnings

coef = 0.d0

d6spi = 6.d0 / pi
dpis6 = pi / 6.d0
d1s3  = 1.d0 / 3.d0
d2s3  = 2.d0 / 3.d0

! --- Si couplage retour thermique :

if ( ltsthe.eq.1 .and. nor.eq.1 ) then

  ! Le couplage thermique n'est plus valable depuis l'implementation du sechage
  write(nfecra,1001) iphyla,ltsthe
  call csexit (1)

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

if ( ippmod(iccoal).ge.0 ) then
  iromf = ipproc(irom1)
else
  iromf = ipproc(irom)
endif

!===============================================================================
! 3. Temperature moyenne Fluide en Kelvin dans le cas ou la phase
!    porteuse est une flamme de charbon pulverise
!===============================================================================

if ( ippmod(iccoal).ge.0 .or. ippmod(icpl3c).ge.0 ) then

   do iel = 1,ncel
     tempf(iel) = propce(iel,ipproc(itemp1))
   enddo

else
  write(nfecra,1000) iphyla, ippmod(icpl3c), ippmod(iccoal)
  call csexit (1)
  !==========
endif

!===============================================================================
! 4. Calcul de la masse d'eau qui s'evapore
!    On suppose pour le calcul de la masse volumique du charbon actif que
!    le sechage a lieu a volume constant
!===============================================================================

! --- Initialisation
do npt = 1,nbpart
  fwat(npt) = 0.d0
enddo

!      Chaleur Latente en J/kg
lv = 2.263d+6
!      Temperature d'ebulition de l'eau
tebl = 100.d0 + tkelvi
!      Temperature limite
tlimit = 302.24d0
!      Temperature mini (apres, la fraction massique d'eau saturante est nulle)
tmini = tlimit*(1.d0-tlimit*rr/(lv*wmole(ih2o)))
!      Nombre de Sherwood (fixé à 2)
sherw=2.0d0

! On prend en compte l'humidité du charbon
do npt = 1,nbpart
! --- Calcul de la fraction massique d'eau saturante
  if (itepa(npt,jisor).gt.0) then
    iel  = itepa(npt,jisor)
    tpk = ettp(npt,jhp)
    if (tpk.ge.tmini) then
      if (tpk.ge.tlimit) then
        aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
        aux2 = aux1 * exp( lv * wmole(ih2o) * (1.0d0/tebl - 1.0d0/tpk) / rr )
      else
        ! On linearise la fraction massique d'eau saturante entre tmini et Tlimit
        ! En Tlimit, la fraction massique d'eau saturante est nulle
        aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
        aux2 = aux1 * exp( lv * wmole(ih2o) * (1.0d0/tebl - 1.0d0/tlimit) / rr ) &
               * (lv*wmole(ih2o) / (rr*tlimit**2)) * (tpk - tmini)
      endif
! --- Calcul du terme source d eau diffusee
      aux3 = max(1.0d0 - aux2, precis)
      aux4 = pi*tepa(npt,jrd0p)*diftl0*sherw*                                    &
             log((1.0d0-propce(iel,ipproc(iym1(ih2o))))/aux3)
    else
      ! Le flux est nul
      aux1 = wmole(ih2o) / propce(iel,ipproc(immel))
      aux4 = 0.0d0
    endif

! --- Verification clipping
    ! Limitation du flux par rapport à qte d'humidite encore presente
    fwat(npt) = min(ettp(npt,jmwat)/dtp,aux4)

    ! Limitation du flux par rapport à la température de saturation
    ! On limite le flux de sechage pour que, à la fin d'un pas de temps,
    ! l'enthalpie de la particule soit suffisament élevée pour que sa pression
    ! saturante en eau soit superieure a la pression partielle d'eau dans l'air
    ! qui l'entoure

    ! Calcul de tsat, temperature saturante à la fraction partielle de l'air
    if (propce(iel,ipproc(iym1(ih2o))) .gt. precis) then
      tsat = 1 / (1/tebl - rr*log(propce(iel,ipproc(iym1(ih2o)))/aux1)         &
                   /(lv * wmole(ih2o)) )
      if (tsat .lt. tlimit) then
        tsat = tmini + propce(iel,ipproc(iym1(ih2o))) / (aux1 *                &
                       exp( lv*wmole(ih2o)*(1.0d0/tebl-1.0d0/tlimit)/rr) *     &
                      (lv*wmole(ih2o)/(rr*tlimit**2)) )
      endif
    else
      tsat = tmini
    endif

    ! On calcule le flux maximum d'evaporation/condensation autorise tel que T(n+1)=Tsat
    icha = itepa(npt,jinch)
    ! Temps caracteristique thermique
    dd2 = ettp(npt,jdp)*ettp(npt,jdp)
    diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)       &
            +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)
    aux1 = tempct(npt,1)*diamp2/dd2
    ! Flux de rayonnement
    aux3 = pi * diamp2                                          &
         * ( propce(iel,ipproc(ilumin))/4.d0                    &
         - stephn*(tpk**4) )
    ! Facteurs exponentiels
    aux5 = dtp/aux1
    aux6 = exp(-aux5)
    ! Flux maximum d'evaporation/condensation
    fwatsat =  ( ( (tpk*aux6-tsat)*ettpa(npt,jmp)*ettpa(npt,jcp) )             &
                    / (aux1*(1.d0-aux6))                                       &
                 + (ettpa(npt,jtf)+tkelvi)*ettpa(npt,jmp)*ettpa(npt,jcp)/aux1  &
                 + aux3 ) / lv

    ! Limitation éventuelle du flux d'evaporation/condensation
    if (fwat(npt) .gt. 0.d0) then
      fwatsat = max(0.d0 , fwatsat)
      fwat(npt) = min( fwatsat , fwat(npt))
    else
      fwatsat = min(0.d0 , fwatsat)
      fwat(npt) = max( fwatsat , fwat(npt))
    endif

  endif
enddo


!===============================================================================
! 5. Calcul des constantes de vitesses SPK1 et SPK2 du transfert
!    de masse par devolatilisation avec des lois d'Arrhenius
!===============================================================================

!     RR --> Constante des gaz parfaits en J/mol/K

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    tpk = ettp(npt,jhp)
    aux1 = 1.d0 / (rr*tpk)
    skp1(npt) = a1ch(icha) * exp( -e1ch(icha) * aux1)
    skp2(npt) = a2ch(icha) * exp( -e2ch(icha) * aux1)
  endif
enddo

!===============================================================================
! 6. Calcul de la masse volumique du coke
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

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    if ( ettpa(npt,jmch).ge.(1.d-3*ettpa(npt,jmp)) ) then
      icha = itepa(npt,jinch)
      mp0  = dpis6 * (tepa(npt,jrd0p)**3) * rho0ch(icha)
      ! mv represente la masse qui a quitté le grain (eau+produits de la devolatilisation)
      mv = mp0 - ettpa(npt,jmch) - ettpa(npt,jmwat) - xashch(icha)*mp0 - ettpa(npt,jmck)
      tepa(npt,jrhock)=rho0ch(icha)- d6spi/(tepa(npt,jrd0p)**3)/(1.d0-xashch(icha))*mv
    endif
  endif
enddo

!===============================================================================
! 7. Calcul du diametre du coeur retrecissant
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    if ( ettpa(npt,jmch).ge.(1.d-3*ettpa(npt,jmp)) ) then
     tepa(npt,jrdck) = tepa(npt,jrd0p)
    else
     tepa(npt,jrdck) =                                             &
             ( (d6spi / ( 1.d0-xashch(icha)) )                    &
              *( ettp(npt,jmch)/rho0ch(icha)                      &
              +ettp(npt,jmck)/tepa(npt,jrhock)))**d1s3
    endif
  endif
enddo

!===============================================================================
! 8. Calcul de la constante globale de reaction
!===============================================================================

! ---  Hypothese Sherwood = 2.

sherw = 2.d0

do npt = 1,nbpart

  if (itepa(npt,jisor).gt.0) then

    icha = itepa(npt,jinch)

    tpk = ettp(npt,jhp)

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
! 9. Calcul de la GAMMAhet , GAMMACH et 0.5(MO2/MC)*(HO2(Tp)-HO2(TF))
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then

    icha = itepa(npt,jinch)
    iel  = itepa(npt,jisor)

    tpk = ettp(npt,jhp)
    tfk = ettp(npt,jtf) + tkelvi

! --- Calcul de la pression partielle en oxygene (atm)
!                                                 ---
!       PO2 = RHO1*RR*T*YO2/MO2

    po2 = propce(iel,iromf) * rr * tempf(iel)                     &
        * propce(iel,ipproc(iym1(io2))) / wmole(io2) / prefth

! --- Calcul de (Surface efficace)/(Mck**2/3) : SE

    se =  ( pi*(1.d0-xashch(icha)) )**d1s3                        &
        * ( 6.d0/   tepa(npt,jrhock)    )**d2s3

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
! 10. Integration Masse d eau
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      aux1 = fwat(npt)*dtp
      ettp(npt,jmwat) = ettpa(npt,jmwat)-aux1

! Clipping
      if ( ettp(npt,jmwat).lt.precis ) then
        ettp(npt,jmwat) = 0.d0
      endif
    endif
  enddo

else if (nor.eq.2 .and. ibord(npt).eq.0) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      aux1 = fwat(npt)*dtp
      ettp(npt,jmwat) = 0.5d0 * ( ettp(npt,jmwat)                   &
                         + ettpa(npt,jmwat)-aux1 )

! Clipping
      if ( ettp(npt,jmwat).lt.precis ) then
        ettp(npt,jmwat) = 0.d0
      endif

    endif
  enddo
endif

!===============================================================================
! 11. Integration Masse de Charbon reactif
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
! 12. Integration Masse de Coke
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      icha = itepa(npt,jinch)

      aux1 = -(skp1(npt) * (1.d0-y1ch(icha))       &
           + skp2(npt) * (1.d0-y2ch(icha)))      &
                       / (skp1(npt)+skp2(npt))

      if ( ettpa(npt,jmck).gt.precis ) then

        aux2 = ettpa(npt,jmck)**d1s3
        aux3 = aux2 + d2s3 *gamhet(npt) * dtp

        ter1 = (aux1*(ettp(npt,jmch)-ettpa(npt,jmch))*aux2-gamhet(npt)*ettpa(npt,jmck)*dtp) /aux3

        tsvar(npt,jmck) = 0.5d0 * ter1

        ettp(npt,jmck) = ettpa(npt,jmck) + ter1

      else

        ter1 = aux1*(ettp(npt,jmch)-ettpa(npt,jmch))

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

      aux1 = -(skp1(npt) * (1.d0-y1ch(icha))      &
             + skp2(npt) * (1.d0-y2ch(icha)))      &
                       / (skp1(npt)+skp2(npt))

      if ( ettpa(npt,jmck).gt.precis ) then

        aux2 = ettpa(npt,jmck)**d1s3
        aux3 = aux2 + d2s3 *gamhet(npt) *dtp


        ter1 = (aux1*(ettp(npt,jmch)-ettpa(npt,jmch))*aux2-gamhet(npt)*ettpa(npt,jmck)*dtp) /aux3

        ettp(npt,jmck) = ettpa(npt,jmck)                          &
                            +tsvar(npt,jmck)+0.5d0*ter1

      else

        ter1 = aux1*(ettp(npt,jmch)-ettpa(npt,jmch))

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
! 13. Integration de la temperature des grains de charbon
!===============================================================================

if (nor.eq.1) then
  do npt = 1,nbpart
    if (itepa(npt,jisor).gt.0) then

      icha = itepa(npt,jinch)
      iel   = itepa(npt,jisor)

      dd2 = ettp(npt,jdp)*ettp(npt,jdp)
      diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)       &
              +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)

      aux1 = tempct(npt,1)*diamp2/dd2

!    Combustion heterogene & sechage
      aux2 = ( -gamhet(npt)*(ettp(npt,jmck)**d2s3)*deltah(npt) )  &
            +( -fwat(npt)                         *lv          )

!    Rayonnement
      aux3 = pi * diamp2                                          &
           * ( propce(iel,ipproc(ilumin))/4.d0                    &
           - stephn*((ettp(npt,jhp))**4) )


      if (ettpa(npt,jmp)*ettpa(npt,jcp) .le. precis) then
        aux4 = (ettpa(npt,jtf)+tkelvi) + aux1*(aux2+aux3)/precis
      else
        aux4 = (ettpa(npt,jtf)+tkelvi) + aux1*(aux2+aux3)/(ettpa(npt,jmp)*ettpa(npt,jcp))
      endif
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

      icha = itepa(npt,jinch)
      iel   = itepa(npt,jisor)

      dd2 = ettp(npt,jdp)*ettp(npt,jdp)
      diamp2 = xashch(icha)*tepa(npt,jrd0p)*tepa(npt,jrd0p)       &
              +(1.d0-xashch(icha))*tepa(npt,jrdck)*tepa(npt,jrdck)

      aux1 = tempct(npt,1)*diamp2/dd2

!    Combustion heterogene & sechage
      aux2 = ( -gamhet(npt)*(ettp(npt,jmck)**d2s3)*deltah(npt) )  &
            +( -fwat(npt)                         *lv          )

!    Rayonnement
      aux3 = pi * diamp2                                          &
           * ( propce(iel,ipproc(ilumin))/4.d0                    &
           - stephn*((ettp(npt,jhp))**4) )


      if (ettpa(npt,jmp)*ettpa(npt,jcp) .le. precis) then
        aux4 = (ettpa(npt,jtf)+tkelvi) + aux1*(aux2+aux3)/precis
      else
        aux4 = (ettpa(npt,jtf)+tkelvi) + aux1*(aux2+aux3)/(ettpa(npt,jmp)*ettpa(npt,jcp))
      endif
      aux5 = dtp/aux1
      aux6 = exp(-aux5)
      ter1 = ettpa(npt,jhp) * aux6
      ter2 = aux4 * ( 1.d0-((1.d0-aux6)/aux5) )
      ettp(npt,jhp) = tsvar(npt,jhp) + 0.5d0*ter1 + ter2
    endif
  enddo
endif

!===============================================================================
! 14 Mise a jour du diametre du coeur retrecissant
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    if ( ettpa(npt,jmch).ge.(1.d-3*ettpa(npt,jmp)) ) then
     tepa(npt,jrdck) = tepa(npt,jrd0p)
    else
     tepa(npt,jrdck) =                                            &
             ( (d6spi / ( 1.d0-xashch(icha)) )                    &
              *( ettp(npt,jmch)/rho0ch(icha)                      &
              +ettp(npt,jmck)/tepa(npt,jrhock)))**d1s3
    endif
  endif
enddo

!===============================================================================
! 15. Calcul du diametre des grains de charbon
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
! 16. Calcul de la masse des grains de charbon
!===============================================================================

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    icha = itepa(npt,jinch)
    mp0  = dpis6 * (tepa(npt,jrd0p)**3) * rho0ch(icha)
    ettp(npt,jmp) = ettp(npt,jmch) + ettp(npt,jmck) + ettp(npt,jmwat)          &
                      + xashch(icha)*mp0
  endif
enddo

!===============================================================================

! Free memory
deallocate(tempf)
deallocate(gamhet, deltah)

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

1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON        ',/,&
'@      EST ACTIVE (LAGICH) AVEC COUPLAGE RETOUR THERMIQUE    ',/,&
'@                                                            ',/,&
'@       IPHYLA = ', I10                                       ,/,&
'@       LTSTHE = ', I10                                       ,/,&
'@                                                            ',/,&
'@  Le transport Lagrangien de particule de charbon ne peut   ',/,&
'@   etre couple avec la phase Eulerienne depuis              ',/,&
'@   l''introduction des termes de sechage au moins           ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de LTSTHE dans la subroutine USLAG1 et ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
