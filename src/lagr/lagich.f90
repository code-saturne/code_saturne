!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                &
   itepa  , propce , ettp   , ettpa  ,                         &
   tepa   , tempct , tsvar  , cpgd1  , cpgd2  ,                &
   cpght                                         )

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
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! cpgd1,cpgd2,     ! tr ! --> ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpmax    !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
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
use field

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp , nvp1 , nvep , nivep

integer          itepa(nbpmax,nivep)

double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision cpgd1(nbpmax), cpgd2(nbpmax), cpght(nbpmax)

! Local variables

integer          npt , iel , icha , mode , iii
integer          ilayer , ilayer_het
double precision aux1 , aux2 , aux3 , aux4 , aux5
double precision volume_couche , rayon(nlayer) , mlayer(nlayer)
double precision mwater(nlayer) , mwat_max, fwat(nlayer), fcoke(nlayer)
double precision rep, prt, xnul, xrkl, sherw
double precision coef , mp0 , d6spi , dpis6 , d1s3 , d2s3, mv
double precision f1mc(ncharm) , f2mc(ncharm)
double precision coefe(ngazem)
double precision phith(nlayer), temp(nlayer)

double precision skp1(nlayer) , skp2(nlayer) , skglob, gamhet, deltah

double precision precis, lv, tebl, tlimit, tmini

double precision, dimension(:), pointer :: cromf
double precision, dimension(:), pointer :: viscl

precis = 1.d-15                   ! Petit nombre (pour la precision numerique)
lv = 2.263d+6                     ! Chaleur Latente en J/kg
tebl = 100.d0 + tkelvi            ! Temperature d'ebulition de l'eau
tlimit = 302.24d0                 ! Temperature limite

! Temperature mini (apres, la fraction massique d'eau saturante est nulle)
tmini = tlimit*(1.d0-tlimit*rr/(lv*wmole(ih2o)))

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Verification de la presence d'une physique
if ( ippmod(iccoal).lt.0 .and. ippmod(icpl3c).lt.0 ) then
  write(nfecra,1000) iphyla, ippmod(icpl3c), ippmod(iccoal)
  call csexit (1)
  !==========
endif

! Initialize variables to avoid compiler warnings
coef = 0.d0
d6spi = 6.d0 / pi
dpis6 = pi / 6.d0
d1s3  = 1.d0 / 3.d0
d2s3  = 2.d0 / 3.d0

! --- Si couplage retour thermique :
if ( ltsthe.eq.1 ) then
  coef = 1.d0 / dble(nordre)
  if (nor.eq.1 ) then
    do npt = 1,nbpart
      cpgd1(npt) = 0.d0
      cpgd2(npt) = 0.d0
      cpght(npt) = 0.d0
    enddo
  endif
endif

!===============================================================================
! 2. Pointeurs vers la masse volumique du fluide porteur
!===============================================================================

if (ippmod(iccoal).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

call field_get_val_s(iprpfl(iviscl), viscl)

!===============================================================================
! 3. Boucle principale sur l'ensemble des particules
!===============================================================================
do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then

    ! Variables generiques
    iel  = itepa(npt,jisor)
    icha = itepa(npt,jinch)
    volume_couche = dpis6 * (tepa(npt,jrd0p)**3) / float(nlayer)

    ! Calcul du Reynolds
    aux1 = sqrt( ( ettp(npt,juf) -ettp(npt,jup) )*                             &
                 ( ettp(npt,juf) -ettp(npt,jup) )                              &
               + ( ettp(npt,jvf) -ettp(npt,jvp) )*                             &
                 ( ettp(npt,jvf) -ettp(npt,jvp) )                              &
               + ( ettp(npt,jwf) -ettp(npt,jwp) )*                             &
                 ( ettp(npt,jwf) -ettp(npt,jwp) )  )
    xnul = viscl(iel) / cromf(iel)
    rep  = aux1 * ettp(npt,jdp) / xnul

    ! Calcul du Prandtl et du Sherwood
    if (ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2) then
      xrkl = diftl0 / cromf(iel)
    else if (ivisls(iscalt).ge.1) then
      xrkl = propce(iel,ipproc(ivisls(iscalt))) / cromf(iel)
    else
      xrkl = visls0(iscalt) / cromf(iel)
    endif
    prt   = xnul / xrkl
    sherw = 2 + 0.55d0 * rep**0.5d0 * prt**(d1s3)

    ! Calcul des rayons de discrétisation
    do ilayer = 1, nlayer
      rayon(ilayer)=(((tepa(npt,jrd0p)/2.0d0)**3)                              &
                      *(float(ilayer)/float(nlayer)))**(1.d0/3.d0)
    enddo
    mp0  = dpis6 * (tepa(npt,jrd0p)**3) * rho0ch(icha)
    mwat_max = xwatch(icha)*mp0/float(nlayer)

    ! Calcul de la quantité d'eau sur chaque couche
    aux1 = ettp(npt,jmwat)
    do ilayer = 1, nlayer
      if (ilayer.eq.nlayer) then
        mwater(ilayer)=max(0.0d0,aux1)
      else
        mwater(ilayer)=max(0.0d0,min(aux1,mwat_max))
      endif
      aux1=aux1-mwater(ilayer)
    enddo

    ! Masse sur chaque couche
    do ilayer = 1, nlayer
      mlayer(ilayer) = xashch(icha)*mp0/float(nlayer)                          &
                      +mwater(ilayer)                                          &
                      +ettp(npt,jmch(ilayer))                                  &
                      +ettp(npt,jmck(ilayer))
    enddo


!===============================================================================
! 4. Calcul de la masse d'eau qui s'evapore
!    On suppose pour le calcul de la masse volumique du charbon actif que
!    le sechage a lieu a volume constant
!===============================================================================

    ! --- Calcul du flux de vapeur pour la particule
    call lagsec                                                                &
    !==========
    ( nbpmax, nvp    , nvp1   , nvep     , nivep ,                             &
      npt   ,                                                                  &
      itepa , propce , ettp   , ettpa    , tepa  , tempct , tsvar ,            &
      rayon , mlayer , mwater , mwat_max , volume_couche  , sherw , fwat   )


!===============================================================================
! 5. Calcul des constantes de vitesses SPK1 et SPK2 du transfert
!    de masse par devolatilisation avec des lois d'Arrhenius
!===============================================================================

    !     RR --> Constante des gaz parfaits en J/mol/K
    do ilayer = 1, nlayer

      aux1 = 1.d0 / (rr*ettp(npt,jhp(ilayer)))
      skp1(ilayer) = a1ch(icha) * exp( -e1ch(icha) * aux1)
      skp2(ilayer) = a2ch(icha) * exp( -e2ch(icha) * aux1)

      aux1 = skp1(ilayer) * y1ch(icha) * ettp(npt,jmch(ilayer))
      aux2 = skp2(ilayer) * y2ch(icha) * ettp(npt,jmch(ilayer))

      ! --- Couplage retour thermique
      if ( ltsthe.eq.1 ) then
        cpgd1(npt) = cpgd1(npt) + coef*aux1
        cpgd2(npt) = cpgd2(npt) + coef*aux2
      endif

    enddo


!===============================================================================
! 6. Calcul de la constante globale de combustion hétérogène
!===============================================================================

    ! --- Repérage de la couche où se déroule la combustion hétérogène
    ! On repere la couche avec du ch la plus externe
    ilayer_het = 1
    do ilayer = 1 , nlayer
      if (ettpa(npt,jmch(ilayer)).gt.0.0d0 ) then
        ilayer_het = ilayer
      endif
    enddo

    ! On verifie cherche s'il reste du ck sur une couche plus externe
    do ilayer = ilayer_het , nlayer
      if (ettpa(npt,jmck(ilayer)).gt.0.0d0 ) then
        ilayer_het = ilayer
      endif
    enddo

    ! --- Coefficient de cinetique chimique de formation de CO
    !       en (kg.m-2.s-1.atm(-n))
    ! Conversion (kcal/mol -> J/mol)
    aux1 = ehetch(icha) * 1.0d3 * xcal2j
    aux2 = ahetch(icha)                                                      &
      * exp(- aux1 / (rr*ettp(npt,jhp(ilayer_het))) )

    ! --- Coefficient de diffusion en  (Kg/m2/s/atm) et constante
    !     globale de reaction
    if ( tepa(npt,jrdck).gt.precis ) then
      ! La constante 2.53d-7 est expliquée dans le tome 5 du rapport sur les
      ! physiques particulières de Code_Saturne (HI-81/04/003/A) équation 80
      aux3 = sherw * 2.53d-7 * (propce(iel,ipproc(itemp1))**0.75d0)           &
                             / tepa(npt,jrdck)
      skglob = (aux2*aux3) / (aux2+aux3)
    else
      skglob = aux2
    endif


!===============================================================================
! 7. Calcul de la GAMMAhet
!===============================================================================

    ! --- Calcul de la pression partielle en oxygene (atm)
    !                                                 ---
    !       PO2 = RHO1*RR*T*YO2/MO2
    aux1 = cromf(iel) * rr * propce(iel,ipproc(itemp1))                 &
         * propce(iel,ipproc(iym1(io2))) / wmole(io2) / prefth

    ! --- Calcul de surface efficace : SE
    aux2 =  pi * (1.0d0-xashch(icha)) * tepa(npt,jrdck)**2

    ! --- Pas de combustion heterogene si Mch/Mp >= 1.D-3
    if ( ettpa(npt,jmch(1)).ge.(1.d-3*mlayer(1)) ) then
      gamhet = 0.d0
    else
      ! --- Calcul de la GamHET
      gamhet = aux1 * aux2 * skglob
    endif

    ! --- Couplage retour thermique
    if ( ltsthe.eq.1 ) then
      cpght(npt) = cpght(npt) + coef * gamhet
    endif


!===============================================================================
! 8. Calcul de la 0.5(MO2/MC)*(HO2(Tp)-HO2(TF))
!===============================================================================
    ! --- Calcul de Hc(Tp)-Mco/Mc Hco2(Tp)+0.5Mo2/Mc Ho2(Tf)

    !        Calcul de Hcoke(TP)
    aux1 = h02ch(icha) + ettp(npt,jcp)*(ettp(npt,jhp(ilayer_het))-trefth)

    !        Calcul de MCO/MC HCO(TP)
    do iii = 1, ngazem
      coefe(iii) = zero
    enddo
    coefe(ico) = wmole(ico) / wmolat(iatc)

    do iii = 1, ncharm
      f1mc(iii) = zero
      f2mc(iii) = zero
    enddo
    mode      = -1
    call cpthp1 ( mode , aux2 , coefe  , f1mc , f2mc ,  ettp(npt,jhp(ilayer_het)) )
    !==========

    !        Calcul de MO2/MC/2. HO2(TF)
    do iii = 1, ngazem
      coefe(iii) = zero
    enddo
    coefe(io2) = wmole(io2) / wmolat(iatc) / 2.d0

    do iii = 1, ncharm
      f1mc(iii) = zero
      f2mc(iii) = zero
    enddo
    mode      = -1
    aux3      = ettp(npt,jtf) + tkelvi
    call cpthp1 ( mode  , aux4 , coefe , f1mc  , f2mc , aux3 )
    !==========

    deltah = aux2 - aux4 - aux1


!===============================================================================
! 9. Integration Masse d eau
!===============================================================================

    if (nor.eq.1) then
      aux1 = 0.0d0
      do ilayer = 1, nlayer
        aux1 = aux1 + fwat(ilayer)*dtp
      enddo
      ettp(npt,jmwat) = ettpa(npt,jmwat)-aux1

      ! Clipping
      if ( ettp(npt,jmwat).lt.precis ) then
        ettp(npt,jmwat) = 0.d0
      endif

    else if (nor.eq.2) then
      aux1 = 0.0d0
      do ilayer = 1, nlayer
        aux1 = aux1 + fwat(ilayer)*dtp
      enddo
      ettp(npt,jmwat) = 0.5d0 * ( ettp(npt,jmwat)+ettpa(npt,jmwat)-aux1 )

      ! Clipping
      if ( ettp(npt,jmwat).lt.precis ) then
        ettp(npt,jmwat) = 0.d0
      endif

    endif


!===============================================================================
! 10. Integration Masse de Charbon reactif
!===============================================================================

    if (nor.eq.1) then
      do ilayer = 1, nlayer
        aux1 = exp(-(skp1(ilayer)+skp2(ilayer))*dtp)
        ettp(npt,jmch(ilayer)) = ettpa(npt,jmch(ilayer))*aux1

        ! Clipping
        if ( ettp(npt,jmch(ilayer)).lt.precis ) then
          ettp(npt,jmch(ilayer)) = 0.d0
        endif
      enddo

    else if (nor.eq.2) then
      do ilayer = 1, nlayer
        aux1 = exp(-(skp1(ilayer)+skp2(ilayer))*dtp)
        ettp(npt,jmch(ilayer)) = 0.5d0 * ( ettp(npt,jmch(ilayer))              &
                                          +ettpa(npt,jmch(ilayer))*aux1 )

        ! Clipping
        if ( ettp(npt,jmch(ilayer)).lt.precis ) then
          ettp(npt,jmch(ilayer)) = 0.d0
        endif
      enddo

    endif


!===============================================================================
! 11. Integration Masse de Coke
!===============================================================================

    if (nor.eq.1) then

      ! On initialise le flux de comb hétérogène effectif
      do ilayer = 1, nlayer
        fcoke(ilayer) = 0.0d0
      enddo

      ! On boucle sur toutes les cellules qui ont du coke ou du charbon reactif
      do ilayer=ilayer_het,1,-1
        aux1 = ( skp1(ilayer) * (1.d0-y1ch(icha))                              &
                +skp2(ilayer) * (1.d0-y2ch(icha)) )                            &
               /( skp1(ilayer)+skp2(ilayer) )

        aux2 = exp(-(skp1(ilayer)+skp2(ilayer))*dtp)
        aux3 = aux1 * ettpa(npt,jmch(ilayer)) * (1.0d0-aux2) / dtp

        if ( ilayer.eq.ilayer_het ) then
          ! Calcul de la masse de coke équivalente
          aux4 = dpis6 * (1.0d0-xashch(icha)) * tepa(npt,jrdck)**3             &
                 * tepa(npt,jrhock(ilayer))

          if ( aux4.gt.precis ) then
            ! On tient compte de la combustion hétérogène
            aux5 = dtp * aux4 * (-gamhet+aux3)/(d2s3*gamhet*dtp+aux4)
            fcoke(ilayer) = gamhet
          else
            ! On néglige la comb hétérogène
            aux5 = dtp * aux3
          endif

        else
          ! On néglige la comb hétérogène
          aux5 = dtp * aux3
        endif

        ettp(npt,jmck(ilayer)) = ettpa(npt,jmck(ilayer)) + aux5
      enddo

      !  Si gamhet est trop important, on le repartit sur plusieurs couches
      do ilayer=ilayer_het,1,-1
        if ( ettp(npt,jmck(ilayer)).lt.0.d0 ) then

          ! On limite la comb hétérogène
          fcoke(ilayer) = fcoke(ilayer)+ettp(npt,jmck(ilayer))

          ! On attaque éventuellement la comb de la couche suivante
          if ( ilayer.gt.2 ) then
            ettp(npt,jmck(ilayer-1)) = ettp(npt,jmck(ilayer-1))+ettp(npt,jmck(ilayer))
          endif

          ! On limite la masse de coke
          ettp(npt,jmck(ilayer)) = 0.d0

        endif
      enddo

    else if (nor.eq.2) then
      ! Pas d'ordre 2 pour le moment
      call csexit(1)
    endif

!===============================================================================
! 12. Integration de la temperature des grains de charbon
!===============================================================================

    do ilayer = 1, nlayer
      ! Terme sources thermiques couche par couche (en W)
      ! Les échanges thermiques avec l'extérieur sont calculés directement par lagtmp
      phith(ilayer) = ( -fcoke(ilayer) * deltah )                              &
                        -fwat (ilayer) * lv
    enddo

    call lagtmp                                                                &
    !==========
    ( nbpmax , nvp    , nvp1  , nvep  , nivep ,                                &
      npt    ,                                                                 &
      itepa  , propce , ettp  , ettpa , tepa  , tempct ,                       &
      rayon  , mlayer , phith , temp  , tsvar , volume_couche )

    do ilayer = 1, nlayer
      ettp(npt,jhp(ilayer)) = temp(ilayer)
    enddo


!===============================================================================
! 13. Mise a jour du diametre de la masse volumique du coke
!===============================================================================

    do ilayer = 1, nlayer

      if (ettpa(npt,jmch(ilayer)).ge.1.d-3*mlayer(ilayer)) then
        ! mv represente la mlayer qui a quitté le grain ( sechage + pyrolyse)
        mv = mp0*(1-xashch(icha))/float(nlayer)                                &
                - ettp(npt,jmch(ilayer))                                       &
                - ettp(npt,jmck(ilayer))                                       &
                - (mwater(ilayer)-fwat(ilayer)*dtp)
        ! masse volumique du coke SEUL
        tepa(npt,jrhock(ilayer)) = rho0ch(icha)                                &
                                  -mv/(volume_couche*(1.d0-xashch(icha)) )
      endif

    enddo


!===============================================================================
! 14. Mise a jour du diametre du coeur retrecissant
!===============================================================================

    ! On repere la couche avec du ch la plus externe
    ilayer_het = 1
    do ilayer = 1 , nlayer
      if (ettpa(npt,jmch(ilayer)).gt.0.0d0 ) then
        ilayer_het = ilayer
      endif
    enddo

    ! On verifie cherche s'il reste du ck sur une couche plus externe
    do ilayer = ilayer_het , nlayer
      if (ettpa(npt,jmck(ilayer)).gt.0.0d0 ) then
        ilayer_het = ilayer
      endif
    enddo

    if (ettp(npt,jmch(ilayer_het)).ge.1.d-3*mlayer(ilayer_het)) then
      ! La pyrolyse n'est pas terminée, le char a le diametre initial
      tepa(npt,jrdck) = 2.0d0*rayon(ilayer_het)
    else
      ! On repartit le char de façon uniforme
      if (ilayer_het.eq.1) then
        tepa(npt,jrdck) =  ( (d6spi/(1.d0-xashch(icha)))                       &
                            *(ettp(npt,jmch(ilayer_het))                       &
                               /rho0ch(icha)                                   &
                             +ettp(npt,jmck(ilayer_het))                       &
                               /tepa(npt,jrhock(ilayer_het)))                  &
                                                               ) **d1s3
        ! Clipping
        if (tepa(npt,jrdck).gt.2.0d0*rayon(ilayer_het)) then
          tepa(npt,jrdck) = 2.0d0*rayon(ilayer_het)
        else if (tepa(npt,jrdck).lt.0.0d0) then
          tepa(npt,jrdck) = 0.0d0
        endif
      else
        tepa(npt,jrdck) = ( (2.0d0*rayon(ilayer_het-1))**3 +                   &
                           ( (d6spi/(1.d0-xashch(icha)))                       &
                            *(ettp(npt,jmch(ilayer_het))                       &
                               /rho0ch(icha)                                   &
                             +ettp(npt,jmck(ilayer_het))                       &
                               /tepa(npt,jrhock(ilayer_het)))                  &
                                                               ))**d1s3
        ! Clipping
        if (tepa(npt,jrdck).gt.2.0d0*rayon(ilayer_het)) then
          tepa(npt,jrdck) = 2.0d0*rayon(ilayer_het)
        else if (tepa(npt,jrdck).lt.2.0d0*rayon(ilayer_het-1)) then
          tepa(npt,jrdck) = 2.0d0*rayon(ilayer_het-1)
        endif
      endif
    endif


!===============================================================================
! 15. Calcul du diametre des grains de charbon
!===============================================================================

    ettp(npt,jdp) = (xashch(icha)*(tepa(npt,jrd0p)**2)            &
                         + (1.d0-xashch(icha))                    &
                          *(tepa(npt,jrdck)**2) )**0.5d0

!===============================================================================
! 16. Calcul de la masse des grains de charbon
!===============================================================================

    aux1 = 0.0d0
    do ilayer = 1, nlayer
      aux1 = aux1 + ettp(npt,jmch(ilayer)) + ettp(npt,jmck(ilayer))
    enddo

    ettp(npt,jmp) = aux1 + ettp(npt,jmwat) + xashch(icha)*mp0


!===============================================================================
! 17. Fin de la boucle principale sur l'ensemble des particules
!===============================================================================
  endif   ! itepa(npt,jisor).gt.0
enddo     !npt = 1,nbpart

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

end subroutine
