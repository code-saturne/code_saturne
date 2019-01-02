!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine cs_coal_readata

!===============================================================================
!  FONCTION  :
!  ---------

! LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
!      RELATIF A LA COMBUSTION CHARBON PULVERISE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use ihmpre

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=150) :: chain1,chain2

integer          it     , ice    , iat    , ii ,jj
integer          ncoel  , inicoe , inicha , ierror
integer          icla   , icha   , is
integer          idebch , ifinch , lonch  , ichai  , ichcoe
integer          atcoel(ngazem,natom)
integer          ierr   , ndim

double precision fcor   , pcisec , pcipur , xashpc
double precision hco20  , ho20   , hh2o0 , hso20
double precision den1   , den2
double precision tmin   , tmax
double precision wmolce(ngazem), ehcoel(ngazem,npot)
double precision cpcoel(ngazem,npot)
double precision dhvol1,dhvol2,hcoke,hchar,ehvol1,ehvol2
double precision wmco,wmo2,wmco2,wmh2o,wmn2,wmc
double precision dmf3,dmf4,dmf5,som1,som2
double precision sm(8),solu(8),mat(8,8)

! PCI-PCS
!
double precision pcibrut,pcssec,pcsbrut,pcspur,xwatpc
!
!MODEL DE NOx
!============
!
!Terminologie
!FM: Fraction massique
!MV: Matiere volatile

!CHARBON et CHAR
!Nombre de mol de C du charbon et du char
double precision molcch,molcck
!Nombre de mol de H, O, N et S du charbon
double precision nalpha(ncharm), nbeta(ncharm), nteta(ncharm),nomega(ncharm)
!Nombre de mol de H, O, N et S du char1 (devolatilisation des MVs legeres)
double precision ngama1(ncharm),ndelt1(ncharm),nkapp1(ncharm),nzeta1(ncharm)
!Nombre de mol de H, O, N et S du char2 (devolatilisation des MVs lourdes)
double precision ngama2(ncharm),ndelt2(ncharm),nkapp2(ncharm),nzeta2(ncharm)
!
!COMPOSITION des MVs
!Rapport H/C des MVs legeres, FM des MVs legeres, Rapport H/C des MVs lourdes,
!FM des MVs lourdes
double precision nchx1(ncharm),ny1ch(ncharm),nchx2(ncharm),ny2ch(ncharm)
!Nombre de mol de CHx1, CO, H2O, H2S, HCN, NH3 des MV legeres
double precision noxa1(ncharm),noxb1(ncharm),noxc1(ncharm),noxd1(ncharm),      &
                 noxe1(ncharm),noxf1(ncharm)
!Nombre de mol de CHx1, CO, H2O, H2S, HCN, NH3 des MV lourdes
double precision noxa2(ncharm),noxb2(ncharm),noxc2(ncharm),noxd2(ncharm),      &
                 noxe2(ncharm),noxf2(ncharm)
!Masse molaire de CHx1 et CHx2
double precision ychx1t,ychx2t
!COMPOSITION des PRODUITS de la REACTION HETEROGENE
!Nombre de mol de O2, CO, HCN et CO (combustion du char1)
double precision noxh1(ncharm),noxi1(ncharm),noxj1(ncharm),noxk1(ncharm)
!Nombre de mol de O2, CO, HCN et CO (combustion du char2)
double precision noxh2(ncharm),noxi2(ncharm),noxj2(ncharm),noxk2(ncharm)
!
!Modele de REBURNING
! Constantes cinetiques tabulees en fonction de la temperature
double precision kf1(7),kf2(7),kf3(7),kf4(7),kr4(7),kf5(7),kr5(7),kf6(7),      &
                 kr6(7),kf7(7),kr7(7)
! Contante selon la doc theorique de FLUENT
double precision pflue
!===============================================================================

!===============================================================================
! 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
!===============================================================================

ncoel = 13

if ( ncoel.gt.ngazem ) then
  write(nfecra,9991) ngazem,ncoel
  call csexit (1)
endif

! ---- Nb de points de tabulation ENTH-TEMP

npo = 8
if ( npo.gt.npot ) then
  write(nfecra,9992) npot,npo
  call csexit (1)
endif

! --- Lecture des noms des constituants elementaires

do ice=1,ncoel
  do inicoe=1,len(nomcoe(ice))
    NOMCOE(ICE)(INICOE:INICOE)=' '
  enddo
enddo

do inicha=1,len(chain1)
  CHAIN1(INICHA:INICHA)=' '
enddo

do inicha=1,len(chain2)
  CHAIN2(INICHA:INICHA)=' '
enddo

chain1 = 'CH4 C2H4 CO H2S H2 HCN NH3 O2 CO2 H2O SO2 N2 C(S)'
call verlon (chain1, idebch, ifinch, lonch)
chain2(1:lonch)=chain1(idebch:ifinch)

ice=1
ichcoe=0
do ichai=1,lonch
  IF (CHAIN2(ICHAI:ICHAI).NE.' ') THEN
    ichcoe=ichcoe+1
    nomcoe(ice)(ichcoe:ichcoe) =chain2(ichai:ichai)
  else
    if (ichcoe.ne.0) then
      ice=ice+1
      ichcoe=0
    endif
  endif
enddo

! --- Temperature Min et Max

tmin = 300.d0
tmax = 2400.d0

! ---- Nb especes atomiques (C, H, O, S et N)

nato = 5
if ( nato.gt.natom ) then
  write(nfecra,9993) natom,nato
  call csexit (1)
  !==========
endif

! ---- Masse molaire especes atomiques
!      Composition des constituants elementaires en fonction
!        des especes elementaires

wmolat(1) = 0.012d0
wmolat(2) = 0.001d0
wmolat(3) = 0.016d0
wmolat(4) = 0.014d0
wmolat(5) = 0.032d0

! CH4
atcoel(1,1) = 1
atcoel(1,2) = 4
atcoel(1,3) = 0
atcoel(1,4) = 0
atcoel(1,5) = 0

! C2H4
atcoel(2,1) = 2
atcoel(2,2) = 4
atcoel(2,3) = 0
atcoel(2,4) = 0
atcoel(2,5) = 0

! CO
atcoel(3,1) = 1
atcoel(3,2) = 0
atcoel(3,3) = 1
atcoel(3,4) = 0
atcoel(3,5) = 0

! H2S
atcoel(4,1) = 0
atcoel(4,2) = 2
atcoel(4,3) = 0
atcoel(4,4) = 0
atcoel(4,5) = 1

! H2
atcoel(5,1) = 0
atcoel(5,2) = 2
atcoel(5,3) = 0
atcoel(5,4) = 0
atcoel(5,5) = 0

! 'HCN
atcoel(6,1) = 1
atcoel(6,2) = 1
atcoel(6,3) = 0
atcoel(6,4) = 1
atcoel(6,5) = 0

! NH3
atcoel(7,1) = 0
atcoel(7,2) = 3
atcoel(7,3) = 0
atcoel(7,4) = 1
atcoel(7,5) = 0

! O2
atcoel(8,1) = 0
atcoel(8,2) = 0
atcoel(8,3) = 2
atcoel(8,4) = 0
atcoel(8,5) = 0

! CO2
atcoel(9,1) = 1
atcoel(9,2) = 0
atcoel(9,3) = 2
atcoel(9,4) = 0
atcoel(9,5) = 0

! H2O
atcoel(10,1) = 0
atcoel(10,2) = 2
atcoel(10,3) = 1
atcoel(10,4) = 0
atcoel(10,5) = 0

! SO2
atcoel(11,1) = 0
atcoel(11,2) = 0
atcoel(11,3) = 2
atcoel(11,4) = 0
atcoel(11,5) = 1

! N2
atcoel(12,1) = 0
atcoel(12,2) = 0
atcoel(12,3) = 0
atcoel(12,4) = 2
atcoel(12,5) = 0

! C(S)
atcoel(13,1) = 1
atcoel(13,2) = 0
atcoel(13,3) = 0
atcoel(13,4) = 0
atcoel(13,5) = 0

! ---- Calcul des masses molaires des constituants elementaires

do ice = 1, ncoel
  wmolce(ice) = 0.d0
  do iat = 1, nato
    wmolce(ice)= wmolce(ice) + atcoel(ice,iat)*wmolat(iat)
  enddo
enddo

!===============================================================================
! 2. CALCULS DE DONNEES COMPLEMENTAIRES
!===============================================================================

! --> Discretisation de la temperature

do it = 1, npo
  th(it) = dble(it-1)*(tmax-tmin)/dble(npo-1) + tmin
enddo

! --> Calcul des enthalpies pour les differentes especes courantes

call pptbht                                                       &
!==========
 ( ncoel  ,                                                       &
   nomcoe , ehcoel , cpcoel , wmolce )

! --> Calcul tabulation enthalpie - temperature pour le melange gazeux

! ---- Nb de constituants gazeux
!     ATTENTION ON COMPTE EGALEMENT CH4 et le monomere CH2
!  on rajoute : H2S, H2 et SO2 (on passe de 5 a 8)
!

ngaze = 2 + 10 + 2*ncharb
ngazg = 2 + 10

! ---- Definition des pointeurs pour les tableaux WMOLE et EHGAZE
!      REMARQUE : Cette position de pointeurs va egalement servir
!                 pour le tableau de pointeurs IYM1
!                 (indices de propriétés)
!                 ON BALAYE JUSTE DE 1 A (NGAZE-2*NCHARB=ngazg)

is    = 0
is    = is + 1
ichx1 = is
is    = is + 1
ichx2 = is
is    = is + 1
ico   = is
is    = is + 1
ih2s  = is
is    = is + 1
ihy   = is
is    = is + 1
ihcn  = is
is    = is + 1
inh3  = is
is    = is + 1
io2   = is
is    = is + 1
ico2  = is
is    = is + 1
ih2o  = is
is    = is + 1
iso2  = is
is    = is + 1
in2   = is
do icha = 1, ncharb
  ichx1c(icha) = is + icha
  ichx2c(icha) = is + icha + ncharb
enddo

! ---- Remplissage de EHGAZE et WMOLE
!         a partir de EHCOEL et WMOLCE

! ------ ATTENTION :
!        On prend par defaut
!          du CH4 pour CHX1m
!          et le monomere CH2 (jouant le role du C2H4) pour CHX2m

do it = 1, npo
  ehgaze(ichx1,it) = ehcoel( 1,it)
  ehgaze(ichx2,it) = ehcoel( 2,it)
  ehgaze(ico  ,it) = ehcoel( 3,it)
  ehgaze(ih2s ,it) = ehcoel( 4,it)
  ehgaze(ihy  ,it) = ehcoel( 5,it)
  ehgaze(ihcn ,it) = ehcoel( 6,it)
  ehgaze(inh3 ,it) = ehcoel( 7,it)
  ehgaze(io2  ,it) = ehcoel( 8,it)
  ehgaze(ico2 ,it) = ehcoel( 9,it)
  ehgaze(ih2o ,it) = ehcoel(10,it)
  ehgaze(iso2 ,it) = ehcoel(11,it)
  ehgaze(in2  ,it) = ehcoel(12,it)
enddo
wmole(ichx1) = wmolce( 1)
! ------ ATTENTION : Il faut prendre la masse molaire du monomere CH2
!                    et non celle du C2H4
wmole(ichx2) = wmolat(iatc)+wmolat(iath)*2.d0
wmole(ico  ) = wmolce( 3)
wmole(ih2s ) = wmolce( 4)
wmole(ihy  ) = wmolce( 5)
wmole(ihcn ) = wmolce( 6)
wmole(inh3 ) = wmolce( 7)
wmole(io2  ) = wmolce( 8)
wmole(ico2 ) = wmolce( 9)
wmole(ih2o ) = wmolce(10)
wmole(iso2 ) = wmolce(11)
wmole(in2  ) = wmolce(12)


! --> Calcul tabulation enthalpie - temperature pour le solide
!     Charbon reactif, Coke et Cendres

! ---- Nb de constituants solide

nsolid = 4*ncharb

! ---- Definition des pointeurs ICH, ICK et IASH

is = 0
do icha = 1, ncharb
  is         = is + 1
  ich(icha)  = is
  is         = is + 1
  ick(icha)  = is
  is         = is + 1
  iash(icha) = is
  is         = is + 1
  iwat(icha) = is
enddo

! ---- Calcul relatif au charbon reactif

do icha = 1, ncharb

! ------ Composition du charbon reactif
!        sous la forme CH(ALPHA)O(BETA)N(TETA)S(OMEGA)
!        et calcul de la masse molaire

  alpha(icha) = (hch(icha)/wmolat(iath))           &
              / (cch(icha)/wmolat(iatc))
  beta(icha)  = (och(icha)/wmolat(iato))           &
              / (cch(icha)/wmolat(iatc))
!==============================================================================
! Actuellement, les especes azotees ne sont pas considerees au cours de la
! combustion de la phase gaz. Ainsi, il faut assurer que n'aucune espece
! azotee est liberee dans la phase gaz.
!==============================================================================
  teta(icha)  = ((nch(icha)-nch(icha))/wmolat(iatn))                          &
              / (cch(icha)/wmolat(iatc))

  omega(icha) = (sch(icha)/wmolat(iats))                                      &
              / (cch(icha)/wmolat(iatc))
!==============================================================================
!Correction des correlations teta - Soufre et omega - Azote.
  wmols(ich(icha)) = wmolat(iatc)                                             &
                   + alpha(icha)*wmolat(iath)                                 &
                   + beta (icha)*wmolat(iato)                                 &
                   + teta (icha)*wmolat(iatn)                                 &
                   + omega(icha)*wmolat(iats)
!==============================================================================
enddo

!
! Transformation des coef de repartition de l'azote en HCN etNH3
!
if ( ieqnox .eq. 1 ) then
  do icha = 1, ncharb
    if ( nch(icha) .gt. 0.D0 ) then
      som1 = crepn1(1,icha)+crepn1(2,icha)
      som2 = crepn2(1,icha)+crepn2(2,icha)
      if ( som1 .lt. 0.D0 .or. som2 .lt. 0.D0 ) then
        write(nfecra,9971) ICHA
        call csexit(1)
      endif
      crepn1(1,icha)= crepn1(1,icha)/som1
      crepn1(2,icha)= crepn1(2,icha)/som1
      crepn2(1,icha)= crepn2(1,icha)/som2
      crepn2(2,icha)= crepn2(2,icha)/som2
    endif
  enddo
endif

! ------ Transformation par la formule de Schaff du
!        PCI sur sec en PCI sur pur LORSQUE IPCI = 1
!
! ------ Calcul du PCI(pur)

!    si IPCI = 1 : PCI(sec )---> PCI(pur)
!    si IPCI = 2 : PCI(brut)---> PCI(pur)
!    si IPCI = 3 : PCS(pur )---> PCI(pur)
!    si IPCI = 4 : PCS(sec )---> PCI(pur)
!    si IPCI = 5 : PCS(brut)---> PCI(pur)
!    si IPCI = 6 : correlation IGT

fcor = 1.2d0
do icha = 1, ncharb

  xwatpc = xwatch(icha)*100.d0

  if (ipci(icha).eq.1) then
    ! On transforme le PCI sur sec J/kg ---> kcal/kg
    pcisec = pcich(icha)/1000.d0/xcal2j
    ! Calcul du PCI sur pur
    ! DS pas de passage de xashsec au lieu de xashpc
    xashpc = xashch(icha)*100.d0
    pcipur = ( pcisec+5.95d0*(zero+(1.d0-fcor)*xashsec(icha)) )*100.d0   &
           / ( 100.d0-zero-fcor*xashsec(icha) )
    !Conversion de PCI sur pur de kcal/kg ----> J/kg
    pcipur = pcipur*1000.d0*xcal2j
    !On ecrase PCICH sur sec par PCIPUR
    pcich(icha) = pcipur

  else if ( ipci(icha).eq.2 ) then

    !On transforme le PCI sur brut J/kg ---> kcal/kg
    pcibrut = pcich(icha)/1000.d0/xcal2j
    !Calcul du PCI sur pur
    pcipur = ( pcibrut + 5.95d0*(xwatpc+(1.d0-fcor)       &
                               *xashsec(icha)) )*100.d0   &
             / ( 100.d0-xwatpc-fcor*xashsec(icha) )
    !Conversion de PCI sur pur de kcal/kg ----> J/kg
    pcipur = pcipur*1000.d0*xcal2j
    !On ecrase PCICH sur sec par PCIPUR
    pcich(icha) = pcipur

  else if ( ipci(icha).ge.3 ) then

    if ( ipci(icha).eq.3 ) then

      ! On transforme le PCS(pur) de J/kg ---> KJ/kg

      pcspur = pcich(icha)/1000.d0

      ! Calcul du PCS(sec)

      pcssec= pcspur*(100.D0-xashsec(icha))/100.D0

    else if ( ipci(icha).eq.4 ) then

      ! On transforme le PCS(psec) de J/kg ---> KJ/kg

      pcisec = pcich(icha)/1000.d0

    else if ( ipci(icha).eq.5 ) then

      ! On transforme le PCS(brut) de J/kg ---> KJ/kg

      pcsbrut = pcich(icha)/1000.d0

     ! Calcul du PCS(sec) en KJ/kg

      pcssec = pcsbrut*100.D0/(100.D0-xwatpc)

    else if ( ipci(icha).eq.6 ) then

      ! Calcul du PCS(sec) en KJ/kg

      pcssec = 340.94d0*cch(icha) + 1322.98d0*hch(icha)            &
              + 68.3844d0*sch(icha)-119.86d0*(och(icha)+nch(icha)) &
              - 15.305d0 *xashsec(icha)
    endif
!
!   Calcul du PCI(sec) a partir du PCS(sec)
!
    pcisec = pcssec -226.d0*hch(icha)
!
!   On transforme le PCI sur sec KJ/kg ---> kcal/kg
!
    pcisec = pcisec/xcal2j
!
!   Calcul du PCI sur pur
!
    pcipur = ( pcisec+5.95d0*(zero+(1.d0-fcor)*xashsec(icha)) )*100.d0   &
            /( 100.d0-zero-fcor*xashsec(icha) )
!
!   Conversion de PCI sur pur de kcal/kg ----> J/kg
!
    pcipur = pcipur*1000.d0*xcal2j
!
!   on ecrase PCICH par PCIPUR
!
    pcich(icha) = pcipur
!
  endif

enddo


! ------ Calcul de : H02 pour CH , CK , ASH et WT

do icha = 1, ncharb
  hco20 = ehgaze(ico2,1)*wmole(ico2)
  hh2o0 = ehgaze(ih2o,1)*wmole(ih2o)
  ho20  = ehgaze(io2 ,1)*wmole(io2 )
  hso20 = ehgaze(iso2,1)*wmole(iso2)
!=============================================================================
! Correction de la correlation teta - soufre.
!
  h02ch(icha) = pcich(icha) +                                                &
     (  hco20 + alpha(icha)/2.d0*hh2o0                                       &
              + omega(icha)     *hso20                                       &
      -( 1.d0 + alpha(icha)/4.d0                                             &
              - beta (icha)/2.d0                                             &
              + omega(icha)      )*ho20  ) / wmols(ich(icha))
!=============================================================================
  eh0sol(ich(icha)) = h02ch(icha)

enddo

!   Pour l'instant H02CK et H02ASH sont nulles

do icha = 1, ncharb
  eh0sol(ick (icha)) = 0.d0
  eh0sol(iash(icha)) = 0.d0
enddo

!   Calcul de H02 pour l'eau liquide : d'apres JANAF
!                   H0 = -15871802.2
!                   CP = 4181.35

do icha = 1, ncharb
  eh0sol(iwat(icha)) = -15871802.2d0
enddo

! Construction table enthalpie/temperature Charbon
!   ---> Pour l'instant meme table que le Gaz

!  Nombre de tabulation : NPOC

npoc = npo

do ii = 1, npoc
  thc(ii) = th(ii)
enddo


! ------ Calcul de EHSOLI pour le charbon reactif
!        Si CP2CH > 0 : HCH = H02CH + CP2CH(T2-TREFTH)
!        Sinon       : On considere le PCI constant qqs T2

do icha = 1, ncharb
  if ( cp2ch(icha).gt.epsicp ) then
    do it = 1, npoc
      ehsoli(ich(icha),it) = eh0sol(ich(icha))                  &
                            + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      hco20 = ehgaze(ico2,it)*wmole(ico2)
      hh2o0 = ehgaze(ih2o,it)*wmole(ih2o)
      ho20  = ehgaze(io2 ,it)*wmole(io2 )
      hso20 = ehgaze(iso2,it)*wmole(iso2)
!==============================================================================
! Correction de la correlation teta - soufre.
!
      ehsoli(ich(icha),it) = pcich(icha)                                      &
           +( hco20 + alpha(icha)/2.d0*hh2o0                                      &
           + omega(icha)     *hso20                                      &
           -( 1.d0 + alpha(icha)/4.d0                                            &
           - beta (icha)/2.d0                                            &
           + omega(icha)      ) *ho20 )  / wmols(ich(icha))
!==============================================================================
    enddo
  endif
  gamma(icha) = zero
  delta(icha) = zero
  kappa(icha) = zero
  zeta (icha) = zero
enddo

! ---- Calcul relatif au coke

! ------ Par defaut
!        Coke = Carbone solide
!        GAMMA = 0 , DELTA = 0 , KAPPA = 0, ZETA = 0
!        Si CP2CH > 0 : HCK = H02CH + CP2CH(T2-TREFTH)
!        Sinon        : HCK = Enthalpie du carbone pur

do icha = 1, ncharb
  wmols(ick(icha)) = wmolce(ncoel)
  if (cp2ch(icha).gt.epsicp) then
    do it = 1, npoc
      ehsoli(ick(icha),it) = eh0sol(ick(icha))                  &
                           + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      ehsoli(ick(icha),it) = ehcoel(ncoel,it)
    enddo
  endif
enddo

! ------ Coke = CH(GAMMA)O(DELTA)N(KAPPA)S(ZETA)
do icha = 1, ncharb
  if ( pcick(icha).gt.epsicp ) then
    gamma(icha) = hck(icha)/cck(icha)
    delta(icha) = ock(icha)/cck(icha)
!==============================================================================
! Correction de la correlation kappa - soufre.
     kappa(icha) = nck(icha)/cck(icha)

! Correction de la correlation zeta - azote.
     zeta (icha) = sck(icha)/cck(icha)

! Correction des correlations kappa - soufre et zeta - azote.
    wmols(ick(icha)) = wmolat(iatc)                                           &
                      + gamma(icha)*wmolat(iath)                              &
                      + delta(icha)*wmolat(iato)                              &
                      + kappa(icha)*wmolat(iatn)                              &
                      + zeta (icha)*wmolat(iats)
!==============================================================================
!
!        On considere le PCI constant qqs T
    do it = 1, npoc
      hco20 = ehgaze(ico2,it)*wmole(ico2)
      hh2o0 = ehgaze(ih2o,it)*wmole(ih2o)
      hso20 = ehgaze(iso2,it)*wmole(iso2)
      ho20  = ehgaze(io2 ,it)*wmole(io2 )
!==============================================================================
! Correction de la correlation kappa - souffre.
      ehsoli(ick(icha),it) = pcick(icha)                                      &
       + ( hco20  + gamma(icha)/2.d0*hh2o0                                    &
                  + zeta(icha)     *hso20                                     &
         - ( 1.d0 + gamma(icha)/4.d0                                          &
                    - delta(icha)/2.d0                                        &
                    + zeta(icha)      )*ho20  )  / wmols(ick(icha))
!==============================================================================
    enddo
  endif
enddo

! ---- Calcul relatif aux cendres

do icha = 1, ncharb
  if (cp2ch(icha).gt.epsicp) then
    do it = 1, npoc
      ehsoli(iash(icha),it) = eh0sol(iash(icha))                &
                            + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      ehsoli(iash(icha),it) = h0ashc(icha)                      &
                            + cpashc(icha)*(thc(it)-trefth)
    enddo
  endif
  wmols(iash(icha)) = zero
enddo

! ---- Calcul relatif a l'eau

do icha = 1, ncharb
  cp2wat(icha) = 4181.35d0

  do it = 1, npoc
    ehsoli(iwat(icha),it) = eh0sol(iwat(icha))                  &
                           + cp2wat(icha)*(thc(it)-trefth)
  enddo
enddo

! --> Calcul relatifs aux matieres volatiles

!  On test que l'on utilise la meme option pour Y1 et Y2
!  pour chaque charbon

do icha = 1, ncharb
  if ( iy1ch(icha).ne. iy2ch(icha) ) then

    write(nfecra,9980) icha,iy1ch(icha),iy2ch(icha)
    call csexit(1)

  endif
  if ( iy1ch(icha).lt. 0 .or.  iy1ch(icha).gt.2 .or.            &
       iy2ch(icha).lt. 0 .or.  iy2ch(icha).gt.2     ) then
    write(nfecra,9981) icha,iy1ch(icha),iy1ch(icha)
    call csexit(1)
  endif
enddo

! Matieres volatiles legeres : [CH(CHX1); CO]
!     CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A1 CH(CHX1) + B1 CO + D1 H2S
!                              + E1 HCN + F1 NH3
!                              + (1-A1-B1-E1) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!       Si IY1CH = 0 CHX1 fixe , Y1CH calcule
!       Si IY1CH = 1 Y1CH fixe , CHX1 calcule
!       Si IY1CH = 2 Y1CH fixe , CHX1 fixe, on ajoute de l'eau
!     CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A1 CH(CHX1) + B1 CO + C1 H2O
!                              + D1 H2S + E1 HCN + F1 NH3
!                              + (1-A1-B1-E1) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!
do icha = 1, ncharb
!
!matrice 8x8
!
! Ligne 1
  mat(1,1) = -gamma(icha)
  mat(1,2) = -gamma(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.d0
  mat(1,5) = 1.d0-gamma(icha)
  mat(1,6) = 3.d0
  mat(1,7) = 1.d0
  mat(1,8) = 0.d0
! Ligne 2
  mat(2,1) = -delta(icha)
  mat(2,2) = 1.d0-delta(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
  mat(2,5) = -delta(icha)
  mat(2,6) = 0.d0
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0

  mat(3,1) = -kappa(icha)
  mat(3,2) = -kappa(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-kappa(icha)
  mat(3,6) = 1.d0
  mat(3,7) = 0.d0
  mat(3,8) = 0.d0

  mat(4,1) = 0.d0
  mat(4,2) = 0.d0
  mat(4,3) = 0.d0
  mat(4,4) = 0.d0
  mat(4,5) = crepn1(2,icha)
  mat(4,6) =-crepn1(1,icha)
  mat(4,7) = 0.d0
  mat(4,8) = 0.d0

  mat(5,1) = -zeta(icha)
  mat(5,2) = -zeta(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
  mat(5,5) = -zeta(icha)
  mat(5,6) = 0.d0
  mat(5,7) = 0.d0
  mat(5,8) = 0.d0

  mat(6,1) = wmolat(iatc)
  mat(6,2) = wmole(ico)
  mat(6,3) = wmole(ih2o)
  mat(6,4) = wmole(ih2s)
  mat(6,5) = wmole(ihcn)
  mat(6,6) = wmole(inh3)
  mat(6,7) = wmolat(iath)

!==============================================================================
!Correction des correlations teta - soufre et omega - azote.
  mat(6,8) = -(wmolat(iatc)+alpha(icha)*wmolat(iath)                          &
                           +beta (icha)*wmolat(iato)                          &
                           +teta (icha)*wmolat(iatn)                          &
                           +omega(icha)*wmolat(iats))
!==============================================================================
! Second membre
!
  sm(1) = alpha(icha)-gamma(icha)
  sm(2) = beta (icha)-delta(icha)
  sm(3) = teta (icha)-kappa(icha)
  sm(4) = 0.d0
  sm(5) = omega(icha)-zeta(icha)
  sm(6) = 0.d0

! On complete la matrice et le second menbre en fonction de
! la valeur de IY1CH

  if ( iy1ch(icha).eq.0 ) then
    ! Valeur de X1
    chx1(icha) = 4.D0
    ! Matrice
    mat(7,1) = -chx1(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = 0.d0
    sm(8) = 0.d0

  else if ( iy1ch(icha).eq.1 ) then

    ! Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0

    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = y1ch(icha)
    sm(8) = 0.D0

  else if ( iy1ch(icha).eq.2 ) then
    ! Valeur de X1
    chx1(icha) = 4.D0

    ! Matrice
    mat(7,1) = -chx1(icha)
    mat(7,2) = 0.D0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
    ! Second membre
    sm(7) = 0.D0
    sm(8) = y1ch(icha)

  endif

  ndim = 8
  call coal_resol_matrice( ndim, mat , sm , solu , ierr)

  if ( ierr .eq. 0 ) then
    a1(icha) = solu(1)
    b1(icha) = solu(2)
    c1(icha) = solu(3)
    d1(icha) = solu(4)
    e1(icha) = solu(5)
    f1(icha) = solu(6)
    if ( iy1ch(icha).eq.0 ) then
      y1ch(icha) = solu(8)
    else if ( iy1ch(icha).eq.1 ) then
      chx1(icha) = solu(7)/solu(1)
    endif
  else
    write(nfecra,9982) icha
    call csexit(1)
  endif
enddo

! Matieres volatiles lourdes : [CH(CHX2); CO]
!     CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A2 CH(CHX1) + B2 CO + D2 H2S
!                              + E2 HCN + F2 NH3
!                              + (1-A2-B2-E2) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!       Si IY2CH = 0 CHX2 fixe, Y2CH calcule
!       Si IY2CH = 1 Y2CH fixe, CHX2 calcule
!       Si IY2CH = 2 Y2CH fixe, CHX2 fixe, on ajoute de l'eau
!     CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A2 CH(CHX1) + B2 CO + C2 H2O
!                              + D2 H2S + E2 HCN + F2 NH3
!                              + (1-A2-B2-E2) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!
do icha = 1, ncharb
  mat(1,1) = -gamma(icha)
  mat(1,2) = -gamma(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.d0
!==============================================================================
  mat(1,5) = 1.d0-gamma(icha)
!==============================================================================
  mat(1,6) = 3.d0
!==============================================================================
  mat(1,7) = 1.D0
  mat(1,8) = 0.D0
!Ligne 2
  mat(2,1) = -delta(icha)
  mat(2,2) = 1.d0-delta(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
!==============================================================================
  mat(2,5) = -delta(icha)
!==============================================================================
  mat(2,6) = 0.d0
!==============================================================================
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0

  mat(3,1) = -kappa(icha)
  mat(3,2) = -kappa(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-kappa(icha)
  mat(3,6) = 1.d0
  mat(3,7) = 0.d0
  mat(3,8) = 0.d0

  mat(4,1) = 0.d0
  mat(4,2) = 0.d0
  mat(4,3) = 0.d0
  mat(4,4) = 0.d0
  mat(4,5) = crepn2(2,icha)
  mat(4,6) =-crepn2(1,icha)
  mat(4,7) = 0.d0
  mat(4,8) = 0.d0

  mat(5,1) = -zeta(icha)
  mat(5,2) = -zeta(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
!==============================================================================
  mat(5,5) = -zeta(icha)
!==============================================================================
  mat(5,6) = 0.d0
  mat(5,7) = 0.d0
  mat(5,8) = 0.d0
!
  mat(6,1) = wmolat(iatc)
  mat(6,2) = wmole(ico)
  mat(6,3) = wmole(ih2o)
  mat(6,4) = wmole(ih2s)
  mat(6,5) = wmole(ihcn)
  mat(6,6) = wmole(inh3)
  mat(6,7) = wmolat(iath)
!
!==============================================================================
! Correction des correlations teta - soufre et omega - azote.
  mat(6,8) = -(wmolat(iatc)+alpha(icha)*wmolat(iath)    &
                           +beta (icha)*wmolat(iato)    &
                           +teta (icha)*wmolat(iatn)    &
                           +omega(icha)*wmolat(iats) )
!==============================================================================
! Second membre

  sm(1) = alpha(icha)-gamma(icha)
  sm(2) = beta (icha)-delta(icha)
  sm(3) = teta (icha)-kappa(icha)
  sm(4) = 0.d0
  sm(5) = omega(icha)-zeta(icha)
  sm(6) = 0.d0

! On complete la matrice et le second menbre en fonction de
! la valeur de IY2CH
!
  if ( iy2ch(icha).eq.0 ) then
! Valeur de X2
    chx2(icha) = 2.d0
! Matrice
    mat(7,1) = -chx2(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
! Second membre
    sm(7) = 0.d0
    sm(8) = 0.d0

  else if ( iy2ch(icha).eq.1 ) then

! Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0

    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
! Second membre
    sm(7) = y2ch(icha)
    sm(8) = 0.d0

  else if ( iy2ch(icha).eq.2 ) then
! Valeur de X2
    chx2(icha) = 2.d0

! Matrice
    mat(7,1) = -chx2(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
! Second membre
    sm(7) = 0.d0
    sm(8) = y2ch(icha)

  endif

  call coal_resol_matrice( ndim, mat , sm , solu , ierr)

  if ( ierr .eq. 0 ) then
    a2(icha) = solu(1)
    b2(icha) = solu(2)
    c2(icha) = solu(3)
    d2(icha) = solu(4)
    e2(icha) = solu(5)
    f2(icha) = solu(6)
    if ( iy2ch(icha).eq.0 ) then
      y2ch(icha) = solu(8)
    else if ( iy1ch(icha).eq.1 ) then
      chx2(icha) = solu(7)/solu(1)
    endif
  else
    write(nfecra,9982) icha
    call csexit(1)
  endif
enddo

! --> Calcul de EHGAZE et de WMOLE
!     pour les especes CH(CHX1) et CH(CHX2)

! ---- Especes CH(CHX1)

do icha = 1, ncharb
  wmole(ichx1c(icha)) = wmolat(iatc)+chx1(icha)*wmolat(iath)
  if ( iy1ch(icha).eq.0 .or. iy1ch(icha).eq.2 ) then
    do it = 1, npo
      ehgaze(ichx1c(icha),it) = ehgaze(ichx1,it)
    enddo
  else
!        On a suppose D(HDEV,1,ICHA) = 0
    do it = 1, npo
      den1 = a1(icha)*wmole(ichx1c(icha))
      ehgaze(ichx1c(icha),it) =                                   &
      ( ( ehsoli(ich(icha),it) -                                  &
          (1.d0-y1ch(icha))*ehsoli(ick(icha),it) )                &
          * wmols(ich(icha))                                      &
        - b1(icha)*wmole(ico)*ehgaze(ico,it) )                    &
      / den1
    enddo
  endif
enddo

! ---- Especes CH(CHX2)

do icha = 1, ncharb
  wmole(ichx2c(icha)) = wmolat(iatc)+chx2(icha)*wmolat(iath)
  if ( iy2ch(icha).eq.0 .or. iy2ch(icha).eq.2 ) then
    do it = 1, npo
      ehgaze(ichx2c(icha),it) = ehgaze(ichx2,it)
    enddo
  else
!        On a suppose D(HDEV,2,ICHA) = 0
    do it = 1, npo
      den2 = a2(icha)*wmole(ichx2c(icha))
      ehgaze(ichx2c(icha),it) =                                   &
      ( ( ehsoli(ich(icha),it) -                                  &
           (1.d0-y2ch(icha))*ehsoli(ick(icha),it) )               &
           * wmols(ich(icha))                                     &
        - b2(icha)*wmole(ico)*ehgaze(ico,it) )                    &
      / den2
    enddo
  endif
enddo


! --> Calcul pour les differentes classes

do icla = 1, nclacp

! ---- Diametre min (m)

  dia2mn(icla) = zero

! ---- Masse volumique (kg/m3)

  rho20(icla)  = rho0ch(ichcor(icla))
  rho2mn(icla) = rho20(icla)*xashch(ichcor(icla))

! ---- Masse initiale de la particule (m)

  xmp0(icla) = rho20(icla)*pi*(diam20(icla)**3)/6.d0

! ---- Masse de cendres de la particule (m)

  xmash(icla) = xmp0(icla)*xashch(ichcor(icla))
!
enddo
!
!==============================================================================
!==============================================================================
! Nouveau modele de NOx (de Marcus et Sandro)
!==============================================================================
!==============================================================================
!
! Initialisation
ychx1t   = zero
ychx2t   = zero
!
! Composition du charbon reactif sous la forme
! CH(nalpha)O(nbeta)N(nteta)S(nomega)
!
do icha = 1, ncharb
!
  nalpha(icha) = (hch(icha)/wmolat(iath))                                     &
                 / (cch(icha)/wmolat(iatc))
  nbeta(icha)  = (och(icha)/wmolat(iato))                                     &
                 / (cch(icha)/wmolat(iatc))
  nteta(icha)  = (nch(icha)/wmolat(iatn))                                     &
                 / (cch(icha)/wmolat(iatc))
  nomega(icha) = (sch(icha)/wmolat(iats))                                     &
                 / (cch(icha)/wmolat(iatc))
!
enddo
!
! Composition du char sous la forme
! CH(NOXGAMMA)O(NOXDELTA)N(NOXKAPPA)S(NOXZETA)
!
do icha = 1, ncharb

  ! Fraction massique d'azote du charbon (un seul calcul au debut)
  nnch(icha) = (nteta(icha)*wmolat(iatn))/                                    &
                 (wmolat(iatc)                                                &
                  +nalpha(icha)*wmolat(iath)                                  &
                  +nbeta (icha)*wmolat(iato)                                  &
                  +nteta (icha)*wmolat(iatn)                                  &
                  +nomega(icha)*wmolat(iats))

! Test sur la composition du char. Si toutes les coefficients valent
! zero, il faut les calculer automatiquement.
  if(cck(icha).eq.0.d0.and.hck(icha).eq.0.d0.and.ock(icha)&
       .eq.0.d0.and.nck(icha).eq.0.d0.and.sck(icha).eq.0.d0) then

    !         Fraction massique d'azote dans le char (en fonction de Y1)
    nnckle(icha) = repnle(icha)*nnch(icha)/                             &
         (1.d0-y1ch(icha))

    !         Fraction massique d'azote dans le char (en fonction de Y2)
    nncklo(icha) = repnlo(icha)*nnch(icha)/                             &
         (1.d0-y2ch(icha))

    !         En cas de liberation de HCN au cours de la combustion heterogene, la
    !         fraction massique d'hydrogene dans le char est (en fonction de Y1):
    nhckle(icha)  = repnck(icha)*nnckle(icha)*wmolat(iath)/wmolat(iatn)

    !         En cas de liberation de HCN au cours de la combustion heterogene, la
    !         fraction massique d'hydrogene dans le char est (en fonction de Y2):
    nhcklo(icha)  = repnck(icha)*nncklo(icha)*wmolat(iath)/wmolat(iatn)

    !         Fraction massique de carbone du char (en fonction de Y1)
    ncckle(icha)  = 1.d0 - nnckle(icha) - nhckle(icha)

    !         Fraction massique de carbone du char (en fonction de Y2)
    nccklo(icha)  = 1.d0 - nncklo(icha) - nhcklo(icha)

    !         Coeffcients (en fonction de Y1)
    nkapp1(icha) = (wmolat(iatc)-wmolat(iatc)*ncckle(icha))/            &
         (ncckle(icha)*(repnck(icha)*wmolat(iath) +        &
         wmolat(iatn)))

    !         Coeffcients (en fonction de Y2)
    nkapp2(icha) = (wmolat(iatc)-wmolat(iatc)*nccklo(icha))/            &
         (nccklo(icha)*(repnck(icha)*wmolat(iath) +        &
         wmolat(iatn)))

    ngama1(icha) = repnck(icha)*nkapp1(icha)


    ngama2(icha) = repnck(icha)*nkapp2(icha)

    ndelt1(icha) = 0.d0
    ndelt2(icha) = 0.d0

    nzeta1 (icha) = 0.d0
    nzeta2 (icha) = 0.d0

    ! Composition donnee par l'utilisateur
  else

    ! Il faut assurer que le coefficient du carbon ne vaut pas zero.
    if(cck(icha).eq.0.d0) then
      write(nfecra,9100)
      call csexit(1)
    endif

    ngama1(icha)    = (hck(icha)/wmolat(iath))                              &
                        / (cck(icha)/wmolat(iatc))
    ngama2(icha)    = ngama1(icha)

    ndelt1(icha)    = (ock(icha)/wmolat(iato))                              &
                        / (cck(icha)/wmolat(iatc))
    ndelt2(icha)    = ndelt1(icha)

    nkapp1(icha)    = (nck(icha)/wmolat(iatn))                              &
                        / (cck(icha)/wmolat(iatc))
    nkapp2(icha)    = nkapp1(icha)

    nzeta1 (icha)    = (sck(icha)/wmolat(iats))                             &
                        / (cck(icha)/wmolat(iatc))
    nzeta2 (icha)    =  nzeta1 (icha)

  endif

  ! Les coefficientes de la premiere reaction de pyrolyse.

  ! Ligne 1
  mat(1,1) = -ngama1(icha)
  mat(1,2) = -ngama1(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.d0
  mat(1,5) = 1.d0-ngama1(icha)
  mat(1,6) = 3.d0
  mat(1,7) = 1.d0
  mat(1,8) = 0.d0
  ! Ligne 2
  mat(2,1) = -ndelt1(icha)
  mat(2,2) = 1.d0-ndelt1(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
  mat(2,5) = -ndelt1(icha)
  mat(2,6) = 0.d0
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0
  ! Ligne 3
  mat(3,1) = -nkapp1(icha)
  mat(3,2) = -nkapp1(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-nkapp1(icha)
  mat(3,6) = 1.d0
  mat(3,7) = 0.d0
  mat(3,8) = 0.d0
  !
  mat(4,1) = 0.d0
  mat(4,2) = 0.d0
  mat(4,3) = 0.d0
  mat(4,4) = 0.d0
  mat(4,5) = crepn1(2,icha)
  mat(4,6) =-crepn1(1,icha)
  mat(4,7) = 0.d0
  mat(4,8) = 0.d0
  !
  mat(5,1) = -nzeta1(icha)
  mat(5,2) = -nzeta1(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
  mat(5,5) = -nzeta1(icha)
  mat(5,6) = 0.d0
  mat(5,7) = 0.d0
  mat(5,8) = 0.d0
  !
  mat(6,1) = wmolat(iatc)
  mat(6,2) = wmole(ico)
  mat(6,3) = wmole(ih2o)
  mat(6,4) = wmole(ih2s)
  mat(6,5) = wmole(ihcn)
  mat(6,6) = wmole(inh3)
  mat(6,7) = wmolat(iath)
  mat(6,8) = -(wmolat(iatc)+nalpha(icha)*wmolat(iath)                        &
       +nbeta (icha)*wmolat(iato)                        &
       +nteta (icha)*wmolat(iatn)                        &
       +nomega(icha)*wmolat(iats))

  ! Second membre

  sm(1) = nalpha(icha)-ngama1(icha)
  sm(2) = nbeta (icha)-ndelt1(icha)
  sm(3) = nteta (icha)-nkapp1(icha)
  sm(4) = 0.d0
  sm(5) = nomega(icha)-nzeta1(icha)
  sm(6) = 0.d0

  ! On complete la matrice et le second membre en fonction de
  ! la valeur de IY1CH

  if ( iy1ch(icha).eq.0 ) then
    ! Valeur de X1
    nchx1(icha) = 4.d0
    ! Matrice
    mat(7,1) = -nchx1(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = 0.d0
    sm(8) = 0.d0

  else if ( iy1ch(icha).eq.1 ) then

    ! Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = y1ch(icha)
    sm(8) = 0.d0

  else if ( iy1ch(icha).eq.2 ) then
    ! Valeur de X1
    nchx1(icha) = 4.d0

    ! Matrice
    mat(7,1) = -nchx1(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
    ! Second membre
    sm(7) = 0.d0
    sm(8) = y1ch(icha)

  endif

  ndim = 8

  call coal_resol_matrice( ndim, mat , sm , solu , ierr)

  if ( ierr .eq. 0 ) then
    noxa1(icha) = solu(1)
    noxb1(icha) = solu(2)
    noxc1(icha) = solu(3)
    noxd1(icha) = solu(4)
    noxe1(icha) = solu(5)
    noxf1(icha) = solu(6)
    if ( iy1ch(icha).eq.0 ) then
      ny1ch(icha) = solu(8)
    else if ( iy1ch(icha).eq.1 ) then
      nchx1(icha) = solu(7)/solu(1)
    endif
  else
    write(nfecra,9982) icha
    call csexit(1)
  endif

  ! Masse molaire de CHX1
  wchx1c(icha) = wmolat(iatc)+nchx1(icha)*wmolat(iath)

  ! Les fractions massiques de HCN,NH3,CHx1 dans les matieres volatiles legeres.
  yhcnle(icha) = (noxe1(icha)*wmole(ihcn)) /                                 &
       (noxa1(icha)*(wmolat(iatc)+nchx1(icha)*wmolat(iath))      &
                    +noxb1(icha)*wmole(ico)                                   &
                    +noxc1(icha)*wmole(ih2o)                                  &
                    +noxd1(icha)*wmole(ih2s)                                  &
                    +noxe1(icha)*wmole(ihcn)                                  &
                    +noxf1(icha)*wmole(inh3))

  ynh3le(icha) = (noxf1(icha)*wmole(inh3)) /                                 &
                    (noxa1(icha)*(wmolat(iatc)+nchx1(icha)*wmolat(iath))      &
                    +noxb1(icha)*wmole(ico)                                   &
                    +noxc1(icha)*wmole(ih2o)                                  &
                    +noxd1(icha)*wmole(ih2s)                                  &
                    +noxe1(icha)*wmole(ihcn)                                  &
                    +noxf1(icha)*wmole(inh3))

  ychxle(icha)   = (noxa1(icha)*(wmolat(iatc)+nchx1(icha)*wmolat(iath)))/    &
                    (noxa1(icha)*(wmolat(iatc)+nchx1(icha)*wmolat(iath))      &
                    +noxb1(icha)*wmole(ico)                                   &
                    +noxc1(icha)*wmole(ih2o)                                  &
                    +noxd1(icha)*wmole(ih2s)                                  &
                    +noxe1(icha)*wmole(ihcn)                                  &
                    +noxf1(icha)*wmole(inh3))

! Les coefficients des produits de la reaction heterogene (Char1)
  noxj1(icha)   =  ngama1(icha)
  noxi1(icha)   =  1.d0 - noxj1(icha)
  noxk1(icha)   =  nkapp1(icha) - noxj1(icha)
  noxh1(icha)   =  (noxi1(icha) + noxk1(icha))/2.d0

! Les fractions massiques de HCN, NO et CO par kg Char1.
  ycoch1(icha) = (noxi1(icha)*wmole(ico))  /                                 &
                  (noxi1(icha)*wmole(ico)+noxj1(icha)*wmole(ihcn)+            &
                   noxk1(icha)*3.d-2)
  yhcnc1(icha) = (noxj1(icha)*wmole(ihcn)) /                                 &
                  (noxi1(icha)*wmole(ico)+noxj1(icha)*wmole(ihcn)+            &
                   noxk1(icha)*3.d-2)
  ynoch1(icha) = (noxk1(icha)*3.d-2)       /                                 &
                  (noxi1(icha)*wmole(ico)+noxj1(icha)*wmole(ihcn)+            &
                   noxk1(icha)*3.d-2)

! Les coefficientes de la deuxieme reaction de pyrolyse.

  mat(1,1) = -ngama2(icha)
  mat(1,2) = -ngama2(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.d0
  mat(1,5) = 1.d0-ngama2(icha)
  mat(1,6) = 3.d0
  mat(1,7) = 1.d0
  mat(1,8) = 0.d0
  ! Ligne 2
  mat(2,1) = -ndelt2(icha)
  mat(2,2) = 1.d0-ndelt2(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
  mat(2,5) = -ndelt2(icha)
  mat(2,6) = 0.d0
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0

  mat(3,1) = -nkapp2(icha)
  mat(3,2) = -nkapp2(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-nkapp2(icha)
  mat(3,6) = 1.d0
  mat(3,7) = 0.d0
  mat(3,8) = 0.d0

  mat(4,1) = 0.d0
  mat(4,2) = 0.d0
  mat(4,3) = 0.d0
  mat(4,4) = 0.d0
  mat(4,5) = crepn2(2,icha)
  mat(4,6) =-crepn2(1,icha)
  mat(4,7) = 0.d0
  mat(4,8) = 0.d0

  mat(5,1) = -nzeta2(icha)
  mat(5,2) = -nzeta2(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
  mat(5,5) = -nzeta2(icha)
  mat(5,6) = 0.d0
  mat(5,7) = 0.d0
  mat(5,8) = 0.d0

  mat(6,1) = wmolat(iatc)
  mat(6,2) = wmole(ico)
  mat(6,3) = wmole(ih2o)
  mat(6,4) = wmole(ih2s)
  mat(6,5) = wmole(ihcn)
  mat(6,6) = wmole(inh3)
  mat(6,7) = wmolat(iath)
  mat(6,8) = -(wmolat(iatc)+nalpha(icha)*wmolat(iath)                        &
                            +nbeta (icha)*wmolat(iato)                        &
                            +nteta (icha)*wmolat(iatn)                        &
                            +nomega(icha)*wmolat(iats) )


  ! Second membre

  sm(1) = nalpha(icha)-ngama2(icha)
  sm(2) = nbeta (icha)-ndelt2(icha)
  sm(3) = nteta (icha)-nkapp2(icha)
  sm(4) = 0.d0
  sm(5) = nomega(icha)-nzeta2(icha)
  sm(6) = 0.d0

  ! On complete la matrice et le second menbre en fonction de
  ! la valeur de IY2CH

  if ( iy2ch(icha).eq.0 ) then
    ! Valeur de X2
    nchx2(icha) = 2.d0
    ! Matrice
    mat(7,1) = -nchx2(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = 0.D0
    sm(8) = 0.D0

  else if ( iy2ch(icha).eq.1 ) then

    ! Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
    ! Second membre
    sm(7) = y2ch(icha)
    sm(8) = 0.d0

  else if ( iy2ch(icha).eq.2 ) then
    ! Valeur de X2
    nchx2(icha) = 2.D0

    ! Matrice
    mat(7,1) = -nchx2(icha)
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0

    mat(8,1) = 0.d0
    mat(8,2) = 0.d0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
    ! Second membre
    sm(7) = 0.d0
    sm(8) = y2ch(icha)

  endif

  call coal_resol_matrice( ndim, mat , sm , solu , ierr)

  if ( ierr .eq. 0 ) then
    noxa2(icha) = solu(1)
    noxb2(icha) = solu(2)
    noxc2(icha) = solu(3)
    noxd2(icha) = solu(4)
    noxe2(icha) = solu(5)
    noxf2(icha) = solu(6)
    if ( iy2ch(icha).eq.0 ) then
      ny2ch(icha) = solu(8)
    else if ( iy1ch(icha).eq.1 ) then
      nchx2(icha) = solu(7)/solu(1)
    endif
  else
    write(nfecra,9982) icha
    call csexit(1)
  endif

  ! Masse molaire de CHX2
  wchx2c(icha) = wmolat(iatc)+nchx2(icha)*wmolat(iath)

  ! Les fractions massiques de HCN,NH3,CHx2 dans les matieres volatiles lourdes.
  yhcnlo(icha) = (noxe2(icha)*wmole(ihcn))/                                  &
                   (noxa2(icha)*(wmolat(iatc)+nchx2(icha)*wmolat(iath))       &
                   +noxb2(icha)*wmole(ico)                                    &
                   +noxc2(icha)*wmole(ih2o)                                   &
                   +noxd2(icha)*wmole(ih2s)                                   &
                   +noxe2(icha)*wmole(ihcn)                                   &
                   +noxf2(icha)*wmole(inh3))

  ynh3lo(icha) = (noxf2(icha)*wmole(inh3))/                                  &
                   (noxa2(icha)*(wmolat(iatc)+nchx2(icha)*wmolat(iath))       &
                   +noxb2(icha)*wmole(ico)                                    &
                   +noxc2(icha)*wmole(ih2o)                                   &
                   +noxd2(icha)*wmole(ih2s)                                   &
                   +noxe2(icha)*wmole(ihcn)                                   &
                   +noxf2(icha)*wmole(inh3))

  ychxlo(icha) = (noxa2(icha)*(wmolat(iatc)+nchx2(icha)*wmolat(iath)))/      &
                   (noxa2(icha)*(wmolat(iatc)+nchx2(icha)*wmolat(iath))       &
                   +noxb2(icha)*wmole(ico)                                    &
                   +noxc2(icha)*wmole(ih2o)                                   &
                   +noxd2(icha)*wmole(ih2s)                                   &
                   +noxe2(icha)*wmole(ihcn)                                   &
                   +noxf2(icha)*wmole(inh3))

  ! Les coefficients des produits de la reaction heterogene (Char2)
  noxj2(icha)   =  ngama2(icha)
  noxi2(icha)   =  1.d0 - noxj2(icha)
  noxk2(icha)   =  nkapp2(icha) - noxj2(icha)
  noxh2(icha)   =  (noxi2(icha) + noxk2(icha))/2.d0

  ! Les fractions massiques de HCN, NO et CO par kg Char2.
  ycoch2(icha) = (noxi2(icha)*wmole(ico))  /                                 &
                  (noxi2(icha)*wmole(ico)+noxj2(icha)*wmole(ihcn)+            &
                   noxk2(icha)*3.d-2)
  yhcnc2(icha) = (noxj2(icha)*wmole(ihcn)) /                                 &
                  (noxi2(icha)*wmole(ico)+noxj2(icha)*wmole(ihcn)+            &
                   noxk2(icha)*3.d-2)
  ynoch2(icha) = (noxk2(icha)*3.d-2)       /                                 &
                  (noxi2(icha)*wmole(ico)+noxj2(icha)*wmole(ihcn)+            &
                   noxk2(icha)*3.d-2)
enddo

 ! Masse molaires des matieres volatiles legeres et lourdes
do icha = 1, ncharb
  ychx1t = ychx1t + ychxle(icha)
  ychx2t = ychx2t + ychxlo(icha)
enddo

do icha = 1, ncharb
  wmchx1 = wmchx1 + ( ychxle(icha)/ychx1t * wchx1c(icha) )
  wmchx2 = wmchx2 + ( ychxlo(icha)/ychx2t * wchx2c(icha) )
enddo

! Discretisation de la temperature

teno(1) =  300.d0
teno(2) =  500.d0
teno(3) = 1000.d0
teno(4) = 1500.d0
teno(5) = 2000.d0
teno(6) = 2500.d0
teno(7) = 3000.d0

kf1(1)  = 1.58d+12
kf1(2)  = 1.85d+12
kf1(3)  = 1.68d+12
kf1(4)  = 1.45d+12
kf1(5)  = 1.26d+12
kf1(6)  = 1.13d+12
kf1(7)  = 1.02d+12

kf2(1)  = 4.10d+13
kf2(2)  = 4.10d+13
kf2(3)  = 4.10d+13
kf2(4)  = 4.10d+13
kf2(5)  = 4.10d+13
kf2(6)  = 4.10d+13
kf2(7)  = 4.10d+13

kf3(1)  = 1.90d+13
kf3(2)  = 1.90d+13
kf3(3)  = 1.90d+13
kf3(4)  = 1.90d+13
kf3(5)  = 1.90d+13
kf3(6)  = 1.90d+13
kf3(7)  = 1.90d+13

kf4(1)  = 8.62d+04
kf4(2)  = 2.84d+08
kf4(3)  = 2.04d+11
kf4(4)  = 2.43d+12
kf4(5)  = 9.61d+12
kf4(6)  = 2.38d+13
kf4(7)  = 4.60d+13

kr4(1)  = 2.05d+04
kr4(2)  = 3.36d+07
kr4(3)  = 1.18d+10
kr4(4)  = 1.20d+11
kr4(5)  = 4.88d+11
kr4(6)  = 1.33d+12
kr4(7)  = 2.92d+12

kf5(1)  = 5.80d+07
kf5(2)  = 4.98d+09
kf5(3)  = 2.31d+11
kf5(4)  = 1.10d+12
kf5(5)  = 2.74d+12
kf5(6)  = 5.14d+12
kf5(7)  = 8.26d+12

kr5(1)  = 1.79d+01
kr5(2)  = 4.76d+05
kr5(3)  = 1.73d+09
kr5(4)  = 3.78d+10
kr5(5)  = 2.11d+11
kr5(6)  = 6.57d+11
kr5(7)  = 1.50d+12

kf6(1)  = 1.86d+14
kf6(2)  = 1.91d+14
kf6(3)  = 1.90d+14
kf6(4)  = 1.81d+14
kf6(5)  = 1.72d+14
kf6(6)  = 1.66d+14
kf6(7)  = 1.63d+14

kr6(1)  = 5.86d+11
kr6(2)  = 4.72d+12
kr6(3)  = 2.26d+13
kr6(4)  = 3.80d+13
kr6(5)  = 4.94d+13
kr6(6)  = 5.78d+13
kr6(7)  = 6.42d+13

kf7(1)  = 1.65d+14
kf7(2)  = 1.65d+14
kf7(3)  = 1.65d+14
kf7(4)  = 1.65d+14
kf7(5)  = 1.65d+14
kf7(6)  = 1.65d+14
kf7(7)  = 1.65d+14

kr7(1)  = 3.23d-03
kr7(2)  = 2.40d+04
kr7(3)  = 3.47d+09
kr7(4)  = 1.89d+11
kr7(5)  = 1.45d+12
kr7(6)  = 5.03d+12
kr7(7)  = 1.17d+13

pflue   = 1.d-4

do ii = 1, 7

  ! Chi2 est calculé selon la cinetique donnee dans le manual de FLUENT.
  chi2(ii) = (4.52e5*(teno(ii)**1.6d0)                          &
              *exp(-80815.d0/cs_physical_constants_r/teno(ii))) &
            /(1.02e5*(teno(ii)**1.6d0)                          &
              *exp(-13802.d0/cs_physical_constants_r/teno(ii)))

  ! JJ indique le rapport H/C du combustible (4=CH4;3=CH3,etc.)
  do jj = 1,4

    if(jj.eq.1) then

      ka(jj,ii) = kf1(ii) * ( kr6(ii) / kf6(ii) ) * 1.d-6 * pflue
      kb(jj,ii) = kf2(ii) * 1.d-6 * pflue
      kc(jj,ii) = kf3(ii) * ( kf7(ii) )/( kr7(ii) ) * 1.d-6 * pflue

    elseif(jj.eq.2) then

      ka(jj,ii) = kf1(ii) * 1.d-6 * pflue
      kb(jj,ii) = kf2(ii) * ( kf6(ii) )/( kr6(ii) ) * 1.d-6 * pflue
      kc(jj,ii) = kf3(ii) * ( kf6(ii) * kf7(ii) )/( kr6(ii) * kr7(ii) )     &
                    * 1.d-6 * pflue

    elseif(jj.eq.3) then

      ka(jj,ii) = kf1(ii) * ( kf5(ii) )/( kr5(ii) ) * 1.d-6 * pflue
      kb(jj,ii) = kf2(ii) * ( kf5(ii) * kf6(ii) )/( kr5(ii) * kr6(ii) )     &
                    * 1.d-6 * pflue
      kc(jj,ii) = kf3(ii) * ( kf5(ii) * kf6(ii) * kf7(ii) )/( kr5(ii)       &
                    * kr6(ii) * kr7(ii) ) * 1.d-6 * pflue

    elseif(jj.eq.4) then

      ka(jj,ii) = kf1(ii) * ( kf4(ii) * kf5(ii) )/( kr4(ii) * kr5(ii) )     &
                    * 1.d-6 * pflue
      kb(jj,ii) = kf2(ii) * ( kf4(ii) * kf5(ii) * kf6(ii) )/( kr4(ii)       &
                    * kr5(ii) * kr6(ii) ) * 1.d-6 * pflue
      kc(jj,ii) = kf3(ii) * ( kf4(ii) * kf5(ii) * kf6(ii) * kf7(ii) )/      &
                    ( kr4(ii) * kr5(ii) * kr6(ii) * kr7(ii) ) * 1.d-6 * pflue

    endif

  enddo

enddo

! AFFICHAGE RECAPITULATIF
! =======================
! Coals data summary
write(nfecra,8000)

!  Constante du modele de devolatilisation

write(nfecra,8100)
do icha=1,ncharb
  ! Coal number
  write(nfecra,8011)  icha
  ! Composition normee du Charbon
  write(nfecra,8029)
  molcch = cch(icha)/cch(icha)
  write(nfecra,8030) molcch, alpha(icha), beta(icha), teta(icha), omega(icha)
  ! Composition normee du Char
  write(nfecra,8031)
  if (cck(icha).eq.0.d0) then
    molcck = 1.d0
  else
    molcck = cck(icha)/cck(icha)
  endif
  write(nfecra,8030) molcck, gamma(icha), delta(icha), kappa(icha), zeta(icha)
  ! Devolatilization model options
  write(nfecra,8100)
  if ( iy1ch(icha) .eq.0 ) then
!   Option 0
    write(nfecra,8102)
  else if ( iy1ch(icha) .eq.1 ) then
    !   Option 1
    write(nfecra,8103)
  else if  ( iy1ch(icha) .eq.2 ) then
    !   Option 2
    write(nfecra,8104)
  endif
! Y1 et Y2
  write(nfecra,8105) y1ch(icha),y2ch(icha)
! Volatile materials Composition
  write(nfecra,8010)
! Coal number, CHx1, CHx2
  write(nfecra,8051)
  ierror = 0
  write(nfecra,8052) icha,chx1(icha),chx2(icha)
  if ( chx1(icha) .le. 0.d0 .or. chx2(icha) .le. 0 ) then
    ierror = 1
  endif
  if ( ierror .eq. 1 ) then
    write(nfecra,9970)
    call csexit(1)
  endif
! Composition MV
  write(nfecra,8012)
  write(nfecra,8013) a1(icha),b1(icha),c1(icha),d1(icha),e1(icha),f1(icha)
  write(nfecra,8014) a2(icha),b2(icha),c2(icha),d2(icha),e2(icha),f2(icha)
enddo

!Model de NOx
write(nfecra,8034)
write(nfecra, 8024) ieqnox, imdnox, irb
! Numero Charbon et Option Y
do icha=1,ncharb
  ! Coal number
  write(nfecra,8011)  icha
  ! Composition du Charbon
  write(nfecra,8029)
  molcch = cch(icha)/cch(icha)
  write(nfecra,8030) molcch,nalpha(icha),nbeta(icha),nteta(icha),nomega(icha)
  ! Composition du Char Reaction 1
  write(nfecra,8036)
  if (cck(icha).eq.0.d0) then
    molcck = 1.d0
  else
    molcck = cck(icha)/cck(icha)
  endif
  write(nfecra,8030) molcck,ngama1(icha),ndelt1(icha),nkapp1(icha),nzeta1(icha)
  ! Composition du Char Reaction 2
  write(nfecra,8039)
  if (cck(icha).eq.0.d0) then
    molcck = 1.d0
  else
    molcck = cck(icha)/cck(icha)
  endif
  write(nfecra,8030) molcck,ngama2(icha),ndelt2(icha),nkapp2(icha),nzeta2(icha)
  ! Devolatilization model options
  write(nfecra,8100)
  if ( iy1ch(icha) .eq.0 ) then
    write(nfecra,8102)
  else if ( iy1ch(icha) .eq.1 ) then
    write(nfecra,8103)
  else if  ( iy1ch(icha) .eq.2 ) then
    write(nfecra,8104)
  endif
!  Y1 et Y2
  write(nfecra,8105) y1ch(icha),y2ch(icha)
  write(nfecra,8010)
  write(nfecra,8051)
  ierror = 0
! Composition MV
  write(nfecra,8052) icha,nchx1(icha),nchx2(icha)
  if ( nchx1(icha) .le. 0.d0 .or. nchx2(icha) .le. 0 ) then
    ierror = 1
  endif
  if ( ierror .eq. 1 ) then
    write(nfecra,9970)
    call csexit(1)
  endif
! Composition des MV
  write(nfecra,8012)
  write(nfecra,8013) noxa1(icha),noxb1(icha),noxc1(icha),noxd1(icha),         &
                     noxe1(icha),noxf1(icha)
  write(nfecra,8014) noxa2(icha),noxb2(icha),noxc2(icha),noxd2(icha),         &
                     noxe2(icha),noxf2(icha)
! Produits reaction heterogene Reaction 1
   write(nfecra,8040)
   write(nfecra,8041)
   write(nfecra,8042) ycoch1(icha),yhcnc1(icha),ynoch1(icha)
   write(nfecra,8043) ycoch2(icha),yhcnc2(icha),ynoch2(icha)
enddo

!  Calcul de DHdev

write(nfecra,8001)
do icha=1,ncharb

  !    Enthalpie du Charbon

  hchar = ehsoli(ich(icha),1)
  hcoke = ehsoli(ick(icha),1)

  !    Enthalpie des matieres volatiles

  ehvol1 = ( a1(icha)*ehgaze(ichx1c(icha),1)*wmole(ichx1c(icha))   &
           +b1(icha)*ehgaze(ico ,1)*wmole(ico)                    &
           +c1(icha)*ehgaze(ih2o,1)*wmole(ih2o)                   &
           +d1(icha)*ehgaze(ih2s,1)*wmole(ih2s) )                 &
         /( a1(icha)*wmole(ichx1c(icha))                          &
           +b1(icha)*wmole(ico)                                   &
           +c1(icha)*wmole(ih2o)                                  &
           +d1(icha)*wmole(ih2s) )

 ehvol2 = ( a2(icha)*ehgaze(ichx2c(icha),1)*wmole(ichx2c(icha))   &
           +b2(icha)*ehgaze(ico ,1)*wmole(ico)                    &
           +c2(icha)*ehgaze(ih2o,1)*wmole(ih2o)                   &
           +d2(icha)*ehgaze(ih2s,1)*wmole(ih2s) )                 &
         /( a2(icha)*wmole(ichx2c(icha))                          &
           +b2(icha)*wmole(ico)                                   &
           +c2(icha)*wmole(ih2o)                                  &
           +d2(icha)*wmole(ih2s) )

 dhvol1 = hchar-y1ch(icha)*ehvol1-(1.d0-y1ch(icha))*hcoke
 dhvol2 = hchar-y2ch(icha)*ehvol2-(1.d0-y2ch(icha))*hcoke

 write(nfecra,8002) icha,dhvol1,dhvol2,pcich(icha)

enddo

!  Loi enthalpie/temperature

write(nfecra,8020)
do icha = 1, ncharb
  write(nfecra,8011) icha
  write(nfecra,8021)
  do ii = 1, npoc
    write(nfecra,8022)thc(ii),ehsoli(ich (icha),ii)               &
                             ,ehsoli(ick (icha),ii)               &
                             ,ehsoli(iash(icha),ii)               &
                             ,ehsoli(iwat(icha),ii)
  enddo
enddo

!  Calcul des AiFj : nbre de mole de i par kg de j a l'origine

wmco   = wmole(ico)
wmo2   = wmole(io2)
wmco2  = wmole(ico2)
wmh2o  = wmole(ih2o)
wmn2   = wmole(in2)
wmc    = wmolat(iatc)

dmf3  = ( oxyo2 (1)*wmo2 +oxyn2 (1)*wmn2     &
         +oxyh2o(1)*wmh2o+oxyco2(1)*wmco2 )

if ( dmf3 .le. 0.d0 ) then
  write(nfecra,9896) oxyo2(1) ,oxyn2(1) ,    &
                     oxyh2o(1),oxyco2(1)
endif

af3(io2)  = oxyo2(1)  / dmf3
af3(in2)  = oxyn2(1)  / dmf3
af3(ih2o) = oxyh2o(1) / dmf3
af3(ico2) = oxyco2(1) / dmf3

if ( noxyd .ge. 2.d0 ) then
  dmf4  = ( oxyo2 (2)*wmo2 +oxyn2 (2)*wmn2                                    &
        +oxyh2o(2)*wmh2o+oxyco2(2)*wmco2 )
  if ( dmf4 .le. 0.d0 ) then
    write(nfecra,9897) oxyo2(2) ,oxyn2(2) ,                                   &
         oxyh2o(2),oxyco2(2)
    call csexit(1)
  endif

  af4(io2)  = oxyo2(2)  / dmf4
  af4(in2)  = oxyn2(2)  / dmf4
  af4(ih2o) = oxyh2o(2) / dmf4
  af4(ico2) = oxyco2(2) / dmf4

endif

if ( noxyd .eq. 3.d0 ) then
  dmf5  = ( oxyo2 (3)*wmo2 +oxyn2 (3)*wmn2                                    &
           +oxyh2o(3)*wmh2o+oxyco2(3)*wmco2 )
  if ( dmf5 .le. 0.d0 ) then
    write(nfecra,9898) oxyo2(3) ,oxyn2(3) ,                                   &
                       oxyh2o(3),oxyco2(3)
    call csexit(1)
  endif

  af5(io2)  = oxyo2(3)  / dmf5
  af5(in2)  = oxyn2(3)  / dmf5
  af5(ih2o) = oxyh2o(3) / dmf5
  af5(ico2) = oxyco2(3) / dmf5

endif
!vapeur
af6(ih2o)  = 1.d0/wmh2o
! coke par o2
af7(ico)   = 1.0d0/wmc
af7(io2)   =-0.5d0/wmc
! coke par co2
af8(ico)   = 2.0d0/wmc
af8(ico2)  =-1.0d0/wmc
! coke par h2o
af9(ih2o)  =-1.0d0/wmc
af9(ico)   = 1.0d0/wmc
af9(ihy)   = 1.0d0/wmc

return

!--------
! Formats
!--------

!================================================================
8000 format                             &
  (/,                                   &
   3X,'** coals data summary **', /,    &
   3X,'------------------------')
!================================================================
8001 format                                                  &
  (/,                                                        &
   '---', '-----------', /,                                  &
      3x, 'Delta Hdev ', /,                                  &
   '---', '-----------', /,                                  &
   '   coal numb. ', 4x,' light', 10x,' heavy',12x, 'PCI' , /, &
   '---------',                                              &
   '----------------------------------------------------' )
!================================================================
8002 format                                                  &
  (/,                                                        &
   3x, i6, 5x, 1e12.4, 5x, 1e12.4, 5x, 1e12.4       )
!================================================================
8010 format                                       &
  (/,                                             &
   3X,'** volatile materials composition **', /,  &
   3X,'------------------------------------')
!================================================================
8011 format                                       &
  (/,                                             &
   '---', '-------------------', /,               &
   3X,'Coal number:', 1X, i3   , /,               &
   '---', '-------------------', / )
!================================================================
8012 format                                                  &
  (/,                                                        &
   '---',                                                    &
   '--------------------------------------',                 &
   '--------------------------------------',              /, &
   12X, 'CHx', 8X,' CO', 8X,' H2O',10X, 'H2S' ,9X,           &
   'HCN', 9X, 'NH3 ', /,                                     &
   '---',                                                    &
   '--------------------------------------',                 &
   '--------------------------------------' )
!================================================================
8013 format                                                  &
  (/,                                                        &
   3X, 'MV1', 6e12.4)
!================================================================
8014 format                                                  &
  (/,                                                        &
   3X, 'MV2', 6e12.4)
!================================================================
8020 format                                  &
  (/,                                        &
   3X,'-------------------------------', /,  &
   3X,'** Enthalpy/Temperature law  **', /,  &
   3X,'-------------------------------')
!================================================================
8021 format                                                  &
  (/,                                                        &
   '---',                                                    &
   '--------------------------------------',                 &
   '--------------------------------------',              /, &
   10X, 'Temp', 8X,' h_ch', 8X,' h_coke',10X,                &
   'h_ash' ,9X, 'h_wat' , /, &
   '---',                                                    &
   '--------------------------------------',                 &
   '--------------------------------------' )
!================================================================
8022 format                                                  &
  (/,                                                        &
   3X, 1e12.4, 3X, 1e12.4, 3X, 1e12.4, 3X, 1e12.4, 3X, 1e12.4)
!================================================================
 8051 format(/,5X,'Coal number',7X,'CHX1',12X,'CHX2')
 8052 format( 10x, i3, 8x,1e12.4,1x,1e12.4)
!================================================================
 8100 format(/,3X,'---------------------------------' , /, &
               3X,' Devolatization model constant'    , /, &
               3X,'---------------------------------')
!================================================================
 8102 format(/,8X,'Option 0: Y1 and Y2 computed  ')
 8103 format(/,8X,'Option 1: Y1 and Y2 fixed     ')
 8104 format(/,8X,'Option 2: CHX1 and CHX2 fixed ')
 8105 format(  8X,'Y1_ch = ', 1e12.4,3X,'Y2_ch = ', 1e12.4)
!====================================================================
 8034 format(/,/,3X,'---------------------------------' , /, &
                 3X,' NOx model                       ' , /, &
                 3X,'---------------------------------')
!====================================================================
 8024 format(/,/,3X,'---------------------------------' , /, &
                 3X,'IEQNOX =', i3                      , /, &
                 3X,'NOx formation Features:'           , /, &
                 3X,'IMDNOX =', i3                      , /, &
                 3X,'Reburning:'                        , /, &
                 3X,'IRB    =', i3                      , /, &
                 3X,'---------------------------------')
!====================================================================
 8029 format(/,3X,'---------------------------------' , /,          &
               3X,' Normalized coal composition'      , /,          &
               3X,'---------------------------------')
 8030 format(8X,'C:',1e12.4,3X,'H:',1e12.4,3X,'O:',1e12.4,3X,       &
                'N:',1e12.4,3X,'S:',1e12.4)
!====================================================================
 8031 format(/,3X,'---------------------------------' , /,          &
               3X,' Normalized char composition'      , /,          &
               3X,'---------------------------------')
!====================================================================
 8036 format(/,3X,'-----------------------------------------' , /,  &
               3X,' Normalized char composition (Reaction 1)' , /,  &
               3X,'-----------------------------------------')
!====================================================================
 8039 format(/,3X,'-----------------------------------------' , /,  &
               3X,' Normalized char composition (Reaction 2)' , /,  &
               3X,'-----------------------------------------')
!====================================================================
 8040 format(/,3X,'---------------------------------' , /, &
               3X,' Heterogeneous reaction products'  , /, &
               3X,'---------------------------------')
 8041 format(8X,'CO',12X,'HCN',12X,'NO')
 8042 format(3X,'R1:',1e12.4,1e12.4,1e12.4)
 8043 format(3X,'R2:',1e12.4,1e12.4,1e12.4)
!====================================================================

 9970 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  Some values of CHX1 and CHX2 are negative or zero.        ',/,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9971 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  Nitrogen was chosen in the coal composition but the       ',/,&
'@  distribution between HCN and NH3 is zero for coal ', i2,    /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9980 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  A different option is used for Y1 and Y2 for coal ', i2,    /,&
'@         IY1CH = ', i2,                                       /,&
'@         IY2CH = ', i2,                                       /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9981 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  An unavailable option is used for Y1 and Y2 for coal ', i2, /,&
'@         IY1CH = ', i2,                                       /,&
'@         IY2CH = ', i2,                                       /,&
'@',                                                            /,&
'@  only values 0 , 1 or 2 are admissible.',                    /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9982 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  Non invertible matrix for volatile matter computation',     /,&
'@  for coal ', i2,                                             /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9991 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  The current number of species must be less than',           /,&
'@                                        or equal to ', i10,   /,&
'@  It is ', i10,' here.',                                      /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9992 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  The number of tabulation points is limited to ', i10,       /,&
'@  It is ', i10,' here.',                                      /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9993 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  The number of elementary species is limited to ', i10,      /,&
'@  It is ', i10,' here.',                                      /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9896 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  LA COMPOSITION DE L''OXYDANT 1 EST ERRONEE                ',/,&
'@     O2  :  ',G15.7,'                                       ',/,&
'@     N2  :  ',G15.7,'                                       ',/,&
'@     H2O :  ',G15.7,'                                       ',/,&
'@     CO2 :  ',G15.7,'                                       ',/,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9897 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  The composition of oxidant 2 is incorrect',                 /,&
'@     O2  :  ', g15.7,                                         /,&
'@     N2  :  ', g15.7,                                         /,&
'@     H2O :  ', g15.7,                                         /,&
'@     CO2 :  ', g15.7,                                         /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9898 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@',                                                            /,&
'@  The composition of oxidant 2 is incorrect',                 /,&
'@     O2  :  ', g15.7,                                         /,&
'@     N2  :  ', g15.7,                                         /,&
'@     H2O :  ', g15.7,                                         /,&
'@     CO2 :  ', g15.7,                                         /,&
'@',                                                            /,&
'@  The computation CAN NOT run',                               /,&
'@',                                                            /,&
'@  Check the setup parameters.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@                Pulverized coal model',                       /,&
'@'                                                             /,&
'@  *** Wrong elementary coke balance ***',                     /,&
'@',                                                            /,&
'@  Check the mass percentages of each element.',               /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!----
! End
!----

end subroutine cs_coal_readata

!===============================================================================
!===============================================================================
!===============================================================================

subroutine coal_resol_matrice &
 ( ndim, aa , bb , xx , ierr)

!===============================================================================
! FONCTION :
! --------
!     finds the solution of the linear system A*XX=BB , dim(A)=NDIM
!     using Gauss elimination method with partial pivoting
!     A,B - modified data ; N - data ; X - result
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! DIM              ! r  ! <-- ! dimension du syteme                            !
! AA               ! r  ! <-->! matrice du systeme                             !
! BB               ! r  ! <-->! second menbre                                  !
! XX               ! r  ! --> ! variable                                       !
! IERR             ! r  ! --> ! gestion des erreurs                            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================
!===============================================================================
! MODULE
!===============================================================================
!
!*************************************************************
!
implicit none
!
!*************************************************************

!===============================================================================
! Arguments
integer                                 , intent(in)    :: ndim
double precision , dimension(ndim,ndim) , intent(inout) :: aa
double precision , dimension(ndim)      , intent(inout) :: bb
double precision , dimension(ndim)      , intent(out)   :: xx
integer                                 , intent(out)   :: ierr
! Variables locales
integer          :: ii,jj,kk,iw
double precision :: epsil,ww,pp,ss

!===============================================================================
! Parametres
EPSIL  = 1.D-10
IERR   = 0
!
! initialisation de la solution et test pour etre sur que les elements
! de la diagonale de la matrice sont non nuls
!
BCLEP: DO II=1,NDIM
         IW = II
         WW = ABS(AA(II,II))
         DO JJ=II,NDIM
           IF ( ABS(AA(JJ,II)) .GT. WW ) THEN
             IW = JJ
             WW = ABS( AA(JJ,II) )
           ENDIF
         ENDDO
         IF ( WW .LE. EPSIL ) THEN
           IERR = 1
           EXIT BCLEP
         ENDIF
!
         DO JJ=II,NDIM
           PP = AA(II,JJ)
           AA(II,JJ) = AA(IW,JJ)
           AA(IW,JJ) = PP
         ENDDO
!
         PP = BB(II)
         BB(II) = BB(IW)
         BB(IW) = PP
!
         DO JJ=II+1,NDIM
           PP = AA(JJ,II)/AA(II,II)
           DO KK=II+1,NDIM
             AA(JJ,KK) = AA(JJ,KK)-PP*AA(II,KK)
           ENDDO
           BB(JJ) = BB(JJ) - PP*BB(II)
        ENDDO
ENDDO BCLEP
!
IF ( IERR .NE. 1 ) THEN
  IF ( ABS(AA(NDIM,NDIM)) .LT. EPSIL ) THEN
    IERR = 1
  ELSE
    XX(NDIM) = BB(NDIM)/AA(NDIM,NDIM)
    DO II=NDIM-1,1,-1
      PP = 1.D0/AA(II,II)
      SS = 0.D0
      DO JJ=II+1,NDIM
        SS = SS+AA(II,JJ)*XX(JJ)
      ENDDO
      XX(II) = PP*(BB(II)-SS)
    ENDDO
  ENDIF
ENDIF
!
!----
! End
!----
end subroutine coal_resol_matrice
