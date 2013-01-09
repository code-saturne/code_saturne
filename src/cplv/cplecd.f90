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

subroutine cplecd
!================
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

!===============================================================================

implicit none

! Arguments

! Local variables

character *150 chain1,chain2
character *12   nomcoe(ngazem)

integer          it     , ice    , iat    , ios , ii , ioxy
integer          ncoel  , inicoe , inicha , ierror
integer          idecal , icla   , iclapc , icha   , is
integer          idebch , ifinch , lonch  , ichai  , ichcoe
integer          atcoel(ngazem,natom)

double precision fcor   , pcisec , pcipur , xashpc
double precision hco20  , ho20   , hh2o0
double precision den1   , den2
double precision tmin   , tmax
double precision ipci(ncharm)
double precision wmolce(ngazem), ehcoel(ngazem,npot)
double precision cpcoel(ngazem,npot),det,matdet
double precision wmv1,wmvch1,wmv2,wmvch2
double precision a11,a12,a13,a21,a22,a23,a31,a32,a33
double precision dhvol1,dhvol2,hcoke,hchar,ehvol1,ehvol2
double precision wmco,wmo2,wmco2,wmh2o,wmn2,wmc
double precision dmf4,dmf6,dmf7


!===============================================================================

!==================================================
! 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
!==================================================

! --> Ouverture du fichier

open ( unit=impfpp, file=ficfpp,                                  &
        STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',      &
                                        iostat=ios, err=99 )
rewind (unit=impfpp,err=99 )

! --> Lecture thermochimie

read (impfpp,*,err=999,end=999 )

! ---- Nb de constituants elementaires (gazeux et solide)

read ( impfpp,*,err=999,end=999 ) ncoel
if ( ncoel.gt.ngazem ) then
  write(nfecra,9991) ngazem,ncoel
  call csexit (1)
endif

! ---- Nb de points de tabulation ENTH-TEMP

read ( impfpp,*,err=999,end=999 ) npo
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

read (impfpp,*,err=999,end=999 )
read (impfpp,1010,err=999,end=999 ) chain1
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

 1010 format(a150)

! --- Temperature Min et Max

read (impfpp,*,err=999,end=999) tmin
read (impfpp,*,err=999,end=999) tmax


! ---- Nb especes atomiques (C, H, O, N, ...)

read (impfpp,*,err=999,end=999 ) nato
if ( nato.gt.natom ) then
  write(nfecra,9993) natom,nato
  call csexit (1)
  !==========
endif

! ---- Masse molaire especes atomiques
!      Composition des constituants elementaires en fonction
!        des especes elementaires

do iat = 1, nato
  read (impfpp,*,err=999,end=999 ) wmolat(iat),                   &
                      ( atcoel(ice,iat),ice=1,ncoel )
enddo

! ---- Calcul des masses molaires des constituants elementaires

do ice = 1, ncoel
  wmolce(ice) = 0.d0
  do iat = 1, nato
    wmolce(ice)= wmolce(ice) + atcoel(ice,iat)*wmolat(iat)
  enddo
enddo


! --> Lecture rayonnement

read (impfpp,*,err=999,end=999 )

! ---- Coefficient d'absorption du melange gazeux

read (impfpp,*,err=999,end=999 ) ckabs1


! --> Lecture caracteristiques charbon

read (impfpp,*,err=999,end=999 )

! ---- Nb de charbons

read ( impfpp,*,err=999,end=999 ) ncharb
if ( ncharb.gt.ncharm ) then
  write(nfecra,9995) ncharm,ncharb
  call csexit (1)
endif

! ---- Nb de classes par charbon

read (impfpp,*,err=999,end=999 )                                  &
                 ( nclpch(icha),icha=1,ncharb )
do icha = 1, ncharb
  if ( nclpch(icha).gt.nclcpm ) then
    write(nfecra,9996) nclcpm,nclpch(icha),icha
    call csexit (1)
  endif
enddo

! ---- Calcul du nb de classes et remplissage de ICHCOR

nclacp = 0
do icha = 1, ncharb
  nclacp = nclacp + nclpch(icha)
enddo
idecal = 0
do icha = 1, ncharb
  do iclapc = 1, nclpch(icha)
    icla = iclapc+idecal
    ichcor(icla) = icha
  enddo
  idecal = nclpch(icha)
enddo

! ---- Diametre initial par classe (m)

read (impfpp,*,err=999,end=999 )                                  &
     ( diam20(icla),icla=1,nclacp )

! ---- Composition elementaire en C, H, O sur sec (% en masse)

read (impfpp,*,err=999,end=999 )                                  &
     ( cch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( hch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( och(icha),icha=1,ncharb )

! ---- PCI sur charbon sec ou pur suivant la valeur de IPCI

read (impfpp,*,err=999,end=999 )                                  &
     ( ipci(icha),pcich(icha),icha=1,ncharb )

! ---- CP moyen du charbon sec (J/kg/K)

read (impfpp,*,err=999,end=999 )                                  &
     ( cp2ch(icha),icha=1,ncharb )

! ---- Masse volumique initiale (kg/m3)

read (impfpp,*,err=999,end=999 )                                  &
     ( rho0ch(icha),icha=1,ncharb )

! ---- Caracteristiques coke

read (impfpp,*,err=999,end=999 )

! ------- Composition elementaire en C, H, O sur sec (%)

read (impfpp,*,err=999,end=999 )                                  &
     ( cck(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( hck(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ock(icha),icha=1,ncharb )

! ------ PCI sur charbon sec

read (impfpp,*,err=999,end=999 )                                  &
     ( pcick(icha),icha=1,ncharb )

! ---- Caracteristiques cendres

read (impfpp,*,err=999,end=999 )

! ------ Taux de cendre (kg/kg) en %

read (impfpp,*,err=999,end=999 )                                  &
     ( xashch(icha),icha=1,ncharb )
!        Transformation en kg/kg
do icha = 1, ncharb
  xashch(icha) = xashch(icha)/100.d0
enddo

! ------ Enthalpie de formation des cendres (J/kg)

read (impfpp,*,err=999,end=999 )                                  &
     ( h0ashc(icha),icha=1,ncharb )

! ------ CP des cendres (J/kg/K)

read (impfpp,*,err=999,end=999 )                                  &
     ( cpashc(icha),icha=1,ncharb )

! ------ Taux d'humidite (kg/kg) en %

read (impfpp,*,err=999,end=999 )                                  &
     ( xwatch(icha),icha=1,ncharb )

!        Transformation en kg/kg
do icha = 1, ncharb
  xwatch(icha) = xwatch(icha)/100.d0
enddo

!      Transformation du taux de cendre de sec
!      sur humide en kg/kg

do icha = 1, ncharb
  xashch(icha) = xashch(icha)*(1.d0-xwatch(icha))
enddo

! ---- Parametres de devolatilisation (modele de Kobayashi)

read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( iy1ch(icha),y1ch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( iy2ch(icha),y2ch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( a1ch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( a2ch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( e1ch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( e2ch(icha),icha=1,ncharb )

! ---- Parametres combustion heterogene pour O2(modele a sphere retrecissante)

read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( ahetch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ehetch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( iochet(icha),icha=1,ncharb)

! ---- Parametres combustion heterogene pour CO2(modele a sphere retrecissante)

read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( ahetc2(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ehetc2(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ioetc2(icha),icha=1,ncharb)

! --> Lecture caracteristiques Oxydants

read (impfpp,*,err=999,end=999 )

! ---- Nb d'oxydants

read ( impfpp,*,err=999,end=999 ) noxyd
if ( noxyd.lt.1 .or. noxyd .gt. 3 ) then
  write(nfecra,9895) noxyd
  call csexit (1)
endif

! ---- Composition en O2,N2,H2O,N2

do ioxy=1,3
  oxyo2 (ioxy) = 0.d0
  oxyn2 (ioxy) = 0.d0
  oxyh2o(ioxy) = 0.d0
  oxyco2(ioxy) = 0.d0
enddo

read (impfpp,*,err=999,end=999 )                                  &
     ( oxyo2(ioxy),ioxy=1,noxyd )
read (impfpp,*,err=999,end=999 )                                  &
     ( oxyn2(ioxy),ioxy=1,noxyd )
read (impfpp,*,err=999,end=999 )                                  &
     ( oxyh2o(ioxy),ioxy=1,noxyd )
read (impfpp,*,err=999,end=999 )                                  &
     ( oxyco2(ioxy),ioxy=1,noxyd )

! --> Fermeture du fichier (ne pas oublier, car l'unite sert pour janaf)

close(impfpp)


!==============================================
! 2. CALCULS DE DONNEES COMPLEMENTAIRES
!==============================================


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

ngaze = 2 + 5 + 2*ncharb

! ---- Definition des pointeurs pour les tableaux WMOLE et EHGAZE
!      REMARQUE : Cette position de pointeurs va egalement servir
!                 pour le tableau de pointeurs IYM1 relatif aux
!                 tableaux PROPCE et PROPFB
!                 ON BALAYE JUSTE DE 1 A (NGAZE-2*NCHARB)

is    = 0
is    = is + 1
ichx1 = is
is    = is + 1
ichx2 = is
is   = is + 1
ico  = is
is   = is + 1
io2  = is
is   = is + 1
ico2 = is
is   = is + 1
ih2o = is
is   = is + 1
in2  = is
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
  ehgaze(ichx1,it) = ehcoel(1,it)
  ehgaze(ichx2,it) = ehcoel(2,it)
  ehgaze(ico  ,it) = ehcoel(3,it)
  ehgaze(io2  ,it) = ehcoel(4,it)
  ehgaze(ico2 ,it) = ehcoel(5,it)
  ehgaze(ih2o ,it) = ehcoel(6,it)
  ehgaze(in2  ,it) = ehcoel(7,it)
enddo
wmole(ichx1) = wmolce(1)
! ------ ATTENTION : Il faut prendre la masse molaire du monomere CH2
!                    et non celle du C2H4
wmole(ichx2) = wmolat(iatc)+wmolat(iath)*2.d0
wmole(ico  ) = wmolce(3)
wmole(io2  ) = wmolce(4)
wmole(ico2 ) = wmolce(5)
wmole(ih2o ) = wmolce(6)
wmole(in2  ) = wmolce(7)


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

alpham = zero
betam  = zero

do icha = 1, ncharb

! ------ Composition du charbon reactif
!        sous la forme CH(ALPHA)O(BETA)
!        et calcul de la masse molaire

  alpha(icha) = (hch(icha)/wmolat(iath))                          &
              / (cch(icha)/wmolat(iatc))
  beta(icha)  = (och(icha)/wmolat(iato))                          &
              / (cch(icha)/wmolat(iatc))
  alpham      = alpham + alpha(icha)
  betam       = betam  + beta(icha)
  wmols(ich(icha)) = wmolat(iatc) + alpha(icha)*wmolat(iath)      &
                   + beta(icha)*wmolat(iato)
enddo
alpham = alpham / dble(ncharb)
betam  = betam  / dble(ncharb)

! ------ Transformation par la formule de Schaff du
!        PCI sur sec en PCI sur pur LORSQUE IPCI = 1

fcor = 1.2
do icha = 1, ncharb

  if (ipci(icha).eq.1) then
!    On donne directement le PCI sur sec J/kg ---> kcal/kg
    pcisec = pcich(icha)/1000.d0/xcal2j
    xashpc = xashch(icha)*100.d0
    pcipur = ( pcisec+5.95d0*(zero+(1.d0-fcor)*xashpc) )*100.d0   &
           / ( 100.d0-zero-fcor*xashpc )
!    Conversion de PCI sur pur de kcal/kg ----> J/kg
    pcipur = pcipur*1000.d0*xcal2j
!    On ecrase PCICH sur sec par PCIPUR
    pcich(icha) = pcipur
  endif

enddo

! ------ Calcul de : H02 pour CH , CK , ASH et WT

do icha = 1, ncharb
  hco20 = ehgaze(ico2,1)*wmole(ico2)
  hh2o0 = ehgaze(ih2o,1)*wmole(ih2o)
  ho20  = ehgaze(io2 ,1)*wmole(io2 )

  h02ch(icha) = pcich(icha) +                                     &
     ( hco20 + alpha(icha)/2.d0*hh2o0 -                           &
       (1.d0+alpha(icha)/4.d0-beta(icha)/2.d0)*ho20 )             &
     / wmols(ich(icha))

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
      ehsoli(ich(icha),it) = eh0sol(ich(icha))                    &
                            + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      hco20 = ehgaze(ico2,it)*wmole(ico2)
      hh2o0 = ehgaze(ih2o,it)*wmole(ih2o)
      ho20  = ehgaze(io2,it) *wmole(io2)
      ehsoli(ich(icha),it) = pcich(icha) +                        &
       ( hco20 + alpha(icha)/2.d0*hh2o0 -                         &
         (1.d0+alpha(icha)/4.d0-beta(icha)/2.d0)*ho20 )           &
       / wmols(ich(icha))
    enddo
  endif
enddo

! ---- Calcul relatif au coke

! ------ Par defaut Coke = Carbone solide
!                   GAMMA = 0 et DELTA = 0
!            Si CP2CH > 0 : HCK = H02CH + CP2CH(T2-TREFTH)
!            Sinon       : HCK = Enthalpie du carbone pur

gamma(icha) = zero
delta(icha) = zero
do icha = 1, ncharb
  wmols(ick(icha)) = wmolce(8)
  if (cp2ch(icha).gt.epsicp) then
    do it = 1, npoc
      ehsoli(ick(icha),it) = eh0sol(ick(icha))                    &
                           + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      ehsoli(ick(icha),it) = ehcoel(8,it)
    enddo
  endif
enddo

! ------ Coke = CH(GAMMA)O(DELTA)

do icha = 1, ncharb
  if ( pcick(icha).gt.epsicp ) then
    gamma(icha) = hck(icha)/cck(icha)
    delta(icha) = ock(icha)/cck(icha)
    wmols(ick(icha)) = wmolat(iatc)+gamma(icha)*wmolat(iath)      &
                     + delta(icha)*wmolat(iato)
!        On considere le PCI constant qqs T
    do it = 1, npoc
      hco20 = ehgaze(ico2,it)*wmole(ico2)
      hh2o0 = ehgaze(ih2o,it)*wmole(ih2o)
      ho20  = ehgaze(io2,it) *wmole(io2)
      ehsoli(ick(icha),it) = pcick(icha) +                        &
       ( hco20 + gamma(icha)/2.d0*hh2o0 -                         &
        (1.d0+gamma(icha)/4.d0-delta(icha)/2.d0)*ho20 )           &
       / wmols(ich(icha))
    enddo
  endif
enddo

! ---- Calcul relatif aux cendres

do icha = 1, ncharb
  if (cp2ch(icha).gt.epsicp) then
    do it = 1, npoc
      ehsoli(iash(icha),it) = eh0sol(iash(icha))                  &
                            + cp2ch(icha)*(thc(it)-trefth)
    enddo
  else
    do it = 1, npoc
      ehsoli(iash(icha),it) = h0ashc(icha)                        &
                            + cpashc(icha)*(thc(it)-trefth)
    enddo
  endif
  wmols(iash(icha)) = zero
enddo

! ---- Calcul relatif a l'eau

do icha = 1, ncharb
  cp2wat(icha) = 4181.35d0

  do it = 1, npoc
    ehsoli(iwat(icha),it) = eh0sol(iwat(icha))                    &
                           +cp2wat(icha)*(thc(it)-trefth)
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
  if ( iy1ch(icha).lt. 0 .or.  iy1ch(icha).le.0 .or.              &
       iy2ch(icha).lt. 0 .or.  iy2ch(icha).le.0     ) then
    write(nfecra,9981) icha,iy1ch(icha),iy1ch(icha)
    call csexit(1)
  endif
enddo

! ---- Matieres volatiles legeres : [CH(CHX1); CO]
!       CH(ALPHA)O(BETA) --> A1 CH(CHX1) + B1 CO
!                            + (1-A1-B1) CH(GAMMA)O(DELTA)
!       Si IY1CH = 0 CHX1 fixe , Y1CH calcule
!       Si IY1CH = 1 Y1CH fixe , CHX1 calcule
!       Si IY1CH = 2 Y1CH fixe, CHX1 fixe, on ajoute de l'eau
!       CH(ALPHA)O(BETA) --> A1 CH(CHX1) + B1 CO + C1 H2O
!                            + (1-A1-B1) CH(GAMMA)O(DELTA)

do icha = 1, ncharb
  if ( iy1ch(icha).eq.0 ) then
    chx1(icha)   = 4.d0
    den1 = ( (1.d0-delta(icha))*chx1(icha)-gamma(icha) )          &
        * wmols(ich(icha))
    y1ch(icha) = ( ( gamma(icha)*(beta(icha)-1.d0) +              &
                     alpha(icha)*(1.d0-delta(icha)) )             &
                     *(wmolat(iatc)+chx1(icha)*wmolat(iath))      &
                 + ( beta(icha)*(chx1(icha)-gamma(icha)) +        &
                     delta(icha)*(alpha(icha)-chx1(icha)) )       &
                     *wmole(ico) )                                &
               / den1
  elseif(iy1ch(icha) .eq. 1) then
    den1 = y1ch(icha)*(1.d0-delta(icha))*wmols(ich(icha))         &
         - wmolat(iath)*(gamma(icha)*(beta(icha)-1.d0) +          &
                         alpha(icha)*(1.d0-delta(icha)))          &
         + wmole(ico)*(delta(icha)-beta(icha))
    chx1(icha) = ( wmolat(iatc)*(alpha(icha)-gamma(icha))         &
                 + wmolat(iato)*(delta(icha)*alpha(icha)-         &
                                 beta(icha) *gamma(icha)  )       &
                 + y1ch(icha)*gamma(icha)*wmols(ich(icha)) )      &
               / den1
  elseif(iy1ch(icha) .eq. 2) then
    chx1(icha)   = 4.d0
  endif

  if ( iy1ch(icha).lt.2) then
   den1     =  (1.d0-delta(icha))*chx1(icha)-gamma(icha)
   a1(icha) = ( gamma(icha)*(beta(icha)-1.d0) +                   &
                alpha(icha)*(1.d0-delta(icha)) )       / den1
   b1(icha) = ( beta(icha)*(chx1(icha)-gamma(icha)) +             &
                delta(icha)*(alpha(icha)-chx1(icha)) ) / den1
   c1(icha) = zero
  else
   wmv1 = y1ch(icha)*( wmolat(iatc)                               &
                      +alpha(icha)*wmolat(iath)                   &
                      +beta(icha)*wmolat(iato) )
   wmvch1 = wmolat(iatc) + chx1(icha)*wmolat(iath)

   a11 = 4.d0-gamma(icha)
   a12 = -gamma(icha)
   a13 = 2.d0
   a21 = -delta(icha)
   a22 = 1.d0-delta(icha)
   a23 = 1.d0
   a31 = wmvch1
   a32 = wmole(ico)
   a33 = wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)

   den1 = 1.d0 / det

   a11 =  (alpha(icha)-gamma(icha))
   a12 = -gamma(icha)
   a13 =  2.d0
   a21 =  (beta(icha)-delta(icha) )
   a22 =  (1.d0-delta(icha))
   a23 =  1.d0
   a31 =  wmv1
   a32 =  wmole(ico)
   a33 =  wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   a1(icha) = den1 * det

   a11 = 4.d0-gamma(icha)
   a12 = alpha(icha)-gamma(icha)
   a13 = 2.d0
   a21 = -delta(icha)
   a22 = beta(icha)-delta(icha)
   a23 = 1.d0
   a31 = wmvch1
   a32 = wmv1
   a33 = wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   b1(icha) = den1 * det

   a11 = 4.d0-gamma(icha)
   a12 = -gamma(icha)
   a13 = alpha(icha)-gamma(icha)
   a21 = -delta(icha)
   a22 = 1.d0-delta(icha)
   a23 = beta(icha)-delta(icha)
   a31 = wmvch1
   a32 = wmole(ico)
   a33 = wmv1
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   c1(icha) = den1 * det
  endif
 enddo

! ---- Matieres volatiles lourdes : [CH(CHX2); CO]
!       CH(ALPHA)O(BETA) --> A2 CH(CHX2) + B2 CO
!                            + (1-A2-B2) CH(GAMMA)O(DELTA)
!       Si IY2CH = 0 CHX2 fixe, Y2CH calcule
!       Si IY2CH = 1 Y2CH fixe, CHX2 calcule
!       Si IY2CH = 2 Y2CH fixe, CHX2 fixe, on ajoute de l'eau
!       CH(ALPHA)O(BETA) --> A2 CH(CHX2) + B2 CO + C2 H2O
!                            + (1-A2-B2) CH(GAMMA)O(DELTA)

do icha = 1, ncharb
  if ( iy2ch(icha).eq.0 ) then
    chx2(icha)   = 2.d0
    den2 = ( (1.d0-delta(icha))*chx2(icha)-gamma(icha) )          &
        * wmols(ich(icha))
    y2ch(icha) = ( ( gamma(icha)*(beta(icha)-1.d0) +              &
                     alpha(icha)*(1.d0-delta(icha)) )             &
                     *(wmolat(iatc)+chx2(icha)*wmolat(iath))      &
                 + ( beta(icha)*(chx2(icha)-gamma(icha)) +        &
                     delta(icha)*(alpha(icha)-chx2(icha)) )       &
                     *wmole(ico) )                                &
               / den2
  elseif(iy2ch(icha).eq.1) then
    den2 = y2ch(icha)*(1.d0-delta(icha))*wmols(ich(icha))         &
         - wmolat(iath)*(gamma(icha)*(beta(icha)-1.d0) +          &
                         alpha(icha)*(1.d0-delta(icha)))          &
         + wmole(ico)*(delta(icha)-beta(icha))
    chx2(icha) = ( wmolat(iatc)*(alpha(icha)-gamma(icha))         &
                 + wmolat(iato)*(delta(icha)*alpha(icha)-         &
                                 beta(icha) *gamma(icha)  )       &
                 + y2ch(icha)*gamma(icha)*wmols(ich(icha)) )      &
               / den2
  else
    y2ch(icha) = min(2.d0*y1ch(icha),(1.d0+y1ch(icha))/2.d0)
    chx2(icha) = 2.d0
  endif

  if( iy2ch(icha).lt.2) then
   den2     = (1.d0-delta(icha))*chx2(icha)-gamma(icha)
   a2(icha) = ( gamma(icha)*(beta(icha)-1.d0) +                   &
                alpha(icha)*(1.d0-delta(icha)) )       / den2
   b2(icha) = ( beta(icha)*(chx2(icha)-gamma(icha)) +             &
                delta(icha)*(alpha(icha)-chx2(icha)) ) / den2
   c2(icha) = zero
  else
   wmv2 = y2ch(icha)*( wmolat(iatc)                               &
                      +alpha(icha)*wmolat(iath)                   &
                      +beta(icha)*wmolat(iato) )
   wmvch2 = wmolat(iatc) + chx2(icha)*wmolat(iath)

   a11 = 2.d0-gamma(icha)
   a12 = - gamma(icha)
   a13 = 2.d0
   a21 = -delta(icha)
   a22 = 1.d0-delta(icha)
   a23 = 1.d0
   a31 = wmvch2
   a32 = wmole(ico)
   a33 = wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   den2 = 1.d0 / det

   a11 = (alpha(icha)-gamma(icha))
   a12 = - gamma(icha)
   a13 = 2.d0
   a21 = (beta(icha)-delta(icha) )
   a22 = (1.d0-delta(icha))
   a23 = 1.d0
   a31 = wmv2
   a32 = wmole(ico)
   a33 = wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   a2(icha) = den2 * det

   a11 = (2.d0-gamma(icha))
   a12 = (alpha(icha)-gamma(icha))
   a13 = 2.d0
   a21 = -delta(icha)
   a22 = (beta(icha)-delta(icha))
   a23 = 1.d0
   a31 = wmvch2
   a32 = wmv2
   a33 = wmole(ih2o)
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   b2(icha) = den2 * det

   a11 = 2.d0-gamma(icha)
   a12 = - gamma(icha)
   a13 = alpha(icha)-gamma(icha)
   a21 = -delta(icha)
   a22 = 1.d0-delta(icha)
   a23 = (beta(icha)-delta(icha))
   a31 = wmvch2
   a32 = wmole(ico)
   a33 = wmv2
   det = matdet(a11, a12, a13, a21, a22, a23, a31, a32, a33)
   c2(icha) = den2 * det
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


enddo

! AFFICHAGE RECAPITULATIF
! =======================

write(nfecra,8000)

!  Constante du modele de devolatilisation

write(nfecra,8100)
do icha=1,ncharb
  write(nfecra,8101) icha
  if ( iy1ch(icha) .eq.0 ) then
    write(nfecra,8102)
  else if ( iy1ch(icha) .eq.1 ) then
    write(nfecra,8103)
  else if  ( iy1ch(icha) .eq.2 ) then
    write(nfecra,8104)
  endif

  write(nfecra,8105) y1ch(icha),y2ch(icha)
enddo


!  Composition matieres volatiles

write(nfecra,8010)

write(nfecra,8051)
ierror = 0
do icha=1,ncharb
  write(nfecra,8052) icha,chx1(icha),chx2(icha)
  if ( chx1(icha) .le. 0.d0 .or. chx2(icha) .le. 0 ) then
    ierror = 1
  endif
enddo
if ( ierror .eq. 1 ) then
  write(nfecra,9970)
  call csexit(1)
endif

do icha=1,ncharb
  write(nfecra,8011)  icha
  write(nfecra,8012)
  write(nfecra,8013) a1(icha),b1(icha),c1(icha)
  write(nfecra,8014) a2(icha),b2(icha),c2(icha)
enddo

!  Calcul de DHdev

write(nfecra,8001)
do icha=1,ncharb

!    Enthalie du Charbon

 hchar = ehsoli(ich(icha),1)
 hcoke = ehsoli(ick(icha),1)

!    Enthalpie des matieres volatiles

 ehvol1 = ( a1(icha)*ehgaze(ichx1c(icha),1)*wmole(ichx1c(icha))   &
           +b1(icha)*ehgaze(ico ,1)*wmole(ico)                    &
           +c1(icha)*ehgaze(ih2o,1)*wmole(ih2o) )                 &
         /( a1(icha)*wmole(ichx1c(icha))                          &
           +b1(icha)*wmole(ico)                                   &
           +c1(icha)*wmole(ih2o) )

 ehvol2 = ( a2(icha)*ehgaze(ichx2c(icha),1)*wmole(ichx2c(icha))   &
           +b2(icha)*ehgaze(ico ,1)*wmole(ico)                    &
           +c2(icha)*ehgaze(ih2o,1)*wmole(ih2o) )                 &
         /( a2(icha)*wmole(ichx2c(icha))                          &
           +b2(icha)*wmole(ico)                                   &
           +c2(icha)*wmole(ih2o) )

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

acof3 =  1.d0/wmc
ao2f3 = -.5d0/wmc

dmf4  = ( oxyo2 (1)*wmo2 +oxyn2 (1)*wmn2                          &
         +oxyh2o(1)*wmh2o+oxyco2(1)*wmco2 )

if ( dmf4 .le. 0.d0 ) then
  write(nfecra,9896) oxyo2(1) ,oxyn2(1) ,                         &
                     oxyh2o(1),oxyco2(1)
  call csexit(1)
endif

ao2f4  = oxyo2(1)  / dmf4
an2f4  = oxyn2(1)  / dmf4
ah2of4 = oxyh2o(1) / dmf4
aco2f4 = oxyco2(1) / dmf4

if ( noxyd .ge. 2.d0 ) then
  dmf6  = ( oxyo2 (2)*wmo2 +oxyn2 (2)*wmn2                        &
           +oxyh2o(2)*wmh2o+oxyco2(2)*wmco2 )
  if ( dmf6 .le. 0.d0 ) then
    write(nfecra,9897) oxyo2(2) ,oxyn2(2) ,                       &
                       oxyh2o(2),oxyco2(2)
    call csexit(1)
  endif

  ao2f6  = oxyo2(2)  / dmf6
  an2f6  = oxyn2(2)  / dmf6
  ah2of6 = oxyh2o(2) / dmf6
  aco2f6 = oxyco2(2) / dmf6

else
  ao2f6  = 0.d0
  an2f6  = 0.d0
  ah2of6 = 0.d0
  aco2f6 = 0.d0
endif

if ( noxyd .eq. 3.d0 ) then
  dmf7  = ( oxyo2 (3)*wmo2 +oxyn2 (3)*wmn2                        &
           +oxyh2o(3)*wmh2o+oxyco2(3)*wmco2 )
  if ( dmf7 .le. 0.d0 ) then
    write(nfecra,9898) oxyo2(3) ,oxyn2(3) ,                       &
                       oxyh2o(3),oxyco2(3)
    call csexit(1)
  endif

  ao2f7  = oxyo2(3)  / dmf7
  an2f7  = oxyn2(3)  / dmf7
  ah2of7 = oxyh2o(3) / dmf7
  aco2f7 = oxyco2(3) / dmf7

else
  ao2f7  = 0.d0
  an2f7  = 0.d0
  ah2of7 = 0.d0
  aco2f7 = 0.d0
endif

ah2of5  = 1.d0/wmh2o

return


!============================
! 3. SORTIE EN ERREUR
!============================

  99  continue
write ( nfecra,9998 )
call csexit (1)
!==========

  999 continue
write ( nfecra,9999 )
call csexit (1)
!==========


!--------
! FORMATS
!--------


 8000 format(1X,' RECAPITULATIF SUR LES CHARBONS :'/,             &
       1X,' ==============================  '  )
 8001 format(/,3X,'Delta Hdev  ',/,                               &
       5X,' Numero Charbon',6X,'Legeres',7X,'Lourdes',9X,'PCI')
 8002 format(10x,i3,8x,g15.7,1x,g15.7,g15.7)

 8010 format(/,3X,'COMPOSITION MATIERES VOLATILES ')
 8011 format(/,5X,' Numero Charbon',I6)
 8012 format(18X,'CHX1 ',12X,'CO',12X,'H2O')
 8013 format( 8X,' MV1 ',G15.7,G15.7,G15.7)
 8014 format( 8X,' MV2 ',G15.7,G15.7,G15.7)
 8020 format(/,3X,'LOI ENTHALPIE/TEMERATURE ')
 8021 format( 10X,'TEMPERATURE',9X,'Hch ',10X,'Hcoke',10X,              &
                             'Hash',10X,'Hwat')
 8022 format( 8x,5(g15.7))

 8050 format(/,3X,'COMPOSITION DU CHARBON ')
 8051 format(/,5X,'Numero Charbon',7X,'CHX1',12X,'CHX2')
 8052 format(10x,i3,8x,g15.7,1x,g15.7)

 8100 format(/,3X,'CONSTANTE DU MODELE DE DEVOLATILISATION')
 8101 format(/,8X,'CHARBON ',I2)
 8102 format(/,12X,'OPTION 0 : Y1 et Y2 CALCULES     ')
 8103 format(/,12X,'OPTION 1 : Y1 et Y2 FIXES        ')
 8104 format(/,12X,'OPTION 2 : CHX1 et CHX2 FIXES    ')

 8105 format(12X,'Y1CH = ',G15.7,3X,'Y2CH = ',G15.7)

 9970 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  CERTAINES VALEURS DE CHX1 ET CHX2 SONT NEGATIVES          ',/,&
'@  OU NULLES                                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9980 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Vous utilisez une option differente pour Y1 et Y2         ',/,&
'@  pour le charbon ',I2,'                                    ',/,&
'@         IY1CH = ',I2,'                                     ',/,&
'@         IY2CH = ',I2,'                                     ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9981 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Vous utilisez une option  pour Y1 et Y2 non disponible    ',/,&
'@  pour le charbon ',I2,'                                    ',/,&
'@         IY1CH = ',I2,'                                     ',/,&
'@         IY2CH = ',I2,'                                     ',/,&
'@                                                            ',/,&
'@  seules les valeur : 0 , 1 ou 2 sont admissibles.          ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9991 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes courantes doit etre inferieur        ',/,&
'@                                  ou egal a',I10             ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9992 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre de points de tabulation est limite a ',I10       ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9993 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes elementaires est limite a ',I10       ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9995 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre de charbons est limite a      ',I10              ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9996 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre de classes par charbon est limite a ',I10        ,/,&
'@   Il vaut ',I10   ,' pour le charbon ',I10                  ,/,&
'@                      dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Erreur a l''ouverture du fichier parametrique.            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Erreur a la lecture du fichier parametrique.              ',/,&
'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,&
'@    ou son format inadapte.                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9895 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Le nombre d''Oxydants doit etre compris entre 1 et 3      ',/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9896 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  LA COMPOSITION DE L''OXYDANT 1 EST ERRONEE                ',/,&
'@     O2  :  ',G15.7,'                                       ',/,&
'@     N2  :  ',G15.7,'                                       ',/,&
'@     H2O :  ',G15.7,'                                       ',/,&
'@     CO2 :  ',G15.7,'                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9897 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  LA COMPOSITION DE L''OXYDANT 2 EST ERRONEE                ',/,&
'@     O2  :  ',G15.7,'                                       ',/,&
'@     N2  :  ',G15.7,'                                       ',/,&
'@     H2O :  ',G15.7,'                                       ',/,&
'@     CO2 :  ',G15.7,'                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9898 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  LA COMPOSITION DE L''OXYDANT 3 EST ERRONEE                ',/,&
'@     O2  :  ',G15.7,'                                       ',/,&
'@     N2  :  ',G15.7,'                                       ',/,&
'@     H2O :  ',G15.7,'                                       ',/,&
'@     CO2 :  ',G15.7,'                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9899 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  LA COMPOSITION DE L''OXYDANT 1 EST ERRONEE                ',/,&
'@     O2  :  ',G15.7,'                                       ',/,&
'@     N2  :  ',G15.7,'                                       ',/,&
'@     H2O :  ',G15.7,'                                       ',/,&
'@     CO2 :  ',G15.7,'                                       ',/,&
'@                                                            ',/,&
'@  ALORS QUE L''OXYDANT 1 DOIT OBLIGATOIREMENT CONTENIR      ',/,&
'@  DE L''OXYGENE                                             ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


end subroutine


function matdet &
!==============

 ( a11 , a12 , a13 , a21 , a22 , a23 , a31 , a32 , a33 )

!===============================================================================
! FONCTION :
! --------

! CALCUL DU DETERMINANT DE LA MATRICE 3x3

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! aa1 a a33        ! r  ! <-- ! coefficient de la matrice                      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

!===============================================================================

implicit none

! Arguments

double precision a11, a12, a13, a21, a22, a23, a31, a32, a33
double precision matdet

!===============================================================================

!===============================================================================
! 1. CALCUL DU DETERMINANT
!===============================================================================

matdet = a11*a22*a33 + a21*a32*a13 + a31*a12*a23                  &
       - a11*a32*a23 - a21*a12*a33 - a31*a22*a13

end function
