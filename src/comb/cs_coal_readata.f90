!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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
!=========================
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

!===============================================================================

implicit none

! Arguments

! Local variables

character *150  chain1,chain2
character *12   nomcoe(ngazem)

integer          it     , ice    , iat    , ios , ii ,jj, ioxy
integer          ncoel  , inicoe , inicha , ierror
integer          idecal , icla   , iclapc , icha   , is
integer          idebch , ifinch , lonch  , ichai  , ichcoe
integer          atcoel(ngazem,natom)
integer          ierr   , ndim

double precision fcor   , pcisec , pcipur , xashpc
double precision hco20  , ho20   , hh2o0 , hso20, hnh30
double precision den1   , den2
double precision tmin   , tmax
double precision ipci(ncharm)
double precision wmolce(ngazem), ehcoel(ngazem,npot)
double precision cpcoel(ngazem,npot),det,matdet
double precision wmv1,wmvch1,wmv2,wmvch2
double precision a11,a12,a13,a21,a22,a23,a31,a32,a33
double precision dhvol1,dhvol2,hcoke,hchar,ehvol1,ehvol2
double precision wmco,wmo2,wmco2,wmh2o,wmn2,wmc
double precision dmf3,dmf4,dmf5,som1,som2
double precision sm(8),solu(8),mat(8,8)

!===============================================================================

!===============================================================================
! 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
!===============================================================================

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

! ---- Nb especes atomiques (C, H, O, S et N)

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

! ---- Composition elementaire en C, H , O , N , S
!                 sur sec (% en masse)

read (impfpp,*,err=999,end=999 )                                  &
     ( cch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( hch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( och(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( nch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( sch(icha),icha=1,ncharb )

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

! ------- Composition elementaire en C , H , O , N , S sur sec (%)

read (impfpp,*,err=999,end=999 )                                  &
     ( cck(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( hck(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ock(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( nck(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( sck(icha),icha=1,ncharb )

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
! ---- Repartition de l'azote entre HCN et NH3
read (impfpp,*,err=999,end=999 )                                  &
     ( crepn1(1,icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( crepn1(2,icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( crepn2(1,icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( crepn2(2,icha),icha=1,ncharb )
! ---- Parametres combustion heterogene pour O2
!      (modele a sphere retrecissante)
read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( ahetch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ehetch(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( iochet(icha),icha=1,ncharb)
! ---- Parametres combustion heterogene pour CO2
!      (modele a sphere retrecissante)
read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( ahetc2(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ehetc2(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ioetc2(icha),icha=1,ncharb)

! ---- Parametres combustion heterogene pour H2O
!      (modele a sphere retrecissante)

read (impfpp,*,err=999,end=999 )

read (impfpp,*,err=999,end=999 )                                  &
     ( ahetwt(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ehetwt(icha),icha=1,ncharb )
read (impfpp,*,err=999,end=999 )                                  &
     ( ioetwt(icha),icha=1,ncharb)

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
!                 pour le tableau de pointeurs IYM1 relatif aux
!                 tableaux PROPCE et PROPFB
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
!        sous la forme CH(ALPHA)O(BETA)S(TETA)N(OMEGA)
!        et calcul de la masse molaire

  alpha(icha) = (hch(icha)/wmolat(iath))           &
              / (cch(icha)/wmolat(iatc))
  beta(icha)  = (och(icha)/wmolat(iato))           &
              / (cch(icha)/wmolat(iatc))
  teta(icha)  = (sch(icha)/wmolat(iats))           &
              / (cch(icha)/wmolat(iatc))
  omega(icha) = (nch(icha)/wmolat(iatn))           &
              / (cch(icha)/wmolat(iatc))

  wmols(ich(icha)) = wmolat(iatc)                  &
                   + alpha(icha)*wmolat(iath)      &
                   + beta (icha)*wmolat(iato)      &
                   + teta (icha)*wmolat(iats)      &
                   + omega(icha)*wmolat(iatn)
enddo

!
! Transformation des coef de repartition de l'azote en HCN etNH3
!
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
  hso20 = ehgaze(iso2,1)*wmole(iso2)

  h02ch(icha) = pcich(icha) +                                  &
     (  hco20 + alpha(icha)/2.d0*hh2o0                         &
              + teta (icha)     *hso20                         &
      -( 1.d0 + alpha(icha)/4.d0                               &
              - beta (icha)/2.d0                               &
              + teta (icha)      )*ho20  ) / wmols(ich(icha))

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
      ehsoli(ich(icha),it) = pcich(icha)                        &
       +( hco20 + alpha(icha)/2.d0*hh2o0                        &
                + teta (icha)     *hso20                        &
        -( 1.d0 + alpha(icha)/4.d0                              &
                - beta (icha)/2.d0                              &
                + teta (icha)      ) *ho20 )  / wmols(ich(icha))
    enddo
  endif
enddo

! ---- Calcul relatif au coke

! ------ Par defaut Coke = Carbone solide
!                   GAMMA = 0 , DELTA = 0 , kappa = 0
!            Si CP2CH > 0 : HCK = H02CH + CP2CH(T2-TREFTH)
!            Sinon       : HCK = Enthalpie du carbone pur

gamma(icha) = zero
delta(icha) = zero
kappa(icha) = zero
zeta (icha) = zero
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

! ------ Coke = CH(GAMMA)O(DELTA)S(KAPPA)N(ZETA)

do icha = 1, ncharb
  if ( pcick(icha).gt.epsicp ) then
    gamma(icha) = hck(icha)/cck(icha)
    delta(icha) = ock(icha)/cck(icha)
    kappa(icha) = sck(icha)/cck(icha)
    zeta (icha) = nck(icha)/cck(icha)
!
    wmols(ick(icha)) = wmolat(iatc)                            &
                      + gamma(icha)*wmolat(iath)               &
                      + delta(icha)*wmolat(iato)               &
                      + kappa(icha)*wmolat(iats)               &
                      + zeta (icha)*wmolat(iatn)

!        On considere le PCI constant qqs T
    do it = 1, npoc
      hco20 = ehgaze(ico2,it)*wmole(ico2)
      hh2o0 = ehgaze(ih2o,it)*wmole(ih2o)
      hso20 = ehgaze(iso2,it)*wmole(iso2)
      ho20  = ehgaze(io2 ,it)*wmole(io2 )
      ehsoli(ick(icha),it) = pcick(icha)                       &
       + ( hco20  + gamma(icha)/2.d0*hh2o0                     &
                  + kappa(icha)     *hso20                     &
         - ( 1.d0 + gamma(icha)/4.d0                           &
                    - delta(icha)/2.d0                         &
                    + kappa(icha)      )*ho20  )  / wmols(ick(icha))
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
  if ( iy1ch(icha).lt. 0 .or.  iy1ch(icha).gt.2 .or.            &
       iy2ch(icha).lt. 0 .or.  iy2ch(icha).gt.2     ) then
    write(nfecra,9981) icha,iy1ch(icha),iy1ch(icha)
    call csexit(1)
  endif
enddo

! ---- Matieres volatiles legeres : [CH(CHX1); CO]
!       CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A1 CH(CHX1) + B1 CO + D1 H2S
!                                  + E1 HCN + F1 NH3
!                                  + (1-A1-B1-E1) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!       Si IY1CH = 0 CHX1 fixe , Y1CH calcule
!       Si IY1CH = 1 Y1CH fixe , CHX1 calcule
!       Si IY1CH = 2 Y1CH fixe , CHX1 fixe, on ajoute de l'eau
!       CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A1 CH(CHX1) + B1 CO + C1 H2O + D1 H2S
!                                  + E1 HCN + F1 NH3
!                                  + (1-A1-B1-E1) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!
do icha = 1, ncharb
!
!matrice 8x8
!
!Ligne 1
  mat(1,1) = -gamma(icha)
  mat(1,2) = -gamma(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.D0
  mat(1,5) = 1.D0-gamma(icha)
  mat(1,6) = 3.D0
  mat(1,7) = 1.D0
  mat(1,8) = 0.D0
!Ligne 2
  mat(2,1) = -delta(icha)
  mat(2,2) = 1.d0-delta(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
  mat(2,5) = -delta(icha)
  mat(2,6) = 0.d0
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0
!
  mat(3,1) = -kappa(icha)
  mat(3,2) = -kappa(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-kappa(icha)
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
  mat(5,1) = -zeta(icha)
  mat(5,2) = -zeta(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
  mat(5,5) = -zeta(icha)
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
  mat(6,8) = -(wmolat(iatc)+alpha(icha)*wmolat(iath)    &
                           +beta (icha)*wmolat(iato)    &
                           +teta (icha)*wmolat(iats)    &
                           +omega(icha)*wmolat(iatn))
!
!Second membre
!
  sm(1) = alpha(icha)-gamma(icha)
  sm(2) = beta (icha)-delta(icha)
  sm(3) = teta (icha)-kappa(icha)
  sm(4) = 0.d0
  sm(5) = omega(icha)-zeta(icha)
  sm(6) = 0.d0
!
! On complete la matrice et le second menbre en fonction de
! la valeur de IY1CH
!
  if ( iy1ch(icha).eq.0 ) then
! Valeur de X1
    chx1(icha) = 4.D0
!Matrice
    mat(7,1) = -chx1(icha)
    mat(7,2) = 0.D0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
!Second membre
    sm(7) = 0.D0
    sm(8) = 0.D0
!
  else if ( iy1ch(icha).eq.1 ) then
!
!Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
!Second membre
    sm(7) = y1ch(icha)
    sm(8) = 0.D0
!
  else if ( iy1ch(icha).eq.2 ) then
! Valeur de X1
    chx1(icha) = 4.D0
!
!Matrice
    mat(7,1) = -chx1(icha)
    mat(7,2) = 0.D0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
!Second membre
    sm(7) = 0.D0
    sm(8) = y1ch(icha)
!
  endif
!
  ndim = 8
  call coal_resol_matrice( ndim, mat , sm , solu , ierr)
!
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
!
! ---- Matieres volatiles lourdes : [CH(CHX2); CO]
!       CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A2 CH(CHX1) + B2 CO + D2 H2S
!                                  + E2 HCN + F2 NH3
!                                  + (1-A2-B2-E2) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!       Si IY2CH = 0 CHX2 fixe, Y2CH calcule
!       Si IY2CH = 1 Y2CH fixe, CHX2 calcule
!       Si IY2CH = 2 Y2CH fixe, CHX2 fixe, on ajoute de l'eau
!       CH(ALPHA)O(BETA)N(OMEGA)S(TETA) --> A2 CH(CHX1) + B2 CO + C2 H2O + D2 H2S
!                                  + E2 HCN + F2 NH3
!                                  + (1-A2-B2-E2) CH(GAMMA)O(DELTA)N(ZETA)S(KAPPA)
!
do icha = 1, ncharb
  mat(1,1) = -gamma(icha)
  mat(1,2) = -gamma(icha)
  mat(1,3) = 2.d0
  mat(1,4) = 2.D0
  mat(1,5) = 1.D0
  mat(1,6) = 0.D0
  mat(1,7) = 1.D0
  mat(1,8) = 0.D0
!Ligne 2
  mat(2,1) = -delta(icha)
  mat(2,2) = 1.d0-delta(icha)
  mat(2,3) = 1.d0
  mat(2,4) = 0.d0
  mat(2,5) = 0.d0
  mat(2,6) = 1.d0
  mat(2,7) = 0.d0
  mat(2,8) = 0.d0
!
  mat(3,1) = -kappa(icha)
  mat(3,2) = -kappa(icha)
  mat(3,3) = 0.d0
  mat(3,4) = 0.d0
  mat(3,5) = 1.d0-kappa(icha)
  mat(3,6) = 1.d0
  mat(3,7) = 0.d0
  mat(3,8) = 0.d0
!
  mat(4,1) = 0.d0
  mat(4,2) = 0.d0
  mat(4,3) = 0.d0
  mat(4,4) = 0.d0
  mat(4,5) = crepn2(2,icha)
  mat(4,6) =-crepn2(1,icha)
  mat(4,7) = 0.d0
  mat(4,8) = 0.d0
!
  mat(5,1) = -zeta(icha)
  mat(5,2) = -zeta(icha)
  mat(5,3) = 0.d0
  mat(5,4) = 1.d0
  mat(5,5) = 0.d0
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
  mat(6,8) = -(wmolat(iatc)+alpha(icha)*wmolat(iath)    &
                           +beta (icha)*wmolat(iato)    &
                           +teta (icha)*wmolat(iats)    &
                           +omega(icha)*wmolat(iatn) )
!
!Second membre
!
  sm(1) = alpha(icha)-gamma(icha)
  sm(2) = beta (icha)-delta(icha)
  sm(3) = teta (icha)-kappa(icha)
  sm(4) = 0.d0
  sm(5) = omega(icha)-zeta(icha)
  sm(6) = 0.d0
!
! On complete la matrice et le second menbre en fonction de
! la valeur de IY2CH
!
  if ( iy2ch(icha).eq.0 ) then
! Valeur de X2
    chx2(icha) = 2.D0
!Matrice
    mat(7,1) = -chx2(icha)
    mat(7,2) = 0.D0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
!Second membre
    sm(7) = 0.D0
    sm(8) = 0.D0
!
  else if ( iy2ch(icha).eq.1 ) then
!
!Matrice
    mat(7,1) = 0.d0
    mat(7,2) = 0.d0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 0.d0
    mat(7,8) = 1.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 1.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 0.d0
!Second membre
    sm(7) = y2ch(icha)
    sm(8) = 0.D0
!
  else if ( iy2ch(icha).eq.2 ) then
! Valeur de X2
    chx2(icha) = 2.D0
!
!Matrice
    mat(7,1) = -chx2(icha)
    mat(7,2) = 0.D0
    mat(7,3) = 0.d0
    mat(7,4) = 0.d0
    mat(7,5) = 0.d0
    mat(7,6) = 0.d0
    mat(7,7) = 1.d0
    mat(7,8) = 0.d0
!
    mat(8,1) = 0.D0
    mat(8,2) = 0.D0
    mat(8,3) = 0.d0
    mat(8,4) = 0.d0
    mat(8,5) = 0.d0
    mat(8,6) = 0.d0
    mat(8,7) = 0.d0
    mat(8,8) = 1.d0
!Second membre
    sm(7) = 0.D0
    sm(8) = y2ch(icha)
!
  endif
!
  call coal_resol_matrice( ndim, mat , sm , solu , ierr)
!
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
  write(nfecra,8013) a1(icha),b1(icha),c1(icha),d1(icha),e1(icha),f1(icha)
  write(nfecra,8014) a2(icha),b2(icha),c2(icha),d2(icha),e2(icha),f2(icha)
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
  dmf4  = ( oxyo2 (2)*wmo2 +oxyn2 (2)*wmn2     &
           +oxyh2o(2)*wmh2o+oxyco2(2)*wmco2 )
  if ( dmf4 .le. 0.d0 ) then
    write(nfecra,9897) oxyo2(2) ,oxyn2(2) ,    &
                       oxyh2o(2),oxyco2(2)
    call csexit(1)
  endif

  af4(io2)  = oxyo2(2)  / dmf4
  af4(in2)  = oxyn2(2)  / dmf4
  af4(ih2o) = oxyh2o(2) / dmf4
  af4(ico2) = oxyco2(2) / dmf4

endif

if ( noxyd .eq. 3.d0 ) then
  dmf5  = ( oxyo2 (3)*wmo2 +oxyn2 (3)*wmn2    &
           +oxyh2o(3)*wmh2o+oxyco2(3)*wmco2 )
  if ( dmf5 .le. 0.d0 ) then
    write(nfecra,9898) oxyo2(3) ,oxyn2(3) ,   &
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
af8(io2)   =-1.0d0/wmc
! coke par h2o
af9(ih2o)  =-1.0d0/wmc
af9(ico)   = 1.0d0/wmc
af9(ihy)   = 1.0d0/wmc

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


!===============================================================================
 8000 format(1X,' RECAPITULATIF SUR LES CHARBONS :'/,             &
       1X,' ==============================  '  )
 8001 format(/,3X,'Delta Hdev  ',/,                               &
       5X,' Numero Charbon',6X,'Legeres',7X,'Lourdes',9X,'PCI')
 8002 format(10x,i3,8x,g15.7,1x,g15.7,g15.7)

 8010 format(/,3X,'COMPOSITION MATIERES VOLATILES ')
 8011 format(/,5X,' Numero Charbon',I6)
 8012 format(18X,'CHX  ',12X,'CO',12X,'H2O',12X,'H2S',12X,'HCN',12X,'NH3')
 8013 format( 8X,' MV1 ',G15.7,G15.7,G15.7,G15.7,G15.7,G15.7)
 8014 format( 8X,' MV2 ',G15.7,G15.7,G15.7,G15.7,G15.7,G15.7)
 8020 format(/,3X,'LOI ENTHALPIE/TEMERATURE ')
 8021 format( 10X,'TEMPERATURE',9X,'Hch ',10X,'Hcoke',10X,        &
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
!===============================================================================

 9970 format(                                                      &
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
 9971 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@  Vous avec pris l''azote dans la composition du charbon    ',/,&
'@  mais votre repartition entre HCN et NH3 est nulle         ',/,&
'@  pour le charbon ',I2,'                                    ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9980 format(                                                     &
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
 9981 format(                                                     &
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
 9982 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,&
'@                                                            ',/,&
'@       MATRICE NON INVERSIBLE POUR LE CALCUL DES            ',/,&
'@                MATIERES VOLATILES                          ',/,&
'@  pour le charbon ',I2,'                                    ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9991 format(                                                     &
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
 9992 format(                                                     &
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
 9993 format(                                                     &
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
 9995 format(                                                     &
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
 9996 format(                                                     &
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
 9998 format(                                                     &
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
 9999 format(                                                     &
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

 9895 format(                                                     &
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

 9896 format(                                                     &
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
 9897 format(                                                     &
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
 9898 format(                                                     &
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
 9899 format(                                                     &
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

!----
! End
!----
end subroutine

!===============================================================================
!===============================================================================
!===============================================================================

subroutine coal_resol_matrice &
!============================
 ( NDIM, AA , BB , XX , IERR)
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
integer                                 , intent(in)    :: NDIM
double precision , dimension(NDIM,NDIM) , intent(inout) :: AA
double precision , dimension(NDIM)      , intent(inout) :: BB
double precision , dimension(NDIM)      , intent(out)   :: XX
integer                                 , intent(out)   :: IERR
! Variables locales
integer          :: ITERMX,IT,II,JJ,KK,IW
double precision :: EPSIL,ERR,X0,WW,PP,SS

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
