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

subroutine cpini1
!================

!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES OPTIONS DES VARIABLES TRANSPORTEES
!                     ET DES VARIABLES ALGEBRIQUES
!               COMBUSTION CHARBON PULVERISE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
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
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "cstnum.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "ihmpre.h"
include "ppcpfu.h"

!===============================================================================

integer          ipp , icla , ii , jj , iok
integer          icha , isc , is, iphas
double precision wmolme

!===============================================================================
!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

do iphas = 1, nphas
  if(iscalt(iphas).ne.-1) then
    write(nfecra,1000)iphas,iscalt(iphas)
    call csexit (1)
    !==========
  endif
enddo
do ii = 1, nscapp
  if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscapp(ii),iscsth(iscapp(ii))
    call csexit (1)
    !==========
  endif
enddo

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! --> Definition des scamin et des scamax des variables transportees

! ---- Variables propres a la suspension gaz - particules

scamin(ihm)   = -grand
scamax(ihm)   = +grand

! ---- Variables propres a la phase dispersee

do icla = 1, nclacp


!        SCAMIN(IH2(ICLA)) = EH0SOL(ICH(ICHCOR(ICLA)))
!    Prise en compte de l'humidite :
!        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
!          SCAMIN(IH2(ICLA)) = EH0SOL(IWAT(ICHCOR(ICLA)))
!        ENDIF

  scamin(ih2(icla)) = -grand

  scamax(ih2(icla)) = +grand
  scamin(inp(icla)) = 0.d0
  scamax(inp(icla)) = +rinfin
  scamin(ixch(icla)) = 0.d0
  scamax(ixch(icla)) = 1.d0
  scamin(ixck(icla)) = 0.d0
  scamax(ixck(icla)) = 1.d0
  if ( ippmod(icp3pl) .eq. 1 ) then
    scamin(ixwt(icla)) = 0.d0
    scamax(ixwt(icla)) = 1.d0
  endif
enddo
do icha = 1, ncharb
  scamin(if1m(icha)) = 0.d0
  scamax(if1m(icha)) = 1.d0
  scamin(if2m(icha)) = 0.d0
  scamax(if2m(icha)) = 1.d0
enddo

! ---- Variables propres a la phase continue

scamin(if3m) = 0.d0
scamax(if3m) = 1.d0
if ( ihtco2 .eq. 1 ) then
  scamin(if3mc2) = 0.d0
  scamax(if3mc2) = 1.d0
endif
scamin(if4p2m) = 0.d0
scamax(if4p2m) = 0.25d0
if ( ippmod(icp3pl) .eq. 1 ) then
  scamin(if5m) = 0.d0
  scamax(if5m) = 1.d0
endif

!    Oxycombustion

if ( noxyd .ge. 2 ) then
  scamin(if6m) = 0.d0
  scamax(if6m) = 1.d0
endif
if ( noxyd .eq. 3 ) then
  scamin(if7m) = 0.d0
  scamax(if7m) = 1.d0
endif

!    Modele de CO

if ( ieqco2 .ge. 1 ) then
  scamin(iyco2) = 0.d0
  scamax(iyco2) = 1.d0
endif

! --> Nature des scalaires transportes

do isc = 1, nscapp

! ---- Type de scalaire (0 passif, 1 temperature en K
!                                 -1 temperature en C
!                                  2 enthalpie)
!      La distinction -1/1 sert pour le rayonnement
  iscsth(iscapp(isc)) = 0

enddo


! ---- On resout en enthalpie avec un CP constant (Cf. cpvarp)

iphas = iphsca(ihm)
iscalt(iphas) = ihm
iscsth(ihm) = 2

! --> Donnees physiques ou numeriques propres aux scalaires CP

do isc = 1, nscapp

  jj = iscapp(isc)

  if ( iscavr(jj).le.0 ) then

!        En combustion on considere que la viscosite turbulente domine
!        ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    visls0(jj) = viscl0(iphsca(jj))

  endif

! ------ Schmidt ou Prandtl turbulent

  sigmas(jj) = 0.7d0

! ------ Coeff dissipation des fluctuations

  rvarfl(jj) = 0.8d0

  ii = isca(iscapp(isc))

! ------ Niveau de detail des impressions pour les variables et
!          donc les scalaires (valeurs 0 ou 1)
!          Si = -10000 non modifie par l'utilisateur -> niveau 1
  if(iwarni(ii).eq.-10000) then
    iwarni(ii) = 1
  endif

!   - Interface Code_Saturne:
!     ======================

!   NOMVAR, ICHRVR,ILISVR, IHISVR are
!   already filled in UINUM1 routine

  if (iihmpr.ne.1) then

! ------ Informations relatives a la resolution des scalaires

!         - Facteur multiplicatif du pas de temps
  cdtvar(ii) = 1.d0

!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
  blencv(ii) = 0.d0

!         - Type de schema convetif second ordre (utile si BLENCV > 0)
!           = 0 : Second Order Linear Upwind
!           = 1 : Centre
  ischcv(ii) = 1

!         - Test de pente pour basculer d'un schema centre vers l'upwind
!           = 0 : utlisation automatique du test de pente
!           = 1 : calcul sans test de pente
  isstpc(ii) = 0

!         - Reconstruction des flux de convetion et de diffusion aux faces
!           = 0 : pas de reconstruction
  ircflu(ii) = 0

  endif

enddo

!   - Interface Code_Saturne:
!     ======================

!   NOMVAR, ICHRVR,ILISVR, IHISVR are
!   already filled in CSENSO routine

if (iihmpr.ne.1) then

!---> Variable courante : nom, sortie chrono, suivi listing, sortie histo

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     NOMVAR( ) = nom de la variable
!     ICHRVR( ) = sortie chono (oui 1/non 0)
!     ILISVR( ) = suivi listing (oui 1/non 0)
!     IHISVR( ) = sortie historique (nombre de sondes et numeros)
!     si IHISVR(.,1)  = -1 sortie sur toutes les sondes

!     NB : Les 8 premiers caracteres du noms seront repris dans le
!          listing 'developpeur'


! ---- Variables propres a la suspension gaz - particules

ipp = ipprtp(isca(ihm))
NOMVAR(IPP)  = 'Enthalpy'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Variables propres a la phase dispersee

do icla = 1, nclacp
  ipp = ipprtp(isca(ixck(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'XCK_CP' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ixch(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'XCH_CP' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(inp(icla)))
  WRITE(NOMVAR(IPP),'(A5,I2.2)')'NP_CP' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ih2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'ENT_CP' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  if ( ippmod(icp3pl) .eq. 1 ) then
    ipp = ipprtp(isca(ixwt(icla)))
    WRITE(NOMVAR(IPP),'(A6,I2.2)')'XWT_CP' ,ICLA
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
enddo
do icha = 1, ncharb
  ipp = ipprtp(isca(if1m(icha)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Fr_MV1' ,ICHA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(if2m(icha)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Fr_MV2' ,ICHA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo

! ---- Variables propres a la phase continue

ipp = ipprtp(isca(if3m))
NOMVAR(IPP)  = 'Fr_HET_O2'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1
if ( ihtco2 .eq. 1 ) then
  ipp = ipprtp(isca(if3mc2))
  NOMVAR(IPP)  = 'Fr_HET_CO2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif
ipp = ipprtp(isca(if4p2m))
NOMVAR(IPP)  = 'Var_AIR'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1
if ( ippmod(icp3pl) .eq. 1 ) then
  ipp = ipprtp(isca(if5m))
  NOMVAR(IPP)  = 'FR_H2O'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if ( noxyd .ge. 2 ) then
  ipp = ipprtp(isca(if6m))
  NOMVAR(IPP)  = 'FR_OXYD2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif
if ( noxyd .eq. 3 ) then
  ipp = ipprtp(isca(if7m))
  NOMVAR(IPP)  = 'FR_OXYD3'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if ( ieqco2 .ge. 1 ) then
  ipp = ipprtp(isca(iyco2))
  NOMVAR(IPP)  = 'FR_CO2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================


! ---> Variables algebriques propres a la suspension gaz - particules

ipp = ipppro(ipproc(immel))
NOMVAR(IPP)   = 'XM'
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! ---> Variables algebriques propres a la phase dispersee

do icla = 1, nclacp
  ipp = ipppro(ipproc(itemp2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Tem_CP' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(irom2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Rho_CP' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(idiam2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Dia_CK' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmdch(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_DCH' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmdv1(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_DV1' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmdv2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_DV2' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmhet(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_HET_O2' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  if ( ihtco2 .eq. 1 ) then
    ipp = ipppro(ipproc(ighco2(icla)))
    WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_HET_CO2' ,ICLA
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  if ( ippmod(icp3pl) .eq. 1 ) then
    ipp = ipppro(ipproc(igmsec(icla)))
    WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_SEC' ,ICLA
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ipp = ipppro(ipproc(ix2(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'Frm_CP' ,ICLA
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
enddo

! ---> Variables algebriques propres a la phase continue

ipp = ipppro(ipproc(itemp1))
NOMVAR(IPP)   = 'Temp_GAZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(irom1))
NOMVAR(IPP)   = 'ROM_GAZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(1)))
NOMVAR(IPP)   = 'YM_CHx1m'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(2)))
NOMVAR(IPP)   = 'YM_CHx2m'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(3)))
NOMVAR(IPP)   = 'YM_CO'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(4)))
NOMVAR(IPP)   = 'YM_O2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(5)))
NOMVAR(IPP)   = 'YM_CO2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(6)))
NOMVAR(IPP)   = 'YM_H2O'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(7)))
NOMVAR(IPP)   = 'YM_N2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

endif



!===============================================================================
! 3. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! ---> Definition des pointeurs du tableau TBMCR utilise dans cpphy1.F
!      et dans les sous-programmes appeles

is = 0
do icha = 1, ncharb
  is          = is + 1
  if1mc(icha) = is
  is          = is + 1
  if2mc(icha) = is
enddo
is      = is + 1
ix1mc   = is
is      = is + 1
ix2mc   = is
is      = is + 1
ichx1f1 = is
is      = is + 1
ichx2f2 = is
is      = is + 1
icof1   = is
is      = is + 1
icof2   = is
is      = is + 1
ih2of1  = is
is      = is + 1
ih2of2  = is


! ---> Initialisation

! ---- Calcul de RO0 a partir de T0 et P0
!        (loi des gaz parfaits applliquee a l'air)

!    On initialise RO0 avec l'oxydant 1 qui est sense etre
!    l'oxydant majoritaire

wmolme = ( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)             &
          +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))           &
        /(oxyo2(1)+oxyn2(1)+oxyh2o(1)+oxyco2(1))

ro0(iphas) = p0(iphas)*wmolme / (rr*t0(iphas))

! ---- Initialisation pour la masse volumique du coke

do icha = 1, ncharb
  rhock(icha) = rho0ch(icha)
enddo

! On met SRROM et DIFTL0 a -GRAND pour forcer l'utilisateur a les
! definir dans uscpi1
! ---> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = -grand

! ---> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
!     C'est cette valeur par defaut qui est TOUJOURS utilisee dans les
!       calculs charbon (un peu etonnant de ne pas prendre en
!       compte les variations de ce parametre physique si on
!       recherche des informations sur les flux thermiques aux parois)

diftl0 =-grand

! ---> Masse volumique variable et viscosite constante (pour les suites)
irovar(iphas) = 1
ivivar(iphas) = 0

!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call uicpi1(srrom)
  !==========

  diftl0 = 4.25d-5

endif

call uscpi1
!==========

!===============================================================================
! 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
call cpveri (iok)
!==========

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
  !==========
else
  write(nfecra,9998)
endif

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (C.P.) DEMANDEE                   ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (C.P.) DEMANDEE                   ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                           &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (uscpi1).',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier uscpi1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
