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

subroutine cs_fuel_param
!=======================

!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES OPTIONS DES VARIABLES TRANSPORTEES
!                     ET DES VARIABLES ALGEBRIQUES
!               COMBUSTION FUEL

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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

integer          ipp , ii , jj , iok , icla
integer          isc
double precision wmolme

!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

if(iscalt.ne.-1) then
  write(nfecra,1000)iscalt
  call csexit (1)
  !==========
endif
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

! ---- Variables propres a la phase dispersee

scamin(ihm) = -grand
scamax(ihm) = +grand

do icla=1,nclafu
  scamin(ing(icla))   = 0.d0
  scamax(ing(icla))   = +rinfin
  scamin(iyfol(icla)) = 0.d0
  scamax(iyfol(icla)) = 4.d-1
  scamin(ih2(icla))   = h02fol
  scamax(ih2(icla))   = +grand
enddo

! ---- Variables propres a la phase continue

scamin(ifvap)  = 0.d0
scamax(ifvap)  = 1.d0
!
! Oxydant
if ( noxyd .ge. 2 ) then
  scamin(if4m) = 0.d0
  scamax(if4m) = 1.d0
endif
!
if ( noxyd .eq. 3 ) then
  scamin(if5m) = 0.d0
  scamax(if5m) = 1.d0
endif
! Combustion heterogene
scamin(if7m)   = 0.d0
scamax(if7m)   = 1.d0
! Variance
scamin(ifvp2m)  = 0.d0
scamax(ifvp2m)  = 0.25d0
! Modele de CO
if ( ieqco2 .ge. 1 ) then
  scamin(iyco2) = 0.d0
  scamax(iyco2) = 1.d0
endif
! Modele de Nox
if ( ieqnox .eq. 1 ) then
  scamin(iyhcn)  = 0.d0
  scamax(iyhcn)  = 1.d0
  scamin(iyno)   = 0.d0
  scamax(iyno)   = 1.d0
  scamin(ihox) = -grand
  scamax(ihox) = +grand
  scamin(iynh3) = 0.d0
  scamax(iynh3) = 1.d0
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

iscalt = ihm
iscsth(ihm)   = 2

! --> Donnees physiques ou numeriques propres aux scalaires CP

do isc = 1, nscapp

  jj = iscapp(isc)

  if ( iscavr(jj) .le. 0 ) then

!        En combustion on considere que la viscosite turbulente domine
!        ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    visls0(jj) = viscl0

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

enddo


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

do icla = 1, nclafu
  ipp = ipprtp(isca(iyfol(icla)))
  WRITE(NOMVAR(IPP),'(A8,I2.2)')'YFOL_FOL' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ih2(icla)))
  WRITE(NOMVAR(IPP),'(A7,I2.2)')'HLF_FOL' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ing(icla)))
  WRITE(NOMVAR(IPP),'(A6,I2.2)')'NG_FOL' ,ICLA
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo

! ---- Variables propres a la phase gaz

ipp = ipprtp(isca(ifvap))
NOMVAR(IPP)  = 'Fr_VAP'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1
ipp = ipprtp(isca(if7m))
NOMVAR(IPP)  = 'Fr_HET'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1
ipp = ipprtp(isca(ifvp2m))
NOMVAR(IPP)  = 'Var_CB'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

if ( ieqco2 .ge. 1 ) then
  ipp = ipprtp(isca(iyco2))
  NOMVAR(IPP)  = 'FR_CO2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if ( ieqnox .eq. 1 ) then
  ipp = ipprtp(isca(iyhcn))
  NOMVAR(IPP)  = 'FR_HCN'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(iyno))
  NOMVAR(IPP)  = 'FR_NO'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ihox))
  NOMVAR(IPP)  = 'Enth_Ox'
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

do icla = 1, nclafu
  ipp = ipppro(ipproc(itemp2(icla)))
  write(nomvar(ipp),'(a7,i2.2)')'Tem_FOL' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(irom2(icla)))
  write(nomvar(ipp),'(a7,i2.2)')'Rho_FOL' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(idiam2(icla)))
  write(nomvar(ipp),'(a6,i2.2)')'Dia_gt' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ih1hlf(icla)))
  write(nomvar(ipp),'(a6,i2.2)')'H1-Hlf' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmeva(icla)))
  write(nomvar(ipp),'(a6,i2.2)')'Ga_EVA' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmhtf(icla)))
  write(nomvar(ipp),'(a6,i2.2)')'Ga_HET' ,icla
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
NOMVAR(IPP)   = 'YM_FO0'
ichrvr(ipp)   =  0
ilisvr(ipp)   =  0
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(2)))
NOMVAR(IPP)   = 'YM_FOV'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(3)))
NOMVAR(IPP)   = 'YM_CO'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(4)))
NOMVAR(IPP)   = 'YM_H2S'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(5)))
NOMVAR(IPP)   = 'YM_H2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(6)))
NOMVAR(IPP)   = 'YM_HCN'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(7)))
NOMVAR(IPP)   = 'YM_NH3'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(8)))
NOMVAR(IPP)   = 'YM_O2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(9)))
NOMVAR(IPP)   = 'YM_CO2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(10)))
NOMVAR(IPP)   = 'YM_H2O'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(11)))
NOMVAR(IPP)   = 'YM_SO2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(12)))
NOMVAR(IPP)   = 'YM_N2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
!
if ( ieqnox .eq. 1 ) then
  ipp = ipppro(ipproc(ighcn1))
  NOMVAR(IPP)   = 'EXP1'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ighcn2))
  NOMVAR(IPP)   = 'EXP2'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ignoth))
  NOMVAR(IPP)   = 'EXP3'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif
!
ipp = ipppro(ipproc(ibcarbone))
NOMVAR(IPP)   = 'Bilan_C'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iboxygen))
NOMVAR(IPP)   = 'Bilan_O'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(ibhydrogen))
NOMVAR(IPP)   = 'Bilan_H'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
!
!===============================================================================
! 3. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! ---> Initialisation

! ---- Calcul de RO0 a partir de T0 et P0
!        (loi des gaz parfaits applliquee a l'air)
!    On initialise RO0 avec l'oxydant 1 qui est sense etre
!    l'oxydant majoritaire

wmolme = ( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)             &
          +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))           &
        /(oxyo2(1)+oxyn2(1)+oxyh2o(1)+oxyco2(1))
ro0 = p0*wmolme / (rr*t0)

! ---- Initialisation pour la masse volumique du coke
!    Maintenant c'est fait dans FULECD
!     RHOKF = RHO0FL

! ---> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.90d0

! ---> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
!     C'est cette valeur par defaut qui est TOUJOURS utilisee dans les
!       calculs charbon (un peu etonnant de ne pas prendre en
!       compte les variations de ce parametre physique si on
!       recherche des informations sur les flux thermiques aux parois)

diftl0      = 4.25d-5
visls0(ihm) = diftl0

! ---> Masse volumique variable et viscosite constante (pour les suites)
irovar = 1
ivivar = 0

!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

call user_fuel_ini1
!==================

!===============================================================================
! 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
call cs_fuel_verify(iok)
!==================

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
  !==========
else
  write(nfecra,9998)
endif


 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (FUEL) DEMANDEE                   ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner, or              ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT = ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (FUEL) DEMANDEE                   ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner, or             ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9998 format(                                                     &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                            (user_fuel_ini1).',/)
 9999 format(                                                     &
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
'@  Verifier user_fuel_ini1.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine
