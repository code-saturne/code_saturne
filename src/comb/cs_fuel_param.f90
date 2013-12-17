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
use field

!===============================================================================

implicit none

integer          ipp , ii , jj , iok , icla
integer          isc
integer          kscmin
double precision wmolme

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! Define clippings for transported variables

do icla = 1,nclafu
  call field_set_key_double(ivarfl(ih2(icla)), kscmin, h02fol)
enddo

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

!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================


! ---> Variables algebriques propres a la suspension gaz - particules

ipp = ipppro(ipproc(immel))
nomprp(ipproc(immel))   = 'XM'
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! ---> Variables algebriques propres a la phase dispersee

do icla = 1, nclafu
  ipp = ipppro(ipproc(itemp2(icla)))
  write(nomprp(ipproc(itemp2(icla))),'(a7,i2.2)')'Tem_FOL' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(irom2(icla)))
  write(nomprp(ipproc(irom2(icla))),'(a7,i2.2)')'Rho_FOL' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(idiam2(icla)))
  write(nomprp(ipproc(idiam2(icla))),'(a6,i2.2)')'Dia_gt' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ih1hlf(icla)))
  write(nomprp(ipproc(ih1hlf(icla))),'(a6,i2.2)')'H1-Hlf' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmeva(icla)))
  write(nomprp(ipproc(igmeva(icla))),'(a6,i2.2)')'Ga_EVA' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(igmhtf(icla)))
  write(nomprp(ipproc(igmhtf(icla))),'(a6,i2.2)')'Ga_HET' ,icla
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
enddo

! ---> Variables algebriques propres a la phase continue

ipp = ipppro(ipproc(itemp1))
nomprp(ipproc(itemp1))   = 'Temp_GAZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(irom1))
nomprp(ipproc(irom1))   = 'ROM_GAZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(1)))
nomprp(ipproc(iym1(1)))   = 'YM_FO0'
ichrvr(ipp)   =  0
ilisvr(ipp)   =  0
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(2)))
nomprp(ipproc(iym1(2)))   = 'YM_FOV'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(3)))
nomprp(ipproc(iym1(3)))   = 'YM_CO'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(4)))
nomprp(ipproc(iym1(4)))   = 'YM_H2S'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(5)))
nomprp(ipproc(iym1(5)))   = 'YM_H2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(6)))
nomprp(ipproc(iym1(6)))   = 'YM_HCN'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(7)))
nomprp(ipproc(iym1(7)))   = 'YM_NH3'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(8)))
nomprp(ipproc(iym1(8)))   = 'YM_O2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(9)))
nomprp(ipproc(iym1(9)))   = 'YM_CO2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(10)))
nomprp(IPP)   = 'YM_H2O'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(11)))
nomprp(ipproc(iym1(11)))   = 'YM_SO2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iym1(12)))
nomprp(ipproc(iym1(12)))   = 'YM_N2'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
!
if ( ieqnox .eq. 1 ) then
  ipp = ipppro(ipproc(ighcn1))
  nomprp(ipproc(ighcn1))   = 'EXP1'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ighcn2))
  nomprp(ipproc(ighcn2))   = 'EXP2'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
  ipp = ipppro(ipproc(ignoth))
  nomprp(ipproc(ignoth))   = 'EXP3'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif
!
ipp = ipppro(ipproc(ibcarbone))
nomprp(ipproc(ibcarbone))   = 'Bilan_C'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(iboxygen))
nomprp(ipproc(iboxygen))   = 'Bilan_O'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
ipp = ipppro(ipproc(ibhydrogen))
nomprp(ipproc(ibhydrogen))   = 'Bilan_H'
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
visls0(iscalt) = diftl0

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
