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

subroutine cs_coal_param
!=======================

!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES OPTIONS DES VARIABLES TRANSPORTEES
!                     ET DES VARIABLES ALGEBRIQUES
!               COMBUSTION CHARBON PULVERISE

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
use ppincl
use ihmpre
use ppcpfu
use cs_coal_incl
use field

!===============================================================================

implicit none

integer          ipp , icla , ii , jj , iok
integer          icha , isc

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

! ---- Variables propres a la suspension gaz - particules

scamin(ihm)   = -grand
scamax(ihm)   = +grand

! ---- Variables propres a la phase dispersee

do icla = 1, nclacp
  scamin(ih2(icla)) = -grand
  scamax(ih2(icla)) = +grand
  scamin(inp(icla)) = 0.d0
  scamax(inp(icla)) = +rinfin
  scamin(ixch(icla)) = 0.d0
  scamax(ixch(icla)) = 1.d0
  scamin(ixck(icla)) = 0.d0
  scamax(ixck(icla)) = 1.d0
  if ( ippmod(iccoal) .eq. 1 ) then
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
if ( noxyd .ge. 2 ) then
  ! Oxydant 2
  scamin(if4m) = 0.d0
  scamax(if4m) = 1.d0
endif
!
if ( noxyd .eq. 3 ) then
  ! Oxydant 3
  scamin(if5m) = 0.d0
  scamax(if5m) = 1.d0
endif
! eau libre
if ( ippmod(iccoal) .ge. 1 ) then
  scamin(if6m) = 0.d0
  scamax(if6m) = 1.d0
endif
! Produits de la combustion du coke par O2
scamin(if7m) = 0.d0
scamax(if7m) = 1.d0
! Produits de la combustion du coke par CO2
if ( ihtco2 .eq. 1 ) then
  scamin(if8m) = 0.d0
  scamax(if8m) = 1.d0
endif
! Produits de la combustion du coke par H2O
if ( ihth2o .eq. 1 ) then
  scamin(if9m) = 0.d0
  scamax(if9m) = 1.d0
endif
! Variance
scamin(ifvp2m) = 0.d0
scamax(ifvp2m) = 0.25d0
!    Modele de CO
if ( ieqco2 .eq. 1 ) then
  scamin(iyco2) = 0.d0
  scamax(iyco2) = 1.d0
endif
if ( ieqnox .ge. 1 ) then
  scamin(iyhcn) = 0.d0
  scamax(iyhcn) = 1.d0
! On donne les limites de la variable transportee qui represente le NH3.
  scamin(iynh3) = 0.d0
  scamax(iynh3) = 1.d0

  scamin(iyno) = 0.d0
  scamax(iyno) = 1.d0
  scamin(ihox) = -grand
  scamax(ihox) = +grand
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
iscsth(ihm) = 2

! --> Donnees physiques ou numeriques propres aux scalaires CP

do isc = 1, nscapp

  jj = iscapp(isc)

  if ( iscavr(jj).le.0 ) then

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
  nomvar(ipp)  = 'Enthalpy'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!TODO remove once the field are totally deployed, up to now
!     it is used to not create twice the same field.
! ---- Variables propres a la phase dispersee

  do icla = 1, nclacp
    ipp = ipprtp(isca(ixck(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'xck_cp' ,icla
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(ixch(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Xch_CP' ,icla
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(inp(icla)))
    write(nomvar(ipp),'(a5,i2.2)')'Np_CP' ,icla
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(ih2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Ent_CP' ,icla
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    if (ippmod(iccoal) .eq. 1) then
      ipp = ipprtp(isca(ixwt(icla)))
      write(nomvar(ipp),'(a6,i2.2)')'Xwt_CP' ,icla
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    endif
    if (i_coal_drift.eq.1) then
      ipp = ipprtp(isca(iagecp_temp(icla)))
      write(nomvar(ipp),'(a8,i2.2)')'X_Age_CP' ,icla
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)=-1
    endif
  enddo
  do icha = 1, ncharb
    ipp = ipprtp(isca(if1m(icha)))
    write(nomvar(ipp),'(a6,i2.2)')'Fr_mv1' ,icha
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(if2m(icha)))
    write(nomvar(ipp),'(a6,i2.2)')'Fr_mv2' ,icha
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
! ---- Variables propres a la phase continue
  if (i_coal_drift.eq.1) then
    ipp = ipprtp(isca(iaggas_temp))
    nomvar(ipp)  = 'X_Age_Gas'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( noxyd .ge. 2 ) then
    ipp = ipprtp(isca(if4m))
    nomvar(ipp)  = 'FR_OXYD2'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( noxyd .eq. 3 ) then
    ipp = ipprtp(isca(if5m))
    nomvar(ipp)  = 'FR_OXYD3'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ippmod(iccoal) .ge. 1 ) then
    ipp = ipprtp(isca(if6m))
    nomvar(ipp)  = 'FR_H2O'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  ipp = ipprtp(isca(if7m))
  nomvar(ipp)  = 'FR_HET_O2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  if ( ihtco2 .eq. 1 ) then
    ipp = ipprtp(isca(if8m))
    nomvar(ipp)  = 'FR_HET_CO2'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ihth2o .eq. 1 ) then
    ipp = ipprtp(isca(if9m))
    nomvar(ipp)  = 'FR_HET_H2O'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
!
  ipp = ipprtp(isca(ifvp2m))
  NOMVAR(IPP)  = 'Var_F1F2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
!
  if ( ieqco2 .eq. 1 ) then
    ipp = ipprtp(isca(iyco2))
    NOMVAR(IPP)  = 'FR_CO2'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if ( ieqnox .ge. 1 ) then
    ipp = ipprtp(isca(iyhcn))
    NOMVAR(IPP)  = 'FR_HCN'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

    ipp = ipprtp(isca(iynh3))
    NOMVAR(IPP)  = 'FR_NH3'
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
  do icla = 1, nclacp
    ipp = ipppro(ipproc(itemp2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Tem_CP' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(irom2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Rho_CP' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(idiam2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Dia_CK' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(igmdch(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Ga_DCH' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(igmdv1(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Ga_DV1' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(igmdv2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Ga_DV2' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(igmhet(icla)))
    write(nomvar(ipp),'(a9,i2.2)')'Ga_HET_O2' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    if ( ihtco2 .eq. 1 ) then
      ipp = ipppro(ipproc(ighco2(icla)))
      write(nomvar(ipp),'(a10,i2.2)')'Ga_HET_CO2' ,icla
      ichrvr(ipp)   = 1
      ilisvr(ipp)   = 1
      ihisvr(ipp,1) = -1
    endif
    if ( ihth2o .eq. 1 ) then
      ipp = ipppro(ipproc(ighh2o(icla)))
      write(nomvar(ipp),'(a10,i2.2)')'Ga_HET_H2O' ,icla
      ichrvr(ipp)   = 1
      ilisvr(ipp)   = 1
      ihisvr(ipp,1) = -1
    endif
    if ( ippmod(iccoal) .eq. 1 ) then
      ipp = ipppro(ipproc(igmsec(icla)))
      write(nomvar(ipp),'(a6,i2.2)')'Ga_SEC' ,icla
      ichrvr(ipp)   = 1
      ilisvr(ipp)   = 1
      ihisvr(ipp,1) = -1
    endif
    ipp = ipppro(ipproc(ix2(icla)))
    write(nomvar(ipp),'(a6,i2.2)')'Frm_CP' ,icla
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo
! ---> Variables algebriques propres a la phase continue
  ipp = ipppro(ipproc(itemp1))
  nomvar(IPP)   = 'Temp_GAZ'
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
  if ( ieqnox .eq. 1 ) then
    ipp = ipppro(ipproc(ighcn1))
    nomvar(IPP)   = 'EXP1'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ighcn2))
    nomvar(IPP)   = 'EXP2'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ignoth))
    nomvar(IPP)   = 'EXP3'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ignh31))
    nomvar(IPP)   = 'EXP4'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ignh32))
    nomvar(IPP)   = 'EXP5'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifhcnd))
    nomvar(IPP)   = 'F_HCN_DEV'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifhcnc))
    nomvar(IPP)   = 'F_HCN_HET'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnh3d))
    nomvar(IPP)   = 'F_NH3_DEV'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnh3c))
    nomvar(IPP)   = 'F_NH3_HET'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnohc))
    nomvar(IPP)   = 'F_NO_HCN'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnonh))
    nomvar(IPP)   = 'F_NO_NH3'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnoch))
    nomvar(IPP)   = 'F_NO_HET'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifnoth))
    nomvar(IPP)   = 'F_NO_THE'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(icnohc))
    nomvar(IPP)   = 'C_NO_HCN'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(icnonh))
    nomvar(IPP)   = 'C_NO_NH3'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(ifhcnr))
    nomvar(IPP)   = 'F_HCN_RB'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(icnorb))
    nomvar(IPP)   = 'C_NO_RB'
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
    ipp = ipppro(ipproc(igrb))
    nomvar(IPP)   = 'EXP_RB'
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
endif
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

do icha = 1, ncharb
  rhock(icha) = rho0ch(icha)
enddo

! On met SRROM et DIFTL0 a -GRAND pour forcer l'utilisateur a les
! definir dans user_coal_ini1
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
irovar = 1
ivivar = 0

!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call uicpi1(srrom, diftl0)
  !==========

  diftl0 = 4.25d-5

endif

call user_coal_ini1
!==================

!===============================================================================
! 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
call cs_coal_verify (iok)
!=====================

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
  !==========
else
  write(nfecra,9998)
endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (C.P.) DEMANDEE                   ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner, or              ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT(',I10   ,') = ',I10                               ,/,&
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
'@    PHYSIQUE PARTICULIERE (C.P.) DEMANDEE                   ',/,&
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
'                                            (user_coal_ini1).',/)
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
'@  Verifier user_coal_ini1.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! End
!----

return
end subroutine
