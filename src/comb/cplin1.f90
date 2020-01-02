!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine cplin1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!       INIT DES OPTIONS DES VARIABLES TRANSPORTEES
!                     ET DES VARIABLES ALGEBRIQUES

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
use ppcpfu
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ii , jj, iok
integer          icha , isc , is
double precision wmolme, turb_schmidt

type(var_cal_opt) :: vcopt

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! --> Nature des scalaires transportes

do isc = 1, nscapp

! ---- Type de scalaire (0 passif, 1 temperature en K
!                                 -1 temperature en C
!                                  2 enthalpie)
!      La distinction -1/1 sert pour le rayonnement
  iscacp(iscapp(isc)) = 0

enddo

! --> Donnees physiques ou numeriques propres aux scalaires CP

do isc = 1, nscapp

  jj = iscapp(isc)

  if ( iscavr(jj).le.0 ) then

! ---- En combustion on considere que la viscosite turbulente domine
!      ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    visls0(jj) = viscl0

  endif

! ---- Schmidt ou Prandtl turbulent

  turb_schmidt = 0.7d0
  call field_set_key_double(ivarfl(isca(jj)), ksigmas, turb_schmidt)

! ---- Coeff dissipation des fluctuations

  rvarfl(jj) = 0.8d0

  ii = isca(iscapp(isc))

  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)

! ---- Informations relatives a la resolution des scalaires

!       - Facteur multiplicatif du pas de temps
  cdtvar(ii) = 1.d0

!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
  vcopt%blencv = 0.d0

!         - Type de schema convetif second ordre (utile si BLENCV > 0)
!           = 0 : Second Order Linear Upwind
!           = 1 : Centre
  vcopt%ischcv = 1

!         - Test de pente pour basculer d'un schema centre vers l'upwind
!           = 0 : utilisation automatique du test de pente
!           = 1 : calcul sans test de pente
  vcopt%isstpc = 0

!         - Reconstruction des flux de convetion et de diffusion aux faces
!           = 0 : pas de reconstruction
  vcopt%ircflu = 0

  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)

enddo

!===============================================================================
! 2. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! Definition des pointeurs du tableau TBMCR utilise dans cpphy1.f90
! et dans les sous-programmes appeles

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

! ---> Initialisation

! ---- Calcul de RO0 a partir de T0 et P0
!        (loi des gaz parfaits appliquee a l'air)

wmolme = (wmole(io2)+xsi*wmole(in2)) / (1.d0+xsi)
ro0 = p0*wmolme / (cs_physical_constants_r*t0)

! ---- Initialisation pour la masse volumique du coke

do icha = 1, ncharb
  rhock(icha) = rho0ch(icha)
enddo

! On met SRROM et DIFTL0 a -GRAND pour forcer l'utilisateur a les
! definir dans cplin1
! ---> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom =-grand

! ---> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0      = -grand

! ---> Masse volumique variable et viscosite constante (pour les suites)
irovar = 1
ivivar = 0

!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

call cs_user_combustion

!===============================================================================
! 4. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
call cplver (iok)
!==========

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
'                                        (cs_user_combustion).',/)
 9999 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (C.P. COUPLE LAGRANGIEN) DEMANDEE ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier cs_user_combustion.'                              ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
