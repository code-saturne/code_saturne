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

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use cs_c_bindings

!===============================================================================

implicit none

integer          ii , jj , iok , icla
integer          isc
integer          kscmin
double precision wmolme, turb_schmidt

type(var_cal_opt) :: vcopt

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

  turb_schmidt = 0.7d0
  call field_set_key_double(ivarfl(isca(jj)), ksigmas, turb_schmidt)

! ------ Coeff dissipation des fluctuations

  rvarfl(jj) = 0.8d0

  ii = isca(iscapp(isc))

  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)

! ------ Informations relatives a la resolution des scalaires

!         - Facteur multiplicatif du pas de temps
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

enddo

!===============================================================================
! 2. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! ---> Initialisation

! ---- Calcul de RO0 a partir de T0 et P0
!        (loi des gaz parfaits applliquee a l'air)
!    On initialise RO0 avec l'oxydant 1 qui est sense etre
!    l'oxydant majoritaire

wmolme = ( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)             &
          +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))           &
        /(oxyo2(1)+oxyn2(1)+oxyh2o(1)+oxyco2(1))
ro0 = p0*wmolme / (cs_physical_constants_r*t0)

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
! 3. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

call cs_user_combustion

!===============================================================================
! 4. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
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
'                                        (cs_user_combustion).',/)
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
'@  Verifier cs_user_combustion.'                              ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine
