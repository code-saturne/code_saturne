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

subroutine cfini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!              LE COMPRESSIBLE SANS CHOC
!   EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS cs_user_parameters.f90

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
use ppincl
use ihmpre
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ii
integer          iok, kscmin, kscmax
double precision  scaclp(4)

type(var_cal_opt) :: vcopt

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! 1.1 Definition des scamin et des scamax des variables transportees
! ==================================================================

! Key id for scamin and scamax
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

call field_get_key_double(ivarfl(isca(ienerg)), kscmin, scaclp(1))
call field_get_key_double(ivarfl(isca(itempk)), kscmin, scaclp(2))
call field_get_key_double(ivarfl(isca(ienerg)), kscmax, scaclp(3))
call field_get_key_double(ivarfl(isca(itempk)), kscmax, scaclp(4))

if ( (abs(scaclp(1)+grand).gt.epzero).or.           &
     (abs(scaclp(2)+grand).gt.epzero).or.           &
     (abs(scaclp(3)-grand).gt.epzero).or.           &
     (abs(scaclp(4)-grand).gt.epzero) ) then
  write(nfecra,2000) scaclp(1), scaclp(3), scaclp(2), scaclp(4)
  call csexit (1)
endif

! 1.2 Nature des scalaires transportes
! ====================================

! Does scalar itempk behave like a temperature ?
! TODO check this; should be 1 for temperature unless handled in
!      another manner

iscacp(itempk) = 0

!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
do ii = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  vcopt%blencv = 0.d0
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
enddo

!===============================================================================
! 2. PARAMETRES GLOBAUX
!===============================================================================

! --- Couplage vitesse/pression (0 : algorithme classique,
!                                1 : couplage instationnaire)
!     Uniquement en monophasique et en incompressible

if( ipucou.ne.0 ) then
  write(nfecra,3000) ipucou
  call csexit (1)
endif

! --- Estimateurs pour Navier-Stokes

!     Interdits en compressible

if( (iescal(iespre).ne.0) .or.                            &
    (iescal(iesder).ne.0) .or.                            &
    (iescal(iescor).ne.0) .or.                            &
    (iescal(iestot).ne.0) ) then
  write(nfecra,4000)
  call csexit (1)
endif

!===============================================================================
! 3. OPTIONS DE CALCUL PAR DEFAUT
!===============================================================================

! --> Conditions aux limites prenant en compte l'equilibre hydrostatique
!     (oui = 1 , non = 0)

icfgrp = 1


! ---> Masse volumique variable (pour les suites)
irovar = 1

!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

if (iihmpr.eq.1) then
  call cscfgp(icfgrp)
endif

call uscfx2
!==========

!===============================================================================
! 5. OPTIONS DE CALCUL OBLIGATOIRES
!     qui pourront etre remontees au dessus de uscfx1
!     selon les developpements
!===============================================================================

! --> Prise en compte de la pression predite pour resoudre Navier-Stokes
!     (oui = 1 , non = 0)

igrdpp = 1

! --> Prediction de pression par une equation d'evolution

!     ATTENTION   PAS ENCORE IMPLEMENTE
!========   LAISSER IPPRED = 0

ippred = 0


!===============================================================================
! 6. VERIFICATIONS
!===============================================================================

iok = 0
if(icfgrp.ne.0.and.icfgrp.ne.1) then
  write(nfecra,5000)'ICFGRP',icfgrp
  iok = 1
endif

if (iok.ne.0) then
  call csexit (1)
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 2000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  Les bornes des variables energie ou temperature           ',/,&
'@    ont ete modifiees :                                     ',/,&
'@                                                            ',/,&
'@                      SCAMIN        SCAMAX                  ',/,&
'@  energie     ',2E14.5                                       ,/,&
'@  temperature ',2E14.5                                       ,/,&
'@                                                            ',/,&
'@  Les bornes de ces variables ne doivent pas etre modifiees.',/,&
'@  On peut modifier les bornes des variables rho et energie  ',/,&
'@  dans uscfx1, mais ce n''est pas conseille.                ',/,&
'@  Il est preferable de gerer les depassements éventuels     ',/,&
'@  au moyen des fonctions contenues dans le fichier          ',/,&
'@  cfther.f90: cf_check_internal_energy, cf_check_temperature',/,&
'@  (arret du calcul en fin de pas de temps en cas de         ',/,&
'@   depassement).                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  L''option IPUCOU = ',I10                                   ,/,&
'@    n''est pas compatible avec le module compressible       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Imposer IPUCOU = 0.                                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  Les estimateurs ne sont pas compatibles avec le module    ',/,&
'@    compressible.                                           ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Imposer IESCAL(.) = 0.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 2000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATAS                ',/,&
'@    =========                                               ',/,&
'@    SPECIFIC PHYSICS MODULES (COMPRESSIBLE) SET             ',/,&
'@                                                            ',/,&
'@  The bounds of the variables energy or temperature         ',/,&
'@    have been modified :                                    ',/,&
'@                                                            ',/,&
'@                      SCAMIN        SCAMAX                  ',/,&
'@  energy      ',2E14.5                                       ,/,&
'@  temperature ',2E14.5                                       ,/,&
'@                                                            ',/,&
'@  The bounds of these variables should not be modified.     ',/,&
'@  It is possible to modify the bounds of the variables      ',/,&
'@  density or energy in uscfx2, but it is not recommended.   ',/,&
'@  It is advised to manage the possible overshoot by the     ',/,&
'@  use of the functions defined in the file cfther.f90:      ',/,&
'@  cf_check_internal_energy, cf_check_temperature (stop of   ',/,&
'@  the calculation at the end of the time step in case of an ',/,&
'@  overshoot).                                               ',/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  Check parameters.                                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATAS                ',/,&
'@    =========                                               ',/,&
'@    SPECIFIC PHYSICS MODULES (COMPRESSIBLE) SET             ',/,&
'@                                                            ',/,&
'@  The option IPUCOU = ',I10                                  ,/,&
'@    is not compatible with the compressible module          ',/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  Impose IPUCOU = 0.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATAS                ',/,&
'@    =========                                               ',/,&
'@    SPECIFIC PHYSICS MODULES (COMPRESSIBLE) SET             ',/,&
'@                                                            ',/,&
'@  The error estimators are not compatible with the          ',/,&
'@    compressible module.                                    ',/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  Impose IESCAL(.) = 0.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP WHILE READING INPUT DATAS                ',/,&
'@    =========                                               ',/,&
'@    SPECIFIC PHYSICS MODULES (COMPRESSIBLE) SET             ',/,&
'@                                                            ',/,&
'@    ',A6,' MUST BE AN INTEGER EGAL TO 0 OR 1                ',/,&
'@    IT HAS VALUE',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
