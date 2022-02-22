!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ii, iok

type(var_cal_opt) :: vcopt

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! Does scalar itempk behave like a temperature ?
! TODO check this; should be 1 for temperature unless handled in
!      another manner

call field_set_key_int(ivarfl(isca(itempk)), kscacp, 0)

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

call cscfgp(icfgrp)

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

return
end subroutine
