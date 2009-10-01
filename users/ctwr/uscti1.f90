!-------------------------------------------------------------------------------

!VERS


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

subroutine uscti1
!================


!===============================================================================
!  FONCTION  :
!  ---------


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
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "ctincl.h"

!===============================================================================

integer          iphas

double precision cpa,cpe,cpv,hv0,rhoe,visc,conduc

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@                      MODULE AEROREFRIGERANTS               ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usctin DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       options generales. Il est indispensable.             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  PARAMETRES POUR L'ECART DE TEMPERATURE IMPOSE
!===============================================================================

!     ACTIVATION
iaeeri = 0

!     ECART DE REFRIGERATION A IMPOSER
vaeeri = 13.d0

!     FREQUENCE DE MODIFICATION DE LA TEMPERATURE
iaeerp = 5

!     PAS DE TEMPERATURE POUR LE CALCUL DE LA PENTE DE ECARTREF(TEAU)
paseri = 0.015d0

!     MAXIMUM DE LA TEMPERATURE D'EAU CHAUDE MOYENNE PONDEREE
aetemx = 80.d0

!     MINIMUM DE LA TEMPERATURE D'EAU REFROIDIE MOYENNE PONDEREE
aetemn = 10.d0

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE ENTREE EAU

nbzsup = 2

!     LISTE DES NBZSUP ZONES D'ECHANGES EN BORD DE L'ENTREE EAU

lizsup(1) = 1
lizsup(2) = 2

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE SORTIE EAU
nbzinf = 2

!     LISTE DES NBZINF ZONES D'ECHANGES EN BORD DE LA SORTIE EAU

lizinf(1) = 1
lizinf(2) = 2

!     INSTANT ACTIVATION ECART IMPOSE

inbaei = 1000.D0

!===============================================================================
! 2.  POST-PROCESSING DES ZONES D'ECHANGES
!===============================================================================

ichrze = 1

!===============================================================================
! 3.  CALCUL SUITE AEROREFRIGERANT
!===============================================================================

isuict = isuite

!===============================================================================
! 4.  PROPRIETES DE L'AIR
!===============================================================================

! Il est deconseille de modifier ici ces proprietes

cpa    = 1006.0d0
cpv    = 1831.0d0
cpe    = 4179.0d0
hv0    = 2501600.0d0
rhoe   = 997.85615d0
visc   = 1.765d-5
conduc = 0.02493d0

call ctprof &
!==========
( cpa, cpv, cpe, hv0, rhoe, visc, conduc, gx, gy, gz )

!----
! FIN
!----

return
end subroutine
