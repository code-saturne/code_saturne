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

subroutine usfui1
!================

!===============================================================================
!  FONCTION  :
!  ---------

!  ROUTINE UTILISATEUR POUR ENTREE DES PARAMETRES DE CALCUL
!  RELATIFS A LA COMBUSTION DU FUEL
!    (COMMONS)

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

implicit none

!===============================================================================
! Common blocks
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
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"
include "ppcpfu.h"

!===============================================================================

integer          jpp , icla

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
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     MODULE COMBUSTION FUEL :                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usfui1 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

!  Sortie chrono, suivi listing, sortie histo
!     Si l'on n'affecte pas les tableaux suivants,
!     les valeurs par defaut seront utilisees

!       ICHRVR( ) = sortie chono (oui 1/non 0)
!       ILISVR( ) = suivi listing (oui 1/non 0)
!       IHISVR( ) = sortie historique (nombre de sondes et numeros)
!       si IHISVR(.,1)  = -1 sortie sur toutes les sondes definies
!                            dans usini1


! --> Variables propres a la suspension gaz - particules

!      - Enthalpie de la suspension
jpp = ipprtp(isca(ihm))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

! --> Variables propres a la phase dispersee

do icla = 1, nclafu
!       - Fraction massique de Fuel
  jpp = ipprtp(isca(iyfol(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1

!       - Nb de particules par kg de melange
  jpp = ipprtp(isca(ing(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1

!       - Enthalpie du Fuel
  jpp = ipprtp(isca(ihlf(icla)))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
enddo


! --> Variables propres a la phase continue

!       - Moyenne de F1 (vapeur )
jpp = ipprtp(isca(ifvap))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - Moyenne du F3 (representatif du C libere sous forme de CO
!       lors de la combustion heterogene)
jpp = ipprtp(isca(ifhtf))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - Variance associe au traceur 4 (air)
jpp = ipprtp(isca(if4p2m))
ichrvr(jpp)  = 1
ilisvr(jpp)  = 1
ihisvr(jpp,1)= -1

!     - YCO2

if ( ieqco2 .ge. 1 ) then
  jpp = ipprtp(isca(iyco2))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
endif

!     - HCN et NO

if ( ieqnox .eq. 1 ) then
  jpp = ipprtp(isca(iyhcn))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
  jpp = ipprtp(isca(iyno))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
  jpp = ipprtp(isca(itaire))
  ichrvr(jpp)  = 1
  ilisvr(jpp)  = 1
  ihisvr(jpp,1)= -1
endif

!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================

!  Sortie chrono, suivi listing, sortie histo
!     Si l'on n'affecte pas les tableaux suivants,
!     les valeurs par defaut seront utilisees

!       ICHRVR( ) = sortie chono (oui 1/non 0)
!       ILISVR( ) = suivi listing (oui 1/non 0)
!       IHISVR( ) = sortie historique (nombre de sondes et numeros)
!       si IHISVR(.,1)  = -1 sortie sur toutes les sondes definies
!                            dans usini1

! --> Variables algebriques propres a la suspension gaz - particules

!     - Masse molaire du melange gazeux
jpp = ipppro(ipproc(immel))
ichrvr(jpp)   = 0
ilisvr(jpp)   = 0
ihisvr(jpp,1) = -1

! --> Variables algebriques propres a la phase dispersee

do icla = 1, nclafu
!       - Temperature des gouttes
  jpp = ipppro(ipproc(itemp3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Masse volumique des gouttes
  jpp = ipppro(ipproc(irom3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Diametre des gouttes
  jpp = ipppro(ipproc(idiam3(icla)))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1

!       - Transfert d'energie par convection-diffusion
  jpp = ipppro(ipproc(ih1hlf(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1

!       - Transfert de masse du a l'evaporation (s-1) < 0
  jpp = ipppro(ipproc(igmeva(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1

!       - Transfert de masse du a la combustion heterogene
  jpp = ipppro(ipproc(igmhtf(icla)))
  ichrvr(jpp)   = 0
  ilisvr(jpp)   = 0
  ihisvr(jpp,1) = -1
enddo

! --> Variables algebriques propres a la phase continue

!     - Temperature du melange gazeux
jpp = ipppro(ipproc(itemp1))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du FOV
jpp = ipppro(ipproc(iym1(ifov)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CO
jpp = ipppro(ipproc(iym1(ico)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du O2
jpp = ipppro(ipproc(iym1(io2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CO2
jpp = ipppro(ipproc(iym1(ico2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du H2O
jpp = ipppro(ipproc(iym1(ih2o)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du N2
jpp = ipppro(ipproc(iym1(in2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du H2S
jpp = ipppro(ipproc(iym1(ih2s)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - Fraction massique (dans le melange gazeux) du SO2
jpp = ipppro(ipproc(iym1(iso2)))
ichrvr(jpp)   = 1
ilisvr(jpp)   = 1
ihisvr(jpp,1) = -1

!     - MODEL NOX :
if ( ieqnox .eq. 1 ) then
  jpp = ipppro(ipproc(ighcn1))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
  jpp = ipppro(ipproc(ighcn2))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
  jpp = ipppro(ipproc(ignoth))
  ichrvr(jpp)   = 1
  ilisvr(jpp)   = 1
  ihisvr(jpp,1) = -1
endif

!===============================================================================
! 3. OPTIONS DE CALCUL
!===============================================================================

! --- Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 4. CONSTANTES PHYSIQUES
!===============================================================================

! ---> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0 = 4.25d-5


!----
! FIN
!----

return

end subroutine
