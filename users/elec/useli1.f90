!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

                  subroutine useli1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!              LE MODULE ELECTRIQUE
!   EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1

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
include "ppincl.h"
include "elincl.h"

!===============================================================================

integer          ipp, iesp , idimve

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
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR useli1 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       options generales. Il est indispensable.             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

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


! --> Variables communes aux versions electriques

! ---- Enthalpie
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Potentiel reel
ipp = ipprtp(isca(ipotr))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!---- Fractions massiques des constituants
if ( ngazg .gt. 1 ) then
  do iesp = 1, ngazg-1
    ipp = ipprtp(isca(iycoel(iesp)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Version effet Joule

if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ipp = ipprtp(isca(ipoti))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! --> Version arc electrique 3D

if ( ippmod(ielarc).ge.2 ) then
  do idimve = 1, ndimve
    ipp = ipprtp(isca(ipotva(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Version conduction ionique
!     indisponible dans la version presente

!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp) )
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Conductivite electrique
ipp = ipppro(ipproc(ivisls(ipotr)))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Puissance volumique dissipee par effet Joule
ipp = ipppro(ipproc(iefjou) )
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- densite de courant reelle
do idimve = 1, ndimve
  ipp = ipppro(ipproc(idjr(idimve)) )
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo

! ---- densite de courant imaginaire
if ( ippmod(ieljou).eq.4 ) then
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(idji(idimve)) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

if ( ippmod(ielarc).ge.1 ) then

! ---- Forces electromagnetiques de Laplace
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(ilapla(idimve)) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo

! ---- Coefficient d'absorption ou TS radiatif
  if ( ixkabe.gt.0 ) then
    ipp = ipppro(ipproc(idrad) )
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
endif

! ---- Charge electrique volumique
if ( ippmod(ielion).ge.1 ) then
  ipp = ipppro(ipproc(iqelec) )
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! 3. OPTIONS DE CALCUL
!===============================================================================

! --> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.d0

! --> Recalage des variables electriques (joule ou arc electrique)
!      IELCOR = 0 : pas de correction
!      IELCOR = 1 : correction
ielcor = 0

!     Intensite de courant imposee (arc electrique) en Ampere
!             et Puissance imposee (effet Joule/verre) en Watt
!       ces valeurs doivent etre positives
!       (et en general elles sont strictement positives)
couimp = 0.d0
puisim = 0.d0

!     Differentiel de potentiel initiale
!       la valeur doit etre strictement positive
dpot = 0.d0


!----
! FIN
!----

return
end
