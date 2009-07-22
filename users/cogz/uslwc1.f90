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

                  subroutine uslwc1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!                  POUR LA COMBUSTION
!           FLAMME DE PREMELANGE : MODELE LWC
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
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "radiat.h"

!===============================================================================

integer          ipp, idirac

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
'@     MODULE COMBUSTION GAZ MODELE LWC :                     ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslwc1 DOIT ETRE COMPLETE',/,&
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


! ---- Fraction de melange
if ( ippmod(icolwc).ge.0 ) then
  ipp = ipprtp(isca(ifm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance de la fraction de melange
  ipp = ipprtp(isca(ifp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Fraction massique de fuel
  ipp = ipprtp(isca(iyfm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance de la fraction massique de fuel
  ipp = ipprtp(isca(iyfp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if (ippmod(icolwc).ge.2) then
    ipp = ipprtp(isca(icoyfp))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
endif

! ---- Enthalpie
if ( ippmod(icolwc).eq.1 .or.                                     &
     ippmod(icolwc).eq.3 .or.                                     &
     ippmod(icolwc).eq.5    ) then
  ipp = ipprtp(isca(ihm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================

! --- Terme source
  ipp = ipppro(ipproc(itsc))
  NOMVAR(IPP)   = 'T.SOURCE'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- temperature
  ipp = ipppro(ipproc(itemp))
  NOMVAR(IPP)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- fraction massique Combustible
  ipp = ipppro(ipproc(iym(1)))
  NOMVAR(IPP)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- fraction massique Oxydant
  ipp = ipppro(ipproc(iym(2)))
  NOMVAR(IPP)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- fraction massique Produit
  ipp = ipppro(ipproc(iym(3)))
  NOMVAR(IPP)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  do idirac = 1, ndirac
    ipp = ipppro(ipproc(irhol(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'RHOL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iteml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TEML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmel(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMEL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmal(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'FMAL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iampl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'AMPL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(itscl(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'TSCL',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(imaml(idirac)))
    WRITE(NOMVAR(IPP),'(A4,I1)') 'MAML',IDIRAC
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo

! ---- Modele LWC AVEC RAYONNEMENT

if ( iirayo.gt.0 ) then

! ---- Coeff d'absorption
  ipp = ipppro(ipproc(ickabs))
  NOMVAR(IPP)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Terme T^4
  ipp = ipppro(ipproc(it4m))
  NOMVAR(IPP)   = 'TEMP4'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Terme T^3
  ipp = ipppro(ipproc(it3m))
  NOMVAR(IPP)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 3. OPTIONS DE CALCUL
!===============================================================================

! --> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.95d0


!===============================================================================
! 4. CONSTANTES PHYSIQUES
!===============================================================================

! --> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0 = 4.25d-5

! --> Constante du modele LWC

! --- Vitesse de reference
 vref = 60.d0
! --- Longueur de reference
 lref = 0.1d0
! --- Temperature d'activation
 ta   = 0.2d5
! --- Temperature de cross-over (combustion du propane)
 tstar= 0.12d4

!----
! FIN
!----

return
end
