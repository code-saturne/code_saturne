!-------------------------------------------------------------------------------

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

subroutine ppcsca
!================


!===============================================================================
!  FONCTION  :
!  ---------

!       CALCUL DE NSCAPP SUIVANT LA PHYSIQUE PARTICULIERE
!       DEMANDEE PAR l'UTILISATEUR

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
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "elincl.h"
include "fuincl.h"
include "ppcpfu.h"
include "atincl.h"

!===============================================================================

! Local variables

integer          iok


!===============================================================================

!===============================================================================
! 1. MODELES DE COMBUSTION GAZ (ICO...)
!===============================================================================

! --> Flamme de diffusion

if ( ippmod(icod3p).eq.0 ) nscapp = 2
if ( ippmod(icod3p).eq.1 ) nscapp = 3

if ( ippmod(icodeq).eq.0 ) nscapp = 2
if ( ippmod(icodeq).eq.1 ) nscapp = 3

! --> Flamme de premelange (modele EBU)

if ( ippmod(icoebu).eq.0 ) nscapp = 1
if ( ippmod(icoebu).eq.1 ) nscapp = 2
if ( ippmod(icoebu).eq.2 ) nscapp = 2
if ( ippmod(icoebu).eq.3 ) nscapp = 3

! ---> Modele BML

if ( ippmod(icobml).eq.0 ) nscapp = 2
if ( ippmod(icobml).eq.1 ) nscapp = 3

! ---> Modele LWC

if ( ippmod(icolwc).eq.0 ) nscapp = 4
if ( ippmod(icolwc).eq.1 ) nscapp = 5
if ( ippmod(icolwc).eq.2 ) nscapp = 5
if ( ippmod(icolwc).eq.3 ) nscapp = 6
if ( ippmod(icolwc).eq.4 ) nscapp = 5
if ( ippmod(icolwc).eq.5 ) nscapp = 6

! ---nombre de Dirac pour le modele LWC

if ( ippmod(icolwc).eq.0 .or.                                     &
     ippmod(icolwc).eq.1           ) then
  ndirac = 2

else if ( ippmod(icolwc).eq.2 .or.                                &
          ippmod(icolwc).eq.3      ) then
  ndirac = 3

else if ( ippmod(icolwc).eq.4 .or.                                &
          ippmod(icolwc).eq.5      ) then
  ndirac = 4
endif

!===============================================================================
! 2. MODELES DE COMBUSTION CHARBON PULVERISE (ICP...)
!===============================================================================

! --> Flamme charbon pulverise

! ------ Transport d'H2
if ( ippmod(icp3pl).eq.0 ) then

  nscapp = 1 + 2*ncharb + 2 + 4*nclacp + (noxyd-1)

  if ( ihtco2.eq. 1) nscapp = nscapp + 1
  if ( ieqco2.ge. 1) nscapp = nscapp + 1

endif

! ------ Transport d'H2 + phase de sechage

if ( ippmod(icp3pl).eq.1 ) then

  nscapp = 1 + 2*ncharb + 3 +5*nclacp + (noxyd-1)

  if ( ihtco2.eq. 1) nscapp = nscapp + 1
  if ( ieqco2.eq. 1) nscapp = nscapp + 1

endif

! --> Flamme charbon pulverise couple Lagrangien

if ( ippmod(icpl3c).eq.0 ) then
  nscapp = 1 + 2*ncharb + 2
endif

!===============================================================================
! 3. MODELE COMPRESSIBLE SANS CHOC : rho, e, s
!===============================================================================

if ( ippmod(icompf).ge.0 ) nscapp = 3*nphas


!===============================================================================
! 4. MODELES ELECTRIQUES : Effet Joule        (IELJOU)
!                          Arc electrique     (IELARC)
!                          Conduction ionique (IELION)
!===============================================================================


! --> Effet Joule

  if ( ippmod(ieljou).eq.1 ) nscapp = 2 + ngazg-1
  if ( ippmod(ieljou).eq.2 ) nscapp = 3 + ngazg-1

  if ( ippmod(ieljou).eq.3 ) nscapp = 2 + ngazg-1
  if ( ippmod(ieljou).eq.4 ) nscapp = 3 + ngazg-1

! --> Arc electrique

  if ( ippmod(ielarc).eq.1 ) nscapp = 2 + ngazg-1
  if ( ippmod(ielarc).eq.2 ) nscapp = 5 + ngazg-1

! --> Conduction ionique

  if ( ippmod(ielion).eq.1 ) nscapp = 2 + ngazg-1

!===============================================================================
! 5. MODELES DE COMBUSTION FUEL (ICFUEL)
!===============================================================================

if ( ippmod(icfuel).ge.0 ) then
  nscapp = 4 + 3*nclafu
  if ( ieqco2.ge. 1) nscapp = nscapp + 1
  if ( ieqnox.eq. 1) nscapp = nscapp + 3
endif

!===============================================================================
! 6. MODELE ATMOSPHERIQUE (IATMOS)
!===============================================================================

if ( ippmod(iatmos).eq.0 ) nscapp = 0
if ( ippmod(iatmos).eq.1 ) nscapp = 1
if ( ippmod(iatmos).eq.2 ) nscapp = 3

!===============================================================================
! 8. MODELISATION DES AEROREFRIGERANTS (IAEROS)
!===============================================================================

if ( ippmod(iaeros).ge.0 ) nscapp = 2

!===============================================================================
! 9. VERIFICATION : UNE SEULE PHASE
!===============================================================================

iok = 0

if( nscapp.gt.0 .and. ippmod(icompf).lt.0 ) then
!                     ^^^^^^^^^^^^^^^^^^^^^^^^^
!     On peut penser faire du multiphasique avec la physique particuliere
!     compressible sans choc... (mettons deux pressions indépendantes
!     par exemple) ... a mettre en place si necessaire ...

  if(nphas.ne.1) then
    write(nfecra,1000)
    iok = iok + 1
  endif
endif

if(iok.ne.0) then
  call csexit (1)
  !==========
endif

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SEULE UNE PHASE PERMISE EN PHYSIQUE PARTICULIERE        ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Arret dans ppcsca.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return
end subroutine
