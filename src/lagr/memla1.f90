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

subroutine memla1 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   lndnod ,                                                       &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   iiitep , iicoce , iityce ,                                     &
   iettp  , iettpa , iitepa , istatc , istatv ,                   &
   itslag , istatf ,                                              &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!         Reservation de la memoire pour les tableaux qui doivent
!         etre conserves en dehors de la boucle en temps, notamment
!         pour le post-processing.

!         Remarque : tous les tableaux ouverts ici doivent
!                    recevoir une initialisation de valeurs
!                    par defaut dans le sous-programme LAGLEC

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ia en entree                             !
! idbra0           ! e  ! <-- ! pointeur de la premiere cas libre              !
!                  !    !     !  dans ra en entree                             !
! ndim             ! e  ! <-- ! dimension (3)                                  !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! iiitep           ! e  ! --> ! pointeur sur itepa                             !
! iicoce           ! e  ! --> ! pointeur sur icocel                            !
! iityce           ! e  ! --> ! pointeur sur itycel                            !
! iettp            ! e  ! --> ! pointeur sur ettp                              !
! iettpa           ! e  ! --> ! pointeur sur ettpa (simple init. ici)          !
! iitepa           ! e  ! --> ! pointeur sur tepa                              !
! istatc           ! e  ! --> ! pointeur sur statis                            !
! itslag           ! e  ! --> ! pointeur sur tslagr                            !
! istatf           ! e  ! --> ! pointeur sur parbor                            !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
!__________________.____._____.________________________________________________.

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
include "entsor.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          lndnod
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          iiitep , iicoce , iityce
integer          iettp  , iettpa
integer          iitepa , istatc , istatv , itslag , istatf
integer          ifinia , ifinra

! Local variables

integer          idebia , idebra

!===============================================================================

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA


!     IETTPA est simplement initialise ici, pour eviter un arret avec les
!       options de compilation check bounds. La reservation memoire
!       correspondante est faite dans memla2.

idebia = idbia0
idebra = idbra0


if (iilagr.eq.0) then

  iiitep =       idebia
  iicoce =       iiitep
  iityce =       iicoce
  iifrla =       iityce
  ifinia =       iifrla

  iettpa =       idebra
  iettp  =       iettpa
  iitepa =       iettp
  istatc =       iitepa
  istatv =       istatc
  itslag =       istatv
  istatf =       itslag
  ifinra =       istatf

else

  iiitep =       idebia
  iicoce =       iiitep + nbpmax * nivep
  iityce =       iicoce + lndnod
  iifrla =       iityce + ncelet + 1
  ifinia =       iifrla + nfabor

  iettpa =       idebra
  iettp  =       iettpa
  iitepa =       iettp  + nbpmax * nvp
  istatc =       iitepa + nbpmax * nvep
  istatv =       istatc + ncelet * nvlsta                         &
                        + ncelet * nvlsta * nbclst
  itslag =       istatv + ncelet * max((nvlsta-1),0)              &
                        + ncelet * max((nvlsta-1),0) * nbclst
  istatf =       itslag + ntersl * ncelet
  ifinra =       istatf + nfabor * nvisbr

endif

!---> VERIFICATION

CALL IASIZE('MEMLA1',IFINIA)
!==========

CALL RASIZE('MEMLA1',IFINRA)
!==========

return
end subroutine
