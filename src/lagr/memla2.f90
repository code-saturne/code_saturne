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

subroutine memla2 &
!================

 ( idbia0 , idbra0 ,                                              &
   nfabor , ncelet , nfac   ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   iindep , iibord , iettpa , iauxl  , iauxl2 ,                   &
   itaup  , iitlag , ipiil  ,                                     &
   ivagau , itsuf  , itsup  , ibx    ,                            &
   igradp , igradv , icroul ,                                     &
   itepct , itsfex , itsvar ,                                     &
   icpgd1 , icpgd2 , icpght ,                                     &
   ibrgau , itebru ,                                              &
   iw1    , iw2    , iw3    ,                                     &
   ifinia , ifinra )


!===============================================================================
!  FONCTION
!  --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

! Reservation de la memoire pour les tableaux qui ne doivent pas
! etre conserves en dehors de la boucle en temps.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! tr ! <-- ! pointeur de la premiere cas libre des          !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ncelet           ! e  ! <-- ! nombre d'elements                              !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! iindep           ! e  ! --> ! pointeur sur indep (num cel de depart          !
! iibord           ! e  ! --> ! pointeur sur ibord                             !
! iettpa           ! e  ! --> ! pointeur sur ettpa                             !
! iauxl            ! e  ! --> ! pointeur sur auxl                              !
! iauxl2           ! e  ! --> ! pointeur sur auxl2                             !
! itaup            ! e  ! --> ! pointeur sur taup                              !
! iitlag           ! e  ! --> ! pointeur sur tlag                              !
! ipiil            ! e  ! --> ! pointeur sur piil                              !
! ivagau           ! e  ! --> ! pointeur sur vagaus                            !
! itsuf            ! e  ! --> ! pointeur sur tsuf                              !
! itsup            ! e  ! --> ! pointeur sur tsup                              !
! ibx              ! e  ! --> ! pointeur sur bx                                !
! igradp           ! e  ! --> ! pointeur sur gradpr                            !
! igradv           ! e  ! --> ! pointeur sur gradvf                            !
! icroul           ! e  ! --> ! pointeur sur croule                            !
! itepct           ! e  ! --> ! pointeur sur tempct                            !
! itsfex           ! e  ! --> ! pointeur sur tsfext                            !
! itsvar           ! e  ! --> ! pointeur sur tsvar                             !
! icpgd1           ! e  ! --> ! pointeur sur cpgd1                             !
! icpgd2           ! e  ! --> ! pointeur sur cpgd2                             !
! icpght           ! e  ! --> ! pointeur sur cpght                             !
! iw1              ! e  ! --> ! pointeur sur w1                                !
! iw2              ! e  ! --> ! pointeur sur w2                                !
! iw3              ! e  ! --> ! pointeur sur w3                                !
!                  !    !     !                                                !
!                  ! tr !     !  tableaux ia/ra                                !
! ifinia           ! tr ! --> ! pointeur de la premiere cas libre dan          !
!                  ! tr !     !  dans ia en sortie                             !
! ifinra           ! tr ! --> ! pointeur de la premiere cas libre dan          !
!                  ! tr !     !  dans ia en sortie                             !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use lagpar
use lagran

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nfabor , ncelet , nfac
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          iindep , iibord , iettpa , iauxl , iauxl2
integer          itaup  , iitlag , ipiil
integer          ivagau , itsuf  , itsup  , ibx
integer          igradp , igradv , icroul
integer          itepct , itsfex , itsvar
integer          icpgd1 , icpgd2 , icpght
integer          ibrgau , itebru
integer          iw1     , iw2    , iw3
integer          ifinia , ifinra

! Local variables

integer          idebia , idebra

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA


iibord =       idebia
ifinia =       iibord + nbpmax

iindep = ifinia
ifinia = iindep + nbpmax

iettpa =       idebra
iauxl  =       iettpa + nbpmax * nvp
itaup  =       iauxl  + nbpmax * 3
iitlag =       itaup  + nbpmax
ipiil  =       iitlag + nbpmax * 3
ivagau =       ipiil  + nbpmax * 3
itsuf  =       ivagau + nbpmax * nvgaus
itsup  =       itsuf  + nbpmax * 3
ibx    =       itsup  + nbpmax * 3
itsvar =       ibx    + nbpmax * 3 * 2
igradp =       itsvar + nbpmax * nvp1
iw1    =       igradp + ncelet * 3
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
ifinra =       iw3    + ncelet

if ( (iphyla.eq.1 .and. itpvar.eq.1) .or.                         &
      iphyla.eq.2                         ) then
  itepct = ifinra
  ifinra = itepct + 2*nbpmax
else
  itepct = 1
endif

if (iilagr.eq.2) then
  itsfex = ifinra
  ifinra = itsfex + nbpmax
else
  itsfex = 1
endif

if ( iilagr.eq.2 .and. iphyla.eq.2                                &
                 .and. ltsthe .eq.1 ) then
  icpgd1 = ifinra
  icpgd2 = icpgd1 + nbpmax
  icpght = icpgd2 + nbpmax
  ifinra = icpght + nbpmax
else
  icpgd1 = 1
  icpgd2 = 1
  icpght = 1
endif

if (modcpl.gt.0) then
  igradv =     ifinra
  ifinra =     igradv + ncelet * 9
else
  igradv =  1
endif

if (iroule.eq.1) then
  icroul = ifinra
  ifinra = icroul + ncelet
else
  icroul = 1
endif

if ( lamvbr .eq. 1 ) then
  ibrgau = ifinra
  itebru = ibrgau + nbpmax * nbrgau
  ifinra = itebru + nbpmax
else
  ibrgau = 1
  itebru = 1
endif

if (nordre.eq.2) then
  iauxl2 = ifinra
  ifinra = iauxl2 + nbpmax*7
else
  iauxl2 = 1
endif

!---> VERIFICATION

call iasize('memla2',ifinia)
!==========

call rasize('memla2',ifinra)
!==========

return
end subroutine
