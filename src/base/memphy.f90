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

subroutine memphy &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , ncelet , ncel   , nfac   , nfabor , nphas  ,          &
   iw1    , iw2    , iw3    , iw4    ,                            &
   iw5    , iw6    , iw7    , iw8    ,                            &
   iw9    , iw10   , iw11   , iw12   , ixmij  ,                   &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE CALCUL DES PROPRIETES PHYSIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nvar             ! e  ! <-- ! nombre de variables                            !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! iw1,2,3,4,5,6    ! e  ! --> ! "pointeur" sur w1,2,3,4,5,6,7,8,9,10,11,12     !
! ixmij            ! e  ! --> ! "pointeur" sur xmij                            !
! ifinia           ! e  ! --> ! pointeur de la premiere cas libre dan          !
!                  !    !     !  dans ia en sortie                             !
! ifinra           ! e  ! --> ! pointeur de la premiere cas libre dan          !
!                  !    !     !  dans ia en sortie                             !
!__________________.____._____.________________________________________________.

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
include "optcal.h"

!===============================================================================


integer          idbia0 ,idbra0
integer          nvar
integer          ncelet , ncel   , nfac   , nfabor, nphas
integer          iw1    , iw2    , iw3    , iw4
integer          iw5    , iw6    , iw7    , iw8
integer          iw9    , iw10   , iw11   , iw12  , ixmij
integer          ifinia , ifinra

integer          idebia , idebra
integer          imemph , iphas

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0


!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia =       idebia

iw1    =       idebra
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
iw4    =       iw3    + ncelet
iw5    =       iw4    + ncelet
iw6    =       iw5    + ncelet
iw7    =       iw6    + ncelet
iw8    =       iw7    + ncelet
iw9    =       iw8    + ncelet

imemph = 0
do iphas =1, nphas
  if(iturb(iphas).eq.41)  imemph = 1
enddo
do iphas =1, nphas
  if(iturb(iphas).eq.42)  imemph = 2
enddo

if(imemph.eq.1) then
  iw10   = iw9  + ncelet
  iw11   = iw10  ! w11 is not used
  iw12   = iw11  ! w12 is not used
  ixmij  = iw12 + ncelet
  ifinra = ixmij+ 6*ncelet
else if(imemph.eq.2) then
  iw10   = iw9  + ncelet
  iw11   = iw10 + ncelet
  iw12   = iw11 + ncelet
  ixmij  = iw12 + ncelet
  ifinra = ixmij ! xmij is not used
else
  iw10   = iw9
  iw11   = iw10
  iw12   = iw11
  ixmij  = iw12
  ifinra = ixmij
endif


!---> VERIFICATION

CALL IASIZE('MEMPHY',IFINIA)
!==========

CALL RASIZE('MEMPHY',IFINRA)
!==========

return
end subroutine
