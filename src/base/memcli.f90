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

subroutine memcli &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                                                                 &
   isvhb  , isvtb  ,                                              &
   iicodc , ircodc ,                                              &
   iw1    , iw2    , iw3    , iw4    , iw5    , iw6    ,          &
   icoefu , irijip , iuetbo , ivsvdr , ihbord , itbord ,          &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE CONDITIONS LIMITES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! isvhb            ! e  ! <-- ! indicateur de stockage de hbord                !
! isvtb            ! e  ! <-- ! indicateur de stockage de tbord                !
! iicodc,ircodc    ! e  ! --> ! "pointeur" sur icodlc rcodcl                   !
! iw1,2,3,4,5,6    ! e  ! --> ! "pointeur" sur w1,2,3,4,5,6                    !
! iuetbo           ! e  ! --> ! "pointeur" sur uetbor                          !
! ihbord           ! e  ! --> ! "pointeur" sur hbord                           !
! itbord           ! e  ! --> ! "pointeur" sur hbord                           !
! irijip           ! e  ! --> ! "pointeur" sur rijipb                          !
! ivsvdr           ! e  ! --> ! "pointeur" sur visvdr                          !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
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
use optcal
use ppppar
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          isvhb  , isvtb
integer          iicodc , ircodc
integer          iw1    , iw2    , iw3    , iw4    , iw5    , iw6
integer          icoefu , irijip , iuetbo , ivsvdr
integer          ihbord , itbord
integer          ifinia , ifinra

integer          idebia , idebra, irij, iiuetb

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0

irij = 0
if(itytur.eq.3) then
  irij = 1
endif

iiuetb = 0
if(itytur.eq.4.and.idries.eq.1)then
  iiuetb = 1
endif

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

iicodc =       idebia
ifinia =       iicodc + nfabor*nvar

ircodc =       idebra
iw1    =       ircodc + nfabor*nvar*3
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
iw4    =       iw3    + ncelet
iw5    =       iw4    + ncelet
iw6    =       iw5    + ncelet
icoefu =       iw6    + ncelet
irijip =       icoefu + nfabor*ndim
iuetbo =       irijip + nfabor*6*irij
ivsvdr =       iuetbo + nfabor*iiuetb
ifinra =       ivsvdr + ncelet*iiuetb

ihbord =       ifinra
if(isvhb.gt.0) then
  ifinra =       ihbord + nfabor
endif
itbord =       ifinra
if(isvtb.gt.0.or.iirayo.gt.0) then
  ifinra =       itbord + nfabor
endif

!---> VERIFICATION

call iasize('memcli',ifinia)
!==========

call rasize('memcli',ifinra)
!==========

return
end subroutine
