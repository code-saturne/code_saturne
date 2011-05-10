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

subroutine memcft &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   iw7    , iw8    , iw9    , iw10   , iw11   , iw12   ,          &
   iviscf , icoefu , ixam   ,                                     &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!             GESTION MEMOIRE DANS CFDTTV POUR MEMCFT

!  (APPELE PAR DTTVAR EN COMPRESSIBLE SANS CHOC)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iw7..12          ! e  ! --> ! "pointeurs" sur w7 a w12                       !
! iviscf           ! e  ! --> ! "pointeur" sur viscf                           !
! icoefu           ! e  ! --> ! "pointeur" sur coefu                           !
! ixam             ! e  ! --> ! "pointeur" sur xam                             !
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

use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          iw7    , iw8    , iw9    , iw10   , iw11   , iw12
integer          iviscf , icoefu , ixam
integer          ifinia , ifinra


! Local variables

integer          idebia, idebra

!===============================================================================


!---> INITIALISATION

idebia = idbia0
idebra = idbra0

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia =       idebia

iw7    =       idebra
iw8    =       iw7    + ncelet
iw9    =       iw8    + ncelet
iw10   =       iw9    + ncelet
iw11   =       iw10   + ncelet
iw12   =       iw11   + ncelet
iviscf =       iw12   + ncelet
icoefu =       iviscf + nfac
ixam   =       icoefu + nfabor*3
ifinra =       ixam   + nfac*2

!---> VERIFICATION

call iasize('memcft',ifinia)
!==========

call rasize('memcft',ifinra)
!==========

return
end subroutine
