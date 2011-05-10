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

subroutine memra2 &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   idtr   , iviscf , iviscb , idam   , ixam   ,                   &
   idrtp  , ismbr  , irovsd ,                                     &

   icorua , icorub , iflxma , iflxmb , itek   ,                   &

   iw1    , iw2    , iw3    , iw4    , iw5    , iw6    , iw7    , &
   iw8    , iw9    , iw10   ,                                     &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  GESTION MEMOIRE POUR LE CALCUL DES TERMES SOURCES RADIATIFS
!  ET DU FLUX INCIDENT RADIATIF PARIETAL (RAYDOM)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! idtr             ! e  ! --> ! "pointeur" sur dtr                             !
! iviscf, b        ! e  ! --> ! "pointeur" sur viscf, viscb                    !
! idam, ixam       ! e  ! --> ! "pointeur" sur dam, xam                        !
! idrtp            ! e  ! --> ! "pointeur" sur drtp                            !
! ismbr            ! e  ! --> ! "pointeur" sur smbr                            !
! irovsd           ! e  ! --> ! "pointeur" sur rovsdt                          !
! icorua,icorub    ! e  ! --> ! "pointeur" sur cofrua et cofrub                !
! iflxma,iflxmb    ! e  ! --> ! "pointeur" sur flumas et flumab                !
! itek             ! e  ! --> ! "pointeur" sur tempk                           !
! iw1,2,...,11     ! e  ! --> ! "pointeur" sur w1 a w10                       !
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
use entsor
use numvar
use optcal
use ppppar
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          idtr
integer          iviscf , iviscb , idam   , ixam
integer          idrtp  , ismbr  , irovsd
integer          iw1    , iw2    , iw3    , iw4    , iw5    , iw6
integer          iw7    , iw8    , iw9    , iw10
integer          ifinia , ifinra

integer          icorua , icorub
integer          iflxma , iflxmb
integer          itek

! Local variables

integer          idebia, idebra

!===============================================================================
!===============================================================================
!---> INITIALISATION
!===============================================================================
idebia = idbia0
idebra = idbra0

!===============================================================================
!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
!===============================================================================

idtr   =       idebra
iviscf =       idtr   + ncelet
iviscb =       iviscf + nfac
idam   =       iviscb + nfabor
ixam   =       idam   + ncelet
idrtp  =       ixam   + nfac*2
ismbr  =       idrtp  + ncelet
irovsd =       ismbr  + ncelet
iw1    =       irovsd + ncelet
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
iw4    =       iw3    + ncelet
iw5    =       iw4    + ncelet
iw6    =       iw5    + ncelet
iw7    =       iw6    + ncelet
iw8    =       iw7    + ncelet
iw9    =       iw8    + ncelet
iw10   =       iw9    + ncelet
itek   =       iw10   + ncelet
ifinra =       itek   + ncelet * nrphas

icorua =       ifinra
icorub =       icorua + nfabor
iflxma =       icorub + nfabor
iflxmb =       iflxma + nfac
ifinra =       iflxmb + nfabor

!===============================================================================
!---> VERIFICATION
!===============================================================================

call rasize('memra2',ifinra)
!==========

return

end subroutine
