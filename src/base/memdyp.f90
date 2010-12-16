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

subroutine memdyp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   idam   , ixam   , ismbr  , irovsd ,                            &
   irtdp  , idrtdp ,                                              &
   iqfx   , iqfy   , iqfz   , icoefq , irho   , irhob  ,          &
   iflua  , iflub  ,                                              &
   icoax  , icobx  , icoay  , icoby  , icoaz  , icobz  ,          &
   iw1    , iw2    , iw3    , iw4    , iw5    , iw6    , iw7    , &
   iw8    , iw9    ,                                              &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE DISTANCE ADIMENSIONNELLE A LA PAROI

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! idam, ixam       ! e  ! --> ! "pointeur" sur dam, xam                        !
! ismbr            ! e  ! --> ! "pointeur" sur smbr                            !
! irovsd           ! e  ! --> ! "pointeur" sur rovsdt                          !
! idrtp            ! e  ! --> ! "pointeur" sur drtp                            !
! idrtdp           ! e  ! --> ! "pointeur" sur drtdp                           !
! iqx, y, z        ! e  ! --> ! "pointeur" sur qx, qy, qz                      !
! icoefq           ! e  ! --> ! "pointeur" sur coefq                           !
! irho, irhob      ! e  ! --> ! "pointeur" sur rho, rho bord                   !
! iflua, iflub     ! e  ! --> ! "pointeur" sur flumas, flumab                  !
! icoax,y,z, b     ! e  ! --> ! "pointeur" sur coefax,y,z, coefbx,y,z          !
! iw1,2,...,9      ! e  ! --> ! "pointeur" sur w1 a w9                         !
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
use numvar
use optcal
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          idam   , ixam   , ismbr  , irovsd
integer          irtdp  , idrtdp
integer          iqfx   , iqfy   , iqfz   , icoefq
integer          irho   , irhob
integer          iflua  , iflub
integer          icoax  , icobx
integer          icoay  , icoby
integer          icoaz  , icobz
integer          iw1    , iw2    , iw3    , iw4    , iw5    , iw6
integer          iw7    , iw8    , iw9
integer          ifinia , ifinra

integer          idebia, idebra

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0


!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia =       idebia

idam   =       idebra
ixam   =       idam   + ncelet
ismbr  =       ixam   + nfac*2
irovsd =       ismbr  + ncelet
irtdp  =       irovsd + ncelet
idrtdp =       irtdp  + ncelet
iqfx   =       idrtdp + ncelet
iqfy   =       iqfx   + ncelet
iqfz   =       iqfy   + ncelet
icoefq =       iqfz   + ncelet
irho   =       icoefq + nfabor*ndim
irhob  =       irho   + ncelet
iflua  =       irhob  + nfabor
iflub  =       iflua  + nfac
icoax  =       iflub  + nfabor
icobx  =       icoax  + nfabor
icoay  =       icobx  + nfabor
icoby  =       icoay  + nfabor
icoaz  =       icoby  + nfabor
icobz  =       icoaz  + nfabor
iw1    =       icobz  + nfabor
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
iw4    =       iw3    + ncelet
iw5    =       iw4    + ncelet
iw6    =       iw5    + ncelet
iw7    =       iw6    + ncelet
iw8    =       iw7    + ncelet
iw9    =       iw8    + ncelet
ifinra =       iw9    + ncelet

!---> VERIFICATION

call rasize('memdyp',ifinra)
!==========

return
end subroutine
