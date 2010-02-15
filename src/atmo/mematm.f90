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

subroutine mematm &
!================

 ( idbia0 , idbra0 , ifinia , ifinra , ra )
!===============================================================================
!  FONCTION
!  --------

!  MEMORY MANAGEMENT FOR THE ATMOSPHERIC SPECIFIC PHYSICS MODULE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ifinia           ! e  ! --> ! pointeur de la premiere cas libre dan          !
!                  !    !     !  dans ia en sortie                             !
! ifinra           ! e  ! --> ! pointeur de la premiere cas libre dan          !
!                  !    !     !  dans ia en sortie                             !
!__________________.____._______________.________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     COMMON DATA
!===============================================================================

include "paramx.h"
include "entsor.h"
include "ppppar.h"
include "atincl.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0, ifinia, ifinra
double precision ra(*)

! Local variables

integer          idebia, idebra, imode
double precision rvoid(1)

!===============================================================================
!  1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
!  2. MEMORY MANAGEMENT FOR THE METEO PROFILES
!===============================================================================

!  2.1 READING THE METEO FILE IN ORDER TO ESTIMATE THE REQUIRED SIZE OF
!  THE TABLES
!-------------------------------------------------------------------------------
if (imeteo.gt.0) then

  imode = 0

!     Nb les arguments ne sont pas utilises quand IMODE=0
  call atlecm &
  !==========
 ( imode      ,             &
   rvoid , rvoid , rvoid ,  &
   rvoid , rvoid , rvoid ,  &
   rvoid , rvoid ,          &
   rvoid , rvoid ,          &
   rvoid , rvoid ,          &
   rvoid , rvoid , rvoid )

  endif

!  2.2 MEMORY MANAGEMENT :
! ------------------------

! --> Integer

ifinia = idebia

! --> Real

itmmet = idebra
izdmet = itmmet + nbmetm
iztmet = izdmet + nbmetd
iumet  = iztmet + nbmett
ivmet  = iumet  + nbmetd*nbmetm
iwmet  = ivmet  + nbmetd*nbmetm
iekmet = iwmet  + nbmetd*nbmetm
iepmet = iekmet + nbmetd*nbmetm
ittmet = iepmet + nbmetd*nbmetm
iqvmet = ittmet + nbmett*nbmetm
ipmer  = iqvmet + nbmett*nbmetm
ixmet  = ipmer  + nbmetm
iymet  = ixmet  + nbmetm
irmet  = iymet  + nbmetm
itpmet = irmet  + nbmett*nbmetm
iphmet = itpmet + nbmett*nbmetm
ifinra = iphmet + nbmett*nbmetm

! --> Verification

CALL IASIZE('MEMATM',IFINIA)
!==========

CALL RASIZE('MEMATM',IFINRA)
!==========

!----
! FIN
!----

return
end subroutine


