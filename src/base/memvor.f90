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

subroutine memvor &
!================

 ( idbia0 , idbra0 , iappel , nfabor , ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE POUR LA METHODE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
!__________________.____._______________.________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "entsor.f90"
include "vorinc.f90"

!===============================================================================

! Arguments

integer          idbia0, idbra0, ifinia, ifinra
integer          iappel, nfabor

! Local variables

integer          idebia, idebra

idebia = idbia0
idebra = idbra0

if(iappel.eq.1) then

  iirepv = idebia
  ifinia = iirepv + nfabor

  ifinra = idebra

  CALL IASIZE('MEMVOR',IFINIA)

elseif(iappel.eq.2) then

  iifagl = idebia
  ifinia = iifagl + nnent*icvmax

  ixyzv  = idebra
  ivisv  = ixyzv  + nnent*icvmax*3
  iw1x   = ivisv  + nnent*icvmax
  iw1y   = iw1x   + nnent*icvmax
  iw1z   = iw1y   + nnent*icvmax
  iw1v   = iw1z   + nnent*icvmax
  iw2x   = iw1v   + nnent*icvmax
  iw2y   = iw2x   + nnent*icvmax
  iw2z   = iw2y   + nnent*icvmax
  iw2v   = iw2z   + nnent*icvmax
  ifinra = iw2v   + nnent*icvmax

  CALL IASIZE('MEMVOR',IFINIA)
  CALL RASIZE('MEMVOR',IFINRA)

elseif(iappel.eq.3) then

  iivrce = idebia
  ifinia = iivrce + nnent*nvomax

  iyzcel = idebra
  iuvort = iyzcel + nnent*icvmax*2
  ivvort = iuvort + nnent*icvmax
  iwvort = ivvort + nnent*icvmax
  iyzvor = iwvort + nnent*icvmax
  iyzvoa = iyzvor + nnent*nvomax*2
  isignv = iyzvoa + nnent*nvomax*2
  ixsigm = isignv + nnent*nvomax
  ixgamm = ixsigm + nnent*nvomax
  ixtmp  = ixgamm + nnent*nvomax*2
  ixtmpl = ixtmp  + nnent*nvomax
  ifinra = ixtmpl + nnent*nvomax

  CALL IASIZE('MEMVOR',IFINIA)
  CALL RASIZE('MEMVOR',IFINRA)
endif

return
end subroutine
