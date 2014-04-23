!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine cs_fuel_thfieldconv2 &
!==============================
 ( ncelet , ncel   ,                  &
   rtp    , propce )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE LA TEMPERATURE DES PARTICULES
!  EN FONCTION DE L'ENTHALPIE DU FOL ET DES CONCENTRATIONS
!
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! eh0              ! tr ! <-- ! tableau reel de travail                        !
! eh1              ! tr ! <-- ! tableau reel de travail                        !
!__________________!____!_____!________________________________________________!
!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================
use paramx
use numvar
use optcal
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
double precision rtp(ncelet,nflown:nvar), propce(ncelet,*)

! Local variables

integer          icel    , icla
integer          ipcte1 , ipcte2
integer          mode

double precision eh2
double precision xsolid(2), mkfini , diamgt
double precision masgut  , mfgout , mkgout , rhofol

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================
! --- Initialisation de T2 a T1

ipcte1 = ipproc(itemp1)
do icla = 1, nclafu
  ipcte2 = ipproc(itemp2(icla))
  do icel = 1, ncel
    propce(icel,ipcte2) = propce(icel,ipcte1)
  enddo
enddo

!===============================================================================
! 2. CALCUL DE LA TEMPERATURE DES PARTICULES
!===============================================================================
!
do icla=1,nclafu
!
  ipcte2 = ipproc(itemp2(icla))

  mkfini = rho0fl*pi/6.d0*dinikf(icla)**3

  do icel = 1, ncel
!
    rhofol = propce(icel,ipproc(irom2 (icla)))
    diamgt = propce(icel,ipproc(idiam2(icla)))
    masgut = rho0fl*pi/6.d0*(diamgt**3.d0)
    if (diamgt.le.dinikf(icla)) then
      mkgout = masgut
    else
      mkgout = mkfini
    endif
    mfgout = masgut - mkgout
    xsolid(1) = 1.d0-fkc
    xsolid(2) = fkc
    if(masgut.gt.zero) then
      xsolid(1) = mfgout / masgut
      xsolid(2) = mkgout / masgut
    endif
    xsolid(1) = min(1.d0,max(0.d0,xsolid(1)))
    xsolid(2) = min(1.d0,max(0.d0,xsolid(2)))
!
    if ( rtp(icel,isca(iyfol(icla))) .gt. (3.d3*epsifl) ) then
!
      eh2 =  rtp(icel,isca(ih2(icla)))/rtp(icel,isca(iyfol(icla)))
!
      mode = 1
      call cs_fuel_htconvers2 &
!     =======================
      (mode, eh2 , xsolid , propce(icel,ipcte2))
!
    endif

  enddo
!
enddo

!----
! END
!----
return
end subroutine
