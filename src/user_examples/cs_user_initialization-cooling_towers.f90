!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization-cooling_towers.f90
!> \brief Cooling towers example
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________


subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          iel, iutile
integer          ilelt, nlelt

double precision d2s3

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

d2s3 = 2.d0/3.d0

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then

!   --- Initialisation de la temperature de l'air a 11 deg Celsius
!                      de l'humidite de l'air a 0.0063
!       pour toutes les cellules

  do iel = 1, ncel

    rtp(iel,isca(itemp4)) = 11.d0
    rtp(iel,isca(ihumid)) = 0.0063d0

  enddo

!   --- Initialisation de la temperature de l'air a 20 deg Celsius
!                      de l'humidite de l'air a 0.012
!                      de la vitesse
!       uniquement pour les cellules de couleur 6

  call getcel('6', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    iel = lstelt(ilelt)

    rtp(iel,iu) = -0.5d0

    rtp(iel,isca(itemp4)) = 20.d0
    rtp(iel,isca(ihumid)) = 0.012d0

  enddo

endif
!< [init]

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_initialization
