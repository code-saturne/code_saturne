!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_user_initialization-gas_3ptchem.f90
!>
!> \brief 3 PTCHEM gas example
!>
!> See \subpage cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

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
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, igg
double precision coefg(ngazgm)

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m, cvar_scalt
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

! Control output

write(nfecra,9001)

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if ( isuite.eq.0 ) then

  do iel = 1, ncel

! ----- Mean Mixture Fraction
    cvar_fm(iel)   = fs(1)

! ----- Variance of Mixture Fraction
    cvar_fp2m(iel) = zero

! ----- Enthalpy
    if ( ippmod(icod3p).eq.1 ) then
      cvar_scalt(iel) = hinfue*fs(1)+hinoxy*(1.d0-fs(1))
    endif

  enddo

endif
!< [init]


!--------
! Formats
!--------

 9001 format(                                       /,&
'  user defined initialization of variables'       ,/,&
                                                    /)


!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_f_initialization
