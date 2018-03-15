!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------
!> \file ctiniv.f90
!> \brief Initialisation of calculation variables for the cooling tower module,
!> it is the counterpart of usiniv.f90.
!>
!> Initialise for example the meteorological field for each cell of
!> the domain by interpolation of the data from the meteo file
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   nvar        total number of variables
!> \param[in]   nscal       total number of scalars
!> \param[in]   dt          time step value
!-------------------------------------------------------------------------------

subroutine ctiniv(nvar   , nscal  , dt)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ctincl
use ppincl
use field
use mesh
use cs_c_bindings

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel, ifac, f_id
integer          iflmas, iflmab

double precision, dimension(:), pointer :: cvar_temp, cvar_templ, cvar_yml
double precision, dimension(:), pointer :: cvar_ymw
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_val_s(ivarfl(isca(iscalt)), cvar_temp)
call field_get_val_s(ivarfl(isca(iyml)), cvar_yml)
call field_get_val_s(ivarfl(isca(iymw)), cvar_ymw)
call field_get_val_s(itml, cvar_templ)

!===============================================================================
! 2. Standard initialization
!===============================================================================

! Only if the simulation is not a restart from another one
if (isuite.eq.0) then
  do iel = 1, ncel
    cvar_temp(iel) = (t0 - tkelvi)
    cvar_ymw(iel) = humidity0 / ( 1.d0 + humidity0)
    ! The liquid values can be adjusted in the packing regions
    ! using 'cs_user_f_initialization'
    cvar_templ(iel) = cvar_temp(iel)
    cvar_yml(iel) = 0.0

  enddo

  call synsca(cvar_temp)
  call synsca(cvar_ymw)
  call synsca(cvar_templ)
  call synsca(cvar_yml)

  ! Diffusivities of the dry air and the injected liquid
  ! Note: this comes after 'cs_user_cooling_towers' so it will overwrite
  !       what users may have specified there
  visls0(iymw) = 1.d-12
  visls0(iyml) = 1.d-12

  ! initialise:
  !   - the enthalpies, which are the solution variables
  !   - the humidity, which users might have modified if they changed the
  !     mass fraction of dry air in the humid air

  call cs_ctwr_init_field_vars(ro0,t0,p0,molmass_rat)

  ! Reference diffusivity of the injected liquid enthalpy
  ! The diffusivity of the injected liquid enthalpy is variable
  ! and, therefore, updated at each time step in 'ctphyv')
  if (cp_l.le.0.0 .or. lambda_l.le.0.0) then
    !!FIXME - stop the code and publish an error message
  else
    visls0(ihml) = lambda_l/cp_l
  endif

else

   !! Add - NAT
   !! Restarts

   ! Diffusivities of the dry air and the injected liquid
   visls0(iymw) = 1.d-12
   visls0(iyml) = 1.d-12

   !! Restarts - recompute the required properties based on the
   !! saved solution variables: for example, the humidty, liquid
   !! vertical velocity, etc.
   call cs_ctwr_restart_field_vars(ro0,t0,p0,humidity0,molmass_rat)

endif

!===============================================================================
! 3. User initialization
!===============================================================================

call cs_user_f_initialization &
( nvar   , nscal  ,                                            &
  dt     )

!===============================================================================
! 5. Imposed injected liquid mass flux in packing zones
!===============================================================================

f_id = ivarfl(isca(iyml))

! Id of the mass flux
call field_get_key_int(f_id, kimasf, iflmas) ! interior mass flux
! Pointer to the internal mass flux
call field_get_val_s(iflmas, imasfl)

! Id of the mass flux
call field_get_key_int(f_id, kbmasf, iflmab) ! boundary mass flux
! Pointer to the Boundary mass flux
call field_get_val_s(iflmab, bmasfl)

! Initialise the liquid mass flow through the cell faces
! and a band of ghost cells on each side of the packing zone in order
! to impose boundary values
call cs_ctwr_init_flow_vars(imasfl)

call synsca(cvar_temp)
call synsca(cvar_ymw)
call synsca(cvar_templ)
call synsca(cvar_yml)

do ifac = 1, nfabor
  bmasfl(ifac) = 0.d0
enddo


!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
