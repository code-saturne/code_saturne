!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine cs_user_metal_structures_source_terms                  &
 ( nvar   , nscal  ,                                              &
   ncmast ,                                                       &
   ltmast , itypst , izmast ,                                     &
   svcond , tmet  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppincl
use mesh
use field
use cs_c_bindings
use cs_f_interfaces

use cs_tagms
!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncmast
integer          izmet, met_znb

integer          ltmast(ncelet)
integer          izmast(ncelet)

double precision itypst(ncelet,nvar)
double precision svcond(ncelet,nvar)
double precision tmet

! Local variables

!< [loc_var_dec]
integer          icmst
integer          ifac, iel, iscal
integer          ivarh
integer          izone
integer          f_id

double precision hvap, tk

type(gas_mix_species_prop) s_h2o_g

double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: cvar_h
!< [loc_var_dec]

!===============================================================================

!< [init]
call field_get_id_try("y_h2o_g", f_id)
if (f_id.ne.-1) &
  call field_get_key_struct_gas_mix_species_prop(f_id, s_h2o_g)
!< [init]

!< [cells_selection]

!===============================================================================
! Select the cells which are associated to the metal structures volume
! with the subroutine getcel
!===============================================================================

izone = 0
met_znb = 1

do izmet = 1 , met_znb

  if (izmet.eq.1) then
    ! Cells associated to the geometric zone
    call getcel('z > -7.0d0 and z < 53.d0', ncmast, ltmast)

    izone = izone + 1

    do icmst = 1, ncmast
      iel = ltmast(icmst)
      izmast(iel) = izone
    enddo
  endif

enddo

!< [cells_selection]

!< [model_settings]

!===============================================================================
! Parameters padding of the wall thermal model and condensation model
! ------------------------------------------------------------------
! The condensation model can be used alone with a constant temperature
! specified by the user at the cold wall (tmet=tmet0 in this case)
! or together with a 0-D thermal model. In the latter case, the two models are
! coupled.
!===============================================================================

if (icondv.eq.0) then

  ! Choice the way to impose the wall temperature (tmet)
  ! at the solid/fluid interface:
  !
  ! with the parameter itagms defined below:
  ! ----------------------------------------
  !  0 : A constant wall temperature imposed is given by the user
  !     ( tmet = tmet0 used as the wall temperature of the metal structures
  !      by the condensation model )
  !  1 : A variable wall temperature is imposed with a 0-D thermal model
  !     ( tmet = tmet(iel,1) computed by cs_metal_structures_tag.f90
  !      and used as the wall temperature of the metal structures
  !      by the condensation model )

  itagms = 1

  ! Wall temperature computed by a 0-D thermal model
  ! with an explicit scheme and variable over time.
  ! ------------------------------------------------
  ! Remark : the wall temperature unit is [°C].
  if(itagms.eq.1) then
    !----------------------------------------
    ! (xem) thickness of the metal structures
    !----------------------------------------
    xem = 0.024d0
    !-------------------------------------------
    ! Initial condition of the 0-D thermal model
    !-------------------------------------------
    tmet0 = 92.d0
    !--------------------------------------------
    ! Physical properties of the metal structures
    !--------------------------------------------
    ! (xro) density (kg.m-3)
    xro_m   = 8000.d0
    ! (xcond) thermal conductivity (W.m-1.C-1)
    xcond_m = 12.8d0
    ! (xcp)   Specific heat (J.kg-1.C-1)
    xcp_m   = 500.0d0
  else
    ! Wall temperature imposed as constant
    ! with a value specified by the user
    tmet = 92.d0
  endif

endif

!< [model_settings]

!< [source_terms_values]

!===============================================================================
! The user can specify here the values of the following arrays used by the
! modelling of the metal structures condensation:
!         - itypst(:,ivar) to specify the condensation source term type,
!         - svcond(:,ivar) the scalar value to multiply by the sink term array
!                          of the metal structures condensation model.
! These two arrays can be filled for each transported scalar.
!===============================================================================

!-- pointer to the specific heat
if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

!-- pointer to the enthalpy value
ivarh = isca(iscalt)
call field_get_val_s(ivarfl(ivarh), cvar_h)

! loop over the cells associated to the metal structure
! source terms zone
do icmst = 1, ncmast
  iel = ltmast(icmst)

  ! Compute the enthalpy value of vapor gas
  if (ntcabs.le.1) then
    tk = t0
  else
    tk = cvar_h(iel)/cpro_cp(iel)
  endif
  hvap = s_h2o_g%cp*tk

  ! Condensation source terms associated
  ! to the metal structures imposed
  ! for each scalar.
  !----------------------------------------
  if (nscal.gt.0) then
    do iscal = 1, nscal
      if (iscal.eq.iscalt) then

        ! enthalpy value used for
        ! the explicit condensation term
        itypst(iel,isca(iscalt)) = 1
        svcond(iel,isca(iscalt)) = hvap
      else

        ! scalar values used for
        ! the explicit condensation term
        itypst(iel,isca(iscal)) = 1
        svcond(iel,isca(iscal)) = 0.d0
      endif
    enddo
  endif

enddo

!< [source_terms_values]

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_metal_structures_source_terms
