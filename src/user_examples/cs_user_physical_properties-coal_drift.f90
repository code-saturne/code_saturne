!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_user_physical_properties-coal_drift.f90
!>
!> \brief Definition of physical variable laws for scalars with a drift.
!>
!> See \ref physical_properties for examples.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
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
use field
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          ivart, iel, ifac
integer          ivar
integer          f_id
integer          icla
integer          iscdri, keydri, iflid, nfld, keyccl

double precision xvart
double precision aa1, bb1, cc1, dd1
double precision aa2, bb2, cc2, dd2
double precision aa3, bb3, cc3, dd3
double precision aa4, bb4, cc4, dd4
double precision aa5, bb5, cc5, dd5
double precision aa6, bb6, cc6, dd6
double precision aa7, bb7, cc7, dd7
double precision visco_O2, visco_CO, visco_H2, visco_N2
double precision visco_SO2, visco_NH3, visco_CO2

character*80     fname

double precision, allocatable, dimension(:) :: visco

double precision, dimension(:), pointer :: cpro_rom1, cpro_rom2, cpro_diam2
double precision, dimension(:), pointer :: cpro_temp, cpro_x2, cpro_x1
double precision, dimension(:), pointer :: cpro_ym1_3, cpro_ym1_5, cpro_ym1_7
double precision, dimension(:), pointer :: cpro_ym1_8
double precision, dimension(:), pointer :: cpro_ym1_9, cpro_ym1_11, cpro_ym1_12
double precision, dimension(:), pointer :: cpro_taup
double precision, dimension(:), pointer :: cpro_taupg
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! 0. Initializations to keep
!===============================================================================

!< [init]
allocate(visco(ncelet))

call field_get_val_s(iym1(3), cpro_ym1_3)
call field_get_val_s(iym1(5), cpro_ym1_5)
call field_get_val_s(iym1(7), cpro_ym1_7)
call field_get_val_s(iym1(8), cpro_ym1_8)
call field_get_val_s(iym1(9), cpro_ym1_9)
call field_get_val_s(iym1(11), cpro_ym1_11)
call field_get_val_s(iym1(12), cpro_ym1_12)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Number of fields
call field_get_n_fields(nfld)
!< [init]

!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!===============================================================================
!  Example:
!  =======
!===============================================================================

!  ===================================================================

!< [example_1]

! Temperature
call field_get_val_s(itemp, cpro_temp)

! Gas density
call field_get_val_s(irom1, cpro_rom1)

! First initialization
if (ntcabs.le.1) then
  do iel = 1, ncel
    visco(iel) = viscl0
    cpro_rom1(iel) = ro0
  enddo
  do icla = 1, nclacp
    call field_get_val_s(irom2(icla), cpro_rom2)
    call field_get_val_s(idiam2(icla), cpro_diam2)
    do iel = 1, ncel
      cpro_rom2(iel)  = rho20(icla)
      cpro_diam2(iel) = diam20(icla)
    enddo
  enddo
endif

!----------------Gas viscosity function of temperature--------------------------
!
!--------------------------1-O2 2-CO 3-H2 4-N2 5-SO2 6-NH3 7-CO2----------------

aa1 = 4.0495d-6
bb1 = 6.22d-8
cc1 = -2.3032d-11
dd1 = 4.4077d-15

aa2 = 9.9987d-6
bb2 = 5.1578d-8
cc2 = -1.8383d-11
dd2 = 3.33307d-15

aa3 = 2.894d-6
bb3 = 2.22508d-8
cc3 = -8.041d-12
dd3 = 1.4619d-15

aa4 = 4.3093d-6
bb4 = 5.0516d-8
cc4 = -1.7869d-11
dd4 = 3.2136d-15

aa5 = -1.9889d-6
bb5 = 5.365d-8
cc5 = -1.4286d-11
dd5 = 2.1639d-15

aa6 = -1.293d-6
bb6 = 4.1194d-8
cc6 = -1.772d-11
dd6 = 1.8699d-15

aa7 = 4.4822d-7
bb7 = 5.4327d-8
cc7 = -1.7581d-11
dd7 = 2.9979d-15

!-------------------------------------------------------------------------------
!      law                    mu   = a + b T + c T**2 + d T**3
!      so      cpro_viscl(iel) = a +b*xvart+c*xvart**2 + d*xvart**3
!-------------------------------------------------------------------------------

if (ntcabs.gt.1) then
  do iel = 1, ncel

    xvart = cpro_temp(iel)
    visco_O2  = aa1 + xvart*bb1 + cc1*xvart**2 + dd1*xvart**3
    visco_CO  = aa2 + xvart*bb2 + cc2*xvart**2 + dd2*xvart**3
    visco_H2  = aa3 + xvart*bb3 + cc3*xvart**2 + dd3*xvart**3
    visco_N2  = aa4 + xvart*bb4 + cc4*xvart**2 + dd4*xvart**3
    visco_SO2 = aa5 + xvart*bb5 + cc5*xvart**2 + dd5*xvart**3
    visco_NH3 = aa6 + xvart*bb6 + cc6*xvart**2 + dd6*xvart**3
    visco_CO2 = aa7 + xvart*bb7 + cc7*xvart**2 + dd7*xvart**3

    ! Viscosity of the mixing
    visco(iel) = ( cpro_ym1_8(iel) * visco_O2                     &
                 + cpro_ym1_3(iel) * visco_CO                     &
                 + cpro_ym1_5(iel) * visco_H2                     &
                 + cpro_ym1_12(iel)* visco_N2                     &
                 + cpro_ym1_11(iel)* visco_SO2                    &
                 + cpro_ym1_7(iel) * visco_NH3                    &
                 + cpro_ym1_9(iel) * visco_CO2 )/                 &
                 ( cpro_ym1_8(iel) + cpro_ym1_3(iel)              &
                 + cpro_ym1_5(iel) + cpro_ym1_12(iel)             &
                 + cpro_ym1_11(iel)+ cpro_ym1_7(iel)              &
                 + cpro_ym1_9(iel))

  enddo
endif

! get x1 = 1 - sum cpro_x2
call field_get_val_s_by_name("x_c", cpro_x1)

! All gas scalars have the same drift as if1m(1)
!-----------------------------------------------

do iflid = 0, nfld-1

  ! Index of the scalar class (<0 if the scalar belongs to the gas phase)
  call field_get_key_int(iflid, keyccl, icla)

  call field_get_key_int(iflid, keydri, iscdri)

  ! We only handle here one scalar with a drift per gas class
  if (icla.le.-1.and.btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    ! Position of variables, coefficients
    ! -----------------------------------

    ! Name of the drift scalar
    call field_get_name(iflid, fname)

    ! Index of the corresponding "relaxation time" (cpro_taupg) for the gas
    ! WARNING: for the gas, this tau might be negative
    call field_get_id('drift_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, cpro_taupg)

    ! Initialize to 0
    do iel = 1, ncel
      cpro_taupg(iel) = 0.d0
    enddo

  endif
enddo

! Loop over coal particle classes
! We only handle here coal class with a drift
!--------------------------------------------

do iflid = 0, nfld-1

  ! Index of the scalar class (<0 if the scalar belongs to the gas phase)
  call field_get_key_int(iflid, keyccl, icla)

  call field_get_key_int(iflid, keydri, iscdri)

  ! We only handle here one scalar with a drift per particle class
  if (icla.ge.1.and.btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    call field_get_val_s(irom2(icla), cpro_rom2)
    call field_get_val_s(idiam2(icla), cpro_diam2)
    call field_get_val_s(ix2(icla), cpro_x2)

    ! Position of variables, coefficients
    ! -----------------------------------

    ! Name of the drift scalar
    call field_get_name(iflid, fname)

    ! Index of the corresponding relaxation time (cpro_taup)
    call field_get_id('drift_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, cpro_taup)

    ! Computation of the relaxation time of the particles
    ! the drift is therefore v_g = tau_p * g
    !----------------------------------------------------

    do iel = 1, ncel

      ! Simple model for Low Reynolds Numbers
      cpro_taup(iel) = cpro_x1(iel) * cpro_rom2(iel)                &
                     * cpro_diam2(iel)**2                           &
                     / (18.d0*visco(iel))

    enddo

    ! Drift for the gas:
    ! tau_pg = - Sum_i X2_i v_gi
    do iel = 1, ncel

      cpro_taupg(iel) = cpro_taupg(iel)                                  &
                      - ( cpro_taup(iel) * cpro_x2(iel) )

    enddo

  endif ! test icla

enddo ! loop on iflid

!< [example_1]

!Free memory
deallocate(visco)

!===============================================================================

!----
! End
!----

return
end subroutine usphyv
