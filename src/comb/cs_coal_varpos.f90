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

subroutine cs_coal_varpos
!========================

!===============================================================================
!  FONCTION  :
!  ---------
!       INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
!                COMBUSTION CHARBON PULVERISE
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use ihmpre
use cs_coal_incl
use field

!===============================================================================

implicit none

integer          icla,  icha, f_id
integer          keyccl, keydri, kscmin, kscmax
integer          iscdri
integer(c_int) :: n_coals, n_classes
character(len=80) :: f_label, f_name

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_coal_combustion(n_coals, n_classes)  &
    bind(C, name='cs_field_pointer_map_coal_combustion')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value        :: n_coals, n_classes
  end subroutine cs_field_pointer_map_coal_combustion

end interface

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. Definition of fields
!===============================================================================

! Thermal model

itherm = 2
call add_model_scalar_field('enthalpy', 'Enthalpy', ihm) !FIXME enth_bulk?
iscalt = ihm

! Set min and max clipping
f_id = ivarfl(isca(iscalt))
call field_set_key_double(f_id, kscmin, -grand)
call field_set_key_double(f_id, kscmax, grand)

! Activate the drift: 0 (no activation),
!                     1 (transported particle velocity)
!                     2 (limit drop particle velocity)
if (i_comb_drift.ge.1) then
  iscdri = 1
endif

! Dispersed phase variables
!--------------------------

! NB: 'c' stands for continuous <> 'p' stands for particles

! "ibset" function set in the variable "iscdri" the fact that
! that the first element of the class (here Np) creates the convective
! flux of the class. The convective flux of the calls is common to all the
! elements of the class!
! note that the function "ibclr" means that the other element of the class
! DO NOT create any additional convective flux (but use the one of the class)

! Number of particles of the class icla per kg of air-coal mixture (bulk)
do icla = 1, nclacp

  write(f_name,'(a,i2.2)') 'n_p_', icla
  write(f_label,'(a,i2.2)') 'Np_', icla
  call add_model_scalar_field(f_name, f_label, inp(icla))
  f_id = ivarfl(isca(inp(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, rinfin)

  ! Scalar with drift: DO create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Mass fraction of reactive coal of the class icla per kg of bulk
do icla = 1, nclacp

  write(f_name,'(a,i2.2)') 'x_p_coal_', icla
  write(f_label,'(a,i2.2)') 'Xp_Ch_', icla
  call add_model_scalar_field(f_name, f_label, ixch(icla))
  f_id = ivarfl(isca(ixch(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Mass fraction of char (coke in French) of the class icla per kg of bulk
do icla = 1, nclacp

  write(f_name,'(a,i2.2)') 'x_p_char_', icla
  write(f_label,'(a,i2.2)') 'Xp_Ck_', icla
  call add_model_scalar_field(f_name, f_label, ixck(icla))
  f_id = ivarfl(isca(ixck(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Mass fraction of water (within the particle) of the class icla per kg of bulk
if (ippmod(iccoal).eq.1) then

  do icla = 1, nclacp

    write(f_name,'(a,i2.2)') 'x_p_wt_', icla
    write(f_label,'(a,i2.2)') 'Xp_wt_', icla
    call add_model_scalar_field(f_name, f_label, ixwt(icla))
    f_id = ivarfl(isca(ixwt(icla)))

    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! Set min and max clipping
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    if (i_comb_drift.ge.1) then
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)
    endif

  enddo

endif

! Enthalpy of the class icla per kg of bulk
! (Enthalpy of the class is the product of the mass fraction of the class
!  by massic enthalpy of the class).
do icla = 1, nclacp

  write(f_name,'(a,i2.2)') 'x_p_h_', icla
  write(f_label,'(a,i2.2)') 'Xp_Ent_', icla
  call add_model_scalar_field(f_name, f_label, ih2(icla))
  f_id = ivarfl(isca(ih2(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, -grand)
  call field_set_key_double(f_id, kscmax, grand)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Age of the class icla time the Np (number of particle per kg of bulk)
if (i_comb_drift.ge.1) then

  do icla = 1, nclacp

    write(f_name,'(a,i2.2)') 'n_p_age_', icla
    write(f_label,'(a,i2.2)') 'Np_Age_', icla
    call add_model_scalar_field(f_name, f_label, inagecp(icla))
    f_id = ivarfl(isca(inagecp(icla)))

    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! Set min and max clipping
    ! TODO: test on ippmod(icoal) used to reproduce previous
    !       behavior; check if it was actually desired.
    if (ippmod(iccoal).eq.1) then
      call field_set_key_double(f_id, kscmin, 0.d0 )
      call field_set_key_double(f_id, kscmax, grand)
    endif

    ! Scalar with drift: BUT Do NOT create additional mass flux
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

    call field_set_key_int(f_id, keydri, iscdri)

  enddo

endif

! Particles velocities (when they are transported)
if (i_comb_drift .eq. 1) then
  do icla = 1, nclacp

    write(f_name,'(a,i2.2)') 'v_x_p_', icla
    write(f_label,'(a,i2.2)') 'Vp_X_', icla
    call add_model_scalar_field(f_name, f_label, iv_p_x(icla))
    f_id = ivarfl(isca(iv_p_x(icla)))

    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)

    write(f_name,'(a,i2.2)') 'v_y_p_', icla
    write(f_label,'(a,i2.2)') 'Vp_Y_', icla
    call add_model_scalar_field(f_name, f_label, iv_p_y(icla))
    f_id = ivarfl(isca(iv_p_y(icla)))

    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)

    write(f_name,'(a,i2.2)') 'v_z_p_', icla
    write(f_label,'(a,i2.2)') 'Vp_Z_', icla
    call add_model_scalar_field(f_name, f_label, iv_p_z(icla))
    f_id = ivarfl(isca(iv_p_z(icla)))

    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)

  enddo
endif

icla = -1

! Continuous phase variables
!---------------------------

! NB: 'c' stands for continuous <> 'p' stands for particles

! Field gas enthalpy

! Enthalpy of the gas phase per kg of bulk
! (The gas phase is a class with a negative icla, Enthalpy of the class
!  is the product of the mass fraction of the class
!  by massic enthalpy of the class).
f_name = 'x_c_h'
f_label = 'Xc_Ent'
call add_model_scalar_field(f_name, f_label, ihgas)
f_id = ivarfl(isca(ihgas))
! Set the index of the scalar class in the field structure
call field_set_key_int(f_id, keyccl, icla)

! The first gas salar contains the drift flux, the others
! Scalar with drift: DO create additional mass flux
if (i_comb_drift.ge.1) then
  iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)
endif

! Light (F8) and heavy (F9) volatile matter

! Mass of the mean value of the tracer 1 (representing the light
! volatiles released by the coal icha) divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
do icha = 1, ncharb

  write(f_name,'(a,i2.2)') 'fr_mv1_', icha
  write(f_label,'(a,i2.2)') 'Fr_mv1_', icha
  call add_model_scalar_field(f_name, f_label, if1m(icha))
  f_id = ivarfl(isca(if1m(icha)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! The first gas salar (x_c_h) creates the drift flux
  if (i_comb_drift.ge.1) then
    ! Scalar with drift: BUT Do NOT create additional mass flux
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Mass of the mean value of the tracer 2 (representing the heavy
! volatiles released by the coal icha) divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
do icha = 1, ncharb

  write(f_name,'(a,i2.2)') 'fr_mv2_', icha
  write(f_label,'(a,i2.2)') 'Fr_mv2_', icha
  call add_model_scalar_field(f_name, f_label, if2m(icha))
  f_id = ivarfl(isca(if2m(icha)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Mass of the Oxydant 2 divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (noxyd .ge. 2) then

  f_name = 'fr_oxyd2'
  f_label = 'FR_OXYD2'
  call add_model_scalar_field(f_name, f_label, if4m)
  f_id = ivarfl(isca(if4m))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Mass of the Oxydant 3 divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (noxyd .ge. 3) then

  f_name = 'fr_oxyd3'
  f_label = 'FR_OXYD3'
  call add_model_scalar_field(f_name, f_label, if5m)
  f_id = ivarfl(isca(if5m))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Mass of the water from coal drying divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (ippmod(iccoal).eq.1) then

  f_name = 'fr_h2o'
  f_label = 'FR_H2O'
  call add_model_scalar_field(f_name, f_label, if6m)
  f_id = ivarfl(isca(if6m))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Mass of the Carbon from coal oxydized by O2 divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
f_name = 'fr_het_o2'
f_label  = 'FR_HET_O2'
call add_model_scalar_field(f_name, f_label, if7m)
f_id = ivarfl(isca(if7m))

! Set the index of the scalar class in the field structure
call field_set_key_int(f_id, keyccl, icla)

! Set min and max clipping
call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

! Scalar with drift: BUT Do NOT create additional mass flux
if (i_comb_drift.ge.1) then
  iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)
endif

! Mass of the Carbon from coal gasified by CO2 divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (ihtco2.eq.1) then

  f_name = 'fr_het_co2'
  f_label  = 'FR_HET_CO2'
  call add_model_scalar_field(f_name, f_label, if8m)
  f_id = ivarfl(isca(if8m))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Mass of the Carbon from coal gasified by H2O divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (ihth2o.eq.1) then

  f_name = 'fr_het_h2o'
  f_label  = 'FR_HET_H2O'
  call add_model_scalar_field(f_name, f_label, if9m)
  f_id = ivarfl(isca(if9m))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Variance
f_name = 'f1f2_variance'
f_label  = 'Var_F1F2'
call add_model_scalar_field(f_name, f_label, ifvp2m)
f_id = ivarfl(isca(ifvp2m))

! Set the index of the scalar class in the field structure
call field_set_key_int(f_id, keyccl, icla)

! Set min and max clipping
call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 0.25d0)

! Scalar with drift: BUT Do NOT create additional mass flux
if (i_comb_drift.ge.1) then
  iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)
endif

! Mass of the Carbon dioxyde (CO or CO2) divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
!FIXME check for the oxycombustion, it would be more relevant to track CO
if (ieqco2.ge.1) then

  f_name = 'x_c_co2'
  f_label  = 'Xc_CO2'
  call add_model_scalar_field(f_name, f_label, iyco2)
  f_id = ivarfl(isca(iyco2))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Mass of the HCN divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
if (ieqnox.eq.1) then

  f_name = 'x_c_hcn'
  f_label = 'Xc_HCN'
  call add_model_scalar_field(f_name, f_label, iyhcn)
  f_id = ivarfl(isca(iyhcn))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! Mass of the NH3 divided by the mass of bulk
  ! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  f_name =  'x_c_nh3'
  f_label =  'Xc_NH3'
  call add_model_scalar_field(f_name, f_label, iynh3)
  f_id = ivarfl(isca(iynh3))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! Mass of the NO divided by the mass of bulk
  ! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  f_name =  'x_c_no'
  f_label =  'Xc_NO'
  call add_model_scalar_field(f_name, f_label, iyno)
  f_id = ivarfl(isca(iyno))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! Enthalpy of the oxydizer times the fraction of gas divided by the mass of bulk
  ! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
  f_name =  'x_c_h_ox'
  f_label =  'Xc_Ent_Ox'
  call add_model_scalar_field(f_name, f_label, ihox)
  f_id = ivarfl(isca(ihox))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, -grand)
  call field_set_key_double(f_id, kscmax,  grand)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

endif

! Age of bulk
!FIXME give the possibility of having the age seperately
if (i_comb_drift.ge.1) then

  f_name = 'age'
  f_label = 'Age'
  call add_model_scalar_field(f_name, f_label, iage)
  f_id = ivarfl(isca(iage))

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0 )
  call field_set_key_double(f_id, kscmax, grand)
endif

! Map to field pointers

n_coals = ncharb
n_classes = nclacp

call cs_field_pointer_map_coal_combustion(n_coals, n_classes)

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!    - PROPRES A LA SUSPENSION : ICP
!===============================================================================

! Although we are in enthalpy formulation, we keep Cp constant

icp = -1

return
end subroutine
