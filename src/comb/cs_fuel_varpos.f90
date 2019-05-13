!-------------------------------------------------------------------------------

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

subroutine cs_fuel_varpos
!========================

!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
!                COMBUSTION FUEL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!
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
use cs_fuel_incl
use field

!===============================================================================

implicit none

integer          icla, f_id
integer          keyccl, keydri, kscmin, kscmax
integer          iscdri
character(len=80) :: f_label, f_name

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key id of the fuel scalar class
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
call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
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

! Number of particles of the class icla per kg of air-fuel mixture (bulk)
do icla = 1, nclafu

  write(f_name,'(a,i2.2)') 'nd_fuel_', icla
  write(f_label,'(a,i2.2)') 'NG_FOL', icla
  call add_model_scalar_field(f_name, f_label, ing(icla))
  f_id = ivarfl(isca(ing(icla)))

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

! Mass fraction of ???? of the class icla per kg of bulk
do icla = 1, nclafu

  write(f_name,'(a,i2.2)') 'x_p_', icla
  write(f_label,'(a,i2.2)') 'YFOL_FOL', icla
  call add_model_scalar_field(f_name, f_label, iyfol(icla))
  f_id = ivarfl(isca(iyfol(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 4.d-1)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Enthalpy of the class icla per kg of bulk
! (Enthalpy of the class is the product of the mass fraction of the class
!  by massic enthalpy of the class).
do icla = 1, nclafu

  write(f_name,'(a,i2.2)') 'x_p_h_', icla
  write(f_label,'(a,i2.2)') 'Xp_Ent_', icla
  call add_model_scalar_field(f_name, f_label, ih2(icla))
  f_id = ivarfl(isca(ih2(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! min clipping is set later, as h02fol is known after reading fuel data
  ! call field_set_key_double(f_id, kscmin, h02fol)
  call field_set_key_double(f_id, kscmax, grand)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_comb_drift.ge.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

enddo

! Particles velocities (when they are transported)
if (i_comb_drift .eq. 1) then
  do icla = 1, nclafu

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

! Mass of vapor divided by the mass of bulk
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
call add_model_scalar_field('fr_vap', 'Fr_VAP', ifvap)
f_id = ivarfl(isca(ifvap))

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

! Combustion heterogene

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

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      Rq : pas de variance associee a un scalaire dans notre cas
!    - PROPRES A LA SUSPENSION : ICP
!===============================================================================

! Although we are in enthalpy formulation, we keep Cp constant

icp = -1

!----
! End
!----

return
end subroutine
