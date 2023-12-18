!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file fldvar.f90
!> \brief Variables definition initialization, according to calculation type
!> selected by the user.
!
!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[out]    nmodpp        number of activated particle physic models
!______________________________________________________________________________

subroutine fldvar &
( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use albase
use lagran
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings
use cfpoin, only:ieos

!===============================================================================

implicit none

! Arguments

integer       nmodpp

! Local variables

integer       ipp
integer       iok, keycpl, nmodpp_compatibility, vof_mask
integer       key_lim_id, kscmin, kscmax

type(var_cal_opt) :: vcopt

!===============================================================================
! Interfaces
!===============================================================================

procedure() :: add_variable_field, add_cdo_variable_field
procedure() :: pplecd, ppvarp, add_model_scalar_field
procedure() :: add_user_scalar_fields, fldvar_check_nvar
procedure() :: init_var_cal_opt

interface

  ! Interface to C function returning number of user-defined variables

  function cs_parameters_n_added_variables() result(n) &
    bind(C, name='cs_parameters_n_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int)                                           :: n
  end function cs_parameters_n_added_variables

end interface

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

call field_get_key_id('coupled', keycpl)

!===============================================================================
! Calcul de nscapp
! Verification du nombre de scalaires
! Construction de iscapp
! Calcul de nscal

!  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
!  On en profite aussi pour remplir ITYTUR et ITYTURT puisque ITURB et ITURT
!    viennent d'etre definis.
!===============================================================================

! ---> Remplissage de ITYTUR
itytur = iturb/10

! ---> Coherence modele
!     Rq : ATTENTION il faudrait renforcer le blindage

iok   = 0
nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
    if (ippmod(ipp).lt.0 .or. ippmod(ipp).gt.5) then
      write(nfecra,6001)
      iok = iok + 1
    endif
  endif
enddo

nmodpp_compatibility = nmodpp

! Compressible module and gas mix are compatible
if (ippmod(igmix).ne.-1 .and. ippmod(icompf) .ne. -1) then
  nmodpp_compatibility = nmodpp_compatibility - 1
endif

! Atmo in humid atmosphere et Couling tower (iaeros) coupling
if (ippmod(iatmos).eq.2 .and. ippmod(iaeros) .ne. -1) then
  nmodpp_compatibility = nmodpp_compatibility - 1
endif

if (nmodpp_compatibility.gt.1) then
  write(nfecra,6000)
  iok = iok + 1
endif

! In case ideal gas mix specific physics was enabled by the user
! together with the compressible module, the equation of state
! indicator is reset to the approprate value automatically (ieos=3)
! and the user is warned.
if (ippmod(igmix).ge.0.and.ippmod(icompf).ge.0.and.ieos.ne.3) then
  ieos = 3
  write(nfecra,6002)
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

! Set global indicator: ippmod(iphpar)
!  0: no specific model
!  1: active specific physical model
ippmod(iphpar) = 0
if (nmodpp.gt.0) then
  ippmod(iphpar) = 1
endif

! Define main variables
!======================

nvar = 0

! Velocity

call add_variable_field('velocity', 'Velocity', 3, iu)
call field_set_key_int(ivarfl(iu), keycpl, 1)

! All components point to same field
iv = iu + 1
iw = iv + 1

! Pressure or hydraulic head for groundwater flow module

call add_variable_field('pressure', 'Pressure', 1, ipr)

! Enabled VoF model if free surface or mass transfer modeling enabled
vof_mask = ior(VOF_FREE_SURFACE, VOF_MERKLE_MASS_TRANSFER)
if (iand(ivofmt, vof_mask).ne.0) &
     ivofmt = ior(VOF_ENABLED, ivofmt)

! Mass balance equation options (pressure)

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! elliptic equation
vcopt%iconv = 0

! compressible algorithm
if ((ippmod(icompf).ge.0).or.(idilat.eq.2.and.ieos.gt.-1)) then
  vcopt%istat = 1
else
  vcopt%istat = 0
endif

! VoF algorithm: activate the weighting for the pressure
if (ivofmt.gt.0) vcopt%iwgrec = 1

call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! void fraction (VoF algorithm)

if (ivofmt.gt.0) then
  call add_variable_field('void_fraction', 'Void Fraction', 1, ivolf2)
  call field_get_key_struct_var_cal_opt(ivarfl(ivolf2), vcopt)
  vcopt%idiff = 0  ! pure convection equation

  ! NVD/TVD scheme
  vcopt%ischcv = 4
  call field_get_key_id("limiter_choice", key_lim_id)
  ! (CICSAM limiter)
  call field_set_key_int(ivarfl(ivolf2), key_lim_id, 11)
  ! Beta Limiter
  vcopt%isstpc = 2

  call field_set_key_struct_var_cal_opt(ivarfl(ivolf2), vcopt)
  ! Bounds for the beta limiter
  call field_get_key_id("min_scalar", kscmin)
  call field_get_key_id("max_scalar", kscmax)
  call field_set_key_double(ivarfl(ivolf2), kscmin, 0.d0)
  call field_set_key_double(ivarfl(ivolf2), kscmax, 1.d0)
endif

! Turbulence

if (itytur.eq.2) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
else if (itytur.eq.3) then
  call add_variable_field('rij', 'Rij', 6, irij)
  call field_set_key_int(ivarfl(irij), keycpl, 1)

  ! All rij components point to same field
  ir11 = irij
  ir22 = ir11 + 1
  ir33 = ir22 + 1
  ir12 = ir33 + 1
  ir23 = ir12 + 1
  ir13 = ir23 + 1

  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
  if (iturb.eq.32) then
    call add_variable_field('alpha', 'Alphap', 1, ial)
    ! Elliptic equation (no convection, no time term)
    call field_get_key_struct_var_cal_opt(ivarfl(ial), vcopt)
    vcopt%istat = 0
    vcopt%iconv = 0
    ! For alpha, we always have a diagonal term, so do not shift the diagonal
    vcopt%idircl = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ial), vcopt)
  endif
else if (itytur.eq.5) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
  call add_variable_field('phi', 'Phi', 1, iphi)
  if (iturb.eq.50) then
    call add_variable_field('f_bar', 'f_bar', 1, ifb)
    call field_get_key_struct_var_cal_opt(ivarfl(ifb), vcopt)
    vcopt%istat = 0
    vcopt%iconv = 0
    ! For fb, we always have a diagonal term, so do not shift the diagonal
    vcopt%idircl = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ifb), vcopt)
  else if (iturb.eq.51) then
    call add_variable_field('alpha', 'Alpha', 1, ial)
    call field_get_key_struct_var_cal_opt(ivarfl(ial), vcopt)
    vcopt%istat = 0
    vcopt%iconv = 0
    ! For alpha, we always have a diagonal term, so do not shift the diagonal
    vcopt%idircl = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ial), vcopt)
  endif
else if (iturb.eq.60) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('omega', 'Omega', 1, iomg)
else if (iturb.eq.70) then
  call add_variable_field('nu_tilda', 'NuTilda', 1, inusa)
endif

! Mesh velocity with ALE
if (iale.ge.1) then

  ! field defined on vertices if CDO-Vb scheme is used
  if (iale.eq.2) then
    call add_cdo_variable_field('mesh_velocity', 'Mesh Velocity',    &
                                3, MESH_LOCATION_VERTICES, 1, iuma)
  else
    call add_variable_field('mesh_velocity', 'Mesh Velocity', 3, iuma)
  endif

  call field_set_key_int(ivarfl(iuma), keycpl, 1)

  ivma = iuma + 1
  iwma = ivma + 1

  call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
  vcopt%istat = 0
  vcopt%iconv = 0
  call field_set_key_struct_var_cal_opt(ivarfl(iuma), vcopt)

endif

! Number of user variables

nscaus = cs_parameters_n_added_variables()

! Specific physics variables
call ppvarp

! Thermal model with no specific physics

if (nmodpp.eq.0) then

  if (itherm .eq. 1) then
    call add_model_scalar_field('temperature', 'Temperature', iscalt)
  else if (itherm .eq. 2) then
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
  else if (itherm .eq. 4) then
    call add_model_scalar_field('internal_energy', 'Eint', iscalt)
  endif

endif

call add_user_scalar_fields

! Map pointers

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

! ---> Verifications

iok = 0

if (nscaus.lt.0) then
  write(nfecra,6010) nscaus
  iok = iok + 1
endif

if (nscaus.gt.0 .or. nscapp.gt.0) then
  if ((nscaus+nscapp).gt.nscamx) then
    if (nscapp.le.0) then
      write(nfecra,6011) nscaus, nscamx, nscamx, nscaus
    else
      write(nfecra,6012) nscaus,nscapp,nscamx,nscamx-nscapp,nscaus+nscapp
    endif
    iok = iok + 1
  endif
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

return

!===============================================================================
! 5. FORMATS
!===============================================================================

 6000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     SEVERAL INCOMPATIBLE MODELS OF SPECIFIC PHYSICS ARE    ',/,&
'@     ARE ENABLED.                                           ',/,&
'@                                                            ',/,&
'@  Only the compressible and gas mix specific physics can be ',/,&
'@    enabled simultaneously.                                 ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     WRONG SELLECTION OF THE MODEL FOR SPECIFIC PHYSICS     ',/,&
'@                                                            ',/,&
'@  The values of the indices of the array IPPMOD are not     ',/,&
'@    admissible                                              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : AT THE INITIAL DATA VERIFICATION            ',/,&
'@    =========                                               ',/,&
'@     EQUATION OF STATE INCOMPATIBLE WITH SELECTED SPECIFIC  ',/,&
'@     PHYSICS                                                ',/,&
'@                                                            ',/,&
'@  The specific physics compressible and gas mix are         ',/,&
'@    simultaneously enabled but the selected equation of     ',/,&
'@    state is not ideal gas mix (ieos different from 3).     ',/,&
'@                                                            ',/,&
'@  The indicator ieos has been reset to 3 and the calculation',/,&
'@  will run.                                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     ERRONEOUS NUMBER OF SCALARS                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars must be an integer either     ',/,&
'@   positive or zero. Here is      NSCAUS  = ',I10            ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@  requested                          is   NSCAUS = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx           is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximmum value allowed of   NSCAUS                    ',/,&
'@                          is in   NSCAMX        = ',I10      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ',I10     ,/,&
'@  The number of scalars necessary for the specific physics'  ,/,&
'@    with the chosen model is              NSCAPP = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximum value allowed for  NSCAUS                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!===============================================================================

!> \brief add field defining a general solved variable, with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name          field name
!> \param[in]  label         field default label, or empty
!> \param[in]  dim           field dimension
!> \param[out] ivar          variable number for defined field
!_______________________________________________________________________________

subroutine add_variable_field &
 ( name, label, dim, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
integer, intent(out)         :: ivar

! Local variables

integer  id, ii

integer, save :: keyvar = -1

! Create field

call variable_field_create(name, label, MESH_LOCATION_CELLS, dim, id)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

ivarfl(ivar) = id

call field_set_key_int(id, keyvar, ivar)

call init_var_cal_opt(id)

if (dim .gt. 1) then
  do ii = 2, dim
    ivarfl(ivar + ii - 1) = id
  enddo
endif

return

end subroutine add_variable_field

!===============================================================================

!> \brief Add a field defining a general solved variable, with default options
!>        This variable is solved with a CDO scheme.
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name          field name
!> \param[in]  label         field default label, or empty
!> \param[in]  dim           field dimension
!> \param[in]  location_id   id of the mesh location where the field is defined
!> \param[in]  has_previous  if greater than 0 then stores previous state
!> \param[out] ivar          variable number for defined field
!_______________________________________________________________________________

subroutine add_cdo_variable_field &
 ( name, label, dim, location_id, has_previous, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim, location_id, has_previous
integer, intent(out)         :: ivar

! Local variables

integer  id, ii

integer, save :: keyvar = -1

! Create field

call variable_cdo_field_create(name, label, location_id, dim, has_previous, id)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

ivarfl(ivar) = id

call field_set_key_int(id, keyvar, ivar)

call init_var_cal_opt(id)

if (dim .gt. 1) then
  do ii = 2, dim
    ivarfl(ivar + ii - 1) = id
  enddo
endif

return

end subroutine add_cdo_variable_field

!===============================================================================

!> \brief add fields defining user solved scalar variables,
!>        with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine add_user_scalar_fields

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt

! Arguments

! Local variables

integer  iscal, nfld1, nfld2
integer  dim, id, ii, ivar, keycpl

integer :: keyvar, keysca

!===============================================================================
! Interfaces
!===============================================================================

interface

  ! Interface to C function building user-defined variables

  subroutine cs_parameters_create_added_variables() &
    bind(C, name='cs_parameters_create_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_create_added_variables

end interface

!===============================================================================

! Create fields

call field_get_n_fields(nfld1)

call cs_parameters_create_added_variables

call field_get_n_fields(nfld2)

! Now map those fields

iscal = 0

call field_get_key_id('coupled', keycpl)
call field_get_key_id("scalar_id", keysca)
call field_get_key_id("variable_id", keyvar)

do id = nfld1, nfld2 - 1

  call field_get_dim(id, dim)

  if (dim.eq.3) then
    call field_set_key_int(id, keycpl, 1)
  else if (dim.ne.1) then
    cycle
  endif

  iscal = iscal + 1

  ivar = nvar + 1
  nvar = nvar + dim
  nscal = nscal + 1

  ! Check we have enough slots
  call fldvar_check_nvar

  isca(iscal) = ivar
  ivarfl(ivar) = id

  call field_set_key_int(id, keyvar, ivar)
  call field_set_key_int(id, keysca, iscal)
  call init_var_cal_opt(id)

  if (dim .gt. 1) then
    do ii = 2, dim
      ivarfl(ivar + ii - 1) = id
    enddo
  endif

enddo

return

!---
! Formats
!---

end subroutine add_user_scalar_fields

!===============================================================================

!> \brief add field defining a non-user solved scalar variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[out] iscal          scalar number for defined field
!_______________________________________________________________________________

subroutine add_model_scalar_field &
 ( name, label, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: add_model_field

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(out)         :: iscal

! Local variables

integer  dim

dim = 1

call add_model_field(name, label, dim, iscal)

return

end subroutine add_model_scalar_field

!===============================================================================
!
!> \brief add field defining a non-user solved variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[in]  dim            field dimension
!> \param[out] iscal          variable number for defined field
!_______________________________________________________________________________

subroutine add_model_field &
 ( name, label, dim, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field
use cs_c_bindings

!===============================================================================

implicit none

procedure() :: add_model_field_indexes

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
integer, intent(out)         :: iscal

! Local variables

integer  id
integer  location_id

location_id = 1 ! variables defined on cells

! Create field

call variable_field_create(name, label, location_id, dim, id)

call add_model_field_indexes(id, iscal)

return

end subroutine add_model_field

!===============================================================================
!
!> \brief add field indexes associated with a new non-user solved
!>        scalar variable, with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  f_id           field id
!> \param[out] ivar           variable number for defined field
!_______________________________________________________________________________

subroutine add_variable_field_indexes &
 ( f_id, ivar )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, init_var_cal_opt
procedure() :: csexit

! Arguments

integer, intent(in)  :: f_id
integer, intent(out) :: ivar

! Local variables

integer  dim, ii

integer, save :: keyvar = -1

! Get field dimension

call field_get_dim(f_id, dim)

if (keyvar.lt.0) then
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_set_key_int(f_id, keyvar, ivar)
call init_var_cal_opt(f_id)

return

end subroutine add_variable_field_indexes

!===============================================================================
!
!> \brief add field indexes associated with a new non-user solved
!>        scalar variable, with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  f_id           field id
!> \param[out] iscal          scalar id for defined field
!_______________________________________________________________________________

subroutine add_model_field_indexes &
 ( f_id, iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

procedure() :: fldvar_check_nvar, fldvar_check_nscapp, init_var_cal_opt
procedure() :: csexit

! Arguments

integer, intent(in)  :: f_id
integer, intent(out) :: iscal

! Local variables

integer  dim, ivar, ii

integer, save :: keyvar = -1
integer, save :: keysca = -1

! Get field dimension

call field_get_dim(f_id, dim)

if (keysca.lt.0) then
  call field_get_key_id("scalar_id", keysca)
  call field_get_key_id("variable_id", keyvar)
endif

ivar = nvar + 1
nvar = nvar + dim
nscal = nscal + 1
iscal = nscaus + nscapp + 1
nscapp = nscapp + 1

! Check we have enough slots
call fldvar_check_nvar
call fldvar_check_nscapp

isca(iscal) = ivar
iscapp(nscapp) = iscal

do ii = 1, dim
  ivarfl(ivar + ii - 1) = f_id
enddo

call field_set_key_int(f_id, keyvar, ivar)
call field_set_key_int(f_id, keysca, iscal)
call init_var_cal_opt(f_id)

return

end subroutine add_model_field_indexes

!===============================================================================

!> \brief Check nvarmx is sufficient for the required number of variables.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nvar

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar

!===============================================================================

implicit none

procedure() :: csexit

! Arguments

! Local variables

if (nvar .gt. nvarmx) then
  write(nfecra,1000) nvar, nvarmx
  call csexit (1)
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined                           ',/,&
'@    corresponds to a number of variables NVAR  >= ', i10     ,/,&
'@  The maximum number of variables allowed                   ',/,&
'@                      in   paramx   is  NVARMX  = ', i10     ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@  If NVARMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine fldvar_check_nvar

!===============================================================================

!> \brief Initialize the given variable calculation option structure with
!>        legacy values (iniini) allowing to later test user modification.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine init_var_cal_opt &
 ( id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use cs_c_bindings
use field
use optcal

!===============================================================================

implicit none

! Arguments
integer id

! Local variables
type(var_cal_opt) :: vcopt

! Most values set by default at in _var_cal_opt default;
! see cs_parameters.c

call field_get_key_struct_var_cal_opt(id, vcopt)

! Undefined values, may be modified by modini
vcopt%isstpc = -999
vcopt%nswrsm = -1
vcopt%thetav = -1.d0
vcopt%blencv = -1.d0
vcopt%epsilo = -1.d0
vcopt%epsrsm = -1.d0
vcopt%relaxv = -1.d0

call field_set_key_struct_var_cal_opt(id, vcopt)

return

end subroutine init_var_cal_opt

!===============================================================================

!> \brief Check nscamx is sufficient for the required number of model scalars.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nscapp

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar

!===============================================================================

implicit none

procedure() :: csexit

! Arguments

! Local variables

if ((nscaus+nscapp).gt.nscamx) then
  write(nfecra,1000) nscaus,nscamx,nscamx-nscaus
  call csexit (1)
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ', i10    ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         est  NSCAMX = ', i10    ,/,&
'@                                                            ',/,&
'@  The maximum value possible for NSCAPP                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAUS = ', i10    ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  If NSCAMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine fldvar_check_nscapp

!===============================================================================
! C bindings (reverse)
!===============================================================================

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new non-user solved
!>        variable, with default options
!
!> \param[in]  f_id    field id

!> \result             variable number for defined field
!-------------------------------------------------------------------------------

function cs_add_variable_field_indexes(f_id) result(ivar) &
  bind(C, name='cs_add_variable_field_indexes')

  use, intrinsic :: iso_c_binding
  use cs_c_bindings

  implicit none

  procedure() :: add_variable_field_indexes

  ! Arguments

  integer(c_int), value :: f_id
  integer(c_int) :: ivar

  ! Local variables

  integer f_id0, ivar0

  f_id0 = f_id

  call add_variable_field_indexes(f_id0, ivar0)

  ivar = ivar0

end function cs_add_variable_field_indexes

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new non-user solved
!>        variable, with default options
!
!> \param[in]  f_id    field id

!> \result             scalar number for defined field
!-------------------------------------------------------------------------------

function cs_c_add_model_field_indexes(f_id) result(iscal) &
  bind(C, name='cs_add_model_field_indexes')

  use, intrinsic :: iso_c_binding
  use cs_c_bindings

  implicit none

  procedure() :: add_model_field_indexes

  ! Arguments

  integer(c_int), value :: f_id
  integer(c_int) :: iscal

  ! Local variables

  integer f_id0, iscal0

  f_id0 = f_id

  call add_model_field_indexes(f_id0, iscal0)

  iscal = iscal0

end function cs_c_add_model_field_indexes

!-------------------------------------------------------------------------------
!> \brief add field indexes associated with a new solved thermal variable,
!>        with default options
!
!> \param[in]  f_id    field id
!-------------------------------------------------------------------------------

subroutine cs_c_add_model_thermal_field_indexes(f_id) &
  bind(C, name='cs_add_model_thermal_field_indexes')

  use, intrinsic :: iso_c_binding
  use optcal
  use cs_c_bindings

  implicit none

  procedure() :: add_model_field_indexes

  ! Arguments

  integer(c_int), value :: f_id

  ! Local variables

  integer f_id0, iscal0

  f_id0 = f_id

  call add_model_field_indexes(f_id0, iscal0)

  iscalt = iscal0

end subroutine cs_c_add_model_thermal_field_indexes

!---------------------------------------------------------------------------
