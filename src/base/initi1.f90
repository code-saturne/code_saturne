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

!> \file initi1.f90
!> \brief Commons initialization.
!>
!------------------------------------------------------------------------------

subroutine initi1

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor
use ppincl, only: ippmod, nmodmx, iatmos
use post
use cs_c_bindings
use field
use lagran
use cpincl
use dimens
use numvar, only : nscaus, nscapp
use radiat
use cdomod

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iok, ipp, nmodpp, have_thermal_model

!===============================================================================

procedure() :: iniini, iniusi, ppini1, addfld, modini, fldini
procedure() ::  usipes, verini, impini

interface

  subroutine cs_gui_boundary_conditions_define(bdy)  &
    bind(C, name='cs_gui_boundary_conditions_define')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: bdy
  end subroutine cs_gui_boundary_conditions_define

  subroutine cs_gui_output()  &
      bind(C, name='cs_gui_output')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_output

  subroutine cs_lagr_add_fields()  &
      bind(C, name='cs_lagr_add_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_lagr_add_fields

  subroutine cs_parameters_global_complete()  &
      bind(C, name='cs_parameters_global_complete')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_global_complete

  subroutine cs_parameters_eqp_complete()  &
      bind(C, name='cs_parameters_eqp_complete')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_eqp_complete

  subroutine cs_parameters_output_complete()  &
      bind(C, name='cs_parameters_output_complete')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_output_complete

  subroutine user_boundary_conditions_setup_wrapper()  &
      bind(C, name='cs_user_boundary_conditions_setup_wrapper')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_boundary_conditions_setup_wrapper

  subroutine user_finalize_setup_wrapper()  &
      bind(C, name='cs_user_finalize_setup_wrapper')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_finalize_setup_wrapper

  subroutine user_syrthes_coupling()  &
      bind(C, name='cs_user_syrthes_coupling')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_syrthes_coupling

end interface

!===============================================================================
! Initialize modules before user has access
!===============================================================================

call iniini

nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
  endif
enddo

!===============================================================================
! User input, variable definitions
!===============================================================================

call iniusi

call ppini1

!===============================================================================
! Map Fortran pointers to C global data
!===============================================================================

call elec_option_init

if (ippmod(iatmos).ge.0) call cs_at_data_assim_initialize

! Lagrangian model options

call lagran_init_map

have_thermal_model = 0
if (iscalt.ge.1) then
  have_thermal_model = 1
endif
call cs_lagr_options_definition(isuite, have_thermal_model, dtref, iccvfg)
call cs_lagr_add_fields

! Additional fields if not in CDO mode only

if (icdo.lt.2) then
  call addfld
endif

!===============================================================================
! Changes after user initialization and additional fields dependent on
! main fields and options.
!===============================================================================

! Do not call this routine if CDO mode only (default variables and properties
! are not defined anymore)
if (icdo.lt.2) then
  call cs_parameters_global_complete
  call modini
  call fldini
endif

call cs_parameters_eqp_complete

nscal = nscaus + nscapp

!===============================================================================
! Time moments called after additionnal creation
!===============================================================================

call cs_gui_time_moments
call cs_user_time_moments

!===============================================================================
! GUI-based boundary condition definitions
!===============================================================================

call cs_gui_boundary_conditions_define(c_null_ptr)

!===============================================================================
! Some final settings
!===============================================================================

! Postprocessing and logging

call cs_gui_output

! Do not call this routine if CDO mode only (default variables and properties
! are not defined anymore)
if (icdo.lt.2) then
   call usipes(nmodpp)
   ! Avoid a second spurious call to this function
   ! Called in the C part if CDO is activated, i.e. when
   ! additional geometric quantities and connectivities are built
   if (icdo.lt.0) then
     call user_boundary_conditions_setup_wrapper
     call user_finalize_setup_wrapper
   endif
endif

call cs_parameters_output_complete

!===============================================================================
! Coherency checks
!===============================================================================

iok = 0

! No verification in CDO mode only. This is done elsewhere
if (icdo.lt.2) then
   call verini (iok)
   call parameters_check
endif

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
else
  write(nfecra,9998)
endif

 9998 format(                                                   /,&
' No error detected during the data verification'              ,/,&
'                          cs_user_parameters.f90 and others).',/)
 9999 format(                                                     &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,&
'@    ========'                                                ,/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run (',i10,' errors).'         ,/,&
'@'                                                            ,/,&
'@  See previous impressions for more informations.'           ,/,&
'@  Verify the provided data in the interface'                 ,/,&
'@    and in user subroutines.'                                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

!===============================================================================
! Output
!===============================================================================

call impini

return
end subroutine
