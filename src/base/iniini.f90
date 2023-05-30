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

!> \file iniini.f90
!> \brief Commons default initialization before handing over the user.
!>
!------------------------------------------------------------------------------

subroutine iniini

!===============================================================================
! Module files
!===============================================================================

use atincl
use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use rotation
use entsor
use pointe
use albase
use alaste
use parall
use period
use cplsat
use ppincl
use ppcpfu
use mesh
use field
use vof
use cavitation
use radiat
use turbomachinery
use cs_nz_condensation, only: init_sizes_pcond
use ctincl
use cfpoin
use vof
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ii, iscal, iest

!===============================================================================

!===============================================================================
! 0. STOCKAGE DES ARGUMENTS ET IMPRESSIONS INITIALES
!===============================================================================

write(nfecra, 900)

 900  format(/,                                                   &
'===============================================================',&
/,/,                                                        &
'                   CALCULATION PREPARATION'                   ,/,&
'                   ======================='                   ,/,&
                                                                /,&
                                                                /,&
' ===========================================================' ,/,&
                                                                /,&
                                                                /)

!===============================================================================
! 0. Global field keys
!===============================================================================

call field_get_key_id("label", keylbl)
call field_get_key_id('log', keylog)
call field_get_key_id('post_vis', keyvis)

call field_get_key_id("inner_mass_flux_id", kimasf)
call field_get_key_id("boundary_mass_flux_id", kbmasf)

call field_get_key_id("diffusivity_id", kivisl)
call field_get_key_id("diffusivity_ref", kvisl0)

call field_get_key_id("is_temperature", kscacp)

call field_get_key_id("density_id", kromsl)

call field_get_key_id("gradient_weighting_id", kwgrec)

call field_get_key_id("source_term_prev_id", kstprv)
call field_get_key_id("source_term_id", kst)

call field_get_key_id("turbulent_schmidt", ksigmas)

call field_get_key_id("turbulent_flux_ctheta", kctheta)

icrom = -1
ibrom = -1

ipori = -1
iporf = -1

!===============================================================================
! Map Fortran pointers to C global data
!===============================================================================

call atmo_init
call time_step_init
call time_step_options_init
call thermal_model_init
call turb_model_init
call turb_rans_model_init
call turb_les_model_init
call turb_hybrid_model_init
call turb_model_constants_init
call wall_functions_init
call physical_constants_init
call porosity_ibm_init
call porosity_from_scan_init
call fluid_properties_init
call space_disc_options_init
call time_scheme_options_init
call velocity_pressure_options_init
call restart_auxiliary_options_init
call turb_reference_values_init
call listing_writing_period_init
call radiat_init
call gas_mix_options_init
call ctwr_properties_init
call map_ale
call cf_model_init
call vof_model_init
call cavitation_model_init

call map_turbomachinery_model(iturbo, ityint)
call init_sizes_pcond()

!===============================================================================
! I/O: entsor.f90
!===============================================================================

! ---> NFECRA vaut 6 par defaut ou 9 en parallele (CSINIT)

! ---> Fichier thermochinie
!        FPP : utilisateur
!        JNF : Janaf
!        Les deux fichiers peuvent partager la meme unite
!          puisqu'ils sont lus l'un a pres l'autre.
!      En prime, INDJON (janaf=1 ou non=0)

impfpp = 25
ficfpp = 'define_ficfpp_in_usppmo'

! ---> Fichiers module atmospherique
call atmo_set_meteo_file_name('meteo')

! ---> Fichiers utilisateurs

do ii = 1, nusrmx
  impusr(ii) = 69+ii
enddo

! Here entsor.f90 is completely initialized

!===============================================================================
! Get mesh metadata.
!===============================================================================

call ledevi(iperio, iperot)

call tstjpe(iperio, iperot)

!===============================================================================
! Position of variables in numvar.f90
!===============================================================================

! Initialize mappings of field ids

do ii = 1, nvarmx
  ivarfl(ii) = -1
enddo

! Scalar to variable mappings

do iscal = 1, nscamx
  isca  (iscal) = 0
  iscapp(iscal) = 0
  iscasp(iscal) = 0
enddo

! Default initialization for specific physical models

call ppinii

!===============================================================================
! Arrays of optcal.f90
!===============================================================================

! Ordering of BC's will be computed after cs_user_boundary_conditions.

do ii = 1, ntypmx
  idebty(ii) = 0
  ifinty(ii) = 0
enddo

! Here, all of optcal.f90 is initialized

!===============================================================================
! Arrays of cstphy.f90
!===============================================================================

! Reset pther to p0
pther = p0

!===============================================================================
! Lagrangian arrays
!===============================================================================

tslagr => rvoid2

return
end subroutine
