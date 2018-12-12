
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

!> \file cs_user_parameters-richards.f90
!>
!> \brief Darcy module parameters example.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

!> \brief User subroutine for selection of specific physics module

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if the XML file from the GUI is used
!>                              (1 : yes, 0 : no)
!______________________________________________________________________________!

subroutine usppmo &
 ( ixmlpu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

!< [richards_activ]
! Set groundwater flow module active (1: active, 0: inactive)
ippmod(idarcy) = 1
!< [richards_activ]


return
end subroutine usppmo

!===============================================================================

!> \brief User subroutine for input of model selection parameters.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      ixmlpu       indicates if the XML file from the GUI is used
!>                              used (1: yes, 0: no
!> \param[in, out] iturb        turbulence model
!> \param[in, out] itherm       thermal model
!> \param[in, out] iale         ale module
!> \param[in, out] ivofmt       vof method
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

subroutine usipph &
!================
 ( ixmlpu, iturb , itherm, iale , ivofmt, icavit )

!===============================================================================

use entsor, only: nfecra ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, ivofmt, icavit

! Local variables

!===============================================================================

!< [richards_warning]
if (ixmlpu.eq.0) then
  ! For groundwater flow module, turbulent model is set to 0 for security reason
  iturb = 0
endif
!< [richards_warning]

return
end subroutine usipph

!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use rotation
use darcy_module
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables
integer       iscal, ifcvsl, ii, key_decay

type(var_cal_opt) :: vcopt

double precision decay_rate

!----------
! Flow part
!----------

!< [richards_iwgrec]

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! Set gradient computation to weighted (1) for high permeability ratio in
! tetrahedral meshes.
! Only works with isotropic permeability and the
! standard least square gradient reconstruction.
if (darcy_anisotropic_permeability.eq.0) then
  vcopt%iwgrec = 1
  imrgra = 1
  call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
endif
!< [richards_iwgrec]

!< [richards_num_flow]
! Set low criteria for gradient reconstruction to obtain a
! smooth velocity field

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

vcopt%nswrsm = 5000
vcopt%epsrsm = 1.d-10
vcopt%nswrgr = 5000
vcopt%epsrgr = 1.d-10
vcopt%epsilo = 1.d-13

call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

!< [richards_num_flow]

! --------------
! Transport part
!---------------

!< [richards_num_trpt]
if (darcy_anisotropic_dispersion.eq.0) then
  do ii = 1, nscal
    ! Set gradient computation to weighted (1) for high permeability ratio in tetrahedral meshes.
    ! Only works with isotropic diffusion.
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    vcopt%iwgrec = 1
    ! Set higher maximum number of iteration for reconstruction to increase the accuracy
    vcopt%nswrsm = 10
    call field_set_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
  enddo
endif
!< [richards_num_trpt]

!< [richards_decay]
! In groundwater module flow, the radioactive decay of solute is treated as a
! source term in the transport equation

! Get radioactive decay rate key
call field_get_key_id("fo_decay_rate", key_decay)
do ii = 1, nscal
  decay_rate = 3.d-1
  ! Set radioactive decay rate
  call field_set_key_double(ivarfl(isca(ii)), key_decay, decay_rate)
enddo
!< [richards_decay]

!---------------

! Set total number of iteration and reference time step
! (can be modified in cs_user_extra_operations.f90)
!< [richards_num_time]
ntmabs = 500
dtref  = 1.d-3
!< [richards_num_time]

! Initialise field id for each scalar
do iscal = 1, nscaus
  if (iscavr(iscal).le.0) then
    ifcvsl = 0
    call field_set_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
  endif
enddo

return
end subroutine usipsu

!===============================================================================

!> \brief User routine for definition of computation parameters dealing
!>        with Darcy module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine user_darcy_ini1

!===============================================================================
! Module files
!===============================================================================

use ihmpre, only: iihmpr
use entsor
use darcy_module
use cs_c_bindings
use numvar

!===============================================================================

implicit none

!===============================================================================

! Local variables

double precision, dimension(:), pointer :: kd, kplus, kminus
type(gwf_soilwater_partition) :: sorption_sca1

!----------
! Flow part
!----------

!< [richards_perm]
! Set permeability to isotropic (0) or anisotropic (1) for all soils
darcy_anisotropic_permeability = 0
!< [richards_perm]

!< [richards_grav]
! Set gravity to pass from H to h. Example for H = h + z:
darcy_gravity = 1
!< [richards_grav]

!< [richards_conv]
! Set convergence criteron of the Newton scheme over pressure (0) or over velocity (1).
! It is recommended to keep the criteron over pressure.
darcy_convergence_criterion = 0
!< [richards_conv]

!---------------
! Transport part
!---------------

!< [richards_disp]
! Set dispersion to isotropic (0) or anisotropic (1) for all solutes
darcy_anisotropic_dispersion = 0
!< [richards_disp]

!< [richards_steady]
! Set if the transport part is based on a steady (0) or unsteady (1) flow field
darcy_unsteady = 0
!< [richards_steady]

!---------------
! Sorption model
!---------------

!< [richards_partition]
! Get soil-water partition structure.
call field_get_key_struct_gwf_soilwater_partition(ivarfl(isca(1)), &
                                                  sorption_sca1)

! Set the sorption model to Kd approach (0) or EK model (1),
! Kd approach is set by default.
sorption_sca1%kinetic = 1

! Enable precipitation model, by default, there is no precipitation.
sorption_sca1%imxsol = 0 ! imxsol will hold the solubility index field id

! Set the modifications in the soil-water partition structure.
call field_set_key_struct_gwf_soilwater_partition(ivarfl(isca(1)), &
                                                  sorption_sca1)
!< [richards_partition]

return

end subroutine user_darcy_ini1

!===============================================================================
