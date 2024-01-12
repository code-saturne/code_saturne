!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \ref f_parameters for examples.
!>
!>   If the code_saturne GUI is used, this file is not required (but may be
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

!> \brief User subroutine for the input of additional user parameters.
!
!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use atincl
use atsoil
use atchem
use atimbr
use sshaerosol
use ctincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer           iiv

!===============================================================================
! Initialize non-standard calculation options for the atmospheric version.
!===============================================================================

!< [usati1]
!  -----------------------------------------------------------------------------
!  Atmospheric imbrication on large scale meteo (atimbr module)
!  -----------------------------------------------------------------------------
!
! --------------------------------------------------------------
! activation flag
! --------------------------------------------------------------
imbrication_flag    = .false.
imbrication_verbose = .false.

! ------------------------------------------------------------------------------
! flags for activating the cressman interpolation for the boundary conditions
! ------------------------------------------------------------------------------
cressman_u     = .true.
cressman_v     = .true.
cressman_tke   = .true.
cressman_eps   = .true.
cressman_theta = .true.
cressman_qw    = .true.
cressman_nc    = .true.

! --------------------------------------------------------------
! numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
horizontal_influence_radius = 8500.d0
vertical_influence_radius = 100.d0

! --------------------------------------------------------------

! dtchemmax: maximal time step (s) for chemistry resolution
dtchemmax = 10.0d0

!< [usati1]

!< [usatsoil]
! Example to modify the soil parameters if activated
if (iatsoil.eq.1) then
  ! Example to modify some Soil constants for minerals (from Wangara test case)
  tab_sol(4)%csol = 1.7e-5
  tab_sol(4)%rugthe = 0.0012
  tab_sol(4)%rugdyn = 0.0012
endif

! Initializing the soil table of each vertical grid
if (iatsoil.ne.0) then
  do iiv = 1, nvert
    soilvert(iiv)%albedo  = 0.25d0
    soilvert(iiv)%emissi  = 0.965d0
    soilvert(iiv)%ttsoil  = 14.77d0
    soilvert(iiv)%totwat  = 0.0043d0
    soilvert(iiv)%pressure = 1023.d0
    soilvert(iiv)%density = 1.23d0
  enddo
endif
!< [usatsoil]

!----
! Formats
!----

return
end subroutine usipsu

!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called for the following models
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!
!===============================================================================

subroutine cs_user_combustion

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
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use radiat

!===============================================================================

implicit none

!< [cs_user_combustion]

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! --- Density Relaxation
!     RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0

!===============================================================================
! 2. Physical Constants
!===============================================================================

! Dynamic Diffusion Coefficient (kg/(m s))
! Should rather be set in cs_user_parameters
call field_set_key_double(ivarfl(isca(iscalt)), kvisl0, 4.25d-5)

! -----------------------------------------------------------------------------
! 2.1 For 3 points combusution model ONLY
! -----------------------------------------------------------------------------

! Reference temperature for fuel and oxydant (K)
tinfue = 436.d0
tinoxy = 353.d0

! -----------------------------------------------------------------------------
! 2.2 For EBU-model ONLY
! -----------------------------------------------------------------------------

! cebu: EBU-model constant
cebu   = 2.5d0

! -----------------------------------------------------------------------------
! 2.3 For Libby-Williams model ONLY
! -----------------------------------------------------------------------------

! Reference velocity
vref = 60.d0
! Reference length scale
lref = 0.1d0
! Activation Temperature
ta   = 0.2d5
! Cross-over Temperature (combustion of propane)
tstar= 0.12d4

!< [cs_user_combustion]

!----
! End
!----

return
end subroutine cs_user_combustion

!===============================================================================

