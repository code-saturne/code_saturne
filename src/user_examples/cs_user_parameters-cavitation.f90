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

!===============================================================================

!> \file cs_user_parameters-cavitation.f90
!>
!> \brief Cavitation parameters example.
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
!> \param[in, out] iale         ALE module
!> \param[in, out] ivofmt       vof method
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

subroutine usipph &
 ( ixmlpu, iturb , itherm, iale , ivofmt, icavit )

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, ivofmt, icavit

! Local variables

!===============================================================================

!>    In this subroutine, only the parameters which already appear may
!>    be set, to the exclusion of any other.
!>
!>    If we are not using the Code_Saturne GUI:
!>    All the parameters which appear in this subroutine must be set.
!>
!>    If we are using the Code_Saturne GUI:
!>    parameters protected by a test of the form:
!>
!>      if (ixmlpu.eq.0) then
!>         ...
!>      endif
!>
!>    should already have been defined using the GUI, so only
!>    experts should consider removing the test and adapting them here.

!===============================================================================

! --- Cavitation module
!    - -1: module not activated
!    -  0: no vaporization/condensation model
!    -  1: Merkle's model
!
!  Specific cavitation module input parameters should be set usipsu
!

icavit = 1

!----
! Formats
!----

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
use ihmpre
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use vof

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, jj, imom, kscmin, kscmax

!===============================================================================

!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

!-------------------------------------------------------------------------------
!>
!> \file cs_user_parameters-cavitation.f90
!>
!>
!>   Cavitation module:
!>   ------------------
!>
!>      Model based on a homogeneous mixture. The physical properties (density
!>        and dynamic viscosity) of the mixture depends on a resolved void
!>        fraction and constant reference properties of the liquid phase and
!>        the gas phase.
!>
!>    User example subroutine for cavitation module parameter definitions
!>
!-------------------------------------------------------------------------------


!===============================================================================
! 1. Definition of the homogeneous mixture physical properties
!===============================================================================

! --- Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
!       liquid phase

!< [phprop_l]
rho1 = 1.d3
mu1 = 1.d-3
!< [phprop_l]

! --- Reference density, in kg/m3, and molecular viscosity, kg/(m s), of the
!       gas phase

!< [phprop_g]
rho2 = 1.d0
mu2 = 1.d-5
!< [phprop_g]

!===============================================================================
! 2. Model parameters of the vaporization term (Merkle model)
!===============================================================================

! ---  Reference saturation pressure in kg/(m s2)

!< [presat]
presat = 2.d3
!< [presat]

! ---  Reference length, in meters, and velocity scales, in m/s, of the flow

!< [scales_inf]
linf = 0.1d0
uinf = 1.d0
!< [scales_inf]

!===============================================================================
! 3. Interaction with turbulence
!===============================================================================

! --- Eddy-viscosity correction (Reboud et al correction)
!      0: deactivated
!      1: activated

!< [reboud_activ]
icvevm = 1
!< [reboud_activ]

!----
! End
!----

end subroutine usipsu
