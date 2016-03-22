!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine usppmo &
!================
 ( ixmlpu )

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

subroutine usipph &
!================
 ( ixmlpu, iturb , itherm, iale , icavit )

!===============================================================================

use entsor, only: nfecra ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, icavit

!< [richards_warning]
if (ixmlpu.eq.0) then
  ! For groundwater flow module, turbulent model is set to 0 for security reason
  iturb = 0
endif
!< [richards_warning]

return
end subroutine usipph

!===============================================================================

subroutine usipsu &
!================

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

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables
integer       iscal, ifcvsl, ii

!----------
! Flow part
!----------

!< [richards_iwgrec]
! Set gradient computation to weighted (1) for high permeability ratio in tetrahedral meshes.
! Only works with isotropic permeability and the standard least square gradient reconstruction.
if (darcy_anisotropic_permeability.eq.0) then
  iwgrec(ipr) = 1
  imrgra = 1
endif
!< [richards_iwgrec]

!< [richards_num_flow]
! Set low criteria for gradient reconstruction to obtain a smooth velocity field
nswrsm(ipr) = 5000
epsrsm(ipr) = 1.d-10
nswrgr(ipr) = 5000
epsrgr(ipr) = 1.d-10
epsilo(ipr) = 1.d-13
!< [richards_num_flow]

! --------------
! Transport part
!---------------

!< [richards_num_trpt]
if (darcy_anisotropic_dispersion.eq.0) then
  do ii = 1, nscal
    ! Set gradient computation to weighted (1) for high permeability ratio in tetrahedral meshes.
    ! Only works with isotropic diffusion.
    iwgrec(isca(ii)) = 1
    ! Set higher maximum number of iteration for reconstruction to increase the accuracy
    nswrsm(isca(ii)) = 10
  enddo
endif
!< [richards_num_trpt]

!---------------

! Set total number of iteration and reference time step (can be modified in cs_user_extra_operations.f90)
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

subroutine user_darcy_ini1
!========================

use ihmpre, only: iihmpr
use entsor
use darcy_module

!===============================================================================

implicit none

!===============================================================================

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

return

end subroutine user_darcy_ini1
