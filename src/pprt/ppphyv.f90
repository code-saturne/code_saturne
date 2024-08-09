!-------------------------------------------------------------------------------

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
! Function:
! ---------

!> \file ppphyv.f90
!>
!> \brief These subroutines fill physical properties which are variable in time
!>        for the dedicated physics modules
!>        (BEFORE and AFTER the user subroutines).
!>
!> \warning:
!>  - it is forbidden to modify the turbulent viscosity here.
!>  - \ref cstphy::icp "icp" must be set to 1 if one wants the specific heat
!>    to be variable in space
!>  - it is necessary to call field_set_key_int(ivarfl(isca(iscal)), kivisl, 0)
!>    if one wants the specific heat to be variable in space
!> \remarks:
!>  - this routine is called at the beginning of each time step,
!>    thus, at the first time step, the only initialized variables are:
!>     - in cs_user_parameters:
!>         - the density (set at ro0)
!>         - the molecular viscosity (set to viscl0)
!>
!> Here can be given:
!>  - the cell density in kg/m3
!>    (and eventually the boundary density)
!>  - the dynamic molecular viscosity in kg/(m s)
!>  - the specific heat Cp in J/(kg degrees)
!>  - scalars diffusivities in kg/(m s)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     mbrom         indicator of prescribed density at the boundary
!_______________________________________________________________________________

subroutine cs_physical_properties1(mbrom) &
  bind(C, name='cs_f_physical_properties1')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_c_bindings
use mesh
use, intrinsic :: iso_c_binding
use field

!===============================================================================

implicit none

! Arguments

integer          mbrom

! Local variables

!===============================================================================
! 1. Fill properties depending on the model
!===============================================================================

! ---> Flamme de diffusion chimie 3 points

if (ippmod(icod3p).ge.0) then
  call d3pphy()
endif

! ---> Diffusion flame steady laminar flamelet approach
! ---> Obtain physical property values from the preprocessed look-up table
if (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_physical_prop(mbrom, izfppp)
endif

! ---> Flamme de premelange : Modele EBU

if (ippmod(icoebu).ge.0) then
  call ebuphy(mbrom, izfppp)
endif

! ---> Flamme de premelange : Modele LWC

if (ippmod(icolwc).ge.0) then
  call lwcphy(mbrom, izfppp)
endif

! ---> Flamme charbon pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_physprop(mbrom)
endif

! ---> Physique particuliere : Versions electriques
!          Effet Joule
!          Arc electrique
!          Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1       ) then
!     En Joule, on impose a l'utilisateur de programmer ses lois
!        sur les proprietes (masse volumique , ...)
!        Des exemples physiques sont fournis dans cs_user_physical_properties.
!     En arc electrique, on lit un fichier de donnees et on interpole.
  call elphyv
endif

! ---> Aerorefrigerants

if (ippmod(iaeros).ge.0) then
  call cs_ctwr_phyvar_update(ro0, t0, p0)

! ---> Atmospheric Flows (except constant density: ippmod(iatmos) = 0)
else if (ippmod(iatmos).ge.1) then
  call atphyv
  call cs_atmo_phyvar_update()
endif

!----
! End
!----

return
end subroutine cs_physical_properties1
