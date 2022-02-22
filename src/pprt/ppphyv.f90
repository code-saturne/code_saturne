!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \brief These subroutineS fill physical properties which are variable in time
!>        for the dedicated physics modules
!>        (BEFORE and AFTER the user surbroutines).
!>
!> \warning:
!>  - it is forbidden to modify the turbulent viscosity here.
!>  - \ref cstphy::icp "icp" must be set to 1 if one wants the specific heat
!>    to be variable in space
!>  - it is necessary to call field_set_key_int(ivarfl(isca(iscal)), kivisl, 0)
!>    if one wants the specific heat to be variable in space
!> \remarks:
!>  - this routine is called at the begining of each time step,
!>    thus, at the first time step, the only initialized variables are:
!>     - in usipsu.f90:
!>         - the density (set at ro0)
!>         - the molecular viscosity (set to viscl0)
!>     - in usppiv.f90:
!>         - the variables (0 by default)
!>
!> Here can be given:
!>  - the cell density in kg/m3
!>    (and eventually the boundary density)
!>  - the dynamic molecular viscosity in kg/(m s)
!>  - the specific heat Cp in J/(kg degres)
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


subroutine cs_physical_properties1 &
 ( mbrom  )

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
use mesh

use field

!===============================================================================

implicit none

! Arguments

integer          mbrom

! Local variables

!===============================================================================

!===============================================================================
! 1. Fill properties depending on the model
!===============================================================================

! If combustion drift activated
if (i_comb_drift.ge.1) then
  call cs_physical_properties_combustion_drift()
endif

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
  call cs_coal_physprop(mbrom, izfppp)
endif

! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

if (ippmod(icpl3c).ge.0) then
  call cplphy(mbrom, izfppp)
endif

! ---> Flamme fuel

if (ippmod(icfuel).ge.0) then
  call cs_fuel_physprop(mbrom, izfppp)
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
  call ctphyv
endif

! ---> Atmospheric Flows (except constant density: ippmod(iatmos) = 0)

if (ippmod(iatmos).ge.1) then
  call atphyv
endif

!----
! End
!----

return
end subroutine cs_physical_properties1

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine cs_physical_properties2

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
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================

interface

  !---------------------------------------------------------------------------

  ! Interface to C function defining groundwater flow physical properties.

  subroutine cs_gwf_delay_update()  &
    bind(C, name='cs_gwf_delay_update')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gwf_delay_update

  !---------------------------------------------------------------------------

end interface

!===============================================================================
! 1. Initializations
!===============================================================================

!===============================================================================
! 2. Fill properties depending on the model
!===============================================================================

! Gas mixture modelling in presence of noncondensable gases or/and
! condensable gas as stream.
if (ippmod(igmix).ge.0) then
  call cs_gas_mix_physical_properties
endif

! Compressible
! Has to be called after the gas mix physical properties. The isobaric specific
! heat (Cp) is computed in cs_gas_mix_physical_properties and used in cfphyv to
! compute the isochoric specific heat, if gas mix specific physics is enabled.
if (ippmod(icompf).ge.0) then
  call cfphyv
endif

! Groundwater flows
! Has to be called after the definition of saturation and soil density in order
! to update the sorbed concentration (in case the kinetic sorption model is
! enabled) and delay.
if (ippmod(idarcy).ge.1) then
  call cs_gwf_delay_update
endif

!----
! End
!----

return
end subroutine cs_physical_properties2
