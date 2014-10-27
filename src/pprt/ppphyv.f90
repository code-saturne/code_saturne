!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
!>  - \ref icp must be set to 1 if one wants the specific heat to be variable
!>    in space
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
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of prescribed density at the boundary
!> \param[in]     dt            time step (per cell)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________


subroutine cs_physical_properties1 &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     , propce )

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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

!===============================================================================
! 2. Fill properties depending on the model
!===============================================================================

! ---> Flamme de diffusion chimie 3 points

  if (ippmod(icod3p).ge.0) then
    call d3pphy(mbrom, izfppp, propce)
  endif

! ---> Flamme de diffusion chimie equilibre

!        if ( ippmod(icodeq).ge.0 )

! ---> Flamme de premelange : Modele EBU

  if (ippmod(icoebu).ge.0) then
    call ebuphy(mbrom, izfppp, propce)
  endif

! ---> Flamme de premelange : Modele BML

!        if ( ippmod(icobml).ge.0 )
!     &     call bmlphy

! ---> Flamme de premelange : Modele LWC

  if (ippmod(icolwc).ge.0) then
    call lwcphy(mbrom, izfppp, propce)
  endif

! ---> Flamme charbon pulverise

   if (ippmod(iccoal).ge.0) then
     call cs_coal_physprop(mbrom, izfppp, propce)
   endif


! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

  if (ippmod(icpl3c).ge.0) then
    call cplphy(mbrom, izfppp, propce)
  endif

! ---> Flamme fuel

  if (ippmod(icfuel).ge.0) then
    call cs_fuel_physprop(mbrom, izfppp, propce)
  endif

! ---> Physique particuliere : Versions electriques
!          Effet Joule
!          Arc electrique
!          Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

!     En Joule, on impose a l'utilisateur de programmer ses lois
!        sur les proprietes (masse volumique , ...)
!        Des exemples physiques sont fournis dans uselph.
!     En arc electrique, on lit un fichier de donnees et on interpole.

  call elphyv &
  !==========
 ( nvar   , nscal  ,                                              &
   mbrom  , izfppp ,                                              &
   dt     )

endif

! ---> Aerorefrigerants

if (ippmod(iaeros).ge.0) then
   call ctphyv
endif

! ---> Atmospheric Flows (except constant density: ippmod(iatmos) = 0)

if (ippmod(iatmos).ge.1) then
   call atphyv(propce)
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
!> \param[in]     nvar          total number of variables
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine cs_physical_properties2 &
 ( propce )

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

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)

! Local variables

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

!===============================================================================
! 2. Fill properties depending on the model
!===============================================================================

! Compressible
if (ippmod(icompf).ge.0) then

  call cfphyv(propce)

endif

! Mixing gases modelling in presence of noncondensable gases or/and
! condensable gas as stream.
if (ippmod(imixg).ge.0) then

  call cs_mixing_gas_physical_properties

endif

!----
! End
!----

return
end subroutine cs_physical_properties2
