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

subroutine covarp

!===============================================================================
!  FONCTION  :
!  ---------

!      INIT DES POSITIONS DES VARIABLES SELON
!              POUR LA COMBUSTION
!        FLAMME DE DIFFUSION ET DE PREMELANGE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use ppppar
use ppthch
use coincl
use ppincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! 1. Define variable fields
!===============================================================================

if (itherm .eq. 2) then
  call field_get_id_try('enthalpy', ihm)
endif

! Flamme de diffusion : chimie 3 points
if (ippmod(icod3p).ge.0) then
  call field_get_id_try('mixture_fraction', ifm)
  call field_get_id_try('mixture_fraction_variance', ifp2m)
endif

! 1.2 Flamme de diffusion : Steady laminar flamelet approach
if (ippmod(islfm).ge.0) then
  call field_get_id_try('mixture_fraction', ifm)
  if (mode_fp2m .eq. 0) then
    call field_get_id_try('mixture_fraction_variance', ifp2m)
  else if(mode_fp2m .eq. 1) then
    call field_get_id_try('mixture_fraction_2nd_moment', ifsqm)
  endif
  ! Flamelet/Progress variable model
  if (ippmod(islfm).ge.2) then
    call field_get_id_try('progress_variable', ipvm)
  endif
endif

! 1.3 Flamme de premelange : modele EBU
if (ippmod(icoebu).ge.0) then
  ! Fraction massique des gaz frais
  call field_get_id_try('fresh_gas_fraction', iygfm)
  if (ippmod(icoebu).eq.2 .or.ippmod(icoebu).eq.3) then
    call field_get_id_try('mixture_fraction', ifm)
  endif
endif

! 1.5 Flamme de premelange : modele LWC
if (ippmod(icolwc).ge.0 ) then
  call field_get_id_try('mixture_fraction', ifm)
  call field_get_id_try('mixture_fraction_variance', ifp2m)
  call field_get_id_try('mass_fraction', iyfm)
  call field_get_id_try('mass_fraction_variance', iyfp2m)
endif

! Soot mass fraction and precursor number
if (isoot.ge.1) then
  call field_get_id_try('soot_mass_fraction', ifsm)
  call field_get_id_try('soot_precursor_number', inpm)
endif

!===============================================================================

return
end subroutine
