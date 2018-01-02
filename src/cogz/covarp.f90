!-------------------------------------------------------------------------------

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

subroutine covarp
!================


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
use cpincl
use ppincl
use ihmpre
use field

!===============================================================================

implicit none

! Local variables

integer        f_id
integer        kscmin, kscmax, kscavr

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_gas_combustion()  &
    bind(C, name='cs_field_pointer_map_gas_combustion')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_field_pointer_map_gas_combustion

  subroutine cs_gui_labels_gas_combustion()  &
    bind(C, name='cs_gui_labels_gas_combustion')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_labels_gas_combustion

end interface

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)
call field_get_key_id("first_moment_id", kscavr)

!===============================================================================
! 1. Define variable fields
!===============================================================================

! 1.1 Flamme de diffusion : chimie 3 points
! =========================================

if (ippmod(icod3p).ge.0) then

  ! Mixture fraction and its variance

  call add_model_scalar_field('mixture_fraction', 'Fra_MEL', ifm)
  f_id = ivarfl(isca(ifm))
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  call add_model_scalar_field('mixture_fraction_variance', 'Var_FrMe', ifp2m)
  f_id = ivarfl(isca(ifp2m))
  call field_set_key_int(f_id, kscavr, ivarfl(isca(ifm)))

  ! Enthalpy

  if (ippmod(icod3p).eq.1) then
    itherm = 2
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
    ! Set min and max clipping
    f_id = ivarfl(isca(iscalt))
    call field_set_key_double(f_id, kscmin, -grand)
    call field_set_key_double(f_id, kscmax, grand)
  endif

  ! Soot mass fraction and precursor number
  if (isoot.ge.1) then

    call add_model_scalar_field('soot_mass_fraction', 'Fra_Soot', ifsm)
    f_id = ivarfl(isca(ifsm))
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)

    call add_model_scalar_field('soot_precursor_number', 'NPr_Soot', inpm)
    f_id = ivarfl(isca(inpm))
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)

  endif

endif

! 1.2 Flamme de premelange : modele EBU
! =====================================

if (ippmod(icoebu).ge.0) then

  ! Fraction massique des gaz frais
  call add_model_scalar_field('fresh_gas_fraction', 'Fra_GF', iygfm)

  f_id = ivarfl(isca(iygfm))
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  if (ippmod(icoebu).eq.2 .or.ippmod(icoebu).eq.3) then
    ! Taux de melange
    call add_model_scalar_field('mixture_fraction', 'Fra_MEL', ifm)
    f_id = ivarfl(isca(ifm))
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)
  endif
  if (ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3) then
    itherm = 2
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
    f_id = ivarfl(isca(iscalt))
    call field_set_key_double(f_id, kscmin, -grand)
    call field_set_key_double(f_id, kscmax,  grand)
  endif
endif

! 1.3 Flamme de premelange : modele BML A DEVELOPPER
! ==================================================

! 1.4 Flamme de premelange : modele LWC
! =====================================

if (ippmod(icolwc).ge.0 ) then

  call add_model_scalar_field('mixture_fraction', 'Fra_MEL', ifm)
  f_id = ivarfl(isca(ifm))
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  call add_model_scalar_field('mixture_fraction_variance', 'Var_FrMe', ifp2m)
  f_id = ivarfl(isca(ifp2m))
  call field_set_key_int(f_id, kscavr, ivarfl(isca(ifm)))


  call add_model_scalar_field('mass_fraction', 'Fra_Mas', iyfm)
  f_id = ivarfl(isca(iyfm))
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  call add_model_scalar_field('mass_fraction_variance', 'Var_FMa', iyfp2m)
  f_id = ivarfl(isca(iyfp2m))
  call field_set_key_int(f_id, kscavr, ivarfl(isca(iyfm)))

  if (ippmod(icolwc).ge.2 ) then
    call add_model_scalar_field('mass_fraction_covariance', 'COYF_PP4', icoyfp)
    f_id = ivarfl(isca(icoyfp))
    call field_set_key_double(f_id, kscmin, -0.25d0)
    call field_set_key_double(f_id, kscmax, 0.25d0)
  endif

  if (ippmod(icolwc).eq.1 .or. &
      ippmod(icolwc).eq.3 .or. &
      ippmod(icolwc).eq.5) then
    itherm = 2
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
  endif

endif

! MAP to C API
call cs_field_pointer_map_gas_combustion

! Mapping for GUI
if (iihmpr.eq.1) then
  call cs_gui_labels_gas_combustion
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      ICP
!===============================================================================

if (ippmod(icod3p).eq.1 .or.   &
    ippmod(icoebu).eq.1 .or.   &
    ippmod(icoebu).eq.3 .or.   &
    ippmod(icolwc).eq.1 .or.   &
    ippmod(icolwc).eq.3 .or.   &
    ippmod(icolwc).eq.5) then

  ! Although we are in enthalpy formulation, we keep Cp constant

  icp = -1

endif

return
end subroutine
