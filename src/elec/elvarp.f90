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

subroutine elvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!      INIT DES POSITIONS DES VARIABLES POUR LE MODULE ELECTRIQUE
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
use ppincl
use elincl
use ihmpre
use field

!===============================================================================

implicit none

! Local variables

integer        iesp , idimve, isc
integer        f_id
integer        kscmin, kscmax

character(len=80) :: f_name, f_label
integer(c_int) :: n_gasses

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_electric_arcs(n_gasses)  &
    bind(C, name='cs_field_pointer_map_electric_arcs')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value        :: n_gasses
  end subroutine cs_field_pointer_map_electric_arcs

  subroutine cs_gui_labels_electric_arcs(n_gasses)  &
    bind(C, name='cs_gui_labels_electric_arcs')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value        :: n_gasses
  end subroutine cs_gui_labels_electric_arcs

end interface

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! 1.0 Dans toutes les versions electriques
! ========================================

! Thermal model

itherm = 2
call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
iscalt = ihm
!iscacp(iscalt)   = 0

! Real potential
call add_model_scalar_field('elec_pot_r', 'POT_EL_R', ipotr)

f_id = ivarfl(isca(ipotr))
call field_set_key_double(f_id, kscmin, -grand)
call field_set_key_double(f_id, kscmax, +grand)

! 1.1 Effet Joule (cas potentiel imaginaire)
! ==========================================

if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ! Imaginary potential
  call add_model_scalar_field('elec_pot_i', 'POT_EL_I', ipoti)
  f_id = ivarfl(isca(ipoti))
  call field_set_key_double(f_id, kscmin, -grand)
  call field_set_key_double(f_id, kscmax, +grand)
endif

! 1.2 Arc electrique
! ==================

if (ippmod(ielarc).ge.2) then

  ! Vector potential
  call add_model_field('vec_potential', 'POT_VEC', ndimve, ipotva)
  f_id = ivarfl(isca(ipotva))
  call field_set_key_double(f_id, kscmin, -grand)
  call field_set_key_double(f_id, kscmax, +grand)
endif

! 1.3 Conduction ionique
! ======================


! 1.4 Dans toutes les versions electriques
! ========================================

! ---- Fractions massiques des constituants

if (ngazg .gt. 1) then
  do iesp = 1, ngazg-1
    write(f_name,'(a13,i2.2)') 'esl_fraction_',iesp
    write(f_label,'(a6,i2.2)') 'YM_ESL',iesp
    call add_model_scalar_field(f_name, f_label, iycoel(iesp))
    f_id = ivarfl(isca(iycoel(iesp)))
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)
  enddo
endif

! Map to field pointers

n_gasses = ngazg

call cs_field_pointer_map_electric_arcs(n_gasses)

! Map labels for GUI

if (iihmpr.eq.1) then
  call cs_gui_labels_electric_arcs(n_gasses)
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!      IVISLS, ICP
!===============================================================================

do isc = 1, nscapp

  if (iscavr(iscapp(isc)).le.0) then

! ---- Viscosite dynamique moleculaire variable pour les
!                                              scalaires ISCAPP(ISC)
!        Pour l'enthalpie en particulier.
!        Pour le potentiel vecteur, voir plus bas
    ivisls(iscapp(isc)) = 1

  endif

enddo

! ---- "Viscosite dynamique moleculaire" = 1
!                                  pour le potentiel vecteur en Arc
if (ippmod(ielarc).ge.2) then
  ivisls(ipotva) = 0
endif

! ---- Cp est variable ; pas sur que ce soit indispensable pour le verre
!                              mais pour le moment c'est comme ca.
icp = 1

return
end subroutine

