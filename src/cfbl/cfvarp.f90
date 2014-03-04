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

subroutine cfvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!              INIT DES POSITIONS DES VARIABLES
!            POUR LE COMPRESSIBLE SANS CHOC SELON
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
use ihmpre

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_compressible()  &
    bind(C, name='cs_field_pointer_map_compressible')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_field_pointer_map_compressible

  subroutine cs_gui_labels_compressible()  &
    bind(C, name='cs_gui_labels_compressible')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_labels_compressible

end interface

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================


if (ippmod(icompf).ge.0) then

  ! Total energy

  itherm = 3
  call add_model_scalar_field('total_energy', 'EnergieT', ienerg)
  iscalt = ienerg

  ! Alias for B.C.
  irunh = ienerg

  ! Temperature (post)
  call add_model_scalar_field('temperature', 'TempK', itempk)

! ---- Viscosite dynamique de reference relative au scalaire ITEMPK
  ivisls(itempk) = 0
  visls0(itempk) = epzero

! ---- Viscosite dynamique de reference relative au scalaire IENERG
  ivisls(ienerg) = 0
  visls0(ienerg) = epzero

! ---- Initialisation par defaut de la viscosite en volume (cste)
  iviscv = 0
  viscv0 = 0.d0

! MAP to C API
call cs_field_pointer_map_compressible

! Mapping for GUI
if (iihmpr.eq.1) then
  call cs_gui_labels_compressible
endif

!===============================================================================
! 2. OPTIONS DE CALCUL
!===============================================================================

! --> Cv constant ou variable (par defaut : constant)
  icv = 0
  cv0 = 0.d0

  call cf_set_thermo_options

!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTILISATEUR
!===============================================================================

  if (iihmpr.eq.1) then
    call csvvva(iviscv)
  endif

endif

!--------
! Formats
!--------


return
end subroutine

