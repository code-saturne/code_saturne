!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file rayprp.f90
!> \brief Properties definitions for radiative model.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine rayprp

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use albase
use lagpar
use lagdim
use lagran
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use ihmpre
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=80) :: f_label, f_name
character(len=16) :: e_label, e_name
integer       irphas
integer       ippok
integer       nprayc, nprprv
integer       itycat, ityloc
logical       ilved, inoprv

integer(c_int) :: n_r_phases

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_radiation(n_r_phases)  &
    bind(C, name='cs_field_pointer_map_radiation')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: n_r_phases
  end subroutine cs_field_pointer_map_radiation

  subroutine cs_gui_labels_radiation(n_r_phases)  &
    bind(C, name='cs_gui_labels_radiation')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: n_r_phases
  end subroutine cs_gui_labels_radiation

end interface

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings
ippok = 0

ilved  = .true.   ! not interleaved by default
inoprv = .false.  ! variables have no previous value

!===============================================================================
! 1. RESERVATION D'UNE PLACE DANS PROPCE SI RAYONNEMENT
!    ET LAGRANGIEN AVEC THERMIQUE DES PARTICULES
!===============================================================================

if (iirayo.gt.0) then

  nprprv = nproce

  ! --- Numeros de propriete

  call add_property_field('luminance', 'Luminance', ilumin)
  call hide_property(ilumin)
  ihisvr(ipppro(ilumin),1) = -1

  call add_property_field_nd('radiative_flux', 'Qrad', 3, iqx)
  iqy = iqx + 1
  iqz = iqy + 1
  call hide_property(iqx)
  call hide_property(iqy)
  call hide_property(iqz)
  ihisvr(ipppro(iqx),1) = -1
  ihisvr(ipppro(iqy),1) = -1
  ihisvr(ipppro(iqz),1) = -1

  do irphas = 1, nrphas

    if (irphas.gt.1) then
      write(e_name,  '("_", i2.2)') irphas
      write(e_label,  '("_", i2.2)') irphas
    else
      e_name = ""
      e_label = ""
    endif

    f_name = 'rad_st' // e_name
    f_label = 'Srad' // e_label
    call add_property_field(f_name, f_label, itsre(irphas))
    call hide_property(itsre(irphas))
    ihisvr(ipppro(itsre(irphas)),1) = -1

    f_name = 'rad_st_implicit' // e_name
    f_label = 'ITSRI' // e_label
    call add_property_field(f_name, f_label, itsri(irphas))
    call hide_property(itsri(irphas))
    ihisvr(ipppro(itsri(irphas)),1) = -1

    f_name = 'rad_absorption' // e_name
    f_label = 'Absorp' // e_label
    call add_property_field(f_name, f_label, iabso(irphas))
    call hide_property(iabso(irphas))
    ihisvr(ipppro(iabso(irphas)),1) = -1

    f_name = 'rad_emission' // e_name
    f_label = 'Emiss' // e_label
    call add_property_field(f_name, f_label, iemi(irphas))
    call hide_property(iemi(irphas))
    ihisvr(ipppro(iemi(irphas)),1) = -1

    f_name = 'rad_absorption_coeff' // e_name
    f_label = 'CoefAb' // e_label
    call add_property_field(f_name, f_label, icak(irphas))
    call hide_property(icak(irphas))
    ihisvr(ipppro(icak(irphas)),1) = -1

  enddo

  nprayc = nproce - nprprv

  ! Boundary face fields

  ilved  = .true.   ! not interleaved by default
  inoprv = .false.  ! variables have no previous value

  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 3 ! boundary faces

  call field_create('wall_temperature',  &
                    itycat, ityloc, 1, ilved, inoprv, itparo)
  call field_set_key_str(itparo, keylbl, 'Wall_temp')

  call field_create('rad_incident_flux',  &
                    itycat, ityloc, 1, ilved, inoprv, iqinci)
  call field_set_key_str(iqinci, keylbl, 'Incident_flux')

  if (imoadf.ge.1) then
    call field_create('Spectral_rad_incident_flux',  &
                      itycat, ityloc, nwsgg, ilved, inoprv, iqinsp)
    call field_set_key_str(iqinsp, keylbl, 'Spectral_incident_flux')
  endif

  call field_create('wall_thermal_conductivity',  &
                    itycat, ityloc, 1, ilved, inoprv, ixlam)
  call field_set_key_str(ixlam, keylbl, 'Th_conductivity')

  call field_create('wall_thickness',  &
                    itycat, ityloc, 1, ilved, inoprv, iepa)
  call field_set_key_str(iepa, keylbl, 'Thickness')

  call field_create('emissivity',  &
                    itycat, ityloc, 1, ilved, inoprv, ieps)
  call field_set_key_str(ieps, keylbl, 'Emissivity')

  call field_create('rad_net_flux',  &
                    itycat, ityloc, 1, ilved, inoprv, ifnet)
  call field_set_key_str(ifnet, keylbl, 'Net_flux')

  call field_create('rad_convective_flux',  &
                    itycat, ityloc, 1, ilved, inoprv, ifconv)
  call field_set_key_str(ifconv, keylbl, 'Convective_flux')

  call field_create('rad_exchange_coefficient',  &
                    itycat, ityloc, 1, ilved, inoprv, ihconv)
  call field_set_key_str(ihconv, keylbl, 'Convective_exch_coef')

! Map to field pointers

n_r_phases = nrphas
call cs_field_pointer_map_radiation(n_r_phases)

endif

!===============================================================================

return
end subroutine rayprp
