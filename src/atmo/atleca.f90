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

!> \file atleca.f90
!> \brief Reads initial aerosol concentration and number
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine atleca
!================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use atincl
use numvar
use atchem
use sshaerosol
use field

implicit none

!===============================================================================

! Arguments


! Local variables

integer    isc, f_id
integer    jsp, jb
integer    impmea
character  label*80
character  ficmea*80

!================================================================================
! ALLOCATE
!================================================================================

if (allocated(dlconc0)) deallocate(dlconc0)
allocate(dlconc0(n_aer*(1+nlayer_aer)))

!================================================================================
! READINGS
!================================================================================

write(nfecra,*) ''
write(nfecra,*) 'reading of aerosols numbers and concentrations'

if (init_aero_with_lib) then

  ! The external library provides the concentrations / numbers
  call sshaerosol_get_aero(dlconc0)

  ! Conversion from microg / m^3 to ppm
  do jb = 1, n_aer * nlayer_aer
    dlconc0(jb) = dlconc0(jb) / (1.d-3 * ro0)
  enddo

  ! Conversion from molecules / m^3 to molecules / kg
  do jb = n_aer * nlayer_aer + 1, n_aer * nlayer_aer + n_aer
    dlconc0(jb) = dlconc0(jb) / ro0
  enddo

else

  ! Read from file
  call atmo_get_aero_conc_file_name(ficmea)
  open(newunit=impmea,file=ficmea,status='old')
  ! Reading aerosol numbers
  do jb = 1, n_aer
    read(impmea,*) dlconc0(nlayer_aer*n_aer+jb)
  enddo
  ! Reading aerosol concentrations
  do jb = 1, n_aer
    do jsp = 1, nlayer_aer
      read(impmea,*) dlconc0(jb+(jsp-1)*n_aer)
    enddo
  enddo
  ! Close
  close(impmea)

endif

!================================================================================
! PRINTINGS
!================================================================================

write(nfecra, *) ''
write(nfecra, *) '==================================================='
write(nfecra, *) 'printing aerosol numbers'
do jb = 1, n_aer
  isc = isca_chem(nespg + n_aer*nlayer_aer + jb)
  f_id = ivarfl(isca(isc))
  call field_get_label(f_id, label)
  write(nfecra,1001) label, dlconc0(n_aer*nlayer_aer+jb)
enddo

write(nfecra, *) ''
write(nfecra, *) '==================================================='
write(nfecra, *) 'printing aerosol concentrations'

do jb = 1, n_aer
  write(nfecra,*) "Size bin number ",jb
  do jsp = 1, nlayer_aer
    isc = isca_chem(nespg + jb + (jsp-1)*n_aer)
    f_id = ivarfl(isca(isc))
    call field_get_label(f_id, label)
    write(nfecra,1001) label, dlconc0(jb + (jsp-1)*n_aer)
  enddo
enddo
1001 format(A20," : ",ES10.2)

end subroutine atleca

