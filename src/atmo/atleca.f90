!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
use siream
use field

implicit none

!===============================================================================

! Arguments


! Local variables

integer    isc, f_id
integer    jsp, jb
character  label*80

!================================================================================
! READINGS
!================================================================================

write(nfecra,*) ''
write(nfecra,*) 'reading of aerosols numbers and concentrations'

open(impmea,file=ficmea,status='old')
! Reading aerosol numbers
do jb = 1, nbin_aer
  read(impmea,*) dlconc0(nesp_aer*nbin_aer+jb)
enddo
! Reading aerosol concentrations
do jb = 1, nbin_aer
  do jsp = 1, nesp_aer
    read(impmea,*) dlconc0(jb+(jsp-1)*nbin_aer)
  enddo
enddo
close(impmea)

!================================================================================
! PRINTINGS
!================================================================================

write(nfecra, *)
write(nfecra, *) '==================================================='
write(nfecra, *) 'printing aerosol numbers'
do jb = 1, nbin_aer
  write(nfecra,1000) jb, dlconc0(nesp_aer*nbin_aer+jb)
enddo
1000 format("Bin ",I2," : ",ES10.2)
write(nfecra, *)
write(nfecra, *) '==================================================='
write(nfecra, *) 'printing aerosol concentrations'

do jb = 1, nbin_aer
  write(nfecra,*) "Bin ",jb
  do jsp = 1, nesp_aer
    isc = (isca_chem(1) - 1) + nespg_siream+jb+(jsp-1)*nbin_aer
    f_id = ivarfl(isca(isc))
    call field_get_label(f_id, label)
    write(nfecra,1001) label, dlconc0(jb+(jsp-1)*nbin_aer)
  enddo
enddo
1001 format(A10," : ",ES10.2)

end subroutine

