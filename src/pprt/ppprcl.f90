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

subroutine ppprcl()  &
 bind(C, name='cs_f_ppprcl')

!===============================================================================
! FONCTION :
! --------

!    PREPARATION DU REMPLISSAGE DES CONDITIONS AUX LIMITES

!           AIGUILLAGE SPECIFIQUE AUX PHYSIQUES PARTICULIERES


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use cfpoin
use atincl
use pointe, only: izfppp
use dimens, only: nvar
use cs_c_bindings
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ifac, izone, ivar

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

! ---> Combustion gaz USEBUC
!      Flamme de diffusion : chimie 3 points

if (ippmod(icod3p).ge.0 .or. ippmod(islfm).ge.0) then

  do izone = 1, nozppm
    qimp(izone)   = zero
    iqimp(izone)  = 0
    ientox(izone) = 0
    ientfu(izone) = 0
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion gaz USEBUC
!      Flamme de premelange : modele EBU

elseif ( ippmod(icoebu).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    qimp(izone)   = zero
    icalke(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    fment(izone)  = zero
    tkent(izone)  = zero
    ientgf(izone) = 0
    ientgb(izone) = 0
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Compressible

elseif (ippmod(icompf).ge.0) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

elseif (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
