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

!===============================================================================
! Function:
! --------
!> \file cs_fuel_physprop2.f90
!>
!> \brief Calculation of the physical properties of the dispersed phase
!>
!>
!> Cell values
!> -----------
!> - Mass fraction of liquid
!>   and eventual clippings
!> - Diameter
!> - Mass density
!>   and eventual clippings
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!______________________________________________________________________________!

subroutine cs_fuel_physprop2 &
 ( ncelet , ncel )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

! Local variables

integer          iel
integer          n1     , n2     , icla

double precision xnp    ,  d1s3
double precision diam2m , diam2x

double precision, dimension(:), pointer :: cvar_yfolcl, cvar_ngcl
double precision, dimension(:), pointer :: cpro_diam2, cpro_rom2

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

d1s3 = 1.d0/3.d0

!===============================================================================
! 2. Calculation for each class
!    - of mass density of fol
!    - of mass fraction of fol
!    - of the coke diameter
!===============================================================================

do icla = 1, nclafu

  n1 = 0
  n2 = 0
  diam2m =  1.d0
  diam2x =  0.d0

  call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl)
  call field_get_val_s(ivarfl(isca(ing(icla))), cvar_ngcl)
  call field_get_val_s(idiam2(icla), cpro_diam2)
  call field_get_val_s(irom2(icla), cpro_rom2)

  do iel = 1, ncel

    !  Mass density
    cpro_rom2(iel) = rho0fl

    ! --- Calculation of the particle diameter

    yfol   = cvar_yfolcl(iel)
    xnp    = cvar_ngcl(iel)
    if ( yfol .gt. epsifl .and. (xnp*yfol) .gt. 0.d0) then

      cpro_diam2(iel) = ( (yfol / cpro_rom2(iel) )           &
                            /(pi/6.d0 * xnp) ) ** d1s3

      if ( cpro_diam2(iel) .gt. dinifl(icla) ) then
        n1 = n1+1
        diam2x = max(diam2x,cpro_diam2(iel))
        cpro_diam2(iel) = dinifl(icla)
      endif

      if ( cpro_diam2(iel) .lt. diniin(icla) ) then
        n2 = n2+1
        diam2m = min(diam2m,cpro_diam2(iel))
        cpro_diam2(iel) = diniin(icla)
      endif

    else
      cpro_diam2(iel) = dinifl(icla)
    endif

  enddo

  if (irangp.ge.0) then

    call parcpt (n1)
    !==========
    call parcpt (n2)
    !==========

    call parmax (diam2x)
    !==========
    call parmin (diam2m)
    !==========
  endif

  if ( n1 .gt. 0 ) then
    write(nfecra,1001) icla, n1, diam2x
  endif
  if ( n2 .gt. 0 ) then
    write(nfecra,1002)  icla, n2, diam2m
  endif

enddo

!--------
! Formats
!--------
!----

 1001 format(/,1X,' clipping in max of class diameter:',I2,           &
       /,10X,' Number of points: ',I8,                           &
       /,10X,' Max value: ',G15.7)
 1002 format(/,1X,' clipping in min of class diametre:',I2,           &
       /,10X,' Number of points: ',I8,                           &
       /,10X,' Min value: ',G15.7)


!----
! End
!----

return
end subroutine
