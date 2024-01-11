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

!===============================================================================
! Function:
! --------
!>  \file  cs_coal_htconvers2.f90
!===============================================================================

!===============================================================================
!> \brief Compute particles enthalpy
!>        Function with temperature and concentrations

!> \param[in]     class_id      class id (0 to n-1)
!> \param[in]     xsolid        mass fraction of components
!> \param[in,out] temper        temperature (in kelvin)

!> \return   mass enthalpy (in \f$ j . kg^{-1}) \f$
!===============================================================================

function cs_coal_thconvers2(class_id, xsolid, temper) result(enthal) &
  bind(C, name='cs_coal_thconvers2')

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer(c_int), value :: class_id

real(c_double) :: xsolid(*)
real(c_double), value :: temper
real(c_double) :: enthal

! Local variables

integer          icla, icha, it , isol , ihflt2

double precision eh1 , eh0

!===============================================================================

icla = class_id + 1

! Conversion mode
ihflt2 = 1

if (ihflt2.eq.0) then

  ! 1. H2 linear function
  !======================

  icha = ichcor(icla)

  enthal = h02ch(icha) + cp2ch(icha)*(temper-trefth)

elseif (ihflt2.ne.0) then

  ! 2. H2 tabulated
  !================

  it = npoc
  if (temper.ge.thc(it)) then
    enthal = zero
    do isol = 1, nsolid
      enthal = enthal + xsolid(isol)*ehsoli(isol,it)
    enddo
    go to 11
  endif

  it = 1
  if (temper.le.thc(it)) then
    enthal = zero
    do isol = 1, nsolid
      enthal = enthal + xsolid(isol)*ehsoli(isol,it)
    enddo
    go to 11
  endif
  it = 1
10 continue

  it = it + 1
  if (temper.le.thc(it)) then
    eh0 = zero
    eh1 = zero
    do isol = 1, nsolid
      eh0 = eh0 + xsolid(isol)*ehsoli(isol,it-1)
      eh1 = eh1 + xsolid(isol)*ehsoli(isol,it)
    enddo
    enthal =   eh0                                                &
             + (eh1-eh0)*(temper-thc(it-1))                       &
                        /(thc(it)-thc(it-1))
    goto 11
  endif
  goto 10
11 continue

endif

!----
! End
!----

end function cs_coal_thconvers2

!===============================================================================

!> \brief Calculating temperature of particles
!>        Function with enthalpy and concentrations

!> \param[in]     icla          class number
!> \param[in]     enthal        mass enthalpy (in \f$ j . kg^{-1}) \f$
!> \param[in]     xsolid        mass fraction of components
!> \param[in,out] temper        temperature (in kelvin)
!> \param[in]     t1            coal inlet temperature

!===============================================================================

subroutine cs_coal_htconvers2 &
 (icla , enthal , xsolid , temper , t1)

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          icla , icha

double precision xsolid(nsolim)
double precision temper , enthal , t1

! Local variables

integer          it , isol , ihflt2

double precision eh1 , eh0 , x2

!===============================================================================
! Conversion mode
ihflt2 = 1

if (ihflt2.eq.0) then

  ! 1. H2 linear function
  !======================

  icha = ichcor(icla)

  temper =  (enthal-h02ch(icha))/cp2ch(icha) + trefth

elseif (ihflt2.ne.0) then

  ! 2. H2 tabulated
  !================

  ! --> Enthalpy law -> temperature (mode = 1)

  x2 = 0.d0
  do isol = 1, nsolid
    x2 = x2 + xsolid(isol)
  enddo

  if (x2 .gt. epsicp) then
    it  = npoc-1
    eh1 = zero
    do isol = 1, nsolid
      eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
    enddo
    if (enthal.ge.eh1) temper = thc(it+1)

    it  = 1
    eh0 = zero
    do isol = 1, nsolid
      eh0 = eh0 + xsolid(isol)*ehsoli(isol,it)
    enddo
    if (enthal.le.eh0) temper = thc(it)

    do it = 1, npoc-1
      eh0 = zero
      eh1 = zero
      do isol = 1, nsolid
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it)
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
      enddo
      if (enthal.ge.eh0 .and. enthal.le.eh1)                     &
        temper =  thc(it)                                        &
                 + (enthal-eh0)*(thc(it+1)-thc(it))/(eh1-eh0)
    enddo

  else
    temper = t1
  endif

endif

!----
! End
!----

return
end subroutine cs_coal_htconvers2
