!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> \file usthht.f90
!>
!> \brief Enthalpy-temperature conversion law definition.
!>
!>  Enthalpy    -> Temperature law (mode =  1) \n
!>  Temperature -> Enthalpy    law (mode = -1)
!>
!> See \subpage us_thht for examples.

subroutine usthht &
!================

 ( mode   , enthal , temper  )

!===============================================================================
!> \brief Enthalpy-temperature conversion law definition.
!>
!>  Enthalpy    -> Temperature law (mode =  1) \n
!>  Temperature -> Enthalpy    law (mode = -1)
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     mode          -1 : t -> h  ;   1 : h -> t
!> \param[in]     enthal        enthalpie
!> \param[in]     temper        temperature
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          mode

double precision enthal, temper

! Local variables

integer         it

!< [loc_var_dec]

! For the second example below,
!   input of ntab > 1 (fictitious) tabulated values of H (=HT)
!   as a function of NTAB values of T (=TH) (be careful with unit K or C)

integer          ntab
parameter       (ntab=5)
double precision ht(ntab), th(ntab)
data             ht /100000.d0,200000.d0,300000.d0,               &
                               400000.d0,500000.d0 /
data             th /   100.d0,   200.d0,   300.d0,               &
                                  400.d0,   500.d0 /

!< [loc_var_dec]

!===============================================================================


!     WARNING, the subroutine is called in loops:
!     =======                    ===============

!       Avoid parallel (collective) operations
!       =====

!===============================================================================
! Examples
!===============================================================================

!< [example_1]

! First example, corresponding to H = CpT with Cp = 4000
! ======================================================

! --- Mode H -> T
if (mode .eq.  1) then
  temper = enthal / 4000.d0

! --- Mode T -> H
else
  enthal = temper * 4000.d0

endif

!< [example_1]

return

!< [example_2]

! Second example, corresponding to a simple interpolation
! based on a tabulation H = f(T) entered as DATA
! ======================================================

! --- Mode H -> T
if (mode .eq.  1) then

  ! Default initialization
  temper = 0.d0

  ! If H is smaller than the smallest tabulated value,
  ! set T to smallest temperature.
  if (enthal.le.ht(1)) then
    temper = th(1)

  ! If H is larger than the largest tabulated value,
  ! set T to largest temperature.
  else if (enthal.ge.ht(ntab)) then
    temper = th(ntab)

  ! Otherwise, use piecewise linear interpolation
  else
    do it = 2, ntab
      if(enthal.le.ht(it)) then
        temper = th(it-1)                                         &
          +(enthal-ht(it-1))*(th(it)-th(it-1))/(ht(it)-ht(it-1))
      endif
    enddo
  endif

! --- Mode T -> H
else

  ! Default initialization
  enthal = 0.d0

  ! If T is smaller than the smallest tabulated value,
  ! set H to smallest enthalpy.
  if (temper.le.th(1)) then
    enthal = ht(1)

  ! If T is larger than the largest tabulated value,
  ! set H to largest enthalpy.
  else if (temper.ge.th(ntab)) then
    enthal = ht(ntab)

  ! Otherwise, use piecewise linear interpolation
  else
    do it = 2, ntab
      if(temper.le.th(it)) then
        enthal = ht(it-1)                                         &
          +(temper-th(it-1))*(ht(it)-ht(it-1))/(th(it)-th(it-1))
      endif
    enddo
  endif
endif

!< [example_2]

return

!----
! End
!----

end subroutine usthht