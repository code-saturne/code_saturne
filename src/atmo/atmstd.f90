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
!> \file atmstd.f90
!> \brief Compute standard atmospheric profile

!> \brief compute standard atmospheric profile (Holton p 374)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]       z          absolute altitude in m
!> \param[out]      p          pressure in pa
!> \param[out]      t          temperature in k
!> \param[out]      r          density in kg/m3
!-------------------------------------------------------------------------------
subroutine atmstd(z,p,t,r)

use cstphy, only:rair

implicit none

double precision z,t,p,r
double precision p0,t0,zt,g,a,t11,p11

p0 = 101325.d0
t0 = 288.15d0
zt = 11000.d0
g = 9.81d0
a = -6.5d-03

if(z.le.zt)then
  t = t0+a*z
  p = p0*(t/t0)**(-g/rair/a)
  r = p/rair/t
else
  t11 = t0 + a*zt
  p11 = p0*(t11/t0)**(-g/rair/a)
  t = t11
  p = p11*exp(-g/rair/t11*(z - zt))
  r = p/rair/t
endif

return
end subroutine atmstd
