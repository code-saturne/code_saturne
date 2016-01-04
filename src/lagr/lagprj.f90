!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine lagprj &
!================

 ( isens ,                                                        &
   vrgx   , vrgy , vrgz ,                                         &
   vrlx   , vrly , vrlz ,                                         &
   a11 , a21 , a31 , a12 , a22 , a32, a13 , a23 , a33 )

!===============================================================================

! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   ------------------------------------------------------

!   Modification of the coordinate system.
!    -If isens = 1: from global to local
!    -If isens = 2: from local to global
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! isens            ! i  ! <-- ! "direction" of the modification                !
! vrgx             ! r  ! <-> ! x-coordinate in the global coordinate system   !
! vrgy             ! r  ! <-> ! y-coordinate in the global coordinate system   !
! vrgz             ! r  ! <-> ! z-coordinate in the global coordinate system   !
! vrlx             ! r  ! <-> ! x-coordinate in the local coordinate system    !
! vrly             ! r  ! <-> ! y-coordinate in the local coordinate system    !
! vrlz             ! r  ! <-> ! z-coordinate in the local coordinate system    !
! a_(ij)           ! r  ! <-- ! coefficients of the transformation matrix      !
!------------------------------------------------------------------------------!
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array


!===============================================================================

!===============================================================================
!     Common blocks
!===============================================================================

use paramx
use entsor

!===============================================================================

implicit none

!-------------------------------------------------------------------------------

! Arguments
integer          isens

double precision vrgx,vrgy,vrgz
double precision vrlx,vrly,vrlz
double precision a11,a21,a31
double precision a12,a22,a32
double precision a13,a23,a33

! Local variables


!===============================================================================

!===============================================================================
! 1. Projection from the global reference frame to the local reference frame
!                             ( global ---> local )
!===============================================================================

if (isens.eq.1) then

  vrlx = a11*vrgx + a21*vrgy + a31*vrgz
  vrly = a12*vrgx + a22*vrgy + a32*vrgz
  vrlz = a13*vrgx + a23*vrgy + a33*vrgz

!===============================================================================
! 2. Projection from the local reference frame to the local reference frame
!                              ( local ---> global )
!===============================================================================

else if (isens.eq.2) then

  vrgx = a11*vrlx + a12*vrly + a13*vrlz
  vrgy = a21*vrlx + a22*vrly + a23*vrlz
  vrgz = a31*vrlx + a32*vrly + a33*vrlz

else

  write(nfecra,1000)
  call csexit(1)

endif

!=======
! FORMAT
!=======

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    VALEUR DE ISENS DANS LAGPRJ INCOHERENTE                ',/, &
'@                                                           ',/, &
'@       VALEURS ADMISES : 1 OU 2                             ',/,&
'@       VALEUR DETECTEE : ',I6                               ,/, &
'@                                                           ',/, &
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
