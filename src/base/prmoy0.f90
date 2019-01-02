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

subroutine prmoy0 &
!================

 ( ncelet , ncel   , cell_f_vol , pvar   )

!===============================================================================
! FONCTION :
! ----------

! RECALAGE DE LA VARIABLE PVAR (PRESSION)POUR OBTENIR
! QU'ELLE SOIT A MOYENNE NULLE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! cell_f_vol       ! tr ! <-- ! fluid volume of cells                          !
! pvar             ! tr ! <-- ! tableau de valeurs au cellules                 !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

double precision cell_f_vol(ncelet), pvar(ncelet)

! Local variables

integer          iel
double precision pmoy

!===============================================================================

pmoy = 0.0d0
do iel = 1, ncel
  pmoy = pmoy + cell_f_vol(iel) * pvar(iel)
enddo
if (irangp.ge.0) then
  call parsom (pmoy)
  !==========
endif

pmoy = pmoy / voltot
do iel = 1, ncel
  pvar(iel) = pvar(iel) - pmoy + pred0
enddo

!----
! FIN
!----

return

end subroutine
