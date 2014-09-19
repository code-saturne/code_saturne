!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine clpalp &
!================

 ( ncelet , ncel   , nvar   )

!===============================================================================
! FONCTION :
! ----------

! Clipping of Alpha

!-------------------------------------------------------------------------------
! Arguments
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! nvar             ! e  ! <-- ! nombre de variables                            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!-------------------------------------------------------------------------------
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use numvar
use cstnum
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel

! VARIABLES LOCALES

integer          iel
integer          iclpmn, iclpmx
double precision vmin(1), vmax(1), var

double precision, dimension(:), pointer :: cvar_al

!===============================================================================

!===============================================================================
!  ---> Stockage Min et Max pour listing
!===============================================================================

call field_get_val_s(ivarfl(ial), cvar_al)

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = cvar_al(iel)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

iclpmn = 0
iclpmx = 0
do iel = 1, ncel
  if (cvar_al(iel).lt.0.d0) then
    iclpmn = iclpmn + 1
    cvar_al(iel) = 0.d0
  elseif(cvar_al(iel).gt.1.d0) then
    iclpmx = iclpmx + 1
    cvar_al(iel) = 1.d0
  endif
enddo

call log_iteration_clipping_field(ivarfl(ial), iclpmn, iclpmx, vmin, vmax)

return

end
