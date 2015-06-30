!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine raysca &
 ( iisca  ,                                                       &
   ncelet , ncel   ,                                              &
   smbrs  , rovsdt , volume )

!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!       PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
!       IMPLICITE ET EXPLICITE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iisca            ! e  ! <-- ! num scalaire temperature ou enthalpie          !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! smbrs(ncelet)    ! tr ! <-- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! <-- ! tableau de travail pour terme instat           !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use cstphy
use optcal
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer          iisca , ncelet , ncel

double precision volume(ncelet)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

! Local variables

integer          iel

double precision, dimension(:), pointer :: cpro_tsri, cpro_tsre

!===============================================================================

!===============================================================================
! Radiative source terms (thermal scalar only)
!===============================================================================

if (iisca.eq.iscalt .and. (itherm.eq.1 .or. itherm.eq.2)) then

  call field_get_val_s(iprpfl(itsri(1)), cpro_tsri)
  call field_get_val_s(iprpfl(itsre(1)), cpro_tsre)

  ! Implicit part

  do iel = 1,ncel
    cpro_tsri(iel) = max(-cpro_tsri(iel),zero)
    rovsdt(iel) = rovsdt(iel) + cpro_tsri(iel)*volume(iel)
  enddo

  ! Explicit part

  do iel = 1,ncel
    smbrs(iel) = smbrs(iel) + cpro_tsre(iel)*volume(iel)
  enddo

endif

return
end subroutine
