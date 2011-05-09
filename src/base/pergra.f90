!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine pergra &
!================

 ( nphjmx , nphas  ,                                              &
   ju     , jv     , jw     ,                                     &
   jtytur ,                                                       &
   jr11   , jr22   , jr33   , jr12   , jr13   , jr23   )

!===============================================================================
! FONCTION :
! --------

! Recuperation de certains COMMON necessaires a PERING


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nphjmx           ! e  ! <-- ! nombre de phases max                           !
! nphas            ! i  ! <-- ! number of phases                               !
! ju, jv, jw       ! te ! --> ! numero de variable pour u, v, w                !
! jtytur           ! te ! --> ! indicateur modele de turbulence                !
! jr11...jr23      ! te ! --> ! numero de variable pour rij                    !
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
use numvar
use optcal

!===============================================================================

implicit none

! Arguments

integer          nphjmx, nphas
integer          ju,jv,jw
integer          jtytur
integer          jr11,jr22,jr33
integer          jr12,jr13,jr23

! Local variables

integer          iphas

!===============================================================================


! Recuperation des COMMON de "optcal"

jtytur = itytur

! Recuperation des COMMON de "numvar"

ju   = iu
jv   = iv
jw   = iw
if(itytur.eq.3) then
  jr11 = ir11
  jr22 = ir22
  jr33 = ir33
  jr12 = ir12
  jr13 = ir13
  jr23 = ir23
else
  jr11 = 0
  jr22 = 0
  jr33 = 0
  jr12 = 0
  jr13 = 0
  jr23 = 0
endif

!----
! FIN
!----

return
end subroutine
