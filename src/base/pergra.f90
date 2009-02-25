!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nphjmx           ! e  ! <-- ! nombre de phases max                           !
! nphas            ! e  ! <-- ! nombre de phases                               !
! ju, jv, jw       ! te ! --> ! numero de variable pour u, v, w                !
! jtytur           ! te ! --> ! indicateur modele de turbulence                !
! jr11...jr23      ! te ! --> ! numero de variable pour rij                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"

!===============================================================================

! Arguments

integer          nphjmx, nphas
integer          ju(nphjmx),jv(nphjmx),jw(nphjmx)
integer          jtytur(nphjmx)
integer          jr11(nphjmx),jr22(nphjmx),jr33(nphjmx)
integer          jr12(nphjmx),jr13(nphjmx),jr23(nphjmx)

! VARIABLES LOCALES

integer          iphas

!===============================================================================


! Recuperation des COMMON de "optcal"

do iphas = 1, nphas
  jtytur(iphas) = itytur(iphas)
enddo

! Recuperation des COMMON de "numvar.h"

do iphas = 1, nphas
  ju  (iphas) = iu  (iphas)
  jv  (iphas) = iv  (iphas)
  jw  (iphas) = iw  (iphas)
  if(itytur(iphas).eq.3) then
    jr11(iphas) = ir11(iphas)
    jr22(iphas) = ir22(iphas)
    jr33(iphas) = ir33(iphas)
    jr12(iphas) = ir12(iphas)
    jr13(iphas) = ir13(iphas)
    jr23(iphas) = ir23(iphas)
  else
    jr11(iphas) = 0
    jr22(iphas) = 0
    jr33(iphas) = 0
    jr12(iphas) = 0
    jr13(iphas) = 0
    jr23(iphas) = 0
  endif
enddo

!----
! FIN
!----

return
end
