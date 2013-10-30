!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine lecamo &
!================

 ( ncelet , ncel   , nfabor , nvar   , nscal  ,                   &
   dt     , rtp    , propce ,                                     &
   frcxt  , prhyd  )

!===============================================================================

! FONCTION :
! ----------
! LECTURE DES FICHIERS SUITE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! --> ! pas de temps                                   !
! rtp              ! tr ! --> ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant        )          !
! propce           ! tr ! --> ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! frcxt(3,ncelet)  ! tr ! --> ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! prhyd(ncelet)    ! ra ! --> ! pression hydrostatic predite                   !
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
use cstphy
use cstnum
use entsor
use optcal
use pointe
use numvar
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfabor
integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision frcxt(3,ncelet), prhyd(ncelet)

! Local variables


!===============================================================================



!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

write(nfecra,1000)


!===============================================================================
! 2. LECTURE DU FICHIER SUITE PRINCIPAL
!===============================================================================

call lecamp &
!==========
( ncelet , ncel   ,                                    &
  nvar   , nscal  ,                                    &
  rtp    )


!===============================================================================
! 3. LECTURE DU FICHIER SUITE AUXILIAIRE
!===============================================================================

if (ileaux.eq.1) then

  call lecamx &
  !==========
( ncelet , ncel   , nfabor , nvar   , nscal  ,       &
  dt     , propce ,                                  &
  frcxt  , prhyd  )

endif


!===============================================================================
! 4. SORTIE
!===============================================================================

write(nfecra,2000)

!===============================================================================
! 5. FORMATS
!===============================================================================


#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
' ----------------------------------------------------------- ',/,&
                                                                /,&
     3X,'** LECTURE DES FICHIERS SUITE PRINCIPAL ET AUXILIAIRE',/,&
     3X,'   ------------------------------------------------- ',/)
 2000 format(/,                                                   &
' ----------------------------------------------------------- ',/)

#else

 1000 format(/,                                                   &
' ----------------------------------------------------------- ',/,&
                                                                /,&
     3X,'** READING MAIN AND AUXILIARY RESTART FILES'          ,/,&
     3X,'   ----------------------------------------'          ,/)
 2000 format(/,                                                   &
' ----------------------------------------------------------- ',/)

#endif


return
end subroutine
