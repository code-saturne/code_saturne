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

subroutine lecamo &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,          &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , propce , propfa , propfb ,                   &
   coefa  , coefb  , frcxt  ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================

! FONCTION :
! ----------
! LECTURE DES FICHIERS SUITE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension du calcul                            !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nnod             ! e  ! <-- ! nombre de noeuds                               !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! dt(ncelet)       ! tr ! --> ! pas de temps                                   !
! rtp              ! tr ! --> ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant        )          !
! propce           ! tr ! --> ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! --> ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! --> ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! --> ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! frcxt(ncelet,    ! tr ! --> ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "optcal.h"
include "pointe.h"
include "numvar.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor, nnod
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision frcxt(ncelet,3,nphas)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra


!===============================================================================



!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

write(nfecra,1000)

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. LECTURE DU FICHIER SUITE PRINCIPAL
!===============================================================================

call lecamp (idebia , idebra ,                                    &
!     ==========
             ncelet , ncel   ,                                    &
             nvar   , nscal  , nphas  ,                           &
             nideve , nrdeve , nituse , nrtuse ,                  &
             idevel , ituser , ia     ,                           &
             rtp    ,                                             &
             rdevel , rtuser , ra     )


!===============================================================================
! 3. LECTURE DU FICHIER SUITE AUXILIAIRE
!===============================================================================

if (ileaux.eq.1) then

  call lecamx (idebia , idebra ,                                  &
!       ==========
               ndim   , ncelet , ncel   , nfac   , nfabor ,       &
               nnod   , nvar   , nscal  , nphas  ,                &
               nideve , nrdeve , nituse , nrtuse ,                &
               idevel , ituser , ia     ,                         &
               dt     , rtp    , propce , propfa , propfb ,       &
               coefa  , coefb  , frcxt  ,                         &
               rdevel , rtuser , ra     )

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
