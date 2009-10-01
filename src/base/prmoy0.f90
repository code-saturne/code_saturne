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

subroutine prmoy0 &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   , nfac   , nfabor ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , idevel , ituser , ia     ,                            &
   volume , pvar   ,                                              &
   rdevel , rtuser , ra    )

!===============================================================================
! FONCTION :
! ----------

! RECALAGE DE LA VARIABLE PVAR (PRESSION)POUR OBTENIR
! QU'ELLE SOIT A MOYENNE NULLE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de la phase traitee                     !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! volume(ncelet    ! tr ! <-- ! volume des elements                            !
! pvar             ! tr ! <-- ! tableau de valeurs au cellules                 !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "parall.h"
include "cstphy.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ncelet , ncel   , nfac , nfabor
integer          nideve , nrdeve , nituse , nrtuse

integer          iphas, idevel(nideve), ituser(nituse)
integer          ia(*)
double precision volume(ncelet), pvar(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          iel
double precision pmoy, pr0iph

!===============================================================================

pmoy = 0.0d0
do iel = 1, ncel
  pmoy = pmoy + volume(iel) * pvar(iel)
enddo
if (irangp.ge.0) then
  call parsom (pmoy)
  !==========
endif

pmoy = pmoy / voltot
pr0iph = pred0(iphas)
do iel = 1, ncel
  pvar(iel) = pvar(iel) - pmoy + pr0iph
enddo

!----
! FIN
!----

return

end subroutine
