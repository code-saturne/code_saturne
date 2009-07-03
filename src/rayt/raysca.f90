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

subroutine raysca &
!================

 ( iisca  ,                                                       &
   ncelet , ncel   ,                                              &
   smbrs  , rovsdt , volume , propce  )

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
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iisca            ! e  ! <-- ! num scalaire temperature ou enthalpie          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! smbrs(ncelet)    ! tr ! <-- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! <-- ! tableau de travail pour terme instat           !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
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
include "cstnum.h"
include "cstphy.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"
include "ppincl.h"
include "radiat.h"
include "entsor.h"
include "numvar.h"

!===============================================================================

! Arguments

integer          iisca , ncelet , ncel

double precision volume(ncelet)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)
double precision propce(ncelet,*)

! VARIABLES LOCALES

integer          iel

!===============================================================================

!===============================================================================
! 1. PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
!===============================================================================


if (abs(iscsth(iisca)).eq.1 .or. iscsth(iisca).eq.2) then

  do iel = 1,ncel
    propce(iel,ipproc(itsri(1))) = max(-propce(iel,ipproc(itsri(1))),zero)
  enddo

  do iel = 1,ncel

!--> PARTIE EXPLICITE

    smbrs(iel) = smbrs(iel)  + propce(iel,ipproc(itsre(1)))*volume(iel)

!--> PARTIE IMPLICITE

    rovsdt(iel)= rovsdt(iel) + propce(iel,ipproc(itsri(1)))*volume(iel)

  enddo

endif
!

!----
! FIN
!----

return

end
