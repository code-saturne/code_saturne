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

subroutine cprays &
!================

  ( ivar   , ncelet , ncel   ,                                    &
    volume , propce , smbrs  , rovsdt )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
!   PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
!   IMPLICITE ET EXPLICITE DANS L'EQUATION DES PARTICULES
!   DE LA CLASSE ICLA

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ivar             ! e  ! <-- ! numero de la variable scalaire                 !
!                  !    !     !   energie (enthalpie h2) pour le               !
!                  !    !     !   charbon                                      !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! volume(ncelet    ! tr ! <-- ! volume des cellules                            !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! smbrs(ncelet     ! tr ! <-- ! second membre du systeme                       !
! rovsdt(ncelet    ! tr ! <-- ! diagonale du systeme                           !
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
include "entsor.h"
include "numvar.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          ivar , ncelet, ncel

double precision volume(ncelet)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)
double precision propce(ncelet,*)

! VARIABLES LOCALES

integer          iel , numcla , ipcl

!===============================================================================

!===============================================================================
! 1. RECHERCHE DE LA ZONE MEMOIRE (IPH) EN FONCTION DU NUMERO DE PHASE
!    COURANT IPHAS POUR TROUVER LES BONS TERMES SOURCES
!===============================================================================

numcla = ivar-isca(ih2(1))+1
ipcl   = 1+numcla

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
!===============================================================================


do iel = 1,ncel
  propce(iel,ipproc(itsri(ipcl))) = max(-propce(iel,ipproc(itsri(ipcl))),zero)
enddo

do iel = 1,ncel
  if ( propce(iel,ipproc(ix2(numcla))) .gt. epzero ) then

!--> PARTIE EXPLICITE

    smbrs(iel)  = smbrs(iel) +  propce(iel,ipproc(itsre(ipcl)))*volume(iel)        &
                               *propce(iel,ipproc(ix2(numcla)))

!--> PARTIE IMPLICITE

    rovsdt(iel) = rovsdt(iel) + propce(iel,ipproc(itsri(ipcl)))*volume(iel)
  endif

enddo


!----
! FIN
!----

return

end subroutine
