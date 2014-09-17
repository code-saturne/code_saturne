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

subroutine lagimp &
!================

 ( tsvar  , auxl1  , auxl2  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LA MASSE DES PARTICULES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! auxl1(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl2(nbpmax)    ! tr ! --- ! tableau de travail                             !
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
use cstphy
use cstnum
use optcal
use entsor
use lagdim, only: nbpmax, nvp1
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax) , auxl2(nbpmax)

! Local variables

integer          npt

!===============================================================================


do npt = 1,nbpart

  if ( ipepa(jisor,npt).gt.0 ) then

    auxl1(npt)     = 1.d0
    auxl2(npt)     = eptp(jmp,npt)

  endif
enddo

call lagitg                                                       &
!==========
 ( jmp    ,                                                       &
   auxl1  , auxl2  , tsvar  )

!===============================================================================

!----
! FIN
!----

end subroutine
