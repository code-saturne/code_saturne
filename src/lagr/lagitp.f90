!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine lagitp &
!================

 ( propce , tempct )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LA TEMPERATURE DES PARTICULES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpart,2)      !    !     !                                                !
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
use lagpar
use lagran
use ppppar
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)
double precision tempct(nbpart,2)

! Local variables

integer          npt , iel
double precision srad
double precision, dimension(:), allocatable :: tcarac, pip

!===============================================================================

allocate(tcarac(nbpart), pip(nbpart))

!===============================================================================
!     REMPLISSAGE DU TEMPS CARACTERISTIQUE ET DU "PSEUDO SECOND MEMBRE"
!===============================================================================

do npt = 1,nbpart

  if (ipepa(jisor,npt).gt.0) then

    tcarac(npt) = tempct(npt,1)

    if (nor.eq.1) then
      pip(npt) = eptpa(jtf,npt)
    else
      pip(npt) = eptp(jtf,npt)
    endif

  endif

enddo

!===============================================================================
!     PRISE EN COMPTE DU RAYONNEMENT S'IL Y A LIEU
!===============================================================================

if (iirayo.gt.0) then

  do npt = 1,nbpart

    iel = ipepa(jisor,npt)

    if (iel.gt.0) then

      if (nor.eq.1) then

        srad = pi *eptpa(jdp,npt) *eptpa(jdp,npt)                 &
                  *pepa(jreps,npt) *(propce(iel,ipproc(ilumin))   &
                        -4.d0 *stephn *eptpa(jtp,npt)**4 )
        pip(npt) = eptpa(jtf,npt)                                 &
               +tcarac(npt) *srad /eptpa(jcp,npt) /eptpa(jmp,npt)
      else

        srad = pi *eptp(jdp,npt) *eptp(jdp,npt) *pepa(jreps,npt)  &
                *(propce(iel,ipproc(ilumin))                      &
                -4.d0 *stephn *eptp(jtp,npt)**4 )
        pip(npt) = eptp(jtf,npt)                                  &
                +tcarac(npt) *srad /eptp(jcp,npt) /eptp(jmp,npt)

      endif

    endif

  enddo

endif

!===============================================================================
!     INTEGRATION
!===============================================================================

call lagitg(jtp, tcarac, pip)

!===============================================================================

deallocate(tcarac, pip)

!----
! FIN
!----

end subroutine
