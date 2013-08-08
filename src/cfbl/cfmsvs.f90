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

subroutine cfmsvs &
!================

 ( nvar   , nscal  ,                                              &
   iscal  ,                                                       &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
   w1     , w2     , w3     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA "VISCOSITE" AUX FACES (Delta t c2)
!   POUR LA RESOLUTION DE LA MASSE VOLUMIQUE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iscal            ! i  ! <-- ! scalar number                                  !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! viscf(nfac)      ! tr ! --> ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --> ! visc*surface/dist aux faces de bord            !
! w1..3(ncelet)    ! tr ! --- ! tableau de travail                             !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iscal


double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision viscf(nfac), viscb(nfabor)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)

! Local variables

integer          ifac  , iel

integer          imvis1, iccfth, imodif

double precision rvoid(1)

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================


do ifac = 1, nfac
  viscf(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  viscb(ifac) = 0.d0
enddo

!===============================================================================
! 2. VISCOSITE AUX FACES
!===============================================================================

! --- Calcul de c2 et affectation a W1
iccfth = 126
imodif = 0
call cfther                                                       &
!==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce ,                            &
   w1     , rvoid  , w2     , w3     )

! --- "Vitesse" de diffusion de RHO = dt*c2
do iel = 1, ncel

  w1(iel) = dt(iel)*w1(iel)

enddo

! --- Calcul de (c2)ij par une moyenne harmonique
imvis1 = 1

call viscfa                                                       &
!==========
 ( imvis1 ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!     Au bord, voir le sous-programme appelant.

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine
