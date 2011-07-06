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

subroutine csc2ts &
!================

 ( nvar   , nscal  ,                                              &
   ncecpl,                                                        &
   ivar   ,                                                       &
   lcecpl ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   crvexp , crvimp ,                                              &
   rvcpce )

!===============================================================================
! FONCTION :
! ----------

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! i  ! <-- ! variable number                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant            prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour part implicit          !
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
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ivar
integer          ncecpl

integer          lcecpl(ncecpl)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision crvexp(ncelet), crvimp(ncelet)
double precision rvcpce(ncecpl)

! Local variables

integer          iel    , ipcrom
integer          ipt    , ielloc
double precision xdis   , xloc   , xtau

!----------------------------------------------------------------------------------


ipcrom = ipproc(irom)

xtau = 100.d0*dtref

do ipt = 1, ncecpl

  ielloc = lcecpl(ipt)

  xdis = rvcpce(ipt)
  xloc = rtpa(ielloc,ivar)

  crvexp(ielloc) = crvexp(ielloc) &
                 + volume(ielloc)*propce(ielloc,ipcrom)*(xdis-xloc)/xtau

enddo

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
