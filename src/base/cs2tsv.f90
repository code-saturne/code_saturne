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

subroutine cs2tsv &
!================

 ( ncecpl,  lcecpl ,                                              &
   vela   , propce , crvexp , rvcpce )

!===============================================================================
! FONCTION :
! ----------

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! crvexp(3,ncelet) ! tr ! --> ! tableau de travail pour part explicit          !
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

integer          ncecpl

integer          lcecpl(ncecpl)

double precision propce(ncelet,*)
double precision crvexp(3,ncelet)
double precision rvcpce(3,ncecpl)
double precision vela(3,ncelet)

! Local variables

integer          ipcrom , isou
integer          ipt    , ielloc
double precision xdis   , xloc   , xtau   , rovtau

!----------------------------------------------------------------------------------


ipcrom = ipproc(irom)

xtau = 100.d0*dtref

do ipt = 1, ncecpl

  ielloc = lcecpl(ipt)

  rovtau = volume(ielloc)*propce(ielloc,ipcrom)/xtau

  do isou = 1, 3
    xdis = rvcpce(isou,ipt)
    xloc = vela(isou,ielloc)
    crvexp(isou,ielloc) = crvexp(isou,ielloc) + rovtau*(xdis-xloc)
  enddo

enddo

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
