!-------------------------------------------------------------------------------

!VERS


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

subroutine usray4 &
!================

 ( nvar   , nscal  ,                                              &
   mode   ,                                                       &
   itypfb ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   tparop , hparop , tempk  ,                                     &
   ra     )

!===============================================================================
! Purpose:
! --------

! User subroutine for input of radiative transfer parameters:

!   Temperature <--> enthalpy convertion
!   Usefull if the thermal scalar is an enthalpy.

!   PRECAUTIONS: ENTHALPY MUST BE CONVERTED IN KELVIN TEMPERATURE

!   Warning: it is forbidden to modify MODE in this subroutine

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! mode             ! i  ! <-- ! convertion mode                                !
!                  !    !     ! mode = 1 enthaly -> temperature                !
!                  !    !     ! mode =-1 temperature -> enthaly                !
! itypfb           ! ia ! <-- ! boundary face types                            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! tparop(nfabor)   ! i  ! <-- ! temperature in kelvin for wall boundary faces  !
! hparop(nfabor)   ! i  ! --> ! enthalpy for wall boundary faces               !
! tempk(ncelet)    ! i  ! --> ! temperature in kelvin                          !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use cpincl
use ppincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mode

integer          itypfb(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision tempk(ncelet)
double precision tparop(nfabor), hparop(nfabor)

double precision ra(*)


! Local variables

integer          iel , ifac , iscal

!===============================================================================

!===============================================================================
! 1 - INITIALISATIONS GENERALES
!===============================================================================


iscal = iscalt

!===============================================================================
!  2.1 - Tempearature in kelvin for cells
!===============================================================================

!---> enthalpy -> temperature convertion (MODE =  1)
!     -----------------------------------------------


if (mode.eq.1) then

  do iel = 1,ncel
    call usthht(mode,rtpa(iel,isca(iscal)),tempk(iel))
  enddo

endif


!===============================================================================
!  2.2 - Enthalpy for wall boundary faces
!===============================================================================

!---> Temperature -> enthalpy (MODE = -1)
!     -----------------------------------


if (mode.eq.-1) then

  do ifac = 1,nfabor

    if (itypfb(ifac).eq.iparoi .or.                               &
        itypfb(ifac).eq.iparug )then

      call usthht(mode,hparop(ifac),tparop(ifac))

    endif

  enddo

endif

!----
! END
!----

return
end subroutine
