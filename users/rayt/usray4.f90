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

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   mode   ,                                                       &
   itypfb ,                                                       &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   tparop , hparop , tempk  ,                                     &
   rdevel , rtuser ,                                              &
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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iphas            ! i  ! <-- ! current phase number                           !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! mode             ! i  ! <-- ! convertion mode                                !
!                  !    !     ! mode = 1 enthaly -> temperature                !
!                  !    !     ! mode =-1 temperature -> enthaly                !
! itypfb           ! ia ! <-- ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary developpement !
! ituser(nituse    ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !                                                !
! tparop(nfabor)   ! i  ! <-- ! temperature in kelvin for wall boundary faces  !
! hparop(nfabor)   ! i  ! --> ! enthalpy for wall boundary faces               !
! tempk(ncelet)    ! i  ! --> ! temperature in kelvin                          !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary developpement    !
! rtuser(nituse)   ! ra ! <-- ! user-reserved real work array                  !
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

integer          idbia0 , idbra0
integer          nvar   , nscal  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          mode

integer          itypfb(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)

double precision tempk(ncelet)
double precision tparop(nfabor), hparop(nfabor)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! Local variables

integer          idebia , idebra
integer          iel , ifac , iscal

!===============================================================================

!===============================================================================
! 1 - INITIALISATIONS GENERALES
!===============================================================================

idebia = idbia0
idebra = idbra0

iscal = iscalt(iphas)

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
