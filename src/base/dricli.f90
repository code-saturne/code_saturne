!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine dricli &
!=====================================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! Purpose:
! -------

!    Fill boundary conditions arrays (icodcl, rcodcl) for scalars with
!    drift velocity in the framework of the Nerisson model



!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor)   ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb(nfabor)   ! ia ! --> ! boundary face types                            !
! izfppp(nfabor)   ! ia ! --> ! boundary face zone number                      !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
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
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use ctincl
use elincl
use pointe, only : uetbor
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal, iscal, iel, ifac, ilelt

integer          icodcl(nfabor,nvar), nlelt
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor), ipcvis

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)

! Local variables

double precision taupplu, vsedim

integer, allocatable, dimension(:) :: lstelt

double precision, allocatable, dimension(:) :: nufluid

!===============================================================================


!===============================================================================
! Initialization
!===============================================================================


!===============================================================================
! Specific treatment for aerosol deposition
!===============================================================================

ipcvis = ipproc(iviscl)

!
   if (.not.(allocated(taupae))) then
      allocate(taupae(nscaus+nscadr,ncel))
   endif
!
allocate(nufluid(ncel))



! Retrieve of the viscosity and of the relaxation time scale
! taupae
!


do iel = 1, ncel

   if (ipcvis.eq.0) then
      nufluid(iel) = viscl0  / propce(iel, ipproc(irom))

      do iscal = nscaus + 1, nscaus + nscadr

         taupae(iscal,iel) = (diapart(iscal)**2)*rhopart(iscal)        &
              / 18.d0 / viscl0
      enddo

   else
      nufluid(iel) = propce(iel,ipcvis) / propce(iel, ipproc(irom))

      do iscal = nscaus + 1, nscaus + nscadr

         taupae(iscal,iel) = (diapart(iscal)**2)*rhopart(iscal)          &
              / 18.d0 / propce(iel,ipcvis)

      enddo

   endif

enddo

do ifac = 1, nfabor

   iel = ifabor(ifac)

   if ((itypfb(ifac).eq.iparoi).or.(itypfb(ifac).eq.iparug)) then


      do iscal = nscaus + 1, nscaus + nscadr

         !!     calculation of tau_p^+
         !!     the dimensionless particle relaxation time

         taupplu = taupae(iscal,iel) * uetbor(ifac) ** 2 / nufluid(iel)

         !!     Calculation of the deposition velocity
         !!
         !!
         vsedim = abs((taupae(iscal,iel)/surfbn(ifac))*(gx*surfbo(1,ifac)     &
              + gy*surfbo(2,ifac)                        &
              + gz*surfbo(3,ifac)))


         icodcl(ifac,isca(iscal)) = 3
         rcodcl(ifac,isca(iscal),3) = rtp(iel, isca(iscal))

      enddo

   endif

enddo

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----
deallocate(nufluid)

return
end subroutine
