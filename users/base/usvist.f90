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

subroutine usvist &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

!===============================================================================
! Purpose:
! -------

! User subroutine.

! Modify turbulent viscosity

! This subroutine is called at beginning of each time step
! after the computation of the turbulent viscosity
! (physical quantities have already been computed in usphyv)

! Turbulent viscosity VISCT (kg/(m s)) can be modified

! A modification of the turbulent viscosity can lead to very
! significant differences betwwen solutions and even give wrong
! results

! This subroutine is therefore reserved to expert users

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss
! ncesmp           ! i  ! <-- ! number of cells with mass source term
! iphas            ! i  ! <-- ! phase number
! icepdc(ncelet    ! te ! <-- ! head loss cell numbering                       !
! icetsm(ncesmp    ! te ! <-- ! numbering of cells with mass source term       !
! itypsm           ! te ! <-- ! kind of mass source for each variable          !
! (ncesmp,nvar)    !    !     !  (cf. ustsma)                                  !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! ra ! <-- ! work array for head loss terms                 !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! ra ! <-- ! values of variables related to mass source     !
! (ncesmp,*   )    !    !     ! term. If ivar=ipr, smacel=mass flux            !
! w1...8(ncelet    ! ra ! --- ! work arrays                                    !
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
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel, iccocg, inc
integer          iuiph, iviph, iwiph
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvst, iphydp
double precision dudx, dudy, dudz, sqdu, visct, rom

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  Example :
!       For phase 1:
!                visct = max(visct, rom * sqrt(dudx**2 + dudy**2 + dudz**2)
!                (intentionally fancyful relation)
!                Remark: incomming viscosity is consistent with the selected
!                turbulence modelling
!       For other phases:
!                We keep the viscosity computed by the selected turbulence
!                modelling

!===============================================================================

!===============================================================================
! 1.1 Selection of the phase to deal with
!===============================================================================

if(iphas.ne.1) then
  return
endif

!===============================================================================
! 1.2 Initialization
!===============================================================================

! --- Memory
idebia = idbia0
idebra = idbra0

! --- Number associated to variables (in RTP)
iuiph = iu(iphas)
iviph = iv(iphas)
iwiph = iw(iphas)

! --- Physical quantity numbers in PROPCE (physical quantities defined
!     at each cell center)
ipcvst = ipproc(ivisct(iphas))
ipcrom = ipproc(irom  (iphas))

! --- Boundary condition number associated to variables in COEFA and COEFB
!      JB=>?  (c.l. std, i.e. non flux)
ipcliu = iclrtp(iuiph,icoef)
ipcliv = iclrtp(iviph,icoef)
ipcliw = iclrtp(iwiph,icoef)

!===============================================================================
! 1.3 Compute velocity gradient
!===============================================================================

iccocg = 1
inc = 1
iphydp = 0

! W1 = DUDX, W2 = DUDY, W3=DUDZ

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iuiph  , imrgra , inc    , iccocg , iphydp ,                   &
   nswrgr(iuiph) , imligr(iuiph) ,                                &
   iwarni(iuiph) , nfecra ,                                       &
   epsrgr(iuiph) , climgr(iuiph) , extrag(iuiph) ,                &
   ia     ,                                                       &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iuiph) , coefa(1,ipcliu) , coefb(1,ipcliu) ,            &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   ra     )

!===============================================================================
! 1.4 Computation of the dynamic viscosity
!===============================================================================

do iel = 1, ncel

! --- Current dynamic viscosity and fluid density
  visct = propce(iel,ipcvst)
  rom   = propce(iel,ipcrom)

! --- Various computations
  dudx = w1(iel)
  dudy = w2(iel)
  dudz = w3(iel)
  sqdu = sqrt(dudx**2+dudy**2+dudz**2)

! --- Computation of the new dynamic viscosity
  visct = max (visct,rom*sqdu)

! --- Store the new computed dynamic viscosity
  propce(iel,ipcvst) = visct

enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
