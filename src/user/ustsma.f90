!-------------------------------------------------------------------------------

!VERS

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

subroutine ustsma &
!================

 ( nvar   , nscal  , ncepdp ,                                     &
   ncesmp , iappel ,                                              &
   icepdc , icetsm , itypsm , izctsm ,                            &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Mass source term

! The subroutine ustsma is called at three different stages in the code
!  (iappel = 1, 2 or 3)

! iappel = 1
!    Calculation of the number of cells where a mass source term is
!    imposed: ncesmp
!    Called once at the beginning of the calculation

! iappel = 2
!    Identification of the cells where a mass source term is imposed:
!    array icesmp(ncesmp)
!    Called once at the beginning of the calculation

! iappel = 3
!    Calculation of the values of the mass source term
!    Called at each time step



! The equation for mass conservation becomes

!           d(rho)/dt + div(rho u) = gamma

! The equation for a variable f becomes

!           d(f)/dt = ..... + gamma*(f_i - f)

!   discretized as

!           rho*(f^(n+1) - f^(n))/dt = .....
!                                    + gamma*(f_i - f^(n+1))

! f_i is the value of f associated to the injected mass.
! Two options are available:
!   - the mass flux is injected with the local value of variable f
!                           --> f_i = f^(n+1)
!                   (the equation for f is therefore not modified)
!
!   - the mass flux is injected with a specific value for f
!                           --> f_i is specified by the user


! Variables to be specified by the user
! =====================================

!  ncesmp: number of cells where a mass source term is imposed

!  icetsm(ieltsm): identification of the cells where a mass source
!                  term is imposed.
!                  For each cell where a mass source term is imposed
!                  (ielstm in [1;ncesmp]), icetsm(ieltsm) is the
!                  global index number of the corresponding cell
!                  (icestm(ieltsm) in [1;ncel])

!  smacel(ieltsm,ipr): value of the injection mass rate gamma (kg/m3/s)
!                             in the ieltsm cell with mass source term

!  itypsm(ieltsm,ivar): type of treatment for variable ivar in the
!                       ieltsm cell with mass source term.
!                     * itypsm = 0 --> injection of ivar at local value
!                     * itypsm = 1 --> injection of ivar at user
!                                      specified value

!  smacel(ieltsm,ivar): specified value for variable ivar associated
!                       to the injected mass in the ieltsm cell with
!                       a mass source term
!                                  except for ivar=ipr

!
! Remarks
! =======
!
! - if itypsm(ieltsm,ivar)=0, smacel(ieltsm,ivar) is not used

! - if smacel(ieltsm,ipr)<0, mass is removed from the system,
!     therefore Code_Saturna automatically considers f_i=f^(n+1),
!     whatever the values of itypsm or smacel specified by the user

! - if a value ivar is not linked for a mass source
!     term is imposed, no source term will be taen into account.

! - if a scalar doesn't evolve following the standard equation
!     d(rho f)/dt + d(rho U f)/dx = ...
!     (alternate convective field for instance), the source term
!     set by this routine will nto be correct (except in case of
!     injection at the local value of the variable). The proper source
!     term should be added directly in ustssc.


! Identification of cells
! =======================
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! iappel           ! i  ! <-- ! indicates which at which stage the routine is  !
!                  !    !     !  is called                                     !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
!                  !    !     !  (usable only for iappel > 1)                  !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see uttsma.f90)                              !
! izctsm(ncelet)   ! ia ! <-- ! cells zone for mass source terms definition    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate                     !
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
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iappel

integer          icepdc(*)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izctsm(ncel)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6)
double precision smacel(ncesmp,nvar)

! Local variables

integer          ieltsm
integer          ifac, ii
integer          ilelt, nlelt
integer          izone

double precision vent, vent2
double precision dh, ustar2
double precision xkent, xeent
double precision flucel
double precision vtot  , gamma

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

if (iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! 1. One or two calls

!   First call:
!
!       iappel = 1: ncesmp: calculation of the number of cells with
!                             mass source term


!   Second call (if ncesmp>0):
!       iappel = 2: icetsm: index number of cells with mass source terms

! WARNINGS
! ========
!   Do not use smacel in this section (it is set on the third call, iappel=3)

!   Do not use icetsm in this section on the first call (iappel=1)

!   This section (iappel=1 or 2) is only accessed at the beginning of a
!     calculation. Should the localization of the mass source terms evolve
!     in time, the user must identify at the beginning all cells that can
!     potentially becomea mass source term.

!===============================================================================


!  1.1 To be completed by the user: cell selection
!  -----------------------------------------------

! Example 1: No mass source term (default)
  ieltsm = 0


! Example 2 : Mass source term one in the cells that
!              have a boundary face of color 3 and the cells
!              with a coordinate X between 2.5 and 5.
!
!     In this test in two parts, one mut pay attention not to count
!      the cells twice (a cell with a boundary face of color 3 can
!      also have a coordinate X between 2.5 and 5).
!     One should also pay attention that, on the first call, the
!      array icetsm doesn't exist yet. It mustn't be used outside
!      of tests (iappel.eq.2).

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

  if (.false.) then

    izone = 0
    ieltsm = 0

    !     Cells with coordinate X between 2.5 and 5.

    call getcel('X > 2.5 and X < 5.0',nlelt,lstelt)

    izone = izone + 1

    do ilelt = 1, nlelt
      ii = lstelt(ilelt)
      izctsm(ii) = izone
      ieltsm = ieltsm + 1
      if (iappel.eq.2) icetsm(ieltsm) = ii
    enddo


    !     Cells with a boundary face of color 3

    call getfbr('3',nlelt,lstelt)

    izone = izone + 1

    do ilelt = 1, nlelt
      ifac = lstelt(ilelt)
      ii   = ifabor(ifac)
      !       The cells that have already been counted above are not
      !        counted again.
      if (.not.(xyzcen(1,ii).lt.500.d0.and.                     &
           xyzcen(1,ii).gt.250.d0)    )then
        ieltsm = ieltsm + 1
        izctsm(ii) = izone
        if (iappel.eq.2) icetsm(ieltsm) = ii
      endif
    enddo

 endif


!  1.2 Generic subsection: do not modify
!  -------------------------------------

! --- For iappel = 1,
!      Specification of ncesmp. This block is valid for both examples.

  if (iappel.eq.1) then
    ncesmp = ieltsm
  endif

!-------------------------------------------------------------------------------

elseif (iappel.eq.3) then

!===============================================================================

! 2. For ncesmp > 0 , third call

!       iappel = 3 : itypsm : type of mass source term
!                    smacel : mass source term


! Remark
! ======
! If itypsm(ieltsm,ivar) is set to 1, smacel(ieltsm,ivar) must be set.

!===============================================================================



!  2.1 To be completed by the user: itypsm and smacel
!  --------------------------------------------------

! Example 1: simulation of an inlet condition by mass source terms
!            and printing of the total mass rate.

  vent = 0.1d0
  vent2 = vent**2
  dh     = 0.5d0
  !
  ! Calculation of the inlet conditions for k and epsilon with standard
  !   laws in a circular pipe.
  ustar2 = 0.d0
  xkent  = epzero
  xeent  = epzero

  call keendb                                                   &
  !==========
( vent2, dh, ro0, viscl0, cmu, xkappa,        &
  ustar2, xkent, xeent )

  flucel = 0.d0
  do ieltsm = 1, ncesmp
    smacel(ieltsm,ipr) = 30000.d0
    itypsm(ieltsm,iv) = 1
    smacel(ieltsm,iv) = vent
    if (itytur.eq.2) then
      itypsm(ieltsm,ik) = 1
      smacel(ieltsm,ik) = xkent
      itypsm(ieltsm,iep) = 1
      smacel(ieltsm,iep) = xeent
    else if (itytur.eq.3) then
      itypsm(ieltsm,ir11) = 1
      itypsm(ieltsm,ir12) = 1
      itypsm(ieltsm,ir13) = 1
      itypsm(ieltsm,ir22) = 1
      itypsm(ieltsm,ir23) = 1
      itypsm(ieltsm,ir33) = 1
      smacel(ieltsm,ir11) = 2.d0/3.d0*xkent
      smacel(ieltsm,ir12) = 0.d0
      smacel(ieltsm,ir13) = 0.d0
      smacel(ieltsm,ir22) = 2.d0/3.d0*xkent
      smacel(ieltsm,ir23) = 0.d0
      smacel(ieltsm,ir33) = 2.d0/3.d0*xkent
      itypsm(ieltsm,iep) = 1
      smacel(ieltsm,iep) = xeent
    else if (iturb.eq.50) then
      itypsm(ieltsm,ik) = 1
      smacel(ieltsm,ik) = xkent
      itypsm(ieltsm,iep) = 1
      smacel(ieltsm,iep) = xeent
      itypsm(ieltsm,iphi) = 1
      smacel(ieltsm,iphi) = 2.d0/3.d0
      ! There is no mass source term in the equation for f_bar
    else if (iturb.eq.60) then
      itypsm(ieltsm,ik) = 1
      smacel(ieltsm,ik) = xkent
      itypsm(ieltsm,iomg)= 1
      smacel(ieltsm,iomg)= xeent/cmu/xkent
    endif
    if (nscal.gt.0) then
      do ii = 1, nscal
        itypsm(ieltsm,isca(ii)) = 1
        smacel(ieltsm,isca(ii)) = 1.d0
      enddo
    endif
    flucel = flucel+                                            &
         volume(icetsm(ieltsm))*smacel(ieltsm,ipr)
  enddo

  if (irangp.ge.0) then
    call parsom (flucel)
  endif

  if (iwarni(ipr).ge.1) then
    write(nfecra,1000) flucel
  endif

!-------------------------------------------------------------------------------

! Example 2 : simulation of a suction (by a pump for instance) with a
!             total rate of 80 000 kg/s.
!             The suction rate is supposed to be uniformly distributed
!             on all the cells selected above.

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

  if (.false.) then

    ! Calculation of the total volume of the area where the mass source
    !   term is imposed (the case of parallel computing is taken into
    !   account with the call to parsom).
    vtot = 0.d0
    do ieltsm = 1, ncesmp
      vtot = vtot + volume(icetsm(ieltsm))
    enddo
    if (irangp.ge.0) then
      call parsom (vtot)
    endif

    ! The mass suction rate is gamma = -80000/vtot (in kg/m3/s)
    ! It is set below, with a test for cases where vtot=0. The total
    ! mass rate is calculated for verification.

    if (vtot.gt.0.d0) then
      gamma = -80000.d0/vtot
    else
      write(nfecra,9000) vtot
      call csexit (1)
    endif

    flucel = 0.d0
    do ieltsm = 1, ncesmp
      smacel(ieltsm,ipr) = gamma
      flucel = flucel+                                          &
           volume(icetsm(ieltsm))*smacel(ieltsm,ipr)
    enddo

    if (irangp.ge.0) then
      call parsom (flucel)
    endif

    if (iwarni(ipr).ge.1) then
      write(nfecra,2000) flucel, vtot
    endif

  endif

!-------------------------------------------------------------------------------

endif

!--------
! Formats
!--------

 1000 format(/,'Mass rate generated in the domain: ',E14.5,/)

 2000 format(/,'Mass flux rate generated in the domain: ',E14.5,/,         &
               '                         distributed on the volume: ',E14.5)

 9000 format(/,'Error in ustsma                ',/,                        &
               '   the volume of the mass suction area is = ',E14.5,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsma
