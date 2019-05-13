!-------------------------------------------------------------------------------

!                      Code_Saturne version 4.0-alpha
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_mass_source_terms.f90
!>
!> \brief Mass source term example.
!>
!> See \subpage cs_user_mass_source_terms for examples.
!>

!-------------------------------------------------------------------------------
!>           Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iappel        indicates which at which stage the routine is
!>                              is called
!> \param[in]     icepdc        index number of cells with head loss terms
!>                              (usable only for iappel > 1)
!> \param[in,out] icetsm        index number of cells with mass source terms
!> \param[in,out] itypsm        type of mass source term for each variable
!>                               (see uttsma.f90)
!> \param[in]     izctsm        cells zone for mass source terms definition
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in,out] smacel        value associated to each variable in the mass
!>                              source terms or mass rate
!______________________________________________________________________________!

subroutine cs_user_mass_source_terms &
 ( nvar   , nscal  , ncepdp ,                                     &
   ncesmp , iappel ,                                              &
   icepdc , icetsm , itypsm , izctsm ,                            &
   dt     ,                                                       &
   ckupdc , smacel )

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iappel

integer          icepdc(*)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izctsm(ncel)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp)
double precision smacel(ncesmp,nvar)

!< [loc_var]
! Local variables

integer          ieltsm
integer          ifac, ii
integer          ilelt, nlelt
integer          izone

double precision wind, wind2
double precision dh, ustar2
double precision xkent, xeent
double precision flucel
double precision vtot  , gamma

integer, allocatable, dimension(:) :: lstelt

type(var_cal_opt) :: vcopt

!< [loc_var]

!===============================================================================

!< [allocate]
! Allocate a temporary array for cells selection
allocate(lstelt(ncel))
!< [allocate]

!< [one_or_two]
if (iappel.eq.1.or.iappel.eq.2) then
!< [one_or_two]

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
!     potentially become a mass source term.

!===============================================================================


!  1.1 To be completed by the user: cell selection
!  -----------------------------------------------

! Example 1: No mass source term (default)

!< [example_1_1]
  ieltsm = 0
!< [example_1_1]

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

!< [example_1_2]
  izone = 0
  ieltsm = 0

  ! Cells with coordinate X between 2.5 and 5.

  call getcel('X > 2.5 and X < 5.0',nlelt,lstelt)

  izone = izone + 1

  do ilelt = 1, nlelt
    ii = lstelt(ilelt)
    izctsm(ii) = izone
    ieltsm = ieltsm + 1
    if (iappel.eq.2) icetsm(ieltsm) = ii
  enddo


  ! Cells with a boundary face of color 3

  call getfbr('3',nlelt,lstelt)

  izone = izone + 1

  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    ii   = ifabor(ifac)
    ! The cells that have already been counted above are not
    ! counted again.
    if (.not.(xyzcen(1,ii).lt.500.d0.and.                     &
         xyzcen(1,ii).gt.250.d0)    )then
      ieltsm = ieltsm + 1
      izctsm(ii) = izone
      if (iappel.eq.2) icetsm(ieltsm) = ii
    endif
  enddo
!< [example_1_2]

!  1.2 Generic subsection: do not modify
!  -------------------------------------

! --- For iappel = 1,
!      Specification of ncesmp. This block is valid for both examples.
!< [generic_sub]
  if (iappel.eq.1) then
    ncesmp = ieltsm
  endif
!< [generic_sub]
!-------------------------------------------------------------------------------

!< [call_3]
elseif (iappel.eq.3) then
!< [call_3]

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

!< [example_2_1]
  wind = 0.1d0
  wind2 = wind**2
  dh     = 0.5d0
!< [example_2_1]


  ! Calculation of the inlet conditions for k and epsilon with standard
  !   laws in a circular pipe.

!< [inlet_cal]
  ustar2 = 0.d0
  xkent  = epzero
  xeent  = epzero

  call turbulence_bc_ke_hyd_diam(wind2, dh, ro0, viscl0,  &
                                 ustar2, xkent, xeent )

  flucel = 0.d0
  do ieltsm = 1, ncesmp
    smacel(ieltsm,ipr) = 30000.d0
    itypsm(ieltsm,iv) = 1
    smacel(ieltsm,iv) = wind
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

  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

  if (vcopt%iwarni.ge.1) then
    write(nfecra,1000) flucel
  endif
!< [inlet_cal]

!-------------------------------------------------------------------------------

! Example 2 : simulation of a suction (by a pump for instance) with a
!             total rate of 80 000 kg/s.
!             The suction rate is supposed to be uniformly distributed
!             on all the cells selected above.

  ! Calculation of the total volume of the area where the mass source
  !   term is imposed (the case of parallel computing is taken into
  !   account with the call to parsom).

!< [calcul_total]
  vtot = 0.d0
  do ieltsm = 1, ncesmp
    vtot = vtot + volume(icetsm(ieltsm))
  enddo
  if (irangp.ge.0) then
    call parsom (vtot)
  endif
!< [calcul_total]

  ! The mass suction rate is gamma = -80000/vtot (in kg/m3/s)
  ! It is set below, with a test for cases where vtot=0. The total
  ! mass rate is calculated for verification.

!< [mass_suction]
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

  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

  if (vcopt%iwarni.ge.1) then
    write(nfecra,2000) flucel, vtot
  endif
!< [mass_suction]

!-------------------------------------------------------------------------------
!< [end_call_3]
endif
!< [end_call_3]
!--------
! Formats
!--------
!< [format]
 1000 format(/,'Mass rate generated in the domain: ',E14.5,/)

 2000 format(/,'Mass flux rate generated in the domain: ',E14.5,/,         &
               '                         distributed on the volume: ',E14.5)

 9000 format(/,'Error in cs_user_mass_source_terms',/,                     &
               '   the volume of the mass suction area is = ',E14.5,/)
!< [format]
!----
! End
!----
!< [deallocate]
! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_mass_source_terms
!< [deallocate]
