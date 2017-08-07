!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \brief Mass source term user subroutine.
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

! Local variables

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


endif

!--------
! Formats
!--------

 1000 format(/,'Mass rate generated in the domain: ',E14.5,/)

 2000 format(/,'Mass flux rate generated in the domain: ',E14.5,/,         &
               '                         distributed on the volume: ',E14.5)

 9000 format(/,'Error in cs_user_mass_source_terms',/,               &
               '   the volume of the mass suction area is = ',E14.5,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_mass_source_terms
