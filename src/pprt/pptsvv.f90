!-------------------------------------------------------------------------------

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
! Function:
! ---------

!> \file pptssc.f90
!>
!> \brief This subroutine defines the source terms for vectors which are part of
!> specific physics models. Source terms are defined over one time step.
!>
!> Warning: rovsdt and smbrs already hold possible user
!> source terms values and thus have to be incremented (and not overwritten).
!>
!> For stability reasons, only positive terms are added to rovsdt, while there
!> are no such constrains on values to be added to smbrs.
!>
!> In the case of a source term of the form cexp, the source term
!> should be implemented as follows:
!> \f[
!>   smbrv(i)  = smbrv(i)  + cexp(i) + \sum_j cimp(i,j)*var(j)
!> \f]
!>
!> rovsdt and smbrs are provided here respectively in kg/s and in kg/s*[scalar].
!> Examples:
!>   velocity \f$ kg m/s^2 \f$
!>   temperature \f$ kg K/s \f$
!>   enthalpy \f$ J/s \f$
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         number in list of additional variables
!> \param[in,out] smbrv         explicit source term part
!______________________________________________________________________________!


subroutine pptsvv(iscal, smbrv)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrv(3,ncelet)

! Local variables

!===============================================================================

! MHD module
! Joule effect
! Electric arcs
! Ionic conduction

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call eltsvv(ivarfl(isca(iscal)), smbrv)
endif

!----
! End
!----

return

end subroutine
