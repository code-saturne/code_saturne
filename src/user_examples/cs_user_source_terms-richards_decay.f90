!-------------------------------------------------------------------------------

!                      Code_Saturne version 4.0-alpha
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine ustssc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use darcy_module

!===============================================================================

implicit none

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

character*80     chaine
integer          iiscvr,  iel
integer          ilelt, nlelt

double precision lambda

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom

integer delay_id
double precision, dimension(:), pointer :: delay, saturation
character*80     fname

!===============================================================================

! In groundwater module flow, the radioactive decay of solute is treated as a
! source term in the transport equation

allocate(lstelt(ncel))

call field_get_label(ivarfl(isca(iscal)), chaine)

iiscvr = iscavr(iscal)

if (iwarni(isca(iscal)).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!< [richards_decay]
! Set the first order decay coefficient
lambda = 1.d-2

! Set radioactive decay for the first solute
if (isca(iscal).eq.1) then
  do iel = 1, ncel
    crvimp(iel) = -volume(iel)*lambda
  enddo
endif
!< [richards_decay]

 1000 format(' User source terms for variable ',A8,/)

deallocate(lstelt)

return
end subroutine ustssc
