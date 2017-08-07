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

!> \file cs_user_source_terms-richards.f90
!>
!> \brief User source terms example.
!>
!> See \subpage cs_user_source_terms and \subpage cs_user_source_terms-scalar_in_a_channel
!> for examples.
!>
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \section ustsnv_use  Usage
!>
!> The additional source term is decomposed into an explicit part (\c crvexp) and
!> an implicit part (\c crvimp) that must be provided here.
!> The resulting equation solved by the code for a velocity is:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
!>   = \tens{crvimp} \cdot \vect{u} + \vect{crvexp}
!> \f]
!>
!> Note that \c crvexp and \c crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m/s2
!>   - crvimp is expressed in kg/s
!>
!> The \c crvexp and \c crvimp arrays are already initialized to 0
!> before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> \remark The additional force on \f$ x_i \f$ direction is given by
!>  \c crvexp(i, iel) + vel(j, iel)* crvimp(j, i).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a
!> \ref getcel command. For more info on the syntax of the \ref getcel command,
!> refer to the user manual or to the comments on the similar command
!> \ref getfbr in the routine \ref cs_user_boundary_conditions.
!
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
use cs_c_bindings

!===============================================================================

implicit none

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

character(len=80) ::  chaine

integer          iiscvr,  iel
integer          icelt, ncelt
integer, allocatable, dimension(:) :: lstcel

double precision lambda, leaching_time, leaching_volume
double precision, dimension(:), pointer :: delay, saturation

type(var_cal_opt) :: vcopt

!===============================================================================

allocate(lstcel(ncel))

call field_get_label(ivarfl(isca(iscal)), chaine)

iiscvr = iscavr(iscal)

call field_get_key_struct_var_cal_opt(ivarfl(isca(iscal)), vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!======================
! Leaching
!======================

!< [richards_leaching]
! Leaching from volumic zone labelled as 'LEACHING_ZONE'
! Set the duration of leaching
leaching_time = 1.d2

call getcel('LEACHING_ZONE', ncelt, lstcel)

! Compute the volume of the leaching zone
leaching_volume = 0.d0
do icelt = 1, ncelt
  iel = lstcel(icelt)
  leaching_volume = leaching_volume + volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(leaching_volume)
endif

! Compute the source term
do icelt = 1, ncelt
  iel = lstcel(icelt)
  !Progressive leaching
  if (ttcabs.lt.leaching_time) then
    crvexp(iel) = 1. / leaching_time * volume(iel) / leaching_volume &
         * exp(-lambda*ttcabs) / saturation(iel) / delay(iel)
  else
    crvexp(iel) = 0.d0
  endif
enddo
!< [richards_leaching]


 1000 format(' User source terms for variable ',A8,/)

deallocate(lstcel)

return
end subroutine ustssc
