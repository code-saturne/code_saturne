!-------------------------------------------------------------------------------

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

!> \file cscpce.f90
!> \brief Preparation of sending velocity variables for coupling between
!> two instances of Code_Saturne via boundary faces.
!> Received indformation will be transformed into boundary condition
!> in subroutine \ref csc2cl.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nptdis
!> \param[in]     ivar          variable number
!> \param[in]     locpts
!> \param[in]     vela          variable value at time step beginning
!> \param[in]     coefav
!> \param[in]     coefbv
!> \param[in]     coopts
!> \param[out]    rvdis
!______________________________________________________________________________

subroutine cscpce &
 ( nptdis , ivar   ,                                              &
   locpts ,                                                       &
   vela   ,                                                       &
   coefav , coefbv ,                                              &
   coopts , rvdis  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use cplsat
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ivar
integer          nptdis

integer          locpts(nptdis)

double precision coopts(3,nptdis), rvdis(3,nptdis)
double precision coefav(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision vela(3,ncelet)

! Local variables

integer          ipt    , iel    , isou   , f_id
integer          inc    , iccocg , nswrgp
integer          iwarnp , imligp

double precision epsrgp , climgp
double precision dx     , dy     , dz

double precision, dimension(:,:,:), allocatable :: gradv

type(var_cal_opt) :: vcopt

!===============================================================================

! Allocate a temporary array
allocate(gradv(3,3,ncelet))

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

inc    = 1
iccocg = 1
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
iwarnp = vcopt%iwarni
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr

if (ivar.le.0) then
  f_id = -1
else
  f_id = ivarfl(ivar)
endif

call cgdvec &
( f_id   , imrgra , inc    , nswrgp , iwarnp , imligp ,          &
  epsrgp , climgp ,                                              &
  coefav , coefbv , vela   ,                                     &
  gradv)

! --- Interpolation

do ipt = 1, nptdis

  iel = locpts(ipt)

  dx = coopts(1,ipt) - xyzcen(1,iel)
  dy = coopts(2,ipt) - xyzcen(2,iel)
  dz = coopts(3,ipt) - xyzcen(3,iel)

  do isou = 1, 3
    rvdis(isou,ipt) = vela(isou,iel) + gradv(1,isou,iel)*dx       &
                                     + gradv(2,isou,iel)*dy       &
                                     + gradv(3,isou,iel)*dz
  enddo

enddo

! Free memory
deallocate(gradv)

!--------
! FORMATS
!--------
!----
! FIN
!----

return
end subroutine
