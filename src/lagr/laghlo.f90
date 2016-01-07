!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
! --------

!> \file laghlo.f90
!> \brief Define Head losses to take into account deposit in the flow
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncepdc        number of cells in head loss zone
!> \param[in]     icepdc        numbers of ncepdp cells with head loss
!> \param[out]    ckupdc        head loss
!_______________________________________________________________________________


subroutine laghlo &
 ( ncepdc , icepdc , ckupdc )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use parall
use period
use mesh
use field
use ppincl
use entsor

!===============================================================================

implicit none

! Arguments

integer          ncepdc
integer          icepdc(ncepdc)
double precision ckupdc(ncepdc,6)

! Local variables

integer          iel, ielpdc, poro_id

double precision v, ck
double precision romf, visccf, lcell

double precision, dimension(:), allocatable :: mdiam

double precision, dimension(:,:), pointer :: cvara_vel
double precision, dimension(:), pointer :: lporo

double precision, dimension(:), pointer :: cromf
double precision, dimension(:), pointer :: viscl

!===============================================================================

! Map field arrays

call field_get_val_prev_v(ivarfl(iu), cvara_vel)

! Check that head loss zone definitions are consistent

if (ncepdc.ne.ncel) then
  write(nfecra,1001)
  call csexit(1)
endif

!=============================================================================

!    ckupdc: compute head loss coefficients in the calculation coordinates,
!            organized in order k11, k22, k33, k12, k13, k23

! Note:
!
!    - make sure diagonal coefficients are positive. The calculation
!      may crash if this is not the case, and no further check will
!      be done

! Pointer on the density w.r.t the flow

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

call field_get_val_s(iprpfl(iviscl), viscl)

!===============================================================================
! Porosity calculation for the influence of the deposit on the flow
! by head losses
!===============================================================================

allocate(mdiam(ncelet))

call field_get_id_try('clogging_porosity', poro_id)

if (poro_id .lt.0) then
  allocate(lporo(ncelet))
else
  call field_get_val_s(poro_id, lporo)
endif

call porcel(mdiam, lporo)

!===============================================================================

! Calculation of the head loss term with the Ergun law
! mdiam :  mean diameter of deposited particles
! lcell :  characteristic length in the flow direction

do ielpdc = 1, ncepdc
  iel = icepdc(ielpdc)
  if (mdiam(iel) .gt. 0.d0) then
    lcell = (volume(iel))**(1.d0/3.d0)
    romf = cromf(iel)
    visccf = viscl(iel) / romf
    v = sqrt(cvara_vel(1,iel)**2 + cvara_vel(2,iel)**2 + cvara_vel(3,iel)**2)
    ck =     v * 1.75d0 * (1 - lporo(iel)) / lporo(iel)**3.d0         &
           * lcell / mdiam(iel)                                       &
         +   (lcell * 150.d0 * visccf ) /  (romf * mdiam(iel)**2)     &
           * (1 - lporo(iel))**2 / lporo(iel)*3
    ckupdc(iel,1) = ck
    ckupdc(iel,2) = ck
    ckupdc(iel,3) = ck
    ckupdc(iel,4) = 0.d0
    ckupdc(iel,5) = 0.d0
    ckupdc(iel,6) = 0.d0
  endif
enddo

if (poro_id .lt.0) deallocate(lporo)
deallocate(mdiam)

return

 1001 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:',                                                  /,&
'@    ======',                                                  /,&
'@   TO BE COMPATIBLE WITH THE LAGRANGIAN DEPOSITION MODEL,'    /,&
'@     HEAD LOSS ZONES MUST COVER THE WHOLE MESH',              /,&
'@ Head loss coefficiets may be locally zero.'                  /,&
'@ Check your case setup.',                                     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

end subroutine laghlo
