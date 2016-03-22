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
!> \file cs_user_radiative_transfer_bcs.f90
!> \brief User subroutine for input of radiative transfer parameters: boundary
!> conditions

!> See \subpage cs_user_radiative_transfer for examples.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
! _____________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     itypfb        boundary face types
!> \param[in]     icodcl        boundary condition code
!>                                - 1  -> Dirichlet
!>                                - 2  -> convective outelet
!>                                - 3  -> flux density
!>                                - 4  -> sliding wall and u.n=0 (velocity)
!>                                - 5  -> friction and u.n=0 (velocity)
!>                                - 6  -> roughness and u.n=0 (velocity)
!>                                - 9  -> free inlet/outlet (velocity)
!>                                inflowing possibly blocked
!> \param[in]     izfrdp        boundary faces -> zone number
!> \param[in]     isothp        boundary face type for radative transfer
!>                                - itpimp -> Gray wall with fixed inside temp
!>                                - ipgrno -> Gray wall with fixed outside temp
!>                                - iprefl -> Reflecting wall with fixed
!>                                         outside temp
!>                                - ifgrno -> Gray wall with fixed
!>                                      conduction flux
!>                                - ifrefl -> Reflecting wall with fixed
!>                                      conduction flux
!> \param[in]     tmin          min value of the wall temperature
!> \param[in]     tmax          max value of the wall temperature
!> \param[in]     tx            relaxation coefficient (0 < tx < 1)
!> \param[in]     dt            time step (per cell)
!> \param[in]     rcodcl        boundary condition values
!>                                rcodcl(3) = flux density value
!>                                (negative for gain) in W/m2
!> \param[in]     thwall        inside current wall temperature (K)
!> \param[in]     qincid        radiative incident flux  (W/m2)
!> \param[in]     hfcnvp        convective exchange coefficient (W/m2/K)
!> \param[in]     flcnvp        convective flux (W/m2)
!> \param[out]    xlamp         conductivity (W/m/K)
!> \param[out]    epap          thickness (m)
!> \param[out]    epsp          emissivity (>0)
!> \param[out]    textp         outside temperature (K)
!> \param[out]    tintp         initial inside temperature (K)
! _____________________________________________________________________________!

subroutine cs_user_radiative_transfer_bcs &
 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   icodcl , izfrdp , isothp ,                                     &
   tmin   , tmax   , tx     ,                                     &
   dt     , rcodcl ,                                              &
   thwall , qincid , hfcnvp , flcnvp ,                            &
   xlamp  , epap   , epsp   , textp  , tintp  )

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
use ppppar
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)
integer          icodcl(nfabor,nvarcl)
integer          izfrdp(nfabor), isothp(nfabor)

double precision tmin , tmax , tx

double precision dt(ncelet)
double precision rcodcl(nfabor,nvarcl,3)

double precision thwall(nfabor), qincid(nfabor)
double precision hfcnvp(nfabor),flcnvp(nfabor)
double precision xlamp(nfabor), epap(nfabor)
double precision epsp(nfabor)
double precision textp(nfabor), tintp(nfabor)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
  !==========
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in definition of boundary conditions'   ,/,&
'@    ======='                                                 ,/,&
'@  The user subroutine ''cs_user_radiative_transfer_bcs.f90'  ,/,&
'@  must be completed.'                                        ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================
! 0. Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


! -------
! Format
! -------


 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in definition of boundary conditions',/,   &
'@    =======',/,                                                 &
'@   Radiative data are missing for face: ',I10,/,                &
'@',/,                                                            &
'@   The user subroutine ''cs_user_radiative_transfer_bcs.f90' ,/,&
'@   must be completed.'                                       ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! ---
! End
! ---

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_radiative_transfer_bcs