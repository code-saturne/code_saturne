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
!> \brief User subroutine for input of radiative transfer parameters: boundary conditions

!> See \subpage cs_user_radiative_transfer for more information.

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

integer          ifac , ivar, iok
integer          ilelt, nlelt

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


!===============================================================================
! 1. IVAR: number of the thermal variable
!===============================================================================

ivar = isca(iscalt)

!===============================================================================
!  2. Min and max values for the wall temperatures (clipping otherwise)
!   TMIN and TMAX are given in Kelvin.
!===============================================================================

tmin = 0.d0
tmax = grand + tkelvi

!===============================================================================
! 3. Assign boundary conditions to boundary wall
!===============================================================================

!     ZONES DEFINITION
!     ================

!     We define zones of wall boundary, and we assign a type.
!       This allows to apply the boundary conditions and realize
!       balance sheets by treating them separately for each zone.

!     For each boundary face ifac (not just the faces of wall)
!       the user defines his own choice by a number of zone
!       IZFRDP(ifac) from color of the boundary face
!         or more generally, their properties (color, groups ...),
!         or boundary conditions specified in cs_user_boundary_conditions,
!         or even of their coordinates.
!     Warning: it is essential that ALL boundary faces
!       have been assigned to a zone.
!     The number of zones (the value of IZFRDP(ifac)) is
!       arbitrarily chosen by the user, but must be a
!       positive integer and less than or equal to NBZRDM
!       (value set in parameter radiat.h).



!     WALL CARACTERISTICS
!     ===================

!      WARNING: the unity of the temperature is the Kelvin
!      -------

!      Mandatory data:
!      ---------------
!      isothp(ifac) boundary face type
!                  = itpimp -> Gray wall with fixed inside temperature
!                  = ipgrno -> Gray wall with fixed outside temperature
!                  = iprefl -> Reflecting wall with fixed outside temperature
!                  = ifgrno -> Gray wall with fixed conduction flux
!                  = ifrefl -> Reflecting wall with fixed conduction flux

!      tintp(ifac) inside wall temperature (Kelvin)
!                  initialize thwall at the first time step.
!                  If isothp = itpimp, the value of thwall is fixed to tintp
!                  In the other case, tintp is only for initialization.


!      Other data (depend of the isothp):
!      ----------------------------------

!      rcodcl = conduction flux
!      epsp   = emissivity
!      xlamp  = conductivity (W/m/K)
!      epap   = thickness (m)
!      textp  = outside temperature (K)


!     EXAMPLE
!     =======

!        Wall boundary faces (IPAROI and IPARUG), are devided into 5 zones
!          located with IFRFAC(IFAC) in the range of number from 51 to 55.
!          For each location a different radiative boundary condition is applied.
!        For all other boundary that are not wall (i.e. inlet, oulet, symetry)
!          the user can define arbritay new zone using the array IFRFAC(IFAC),
!          for wich a value can be arbitrarily choosen between 1 and NBZRDM.
!
!     Warning: it is forbidden to modify thwall and qincid in this subroutine
!     ========

!    Indicator for forgotten faces.
iok = 0

!                           WARNING

!   -------------------------------------------------------------------
!-->   For all boundary faces that are not wall it is MANDATORY to
!      impose a number of zone in the array izfrdp.
!      For each zone, informations will be displayed in the listing.
!       ------------------------------------

do ifac = 1, nfabor

  if     (itypfb(ifac).eq.isolib) then
    izfrdp(ifac) = 60
  elseif (itypfb(ifac).eq.ifrent) then
    izfrdp(ifac) = 61
  elseif (itypfb(ifac).eq.ientre) then
    izfrdp(ifac) = 62
  elseif (itypfb(ifac).eq.i_convective_inlet) then
    izfrdp(ifac) = 63
  elseif (itypfb(ifac).eq.isymet) then
    izfrdp(ifac) = 64


!   -------------------------------------------------------------------
!-->
!      Verification that all boundary faces have been treated.
!       ------------------------------------

  elseif ( itypfb(ifac).eq.iparoi .or.                      &
           itypfb(ifac).eq.iparug     ) then
    if (izfrdp(ifac) .eq. -1) then
      write(nfecra,1000)ifac
      iok = iok + 1
    endif
  endif

!     End of the loop on the boundary faces
!     -------------------------------------

enddo

! Stop if there are forgotten faces
if(iok.ne.0) then
  call csexit (1)
  !==========
endif

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
