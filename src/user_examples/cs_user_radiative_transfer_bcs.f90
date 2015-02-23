!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!< [loc_var]
! Local variables

integer          ifac , ivar, iok
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt
!< [loc_var]

!===============================================================================
! 0. Initialization
!===============================================================================

!< [allocate]
! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))
!< [allocate]

!===============================================================================
! 1. IVAR: number of the thermal variable
!===============================================================================

!< [ivar]
ivar = isca(iscalt)
!< [ivar]

!===============================================================================
!  2. Min and max values for the wall temperatures (clipping otherwise)
!   TMIN and TMAX are given in Kelvin.
!===============================================================================

!< [temp]
tmin = 0.d0
tmax = grand + tkelvi
!<[temp]

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

!< [forg]
iok = 0
!< [forg]

!   -------------------------------------------------------------------
!-->  Example 1:
!       For wall boundary faces, selection criteria: color 1
!       Gray or black wall with profil of fixed inside temperature
!       ------------------------------------

!< [example_1]
call getfbr('1',nlelt,lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if (itypfb(ifac).eq.iparoi) then

    ! zone number
    izfrdp(ifac) = 51

    ! Type of condition: gray or black wall with fixed inside temperature
    isothp(ifac) = itpimp

    ! Emissivity
    epsp  (ifac) = 0.1d0

     ! Profil of fixed inside temperature
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo
!< [example_1]

!   -------------------------------------------------------------------
!-->  Example 2 :
!       For wall boundary faces, selection criteria: color 2
!       Gray or black wall with fixed outside temperature TEXTP
!       ------------------------------------

!< [example_2]
call getfbr('2',nlelt,lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if (itypfb(ifac).eq.iparug) then

    ! zone number
    izfrdp(ifac) = 52

    ! Type of condition: gray or black wall with fixed outside temperature TEXTP
    isothp(ifac) = ipgrno

    ! Emissivity
    epsp  (ifac) = 0.9d0
    ! Conductivity (W/m/K)
    xlamp (ifac) = 3.0d0
    ! Thickness    (m)
    epap  (ifac) = 0.1d0
    ! Fixed outside temperature: 473.16 K
    textp (ifac) = 200.d0 + tkelvi
    ! Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo
!< [example_2]

!   -------------------------------------------------------------------
!-->  Exemple 3 :
!       For wall boundary faces, selection criteria: color 3
!       Reflecting wall (EPSP = 0) with fixed outside temperature TEXTP
!       ------------------------------------

!< [example_3]
call getfbr('3',nlelt,lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if (itypfb(ifac).eq.iparoi) then

    ! zone number
    izfrdp(ifac) = 53

    ! Type of condition: reflecting wall with fixed outside temperature TEXTP
    isothp(ifac) = iprefl

    ! Conductivity (W/m/K)
    xlamp (ifac) = 3.0d0
    ! Thickness    (m)
    epap  (ifac) = 0.1d0
    ! Fixed outside temperature: 473.16 K
    textp (ifac) = 200.d0 + tkelvi
    ! Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo
!< [example_3]

!   -------------------------------------------------------------------
!-->  Example 4 :
!      For wall boundary faces which have the color 4:
!           gray or black wall and fixed conduction flux through the wall

!        XLAMP
!        -----(Tparop-Textp) = fixed conduction flux     (W/m2)
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       If the conduction flux is zero then the wall is adiabatic.
!       The array RCODCL(IFAC,IVAR,3) has the value of the flux.
!       Flux density (< 0 if gain for the fluid)
!         For temperatures T,    in Watt/m2:
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         For enthalpies H,      in Watt/m2:
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

!< [example_4]
call getfbr('4',nlelt,lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if (itypfb(ifac).eq.iparoi) then

!      zone number
    izfrdp(ifac) = 54

    ! Type of condition: gray or black wall with fixed conduction flux through the wall
    isothp(ifac) = ifgrno

    ! Emissivity
    epsp  (ifac) = 0.9d0
    ! Conduction flux (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
    ! Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo
!< [example_4]

!   -------------------------------------------------------------------
!-->  Example 5 :
!      For wall boundary faces which have the color 5:
!           reflecting wall and fixed conduction flux through the wall

!      Equivalent to impose a Neumann to the fluid

!        XLAMP
!        -----(Tparop-Textp) = fixed conduction flux and EPSP = 0
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       If the conduction flux is zero then the wall is adiabatic.
!       Flux density (< 0 if gain for the fluid)
!         For temperatures T,    in Watt/m2:
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         For enthalpies H,      in Watt/m2:
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

!< [example_5]
call getfbr('5',nlelt,lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if (itypfb(ifac).eq.iparoi) then

    ! zone number
    izfrdp(ifac) = 55

    ! Type of condition: reflecting wall with fixed conduction flux through the wall
    isothp(ifac) = ifrefl

    ! Conduction flux (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
    ! Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo
!< [example_5]

!                           WARNING

!   -------------------------------------------------------------------
!-->   For all boundary faces that are not wall it is MANDATORY to
!      impose a number of zone in the array izfrdp.
!      For each zone, informations will be displayed in the listing.
!       ------------------------------------

!< [w]
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
!<[w]

!   -------------------------------------------------------------------
!-->
!      Verification that all boundary faces have been treated.
!       ------------------------------------

!< [check]
  elseif ( itypfb(ifac).eq.iparoi .or.                      &
           itypfb(ifac).eq.iparug     ) then
    if (izfrdp(ifac) .eq. -1) then
      write(nfecra,1000)ifac
      iok = iok + 1
    endif
  endif
!< [check]


!     End of the loop on the boundary faces
!     -------------------------------------

!< [end_radiative]
enddo

! Stop if there are forgotten faces
if(iok.ne.0) then
  call csexit (1)
  !==========
endif
!< [end_radiative]

! -------
! Format
! -------

!< [format_radiative]
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
!< [format_radiative]

! ---
! End
! ---

!< [deallocate]
! Deallocate the temporary array
deallocate(lstelt)
!< [deallocate]

return
end subroutine cs_user_radiative_transfer_bcs
