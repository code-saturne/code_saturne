!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_metal_structures_tag.f90
!>
!> \brief The 0-D thermal model to compute the temperature at the metal
!> structures wall and pass to the volume condensation modelling to be able to
!> model the metal structures effects.
!> This metal structures temperature computed is passed to the volume condensation
!> model to estimate the heat flux at the metall structures wall where
!> the condensation occurs.
!>
!> This subroutine is used to compute at each cell the
!> \f$T^{v}_{\mbox{metal}} \f$ at metal structures wall.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     dt            time step of the 1D thermal model
!_______________________________________________________________________________

subroutine cs_metal_structures_tag &
 ( ncmast , ltmast ,                          &
   dt     )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use radiat
use parall
use mesh
use field
use pointe, only:svcond, flxmst
use cs_tagms

!===============================================================================

implicit none

! Arguments

integer          ncmast, ltmast(ncelet)

double precision dt(ncelet)

! Local variables

integer          icmst
integer          iel

double precision xlcond   , flux    , vol_metal
double precision tw_metal , t_sym   , unstom
double precision tau_min  , tau_max
double precision tpminf   , tpmaxf  , tpmins , tpmaxs

!===============================================================================

!===============================================================================
! 0 - Initialization
!===============================================================================

! Vaporization latent heat
xlcond  = 2278.0d+3

tau_min = +1.d20 ; tau_max = -1.d20
tpminf  = +1.d20 ; tpmaxf  = -1.d20
tpmins  = +1.d20 ; tpmaxs  = -1.d20

! Metal structures volume
vol_metal = 0.d0
do icmst = 1, ncmast
  iel = ltmast(icmst)
  vol_metal = vol_metal + volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(vol_metal)
endif

!===============================================================================
! 1 - 0-D Thermal model with a explicit scheme
!===============================================================================

do icmst = 1,ncmast
  iel= ltmast(icmst)

  !-------------
  ! Fluid border
  !-------------

  ! Explicit flux recovered at the fluid frontier
  flux = ( flxmst(iel) - svcond(iel, ipr)*xlcond)    &
         /( s_metal*volume(iel)/vol_metal )

  ! temperature at the metal structures wall and the symmetry
  ! on the middle of the metal structure cell
  tw_metal = t_metal(iel,1)
  t_sym    = t_metal(iel,2)

  ! the characteristic time of heat propagation on
  ! a half thick metal structure
  unstom = xcond_m*s_metal/(xem/2.d0*m_metal*xcp_m/2.d0)

  tau_min = min(tau_min, 1.d0/unstom)
  tau_max = max(tau_max, 1.d0/unstom)

  ! Slove a 0-D unsteady conduction problem
  ! at the fluid and the symmetry frontiers
  ! with both equations given t_1 and t_2
  ! respectively  with t_1 the temperature
  ! past to the condensation correlations
  ! for the metal structures modelling.

  ! Compute t_1(n+1) near the mass wall
  t_metal(iel,1) = tw_metal                     &
                 + dt(iel)*unstom*              &
                    (  flux*xem/(2.d0*xcond_m)  &
                     + t_sym - tw_metal )

  ! Compute t_2(n+1) near the symmetry
  t_metal(iel,2) = t_sym                        &
                 + dt(iel)*unstom*(tw_metal-t_sym)
enddo

!===============================================================================
! 2 - Print Min/max values of temperatures and characteristic time scale
!===============================================================================

if( mod(ntcabs,ntlist).eq.0 ) then

  do icmst = 1, ncmast

    iel = ltmast(icmst)

    tpminf =min(tpminf,t_metal(iel,1))
    tpmaxf =max(tpmaxf,t_metal(iel,1))
    tpmins =min(tpmins,t_metal(iel,2))
    tpmaxs =max(tpmaxs,t_metal(iel,2))
  enddo

  if ( irangp .ge. 0 ) then
    call parmin(tpminf)
    call parmax(tpmaxf)
    call parmin(tpmins)
    call parmax(tpmaxs)
    call parmin(tau_min)
    call parmax(tau_max)
  endif

  write(nfecra,1000)
  write(nfecra,1001) ttcabs, tpminf, tpmaxf, &
                             tpmins, tpmaxs, &
                            tau_min, tau_max
  write(nfecra,1002)

endif

!--------
! Formats
!--------

 1000 format(/,&
          3x,'======================================== ',/, &
          3x,'Resolution of the 0-D thermal problem    ',/, &
          3x,' coupled with condensation correlations  ',/, &
          3x,'to model the metal structures effects    ',/, &
          3x,'======================================== ',/, &
             /,&
  3x,'------------------------------------------'   ,   &
     '------------------------------------'         ,/, &
     '------------------------------------'         ,/, &
    3x,' time', 8x,'Tp_fl (min) ',5x,'Tp_fl  (max)',6x, &
                 'Tp_sym(min) ',5x,'Tp_sym (max)'  ,/,  &
                 'tau   (min) ',5x,'tau    (max)'  ,/,  &
  3x,'  (s) ',8x, ' (C)       ' ,5x,' (C)        ',6x,  &
                  ' (C)       ' ,5x,' (C)        '  ,/, &
                  ' (-)       ' ,5x,' (-)        '  ,/, &
  3x,'------------------------------------------',      &
     '------------------------------------',            &
     '------------------------------------' )
 1001 format( 3x, 7(g15.7,1x) )
 1002 format(&
  3X,'------------------------------------------'   ,   &
  3x,'------------------------------------',            &
     '------------------------------------' )

!----
! End
!----

return

end subroutine cs_metal_structures_tag
