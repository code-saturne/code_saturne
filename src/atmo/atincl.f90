!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

! Module for atmospheric module specific physics

module atincl

  !=============================================================================

  use ppppar

  !=============================================================================

  ! 1. Pointers specific to the atmospheric physics

  ! 1.1 Pointers specific to the input meteo profile (in ra)
  !-------------------------------------------------------------------------------
  !   Pointers specific to values read in the input meteo file:
  !                    itmmet ---> time (in sec) of the meteo profile
  !                    izdmet ---> altitudes of the dynamic profiles
  !                    iztmet ---> altitudes of the temperature profile
  !                    iumet, ivmet, iwmet  --> meteo u, v, w profiles
  !                    iekmet --->  meteo turbulent kinetic energy profile
  !                    iepmet ---> meteo turbulent dissipation profile
  !                    ittmet --->  meteo temperature (Celsius) profile
  !                    iqvmet ---> meteo specific humidity profile
  !                    ipmer  ---> Sea level pressure
  !                    ixmet, iymet --> cooordinates of the meteo profile

  !   Pointers specific to values calculated from the meteo file
  !   (cf atlecm.f90):
  !                    irmet -->  density profile
  !                    itpmet -->  potential temperature profile
  !                    iphmet -->  hydro. pressure from Laplace integration

  integer, save ::  itmmet,  &
                    izdmet,  &
                    iztmet,  &
                    iumet,   &
                    ivmet,   &
                    iwmet,   &
                    iekmet,  &
                    iepmet,  &
                    ittmet,  &
                    iqvmet,  &
                    ipmer,   &
                    ixmet,   &
                    iymet,   &
                    irmet,   &
                    itpmet,  &
                    iphmet

  ! 1.2 Pointers for the positions of the variables (in rtp, rtpa)
  !-------------------------------------------------------------------------------

  !   Variables specific to the atmospheric physics:
  !   ippmod(iatmos) = 1 (Dry atmosphere):
  !                    itempp---> potential temperature
  !   ippmod(IATMOS) = 2 (Humid atmosphere):
  !                    itempl---> liquid potential temperature
  !                    itotwt---> total water content
  !                    intdrp---> total number of droplets

  integer, save :: itempp, itempl, itotwt, intdrp

  ! 1.3 Pointers for the positions of the properties for the specific phys.
  !      (ipproc in propce, propfa, propfb)
  !-------------------------------------------------------------------------------

  !   Properties specific to the atmospheric physics:
  !   ippmod(iatmos) = 1 or 2 (Dry or Humid atmosphere):
  !                    itempc---> temperature (in celsius)
  !   ippmod(iatmos) = 2 (Humid atmosphere):
  !                    iliqwt---> liquid water content

  integer, save :: itempc,iliqwt

  !-------------------------------------------------------------------------------

  ! 2. Data specific to the atmospheric physics

  ! 2.1 Constant specific
  !-------------------------------------------------------------------------------
  ! rair --> perfect gas constant for air (mixture) defined in atini1

  double precision, save :: rair

  ! 2.2 Data specific to the input meteo profile
  !-------------------------------------------------------------------------------
  !                   imeteo --> flag for reading the meteo input file
  !                               = 0 -> no reading
  !                              = 1 -> reading
  !                   nbmetd --> numbers of altitudes for the dynamics
  !                   nbmett --> numbers of altitudes for the temperature
  !                                and specific humidity
  !                   nbmetm --> numbers of time steps for the meteo profiles
  !                   iprofm --> read zone boundary conditions from profile

  integer, save :: imeteo, nbmetd, nbmett, nbmetm, iprofm(nozppm)

  !=============================================================================

end module atincl
