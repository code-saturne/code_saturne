!-------------------------------------------------------------------------------

!VERS


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

!                              atincl.h

!===============================================================================

!            INCLUDE FOR THE ATMOSPHERIC MODULE SPECIFIC PHYSICS

!-------------------------------------------------------------------------------

! 1. Pointers specific to the atmospheric physics

! 1.1 Pointers specific to the input meteo profile (in ra)
!-------------------------------------------------------------------------------
!         /IPROM/
!   Pointers specific to values read in the input meteo file:
!                    ITMMET ---> time (in sec) of the meteo profile
!                    IZDMET ---> altitudes of the dynamic profiles
!                    IZTMET ---> altitudes of the temperature profile
!                    IUMET, IVMET, IWMET  --> meteo u, v, w profiles
!                    IEKMET --->  meteo turbulent kinetic energy profile
!                    IEPMET ---> meteo turbulent dissipation profile
!                    ITTMET --->  meteo temperature (Celsius) profile
!                    IQVMET ---> meteo specific humidity profile
!                    IPMER  ---> Sea level pressure
!                    IXMET, IYMET --> cooordinates of the meteo profile

!   Pointers specific to values calculated from the meteo file
!   (cf atlecm.f90):
!                    IRMET -->  density profile
!                    ITPMET -->  potential temperature profile
!                    IPHMET -->  hydro. pressure from Laplace integration

integer        itmmet,                                            &
               izdmet,                                            &
               iztmet,                                            &
               iumet,                                             &
               ivmet,                                             &
               iwmet,                                             &
               iekmet,                                            &
               iepmet,                                            &
               ittmet,                                            &
               iqvmet,                                            &
               ipmer,                                             &
               ixmet,                                             &
               iymet,                                             &
               irmet,                                             &
               itpmet,                                            &
               iphmet

common /iprom/ itmmet,                                            &
               izdmet,                                            &
               iztmet,                                            &
               iumet,                                             &
               ivmet,                                             &
               iwmet,                                             &
               iekmet,                                            &
               iepmet,                                            &
               ittmet,                                            &
               iqvmet,                                            &
               ipmer,                                             &
               ixmet,                                             &
               iymet,                                             &
               irmet,                                             &
               itpmet,                                            &
               iphmet

! 1.2 Pointers for the positions of the variables (in rtp, rtpa)
!-------------------------------------------------------------------------------

!         /IATVAR/
!   Variables specific to the atmospheric physics:
!   IPPMOD(IATMOS) = 1 (Dry atmosphere):
!                    ITEMPP---> potential temperature
!   IPPMOD(IATMOS) = 2 (Humid atmosphere):
!                    ITEMPL---> liquid potential temperature
!                    ITOTWT---> total water content
!                    INTDRP---> total number of droplets

integer        itempp, itempl, itotwt, intdrp

common /iatvar/ itempp, itempl, itotwt, intdrp

! 1.3 Pointers for the positions of the properties for the specific phys.
!      (ipproc in propce, propfa, propfb)
!-------------------------------------------------------------------------------

!         /IATPRO/
!   Properties specific to the atmospheric physics:
!   IPPMOD(IATMOS) = 1 or 2 (Dry or Humid atmosphere):
!                    ITEMPC---> temperature (in celsius)
!   IPPMOD(IATMOS) = 2 (Humid atmosphere):
!                    ILIQWT---> liquid water content

integer        itempc,iliqwt

common /iatpro/ itempc,iliqwt


!-------------------------------------------------------------------------------

! 2. Data specific to the atmospheric physics

! 2.1 Constant specific
!-------------------------------------------------------------------------------
!          /ATMCST/ RAIR --> perfect gaz constant for air (mixture)
!                            defined in atini1

 double precision rair
 common /atmcst/ rair

! 2.2 Data specific to the input meteo profile
!-------------------------------------------------------------------------------
!          /PROMET/
!                   IMETEO --> flag for reading the meteo input file
!                               = 0 -> no reading
!                              = 1 -> reading
!                   NBMETD --> numbers of altitudes for the dynamics
!                   NBMETT --> numbers of altitudes for the temperature
!                                and specific humidity
!                   NBMETM --> numbers of time steps for the meteo profiles
!                   IPROFM --> read zone boundary conditions from profile

  integer   imeteo, nbmetd, nbmett, nbmetm,                       &
            iprofm(nozppm)

  common / iimet / imeteo, nbmetd, nbmett,nbmetm, iprofm
