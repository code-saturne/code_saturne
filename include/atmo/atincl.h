c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C                              atincl.h
C
C***********************************************************************
C
C            INCLUDE FOR THE ATMOSPHERIC MODULE SPECIFIC PHYSICS
C
C-----------------------------------------------------------------------
C
C 1. Pointers specific to the atmospheric physics
C
C ----------------------------------------------------------------------
C
C 1.1 Pointers specfic to the input meteo profile
C-----------------------------------------------------------------------
C         /IPROM/
C   Pointers specific to values read in the input meteo file:
C                    ITMMET ---> time (in sec) of the meteo profile
C                    IZDMET ---> altitudes of the dynamic profiles
C                    IZTMET ---> altitudes of the temperature profile
C                    IUMET, IVMET, IWMET  --> meteo u, v, w profiles
C                    IEKMET --->  meteo turbulent kinetic energy profile
C                    IEPMET ---> meteo turbulent dissipation profile
C                    ITTMET --->  meteo temperature (Celsius) profile
C                    IQVMET ---> meteo specific humidity profile
C                    IPMER  ---> Sea level pressure
C                    IXMET, IYMET --> cooordinates of the meteo profile
C
C   Pointers specific to values calculated from the meteo file
C   (cf ATLECM.F):
C                    IRMET -->  density profile
C                    ITPMET -->  potential temperature profile
C                    IPHMET -->  hydro. pressure from Laplace integration
C
      INTEGER        ITMMET,
     &               IZDMET,
     &               IZTMET,
     &               IUMET,
     &               IVMET,
     &               IWMET,
     &               IEKMET,
     &               IEPMET,
     &               ITTMET,
     &               IQVMET,
     &               IPMER,
     &               IXMET,
     &               IYMET,
     &               IRMET,
     &               ITPMET,
     &               IPHMET
C
      COMMON /IPROM/ ITMMET,
     &               IZDMET,
     &               IZTMET,
     &               IUMET,
     &               IVMET,
     &               IWMET,
     &               IEKMET,
     &               IEPMET,
     &               ITTMET,
     &               IQVMET,
     &               IPMER,
     &               IXMET,
     &               IYMET,
     &               IRMET,
     &               ITPMET,
     &               IPHMET


C-----------------------------------------------------------------------
C
C 2. Data specific to the atmospheric physics
C
C ----------------------------------------------------------------------
C
C 2.1 Data specific to the input meteo profile
C-----------------------------------------------------------------------
C          /PROMET/
C                   IMETEO --> flag for reading the meteo input file
C                               = 0 -> no reading
C                              = 1 -> reading
C                   NBMETD --> numbers of altitudes for the dynamics
C                   NBMETT --> numbers of altitudes for the temperature
C                                and specific humidity
C                   NBMETM --> numbers of time steps for the meteo profiles
C                   IPROFM --> read zone boundary conditions from profile
C
        INTEGER   IMETEO, NBMETD, NBMETT, NBMETM,
     &            IPROFM(NOZPPM)
C
        COMMON / IIMET / IMETEO, NBMETD, NBMETT,NBMETM, IPROFM
