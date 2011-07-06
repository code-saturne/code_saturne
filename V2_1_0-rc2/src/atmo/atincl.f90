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

  ! 1.1 Arrays specific to the input meteo profile
  !-----------------------------------------------

  !   Arrays specific to values read in the input meteo file:
  !                    tmmet ---> time (in sec) of the meteo profile
  !                    zdmet ---> altitudes of the dynamic profiles
  !                    ztmet ---> altitudes of the temperature profile
  !                    umet, vmet, wmet  --> meteo u, v, w profiles
  !                    ekmet --->  meteo turbulent kinetic energy profile
  !                    epmet ---> meteo turbulent dissipation profile
  !                    ttmet --->  meteo temperature (Celsius) profile
  !                    qvmet ---> meteo specific humidity profile
  !                    pmer  ---> Sea level pressure
  !                    xmet, ymet --> cooordinates of the meteo profile

  !   Arrays specific to values calculated from the meteo file
  !   (cf atlecm.f90):
  !                    rmet -->  density profile
  !                    tpmet -->  potential temperature profile
  !                    phmet -->  hydro. pressure from Laplace integration

  double precision, allocatable, dimension(:) :: tmmet, zdmet, ztmet
  double precision, allocatable, dimension(:,:) :: umet, vmet, wmet
  double precision, allocatable, dimension(:,:) :: ekmet, epmet, ttmet, qvmet
  double precision, allocatable, dimension(:) :: pmer
  double precision, allocatable, dimension(:) :: xmet, ymet
  double precision, allocatable, dimension(:,:) :: rmet, tpmet, phmet

  ! 1.2 Pointers for the positions of the variables (in rtp, rtpa)
  !---------------------------------------------------------------

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
  !------------------------------------------------------------------------

  !   Properties specific to the atmospheric physics:
  !   ippmod(iatmos) = 1 or 2 (Dry or Humid atmosphere):
  !                    itempc---> temperature (in celsius)
  !   ippmod(iatmos) = 2 (Humid atmosphere):
  !                    iliqwt---> liquid water content

  integer, save :: itempc, iliqwt

  !----------------------------------------------------------------------------

  ! 2. Data specific to the atmospheric physics

  ! 2.1 Data specific to the input meteo profile
  !----------------------------------------------
  !                   imeteo --> flag for reading the meteo input file
  !                               = 0 -> no reading
  !                              = 1 -> reading
  !                   nbmetd --> numbers of altitudes for the dynamics
  !                   nbmett --> numbers of altitudes for the temperature
  !                                and specific humidity
  !                   nbmetm --> numbers of time steps for the meteo profiles
  !                   iprofm --> read zone boundary conditions from profile

  integer, save :: imeteo, nbmetd, nbmett, nbmetm, iprofm(nozppm)

contains

  !=============================================================================

  subroutine init_meteo

    implicit none

    integer :: imode
    double precision :: rvoid(1)

    if (imeteo.gt.0) then

      imode = 0

      !     Nb les arguments ne sont pas utilises quand IMODE=0
      call atlecm &
      !==========
    ( imode ,                  &
      nbmetd, nbmett, nbmetm,  &
      rvoid , rvoid , rvoid ,  &
      rvoid , rvoid , rvoid ,  &
      rvoid , rvoid ,          &
      rvoid , rvoid ,          &
      rvoid , rvoid ,          &
      rvoid , rvoid , rvoid )

      allocate(tmmet(nbmetm), zdmet(nbmetd), ztmet(nbmett))
      allocate(umet(nbmetd,nbmetm), vmet(nbmetd,nbmetm), wmet(nbmetd,nbmetm))
      allocate(ekmet(nbmetd,nbmetm), epmet(nbmetd,nbmetm))
      allocate(ttmet(nbmett,nbmetm), qvmet(nbmett, nbmetm))
      allocate(pmer(nbmetm))
      allocate(xmet(nbmetm), ymet(nbmetm))
      allocate(rmet(nbmett,nbmetm), tpmet(nbmett,nbmetm), phmet(nbmett,nbmetm))

    endif

  end subroutine init_meteo

  !=============================================================================

  subroutine finalize_meteo

    implicit none

    if (imeteo.gt.0) then

      deallocate(tmmet, zdmet, ztmet)
      deallocate(umet, vmet, wmet)
      deallocate(ekmet, epmet)
      deallocate(ttmet, qvmet)
      deallocate(pmer)
      deallocate(xmet, ymet)
      deallocate(rmet, tpmet, phmet)

    endif

  end subroutine finalize_meteo

end module atincl
