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

!> \file siream.f90
!> \brief Module for atmospheric aerosol chemistry in the atmospheric module

module siream

!=============================================================================

use ppppar, only: nozppm

!=============================================================================
!> \defgroup at_aerosol_chemistry Aerosol chemistry parameters for the
!>                                atmospheric module

!> \addtogroup at_aerosol_chemistry
!> \{

!> Flag to activate or not siream model
integer, save :: iaerosol
!> Flag to desactivate gaseous chemistry which is automatically
!> activated otherwise
integer, save :: inogaseouschemistry

! Number of aerosol species
!> - 2  --> inert
!> - 5  --> inorganics
!> - 13 --> organics and water
!> - 21 --> SIREAM
integer nesp_aer
parameter (nesp_aer = 21)

!> Number of gaseous species passed to SIREAM
integer, save ::  nespg_siream

!> Number of aerosol bins (can vary depending on the user)
integer nbin_aer
parameter (nbin_aer = 5)

!> Number of cycle in aerosol computation between ts and tf
integer ncycle_aer

!> Number of aerosol module options
integer noptions_aer
parameter (noptions_aer = 14)

!> 1D list of aerosol module options
integer, save :: options_aer(noptions_aer)

!> Names of particular species in SIREAM
character(len=10), dimension(nesp_aer) :: esp_siream

!> logical unit of the aerosol concentration profiles file
integer, save         ::  impmea
!> name of the aerosol concentration profiles file
character(len=10), save    ::  ficmea

!> Aerosol diameters at bin bounds
double precision, save :: bin_bound_aer(nbin_aer + 1)
!> Fixed aerosol density ([g/m^3])
double precision, save :: fixed_density_aer
!> Size variable aerosol density ([g/m^3])
double precision, save :: density_aer(nbin_aer)
!> Coagulation couples for each bin
integer, save :: couples_coag(nbin_aer)
!> First bin index of coagulation couples
integer, save :: first_index_coag(nbin_aer, 4 * nbin_aer)
!> Second bin index of coagulation couples
integer, save :: second_index_coag(nbin_aer, 4 * nbin_aer)
!> Coagulation partition coefficient
double precision, save :: coefficient_coag(nbin_aer, nbin_aer, nbin_aer)

!> Initial gaseous and particulate concentrations
!> and aerosol number read in file SIREAM
double precision, save, dimension(nesp_aer*nbin_aer+nbin_aer) :: dlconc0

!> read zone boundary conditions from file
integer, save :: iprofa(nozppm)

!> Flag to activate or not coagulation
integer, save :: icoag_siream
!> Flag to activate or not condensation/evaporation
integer, save :: icond_siream
!> Flag to activate or not nucleation
integer, save :: inucl_siream
!> Flag to consider or not kelvin effect
integer, save :: ikelv_siream
!> Cutting bin between equilibrium and dynamic bins
integer, save :: icut_siream
!> Sulfate condensation computation method
integer, save :: isulfcond_siream
!> Solver for dynamic condensation
integer, save :: kdslv_siream
!> Redistribution method of lagrangian bins
integer, save :: iredist_siream
!> Nucleation model
integer, save :: itern_siream
!> Method used for estimation of wet diameter
integer, save :: ithrm_siream

!> \}

!=============================================================================

contains

!=============================================================================

subroutine init_aerosols

implicit none

! Compute coagulation partition coefficients
call plug_compute_coagulation_coefficient(nbin_aer, bin_bound_aer,     &
     couples_coag, first_index_coag, second_index_coag, coefficient_coag)

end subroutine init_aerosols

!=============================================================================

end module siream
