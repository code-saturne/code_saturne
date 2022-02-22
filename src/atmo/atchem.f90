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

!> \file atchem.f90
!> \brief Module for chemistry in the atmospheric module

module atchem

  !=============================================================================

  use, intrinsic :: iso_c_binding
  use ppppar, only: nozppm

  !=============================================================================

  !> \defgroup at_gaseous_chemistry Gaseous chemistry parameters for the
  !>                                atmospheric module

  !> \addtogroup at_gaseous_chemistry
  !> \{
  ! Useful constants for chemistry
  !> Avogadro constant (molecules/mol)
  double precision :: navo
  !> Molar mass of dry air constant (Kg/mol)
  double precision :: Mair
  parameter (navo = 6.022d+23)         ! Molecules/mol
  parameter (Mair = 28.9d-3)           ! Kg/mol

  !> Choice of chemistry resolution scheme
  !> - 0 --> no atmospheric chemistry
  !> - 1 --> quasi steady equilibrium NOx scheme with 4 species and 5 reactions
  !> - 2 --> scheme with 20 species and 34 reactions
  !> - 3 --> scheme CB05 with 52 species and 155 reactions
  !> - 4 --> user defined schema
  integer(c_int), pointer, save :: ichemistry

  !> ifilechemistry: choice to read (=1,2,3,4, according to the scheme) or not (0)
  !> a concentration profile file
  integer, save :: ifilechemistry
  !> isepchemistry: splitted (=1) or semi-coupled (=2, pu-sun) resolution
  !> of chemistry
  integer, save :: isepchemistry
  !> photolysis: inclusion (true) or not (false) of photolysis reactions
  logical(kind=c_bool), pointer, save :: photolysis
  !> Number of chemical species
  integer(c_int), pointer, save :: nespg
  !> Number of chemical reactions
  integer(c_int), pointer, save :: nrg

  !> scalar id for chemical species
  integer(c_int), dimension(:), pointer, save ::  isca_chem
  !> Molar mass of chemical species (g/mol)
  double precision, dimension(:), pointer ::  dmmk
  !> pointer to deal with different orders of chemical species
  integer(c_int), dimension(:), pointer ::  chempoint
  !> conversion factors for reaction rates Jacobian matrix
  double precision, allocatable, dimension(:) ::  conv_factor_jac
  !> kinetics constants
  double precision, allocatable, dimension(:) ::  reacnum

  !> maximal time step for chemistry resolution
  double precision dtchemmax

  !> number of time steps for the concentration profiles file
  integer, save         ::  nbchim
  !> number of altitudes for the concentration profiles file
  integer, save         ::  nbchmz
  !> number of initialized chemical species in the concentration profiles file
  integer, save         ::  nespgi = 0

  !> indices of chemical species in the concentration profiles file
  integer, allocatable, dimension(:)          :: idespgi
  !> concentration profiles
  double precision, allocatable, dimension(:) :: espnum
  !> altitudes of the concentration profiles
  double precision, allocatable, dimension(:) :: zproc
  !> time steps of the concentration profiles
  double precision, allocatable, dimension(:) :: tchem
  !> X coordinates of concentration profiles
  double precision, allocatable, dimension(:) :: xchem
  !> Y coordinates of concentration profiles
  double precision, allocatable, dimension(:) :: ychem
  !> read zone boundary conditions from profile
  integer, save :: iprofc(nozppm)

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function returning a chemistry concentration file name

    subroutine cs_f_atmo_get_chem_conc_file_name(f_name_max, f_name, f_name_len)  &
      bind(C, name='cs_f_atmo_get_chem_conc_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_name_max
      type(c_ptr), intent(out)    :: f_name
      integer(c_int), intent(out) :: f_name_len
    end subroutine cs_f_atmo_get_chem_conc_file_name

    !---------------------------------------------------------------------------

    ! Interface to C function returning a aerosol concentration file name

    subroutine cs_f_atmo_get_aero_conc_file_name(f_name_max, f_name, f_name_len)  &
      bind(C, name='cs_f_atmo_get_aero_conc_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_name_max
      type(c_ptr), intent(out)    :: f_name
      integer(c_int), intent(out) :: f_name_len
    end subroutine cs_f_atmo_get_aero_conc_file_name

    !---------------------------------------------------------------------------

    !> \brief Deallocate arrays for atmo chemistry

    subroutine cs_f_atmo_chem_finalize() &
      bind(C, name='cs_f_atmo_chem_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_atmo_chem_finalize

    !---------------------------------------------------------------------------


  end interface

  !=============================================================================

contains

  !-------------------------------------------------------------------------------

  !> \brief Return chemistry concentration file name

  !> \param[out]  name   chemistry concentration file name

  subroutine atmo_get_chem_conc_file_name(name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(out) :: name

    ! Local variables

    integer :: i
    integer(c_int) :: name_max, c_name_len
    type(c_ptr) :: c_name_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_name

    name_max = len(name)

    call cs_f_atmo_get_chem_conc_file_name(name_max, c_name_p, c_name_len)
    call c_f_pointer(c_name_p, c_name, [c_name_len])

    do i = 1, c_name_len
      name(i:i) = c_name(i)
    enddo
    do i = c_name_len + 1, name_max
      name(i:i) = ' '
    enddo

    return

  end subroutine atmo_get_chem_conc_file_name

  !=============================================================================

  !> \brief Return aerosol concentration file name

  !> \param[out]  name   aerosol concentration file name

  subroutine atmo_get_aero_conc_file_name(name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(out) :: name

    ! Local variables

    integer :: i
    integer(c_int) :: name_max, c_name_len
    type(c_ptr) :: c_name_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_name

    name_max = len(name)

    call cs_f_atmo_get_aero_conc_file_name(name_max, c_name_p, c_name_len)
    call c_f_pointer(c_name_p, c_name, [c_name_len])

    do i = 1, c_name_len
      name(i:i) = c_name(i)
    enddo
    do i = c_name_len + 1, name_max
      name(i:i) = ' '
    enddo

    return

  end subroutine atmo_get_aero_conc_file_name

  !=============================================================================

  !> \brief Allocate some atmoshperic chemistry arrays
  subroutine init_chemistry

    implicit none

    integer imode

    ! First reading of concentration profiles file
    imode = 0

    call atlecc(imode)

    ! Dynamical allocations

    allocate(conv_factor_jac(nespg*nespg))
    allocate(idespgi(nespgi))
    allocate(espnum(nespg*nbchim*nbchmz))
    allocate(zproc(nbchmz))
    allocate(tchem(nbchim))
    allocate(xchem(nbchim))
    allocate(ychem(nbchim))

  end subroutine init_chemistry

  !=============================================================================

  !> \brief Allocate memory relative to mesh size

  subroutine init_chemistry_reacnum

    use mesh, only: ncel

    implicit none

    ! Dynamical allocations
    allocate(reacnum(ncel*nrg))

  end subroutine init_chemistry_reacnum

  !=============================================================================

  !> \brief Initialize species_to_field_id
  subroutine cs_atmo_chem_init_c_chemistry

    use numvar, only : ivarfl, isca
    use cs_c_bindings

    implicit none

    ! Local variables
    integer i
    integer(c_int), dimension(nespg) :: c_species_to_fid

    do i = 1, nespg
      c_species_to_fid(i) = ivarfl(isca(isca_chem(i)))
    enddo

    call cs_f_atmo_chem_initialize_species_to_fid(c_species_to_fid)

  end subroutine cs_atmo_chem_init_c_chemistry

  !=============================================================================

  !> \brief Map pointers to arrays

  subroutine init_chemistry_pointers

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables

    type(c_ptr) :: c_species_to_scalar_id, c_molar_mass, c_chempoint

    call cs_f_atmo_chem_arrays_get_pointers(c_species_to_scalar_id, &
                                            c_molar_mass, &
                                            c_chempoint)

    call c_f_pointer(c_species_to_scalar_id, isca_chem, [nespg])
    call c_f_pointer(c_molar_mass, dmmk, [nespg])
    call c_f_pointer(c_chempoint, chempoint, [nespg])

  end subroutine init_chemistry_pointers

  !=============================================================================

  !> \brief deallocate the space

  subroutine finalize_chemistry

    use cs_c_bindings

    implicit none

    call cs_f_atmo_chem_finalize()

    deallocate(conv_factor_jac)
    deallocate(reacnum)
    deallocate(idespgi)
    deallocate(espnum)
    deallocate(zproc)
    deallocate(tchem)
    deallocate(xchem)
    deallocate(ychem)

  end subroutine finalize_chemistry

  !=============================================================================

end module atchem
