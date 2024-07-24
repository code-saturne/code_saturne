!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file pointe.f90
!> \brief Module for pointer variables

module pointe

  !=============================================================================

  use, intrinsic :: iso_c_binding
  use paramx

  implicit none

  !=============================================================================
  !> \defgroup pointer_variables Module for pointer variables

  !> \addtogroup pointer_variables
  !> \{

  !> \defgroup fortran_pointer_containers Containers to Fortran array pointers.
  !> An array of one of these derived types can be used to manage a set of
  !> pointers (Fortran not allow arrays of pointers directly, so this
  !> technique is a classical workaround.

  ! Note also that Fortran bounds remapping could in theory be used to
  ! handle different pointer shapes with a single type.

  !> \addtogroup fortran_pointer_containers
  !> \{

  !> container for rank 1 double precision array pointer.
  type pmapper_double_r1
    double precision, dimension(:),  pointer :: p !< rank 1 array pointer
  end type pmapper_double_r1

  !> \}

  !=============================================================================

  !> \defgroup coupled_case Specific arrays for the coupled case

  !> \addtogroup coupled_case
  !> \{

  !> \anchor itypfb
  !> boundary condition type at the boundary face \c ifac
  !> (see \ref cs_user_boundary_conditions)
  integer, dimension(:), pointer, save :: itypfb

  !> indirection array allowing to sort the boundary faces
  !> according to their boundary condition type \c itypfb
  !integer, allocatable, dimension(:) :: itrifb
  integer, dimension(:), pointer, save :: itrifb

  !> to identify boundary zones associated with boundary faces
  !> (specific physics models)
  integer, dimension(:), pointer :: izfppp

  !> \}

  !=============================================================================

  !> \defgroup porosity_ibm Porosity from immersed boundaries parameters

  !> \addtogroup porosity_ibm
  !> \{

  !> Activate the computation
  integer(c_int), pointer, save :: ibm_porosity_mode

  !> \}

  !=============================================================================

  !> \defgroup porosity_from_scan Porosity from scan module parameters

  !> \addtogroup porosity_from_scan
  !> \{

  !> Activate the computation
  logical(c_bool), pointer, save :: compute_porosity_from_scan

  !> \}

  !=============================================================================

  !... Auxiliaires
  !> \addtogroup auxiliary
  !> \{

  !> liquid-vapor mass transfer term for cavitating flows
  !> and its derivative with respect to pressure
  double precision, allocatable, target, dimension(:) :: gamcav, dgdpca

  !> reference point for wall condensation,
  !> used in forced and mixed convection regimes
  double precision, allocatable, dimension(:,:) :: xref_cond

  !> \}

  !> \}

contains

  !=============================================================================

  ! Initialize auxiliary arrays

  subroutine init_aux_arrays() &
    bind(C, name='cs_f_init_aux_arrays')

    use mesh, only: ncelet, nfabor
    use paramx
    use optcal
    use cs_c_bindings

    implicit none

    ! Arguments

    ! Local variables

    ! Boundary-face related arrays

    allocate(itrifb(nfabor))

    ! liquid-vapor mass transfer term for cavitating flows
    ! and its part implicit in pressure
    if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) then
      allocate(gamcav(ncelet), dgdpca(ncelet))
    endif

  end subroutine init_aux_arrays

  !=============================================================================

  ! Resize auxiliary arrays

  subroutine resize_aux_arrays() &
   bind(C, name='cs_fortran_resize_aux_arrays')
   use, intrinsic :: iso_c_binding
   use cs_c_bindings

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    ! Local variables

    integer iel
    double precision, allocatable, dimension(:) :: buffer

    ! Resize/copy arrays

    allocate(buffer(ncelet))

    ! liquid-vapor mass transfer term for cavitating flows
    ! and its part implicit in pressure

    if (allocated(gamcav)) then
      do iel = 1, ncel
        buffer(iel) = gamcav(iel)
      enddo
      deallocate(gamcav)
      call synsca (buffer)
      allocate(gamcav(ncelet))
      do iel = 1, ncelet
        gamcav(iel) = buffer(iel)
      enddo

      do iel = 1, ncel
        buffer(iel) = dgdpca(iel)
      enddo
      deallocate(dgdpca)
      call synsca (buffer)
      allocate(dgdpca(ncelet))
      do iel = 1, ncelet
        dgdpca(iel) = buffer(iel)
      enddo
    endif

    deallocate(buffer)

  end subroutine resize_aux_arrays

  !=============================================================================

  ! Free auxiliary arrays

  subroutine finalize_aux_arrays() &
    bind(C, name='cs_f_finalize_aux_arrays')

    deallocate(itrifb)
    if (allocated(gamcav)) deallocate(gamcav, dgdpca)

  end subroutine finalize_aux_arrays

  !=============================================================================

  subroutine boundary_conditions_init() &
    bind(C, name='cs_f_boundary_conditions_init')

    use, intrinsic :: iso_c_binding
    use mesh
    use cs_c_bindings

    implicit none

    ! Local variables

    type(c_ptr) :: c_itypfb, c_izfppp, c_itrifb

    call cs_f_boundary_conditions_create

    call cs_f_boundary_conditions_get_pointers(c_itypfb, c_izfppp, c_itrifb)

    call c_f_pointer(c_itypfb, itypfb, [nfabor])
    call c_f_pointer(c_izfppp, izfppp, [nfabor])
    call c_f_pointer(c_itrifb, itrifb, [nfabor])

  end subroutine boundary_conditions_init

  !=============================================================================

  !> \brief Allocate the cs_glob_porosity_ibm structure.

  subroutine porosity_ibm_init

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_ibm_porosity_mode

    call cs_f_porosity_ibm_get_pointer(c_ibm_porosity_mode)

    call c_f_pointer(c_ibm_porosity_mode, ibm_porosity_mode)

    return

  end subroutine porosity_ibm_init

  !=============================================================================

  !> \brief Allocate the cs_glob_porosity_from_scan structure.

  subroutine porosity_from_scan_init

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_compute_from_scan

    call cs_f_porosity_from_scan_get_pointer(c_compute_from_scan)

    call c_f_pointer(c_compute_from_scan, compute_porosity_from_scan)

  end subroutine porosity_from_scan_init

  !=============================================================================


  !> \brief Return C pointer to cavitation "dgdpca" array

  !> \return  cav_dgpd

  function cs_get_cavitation_dgdp_st() result(cav_dgpd) &
    bind(C, name='cs_get_cavitation_dgdp_st')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr) :: cav_dgpd
    cav_dgpd = c_loc(dgdpca)
  end function cs_get_cavitation_dgdp_st

  !=============================================================================

  !> \brief Return C pointer to cavitation "gamcav" array

  !> \return  cav_gam

  function cs_get_cavitation_gam() result(cav_gam) &
    bind(C, name='cs_get_cavitation_gam')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr) :: cav_gam
    cav_gam = c_loc(gamcav)
  end function cs_get_cavitation_gam

  !=============================================================================

end module pointe
