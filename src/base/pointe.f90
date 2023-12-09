!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

  !> container for rank 2 double precision array pointer.
  type pmapper_double_r2
    double precision, dimension(:,:),  pointer :: p !< rank 2 array pointer
  end type pmapper_double_r2

  !> \}

  !=============================================================================

  !> \defgroup dummy_arrays Dummy target arrays for null pointers

  !> \addtogroup dummy_arrays
  !> \{

  integer, dimension(1),   target :: ivoid1
  integer, dimension(1,1), target :: ivoid2

  double precision, dimension(1),     target :: rvoid1
  double precision, dimension(1,1),   target :: rvoid2

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

  !> \defgroup thermal_1D Thermal 1D module parameters

  !> \addtogroup thermal_1D
  !> \{

  !> number of boundary faces which are coupled
  !> with a wall 1D thermal module. See the user subroutine
  !> \ref cs_user_1d_wall_thermal
  integer(c_int), pointer, save :: nfpt1d

  !> global number of boundary faces which are coupled with
  !> a wall 1D thermal module. (ie sum over all ranks of nfpt1d)
  integer(c_int), pointer, save :: nfpt1t

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

  !> \anchor ncepdc
  !> number of cells in which a pressure drop is imposed.
  integer, save :: ncepdc

  !> \anchor icepdc
  !> number of the \c ncepdc cells in which a pressure drop is imposed.
  !> See \c {iicepd}
  integer, allocatable, dimension(:) :: icepdc

  !> \anchor ckupdc
  !> value of the coefficients of the pressure drop tensor of the
  !> \c ncepdc cells in which a pressure drop is imposed.
  !> Note the 6 values are interleaved as follows: (k11, k22, k33, k12, k23, k13).
  !> See \c ickpdc
  real(c_double), allocatable, dimension(:,:), target :: ckupdc
  type(c_ptr) :: p_ckupdc = c_null_ptr
  bind(C, name='cs_glob_ckupdc') :: p_ckupdc

  !> \anchor ncetsm
  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm also
  integer, save :: ncetsm

  !> \anchor icetsm
  !> number of the \c ncetsm cells in which a mass injection is imposed.
  !> See \c iicesm and the \c cs_equation_add_volume_mass_injection_* functions
  integer, allocatable, dimension(:), target :: icetsm

  !> \anchor itypsm
  !> type of mass source term for each variable
  !> - 0 for an injection at ambient value,
  !> - 1 for an injection at imposed value.
  integer, allocatable, dimension(:,:), target :: itypsm

  !> \anchor smacel
  !> value of the mass source term for pressure.
  !> For the other variables, eventual imposed injection value.
  !> See the user subroutine \ref cs_user_mass_source_terms
  double precision, allocatable, dimension(:,:), target :: smacel

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

  subroutine init_aux_arrays(ncelet, nfabor)

    use paramx
    use parall
    use period
    use optcal
    use entsor
    use ppincl
    use field
    use cs_c_bindings

    implicit none

    ! Arguments

    integer, intent(in) :: ncelet, nfabor

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

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    ! Local variables

    integer iel
    double precision, allocatable, dimension(:) :: buffer

    procedure() :: synsca

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

  subroutine finalize_aux_arrays

    deallocate(itrifb)
    if (allocated(gamcav)) deallocate(gamcav, dgdpca)

  end subroutine finalize_aux_arrays

  !=============================================================================

  subroutine init_kpdc

    allocate(icepdc(ncepdc))
    allocate(ckupdc(6,ncepdc))
    p_ckupdc = c_loc(ckupdc)

  end subroutine init_kpdc

  !=============================================================================

  subroutine finalize_kpdc

    deallocate(icepdc)
    deallocate(ckupdc)

  end subroutine finalize_kpdc

  !=============================================================================

  subroutine init_tsma(nvar)

    implicit none

    integer :: nvar

    allocate(icetsm(ncetsm))
    allocate(itypsm(ncetsm,nvar))
    allocate(smacel(ncetsm,nvar))

  end subroutine init_tsma

  !=============================================================================

  subroutine finalize_tsma

    ncetsm = 0
    deallocate(icetsm)
    deallocate(itypsm)
    deallocate(smacel)

  end subroutine finalize_tsma

  !=============================================================================

  subroutine boundary_conditions_init

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

  subroutine boundary_conditions_finalize

    use cs_c_bindings

    implicit none

    call cs_f_boundary_conditions_free

  end subroutine boundary_conditions_finalize

  !=============================================================================

  !> \brief Allocate the cs_glob_1d_wall_thermal structure.

  subroutine init_1d_wall_thermal

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_nfpt1d, c_nfpt1t

    call cs_1d_wall_thermal_create

    call cs_f_1d_wall_thermal_get_pointers(c_nfpt1d, c_nfpt1t)

    call c_f_pointer(c_nfpt1d, nfpt1d)
    call c_f_pointer(c_nfpt1t, nfpt1t)

  end subroutine init_1d_wall_thermal

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

  !> \brief Return pointer to the ifpt1d array for the 1D wall thermal module.

  !> \param[out]    ifpt1d         pointer to ifpt1d

  subroutine cs_1d_wall_thermal_get_faces(ifpt1d)

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    integer, dimension(:), pointer, intent(out) :: ifpt1d

    ! Local variables

    type(c_ptr) :: c_ifpt1d

    call cs_f_1d_wall_thermal_get_faces(c_ifpt1d)
    call c_f_pointer(c_ifpt1d, ifpt1d, [nfpt1d])

  end subroutine cs_1d_wall_thermal_get_faces

  !=============================================================================

  !> \brief Return pointer to the tppt1d array for the 1D wall thermal module.

  !> \param[out]    tppt1d         pointer to tppt1d

  subroutine cs_1d_wall_thermal_get_temp(tppt1d)

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    double precision, dimension(:), pointer, intent(out) :: tppt1d

    ! Local variables

    type(c_ptr) :: c_tppt1d

    call cs_f_1d_wall_thermal_get_temp(c_tppt1d)
    call c_f_pointer(c_tppt1d, tppt1d, [nfpt1d])

  end subroutine cs_1d_wall_thermal_get_temp

  !=============================================================================

  !> \brief Return pointers to the mass source term arrays

  !> \param[in]   var_id   id of associated variables
  !> \param[out]  ncesmp   number of cells with mass source terms
  !> \param[out]  icetsm   cell numbers with mass source terms (1 to n)
  !> \param[out]  itpsmp   mass source types (0: ambient value, 1: smacel value)
  !> \param[out]  smacel     mass source values

  subroutine cs_f_volume_mass_injection_get_arrays        &
    (var_id, ncesmp, icetsm_p, itypsm_p, smacel_p)        &
    bind(C, name='cs_f_volume_mass_injection_get_arrays')

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), value :: var_id
    integer(c_int) :: ncesmp
    type(c_ptr) :: icetsm_p, itypsm_p, smacel_p

    ! Local variables

    integer ivar

    ! Get pointer values

    ivar= var_id

    ncesmp = ncetsm
    if (ncetsm.gt.0) then
      icetsm_p = c_loc(icetsm(1))
      itypsm_p = c_loc(itypsm(1, var_id))
      smacel_p = c_loc(smacel(1, var_id))
    else
      icetsm_p = c_null_ptr
      itypsm_p = c_null_ptr
      smacel_p = c_null_ptr
    endif

  end subroutine cs_f_volume_mass_injection_get_arrays

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
