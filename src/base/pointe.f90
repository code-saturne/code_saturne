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
  ! handle different pointer shapes with a single type, but this feature
  ! of Fortran 2003 (extended in 2008) is supported only by very recent
  ! compilers (2010 or later), and even recent compilers such as Intel 12
  ! reportedly have bugs using this feature, so we will need to avoid
  ! depending on it for a few years.

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

  !> container for rank 3 double precision array pointer.
  type pmapper_double_r3
    double precision, dimension(:,:,:),  pointer :: p !< rank 3 array pointer
  end type pmapper_double_r3

  !> \}

  !=============================================================================

  !> \defgroup dummy_arrays Dummy target arrays for null pointers

  !> \addtogroup dummy_arrays
  !> \{

  integer, dimension(1),   target :: ivoid1
  integer, dimension(1,1), target :: ivoid2

  double precision, dimension(1),     target :: rvoid1
  double precision, dimension(1,1),   target :: rvoid2
  double precision, dimension(1,1,1), target :: rvoid3


  !> \}

  !=============================================================================

  !> \defgroup auxiliary Auxiliary variables

  !> \addtogroup auxiliary
  !> \{

  !... Auxiliaires

  !> non-dimensional distance \f$y^+\f$ between a given volume and the closest
  !> wall, when it is necessary (LES with van Driest-wall damping).
  !> The adimensional distance \f$y^+\f$ between the center of the cell \c iel
  !> and the closest wall is therefore \c yplpar(iel1)
  double precision, allocatable, dimension(:)   :: yplpar

  !> \}
  !=============================================================================

  !> \defgroup coupled_case Specific arrays for the coupled case

  !> \addtogroup coupled_case
  !> \{

  !> \anchor itypfb
  !> boundary condition type at the boundary face \c ifac
  !> (see user subroutine \ref cs\_user\_boundary\_conditions)
  integer, dimension(:), pointer, save :: itypfb

  !> indirection array allowing to sort the boundary faces
  !> according to their boundary condition type \c itypfb
  integer, allocatable, dimension(:) :: itrifb

  !> to identify boundary zones associated with boundary faces
  !> (specific physics models)
  integer, dimension(:), pointer :: izfppp

  !> the index of the structure, (\c idfstr(ifac) where \c ifac is the index
  !> of the face), 0 if the face is not coupled to any structure.
  integer, allocatable, dimension(:) :: idfstr

  !> \}

  !=============================================================================

  !... Parametres du module thermique 1D
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

  !... Auxiliaires
  !> \addtogroup auxiliary
  !> \{

  !> \anchor ncepdc
  !> number of cells in which a pressure drop is imposed.
  !> See the user subroutine \ref cs_f_user_head_losses
  integer, save :: ncepdc

  !> \anchor icepdc
  !> number of the \c ncepdc cells in which a pressure drop is imposed.
  !> See \c {iicepd} and the user subroutine \ref cs_f_user_head_losses
  integer, allocatable, dimension(:) :: icepdc

  !> zone with head losses
  integer, allocatable, dimension(:) :: izcpdc

  !> \anchor ckupdc
  !> value of the coefficients of the pressure drop tensor of the
  !> \c ncepdc cells in which a pressure drop is imposed.
  !> Note the 6 values are sorted as follows: (k11, k22, k33, k12, k23, k33).
  !> See \c ickpdc and the user subroutine \ref cs_f_user_head_losses
  double precision, allocatable, dimension(:,:) :: ckupdc

  !> Head loss factor of the fluid outside the domain, between infinity and
  !> the entrance (for \ref ifrent boundary type). The default value is 0,
  !> dimensionless factor. The user may give a value in
  !> \ref cs_user_boundary_conditions in the array
  !> \c rcodcl(ifac, \ref ipr, 2).
  double precision, allocatable, dimension(:) :: b_head_loss

  !> \anchor ncetsm
  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref cs_user_mass_source_terms
  integer, save :: ncetsm

  !> \anchor icetsm
  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref cs_user_mass_source_terms}}
  integer, allocatable, dimension(:) :: icetsm

  !> zone where a mass source term is imposed.
  integer, allocatable, dimension(:) :: izctsm

  !> \anchor itypsm
  !> type of mass source term for each variable
  !> - 0 for an injection at ambient value,
  !> - 1 for an injection at imposed value.
  !> See the user subroutine \ref cs_user_mass_source_terms
  integer, allocatable, dimension(:,:) :: itypsm

  !> \anchor smacel
  !> value of the mass source term for pressure.
  !> For the other variables, eventual imposed injection value.
  !> See the user subroutine \ref cs_user_mass_source_terms
  double precision, allocatable, dimension(:,:) :: smacel

  !> liquid-vapour mass transfer term for cavitating flows
  !> and its derivative with respect to pressure
  double precision, allocatable, dimension(:) :: gamcav, dgdpca

  !> number of the nfbpcd faces in which a condensation source terms is imposed.
  !> See \c ifbpcd and the user subroutine \ref cs_user_boundary_mass_source_terms
  integer, save :: nfbpcd

  !> list on the nfbpcd faces in which a condensation source terms is imposed.
  !> See \c ifbpcd and the user subroutine \ref cs_user_boundary_mass_source_terms
  integer, allocatable, dimension(:) :: ifbpcd

  !> zone where a condensation source terms is imposed.
  integer, allocatable, dimension(:) :: izftcd

  !> type of condensation source terms for each variable
  !> - 0 for an variable at ambient value,
  !> - 1 for an variable at imposed value.
  !> See the user subroutine \ref cs_user_boundary_mass_source_terms
  integer, allocatable, dimension(:,:) :: itypcd

  !> value of the condensation source terms for pressure.
  !> For the other variables, eventual imposed specific value.
  !> See the user subroutine \ref cs_user_boundary_mass_source_terms
  double precision, allocatable, dimension(:,:) :: spcond

  !> value of the thermal flux for the condensation model.
  !> See the user subroutine \ref cs_user_boundary_mass_source_terms
  double precision, allocatable, dimension(:) :: thermal_condensation_flux

  !> value of the thermal exchange coefficient associated to
  !> the condensation model used.
  !> See the user subroutine \ref cs_user_boundary_mass_source_terms
  double precision, allocatable, dimension(:) :: hpcond

  !> Specific 1D thermal model with implicit time scheme (only used
  !> with condensation modelling to the cold wall)
  !> flthr     ! external heat flux used as flux conditions
  !>           ! of the 1d thermal model (in unit \f$W.m^{-2}\f$).
  double precision, allocatable, dimension(:) :: flthr
  !> dflthr    ! external heat flux derivative used as flux conditions
  !>           ! of the 1d thermal model (in unit \f$W.m^{-3}\f$).
  double precision, allocatable, dimension(:) :: dflthr

  !> number of the ncmast cells in which a condensation source terms is imposed.
  !> See \c lstmast list and the subroutine \ref cs_user_metal_structures_source_terms
  integer, save :: ncmast

  !> list on the ncmast cells in which a condensation source terms is imposed.
  !> See  the user subroutine \ref cs_user_metal_structures_source_terms.
  integer, allocatable, dimension(:) :: ltmast

  !> zone type where a condensation source terms is imposed to model
  !> the metal structures condensation on a volumic zone.
  integer, allocatable, dimension(:) :: izmast

  !> type of condensation source terms for each variable
  !> - 0 for a variable at ambient value,
  !> - 1 for a variable at imposed value.
  !> See the user subroutine \ref  cs_user_metal_structures_source_terms.
  integer, allocatable, dimension(:,:) :: itypst

  !> value of the condensation source terms for pressure
  !> associated to the metal structures modelling.
  !> For the other variables, eventual imposed specific value.
  !> See the user subroutine \ref cs_user_metal_structures_source_terms.
  double precision, allocatable, dimension(:,:) :: svcond

  !> value of the thermal flux for the condensation model
  !> associated to the metal structures modelling.
  !> See the user subroutine \ref cs_user_metal_structures_source_terms.
  double precision, allocatable, dimension(:) :: flxmst


  !> \}

  !=============================================================================

  !> \addtogroup lagran
  !> \{
  !> \defgroup lag_arrays Lagrangian arrays

  !> \addtogroup lag_arrays
  !> \{

  !> \anchor tslagr
  double precision, pointer, dimension(:,:), save :: tslagr
  !> \}

  !> \}

  !> \}

contains

  !=============================================================================

  ! Initialize auxiliary arrays

  subroutine init_aux_arrays &

( ncelet , nfabor )

    use paramx
    use numvar, only: ipr, iu, ivarfl
    use parall
    use period
    use optcal
    use entsor
    use ppincl
    use lagran
    use albase
    use ihmpre
    use field
    use cs_c_bindings

    implicit none

    ! Arguments

    integer, intent(in) :: ncelet, nfabor

    ! Local variables

    type(var_cal_opt) :: vcopt

    ! Boundary-face related arrays

    allocate(itrifb(nfabor))

    ! ALE array for structure definition

    if (iale.eq.1) then
      allocate(idfstr(nfabor))
    endif

    ! Also tensorial diffusion for the velocity in case of tensorial porosity
    if (iporos.eq.2) then
      ! Tensorial diffusivity
      call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
      vcopt%idften = ANISOTROPIC_LEFT_DIFFUSION
      call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)
    endif

    ! Diagonal cell tensor for the pressure solving when needed
    if (ncpdct.gt.0.or.ipucou.eq.1.or.iporos.eq.2) then
      call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
      vcopt%idften = ANISOTROPIC_LEFT_DIFFUSION
      call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
    endif

    ! Wall-distance calculation

    if (itytur.eq.4 .and. idries.eq.1) then
      allocate(yplpar(ncelet))
    endif

    ! liquid-vapour mass transfer term for cavitating flows
    ! and its part implicit in pressure
    if (icavit.ge.0) then
      allocate(gamcav(ncelet), dgdpca(ncelet))
    endif

    return

  end subroutine init_aux_arrays

  !=============================================================================

  ! Resize auxiliary arrays

  subroutine resize_aux_arrays

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    ! Local variables

    integer iel
    double precision, allocatable, dimension(:) :: buffer

    ! Resize/copy arrays

    allocate(buffer(ncelet))

    if (allocated(yplpar)) then
      do iel = 1, ncel
        buffer(iel) = yplpar(iel)
      enddo
      deallocate(yplpar)
      call synsca (buffer)
      allocate(yplpar(ncelet))
      do iel = 1, ncelet
        yplpar(iel) = buffer(iel)
      enddo
    endif

    ! liquid-vapour mass transfer term for cavitating flows
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

    return

  end subroutine resize_aux_arrays

  !=============================================================================

  ! Free auxiliary arrays

  subroutine finalize_aux_arrays

    deallocate(itrifb)
    if (allocated(idfstr)) deallocate(idfstr)
    if (allocated(izcpdc)) deallocate(izcpdc)
    if (allocated(izctsm)) deallocate(izctsm)
    if (allocated(yplpar)) deallocate(yplpar)
    if (allocated(b_head_loss)) deallocate(b_head_loss)
    if (allocated(gamcav)) deallocate(gamcav, dgdpca)

    return

  end subroutine finalize_aux_arrays

  !=============================================================================

  subroutine init_kpdc

    allocate(icepdc(ncepdc))
    allocate(ckupdc(ncepdc,6))

  end subroutine init_kpdc

  !=============================================================================

  subroutine finalize_kpdc

    deallocate(icepdc)
    deallocate(ckupdc)

  end subroutine finalize_kpdc

  !=============================================================================

  subroutine init_tsma ( nvar )

    implicit none

    integer :: nvar

    allocate(icetsm(ncetsm))
    allocate(itypsm(ncetsm,nvar))
    allocate(smacel(ncetsm,nvar))

  end subroutine init_tsma

  !=============================================================================

  subroutine finalize_tsma

    deallocate(icetsm)
    deallocate(itypsm)
    deallocate(smacel)

  end subroutine finalize_tsma

  !=============================================================================
  !=============================================================================

  subroutine init_pcond ( nvar )

    implicit none

    integer :: nvar

    allocate(ifbpcd(nfbpcd))
    allocate(itypcd(nfbpcd,nvar))
    allocate(spcond(nfbpcd,nvar))
    allocate(thermal_condensation_flux(nfbpcd))
    allocate(hpcond(nfbpcd))
    allocate(flthr(nfbpcd),dflthr(nfbpcd))

    !---> Array initialization
    flthr(:)  = 0.d0
    dflthr(:) = 0.d0

  end subroutine init_pcond

  !=============================================================================

  subroutine finalize_pcond
    deallocate(ifbpcd)
    deallocate(itypcd)
    deallocate(spcond)
    deallocate(thermal_condensation_flux)
    deallocate(hpcond)
    deallocate(flthr, dflthr)

  end subroutine finalize_pcond
  !=============================================================================

  subroutine init_vcond ( nvar , ncelet)

    implicit none

    integer :: nvar, ncelet

    allocate(ltmast(ncelet))
    allocate(izmast(ncelet))
    allocate(itypst(ncelet, nvar))
    allocate(svcond(ncelet, nvar))
    allocate(flxmst(ncelet))

  end subroutine init_vcond

  !=============================================================================

  subroutine finalize_vcond
    deallocate(ltmast)
    deallocate(itypst)
    deallocate(izmast)
    deallocate(svcond)
    deallocate(flxmst)

  end subroutine finalize_vcond

  !=============================================================================

  subroutine boundary_conditions_init

    use, intrinsic :: iso_c_binding
    use mesh
    use cs_c_bindings

    implicit none

    ! Local variables

    type(c_ptr) :: c_itypfb, c_izfppp

    call cs_f_boundary_conditions_create

    call cs_f_boundary_conditions_get_pointers(c_itypfb, c_izfppp)

    call c_f_pointer(c_itypfb, itypfb, [nfabor])
    call c_f_pointer(c_izfppp, izfppp, [nfabor])

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

    return

  end subroutine init_1d_wall_thermal

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

end module pointe
