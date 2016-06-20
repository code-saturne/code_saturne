!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

  use paramx

  implicit none

  !=============================================================================
  !> \defgroup pointer_variables Module for pointer variables

  !> \addtogroup pointer_variables
  !> \{

  !> \defgroup fortran_pointer_containers Containers to Fortran array pointers.
  !> An array of one of these these derived types can be used to manage a set of
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

  !> distance between the center of a given volume and the closest wall,
  !> when it is necessary (\f$R_{ij}-\varepsilon\f$ with wall echo,
  !> LES with van Driest-wall damping, or \f$k-\omega\f$ (SST) turbulence model).
  !> The distance between the center of the cell
  !> \c iel and the closest wall is \c dispar(iel)
  double precision, allocatable, dimension(:)   :: dispar

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
  integer, dimension(:), pointer :: itypfb(:)

  !> indirection array allowing to sort the boundary faces
  !> according to their boundary condition type \c itypfb
  integer, allocatable, dimension(:) :: itrifb

  !> to identify boundary zones associated with boundary faces
  !> (particular physics)
  integer, allocatable, dimension(:) :: izfppp

  !> to identify boundary zones associated with boundary faces
  !> (radiative transfer)
  integer, allocatable, dimension(:) :: izfrad

  !> the index of the structure, (\c idfstr(ifac) where \c ifac is the index
  !> of the face), 0 if the face is not coupled to any structure.
  integer, allocatable, dimension(:) :: idfstr

  !> square of the norm of the deviatoric part of the deformation
  !> rate tensor (\f$S^2=2S_{ij}^D S_{ij}^D\f$). This array is defined only
  !> with the \f$k-\omega\f$ (SST) turbulence model
  double precision, allocatable, dimension(:)   :: s2kw

  !> divergence of the velocity. More precisely it is the trace of the velocity
  !> gradient (and not a finite volume divergence term). In the
  !> cell \c iel,  \f$div(\vect{u})\f$ is given by \c divukw(iel1).
  !> This array is defined only with the \f$k-\omega\f$ SST turbulence model
  !> (because in this case it may be calculated at the same time as \f$S^2\f$)
  double precision, allocatable, dimension(:)   :: divukw

  !> strain rate tensor at the previous time step
  double precision, allocatable, dimension(:,:) :: straio

  !> \}

  !=============================================================================

  !... Parametres du module thermique 1D
  !> \defgroup thermal_1D Thermal 1D module parameters

  !> \addtogroup thermal_1D
  !> \{


  !> number of boundary faces which are coupled
  !> with a wall 1D thermal module. See the user subroutine \ref uspt1d
  integer, save :: nfpt1d

  ! TODO
  integer, save :: nmxt1d

  !> zones of t1d, dimensioned with nfabor (TODO)
  integer, allocatable, dimension(:) :: izft1d

  !> number of discretisation cells in the 1D wall for the
  !> \c nfpt1d boundary faces which are coupled with a wall 1D thermal module.
  !> The number of cells for these boundary faces is given by
  !> \c nppt1d(ii), with 1 <= ii <= nfpt1d.
  !> See the user subroutine \ref uspt1d
  integer, allocatable, dimension(:) :: nppt1d

  !> array allowing to mark out the numbers of
  !> the \c nfpt1d boundary faces which are coupled with a wall 1D
  !> thermal module. The numbers of these boundary faces are
  !> given by \c ifpt1d(ii), with 1 <= ii <= \c nfpt1d.
  !> See the user subroutine \ref uspt1d
  integer, allocatable, dimension(:) :: ifpt1d

  !> typical boundary condition at the external (pseudo) wall:
  !> Dirichlet condition (\c iclt1d=1) or flux condition (\c iclt1d=3)
  integer, allocatable, dimension(:) :: iclt1d

  !> thickness of the 1D wall for the \c nfpt1d boundary faces
  !> which are coupled with a wall 1D thermal module.
  !> The wall thickness for these boundary faces is therefore given by
  !> \c eppt1d(ii), with 1 <= ii <= \c nfpt1d.
  !> See the user subroutine \ref uspt1d
  double precision, allocatable, dimension(:) :: eppt1d

  !> geometry of the pseudo wall mesh (refined as a fluid
  !> if \c rgt1d is smaller than 1
  double precision, allocatable, dimension(:) :: rgpt1d

  !> initialisation temperature of the wall (uniform in thickness).
  !> In the course of the calculation, the array stores the temperature
  !> of the solid at the fluid/solid interface.
  double precision, allocatable, dimension(:) :: tppt1d

  !> external temperature of the pseudo wall in the Dirichlet case.
  double precision, allocatable, dimension(:) :: tept1d

  !> external coefficient of transfer in the pseudo wall under Dirichlet conditions
  !> (in \f$W.m^{-2}.K^.\f$).
  double precision, allocatable, dimension(:) :: hept1d

  ! fept1d ! nfpt1d                  ! flux thermique exterieur
  !> external heat flux in the pseudo wall under the flux conditions
  !> (in \f$W.m^{-2}\f$, negative value for energy entering the wall).
  double precision, allocatable, dimension(:) :: fept1d

  !> thermal diffusivity
  double precision, allocatable, dimension(:) :: xlmbt1

  !> volumetric heat capacity \f$\rho C_p\f$ of the wall uniform in thickness
  !> (in \f$J.m^{-3}.K^{-1}\f$).
  double precision, allocatable, dimension(:) :: rcpt1d

  !> physical time step associated with the solved 1D equation of the pseudo wall
  !> (which can be different from the time step in the calculation).
  double precision, allocatable, dimension(:) :: dtpt1d

  !> \}

  !=============================================================================

  !... Auxiliaires
  !> \addtogroup auxiliary
  !> \{


  !> number of cells in which a pressure drop is imposed.
  !> See the user subroutine \ref cs_user_head_losses
  integer, save :: ncepdc

  !> number of the \c ncepdc cells in which a pressure drop is imposed.
  !> See \c {iicepd} and the user subroutine \ref cs_user_head_losses
  integer, allocatable, dimension(:) :: icepdc

  !> zone with head losses
  integer, allocatable, dimension(:) :: izcpdc

  !> value of the coefficients of the pressure drop tensor of the
  !> \c ncepdc cells in which a pressure drop is imposed.
  !> Note the 6 values are sorted as follows: (k11, k22, k33, k12, k23, k33).
  !> See \c ickpdc and the user subroutine \ref cs_user_head_losses
  double precision, allocatable, dimension(:,:) :: ckupdc

  !> Head loss factor of the fluid outside the domain, between infinity and
  !> the entrance (for \ref ifrent boundary type). The default value is 0,
  !> dimensionless factor. The user may give a value in
  !> \ref cs_user_boundary_conditions in the array
  !> \c rcodcl(ifac, \ref ipr, 2).
  double precision, allocatable, dimension(:) :: b_head_loss

  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref cs_user_mass_source_terms
  integer, save :: ncetsm

  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref cs_user_mass_source_terms}}
  integer, allocatable, dimension(:) :: icetsm

  !> zone where a mass source term is imposed.
  integer, allocatable, dimension(:) :: izctsm

  !> type of mass source term for each variable
  !> - 0 for an injection at ambient value,
  !> - 1 for an injection at imposed value.
  !> See the user subroutine \ref cs_user_mass_source_terms
  integer, allocatable, dimension(:,:) :: itypsm

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
    use radiat, only: iirayo
    use albase
    use ihmpre
    use field

    implicit none

    ! Arguments

    integer, intent(in) :: ncelet, nfabor

    ! Local variables
    integer                ifac
    integer                kdiftn

    ! Key id for diffusivity tensor
    call field_get_key_id("diffusivity_tensor", kdiftn)

    ! Boundary-face related arrays

    allocate(itrifb(nfabor))
    if (ippmod(iphpar).ge.1 .or. iihmpr.eq.1) then
      allocate(izfppp(nfabor))
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif
    if (iirayo.gt.0) then
      allocate(izfrad(nfabor))
      do ifac = 1, nfabor
        izfrad(ifac) = 0
      enddo
    endif

    ! ALE array for structure definition

    if (iale.eq.1) then
      allocate(idfstr(nfabor))
    endif

    ! Also tensorial diffusion for the velocity in case of tensorial porosity
    if (iporos.eq.2) then
      idften(iu) = 6
      ! Key word: tensorial diffusivity
      call field_set_key_int(ivarfl(iu), kdiftn, idften(iu))
    endif

    ! Diagonal cell tensor for the pressure solving when needed
    if (ncpdct.gt.0.or.ipucou.eq.1.or.iporos.eq.2) then
      idften(ipr) = 6
    endif

    ! Wall-distance calculation

    if (ineedy.eq.1) then
      allocate(dispar(ncelet))
      if (itytur.eq.4 .and. idries.eq.1) then
        allocate(yplpar(ncelet))
      endif
    endif

    ! Temporary storage arrays for k-omega model

    if (iturb.eq.60) then
      allocate(s2kw(ncelet))
      allocate(divukw(ncelet))
    endif

    ! Strain rate tensor at the previous time step
    ! if rotation curvature correction of eddy viscosity
    if (irccor.eq.1) then
      if (idtvar.ge.0) then
        allocate(straio(ncelet,6))
      endif
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

    integer iel, isou
    double precision, allocatable, dimension(:) :: buffer
    double precision, allocatable, dimension(:,:) :: buff2

    ! Resize/copy arrays

    allocate(buffer(ncelet))

    ! Wall-distance calculation

    if (allocated(dispar)) then
      do iel = 1, ncel
        buffer(iel) = dispar(iel)
      enddo
      deallocate(dispar)
      call synsca (buffer)
      allocate(dispar(ncelet))
      do iel = 1, ncelet
        dispar(iel) = buffer(iel)
      enddo
    endif

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

    ! Temporary storage arrays for k-omega model

    if (allocated(s2kw)) then
      do iel = 1, ncel
        buffer(iel) = s2kw(iel)
      enddo
      deallocate(s2kw)
      call synsca (buffer)
      allocate(s2kw(ncelet))
      do iel = 1, ncelet
        s2kw(iel) = buffer(iel)
      enddo

      do iel = 1, ncel
        buffer(iel) = divukw(iel)
      enddo
      deallocate(divukw)
      call synsca (buffer)
      allocate(divukw(ncelet))
      do iel = 1, ncelet
        divukw(iel) = buffer(iel)
      enddo
    endif

    ! Strain rate tensor at the previous time step
    ! for rotation-curvature correction of eddy viscosity

    if (allocated(straio)) then
      allocate(buff2(ncel,6))
      do isou = 1, 6
        do iel = 1, ncel
          buff2(iel,isou) = straio(iel,isou)
        enddo
      enddo
      deallocate(straio)
      allocate(straio(ncelet,6))
      do isou = 1, 6
        do iel = 1, ncel
          straio(iel,isou) = buff2(iel,isou)
        enddo
      enddo
      deallocate(buff2)
      call synten (straio(1,1), straio(1,4), straio(1,5), &
                   straio(1,4), straio(1,2), straio(1,6), &
                   straio(1,5), straio(1,6), straio(1,3))
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
    if (allocated(izfppp)) deallocate(izfppp)
    if (allocated(izfrad)) deallocate(izfrad)
    if (allocated(idfstr)) deallocate(idfstr)
    if (allocated(izcpdc)) deallocate(izcpdc)
    if (allocated(izctsm)) deallocate(izctsm)
    if (allocated(izft1d)) deallocate(izft1d)
    if (allocated(dispar)) deallocate(dispar)
    if (allocated(yplpar)) deallocate(yplpar)
    if (allocated(s2kw)) deallocate(s2kw, divukw)
    if (allocated(straio))  deallocate(straio)
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

  subroutine init_pt1d

    use cstnum

    allocate(nppt1d(nfpt1d), ifpt1d(nfpt1d), iclt1d(nfpt1d))
    allocate(eppt1d(nfpt1d), rgpt1d(nfpt1d), tppt1d(nfpt1d))
    allocate(tept1d(nfpt1d), hept1d(nfpt1d), fept1d(nfpt1d))
    allocate(xlmbt1(nfpt1d), rcpt1d(nfpt1d), dtpt1d(nfpt1d))

    !---> INITIALISATION DES TABLEAUX
    !     a des valeurs sortant en erreur dans vert1d
    !     sauf pour les variables de conditions aux limites
    !     qui sont initialisees a des valeurs standard
    !     (flux nul, Timpose=0, coef d'echange infini)

    ifpt1d(:) = -999
    nppt1d(:) = -999
    iclt1d(:) = 3
    eppt1d(:) = -999.d0
    rgpt1d(:) = -999.d0
    tppt1d(:) = 0.d0
    tept1d(:) = 0.d0
    hept1d(:) = rinfin
    fept1d(:) = 0.d0
    xlmbt1(:) = -999.d0
    rcpt1d(:) = -999.d0
    dtpt1d(:) = -999.d0

  end subroutine init_pt1d

  !=============================================================================

  subroutine finalize_pt1d

    deallocate(nppt1d, ifpt1d, iclt1d)
    deallocate(eppt1d, rgpt1d, tppt1d)
    deallocate(tept1d, hept1d, fept1d)
    deallocate(xlmbt1, rcpt1d, dtpt1d)

  end subroutine finalize_pt1d

  !=============================================================================

  subroutine boundary_conditions_init

    use, intrinsic :: iso_c_binding
    use mesh
    use cs_c_bindings

    implicit none

    ! Local variables

    type(c_ptr) :: c_itypfb

    call cs_f_boundary_conditions_type_create

    call cs_f_boundary_conditions_type_get_pointer(c_itypfb)

    call c_f_pointer(c_itypfb, itypfb, [nfabor])

  end subroutine boundary_conditions_init

  !=============================================================================

  subroutine boundary_conditions_finalize

    use cs_c_bindings

    implicit none

    call cs_f_boundary_conditions_type_free

  end subroutine boundary_conditions_finalize

end module pointe
