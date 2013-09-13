!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

  !> \defgroup auxiliary Auxiliary variables

  !> \addtogroup auxiliary
  !> \{

  !... Auxiliaires

  !> distance between the center of a given volume and the closest wall,
  !> when it is necessary (\f$R_{ij}-\varepsilon\f$ with wall echo,
  !> LES with van Driest-wall damping, or \f$k-\omega\f$ (SST) turbulence model)
  !> and when \c icdpar=1. The distance between the center of the cell
  !> \c iel and the closest wall is \c dispar(iel)
  double precision, allocatable, dimension(:)   :: dispar

  !> non-dimensional distance \f$y^+\f$ between a given volume and the closest wall,
  !> when it is necessary (LES with van Driest-wall damping) and when \c icdpar=1.
  !> The adimensional distance \f$y^+\f$ between the center of the cell \c iel
  !> and the closest wall is therefore \c yplpar(iel1)
  double precision, allocatable, dimension(:)   :: yplpar

  !> \f$y^+\f$ at boundary, if post-processed
  double precision, allocatable, dimension(:)   :: yplbr

  !> friction velocity at the wall, in the case of a LES calculation
  !> with van Driest-wall damping
  double precision, allocatable, dimension(:)   :: uetbor

  !> stresses at boundary (if post-processed)
  double precision, allocatable, dimension(:,:) :: forbr

  !> \}

  !=============================================================================

  !> \defgroup coupled_case Specific arrays for the coupled case

  !> \addtogroup coupled_case
  !> \{

  !> boundary conditions for the velocity vector with
  !> the coupled velocity components algorithm (\c ivelco=1): see \ref note_2
  double precision, dimension(:,:),   allocatable :: coefau

  !> boundary conditions for the velocity diffusion flux with
  !> the coupled velocity components algorithm (\c ivelco=1): see \ref note_2
  double precision, dimension(:,:),   allocatable :: cofafu

  !> boundary conditions for the velocity convective flux (only for
  !> compressible flows).
  double precision, dimension(:,:),   allocatable :: cofacu

  !> boundary conditions for the velocity vector with
  !> the coupled velocity components algorithm (\c ivelco=1): see \ref note_2
  double precision, dimension(:,:,:), allocatable :: coefbu

  !> boundary conditions for the velocity diffusion flux with
  !> the coupled velocity components algorithm (\c ivelco=1): see \ref note_2
  double precision, dimension(:,:,:), allocatable :: cofbfu

  !> boundary conditions for the velocity convective flux (only for
  !> compressible flows).
  double precision, dimension(:,:,:), allocatable :: cofbcu

  !> explicit Boundary conditions for the mesh velocity.
  !> dim = (3,nfabor)
  double precision, dimension(:,:), allocatable ::   cfaale

  !> explicit Boundary conditions for the mesh velocity.
  !> dim = (3,nfabor)
  double precision, dimension(:,:), allocatable ::   claale

  !> implicit Boundary conditions for the mesh velocity.
  !> dim = (3,3,nfabor)
  double precision, dimension(:,:,:), allocatable :: cfbale

  !> implicit Boundary conditions for the mesh velocity.
  !> dim = (3,3,nfabor)
  double precision, dimension(:,:,:), allocatable :: clbale

  !> boundary condition type at the boundary face \c ifac
  !> (see user subroutine \ref cs\_user\_boundary\_conditions)
  integer, allocatable, dimension(:) :: itypfb

  !> indirection array allowing to sort the boundary faces
  !> according to their boundary condition type \c itypfb
  integer, allocatable, dimension(:) :: itrifb

  !> to identify boundary zones associated wiyth boundary faces
  !> (particular physics)
  integer, allocatable, dimension(:) :: izfppp

  !> to identify boundary zones associated wiyth boundary faces
  !> (radiative transfert)
  integer, allocatable, dimension(:) :: izfrad

  !> number of the wall face (type \c itypfb=iparoi or \c iparug)
  !> which is closest to the center of a given volume when necessary
  !> (\f$R_{ij}-\varepsilon\f$ with wall echo, LES with van Driest-wall damping,
  !> or \f$k-\omega\f$ (SST) turbulence model) and when \c icdpar=2.
  !> The number of the wall face which is the closest to
  !> the center of the cell \c iel is \c ifapat(iel1).
  !> This calculation method is not compatible with parallelism and periodicity
  integer, allocatable, dimension(:) :: ifapat

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
  !> See the user subroutine \ref uskpdc
  integer, save :: ncepdc

  !> number of the \c ncepdc cells in which a pressure drop is imposed.
  !> See \c {iicepd} and the user subroutine \ref uskpdc
  integer, allocatable, dimension(:) :: icepdc

  !> zone with head losses
  integer, allocatable, dimension(:) :: izcpdc

  !> value of the coefficients of the pressure drop tensor of the
  !> \c ncepdc cells in which a pressure drop is imposed.
  !> Note the 6 values are sorted as follows: (k11, k22, k33, k12, k23, k33).
  !> See \c ickpdc and the user subroutine ref uskpdc
  double precision, allocatable, dimension(:,:) :: ckupdc

  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref ustsma
  integer, save :: ncetsm

  !> number of the \c ncetsm cells in which a mass source term is imposed.
  !> See \c iicesm and the user subroutine \ref ustsma}}
  integer, allocatable, dimension(:) :: icetsm

  !> zone where a mass source term is imposed.
  integer, allocatable, dimension(:) :: izctsm

  !> type of mass source term for each variable
  !> - 0 for an injection at ambient value,
  !> - 1 for an injection at imposed value.
  !> See the user subroutine \ref ustsma
  integer, allocatable, dimension(:,:) :: itypsm

  !> value of the mass source term for pressure.
  !> For the other variables, eventual imposed injection value.
  !> See the user subroutine \ref ustsma
  double precision, allocatable, dimension(:,:) :: smacel

  !> value of the porosity
  double precision, allocatable, dimension(:) :: porosi

  !> value of the porosity (for convection and diffusion only)
  double precision, allocatable, dimension(:,:) :: porosf

  !> symmetric tensor cell visco
  double precision, allocatable, dimension(:,:) :: visten

  !> diagonal tensor cell tensor for pressure
  double precision, allocatable, dimension(:,:) :: dttens

  !> \}

  !> \}

contains

  !=============================================================================

  ! Initialize auxiliary arrays

  subroutine init_aux_arrays &

( ncelet , nfabor )

    use paramx
    use numvar, only: ipr, iu
    use parall
    use period
    use optcal
    use entsor
    use ppincl
    use lagran
    use radiat
    use albase
    use ihmpre

    implicit none

    ! Arguments

    integer, intent(in) :: ncelet, nfabor

    ! Local variables
    integer                iok, ivar, iscal, iel

    ! Boundary-face related arrays

    allocate(itrifb(nfabor), itypfb(nfabor))
    if (ippmod(iphpar).ge.1 .or. iihmpr.eq.1) then
      allocate(izfppp(nfabor))
    endif
    if (iirayo.gt.0) then
      allocate(izfrad(nfabor))
    endif

    ! ALE array for structure definition

    if (iale.eq.1) then
      allocate(idfstr(nfabor))
      allocate(cfaale(3,nfabor), claale(3,nfabor))
      allocate(cfbale(3,3,nfabor), clbale(3,3,nfabor))
    endif

    ! Boundary condition for the velocity when components are coupled

    allocate(coefau(3,nfabor),cofafu(3,nfabor))
    allocate(coefbu(3,3,nfabor),cofbfu(3,3,nfabor))
    if (ippmod(icompf).ge.0) then
      allocate(cofacu(3,nfabor))
      allocate(cofbcu(3,3,nfabor))
    endif

    ! Porosity array when needed

    if (iporos.ge.1) then
      allocate(porosi(ncelet))
      do iel = 1, ncelet
        porosi(iel) = 1.d0
      enddo
    endif
    if (iporos.eq.2) then
      allocate(porosf(6, ncelet))
      do iel = 1, ncelet
        porosf(1, iel) = 1.d0
        porosf(2, iel) = 1.d0
        porosf(3, iel) = 1.d0
        porosf(4, iel) = 0.d0
        porosf(5, iel) = 0.d0
        porosf(6, iel) = 0.d0
      enddo
    endif

    ! Symmetric cell diffusivity when needed
    iok = 0
    do ivar = 1, nvarmx
      if (idften(ivar).eq.6) iok = 1
    enddo

    do iscal = 1, nscamx
      if (ityturt(iscal).eq.3) iok = 1
    enddo

    ! Also tensorial diffusion for the velocity in case of tensorial porosity
    if (iporos.eq.2) then
      idften(iu) = 6
      iok = 1
    endif

    if (iok.eq.1) then
      allocate(visten(6,ncelet))
    endif

    ! Diagonal cell tensor for the pressure solving when needed
    if (ncpdct.gt.0.or.ipucou.eq.1.or.iporos.eq.2) then
      idften(ipr) = 6
      allocate(dttens(6,ncelet))
    endif

    ! Wall-distance calculation

    if (ineedy.eq.1 .and. abs(icdpar).eq.1) then
      allocate(dispar(ncelet))
      if (     (itytur.eq.4 .and. idries.eq.1) &
          .or. (iilagr.ge.1 .and. iroule.eq.2) ) then
        allocate(yplpar(ncelet))
      endif
    endif
    if (ineedy.eq.1 .and. abs(icdpar).eq.2) then
      allocate(ifapat(ncelet))
    endif

    ! Forces on boundary faces

    if (ineedf.eq.1) then
      allocate(forbr(3,nfabor))
    endif

    ! Friction velocity on boundary faces

    if (     (itytur.eq.4 .and. idries.eq.1) &
        .or. (iilagr.ge.1 .and. idepst.gt.0) ) then
      allocate(uetbor(nfabor))
    endif

    ! Non-dimensional distance to the wall (for post-processing)

    if (ipstdv(ipstyp).ne.0) then
      allocate(yplbr(nfabor))
    endif

    ! Temporary storage arrays for k-omega model

    if (iturb.eq.60) then
      allocate(s2kw(ncelet))
      allocate(divukw(ncelet))
    endif

   ! Strain rate tensor at the previous time step
   ! if rotation curvature correction
   if (irccor.eq.1) then
      if (idtvar.ge.0) then
        allocate(straio(ncelet,6))
      endif
    endif

    return

  end subroutine init_aux_arrays

  !=============================================================================

  ! Free auxiliary arrays

  subroutine finalize_aux_arrays

    deallocate(itrifb, itypfb)
    if (allocated(izfppp)) deallocate(izfppp)
    if (allocated(izfrad)) deallocate(izfrad)
    if (allocated(idfstr)) deallocate(idfstr)
    if (allocated(izcpdc)) deallocate(izcpdc)
    if (allocated(izctsm)) deallocate(izctsm)
    if (allocated(izft1d)) deallocate(izft1d)
    if (allocated(coefau)) deallocate(coefau, cofafu,                   &
                                      coefbu, cofbfu)
    if (allocated(cofacu)) deallocate(cofacu, cofbcu)
    if (allocated(porosi)) deallocate(porosi)
    if (allocated(visten)) deallocate(visten)
    if (allocated(dttens)) deallocate(dttens)
    if (allocated(cfaale)) deallocate(cfaale, cfbale, claale, clbale)
    if (allocated(dispar)) deallocate(dispar)
    if (allocated(yplpar)) deallocate(yplpar)
    if (allocated(ifapat)) deallocate(ifapat)
    if (allocated(forbr)) deallocate(forbr)
    if (allocated(uetbor)) deallocate(uetbor)
    if (allocated(yplbr)) deallocate(yplbr)
    if (allocated(s2kw)) deallocate(s2kw, divukw)
    if (allocated(straio))  deallocate(straio)

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

end module pointe
