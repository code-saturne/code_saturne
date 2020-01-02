!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file numvar.f90
!> \brief Module for variable numbering

module numvar

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup numvar Module for variable numbering

  !> \addtogroup numvar
  !> \{

  !----------------------------------------------------------------------------
  ! Main variables
  !----------------------------------------------------------------------------

  !> \defgroup main_variables Main variables
  !> \brief Main variable field indices (previously stored in rtp, rtpa).

  !> \addtogroup main_variables
  !> \{

  !> \anchor ipr
  !> pressure
  integer, save :: ipr

  !> \anchor iu
  !> velocity component \f$ u_x \f$
  integer, save :: iu

  !> \anchor iv
  !> velocity component \f$ u_y \f$
  integer, save :: iv

  !> \anchor iw
  !> velocity component \f$ u_z \f$
  integer, save :: iw

  !> \anchor ivolf2
  !> void fraction for VOF method
  integer, save :: ivolf2

  !> \anchor ik
  !> turbulent kinetic energy \f$ k \f$
  integer, save :: ik

  !> \anchor iep
  !> turbulent dissipation \f$ \varepsilon \f$
  integer, save :: iep

  !> \anchor ir11
  !> Reynolds stress component \f$ R_{xx} \f$
  integer, save :: ir11

  !> \anchor ir22
  !> Reynolds stress component \f$ R_{yy} \f$
  integer, save :: ir22

  !> \anchor ir33
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir33

  !> \anchor ir12
  !> Reynolds stress component \f$ R_{xy} \f$
  integer, save :: ir12

  !> \anchor ir23
  !> Reynolds stress component \f$ R_{yz} \f$
  integer, save :: ir23

  !> \anchor ir13
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir13

  !> \anchor irij
  !> Reynolds stress tenso \f$ R_{ij} \f$
  integer, save :: irij

  !> \anchor iphi
  !> variable \f$ \phi \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: iphi

  !> \anchor ifb
  !> variable \f$ f_b \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: ifb

  !> \anchor ial
  !> variable \f$ \alpha \f$ of the \f$ Bl-v^2-k \f$ model
  integer, save :: ial

  !> \anchor iomg
  !> variable \f$ \omega \f$ of the \f$ k-\omega \f$ SST
  integer, save :: iomg

  !> \anchor inusa
  !> variable \f$ \widetilde{\nu}_T \f$ of the Spalart Allmaras
  integer, save :: inusa

  !> \anchor isca
  !> isca(i) is the index of the scalar i
  integer, save :: isca(nscamx)

  !> \anchor iscapp
  !> iscapp(i) is the index of the specific physics scalar i
  integer, save :: iscapp(nscamx)

  !> \anchor nscaus
  !> number of user scalars solutions of an advection equation
  integer, save :: nscaus

  !> \anchor nscapp
  !> number of specific physics scalars
  integer, save :: nscapp

  !> \anchor nscasp
  !> number of species scalars
  integer(c_int), pointer, save :: nscasp

  !> \anchor iuma
  !> mesh velocity component \f$ w_x \f$
  integer, save :: iuma

  !> \anchor ivma
  !> mesh velocity component \f$ w_y \f$
  integer, save :: ivma

  !> \anchor iwma
  !> mesh velocity component \f$ w_z \f$
  integer, save :: iwma

  !> \}

  !----------------------------------------------------------------------------
  ! Physical properties
  !----------------------------------------------------------------------------

  !> \defgroup physical_prop Physical properties
  !> \brief Physical properties are using the field API.
  !> See \ref cs_user_boundary_conditions for some examples.

  !> \addtogroup physical_prop
  !> \{

  !> Density at the current time step (equal to icrom, kept for compatibility)
  integer, save :: irom

  !> dynamic molecular viscosity (in kg/(m.s))
  integer, save :: iviscl

  !> dynamic turbulent viscosity
  integer, save :: ivisct

  !> error estimator for Navier-Stokes
  integer, save :: iestim(nestmx)

  !> interior and boundary convective mass flux key ids of the variables
  integer, save :: kimasf, kbmasf

  !> constant diffusivity field id key for scalars
  integer, save :: kvisl0

  !> variable diffusivity field id key for scalars
  integer, save :: kivisl

  !> variable density field id key for scalars
  integer, save :: kromsl

  !> source terms at previous time step for 2nd order
  integer, save :: kstprv

  !> source terms at the current time step (used for limiters)
  integer, save :: kst

  !> turbulent schmidt key for scalars
  integer, save :: ksigmas

  !> convective mass flux of the variables at the previous time-step
  integer, save :: ifluaa(nvarmx)

  !> cell density field ids of the variables
  integer, save :: icrom

  !> boundary density field ids of the variables
  integer, save :: ibrom

  !> field ids of the cell porosity
  integer, save :: ipori, iporf

  !> dynamic constant of Smagorinsky
  integer, save :: ismago

  !> field ids of the anisotropic viscosity
  !> \remark turbulent or Darcy module anisotropic diffusion
  integer, save :: ivsten, ivstes

  !> Courant number
  integer, save :: icour

  !> Fourier number
  integer, save :: ifour

  !> Total pressure at cell centers
  !> \f$ P_{tot} = P^\star +\rho \vect{g} \cdot (\vect{x}-\vect{x}_0) \f$
  integer, save :: iprtot

  !> Mesh velocity viscosity for the ALE module
  !> \remark might be orthotropic
  integer, save :: ivisma

  !> pointer for dilatation source terms
  integer, save :: iustdy(nscamx)

  !> pointer for global dilatation source terms
  integer, save :: itsrho

  !> pointer for thermal expansion coefficient
  integer, save :: ibeta

  !> pointer for deduced mass fraction in case of gas mix
  integer, save :: iddgas

  !> pointer for gas mix molar mass
  integer, save :: igmxml

  !> field id of the stresses at boundary  (if post-processed)
  integer, save :: iforbr

  !>  field id of \f$y^+\f$ at boundary (if post-processed)
  integer, save :: iyplbr

  !>  field id of temperature at boundary
  integer, save ::  itempb

  !> field id of the square of the norm of the deviatoric
  !> part of the deformation rate tensor (\f$S^2=2S_{ij}^D S_{ij}^D\f$).
  !> Field defined only with the \f$k-\omega\f$ (SST) turbulence model
  integer, save :: is2kw

  !> field id of the divergence of the velocity. More precisely,
  !> it is the trace of the velocity gradient (and not a finite volume
  !> divergence term). In the cell \c iel,  \f$div(\vect{u})\f$ is given
  !> by \c divukw(iel1). This array is defined only with the \f$k-\omega\f$ SST
  !> turbulence model (because in this case it may be calculated at the same
  !> time as \f$S^2\f$)
  integer, save :: idivukw

  !> field id of the strain rate tensor at the previous time step
  integer, save :: istraio

  !> \}
  !----------------------------------------------------------------------------
  ! Numerical properties
  !----------------------------------------------------------------------------

  !> \defgroup numerical_prop Numerical properties

  !> \addtogroup numerical_prop
  !> \{

  !> Weighting for gradient calculation on variables
  integer, save :: kwgrec

  !> \}
  !----------------------------------------------------------------------------
  ! Mapping to field structures
  !----------------------------------------------------------------------------

  !> \defgroup field_map Mapping to field structures

  !> \addtogroup field_map
  !> \{

  !> Field id for variable i
  integer, save :: ivarfl(nvarmx)

  !> Field id for the dttens tensor
  integer, save :: idtten

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving the number of species in the gas mix
    ! if gas mix model is enabled (igmix)

    subroutine cs_f_gas_mix_get_pointers(nscasp) &
      bind(C, name='cs_f_gas_mix_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nscasp
    end subroutine cs_f_gas_mix_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \addtogroup field_map
  !> \{

  !> \anchor iprpfl
  !> Identity function for compatibility with deprecated iprpfl array

  elemental pure function iprpfl(f_id) result(r_f_id)

    implicit none

    ! Parameters

    integer, intent(in) :: f_id
    integer             :: r_f_id

    ! Function body

    r_f_id = f_id

  end function iprpfl

  !> \}

  !> \brief Initialize Fortran gas mix API.
  !> This maps Fortran pointers to global C variables.

  subroutine gas_mix_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_nscasp

    call cs_f_gas_mix_get_pointers(c_nscasp)

    call c_f_pointer(c_nscasp, nscasp)

  end subroutine gas_mix_options_init

  !=============================================================================

  !> \}

end module numvar
