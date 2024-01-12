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
  integer, save :: ipr = 0

  !> \anchor iu
  !> velocity component \f$ u_x \f$
  integer, save :: iu = 0

  !> \anchor iv
  !> velocity component \f$ u_y \f$
  integer, save :: iv = 0

  !> \anchor iw
  !> velocity component \f$ u_z \f$
  integer, save :: iw = 0

  !> \anchor ivolf2
  !> void fraction for VOF method
  integer, save :: ivolf2 = 0

  !> \anchor ik
  !> turbulent kinetic energy \f$ k \f$
  integer, save :: ik = 0

  !> \anchor iep
  !> turbulent dissipation \f$ \varepsilon \f$
  integer, save :: iep = 0

  !> \anchor ir11
  !> Reynolds stress component \f$ R_{xx} \f$
  integer, save :: ir11 = 0

  !> \anchor ir22
  !> Reynolds stress component \f$ R_{yy} \f$
  integer, save :: ir22 = 0

  !> \anchor ir33
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir33 = 0

  !> \anchor ir12
  !> Reynolds stress component \f$ R_{xy} \f$
  integer, save :: ir12 = 0

  !> \anchor ir23
  !> Reynolds stress component \f$ R_{yz} \f$
  integer, save :: ir23 = 0

  !> \anchor ir13
  !> Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir13 = 0

  !> \anchor irij
  !> Reynolds stress tenso \f$ R_{ij} \f$
  integer, save :: irij = 0

  !> \anchor iphi
  !> variable \f$ \phi \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: iphi = 0

  !> \anchor ifb
  !> variable \f$ f_b \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: ifb = 0

  !> \anchor ial
  !> variable \f$ \alpha \f$ of the \f$ Bl-v^2-k \f$ model
  integer, save :: ial = 0

  !> \anchor iomg
  !> variable \f$ \omega \f$ of the \f$ k-\omega \f$ SST
  integer, save :: iomg = 0

  !> \anchor inusa
  !> variable \f$ \widetilde{\nu}_T \f$ of the Spalart Allmaras
  integer, save :: inusa = 0

  !> \anchor isca
  !> isca(i) is the index of the scalar i
  integer, save :: isca(nscamx)

  !> \anchor iscapp
  !> iscapp(i) is the index of the specific physics scalar i
  integer, save :: iscapp(nscamx)

  !> \anchor nscaus
  !> number of user scalars solutions of an advection equation
  integer, save :: nscaus = 0

  !> \anchor nscapp
  !> number of specific physics scalars
  integer, save :: nscapp = 0

  !> \anchor nscasp
  !> number of species scalars
  integer(c_int), pointer, save :: nscasp

  !> \anchor iuma
  !> mesh velocity component \f$ w_x \f$
  integer, save :: iuma = 0

  !> \anchor ivma
  !> mesh velocity component \f$ w_y \f$
  integer, save :: ivma = 0

  !> \anchor iwma
  !> mesh velocity component \f$ w_z \f$
  integer, save :: iwma = 0

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
  integer, save :: irom = -1

  !> dynamic molecular viscosity (in kg/(m.s))
  integer, save :: iviscl = -1

  !> dynamic turbulent viscosity
  integer, save :: ivisct = -1

  !> interior and boundary convective mass flux key ids of the variables
  integer, save :: kimasf = -1, kbmasf = -1

  !> constant diffusivity field id key for scalars
  integer, save :: kvisl0 = -1

  !> variable diffusivity field id key for scalars
  integer, save :: kivisl = -1

  !> do scalars behave as a temperature (regarding multiplication by Cp) ?
  integer, save :: kscacp = -1

  !> variable density field id key for scalars
  integer, save :: kromsl = -1

  !> source terms at previous time step for 2nd order
  integer, save :: kstprv = -1

  !> source terms at the current time step (used for limiters)
  integer, save :: kst = -1

  !> turbulent schmidt key for scalars
  integer, save :: ksigmas = -1

  !> turbulent flow models GGDH and AFM key for scalar
  integer, save :: kctheta = -1

  !> convective mass flux of the variables at the previous time-step
  integer, save :: ifluaa(nvarmx)

  !> cell density field ids of the variables
  integer, save :: icrom = -1

  !> boundary density field ids of the variables
  integer, save :: ibrom = -1

  !> field ids of the cell porosity
  integer, save :: ipori = -1, iporf = -1

  !> dynamic constant of Smagorinsky
  integer, save :: ismago = -1

  !> field ids of the anisotropic viscosity
  !> \remark turbulent or Darcy module anisotropic diffusion
  integer, save :: ivsten = -1, ivstes = -1

  !> Courant number
  integer, save :: icour = -1

  !> Fourier number
  integer, save :: ifour = -1

  !> Total pressure at cell centers
  !> \f$ P_{tot} = P^\star +\rho \vect{g} \cdot (\vect{x}-\vect{x}_0) \f$
  integer, save :: iprtot = -1

  !> Mesh velocity viscosity for the ALE module
  !> \remark might be orthotropic
  integer, save :: ivisma = -1

  !> pointer for dilatation source terms
  integer, save :: iustdy(nscamx)

  !> pointer for global dilatation source terms
  integer, save :: itsrho = -1

  !> pointer for thermal expansion coefficient
  integer, save :: ibeta = -1

  !> pointer for deduced mass fraction in case of gas mix
  integer, save :: iddgas = -1

  !> pointer for gas mix molar mass
  integer, save :: igmxml = -1

  !> field id of the stresses at boundary  (if post-processed)
  integer, save :: iforbr = -1

  !>  field id of \f$y^+\f$ at boundary (if post-processed)
  integer, save :: iyplbr = -1

  !>  field id of temperature at boundary
  integer, save ::  itempb = -1

  !> \}
  !----------------------------------------------------------------------------
  ! Numerical properties
  !----------------------------------------------------------------------------

  !> \defgroup numerical_prop Numerical properties

  !> \addtogroup numerical_prop
  !> \{

  !> Weighting for gradient calculation on variables
  integer, save :: kwgrec = -1

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
  integer, save :: idtten = -1

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
