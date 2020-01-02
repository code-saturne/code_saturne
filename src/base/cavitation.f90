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

!> \file cavitation.f90
!> Module for cavitation modeling
!>
!> Please refer to the
!> <a href="../../theory.pdf#cavitation"><b>cavitation model</b></a> section
!> of the theory guide for more informations.

module cavitation

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \addtogroup cavitation
  !> \{

  !----------------------------------------------------------------------------
  ! Vaporization/condensation model
  !----------------------------------------------------------------------------

  !> \addtogroup cav_source_term
  !> \{

  !> reference saturation pressure (kg/(m s2))
  real(c_double), pointer, save :: presat

  !> reference velocity of the flow (m/s)
  real(c_double), pointer, save :: uinf

  !> reference length scale of the flow (m)
  real(c_double), pointer, save :: linf

  !> constant Cdest of the condensation source term (Merkle model)
  real(c_double), pointer, save :: cdest

  !> constant Cprod of the vaporization source term (Merkle model)
  real(c_double), pointer, save :: cprod

  !> \}

  !----------------------------------------------------------------------------
  ! Interaction with turbulence
  !----------------------------------------------------------------------------

  !> \addtogroup cav_turbulence
  !> \{

  !> activation of the eddy-viscosity correction (Reboud correction)
  !>    - 1: activated
  !>    - 0: desactivated
  integer(c_int), pointer, save :: icvevm

  !> constant mcav of the eddy-viscosity correction (Reboud correction)
  real(c_double), pointer, save :: mcav

  !> \}

  !----------------------------------------------------------------------------
  ! Numerical parameters
  !----------------------------------------------------------------------------

  !> \addtogroup cav_numerics
  !> \{

  !> implicitation in pressure of the vaporization/condensation model
  !>    - 1: activated
  !>    - 0: desactivated
  integer(c_int), pointer, save :: itscvi

  !> \}

  !> \}

  !=============================================================================

  interface

     !---------------------------------------------------------------------------

     !> \cond DOXYGEN_SHOULD_SKIP_THIS

     !---------------------------------------------------------------------------

     ! Interface to C function retrieving pointers to cavitation model indicator
     ! and parameters

     subroutine cs_f_cavitation_get_pointers(presat, uinf, linf, cdest, cprod, &
                                             icvevm, mcav, itscvi)             &
       bind(C, name='cs_f_cavitation_get_pointers')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(out) :: presat, uinf, linf, cdest, cprod
       type(c_ptr), intent(out) :: icvevm, mcav, itscvi
     end subroutine cs_f_cavitation_get_pointers

     !---------------------------------------------------------------------------

     !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

     !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran cavitation model API.
  !> This maps Fortran pointers to global C structure members and indicator.

  subroutine cavitation_model_init

    use, intrinsic :: iso_c_binding

    implicit none

    ! Local variables

    type(c_ptr) :: c_presat, c_uinf, c_linf, c_cdest, c_cprod
    type(c_ptr) :: c_icvevm, c_mcav, c_itscvi

    call cs_f_cavitation_get_pointers(c_presat, c_uinf, c_linf,           &
                                      c_cdest, c_cprod, c_icvevm, c_mcav, &
                                      c_itscvi)

    call c_f_pointer(c_presat, presat)
    call c_f_pointer(c_uinf, uinf)
    call c_f_pointer(c_linf, linf)
    call c_f_pointer(c_cdest, cdest)
    call c_f_pointer(c_cprod, cprod)
    call c_f_pointer(c_icvevm, icvevm)
    call c_f_pointer(c_mcav, mcav)
    call c_f_pointer(c_itscvi, itscvi)

  end subroutine cavitation_model_init

  !=============================================================================

  !> \brief Compute the vaporization source term
  !> \f$ \Gamma_V \left(\alpha, p\right) = m^+ + m^- \f$ using the
  !> Merkle model:
  !> \f[
  !> m^+ = -\dfrac{C_{prod} \rho_l \min \left( p-p_V,0 \right)\alpha(1-\alpha)}
  !>              {0.5\rho_lu_\infty^2t_\infty},
  !> \f]
  !> \f[
  !> m^- = -\dfrac{C_{dest} \rho_v \max \left( p-p_V,0 \right)\alpha(1-\alpha)}
  !>              {0.5\rho_lu_\infty^2t_\infty},
  !> \f]
  !> with \f$ C_{prod}, C_{dest} \f$ empirical constants,
  !> \f$ t_\infty=l_\infty/u_\infty \f$ a reference time scale and \f$p_V\f$
  !> the reference saturation pressure.
  !> \f$ l_\infty \f$, \f$ u_\infty \f$ and \f$p_V\f$ may be provided by
  !> the user (user function).
  !> Note that the r.h.s. of the void fraction transport equation is
  !> \f$ \Gamma_V/\rho_v \f$.

  !> \param[in]  pressure  Pressure array
  !> \param[in]  voidf     Void fraction array

  subroutine cavitation_compute_source_term(pressure, voidf)

    use optcal
    use pointe, only: gamcav, dgdpca
    use mesh, only: ncel, ncelet
    use vof

    ! Arguments

    double precision pressure(ncelet), voidf(ncelet)

    ! Local variables

    integer iel
    double precision tinf, cond, cvap, condens, vaporis

    if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) then

      ! Merkle model

      tinf = linf/uinf

      cond = (cdest*rho2)/(0.5d0*rho1*uinf*uinf*tinf)
      cvap = (cprod*rho1)/(0.5d0*rho1*uinf*uinf*tinf)

      do iel = 1, ncel
        condens = -cond*max(0.d0, pressure(iel) - presat) &
             *voidf(iel)*(1.d0 - voidf(iel))
        vaporis = -cvap*min(0.d0, pressure(iel) - presat) &
             *voidf(iel)*(1.d0 - voidf(iel))
        gamcav(iel) = condens + vaporis
        if (gamcav(iel).lt.0) then
          dgdpca(iel) = -cond*voidf(iel)*(1.d0 - voidf(iel))
        else
          dgdpca(iel) = -cvap*voidf(iel)*(1.d0 - voidf(iel))
        endif
      enddo

    endif

  end subroutine cavitation_compute_source_term

  !=============================================================================

  !> \brief Modify eddy viscosity using the Reboud correction:
  !>\f[
  !> \mu_t'= \dfrac{\rho_v + (1-\alpha)^{mcav}(\rho_l-\rho_v)}{\rho}\mu_t.
  !>\f]
  !>

  !> \param[in]      crom   density array
  !> \param[in]      voidf  void fraction array
  !> \param[in,out]  visct  turbulent viscosity

  subroutine cavitation_correct_visc_turb (crom, voidf, visct)

    use mesh
    use numvar
    use field
    use vof

    ! Arguments

    double precision crom(ncelet), voidf(ncelet)
    double precision visct(ncelet)

    ! Local variables

    integer iel
    double precision frho

    do iel = 1, ncel
      frho =  ( rho2 + (1.d0-voidf(iel))**mcav*(rho1 - rho2) ) &
             /max(crom(iel),1.d-12)
      visct(iel) = frho*visct(iel)
    enddo

  end subroutine cavitation_correct_visc_turb

  !=============================================================================

end module cavitation
