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

!> \file cavitation.f90
!> Module for cavitation modeling

module cavitation

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup cavitation Module for cavitating flow modelling

  !> \addtogroup cavitation
  !> \{

  !----------------------------------------------------------------------------
  ! Homogeneous mixture physical properties
  !----------------------------------------------------------------------------

  !> \defgroup cav_mixture_properties Mixture properties

  !> \addtogroup cav_mixture_properties
  !> \{

  !> reference density of the liquid phase in kg/m3
  double precision, save :: rol

  !> reference density of the gas phase in kg/m3
  double precision, save :: rov

  !> reference molecular viscosity of the liquid phase kg/(m s)
  double precision, save :: mul

  !> reference molecular viscosity of the gas phase kg/(m s)
  double precision, save :: muv

  !> \}

  !----------------------------------------------------------------------------
  ! Vaporization/condensation model
  !----------------------------------------------------------------------------

  !> \defgroup cav_source_term Vaporization/condensation model

  !> \addtogroup cav_source_term
  !> \{

  !> reference saturation pressure in kg/(m s2)
  double precision, save :: presat

  !> reference length scale of the flow in meters
  double precision, save :: linf

  !> reference velocity of the flow in m/s
  double precision, save :: uinf

  !> constant Cprod of the vaporization source term (Merkle model)
  double precision, save :: cprod

  !> constant Cdest of the condensation source term (Merkle model)
  double precision, save :: cdest

  !> \}

  !----------------------------------------------------------------------------
  ! Interaction with turbulence
  !----------------------------------------------------------------------------

  !> \defgroup cav_turbulence Interaction with turbulence

  !> \addtogroup cav_turbulence
  !> \{

  !> activation of the eddy-viscosity correction (Reboud correction)
  !>    - 1: activated
  !>    - 0: desactivated
  integer,          save :: icvevm

  !> constant mcav of the eddy-viscosity correction (Reboud correction)
  double precision, save :: mcav

  !> \}

  !----------------------------------------------------------------------------
  ! Numerical parameters
  !----------------------------------------------------------------------------

  !> \defgroup cav_numerics Numerical parameters

  !> \addtogroup cav_numerics
  !> \{

  !> implicitation in pressure of the vaporization/condensation model
  !>    - 1: activated
  !>    - 0: desactivated
  integer,          save :: itscvi

  !> clipping min. for the void fraction
  double precision, save :: clvfmn

  !> clipping max. for the void fraction
  double precision, save :: clvfmx

  !> \}

  !> \}

contains

  !=============================================================================

  !> \brief  Default initialization of the module variables.

  subroutine init_cavitation

    use cstphy
    use field
    use numvar

    ! Mixture physical properties
    !----------------------------

    rol = 1.d3
    rov = 1.d0
    mul = 1.d-3
    muv = 1.d-5

    ! Vaporization/condensation model parameters
    !-------------------------------------------

    ! physical variables
    presat = 2.d3
    uinf = uref
    linf = 0.1d0

    ! Merkle model constants
    cdest = 5.d1
    cprod = 1.d4

    ! Interaction with turbulence
    !----------------------------

    ! Activate the eddy-viscosity correction (Reboud correction)
    icvevm = 1

    ! Reboud correction constant
    mcav = 10.d0

    ! Numercial parameters
    !---------------------

    ! implicitation in pressure
    itscvi = 1

    ! min/max clippings for the void fraction
    clvfmn = 1.d-12
    clvfmx = 1.d0 - clvfmn

  end subroutine init_cavitation

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

    ! Arguments

    double precision pressure(ncelet), voidf(ncelet)

    ! Local variables

    integer iel
    double precision tinf, cond, cvap, condens, vaporis

    if (icavit.eq.0) then

      ! No model

      do iel = 1, ncel
        gamcav(iel) = 0.d0
        dgdpca(iel) = 0.d0
      enddo

    elseif (icavit.eq.1) then

      ! Merkle model

      tinf = linf/uinf

      cond = (cdest*rov)/(0.5d0*rol*uinf*uinf*tinf)
      cvap = (cprod*rol)/(0.5d0*rol*uinf*uinf*tinf)

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

  !> \brief Compute the mixture density, mixture dynamic viscosity and mixture
  !>        mass flux given the volumetric flux, the void fraction and the
  !>        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
  !>        (liquid), \f$ \rho_v, \mu_v \f$ (gas).
  !>        One have:
  !> \f[
  !> \rho_\celli = \alpha_\celli \rho_v + (1-\alpha_\celli) \rho_l,
  !> \f]
  !> \f[
  !> \mu_\celli = \alpha_\celli \mu_v + (1-\alpha_\celli) \mu_l,
  !> \f]
  !> \f[
  !> \left( \rho\vect{u}\cdot\vect{S} \right)_\ij = \\ \left\lbrace
  !> \begin{array}{ll}
  !>   \rho_\celli (\vect{u}\cdot\vect{S})_\ij
  !>  &\text{ if } (\vect{u}\cdot\vect{S})_\ij>0, \\
  !>   \rho_\cellj (\vect{u}\cdot\vect{S})_\ij
  !>  &\text{ otherwise },
  !> \end{array} \right.
  !> \f]
  !> \f[
  !> \left( \rho\vect{u}\cdot\vect{S} \right)_\ib = \\ \left\lbrace
  !> \begin{array}{ll}
  !>   \rho_\celli (\vect{u}\cdot\vect{S})_\ib
  !>  &\text{ if } (\vect{u}\cdot\vect{S})_\ib>0, \\
  !>   \rho_b (\vect{u}\cdot\vect{S})_\ib
  !>  &\text{ otherwise }.
  !> \end{array} \right.
  !> \f]
  !>

  !> \param[in]  voidf  Void fraction array
  !> \param[in]  coavoi Void fraction boundary coefficient array
  !> \param[in]  cobvoi Void fraction boundary coefficient array
  !> \param[in]  ivoifl Volumetric flux at internal faces array
  !> \param[in]  bvoifl Volumetric flux at boundary faces array
  !> \param[out] crom   Density at cell center array
  !> \param[out] brom   Density at boudary faces array
  !> \param[out] viscl  Dynamic viscosity array
  !> \param[out] imasfl Mass flux at internal faces array
  !> \param[out] bmasfl Mass flux at internal boundary faces array

  subroutine cavitation_update_phys_prop  &
            ( voidf, coavoi, cobvoi, ivoifl, bvoifl, &
              crom, brom, viscl, imasfl, bmasfl )

    use paramx
    use pointe, only: itypfb
    use mesh

    ! Arguments

    double precision voidf(ncelet)
    double precision coavoi(nfabor), cobvoi(nfabor)
    double precision ivoifl(nfac), bvoifl(nfabor)
    double precision crom(ncelet), brom(nfabor)
    double precision viscl(ncelet)
    double precision imasfl(nfac), bmasfl(nfabor)

    ! Local variables

    integer iel, ifac, ii, jj
    double precision bvoidf, flui, fluj, flub

    ! Update mixture density

    do iel = 1, ncelet
      crom(iel) = rov*voidf(iel) + rol*(1.d0 - voidf(iel))
    enddo

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      bvoidf = coavoi(ifac) + cobvoi(ifac)*voidf(iel)
      brom(ifac) = rov*bvoidf + rol*(1.d0 - bvoidf)
    enddo

    ! Update mixture viscosity

    do iel = 1 , ncelet
      viscl(iel) = muv*voidf(iel) + mul*(1.d0 - voidf(iel))
    enddo

    ! Update mass flux

    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      flui = 0.5d0*( ivoifl(ifac) + abs (ivoifl(ifac)) )
      fluj = 0.5d0*( ivoifl(ifac) - abs (ivoifl(ifac)) )

      imasfl(ifac) = imasfl(ifac) + ( flui*crom(ii) + fluj*crom(jj) )

    enddo

    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparoi .or. itypfb(ifac) .eq. isymet) then
        bmasfl(ifac) = 0.d0
      else
        iel = ifabor(ifac)
        flui = 0.5d0*( bvoifl(ifac) + abs(bvoifl(ifac)) )
        flub = 0.5d0*( bvoifl(ifac) - abs(bvoifl(ifac)) )

        bmasfl(ifac) = bmasfl(ifac) + ( flui*crom(iel) + flub*brom(ifac) )

      endif
    enddo

  end subroutine cavitation_update_phys_prop

  !=============================================================================

  !> \brief Modify eddy viscosity using the Reboud correction:
  !>\f[
  !> \mu_t'= \dfrac{\rho_v + (1-\alpha)^{mcav}(\rho_l-\rho_v)}{\rho}\mu_t.
  !>\f]
  !>

  !> \param[in]      crom   Density array
  !> \param[in]      voidf  Void fraction array
  !> \param[in,out]  visct  Turbulent viscosity array

  subroutine cavitation_correct_visc_turb (crom, voidf, visct)

    use mesh

    ! Arguments

    double precision crom(ncelet), voidf(ncelet)
    double precision visct(ncelet)

    ! Local variables

    integer iel
    double precision frho

    do iel = 1, ncel
      frho = ( rov + (1.d0-voidf(iel))**mcav*(rol - rov) )/max(crom(iel),1.d-12)
      visct(iel) = frho*visct(iel)
    enddo

  end subroutine cavitation_correct_visc_turb

  !=============================================================================

  !> \brief Print the global mixture mass budget:
  !> \f[
  !> \sum_i\left(
  !> |\Omega_i|\dfrac{\alpha_i^n - \alpha_i^{n-1}}{\Delta t} +
  !> \sum_{j\in\Face{\celli}}\left(\rho\vect{u}\vect{S}\right)_{ij}^n
  !> \right).
  !> \f]
  !>

  !> \param[in]  crom        Density at cell centers at current time step
  !> \param[in]  croma       Density at cell centers at previous time step
  !> \param[in]  brom        Density at boundary faces at current time step
  !> \param[in]  dt          Time step
  !> \param[in]  imasfl_rel  Mass flux at internal faces array
  !> \param[in]  bmasfl_rel  Mass flux at internal boundary faces array

  subroutine cavitation_print_mass_budget  &
            ( crom, croma, brom, dt, imasfl_rel, bmasfl_rel )

    use cstphy
    use entsor
    use mesh
    use parall
    use rotation
    use turbomachinery

    ! Arguments

    double precision crom(ncelet), croma(ncelet)
    double precision brom(nfabor)
    double precision dt(ncelet)
    double precision, dimension(:), target :: imasfl_rel, bmasfl_rel

    ! Local variables

    integer init, iel, iel1, iel2, ifac
    double precision bilglo, rhofac, vr(3), vr1(3), vr2(3)
    double precision, dimension(:), allocatable :: divro, tinsro
    double precision, dimension(:), allocatable, target :: imasfl_abs, bmasfl_abs
    double precision, dimension(:), pointer :: imasfl, bmasfl

    ! Initialization
    allocate(divro(ncelet))
    allocate(tinsro(ncelet))
    do iel = 1, ncel
      divro(iel)  = 0.d0
      tinsro(iel) = 0.d0
    enddo

    if (.not.(icorio.eq.1.or.iturbo.eq.1)) then
      imasfl => imasfl_rel(1:nfac)
      bmasfl => bmasfl_rel(1:nfabor)
    else
      allocate(imasfl_abs(nfac), bmasfl_abs(nfabor))
      !$omp parallel do private(iel1, iel2, rhofac)
      do ifac = 1, nfac
        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)
        if (irotce(iel1).ne.0 .or. irotce(iel2).ne.0) then

          rhofac = 0.5d0*(crom(iel1) + crom(iel2))
          call rotation_velocity(irotce(iel1), cdgfac(:,ifac), vr1)
          call rotation_velocity(irotce(iel2), cdgfac(:,ifac), vr2)

          imasfl_abs(ifac) = imasfl_rel(ifac) + 0.5d0 *rhofac*(         &
                                  surfac(1,ifac)*(vr1(1) + vr2(1))      &
                                + surfac(2,ifac)*(vr1(2) + vr2(2))      &
                                + surfac(3,ifac)*(vr1(3) + vr2(3)) )
        else
          imasfl_abs(ifac) = imasfl_rel(ifac)
        endif
      enddo
      !$omp parallel do private(iel, rhofac) &
      !$omp          if(nfabor > thr_n_min)
      do ifac = 1, nfabor
        iel = ifabor(ifac)
        if (irotce(iel).ne.0) then

          rhofac = brom(ifac)
          call rotation_velocity(irotce(iel), cdgfbo(:,ifac), vr)

          bmasfl_abs(ifac) = bmasfl_rel(ifac) + rhofac*( surfbo(1,ifac)*vr(1) &
                                              + surfbo(2,ifac)*vr(2)     &
                                              + surfbo(3,ifac)*vr(3) )
        else
          bmasfl_abs(ifac) = bmasfl_rel(ifac)
        endif
      enddo
      imasfl => imasfl_abs(1:nfac)
      bmasfl => bmasfl_abs(1:nfabor)
    endif

    ! (Absolute) Mass flux divergence
    init = 1
    call divmas(init, imasfl, bmasfl, divro)

    ! Unsteady term
    do iel = 1, ncel
      tinsro(iel) = volume(iel)*(crom(iel)-croma(iel))/dt(iel)
    enddo

    ! Budget
    bilglo = 0.d0
    do iel = 1, ncel
      bilglo = bilglo + tinsro(iel) + divro(iel)
    enddo
    if (irangp.ge.0)  call parsom(bilglo)

    ! Printings
    write(nfecra,1000) bilglo

    ! Finalization
    deallocate(divro, tinsro)
    if(allocated(imasfl_abs)) deallocate(imasfl_abs, bmasfl_abs)

    !Format
#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** MODELE DE CAVITATION'                                   ,/,&
'      --------------------'                                   ,/,&
'   Bilan de masse global du melange :',e12.4,                  /)

#else

 1000 format(/,                                                   &
'   ** CAVITATION MODELLING'                                   ,/,&
'      --------------------'                                   ,/,&
'   Mixture global mass budget:',e12.4,                         /)

#endif

    ! End

  end subroutine cavitation_print_mass_budget

  !=============================================================================

end module cavitation
