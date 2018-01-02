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

!> \file vof.f90
!> Module for Volume-Of-Fluid method

module vof

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup VOF Module for free surface flow modelling

  !> \addtogroup vof
  !> \{

  !----------------------------------------------------------------------------
  ! Homogeneous mixture physical properties
  !----------------------------------------------------------------------------

  !> \defgroup vof_mixture_properties Mixture properties

  !> \addtogroup vof_mixture_properties
  !> \{

  !> reference density of fluid 1 (kg/m3).
  !> By convention, liquid phase for cavitation model.
  double precision, save :: rho1

  !> reference density of fluid 2 (kg/m3).
  !> By convention, gas phase for cavitation model.
  double precision, save :: rho2

  !> reference molecular viscosity of fluid 1 (kg/(m s))
  double precision, save :: mu1

  !> reference molecular viscosity of fluid 2 (kg/(m s))
  double precision, save :: mu2

  !> \}

  !----------------------------------------------------------------------------
  ! Numerical parameters
  !----------------------------------------------------------------------------

  !> \defgroup vof_numerics Numerical parameters

  !> \addtogroup vof_numerics
  !> \{

  !> clipping min. for the volume fraction
  double precision, save :: clvfmn

  !> clipping max. for the volume fraction
  double precision, save :: clvfmx

  !> \}

  !> \}

contains

  !=============================================================================

  !> \brief  Default initialization of the module variables.

  subroutine init_vof

    use cstnum

    ! Mixture physical properties
    !----------------------------

    rho1 = 1.d3
    rho2  = 1.d0
    mu1  = 1.d-3
    mu2   = 1.d-5

    ! Numerical parameters
    !---------------------

    ! min/max clippings for the volume fraction
    clvfmn = -grand*10.d0
    clvfmx = grand*10.d0

  end subroutine init_vof

  !=============================================================================

  !> \brief Compute the mixture density, mixture dynamic viscosity and mixture
  !>        mass flux given the volumetric flux, the volume fraction and the
  !>        reference density and dynamic viscosity \f$ \rho_l, \mu_l \f$
  !>        (liquid), \f$ \rho_v, \mu_v \f$ (gas) as follows:
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
  !> \param[out] imasfl Mass flux at internal faces array
  !> \param[out] bmasfl Mass flux at internal boundary faces array

  subroutine vof_update_phys_prop  &
            ( voidf, coavoi, cobvoi, ivoifl, bvoifl, &
              crom, brom, imasfl, bmasfl )

    use field
    use paramx
    use pointe, only: itypfb
    use numvar, only: iviscl
    use mesh

    ! Arguments

    double precision voidf(ncelet)
    double precision coavoi(nfabor), cobvoi(nfabor)
    double precision ivoifl(nfac), bvoifl(nfabor)
    double precision crom(ncelet), brom(nfabor)
    double precision imasfl(nfac), bmasfl(nfabor)

    ! Local variables

    integer iel, ifac, ii, jj
    double precision bvoidf, flui, fluj, flub, evof
    double precision, dimension(:), pointer :: viscl

    call field_get_val_s(iviscl, viscl)

    ! Update mixture density and viscocity on cells

    do iel = 1, ncelet
      evof = voidf(iel)
      crom(iel) = rho2*evof + rho1*(1.d0 - evof)
      viscl(iel) = mu2*evof + mu1*(1.d0 - evof)
    enddo

    ! Update mixture density on boundary faces

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      evof = voidf(iel)
      bvoidf = coavoi(ifac) + cobvoi(ifac)*evof
      brom(ifac) = rho2*bvoidf + rho1*(1.d0 - bvoidf)
    enddo

    ! Update mass flux

    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      flui = 0.5d0*(ivoifl(ifac) + abs(ivoifl(ifac)))
      fluj = 0.5d0*(ivoifl(ifac) - abs(ivoifl(ifac)))
      imasfl(ifac) = imasfl(ifac) + (flui*crom(ii) + fluj*crom(jj))
    enddo

    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.isymet) then
        bmasfl(ifac) = 0.d0
      else
        iel = ifabor(ifac)
        flui = 0.5d0*(bvoifl(ifac) + abs(bvoifl(ifac)))
        flub = 0.5d0*(bvoifl(ifac) - abs(bvoifl(ifac)))
        bmasfl(ifac) = bmasfl(ifac) + (flui*crom(iel) + flub*brom(ifac))
      endif
    enddo

  end subroutine vof_update_phys_prop

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

  subroutine vof_print_mass_budget  &
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
'   ** ALGORITHME VOF      '                                   ,/,&
'      --------------------'                                   ,/,&
'   Bilan de masse global du melange :',e12.4,                  /)

#else

 1000 format(/,                                                   &
'   ** VOF ALGORITHM       '                                   ,/,&
'      --------------------'                                   ,/,&
'   Mixture global mass budget:',e12.4,                         /)

#endif

    ! End

  end subroutine vof_print_mass_budget

  !=============================================================================

end module vof
