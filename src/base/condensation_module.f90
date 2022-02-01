module condensation_module

  implicit none

  !> Constants for the correlation of steam saturated pressure
  double precision, parameter :: pr_c = 221.2d+5
  double precision, parameter :: t_c = 647.3d0
  double precision, parameter :: patm = 101320.0d0
  double precision, parameter :: C_k1 = -  7.691234564d0
  double precision, parameter :: C_k2 = - 26.08023696d0
  double precision, parameter :: C_k3 = -168.1706546d0
  double precision, parameter :: C_k4 =   64.23285504d0
  double precision, parameter :: C_k5 = -118.9646225d0
  double precision, parameter :: C_k6 =    4.16711732d0
  double precision, parameter :: C_k7 =   20.9750676d0
  double precision, parameter :: C_k8 = -  1.d+9
  double precision, parameter :: C_k9 =    6.d0

  !> Caracteristic length
  double precision, parameter :: l_cell = 1.0d0

  !> Vaporization latent heat
  double precision, parameter :: lcond = 2278.0d+3

  contains

!>===================================================================================
  subroutine compute_mix_properties(length, mix_mol_mas, mol_mas_ncond, x_h2o_g, diff_coeff)
!>===================================================================================

    use numvar, only: nscasp, ivarfl, isca
    use cstnum, only: zero
    use field, only: field_get_val_s, field_get_val_s_by_name, field_get_id
    use cs_c_bindings, only: gas_mix_species_prop, field_get_key_struct_gas_mix_species_prop
    use optcal, only: idilat, iscasp, iscalt
    use cstphy, only: icp, pther, p0
    use entsor, only: nfecra

    implicit none

    !> inout
    integer, intent(in) :: length
    double precision, intent(inout) :: mix_mol_mas(length)
    double precision, intent(inout) :: mol_mas_ncond(length)
    double precision, intent(inout) :: x_h2o_g(length)
    double precision, intent(inout) :: diff_coeff(length)

    !> Local
    integer :: iel, f_id, iesp
    double precision, dimension(:), pointer :: y_h2o_g
    double precision, dimension(:), pointer :: cvar_yk, cvar_enth, cpro_cp
    type(gas_mix_species_prop) :: s_h2o_g, s_k

    double precision :: temperature, pressure
    double precision, dimension(1:nscasp+1) :: y_k, mol_mas_k, vol_dif_k

    ! Initialize fields
    call field_get_val_s_by_name("y_h2o_g", y_h2o_g)
    call field_get_id("y_h2o_g", f_id)
    call field_get_key_struct_gas_mix_species_prop(f_id, s_h2o_g)
    call field_get_val_s(ivarfl(isca(iscalt)), cvar_enth)
    if (icp.ge.0) then
      call field_get_val_s(icp, cpro_cp)
    else
      write(nfecra,1000) icp
      call csexit (1)
    endif

    do iesp = 1, nscasp
      call field_get_key_struct_gas_mix_species_prop( &
           ivarfl(isca(iscasp(iesp))), s_k)
      mol_mas_k(iesp) = s_k%mol_mas
      vol_dif_k(iesp) = s_k%vol_dif
    enddo
    mol_mas_k(nscasp+1) = s_h2o_g%mol_mas
    vol_dif_k(nscasp+1) = s_h2o_g%vol_dif

    do iel = 1, length
      ! Water
      y_k(nscasp+1) = 1.0d0
      do iesp = 1, nscasp
        call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)
        y_k(iesp) = cvar_yk(iel)
        y_k(nscasp+1) = y_k(nscasp+1) - y_k(iesp) ! steam fraction
      enddo
      ! Get global mixture molecular weight
      call compute_mix_mol_mas(nscasp+1, y_k, mol_mas_k, mix_mol_mas(iel))
      ! Get mole fraction of steam
      call compute_mole_fraction(y_k(nscasp+1), mix_mol_mas(iel), mol_mas_k(nscasp+1), x_h2o_g(iel))
      ! Get non condensable mix molecular weight
      call compute_mix_mol_mas(nscasp, y_k(1:nscasp), mol_mas_k(1:nscasp), mol_mas_ncond(iel))
      mol_mas_ncond(iel) = mol_mas_ncond(iel) * (1.0d0 - y_k(nscasp+1))
      if (idilat == 3) then
        pressure = pther
      else
        pressure = p0
      endif
      temperature = cvar_enth(iel)/cpro_cp(iel)
      call compute_steam_binary_diffusion(nscasp, y_k(1:nscasp), mol_mas_k(1:nscasp), &
                                            vol_dif_k(1:nscasp), mol_mas_k(nscasp+1), &
                                            vol_dif_k(nscasp+1), mix_mol_mas(iel), temperature, &
                                            pressure, diff_coeff(iel))

    enddo


 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      copain model prescribes a variable specific heat.',/,     &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu',/,                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
!----
! End
!----

  end subroutine compute_mix_properties
!>===================================================================================

!>===================================================================================
  subroutine compute_mix_mol_mas(nb_species, y_k, mol_mas_k, mix_mol_mas)
!>===================================================================================

    implicit none
   ! Number of non condensable species
   ! Mass fractions of non condensable species
   ! Molecular weights of non condensable species
   ! Molecular weight of steam
   ! Molecular weight of steam and non condensable mixture
    !> inout
    integer, intent(in) :: nb_species
    double precision, intent(in) :: y_k(nb_species)
    double precision, intent(in) :: mol_mas_k(nb_species)
    double precision, intent(inout) :: mix_mol_mas

    !> local
    double precision :: y_steam
    integer :: k

    mix_mol_mas = 0.0d0

    do k = 1, nb_species
      mix_mol_mas = mix_mol_mas + y_k(k) / mol_mas_k(k)
    enddo

    mix_mol_mas = 1.0d0 / mix_mol_mas

  end subroutine compute_mix_mol_mas
!>===================================================================================

!>===================================================================================
  subroutine compute_mole_fraction(mass_fraction, mix_mol_mas, mol_mas, mole_fraction)
!>===================================================================================

    implicit none

    double precision, intent(in) :: mass_fraction, mix_mol_mas, mol_mas
    double precision, intent(out) :: mole_fraction

    mole_fraction = mass_fraction*mix_mol_mas/mol_mas

  end subroutine compute_mole_fraction
!>===================================================================================

!>===================================================================================
  subroutine compute_steam_binary_diffusion(nb_ncond, y_ncond_k, mol_mas_ncond_k, &
                                            vol_dif_ncond_k, mol_mas_steam, &
                                            vol_dif_steam, mix_mol_mas, temperature, &
                                            pressure, diffusion)
!>===================================================================================

    implicit none

    !> inout
    integer, intent(in) :: nb_ncond
    double precision, intent(in) :: y_ncond_k(nb_ncond)
    double precision, intent(in) :: mol_mas_ncond_k(nb_ncond)
    double precision, intent(in) :: vol_dif_ncond_k(nb_ncond)
    double precision, intent(in) :: mix_mol_mas, mol_mas_steam, vol_dif_steam
    double precision, intent(in) :: temperature
    double precision, intent(in) :: pressure
    double precision, intent(out) :: diffusion

    !> local
    integer :: k
    double precision :: ratio_tkpr, xmab, xvab, a1
    double precision :: x_ncond_tot, x_k

    diffusion = 0.0d0
    x_ncond_tot = 0.0d0
    do k = 1, nb_ncond
      call compute_mole_fraction(y_ncond_k(k),&
                                 mix_mol_mas, &
                                 mol_mas_ncond_k(k),&
                                 x_k)
      ratio_tkpr = temperature**1.75d0 / pressure
      xmab = sqrt(2.d0/( 1.d0/(mol_mas_steam*1000.d0) &
                        +1.d0/(mol_mas_ncond_k(k)*1000.d0) ) )
      xvab = ( vol_dif_steam**(1.d0/3.d0) &
              + vol_dif_ncond_k(k)**(1.d0/3.d0) )**2.d0
      a1   = 1.43d-7/(xmab*xvab)*patm
      diffusion = diffusion + x_k/(a1*ratio_tkpr)
      x_ncond_tot = x_ncond_tot + x_k
    enddo

    diffusion = x_ncond_tot / diffusion

  end subroutine compute_steam_binary_diffusion
!>===================================================================================


!>===================================================================================
  subroutine compute_psat(temperature, psat)
!>===================================================================================

    implicit none

    double precision, intent(in) :: temperature ! in kelvins
    double precision, intent(inout) :: psat

    double precision :: dtheta

    dtheta = temperature / t_c

    psat  = pr_c*exp((1.d0/dtheta)       &
          *( C_k1*(1.d0-dtheta)          &
            +C_k2*(1.d0-dtheta)**2       &
            +C_k3*(1.d0-dtheta)**3       &
            +C_k4*(1.d0-dtheta)**4       &
            +C_k5*(1.d0-dtheta)**5)      &
          / (1.d0+C_k6*(1.d0-dtheta)     &
                 +C_k7*(1.d0-dtheta)**2) &
        -(1.d0-dtheta)/(C_k8*(1.d0-dtheta)**2+C_k9))

  end subroutine compute_psat
!>===================================================================================

!>===================================================================================
subroutine compute_mac_adams(theta, grashof, schmidt, sherwood)
!>===================================================================================
  implicit none
  double precision, intent(in) :: theta, grashof, schmidt
  double precision, intent(out) :: sherwood
  sherwood = theta*0.13d0*(grashof*schmidt)**(1.d0/3.d0)
end subroutine compute_mac_adams
!>===================================================================================

!>===================================================================================
subroutine compute_schlichting(theta, reynolds, schmidt, sherwood)
!>===================================================================================
  implicit none
  double precision, intent(in) :: theta, reynolds, schmidt
  double precision, intent(out) :: sherwood
  sherwood = theta*0.0296d0*(reynolds**0.8d0)*(schmidt**(1.d0/3.d0))
end subroutine compute_schlichting
!>===================================================================================

!>===================================================================================
subroutine compute_incropera_1(theta, reynolds, grashof, schmidt, sherwood)
!>===================================================================================
  implicit none
  !> inout
  double precision, intent(in) :: theta, reynolds, grashof, schmidt
  double precision, intent(out) :: sherwood
  !> local
  double precision :: sherwood_natural, sherwood_forced
  call compute_mac_adams(theta, grashof, schmidt, sherwood_natural)
  call compute_schlichting(theta, reynolds, schmidt, sherwood_forced)
  sherwood = (abs(sherwood_forced**3.0d0 - sherwood_natural**3.0d0))**(1.0d0/3.0d0)
end subroutine compute_incropera_1
!>===================================================================================

!>===================================================================================
subroutine compute_incropera_2(theta, reynolds, grashof, schmidt, sherwood)
!>===================================================================================
  implicit none
  !> inout
  double precision, intent(in) :: theta, reynolds, grashof, schmidt
  double precision, intent(out) :: sherwood
  !> local
  double precision :: sherwood_natural, sherwood_forced
  call compute_mac_adams(theta, grashof, schmidt, sherwood_natural)
  call compute_schlichting(theta, reynolds, schmidt, sherwood_forced)
  sherwood = (abs(sherwood_forced**3.0d0 + sherwood_natural**3.0d0))**(1.0d0/3.0d0)
end subroutine compute_incropera_2
!>===================================================================================

subroutine compute_exchange_adimensional(theta, reynolds, grashof, schmidt_or_prandtl, convection_regime, sherwood_or_nusselt)
  implicit none
  !> inout
  double precision, intent(in) :: theta, reynolds, grashof, schmidt_or_prandtl
  integer, intent(in) :: convection_regime
  double precision, intent(out) :: sherwood_or_nusselt
  !> local
  double precision :: forced_value, natural_value
  if (convection_regime == 1) then
    call compute_mac_adams(theta, grashof, schmidt_or_prandtl, sherwood_or_nusselt)
  else if (convection_regime == 2) then
    call compute_schlichting(theta, reynolds, schmidt_or_prandtl, sherwood_or_nusselt)
  else if (convection_regime == 3) then
    call compute_incropera_1(theta, reynolds, grashof, schmidt_or_prandtl, sherwood_or_nusselt)
  else if (convection_regime == 4) then
    call compute_incropera_2(theta, reynolds, grashof, schmidt_or_prandtl, sherwood_or_nusselt)
  else if (convection_regime == 5) then
    call compute_mac_adams(theta, grashof, schmidt_or_prandtl, natural_value)
    call compute_schlichting(theta, reynolds, schmidt_or_prandtl, forced_value)
    sherwood_or_nusselt = max(forced_value, natural_value)
  endif
end subroutine compute_exchange_adimensional

!>===================================================================================
subroutine compute_schmidt(kin_viscosity, mass_diffusivity, schmidt)
!>===================================================================================
  implicit none
  double precision, intent(in) :: kin_viscosity, mass_diffusivity
  double precision, intent(out) :: schmidt
  schmidt = kin_viscosity/mass_diffusivity
end subroutine compute_schmidt
!>===================================================================================

!>===================================================================================
subroutine compute_nusselt_natural_convection(theta, grashof, prandtl, nusselt)
!>===================================================================================
  implicit none
  double precision, intent(in) :: theta, grashof, prandtl
  double precision, intent(out) :: nusselt
  nusselt = theta*(0.13d0*(grashof*prandtl)**(1.d0/3.d0))
end subroutine compute_nusselt_natural_convection
!>===================================================================================

!>===================================================================================
subroutine compute_prandtl(dyn_viscosity, lambda_over_cp, prandtl)
!>===================================================================================
  implicit none
  double precision, intent(in) :: dyn_viscosity, lambda_over_cp
  double precision, intent(out) :: prandtl
  prandtl = dyn_viscosity / lambda_over_cp
end subroutine compute_prandtl
!>===================================================================================

!>===================================================================================
subroutine get_wall_temperature(izone, iface, wall_temperature)
!>===================================================================================
  use cs_nz_condensation, only: iztag1d, ztpar
  use optcal, only: isuite, ntcabs
  use cs_nz_tagmr, only: ztpar0, ztmur
  use cstphy, only: tkelvi
  implicit none
  integer, intent(in) :: iface, izone
  double precision, intent(out) :: wall_temperature
  double precision :: t_wall
  if(iztag1d(izone).eq.1) then
    if(isuite.eq.0.and.ntcabs.eq.1) then
      t_wall = ztpar0(izone)
    else
      t_wall = ztmur(iface,1)
    endif
  else
    t_wall = ztpar(izone)
  endif

  wall_temperature = t_wall + tkelvi
end subroutine get_wall_temperature
!>===================================================================================

!>===================================================================================
subroutine compute_grashof(gravity, drho, length, kin_viscosity, grashof)
!>===================================================================================
  implicit none
  double precision, intent(in) :: gravity, drho, length, kin_viscosity
  double precision, intent(out) :: grashof
  grashof = gravity*dabs(drho)*length**3/(kin_viscosity**2)
end subroutine compute_grashof
!>===================================================================================

!>===================================================================================
subroutine get_temperature(enthalpy, cp, temperature)
!>===================================================================================
  implicit none
  double precision, intent(in) :: enthalpy, cp
  double precision, intent(out) :: temperature
  temperature = enthalpy / cp
end subroutine get_temperature
!>===================================================================================

subroutine compute_characteristic_length(point, reference, axis, length)
  implicit none
  !> inout
  double precision, intent(in) :: point(3), reference(3), axis(3)
  double precision, intent(out) :: length
  !> local
  integer :: idir
  length = 0.0d0
  do idir = 1, 3
    length = length + (point(idir) - reference(idir)) * axis(idir)
  enddo
end subroutine compute_characteristic_length

! TODO : replace by cs_math_3_normalize when switching from Fortran to C
subroutine normalize_vector(nb_dim, vector)
  use cstnum, only: epzero
  implicit none
  !> inout
  integer, intent(in) :: nb_dim
  double precision, intent(inout) :: vector(nb_dim)
  !> local
  integer :: i
  double precision :: norm
  norm = 0.0d0
  do i = 1, nb_dim
    norm = norm + vector(i)*vector(i)
  enddo
  norm = dsqrt(norm)
  if (norm > epzero) then
    do i = 1, nb_dim
      vector(i) = vector(i) / norm
    enddo
  endif
end subroutine normalize_vector

subroutine compute_tangential_velocity(velocity, normal, coeff, tangential_component)
  implicit none
  !> inout
  double precision, intent(in) :: velocity(3), normal(3), coeff
  double precision, intent(out) :: tangential_component
  !> local
  integer :: i
  double precision :: u_square, u_normal
  u_square = 0.0d0
  u_normal   = 0.0d0
  do i=1,3
    u_normal = u_normal + velocity(i) * normal(i) * coeff
    u_square = u_square + velocity(i) * velocity(i)
  enddo
  tangential_component = dsqrt(u_square - u_normal**(2.d0))
end subroutine compute_tangential_velocity

end module condensation_module
