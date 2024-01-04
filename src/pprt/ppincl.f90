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

!> \file ppincl.f90
!> General module for specific physics

module ppincl

  !===========================================================================

  use, intrinsic :: iso_c_binding

  use ppppar
  use ppthch

  !=============================================================================

  implicit none

  !===========================================================================

  !> \defgroup ppincl General module for specific physics

  !> \addtogroup ppincl
  !> \{

  !> \defgroup choice Indicator table for specific physics

  !> \addtogroup choice
  !> \{

  !----------------------------------------------------------------------------
  !--> TABLEAU INDICATEURS DU CHOIX DE LA PHYSIQUE PARTICULIERE CHOISIE

  !> number of specific physics
  integer   nmodmx
  parameter(nmodmx = 15)

  !> global indicator for speciphic physics
  !> By default, all the indicators ippmod(i.....) are initialized to -1,
  !> which means that no specific physics is activated.
  !>   - Diffusion flame in the framework of “3 points” rapid complete chemistry:
  !>  indicator ippmod(icod3p)
  !>      - ippmod(icod3p) = 0 adiabatic conditions
  !>      - ippmod(icod3p) = 1 permeatic conditions (enthalpy transport)
  !>      - ippmod(icod3p) =-1 module not activated
  !>
  !>   - Diffusion flames in the framework of steady laminar flamelet approach
  !>  indicator ippmod(islfm)
  !>      - ippmod(islfm) = 0 classic steady laminar flamelet model:
  !>                          adiabatic conditions
  !>      - ippmod(islfm) = 1 classic steady laminar flamelet model:
  !>                          non-adiabatic conditions (enthalpy transport with heat loss)
  !>      - ippmod(islfm) = 2 flamelet/progress variable model:
  !>                          adiabatic conditions
  !>      - ippmod(islfm) = 3 flamelet/progress variable model:
  !>                          non-adiabatic conditions (enthalpy transport with heat loss)
  !>      - ippmod(islfm) =-1 module not activated

  !>   - Eddy Break Up pre-mixed flame: indicator ippmod(icoebu)
  !>      - ippmod(icoebu) = 0 adiabatic conditions at constant richness
  !>      - ippmod(icoebu) = 1 permeatic conditions at constant richness
  !>      - ippmod(icoebu) = 2 adiabatic conditions at variable richness
  !>      - ippmod(icoebu) = 3 permeatic conditions at variable richness
  !>      - ippmod(icoebu) =-1 module not activated
  !>   - Libby-Williams pre-mixed flame: indicator ippmod(icolwc)
  !>      - ippmod(icolwc)=0 two peak model with adiabiatic conditions.
  !>      - ippmod(icolwc)=1 two peak model with permeatic conditions.
  !>      - ippmod(icolwc)=2 three peak model with adiabiatic conditions.
  !>      - ippmod(icolwc)=3 three peak model with permeatic conditions.
  !>      - ippmod(icolwc)=4 four peak model with adiabiatic conditions.
  !>      - ippmod(icolwc)=5 four peak model with permeatic conditions.
  !>      - ippmod(icolwc)=-1 module not activated.
  !>   - Multi-coals and multi-classes pulverised coal combustion: indicator ippmod(iccoal)
  !>
  !>  The number of different coals must be inferior or equal to ncharm = 5.
  !>  The number of particle size classes nclpch(icha) for the coal icha, must be
  !>  inferior or equal to ncpcmx = 20.
  !>      - ippmod(iccoal) = 0 imbalance between the temperature of the continuous and the
  !>        solid phases
  !>      - ippmod(iccoal) = 1 otherwise
  !>      - ippmod(iccoal) =-1 module not activated
  !>   - Electric arcs module (Joule effect and Laplace forces): indicator ippmod(ielarc)
  !>      - ippmod(ielarc) = 1 determination of the magnetic field by means of the Ampere’
  !>   theorem  (not available)
  !>      - ippmod(ielarc) = 2 determination of the magnetic field by means of the vector potential
  !>      - ippmod(ielarc) =-1 module not activated
  !>   - Joule effect module (Laplace forces not taken into account): indicator ippmod(ieljou)
  !>      - ippmod(ieljou) = 1 use of a real potential
  !>      - ippmod(ieljou) = 2 use of a complex potential
  !>      - ippmod(ieljou) = 3 use of real potential and specific boundary conditions
  !>  for transformers.
  !>      - ippmod(ieljou) = 4 use of complex potential and specific boundary conditions
  !>  for transformers.
  !>      - ippmod(ieljou) =-1 module not activated
  !>   - compressible flow module: indicator ippmod(icompf)
  !>      - ippmod(icompf) = 2 module activated: homogeneous two phase model
  !>      - ippmod(icompf) = 1 module activated: single phase model
  !>      - ippmod(icompf) = 0 module activated: single phase barotropic model
  !>      - ippmod(icompf) =-1 module not activated
  !>   - atmospheric flow module: indicator ippmod(iatmos)
  !>      - ippmod(iatmos) =-1 module not activated
  !>      - ippmod(iatmos) = 0 standard modelling
  !>      - ippmod(iatmos) = 1 dry atmosphere
  !>      - ippmod(iatmos) = 2 humid atmosphere

  integer(c_int), pointer, save :: ippmod(:)

  !> ippmod(iphpar) is a global indicator for the specific physics:
  !>  - 0: no specific physics
  !>  - 1: switch on the specific physics
  !>  - 2: switch on the specific physics plus radiative transfer
  !>       with a parametric file
  integer :: iphpar

  !> pointer for specific physics
  !> - ippmod(icod3p) = 0 adiabatic conditions
  !> - ippmod(icod3p) = 1 permeatic conditions (enthalpy transport)
  !> - ippmod(icod3p) =-1 module not activated
  integer ::  icod3p

  !> pointer to specify steady laminar flamelet approach
  !> - ippmod(islfm) = 0 classic steady laminar flamelet model:
  !>                     adiabatic conditions
  !> - ippmod(islfm) = 1 classic steady laminar flamelet model:
  !>                     non-adiabatic conditions (enthalpy transport with heat loss)
  !> - ippmod(islfm) = 2 flamelet/progress variable model:
  !>                     adiabatic conditions
  !> - ippmod(islfm) = 3 flamelet/progress variable model:
  !>                     non-adiabatic conditions (enthalpy transport with heat loss)
  !> - ippmod(islfm) =-1 module not activated
  integer :: islfm

  !> pointer to specify Eddy Break Up pre-mixed flame with indicator ippmod(icoebu)
  !> - ippmod(icoebu) = 0 adiabatic conditions at constant richness
  !> - ippmod(icoebu) = 1 permeatic conditions at constant richness
  !> - ippmod(icoebu) = 2 adiabatic conditions at variable richness
  !> - ippmod(icoebu) = 3 permeatic conditions at variable richness
  !> - ippmod(icoebu) =-1 module not activated
  integer ::  icoebu

  !> pointer to specify Libby-Williams pre-mixed flame withy indicator ippmod(icolwc)
  !> - ippmod(icolwc)=0 two peak model with adiabiatic conditions.
  !> - ippmod(icolwc)=1 two peak model with permeatic conditions.
  !> - ippmod(icolwc)=2 three peak model with adiabiatic conditions.
  !> - ippmod(icolwc)=3 three peak model with permeatic conditions.
  !> - ippmod(icolwc)=4 four peak model with adiabiatic conditions.
  !> - ippmod(icolwc)=5 four peak model with permeatic conditions.
  !> - ippmod(icolwc)=-1 module not activated.
  integer ::  icolwc

  ! TODO Modeles propres a la combustion gaz ICO...
  integer(c_int), pointer, save ::  isoot

  !> pointer to specify Joule effect module (Laplace forces not taken into account)
  !> with indicator ippmod(ieljou):
  !> - ippmod(ieljou) = 1 use of a real potential
  !> - ippmod(ieljou) = 2 use of a complex potential
  !> - ippmod(ieljou) = 3 use of real potential and specific boundary conditions
  !>  for transformers.
  !> - ippmod(ieljou) = 4 use of complex potential and specific boundary conditions
  !>  for transformers.
  !> - ippmod(ieljou) =-1 module not activated
  integer ::  ieljou

  !> pointer to specify Electric arcs module (Joule effect and Laplace forces)
  !> with indicator ippmod(ielarc):
  !> - ippmod(ielarc) = 1 determination of the magnetic field by means of the Ampere’
  !> theorem  (not available)
  !> - ippmod(ielarc) = 2 determination of the magnetic field by means of the vector potential
  !> - ippmod(ielarc) =-1 module not activated
  integer ::  ielarc

  !> pointer to specify multi-coals and multi-classes pulverised coal combustion
  !> with indicator ippmod(iccoal).
  !> The number of different coals must be inferior or equal to ncharm = 3.
  !> The number of particle size classes nclpch(icha) for the coal icha, must be
  !> inferior or equal to ncpcmx = 10.
  !> - ippmod(iccoal) = 0 imbalance between the temperature of the continuous and the
  !> solid phases
  !> - ippmod(iccoal) = 1 otherwise
  !> - ippmod(iccoal) =-1 module not activated
  integer ::  iccoal

  !> pointer to specify compressible module with indicator ippmod(icompf)
  !>      - ippmod(icompf) = 2 module activated: homogeneous two phase model
  !>      - ippmod(icompf) = 1 module activated: single phase model
  !>      - ippmod(icompf) = 0 module activated: single phase barotropic model
  !>      - ippmod(icompf) =-1 module not activated
  integer ::  icompf

  !> pointer to specify atmospheric flow module with indicator ippmod(iatmos)
  !> - ippmod(iatmos) =-1 module not activated
  !> - ippmod(iatmos) = 0 standard modelling
  !> - ippmod(iatmos) = 1 dry atmosphere
  !> - ippmod(iatmos) = 2 humid atmosphere
  integer ::  iatmos

  !> pointer to specify cooling towers module with indicator ippmod(iaeros)
  !> - ippmod(iaeros) =-1 module not activated
  !> - ippmod(iaeros) >= 0  activated
  integer ::  iaeros

  !> pointer to specify gas mixture module with indicator ippmod(igmix)
  !> - ippmod(igmix) =-1 module not activated
  !> - ippmod(igmix) = 0  Air/Helium   gas mixtures
  !> - ippmod(igmix) = 1  Air/Hydrogen gas mixtures
  !> - ippmod(igmix) = 2  Air/Steam    gas mixtures
  !> - ippmod(igmix) = 3  Air/Helium/Steam gas mixtures
  !> - ippmod(igmix) = 4  Air/Hydrogen/Steam gas mixtures
  !> - ippmod(igmix) = 5  Air/Helium with O2 from the air deduced

  integer ::  igmix

  parameter       (iphpar = 1,  icod3p = 2,  islfm = 3,             &
                   icoebu = 4,  icolwc = 5,  iccoal = 6,            &
                   ieljou = 7,  ielarc = 8,  icompf = 9,            &
                   iatmos = 10, iaeros = 11,                        &
                   igmix  = 12)

  !> \}

  !--> PARAMETERS ASSOCIATED WITH THE GAS MIXTURE MODELLING

  !> \defgroup modelling chosen by the user

  !> \addtogroup gas_mixture
  !> \{

  !> Specific condensation modelling
  !>      if = -1 module not activated
  !>      if =  0 condensation source terms activated
  integer(c_int), pointer, save ::  icondb

  !> Wall condensation correlation (only if icondb > -1)
  !>      if =  0 legacy copain model
  !>      if =  1 updated copain model (Benteboula and Dabbene 2020)
  !>      if =  2 Uchida model (TODO : add ref)
  !>      if =  3 Dehbi model (Dehbi 2015)
  integer(c_int), pointer, save ::  icondb_model

  !> Specific condensation modelling
  !>      if = -1 module not activated
  !>      if =  0 condensation source terms with metal
  !>                               structures activate
  integer(c_int), pointer, save ::  icondv

  !> \}


  !--> POINTEURS VARIABLES COMBUSTION GAZ

  !> \defgroup gaz_combustion Gaz combustion variables pointers

  !> \addtogroup gaz_combustion
  !> \{

  ! ---- Variables transportees

  !> pointer to specify the mixing rate in isca(ifm)
  integer, save :: ifm

  !> pointer to specify the variance of the mixing rate in isca(ifp2m)
  integer, save :: ifp2m

  !> pointer to specify the second moment of the mixing rate in isca(ifsqm):
  integer, save :: ifsqm

  !> pointer to specify the transported progress variable ippmod(islfm) >= 2:
  integer, save :: ipvm

  !> pointer to specify the fresh gas mass fraction in isca(iygfm)
  integer, save :: iygfm

  !> the intersection computation mode. If its value is:
  !> - 1 (default), the original algorithm is used. Care should be taken to clip
  !> the intersection on an extremity.
  !> - 2, a new intersection algorithm is used. Caution should be used to avoid to clip
  !> the intersection on an extremity.
  integer, save :: icm

  ! TODO
  !> transported variable
  integer, save :: icp2m

  ! TODO
  !> transported variable
  integer, save :: ifpcpm

  ! TODO
  !> transported variable
  integer, save :: iyfm
  ! TODO
  !> transported variable
  integer, save :: iyfp2m
  ! TODO
  !> transported variable
  integer, save :: icoyfp

  ! ---- Variables d'etat

  !> mass fractions :
  !>  - iym(1): is fuel mass fraction
  !>  - iym(2): oxidiser mass fraction
  !>  - iym(3): product mass fraction
  !> ibym() contains the matching field ids at boundary faces
  integer, save :: iym(ngazgm)
  integer, save :: ibym(ngazgm)

  !> state variable (temperature)
  integer, save :: itemp
  !> state variable
  integer, save :: ifmin
  !> state variable
  integer, save :: ifmax

  !> state variable: Pointer to the reconstructed variance in case of mode_fp2m = 1
  integer, save :: irecvr

  !> state variable: Pointer to the total scalar dissipation rate
  integer, save :: itotki

  !> state variable: Pointer to volumetric heat release rate
  integer, save :: ihrr

  !> state variable: Pointer to enthalpy defect
  integer, save :: ixr

  !> state variable: Pointer to enthalpy defect
  integer, save :: iomgc

  !> state variable: absorption coefficient, when the radiation modelling is activated
  integer, save :: ickabs

  !> state variable:  \f$T^2\f$ term
  integer, save :: it2m
  !> state variable:  \f$T^3\f$ term, when the radiation modelling is activated
  integer, save :: it3m
  !> state variable:  \f$T^4\f$ term, when the radiation modelling is activated
  integer, save :: it4m

  ! pointer for source term in combustion
  ! TODO
  integer, save :: itsc

  !> pointer for soot precursor number in isca (isoot = 1)
  integer, save :: inpm

  !> pointer for soot mass fraction in isca (isoot = 1)
  integer, save :: ifsm

  !> \}

  !--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE

  !> \defgroup coal_combustion  Pulverized coal combustion variables

  !> \addtogroup coal_combustion
  !> \{

  ! ---- Variables transportees
  !        Phase continue (melange gazeux)

  !> mean value of the tracer 1 representing the light
  !> volatiles released by the coal \c icha
  integer, save :: if1m(ncharm)

  !> mean value of the tracer 2 representing the heavy
  !> volatiles released by the coal \c icha
  integer, save :: if2m(ncharm)

  !> tracer 4: mass of the oxydant 2 divided by the mass of bulk
  integer, save :: if4m
  !> tracer 5: mass of the oxydant 3 divided by the mass of bulk
  integer, save :: if5m
  !> tracer 6: water coming from drying
  integer, save :: if6m
  !> tracer 7: mass of the carbon from coal oxydized by O2
  !> divided by the mass of bulk
  integer, save :: if7m
  !> tracer 8: mass of the carbon from coal gasified by CO2
  !> divided by the mass of bulk
  integer, save :: if8m
  !> tracer 9: mass of the Carbon from coal gasified by H2O
  !> divided by the mass of bulk
  integer, save :: if9m

  !> f1f2 variance
  integer, save :: ifvp2m

  !        Phase dispersee (classe de particules)
  !> coke mass fraction related to the class icla
  integer, save :: ixck(nclcpm)

  !> reactive coal mass fraction related to the class \c icla
  integer, save :: ixch(nclcpm)

  !> number of particles of the class \c icla per kg of air-coal mixture
  integer, save :: inp(nclcpm)

  !>  mass enthalpy of the coal of class \c icla, if we are in permeatic conditions
  integer, save :: ih2(nclcpm)

  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase (particle class)
  integer, save :: ixwt(nclcpm)

  !> Pointer to Np*age(particles)
  integer, save :: inagecp(nclcpm)

  !>
  integer, save :: iv_p_x(nclcpm)

  !>
  integer, save :: iv_p_y(nclcpm)

  !>
  integer, save :: iv_p_z(nclcpm)

  ! ---- Variables d'etat
  !        Phase continue (melange gazeux)

  !> mass fractions:
  !>  - iym1(1): mass fraction of \f$CH_{X1m}\f$ (light volatiles) in the gas mixture
  !>  - iym1(2): mass fraction of \f$CH_{X2m}\f$ (heavy volatiles) in the gas mixture
  !>  - iym1(3): mass fraction of CO in the gas mixture
  !>  - iym1(4): mass fraction of \f$O_2\f$ in the gas mixture
  !>  - iym1(5): mass fraction of \f$CO_2\f$ in the gas mixture
  !>  - iym1(6): mass fraction of \f$H_2O\f$ in the gas mixture
  !>  - iym1(7): mass fraction of \f$N_2\f$ in the gas mixture
  integer, save :: iym1(ngazem)

  ! TODO absent de la doc utilisateur
  !> State variables of continuous phase (gas mixture)
  integer, save :: irom1

  !>  molar mass of the gas mixture
  integer, save :: immel

  !        Phase dispersee (classes de particules)

  !> temperature of the particles of the class \c icla
  integer, save :: itemp2(nclcpm)

  !> density of the particles of the class \c icla
  integer, save :: irom2(nclcpm)

  !> diameter of the particles of the class \c icla
  integer, save :: idiam2(nclcpm)

  !>  solid mass fraction of the class \c icla
  integer, save :: ix2(nclcpm)

  !> disappearance rate of the reactive coal of the class \c icla
  integer, save :: igmdch(nclcpm)

  !> coke disappearance rate of the coke burnout of the class \c icla
  integer, save :: igmhet(nclcpm)

  !> Implicite part of the exchanges to the gas by molecular distribution
  integer, save :: igmtr(nclcpm)

  ! TODO absent de la doc utilisateur
  !> State variables of dispersed phase (particles class)
  integer, save :: ighco2(nclcpm)

  !>  mass transfer caused by the release of light volatiles  of the class \c icla
  integer, save :: igmdv1(nclcpm)

  !>  mass transfer caused by the release of heavy volatiles  of the class \c icla
  integer, save :: igmdv2(nclcpm)

  ! TODO absent de la doc utilisateur
  !> State variables of dispersed phase (particles class)
  integer, save :: igmsec(nclcpm)

  !> Used for bulk balance of Carbon
  integer, save :: ibcarbone
  !> Used for bulk balance of Oxygen
  integer, save :: iboxygen
  !> Used for bulk balance of Hydrogen
  integer, save :: ibhydrogen

  !> \}

  !--> POINTEURS VARIABLES COMBUSTION FUEL

  !> \defgroup fuel_combustion Fuel combustion variables

  !> \addtogroup fuel_combustion
  !> \{

  !        Phase dispersee
  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase
  integer, save :: ihlf(nclcpm)

  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase
  integer, save :: ixkf(nclcpm)

  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase
  integer, save :: ixfol(nclcpm)

  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase
  integer, save :: ing(nclcpm)

  ! ---- Variables d'etat
  !        Phase continue

  ! TODO absent de la doc utilisateur
  !> state variable of continuous phase
  integer, save :: iyfol(nclcpm)
  !        Phase dispersee

  ! TODO absent de la doc utilisateur
  !> state variable of dispersed phase
  integer, save :: ih1hlf(nclcpm)
  ! TODO absent de la doc utilisateur
  !> state variable of dispersed phase
  integer, save :: igmhtf(nclcpm)
  ! TODO absent de la doc utilisateur
  !> state variable of dispersed phase
  integer, save :: igmeva(nclcpm)

  !> \}


  !> \defgroup compressible Compressible models options

  !> \addtogroup compressible
  !> \{

  !> specific total energy for compressible algorithm
  integer, save :: ienerg

  !> temperature deduced from the specific total energy
  integer, save :: itempk

  !> \defgroup comp_homogeneous Homogeneous two-phase compressible model options

  !> \addtogroup comp_homogeneous
  !> \{

  !> \anchor ifracv
  !> homogeneous model, volume fraction \f$ \alpha \f$
  integer, save :: ifracv

  !> \anchor ifracm
  !> homogeneous model, mass fraction \f$ y \f$
  integer, save :: ifracm

  !> \anchor ifrace
  !> homogeneous model, energy fraction \f$ z \f$
  integer, save :: ifrace

  !> \}

  !> \defgroup common Common

  !> \addtogroup common
  !> \{

  ! ---- Aliases pour les conditions aux limites

  !> alias for boundary conditions
  integer, save :: irun

  !> alias for boundary conditions
  integer, save :: irunh

  !> reference volume viscosity
  real(c_double), pointer, save :: viscv0

  !> \}

  !> \defgroup enthalpy Enthalpic variables pointers

  !> \addtogroup enthalpy
  !> \{

  !> enthalpy, if transported or if deduced
  integer, save :: ihm

  !> \anchor srrom
  !> with gas combustion, or pulverised coal, \ref srrom
  !> is the sub-relaxation coefficient for the density, following the formula:
  !> \f$\rho^{n+1}$\,=\,srrom\,$\rho^n$+(1-srrom)\,$\rho^{n+1}\f$
  !> hence, with a zero value, there is no sub-relaxation.
  !> \ref srrom is initialized to \ref cstnum::grand "-grand" and the user must
  !> specify a proper value through the GUI or the initialization subroutine
  !> (\ref cs_user_combustion).
  !> It is automatically used after the second time-step.
  !>
  !> Always useful with gas combustion or pulverized coal.
  double precision, save :: srrom

  !> \}


  !> \defgroup boundary_conditions Boundary conditions

  !> \addtogroup boundary_conditions
  !> \{

  !> imposed flow zone indicator
  !> in a way which is similar to the process described in the framework of the EBU module,
  !> the user chooses for every inlet face to impose the mass flow or not
  !> (\ref iqimp "iqimp"(izone)=1 or 0). If the mass flow is imposed, the user
  !> must set the air mass flow value \ref coincl::qimp "qimp"(izone) ant its direction
  !> in \ref rcodcl "rcodcl"(ifac,\ref iu), \ref rcodcl "rcodcl"(ifac,\ref iv)
  !> and \ref rcodcl "rcodcl"(ifac,\ref iw).
  !> If the velocity is imposed, he has to set  \ref rcodcl "rcodcl"(ifac,\ref iu),
  !> \ref rcodcl "rcodcl"(ifac,\ref iv), and \ref rcodcl "rcodcl"(ifac,\ref iw).
  integer(c_int), pointer, save :: iqimp(:)

  !> condition type turbulence indicator
  !>  - 0 : given by the user
  !>  - 1 : automatic, from hydraulic diameter and input velocity performed.
  !>  - 2 : automatic, from turbulent intensity and input velocity performed.
  integer(c_int), pointer, save :: icalke(:)

  !> turbulent intensity (k=1.5(uref*xintur)**2)
  real(c_double), pointer, save :: xintur(:)

  !> hydraulic diameter
  real(c_double), pointer, save :: dh(:)

  !> index of maximum reached boundary zone
  integer, save :: nozapm

  !> number of boundary zones on current process
  integer, save :: nzfppp

  !> list of boundary zones index
  integer, save :: ilzppp(nozppm)

  !> \}

  !=============================================================================

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure

    subroutine cs_f_fluid_properties_pp_get_pointers(viscv0)   &
      bind(C, name='cs_f_fluid_properties_pp_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: viscv0
    end subroutine cs_f_fluid_properties_pp_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_physical_model_get_pointers(p_ippmod)    &
      bind(C, name='cs_f_physical_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ippmod
    end subroutine cs_f_physical_model_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_combustion_model_get_pointers(p_isoot)    &
      bind(C, name='cs_f_combustion_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_isoot
    end subroutine cs_f_combustion_model_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_wall_condensation_get_model_pointers(p_icondb, &
                                                         p_icondv, &
                                                         p_icondb_model) &
      bind(C, name='cs_f_wall_condensation_get_model_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_icondb, p_icondv, p_icondb_model
    end subroutine cs_f_wall_condensation_get_model_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving BC zone array pointers

    subroutine cs_f_boundary_conditions_get_ppincl_pointers(p_iqimp, p_icalke,  &
                                                            p_xintur, p_dh)     &
      bind(C, name='cs_f_boundary_conditions_get_ppincl_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_iqimp, p_icalke, p_xintur, p_dh
    end subroutine cs_f_boundary_conditions_get_ppincl_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran physical models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine pp_models_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_viscv0, p_ippmod, p_isoot
    type(c_ptr) :: p_icondb, p_icondv, p_icondb_model

    call cs_f_fluid_properties_pp_get_pointers(p_viscv0)
    call c_f_pointer(p_viscv0, viscv0)

    call cs_f_physical_model_get_pointers(p_ippmod)
    call cs_f_combustion_model_get_pointers(p_isoot)
    call cs_f_wall_condensation_get_model_pointers(p_icondb,  &
                                                   p_icondv,  &
                                                   p_icondb_model)

    call c_f_pointer(p_ippmod, ippmod, [nmodmx])
    call c_f_pointer(p_isoot, isoot)
    call c_f_pointer(p_icondb, icondb)
    call c_f_pointer(p_icondv, icondv)
    call c_f_pointer(p_icondb_model, icondb_model)

  end subroutine pp_models_init

  !=============================================================================

  !> \brief Map Fortran physical models boundary condition info.
  !> This maps Fortran pointers to global C variables.

  subroutine pp_models_bc_map

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_iqimp, p_icalke, p_xintur, p_dh

    call cs_f_boundary_conditions_get_ppincl_pointers(p_iqimp, p_icalke,  &
                                                      p_xintur, p_dh)

    call c_f_pointer(p_iqimp, iqimp, [nozppm])
    call c_f_pointer(p_icalke, icalke, [nozppm])
    call c_f_pointer(p_xintur, xintur, [nozppm])
    call c_f_pointer(p_dh, dh, [nozppm])

  end subroutine pp_models_bc_map

  !=============================================================================

end module ppincl
