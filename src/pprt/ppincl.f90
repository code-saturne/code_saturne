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
  parameter(nmodmx = 14)

  !> global indicator for speciphic physics
  !> By default, all the indicators ippmod(i.....) are initialized to -1,
  !> which means that no specific physics is activated.
  !>   - Diffusion flame in the framework of “3 points” rapid complete chemistry:
  !>  indicator ippmod(icod3p)
  !>      - ippmod(icod3p) = 0 adiabatic conditions
  !>      - ippmod(icod3p) = 1 permeatic conditions (enthalpy transport)
  !>      - ippmod(icod3p) =-1 module not activated
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
  !>  The number of different coals must be inferior or equal to ncharm = 3.
  !>  The number of particle size classes nclpch(icha) for the coal icha, must be
  !>  inferior or equal to ncpcmx = 10.
  !>      - ippmod(iccoal) = 0 imbalance between the temperature of the continuous and the
  !>        solid phases
  !>      - ippmod(iccoal) = 1 otherwise
  !>      - ippmod(iccoal) =-1 module not activated
  !>   - Multi-classes pulverised heavy fuel combustion: indicator ippmod(icfuel)
  !>      - ippmod(icfuel) = 0 module activated
  !>      - ippmod(icfuel) =-1 module not activated
  !>   - Lagrangian modelling of multi-coals and multi-classes pulverised coal combustion:
  !>     indicator ippmod(icpl3c) The number of different coals must be inferior or equal
  !>     to ncharm = 3. The number of particle size classes nclpch(icha) for the coal icha,
  !>     must be inferior or equal to ncpcmx = 10.
  !>      - ippmod(icpl3c) = 1 coupling with the Lagrangian module, with transport of H2
  !>      - ippmod(icpl3c) =-1 module not activated
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

  !> pointer to specify Lagrangian modelling of multi-coals and multi-classes pulverised coal
  !> combustion with indicator ippmod(icpl3c).
  !> The number of different coals must be inferior or equal
  !>  to ncharm = 3. The number of particle size classes nclpch(icha) for the coal icha,
  !>  must be inferior or equal to ncpcmx = 10.
  !>  - ippmod(icpl3c) = 1 coupling with the Lagrangian module, with transport of H2
  !>  - ippmod(icpl3c) =-1 module not activated
  integer ::  icpl3c

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

  !> coal with drift (0: without drift (default), 1: with)
  integer ::  i_comb_drift

  !> pointer to specify multi-classes pulverised heavy fuel combustion
  !> with indicator ippmod(icfuel)
  !> - ippmod(icfuel) = 0 module activated
  !> - ippmod(icfuel) =-1 module not activated
  integer ::  icfuel

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

  integer ::  igmix

  !> pointer to specify richards model
  !> - ippmod(iricha) =-1 module not activated
  !> - ippmod(iricha) = 1 module activated
  integer ::  idarcy

  parameter       (iphpar = 1 , icod3p = 2 ,                        &
                   icoebu = 3 , icolwc = 4 ,                        &
                   icpl3c = 5 , iccoal = 6 , icfuel = 7 ,           &
                   ieljou = 8 , ielarc = 9 , icompf = 10,           &
                   iatmos = 11, iaeros = 12,                        &
                   igmix  = 13, idarcy = 14)

  !> \}

  !--> PARAMETERS ASSOCIATED WITH THE GAS MIXTURE MODELLING

  !> \defgroup modelling chosen by the user

  !> \addtogroup gas_mixture
  !> \{

  !> Specific condensation modelling
  !>      if = -1 module not activated
  !>      if =  0 condensation source terms activated
  integer, save ::  icondb

  !> Specific condensation modelling
  !>      if = -1 module not activated
  !>      if =  0 condensation source terms with metal
  !>                               structures activate
  integer, save ::  icondv

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

  !> state variable: absorption coefficient, when the radiation modelling is activated
  integer, save :: ickabs

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

  !> mean value of the tracer 3 representing
  !> the carbon released as CO during coke burnout
  integer, save :: if3m

  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if4m
  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if5m
  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if6m
  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if7m
  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if8m
  ! TODO absent de la doc utilisateur
  !> transported variable of countinuous phase (gas mixture)
  integer, save :: if9m

  !> the variance associated with the tracer 4 representing the air
  !> (the mean value of this tracer is not transported, it can be deduced directly
  !> from the three others)
  integer, save :: if4p2m

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

  !> Pointer to age of bulk
  integer, save :: iage

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

  !> temperature of the gas mixture
  integer, save :: itemp1

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

  ! ---- Variables transportees
  !        Phase continue
  ! TODO absent de la doc utilisateur
  !> transported variable of continuous phase
  integer, save :: ifvap

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


  !> \defgroup compressible Compressible module

  !> \addtogroup compressible
  !> \{

  !> specific total energy for compressible algorithm
  integer, save :: ienerg

  !> temperature deduced from the specific total energy
  integer, save :: itempk

  !> \defgroup comp_homogeneous Homogeneous two-phase flow model

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

  !> \defgroup comp_properties Physical properties

  !> \addtogroup comp_properties
  !> \{

  !> additional property:
  !>  - 0 indicates that the volume viscosity is constant and equal to
  !> the reference volume viscosity \ref viscv0.
  !>  - 1 indicates that the volume viscosity is variable: its
  !> variation law must be specified in the user subroutine \ref usphyv.
  !>
  !> Always useful.
  !> The volume viscosity \f$\kappa\f$ is defined by the formula expressing the stress:
  !> \f{eqnarray*}{
  !> \tens{\sigma} = -P\,\tens{Id} + \mu (\grad\,\vect{u} +
  !> \ ^{t}\ggrad\,\vect{u})
  !> +(\kappa-\frac{2}{3}\mu)\,\dive(\vect{u})\,\tens{Id}
  !> \f}
  !>
  integer, save :: iviscv

  !> \}
  !> \}


  !> \defgroup common Common

  !> \addtogroup common
  !> \{

  ! ---- Aliases pour les conditions aux limites

  !> alias for boundary conditions
  integer, save :: irun

  !> alias for boundary conditions
  integer, save :: irunh

  !> reference volume viscosity (noted \f$\kappa\f$ in the equation
  !> expressing \f$\tens{\sigma}\f$ in the paragraph dedicated to \ref iviscv)
  !> always useful, it is the used value, unless the user specifies the volume
  !> viscosity in the user subroutine \ref usphyv
  double precision, save :: viscv0

  !> pressure predicion by an evolution equation
  integer, save :: ippred

  !> indicates whether the pressure should be updated (=1) or not (=0) after the
  !> solution of the acoustic equation always usef
  integer, save :: igrdpp

  ! --- Conditions aux limites prenant en compte l'equilibre hydrostatique

  !> indicates if the boundary conditions should take into account (=1)
  !> or not (=0) the hydrostatic balance.
  !>
  !> Always useful.
  !>
  !> In the cases where gravity is predominant, taking into account
  !> the hydrostatic pressure allows to get rid of the disturbances which
  !> may appear near the horizontal walls when the flow is little convective.
  !>
  !> Otherwise, when \ref icfgrp=0, the pressure condition is calculated
  !> from the solution of the unidimensional Euler equations for a perfect
  !> gas near a wall, for the variables "normal velocity", "density" and
  !> "pressure":
  !>
  !> Case of an expansion (M <= 0):
  !> \f{align*}{
  !>    P_p &= 0 & \textrm{if } 1 + \displaystyle\frac{\gamma-1}{2}M<0
  !> \\ P_p &= P_i \left(1 + \displaystyle\frac{\gamma-1}{2}M\right)
  !> ^{\frac{2\gamma}{\gamma-1}} & \textrm{otherwise}
  !> \f}
  !>
  !> Case of a schock (M > 0):
  !> \f{eqnarray*}{
  !> P_p = P_i \left(1 + \displaystyle\frac{\gamma(\gamma+1)}{4}M^2
  !> +\gamma M \displaystyle\sqrt{1+\displaystyle\frac{(\gamma+1)^2}{16}M^2}\right)
  !> \f}
  !>
  !> with \f$M = \displaystyle\frac{\vect{u}_i \cdot \vect{n}}{c_i}\f$,
  !> internal Mach number calculated with the variables taken in the cell.
  integer, save :: icfgrp

  !> \}


  !> \defgroup cooling Cooling Towers Pointers

  !> \addtogroup cooling
  !> \{

  !> Molar mass ratio = Ratio of the molar mass (H2O) over the molar mass (air)
  double precision :: molmass_rat
  parameter(molmass_rat = 0.622d0)

  !> \defgroup cool_transported  Transported variables

  !> \addtogroup cool_transported
  !> \{

  !> Pointer to define the mass fraction of the water in array isca(iymw)
  integer, save :: iymw

  !> Pointer to define mass fraction of liquid water which is injected in the packing
  integer, save :: iyml

  !> Pointer to define mass fraction of liquid water which is injected in the rain zones
  integer, save :: iy_p_l

  !> Pointer to define the temperature of liquide water in rain zone (in fact Y_l.T_l)
  integer, save :: it_p_l

  !> Pointer to define the temperature (property, deduced from enthalpy)
  !> of liquid water which is injected in the packing
  integer, save :: itml

  !> Pointer to define the field enthalpy of liquid water
  integer, save :: ihml

  !> Pointer to define the air humidity
  integer, save :: ihumid

  !> Pointer to define the vertical velocity of liquid water
  integer, save :: ivertvel

  !> \}
  !> \}

  !> \defgroup enthalpy Enthalpic variables pointers

  !> \addtogroup enthalpy
  !> \{

  !> enthalpy, if transported or if deduced
  integer, save :: ihm

  !> \anchor srrom
  !> with gas combustion, pulverised coal or the electric module, \ref srrom
  !> is the sub-relaxation coefficient for the density, following the formula:
  !> \f$\rho^{n+1}$\,=\,srrom\,$\rho^n$+(1-srrom)\,$\rho^{n+1}\f$
  !> hence, with a zero value, there is no sub-relaxation.
  !> With combustion and pulverized coal, \ref srrom is initialized to
  !> \ref cstnum::grand "-grand" and the user must specify a proper value through
  !> the Interface or the initialization subroutine (\ref cs_user_combustion).
  !> With gas combustion, pulverised coal or electric arcs, \ref srrom is
  !> automatically used after the second time-step. With Joule effect,
  !> the user decides whether or not it will be used in \ref cs_user_physical_properties
  !> from the coding law giving the density.
  !>
  !> Always useful with gas combustion, pulverized coal or the electric module.
  double precision, save :: srrom

  !> \}


  !> \defgroup boundary_conditions Boundary conditions

  !> \addtogroup boundary_conditions
  !> \{

  !> imposed flow zone indicator
  !> in a way which is similar to the process described in the framework of the EBU module,
  !> the user chooses for every inlet face to impose the mass flow or not
  !> (\ref iqimp "iqimp"(izone)=1 or 0). If the mass flow is imposed, the user
  !> must set the air mass flow value \ref cpincl::qimpat "qimpat"(izone), its direction in
  !> \ref rcodcl "rcodcl"(ifac,\ref iu), \ref rcodcl "rcodcl"(ifac,\ref iv)
  !> and \ref rcodcl "rcodcl"(ifac,\ref iw) and the incoming
  !> air temperature \ref cpincl::timpat "timpat"(izone) in Kelvin.
  !> If the velocity is imposed, he has to set  \ref rcodcl "rcodcl"(ifac,\ref iu),
  !> \ref rcodcl "rcodcl"(ifac,\ref iv), and \ref rcodcl "rcodcl"(ifac,\ref iw).
  integer, save ::          iqimp(nozppm)

  !> condition type turbulence indicator
  !>  - 0 : given by the user
  !>  - 1 : automatic, from hydraulic diameter and input velocity performed.
  !>  - 2 : automatic, from turbulent intensity and input velocity performed.
  integer, save ::          icalke(nozppm)

  !> turbulent intensity (k=1.5(uref*xintur)**2)
  double precision, save :: xintur(nozppm)

  !> hydraulic diameter
  double precision, save :: dh(nozppm)

  !> index of maximum reached boundary zone
  integer, save :: nozapm

  !> number of boundary zones on current process
  integer, save :: nzfppp

  !> list of boundary zones index
  integer, save :: ilzppp(nbzppm)

  !> \}

  !=============================================================================

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

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

    type(c_ptr) :: p_ippmod, p_isoot

    call cs_f_physical_model_get_pointers(p_ippmod)
    call cs_f_combustion_model_get_pointers(p_isoot)

    call c_f_pointer(p_ippmod, ippmod, [nmodmx])
    call c_f_pointer(p_isoot, isoot)

  end subroutine pp_models_init

  !=============================================================================

end module ppincl
