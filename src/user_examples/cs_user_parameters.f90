!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

!> \brief User subroutine for selection of specific physics module

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if the XML file from the GUI is used
!>                              (1 : yes, 0 : no)
!______________________________________________________________________________!

subroutine usppmo &
 ( ixmlpu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

!< [usppmo]

!===============================================================================
! 1.  Choice for a specific physics
!===============================================================================

! --- cod3p: Diffusion flame with complete fast chemistry (3 points)
! ==========

!        if = -1   module not activated
!        if =  0   adiabatic model
!        if =  1   extended model with enthalpy source term

if (ixmlpu.eq.0) then

  ippmod(icod3p) = -1

endif

! --- coebu: Eddy-Break Up pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference Spalding model
!                   (adiabatic, homogeneous mixture fraction)
!        if =  1   extended model with enthalpy source term
!                   (homogeneous mixture fraction : perfect premix)
!        if =  2   extended model with mixture fraction transport
!                   (adiabatic, no variance of mixture fraction)
!        if =  3   extended model with enthalpy and mixture fraction transport
!                   (dilution, thermal losses, etc.)

if (ixmlpu.eq.0) then

  ippmod(icoebu) = -1

endif

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference two-peak model with adiabatic condition
!        if =  1   extended two-peak model with enthapy source terms
!        if =  2   extended three-peak model, adiabatic
!        if =  3   extended three-peak model with enthalpy source terms
!        if =  4   extended four-peak model, adiabatic
!        if =  5   extended four-peak model with enthalpy source terms

if (ixmlpu.eq.0) then

  ippmod(icolwc) = -1

endif


! --- Soot model
! =================

!        if = -1   module not activated
!        if =  0   constant fraction of fuel Xsoot
!        if =  1   2 equations model of Moss et al.

isoot = 0

xsoot  = 0.1d0 ! ( if isoot = 0 )
rosoot = 2000.d0 ! kg/m3

! --- cfuel: Heavy fuel oil combustion
! ==========

!        Progressive evaporation (temperature gap)
!        Char residue
!        Sulphur tracking

!        if = -1   module not activated
!        if = 0    module activated

if (ixmlpu.eq.0) then

  ippmod(icfuel) = -1

endif

! --- coal :
! ==========
!
!     Pulverized coal combustion
!        Description of granulometry
!        Assumption of diffusion flame around particles
!         (extension of 3-point fast chemistry "D3P")
!        Between a mixture of gaseous fuels (volatiles matters, CO from char
!                                            oxydation)
!            and a mixture of oxidisers (air and water vapor)
!        Enthalpy for both mix and solid phase are solved
!
!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying

if (ixmlpu.eq.0) then

  ippmod(iccoal) = -1

endif

! Activate the drift: 0 (no activation),
!                     1 (transported particle velocity)
!                     2 (limit drop particle velocity)

i_comb_drift = 1


! --- cpl3c: Pulverized coal with Lagrangian reciprocal approach
! ==========

!        Not recently tested... at least outdated, may be obsolete

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying (NOT functional)

if (ixmlpu.eq.0) then

  ippmod(icpl3c) = -1

endif

! --- compf: Compressible flows
! ==========

!        if = -1   module not activated
!        if =  0   module activated
!        if =  1   barotropic version
!        if =  2   homogeneous two phase model

if (ixmlpu.eq.0) then

  ippmod(icompf) = -1

endif

! --- eljou: Joule effect
! ==========

!        if = -1   module not activated
!        if = 1    Potentiel reel
!        if = 2    Potentiel complexe
!        if = 3    Potentiel reel     + CDL Transfo
!        if = 4    Potentiel complexe + CDL Transfo

if (ixmlpu.eq.0) then

  ippmod(ieljou) = -1

endif

! --- elarc: Electric arcs
! ==========

!        if = -1   module not activated
!        if = 1    electric potential
!        if = 2    electric potential and vector potential (hence 3D modelling)

if (ixmlpu.eq.0) then

  ippmod(ielarc) = -1

endif

! --- atmos: Atmospheric flows
! ==========

!        if = -1   module not activated
!        if = 0    standard modelling
!        if = 1    dry atmosphere
!        if = 2    humid atmosphere (experimental)

if (ixmlpu.eq.0) then

  ippmod(iatmos) = -1

endif

! --- aeros: Cooling towers
! ==========

!        if = -1   module not activated
!        if = 0    no model (NOT functional)
!        if = 1    Poppe's model
!        if = 2    Merkel's model

if (ixmlpu.eq.0) then

  ippmod(iaeros) = -1

endif

! --- igmix: Gas mixtures modelling
! ==========
!        if =-1 module not activated
!        if = 0  Air/Helium   gas mixtures
!        if = 1  Air/Hydrogen gas mixtures
!        if = 2  Air/Steam    gas mixtures
!        if = 3  Air/Helium/Steam gas mixtures
!        if = 4  Air/Hydrogen/Steam gas mixtures


ippmod(igmix) = 0


! Radiative transfer module (iirayo)
!--------------------------
!        if = 0: not activated (Default)
!        if = 1: DOM
!        if = 2: approximation P1 method

iirayo = 1

! --- richards model
! ==========

!        if = -1   module not activated
!        if =  1   module activated

ippmod(idarcy) = -1

!===============================================================================
! 2.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)

!       if = 0   user-specified
!       if = 1   tabulated by JANAF (default)

if (ixmlpu.eq.0) then

  indjon = 1

endif

!===============================================================================
! 2.  Data file related to modules above
!===============================================================================

if (ixmlpu.eq.0) then

  ! Combustion

  if (     ippmod(icod3p).ge.0                                          &
      .or. ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0) then

    if (indjon.eq.1) then
      ficfpp = 'dp_C3P'
    else
      ficfpp = 'dp_C3PSJ'
    endif

  endif

  ! Fuel combustion

  if (ippmod(icfuel).ge.0) then
    ficfpp = 'dp_FUE'
  endif

  ! Atmospheric flows

  if (ippmod(iatmos).ge.0) then
    ficmet = 'meteo'
  endif

 if (ippmod(igmix).ge.0) then
   ! Specific condensation modelling

   ! wall condensation
   !      if = -1 module not activated
   !      if =  0 condensation source terms activated
   icondb = -1

   ! internal condensation
   !      if = -1 module not activated
   !      if =  0 condensation source terms with metal
   !                               structures activate
   icondv = -1
 endif

endif

!< [usppmo]

!----
! End
!----

return
end subroutine usppmo


!===============================================================================

!> \brief User subroutine for input of model selection parameters.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]      ixmlpu       indicates if the XML file from the GUI is used
!>                              used (1: yes, 0: no
!> \param[in, out] iturb        turbulence model
!> \param[in, out] itherm       thermal model
!> \param[in, out] iale         ale module
!> \param[in, out] ivofmt       vof method
!> \param[in, out] icavit       cavitation model
!______________________________________________________________________________!

subroutine usipph &
 ( ixmlpu, iturb , itherm, iale , ivofmt, icavit )

!===============================================================================
! Module files
!===============================================================================

use entsor, only: nfecra ! No other module should appear here
use optcal, only: irijco ! No other module should appear here

!===============================================================================

implicit none

! Arguments

integer ixmlpu
integer iturb, itherm, iale, ivofmt, icavit

! Local variables

!===============================================================================

!>    In this subroutine, only the parameters which already appear may
!>    be set, to the exclusion of any other.
!>
!>    If we are not using the Code_Saturne GUI:
!>    All the parameters which appear in this subroutine must be set.
!>
!>    If we are using the Code_Saturne GUI:
!>    parameters protected by a test of the form:
!>
!>      if (ixmlpu.eq.0) then
!>         ...
!>      endif
!>
!>    should already have been defined using the GUI, so only
!>    experts should consider removing the test and adapting them here.

!===============================================================================

!< [usipph]

! --- Cavitation module
!    - -1: module not activated
!    -  0: no vaporization/condensation model
!    -  1: Merkle's model
!
!  Specific cavitation module input parameters should be set usipsu
!  (see example in cs_user_parameters-cavitation.f90)
!

icavit = -1

! --- Enable two phase homogeneous model
!     (compressible module should be enabled first)


icfhgn = 1


!< [usipph]

!----
! Formats
!----


return
end subroutine usipph


!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

logical       inoprv
integer       ii, jj, ivar, kscmin, kscmax, keydri, kbfid, kccmin, kccmax
integer       f_id, idim1, itycat, ityloc, iscdri, iscal, ifcvsl, b_f_id

type(var_cal_opt) :: vcopt

!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================

!< [usipsu]

! Calculation options (optcal)
! ============================

! In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! By default, this file is read, but it may be useful to deactivate
! its use when restarting after a preprocessing stage possibly leading
! to a different number of faces (such as simply joining meshes on
! a different architecture or optimization level or with different options).

! Writing of auxiliary restart files may also be deactivated using: iecaux = 0

ileaux = 0


! --- Time stepping  (0 : uniform and constant
!                     1 : variable in time, uniform in space
!                     2 : variable in time and space
!                    -1 : steady algorithm)

idtvar = 0


! --- Duration
!       ntmabs = absolute number of the last time step required
!         if we have already run 10 time steps and want to
!         run 10 more, ntmabs must be set to 10 + 10 = 20

ntmabs = 10


! --- Reference time step
!     The example given below is probably not adapted to your case.


dtref  = 0.01d0

! --- Maximum time step: dtmax
!     Set a value base on characteristic values of your case.
!      otherwise, the code will use a multiple of dtref by default.
!     Example with
!        Ld: "dynamic" length (for example, the domain length)
!        Ud: characteristic flow velocity
!        Lt: thermal length (for example, the domain height gravity-wise)
!        Delta_rho/rho: relative density difference
!        g: gravity acceleration

!     dtmax = min(Ld/Ud, sqrt(Lt/(g.Delta_rho/rho)))


! --- Algorithm to take into account the thermodynamical pressure variation in time
!     (not used by default except if idilat = 3)

!     by default:
!     ----------
!      - the thermodynamic pressure (pther) is initialized with p0 = p_atmos
!      - the maximum thermodynamic pressure (pthermax) is initialized with -1
!        (no maximum by default, this term is used to model a venting effect when
!         a positive value is given by the user)
!      - a global leak can be set through a leakage surface sleak with a head
!      loss kleak of 2.9 (Idelcick)

ipthrm = 0

pthermax= -1.d0

sleak = 0.d0
kleak = 2.9d0


! --- Temperature or enthalpy

!   When used without specific physics, if we have chosen to solve in temperature
!     (that is if itherm = 1), the fluid temperature is considered to be in
!     degrees Kelvin by default (be careful for boundary conditions an expression
!     of physical properties depending on temperature)t.

!     If we wish for the fluid solver to work with a temperature in degrees Celsius,
!     we must set itpscl = 2.

!     This is recommended for Syrthes Coupling, but not recommended for the
!     radiative model, as it is a source of user errors in this case:
!     Indeed, the boundary conditions for the fluid temperature will then be
!     in degrees Celsius, while the boundary conditions for radiation in
!     cs_user_radiative_transfer_bcs must still be in Kelvin.

if (nmodpp.eq.0) then
  itpscl = 2
endif


!   If a USER scalar behaves like a temperature (relative to Cp):
!     we set iscacp(isca) = 1.
!
!   Otherwise, we do not modify iscacp(isca)


if (nscaus.gt.0) then
  do ii = 1, nscaus
    iscacp(isca(ii)) = 1
  enddo
endif

! --- Calculation (restart) with frozen velocity field (1 yes, 0 no)


iccvfg = 1


! --- Vortex method for inlet conditions in L.E.S.
!       (0: not activated,  1: activated)
!     The vortex method only regards the L.E.S. models
!     To use the vortex method, edit the 'usvort.f90' user file.


if (itytur.eq.4) then
  ivrtex = 1
endif

! --- Convective scheme

!     blencv = 0 for upwind (order 1 in space, "stable but diffusive")
!            = 1 for centered/second order (order 2 in space)
!       we may use intermediate real values.
!       Here we choose:
!         for the velocity and user scalars:
!           an upwind-centered scheme with 100% centering (blencv=1)
!         for other variables
!           the default code value (upwind standard, centered in LES)

!     Specifically, for user scalars
!       if we suspect an excessive level of numerical diffusion on
!         a variable ivar representing a user scalar
!         iscal (with ivar=isca(iscal)), it may be useful to set
!         blencv = 1.0d0 to use a second-order scheme in space for
!         convection. For temperature or enthalpy in particular, we
!         may thus choose in this case:
!
!         call field_get_key_struct_var_cal_opt(ivarfl(isca(iscalt)), vcopt)
!         vcopt%blencv = 1.0d0
!         call field_set_key_struct_var_cal_opt(ivarfl(isca(iscalt)), vcopt)

!       For non-user scalars relative to specific physics (coal, combustion,
!         electric arcs: see usppmo) implicitly defined by the model,
!         the corresponding information is set automatically elsewhere:
!         we do not modify blencv here.

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
vcopt%blencv = 1.0d0
call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)

if (nscaus.ge.1) then
  do ii = 1, nscaus
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    vcopt%blencv = 1.0d0
    call field_set_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
  enddo
endif


! --- Linear solver parameters (for each unknown)

!     epsilo: relative precision for the solution of the linear system.


if (nscaus.ge.1) then
  do ii = 1, nscaus
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    vcopt%epsilo = 1.d-6
    call field_set_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
  enddo
endif


! --- Dynamic reconstruction sweeps to handle non-orthogonlaities
!     This parameter computes automatically a dynamic relax factor,
!     and can be activated for any variable.
!      - iswdyn = 1: means that the last increment is relaxed
!      - iswdyn = 2: means that the last two increments are used to
!                         relax
!     NB: when iswdyn is greater than 1, then the number of
!         non-orthogonality sweeps is increased to 20.

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
vcopt%iswdyn = 1
call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! --- Rotation/curvature correction for eddy-viscosity turbulence models
!      0: deactivated
!      1: activated


irccor = 1


! --- Stabilization in turbulent regime

!     For difficult cases, a stabilization may be obtained by not
!     reconstructing the convective and diffusive flux for variables
!     of the turbulence model, that is
!       in k-epsilon: if (itytur.eq.2) then
!          ircflu(ik)   = 0 and ircflu(iep)  = 0
!       in Rij-epsilon: if (itytur.eq.3) then
!          ircflu(ir11) = 0,    ircflu(ir22) = 0,
!          ircflu(ir33) = 0,
!          ircflu(ir12) = 0,    ircflu(ir23) = 0,
!          ircflu(ir23) = 0,
!                                  and ircflu(iep)  = 0
!     (note that variable itytur is equal to iturb/10)


if (itytur.eq.2) then
  call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  vcopt%ircflu = 0
  call field_set_key_struct_var_cal_opt(ivarfl(ik), vcopt)

  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt)
  vcopt%ircflu = 0
  call field_set_key_struct_var_cal_opt(ivarfl(iep), vcopt)
endif


! --- Advanced re-initialization for EBRSM or k-omega models

!     - 0: switch off (default)
!     - 1: switch on

reinit_turb = 1


! --- Turbulent diffusion model for second moment closure (iturb = 3x)
!      0: scalar diffusivity (Shir model)
!      1: tensorial diffusivity (Daly and Harlow model, default model)


if (itytur.eq.3) then
  idirsm = 1
endif


! --- Advanced choice of Wall function


iwallf = 5


! Physical constants (cstphy)
! ===========================

! --- gravity (g in m/s2, with the sign in the calculation coordinate axes).

gx = 0.d0
gy = 0.d0
gz = 0.d0


! --- rotation of the reference frame (omega in rad/s)

!       If the rotation is not nul, then
!          icorio = 0: rotation is taken into account by rotating the mesh
!                      (simulation in the absolute frame)
!                 = 1: rotation is taken into account by Coriolis source terms
!                      (simulation in the relative frame)


icorio = 0

call rotation_define(0.d0, 0.d0, 0.d0,  &    ! rotation vector
                     0.d0, 0.d0, 0.d0)       ! invariant point



! --- Reference fluid properties

!       ro0        : density in kg/m3
!       viscl0     : dynamic viscosity in kg/(m s)
!       cp0        : specific heat in J/(Kelvin kg)
!       t0         : reference temperature in Kelvin
!       p0         : total reference pressure in Pascal
!                    the calculation is based on a
!                    reduced pressure P*=Ptot-ro0*g.(x-xref)
!                    (except in compressible case)
!       xyzp0(3)   : coordinates of the reference point for
!                    the total pressure (where it is equal to p0)

!     In general, it is not necessary to furnish a reference point xyz0.
!       If there are outlets, the code will take the center of the
!       reference outlet face.
!       On the other hand, if we plan to explicitly fix Dirichlet conditions
!       for pressure, it is better to indicate to which reference the
!       values relate (for a better resolution of reduced pressure).


!     Other properties are given by default in all cases.

!     Nonetheless, we may note that:

!       In the standard case (no gas combustion, coal, electric arcs,
!                             compressibility):
!       ---------------------
!         ro0, viscl0 and cp0
!             are useful and represent either the fluid properties if they
!             are constant, either simple mean values for the initialization
!             if properties are variable and defined in usphyv.
!         t0  is not useful
!         p0  is useful but is not used in an equation of state. p0
!             is a reference value for the incompressible solver
!             which will serve to set the (possible) domain outlet pressure.
!             We may also take it as 0 or as a physical value in Pascals.

!       With the electric module:
!       ------------------------
!         ro0, viscl0 and cp0
!             are useful but simply represent mean initial values;
!             the density, molecular dynamic viscosity, and specific
!             heat are necessarily defined as fields (whether they are
!             physically variable or not): see cs_user_physical_properties
!             for the Joule effect
!             module and the electric arcs dp_ELE data file.
!         t0  is useful an must be in Kelvin (> 0) but represents a simple
!             initialization value.
!         p0  is useful bu is not used in the equation of state. p0
!             is a reference value for the incompressible solver which
!             will be used to calibrate the (possible) outlet pressure
!             of the domain. We may take it as zero or as a physical
!             value in Pascals.

!       With gas combustion:
!       --------------------
!         ro0 is not useful (it is automatically recalculated by the
!             law of ideal gases from t0 and p0).
!         viscl0 is indispensable: it is the molecular dynamic viscosity,
!             assumed constant for the fluid.
!         cp0 is indispensable: it is the heat capacity, assumed constant,
!             (modelization of source terms involving a local Nusselt in
!             the Lagrangian module, reference value allowing the
!             calculation of a radiative
!             (temperature, exchange coefficient) couple).
!         t0  is indispensible and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).

!       With pulverized coal:
!       ---------------------
!         ro0 is not useful (it is automatically recalculated by the
!             law of ideal gases from t0 and p0).
!         viscl0 is indispensable: it is the molecular dynamic viscosity,
!             assumed constant for the fluid (its effect is expected to
!             be small compared to turbulent effects).
!         cp0 is indispensable: it is the heat capacity, assumed constant,
!             (modelization of source terms involving a local Nusselt in
!             the coal or Lagrangian module, reference value allowing the
!             calculation of a radiative
!             (temperature, exchange coefficient) couple).
!         t0  is indispensable and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).

!       With compressibility:
!       ---------------------
!         ro0 is not useful, stricto sensu; nonetheless, as experience
!             shows that users often use this variable, it is required
!             to assign to it a strictly positive value (for example,
!             an initial value).
!         viscl0 is useful and represents the molecular dynamic viscosity,
!             when it is constant, or a value which will be used during
!             initializations (or in inlet turbulence conditions,
!             depending on the user choice.
!         cp0 is indispensable: it is the heat capacity, assumed constant
!             in the thermodynamics available by default
!         t0  is indispensable and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).
!             With the thermodynamic law available by default,
!             t0 and p0 are used for the initialization of the density.
!         xyzp0 is not useful because the pressure variable directly
!             represents the total pressure.

ro0    = 1.17862d0
viscl0 = 1.83337d-5
cp0    = 1017.24d0

t0 = 20.d0 + 273.15d0
p0 = 1.01325d5


! --- irovar, ivivar, icp: constant or variable density,
!                          viscosity/diffusivity, and specific heat

!     When a specific physics module is active
!       (coal, combustion, electric arcs, compressible: see usppmo)
!       we MUST NOT set variables 'irovar', 'ivivar', and 'icp' here, as
!       they are defined automatically.
!     Nonetheless, for the compressible case, ivivar may be modified
!       in the uscfx2 user subroutine.

!     When no specific physics module is active, we may specify if the
!       density, specific heat, and the molecular viscosity
!       are constant (irovar=0, ivivar=0, icp=-1), which is the default
!       or variable (irovar=1, ivivar=1, icp=0)

!     For those properties we choose as variable, the corresponding law
!       must be defined in usphyv
!       (incs_user_physical_properties.f90);
!       if they are constant, they take values ro0, viscl0, and cp0.

irovar = 1
ivivar = 1
icp = -1

! We only specify XYZ0 if we explicitely fix Dirichlet conditions
! for the pressure.

xyzp0(1) = 0.d0
xyzp0(2) = 0.d0
xyzp0(3) = 0.d0

! --- Variable diffusivity field id (ifcvsl>=0) or constant
!     diffusivity (ifcvsl=-1) for the thermal scalar and USER scalars.

!     With ifcvsl = 0, the field will be added automatically, and later calls to
!       field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
!       will return its id.
!     With ifcvsl > 0, the id of an existing, predifined field is given. This
!       may allow sharing a diffusivity between multiple scalars.

!     For user scalars iscal which represent the variance of another user
!       scalar, the diffusivity of the variance of a scalar is assumed to
!       have the same behavior as the diffusivity of this scalar,
!       so values set here will be ignored.

!     For non-user scalars relative to specific physics (coal, combustion,
!       electric arcs: see usppmo) implicitly defined in the model,
!       the diffusivity should not be modified here.

!     Caution:   complete usphyv with the law defining the diffusivity
!     ========   if and only if ifcvsl = 0 has been set here.


! For thermal scalar
if (ippmod(icompf).ge.0) then
  ifcvsl = -1
  call field_set_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
else if (iscalt.gt.0) then
  ifcvsl = -1
  call field_set_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)
endif

do iscal = 1, nscaus
  if (iscavr(iscal).le.0) then
    ifcvsl = -1
    call field_set_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
  endif
enddo


! --- Variable density field id (ifcvsl>=0) or bulk
!     density (ifcvsl=-1) for USER scalars.

!     With ifcvsl = 0, the field will be added automatically, and later calls to
!       field_get_key_int(ivarfl(isca(iscal)), kromsl, ifcvsl)
!       will return its id.
!     With ifcvsl > 0, the id of an existing, predifined field is given. This
!       may allow sharing a density between multiple scalars.

!     For user scalars iscal which represent the variance of another user
!       scalar, the density of the variance of a scalar is assumed to
!       have the same behavior as the density of this scalar,
!       so values set here will be ignored.

!     Caution:   complete usphyv with the law defining the density
!     ========   if and only if ifcvsl = 0 has been set here.

do iscal = 1, nscaus
  if (iscavr(iscal).le.0) then
    ifcvsl = -1
    call field_set_key_int(ivarfl(isca(iscal)), kromsl, ifcvsl)
  endif
enddo


! --- Turbulent flux model u'T' for the scalar T
!     Algebraic Model
!      0  SGDH
!      10 GGDH
!      11 EB-GGDH (Elliptic Blending)
!      20 AFM
!      21 EB-AFM (Elliptic Blending)
!     Model with transport equations
!      30 DFM
!      31 EB-DFM (Elliptic Blending)

! GGDH for thermal scalar:
if (iscalt.gt.0) iturt(iscalt) = 10

! GGDH for all the scalars:
do jj = 1, nscaus
  iturt(jj) = 10
enddo


! --- Reference velocity for turbulence initialization (m2/s)
!       (useful only with turbulence)

uref = 1.d0


! --- Reference length scale in meters for initialization
!       of epsilon (and specific clipping of turbulence, but
!       this is not the default option)
!       Assign a value of the order of the largest dimension of the
!       physical domain in which the flow may develop.
!       If a negative value is set here, or no value set and the GUI not
!       used, the cubic root of the domain will be used.
!       (useful only for turbulence).

almax = 0.5

! Error estimators for Navier-Stokes (non-frozen velocity field)

! We recommend running a calculation restart on a few time steps
! with the activation of the most interesting of those.
! (=2 to activate, =0 to deactivate).

iescal(iescor) = 2   ! div(rho u) -Gamma
iescal(iestot) = 2   ! resolution precision for the momentum

! ALE (Arbitrary Lagrangian Eulerian) related options
!====================================================

! Number of iterations for fluid initialization. Contrary to ntmabs,
! nalinf is not an absolute iteration number, meaning that in case of
! restart calculation nalinf corresponds to the number of iterations
! for fuid initialization beginning from the first current iteration of
! the calculation restart. In general nalinf = 0 in that case.

nalinf = 75


! Maximum number of iterations in case of implicit Fluid Structure Coupling
! with structural calculations (internal and/or external
! (i.e. using code_aster)).
! nalimx = 1, in case of explicit FSI algorithm.

nalimx = 15


! Relative precision of sub-cycling Fluid Structure Coupling algorithm.

epalim = 1.d-5


!< [usipsu]

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

subroutine usipes &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use field
use parall
use period
use ihmpre
use post
use ppppar
use ppthch
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, ipp, f_id

type(var_cal_opt) :: vcopt

!===============================================================================

!>     This subroutine allows setting parameters
!>     which do not already appear in the other subroutines of this file.
!>
!>     It is possible to add or remove parameters.
!>     The number of physical properties and variables is known here.

!===============================================================================

!===============================================================================
! 1. Logging
!===============================================================================

! Frequency of log output

ntlist = 1


! Log verbosity

do ii = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  vcopt%iwarni = 1
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
enddo

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
vcopt%iwarni = 2
call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
vcopt%iwarni = 2
call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)

!===============================================================================
! 2. Definition of deformable structure time plots
!===============================================================================

! structures output step

nthist = 1
frhist = -1.d0

tplfmt = 1 ! time plot format (1: .dat, 2: .csv, 3: both)

!===============================================================================
! 3. Fine control of variables output
!===============================================================================

! Per variable output control.
! More examples are provided in cs_user_parameters-output.f90

! User scalar variables.

if (isca(1).gt.0.and.nscaus.ge.1) then
  f_id = ivarfl(isca(1))
  call field_set_key_str(f_id, keylbl, 'Scalar 1')
  call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
  call field_set_key_int(f_id, keylog, 1)
endif

if (isca(2).gt.0.and.nscaus.ge.2) then
  f_id = ivarfl(isca(2))
  call field_set_key_str(f_id, keylbl, 'Scalar 2')
  call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
  call field_set_key_int(f_id, keylog, 1)
endif


return
end subroutine usipes


!===============================================================================


!> \brief Initialize non-standard calculation options for the atmospheric version.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine usati1

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use atincl
use atsoil
use atchem
use atimbr
use siream

!===============================================================================

implicit none

!===============================================================================

!< [usati1]

!===============================================================================
! 1. Example of calculation options to modify
!===============================================================================

!!! Reading the meteo file

imeteo = 1

!!! For radiative model or chemistry

! Time of the simulation
!  syear  --> starting year
!  squant --> starting quantile
!  shour  --> starting hour (UTC)
!  smin   --> starting minute
!  ssec   --> starting second

syear = 1994
squant = 1
shour = 1
smin = 0
ssec = 0.d0

! Geographic position
! xlon --> longitude of the domain origin
! xlat --> latitude of the domain origin

xlon = 0.d0
xlat = 45.d0

!  -----------------------------------------------------------------------------
!  Atmospheric imbrication on large scale meteo (atimbr module)
!  -----------------------------------------------------------------------------
!
! --------------------------------------------------------------
! activation flag
! --------------------------------------------------------------
imbrication_flag    = .false.
imbrication_verbose = .false.

! ------------------------------------------------------------------------------
! flags for activating the cressman interpolation for the boundary conditions
! ------------------------------------------------------------------------------
cressman_u     = .true.
cressman_v     = .true.
cressman_tke   = .true.
cressman_eps   = .true.
cressman_theta = .true.
cressman_qw    = .true.
cressman_nc    = .true.

! --------------------------------------------------------------
! numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
horizontal_influence_radius = 8500.d0
vertical_influence_radius = 100.d0

! --------------------------------------------------------------

!!! Gaseous chemistry

! ichemistry: choice of chemistry resolution scheme
!0 --> no atmospheric chemistry
!1 --> quasi steady equilibrium NOx scheme with 4 species and 5 reactions
!2 --> scheme with 20 species and 34 reactions
!3 --> scheme CB05 with 52 species and 155 reactions
!4 --> user defined schema
ichemistry = 0

! ificchemistry: choice to read (=1,2,3,4, according to the scheme)
! or not (0) a concentration profile file
! if ichemistry>0 ifilechemistry is automaticaly set to ichemistry
ifilechemistry = 0

! isepchemistry: splitted (=1) or semi-coupled (=2, pu-sun)
! resolution of chemistry
isepchemistry = 1

! iphotolysis: inclusion (=1) or not (=2) of photolysis reactions
iphotolysis = 1

! dtchemmax: maximal time step (s) for chemistry resolution
dtchemmax = 10.0d0

!!! Aerosol chemistry

! iaerosol: flag to activate aerosol chemistry
! if iaerosol = 1, ichemistry is automatically set to 3 (scheme 3)
iaerosol = 1

! inogaseouschemistry: flag to prevent automatic resolution (=1)
! of gaseous chemistry (scheme 3)
inogaseouschemistry = 0

! ncycle_aer: number of iterations for time splitting
ncycle_aer = 1

! icoag_siream: flag to activate (=1) or not (=0) coagulation
icoag_siream = 1

! icond_siream: flag to activate (=1) or not (=0) condensation/evaporation
icond_siream = 1

! inucl_siream: flag to activate (=1) or not (=0) nucleation
inucl_siream = 1

! icut_siream: cutting bin between equilibrium (1 to icut_siream)
! and dynamic bins (icut_siream to nbin_aer)
icut_siream = nbin_aer

!< [usati1]

!----
! End
!----

return
end subroutine usati1


!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called:
!>
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!>  - Lagrangian module coupled with pulverized coal:
!>    Eulerian combustion of pulverized coal and
!>    Lagrangian transport of coal particles
!>  - Pulverised coal combustion
!>  - Fuel (oil) combustion
!
!===============================================================================

subroutine cs_user_combustion

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ihmpre
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use radiat

!===============================================================================

implicit none

!< [cs_user_combustion]

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! --- Kinetic model for CO <=> CO2

!     if = 0  unused (maximal conversion in turbulent model)
!     if = 1  transport of CO2 mass fraction
!     if = 2  transport of CO mass fraction (coal and fuel only)

ieqco2 = 0

! --- Density Relaxation
!     RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0

!===============================================================================
! 2. Physical Constants
!===============================================================================

! diftl0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5


! -----------------------------------------------------------------------------
! 2.1 For 3 points combusution model ONLY
! -----------------------------------------------------------------------------

! Reference temperature for fuel and oxydant (K)
tinfue = 436.d0
tinoxy = 353.d0


! -----------------------------------------------------------------------------
! 2.2 For EBU-model ONLY
! -----------------------------------------------------------------------------

! cebu: EBU-model constant
cebu   = 2.5d0


! -----------------------------------------------------------------------------
! 2.3 For Libby-Williams model ONLY
! -----------------------------------------------------------------------------

! Reference velocity
vref = 60.d0
! Reference length scale
lref = 0.1d0
! Activation Temperature
ta   = 0.2d5
! Cross-over Temperature (combustion of propane)
tstar= 0.12d4

!< [cs_user_combustion]

!----
! End
!----

return
end subroutine cs_user_combustion


!===============================================================================

!> \brief User subroutine.

!> Initialize non standard options for the compressible flow scheme such
!> as the variability of the thermal conductivity and the volume viscosity.
!> Their values can be given in the subroutine \ref uscfx2 .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx1

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Arguments


! Local variables

integer :: ifcvsl, nphases

double precision :: cv(2), gamma(2), pinf(2), qprim(2)

!===============================================================================

!===============================================================================

!< [uscfx1]

!===============================================================================
! 1. Properties options
!===============================================================================

if (iihmpr.eq.0) then   !  Remove test to set values here when also using GUI.

  ! --> Molecular thermal conductivity
  !       constant  : ifcvsl = -1
  !       variable  : ifcvsl = 0

  ifcvsl = -1
  call field_set_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)

  ! --> Volumetric molecular viscosity
  !       iviscv = -1 : uniform  in space and constant in time
  !              =  0 : variable in space and time

  iviscv = -1

endif

!< [uscfx1]

!----
! End
!----

return
end subroutine uscfx1


!===============================================================================


!> \brief User subroutine.
!>
!> Set values for the reference volumic viscosity, the reference
!> conductivity and the molar mass for compressible flow.
!>
!> Initialize non standard options for the compressible flow scheme such
!> as the hydrostatic equilibrium.
!>
!> In addition to options set in the user subroutine \ref uscfx1 (or in
!> the GUI): this subroutine allows to set a switch to indicate if the
!> molecular viscosity is constant, its values being given in the user
!> subroutine \ref usipsu .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx2

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

! Local variables

!< [uscfx2]

!===============================================================================
! 1. Physical properties
!===============================================================================

if (iihmpr.eq.0) then   !  Remove test to set values here when also using GUI.

! --> Molecular viscosity
!       constant  : ivivar = 0
!       variable  : ivivar = 1

  ivivar = 0

! --> Reference molecular thermal conductivity
!       visls0 = lambda0  (molecular thermal conductivity, W/(m K))

!       WARNING: visls0 must be strictly positive
!         (set a realistic value here even if conductivity is variable)

  visls0(itempk) = 3.d-2

!       If the molecular thermal conductivity is variable, its values
!         must be provided in the user subroutine 'usphyv'

! --> Volumetric molecular viscosity

!       Reference volumetric molecular viscosity

!       viscv0 = kappa0  (volumetric molecular viscosity, kg/(m s))

  viscv0 = 0.d0

!       If the volumetric molecular viscosity is variable, its values
!         must be provided in the user subroutine 'usphyv'

! --> Molar mass of the gas (kg/mol)

!       For example with dry air, xmasml is around 28.8d-3 kg/mol

  xmasmr = 0.028966

! --> Hydrostatic equilibrium

!       Specify if the hydrostatic equilibrium must be accounted for
!         (yes = 1 , no = 0)

  icfgrp = 1

endif

!< [uscfx2]

!----
! End
!----

return
end subroutine uscfx2


!===============================================================================

!> \brief Definition of cooling tower model and exchange zones

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine cs_user_cooling_towers

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! 1. Parameters for prescibed temperature difference
!===============================================================================

! Air and liquid properties
!< [cs_user_cooling_towers]
cp_a    = 1006.0d0
cp_v    = 1831.0d0
cp_l    = 4179.0d0
hv0    = 2501600.0d0
rho_l   = 997.85615d0
viscl0 = 1.765d-5
lambda_l = 0.02493d0
humidity0 = 0.d0
droplet_diam = 0.005d0
!< [cs_user_cooling_towers]

!----
! End
!----

return
end subroutine cs_user_cooling_towers


!===============================================================================

!> \brief User routine for definition of computation parameters dealing
!>        with Darcy module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine user_darcy_ini1

!===============================================================================
! Module files
!===============================================================================

use ihmpre, only: iihmpr
use entsor
use darcy_module

!===============================================================================

implicit none

!===============================================================================

!< [user_darcy_ini1]

darcy_anisotropic_permeability = 0 ! permeability : 0 isotrop, 1 anisotrop

darcy_anisotropic_dispersion = 0 ! dispersion : 0 isotrop, 1 anisotrop

darcy_unsteady = 0 ! 0 steady flow, 1 unsteady flow

darcy_convergence_criterion = 0 ! convergence criterion of Newton scheme : 0, over pressure, 1, over velocity

darcy_gravity = 0 ! gravity is taken into account : 0 no, 1 yes

!< [user_darcy_ini1]

!----
! End
!----

return

end subroutine user_darcy_ini1


!===============================================================================
