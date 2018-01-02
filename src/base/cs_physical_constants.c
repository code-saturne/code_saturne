/*============================================================================
 * Base physical constants data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_physical_constants.c
        Base physical constants and fluid properties data.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_physical_constants_t

  \brief Physical constants descriptor.

  Members of this structure are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_physical_constants_t::gx
        x component of the gravity vector
  \var  cs_physical_constants_t::gy
        y component of the gravity vector
  \var  cs_physical_constants_t::gz
        z component of the gravity vector
  \var  cs_physical_constants_t::icorio
        Coriolis source terms
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_fluid_properties_t

  \brief Fluid properties descriptor.

  Members of these fluid properties are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_fluid_properties_t::ixyzp0
        filling \ref xyzp0 indicator
  \var  cs_fluid_properties_t::ieos
        indicator of equation of state

        useful only for the compressible module. Ideal gas with a constant
        adiabatic coefficient (\ref ieos=1), stiffened gas (\ref ieos=2) and
        mix of ideal gas (\ref ieos=3) are available, but the user can complete
        the file \ref cs_cf_thermo.h, which is not a user source, to add new
        equations of state.
  \var  cs_fluid_properties_t::icp
        indicates if the isobaric specific heat \f$C_p\f$ is variable
        - -1: uniform, no property field is declared
        -  0: variable, \f$C_p\f$ is declared as a property field\n
        When gas or coal combustion is activated, \ref icp is automatically
        set to 0 (constant \f$C_p\f$). With the electric module, it is
        automatically set to 1. The user is not allowed to modify these
        default choices.\n When \ref icp = 1 is specified, the code
        automatically modifies this value to make \ref icp designate the
        effective index-number of the property "specific heat". For each cell
        iel, the value of \f$C_p\f$ is then specified by the user in the
        appropriate subroutine (\ref cs_user_physical_properties for the
        standard physics).\n Useful if there is
        1\f$\leqslant\f$N\f$\leqslant\f$\ref dimens::nscal "nscal" so that
        iscsth(n)=1 (there is a scalar temperature) or with the compressible
        module for non perfect gases.
  \var  cs_fluid_properties_t::icv
        property index of the isochoric specific heat
        - -1: uniform isochoric specific heat (no property field defined)
  \var  cs_fluid_properties_t::irovar
        variable density field \f$ \rho \f$:
        - 1: true, its variation law be given either
        in the GUI, or in the user subroutine
        \ref cs_user_physical_properties .\n
        See \subpage physical_properties for more informations.
        - 0: false, its value is the reference density
        \ref ro0.
  \var  cs_fluid_properties_t::ivivar
        variable viscosity field \f$ \mu \f$:
           - 1: true, its variation law be given either
        in the GUI, or in the user subroutine
        \ref cs_user_physical_properties .\n
        See \subpage physical_properties for more informations.
           - 0: false, its value is the reference molecular
        dynamic viscosity \ref viscl0
  \var  cs_fluid_properties_t::ivsuth
        Sutherland law for laminar viscosity and thermal conductivity
        Only useful in gas mix (igmix) specific physics
        - 1: Sutherland law
        - 0: low temperature law (linear except for helium)
  \var  cs_fluid_properties_t::ro0
        reference density

        Negative value: not initialized.
        Its value is not used in gas or coal combustion modelling (it will be
        calculated following the perfect gas law, with \f$P_0\f$ and \f$T_0\f$).
        With the compressible module, it is also not used by the code, but it
        may be (and often is) referenced by the user in user subroutines; it is
        therefore better to specify its value.

        Always useful otherwise, even if a law defining the density is given by
        the user subroutines \ref usphyv or \ref cs_user_physical_properties.
        Indeed, except with the compressible module, CS  does not use the total
        pressure \f$P\f$ when solving the Navier-Stokes equation, but a reduced
        pressure \f$ P^*=P-\rho_0\vect{g}.(\vect{x}-\vect{x}_0)+ P^*_0-P_0 \f$,
        where \f$\vect{x_0}\f$ is a reference point (see \ref xyzp0) and \f$
        P^*_0 \f$ and \f$ P_0 \f$ are reference values (see \ref pred0 and
        \ref p0). Hence, the term \f$-\grad{P}+\rho\vect{g}\f$ in the equation
        is treated as \f$-\grad{P^*}+(\rho-\rho_0)\vect{g}\f$. The closer
        \ref ro0 is to the value of \f$
        \rho\f$, the more \f$P^*\f$ will tend to
        represent only the dynamic part of the pressure and the faster and more
        precise its solution will be. Whatever the value of \ref ro0, both \f$
        P\f$ and \f$P^*\f$ appear in the listing and the post-processing
        outputs with the compressible module, the calculation is made directly
        on the total pressure.
  \var  cs_fluid_properties_t::viscl0
        reference molecular dynamic viscosity

        Negative value: not initialized.

        Always useful, it is the used value unless the user specifies the
        viscosity in the subroutine \ref usphyv.
  \var  cs_fluid_properties_t::p0
        reference pressure for the total pressure

        Except with the compressible module, the total pressure \f$P\f$ is
        evaluated from the reduced pressure \f$P^*\f$ so that \f$P\f$ is equal
        to \ref p0 at the reference position \f$\vect{x}_0\f$ (given by
        \ref xyzp0).
        With the compressible module, the total pressure is solved directly.
        Always useful.
  \var  cs_fluid_properties_t::pred0
        reference value for the reduced pressure \f$P^*\f$ (see \ref ro0)

        It is especially used to initialise the reduced pressure and as a
        reference value for the outlet boundary conditions.
        For an optimised precision in the resolution of \f$P^*\f$, it is wiser
        to keep \ref pred0 to 0.
        With the compressible module, the "pressure" variable appearing in the
        equations directly represents the total pressure.
        It is therefore initialized to \ref p0 and not \ref pred0 (see
        \ref ro0).
        Always useful, except with the compressible module.
  \var  cs_fluid_properties_t::xyzp0[3]
        coordinates of the reference point \f$\vect{x}_0\f$ for the total
        pressure

        - When there are no Dirichlet conditions for the pressure (closed
        domain), \ref xyzp0 does not need to be specified (unless the total
        pressure has a clear physical meaning in the configuration treated).
        - When Dirichlet conditions on the pressure are specified but only
        through stantard outlet conditions (as it is in most configurations),
        \ref xyzp0 does not need to be specified by the user, since it will be
        set to the coordinates of the reference outlet face (\em i.e. the code
        will automatically select a reference outlet boundary face and set
        \ref xyzp0 so that \f$P\f$ equals \ref p0 at this face). Nonetheless, if
        \ref xyzp0 is specified by the user, the calculation will remain
        correct.
        - When direct Dirichlet conditions are specified by the user (specific
        value set on specific boundary faces), it is better to specify the
        corresponding reference point (\em i.e. specify where the total pressure
        is \ref p0). This way, the boundary conditions for the reduced pressure
        will be close to \ref pred0, ensuring an optimal precision in the
        resolution. If \ref xyzp0 is not specified, the reduced pressure will be
        shifted, but the calculations will remain correct.
        - With the compressible module, the "pressure" variable appearing in the
        equations directly represents the total pressure. \ref xyzp0 is
        therefore not used.

        Always useful, except with the compressible module.
  \var  cs_fluid_properties_t::t0
        reference temperature

        Useful for the specific physics gas or coal combustion (initialization
        of the density), for the electricity modules to initialize the domain
        temperature and for the compressible module (initializations).
        It must be given in Kelvin.
  \var  cs_fluid_properties_t::cp0
        reference specific heat

        Useful if there is 1 <= n <= nscaus
        so that \ref cs_thermal_model_t::iscalt "cs_glob_thermal_model->iscalt" = n
        and \ref cs_thermal_model_t::itherm "cs_glob_thermal_model->itherm" =  1
       (there is a scalar "temperature"), unless the
        user specifies the specific heat in the user subroutine \ref usphyv
        (\ref cstphy::icp "icp" > 0) with the compressible module or coal combustion,
        \ref cp0 is also needed even when there is no user scalar. \note
        None of the scalars from the specific physics is a temperature. \note
        When using the Graphical Interface, \ref cp0 is also used to
        calculate the diffusivity of the thermal scalars, based on their
        conductivity; it is therefore needed, unless the diffusivity is also
        specified in \ref usphyv.
  \var  cs_fluid_properties_t::cv0
        reference isochoric specific heat (J/kg/K)

        Useful for the compressible module
  \var  cs_fluid_properties_t::xmasmr
        molar mass of the perfect gas in \f$ kg/mol \f$
        (if \ref cstphy::ieos "ieos"=1)

        Always useful.
  \var  cs_fluid_properties_t::psginf
        stiffened gas (\ref cs_fluid_properties_t::ieos "ieos"=2) limit pressure
        (zero in perfect gas)
  \var  cs_fluid_properties_t::gammasg
        stiffened gas (\ref cs_fluid_properties_t::ieos "ieos"=2) polytropic coefficient
        (dimensionless)
  \var  cs_fluid_properties_t::ipthrm
        uniform variable thermodynamic pressure:
        - 0: false
        - 1: true
  \var  cs_fluid_properties_t::pther
        Thermodynamic pressure for the current time step.
  \var  cs_fluid_properties_t::pthera
        thermodynamic pressure for the previous time step
  \var  cs_fluid_properties_t::pthermax
        thermodynamic maximum pressure for user clipping, used to model a
        venting effect
  \var  cs_fluid_properties_t::sleak
        Leak surface
  \var  cs_fluid_properties_t::kleak
        Leak head loss (2.9 by default, from Idelcick)
  \var  cs_fluid_properties_t::roref
        Initial reference density
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main physical constants structure and associated pointer */

static cs_physical_constants_t _physical_constants = {
  .gravity[0] = 0.,
  .gravity[1] = 0.,
  .gravity[2] = 0.,
  .icorio = 0};

/* main fluid properties structure and associated pointer */

static cs_fluid_properties_t  _fluid_properties = {
  .ixyzp0   = -1,
  .ieos     = -1,
  .icp      = -1,
  .icv      = -1,
  .irovar   = 0,
  .ivivar   = 0,
  .ivsuth   = 0,
  .ro0      = 1.17862,
  .viscl0   = 1.83337e-5,
  .p0       = 1.01325e5,
  .pred0    = 0.,
  .xyzp0    = {-1.e30, -1.e30, -1.e30},
  .t0       = 293.15,
  .cp0      = 1017.24,
  .cv0      = 0.,
  .xmasmr   = 0.028966, /* air molar mass */
  .psginf   = 0.,
  .gammasg  = 1.4,
  .ipthrm   = 0,
  .pther    = 1.013e5,
  .pthera   = 0.,
  .pthermax = -1.,
  .sleak    = 0.,
  .kleak    = 2.9,
  .roref    = 1.17862};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Ideal gas constant (\f$J.mol^{-1}.K^{-1}\f$) */

const double cs_physical_constants_r = 8.31446;

/*! Boltzmann constant (\f$J.K^{-1}\f$) */

const double cs_physical_constants_kb = 1.38e-23;

/*! Conversion from Celsius to Kelvin: 275.15   */

const double cs_physical_constants_celsius_to_kelvin = 273.15;

/* !Stephan constant for the radiative module \f$\sigma\f$
   in \f$W.m^{-2}.K^{-4}\f$ */

const double cs_physical_constants_stephn = 5.6703e-8;

const cs_physical_constants_t  *cs_glob_physical_constants =
                                  &_physical_constants;

const cs_fluid_properties_t  *cs_glob_fluid_properties = &_fluid_properties;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_physical_constants_get_pointers(double  **gx,
                                     double  **gy,
                                     double  **gz,
                                     int     **icorio);

void
cs_f_fluid_properties_get_pointers(int     **ixyzp0,
                                   int     **ieos,
                                   int     **icp,
                                   int     **icv,
                                   int     **irovar,
                                   int     **ivivar,
                                   int     **ivsuth,
                                   double  **ro0,
                                   double  **viscl0,
                                   double  **p0,
                                   double  **pred0,
                                   double  **xyzp0,
                                   double  **t0,
                                   double  **cp0,
                                   double  **cv0,
                                   double  **xmasmr,
                                   double  **psginf,
                                   double  **gammasg,
                                   int     **ipthrm,
                                   double  **pther,
                                   double  **pthera,
                                   double  **pthermax,
                                   double  **sleak,
                                   double  **kleak,
                                   double  **roref);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global physical constants structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   gx     --> pointer to cs_glob_physical_constants->gravity[0]
 *   gy     --> pointer to cs_glob_physical_constants->gravity[1]
 *   gz     --> pointer to cs_glob_physical_constants->gravity[2]
 *   icorio --> pointer to cs_glob_physical_constants->icorio
 *----------------------------------------------------------------------------*/

void
cs_f_physical_constants_get_pointers(double  **gx,
                                     double  **gy,
                                     double  **gz,
                                     int     **icorio)
{
  *gx     = &(_physical_constants.gravity[0]);
  *gy     = &(_physical_constants.gravity[1]);
  *gz     = &(_physical_constants.gravity[2]);
  *icorio = &(_physical_constants.icorio);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global fluid properties structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ixyzp0   --> pointer to cs_glob_fluid_properties->ixyzp0
 *   ieos     --> pointer to cs_glob_fluid_properties->ieos
 *   icp      --> pointer to cs_glob_fluid_properties->icp
 *   icv      --> pointer to cs_glob_fluid_properties->icv
 *   irovar   --> pointer to cs_glob_fluid_properties->irovar
 *   ivivar   --> pointer to cs_glob_fluid_properties->ivivar
 *   ivsuth   --> pointer to cs_glob_fluid_properties->ivsuth
 *   ro0      --> pointer to cs_glob_fluid_properties->ro0
 *   viscl0   --> pointer to cs_glob_fluid_properties->viscl0
 *   p0       --> pointer to cs_glob_fluid_properties->p0
 *   pred0    --> pointer to cs_glob_fluid_properties->pred0
 *   xyzp0    --> pointer to cs_glob_fluid_properties->xyzp0
 *   t0       --> pointer to cs_glob_fluid_properties->t0
 *   cp0      --> pointer to cs_glob_fluid_properties->cp0
 *   cv0      --> pointer to cs_glob_fluid_properties->cv0
 *   xmasmr   --> pointer to cs_glob_fluid_properties->xmasmr
 *   psginf   --> pointer to cs_glob_fluid_properties->psginf
 *   gammasg  --> pointer to cs_glob_fluid_properties->gammasg
 *   ipthrm   --> pointer to cs_glob_fluid_properties->ipthrm
 *   pther    --> pointer to cs_glob_fluid_properties->pther
 *   pthera   --> pointer to cs_glob_fluid_properties->pthera
 *   pthermax --> pointer to cs_glob_fluid_properties->pthermax
 *   sleak    --> pointer to cs_glob_fluid_properties->sleak
 *   kleak    --> pointer to cs_glob_fluid_properties->kleak
 *   roref    --> pointer to cs_glob_fluid_properties->roref
 *----------------------------------------------------------------------------*/

void
cs_f_fluid_properties_get_pointers(int     **ixyzp0,
                                   int     **ieos,
                                   int     **icp,
                                   int     **icv,
                                   int     **irovar,
                                   int     **ivivar,
                                   int     **ivsuth,
                                   double  **ro0,
                                   double  **viscl0,
                                   double  **p0,
                                   double  **pred0,
                                   double  **xyzp0,
                                   double  **t0,
                                   double  **cp0,
                                   double  **cv0,
                                   double  **xmasmr,
                                   double  **psginf,
                                   double  **gammasg,
                                   int     **ipthrm,
                                   double  **pther,
                                   double  **pthera,
                                   double  **pthermax,
                                   double  **sleak,
                                   double  **kleak,
                                   double  **roref)
{
  *ixyzp0   = &(_fluid_properties.ixyzp0);
  *ieos     = &(_fluid_properties.ieos);
  *icp      = &(_fluid_properties.icp);
  *icv      = &(_fluid_properties.icv);
  *irovar   = &(_fluid_properties.irovar);
  *ivivar   = &(_fluid_properties.ivivar);
  *ivsuth   = &(_fluid_properties.ivsuth);
  *ro0      = &(_fluid_properties.ro0);
  *viscl0   = &(_fluid_properties.viscl0);
  *p0       = &(_fluid_properties.p0);
  *pred0    = &(_fluid_properties.pred0);
  *xyzp0    =  (_fluid_properties.xyzp0);
  *t0       = &(_fluid_properties.t0);
  *cp0      = &(_fluid_properties.cp0);
  *cv0      = &(_fluid_properties.cv0);
  *xmasmr   = &(_fluid_properties.xmasmr);
  *psginf   = &(_fluid_properties.psginf);
  *gammasg  = &(_fluid_properties.gammasg);
  *ipthrm   = &(_fluid_properties.ipthrm);
  *pther    = &(_fluid_properties.pther);
  *pthera   = &(_fluid_properties.pthera);
  *pthermax = &(_fluid_properties.pthermax);
  *sleak    = &(_fluid_properties.sleak);
  *kleak    = &(_fluid_properties.kleak);
  *roref    = &(_fluid_properties.roref);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_physical_constants
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_physical_constants_t *
cs_get_glob_physical_constants(void)
{
  return &_physical_constants;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_fluid_properties
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_fluid_properties_t *
cs_get_glob_fluid_properties(void)
{
  return &_fluid_properties;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Print the physical constants structure to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_physical_constants_log_setup(void)
{

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Physical properties\n"
       "-------------------\n\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("    gravity:    [%14.5e, (Gravity x component)\n"
       "                 %14.5e, (Gravity y component)\n"
       "                 %14.5e] (Gravity z component)\n"
       "    icorio:      %14d (Coriolis source terms)\n\n"),
       cs_glob_physical_constants->gravity[0],
       cs_glob_physical_constants->gravity[1],
       cs_glob_physical_constants->gravity[2],
       cs_glob_physical_constants->icorio);
}

/*----------------------------------------------------------------------------
 *!
 * \brief Print the fluid properties structure to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_fluid_properties_log_setup(void)
{
  cs_log_printf
    (CS_LOG_SETUP,
     _("   Continuous phase:\n"
       "    ro0:         %14.5e (Reference density)\n"
       "    viscl0:      %14.5e (Ref. molecular dyn. visc.)\n"
       "    cp0:         %14.5e (Ref. specific heat)\n"
       "    icp:         %14d (> 0: Variable cp (cs_user_physical_properties))\n"
       "    p0:          %14.5e (Ref. total pressure)\n"
       "    pred0:       %14.5e (Ref. reduced pressure)\n"
       "    t0:          %14.5e (Ref. temperature)\n\n"
       "    irovar:      %14d (Density constant(0) or not(1))\n"
       "    ivivar:      %14d (Molec. visc cst.(0) or not(1))\n\n"
       "    Initial reference point for pressure\n"
       "    xyzp0:       %14.5e %14.5e %14.5e\n"),
       cs_glob_fluid_properties->ro0,
       cs_glob_fluid_properties->viscl0,
       cs_glob_fluid_properties->cp0,
       cs_glob_fluid_properties->icp,
       cs_glob_fluid_properties->p0,
       cs_glob_fluid_properties->pred0,
       cs_glob_fluid_properties->t0,
       cs_glob_fluid_properties->irovar,
       cs_glob_fluid_properties->ivivar,
       cs_glob_fluid_properties->xyzp0[0],
       cs_glob_fluid_properties->xyzp0[1],
       cs_glob_fluid_properties->xyzp0[2]);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
