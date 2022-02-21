/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*
===============================================================================
 Function:
 ---------

> \file cs_user_wall_condensation.c
>
> \brief Source terms associated at the boundary faces and the neighboring
> cells with surface condensation.
>
> This subroutine fills the condensation source terms for each variable at
> the cell center associated to the boundary faces identifed in the mesh.
> The fluid exchange coefficient is computed with a empiric law to be
> imposed at the boundary face where the condensation phenomenon occurs.
>
> This user subroutine is called which allows the setting of
> \f$ \gamma_{\mbox{cond}} \f$ the condensation source term.
>
> This function fills the condensation source term array gamma_cond adding
> to the following equations:
>
> - The equation for mass conservation:
> \f[ D\frac{rho}{dt} + divs \left( \rho \vect{u}^n\right) = \Gamma _{cond}
> \f]
>
> - The equation for a variable \f$\Phi \f$:
> \f[ D\frac{\phi}{dt} = ... + \Gamma _{cond}*(\Phi _i - \Phi)
> \f]
>
> discretized as below:
>
> \f[ \rho*\dfrac{\Phi^{n+1}-\Phi^{n}}/dt = ...
>                            + \Gamma _{cond}*(\Phi _i - \Phi^{n+1})
> \f]
>
> \remarks
>  - \f$ \Phi _i \f$ is the value of \f$ \Phi \f$ associated to the
>    injected condensation rate.
>
>    With 2 options are available:
>       - the condensation rate is injected with the local value
>         of variable \f$ \Phi = \Phi ^{n+1}\f$
>         in this case the \f$ \Phi \f$ variable is not modified.
>
>       - the condensation rate is injected with a specific value
>         for \f$ \Phi = \Phi _i \f$ the specified value given by the
>         user.
>
> \section use Usage
>
> The three stages in the code where this User subroutine
> is called (with \code iappel = 1, 2 and 3\endcode)
>
> \code iappel = 1 \endcode
>  - Calculation of the number of cells where a mass source term is
>    imposed: ncesmp
>    Called once at the beginning of the calculation
>
> \code iappel = 2 \endcode
>   - Identification of the cells where a mass source term is imposed:
>     array icesmp(ncesmp)
>     Called once at the beginning of the calculation
>
> \code iappel = 3 \endcode
>   - Calculation of the values of the mass source term
>     Called at each time step
>
> \section the specific variables to define with is user subroutine
>
>  - ifbpcd(ieltcd): identification of the faces where a condensation
>                    source term is imposed.
>
>  - itypcd(ieltcd,ivar): type of treatment for variable ivar in the
>                       ieltcd cell with condensation source term.
>                     - itypcd = 0 --> injection of ivar at local value
>                     - itypcd = 1 --> injection of ivar at user
>                                      specified value.
>
>  - spcond(ielscd,ipr): value of the injection condensation rate
>                       gamma_cond (kg/m3/s) in the ieltcd cell
>                       with condensation source term.
>
>  - spcond(ieltcd,ivar): specified value for variable ivar associated
>                        to the injected condensation in the ieltcd
>                        cell with a condensation source term except
>                        for ivar=ipr.
>
> \remarks
>  - For each face where a condensation source terms is imposed ielscd
>    in [1;nfbpcd]), ifbpcd(ielscd) is the global index number of the
>    corresponding face (ifbpcd(ieltcd) in [1;ncel]).
>  - if itypcd(ieltcd,ivar)=0, spcond(ielpcd,ivar) is not used.
>  - if spcond(ieltcd,ipr)<0, mass is removed from the system,
>     therefore Code_Saturna automatically considers f_i=f^(n+1),
>     whatever the values of itypcd or smacel specified by the user
>
>   \par Examples of settings for boundary condensation mass source terms
>        Examples are available
>        \ref condens_h_boundary "here".
>
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
 Arguments
______________________________________________________________________________.
  mode           name          role                                           !
______________________________________________________________________________!
> \param[in]     nvar          total number of variables
> \param[in]     nscal         total number of scalars
> \param[in]     iappel        indicates which at which stage the routine is
_______________________________________________________________________________
*/

void
cs_user_wall_condensation(int nvar, int nscal, int iappel)
{
  cs_lnum_t  nlelt;
  cs_lnum_t *lstelt;
  cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;
  cs_lnum_t  nfabor = cs_glob_mesh->n_b_faces;

  cs_wall_cond_t *           wall_cond    = cs_get_glob_wall_cond();
  cs_wall_cond_1d_thermal_t *wall_thermal = cs_get_glob_wall_cond_1d_thermal();

  BFT_MALLOC(lstelt, nfabor, cs_lnum_t);

  int ieltcd = 0;
  int izone  = 0;

  /*==========================================================================
   1. One or two calls
   -------------------
    - iappel = 1: nfbpcd: calculation of the number of faces with
                               condensation source term
    - iappel = 2: ifbpcd: index number of faces with condensation source terms

   Remarks
   =======
    - Do not use spcond in this section (it is set on the third call, iappel=3)
    - Do not use ifbpcd in this section on the first call (iappel=1)
    - This section (iappel=1 or 2) is only accessed at the beginning of a
       calculation. Should the localization of the condensation source terms
       evolve in time, the user must identify at the beginning all cells that
       can potentially become a condensation  source term.
  ===========================================================================*/

  /*! < [zones_definition] */
  if (iappel == 1 || iappel == 2) {

    cs_zone_t *zone = cs_boundary_zone_by_name_try("cold_wall");

    for (ieltcd = 0; ieltcd < zone->n_elts; ieltcd++) {
      if (iappel == 2) {
        cs_lnum_t ifac             = zone->elt_ids[ieltcd];
        wall_cond->ifbpcd[ieltcd]  = ifac;
        wall_cond->izzftcd[ieltcd] = izone;
      }
    }
  }

  if (iappel == 1) {
    wall_cond->nfbpcd = ieltcd + 1;
    wall_cond->nzones = izone + 1;
  }

  /*! < [zones_definition] */

  int iz = 0; // Monozone

  /*! [model_settings] */

  /*
  ===============================================================================
   Parameters of the 1-D thermal model and condensation model
   ------------------------------------------------------------------
   Both models can be activated and coupled together or the condensation model
   can be used with a constant wall temperature specified by the user
   (at iappel=3 tpar=ztpar0(iz) in this case).
  ===============================================================================
  */

  if (iappel == 2) {
    if (wall_cond->icondb == 0) {

      /*
       *    izcophc = model for the mass transfer (condensation) coefficient
       *    ----------------------------------------------------------------
       *    Integer.
       *    1 : Turbulent wall law
       *    2 : Natural convection correlation
       *    3 : Maximum of the two previous (for mixed regime)
       *    */
      wall_cond->izcophc[iz] = 3;

      /*
       *    izcophg = model for the thermal exchange coefficient
       *    ----------------------------------------------------------------
       *    Integer.
       *    1 : Turbulent wall law
       *    2 : Natural convection correlation
       *    3 : Maximum of the two previous (for mixed regime)
       *    */
      wall_cond->izcophg[iz] = 3;

      /*
       *    iztag1d = on/off switch for 1D thermal module
       *    ----------------------------------------------------------------
       *    Integer.
       *    0 : Constant wall temperature (equal to ztpar0(iz))
       *    1 : Variable wall temperature computed with a 1D model
       *    */
      wall_cond->iztag1d[iz] = 1;

      if (wall_cond->iztag1d[iz] == 1) {

        /*
         *      ztheta = proportion of implicitation in the space discretization
         * scheme
         *      ----------------------------------------------------------------
         *      Float in the range [0, 1].
         *      Special values:
         *        0 : explicit scheme
         *        1 : implicit scheme
         *      */
        wall_thermal->ztheta[iz] = 1.0;

        /*
         *      zdxmin = Wall cell size parameters
         *      ----------------------------------------------------------------
         *      Float
         *      Special values:
         *        <=0 : Constant cell size
         *        > 0 : Variable cell size. In this case, the first cell size
         * (fluid side) is set to zdxmin [meters].
         *      */
        wall_thermal->zdxmin[iz] = 0.0;

        /*
         *      znmur = Number of cells in the wall mesh
         *      ----------------------------------------------------------------
         *      Positive integer
         *      */
        wall_thermal->znmur[iz] = 10;

        /*
         *      zepais = Total thickness of the solid wall [meters]
         *      ----------------------------------------------------------------
         *      Positive float
         *      */
        wall_thermal->zepais[iz] = 0.024;

        /*
         *      ztpar0 = Initial temperature in the solid wall [celsius]
         *      ----------------------------------------------------------------
         *      Float.
         *      */
        wall_thermal->ztpar0[iz] = 26.57;
      }
    }
  }
  /*! [model_settings] */

  /*
    ===============================================================================
     2. For nfbpcd > 0 , third call
        iappel = 3 : itypcd: type of condensation source term
                      spcond: condensation source term
     Remark
     ======
     If itypcd(ieltcd,ivar) is set to 1, spcond(ieltcd,ivar) must be set.
    ===============================================================================
  */

  else if (iappel == 3) {

    /*! [solid_wall] */

    // Fill in parameters for wall thermal mode
    // TODO : move to iappel == 2 ?
    if (wall_cond->icondb == 0) {
      if (wall_cond->iztag1d[iz] == 1) {

        /*
         *      zhext = External exchange coefficient
         * [watt.meter^(-2).kelvin^(-1)]
         *      ----------------------------------------------------------------
         *      Positive float.
         *      */
        wall_thermal->zhext[iz] = 1.e8;

        /*
         *      zhext = External temperature [celsius]
         *      ----------------------------------------------------------------
         *      Float.
         *      */
        wall_thermal->ztext[iz] = 26.57;

        /*
         *      zrob = Solid wall density [kilogram.meter^(-3)]
         *      ----------------------------------------------------------------
         *      Positive float.
         *      */
        wall_thermal->zrob[iz] = 8000.0;

        /*
         *      zcondb = Solid wall thermal conductivity
         * [watt.meter^(-1).celsius^(-1)]
         *      ----------------------------------------------------------------
         *      Positive float.
         *      */
        wall_thermal->zcondb[iz] = 12.8;

        /*
         *      zcpb = Solid wall specific heat
         * [joule.kilogram^(-1).celsius^(-1)]
         *      ----------------------------------------------------------------
         *      Positive float.
         *      */
        wall_thermal->zcpb[iz] = 500.0;
      }
      else {

        /*
         *      ztpar = Constant wall temperature [celsius]
         *      ----------------------------------------------------------------
         *      Float.
         *      */
        wall_cond->ztpar[iz] = 26.57;
      }
    }
    /*! [solid_wall] */

    /*
     * From here on to the end : fill in
     *
     * itypcd = for all variables except mass, type of source term
     * ---------------------------------------------------------------------
     *  Array of integers in [0, 1].
     *  0 : use ambient value
     *  1 : impose user value (in this case, corresponding value of spcond must
     *      be provided
     *
     * spcond = user values for condensation source term (if itypcd = 1)
     * ---------------------------------------------------------------------
     *  Array of floats.
     * */

    if (CS_F_(cp) == NULL)
      bft_error(__FILE__, __LINE__, 0, _("error lambda not variable\n"));

    if (CS_F_(h) == NULL)
      bft_error(__FILE__, __LINE__, 0, _("error lambda not variable\n"));

    cs_real_t *cpro_cp = CS_F_(cp)->val;
    cs_real_t *cvar_h  = CS_F_(h)->val;

    // Get specific heat of steam
    int                       k_id = cs_gas_mix_get_field_key();
    cs_field_t *              f    = cs_field_by_name("y_h2o_g");
    cs_gas_mix_species_prop_t gmp;
    cs_field_get_key_struct(f, k_id, &gmp);
    cs_real_t cp_vap = gmp.cp;

    // Get variable ids of quantities of interest
    const int var_id_key = cs_field_key_id("variable_id");

    f      = cs_field_by_name("velocity");
    int iu = cs_field_get_key_int(f, var_id_key) - 1;
    int iv = iu + 1;
    int iw = iv + 1;

    const cs_turb_model_t *turb_mdl = cs_glob_turb_model;
    int                    ik, iep;
    if (turb_mdl->itytur == 2) {
      f   = cs_field_by_name("k");
      ik  = cs_field_get_key_int(f, var_id_key) - 1;
      f   = cs_field_by_name("epsilon");
      iep = cs_field_get_key_int(f, var_id_key) - 1;
    }

    const int keysca   = cs_field_key_id("scalar_id");
    const int n_fields = cs_field_n_fields();

    /*! [source_term_values] */
    const int nfbpcd = wall_cond->nfbpcd;
    for (ieltcd = 0; ieltcd < nfbpcd; ieltcd++) {

      // Enthalpy of steam
      cs_lnum_t ifac   = wall_cond->ifbpcd[ieltcd];
      cs_lnum_t iel    = ifabor[ifac];
      cs_real_t tk     = cs_glob_fluid_properties->t0;
      cs_lnum_t ntcabs = cs_get_glob_time_step()->nt_cur;
      if (ntcabs >= 1)
        tk = cvar_h[iel] / cpro_cp[iel];
      cs_real_t hvap = (cp_vap * tk);

      // Source term for momentum
      wall_cond->itypcd[iu * nfbpcd + ieltcd] = 0;
      wall_cond->itypcd[iv * nfbpcd + ieltcd] = 0;
      wall_cond->itypcd[iw * nfbpcd + ieltcd] = 0;

      wall_cond->spcond[iu * nfbpcd + ieltcd] = 0.0;
      wall_cond->spcond[iv * nfbpcd + ieltcd] = 0.0;
      wall_cond->spcond[iw * nfbpcd + ieltcd] = 0.0;

      // Source term for turbulence
      // TODO generalize to all turbulence models
      if (turb_mdl->itytur == 2) {
        wall_cond->itypcd[ik * nfbpcd + ieltcd]  = 0;
        wall_cond->itypcd[iep * nfbpcd + ieltcd] = 0;
        wall_cond->spcond[ik * nfbpcd + ieltcd]  = 0.0;
        wall_cond->spcond[iep * nfbpcd + ieltcd] = 0.0;
      }

      // Source term for scalars
      for (int f_id = 0; f_id < n_fields; f_id++) {
        f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_VARIABLE) {
          int iscal = cs_field_get_key_int(f, keysca);
          if (iscal > 0) {
            int ivar = cs_field_get_key_int(f, var_id_key) - 1;
            wall_cond->itypcd[ivar * nfbpcd + ieltcd] = 1;
            if (f == cs_thermal_model_field()) {
              wall_cond->spcond[ivar * nfbpcd + ieltcd] = hvap; // enthalpy
            }
            else {
              wall_cond->spcond[ivar * nfbpcd + ieltcd]
                = 0.0; // non-condensable species
            }
          }
        }
      }
      /*! [source_term_values] */
    }
  }

  BFT_FREE(lstelt);
};

/*----------------------------------------------------------------------------*/

END_C_DECLS
