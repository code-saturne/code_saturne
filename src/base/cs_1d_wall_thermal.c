/*============================================================================
 * Modelling the thermal wall with a 1D approach.
 *============================================================================*/

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_math.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_lagr.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_1d_wall_thermal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_1d_wall_thermal.c
        Modelling the thermal wall with a 1D approach.

  \struct cs_1d_wall_thermal_local_model_t
          1D wall thermal module parameters and variables
          for each coupled face.

  \var  cs_1d_wall_thermal_local_model_t::nppt1d
        Number of discretisation cells in the 1D wall
        for the nfpt1d boundary faces which are coupled
        with a 1D wall thermal module.
  \var  cs_1d_wall_thermal_local_model_t::iclt1d
        Boundary condition type at the external (pseudo)
        wall:
        - 1: Dirichlet,
        - 2: Flux condition.
  \var  cs_1d_wall_thermal_local_model_t::eppt1d
        Thickness of the 1D wall for the nfpt1d boundary faces
        which are coupled with a 1D wall thermal module.
  \var  cs_1d_wall_thermal_local_model_t::rgpt1d
        Geometry of the pseudo wall mesh (refined as a fluid
        if rgt1d is smaller than 1.
  \var  cs_1d_wall_thermal_local_model_t::tept1d
        External temperature of the pseudo wall in the Dirichlet case.
  \var  cs_1d_wall_thermal_local_model_t::hept1d
        External coefficient of transfer in the pseudo wall
        under Dirichlet conditions, (in \f$W.m^{-2}.K\f$).
  \var  cs_1d_wall_thermal_local_model_t::fept1d
        External heat flux in the pseudo wall under the flux
        conditions (in \f$W.m^{-2}\f$, negative value for energy
        entering the wall).
  \var  cs_1d_wall_thermal_local_model_t::xlmbt1
        Thermal diffusivity.
  \var  cs_1d_wall_thermal_local_model_t::rcpt1d
        Volumetric heat capacity rho C_p of the wall uniform
        in thickness (\f$J.m^{-3}.K^{-1}\f$).
  \var  cs_1d_wall_thermal_local_model_t::dtpt1d
        Physical time step associated with the solved 1D equation
        of the pseudo wall (which can be different from the time
        step of the calculation).
  \var  cs_1d_wall_thermal_local_model_t::z
        Discretization points coordinates.
  \var  cs_1d_wall_thermal_local_model_t::t
        Temperature at each point of discretization.

  \struct cs_1d_wall_thermal_t

  \brief 1D wall thermal module descriptor.

  \var  cs_1d_wall_thermal_t::nfpt1d
        Number of boundary faces which are coupled
        with a 1D wall thermal module
  Zones of t1d, dimensioned with nfabor
        Global number of boundary faces which are coupled with
        a 1D wall thermal module, i.e. sum of nfpt1d over all
        ranks
  \var  cs_1d_wall_thermal_t::nmxt1d
  \var  cs_1d_wall_thermal_t::izft1d
        Zones of t1d, dimensioned with nfabor
  \var  cs_1d_wall_thermal_t::ifpt1d
        Array allowing to mark out the numbers of
        the nfpt1d boundary faces which are coupled with
        a 1D wall. The test on \ref ifpt1d implicitly assumes
        that the array is completed in ascending order
        (i.e ifpt1d[ii]\f$>\f$ifpt1d[jj] if ii\f$>\f$jj.
        This will be the case if the coupled faces are defined
        starting from the unique loop on the boundary faces (as
        in the example, see \ref us_pt1d). If this is not
        the case, contact the development team to short circuit
        the test.
  \var  cs_1d_wall_thermal_t::tppt1d
        Initialization temperature of the wall (uniform in thickness).
        During the calculation, the array stores the temperature
        of the solid at the fluid/solid interface.
        Other than for re-reading a file, \ref tppt1d is not used.
        \ref cs_1d_wall_thermal_local_model_t::nppt1d "nppt1d" ,
        \ref ifpt1d, \ref cs_1d_wall_thermal_local_model_t::rgpt1d "rgpt1d"
        and \ref cs_1d_wall_thermal_local_model_t::eppt1d "eppt1d" are
        compared to data from the follow-up file and they must be identical.
  \var  cs_1d_wall_thermal_t::local_models
        Array of structures containing the 1D wall thermal
        module parameters and variables for each coupled face.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variable
 *============================================================================*/

static cs_1d_wall_thermal_t _1d_wall_thermal =
{
  .nfpt1d = 0,
  .nfpt1t = 0,
  .nmxt1d = 0,
  .izft1d = NULL,
  .ifpt1d = NULL,
  .tppt1d = NULL,
  .local_models = NULL
};

const cs_1d_wall_thermal_t *cs_glob_1d_wall_thermal = &_1d_wall_thermal;

static cs_restart_t *cs_glob_tpar1d_suite = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_1d_wall_thermal_get_pointers(cs_lnum_t     **nfpt1d,
                                  cs_gnum_t     **nfpt1t);

void
cs_f_1d_wall_thermal_get_faces(cs_lnum_t **ifpt1d);

void
cs_f_1d_wall_thermal_get_temp(cs_real_t **tppt1d);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the discretization points coordinates array and
          the temperature at each point of discretization.
 */
/*----------------------------------------------------------------------------*/

static void
_1d_wall_thermal_local_models_init(void)
{
  cs_lnum_t ii;

  /* Computation of nmxt1d */
  for (ii = 0; ii < _1d_wall_thermal.nfpt1d ; ii++) {
    _1d_wall_thermal.nmxt1d = CS_MAX(_1d_wall_thermal.local_models[ii].nppt1d,
                                     _1d_wall_thermal.nmxt1d);
  }

  /* if necessary, sum over all the processors */
  cs_parall_max(1, CS_INT_TYPE, &_1d_wall_thermal.nmxt1d);

  /* Initialization of the number of discretization points in each structure
     Computation of the total number of discretization points */
  cs_lnum_t nb_pts_tot = 0;

  for (ii = 0; ii < _1d_wall_thermal.nfpt1d ; ii++)
    nb_pts_tot += _1d_wall_thermal.local_models[ii].nppt1d;

  /* Allocate the "t" arrays: Temperature in each point of discretization
          and the "z" arrays: Coordonnates of each point of discretization */

  if (_1d_wall_thermal.nfpt1d > 0) {
    BFT_MALLOC(_1d_wall_thermal.local_models->z, 2 * nb_pts_tot, cs_real_t);
    _1d_wall_thermal.local_models->t =   _1d_wall_thermal.local_models->z
                                       + nb_pts_tot;
  }

  for (ii = 1 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
    _1d_wall_thermal.local_models[ii].z
      =   _1d_wall_thermal.local_models[ii-1].z
        + _1d_wall_thermal.local_models[ii-1].nppt1d;
    _1d_wall_thermal.local_models[ii].t
      =   _1d_wall_thermal.local_models[ii-1].t
        + _1d_wall_thermal.local_models[ii-1].nppt1d;
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Get pointers to members of the global 1d wall thermal structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 * \param[out]   nfpt1d   pointer to cs_glob_1d_wall_thermal->nfpt1d
 * \param[out]   nfpt1t   pointer to cs_glob_1d_wall_thermal->nfpt1t
 */
/*----------------------------------------------------------------------------*/

void
cs_f_1d_wall_thermal_get_pointers(cs_lnum_t     **nfpt1d,
                                  cs_gnum_t     **nfpt1t)
{
  *nfpt1d = &(_1d_wall_thermal.nfpt1d);
  *nfpt1t = &(_1d_wall_thermal.nfpt1t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the ifpt1d array for the 1D wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_f_1d_wall_thermal_get_faces(cs_lnum_t **ifpt1d)
{
  *ifpt1d = _1d_wall_thermal.ifpt1d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the tppt1d array for the 1D wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_f_1d_wall_thermal_get_temp(cs_real_t **tppt1d)
{
  *tppt1d = _1d_wall_thermal.tppt1d;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the cs_glob_1d_wall_thermal structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_create(void)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  _1d_wall_thermal.nfpt1d = 0;
  _1d_wall_thermal.nfpt1t = 0;
  _1d_wall_thermal.nmxt1d = 0;

  /* Allocate the izft1d array */
  BFT_MALLOC(_1d_wall_thermal.izft1d, n_b_faces, cs_lnum_t);

  for (cs_lnum_t ifac = 0 ; ifac < n_b_faces ; ifac++) {
    _1d_wall_thermal.izft1d[ifac] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate the array of structures local_models.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_local_models_create(void)
{
  cs_lnum_t ii;

  /* Allocate the ifpt1d array */
  BFT_MALLOC(_1d_wall_thermal.ifpt1d, _1d_wall_thermal.nfpt1d, cs_lnum_t);

  /* Allocate the tppt1d array */
  BFT_MALLOC(_1d_wall_thermal.tppt1d, _1d_wall_thermal.nfpt1d, cs_real_t);

  /* Allocate the local_models member of the cs_glob_1d_wall_thermal
   * structure */
  BFT_MALLOC(_1d_wall_thermal.local_models,
             _1d_wall_thermal.nfpt1d,
             cs_1d_wall_thermal_local_model_t);

  for (ii = 0; ii < _1d_wall_thermal.nfpt1d ; ii++) {
    _1d_wall_thermal.local_models[ii].nppt1d = -999;
    _1d_wall_thermal.local_models[ii].iclt1d = 3;
    _1d_wall_thermal.ifpt1d[ii] = -999;
    _1d_wall_thermal.local_models[ii].eppt1d = -999.;
    _1d_wall_thermal.local_models[ii].rgpt1d = -999.;
    _1d_wall_thermal.tppt1d[ii] = 0.;
    _1d_wall_thermal.local_models[ii].tept1d = 0.;
    _1d_wall_thermal.local_models[ii].hept1d = 1.e30;
    _1d_wall_thermal.local_models[ii].fept1d = 0.;
    _1d_wall_thermal.local_models[ii].xlmbt1 = -999.;
    _1d_wall_thermal.local_models[ii].rcpt1d = -999.;
    _1d_wall_thermal.local_models[ii].dtpt1d = -999.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the 1D mesh for each face and initialize the temperature.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_mesh_create(void)
{
  cs_lnum_t ii, kk;
  cs_real_t m, rr, e, n;
  cs_real_t *zz;

  /* Allocate the global structure: cs_glob_par1d and the number of
     discretization points on each face */
  if (_1d_wall_thermal.nfpt1t > 0)
   _1d_wall_thermal_local_models_init();

  for (ii = 0; ii < _1d_wall_thermal.nfpt1d; ii++) {

    n = _1d_wall_thermal.local_models[ii].nppt1d;
    e = _1d_wall_thermal.local_models[ii].eppt1d;

    /* Initialization of the Temperature */
    for (kk = 0; kk < n; kk++) {
      (_1d_wall_thermal.local_models[ii].t)[kk] = _1d_wall_thermal.tppt1d[ii];
    }

    /* Mesh */
    zz = _1d_wall_thermal.local_models[ii].z;
    rr = _1d_wall_thermal.local_models[ii].rgpt1d;

    /* Regular */
    if (fabs(rr-1.0) <= 1.e-6) {
      zz[0] = e/n/2.;
      for (kk = 1 ; kk < n ; kk++) {
        zz[kk] = zz[kk-1] + e/n;
      }
    }

    /* Geometric */
    else {
      m = e*(1.-rr)/(1.-pow(rr,n));
      *zz = m/2.;
      for (kk = 1 ; kk < n ; kk++) {
        zz[kk] = zz[kk-1]+m/2.;
        m = m*rr;
        zz[kk] = zz[kk]+m/2.;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the 1D equation for a given face
 *
 * \param[in]   ii   face number
 * \param[in]   tf   fluid temperature at the boundarys
 * \param[in]   hf   exchange coefficient for the fluid
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_solve(cs_lnum_t ii,
                         cs_real_t tf,
                         cs_real_t hf)
{
  cs_lnum_t kk;
  cs_real_t qinc, eps;
  cs_lnum_t ifac = _1d_wall_thermal.ifpt1d[ii] - 1;

  if (cs_glob_lagr_extra_module->radiative_model >= 1) {
    /* coupling with radiative module, qinc and qeps != 0 */
    /* incident radiative flux at the boundary at the boundary */
    qinc = CS_F_(qinci)->val[ifac];

    /* emissivity */
    eps = CS_F_(emissivity)->val[ifac];
  } else {
    /* without coupling with radiative module */
    qinc = 0.;
    eps = 0.;
  }

  /* thermal diffusivity */
  cs_real_t xlmbt1 = _1d_wall_thermal.local_models[ii].xlmbt1;

  /* volumetric heat capacity of the wall */
  cs_real_t rcp = _1d_wall_thermal.local_models[ii].rcpt1d;

  /* exchange coefficient on the exterior wall */
  cs_real_t hept1d = _1d_wall_thermal.local_models[ii].hept1d;

  /* flux on the exterior wall */
  cs_real_t fept1d = _1d_wall_thermal.local_models[ii].fept1d;

  /* thickness of the 1D wall */
  cs_real_t eppt1d = _1d_wall_thermal.local_models[ii].eppt1d;

  /* temperature on the exterior boundary */
  cs_real_t tept1d = _1d_wall_thermal.local_models[ii].tept1d;

  /* time-step for the solid resolution */
  cs_real_t dtpt1d = _1d_wall_thermal.local_models[ii].dtpt1d;

  /* type of exterior boundary condition */
  int iclt1d = _1d_wall_thermal.local_models[ii].iclt1d;

  cs_real_t a1; /* extrapolation coefficient for temperature1 */
  cs_real_t h2; /* thermal exchange coefficient on T(1) */
  cs_real_t f3; /* thermal flux on Tfluide */
  cs_real_t a4; /* extrapolation coefficient for temperature4 */

  cs_real_t h5 = 0.; /* thermal exchange coefficient on T(n) */
  cs_real_t f6 = 0.; /* thermal flux on Text */

  cs_real_t m;

  cs_real_t _al[128];
  cs_real_t *al, *bl, *cl, *dl;
  cs_real_t *zz;
  cs_lnum_t n;

  n = _1d_wall_thermal.local_models[ii].nppt1d;

  if (n > 32)
    BFT_MALLOC(al, 4*n, cs_real_t);
  else
    al = _al;

  bl = al+n;
  cl = bl+n;
  dl = cl+n;

  zz = _1d_wall_thermal.local_models[ii].z;

  /* Build the tri-diagonal matrix */

  /* Boundary conditions on the fluid side: flux conservation */
  /*   flux in the fluid = flux in the solid = f3 + h2*T1 */

  a1 = 1./hf + zz[0]/xlmbt1;
  h2 = -1./a1; // TAKE CARE TO THE MINUS !
  f3 = -h2*tf + qinc;

  /* Boundary conditions on the exterior */
  /*   flux in the fluid = flux in the solid = f6 + h5*T(n-1) */

  /* Dirichlet condition */
  if (iclt1d == 1) {
    a4 = 1./hept1d + (eppt1d - zz[n-1])/xlmbt1;
    h5 = -1./a4;
    f6 = -h5*tept1d;
  }
  /* Forced flux condition */
  else if (iclt1d == 3) {
    h5 = 0.;
    f6 = fept1d;
  }

  /* Mesh interior points */
  for (kk = 1 ; kk <= n-1; kk++) {
    al[kk] = -xlmbt1/(zz[kk]-zz[kk-1]);
  }

  m = 2*zz[0];
  for (kk = 1 ; kk <= n-2 ; kk++) {
    m = 2*(zz[kk]-zz[kk-1])-m;
    bl[kk] = rcp/dtpt1d*m + xlmbt1/(zz[kk+1]-zz[kk]) + xlmbt1/(zz[kk]-zz[kk-1]);
  }

  for (kk = 0; kk <= n-2; kk++) {
    cl[kk] =  -xlmbt1/(zz[kk+1]-zz[kk]);
  }

  m = 2*zz[0];
  dl[0] = rcp/dtpt1d*m*(_1d_wall_thermal.local_models[ii].t)[0];

  for (kk = 1; kk <= n-1; kk++) {
    m = 2*(zz[kk]-zz[kk-1])-m;
    dl[kk] = rcp/dtpt1d*m*(_1d_wall_thermal.local_models[ii].t)[kk];
  }

  /* Boundary points */
  /* bl[0] and bl[n-1] are initialized here and set later,
     in the case where 0 = n-1 */
  bl[0] = 0.;
  bl[n-1] = 0.;
  al[0] = 0.;
  bl[0] += rcp/dtpt1d*2*zz[0] + xlmbt1/(zz[1]-zz[0]) - h2
         + eps*cs_physical_constants_stephan
         * pow(_1d_wall_thermal.local_models[ii].t[0],3.);
  dl[0] += f3;
  bl[n-1] += rcp/dtpt1d*2*(_1d_wall_thermal.local_models[ii].eppt1d-zz[n-1])
           + xlmbt1/(zz[n-1]-zz[n-2]) - h5;
  cl[n-1] = 0.;
  dl[n-1] += f6;

  /* System resolution by a Cholesky method ("dual-scan") */
  for (kk = 1 ; kk <= n-1 ; kk++) {
    bl[kk] -= al[kk]*cl[kk-1]/bl[kk-1];
    dl[kk] -= al[kk]*dl[kk-1]/bl[kk-1];
  }

  _1d_wall_thermal.local_models[ii].t[n-1] = dl[n-1]/bl[n-1];

  for (kk=n-2; kk>=0; kk--) {
    _1d_wall_thermal.local_models[ii].t[kk] =
       (dl[kk] - cl[kk]*_1d_wall_thermal.local_models[ii].t[kk+1])/bl[kk];
  }

  /* Compute the new value of tp */
  _1d_wall_thermal.tppt1d[ii] = hf + xlmbt1/zz[0];
  _1d_wall_thermal.tppt1d[ii]
    = 1./_1d_wall_thermal.tppt1d[ii]
         *(xlmbt1*_1d_wall_thermal.local_models[ii].t[0]/zz[0]
           + hf*tf);

  if (al != _al)
    BFT_FREE(al);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read the restart file of the 1D-wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_read(void)
{
  bool                corresp_cel, corresp_fac, corresp_fbr, corresp_som;
  cs_lnum_t           nbvent;
  cs_lnum_t           ii, jj, ifac, indfac, ierror;
  cs_lnum_t           n_b_faces = cs_glob_mesh->n_b_faces;

  char nomsui[] = "1dwall_module.csc";

  cs_restart_t             *suite;
  cs_mesh_location_type_t   support;
  cs_restart_val_type_t     typ_val;

  ierror = CS_RESTART_SUCCESS;

  /* Computation of nmxt1d */
  for (ii = 0; ii < _1d_wall_thermal.nfpt1d ; ii++) {
    _1d_wall_thermal.nmxt1d = CS_MAX(_1d_wall_thermal.local_models[ii].nppt1d,
                                     _1d_wall_thermal.nmxt1d);
  }

  /* if necessary, sum over all the processors */
  cs_parall_max(1, CS_INT_TYPE, &_1d_wall_thermal.nmxt1d);

  /* Open the restart file */
  cs_glob_tpar1d_suite = cs_restart_create(nomsui,
                                           NULL,
                                           CS_RESTART_MODE_READ);

  if (cs_glob_tpar1d_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the 1D-wall thermal module restart file "
                "in read mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              nomsui);

  /* Pointer to the global restart structure */
  suite = cs_glob_tpar1d_suite;

  /* Verification of the associated "support" to the restart file */
  cs_restart_check_base_location(suite, &corresp_cel, &corresp_fac,
                                 &corresp_fbr, &corresp_som);

  /* Only boundary faces are of interest */
  indfac = (corresp_fbr == true ? 1 : 0);
  if (indfac == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while reading the 1D-wall thermal module restart file.\n"
                "The number of boundary faces has been modified\n"
                "Verify that the restart file corresponds to "
                "the present study.\n"));


  {
    /* Read the header */
    cs_lnum_t  *tabvar;

    BFT_MALLOC(tabvar, 1, cs_lnum_t);

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_int;

    ierror = cs_restart_read_section(suite,
                                     "version_fichier_suite_module_1d",
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       INCORRECT FILE TYPE\n"
                  "\n"
                  "The file %s does not seem to be a restart file\n"
                  "for the 1D-wall thermal module.\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to a\n"
                  "restart file for the 1D-wall thermal module.\n"),
                nomsui);

    BFT_FREE(tabvar);
  }

  {
    /* Read the number of discretiaztion points and coherency checks
       with the input data of USPT1D */
    cs_lnum_t  *tabvar;
    cs_lnum_t  mfpt1d;
    cs_gnum_t  mfpt1t;
    int  iok;

    BFT_MALLOC(tabvar, n_b_faces, cs_lnum_t);

    const char nomrub[] = "nb_pts_discretis";
    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_int;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Coherency checks between the read NFPT1T and the one from USPT1D */
    mfpt1d = 0;
    for (ifac = 0; ifac < n_b_faces; ifac++) {
      if (tabvar[ifac] > 0) mfpt1d++;
    }
    mfpt1t = mfpt1d;
    /* if necessary, sum over all the processors */
    cs_parall_counter(&mfpt1t, 1);
    if (mfpt1t != _1d_wall_thermal.nfpt1t)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "The number of faces with 1D thermal module has\n"
                  "been modified.\n"
                  "PREVIOUS: %lu boundary faces (total)\n"
                  "CURRENT:  %lu boundary faces (total)\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to a\n"
                  "restart file for the 1D-wall thermal module.\n"
                  "Verify uspt1d.\n"), mfpt1t, _1d_wall_thermal.nfpt1t);

    /* Coherency check between read NFPT1D/IFPT1D and the ones from USPT1D
       One already knows that the number of faces are equal, it is then
       enough to check that all selected faces in USPT1D were also
       selected in the previous calculation */
    iok = 0;
    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d; ii++){
        ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
        if (tabvar[ifac] != _1d_wall_thermal.local_models[ii].nppt1d)
          iok++;
    }
    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "IFPT1D or NPPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file correspond to\n"
                  "the present study.\n"
                  "Verify uspt1d\n"
                  "(refer to the user manual for the specificities\n"
                  "of the test on IFPT1D)"));

    /* Allocate the cs_glob_par1d structure */
    _1d_wall_thermal_local_models_init();

    BFT_FREE(tabvar);
  }

  {
    /* Read the wall thickness and check the coherency with USPT1D*/
    cs_real_t   *tabvar;
    cs_lnum_t   iok;
    cs_real_t eppt1d;

    BFT_MALLOC(tabvar, n_b_faces, cs_real_t);

    const char nomrub[] = "epaisseur_paroi";
    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Coherence check between the read EPPT1D and the one from USPT1D */
    iok = 0;
    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      eppt1d = _1d_wall_thermal.local_models[ii].eppt1d;
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      if (fabs(tabvar[ifac] - eppt1d)/eppt1d > 1.e-10) iok++;
    }
    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "The parameter EPPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to\n"
                  "the present study.\n"
                  "Verify uspt1d\n"));

    for (ii = 0; ii < _1d_wall_thermal.nfpt1d; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      _1d_wall_thermal.local_models[ii].eppt1d = tabvar[ifac];
    }

    BFT_FREE(tabvar);
  }

  {
    /* Read the interior boundary temperature */
    cs_real_t  *tabvar;

    BFT_MALLOC(tabvar, n_b_faces, cs_real_t);

    const char nomrub[] = "temperature_bord_int";
    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      _1d_wall_thermal.tppt1d[ii] = tabvar[ifac];
    }

    BFT_FREE(tabvar);
  }

  { /* Read the 1D-mesh coordinates */
    cs_lnum_t   nptmx;
    cs_lnum_t   iok;
    cs_real_t   *tabvar;
    cs_real_t   zz1, zz2, rrgpt1, rgpt1d;

    nptmx = n_b_faces * _1d_wall_thermal.nmxt1d;
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    const char nomrub[] = "coords_maillages_1d";
    nbvent  = _1d_wall_thermal.nmxt1d;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Now one have the cell centers, RGPT1D can be tested */
    iok = 0;
    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      if (_1d_wall_thermal.local_models[ii].nppt1d > 1) {
        zz1 = tabvar[0 + _1d_wall_thermal.nmxt1d*ifac];
        zz2 = tabvar[1 + _1d_wall_thermal.nmxt1d*ifac];
        rrgpt1 = (zz2-2.*zz1)/zz1;
        rgpt1d = _1d_wall_thermal.local_models[ii].rgpt1d;
        if (fabs(rrgpt1-rgpt1d)/rgpt1d > 1.e-10) iok++;
      }
    }

    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND OLD DATA ARE DIFFERENT\n"
                  "\n"
                  "The parameter RGPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                   "Verify that the restart file correspond to\n"
                  "the present study\n"
                  "Verify uspt1d\n"));

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      /* The array is filled until the number of discretization points of
         the given face is reached */
      for (jj = 0; jj < _1d_wall_thermal.local_models[ii].nppt1d ; jj++)
        _1d_wall_thermal.local_models[ii].z[jj]
          = tabvar[jj + _1d_wall_thermal.nmxt1d*ifac];
    }

    BFT_FREE(tabvar);
  }

  { /* Read the wall temperature */
    cs_lnum_t   nptmx;
    cs_real_t   *tabvar;

    nptmx = n_b_faces * _1d_wall_thermal.nmxt1d;
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    const char nomrub[] = "temperature_interne";
    nbvent  = _1d_wall_thermal.nmxt1d;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_restart_read_section(suite,
                                     nomrub,
                                     support,
                                     nbvent,
                                     typ_val,
                                     tabvar);

    if (ierror < CS_RESTART_SUCCESS) {
      cs_base_warn(__FILE__,__LINE__);
      bft_printf(_("Problem while reading the section in the restart file\n"
                   "for the 1D-wall thermal module:\n"
                   "<%s>\n"), nomrub);
    }

    for (ii = 0; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached */
      for (jj = 0 ; jj < _1d_wall_thermal.local_models[ii].nppt1d ; jj++)
        _1d_wall_thermal.local_models[ii].t[jj]
          = tabvar[jj + _1d_wall_thermal.nmxt1d*ifac];

    }

    BFT_FREE(tabvar);
  }

  cs_restart_read_fields(suite, CS_RESTART_1D_WALL_THERMAL);

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_tpar1d_suite);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write the restart file of the 1D-wall thermal module.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_write(void)
{
  cs_lnum_t            nbvent;
  cs_lnum_t            ii, jj, ifac;

  char nomsui[] = "1dwall_module.csc";

  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_restart_t             *suite;
  cs_mesh_location_type_t   support;
  cs_restart_val_type_t     typ_val;

  /* Open the restart file */
  cs_glob_tpar1d_suite = cs_restart_create(nomsui,
                                           NULL,
                                           CS_RESTART_MODE_WRITE);

  if (cs_glob_tpar1d_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the 1D-wall thermal module restart "
                "file in write mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              nomsui);

  /* Pointer to the global restart structure */
  suite = cs_glob_tpar1d_suite;

  {
    /* Write the header */
    cs_lnum_t  tabvar[1] = {120};

    nbvent  = 1;
    support = CS_MESH_LOCATION_NONE;
    typ_val = CS_TYPE_int;

    cs_restart_write_section(suite,
                             "version_fichier_suite_module_1d",
                             support,
                             nbvent,
                             typ_val,
                             tabvar);
  }

  {
    /* Write the number of discretization points */
    cs_lnum_t  *tabvar;

    BFT_MALLOC(tabvar, n_b_faces, cs_lnum_t);

    for (ii = 0 ; ii < n_b_faces ; ii++)
      tabvar[ii] = 0;

    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_int;

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      tabvar[ifac] =_1d_wall_thermal.local_models[ii].nppt1d;
    }

    cs_restart_write_section(suite,
                             "nb_pts_discretis",
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  {
    /* Write the wall thickness */
    cs_real_t   *tabvar;

    BFT_MALLOC(tabvar, n_b_faces, cs_real_t);

    for (ii = 0 ; ii < n_b_faces ; ii++)
      tabvar[ii] = 0.;

    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      tabvar[ifac] = _1d_wall_thermal.local_models[ii].eppt1d;
    }

    cs_restart_write_section(suite,
                             "epaisseur_paroi",
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  {
    /* Write the internal wall-boundary temperature */
    cs_real_t  *tabvar;

    BFT_MALLOC(tabvar, n_b_faces, cs_real_t);

    for (ii = 0 ; ii < n_b_faces ; ii++)
      tabvar[ii] = 0.0;

    nbvent  = 1;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;
      tabvar[ifac] = _1d_wall_thermal.tppt1d[ii];
    }

    cs_restart_write_section(suite,
                             "temperature_bord_int",
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  {
    /* Write the 1D-mesh coordinates */
    cs_lnum_t   nptmx;
    cs_real_t   *tabvar;

    nptmx = n_b_faces * _1d_wall_thermal.nmxt1d;
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    for (ii = 0 ; ii < nptmx ; ii++)
      tabvar[ii] = 0.;

    nbvent  = _1d_wall_thermal.nmxt1d;
    support = CS_MESH_LOCATION_BOUNDARY_FACES;
    typ_val = CS_TYPE_cs_real_t;

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached (the following cases, up to nmxt1d, are
         already set to 0 during the initalization of tabvar */
      for (jj = 0 ; jj < _1d_wall_thermal.local_models[ii].nppt1d ; jj++)
        tabvar[jj + _1d_wall_thermal.nmxt1d*ifac]
          = _1d_wall_thermal.local_models[ii].z[jj];
    }

    cs_restart_write_section(suite,
                             "coords_maillages_1d",
                             support,
                             nbvent,
                             typ_val,
                             tabvar);

    BFT_FREE(tabvar);
  }

  {
    /* Write the wall-interior temperature */
    cs_lnum_t   nptmx;
    cs_real_t   *tabvar;

    nptmx = n_b_faces * _1d_wall_thermal.nmxt1d;
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    for (ii = 0 ; ii < nptmx ; ii++)
      tabvar[ii] = 0.;

    for (ii = 0 ; ii < _1d_wall_thermal.nfpt1d ; ii++) {
      ifac = _1d_wall_thermal.ifpt1d[ii] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached (the following cases, up to nmxt1d, are
         already set to 0 during the initalization of tabvar */
      for (jj = 0 ; jj < _1d_wall_thermal.local_models[ii].nppt1d ; jj++)
        tabvar[jj + _1d_wall_thermal.nmxt1d*ifac]
          = _1d_wall_thermal.local_models[ii].t[jj];

    }

    cs_restart_write_section(suite,
                             "temperature_interne",
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             _1d_wall_thermal.nmxt1d,
                             CS_TYPE_cs_real_t,
                             tabvar);

    BFT_FREE(tabvar);
  }

  cs_restart_write_fields(suite, CS_RESTART_1D_WALL_THERMAL);

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_tpar1d_suite);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the array of structures local_models.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_free(void)
{
  if (_1d_wall_thermal.local_models != NULL)
    BFT_FREE(_1d_wall_thermal.local_models->z);
  BFT_FREE(_1d_wall_thermal.local_models);
  BFT_FREE(_1d_wall_thermal.ifpt1d);
  BFT_FREE(_1d_wall_thermal.tppt1d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy the global 1d wall thermal structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_finalize(void)
{
  BFT_FREE(_1d_wall_thermal.izft1d);
  cs_glob_1d_wall_thermal = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_1d_wall_thermal.
 */
/*----------------------------------------------------------------------------*/

cs_1d_wall_thermal_t *
cs_get_glob_1d_wall_thermal(void)
{
  return &_1d_wall_thermal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information about the 1d wall thermal computation
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_log(void)
{
  // TODO separate min and max search per zone
  cs_real_t Tp_f_min =  cs_math_big_r;
  cs_real_t Tp_f_max = -cs_math_big_r;
  cs_real_t Tp_ext_min =  cs_math_big_r;
  cs_real_t Tp_ext_max = -cs_math_big_r;

  for (cs_lnum_t ii = 0; ii < _1d_wall_thermal.nfpt1d; ii++) {
    cs_real_t Tp_f = (_1d_wall_thermal.local_models[ii].t)[0];
    cs_lnum_t nppt1d = _1d_wall_thermal.local_models[ii].nppt1d;
    cs_real_t Tp_ext  = (_1d_wall_thermal.local_models[ii].t)[nppt1d-1];
    Tp_f_min = cs_math_fmin(Tp_f_min, Tp_f);
    Tp_f_max = cs_math_fmax(Tp_f_max, Tp_f);
    Tp_ext_min = cs_math_fmin(Tp_ext_min, Tp_ext);
    Tp_ext_max = cs_math_fmax(Tp_ext_max, Tp_ext);
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_min(1, CS_DOUBLE, &Tp_f_min);
    cs_parall_max(1, CS_DOUBLE, &Tp_f_max);
    cs_parall_min(1, CS_DOUBLE, &Tp_ext_min);
    cs_parall_max(1, CS_DOUBLE, &Tp_ext_max);
  }

  bft_printf("   ================================\n");
  bft_printf("    1-D wall thermal resolution\n");
  bft_printf("   ================================\n");
  bft_printf("   Minmax temperature at fluid side    : %15.12e    %15.12e\n",
      Tp_f_min, Tp_f_max);
  bft_printf("   Minmax temperature at external side : %15.12e    %15.12e\n",
      Tp_ext_min, Tp_ext_max);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
