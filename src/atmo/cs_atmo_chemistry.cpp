/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "atmo/cs_air_props.h"
#include "base/cs_base.h"
#include "base/cs_boundary_conditions.h"
#include "cdo/cs_domain.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters_check.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_thermal_model.h"
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo.h"
#include "atmo/cs_atmo_aerosol_ssh.h"
#include "atmo/cs_atmo_chemistry.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* atmo chemistry options structure */
static cs_atmo_chemistry_t _atmo_chem = {
  .model = 0,
  .n_species = 0,
  .n_reactions = 0,
  .chemistry_sep_mode = 2,
  .chemistry_with_photolysis = true,
  .aerosol_model = CS_ATMO_AEROSOL_OFF,
  .frozen_gas_chem = false,
  .init_gas_with_lib = false,
  .init_aero_with_lib = false,
  .n_layer = 0,
  .n_size = 0,
  .spack_file_name = nullptr,
  .species_to_scalar_id = nullptr,
  .species_to_field_id = nullptr,
  .species_profiles_to_field_id = nullptr,
  .molar_mass = nullptr,
  .chempoint = nullptr,
  .conv_factor_jac = nullptr,
  .reacnum = nullptr,
  .dlconc0 = nullptr,
  .aero_file_name = nullptr,
  .chem_conc_file_name = nullptr,
  .aero_conc_file_name = nullptr,
  .nt_step_profiles   = 0,
  .n_z_profiles       = 0,
  .n_species_profiles = 0,
  .conc_profiles   = nullptr,
  .z_conc_profiles = nullptr,
  .t_conc_profiles = nullptr,
  .x_conc_profiles = nullptr,
  .y_conc_profiles = nullptr
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_chemistry_t *cs_glob_atmo_chemistry = &_atmo_chem;

static const char *cs_atmo_aerosol_type_enum_name[]
  = {"CS_ATMO_AEROSOL_OFF",
     "CS_ATMO_AEROSOL_SSH"};

static const char *cs_atmo_aerosol_type_name[]
  = {N_("No atmospheric aerosol"),
     N_("Atmospheric aerosol using external code SSH-aerosol")};

static int _init_atmo_chemistry = 1;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_atmo_get_chem_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len);

void
cs_f_atmo_get_aero_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len);

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint);

void
cs_f_atmo_chem_get_pointers(int                    **model,
                            int                    **isepchemistry,
                            int                    **n_species,
                            int                    **n_reactions,
                            bool                   **chemistry_with_photolysis,
                            cs_atmo_aerosol_type_t **aerosol_model,
                            bool                   **frozen_gas_chem,
                            bool                   **init_gas_with_lib,
                            bool                   **init_aero_with_lib,
                            int                    **n_layer,
                            int                    **n_size);
void
cs_f_atmo_chem_initialize_reacnum(cs_real_t **reacnum);

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid);

void
cs_f_atmo_chem_initialize_species_profiles_to_fid(int *species_profiles_fid);

void
cs_f_atmo_chem_finalize(void);

void
cs_f_atmo_get_chem_conc_profiles(int **nbchim,
                                 int **nbchmz,
                                 int **nespgi);
void
cs_f_atmo_get_arrays_chem_conc_profiles(cs_real_t **espnum,
                                        cs_real_t **zproc,
                                        cs_real_t **tchem,
                                        cs_real_t **xchem,
                                        cs_real_t **ychem,
                                        cs_real_t **conv_factor_jac);

void
cs_f_read_aerosol(void);

void
cs_f_ssh_dimensions(int  *spack_n_species,
                    int  *n_reactions,
                    int  *n_photolysis);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_atmo_chem_get_pointers(int                    **model,
                            int                    **isepchemistry,
                            int                    **n_species,
                            int                    **n_reactions,
                            bool                   **chemistry_with_photolysis,
                            cs_atmo_aerosol_type_t **aerosol_model,
                            bool                   **frozen_gas_chem,
                            bool                   **init_gas_with_lib,
                            bool                   **init_aero_with_lib,
                            int                    **n_layer,
                            int                    **n_size)
{
  *model = &(_atmo_chem.model);
  *isepchemistry = &(_atmo_chem.chemistry_sep_mode);
  *n_species = &(_atmo_chem.n_species);
  *n_reactions = &(_atmo_chem.n_reactions);
  *chemistry_with_photolysis = &(_atmo_chem.chemistry_with_photolysis);
  *aerosol_model = &(_atmo_chem.aerosol_model);
  *frozen_gas_chem = &(_atmo_chem.frozen_gas_chem);
  *init_gas_with_lib = &(_atmo_chem.init_gas_with_lib);
  *init_aero_with_lib = &(_atmo_chem.init_aero_with_lib);
  *n_layer = &(_atmo_chem.n_layer);
  *n_size = &(_atmo_chem.n_size);
}

/*----------------------------------------------------------------------------
 * Return the name of the chemistry concentration file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_chem_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len)
{
  *name = _atmo_chem.chem_conc_file_name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving chemistry concentration file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       _atmo_chem.chem_conc_file_name, name_max, *name, *name_len);
  }
}

/*----------------------------------------------------------------------------
 * Return the name of the aerosol concentration file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_aero_conc_file_name(int           name_max,
                                  const char  **name,
                                  int          *name_len)
{
  *name = _atmo_chem.aero_conc_file_name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving chemistry concentration file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       _atmo_chem.aero_conc_file_name, name_max, *name, *name_len);
  }
}

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint)
{
  *species_to_scalar_id = (_atmo_chem.species_to_scalar_id);
  *molar_mass = (_atmo_chem.molar_mass);
  *chempoint = (_atmo_chem.chempoint);
}

void
cs_f_atmo_chem_initialize_reacnum(cs_real_t **reacnum)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  if (_atmo_chem.reacnum == nullptr)
    BFT_MALLOC(_atmo_chem.reacnum, _atmo_chem.n_reactions*n_cells, cs_real_t);

  *reacnum = _atmo_chem.reacnum;
}

void
cs_f_atmo_chem_initialize_species_to_fid(int *species_fid)
{
  assert(species_fid != nullptr);
  assert(_atmo_chem.species_to_field_id != nullptr);

  for (int i = 0; i < _atmo_chem.n_species; i++)
    _atmo_chem.species_to_field_id[i] = species_fid[i];
}

void
cs_f_atmo_chem_initialize_species_profiles_to_fid(int *species_profiles_fid)
{
  if (_atmo_chem.species_profiles_to_field_id == nullptr)
    BFT_MALLOC(_atmo_chem.species_profiles_to_field_id,
               _atmo_chem.n_species_profiles, int);

  for (int i = 0; i < _atmo_chem.n_species_profiles; i++)
    _atmo_chem.species_profiles_to_field_id[i] = species_profiles_fid[i];
}

void
cs_f_atmo_get_chem_conc_profiles(int **nbchim,
                                 int **nbchmz,
                                 int **nespgi)
{
  *nbchmz = &(_atmo_chem.n_z_profiles);
  *nbchim = &(_atmo_chem.nt_step_profiles);
  *nespgi = &(_atmo_chem.n_species_profiles);
}

void
cs_f_atmo_get_arrays_chem_conc_profiles(cs_real_t **espnum,
                                        cs_real_t **zproc,
                                        cs_real_t **tchem,
                                        cs_real_t **xchem,
                                        cs_real_t **ychem,
                                        cs_real_t **conv_factor_jac)
{
  const int nespg = _atmo_chem.n_species;
  const int _nbchmz = _atmo_chem.n_z_profiles;
  const int _nbchim = _atmo_chem.nt_step_profiles;
  const int size = _atmo_chem.n_species*_nbchmz*_nbchim;

  if (_atmo_chem.conc_profiles == nullptr)
    BFT_MALLOC(_atmo_chem.conc_profiles, size, cs_real_t);

  if (_atmo_chem.z_conc_profiles == nullptr)
    BFT_MALLOC(_atmo_chem.z_conc_profiles, _nbchmz, cs_real_t);

  if (_atmo_chem.t_conc_profiles == nullptr)
    BFT_MALLOC(_atmo_chem.t_conc_profiles, _nbchim, cs_real_t);

  if (_atmo_chem.x_conc_profiles == nullptr)
    BFT_MALLOC(_atmo_chem.x_conc_profiles, _nbchim, cs_real_t);

  if (_atmo_chem.y_conc_profiles == nullptr)
    BFT_MALLOC(_atmo_chem.y_conc_profiles, _nbchim, cs_real_t);

  if (_atmo_chem.conv_factor_jac == nullptr)
    BFT_MALLOC(_atmo_chem.conv_factor_jac, nespg*nespg, cs_real_t);

  *espnum = _atmo_chem.conc_profiles;
  *zproc  = _atmo_chem.z_conc_profiles;
  *tchem  = _atmo_chem.t_conc_profiles;
  *xchem  = _atmo_chem.x_conc_profiles;
  *ychem  = _atmo_chem.y_conc_profiles;
  *conv_factor_jac = _atmo_chem.conv_factor_jac;
}

void
cs_f_atmo_chem_finalize(void)
{
  if (_atmo_chem.aerosol_model != CS_ATMO_AEROSOL_OFF)
    cs_atmo_aerosol_finalize();

  BFT_FREE(_atmo_chem.reacnum);
  BFT_FREE(_atmo_chem.species_to_scalar_id);
  BFT_FREE(_atmo_chem.species_to_field_id);
  BFT_FREE(_atmo_chem.molar_mass);
  BFT_FREE(_atmo_chem.chempoint);
  BFT_FREE(_atmo_chem.conv_factor_jac);
  BFT_FREE(_atmo_chem.spack_file_name);
  BFT_FREE(_atmo_chem.aero_file_name);
  BFT_FREE(_atmo_chem.chem_conc_file_name);
  BFT_FREE(_atmo_chem.conc_profiles);
  BFT_FREE(_atmo_chem.z_conc_profiles);
  BFT_FREE(_atmo_chem.t_conc_profiles);
  BFT_FREE(_atmo_chem.x_conc_profiles);
  BFT_FREE(_atmo_chem.y_conc_profiles);
  BFT_FREE(_atmo_chem.dlconc0);
  BFT_FREE(_atmo_chem.species_profiles_to_field_id);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Beta function: beta(x,y) = gamma(x)*gamma(y)/gamma(x+y)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_beta(const cs_real_t x,
      const cs_real_t y)

{
  cs_real_t beta;

  beta = tgamma(x)*tgamma(y)/tgamma(x+y);

  return beta;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute hypergeometric function for |x| < 1 using a series
 *        (cf. for the definition of this function, see for example:
 *        http://mathworld.wolfram.com/hypergeometricfunction.html)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_hypser(cs_real_t  a,
        cs_real_t  b,
        cs_real_t  c,
        cs_real_t  x)
{
  cs_real_t hypser;

  const int maxiter = 10000;
  const cs_real_t error = 1.e-08;

  if (cs::abs(x) >= 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s\n"
              "The x parameter should verify |x| < 1,  x = %12.5e",
              __func__, x);

  cs_real_t fac = 1;
  cs_real_t temp = fac;
  cs_real_t aa  = a;
  cs_real_t bb  = b;
  cs_real_t cc  = c;

  for (int nn = 1; nn < maxiter; nn++) {
    fac = ((aa*bb)/cc)*fac;
    fac = fac*x/nn;
    hypser = fac + temp;

    if (cs::abs(hypser - temp) <= error)
      return hypser;

    aa += 1;
    bb += 1;
    cc += 1;
    temp = hypser;
  }

  return hypser;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Hypergeometric function
 *        (see http://mathworld.wolfram.com/hypergeometricfunction.html
 *        for definition)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_hypgeo(cs_real_t  a,
        cs_real_t  b,
        cs_real_t  c,
        cs_real_t  x)
{
  cs_real_t hypgeo;

  const cs_real_t pp = 0.1;

  /* Initialization */

  const cs_real_t gammaa   = tgamma(a);
  const cs_real_t gammab   = tgamma(b);
  const cs_real_t gammac   = tgamma(c);
  const cs_real_t gammabma = tgamma(b-a);
  const cs_real_t gammacma = tgamma(c-a);
  const cs_real_t gammaamb = tgamma(a-b);
  const cs_real_t gammacmb = tgamma(c-b);

  /* Compute hypergeometric function by convergent series for |x|<1 */

  if (x >= (pp - 1)) {
    hypgeo = _hypser(a, b, c, x);
                     }
  else if (x >= -1-pp) {
    cs_real_t y1 = _hypser(a, a+1.-c, a+1.-b, 1./x);
    cs_real_t y2 = _hypser(b, b+1.-c, b+1.-a, 1./x);
    hypgeo = (gammac*gammabma*y1*pow(-1./x, a))/(gammab*gammacma)
           + (gammac*gammaamb*y2*pow(-1./x, b))/(gammaa*gammacmb);
  }
  else {
    cs_real_t y1 = _hypser(a, a+1.-c, a+1.-b, 1./(-1.-pp));
    cs_real_t y2 = _hypser(b, b+1.-c, b+1.-a, 1./(-1.-pp));
    cs_real_t hyp1 = (gammac*gammabma*y1*pow(-1./(-1.-pp), a))/(gammab*gammacma)
                   + (gammac*gammaamb*y2*pow(-1./(-1.-pp), b))/(gammaa*gammacmb);
    cs_real_t hyp2 = _hypser(a, b, c, -1.+pp);
    hypgeo = hyp1 + (x - (-1.-pp))*(hyp2 - hyp1)/(2.*pp);
  }

  return hypgeo;
}

/*----------------------------------------------------------------------------
 * Convert string to lower case
 *----------------------------------------------------------------------------*/

static void
_strtolower(char        *dest,
            const char  *src)
{
  char *_dest = dest;
  while (*src) {
    *_dest = tolower(*src);
    src++;
    _dest++;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_initialize();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_finalize();

  BFT_FREE(cs_glob_atmo_chemistry->aero_conc_file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 *
 * \param[out]  array  gas concentrations
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t  *array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_get_gas(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF);

  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH)
    cs_atmo_aerosol_ssh_time_advance();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute aerosol cloud droplets nucleation when using the atmospheric
 * humid model using a microphysical model.
 *
 * It is taken into account as an additional step split from advection-diffusion
 * equation, hence the droplet number is first clipped if necessary.
 *
 * \param[out]  nc      droplet number (scalar) in 1/cm**3
 * \param[in]   rom     density of air in kg/m**3
 * \param[in]   qldia   mass fraction of liquid water in kg/kg
 * \param[in]   pphy    true pressure in pascals
 * \param[in]   refrad  radiative cooling
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_nuclea(cs_real_t         *nc,
                       const cs_real_t   *rom,
                       const cs_real_t   *qldia,
                       const cs_real_t   *pphy,
                       const cs_real_t   *refrad)
{
  /* Initialization
   * -------------- */

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;

  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;
  const cs_real_t *tempc = cs_field_by_name("real_temperature")->val;

  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t cp = phys_pro->cp0;
  cs_real_t rvsra = phys_pro->rvsra;
  cs_real_t clatev = phys_pro->clatev;
  cs_real_t rair = phys_pro->r_pg_cnst;

  const cs_real_t rhowater = 1000; // kg/m^3
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t nuc = 0.;
  cs_real_t fbeta = 0.;
  cs_real_t constc = 0.;
  cs_real_t constk = 0.;
  cs_real_t constmu = 0;
  cs_real_t constbeta = 0;
  cs_real_t rvap = rair*rvsra;
  const int modnuc = cs_glob_atmo_option->nucleation_model;

  if (modnuc == 1) {
    /* Constants for the model of Pruppacher and Klett (1997)
     * Case of Buffalo, New-York, Kocmond [1965] */
    constc = 3500.;
    constk = 0.9;

    /* fbeta : beta function (constk/2, 3/2) */
    fbeta = _beta(constk/2., 1.5);
  }
  else if (modnuc == 2) {
    /* Constants for model of Cohard and Pinty (1998)
     * (general case) */
    constc    = 3270.;
    constk    = 1.56;
    constmu   = 0.7;
    constbeta = 136.;

    /* Polluted Air Case. See Hudson and Li (1995)
     * constc    = 1865. ! 8700.d0 if Cohard and Pinty used on PARISFOG
     * constk    = 0.86
     * constmu   = 1.50
     * constbeta = 6.80 */
    fbeta = _beta(constk/2., 1.5);
  }

  /* Computation new Nc field (in cm^-3)
   * ----------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (qldia[c_id] <= cs_math_epzero) {

      nc[c_id] = 0.;

    }
    else { /* qldia[c_id] > cs_math_epzero */

      /* New droplets are created by nucleation if vel_z > 0 or refrad < 0,
       * then if number of new droplets is superior to nc,
       * the difference is added */

      if ((vel[c_id][2] > cs_math_epzero) || (refrad[c_id] < cs_math_epzero)) {

        cs_real_t tempk = tempc[c_id] + tkelvi;
        cs_real_t esat = cs_air_pwv_sat(tempc[c_id]);
        cs_real_t kka = ((5.69 + 0.017*tempc[c_id])/0.239)*1.e-3;
        cs_real_t ddv = 0.211*pow(tempk/tkelvi, 1.94)*(101325./pphy[c_id])*1.e-4;
        cs_real_t aa1 =   0.622*clatev*9.81/(rair*cp*cs_math_pow2(tempk))
                        - 9.81/(rair*tempk);
        cs_real_t aa2 =   rair*tempk/(0.622*esat)
                        + (0.622*cs_math_pow2(clatev))/(tempk*pphy[c_id]*cp);
        cs_real_t aa3 = 1.0/((1000.*rvap*tempk)/(esat*ddv)
                        + clatev*1000.*(clatev/(tempk*rvap) -1.0)/(kka*tempk));
        cs_real_t aa4 = -0.622*clatev/(rair*cs_math_pow2(tempk));

        if ((aa1*vel[c_id][2]+aa4*refrad[c_id]) > cs_math_epzero) {

          // 1.1  Model of Pruppacher and Klett (1997)
          if (modnuc == 1) {
            const cs_real_t cons2 = pow(constc, 2./(constk+2.));
            const cs_real_t cons3 = pow(0.01*(  aa1*vel[c_id][2]
                                              + aa4*refrad[c_id]), 3./2.);
            const cs_real_t cons4
              = (2.*cs_math_pi*rhowater*aa2*pow(aa3, 3./2.)*constk*fbeta);
            nuc = cons2*pow(cons3/cons4, constk/(constk + 2.));
          }

          // 1.2  Modele of Cohard and Pinty (1998)
          else if (modnuc == 2) {
            /* compute max supersaturation: sursat */
            cs_real_t yy = 1.0;
            cs_real_t sursat = 0;
            cs_real_t tmpsur = 0;
            for (int ii = 0; ii < 20; ii++) {
              tmpsur = sursat;
              yy = _hypgeo(constmu, constk/2, constk/2 + 3./2,
                        -constbeta*cs_math_pow2(sursat));
              const cs_real_t cons1
                = pow(0.01*(aa1*vel[c_id][2] + aa4*refrad[c_id]), 3./2);
              const cs_real_t cons2
                = (2*constk*constc*cs_math_pi*rhowater*aa2*fbeta
                   *pow(aa3, 3.0/2));

              sursat = pow((cons1/cons2)/yy, 1./(constk + 2.));
            }

            if (cs::abs(tmpsur - sursat) > 1.e-2)
              bft_error(__FILE__, __LINE__, 0,
                        "Maximum saturation has not converged\n"
                        "Residue = %12.5e", cs::abs(tmpsur - sursat));

            nuc = constc*pow(sursat, constk)*_hypgeo(constmu,
                                                     constk/2,
                                                     constk/2 + 1.,
                                                     -constbeta*
                                                     cs_math_pow2(sursat));
            if (nuc < 0.)
              bft_error(__FILE__, __LINE__, 0,
                        _("Cohard and Pindy model (1998).\n"
                          "The nucleation is negative."));
          }

          /* 1.3  Model of Abdul-Razzak and al. (1998)
           * This model requires a fine knowledge of aerosols chemistry. */
          else if (modnuc == 3) {

            /* Warning: this set of constants fits for PARISFOG case.
             * They were determined using fine data on aerosols measured
             * during PARISFOG POI 13 */
            cs_real_t sigaero3 = 1.3;
            cs_real_t sigaero2 = 1.59;
            cs_real_t sigaero1 = 1.570;
            cs_real_t ntotal1 = 8700.;
            cs_real_t ntotal2 = 8300.;
            cs_real_t ntotal3 = 1000.;
            cs_real_t rayonm3 = 0.4e-6;
            cs_real_t rayonm2 = 0.055e-6;
            cs_real_t rayonm1 = 0.0165e-6;

            /* surface tension water-steam [N/m], tempk [K] */
            const cs_real_t tauvw = 0.0750;

            /* fractions in %
             * fisrt component: sulfate
             * second component: nitrate
             * third component: black carbon
             * fourth component: organic carbon */
            const cs_real_t fraero[4] = {0.25, 0.2, 0.164, 0.386};

            /* number of ions produced by dissociation of a salt molecule
               in water (if NaCl, nbion=2) */
            const cs_real_t nbion[4] = {3., 2., 0., 1.};

            /* osmotic coefficient */
            const cs_real_t coefosm[4] = {1., 1., 1., 1.};

            /* mass fraction of of soluble
               substance in aerosol mix: msolu/mmix */
            const cs_real_t frasolu[4] = {1., 1., 0., 0.1};

            /* molecular mass of water */
            const cs_real_t mmh2o = 18.e-3;

            /* molecular mass of aerosol components */
            const cs_real_t mmaero[4] = {132e-3, 80e-3, 250e-3, 250e-3};

            /* density of aerosol !FIXME to check */
            const cs_real_t rhoaero[4] = {1.77e3, 1.77e3, 1.77e3, 1.77e3};

            /* coefficient for Kelvin effect: coefa [/m] */
            const cs_real_t coefa = 2*tauvw/(rhowater*rvap*tempk);

            /* FIXME
             * cf Pruppacher and Klett 1997 (reprinted correction 2000) (6-28)
             * coefa = 3.3d-7 / tempk */

            /* coefficient for Raoult effect: coefb [-] */
            cs_real_t numcb = 0;
            cs_real_t dencb = 0;
            for (int ii = 0; ii < 4; ii++) {
              numcb += fraero[ii]*nbion[ii]*coefosm[ii]*frasolu[ii]/mmaero[ii];
              dencb += fraero[ii]/rhoaero[ii];
            }
            const cs_real_t coefb = mmh2o*numcb/(dencb*rhowater);

            /* supersaturation [-] */
            cs_real_t sursat1 = (2./sqrt(coefb))*(pow(coefa/3./rayonm1, 1.5));
            cs_real_t sursat2 = (2./sqrt(coefb))*(pow(coefa/3./rayonm2, 1.5));
            cs_real_t sursat3 = (2./sqrt(coefb))*(pow(coefa/3./rayonm3, 1.5));

            /* standard deviation function */
            cs_real_t foncf1 = 0.5*exp(2.5*cs_math_pow2(log(sigaero1)));
            cs_real_t foncf2 = 0.5*exp(2.5*cs_math_pow2(log(sigaero2)));
            cs_real_t foncf3 = 0.5*exp(2.5*cs_math_pow2(log(sigaero3)));
            cs_real_t foncg1 = 1.0 + 0.25*(log(sigaero1));
            cs_real_t foncg2 = 1.0 + 0.25*(log(sigaero2));
            cs_real_t foncg3 = 1.0 + 0.25*(log(sigaero3));

            /* coefficient corresponding to vertical velocity |aa1| */
            aa1 =   0.622*clatev*9.81/(rair*cp*cs_math_pow2(tempk))
                  - 9.81/(rair*tempk);

            /* coefficient corresponding to liquid water condensation |aa2| */
            aa2 =   rair*tempk/(0.622*esat)+(0.622*cs_math_pow2(clatev))
                  / (tempk*pphy[c_id]*cp);

            /* coefficient corresponding to the droplet growth |aa3| */
            ddv = 0.211*pow(tempk/tkelvi, 1.94)*(101325./pphy[c_id])*1.e-4;
            rvap = rvsra*rair;
            aa3 =   1/((1000*rvap*tempk)/(esat*ddv)
                  + clatev*1000*(clatev/(tempk*rvap)-1.0)/(kka*tempk));

            /* coefficient corresponding to infrared radiation |aa4| */
            aa4 = -0.622*clatev/(rair*cs_math_pow2(tempk));

            /* alphaV/G */
            cs_real_t coefavg = (aa1*vel[c_id][2] + aa4*refrad[c_id])/aa3;

            /* coefficient zeta */
            cs_real_t coefzeta = 2*coefa/3*sqrt(coefavg);

            /* coefficient eta - Ntot [m^(-3)] */
            const cs_real_t coefeta1
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal1*1.e6);
            const cs_real_t coefeta2
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal2*1.e6);
            const cs_real_t coefeta3
              = pow(coefavg, 1.5)/(2*cs_math_pi*rhowater*aa2*ntotal3*1.e6);

            /* max supersaturation (smax) */
            cs_real_t x_tp1 = foncf1*pow(coefzeta/coefeta1, 1.5);
            cs_real_t x_tp2 = pow((sursat1*sursat1)/(coefeta1+3*coefzeta), 0.75);
            cs_real_t smax1 = (x_tp1 + foncg1*x_tp2)/sursat1/sursat1;

            x_tp1 = foncf2*pow(coefzeta/coefeta2, 1.5);
            x_tp2 = pow((sursat2*sursat2)/(coefeta2+3*coefzeta), 0.75);
            cs_real_t smax2 = (x_tp1 + foncg2*x_tp2)/sursat2/sursat2;

            x_tp1 = foncf3*(pow(coefzeta/coefeta3, 1.5));
            x_tp2 = pow((sursat3*sursat3)/(coefeta3+3*coefzeta), 0.75);
            cs_real_t smax3 = (x_tp1 + foncg3*x_tp2)/sursat3/sursat3;

            cs_real_t smax = 1.0/sqrt(smax3 + smax2 + smax1);

            if ((smax > 1.0) || (sursat1 > 1.0) ||
                (sursat2 > 1.0) || (sursat3 > 1.0)) {
              bft_error(__FILE__, __LINE__, 0,
                        _("%s: Abdul-Razzak and al. model (1998).\n"
                          " Negative supersaturation."), __func__);
            }

            cs_real_t f3_ov_sqrt2 = 3.*sqrt(2.);
            cs_real_t ugauss1 =   (2.*log(sursat1/smax))
                                / (f3_ov_sqrt2*log(sigaero1));
            cs_real_t ugauss2 =   (2.*log(sursat2/smax))
                                / (f3_ov_sqrt2*log(sigaero2));
            cs_real_t ugauss3 =   (2.*log(sursat3/smax))
                                / (f3_ov_sqrt2*log(sigaero3));

            const cs_real_t nuc1 = 0.5 * ntotal1 * (1.0-erf(ugauss1));
            const cs_real_t nuc2 = 0.5 * ntotal2 * (1.0-erf(ugauss2));
            const cs_real_t nuc3 = 0.5 * ntotal3 * (1.0-erf(ugauss3));

            nuc = nuc1 + nuc2 + nuc3;

          } /* end nucleation model */

          /* 2. Compute difference */
          nuc = cs::max(nuc - nc[c_id], 0.0);
        }
        else {
          nuc = 0.0;
        }

        nc[c_id] += nuc;
      }

      /* 3. if qc > 0, w <= 0 and nc = 0,
       * we impose nc > 0 so that the mean volume radius = 10 microns */
      else if (nc[c_id] < cs_math_epzero) {

        nc[c_id] =   1.e-6*(3.*rom[c_id]*qldia[c_id])
                   / (4.*cs_math_pi*rhowater*cs_math_pow3(10e-6));

      } /* end Vel_z > cs_math_epzero*/

    } /* end qldia[c_id] > cs_math_epzero */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reads initial aerosol concentration and number
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_read_aerosol(void)
{
  const cs_atmo_chemistry_t *at_chem = cs_glob_atmo_chemistry;
  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  const int n_aer = at_chem->n_size;
  const int nlayer_aer = at_chem->n_layer;
  const int size = n_aer * nlayer_aer;

  cs_atmo_chem_initialize_dlconc0();

  bft_printf("reading of aerosols numbers and concentrations\n");

  if (at_chem->init_aero_with_lib) {

    // The external library provides the concentrations / numbers
    if (at_chem->aerosol_model == CS_ATMO_AEROSOL_SSH)
      cs_atmo_aerosol_ssh_get_aero(at_chem->dlconc0);

    // Conversion from molecules / m^3 to ppm
    const cs_real_t ro0 = 1e3*phys_pro->ro0;
    for (int jb = 0; jb < size; jb++)
      at_chem->dlconc0[jb] = at_chem->dlconc0[jb] / ro0;

    // Conversion from molecules / m^3 to molecules / kg
    const int end_size = size + n_aer;
    for (int jb = size; jb < end_size; jb++)
      at_chem->dlconc0[jb] = at_chem->dlconc0[jb] / phys_pro->ro0;

  }
  else {

    // Read from file
    FILE *file = fopen(at_chem->aero_conc_file_name, "r");

    // Reading aerosol numbers
    for (int jb = 0; jb < n_aer; jb++)
      fscanf(file,"%le", &at_chem->dlconc0[jb + size]);

    // Reading aerosol concentrations
    for (int jb = 0; jb < n_aer; jb++)
      for (int jsp = 0; jsp < nlayer_aer; jsp++)
        fscanf(file,"%le", &at_chem->dlconc0[jb + jsp*n_aer]);

    fclose(file);

  }

  /* Logging
     ------- */

  bft_printf("\n===================================================\n");
  bft_printf("printing aerosol numbers\n");
  for (int jb = 0; jb < n_aer; jb++) {
    const int f_id
      = at_chem->species_to_scalar_id[at_chem->n_species + size + jb];
    const cs_field_t *f = cs_field_by_id(f_id);
    bft_printf("%s : %10.2le\n ",
               cs_field_get_label(f),
               at_chem->dlconc0[jb + size]);
  }

  bft_printf("\n===================================================\n");
  bft_printf("printing aerosol concentrations\n");
  for (int jb = 0; jb < n_aer; jb++) {
    bft_printf("Size bin id %d\n", jb);
    for (int jsp = 0; jsp < nlayer_aer; jsp++) {
      const int f_id
        = at_chem->species_to_scalar_id[at_chem->n_species + jsp*n_aer];
      const cs_field_t *f = cs_field_by_id(f_id);
      bft_printf("%s : %10.2le\n ",
                 cs_field_get_label(f),
                 at_chem->dlconc0[jb + jsp*n_aer]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the chemistry file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_chem_conc_file_name(const char *file_name)
{
  if (file_name == nullptr) {
    return;
  }

  if (_atmo_chem.chem_conc_file_name == nullptr) {
    BFT_MALLOC(_atmo_chem.chem_conc_file_name,
               strlen(file_name) + 1,
               char);
  }
  else {
    BFT_REALLOC(_atmo_chem.chem_conc_file_name,
                strlen(file_name) + 1,
                char);
  }

  sprintf(_atmo_chem.chem_conc_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the aerosol file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_set_aero_conc_file_name(const char *file_name)
{
  if (file_name == nullptr) {
    return;
  }
  if (cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_OFF) {
    return;
  }

  if (_atmo_chem.aero_conc_file_name == nullptr) {
    BFT_MALLOC(_atmo_chem.aero_conc_file_name,
               strlen(file_name) + 1,
               char);
  }
  else {
    BFT_REALLOC(_atmo_chem.aero_conc_file_name,
                strlen(file_name) + 1,
                char);
  }

  sprintf(_atmo_chem.aero_conc_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize chemistry array.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_init_chemistry(void)
{
  /* Initialization of the chemical scheme
     quasi steady equilibrium NOx scheme */
  if (_atmo_chem.model == 1) {

    _atmo_chem.n_species = 4;
    _atmo_chem.n_reactions = 5;

    if (_atmo_chem.chempoint == nullptr)
      BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);

    if (_atmo_chem.molar_mass == nullptr)
      BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);

    _atmo_chem.chempoint[0] = 4;
    _atmo_chem.chempoint[1] = 3;
    _atmo_chem.chempoint[2] = 2;
    _atmo_chem.chempoint[3] = 1;
    _atmo_chem.molar_mass[0] = 30.0;
    _atmo_chem.molar_mass[1] = 46.0;
    _atmo_chem.molar_mass[2] = 48.0;
    _atmo_chem.molar_mass[3] = 16.0;

    if (_atmo_chem.species_to_field_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);

    if (_atmo_chem.species_to_scalar_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);

    const char *label[] = {"NO", "NO2", "O3", "O3P"};

    const char *name[]
      = {"species_no", "species_no2", "species_o3", "species_o3p"};

    const int kscal = cs_field_key_id_try("scalar_id");
    for (int ii = 0; ii < _atmo_chem.n_species; ii++) {
      int f_id = cs_variable_field_create(name[ii],
                                          label[ii],
                                          CS_MESH_LOCATION_CELLS,
                                          1);
      cs_field_t *f = cs_field_by_id(f_id);
      cs_add_model_field_indexes(f->id);
      _atmo_chem.species_to_field_id[ii] = f->id;
      _atmo_chem.species_to_scalar_id[ii] = cs_field_get_key_int(f, kscal);
    }

  }
  else if (_atmo_chem.model == 2) {

    _atmo_chem.n_species = 20;
    _atmo_chem.n_reactions = 34;

    if (_atmo_chem.chempoint == nullptr)
      BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);
    if (_atmo_chem.molar_mass == nullptr)
      BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);

    _atmo_chem.chempoint[0]  = 20, _atmo_chem.chempoint[1]  = 19;
    _atmo_chem.chempoint[2]  = 16, _atmo_chem.chempoint[3]  = 17;
    _atmo_chem.chempoint[4]  = 2,  _atmo_chem.chempoint[5]  = 15;
    _atmo_chem.chempoint[6]  = 14, _atmo_chem.chempoint[7]  = 3;
    _atmo_chem.chempoint[8]  = 18, _atmo_chem.chempoint[9]  = 7;
    _atmo_chem.chempoint[10] = 8,  _atmo_chem.chempoint[11] = 9;
    _atmo_chem.chempoint[12] = 4,  _atmo_chem.chempoint[13] = 10;
    _atmo_chem.chempoint[14] = 1,  _atmo_chem.chempoint[15] = 12;
    _atmo_chem.chempoint[16] = 11, _atmo_chem.chempoint[17] = 13;
    _atmo_chem.chempoint[18] = 5,  _atmo_chem.chempoint[19] = 6;

    _atmo_chem.molar_mass[0]  = 30.000;  // Molar mass (g/mol) NO
    _atmo_chem.molar_mass[1]  = 46.000;  // Molar mass (g/mol) NO2
    _atmo_chem.molar_mass[2]  = 48.000;  // Molar mass (g/mol) O3
    _atmo_chem.molar_mass[3]  = 16.000;  // Molar mass (g/mol) O3P
    _atmo_chem.molar_mass[4]  = 16.000;  // Molar mass (g/mol) O1D
    _atmo_chem.molar_mass[5]  = 17.010;  // Molar mass (g/mol) OH
    _atmo_chem.molar_mass[6]  = 33.010;  // Molar mass (g/mol) HO2
    _atmo_chem.molar_mass[7]  = 34.010;  // Molar mass (g/mol) H2O2
    _atmo_chem.molar_mass[8]  = 62.010;  // Molar mass (g/mol) NO3
    _atmo_chem.molar_mass[9]  = 108.01;  // Molar mass (g/mol) N2O5
    _atmo_chem.molar_mass[10] = 47.010;  // Molar mass (g/mol) HONO
    _atmo_chem.molar_mass[11] = 63.010;  // Molar mass (g/mol) HNO3
    _atmo_chem.molar_mass[12] = 28.010;  // Molar mass (g/mol) CO
    _atmo_chem.molar_mass[13] = 30.030;  // Molar mass (g/mol) HCHO
    _atmo_chem.molar_mass[14] = 44.050;  // Molar mass (g/mol) ALD2
    _atmo_chem.molar_mass[15] = 75.040;  // Molar mass (g/mol) C2O3
    _atmo_chem.molar_mass[16] = 121.05;  // Molar mass (g/mol) PAN
    _atmo_chem.molar_mass[17] = 47.030;  // Molar mass (g/mol) XO2
    _atmo_chem.molar_mass[18] = 64.060;  // Molar mass (g/mol) SO2
    _atmo_chem.molar_mass[19] = 98.080;  // Molar mass (g/mol) H2SO4

    const char *label[] = {"NO",   "NO2",  "O3",   "O3P",  "O1D",
                           "OH",   "HO2",  "H2O2", "NO3",  "N2O5",
                           "HONO", "HNO3", "CO",   "HCHO", "ALD2",
                           "C2O3", "PAN",  "XO2",  "SO2",  "H2SO4"};

    const char *name[]
      = {"species_no",  "species_no2",  "species_o3",  "species_o3p",
         "species_o1d", "species_oh",   "species_ho2", "species_h2o2",
         "species_no3", "species_n2o5", "species_hono", "species_hno3",
         "species_co",  "species_hcho", "species_ald2","species_c2o3",
         "species_pan", "species_xo2",  "species_so2", "species_h2so4"};

    if (_atmo_chem.species_to_field_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
    if (_atmo_chem.species_to_scalar_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);

    const int kscal = cs_field_key_id_try("scalar_id");
    for (int ii = 0; ii < _atmo_chem.n_species; ii++) {
      int f_id = cs_variable_field_create(name[ii],
                                          label[ii],
                                          CS_MESH_LOCATION_CELLS,
                                          1);
      cs_field_t *f = cs_field_by_id(f_id);
      cs_add_model_field_indexes(f->id);
      _atmo_chem.species_to_field_id[ii] = f->id;
      _atmo_chem.species_to_scalar_id[ii] = cs_field_get_key_int(f, kscal);
    }

  }
  else if (_atmo_chem.model == 3) {
    _atmo_chem.n_species = 52;
    _atmo_chem.n_reactions = 155;

    if (_atmo_chem.chempoint == nullptr)
      BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);
    if (_atmo_chem.molar_mass == nullptr)
      BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);

    _atmo_chem.chempoint[0]  = 48, _atmo_chem.chempoint[1]  = 52;
    _atmo_chem.chempoint[2]  = 47, _atmo_chem.chempoint[3]  = 43;
    _atmo_chem.chempoint[4]  = 1,  _atmo_chem.chempoint[5]  = 42;
    _atmo_chem.chempoint[6]  = 50, _atmo_chem.chempoint[7]  = 17;
    _atmo_chem.chempoint[8]  = 44, _atmo_chem.chempoint[9]  = 9;
    _atmo_chem.chempoint[10] = 15, _atmo_chem.chempoint[11] = 38;
    _atmo_chem.chempoint[12] = 13, _atmo_chem.chempoint[13] = 37;
    _atmo_chem.chempoint[14] = 41, _atmo_chem.chempoint[15] = 45;
    _atmo_chem.chempoint[16] = 51, _atmo_chem.chempoint[17] = 10;
    _atmo_chem.chempoint[18] = 35, _atmo_chem.chempoint[19] = 46;
    _atmo_chem.chempoint[20] = 14, _atmo_chem.chempoint[21] = 49;
    _atmo_chem.chempoint[22] = 39, _atmo_chem.chempoint[23] = 33;
    _atmo_chem.chempoint[24] = 2,  _atmo_chem.chempoint[25] = 3;
    _atmo_chem.chempoint[26] = 40, _atmo_chem.chempoint[27] = 11;
    _atmo_chem.chempoint[28] = 19, _atmo_chem.chempoint[29] = 20;
    _atmo_chem.chempoint[30] = 4,  _atmo_chem.chempoint[31] = 21;
    _atmo_chem.chempoint[32] = 36, _atmo_chem.chempoint[33] = 22;
    _atmo_chem.chempoint[34] = 34, _atmo_chem.chempoint[35] = 16;
    _atmo_chem.chempoint[36] = 23, _atmo_chem.chempoint[37] = 24;
    _atmo_chem.chempoint[38] = 25, _atmo_chem.chempoint[39] = 31;
    _atmo_chem.chempoint[40] = 32, _atmo_chem.chempoint[41] = 26;
    _atmo_chem.chempoint[42] = 5,  _atmo_chem.chempoint[43] = 6;
    _atmo_chem.chempoint[44] = 27, _atmo_chem.chempoint[45] = 12;
    _atmo_chem.chempoint[46] = 28, _atmo_chem.chempoint[47] = 30;
    _atmo_chem.chempoint[48] = 29, _atmo_chem.chempoint[49] = 7;
    _atmo_chem.chempoint[50] = 8,  _atmo_chem.chempoint[51] = 18;

    _atmo_chem.molar_mass[0]  = 30.000, _atmo_chem.molar_mass[1]  = 46.000;
    _atmo_chem.molar_mass[2]  = 48.000, _atmo_chem.molar_mass[3]  = 16.000;
    _atmo_chem.molar_mass[4]  = 16.000, _atmo_chem.molar_mass[5]  = 17.010;
    _atmo_chem.molar_mass[6]  = 33.010, _atmo_chem.molar_mass[7]  = 34.010;
    _atmo_chem.molar_mass[8]  = 62.010, _atmo_chem.molar_mass[9]  = 108.01;
    _atmo_chem.molar_mass[10] = 47.010, _atmo_chem.molar_mass[11] = 63.010;
    _atmo_chem.molar_mass[12] = 79.010, _atmo_chem.molar_mass[13] = 28.010;
    _atmo_chem.molar_mass[14] = 30.030, _atmo_chem.molar_mass[15] = 44.050;
    _atmo_chem.molar_mass[16] = 75.040, _atmo_chem.molar_mass[17] = 121.05;
    _atmo_chem.molar_mass[18] = 43.040, _atmo_chem.molar_mass[19] = 74.040;
    _atmo_chem.molar_mass[20] = 120.04, _atmo_chem.molar_mass[21] = 47.030;
    _atmo_chem.molar_mass[22] = 47.030, _atmo_chem.molar_mass[23] = 77.040;
    _atmo_chem.molar_mass[24] = 46.070, _atmo_chem.molar_mass[25] = 16.040;
    _atmo_chem.molar_mass[26] = 47.030, _atmo_chem.molar_mass[27] = 32.040;
    _atmo_chem.molar_mass[28] = 48.040, _atmo_chem.molar_mass[29] = 46.030;
    _atmo_chem.molar_mass[30] = 30.070, _atmo_chem.molar_mass[31] = 47.030;
    _atmo_chem.molar_mass[32] = 60.050, _atmo_chem.molar_mass[33] = 76.050;
    _atmo_chem.molar_mass[34] = 15.030, _atmo_chem.molar_mass[35] = 16.000;
    _atmo_chem.molar_mass[36] = 28.050, _atmo_chem.molar_mass[37] = 27.050;
    _atmo_chem.molar_mass[38] = 56.110, _atmo_chem.molar_mass[39] = 68.180;
    _atmo_chem.molar_mass[40] = 70.090, _atmo_chem.molar_mass[41] = 136.24;
    _atmo_chem.molar_mass[42] = 92.140, _atmo_chem.molar_mass[43] = 106.16;
    _atmo_chem.molar_mass[44] = 108.14, _atmo_chem.molar_mass[45] = 141.15;
    _atmo_chem.molar_mass[46] = 48.040, _atmo_chem.molar_mass[47] = 108.14;
    _atmo_chem.molar_mass[48] = 72.060, _atmo_chem.molar_mass[49] = 64.060;
    _atmo_chem.molar_mass[50] = 98.080, _atmo_chem.molar_mass[51] = 63.030;

    const char *label[] = {"NO",   "NO2",  "O3",    "O3P",
                           "O1D",  "OH",   "HO2",   "H2O2",
                           "NO3",  "N2O5", "HONO",  "HNO3",
                           "HNO4", "CO",   "HCHO" , "ALD2",
                           "C2O3", "PAN",  "ALDX",  "CXO3",
                           "PANX", "XO2",  "XO2N",  "NTR",
                           "ETOH", "CH4",  "MEO2",  "MEOH",
                           "MEPX", "FACD", "ETHA",  "ROOH",
                           "AACD", "PACD", "PAR",   "ROR",
                           "ETH",  "OLE",  "IOLE",  "ISOP",
                           "ISPD", "TERP", "TOL",   "XYL",
                           "CRES", "TO2",  "OPEN",   "CRO",
                           "MGLY", "SO2",  "H2SO4", "HCO3"};
    const char *name[]
      = {"species_no",   "species_no2",  "species_o3",    "species_o3p",
         "species_o1d",  "species_oh",   "species_ho2",   "species_h2o2",
         "species_no3",  "species_n2o5", "species_hono",  "species_hno3",
         "species_hno4", "species_co",   "species_hcho",  "species_ald2",
         "species_c2o3", "species_pan",  "species_aldx",  "species_cxo3",
         "species_panx", "species_xo2",  "species_xo2n",  "species_ntr",
         "species_etoh", "species_ch4",  "species_meo2",  "species_meoh",
         "species_mepx", "species_facd", "species_etha",  "species_rooh",
         "species_aacd", "species_pacd", "species_par",   "species_ror",
         "species_eth",  "species_ole",  "species_iole",  "species_isop",
         "species_ispd", "species_terp", "species_tol",   "species_xyl",
         "species_cres", "species_to2",  "species_open",  "species_cro",
         "species_mgly", "species_so2",  "species_h2so4", "species_hco3"};

    if (_atmo_chem.species_to_field_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
    if (_atmo_chem.species_to_scalar_id == nullptr)
      BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);

    const int kscal = cs_field_key_id_try("scalar_id");
    for (int ii = 0; ii < _atmo_chem.n_species; ii++) {
      int f_id = cs_variable_field_create(name[ii],
                                          label[ii],
                                          CS_MESH_LOCATION_CELLS,
                                          1);
      cs_field_t *f = cs_field_by_id(f_id);
      cs_add_model_field_indexes(f->id);
      _atmo_chem.species_to_field_id[ii] = f->id;
      _atmo_chem.species_to_scalar_id[ii] = cs_field_get_key_int(f, kscal);
    }

  }
  //User defined chemistry using SPACK file and routines
  else if (_atmo_chem.model == 4)  {

    /* This function read the number of species, their molar mass
     * and creates variables */
    cs_atmo_declare_chem_from_spack();

    int  spack_n_species, n_photolysis;
    cs_f_ssh_dimensions(&spack_n_species,
                        &_atmo_chem.n_reactions,
                        &n_photolysis);

    if (spack_n_species != _atmo_chem.n_species)
      cs_parameters_error
        (CS_ABORT_IMMEDIATE,
         _("In atmospheric chemistry inpu data"),
         _("The number of gaseous species read from the SPACK file\n"
           "is not equal to the one read in the SPACK source file."));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function set the file name of the SPACK file.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_spack_file_name(const char *file_name)
{
  if (file_name == nullptr) {
    _atmo_chem.model = 0;
    return;
  }

  _atmo_chem.model = 4;

  BFT_MALLOC(_atmo_chem.spack_file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_atmo_chem.spack_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function sets the file name to initialize the aerosol library.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_set_aerosol_file_name(const char *file_name)
{
  if (file_name == nullptr) {
    _atmo_chem.aerosol_model = CS_ATMO_AEROSOL_OFF;
    return;
  }

  _atmo_chem.aerosol_model = CS_ATMO_AEROSOL_SSH;

  BFT_MALLOC(_atmo_chem.aero_file_name,
             strlen(file_name) + 1,
             char);

  sprintf(_atmo_chem.aero_file_name, "%s", file_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function declares additional transported variables for
 *        atmospheric module for the chemistry defined from SPACK.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_declare_chem_from_spack(void)
{
  assert(_atmo_chem.model == 4);

  if (_atmo_chem.spack_file_name == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Atmo chemistry from SPACK file: requires a SPACK file."));

  char line[512] = "";

  /* Open file */
  bft_printf("SPACK file for atmo chemistry:\n    %s \n",
             _atmo_chem.spack_file_name);

  FILE *file = fopen(_atmo_chem.spack_file_name, "rt");
  if (file == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Atmo chemistry from SPACK file: Could not open file."));

  int line_num = 0;

  /* Read "[species]" */
  while (true) {
    line_num++;
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not skip header."));

    if (strncmp("[species]", line, 512) == 0)
      break;
  }

  /* Read SPACK: first loop count the number of species */
  for (int i = 1; true; i++ ) {
    /* Read species */
    line_num++;
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not read line %d."),
                line_num);

    /* When reach [molecular_waight]: break */
    if (strncmp("[molecular_weight]", line, 512) == 0)
      break;
    else {
      _atmo_chem.n_species = i;
    }
  }

  /* Now allocate arrays */
  BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
  BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);
  BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);
  BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);

  /* Read SPACK: second loop Create variables and read molar mass */
  for (int i = 0; i < _atmo_chem.n_species; i++ ) {
    char name[512] = "";
    char label[512] = "";

    /* Read species name and molar mass */
    line_num++;
    if (fscanf(file, "%s %lf\n", line, &(_atmo_chem.molar_mass[i])) != 2)
      bft_error(__FILE__, __LINE__, 0,
                _("Atmospheric chemistry from SPACK file:\n"
                  "  could not read species name and molar mass, line %d."),
                line_num);

    /* The order is already ok */
    _atmo_chem.chempoint[i] = i+1; //FIXME ?

    /* Build name of the field:
     * species_name in lower case */
    strcpy(name, "species_");
    _strtolower(label, line);
    strcat(name, label);

    /* Field of dimension 1 */
    /* Give the original name as label */
    _atmo_chem.species_to_field_id[i]
      = cs_variable_field_create(name, line, CS_MESH_LOCATION_CELLS, 1);

    /* Scalar field, store in isca_chem/species_to_scalar_id (FORTRAN/C) array */
    _atmo_chem.species_to_scalar_id[i]
      = cs_add_model_field_indexes(_atmo_chem.species_to_field_id[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deactivate chemistry initialization procedure
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_initialization_deactivate(void)
{
  _init_atmo_chemistry = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the chemistry module needs initialization
 *
 * \return int value : 1 if needed, 0 if not
 */
/*----------------------------------------------------------------------------*/

int
cs_atmo_chemistry_need_initialization(void)
{
  return _init_atmo_chemistry;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric chemistry options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chemistry_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 0)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Atmospheric chemistry options\n"
                  "-----------------------------\n\n"));

  if (cs_glob_atmo_chemistry->model == 0) {

    /* No atmospheric chemistry */
    cs_log_printf(CS_LOG_SETUP,
                  _("  No atmospheric chemistry\n\n"));

  }
  else if (   cs_glob_atmo_chemistry->model == 1
           || cs_glob_atmo_chemistry->model == 2
           || cs_glob_atmo_chemistry->model == 3) {

    /* Pre-defined schemes */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric chemistry activated\n\n"
         "    Pre-defined scheme %12d\n\n"
         "      n_species: %18d (Number of species)\n"
         "      n_reactions: %16d (Number of reactions)\n"
         "      photolysis: %17s\n"
         "      frozen_gas_chem: %12s\n\n"),
       cs_glob_atmo_chemistry->model,
       cs_glob_atmo_chemistry->n_species,
       cs_glob_atmo_chemistry->n_reactions,
       cs_glob_atmo_chemistry->chemistry_with_photolysis ? "Yes": "No",
       cs_glob_atmo_chemistry->frozen_gas_chem ? "Yes": "No");

  }
  else if (cs_glob_atmo_chemistry->model == 4) {

    /* User defined SPACK chemistry */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric chemistry activated\n\n"
         "    User-defined SPACK scheme\n\n"
         "      n_species: %18d (Number of species)\n"
         "      n_reactions: %16d (Number of reactions)\n"
         "      photolysis: %17s\n"
         "      frozen_gas_chem: %12s\n"
         "      Spack file: %17s\n"),
       cs_glob_atmo_chemistry->n_species,
       cs_glob_atmo_chemistry->n_reactions,
       cs_glob_atmo_chemistry->chemistry_with_photolysis ? "Yes": "No",
       cs_glob_atmo_chemistry->frozen_gas_chem ? "Yes": "No",
       cs_glob_atmo_chemistry->spack_file_name);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the atmospheric aerosols options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 0)
    return;

  cs_log_printf
    (CS_LOG_SETUP,
     _("\nAtmospheric aerosols options\n"
       "----------------------------\n\n"));

  cs_atmo_aerosol_type_t atm_aer_type = cs_glob_atmo_chemistry->aerosol_model;

  if (atm_aer_type == CS_ATMO_AEROSOL_OFF)
    cs_log_printf(CS_LOG_SETUP, _("  %s\n"),
                  cs_atmo_aerosol_type_name[atm_aer_type]);

  else if (atm_aer_type == CS_ATMO_AEROSOL_SSH) {

    /* Atmospheric chemistry with SSH */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Atmospheric aerosols activated\n\n"
         "    Global parameters\n\n"
         "      model: %22s (%s)\n"
         "      n_layer: %20d (Number of layers inside aerosols)\n"
         "      n_size:  %20d (Number of resolved aerosols)\n"
         "      init_gas_with_lib: %10s\n"
         "      init_aero_with_lib: %9s\n"
         "      namelist: %s\n\n"),
       cs_atmo_aerosol_type_enum_name[atm_aer_type],
       cs_atmo_aerosol_type_name[atm_aer_type],
       cs_glob_atmo_chemistry->n_layer,
       cs_glob_atmo_chemistry->n_size,
       cs_glob_atmo_chemistry->init_gas_with_lib ? "Yes": "No",
       cs_glob_atmo_chemistry->init_aero_with_lib ? "Yes": "No",
       cs_glob_atmo_chemistry->aero_file_name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief initialize gaseous and particulate concentrations and aerosol number
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_chem_initialize_dlconc0(void)
{
  if (_atmo_chem.dlconc0 != nullptr)
    return;

  const int n_aer = _atmo_chem.n_size;
  const int nlayer_aer = _atmo_chem.n_layer;
  const int size = n_aer*(1+nlayer_aer);

  BFT_MALLOC(_atmo_chem.dlconc0, size, cs_real_t);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
