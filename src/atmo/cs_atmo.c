/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* global atmo options structure */
static cs_atmo_option_t  _atmo_option = {
  .syear = -1,
  .squant = -1,
  .shour = -1,
  .smin = -1,
  .ssec = -1.,
  .longitute = 1e12, // TODO use cs_math_big_r
  .latitude = 1e12,
  .compute_z_ground = false
};

/* atmo chemistry options structure */
static cs_atmo_chemistry_t _atmo_chem = {
  .model = 0,
  .n_species = 0,
  .n_reactions = 0,
  .spack_file_name = NULL,
  .species_to_scalar_id = NULL,
  .species_to_field_id = NULL,
  .molar_mass = NULL,
  .chempoint = NULL
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_option_t *cs_glob_atmo_option = &_atmo_option;
cs_atmo_chemistry_t *cs_glob_atmo_chemistry = &_atmo_chem;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_atmo_get_pointers(int       **syear,
                       int       **squant,
                       int       **shour,
                       int       **smin,
                       cs_real_t **ssec,
                       cs_real_t **longitute,
                       cs_real_t **latitude,
                       bool      **compute_z_ground,
                       int       **model,
                       int       **n_species,
                       int       **n_reactions);

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint);
void
cs_f_atmo_chem_finalize(void);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointer
 *----------------------------------------------------------------------------*/

void
cs_f_atmo_get_pointers(int       **syear,
                       int       **squant,
                       int       **shour,
                       int       **smin,
                       cs_real_t **ssec,
                       cs_real_t **longitute,
                       cs_real_t **latitude,
                       bool      **compute_z_ground,
                       int       **model,
                       int       **n_species,
                       int       **n_reactions)
{
  *syear     = &(_atmo_option.syear);
  *squant    = &(_atmo_option.squant);
  *shour     = &(_atmo_option.shour);
  *smin      = &(_atmo_option.smin);
  *ssec      = &(_atmo_option.ssec);
  *longitute = &(_atmo_option.longitute);
  *latitude  = &(_atmo_option.latitude);
  *compute_z_ground = &(_atmo_option.compute_z_ground);
  *model = &(_atmo_chem.model);
  *n_species = &(_atmo_chem.n_species);
  *n_reactions = &(_atmo_chem.n_reactions);
}

void
cs_f_atmo_chem_arrays_get_pointers(int       **species_to_scalar_id,
                                   cs_real_t **molar_mass,
                                   int       **chempoint)
{

  if (_atmo_chem.species_to_scalar_id == NULL)
    BFT_MALLOC(_atmo_chem.species_to_scalar_id, _atmo_chem.n_species, int);
  if (_atmo_chem.species_to_field_id == NULL)
    BFT_MALLOC(_atmo_chem.species_to_field_id, _atmo_chem.n_species, int);
  if (_atmo_chem.molar_mass == NULL)
    BFT_MALLOC(_atmo_chem.molar_mass, _atmo_chem.n_species, cs_real_t);
  if (_atmo_chem.chempoint == NULL)
    BFT_MALLOC(_atmo_chem.chempoint, _atmo_chem.n_species, int);

  *species_to_scalar_id = (_atmo_chem.species_to_scalar_id);
  *molar_mass = (_atmo_chem.molar_mass);
  *chempoint = (_atmo_chem.chempoint);
}

void
cs_f_atmo_chem_finalize(void)
{
  BFT_FREE(_atmo_chem.species_to_scalar_id);
  BFT_FREE(_atmo_chem.species_to_field_id);
  BFT_FREE(_atmo_chem.molar_mass);
  BFT_FREE(_atmo_chem.chempoint);
  BFT_FREE(_atmo_chem.spack_file_name);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
 * \brief This function computes the ground elevation
 *
 *  This function solves the following transport equation on \f$ \varia \f$:
 *  \f[
 *  \dfrac{\partial \varia}{\partial t} + \divs \left( \varia \vect{g} \right)
 *      - \divs \left( \vect{V} \right) \varia = 0
 *  \f]
 *  where \f$ \vect{g} \f$ is the gravity field
 *
 *  The boundary conditions on \f$ \varia \f$ read:
 *  \f[
 *   \varia = z \textrm{ on walls}
 *  \f]
 *  \f[
 *   \dfrac{\partial \varia}{\partial n} = 0 \textrm{ elsewhere}
 *  \f]
 *
 *  Remarks:
 *  - a steady state is looked for.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_z_ground_compute(void)
{

  /* Initialization
   *===============*/

  if (!_atmo_option.compute_z_ground)
    return;

  const cs_domain_t *domain = cs_glob_domain;
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_real_3_t *restrict i_face_normal =
     (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal =
     (const cs_real_3_t *restrict)mq->b_face_normal;
  const cs_real_3_t *restrict b_face_cog =
     (const cs_real_3_t *restrict)mq->b_face_cog;

  const int *bc_type = cs_glob_bc_type;

  /* Pointer to z_ground field */
  cs_field_t *f = cs_field_by_name_try("z_ground");

  cs_real_t *restrict i_massflux =
    cs_field_by_id(
        cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id")))->val;
  cs_real_t *restrict b_massflux =
    cs_field_by_id(
        cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id")))->val;

  cs_var_cal_opt_t vcopt;
  cs_field_get_key_struct(f, cs_field_key_id("var_cal_opt"), &vcopt);

  cs_real_3_t normal;
  /* Normal direction is given by the gravity */
  cs_math_3_normalise(
      (const cs_real_t *)(cs_glob_physical_constants->gravity),
      normal);

  for (int i = 0; i < 3; i++)
    normal[i] *= -1.;

  /* Compute the mass flux due to V = - g / ||g||
   *=============================================*/

  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
    i_massflux[face_id] = cs_math_3_dot_product(normal, i_face_normal[face_id]);

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++)
    b_massflux[face_id] = cs_math_3_dot_product(normal, b_face_normal[face_id]);


  /* Boundary conditions
   *====================*/

  cs_real_t norm = 0.;
  cs_real_t ground_surf = 0.;

  /* Dirichlet at walls, homogeneous Neumann elsewhere */

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {
    /* Dirichlet BCs */
    if ((bc_type[face_id] == CS_SMOOTHWALL || bc_type[face_id] == CS_ROUGHWALL)
        && b_massflux[face_id] <= 0.) {

      vcopt.ndircl = 1;
      cs_real_t hint = 1. / mq->b_dist[face_id];
      cs_real_t pimp = cs_math_3_dot_product(b_face_cog[face_id], normal);

      cs_boundary_conditions_set_dirichlet_scalar(&(f->bc_coeffs->a[face_id]),
                                                  &(f->bc_coeffs->af[face_id]),
                                                  &(f->bc_coeffs->b[face_id]),
                                                  &(f->bc_coeffs->bf[face_id]),
                                                  pimp,
                                                  hint,
                                                  cs_math_infinite_r);
      norm += cs_math_pow2(f->bc_coeffs->a[face_id]) * mq->b_face_surf[face_id];
      ground_surf += mq->b_face_surf[face_id];
    }
    /* Neumann Boundary Conditions */
    else {

      cs_real_t hint = 1. / mq->b_dist[face_id];
      cs_real_t qimp = 0.;

      cs_boundary_conditions_set_neumann_scalar(&(f->bc_coeffs->a[face_id]),
                                                &(f->bc_coeffs->af[face_id]),
                                                &(f->bc_coeffs->b[face_id]),
                                                &(f->bc_coeffs->bf[face_id]),
                                                qimp,
                                                hint);

    }
  }

  cs_parall_max(1, CS_INT_TYPE, &(vcopt.ndircl));

  /* Matrix
   *=======*/

  cs_real_t *rovsdt, *dpvar;
  BFT_MALLOC(rovsdt, m->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(dpvar, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    rovsdt[cell_id] = 0.;

  /* Right hand side
   *================*/

  cs_real_t *rhs;
  BFT_MALLOC(rhs, m->n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id< m->n_cells_with_ghosts; cell_id++)
    rhs[cell_id] = 0.;

  /* Norm
   *======*/

  cs_parall_sum(1, CS_REAL_TYPE, &norm);
  cs_parall_sum(1, CS_REAL_TYPE, &ground_surf);

  if (ground_surf > 0.)
    norm = sqrt(norm / ground_surf) * mq->tot_vol;
  else {
    bft_printf("No ground BC or no gravity: no computation of ground elevation.\n");
    return;
  }

  /* Solving
   *=========*/

  /* In case of a theta-scheme, set theta = 1;
     no relaxation in steady case either */

  cs_equation_iterative_solve_scalar(0,   /* idtvar: no steady state algo */
                                     -1,  /* no over loops */
                                     f->id,
                                     f->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     norm,
                                     &vcopt,
                                     f->val_pre,
                                     f->val,
                                     f->bc_coeffs->a,
                                     f->bc_coeffs->b,
                                     f->bc_coeffs->af,
                                     f->bc_coeffs->bf,
                                     i_massflux,
                                     b_massflux,
                                     i_massflux, /* viscosity, not used */
                                     b_massflux, /* viscosity, not used */
                                     i_massflux, /* viscosity, not used */
                                     b_massflux, /* viscosity, not used */
                                     NULL,
                                     NULL,
                                     NULL,
                                     0, /* icvflb (upwind) */
                                     NULL,
                                     rovsdt,
                                     rhs,
                                     f->val,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Free memory */
  BFT_FREE(dpvar);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);

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
  if (file_name == NULL) {
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
 * \brief This function declares additional transported variables for
 *        atmospheric module for the chemistry defined from SPACK.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_declare_chem_from_spack(void)
{
  assert(_atmo_chem.model == 4);

  if (_atmo_chem.spack_file_name == NULL)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo chemistry from SPACK file: requires a SPACK file."));

  char line[512] = "";

  /* Open file */
  bft_printf("SPACK file for atmo chemistry:\n    %s \n",
      _atmo_chem.spack_file_name);

  FILE* file = fopen(_atmo_chem.spack_file_name, "rt");
  if (file == NULL)
    bft_error(__FILE__,__LINE__, 0,
              _("Atmo chemistry from SPACK file: Could not open file."));

  bool loop = true;

  /* Read "[species]" */
  for (int i = 1; loop; i++ ) {
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not skip header."));

    if (strcmp("[species]", line) == 0)
      loop = false;
  }

  loop = true;
  /* Read SPACK: first loop count the number of species */
  for (int i = 1; loop; i++ ) {
    /* Read species */
    if (fscanf(file, "%s\n", line) != 1)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmo chemistry from SPACK file: Could not read line %d."),
                i);

    /* When reach [molecular_waight]: break */
    if (strcmp("[molecular_weight]", line) == 0)
      loop = false;
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
    if (fscanf(file, "%s %lf\n", line, &(_atmo_chem.molar_mass[i])) != 2)
      bft_error(__FILE__,__LINE__, 0,
                _("Atmo chemistry from SPACK file: could not read species name and molar mass, line %d."),
                i);

    /* The order is already ok */
    _atmo_chem.chempoint[i] = i+1;//FIXME ?

    /* Build name of the field:
     * species_name in lower case */
    strcpy(name, "species_");
    _strtolower(label, line);
    strcat(name, label);

    /* Field of dimension 1 */
    _atmo_chem.species_to_field_id[i] = cs_variable_field_create(name, line, CS_MESH_LOCATION_CELLS, 1);

    /* Scalar field, store in isca_chem/species_to_scalar_id (FORTRAN/C) array */
    _atmo_chem.species_to_scalar_id[i] = cs_add_model_field_indexes(_atmo_chem.species_to_field_id[i]);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
