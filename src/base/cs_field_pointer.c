/*============================================================================
 * Field pointers and ids for standard and model fields
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

#include "cs_field.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field_pointer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field_pointer.c
        Field pointers and ids for standard and model fields.
*/

/*! \fn CS_ENUMF_(e)
 * \brief Macro used for scoping of field pointer enums.
 *
 * This macro replaces CS_ENUMF_ by CS_FIELD_POINTER_ and allows
 * to rebuild a full enumerated field pointer id.
 *
 * \param [in] e suffix of enumerated field pointer id.
 */

/*! \fn CS_F_(e)
 * \brief Macro used to return a field pointer by its enumerated value.
 *
 * This macro replaces CS_F_ by an access to the global array of field pointers
 * \ref cs_glob_field_pointers using a rebuilt enumerated field pointer id.
 *
 * \param [in] e suffix of enumerated field pointer id.
 */

/*! \fn CS_FI_(e, i)
 * \brief Macro used to return a field pointer by its enumerated value.
 *
 * This macro replaces CS_FI_ by an access to the global array of field pointers
 * \ref cs_glob_field_pointers using a rebuilt enumerated field pointer id and
 * its field sublist index.
 *
 * \param [in] e suffix of enumerated field pointer id.
 * \param [in] i field enumerator value.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of pointers (initially fixed, but extensible in case
   user fields should be added after the model fields) */

static cs_field_pointer_id_t _n_pointers = 0;
static union cs_field_pointer_val_t  *_field_pointer = NULL;

/* Handling of sublists */

static bool  *_is_sublist = NULL;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointers */

union cs_field_pointer_val_t  *cs_glob_field_pointers = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize field pointers
 *----------------------------------------------------------------------------*/

static void
_init_pointers(void)
{
  assert(_field_pointer == NULL);

  _n_pointers = CS_FIELD_N_POINTERS;
  BFT_MALLOC(_field_pointer, _n_pointers, union cs_field_pointer_val_t);
  BFT_MALLOC(_is_sublist, _n_pointers, bool);

  for (cs_field_pointer_id_t i = 0; i < _n_pointers; i++) {
    _field_pointer[i].f = NULL;
    _is_sublist[i] = false;
  }

  cs_glob_field_pointers = _field_pointer;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all field pointer data.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_destroy_all(void)
{
  for (cs_field_pointer_id_t i = 0; i < _n_pointers; i++) {
    if (_is_sublist[i])
      BFT_FREE(_field_pointer[i].a);
  }
  BFT_FREE(_field_pointer);
  BFT_FREE(_is_sublist);

  cs_glob_field_pointers = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map a simple field to an enumerated pointer.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e).
 *
 * \param[in]  e   field enumerator value
 * \param[in]  f   pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map(cs_field_pointer_id_t   e,
                     cs_field_t             *f)
{
  if (_field_pointer == NULL)
    _init_pointers();

  assert(e < _n_pointers);

  union cs_field_pointer_val_t v;
  v.f = f;
  _field_pointer[e] = v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map a field to an (enumerated pointer, index) couple.
 *
 * This sort of mapping may be used for sets of fields whose size
 * is not known in advance.
 *
 * The associated field pointer may then be retreived using \ref CS_FI_(e, i).
 *
 * \param[in]  e      field enumerator value
 * \param[in]  index  field enumerator index
 * \param[in]  f      pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_indexed(cs_field_pointer_id_t   e,
                             int                     index,
                             cs_field_t             *f)
{
  assert(index >= 0);

  if (_field_pointer == NULL)
    _init_pointers();

  struct cs_field_pointer_array_t *a;
  union cs_field_pointer_val_t v;

  int i;
  int _sub_size = index + 1;
  int _sub_size_prev = 0;

  /* Check for previous size */

  assert(e < _n_pointers);

  v = _field_pointer[e];

  /* Check also that we did not already use in incompatible mapping */

  if (v.f != NULL && ! _is_sublist[e]) {
    cs_field_t *_f = v.f;
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field enum %d is already mapped as non-indexed\n"
                "to field id %d (%s), so it cannot be mapped as indexed."),
              __func__, (int)e, _f->id, _f->name);
  }

  a = v.a;

  if (a != NULL)
    _sub_size_prev = a->n;

  if (_sub_size_prev < _sub_size) {
    /* BFT_MALLOC does not directly handle C flexible array members,
       (which we use here to minimize type width), so use bytes */
    void *p = a;
    size_t _s =   sizeof(struct cs_field_pointer_array_t)
                + sizeof(cs_field_t *) * _sub_size;
    BFT_REALLOC(p, _s, unsigned char);
    a = p;
    a->n = _sub_size;
    for (i = _sub_size_prev; i < index; i++)
      a->p[i] = NULL;
    v.a = a;
  }

  _is_sublist[e] = true;

  v.a->p[index] = f;
  _field_pointer[e] = v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_base(void)
{
  cs_field_pointer_map(CS_ENUMF_(dt),
                       cs_field_by_name_try("dt"));

  cs_field_pointer_map(CS_ENUMF_(p),
                       cs_field_by_name_try("pressure"));
  cs_field_pointer_map(CS_ENUMF_(vel),
                       cs_field_by_name_try("velocity"));

  cs_field_pointer_map(CS_ENUMF_(k),
                       cs_field_by_name_try("k"));
  cs_field_pointer_map(CS_ENUMF_(eps),
                       cs_field_by_name_try("epsilon"));

  cs_field_pointer_map(CS_ENUMF_(r11), cs_field_by_name_try("r11"));
  cs_field_pointer_map(CS_ENUMF_(r22), cs_field_by_name_try("r22"));
  cs_field_pointer_map(CS_ENUMF_(r33), cs_field_by_name_try("r33"));
  cs_field_pointer_map(CS_ENUMF_(r12), cs_field_by_name_try("r12"));
  cs_field_pointer_map(CS_ENUMF_(r23), cs_field_by_name_try("r23"));
  cs_field_pointer_map(CS_ENUMF_(r13), cs_field_by_name_try("r13"));
  cs_field_pointer_map(CS_ENUMF_(rij), cs_field_by_name_try("rij"));

  cs_field_pointer_map(CS_ENUMF_(phi), cs_field_by_name_try("phi"));
  cs_field_pointer_map(CS_ENUMF_(f_bar), cs_field_by_name_try("f_bar"));
  cs_field_pointer_map(CS_ENUMF_(alp_bl), cs_field_by_name_try("alpha"));

  cs_field_pointer_map(CS_ENUMF_(omg), cs_field_by_name_try("omega"));
  cs_field_pointer_map(CS_ENUMF_(nusa), cs_field_by_name_try("nu_tilda"));

  cs_field_pointer_map(CS_ENUMF_(hybrid_blend), cs_field_by_name_try("hybrid_blend"));

  cs_field_pointer_map(CS_ENUMF_(mesh_u),
                       cs_field_by_name_try("mesh_velocity"));

  cs_field_pointer_map(CS_ENUMF_(void_f),
                       cs_field_by_name_try("void_fraction"));

  cs_field_pointer_map(CS_ENUMF_(h),
                       cs_field_by_name_try("enthalpy"));
  cs_field_pointer_map(CS_ENUMF_(t),
                       cs_field_by_name_try("temperature"));

  cs_field_pointer_map(CS_ENUMF_(rho),
                       cs_field_by_name_try("density"));

  cs_field_pointer_map(CS_ENUMF_(cp),
                       cs_field_by_name_try("specific_heat"));

  cs_field_pointer_map(CS_ENUMF_(mu),
                       cs_field_by_name_try("molecular_viscosity"));
  cs_field_pointer_map(CS_ENUMF_(mu_t),
                       cs_field_by_name_try("turbulent_viscosity"));

  cs_field_pointer_map(CS_ENUMF_(vism),
                       cs_field_by_name_try("mesh_viscosity"));

  cs_field_pointer_map(CS_ENUMF_(poro),
                       cs_field_by_name_try("porosity"));
  cs_field_pointer_map(CS_ENUMF_(t_poro),
                       cs_field_by_name_try("tensorial_porosity"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map some boundary fields to enumerated pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_boundary(void)
{
  cs_field_pointer_map(CS_ENUMF_(t_b),
                       cs_field_by_name_try("boundary_temperature"));

  cs_field_pointer_map(CS_ENUMF_(rho_b),
                       cs_field_by_name_try("boundary_density"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for atmospheric models
 *
 * \param[in]  n_chem_species  number of chemical species
 * \param[in]  species_f_id    field id for each chemical species
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_atmospheric(int        n_chem_species,
                                 const int  species_f_id[])
{
  cs_field_pointer_map(CS_ENUMF_(pot_t),
                       cs_field_by_name_try("temperature"));

  cs_field_pointer_map(CS_ENUMF_(totwt),
                       cs_field_by_name_try("total_water"));
  cs_field_pointer_map(CS_ENUMF_(ntdrp),
                       cs_field_by_name_try("number_of_droplets"));

  for (int i = 0; i < n_chem_species; i++)
    cs_field_pointer_map_indexed(CS_ENUMF_(chemistry),
                                 i,
                                 cs_field_by_id(species_f_id[i]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for coal combustion
 *
 * \param[in]  n_coals    number of coals
 * \param[in]  n_classes  number of coal classes
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_coal_combustion(int  n_coals,
                                     int  n_classes)
{
  char s[64];

  cs_field_pointer_map(CS_ENUMF_(h),
                       cs_field_by_name_try("enthalpy"));

  for (int i = 0; i < n_classes; i++) {
    snprintf(s, 63, "n_p_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(np), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_classes; i++) {
    snprintf(s, 63, "x_p_coal_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(xch), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_classes; i++) {
    snprintf(s, 63, "x_p_char_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(xck), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_classes; i++) {
    snprintf(s, 63, "x_p_wt_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(xwt), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_classes; i++) {
    snprintf(s, 63, "x_p_h_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(h2), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_coals; i++) {
    snprintf(s, 63, "fr_mv1_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(f1m), i, cs_field_by_name_try(s));
  }

  for (int i = 0; i < n_coals; i++) {
    snprintf(s, 63, "fr_mv2_%02d", i+1); s[63] = '\0';
    cs_field_pointer_map_indexed(CS_ENUMF_(f2m), i, cs_field_by_name_try(s));
  }

  cs_field_pointer_map(CS_ENUMF_(f4m), cs_field_by_name_try("fr_oxyd2"));
  cs_field_pointer_map(CS_ENUMF_(f5m), cs_field_by_name_try("fr_oxyd3"));
  cs_field_pointer_map(CS_ENUMF_(f6m), cs_field_by_name_try("fr_h2o"));
  cs_field_pointer_map(CS_ENUMF_(f7m), cs_field_by_name_try("fr_het_o2"));
  cs_field_pointer_map(CS_ENUMF_(f8m), cs_field_by_name_try("fr_het_co2"));
  cs_field_pointer_map(CS_ENUMF_(f9m), cs_field_by_name_try("fr_het_h2o"));

  cs_field_pointer_map(CS_ENUMF_(fvp2m), cs_field_by_name_try("f1f2_variance"));

  cs_field_pointer_map(CS_ENUMF_(yco2), cs_field_by_name_try("x_c_co2"));
  cs_field_pointer_map(CS_ENUMF_(yhcn), cs_field_by_name_try("x_c_hcn"));
  cs_field_pointer_map(CS_ENUMF_(yno), cs_field_by_name_try("x_c_no"));
  cs_field_pointer_map(CS_ENUMF_(ynh3), cs_field_by_name_try("x_c_nh3"));

  cs_field_pointer_map(CS_ENUMF_(hox), cs_field_by_name_try("x_c_h_ox"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for compressible model
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_compressible(void)
{
  cs_field_pointer_map(CS_ENUMF_(e_tot),
                       cs_field_by_name_try("total_energy"));

  cs_field_pointer_map(CS_ENUMF_(t_kelvin),
                       cs_field_by_name_try("temperature"));

  /* Also map to main temperature pointer */

  cs_field_pointer_map(CS_ENUMF_(t),
                       cs_field_by_name_try("temperature"));

  /* map volume specific heat if it is non constant */

  cs_field_pointer_map(CS_ENUMF_(cv),
                       cs_field_by_name_try("specific_heat_const_vol"));

  /* map fractions for homogeneous two phase model */

  cs_field_pointer_map(CS_ENUMF_(volume_f),
                       cs_field_by_name_try("volume_fraction"));

  cs_field_pointer_map(CS_ENUMF_(mass_f),
                       cs_field_by_name_try("mass_fraction"));

  cs_field_pointer_map(CS_ENUMF_(energy_f),
                       cs_field_by_name_try("energy_fraction"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for gas mix model
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_gas_mix(void)
{
  cs_field_pointer_map(CS_ENUMF_(mol_mass),
                       cs_field_by_name_try("mix_mol_mas"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for gas combustion.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_gas_combustion(void)
{
  cs_field_pointer_map(CS_ENUMF_(h),
                       cs_field_by_name_try("enthalpy"));

  cs_field_pointer_map(CS_ENUMF_(fm),
                       cs_field_by_name_try("mixture_fraction"));
  cs_field_pointer_map(CS_ENUMF_(fp2m),
                       cs_field_by_name_try("mixture_fraction_variance"));

  cs_field_pointer_map(CS_ENUMF_(fsm),
                       cs_field_by_name_try("soot_mass_fraction"));

  cs_field_pointer_map(CS_ENUMF_(npm),
                       cs_field_by_name_try("soot_precursor_number"));

  cs_field_pointer_map(CS_ENUMF_(ygfm),
                       cs_field_by_name_try("fresh_gas_fraction"));

  cs_field_pointer_map(CS_ENUMF_(yfm),
                       cs_field_by_name_try("mass_fraction"));
  cs_field_pointer_map(CS_ENUMF_(yfp2m),
                       cs_field_by_name_try("mass_fraction_variance"));
  cs_field_pointer_map(CS_ENUMF_(coyfp),
                       cs_field_by_name_try("mass_fraction_covariance"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers for groundwater flows
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_groundwater(void)
{
  cs_field_pointer_map(CS_ENUMF_(head),
                       cs_field_by_name_try("hydraulic_head"));
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
