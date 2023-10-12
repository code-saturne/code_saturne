#ifndef __CS_FIELD_POINTER_H__
#define __CS_FIELD_POINTER_H__

/*============================================================================
 * Field pointers and ids for standard and model fields
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Macro used for scoping of field pointer enums */

#define CS_ENUMF_(e) CS_FIELD_POINTER_ ## e

/* Macro used to return a field pointer by its enumerated value */

#define CS_F_(e) cs_glob_field_pointers[CS_FIELD_POINTER_ ## e].f

#define CS_FI_(e, i) cs_glob_field_pointers[CS_FIELD_POINTER_ ## e].p[i]

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Suffixes of enumerated field pointer ids, arguments of the macro CS_F_(). */

typedef enum {

  /* Base variables and properties */

  CS_ENUMF_(dt),           /*!< local time step */

  CS_ENUMF_(p),            /*!< pressure */
  CS_ENUMF_(vel),          /*!< velocity */

  CS_ENUMF_(k),            /*!< turbulent kinetic energy \f$ k \f$ */
  CS_ENUMF_(eps),          /*!< turbulent dissipation \f$ \varepsilon \f$ */

  CS_ENUMF_(rij),          /*!< Reynolds stress tensor \f$ R_{ij} \f$ */

  CS_ENUMF_(phi),          /*!< \f$ \phi \f$ for \f$ \phi-f_b \f$ model */
  CS_ENUMF_(f_bar),        /*!< \f$ f_b \f$ for \f$ \phi-f_b \f$ model */
  CS_ENUMF_(alp_bl),        /*!< \f$ \alpha \f$ for \f$ Bl-v^2-k \f$
                                 or EBRSM model */

  CS_ENUMF_(omg),          /*!< \f$ \omega \f$ for \f$ k-\omega \f$ SST model */
  CS_ENUMF_(nusa),         /*!< \f$ \widetilde{\nu}_T \f$ for Spalart-
                                Allmaras */

  CS_ENUMF_(hybrid_blend), /*!< Blending factor for DDES*/

  CS_ENUMF_(mesh_u),       /*!< mesh velocity */

  CS_ENUMF_(void_f),       /*!< void fraction */
  CS_ENUMF_(vol_f),        /*!< volume fraction */

  CS_ENUMF_(h),            /*!< enthalpy */
  CS_ENUMF_(t),            /*!< temperature*/
  CS_ENUMF_(t_b),          /*!< temperature (at boundary faces)*/
  CS_ENUMF_(e_tot),        /*!< total energy */
  CS_ENUMF_(h_tot),        /*!< total enthalpy */

  CS_ENUMF_(rho),          /*!< density (at cells) */
  CS_ENUMF_(rho_b),        /*!< density (at boundary faces) */

  CS_ENUMF_(cp),           /*!< isobaric specific heat */
  CS_ENUMF_(cv),           /*!< isochoric specific heat */

  CS_ENUMF_(mu),           /*!< molecular viscosity */
  CS_ENUMF_(mu_t),         /*!< turbulent dynamic viscosity */

  CS_ENUMF_(lambda),       /*!< Thermal conductivity */
  CS_ENUMF_(th_diff),      /*!< Thermal diffusivity */

  CS_ENUMF_(poro),         /*!< porosity */
  CS_ENUMF_(if_poro),      /*!< internal faces porosity */
  CS_ENUMF_(t_poro),       /*!< tensorial porosity */

  /* Specific physics variables and properties */

  CS_ENUMF_(t_kelvin),     /*!< temperature, in Kelvin */

  CS_ENUMF_(vism),         /*!< mesh viscosity */

  CS_ENUMF_(volume_f),     /*!< homogeneous model volume fraction */
  CS_ENUMF_(mass_f),       /*!< homogeneous model mass fraction */
  CS_ENUMF_(energy_f),     /*!< homogeneous model energy fraction */

  CS_ENUMF_(pot_t),        /*!< potential temperature */
  CS_ENUMF_(ntdrp),        /*!< total number of droplets */
  CS_ENUMF_(chemistry),    /*!< chemistry species (indexed) */

  CS_ENUMF_(fm),           /*!< mixture fraction */
  CS_ENUMF_(fp2m),         /*!< mixture fraction variance */

  CS_ENUMF_(fsm),          /*!< soot mass fraction */
  CS_ENUMF_(npm),          /*!< soot precursor number */
  CS_ENUMF_(ygfm),         /*!< fresh gas fraction */

  CS_ENUMF_(yfm),          /*!< mass fraction */
  CS_ENUMF_(yfp2m),        /*!< mass fraction variance */
  CS_ENUMF_(coyfp),        /*!< mass fraction covariance */

  CS_ENUMF_(np),           /*!< particles per kg for coal class */
  CS_ENUMF_(xch),          /*!< reactive coal mass fraction for coal class */
  CS_ENUMF_(xck),          /*!< coke mass fraction for coal class */
  CS_ENUMF_(xwt),          /*!< water mass fraction for coal class */
  CS_ENUMF_(h2),           /*!< mass enthalpy for coal class (permeatic case) */
  CS_ENUMF_(f1m),          /*!< mean value light volatiles for coal class */
  CS_ENUMF_(f2m),          /*!< mean value heavy volatiles for coal class */
  CS_ENUMF_(f4m),          /*!< oxydant 2 mass fraction */
  CS_ENUMF_(f5m),          /*!< oxydant 3 mass fraction */
  CS_ENUMF_(f6m),          /*!< water from coal drying mass fraction */
  CS_ENUMF_(f7m),          /*!< carbon from coal oxidyzed by O2 mass fraction */
  CS_ENUMF_(f8m),          /*!< carbon from coal gasified by CO2 mass fraction */
  CS_ENUMF_(f9m),          /*!< carbon from coal gasified by H2O mass fraction */
  CS_ENUMF_(fvp2m),        /*!< f1f2 variance */
  CS_ENUMF_(yco2),         /*!< CO2 fraction */
  CS_ENUMF_(yhcn),         /*!< HCN fraction */
  CS_ENUMF_(yno),          /*!< NO fraction */
  CS_ENUMF_(ynh3),         /*!< NH3 enthalpy */
  CS_ENUMF_(hox),          /*!< Ox enthalpy */

  CS_ENUMF_(potr),         /*!< Electric potential, real part */
  CS_ENUMF_(poti),         /*!< Electric potential, imaginary part */
  CS_ENUMF_(potva),        /*!< Vector potential */
  CS_ENUMF_(ycoel),        /*!< Constituent mass fraction */
  CS_ENUMF_(joulp),        /*!< Joule power */
  CS_ENUMF_(radsc),        /*!< radiation source */
  CS_ENUMF_(elech),        /*!< electric charge */
  CS_ENUMF_(curre),        /*!< current real */
  CS_ENUMF_(curim),        /*!< current imaginary */
  CS_ENUMF_(laplf),        /*!< laplace forces */
  CS_ENUMF_(magfl),        /*!< magnetic field */
  CS_ENUMF_(elefl),        /*!< electric field */

  CS_ENUMF_(rad_energy),   /*!< Radiative energy */
  CS_ENUMF_(rad_q),        /*!< Radiative flux */
  CS_ENUMF_(radiance),     /*!< (spectral) radiance field */

  CS_ENUMF_(rad_est),      /*!< Radiative flux explicit source term */
  CS_ENUMF_(rad_ist),      /*!< Radiative flux implicit source term */
  CS_ENUMF_(rad_abs),      /*!< Radiative absorption */
  CS_ENUMF_(rad_emi),      /*!< Radiative emission */
  CS_ENUMF_(rad_cak),      /*!< Radiative absorption coefficient */

  CS_ENUMF_(qinci),        /*!< Radiative incident radiative flux density */
  CS_ENUMF_(qinsp),        /*!< Spectral radiative incident flux */
  CS_ENUMF_(xlam),         /*!< Wall thermal conductivity */
  CS_ENUMF_(epa),          /*!< Wall thickness */
  CS_ENUMF_(emissivity),   /*!< Wall emissivity */
  CS_ENUMF_(fnet),         /*!< Boundary radiative flux */
  CS_ENUMF_(fconv),        /*!< Boundary radiative convective flux */
  CS_ENUMF_(hconv),        /*!< radiative exchange coefficient */
  CS_ENUMF_(fup),          /*!< Spectral upward radiative flux */
  CS_ENUMF_(fdown),        /*!< Spectral downward radiative flux */
  CS_ENUMF_(rad_ck_up),    /*!< Spectral upward Radiative absorption coefficient */
  CS_ENUMF_(rad_ck_down),  /*!< Spectral downward Radiative absorption
                             coefficient */

  CS_ENUMF_(mol_mass),     /*!< gas mix molar max */

  CS_ENUMF_(head),         /*!< hydraulic head */

  /* Cooling tower fields */
  CS_ENUMF_(humid),          /*!< Humidity */
  CS_ENUMF_(ym_w),           /*!< Mass fraction of dry air in humid air */
  CS_ENUMF_(t_l),            /*!< Injected liquid water temperature */
  CS_ENUMF_(h_l),            /*!< Injected liquid water enthalpy */
  CS_ENUMF_(y_l_pack),       /*!< Mass of liquid per unit volume of cell */
  CS_ENUMF_(thermal_diff_h), /*!< Humid air thermal diffusivity  */
  CS_ENUMF_(thermal_diff_l), /*!< Injected liquid water thermal diffusivity */
  CS_ENUMF_(pack_zone_id),   /*!< Id of the packing zone */

  /* NCFD fields */
  CS_ENUMF_(yf_ncond),       /*!< non-condensable mass fraction */
  CS_ENUMF_(qp),             /*!< Turbulent Kinetic Energy q2 */
  CS_ENUMF_(qfp),            /*!< Covariance q12 */
  CS_ENUMF_(qfpxx),          /*!< XX component of qfp */
  CS_ENUMF_(qfpxy),          /*!< XY component of qfp */
  CS_ENUMF_(qfpxz),          /*!< XZ component of qfp */
  CS_ENUMF_(qfpyx),          /*!< YX component of qfp */
  CS_ENUMF_(qfpyy),          /*!< YY component of qfp */
  CS_ENUMF_(qfpyz),          /*!< YZ component of qfp */
  CS_ENUMF_(qfpzx),          /*!< ZX component of qfp */
  CS_ENUMF_(qfpzy),          /*!< ZY component of qfp */
  CS_ENUMF_(qfpzz),          /*!< ZZ component of qfp */
  CS_ENUMF_(gamma),          /*!< Interfacial mass transfer */
  CS_ENUMF_(ia),             /*!< Interfacial area */
  CS_ENUMF_(x2),             /*!< x2 for droplets */
  CS_ENUMF_(d32),            /*!< Sauter diameter */
  CS_ENUMF_(drag),           /*!< Phases drag */
  CS_ENUMF_(ad_mass),        /*!< Added mass */
  CS_ENUMF_(wlubr),          /*!< Wall lubrication */
  CS_ENUMF_(th_diff_t),      /*!< Turbulent thermal diffusivity */
  CS_ENUMF_(drho_dp),        /*!< drho over dp */
  CS_ENUMF_(drho_dh),        /*!< drho over dh */
  CS_ENUMF_(tau12_t),        /*!< turbulent tau12 */
  CS_ENUMF_(lift),           /*!< Lift coefficient */
  CS_ENUMF_(disp_t),         /*!< Turbulent dispersion */
  CS_ENUMF_(surf_tens),      /*!< Surface tension */
  CS_ENUMF_(sl_corr),        /*!< Free surface correction weight for GLIM */
  CS_ENUMF_(fi),             /*!< field_a vol_f weight for GLIM */
  CS_ENUMF_(fj),             /*!< field_b vol_f weight for GLIM */
  CS_ENUMF_(drift_vel),      /*!< Particles drift velocity */
  CS_ENUMF_(yplus),          /*!< Wall distance: y+ */
  CS_ENUMF_(vel_mean),       /*!< Mean velocity (for dispersed phases) */
  CS_ENUMF_(vel_rel),        /*!< Relative velocity (for dispersed phases) */
  CS_ENUMF_(dt_dp),          /*!< dtemp/dpress derivative */
  CS_ENUMF_(kindiff),        /*!< Particles kinetic diffusivity */
  CS_ENUMF_(coldiff),        /*!< Particles collisional diffusivity */
  CS_ENUMF_(elast),          /*!< Particles restitution coefficient */
  CS_ENUMF_(c_alpha),        /*!< Weighting coefficient for GEMMA model */

  /* Added variables (scalars) */
  CS_ENUMF_(add_var),        /*!< User added variables */

  /* User-defined arrays */
  CS_ENUMF_(user),

  /* End of attributes */

  CS_FIELD_N_POINTERS

} cs_field_pointer_id_t;

/*! Field pointer array type */

struct cs_field_pointer_array_t {
  cs_field_t   *f;   /*!< pointer to single (first) field */
  cs_field_t  **p;   /*!< array of field pointers */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointers */

extern struct cs_field_pointer_array_t  *cs_glob_field_pointers;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ensure field pointer array is initialized.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_ensure_init(void);

/*----------------------------------------------------------------------------
 * Free all field pointer data.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_destroy_all(void);

/*----------------------------------------------------------------------------
 * Map a simple field to an enumerated pointer.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e).
 *
 * parameters:
 *   e <--  field enumerator value
 *   f <--  pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map(cs_field_pointer_id_t   e,
                     cs_field_t             *f);

/*----------------------------------------------------------------------------
 * Map a field to an (enumerated pointer, index) couple.
 *
 * This sort of mapping may be used for sets of fields whose size
 * is not known in advance.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e, i).
 *
 * parameters:
 *   e     <-- field enumerator value
 *   index <-- field enumerator index
 *   f     <-- pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_indexed(cs_field_pointer_id_t   e,
                             int                     index,
                             cs_field_t             *f);

/*----------------------------------------------------------------------------
 * Map base fields to enumerated pointers.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_base(void);

/*----------------------------------------------------------------------------
 * Map some boundary fields to enumerated pointers.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_boundary(void);

/*----------------------------------------------------------------------------
 * Map base fields to enumerated pointers for atmospheric models
 *
 * parameters:
 *   n_chem_species <-- number of chemical species
 *   species_f_if   <-- field id for each chemical species
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_atmospheric(int        n_chem_species,
                                 const int  species_f_id[]);

/*----------------------------------------------------------------------------
 * Map base fields to enumerated pointers for atmospheric models
 *
 * parameters:
 *   n_coals   <-- number of coals
 *   n_classes <-- number of coal classes
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_coal_combustion(int  n_coals,
                                     int  n_classes);

/*----------------------------------------------------------------------------*
 * Map base fields to enumerated pointers for compressible model
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_compressible(void);

/*----------------------------------------------------------------------------*
 * Map base fields to enumerated pointers for gas mix model
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_gas_mix(void);

/*----------------------------------------------------------------------------
 * Map base fields to enumerated pointers for gas combustion.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_gas_combustion(void);

/*----------------------------------------------------------------------------*/
/*
 * Map base fields to enumerated pointers for groundwater flows
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_groundwater(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_POINTER_H__ */
