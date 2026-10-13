/*============================================================================
 * Physical units management and catalogue
 *============================================================================*/

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

#include "bft/bft_mem.h"

#include "base/cs_defs.h"
#include "base/cs_map.h"
#include "base/cs_thermal_model.h"
#include "base/cs_physical_units.h"

/*----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

typedef struct {
  const char *symbol;
  const char *description;
} cs_unit_def_t;

/*----------------------------------------------------------------------------
 * Static variables
 *----------------------------------------------------------------------------*/

static const cs_unit_def_t _standard_units[] = {
  {"-",        "dimensionless"},
  {"%",        "percentage"},
  {"rad",      "radian"},
  {"deg",      "degree"},
  {"s",        "second"},
  {"s-1",      "per second"},
  {"s-2",      "per second squared"},
  {"Hz",       "Hertz"},
  {"rad/s",    "radian per second"},
  {"m",        "meter"},
  {"m2",       "square meter"},
  {"m3",       "cubic meter"},
  {"m/s",      "meter per second"},
  {"m/s2",     "meter per second squared"},
  {"m2/s",     "square meter per second"},
  {"m2/s2",    "square meter per second squared"},
  {"m2/s3",    "square meter per cubic second"},
  {"m2/s4",    "square meter per second to the fourth"},
  {"kg",       "kilogram"},
  {"g",        "gram"},
  {"kg/m3",    "kilogram per cubic meter"},
  {"kg/s",     "kilogram per second"},
  {"m3/s",     "cubic meter per second"},
  {"Pa",       "Pascal"},
  {"bar",      "bar"},
  {"Pa.s",     "Pascal-second"},
  {"Pa/m",     "Pascal per meter"},
  {"N",        "Newton"},
  {"N.m",      "Newton-meter"},
  {"J",        "Joule"},
  {"J/m3",     "Joule per cubic meter"},
  {"J/kg",     "Joule per kilogram"},
  {"J/(kg.K)", "Joule per kilogram-Kelvin"},
  {"1/K",      "per Kelvin"},
  {"W",        "Watt"},
  {"W/m2",     "Watt per square meter"},
  {"W/m3",     "Watt per cubic meter"},
  {"W/(m.K)",  "Watt per meter-Kelvin"},
  {"K",        "Kelvin"},
  {"degC",     "degree Celsius"},
  {"mol",      "mole"},
  {"mol/m3",   "mole per cubic meter"},
  {"mol/s",    "mole per second"},
  {"kg/mol",   "kilogram per mole"},
  {"J/mol",    "Joule per mole"},
  {"J/(mol.K)","Joule per mole-Kelvin"},
  {"A",        "Ampere"},
  {"V",        "Volt"},
  {"T",        "Tesla"},
  {"S/m",      "Siemens per meter"},
  {"Ohm.m",    "Ohm-meter"}
};

static cs_map_name_to_id_t  *_units_map = nullptr;
static char                 **_units_desc = nullptr;
static size_t                 _units_desc_size = 0;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize physical units map and populate with standard catalogue.
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_units_init(void)
{
  if (_units_map != nullptr)
    return;

  _units_map = cs_map_name_to_id_create();
  size_t n_std = sizeof(_standard_units) / sizeof(_standard_units[0]);
  _units_desc_size = n_std;
  CS_MALLOC(_units_desc, _units_desc_size, char *);
  for (size_t i = 0; i < _units_desc_size; i++)
    _units_desc[i] = nullptr;

  for (size_t i = 0; i < n_std; i++) {
    int id = cs_map_name_to_id(_units_map, _standard_units[i].symbol);
    if (_standard_units[i].description != nullptr) {
      size_t l = strlen(_standard_units[i].description);
      CS_MALLOC(_units_desc[id], l + 1, char);
      strcpy(_units_desc[id], _standard_units[i].description);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or register a unit symbol with an optional description.
 *
 * If the unit symbol already exists, its description is updated if non-null.
 * If symbol is nullptr or empty, nullptr is returned.
 *
 * \param[in] symbol       unit symbol (e.g. "m/s", "Pa", "kg/m3")
 * \param[in] description  optional unit description, or nullptr
 *
 * \return pointer to persistent interned unit symbol string, or nullptr
 */
/*----------------------------------------------------------------------------*/

const char *
cs_physical_unit_define(const char  *symbol,
                        const char  *description)
{
  if (symbol == nullptr || symbol[0] == '\0')
    return nullptr;

  cs_physical_units_init();

  int id = cs_map_name_to_id(_units_map, symbol);

  if (id >= (int)_units_desc_size) {
    size_t new_size = _units_desc_size == 0 ? 16 : _units_desc_size * 2;
    while ((int)new_size <= id)
      new_size *= 2;
    CS_REALLOC(_units_desc, new_size, char *);
    for (size_t i = _units_desc_size; i < new_size; i++)
      _units_desc[i] = nullptr;
    _units_desc_size = new_size;
  }

  if (description != nullptr) {
    CS_FREE(_units_desc[id]);
    size_t l = strlen(description);
    CS_MALLOC(_units_desc[id], l + 1, char);
    strcpy(_units_desc[id], description);
  }

  return cs_map_name_to_id_reverse(_units_map, id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the description associated with a unit symbol.
 *
 * \param[in] symbol  unit symbol
 *
 * \return pointer to unit description, or nullptr if none registered
 */
/*----------------------------------------------------------------------------*/

const char *
cs_physical_unit_get_description(const char  *symbol)
{
  if (symbol == nullptr || symbol[0] == '\0')
    return nullptr;

  cs_physical_units_init();

  int id = cs_map_name_to_id_try(_units_map, symbol);
  if (id > -1 && id < (int)_units_desc_size)
    return _units_desc[id];

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the default physical unit symbol matching a canonical name.
 *
 * \param[in] name  canonical quantity or property name
 *
 * \return pointer to persistent interned unit symbol string, or nullptr
 */
/*----------------------------------------------------------------------------*/

const char *
cs_physical_unit_get_default(const char  *name)
{
  if (name == nullptr || name[0] == '\0')
    return nullptr;

  /* Kinematic fields */
  if (   strcmp(name, "velocity") == 0
      || strcmp(name, "i_velocity") == 0
      || strcmp(name, "b_velocity") == 0
      || strcmp(name, "mesh_velocity") == 0
      || strcmp(name, "drift_velocity") == 0
      || strcmp(name, "inner_face_velocity") == 0
      || strcmp(name, "boundary_face_velocity") == 0
      || strcmp(name, "meteo_velocity") == 0
      || strcmp(name, "vel_mag_mean") == 0
      || strcmp(name, "velocity_mean") == 0)
    return cs_physical_unit_define("m/s", nullptr);

  /* Acceleration */
  else if (   strcmp(name, "acceleration") == 0
           || strcmp(name, "gravity") == 0)
    return cs_physical_unit_define("m/s2", nullptr);

  /* Pressure */
  else if (   strcmp(name, "pressure") == 0
           || strcmp(name, "total_pressure") == 0
           || strcmp(name, "hydrostatic_pressure") == 0
           || strcmp(name, "meteo_pressure") == 0
           || strcmp(name, "algo:meteo_pressure_corrected") == 0
           || strcmp(name, "p0") == 0
           || strcmp(name, "pred0") == 0
           || strcmp(name, "pther") == 0
           || strcmp(name, "pthera") == 0
           || strcmp(name, "pthermax") == 0)
    return cs_physical_unit_define("Pa", nullptr);

  /* Pressure gradient */
  else if (   strcmp(name, "algo:pressure_gradient") == 0
           || strcmp(name, "algo:pressure_increment_gradient") == 0)
    return cs_physical_unit_define("Pa/m", nullptr);

  /* Density */
  else if (   strcmp(name, "density") == 0
           || strcmp(name, "boundary_density") == 0
           || strcmp(name, "meteo_density") == 0
           || strcmp(name, "ro0") == 0
           || strcmp(name, "roref") == 0)
    return cs_physical_unit_define("kg/m3", nullptr);

  /* Viscosity */
  else if (   strcmp(name, "molecular_viscosity") == 0
           || strcmp(name, "turbulent_viscosity") == 0
           || strcmp(name, "dynamic_viscosity") == 0
           || strcmp(name, "viscl0") == 0
           || strcmp(name, "viscv0") == 0)
    return cs_physical_unit_define("Pa.s", nullptr);

  else if (strcmp(name, "kinematic_viscosity") == 0)
    return cs_physical_unit_define("m2/s", nullptr);

  /* Thermal scalar fields */
  else if (   strcmp(name, "temperature") == 0
           || strcmp(name, "meteo_temperature") == 0
           || strcmp(name, "meteo_pot_temperature") == 0
           || strcmp(name, "t0") == 0) {
    if (   cs_glob_thermal_model != nullptr
        && cs_glob_thermal_model->temperature_scale
           == CS_TEMPERATURE_SCALE_CELSIUS)
      return cs_physical_unit_define("degC", nullptr);
    else
      return cs_physical_unit_define("K", nullptr);
  }

  else if (   strcmp(name, "enthalpy") == 0
           || strcmp(name, "total_energy") == 0
           || strcmp(name, "internal_energy") == 0
           || strcmp(name, "eint0") == 0
           || strcmp(name, "l00") == 0
           || strcmp(name, "clatev") == 0)
    return cs_physical_unit_define("J/kg", nullptr);

  else if (   strcmp(name, "specific_heat") == 0
           || strcmp(name, "isobaric_heat_capacity") == 0
           || strcmp(name, "cp0") == 0
           || strcmp(name, "cv0") == 0
           || strcmp(name, "cpv0") == 0
           || strcmp(name, "cvl") == 0
           || strcmp(name, "r_pg_cnst") == 0
           || strcmp(name, "r_v_cnst") == 0)
    return cs_physical_unit_define("J/(kg.K)", nullptr);

  else if (   strcmp(name, "thermal_conductivity") == 0
           || strcmp(name, "lambda0") == 0)
    return cs_physical_unit_define("W/(m.K)", nullptr);

  else if (strcmp(name, "thermal_expansion") == 0)
    return cs_physical_unit_define("1/K", nullptr);

  else if (strcmp(name, "rho_k_prev") == 0)
    return cs_physical_unit_define("J/m3", nullptr);

  else if (strcmp(name, "kinetic_energy_thermal_st") == 0)
    return cs_physical_unit_define("W/m3", nullptr);

  /* Turbulence fields */
  else if (   strcmp(name, "k") == 0
           || strcmp(name, "rij") == 0
           || strcmp(name, "k_tot") == 0
           || strcmp(name, "k_mod") == 0
           || strcmp(name, "k_res") == 0
           || strcmp(name, "ui2_mean") == 0
           || strcmp(name, "meteo_tke") == 0
           || strcmp(name, "meteo_rij") == 0)
    return cs_physical_unit_define("m2/s2", nullptr);

  else if (   strcmp(name, "epsilon") == 0
           || strcmp(name, "eps_mod") == 0
           || strcmp(name, "meteo_eps") == 0)
    return cs_physical_unit_define("m2/s3", nullptr);

  else if (   strcmp(name, "omega") == 0
           || strcmp(name, "omg_mod") == 0
           || strcmp(name, "f_bar") == 0
           || strcmp(name, "strain_rate_tensor") == 0)
    return cs_physical_unit_define("s-1", nullptr);

  else if (strcmp(name, "s2") == 0)
    return cs_physical_unit_define("s-2", nullptr);

  else if (strcmp(name, "nu_tilda") == 0)
    return cs_physical_unit_define("m2/s", nullptr);

  /* Geometric fields */
  else if (   strcmp(name, "wall_distance") == 0
           || strcmp(name, "z_ground") == 0
           || strcmp(name, "Delta_max") == 0
           || strcmp(name, "xyzp0") == 0)
    return cs_physical_unit_define("m", nullptr);

  else if (strcmp(name, "sleak") == 0)
    return cs_physical_unit_define("m2", nullptr);

  /* Time */
  else if (   strcmp(name, "time") == 0
           || strcmp(name, "age") == 0
           || strcmp(name, "htles_t") == 0)
    return cs_physical_unit_define("s", nullptr);

  /* Molar mass */
  else if (strcmp(name, "xmasmr") == 0)
    return cs_physical_unit_define("kg/mol", nullptr);

  /* Dimensionless fields */
  else if (   strcmp(name, "wall_yplus") == 0
           || strcmp(name, "void_fraction") == 0
           || strcmp(name, "alpha") == 0
           || strcmp(name, "phi") == 0
           || strcmp(name, "courant_number") == 0
           || strcmp(name, "fourier_number") == 0
           || strcmp(name, "volume_courant_number") == 0
           || strcmp(name, "cfl_t") == 0
           || strcmp(name, "cfl_p") == 0
           || strcmp(name, "hybrid_blend") == 0
           || strcmp(name, "meteo_humidity") == 0
           || strcmp(name, "yv") == 0
           || strcmp(name, "htles_psi") == 0
           || strcmp(name, "htles_r") == 0
           || strcmp(name, "htles_icc") == 0
           || strcmp(name, "htles_fs") == 0
           || strcmp(name, "f1_kwsst") == 0
           || strcmp(name, "algo:rij_beta2") == 0
           || strcmp(name, "algo:cmu") == 0
           || strcmp(name, "smagorinsky_constant^2") == 0
           || strcmp(name, "rvsra") == 0
           || strcmp(name, "kleak") == 0)
    return cs_physical_unit_define("-", nullptr);

  /* Mass / Volume fluxes */
  else if (   strcmp(name, "mass_flux") == 0
           || strcmp(name, "boundary_mass_flux") == 0
           || strcmp(name, "inner_mass_flux") == 0)
    return cs_physical_unit_define("kg/s", nullptr);

  else if (   strcmp(name, "volume_flux") == 0
           || strcmp(name, "boundary_volume_flux") == 0)
    return cs_physical_unit_define("m3/s", nullptr);

  else if (strcmp(name, "hybrid_sas_source_term") == 0)
    return cs_physical_unit_define("m2/s4", nullptr);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize physical units dictionary and free allocated resources.
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_units_finalize(void)
{
  if (_units_map != nullptr) {
    for (size_t i = 0; i < _units_desc_size; i++)
      CS_FREE(_units_desc[i]);
    CS_FREE(_units_desc);
    _units_desc = nullptr;
    _units_desc_size = 0;
    cs_map_name_to_id_destroy(&_units_map);
  }
}
