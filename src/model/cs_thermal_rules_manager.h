#ifndef CS_THERMAL_RULES_MANAGER_H
#define CS_THERMAL_RULES_MANAGER_H

/*============================================================================
 * Parse and manage ThermalRules.xml
 * Manages validation rules for the thermal model
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <map>
#include <string>
#include <vector>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_tree.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Structure for required fields based on thermal_variable */
struct cs_thermal_field_config_t {
  std::string  field_name;
  std::string  field_type;        // variable, property, boundary, internal
  std::string  location;          // cells, boundary_faces, interior_faces
  int          dimension;         // 1 for scalae, 3 for vector
  bool         post_vis;          // TODO: use int as several bits are possible
  bool         log;
};

/* Structure for field creation rules */
struct cs_thermal_variable_rule_t {
  std::string thermal_variable;  // temperature, enthalpy, internal_energy
  std::vector<cs_thermal_field_config_t> required_fields;
  int is_temperature_flag;       // 0 or 1
  std::string diffusivity_formula;  // "lambda" or "lambda/cp"
  std::string thermo_plane;      // "PT", "PH", "PS"
};

/* Structure for validation constraints */
struct cs_thermal_constraint_t {
  std::string target;
  std::string               type;         /* range, positive_strict, forbidden,
                                             requires, warning */
  double                    min_value;
  double                    max_value;
  bool                      exclusive_min;
  std::string               error_message;
  std::vector<std::string>  conditions;
  std::vector<std::string>  required_flags;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_thermal_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  /*--------------------------------------------------------------------------
   * Maps for cs_parameters.cpp
   *--------------------------------------------------------------------------*/

  // Field creation rules per thermal_variable
  std::map<std::string, cs_thermal_variable_rule_t> thermal_variable_rules_;

  // Additionsl conditional fields
  std::map<std::string, std::vector<cs_thermal_field_config_t>> conditional_fields_;
  // Keys: "kinetic_st", "moist_air", "compressible_temp_energy"

  // Diffusivity formulas
  std::map<std::string, std::string> diffusivity_formulas_;  // thermal_var → formule

  // Plans thermodynamiques
  std::map<std::string, std::string> thermo_planes_;  // thermal_var → plane

  /*--------------------------------------------------------------------------
   * Maps for cs_parameters_check.cpp
   *--------------------------------------------------------------------------*/

  // Validation constraints
  std::vector<cs_thermal_constraint_t> constraints_;

  // Value ranges
  std::map<std::string, double> min_values_;
  std::map<std::string, double> max_values_;

  /*--------------------------------------------------------------------------
   * Maps for cs_gui.cpp
   *--------------------------------------------------------------------------*/

  // Property calculation methods
  std::map<std::string, std::string> property_methods_;

  /*--------------------------------------------------------------------------
   * Enum mappings
   *--------------------------------------------------------------------------*/

  std::map<std::string, int> thermal_variable_enum_;
  std::map<std::string, int> eos_enum_;
  std::map<std::string, int> temp_scale_enum_;
  std::map<std::string, int> thermo_plane_enum_;

  std::map<int, std::string> enum_to_thermal_var_;
  std::map<int, std::string> enum_to_eos_;
  std::map<int, std::string> enum_to_temp_scale_;

  /*--------------------------------------------------------------------------
   * Defaults
   *--------------------------------------------------------------------------*/

  std::map<std::string, std::string> defaults_;

  // Helper
  static inline bool
  _strcmp(const char *s1, const char *s2);

  // Parsing methods
  void
  parse_definitions_();

  void
  parse_validation_rules_();

  void
  parse_validations_();

  void
  parse_defaults_();

  void
  parse_mappings_();

public:

  // Constructor & Destructor
  cs_thermal_rules_manager(const char *rules_xml_path);
  ~cs_thermal_rules_manager();

  /*--------------------------------------------------------------------------
   * Getters for cs_parameters.cpp
   *--------------------------------------------------------------------------*/

  // Get complete configuration for a thermal variable
  const cs_thermal_variable_rule_t *
  get_thermal_variable_rule(const char  *thermal_var) const;

  // Check if additional fields are required
  bool
  requires_kinetic_st_fields(int  has_kinetic_st) const;

  bool
  requires_moist_air_field(int  ieos) const;

  bool
  requires_compressible_fields(int  ieos,
                               int  thermal_var) const;

  // Get additional fields
  const std::vector<cs_thermal_field_config_t> *
  get_kinetic_st_fields() const;

  const std::vector<cs_thermal_field_config_t> *
  get_moist_air_fields() const;

  const std::vector<cs_thermal_field_config_t> *
  get_compressible_fields() const;

  // Get diffusivity formula
  const char *
  get_diffusivity_formula(const char  *thermal_var) const;

  // Get thermodynamic plane
  const char *
  get_thermo_plane(const char  *thermal_var) const;

  int
  get_thermo_plane_enum(const char  *thermal_var) const;

  /*--------------------------------------------------------------------------
   * Getters for cs_parameters_check.cpp
   *--------------------------------------------------------------------------*/

  // Validation ranges
  double
  get_min_value(const char  *param_name) const;

  double
  get_max_value(const char  *param_name) const;

  // Check constraints
  bool
  check_constraint(const char  *constraint_name,
                   int          thermal_var,
                   int          temp_scale,
                   int          has_gravity,
                   int          density_variable) const;

  // Get error message
  const char *
  get_constraint_error_message(const char  *constraint_name) const;

  // Check if Lagrangian order 2 is forbidden
  bool
  is_lagrangian_second_order_forbidden(int  thermal_var,
                                       int  temp_scale) const;

  /*--------------------------------------------------------------------------
   * Getters for cs_gui.cpp
   *--------------------------------------------------------------------------*/

  // Property calculation methods
  const char* get_property_method(const char *property_name) const;

  /*--------------------------------------------------------------------------
   * Enum conversions
   *--------------------------------------------------------------------------*/

  int
  get_thermal_variable_enum(const char  *name) const;

  int
  get_eos_enum(const char  *name) const;

  int
  get_temp_scale_enum(const char  *name) const;

  int
  get_thermo_plane_enum_by_name(const char  *name) const;

  const char *
  get_thermal_variable_name(int  enum_val) const;

  const char *
  get_eos_name(int  enum_val) const;

  const char *
  get_temp_scale_name(int  enum_val) const;

  // Get C++ constant name
  const char *
  get_thermal_variable_constant(int  enum_val) const;

  const char *
  get_eos_constant(int  enum_val) const;

  const char *
  get_temp_scale_constant(int  enum_val) const;

  /*--------------------------------------------------------------------------
   * Defauls
   *--------------------------------------------------------------------------*/

  const char *
  get_default(const char  *key) const;

  double
  get_default_double(const char  *key) const;
};

/*============================================================================
 * Singleton accessor (thread-safe)
 *============================================================================*/

cs_thermal_rules_manager *
cs_get_thermal_rules_manager(bool  no_instanciate = false);

/*----------------------------------------------------------------------------*/

#endif /* CS_THERMAL_RULES_MANAGER_H */
