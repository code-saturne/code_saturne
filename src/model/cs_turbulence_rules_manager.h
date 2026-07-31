#ifndef CS_TURBULENCE_RULES_MANAGER_H
#define CS_TURBULENCE_RULES_MANAGER_H

/*============================================================================
 * Class to parse and manage TurbulenceRules.xml
 * Includes validation rules management for cs_parameters_check.cpp
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

#include <map>
#include <string>
#include <vector>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_tree.h"

/*============================================================================
 * Structure definitions for Output module
 *============================================================================*/

/* Structure for Output field default values */
struct cs_output_field_defaults_t {
  bool listing;
  bool postprocessing;
  bool probes;
};

/* Structure pour la config des variables surfaciques */
struct cs_boundary_var_config_t {
  std::string label;
  bool default_status;
  std::vector<std::string> forbidden_models;
  bool requires_thermal;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_turbulence_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  /*-------------------------------------------------------------------------
   * Maps for cs_gui.cpp
   *-------------------------------------------------------------------------*/
  std::map<std::string, int> model_enum_map_;
  std::map<std::string, bool> requires_wall_function_;
  std::map<std::string, bool> requires_gravity_terms_;
  std::map<std::string, bool> requires_mixing_length_;
  std::map<std::string, std::string> forced_gravity_;
  std::map<std::string, bool> is_les_model_;
  std::map<std::string, bool> is_rij_model_;

  /*-------------------------------------------------------------------------
   * Maps for cs_parameters.cpp
   *-------------------------------------------------------------------------*/
  std::map<std::string, std::string> default_wall_function_type_;
  std::map<std::string, double> yplus_limit_;
  std::map<std::string, std::string> diffusion_model_;
  std::map<int, std::string> enum_to_name_;

  /*-------------------------------------------------------------------------
   * Maps for cs_parameters_check.cpp
   *-------------------------------------------------------------------------*/
  std::map<std::string, int> validation_min_;
  std::map<std::string, int> validation_max_;
  std::map<std::string, std::vector<int>> compatible_itytur_values_;
  std::map<std::string, std::vector<std::string>> models_requiring_zero_;
  std::map<std::string, std::vector<std::string>> not_available_for_models_;

  /*-------------------------------------------------------------------------
   * Output data
   *-------------------------------------------------------------------------*/
  std::map<std::string, cs_output_field_defaults_t> output_field_defaults_;
  std::map<std::string, cs_boundary_var_config_t> output_boundary_var_config_;
  std::map<std::string, std::vector<std::string>> output_tensor_component_names_;
  std::map<std::string, std::string> output_defaults_;

  // Helper
  static inline bool _strcmp(const char *s1, const char *s2);

  // Parsing methods
  void parse_model_groups_();
  void parse_validation_rules_();
  void parse_model_requirements_();
  void build_model_enum_map_();
  void parse_numerical_parameters_();
  void parse_validation_rules_check_();

  // Output parsing (NOUVEAU)
  void parse_output_module_();
  void parse_output_field_defaults_();
  void parse_output_boundary_variables_();
  void parse_output_mappings_();

public:

  // Constructor & Destructor
  cs_turbulence_rules_manager(const char *rules_xml_path);
  ~cs_turbulence_rules_manager();

  /*-------------------------------------------------------------------------
   * Getters for cs_gui.cpp
   *-------------------------------------------------------------------------*/

  int
  get_model_enum(const char  *model_name) const;

  bool
  requires_wall_function(const char  *model_name) const;

  bool
  requires_gravity_terms(const char  *model_name) const;

  bool
  requires_mixing_length(const char  *model_name) const;

  const char *
  get_forced_gravity(const char  *model_name) const;

  bool
  is_les_model(const char  *model_name) const;

  bool
  is_rij_model(const char  *model_name) const;

  bool
  supports_coupled_rij(const char  *model_name) const;

  /*-------------------------------------------------------------------------
   * Getters for cs_parameters.cpp
   *-------------------------------------------------------------------------*/

  const char*
  get_model_name_from_enum(int  model_enum) const;

  const char*
  get_default_wall_function_type(const char  *model_name) const;

  double
  get_yplus_limit(const char  *model_name) const;

  const char*
  get_diffusion_model(const char  *model_name) const;

  int
  get_wall_function_enum(const char  *wall_function_type) const;

  /*-------------------------------------------------------------------------
   * Getters for cs_parameters_check.cpp
   *-------------------------------------------------------------------------*/

  // Get validation ranges
  int
  get_validation_min(const char *param_name) const;
  int
  get_validation_max(const char  *param_name) const;

  // Get compatible itytur values
  const std::vector<int> *
  get_compatible_itytur_values(const char  *param_name) const;

  // Check if a parameter must be 0 for a given model
  bool
  requires_zero_value(const char  *param_name,
                      int          model_enum) const;

  // Check if a parameter is available for a given model
  bool
  is_available_for_model(const char  *param_name,
                         int          model_enum) const;

  // Helper: convert enum to cpnstant name
  std::string
  get_model_constant_name(int model_enum) const;

  /*-------------------------------------------------------------------------
   * Output getters
   *-------------------------------------------------------------------------*/

  cs_output_field_defaults_t
  get_output_field_defaults(const char  *field_type) const;

  cs_boundary_var_config_t
  get_output_boundary_var_config(const char  *var_name) const;
  bool
  is_output_boundary_var_forbidden(const char  *var_name,
                                   const char  *physics_model) const;
  const std::vector<std::string> *
  get_output_tensor_component_names(const char  *tensor_name) const;

  const char*
  get_output_default(const char *key) const;

  bool
  should_disable_yplus(const char *turb_model) const;
};

/*============================================================================
 * Singleton accessor (thread-safe)
 *============================================================================*/

cs_turbulence_rules_manager *
cs_get_turbulence_rules_manager(bool  no_instanciate = false);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_RULES_MANAGER_H */
