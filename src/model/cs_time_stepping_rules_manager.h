#ifndef CS_TIME_STEPPING_RULES_MANAGER_H
#define CS_TIME_STEPPING_RULES_MANAGER_H

/*============================================================================
 * Parse and manage TimeSteppingRules.xml
 * Manages time stepping schemes and theta schemes
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

/* Theta configuration for a property */
struct cs_theta_config_t {
  std::string property_name;
  std::map<int, double> extrap_method_to_theta;  // extrapolation method → theta
  std::map<int, double> source_order_to_theta;   // source term order → theta
  std::string description;
};

/* Time step limitss */
struct cs_time_step_limits_t {
  double dt_ref;
  double dt_min;
  double dt_max;
  double cfl_max;
  double courant_max;
  double fourier_max;
};

/* Validation constraints */
struct cs_time_stepping_constraint_t {
  std::string target;
  std::string type;
  double min_value;
  double max_value;
  std::string error_message;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_time_stepping_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  /*--------------------------------------------------------------------------
   * Maps for cs_parameters.cpp
   *--------------------------------------------------------------------------*/

  // Theta configurations per property
  std::map<std::string, cs_theta_config_t> theta_configs_;

  // Default values
  std::map<std::string, std::string> defaults_;

  // Validation constraints
  std::vector<cs_time_stepping_constraint_t> constraints_;

  /*--------------------------------------------------------------------------
   * Maps for enums
   *--------------------------------------------------------------------------*/

  std::map<std::string, int> extrap_method_enum_;
  std::map<std::string, int> source_order_enum_;

  std::map<int, std::string> enum_to_extrap_method_;
  std::map<int, std::string> enum_to_source_order_;

  // Helper
  static inline bool
  _strcmp(const char  *s1,
          const char  *s2);

  // Parsing methods
  void parse_definitions_();
  void parse_theta_schemes_();
  void parse_defaults_();
  void parse_validation_rules_();
  void parse_mappings_();
  void parse_automatic_settings_();

public:

  // Constructor & Destructor
  cs_time_stepping_rules_manager(const char *rules_xml_path);
  ~cs_time_stepping_rules_manager();

  /*--------------------------------------------------------------------------
   * Getters for cs_parameters.cpp
   *--------------------------------------------------------------------------*/

  // Get theta according to extrapolation method
  double
  get_theta_for_extrapolation(const char  *property_name,
                              int          extrap_method) const;

  // Get theta according to source term order
  double
  get_theta_for_source_term(const char  *source_term_name,
                            int          source_order) const;

  // Get time step limits
  cs_time_step_limits_t
  get_time_step_limits() const;

  // Check if theta should be auto-initialized
  bool
  should_auto_init_theta(const char  *property_name,
                         double       current_theta) const;

  /*--------------------------------------------------------------------------
   * Getters for cs_gui.cpp
   *--------------------------------------------------------------------------*/

  // Get default value
  const char* get_default(const char *key) const;
  double get_default_double(const char *key) const;
  int get_default_int(const char *key) const;

  /*--------------------------------------------------------------------------
   * Getters for cs_parameters_check.cpp
   *--------------------------------------------------------------------------*/

  // Validation ranges
  double
  get_min_value(const char  *param_name) const;

  double
  get_max_value(const char  *param_name) const;

  // Check constraint
  bool validate_dt(double  dt_ref,
                   double  dt_min,
                   double  dt_max) const;
  bool validate_theta(double  theta) const;
  bool validate_cfl(double  cfl_max) const;

  // Error messages
  const char* get_constraint_error_message(const char *constraint_name) const;

  /*--------------------------------------------------------------------------
   * Conversion enums
   *--------------------------------------------------------------------------*/

  int
  get_extrapolation_method_enum(const char  *name) const;

  int
  get_source_order_enum(const char  *name) const;

  const char *
  get_extrapolation_method_name(int  enum_val) const;

  const char *
  get_source_order_name(int  enum_val) const;

  // Get C++ constant name
  const char *
  get_extrapolation_method_constant(int  enum_val) const;

  const char *
  get_source_order_constant(int enum_val) const;
};

/*============================================================================
 * Singleton accessor (thread-safe)
 *============================================================================*/

cs_time_stepping_rules_manager* cs_get_time_stepping_rules_manager();

/*----------------------------------------------------------------------------*/

#endif /* CS_TIME_STEPPING_RULES_MANAGER_H */
