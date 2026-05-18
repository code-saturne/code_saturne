/*============================================================================
 * cs_thermal_rules_manager.cpp
 * 
 * Implementation of the parser for ThermalRules.xml
 *============================================================================*/

#include "cs_thermal_rules_manager.h"
#include "base/cs_log.h"
#include "base/cs_file.h"
#include <cmath>
#include "base/cs_base.h"
#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include <cstring>
#include <cstdlib>
#include "gui/cs_tree_xml.h"
#include "gui/cs_tree_xml.h"
#include "gui/cs_gui_util.h"

#include "base/cs_field.h"
#include "base/cs_wall_functions.h"

/*============================================================================
 * Static helpers
 *============================================================================*/

inline bool
cs_thermal_rules_manager::_strcmp(const char *s1, const char *s2)
{
  if (s1 == nullptr || s2 == nullptr)
    return false;
  return (std::strcmp(s1, s2) == 0);
}

/*============================================================================
 * Constructor
 *============================================================================*/

cs_thermal_rules_manager::cs_thermal_rules_manager(const char *rules_xml_path)
  : rules_tree_(nullptr)
{
  bft_printf("\n=== Initializing cs_thermal_rules_manager ===\n");
  bft_printf("Loading: %s\n", rules_xml_path);
  
  // Charger le fichier XML (comme cs_turbulence_rules_manager)
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  
  if (rules_tree_ == nullptr) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Could not load ThermalRules.xml from path: %s\n"),
              rules_xml_path);
  }
  
  // Parse the different sections
  parse_definitions_();
  parse_validation_rules_();
  parse_validations_();
  parse_defaults_();
  parse_mappings_();
  
  bft_printf("cs_thermal_rules_manager initialized successfully.\n\n");
}

/*============================================================================
 * Destructor
 *============================================================================*/

cs_thermal_rules_manager::~cs_thermal_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Parse Definitions
 *============================================================================*/

void
cs_thermal_rules_manager::parse_definitions_()
{
  bft_printf("  Parsing Definitions...\n");
  
  // Nothing special to parse for Definitions
  // Types are just documented in the XML
  
  bft_printf("  Definitions parsed.\n");
}

/*============================================================================
 * Parse ValidationRules (field creation)
 *============================================================================*/

void
cs_thermal_rules_manager::parse_validation_rules_()
{
  bft_printf("  Parsing ValidationRules...\n");
  
  cs_tree_node_t *val_rules = cs_tree_find_node(rules_tree_,
                                                  "ValidationRules");
  if (val_rules == nullptr)
    return;
  
  // Parse RequiredFields rules
  cs_tree_node_t *rule = cs_tree_node_get_child(val_rules, "Rule");
  while (rule != nullptr) {
    
    const char *rule_type = cs_tree_node_get_tag(rule, "Type");
    
    if (rule_type != nullptr && _strcmp(rule_type, "RequiredFields")) {
      const char *thermal_var = cs_tree_node_get_tag(rule, "ThermalVariable");
      
      if (thermal_var != nullptr) {
        cs_thermal_variable_rule_t var_rule;
        var_rule.thermal_variable = std::string(thermal_var);
        
        // Parse fields to create
        cs_tree_node_t *create_field = cs_tree_node_get_child(rule, "CreateField");
        while (create_field != nullptr) {
          
          const char *field_name = cs_tree_node_get_tag(create_field, "Name");
          const char *field_type = cs_tree_node_get_tag(create_field, "FieldType");
          const char *location = cs_tree_node_get_tag(create_field, "Location");
          const char *post_vis = cs_tree_node_get_tag(create_field, "PostVis");
          const char *log = cs_tree_node_get_tag(create_field, "Log");
          
          if (field_name != nullptr && field_type != nullptr && location != nullptr) {
            cs_thermal_field_config_t field_cfg;
            field_cfg.field_name = std::string(field_name);
            field_cfg.field_type = std::string(field_type);
            field_cfg.location = std::string(location);
            field_cfg.dimension = 1;  // Default scalar
            field_cfg.post_vis = (post_vis != nullptr && _strcmp(post_vis, "on"));
            field_cfg.log = (log != nullptr && _strcmp(log, "on"));
            
            var_rule.required_fields.push_back(field_cfg);
          }
          
          create_field = cs_tree_node_get_next_of_name(create_field);
        }
        
        // Parser SetFlag
        cs_tree_node_t *set_flag = cs_tree_find_node(rule, "SetFlag");
        if (set_flag != nullptr) {
          const char *flag_name = cs_tree_node_get_tag(set_flag, "Name");
          const char *flag_value = cs_tree_node_get_tag(set_flag, "Value");
          
          if (flag_name != nullptr && _strcmp(flag_name, "is_temperature")) {
            var_rule.is_temperature_flag = (flag_value != nullptr) ? atoi(flag_value) : 0;
          }
        }
        
        // Parser DiffusivityFormula
        cs_tree_node_t *diff_formula = cs_tree_find_node(rule, "DiffusivityFormula");
        if (diff_formula != nullptr) {
          const char *formula = cs_tree_node_get_value_str(diff_formula);
          var_rule.diffusivity_formula = (formula != nullptr) ? std::string(formula) : "";
        }
        
        // Parser ThermoPlane
        cs_tree_node_t *thermo_plane = cs_tree_find_node(rule, "ThermoPlane");
        if (thermo_plane != nullptr) {
          const char *plane = cs_tree_node_get_value_str(thermo_plane);
          var_rule.thermo_plane = (plane != nullptr) ? std::string(plane) : "PT";
        }
        
        // Store the rule
        thermal_variable_rules_[std::string(thermal_var)] = var_rule;
        
        bft_printf("    RequiredFields rule for '%s': %lu fields\n",
                   thermal_var, var_rule.required_fields.size());
      }
    }
    
    // Parse conditional field rules
    else if (rule_type != nullptr && _strcmp(rule_type, "KineticSourceTerm")) {
      std::vector<cs_thermal_field_config_t> fields;
      
      cs_tree_node_t *create_field = cs_tree_node_get_child(rule, "CreateField");
      while (create_field != nullptr) {
        const char *field_name = cs_tree_node_get_tag(create_field, "Name");
        const char *field_type = cs_tree_node_get_tag(create_field, "FieldType");
        const char *location = cs_tree_node_get_tag(create_field, "Location");
        const char *dim = cs_tree_node_get_tag(create_field, "Dim");
        
        if (field_name != nullptr) {
          cs_thermal_field_config_t field_cfg;
          field_cfg.field_name = std::string(field_name);
          field_cfg.field_type = (field_type != nullptr) ? std::string(field_type) : "property";
          field_cfg.location = (location != nullptr) ? std::string(location) : "cells";
          field_cfg.dimension = (dim != nullptr) ? atoi(dim) : 1;
          field_cfg.post_vis = true;
          field_cfg.log = true;
          
          fields.push_back(field_cfg);
        }
        
        create_field = cs_tree_node_get_next_of_name(create_field);
      }
      
      conditional_fields_["kinetic_st"] = fields;
      bft_printf("    KineticSourceTerm: %lu fields\n", fields.size());
    }
    
    else if (rule_type != nullptr && _strcmp(rule_type, "MoistAir")) {
      std::vector<cs_thermal_field_config_t> fields;
      
      cs_tree_node_t *create_field = cs_tree_node_get_child(rule, "CreateField");
      while (create_field != nullptr) {
        const char *field_name = cs_tree_node_get_tag(create_field, "Name");
        
        if (field_name != nullptr) {
          cs_thermal_field_config_t field_cfg;
          field_cfg.field_name = std::string(field_name);
          field_cfg.field_type = "property";
          field_cfg.location = "cells";
          field_cfg.dimension = 1;
          field_cfg.post_vis = true;
          field_cfg.log = true;
          
          fields.push_back(field_cfg);
        }
        
        create_field = cs_tree_node_get_next_of_name(create_field);
      }
      
      conditional_fields_["moist_air"] = fields;
      bft_printf("    MoistAir: %lu fields\n", fields.size());
    }
    
    else if (rule_type != nullptr && _strcmp(rule_type, "CompressibleFlow")) {
      std::vector<cs_thermal_field_config_t> fields;
      
      cs_tree_node_t *create_field = cs_tree_node_get_child(rule, "CreateField");
      while (create_field != nullptr) {
        const char *field_name = cs_tree_node_get_tag(create_field, "Name");
        const char *dim = cs_tree_node_get_tag(create_field, "Dim");
        
        if (field_name != nullptr) {
          cs_thermal_field_config_t field_cfg;
          field_cfg.field_name = std::string(field_name);
          field_cfg.field_type = "internal";
          field_cfg.location = "cells";
          field_cfg.dimension = (dim != nullptr) ? atoi(dim) : 1;
          field_cfg.post_vis = false;
          field_cfg.log = false;
          
          fields.push_back(field_cfg);
        }
        
        create_field = cs_tree_node_get_next_of_name(create_field);
      }
      
      conditional_fields_["compressible"] = fields;
      bft_printf("    CompressibleFlow: %lu fields\n", fields.size());
    }
    
    rule = cs_tree_node_get_next_of_name(rule);
  }
  
  bft_printf("  ValidationRules parsed: %lu thermal variable rules\n",
             thermal_variable_rules_.size());
}

/*============================================================================
 * Parse Validations (contraintes)
 *============================================================================*/

void
cs_thermal_rules_manager::parse_validations_()
{
  bft_printf("  Parsing Validations...\n");
  
  cs_tree_node_t *validations = cs_tree_find_node(rules_tree_, "Validations");
  if (validations == nullptr)
    return;
  
  cs_tree_node_t *constraint = cs_tree_node_get_child(validations, "Constraint");
  while (constraint != nullptr) {
    
    const char *target = cs_tree_node_get_tag(constraint, "Target");
    const char *type = cs_tree_node_get_tag(constraint, "Type");
    
    if (target != nullptr && type != nullptr) {
      cs_thermal_constraint_t tc;
      tc.target = std::string(target);
      tc.type = std::string(type);
      tc.min_value = 0.0;
      tc.max_value = 0.0;
      tc.exclusive_min = false;
      
      // Parser MinValue/MaxValue
      cs_tree_node_t *min_node = cs_tree_find_node(constraint, "MinValue");
      if (min_node != nullptr) {
        const char *min_val = cs_tree_node_get_value_str(min_node);
        tc.min_value = (min_val != nullptr) ? atof(min_val) : 0.0;
      }
      
      cs_tree_node_t *max_node = cs_tree_find_node(constraint, "MaxValue");
      if (max_node != nullptr) {
        const char *max_val = cs_tree_node_get_value_str(max_node);
        tc.max_value = (max_val != nullptr) ? atof(max_val) : 0.0;
      }
      
      // Parser ErrorMessage
      cs_tree_node_t *err_msg = cs_tree_find_node(constraint, "ErrorMessage");
      if (err_msg != nullptr) {
        const char *msg = cs_tree_node_get_value_str(err_msg);
        tc.error_message = (msg != nullptr) ? std::string(msg) : "";
      }
      
      constraints_.push_back(tc);
      
      // Store in maps for fast access
      if (_strcmp(type, "range")) {
        min_values_[std::string(target)] = tc.min_value;
        max_values_[std::string(target)] = tc.max_value;
      }
    }
    
    constraint = cs_tree_node_get_next_of_name(constraint);
  }
  
  bft_printf("  Validations parsed: %lu constraints\n", constraints_.size());
}

/*============================================================================
 * Parse Defaults
 *============================================================================*/

void
cs_thermal_rules_manager::parse_defaults_()
{
  bft_printf("  Parsing Defaults...\n");
  
  cs_tree_node_t *defaults = cs_tree_find_node(rules_tree_, "Defaults");
  if (defaults == nullptr)
    return;
  
  cs_tree_node_t *def = cs_tree_node_get_child(defaults, "Default");
  while (def != nullptr) {
    const char *key = cs_tree_node_get_tag(def, "Key");
    const char *value = cs_tree_node_get_tag(def, "Value");
    
    if (key != nullptr && value != nullptr) {
      defaults_[std::string(key)] = std::string(value);
    }
    
    def = cs_tree_node_get_next_of_name(def);
  }
  
  bft_printf("  Defaults parsed: %lu entries\n", defaults_.size());
}

/*============================================================================
 * Parse Mappings (enums)
 *============================================================================*/

void
cs_thermal_rules_manager::parse_mappings_()
{
  bft_printf("  Parsing Mappings...\n");
  
  cs_tree_node_t *mappings = cs_tree_find_node(rules_tree_, "Mappings");
  if (mappings == nullptr)
    return;
  
  // Parser ThermalVariable enum
  cs_tree_node_t *enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='ThermalVariable']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        thermal_variable_enum_[std::string(value)] = enum_int;
        enum_to_thermal_var_[enum_int] = std::string(value);
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  // Parser EquationOfState enum
  enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='EquationOfState']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        eos_enum_[std::string(value)] = enum_int;
        enum_to_eos_[enum_int] = std::string(value);
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  // Parser TemperatureScale enum
  enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='TemperatureScale']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        temp_scale_enum_[std::string(value)] = enum_int;
        enum_to_temp_scale_[enum_int] = std::string(value);
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  // Parser ThermoPlane enum
  enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='ThermoPlane']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        thermo_plane_enum_[std::string(value)] = enum_int;
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  bft_printf("  Mappings parsed.\n");
}

/*============================================================================
 * GETTERS POUR cs_parameters.cpp
 *============================================================================*/

const cs_thermal_variable_rule_t*
cs_thermal_rules_manager::get_thermal_variable_rule(const char *thermal_var) const
{
  auto it = thermal_variable_rules_.find(std::string(thermal_var));
  if (it != thermal_variable_rules_.end())
    return &(it->second);
  return nullptr;
}

bool
cs_thermal_rules_manager::requires_kinetic_st_fields(int has_kinetic_st) const
{
  return (has_kinetic_st == 1 && conditional_fields_.count("kinetic_st") > 0);
}

bool
cs_thermal_rules_manager::requires_moist_air_field(int ieos) const
{
  // CS_EOS_MOIST_AIR = 2
  return (ieos == 2 && conditional_fields_.count("moist_air") > 0);
}

bool
cs_thermal_rules_manager::requires_compressible_fields(int ieos, int thermal_var) const
{
  // ieos != CS_EOS_NONE (0) ET (thermal_var == TEMPERATURE (1) OU INTERNAL_ENERGY (3))
  return (ieos != 0 && (thermal_var == 1 || thermal_var == 3) &&
          conditional_fields_.count("compressible") > 0);
}

const std::vector<cs_thermal_field_config_t>*
cs_thermal_rules_manager::get_kinetic_st_fields() const
{
  auto it = conditional_fields_.find("kinetic_st");
  if (it != conditional_fields_.end())
    return &(it->second);
  return nullptr;
}

const std::vector<cs_thermal_field_config_t>*
cs_thermal_rules_manager::get_moist_air_fields() const
{
  auto it = conditional_fields_.find("moist_air");
  if (it != conditional_fields_.end())
    return &(it->second);
  return nullptr;
}

const std::vector<cs_thermal_field_config_t>*
cs_thermal_rules_manager::get_compressible_fields() const
{
  auto it = conditional_fields_.find("compressible");
  if (it != conditional_fields_.end())
    return &(it->second);
  return nullptr;
}

const char*
cs_thermal_rules_manager::get_diffusivity_formula(const char *thermal_var) const
{
  auto it = thermal_variable_rules_.find(std::string(thermal_var));
  if (it != thermal_variable_rules_.end())
    return it->second.diffusivity_formula.c_str();
  return nullptr;
}

const char*
cs_thermal_rules_manager::get_thermo_plane(const char *thermal_var) const
{
  auto it = thermal_variable_rules_.find(std::string(thermal_var));
  if (it != thermal_variable_rules_.end())
    return it->second.thermo_plane.c_str();
  return "PT";
}

int
cs_thermal_rules_manager::get_thermo_plane_enum(const char *thermal_var) const
{
  const char *plane = get_thermo_plane(thermal_var);
  if (plane != nullptr) {
    auto it = thermo_plane_enum_.find(std::string(plane));
    if (it != thermo_plane_enum_.end())
      return it->second;
  }
  return 0;  // CS_PHYS_PROP_PLANE_PT by default
}

/*============================================================================
 * GETTERS POUR cs_parameters_check.cpp
 *============================================================================*/

double
cs_thermal_rules_manager::get_min_value(const char *param_name) const
{
  auto it = min_values_.find(std::string(param_name));
  if (it != min_values_.end())
    return it->second;
  return 0.0;
}

double
cs_thermal_rules_manager::get_max_value(const char *param_name) const
{
  auto it = max_values_.find(std::string(param_name));
  if (it != max_values_.end())
    return it->second;
  return 0.0;
}

bool
cs_thermal_rules_manager::check_constraint(const char *constraint_name,
                                       int thermal_var,
                                       int temp_scale,
                                       int has_gravity,
                                       int density_variable) const
{
  // Rechercher la contrainte
  for (const auto &c : constraints_) {
    if (c.target == std::string(constraint_name)) {
      // Implement checks according to type
      // TODO: To be completed as needed
      return true;
    }
  }
  return false;
}

const char*
cs_thermal_rules_manager::get_constraint_error_message(const char *constraint_name) const
{
  for (const auto &c : constraints_) {
    if (c.target == std::string(constraint_name))
      return c.error_message.c_str();
  }
  return nullptr;
}

bool
cs_thermal_rules_manager::is_lagrangian_second_order_forbidden(int thermal_var,
                                                            int temp_scale) const
{
  // Interdit si (TEMPERATURE + KELVIN) OU (ENTHALPY)
  // thermal_var: 1=TEMPERATURE, 2=ENTHALPY
  // temp_scale: 1=KELVIN
  
  if (thermal_var == 1 && temp_scale == 1)  // TEMPERATURE + KELVIN
    return true;
  if (thermal_var == 2)  // ENTHALPY
    return true;
  
  return false;
}

/*============================================================================
 * GETTERS POUR cs_gui.cpp
 *============================================================================*/

const char*
cs_thermal_rules_manager::get_property_method(const char *property_name) const
{
  auto it = property_methods_.find(std::string(property_name));
  if (it != property_methods_.end())
    return it->second.c_str();
  return "constant";
}

/*============================================================================
 * ENUM CONVERSIONS
 *============================================================================*/

int
cs_thermal_rules_manager::get_thermal_variable_enum(const char *name) const
{
  auto it = thermal_variable_enum_.find(std::string(name));
  if (it != thermal_variable_enum_.end())
    return it->second;
  return 0;  // CS_THERMAL_MODEL_NONE
}

int
cs_thermal_rules_manager::get_eos_enum(const char *name) const
{
  auto it = eos_enum_.find(std::string(name));
  if (it != eos_enum_.end())
    return it->second;
  return 0;  // CS_EOS_NONE
}

int
cs_thermal_rules_manager::get_temp_scale_enum(const char *name) const
{
  auto it = temp_scale_enum_.find(std::string(name));
  if (it != temp_scale_enum_.end())
    return it->second;
  return 1;  // CS_TEMPERATURE_SCALE_KELVIN by default
}

int
cs_thermal_rules_manager::get_thermo_plane_enum_by_name(const char *name) const
{
  auto it = thermo_plane_enum_.find(std::string(name));
  if (it != thermo_plane_enum_.end())
    return it->second;
  return 0;
}

const char*
cs_thermal_rules_manager::get_thermal_variable_name(int enum_val) const
{
  auto it = enum_to_thermal_var_.find(enum_val);
  if (it != enum_to_thermal_var_.end())
    return it->second.c_str();
  return "none";
}

const char*
cs_thermal_rules_manager::get_eos_name(int enum_val) const
{
  auto it = enum_to_eos_.find(enum_val);
  if (it != enum_to_eos_.end())
    return it->second.c_str();
  return "none";
}

const char*
cs_thermal_rules_manager::get_temp_scale_name(int enum_val) const
{
  auto it = enum_to_temp_scale_.find(enum_val);
  if (it != enum_to_temp_scale_.end())
    return it->second.c_str();
  return "none";
}

const char*
cs_thermal_rules_manager::get_thermal_variable_constant(int enum_val) const
{
  // Mapper enum → nom de constante C++
  switch (enum_val) {
    case 0: return "CS_THERMAL_MODEL_NONE";
    case 1: return "CS_THERMAL_MODEL_TEMPERATURE";
    case 2: return "CS_THERMAL_MODEL_ENTHALPY";
    case 3: return "CS_THERMAL_MODEL_INTERNAL_ENERGY";
    default: return "CS_THERMAL_MODEL_NONE";
  }
}

const char*
cs_thermal_rules_manager::get_eos_constant(int enum_val) const
{
  switch (enum_val) {
    case 0: return "CS_EOS_NONE";
    case 1: return "CS_EOS_IDEAL_GAS";
    case 2: return "CS_EOS_MOIST_AIR";
    case 3: return "CS_EOS_REAL_GAS";
    default: return "CS_EOS_NONE";
  }
}

const char*
cs_thermal_rules_manager::get_temp_scale_constant(int enum_val) const
{
  switch (enum_val) {
    case 0: return "CS_TEMPERATURE_SCALE_NONE";
    case 1: return "CS_TEMPERATURE_SCALE_KELVIN";
    case 2: return "CS_TEMPERATURE_SCALE_CELSIUS";
    default: return "CS_TEMPERATURE_SCALE_KELVIN";
  }
}

/*============================================================================
 * DEFAULTS
 *============================================================================*/

const char*
cs_thermal_rules_manager::get_default(const char *key) const
{
  auto it = defaults_.find(std::string(key));
  if (it != defaults_.end())
    return it->second.c_str();
  return nullptr;
}

double
cs_thermal_rules_manager::get_default_double(const char *key) const
{
  const char *val = get_default(key);
  if (val != nullptr)
    return atof(val);
  return 0.0;
}

/*============================================================================
 * Singleton
 *============================================================================*/

static cs_thermal_rules_manager *g_thermal_rules_manager = nullptr;

cs_thermal_rules_manager*
cs_get_thermal_rules_manager()
{
  if (g_thermal_rules_manager == nullptr) {
    // Search for ThermalRules.xml in the installation directory
    char rules_path[1024];
    const char *install_prefix = cs_base_get_pkgdatadir();
    snprintf(rules_path, 1024, "%s/model/ThermalRules.xml", install_prefix);
    
    // Check if the file exists
    if (!cs_file_isreg(rules_path)) {
      // Otherwise, search in the current directory
      snprintf(rules_path, 1024, "ThermalRules.xml");
    }
    
    g_thermal_rules_manager = new cs_thermal_rules_manager(rules_path);
  }
  return g_thermal_rules_manager;
}