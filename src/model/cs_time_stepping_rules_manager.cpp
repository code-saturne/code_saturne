/*============================================================================
 * cs_time_stepping_rules_manager.cpp
 * 
 * Implementation of the parser for TimeSteppingRules.xml
 *============================================================================*/

#include "cs_time_stepping_rules_manager.h"
#include "base/cs_log.h"
#include "base/cs_file.h"
#include "base/cs_base.h"
#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "gui/cs_tree_xml.h"
#include <cstring>
#include <cstdlib>
#include <cmath>


/*============================================================================
 * Static helpers
 *============================================================================*/

inline bool
cs_time_stepping_rules_manager::_strcmp(const char *s1, const char *s2)
{
  if (s1 == nullptr || s2 == nullptr)
    return false;
  return (std::strcmp(s1, s2) == 0);
}

/*============================================================================
 * Constructor
 *============================================================================*/

cs_time_stepping_rules_manager::cs_time_stepping_rules_manager(const char *rules_xml_path)
  : rules_tree_(nullptr)
{
  bft_printf("\n=== Initializing cs_time_stepping_rules_manager ===\n");
  bft_printf("Loading: %s\n", rules_xml_path);
  
  // Charger le fichier XML
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  
  if (rules_tree_ == nullptr) {
    bft_error(__FILE__, __LINE__, 0,
              _("Error: Could not load TimeSteppingRules.xml from path: %s\n"),
              rules_xml_path);
  }
  
  // Parse the different sections
  parse_definitions_();
  parse_theta_schemes_();
  parse_defaults_();
  parse_validation_rules_();
  parse_mappings_();
  parse_automatic_settings_();
  
  bft_printf("cs_time_stepping_rules_manager initialized successfully.\n\n");
}

/*============================================================================
 * Destructor
 *============================================================================*/

cs_time_stepping_rules_manager::~cs_time_stepping_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Parse Definitions
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_definitions_()
{
  bft_printf("  Parsing Definitions...\n");
  
  // Nothing special to parse for Definitions
  // Types are just documented in the XML
  
  bft_printf("  Definitions parsed.\n");
}

/*============================================================================
 * Parse Theta Schemes
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_theta_schemes_()
{
  bft_printf("  Parsing ThetaSchemes...\n");
  
  cs_tree_node_t *theta_schemes = cs_tree_find_node(rules_tree_, "ThetaSchemes");
  if (theta_schemes == nullptr)
    return;
  
  cs_tree_node_t *config = cs_tree_node_get_child(theta_schemes, "cs_theta_config_t");
  while (config != nullptr) {
    const char *property = cs_tree_node_get_tag(config, "Property");
    
    if (property != nullptr) {
      cs_theta_config_t tc;
      tc.property_name = std::string(property);
      
      // Parser les extrapolations
      cs_tree_node_t *extrap = cs_tree_node_get_child(config, "Extrapolation");
      while (extrap != nullptr) {
        const char *method = cs_tree_node_get_tag(extrap, "Method");
        const char *theta = cs_tree_node_get_tag(extrap, "Theta");
        
        if (method != nullptr && theta != nullptr) {
          int method_enum = -1;
          if (_strcmp(method, "none")) method_enum = 0;
          else if (_strcmp(method, "linear")) method_enum = 1;
          else if (_strcmp(method, "second_order")) method_enum = 2;
          
          if (method_enum >= 0) {
            tc.extrap_method_to_theta[method_enum] = atof(theta);
          }
        }
        
        extrap = cs_tree_node_get_next_of_name(extrap);
      }
      
      // Parser les ordres de termes sources
      cs_tree_node_t *source_order = cs_tree_node_get_child(config, "SourceTermOrder");
      while (source_order != nullptr) {
        const char *order = cs_tree_node_get_tag(source_order, "Order");
        const char *theta = cs_tree_node_get_tag(source_order, "Theta");
        
        if (order != nullptr && theta != nullptr) {
          int order_enum = -1;
          if (_strcmp(order, "explicit")) order_enum = 0;
          else if (_strcmp(order, "semi_implicit_1")) order_enum = 1;
          else if (_strcmp(order, "semi_implicit_2")) order_enum = 2;
          
          if (order_enum >= 0) {
            tc.source_order_to_theta[order_enum] = atof(theta);
          }
        }
        
        source_order = cs_tree_node_get_next_of_name(source_order);
      }
      
      // Parser description
      cs_tree_node_t *desc = cs_tree_find_node(config, "Description");
      if (desc != nullptr) {
        const char *desc_text = cs_tree_node_get_value_str(desc);
        tc.description = (desc_text != nullptr) ? std::string(desc_text) : "";
      }
      
      theta_configs_[std::string(property)] = tc;
      
      bft_printf("    cs_theta_config_t for '%s': %lu extrapolations, %lu source orders\n",
                 property, 
                 tc.extrap_method_to_theta.size(),
                 tc.source_order_to_theta.size());
    }
    
    config = cs_tree_node_get_next_of_name(config);
  }
  
  bft_printf("  ThetaSchemes parsed: %lu configurations\n", theta_configs_.size());
}

/*============================================================================
 * Parse Defaults
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_defaults_()
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
 * Parse Validation Rules
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_validation_rules_()
{
  bft_printf("  Parsing ValidationRules...\n");
  
  cs_tree_node_t *validation = cs_tree_find_node(rules_tree_, "ValidationRules");
  if (validation == nullptr)
    return;
  
  cs_tree_node_t *rule = cs_tree_node_get_child(validation, "Rule");
  while (rule != nullptr) {
    
    cs_tree_node_t *constraint = cs_tree_node_get_child(rule, "Constraint");
    while (constraint != nullptr) {
      const char *target = cs_tree_node_get_tag(constraint, "Target");
      const char *type = cs_tree_node_get_tag(constraint, "Type");
      
      if (target != nullptr && type != nullptr) {
        cs_time_stepping_constraint_t tsc;
        tsc.target = std::string(target);
        tsc.type = std::string(type);
        tsc.min_value = 0.0;
        tsc.max_value = 0.0;
        
        // Parser MinValue/MaxValue
        cs_tree_node_t *min_node = cs_tree_find_node(constraint, "MinValue");
        if (min_node != nullptr) {
          const char *min_val = cs_tree_node_get_value_str(min_node);
          tsc.min_value = (min_val != nullptr) ? atof(min_val) : 0.0;
        }
        
        cs_tree_node_t *max_node = cs_tree_find_node(constraint, "MaxValue");
        if (max_node != nullptr) {
          const char *max_val = cs_tree_node_get_value_str(max_node);
          tsc.max_value = (max_val != nullptr) ? atof(max_val) : 0.0;
        }
        
        // Parser ErrorMessage
        cs_tree_node_t *err_msg = cs_tree_find_node(constraint, "ErrorMessage");
        if (err_msg != nullptr) {
          const char *msg = cs_tree_node_get_value_str(err_msg);
          tsc.error_message = (msg != nullptr) ? std::string(msg) : "";
        }
        
        constraints_.push_back(tsc);
      }
      
      constraint = cs_tree_node_get_next_of_name(constraint);
    }
    
    rule = cs_tree_node_get_next_of_name(rule);
  }
  
  bft_printf("  ValidationRules parsed: %lu constraints\n", constraints_.size());
}

/*============================================================================
 * Parse Mappings
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_mappings_()
{
  bft_printf("  Parsing Mappings...\n");
  
  cs_tree_node_t *mappings = cs_tree_find_node(rules_tree_, "Mappings");
  if (mappings == nullptr)
    return;
  
  // Parser ExtrapolationMethod enum
  cs_tree_node_t *enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='ExtrapolationMethod']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        extrap_method_enum_[std::string(value)] = enum_int;
        enum_to_extrap_method_[enum_int] = std::string(value);
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  // Parser SourceTermOrder enum
  enum_map = cs_tree_find_node(mappings, "EnumMapping[@Name='SourceTermOrder']");
  if (enum_map != nullptr) {
    cs_tree_node_t *entry = cs_tree_node_get_child(enum_map, "Entry");
    while (entry != nullptr) {
      const char *value = cs_tree_node_get_tag(entry, "Value");
      const char *enum_val = cs_tree_node_get_tag(entry, "Enum");
      
      if (value != nullptr && enum_val != nullptr) {
        int enum_int = atoi(enum_val);
        source_order_enum_[std::string(value)] = enum_int;
        enum_to_source_order_[enum_int] = std::string(value);
      }
      
      entry = cs_tree_node_get_next_of_name(entry);
    }
  }
  
  bft_printf("  Mappings parsed.\n");
}

/*============================================================================
 * Parse Automatic Settings
 *============================================================================*/

void
cs_time_stepping_rules_manager::parse_automatic_settings_()
{
  bft_printf("  Parsing AutomaticSettings...\n");
  
  // These settings are already integrated in cs_theta_config_ts
  // This function is for future additions if needed
  
  bft_printf("  AutomaticSettings parsed.\n");
}

/*============================================================================
 * GETTERS POUR cs_parameters.cpp
 *============================================================================*/

double
cs_time_stepping_rules_manager::get_theta_for_extrapolation(const char *property_name, 
                                                        int extrap_method) const
{
  auto it = theta_configs_.find(std::string(property_name));
  if (it != theta_configs_.end()) {
    auto theta_it = it->second.extrap_method_to_theta.find(extrap_method);
    if (theta_it != it->second.extrap_method_to_theta.end())
      return theta_it->second;
  }
  return -999.0;  // Valeur invalide
}

double
cs_time_stepping_rules_manager::get_theta_for_source_term(const char *source_term_name, 
                                                      int source_order) const
{
  auto it = theta_configs_.find(std::string(source_term_name));
  if (it != theta_configs_.end()) {
    auto theta_it = it->second.source_order_to_theta.find(source_order);
    if (theta_it != it->second.source_order_to_theta.end())
      return theta_it->second;
  }
  return -999.0;  // Valeur invalide
}

cs_time_step_limits_t
cs_time_stepping_rules_manager::get_time_step_limits() const
{
  cs_time_step_limits_t limits;
  limits.dt_ref = get_default_double("time_step_ref");
  limits.dt_min = get_default_double("time_step_min");
  limits.dt_max = get_default_double("time_step_max");
  limits.cfl_max = get_default_double("cfl_max");
  limits.courant_max = get_default_double("courant_max");
  limits.fourier_max = get_default_double("fourier_max");
  return limits;
}

bool
cs_time_stepping_rules_manager::should_auto_init_theta(const char *property_name, 
                                                  double current_theta) const
{
  // If theta is -999, it must be auto-initialized
  return (fabs(current_theta + 999.0) < 1e-6);
}

/*============================================================================
 * GETTERS POUR cs_gui.cpp
 *============================================================================*/

const char*
cs_time_stepping_rules_manager::get_default(const char *key) const
{
  auto it = defaults_.find(std::string(key));
  if (it != defaults_.end())
    return it->second.c_str();
  return nullptr;
}

double
cs_time_stepping_rules_manager::get_default_double(const char *key) const
{
  const char *val = get_default(key);
  if (val != nullptr)
    return atof(val);
  return 0.0;
}

int
cs_time_stepping_rules_manager::get_default_int(const char *key) const
{
  const char *val = get_default(key);
  if (val != nullptr)
    return atoi(val);
  return 0;
}

/*============================================================================
 * GETTERS POUR cs_parameters_check.cpp
 *============================================================================*/

double
cs_time_stepping_rules_manager::get_min_value(const char *param_name) const
{
  for (const auto &c : constraints_) {
    if (c.target == std::string(param_name) && c.type == "range")
      return c.min_value;
  }
  return 0.0;
}

double
cs_time_stepping_rules_manager::get_max_value(const char *param_name) const
{
  for (const auto &c : constraints_) {
    if (c.target == std::string(param_name) && c.type == "range")
      return c.max_value;
  }
  return 0.0;
}

bool
cs_time_stepping_rules_manager::validate_dt(double dt_ref, double dt_min, double dt_max) const
{
  if (dt_ref <= 0.0) return false;
  if (dt_min < 0.0) return false;
  if (dt_max <= 0.0) return false;
  if (dt_max <= dt_min) return false;
  return true;
}

bool
cs_time_stepping_rules_manager::validate_theta(double theta) const
{
  return (theta >= 0.0 && theta <= 1.0);
}

bool
cs_time_stepping_rules_manager::validate_cfl(double cfl_max) const
{
  return (cfl_max > 0.0);
}

const char*
cs_time_stepping_rules_manager::get_constraint_error_message(const char *constraint_name) const
{
  for (const auto &c : constraints_) {
    if (c.target == std::string(constraint_name))
      return c.error_message.c_str();
  }
  return nullptr;
}

/*============================================================================
 * ENUM CONVERSIONS
 *============================================================================*/

int
cs_time_stepping_rules_manager::get_extrapolation_method_enum(const char *name) const
{
  auto it = extrap_method_enum_.find(std::string(name));
  if (it != extrap_method_enum_.end())
    return it->second;
  return 0;
}

int
cs_time_stepping_rules_manager::get_source_order_enum(const char *name) const
{
  auto it = source_order_enum_.find(std::string(name));
  if (it != source_order_enum_.end())
    return it->second;
  return 0;
}

const char*
cs_time_stepping_rules_manager::get_extrapolation_method_name(int enum_val) const
{
  auto it = enum_to_extrap_method_.find(enum_val);
  if (it != enum_to_extrap_method_.end())
    return it->second.c_str();
  return "none";
}

const char*
cs_time_stepping_rules_manager::get_source_order_name(int enum_val) const
{
  auto it = enum_to_source_order_.find(enum_val);
  if (it != enum_to_source_order_.end())
    return it->second.c_str();
  return "explicit";
}

const char*
cs_time_stepping_rules_manager::get_extrapolation_method_constant(int enum_val) const
{
  switch (enum_val) {
    case 0: return "CS_TIME_EXTRAP_NONE";
    case 1: return "CS_TIME_EXTRAP_LINEAR";
    case 2: return "CS_TIME_EXTRAP_SECOND_ORDER";
    default: return "CS_TIME_EXTRAP_NONE";
  }
}

const char*
cs_time_stepping_rules_manager::get_source_order_constant(int enum_val) const
{
  switch (enum_val) {
    case 0: return "CS_SOURCE_TERM_EXPLICIT";
    case 1: return "CS_SOURCE_TERM_SEMI_IMPLICIT_1";
    case 2: return "CS_SOURCE_TERM_SEMI_IMPLICIT_2";
    default: return "CS_SOURCE_TERM_EXPLICIT";
  }
}

/*============================================================================
 * Singleton
 *============================================================================*/

static cs_time_stepping_rules_manager *g_timestep_rules_manager = nullptr;

cs_time_stepping_rules_manager*
cs_get_time_stepping_rules_manager()
{
  if (g_timestep_rules_manager == nullptr) {
    // Search for TimeSteppingRules.xml in the installation directory
    char rules_path[1024];
    const char *install_prefix = cs_base_get_pkgdatadir();
    snprintf(rules_path, 1024, "%s/model/TimeSteppingRules.xml", install_prefix);
    
    // Check if the file exists
    if (!cs_file_isreg(rules_path)) {
      // Otherwise, search in the current directory
      snprintf(rules_path, 1024, "TimeSteppingRules.xml");
    }
    
    g_timestep_rules_manager = new cs_time_stepping_rules_manager(rules_path);
  }
  return g_timestep_rules_manager;
}