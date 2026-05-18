/*============================================================================
 * cs_turbulence_rules_manager.cpp
 * 
 * Complete implementation with validation rules
 * 
 * 
 * Date: February 2026
 *============================================================================*/

#include "cs_turbulence_rules_manager.h"
#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "gui/cs_tree_xml.h"
#include "gui/cs_gui_util.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_field.h"
#include "base/cs_wall_functions.h"
#include <cstring>
#include <cstdlib>

/*============================================================================
 * Static global variable (singleton instance)
 *============================================================================*/

static cs_turbulence_rules_manager *g_rules_manager = nullptr;

/*============================================================================
 * Singleton accessor
 *============================================================================*/

cs_turbulence_rules_manager* cs_get_turbulence_rules_manager()
{
  if (g_rules_manager == nullptr) {
    char rules_path[1024];
    
    const char *datadir = cs_base_get_pkgdatadir();
    snprintf(rules_path, 1024, "%s/model/TurbulenceRules.xml", datadir);
    FILE *test_file = fopen(rules_path, "r");
    
    
    
    if (test_file != nullptr) {
      fclose(test_file);
    } else {
      bft_printf("WARNING: TurbulenceRules.xml not found\n");
      snprintf(rules_path, 1024, "TurbulenceRules.xml");
    }
    
    bft_printf("\n==========================================================\n");
    bft_printf("  CHARGEMENT TURBULENCE RULES DEPUIS XML\n");
    bft_printf("==========================================================\n");
    bft_printf("Chemin: %s\n", rules_path);
    
    g_rules_manager = new cs_turbulence_rules_manager(rules_path);
  }
  
  return g_rules_manager;
}

/*============================================================================
 * Private methods
 *============================================================================*/

bool cs_turbulence_rules_manager::_strcmp(const char *s1, const char *s2)
{
  if (s1 == nullptr || s2 == nullptr) return false;
  return strcmp(s1, s2) == 0;
}

void cs_turbulence_rules_manager::parse_model_groups_()
{
  cs_tree_node_t *groups_node = cs_tree_find_node(rules_tree_, "ModelGroups");
  if (groups_node == nullptr) return;
  
  cs_tree_node_t *group = cs_tree_node_get_child(groups_node, "Group");
  while (group != nullptr) {
    const char *group_name = cs_tree_node_get_tag(group, "Name");
    if (group_name != nullptr) {
      cs_tree_node_t *model = cs_tree_node_get_child(group, "Model");
      while (model != nullptr) {
        const char *model_name = cs_tree_node_get_value_str(model);
        if (model_name != nullptr) {
          if (_strcmp(group_name, "LES")) {
            is_les_model_[model_name] = true;
          }
          else if (_strcmp(group_name, "RijModels")) {
            is_rij_model_[model_name] = true;
          }
        }
        model = cs_tree_node_get_next_of_name(model);
      }
    }
    group = cs_tree_node_get_next_of_name(group);
  }
}

void cs_turbulence_rules_manager::parse_validation_rules_()
{
  cs_tree_node_t *rules_node = cs_tree_find_node(rules_tree_, "ValidationRules");
  if (rules_node == nullptr) return;
  
  cs_tree_node_t *rule = cs_tree_node_get_child(rules_node, "Rule");
  while (rule != nullptr) {
    const char *rule_type = cs_tree_node_get_tag(rule, "Type");
    const char *model_name = cs_tree_node_get_tag(rule, "Model");
    
    if (model_name != nullptr && rule_type != nullptr) {
      if (_strcmp(rule_type, "GravityConstraint")) {
        const char *force_value = 
          cs_tree_node_get_child_value_str(rule, "ForceValue");
        if (force_value != nullptr) {
          forced_gravity_[model_name] = force_value;
        }
      }
    }
    rule = cs_tree_node_get_next_of_name(rule);
  }
}

void cs_turbulence_rules_manager::parse_model_requirements_()
{
  cs_tree_node_t *rules_node = cs_tree_find_node(rules_tree_, "ValidationRules");
  if (rules_node == nullptr) return;
  
  cs_tree_node_t *rule = cs_tree_node_get_child(rules_node, "Rule");
  while (rule != nullptr) {
    const char *rule_type = cs_tree_node_get_tag(rule, "Type");
    const char *model_name = cs_tree_node_get_tag(rule, "Model");
    
    if (model_name != nullptr && rule_type != nullptr) {
      if (_strcmp(rule_type, "RequiredComponents")) {
        cs_tree_node_t *mixing = cs_tree_find_node_simple(rule, "RequireMixingLength");
        if (mixing != nullptr) {
          requires_mixing_length_[model_name] = true;
        }
        
        cs_tree_node_t *var = cs_tree_find_node_simple(rule, "RequireVariable");
        if (var != nullptr) {
          if (!is_les_model_[model_name]) {
            requires_wall_function_[model_name] = true;
            requires_gravity_terms_[model_name] = true;
          }
        }
      }
    }
    rule = cs_tree_node_get_next_of_name(rule);
  }
}

void cs_turbulence_rules_manager::build_model_enum_map_()
{
  model_enum_map_["off"]               = CS_TURB_NONE;
  model_enum_map_["mixing_length"]     = CS_TURB_MIXING_LENGTH;
  model_enum_map_["k-epsilon"]         = CS_TURB_K_EPSILON;
  model_enum_map_["k-epsilon-PL"]      = CS_TURB_K_EPSILON_LIN_PROD;
  model_enum_map_["Rij-epsilon"]       = CS_TURB_RIJ_EPSILON_LRR;
  model_enum_map_["Rij-SSG"]           = CS_TURB_RIJ_EPSILON_SSG;
  model_enum_map_["Rij-EBRSM"]         = CS_TURB_RIJ_EPSILON_EBRSM;
  model_enum_map_["LES_Smagorinsky"]   = CS_TURB_LES_SMAGO_CONST;
  model_enum_map_["LES_dynamique"]     = CS_TURB_LES_SMAGO_DYN;
  model_enum_map_["LES_WALE"]          = CS_TURB_LES_WALE;
  model_enum_map_["v2f-phi"]           = CS_TURB_V2F_PHI;
  model_enum_map_["v2f-BL-v2/k"]       = CS_TURB_V2F_BL_V2K;
  model_enum_map_["k-omega-SST"]       = CS_TURB_K_OMEGA;
  model_enum_map_["Spalart-Allmaras"]  = CS_TURB_SPALART_ALLMARAS;
}

void cs_turbulence_rules_manager::parse_numerical_parameters_()
{
  cs_tree_node_t *models_node = cs_tree_find_node(rules_tree_, "TurbulenceModels");
  if (models_node == nullptr) {
    bft_printf("WARNING: TurbulenceModels node not found in XML\n");
    return;
  }
  
  cs_tree_node_t *model = cs_tree_node_get_child(models_node, "TurbulenceModel");
  
  while (model != nullptr) {
    const char *model_name = cs_tree_node_get_tag(model, "Name");
    
    if (model_name != nullptr) {
      int model_enum = 0;
      cs_tree_node_t *enum_node = cs_tree_get_node(model, "Enum");
      if (enum_node != nullptr) {
        const char *enum_str = cs_tree_node_get_value_str(enum_node);
        if (enum_str != nullptr) {
          model_enum = atoi(enum_str);
          enum_to_name_[model_enum] = model_name;
        }
      }
      
      cs_tree_node_t *num_params = cs_tree_get_node(model, "NumericalParameters");
      if (num_params != nullptr) {
        
        cs_tree_node_t *wf = cs_tree_get_node(num_params, "WallFunction");
        if (wf != nullptr) {
          const char *default_type = 
            cs_tree_node_get_child_value_str(wf, "DefaultType");
          if (default_type != nullptr) {
            default_wall_function_type_[model_name] = default_type;
          }
          
          cs_tree_node_t *yplus_node = cs_tree_get_node(wf, "YPlusLimit");
          if (yplus_node != nullptr) {
            const char *yplus_str = cs_tree_node_get_value_str(yplus_node);
            if (yplus_str != nullptr) {
              yplus_limit_[model_name] = atof(yplus_str);
            }
          }
        }
        
        const char *diff_model = 
          cs_tree_node_get_child_value_str(num_params, "DiffusionModel");
        if (diff_model != nullptr) {
          diffusion_model_[model_name] = diff_model;
        }
      }
    }
    
    model = cs_tree_node_get_next_of_name(model);
  }
}

void cs_turbulence_rules_manager::parse_validation_rules_check_()
{
  // SECTION 1: Parser <Validations> pour les ranges (MinValue/MaxValue)
  cs_tree_node_t *validations = cs_tree_find_node(rules_tree_, "Validations");
  if (validations != nullptr) {
    cs_tree_node_t *constraint = cs_tree_node_get_child(validations, "Constraint");
    while (constraint != nullptr) {
      const char *target = cs_tree_node_get_tag(constraint, "Target");
      const char *type = cs_tree_node_get_tag(constraint, "Type");
      
      if (target != nullptr && type != nullptr && _strcmp(type, "range")) {
        const char *min_str = cs_tree_node_get_child_value_str(constraint, "MinValue");
        if (min_str != nullptr) {
          validation_min_[target] = atoi(min_str);
        }
        
        const char *max_str = cs_tree_node_get_child_value_str(constraint, "MaxValue");
        if (max_str != nullptr) {
          validation_max_[target] = atoi(max_str);
        }
      }
      
      constraint = cs_tree_node_get_next_of_name(constraint);
    }
  }
  
  // SECTION 2: Parse <ValidationRules> for compatibility rules
  cs_tree_node_t *val_rules = cs_tree_find_node(rules_tree_, "ValidationRules");
  if (val_rules == nullptr) {
    return;
  }
  
  cs_tree_node_t *compat_rules = cs_tree_get_node(val_rules, "CompatibilityRules");
  if (compat_rules != nullptr) {
    cs_tree_node_t *rule = cs_tree_node_get_child(compat_rules, "Rule");
    while (rule != nullptr) {
      const char *param_name = cs_tree_node_get_tag(rule, "Parameter");
      if (param_name != nullptr) {
        cs_tree_node_t *must_be_in = cs_tree_get_node(rule, "MustBeIn");
        if (must_be_in != nullptr) {
          std::vector<int> values;
          cs_tree_node_t *value = cs_tree_node_get_child(must_be_in, "Value");
          while (value != nullptr) {
            const char *val_str = cs_tree_node_get_value_str(value);
            if (val_str != nullptr) {
              values.push_back(atoi(val_str));
            }
            value = cs_tree_node_get_next_of_name(value);
          }
          if (!values.empty()) {
            compatible_itytur_values_[param_name] = values;
          }
        }
      }
      rule = cs_tree_node_get_next_of_name(rule);
    }
  }
  
  cs_tree_node_t *constraints = cs_tree_get_node(val_rules, "Constraints");
  if (constraints != nullptr) {
    cs_tree_node_t *constraint = cs_tree_node_get_child(constraints, "Constraint");
    while (constraint != nullptr) {
      const char *param_name = cs_tree_node_get_tag(constraint, "Parameter");
      if (param_name != nullptr) {
        cs_tree_node_t *for_models = cs_tree_get_node(constraint, "ForModels");
        if (for_models != nullptr) {
          std::vector<std::string> models;
          cs_tree_node_t *model = cs_tree_node_get_child(for_models, "Model");
          while (model != nullptr) {
            const char *model_name = cs_tree_node_get_value_str(model);
            if (model_name != nullptr) {
              models.push_back(model_name);
            }
            model = cs_tree_node_get_next_of_name(model);
          }
          if (!models.empty()) {
            models_requiring_zero_[param_name] = models;
          }
        }
      }
      constraint = cs_tree_node_get_next_of_name(constraint);
    }
  }
  
  cs_tree_node_t *exclusions = cs_tree_get_node(val_rules, "Exclusions");
  if (exclusions != nullptr) {
    cs_tree_node_t *exclusion = cs_tree_node_get_child(exclusions, "Exclusion");
    while (exclusion != nullptr) {
      const char *param_name = cs_tree_node_get_tag(exclusion, "Parameter");
      if (param_name != nullptr) {
        cs_tree_node_t *not_available = cs_tree_get_node(exclusion, "NotAvailableFor");
        if (not_available != nullptr) {
          std::vector<std::string> models;
          cs_tree_node_t *model = cs_tree_node_get_child(not_available, "Model");
          while (model != nullptr) {
            const char *model_name = cs_tree_node_get_value_str(model);
            if (model_name != nullptr) {
              models.push_back(model_name);
            }
            model = cs_tree_node_get_next_of_name(model);
          }
          if (!models.empty()) {
            not_available_for_models_[param_name] = models;
          }
        }
      }
      exclusion = cs_tree_node_get_next_of_name(exclusion);
    }
  }
}

/*============================================================================
 * Constructor & Destructor
 *============================================================================*/

cs_turbulence_rules_manager::cs_turbulence_rules_manager(const char *rules_xml_path)
  : rules_tree_(nullptr)
{
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  
  if (rules_tree_ == nullptr) {
    bft_error(__FILE__, __LINE__, 0,
              "Impossible de charger TurbulenceRules.xml: %s", rules_xml_path);
  }
  
  parse_model_groups_();
  parse_validation_rules_();
  parse_model_requirements_();
  build_model_enum_map_();
  parse_numerical_parameters_();
  parse_validation_rules_check_();
  
  bft_printf("[cs_turbulence_rules_manager] Chargement complet depuis %s\n", 
             rules_xml_path);
  bft_printf("  - %zu models with mapped enum\n", model_enum_map_.size());
  bft_printf("  - %zu LES models detected\n", is_les_model_.size());
  bft_printf("  - %zu Rij models detected\n", is_rij_model_.size());
  bft_printf("  - %zu models with configured wall function\n", 
             default_wall_function_type_.size());
  bft_printf("  - %zu models with configured ypluli\n", 
             yplus_limit_.size());
  bft_printf("  - %zu models with configured diffusion\n", 
             diffusion_model_.size());
  bft_printf("  - %zu parameters with validation ranges\n", 
             validation_min_.size());
  bft_printf("  - %zu compatibility rules loaded\n\n", 
             compatible_itytur_values_.size());
}

cs_turbulence_rules_manager::~cs_turbulence_rules_manager()
{
  if (rules_tree_ != nullptr) {
    cs_tree_node_free(&rules_tree_);
  }
}

/*============================================================================
 * Getters pour cs_gui.cpp
 *============================================================================*/

int cs_turbulence_rules_manager::get_model_enum(const char *model_name) const
{
  if (model_name == nullptr) return CS_TURB_NONE;
  auto it = model_enum_map_.find(model_name);
  if (it != model_enum_map_.end()) {
    return it->second;
  }
  bft_error(__FILE__, __LINE__, 0,
            "Unknown turbulence model: %s", model_name);
  return CS_TURB_NONE;
}

bool cs_turbulence_rules_manager::requires_wall_function(const char *model_name) const
{
  if (model_name == nullptr) return false;
  auto it = requires_wall_function_.find(model_name);
  return (it != requires_wall_function_.end() && it->second);
}

bool cs_turbulence_rules_manager::requires_gravity_terms(const char *model_name) const
{
  if (model_name == nullptr) return false;
  auto it = requires_gravity_terms_.find(model_name);
  return (it != requires_gravity_terms_.end() && it->second);
}

bool cs_turbulence_rules_manager::requires_mixing_length(const char *model_name) const
{
  if (model_name == nullptr) return false;
  auto it = requires_mixing_length_.find(model_name);
  return (it != requires_mixing_length_.end() && it->second);
}

const char* cs_turbulence_rules_manager::get_forced_gravity(const char *model_name) const
{
  if (model_name == nullptr) return nullptr;
  auto it = forced_gravity_.find(model_name);
  if (it != forced_gravity_.end()) {
    return it->second.c_str();
  }
  return nullptr;
}

bool cs_turbulence_rules_manager::is_les_model(const char *model_name) const
{
  if (model_name == nullptr) return false;
  auto it = is_les_model_.find(model_name);
  return (it != is_les_model_.end() && it->second);
}

bool cs_turbulence_rules_manager::is_rij_model(const char *model_name) const
{
  if (model_name == nullptr) return false;
  auto it = is_rij_model_.find(model_name);
  return (it != is_rij_model_.end() && it->second);
}

bool cs_turbulence_rules_manager::supports_coupled_rij(const char *model_name) const
{
  if (model_name == nullptr) return false;
  return (_strcmp(model_name, "Rij-SSG") || _strcmp(model_name, "Rij-EBRSM"));
}

/*============================================================================
 * Getters pour cs_parameters.cpp
 *============================================================================*/

const char* cs_turbulence_rules_manager::get_model_name_from_enum(int model_enum) const
{
  auto it = enum_to_name_.find(model_enum);
  if (it != enum_to_name_.end()) {
    return it->second.c_str();
  }
  return nullptr;
}

const char* cs_turbulence_rules_manager::get_default_wall_function_type(const char *model_name) const
{
  if (model_name == nullptr) return nullptr;
  auto it = default_wall_function_type_.find(model_name);
  if (it != default_wall_function_type_.end()) {
    return it->second.c_str();
  }
  return "2SCALES_LOG";
}

double cs_turbulence_rules_manager::get_yplus_limit(const char *model_name) const
{
  if (model_name == nullptr) return 2.38;
  auto it = yplus_limit_.find(model_name);
  if (it != yplus_limit_.end()) {
    return it->second;
  }
  return 2.38;
}

const char* cs_turbulence_rules_manager::get_diffusion_model(const char *model_name) const
{
  if (model_name == nullptr) return nullptr;
  auto it = diffusion_model_.find(model_name);
  if (it != diffusion_model_.end()) {
    return it->second.c_str();
  }
  return "ISOTROPIC";
}

int cs_turbulence_rules_manager::get_wall_function_enum(const char *wall_function_type) const
{
  if (wall_function_type == nullptr) return CS_WALL_F_2SCALES_LOG;
  
  if (_strcmp(wall_function_type, "DISABLED"))
    return CS_WALL_F_DISABLED;
  else if (_strcmp(wall_function_type, "1SCALE_LOG"))
    return CS_WALL_F_1SCALE_LOG;
  else if (_strcmp(wall_function_type, "2SCALES_LOG"))
    return CS_WALL_F_2SCALES_LOG;
  else if (_strcmp(wall_function_type, "2SCALES_CONTINUOUS"))
    return CS_WALL_F_2SCALES_CONTINUOUS;
  else if (_strcmp(wall_function_type, "SCALABLE_2SCALES_LOG"))
    return CS_WALL_F_SCALABLE_2SCALES_LOG;
  else if (_strcmp(wall_function_type, "2SCALES_VDRIEST"))
    return CS_WALL_F_2SCALES_VDRIEST;
  else if (_strcmp(wall_function_type, "2SCALES_SMOOTH_ROUGH"))
    return CS_WALL_F_2SCALES_SMOOTH_ROUGH;
  else
    return CS_WALL_F_2SCALES_LOG;
}

/*============================================================================
 * Getters pour cs_parameters_check.cpp
 *============================================================================*/

int cs_turbulence_rules_manager::get_validation_min(const char *param_name) const
{
  if (param_name == nullptr) return 0;
  auto it = validation_min_.find(param_name);
  if (it != validation_min_.end()) {
    return it->second;
  }
  return 0;
}

int cs_turbulence_rules_manager::get_validation_max(const char *param_name) const
{
  if (param_name == nullptr) return 0;
  auto it = validation_max_.find(param_name);
  if (it != validation_max_.end()) {
    return it->second;
  }
  return 0;
}

const std::vector<int>* cs_turbulence_rules_manager::get_compatible_itytur_values(const char *param_name) const
{
  if (param_name == nullptr) return nullptr;
  auto it = compatible_itytur_values_.find(param_name);
  if (it != compatible_itytur_values_.end()) {
    return &(it->second);
  }
  return nullptr;
}

bool cs_turbulence_rules_manager::requires_zero_value(const char *param_name, int model_enum) const
{
  if (param_name == nullptr) return false;
  
  auto it = models_requiring_zero_.find(param_name);
  if (it == models_requiring_zero_.end()) return false;
  
  std::string model_const = get_model_constant_name(model_enum);
  if (model_const.empty()) return false;
  
  for (const auto& model : it->second) {
    if (model == model_const) {
      return true;
    }
  }
  return false;
}

bool cs_turbulence_rules_manager::is_available_for_model(const char *param_name, int model_enum) const
{
  if (param_name == nullptr) return true;
  
  auto it = not_available_for_models_.find(param_name);
  if (it == not_available_for_models_.end()) return true;
  
  std::string model_const = get_model_constant_name(model_enum);
  if (model_const.empty()) return true;
  
  for (const auto& model : it->second) {
    if (model == model_const) {
      return false;
    }
  }
  return true;
}

std::string cs_turbulence_rules_manager::get_model_constant_name(int model_enum) const
{
  switch (model_enum) {
    case 0:  return "CS_TURB_NONE";
    case 10: return "CS_TURB_MIXING_LENGTH";
    case 20: return "CS_TURB_K_EPSILON";
    case 21: return "CS_TURB_K_EPSILON_LIN_PROD";
    case 22: return "CS_TURB_K_EPSILON_LS";
    case 23: return "CS_TURB_K_EPSILON_QUAD";
    case 30: return "CS_TURB_RIJ_EPSILON_LRR";
    case 31: return "CS_TURB_RIJ_EPSILON_SSG";
    case 32: return "CS_TURB_RIJ_EPSILON_EBRSM";
    case 40: return "CS_TURB_LES_SMAGO_CONST";
    case 41: return "CS_TURB_LES_SMAGO_DYN";
    case 42: return "CS_TURB_LES_WALE";
    case 50: return "CS_TURB_V2F_PHI";
    case 51: return "CS_TURB_V2F_BL_V2K";
    case 60: return "CS_TURB_K_OMEGA";
    case 70: return "CS_TURB_SPALART_ALLMARAS";
    default: return "";
  }
}

/*============================================================================
 * PARSING DU MODULE OUTPUT
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Parse le module Output depuis le XML
 *----------------------------------------------------------------------------*/

void
cs_turbulence_rules_manager::parse_output_module_()
{
  bft_printf("\n--- Parsing Output module ---\n");
  
  cs_tree_node_t *output_module = cs_tree_find_node(rules_tree_, 
                                                     "Module[@name='Output']");
  if (output_module == nullptr) {
    bft_printf("WARNING: Module 'Output' not found in CodeSaturneRules.xml\n");
    return;
  }
  
  parse_output_field_defaults_();
  parse_output_boundary_variables_();
  parse_output_mappings_();
  
  bft_printf("Output module parsed successfully\n");
}

/*----------------------------------------------------------------------------
 * Parse default field values
 *----------------------------------------------------------------------------*/

void
cs_turbulence_rules_manager::parse_output_field_defaults_()
{
  cs_tree_node_t *output_module = cs_tree_find_node(rules_tree_,
                                                     "Module[@name='Output']");
  if (output_module == nullptr)
    return;
  
  cs_tree_node_t *validation_rules = cs_tree_find_node(output_module,
                                                        "ValidationRules");
  if (validation_rules == nullptr)
    return;
  
  cs_tree_node_t *rule = cs_tree_node_get_child(validation_rules, "Rule");
  while (rule != nullptr) {
    const char *rule_type = cs_tree_node_get_tag(rule, "Type");
    if (rule_type != nullptr && _strcmp(rule_type, "FieldOutputDefaults")) {
      
      const char *target = cs_tree_node_get_tag(rule, "Target");
      if (target != nullptr) {
        cs_output_field_defaults_t defaults;
        
        cs_tree_node_t *listing = cs_tree_find_node(rule, "ListingPrinting");
        const char *listing_val = (listing != nullptr) ? 
                                   cs_tree_node_get_value_str(listing) : "on";
        defaults.listing = _strcmp(listing_val, "on");
        
        cs_tree_node_t *post = cs_tree_find_node(rule, "PostprocessingRecording");
        const char *post_val = (post != nullptr) ? 
                               cs_tree_node_get_value_str(post) : "on";
        defaults.postprocessing = _strcmp(post_val, "on");
        
        cs_tree_node_t *probes = cs_tree_find_node(rule, "ProbesRecording");
        const char *probes_val = (probes != nullptr) ? 
                                 cs_tree_node_get_value_str(probes) : "on";
        defaults.probes = _strcmp(probes_val, "on");
        
        output_field_defaults_[std::string(target)] = defaults;
        
        bft_printf("  Field defaults for '%s': listing=%d, post=%d, probes=%d\n",
                   target, defaults.listing, defaults.postprocessing, defaults.probes);
      }
    }
    
    rule = cs_tree_node_get_next_of_name(rule);
  }
}

/*----------------------------------------------------------------------------
 * Parse les variables surfaciques
 *----------------------------------------------------------------------------*/

void
cs_turbulence_rules_manager::parse_output_boundary_variables_()
{
  cs_tree_node_t *output_module = cs_tree_find_node(rules_tree_,
                                                     "Module[@name='Output']");
  if (output_module == nullptr)
    return;
  
  cs_tree_node_t *validation_rules = cs_tree_find_node(output_module,
                                                        "ValidationRules");
  if (validation_rules == nullptr)
    return;
  
  cs_tree_node_t *rule = cs_tree_node_get_child(validation_rules, "Rule");
  while (rule != nullptr) {
    const char *rule_type = cs_tree_node_get_tag(rule, "Type");
    if (rule_type != nullptr && _strcmp(rule_type, "BoundaryVariableConfig")) {
      
      const char *target = cs_tree_node_get_tag(rule, "Target");
      if (target != nullptr) {
        cs_boundary_var_config_t config;
        
        // Label
        cs_tree_node_t *label_node = cs_tree_find_node(rule, "Label");
        if (label_node != nullptr) {
          const char *label_val = cs_tree_node_get_value_str(label_node);
          config.label = (label_val != nullptr) ? std::string(label_val) : std::string(target);
        } else {
          config.label = std::string(target);
        }
        
        // Default status
        cs_tree_node_t *default_node = cs_tree_find_node(rule, "DefaultStatus");
        if (default_node != nullptr) {
          const char *default_val = cs_tree_node_get_value_str(default_node);
          config.default_status = (default_val != nullptr && _strcmp(default_val, "on"));
        } else {
          config.default_status = false;
        }
        
        // Forbidden models
        cs_tree_node_t *forbidden_node = cs_tree_find_node(rule, "ForbiddenWhen");
        if (forbidden_node != nullptr) {
          // Parser PhysicalModel
          cs_tree_node_t *model = cs_tree_node_get_child(forbidden_node, "PhysicalModel");
          while (model != nullptr) {
            const char *model_name = cs_tree_node_get_value_str(model);
            if (model_name != nullptr) {
              config.forbidden_models.push_back(std::string(model_name));
            }
            model = cs_tree_node_get_next_of_name(model);
          }
          
          // Parser FieldType
          cs_tree_node_t *field = cs_tree_node_get_child(forbidden_node, "FieldType");
          while (field != nullptr) {
            const char *field_name = cs_tree_node_get_value_str(field);
            if (field_name != nullptr) {
              config.forbidden_models.push_back(std::string(field_name));
            }
            field = cs_tree_node_get_next_of_name(field);
          }
        }
        
        // RequiresThermalModel
        cs_tree_node_t *thermal_node = cs_tree_find_node(rule, "RequiresThermalModel");
        config.requires_thermal = (thermal_node != nullptr);
        
        output_boundary_var_config_[std::string(target)] = config;
        
        bft_printf("  Boundary var '%s': label='%s', default=%d\n",
                   target, config.label.c_str(), config.default_status);
      }
    }
    
    rule = cs_tree_node_get_next_of_name(rule);
  }
}

/*----------------------------------------------------------------------------
 * Parse les mappings (TensorComponents, etc.)
 *----------------------------------------------------------------------------*/

void
cs_turbulence_rules_manager::parse_output_mappings_()
{
  cs_tree_node_t *output_module = cs_tree_find_node(rules_tree_,
                                                     "Module[@name='Output']");
  if (output_module == nullptr)
    return;
  
  cs_tree_node_t *mappings = cs_tree_find_node(output_module, "Mappings");
  if (mappings == nullptr)
    return;
  
  // Parser TensorComponents
  cs_tree_node_t *tensor = cs_tree_node_get_child(mappings, "TensorComponents");
  while (tensor != nullptr) {
    const char *tensor_name = cs_tree_node_get_tag(tensor, "Name");
    if (tensor_name != nullptr) {
      std::vector<std::string> component_names;
      
      cs_tree_node_t *comp = cs_tree_node_get_child(tensor, "Component");
      while (comp != nullptr) {
        const char *comp_name = cs_tree_node_get_tag(comp, "Name");
        if (comp_name != nullptr) {
          component_names.push_back(std::string(comp_name));
        }
        comp = cs_tree_node_get_next_of_name(comp);
      }
      
      output_tensor_component_names_[std::string(tensor_name)] = component_names;
      
      bft_printf("  Tensor '%s': %lu components\n", 
                 tensor_name, component_names.size());
    }
    
    tensor = cs_tree_node_get_next_of_name(tensor);
  }
  
  // Parser Defaults
  cs_tree_node_t *defaults_node = cs_tree_find_node(output_module, "Defaults");
  if (defaults_node != nullptr) {
    cs_tree_node_t *def = cs_tree_node_get_child(defaults_node, "Default");
    while (def != nullptr) {
      const char *key = cs_tree_node_get_tag(def, "Key");
      const char *value = cs_tree_node_get_tag(def, "Value");
      
      if (key != nullptr && value != nullptr) {
        output_defaults_[std::string(key)] = std::string(value);
      }
      
      def = cs_tree_node_get_next_of_name(def);
    }
  }
}

/*============================================================================
 * OUTPUT GETTERS (publics)
 *============================================================================*/

cs_output_field_defaults_t
cs_turbulence_rules_manager::get_output_field_defaults(const char *field_type) const
{
  auto it = output_field_defaults_.find(std::string(field_type));
  
  if (it != output_field_defaults_.end()) {
    return it->second;
  }
  
  // Default values if not found
  cs_output_field_defaults_t defaults;
  defaults.listing = true;
  defaults.postprocessing = true;
  defaults.probes = true;
  return defaults;
}

cs_boundary_var_config_t
cs_turbulence_rules_manager::get_output_boundary_var_config(const char *var_name) const
{
  auto it = output_boundary_var_config_.find(std::string(var_name));
  
  if (it != output_boundary_var_config_.end()) {
    return it->second;
  }
  
  // Empty config if not found
  cs_boundary_var_config_t empty_config;
  empty_config.label = std::string(var_name);
  empty_config.default_status = false;
  empty_config.requires_thermal = false;
  return empty_config;
}

bool
cs_turbulence_rules_manager::is_output_boundary_var_forbidden(
    const char *var_name, 
    const char *physics_model) const
{
  cs_boundary_var_config_t config = get_output_boundary_var_config(var_name);
  
  std::string model_str(physics_model);
  for (const auto& forbidden : config.forbidden_models) {
    if (forbidden == model_str) {
      return true;
    }
  }
  
  return false;
}

const std::vector<std::string>*
cs_turbulence_rules_manager::get_output_tensor_component_names(
    const char *tensor_name) const
{
  auto it = output_tensor_component_names_.find(std::string(tensor_name));
  
  if (it != output_tensor_component_names_.end()) {
    return &(it->second);
  }
  
  return nullptr;
}

const char*
cs_turbulence_rules_manager::get_output_default(const char *key) const
{
  auto it = output_defaults_.find(std::string(key));
  
  if (it != output_defaults_.end()) {
    return it->second.c_str();
  }
  
  return nullptr;
}

bool
cs_turbulence_rules_manager::should_disable_yplus(const char *turb_model) const
{
  return _strcmp(turb_model, "off");
}

/*============================================================================
 * FIN DU CODE À AJOUTER
 *============================================================================*/