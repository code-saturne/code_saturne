/*============================================================================
 * cs_fields_rules_manager.cpp
 *
 * Implementation du manager pour FieldsRules.xml.
 *
 * 
 *============================================================================*/

#include "cs_fields_rules_manager.h"
#include "base/cs_base.h"
#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "gui/cs_tree_xml.h"
#include <cstring>
#include <cstdlib>

/*============================================================================
 * Static helpers
 *============================================================================*/

inline bool
cs_fields_rules_manager::_strcmp(const char *s1, const char *s2)
{
  if (s1 == nullptr || s2 == nullptr) return false;
  return (std::strcmp(s1, s2) == 0);
}

inline int
cs_fields_rules_manager::_atoi_safe(const char *s, int def)
{
  if (s == nullptr) return def;
  return atoi(s);
}

inline double
cs_fields_rules_manager::_atof_safe(const char *s, double def)
{
  if (s == nullptr) return def;
  return atof(s);
}

/*============================================================================
 * Constructor
 *============================================================================*/

cs_fields_rules_manager::cs_fields_rules_manager(const char *rules_xml_path)
  : rules_tree_(nullptr)
{
  bft_printf("\n=== Initializing cs_fields_rules_manager ===\n");
  bft_printf("Loading: %s\n", rules_xml_path);

  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);

  if (rules_tree_ == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Error: Could not load FieldsRules.xml from: %s\n",
              rules_xml_path);

  parse_modules_();

  bft_printf("cs_fields_rules_manager initialized: %lu modules\n\n",
             modules_.size());
}

/*============================================================================
 * Destructor
 *============================================================================*/

cs_fields_rules_manager::~cs_fields_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Parse a single Field node
 *============================================================================*/

cs_field_creation_rule_t
cs_fields_rules_manager::parse_field_(cs_tree_node_t *field_node)
{
  cs_field_creation_rule_t rule;

  /* Attributs de base */
  auto get_attr = [&](const char *attr) -> const char * {
    return cs_tree_node_get_tag(field_node, attr);
  };

  const char *name = get_attr("Name");
  if (name) rule.name = std::string(name);

  const char *label = get_attr("Label");
  if (label) rule.label = std::string(label);

  rule.dimension = _atoi_safe(get_attr("Dimension"), 1);

  const char *type = get_attr("Type");
  rule.type = type ? std::string(type) : "variable";

  const char *loc = get_attr("Location");
  rule.location = loc ? std::string(loc) : "cells";

  const char *pointer = get_attr("Pointer");
  if (pointer) rule.pointer = std::string(pointer);

  /* Boolean flags */
  rule.coupled        = _strcmp(get_attr("CoupledKey"),    "1");
  rule.no_convection  = _strcmp(get_attr("NoConvection"),  "1");
  rule.no_time_term   = _strcmp(get_attr("NoTimeTerm"),    "1");
  rule.no_diag_shift  = _strcmp(get_attr("NoDiagShift"),   "1");
  rule.no_diffusion   = _strcmp(get_attr("NoDiffusion"),   "1");
  rule.hide_if_steady = _strcmp(get_attr("HideIfSteady"),  "1");
  rule.post_process   = _strcmp(get_attr("PostProcess"),   "1");
  rule.log            = _strcmp(get_attr("Log"),           "1");

  /* Attributs optionnels */
  const char *him = get_attr("HideIfModel");
  if (him) rule.hide_if_model = std::string(him);

  const char *rf = get_attr("RestartFile");
  if (rf) rule.restart_file = std::string(rf);

  const char *em = get_attr("ExcludeModel");
  if (em) rule.exclude_model = std::string(em);

  rule.amr_scheme        = _atoi_safe(get_attr("AMRScheme"), 0);
  rule.scalar_min        = _atof_safe(get_attr("ScalarMin"), 0.0);
  rule.scalar_max        = _atof_safe(get_attr("ScalarMax"), 0.0);
  rule.convection_scheme = _atoi_safe(get_attr("ConvectionScheme"), 0);
  rule.limiter           = _atoi_safe(get_attr("Limiter"), -1);
  rule.slope_test        = _atoi_safe(get_attr("SlopeTest"), 0);

  return rule;
}

/*============================================================================
 * Parse all modules
 *============================================================================*/

void
cs_fields_rules_manager::parse_modules_()
{
  cs_tree_node_t *fcr = cs_tree_find_node(rules_tree_, "FieldCreationRules");
  if (fcr == nullptr) {
    bft_printf("WARNING: No FieldCreationRules found in FieldsRules.xml\n");
    return;
  }

  for (cs_tree_node_t *module = cs_tree_node_get_child(fcr, "Module");
       module != nullptr;
       module = cs_tree_node_get_next_of_name(module)) {

    const char *module_name = cs_tree_node_get_tag(module, "Name");
    if (module_name == nullptr) continue;

    std::vector<cs_field_creation_entry_t> entries;

    for (cs_tree_node_t *rule = cs_tree_node_get_child(module, "Rule");
         rule != nullptr;
         rule = cs_tree_node_get_next_of_name(rule)) {

      const char *condition = cs_tree_node_get_tag(rule, "Condition");
      cs_field_creation_entry_t entry;
      entry.condition = condition ? std::string(condition) : "always";

      /* Description */
      cs_tree_node_t *desc = cs_tree_find_node(rule, "Description");
      if (desc != nullptr) {
        const char *d = cs_tree_node_get_value_str(desc);
        if (d) entry.description = std::string(d);
      }

      /* Champs */
      for (cs_tree_node_t *field = cs_tree_node_get_child(rule, "Field");
           field != nullptr;
           field = cs_tree_node_get_next_of_name(field)) {
        entry.fields.push_back(parse_field_(field));
      }

      entries.push_back(entry);
    }

    modules_[std::string(module_name)] = entries;
    bft_printf("  Module '%s': %lu rules\n", module_name, entries.size());
  }
}

/*============================================================================
 * GETTERS
 *============================================================================*/

const std::vector<cs_field_creation_rule_t>*
cs_fields_rules_manager::get_fields(const char *module_name,
                                    const char *condition) const
{
  if (module_name == nullptr || condition == nullptr) return nullptr;

  auto it = modules_.find(std::string(module_name));
  if (it == modules_.end()) return nullptr;

  /* Chercher la premiere rule matching - compatibilite ancienne API */
  for (const auto &entry : it->second) {
    if (entry.condition == std::string(condition))
      return &entry.fields;
  }
  return nullptr;
}

std::vector<cs_field_creation_rule_t>
cs_fields_rules_manager::get_all_fields(const char *module_name,
                                        const char *condition) const
{
  std::vector<cs_field_creation_rule_t> result;
  if (module_name == nullptr || condition == nullptr) return result;

  auto it = modules_.find(std::string(module_name));
  if (it == modules_.end()) return result;

  /* Concatener les champs de TOUTES les Rules avec cette condition */
  for (const auto &entry : it->second) {
    if (entry.condition == std::string(condition)) {
      for (const auto &field : entry.fields)
        result.push_back(field);
    }
  }
  return result;
}

bool
cs_fields_rules_manager::has_module(const char *module_name) const
{
  if (module_name == nullptr) return false;
  return modules_.find(std::string(module_name)) != modules_.end();
}

const std::vector<cs_field_creation_entry_t>*
cs_fields_rules_manager::get_module_rules(const char *module_name) const
{
  if (module_name == nullptr) return nullptr;
  auto it = modules_.find(std::string(module_name));
  if (it == modules_.end()) return nullptr;
  return &it->second;
}

/*============================================================================
 * Singleton
 *============================================================================*/

static cs_fields_rules_manager *g_fields_manager = nullptr;

cs_fields_rules_manager*
cs_get_fields_rules_manager()
{
  if (g_fields_manager == nullptr) {
    char rules_path[1024];
    const char *datadir = cs_base_get_pkgdatadir();
    snprintf(rules_path, 1024, "%s/model/FieldsRules.xml", datadir);

    FILE *test_file = fopen(rules_path, "r");
    if (test_file != nullptr) {
      fclose(test_file);
    }
    else {
      bft_printf("WARNING: FieldsRules.xml not found in %s\n", rules_path);
      snprintf(rules_path, 1024, "FieldsRules.xml");
    }

    g_fields_manager = new cs_fields_rules_manager(rules_path);
  }
  return g_fields_manager;
}
