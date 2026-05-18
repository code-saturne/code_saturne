/*============================================================================
 * cs_elec_rules_manager.h
 *
 * Class to parse and manage ElectricalRules.xml
 *============================================================================*/
#ifndef CS_ELEC_RULES_MANAGER_H
#define CS_ELEC_RULES_MANAGER_H

#include "base/cs_defs.h"
#include "base/cs_tree.h"
#include <map>
#include <string>
#include <vector>

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Joule model entry */
struct cs_elec_joule_model_t {
  std::string name;         /* "AC/DC", "three-phase"... */
  int         ieljou;       /* -1, 1, 2, 3, 4 */
  std::string description;
};

/* Electric arc model entry */
struct cs_elec_arc_model_t {
  std::string name;         /* "off", "on" */
  int         ielarc;       /* -1, 2 */
  std::string description;
};

/* Physical property entry */
struct cs_elec_property_t {
  std::string name;
  std::string label;
  std::string unit;
  std::string model;  /* "arc", "joule_transformer", or "" for all */
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_elec_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_elec_joule_model_t> joule_models_;
  std::map<std::string, cs_elec_arc_model_t>   arc_models_;
  std::map<std::string, std::string>            defaults_;
  std::vector<cs_elec_property_t>               properties_;
  std::vector<int>                              allowed_ieljou_;
  std::vector<int>                              allowed_ielarc_;

  void parse_rules_();

public:
  cs_elec_rules_manager(const char *rules_xml_path);
  ~cs_elec_rules_manager();

  /* Get ieljou value from joule model name */
  int get_ieljou(const std::string &name) const;

  /* Get ielarc value from arc model name */
  int get_ielarc(const std::string &name) const;

  /* Check if ieljou value is valid */
  bool is_valid_ieljou(int ieljou) const;

  /* Check if ielarc value is valid */
  bool is_valid_ielarc(int ielarc) const;

  /* Get default value */
  std::string get_default(const std::string &key) const;
  double      get_default_double(const std::string &key) const;

  /* Get properties for a given model */
  std::vector<cs_elec_property_t>
  get_properties(const std::string &model) const;
};

cs_elec_rules_manager *cs_get_elec_rules_manager(void);

#endif /* CS_ELEC_RULES_MANAGER_H */
