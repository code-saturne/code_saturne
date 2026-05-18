/*============================================================================
 * cs_lagr_rules_manager.h
 *
 * Class to parse and manage LagrangianRules.xml
 *============================================================================*/
#ifndef CS_LAGR_RULES_MANAGER_H
#define CS_LAGR_RULES_MANAGER_H

#include "base/cs_defs.h"
#include "base/cs_tree.h"
#include <map>
#include <string>

/*============================================================================
 * Structure definitions
 *============================================================================*/

struct cs_lagr_coupling_mode_t {
  std::string name;         /* "off", "one_way", "two_way", "frozen" */
  std::string enum_name;    /* "CS_LAGR_OFF", "CS_LAGR_ONEWAY_COUPLING"... */
  int         int_value;    /* 0, 1, 2, 3 */
  std::string description;
};

struct cs_lagr_phys_model_t {
  std::string name;         /* "off", "thermal", "coal", "ctwr" */
  std::string enum_name;    /* "CS_LAGR_PHYS_OFF"... */
  int         int_value;    /* 0, 1, 2, 3 */
  std::string description;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_lagr_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_lagr_coupling_mode_t> coupling_modes_;
  std::map<std::string, cs_lagr_phys_model_t>    phys_models_;
  std::map<std::string, std::string>              defaults_;

  void parse_rules_();

public:
  cs_lagr_rules_manager(const char *rules_xml_path);
  ~cs_lagr_rules_manager();

  /* Get coupling mode int value from name */
  int get_coupling_mode(const std::string &name) const;

  /* Get physical model int value from name */
  int get_phys_model(const std::string &name) const;

  /* Get default value */
  std::string get_default(const std::string &key) const;
  double      get_default_double(const std::string &key) const;
  int         get_default_int(const std::string &key) const;

  /* Check if model requires restart */
  bool requires_restart(const std::string &coupling_mode) const;
};

cs_lagr_rules_manager *cs_get_lagr_rules_manager(void);

#endif /* CS_LAGR_RULES_MANAGER_H */
