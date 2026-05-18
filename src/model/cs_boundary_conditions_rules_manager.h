/*============================================================================
 * cs_boundary_conditions_rules_manager.h
 *
 * Class to parse and manage BoundaryConditionsRules.xml
 *============================================================================*/
#ifndef CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H
#define CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H

#include "base/cs_defs.h"
#include "base/cs_tree.h"
#include <map>
#include <string>

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* BC code mapping: boundary_type -> icodcl */
struct cs_bc_code_entry_t {
  std::string boundary_type;  /* e.g. "CS_BOUNDARY_WALL" */
  int         icodcl;         /* legacy solver code: 1, 5, 6, 13... */
  std::string label;          /* human-readable label */
};

/* Legacy type mapping: nature -> legacy code */
struct cs_bc_legacy_entry_t {
  std::string nature;       /* "wall", "symmetry", "default" */
  std::string legacy_code;  /* "CS_SMOOTHWALL", "CS_SYMMETRY" */
  int         value;        /* numeric value */
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_boundary_conditions_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_bc_code_entry_t>  bc_code_map_;
  std::map<std::string, cs_bc_legacy_entry_t> legacy_map_;
  int default_icodcl_;

  void parse_rules_();

public:
  cs_boundary_conditions_rules_manager(const char *rules_xml_path);
  ~cs_boundary_conditions_rules_manager();

  /* Get icodcl for a given boundary type flag name */
  int get_icodcl(const std::string &boundary_type) const;

  /* Get legacy type code for a given nature */
  int get_legacy_type(const std::string &nature) const;
};

cs_boundary_conditions_rules_manager *
cs_get_boundary_conditions_rules_manager(void);

#endif /* CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H */
