/*============================================================================
 * cs_turbulent_flux_rules_manager.h
 *
 * Manager pour parser TurbulentFluxRules.xml et exposer les regles de
 * creation des champs de flux turbulent.
 *
 *
 *============================================================================*/
#ifndef CS_TURBULENT_FLUX_RULES_MANAGER_H
#define CS_TURBULENT_FLUX_RULES_MANAGER_H

#include "base/cs_defs.h"
#include "base/cs_tree.h"

#include <map>
#include <string>
#include <vector>

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Description d'un champ de flux turbulent a creer */
struct cs_turb_flux_field_t {
  std::string name_template;   /* ex: "{scalar}_turbulent_flux" */
  std::string label_template;  /* ex: "{scalar} Turbulent Flux" */
  int         dimension;       /* 1 ou 3 */
  std::string type;            /* "variable" ou "property" */
  bool        coupled;         /* Coupled=1 pour DFM */
  bool        clipping;        /* Clipping=1 pour DFM */
  bool        post_process;    /* PostProcess=1 */
  bool        log;             /* Log=1 */
};

/* Regle pour un modele de flux turbulent */
struct cs_turb_flux_rule_t {
  std::string model_name;      /* "SGDH", "GGDH", "DFM", etc. */
  int         numeric_value;   /* 0, 10, 11, 20, 21, 30, 31 */
  std::string description;
  std::vector<cs_turb_flux_field_t> fields;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_turbulent_flux_rules_manager {
private:
  cs_tree_node_t *rules_tree_;
  std::map<std::string, cs_turb_flux_rule_t> rules_by_model_;
  std::map<int, cs_turb_flux_rule_t>         rules_by_value_;

  void parse_rules_();
  cs_turb_flux_field_t parse_field_(cs_tree_node_t *field_node);

  static inline int _atoi_safe(const char *s, int def = 0);

public:
  cs_turbulent_flux_rules_manager(const char *rules_xml_path);
  ~cs_turbulent_flux_rules_manager();

  /* Obtenir la regle pour un modele donne (ex: "DFM") */
  const cs_turb_flux_rule_t *get_rule_by_model(const std::string &model) const;

  /* Obtenir la regle pour une valeur numerique (ex: 30 pour DFM) */
  const cs_turb_flux_rule_t *get_rule_by_value(int numeric_value) const;

  /* Obtenir le nom du champ en substituant {scalar} */
  static std::string resolve_name(const std::string &name_template,
                                  const std::string &scalar_name);

  /* Obtenir la valeur numerique pour un modele GUI (ex: "DFM" -> 30) */
  int get_numeric_value(const std::string &model_name) const;

  /* Verifier si un modele cree des champs supplementaires */
  bool creates_fields(const std::string &model_name) const;

  /* Verifier si un modele cree un champ alpha (EB-*) */
  bool creates_alpha_field(const std::string &model_name) const;
};

/*============================================================================
 * Singleton access
 *============================================================================*/

cs_turbulent_flux_rules_manager *cs_get_turbulent_flux_rules_manager(void);

#endif /* CS_TURBULENT_FLUX_RULES_MANAGER_H */
