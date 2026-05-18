/*============================================================================
 * cs_fields_rules_manager.h
 *
 * Manager pour parser FieldsRules.xml et exposer les regles de creation
 * de champs a cs_setup.cpp et cs_parameters.cpp.
 *
 *
 *============================================================================*/

#ifndef CS_FIELDS_RULES_MANAGER_H
#define CS_FIELDS_RULES_MANAGER_H

#include "base/cs_defs.h"
#include "base/cs_tree.h"
#include <map>
#include <string>
#include <vector>

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Description d'un champ a creer */
struct cs_field_creation_rule_t {
  std::string name;           /* nom du champ ex: "k" */
  std::string label;          /* label ex: "Turb Kinetic Energy" */
  int         dimension;      /* 1, 3 ou 6 */
  std::string type;           /* "variable", "property", "postprocess", "internal" */
  std::string location;       /* "cells", "boundary_faces", "interior_faces" */
  std::string pointer;        /* ex: "CS_ENUMF_(k)" */
  bool        coupled;        /* CoupledKey=1 */
  bool        no_convection;  /* iconv=0 */
  bool        no_time_term;   /* istat=0 */
  bool        no_diag_shift;  /* idircl=0 */
  bool        no_diffusion;   /* idiff=0 */
  bool        hide_if_steady; /* cache si regime permanent */
  std::string hide_if_model;  /* cache si modele == X */
  bool        post_process;   /* CS_POST_ON_LOCATION */
  bool        log;            /* log=1 */
  std::string restart_file;   /* "auxiliary" */
  std::string exclude_model;  /* exclure pour ce modele */
  int         amr_scheme;     /* schema AMR */
  double      scalar_min;     /* borne min */
  double      scalar_max;     /* borne max */
  int         convection_scheme; /* ischcv */
  int         limiter;        /* limiteur */
  int         slope_test;     /* isstpc */
};

/* Regle de creation : condition + liste de champs */
struct cs_field_creation_entry_t {
  std::string condition;
  std::string description;
  std::vector<cs_field_creation_rule_t> fields;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_fields_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  /* Champs par module et condition */
  /* module_name -> liste de regles */
  std::map<std::string, std::vector<cs_field_creation_entry_t>> modules_;

  /* Parsing */
  void parse_modules_();
  cs_field_creation_rule_t parse_field_(cs_tree_node_t *field_node);

  static inline bool _strcmp(const char *s1, const char *s2);
  static inline int  _atoi_safe(const char *s, int def = 0);
  static inline double _atof_safe(const char *s, double def = 0.0);

public:
  cs_fields_rules_manager(const char *rules_xml_path);
  ~cs_fields_rules_manager();

  /* =========================================================
   * GETTERS POUR cs_setup.cpp et cs_parameters.cpp
   * ========================================================= */

  /* Obtenir tous les champs a creer pour une condition donnee
   * dans un module donne.
   * Ex: get_fields("Turbulence", "itytur_eq_2")
   */
  const std::vector<cs_field_creation_rule_t>*
  get_fields(const char *module_name, const char *condition) const;

  /* Retourne TOUS les champs de toutes les Rules avec cette condition */
  std::vector<cs_field_creation_rule_t>
  get_all_fields(const char *module_name, const char *condition) const;

  /* Verifier si un module existe */
  bool has_module(const char *module_name) const;

  /* Obtenir toutes les regles d'un module */
  const std::vector<cs_field_creation_entry_t>*
  get_module_rules(const char *module_name) const;

};

/*============================================================================
 * Singleton accessor
 *============================================================================*/

cs_fields_rules_manager* cs_get_fields_rules_manager();

#endif /* CS_FIELDS_RULES_MANAGER_H */
