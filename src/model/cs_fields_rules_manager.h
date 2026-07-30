#ifndef CS_FIELDS_RULES_MANAGER_H
#define CS_FIELDS_RULES_MANAGER_H

/*============================================================================
 * Parse FieldsRules.xml and expose field creation rules
 * for cs_setup.cpp and cs_parameters.cpp.
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

/* Description of a field to creatr */
struct cs_field_creation_rule_t {
  std::string  name;           /* field name (ex: "k") */
  std::string  label;          /* field label (ex: "Turb Kinetic Energy") */
  int          dimension;      /* 1, 3 or 6 */
  std::string  type;           /* "variable", "property", "postprocess",
                                  "internal" */
  std::string  location;       /* "cells", "boundary_faces", "interior_faces" */
  std::string  pointer;        /* ex: "CS_ENUMF_(k)" */
  bool         coupled;        /* CoupledKey=1 */
  bool         no_convection;  /* iconv=0 */
  bool         no_time_term;   /* istat=0 */
  bool         no_diag_shift;  /* idircl=0 */
  bool         no_diffusion;   /* idiff=0 */
  bool         hide_if_steady; /* cache si regime permanent */
  std::string  hide_if_model;  /* cache si modele == X */
  bool         post_process;   /* CS_POST_ON_LOCATION */
  bool         log;            /* log=1 */
  std::string  restart_file;   /* "auxiliary" */
  std::string  exclude_model;  /* exclure pour ce modele */
  int          amr_scheme;     /* schema AMR */
  double       scalar_min;     /* borne min */
  double       scalar_max;     /* borne max */
  int          convection_scheme; /* ischcv */
  int          limiter;        /* limiteur */
  int          slope_test;     /* isstpc */
};

/* Creation rule: condition + list of fields */
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

  /* Return fields to create for a given condition in a given module.
   * Ex: get_fields("Turbulence", "itytur_eq_2") */

  const std::vector<cs_field_creation_rule_t> *
  get_fields(const char  *module_name,
             const char  *condition) const;

  /* Return ALL fields and all rules with this condition */
  std::vector<cs_field_creation_rule_t>
  get_all_fields(const char  *module_name,
                 const char  *condition) const;

  /* Check if a module exists */
  bool has_module(const char  *module_name) const;

  /* Obtain all rules of a module */
  const std::vector<cs_field_creation_entry_t>*
  get_module_rules(const char  *module_name) const;
};

/*============================================================================
 * Singleton accessor
 *============================================================================*/

cs_fields_rules_manager* cs_get_fields_rules_manager();

/*----------------------------------------------------------------------------*/

#endif /* CS_FIELDS_RULES_MANAGER_H */
