#ifndef CS_TURBULENT_FLUX_RULES_MANAGER_H
#define CS_TURBULENT_FLUX_RULES_MANAGER_H

/*============================================================================
 * Parse TurbulentFluxRules.xml and expose rules for creation of
 * turbulente flux fields.
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

  /* Obtain rule fro a given model (ex: "DFM") */
  const cs_turb_flux_rule_t *
  get_rule_by_model(const std::string &model) const;

  /* Obtain rule for a numerical value (ex: 30 for DFM) */
  const cs_turb_flux_rule_t *
  get_rule_by_value(int numeric_value) const;

  /* Obtain name of a field by substituting {scalar} */
  static std::string
  resolve_name(const std::string  &name_template,
               const std::string  &scalar_name);

  /* Obtain numerical value for a GUI model (ex: "DFM" -> 30) */
  int
  get_numeric_value(const std::string  &model_name) const;

  /* Check if a model created additional fields */
  bool
  creates_fields(const std::string  &model_name) const;

  /* Check if a model created an alpha field (EB-*) */
  bool
  creates_alpha_field(const std::string  &model_name) const;
};

/*============================================================================
 * Singleton access
 *============================================================================*/

cs_turbulent_flux_rules_manager *
cs_get_turbulent_flux_rules_manager(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENT_FLUX_RULES_MANAGER_H */
