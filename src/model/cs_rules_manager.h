#ifndef CS_RULES_MANAGER_H
#define CS_RULES_MANAGER_H

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
#include "base/cs_file.h"
#include "base/cs_tree.h"
#include "gui/cs_tree_xml.h"

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_rules_manager {
public:

  cs_rules_manager
  (
    const char *rules_xml_name
  )
  {
    /* Create tree */
    rules_tree_ = cs_tree_node_create("");

    /* Build path */
    char rules_path[1024];
    const char *install_prefix = cs_base_get_pkgdatadir();
    snprintf(rules_path, 1024, "%s/model/%s",
             install_prefix, rules_xml_name);

    if (!cs_file_isreg(rules_path)) {
      // If file does not exist, search in current directory
      snprintf(rules_path, 1024, "%s", rules_xml_name);
    }
    /* Read file */
    cs_tree_xml_read(rules_tree_, rules_path);

    /* Error if file not found... */
    if (rules_tree_ == nullptr)
      bft_error(__FILE__, __LINE__, 0,
                "Cannot load %s: %s\n", rules_xml_name, rules_path);
  }

  ~cs_rules_manager()
  {
    if (rules_tree_ != nullptr)
      cs_tree_node_free(&rules_tree_);
  }

protected:
  /* Members */
  cs_tree_node_t *rules_tree_{nullptr};

  /* Methods */
  static inline int
  _atoi_safe
  (
    const char *s,
    int def = 0
  )
  {
    if (s == nullptr || s[0] == '\0')
      return def;

    return atoi(s);
  }

  static inline double
  _atof_safe
  (
    const char *s,
    double def = 0.0
  )
  {
    if (s == nullptr || s[0] == '\0')
      return def;

    return atof(s);
  }
  static inline bool
  _strcmp
  (
    const char *s1,
    const char *s2
  )
  {
    if (s1 == nullptr || s2 == nullptr)
      return false;

    return (strcmp(s1, s2) == 0);
  }
};

/*----------------------------------------------------------------------------*/

#endif /* CS_RULES_MANAGER_H */
