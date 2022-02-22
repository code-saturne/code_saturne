#ifndef __CS_GUI_UTIL_H__
#define __CS_GUI_UTIL_H__

/*============================================================================
 * Management of the GUI parameters file: xpath request and utilities
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_tree.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Load the XML file in memory.
 *
 * parameter:
 *   filename <-- XML file containing the parameters
 *
 * returns:
 *   error code (0 in case of success)
 *----------------------------------------------------------------------------*/

int
cs_gui_load_file(const char  *filename);

/*-----------------------------------------------------------------------------
 * Check the xml file version.
 *----------------------------------------------------------------------------*/

void
cs_gui_check_version(void);

/*-----------------------------------------------------------------------------
 * Return the number of characters needed to print an integer number
 *
 * parameters:
 *   num <-- integer number
 *
 * returns:
 *   number of characters required
 *----------------------------------------------------------------------------*/

int
cs_gui_characters_number(int num);

/*-----------------------------------------------------------------------------
 * Compare two strings.
 *
 * parameters:
 *   s1 <-- first string
 *   s2 <-- second string
 *
 * returns:
 *   1 if the strings are equal, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_gui_strcmp(const char  *s1,
              const char  *s2);

/*-----------------------------------------------------------------------------
 * Test if 2 real values are equal (avoiding compiler warnings)
 *
 * parameters:
 *   v1 <-- first value to compare
 *   v2 <-- second value to compare
 *
 * returns:
 *   1 if values are equal, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_gui_is_equal_real(cs_real_t v1,
                     cs_real_t v2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returned.
 *
 * \param[in]       node    node whose value is queried
 * \param[in, out]  value   queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_int(cs_tree_node_t  *node,
                    int             *value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an real value based on a tree node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returned.
 *
 * \param[in]       node    node whose value is queried
 * \param[in, out]  value   queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_real(cs_tree_node_t  *node,
                     cs_real_t       *value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer-valued status value based on a node's status tag.
 *
 * The status is defined in a string-valued child (tag) node. If no such
 * tag is present, the initial status is unchanged.
 *
 * \param[in]       node    node whose status is queried
 * \param[in, out]  status  status value (0 or 1)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_status_int(cs_tree_node_t  *node,
                           int             *status);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an bool-valued status value based on a node's status tag.
 *
 * The status is defined in a string-valued child (tag) node. If no such
 * tag is present, the initial status is unchanged.
 *
 * \param[in]       node    node whose status is queried
 * \param[in, out]  status  status value (0 or 1)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_status_bool(cs_tree_node_t  *node,
                            bool            *status);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string value associated with a "tag" child node and
 *         whose presence should be guaranteed.
 *
 * If the matching child node is not present, an error is produced
 *
 * \param[in]  node      tree node which should have a "tag" child
 * \param[in]  tag_name  name of tag child node
 *
 * \return  pointer to matching child string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_gui_node_get_tag(cs_tree_node_t  *node,
                    const char      *tag_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree's child node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returned.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_int(cs_tree_node_t  *node,
                          const char      *child_name,
                          int             *value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree's child node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returned.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_real(cs_tree_node_t  *node,
                           const char      *child_name,
                           cs_real_t       *value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer-valued status value based on a node child's
 *         status tag.
 *
 * The status is defined in a string-valued child (tag) node. If no such
 * child and tag is present, the initial status is unchanged.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  status      queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_status_int(cs_tree_node_t  *node,
                                 const char      *child_name,
                                 int             *status);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update a bool-valued status value based on a node child's
 *         status tag.
 *
 * The status is defined in a string-valued child (tag) node. If no such
 * child and tag is present, the initial status is unchanged.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  status      queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_status_bool(cs_tree_node_t  *node,
                                  const char      *child_name,
                                  bool            *status);

/*-----------------------------------------------------------------------------
 * Add timing increment to global MEI time counter.
 *
 * parameters:
 *   t <-- timing increment to add
 *----------------------------------------------------------------------------*/

void
cs_gui_add_mei_time(double t);

/*-----------------------------------------------------------------------------
 * Get cumulative global MEI time counter.
 *
 * returns:
 *   cumulative global MEI time counter
 *----------------------------------------------------------------------------*/

double
cs_gui_get_mei_times(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_UTIL_H__ */
