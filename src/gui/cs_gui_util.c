/*============================================================================
 * Management of the GUI parameters file: utility functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_parameters.h"
#include "cs_tree.h"
#include "cs_tree_xml.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_util.h"


/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define XML_READER_VERSION 2.0

/*=============================================================================
 * Global variables
 *============================================================================*/

double _cs_gui_mei_time = 0.;

static bool _setup_read = false;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate if an XML file has been loaded.
 *
 * \return  1 if an XML file has been loaded, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_gui_file_is_loaded(void)
{
  int retval = 0;
  if (_setup_read)
    retval = 1;
  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Load the XML file in memory.
 *
 * \param[in]  filename  XML file containing the parameters
 *
 * \return  error code (0 in case of success)
 */
/*----------------------------------------------------------------------------*/

int
cs_gui_load_file(const char  *filename)
{
  int argerr = 0;

  if (cs_glob_tree == NULL)
    cs_glob_tree = cs_tree_node_create(NULL);

  cs_tree_xml_read(cs_glob_tree, filename);

  _setup_read = true;
  return argerr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the XML file version.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_check_version(void)
{
  double version_sat = XML_READER_VERSION;
  double major, maj_sat;

  cs_tree_node_t *tn = cs_glob_tree;
  tn = cs_tree_get_node(tn, "Code_Saturne_GUI");
  if (tn == NULL)
    tn = cs_tree_get_node(cs_glob_tree, "NEPTUNE_CFD_GUI");

  const char *version = cs_tree_node_get_tag(tn, "version");

  double version_number = (version != NULL) ? atof(version) : 0.0;

  double minus = modf(version_number, &major);
  double min_sat = modf(version_sat, &maj_sat);

  if (!cs_gui_is_equal_real(major, maj_sat))
    bft_error(__FILE__, __LINE__, 0,
              _("========================================================\n"
                "   ** Invalid version of the XML file\n"
                "      -------------------------------------- \n"
                "      XML file version: %.1f  \n"
                "      XML reader version: %.1f \n"
                "========================================================\n"),
              version_number, version_sat);

  if (!cs_gui_is_equal_real(minus, min_sat)) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("========================================================\n"
                 "   ** Unexpected version XML file version\n"
                 "      -----------------------------------\n"
                 "      XML file version: %.1f  \n"
                 "      XML reader version: %.1f \n"
                 "\n"
                 "      It is recommended to rebuild a new XML file.\n"
                 "========================================================\n"),
               version_number, version_sat);
  }
}

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
cs_gui_characters_number(int num)
{
  int i      = 1;
  int number = 0;

  assert(num >= 0);

  if (num == 0)
    number ++;
  else
    for (i=1; i <= num; i *= 10)
      number++;

  return number;
}

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
              const char  *s2)
{
  if (s1 == NULL || s2 == NULL) return 0;
  if (strlen(s1) != strlen(s2)) return 0;
  if (!strncmp(s1, s2, strlen(s1))) return 1;
  return 0;
}

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
                     cs_real_t v2)
{
  size_t i;
  int retval = 1;

  const unsigned char *_v1 = (const unsigned char *)(&v1);
  const unsigned char *_v2 = (const unsigned char *)(&v2);

  for (i = 0; i < sizeof(cs_real_t); i++) {
    if (_v1[i] != _v2[i])
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returne.
 *
 * \param[in]       node    node whose value is queried
 * \param[in, out]  value   queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_int(cs_tree_node_t  *node,
                    int             *value)
{
  if (node != NULL) {

    const int *v_i = cs_tree_node_get_values_int(node);

    if (node->size != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Expected 1 value for node %s, not %d"),
                node->name, node->size);

    if (v_i != NULL)
      *value = v_i[0];
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Missing values for node %s"), node->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an real value based on a tree node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returne.
 *
 * \param[in]       node    node whose value is queried
 * \param[in, out]  value   queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_real(cs_tree_node_t  *node,
                     cs_real_t       *value)
{
  if (node != NULL) {

    const cs_real_t *v_r = cs_tree_node_get_values_real(node);

    if (node->size != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Expected 1 value for node %s, not %d"),
                node->name, node->size);

    if (v_r != NULL)
      *value = v_r[0];
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Missing values for node %s"), node->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer-valued status value based on a node's status tag.
 *
 * The status is defined in a string-valued child (tag) node. If no such
 * tag is present, the initial status is unchanged.
 *
 * \param[in]       tn      node whose status is queried
 * \param[in, out]  status  status value (0 or 1)
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_status_int(cs_tree_node_t  *tn,
                           int             *status)
{
  const char  *value = cs_tree_node_get_tag(tn, "status");

  if (cs_gui_strcmp(value, "on"))
    *status = 1;
  else if (cs_gui_strcmp(value, "off"))
    *status = 0;
  else if (value != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid status value: %s"), value);
}

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
                            bool            *status)
{
  const char  *value = cs_tree_node_get_tag(node, "status");

  if (cs_gui_strcmp(value, "on"))
    *status = true;
  else if (cs_gui_strcmp(value, "off"))
    *status = false;
  else if (value != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid status value: %s"), value);
}

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
                    const char      *tag_name)
{
  const char *name = cs_tree_node_get_tag(node, tag_name);

  if (name == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Incorrect setup tree definition for the following node:\n"));
    cs_tree_dump(CS_LOG_DEFAULT, 2, node);
    bft_error(__FILE__, __LINE__, 0,
              _("Missing child (tag) node: %s"), tag_name);
  }

  return name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree's child node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returne.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_int(cs_tree_node_t  *node,
                          const char      *child_name,
                          int             *value)
{
  cs_tree_node_t *tn_c = cs_tree_node_get_child(node, child_name);

  if (tn_c != NULL) {
    const int *v_i = cs_tree_node_get_values_int(tn_c);

    if (tn_c->size != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Expected 1 value for node %s, not %d"),
                tn_c->name, tn_c->size);


    if (v_i != NULL)
      *value = v_i[0];
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Missing values for node %s"), tn_c->name);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update an integer value based on a tree's child node
 *
 * If no node is present, the initial value is unchanged.
 * If the node is present but the value missing, an error is returne.
 *
 * \param[in]       node        node whose value is queried
 * \param[in]       child_name  name of child node
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_real(cs_tree_node_t  *node,
                           const char      *child_name,
                           cs_real_t       *value)
{
  cs_tree_node_t *tn_c = cs_tree_node_get_child(node, child_name);

  if (tn_c != NULL) {
    const cs_real_t *v_r = cs_tree_node_get_values_real(tn_c);

    if (tn_c->size != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Expected 1 value for node %s, not %d"),
                tn_c->name, tn_c->size);


    if (v_r != NULL)
      *value = v_r[0];
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Missing values for node %s"), tn_c->name);
  }
}

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
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_status_int(cs_tree_node_t  *tn,
                                 const char      *child_name,
                                 int             *status)
{
  cs_tree_node_t *tn_c = cs_tree_node_get_child(tn, child_name);

  const char  *value = cs_tree_node_get_tag(tn_c, "status");

  if (value != NULL) {
    if (! strcmp(value, "on"))
      *status = 1;
    else if (! strcmp(value, "off"))
      *status = 0;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid status value: %s"), value);
  }
}

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
 * \param[in, out]  value       queried value
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_node_get_child_status_bool(cs_tree_node_t  *node,
                                  const char      *child_name,
                                  bool            *status)
{
  cs_tree_node_t *tn_c = cs_tree_node_get_child(node, child_name);

  const char  *value = cs_tree_node_get_tag(tn_c, "status");

  if (value != NULL) {
    if (! strcmp(value, "on"))
      *status = true;
    else if (! strcmp(value, "off"))
      *status = false;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid status value: %s"), value);
  }
}

/*-----------------------------------------------------------------------------
 * Add timing increment to global MEI time counter.
 *
 * parameters:
 *   t <-- timing increment to add
 *----------------------------------------------------------------------------*/

void
cs_gui_add_mei_time(double t)
{
  _cs_gui_mei_time += t;
}

/*-----------------------------------------------------------------------------
 * Get cumulative global MEI time counter.
 *
 * returns:
 *   cumulative global MEI time counter
 *----------------------------------------------------------------------------*/

double
cs_gui_get_mei_times(void)
{
  return _cs_gui_mei_time;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
