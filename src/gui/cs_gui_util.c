/*============================================================================
 * Management of the GUI parameters file: utility functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * libxml2 library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_LIBXML2)

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#endif

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

#if defined(HAVE_LIBXML2)
xmlDocPtr docxml            = NULL;   /* Pointer on the XML document  */
xmlXPathContextPtr xpathCtx = NULL;   /* Pointer on the XPath Context */
xmlNodePtr node             = NULL;   /* Pointer on the root node     */
const char *xmlRootName     = NULL;   /* Name of the root node        */
#endif

double _cs_gui_mei_time = 0.;

static bool _setup_read = false;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Indicate if an XML file has been loaded
 *
 * returns:
 *    1 if an XML file has been loaded, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_gui_file_is_loaded(void)
{
  int retval = 0;
  if (_setup_read)
    retval = 1;
  return retval;
}

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
cs_gui_load_file(const char  *filename)
{
  int argerr = 0;

  if (cs_glob_tree == NULL)
    cs_glob_tree = cs_tree_node_create(NULL);

  cs_tree_xml_read(cs_glob_tree, filename);

#if defined(HAVE_LIBXML2)

  int file_descriptor = 0;

  assert(filename);

  /* Verification of the existence of the file while opening it */

  file_descriptor = open(filename, O_RDONLY);

  if (file_descriptor ==  -1) {

    bft_error(__FILE__, __LINE__, 0,
              _("Unable to open the file: %s\n"), filename);
    argerr = 2;
    return argerr;

  }
  else {

    /* If file exists, close it. It will be reopen by xmlParseFile */
    close(file_descriptor);

  }

  /* libxml initialization */
  xmlInitParser();
  LIBXML_TEST_VERSION

  /* Loading the xml file */
  docxml = xmlParseFile(filename);

  if (docxml == NULL) {

    bft_error(__FILE__, __LINE__, 0,
              _("Unable to parse the file: %s\n"), filename);
    argerr = 2;

  }
  else {

    /* Contexte definition */
    xpathCtx = xmlXPathNewContext(docxml);

    /* Get the root node of the xml document and more particularly
       of its label */
    node = xmlDocGetRootElement(docxml);
    xmlRootName = (const char*) node->name;

  }

  /* Check the Interface version */
  cs_gui_check_version();

#else

  bft_error(__FILE__, __LINE__, 0,
            _("%s was built without XML support,\n"
              "so parameter file \"%s\" may not be loaded.\n"),
            "Code_Saturne", filename);

  argerr = 1;

#endif

  _setup_read = true;
  return argerr;
}

/*-----------------------------------------------------------------------------
 * Check the xml file version.
 *----------------------------------------------------------------------------*/

void
cs_gui_check_version(void)
{
  double version_sat = XML_READER_VERSION;
  double major, maj_sat;

  cs_tree_node_t *tn = cs_glob_tree;
  tn = cs_tree_get_node(tn, "/Code_Saturne_GUI");
  if (tn == NULL)
    tn = cs_tree_get_node(cs_glob_tree, "/NEPTUNE_CFD_GUI");

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

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with the root node.
 *
 * returns:
 *   the root path
 *----------------------------------------------------------------------------*/

char*
cs_xpath_init_path(void)
{
#if defined(HAVE_LIBXML2)

  char *path = NULL;

  BFT_MALLOC(path, strlen("/") + strlen(xmlRootName) + 1, char);
  strcpy(path, "/");
  strcat(path, xmlRootName);

  return path;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return NULL;

#endif
}

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with a short way.
 *
 * returns:
 *   the short path.
 *----------------------------------------------------------------------------*/

char*
cs_xpath_short_path(void)
{
  char *path = NULL;

  BFT_MALLOC(path, strlen("/")+1, char);
  strcpy(path, "/");

  return path;
}

/*----------------------------------------------------------------------------
 * Add all elements (*) to the path.
 *
 * parameter:
 *   path <--> path for the xpath request
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_all_elements(char **path)
{
  assert(path);

  BFT_REALLOC(*path,
              strlen(*path) +strlen("/*") +1,
              char);

  strcat(*path, "/*");
}

/*----------------------------------------------------------------------------
 * Add an element (i.e. markup's label) to the path.
 *
 * parameters:
 *   path    <-> path for the xpath request
 *   element <-- label of the new element in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_element(char        **path,
                     const char   *element)
{
  assert(path);

  if (element != NULL) {

    BFT_REALLOC(*path,
                strlen(*path)+ strlen(element)+ strlen("/") +1,
                char);

    strcat(*path, "/");
    strcat(*path, element);

  }
}

/*----------------------------------------------------------------------------
 * Add a list of elements (i.e. markup's label) to the path.
 *
 * parameters:
 *   path <->  path for the xpath request
 *   nbr  <--  size of the labels list
 *   ...  <--  list of labels of new elements in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_elements(char  **path,
                      int     nbr,
                      ...)
{
  va_list list;
  char *elt = NULL;
  int i;

  assert(path);

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(*path,
                  strlen(*path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(*path, "/");
      strcat(*path, elt);
    }
  }

  va_end(list);
}

/*----------------------------------------------------------------------------
 * Add an element's attribute to the path.
 *
 * parameters:
 *   path           <-> path for the xpath request
 *   attribute_name <-- label of the new attribute in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_attribute(char        **path,
                       const char   *attribute_name)
{
  assert(path);
  assert(attribute_name);

  BFT_REALLOC(*path,
              strlen(*path)+ strlen(attribute_name)+ strlen("/@")+1,
              char);

  strcat(*path, "/@");
  strcat(*path, attribute_name);
}

/*----------------------------------------------------------------------------
 * Add the i'th element to the path.
 *
 * parameters:
 *   path    <-> path for the xpath request
 *   element <-- label of the new element in the path
 *   num     <-- number of the element's markup
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_element_num(char        **path,
                         const char   *element,
                         int           num)
{
  int   nfigures = 0;
  char *strnum = NULL;

  assert(path);
  assert(element);

  nfigures = cs_gui_characters_number(num);

  BFT_MALLOC(strnum, nfigures+1, char);

  BFT_REALLOC(*path,
              strlen(*path)+
              strlen("/")+
              strlen(element)+
              strlen("[")+
              nfigures+
              strlen("]")+1,
              char);

  strcat(*path, "/");
  strcat(*path, element);
  sprintf(strnum,"%d", num);
  strcat(*path, "[");
  strcat(*path, strnum);
  strcat(*path, "]");

  BFT_FREE(strnum);
}

/*----------------------------------------------------------------------------
 * Add a test on a value associated to an attribute to the path.
 *
 * parameters:
 *   path            <-> path for the xpath request
 *   attribute_type  <-- label of the attribute for the test in the path
 *   attribute_value <-- value of the attribute for the test in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_test_attribute(char        **path,
                            const char   *attribute_type,
                            const char   *attribute_value)
{
  assert(path);
  assert(attribute_type);
  assert(attribute_value);

  BFT_REALLOC(*path, strlen(*path)+
              strlen("[@")+
              strlen(attribute_type)+
              strlen("='")+
              strlen(attribute_value)+
              strlen("']")+1,
              char);

  strcat(*path, "[@");
  strcat(*path, attribute_type);
  strcat(*path, "='");
  strcat(*path, attribute_value);
  strcat(*path, "']");
}

/*----------------------------------------------------------------------------
 * Add the 'text()' xpath function to the path.
 *
 * parameters:
 *   path <->  path for the xpath request
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_function_text(char  **path)
{
  assert(path);

  BFT_REALLOC(*path,
              strlen(*path)+ strlen("/text()")+1,
              char);

  strcat(*path, "/text()");
}

/*----------------------------------------------------------------------------
 * Return a list of attribute node names from the xpath request in an array.
 *
 * Example: from <a attr="c"/><b attr="d"/> return {c,d}
 *
 * parameters:
 *   path <-- path for the xpath request
 *   size --> array size
 *----------------------------------------------------------------------------*/

char **
cs_gui_get_attribute_values(char  *path,
                            int   *size)
{
#if defined(HAVE_LIBXML2)

  char             **nodes_name = NULL;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         cur;
  xmlXPathObjectPtr  xpathObj;
  int                i;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nodes = xpathObj->nodesetval;
  *size  = (nodes) ? nodes->nodeNr : 0;

  if (*size != 0) {

    BFT_MALLOC(nodes_name, *size, char*);

    for (i =0; i < *size; i++) {
      assert(nodes->nodeTab[0]);

      if (nodes->nodeTab[i]->type == XML_ATTRIBUTE_NODE) {
        cur = nodes->nodeTab[i];
        BFT_MALLOC(nodes_name[i],
                   strlen((char *) cur->children->content)+1,
                   char);
        strcpy(nodes_name[i], (char*) cur->children->content);
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_ATTRIBUTE_NODE.\nXpath: %s\n"),
                  path);
    }
  }

  xmlXPathFreeObject(xpathObj);
  return nodes_name;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return NULL;

#endif
}

/*----------------------------------------------------------------------------
 * Return the value of an element's attribute.
 *
 * Example: from <a b="c"/> return c
 *
 * parameters:
 *   path <-- path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_attribute_value(char  *path)
{
  char **array = NULL;
  char  *attr = NULL;
  int    size;

  array = cs_gui_get_attribute_values(path, &size);

  if ((array == NULL) || (size == 0))
    return NULL;

  if (size > 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Several attributes found: %i \n"
                "The first one is %s \nXpath: %s\n"),
                size, array[0], path);

  BFT_MALLOC(attr, strlen(array[0])+1, char);
  strcpy(attr, array[0]);

  for (int i = 0; i < size; i++)
    BFT_FREE(array[i]);
  BFT_FREE(array);

  return attr;
}

/*----------------------------------------------------------------------------
 * Return a list of children nodes name from the xpath request in an array.
 *
 * Example: from <a>3<\a><b>4<\b> return {a,b}
 *
 * parameters:
 *   path <-- path for the xpath request
 *   size --> array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_nodes_name(char  *path,
                      int   *size)
{
#if defined(HAVE_LIBXML2)

  char             **nodes_name = NULL;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         cur;
  int                i;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nodes = xpathObj->nodesetval;
  *size = (nodes) ? nodes->nodeNr : 0;

  if (*size != 0) {

    BFT_MALLOC(nodes_name, *size, char*);

    for (i =0; i < *size; i++) {
      assert(nodes->nodeTab[0]);

      if (nodes->nodeTab[i]->type == XML_ELEMENT_NODE) {
        cur = nodes->nodeTab[i];
        BFT_MALLOC(nodes_name[i], strlen((const char *) cur->name)+1, char);
        strcpy(nodes_name[i], (const char*) cur->name);
      } else
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_ELEMENT_NODE.\nXpath: %s\n"),
                  path);
    }
  }

  xmlXPathFreeObject(xpathObj);
  return nodes_name;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return NULL;

#endif
}


/*----------------------------------------------------------------------------
 * Return a single node's name from the xpath request.
 *
 * parameter:
 *   path <-- path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_node_name(char  *path)
{
  char **array = NULL;
  char  *name = NULL;
  int    size;

  array = cs_gui_get_nodes_name(path, &size);

  if ((array == NULL) || (size == 0))
    return NULL;

  if (size > 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Several nodes name found: %i \n"
                "The first one is %s \nXpath: %s\n"),
                size, array[0], path);

  BFT_MALLOC(name, strlen(array[0])+1, char);
  strcpy(name, array[0]);

  BFT_FREE(array[0]);
  BFT_FREE(array);

  return name;
}

/*----------------------------------------------------------------------------
 * Return a list of children text nodes from the xpath request in an array.
 *
 * Example: from <a>3<\a><a>4<\a> return {3,4}
 *
 * parameters:
 *   path <-- path for the xpath request
 *   size --> array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_text_values(char  *path,
                       int   *size)
{
#if defined(HAVE_LIBXML2)

  char             **text_name = NULL;
  xmlXPathObjectPtr  xpathObj;
  xmlNodeSetPtr      nodes;
  xmlNodePtr         cur;
  int                i;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nodes = xpathObj->nodesetval;
  *size = (nodes) ? nodes->nodeNr : 0;

  if (*size != 0) {

    BFT_MALLOC(text_name, *size, char*);

    for (i =0; i < *size; i++) {

      assert(nodes->nodeTab[0]);

      if (nodes->nodeTab[i]->type == XML_TEXT_NODE) {
        cur = nodes->nodeTab[i];
        BFT_MALLOC(text_name[i], strlen((char *) cur->content)+1, char);
        strcpy(text_name[i], (char*) cur->content);
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_TEXT_NODE.\nXpath: %s\n"),
                  path);
    }
  }

  xmlXPathFreeObject(xpathObj);
  return text_name;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return NULL;

#endif
}

/*----------------------------------------------------------------------------
 * Return a single child text node from the xpath request.
 *
 * parameter:
 *   path <-- path for the xpath request
 *
 * returns:
 *   child node based on request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_text_value(char  *path)
{
  char **array = NULL;
  char  *text = NULL;
  int    size;

  array = cs_gui_get_text_values(path, &size);

  if ((array == NULL) || (size == 0))
    return NULL;

  if (size > 1)
    bft_error(__FILE__, __LINE__, 0,
            _("Several text node found: %i \n"
              "The first one is %s \nXpath: %s\n"),
              size, array[0], path);

  BFT_MALLOC(text, strlen(array[0])+1, char);
  strcpy(text, array[0]);

  BFT_FREE(array[0]);
  BFT_FREE(array);

  return text;
}

/*----------------------------------------------------------------------------
 * Query a double value parameter.
 *
 * parameters:
 *   path  <-- path for the xpath request
 *   value --> double result of the xpath request
 *
 * returns:
 *   1 if the xpath request succeeded, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_gui_get_double(char    *path,
                  double  *value)
{
  char *text_name = NULL;
  int   test;

  text_name = cs_gui_get_text_value(path);

  if (text_name == NULL)
    test = 0;
  else {
    *value = atof(text_name);
    BFT_FREE(text_name);
    test = 1;
  }
  return test;
}

/*----------------------------------------------------------------------------
 * Query an integer value parameter.
 *
 * parameters:
 *   path  <-- path for the xpath request
 *   value --> double result of the xpath request
 *
 * returns:
 *   1 if the xpath request succeeded, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_gui_get_int(char  *path,
               int   *value)
{
  char *text_name = NULL;
  int   test;

  text_name = cs_gui_get_text_value(path);

  if (text_name == NULL)
    test = 0;
  else {
    *value = atoi(text_name);
    BFT_FREE(text_name);
    test = 1;
  }
  return test;
}

/*----------------------------------------------------------------------------
 * Query the number of elements (i.e. the number of xml markups)
 * from a xpath request.
 *
 * Example: from <a>3<\a><a>4<\a> return 2
 *
 * parameters:
 *   path <-- path for the xpath request
 *
 * returns:
 *   the number of elements in xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_nb_element(char  *path)
{
#if defined(HAVE_LIBXML2)

  xmlXPathObjectPtr xpathObj;
  int nb;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nb = (xpathObj->nodesetval) ? xpathObj->nodesetval->nodeNr : 0;

  xmlXPathFreeObject(xpathObj);

  return nb;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return 0;

#endif
}

/*----------------------------------------------------------------------------
 * Query the maximum integer value from an xpath request result list.
 *
 * Example: from <a>3<\a><a>4<\a> return 4
 *
 * parameters:
 *   path <-- path for the xpath request
 *
 * returns:
 *   the maximum integer value in the list
 *----------------------------------------------------------------------------*/

int
cs_gui_get_max_value(char  *path)
{
#if defined(HAVE_LIBXML2)

  xmlXPathObjectPtr xpathObj;
  xmlNodeSetPtr     nodes;
  xmlNodePtr        cur;
  int               max_val=0;
  int               size;
  int               i;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nodes=xpathObj->nodesetval;
  size = (nodes) ? nodes->nodeNr : 0;

  if (size == 0)
    bft_error (__FILE__, __LINE__, 0, _("No markup found: %s \n"), path);

  else {
    for (i=0; i <size; i++) {
      assert(nodes->nodeTab[i]);
      if(nodes->nodeTab[i]->type == XML_TEXT_NODE) {
        cur = nodes->nodeTab[i];
        max_val = CS_MAX(max_val, atoi((char*) cur->content));
      }
      else {
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_TEXT_NODE.\nXpath: %s\n"), path);
      }
    }
  }

  xmlXPathFreeObject(xpathObj);

  return max_val;

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Code_Saturne has been compiled without XML support."));

  return 0;

#endif
}

/*-----------------------------------------------------------------------------
 * Evaluate the "status" attribute value.
 *
 * parameter:
 *   path   <-- path for the xpath request
 *   result --> status="on" return 1, status="off" return 0
 *
 * returns:
 *   1 if the xpath request has succeeded, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_gui_get_status(char  *path,
                  int   *result)
{
  char *status;
  int   istatus;

  status = cs_gui_get_attribute_value(path);

  if (status == NULL)
    istatus = 0;
  else {
    istatus = 1;

    if (cs_gui_strcmp(status, "on"))
      *result = 1;
    else if (cs_gui_strcmp(status, "off"))
      *result = 0;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid attribute value: %s \nXpath: %s\n"), status, path);

    BFT_FREE(status);
  }

  return istatus;
}

/*-----------------------------------------------------------------------------
 * Return the xml markup quantity.
 *
 * parameters:
 *   markup <--  path for the markup
 *   flag   <--  1: initialize the path with the root node;
 *               0: initialize the path with a short way
 *
 * returns:
 *   XML markup quantity
 *----------------------------------------------------------------------------*/

int
cs_gui_get_tag_count(const char  *markup,
                     int          flag)
{
  char *path = NULL;
  int   number = 0;

  if (flag) {
    path = cs_xpath_init_path();
  } else {
    BFT_MALLOC(path, strlen("/")+1, char);
    strcpy(path, "/");
  }
  cs_xpath_add_element(&path, markup);
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number;
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
  if ( strlen(s1) != strlen(s2)) return 0;
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
