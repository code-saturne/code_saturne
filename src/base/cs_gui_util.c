/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Reader of the parameters file: xpath request and utilities
 *============================================================================*/


#if defined(_CS_HAVE_XML)


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
 * BFT library headers
 *----------------------------------------------------------------------------*/


#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>


/*----------------------------------------------------------------------------
 * libxml2 library headers
 *----------------------------------------------------------------------------*/


#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/


#include "cs_base.h"


/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/


#include "cs_gui_util.h"


/*----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/


#define XML_READER_VERSION 0.0


/*----------------------------------------------------------------------------
 * Global variables
 *----------------------------------------------------------------------------*/

xmlDocPtr docxml            = NULL;   /* Pointer on the XML document  */
xmlXPathContextPtr xpathCtx = NULL;   /* Pointer on the XPath Context */
xmlNodePtr node             = NULL;   /* Pointer on the root node     */
const char *xmlRootName     = NULL;   /* Name of the root node        */

/*============================================================================
 * Fortran API public functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the information if the requested xml file is missing or not.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIHMP (ITURB, IDEUCH, IGRAKE, IGRAKI)
 * *****************
 *
 * INTEGER          IIHMPR   <--   1 if the file exists, 0 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF (csihmp, CSIHMP) (int *const iihmpr)
{
  if (docxml == NULL)
    *iihmpr = 0;
  else
    *iihmpr = 1;
}

/*============================================================================
 * C API public functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Load the xml file in memory. Return an error code for the main program.
 *
 * parameter:
 *   filename            -->  xml file containing the parameters
 *----------------------------------------------------------------------------*/

int
cs_gui_file_loading(const char *const filename)
{
  int file_descriptor = 0;
  int argerr = 0;

  assert(filename);

  /* printf("numero rang proc = %i\n", (int)cs_glob_base_rang); */

  /* Vérification de l'existence du fichier par son ouverture */
  file_descriptor = open(filename, O_RDONLY);

  if (file_descriptor ==  -1) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf( _("Unable to open the file: %s\n"), filename);
    argerr = 2;
    return argerr;

  } else {

  /* Si le fichier existe, on le referme. Il sera réouvert par xmlParseFile */
    close(file_descriptor);

  }

  /* libxml initialization */
  xmlInitParser();
  LIBXML_TEST_VERSION

  /* Loading the xml file */
  docxml = xmlParseFile(filename);


  if (docxml == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf (_("Unable to parse the file: %s\n"), filename);
    argerr = 2;

  } else {

   /* Contexte definition */
   xpathCtx = xmlXPathNewContext(docxml);

   /* Récupération du noeud racine du document xml et
      plus particulièrement de son label */
   node = xmlDocGetRootElement(docxml);
   xmlRootName = (const char*) node->name;
   }
   /* Verification de la version de l'interface */
   cs_gui_get_version();

   return argerr;
}

/*-----------------------------------------------------------------------------
 * Check the xml file version.
 *----------------------------------------------------------------------------*/

void
cs_gui_get_version(void)
{
  char *path;
  char *version;
  double version_number;
  double version_sat = XML_READER_VERSION;
  double major;
  double maj_sat;
  double min_sat;
  double minus;

  path = cs_xpath_init_path();
  cs_xpath_add_attribute(&path, "version");

  version = cs_gui_get_attribute_value(path);
  version_number = atof(version);

  minus = modf(version_number, &major);
  min_sat = modf(version_sat, &maj_sat);

  if(major != maj_sat)
     bft_error(__FILE__, __LINE__, 0,
               _("========================================================\n"
                 "   ** INVALID VERSION OF THE XML FILE\n"
                 "      -------------------------------------- \n"
                 "      XML FILE VERSION: %.1f  \n"
                 "      XML READER VERSION: %.1f \n"
                 "========================================================\n"),
                  version_number, version_sat);

  if(minus != min_sat) {
     cs_base_warn(__FILE__, __LINE__);
     bft_printf(_("========================================================\n"
                 "   ** INCOMPATIBLE VERSION OF THE XML FILE\n"
                 "      -------------------------------------- \n"
                 "      XML FILE VERSION: %.1f  \n"
                 "      XML READER VERSION: %.1f \n"
                 "\n"
                 "      YOU SHOULD RESTART YOUR CALCUL WITH A NEW XML FILE\n"
                 "========================================================\n"),
                 version_number, version_sat);
  }
  BFT_FREE(version);
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with the root node.
 * Return the root path.
 *----------------------------------------------------------------------------*/

char*
cs_xpath_init_path(void)
{
  char *path = NULL;

  BFT_MALLOC(path, strlen("/") + strlen(xmlRootName) + 1, char);
  strcpy(path, "/");
  strcat(path, xmlRootName);

  return path;
}

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with a short way.
 * Return the short path.
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
 * Add all element (*) to the path.
 *
 * parameter:
 *   path               <-->  path for the xpath request
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
 *   path               <-->  path for the xpath request
 *   element             -->  label of the new element in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_element(      char **      path,
                     const char  *const element)
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
 *   path               <-->  path for the xpath request
 *   nbr                 -->  size of the labels list
 *   ...                 -->  list of labels of new elements in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_elements(      char **path,
                      const int    nbr, ...)
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
 *   path               <-->  path for the xpath request
 *   attribute_name      -->  label of the new attribute in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_attribute(      char      **path,
                       const char *const attribute_name)
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
 * Add the i'st element to the path.
 *
 * parameters:
 *   path               <-->  path for the xpath request
 *   element             -->  label of the new element in the path
 *   num                 -->  number of the element's markup
 *----------------------------------------------------------------------------*/

void cs_xpath_add_element_num(      char  **     path,
                              const char  *const element,
                              const int          num)
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
 *   path               <-->  path for the xpath request
 *   attribute_type      -->  label of the attribute for the test in the path
 *   attribute_value     -->  value of the attribute for the test in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_test_attribute(      char **      path,
                            const char  *const attribute_type,
                            const char  *const attribute_value)
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
 * parameter:
 *   path               <-->  path for the xpath request
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_function_text(char **path)
{
  assert(path);

  BFT_REALLOC(*path,
              strlen(*path)+ strlen("/text()")+1,
              char);

  strcat(*path, "/text()");
}

/*----------------------------------------------------------------------------
 * Return a list of attributes nodes name from the xpath request in an array.
 * Example: from <a attr="c"/><b attr="d"/> return {c,d}
 *
 * parameter:
 *   path                -->  path for the xpath request
 *   size               <--   array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_attribute_values(char *const path,
                            int  *const size)
{
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
        BFT_MALLOC(nodes_name[i], strlen((char *) cur->children->content)+1, char);
        strcpy(nodes_name[i], (char*) cur->children->content);
      } else
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_ATTRIBUTE_NODE.\nXpath: %s\n"),
                  path);
    }
  }

  xmlXPathFreeObject(xpathObj);
  return nodes_name;
}

/*----------------------------------------------------------------------------
 * Return the value of an element's attribute.
 * Example: from <a b="c"/> return c
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_attribute_value(char *const path)
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
                size, attr[0], path);

  BFT_MALLOC(attr, strlen(array[0])+1, char);
  strcpy(attr, array[0]);

  BFT_FREE(array[0]);
  BFT_FREE(array);

  return attr;
}

/*----------------------------------------------------------------------------
 * Return a list of children nodes name from the xpath request in an array.
 * Example: from <a>3<\a><b>4<\b> return {a,b}
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   size               <--   array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_nodes_name(char *const path,
                      int  *const size)
{
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
}


/*----------------------------------------------------------------------------
 * Return a single node's name from the xpath request.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_node_name(char *const path)
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
 * Example: from <a>3<\a><a>4<\a> return {3,4}
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   size               <--   array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_text_values(char   *const path,
                       int    *const size)
{
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
      } else
        bft_error(__FILE__, __LINE__, 0,
                  _("The node type is not XML_TEXT_NODE.\nXpath: %s\n"),
                  path);
    }
  }

  xmlXPathFreeObject(xpathObj);
  return text_name;
}

/*----------------------------------------------------------------------------
 * Return a single children text node from the xpath request.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_text_value(char *const path)
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
 * Modify the value parameter and return 1 if the xpath request succeeded,
 * otherwise just return 0.
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   value              <--   double result of the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_double(char   *const path,
                  double *const value)
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
 * Modify the value parameter and return 1 if the xpath request succeeded,
 * otherwise just return 0.
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   value              <--   integer result of the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_int(char *const path,
               int  *const value)
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
 * Return the number of elements (i.e. the number of xml markups)
 * from a xpath request.
 * Example: from <a>3<\a><a>4<\a> return 2
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_nb_element(char *const path)
{
  xmlXPathObjectPtr xpathObj;
  int nb;

  assert(path);

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpathCtx);

  if (xpathObj == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  nb = (xpathObj->nodesetval) ? xpathObj->nodesetval->nodeNr : 0;

  xmlXPathFreeObject(xpathObj);

  return nb;
}

/*----------------------------------------------------------------------------
 * Return the integer max value from a list, which is a xpath request result.
 * Example: from <a>3<\a><a>4<\a> return 4
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

int cs_gui_get_max_value(char *const path)
{
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
}

/*-----------------------------------------------------------------------------
 * Evaluate the "status" attribute value.
 * Return 1 if the xpath request has succeeded, 0 otherwise.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *   result             <--   status="on" return 1, status="off" return 0
 *----------------------------------------------------------------------------*/

int
cs_gui_get_status(char *const path,
                  int  *const result)
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
 * Return the number of characters needed to write an integer number
 *
 * parameter:
 *   num                -->  integer number
 *----------------------------------------------------------------------------*/

int
cs_gui_characters_number(const int num)
{
  int i      = 1;
  int number = 0;

  assert(num > 0);

  if (num == 0)
    number ++;
  else
    for (i=1; i <= num; i *= 10)
      number++;

  return number;
}

/*-----------------------------------------------------------------------------
 * Comparison between two string: return 1 if the two string are equal, 0
 * otherwise.
 *
 * parameters:
 *   s1                -->  first string
 *   s2                -->  second string
 *----------------------------------------------------------------------------*/

int
cs_gui_strcmp(const char *const s1,
              const char *const s2)
{
  if (s1 == NULL || s2 == NULL) return 0;
  if ( strlen(s1) != strlen(s2)) return 0;
  if (!strncmp(s1, s2, strlen(s1))) return 1;
  return 0;
}

/*-----------------------------------------------------------------------------
 * Copy a C string into a Fortran string.
 *
 * parameters:
 *   chainef          <--> Fortran string
 *   chainc            -->  C string
 *   lstrF             -->  maximum length of the Fortran string
 *----------------------------------------------------------------------------*/

void
cs_gui_strcpy_c2f(      char *const chainef,
                  const char *const chainec,
                  const int         lstrF)
{
  int i;

  assert(chainec != NULL);
  assert(lstrF > 0);

  strncpy(chainef, chainec, strlen(chainec));

  for (i = strlen(chainec); i < lstrF ; i++)
    chainef[i] = ' ';
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Stocke dans values le resultat de la requete path
 * Exemple pour :
 *                             <balise>3</balise>
 *                             <balise>4</balise>
 * stocke dans le tableau de flottant values les valeurs 3 et 4
 *
 * retourne 0 ou le nombre de balises
 *----------------------------------------------------------------------------*/

/*
int
cs_gui_get_double_values(char    *const path,
                         double **const values)
{
  char **text_name;
  int    size;
  int    i;
  int    test = 0;

  cs_gui_get_text_values(path, &text_name, &size);

  BFT_MALLOC(*values, size, double);

  if (text_name == NULL)
    test= 0;
  else {
    for (i=0 ; i < size ; i++){
      (*values)[i]=atof(text_name[i]);
      BFT_FREE(text_name[i]);
    }
    BFT_FREE(text_name);
    test = size;
  }
  return test;
}
*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_HAVE_XML */
