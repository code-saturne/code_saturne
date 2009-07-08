/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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

#ifndef __CS_GUI_UTIL_H__
#define __CS_GUI_UTIL_H__

/*============================================================================
 * Management of the GUI parameters file: xpath request and utilities
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the information if the requested xml file is missing
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIHMP (IIHMPR)
 * *****************
 *
 * INTEGER          IIHMPR   <--   1 if the file exists, 0 otherwise
 *----------------------------------------------------------------------------*/

void
CS_PROCF (csihmp, CSIHMP) (int *const iihmpr);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Load the xml file in memory. Return an error code for the main programme.
 *
 * parameter:
 *   filename            -->  xml file containing the parameters
 *----------------------------------------------------------------------------*/

int
cs_gui_file_loading(const char *const filename);

/*----------------------------------------------------------------------------
 * Check the xml file version.
 *----------------------------------------------------------------------------*/

void
cs_gui_get_version(void);

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with the root node.
 * Return the root path.
 *----------------------------------------------------------------------------*/

char*
cs_xpath_init_path(void);

/*----------------------------------------------------------------------------
 * Initialize the path for the xpath request with a short way.
 * Return the short path.
 *----------------------------------------------------------------------------*/

char*
cs_xpath_short_path(void);

/*----------------------------------------------------------------------------
 * Add all element (*) to the path.
 *
 * parameter:
 *   path               <-->  path for the xpath request
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_all_elements(char **path);

/*----------------------------------------------------------------------------
 * Add an element (i.e. markup's label) to the path.
 *
 * parameter:
 *   path               <-->  path for the xpath request
 *   element             -->  label of the new element in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_element(      char **     path,
                     const char *const element);

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
                      const int    nbr, ...);

/*----------------------------------------------------------------------------
 * Add an element's attribute to the path.
 *
 * parameter:
 *   path               <-->  path for the xpath request
 *   attribute_name      -->  label of the new attribute in the path
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_attribute(      char **     path,
                       const char *const attribute_name);

/*----------------------------------------------------------------------------
 * Add the i'st element to the path.
 *
 * parameters:
 *   path               <-->  path for the xpath request
 *   element             -->  label of the new element in the path
 *   num                 -->  number of the element's markup
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_element_num(      char **     path,
                         const char *const element,
                         const int         num);

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
                            const char  *const attribute_value);

/*----------------------------------------------------------------------------
 * Add the 'text()' xpath function to the path.
 *
 * parameter:
 *   path               <-->  path for the xpath request
 *----------------------------------------------------------------------------*/

void
cs_xpath_add_function_text(char **path);

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
                            int  *const size);

/*----------------------------------------------------------------------------
 * Return the value of an element's attribute.
 * Example: from <a b="c"/> return c
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_attribute_value(char *const path);

/*----------------------------------------------------------------------------
 * Return a list of children nodes name from the xpath request in an array.
 * Example: from <a>3<\a><b>4<\b> return {a,b}
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   size               <--   array size
 *----------------------------------------------------------------------------*/

char**
cs_gui_get_nodes_name(char   *const path,
                      int    *const size);

/*----------------------------------------------------------------------------
 * Return a single node's name from the xpath request.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_node_name(char *const path);

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
                       int    *const size);

/*----------------------------------------------------------------------------
 * Return a single children text node from the xpath request.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_text_value(char *const path);

/*----------------------------------------------------------------------------
 * Modify the value parameter and return 1 if the xpath request succeeded,
 * otherwise just return 0.
 * Example: from <a>3<\a><a>4<\a> return {3,4}
 *
 * parameters:
 *   path                -->  path for the xpath request
 *   value              <--   double result of the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_double(char   *const path,
                  double *const value);

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
               int  *const value);

/*----------------------------------------------------------------------------
 * Return the number of elements (i.e. the number of xml markups)
 * from a xpath request.
 * Example: from <a>3<\a><a>4<\a> return 2
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_nb_element(char *const path);

/*----------------------------------------------------------------------------
 * Return the integer max value from a list, which is a xpath request result.
 * Example: from <a>3<\a><a>4<\a> return 4
 *
 * parameter:
 *   path                -->  path for the xpath request
 *----------------------------------------------------------------------------*/

int
cs_gui_get_max_value(char *const path);

/*-----------------------------------------------------------------------------
 * Evaluate the "status" attribute value.
 * Return 1 if the xpath request has succeeded, 0 otherwise.
 *
 * parameter:
 *   path                -->  path for the xpath request
 *   result             <--   status="on" return 1,  status="off" return 0
 *----------------------------------------------------------------------------*/


int
cs_gui_get_status(char *const path,
                  int  *const result);

/*-----------------------------------------------------------------------------
 * Returns the number of markup described with a path
 *
 * parameters:
 *   markup             -->  path for the markup
 *   flag               -->  1: initialize the path with the root node;
 *                           0: initialize the path with a short way
 *----------------------------------------------------------------------------*/

int
cs_gui_get_tag_number(const char *const markup, const int flag);

/*-----------------------------------------------------------------------------
 *  Return the number of sign needed to write an integer number
 *
 * parameter:
 *   num                -->  integer number
 *----------------------------------------------------------------------------*/

int
cs_gui_characters_number(const int num);

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
              const char *const s2);

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
                  const int         lstrF);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_UTIL_H__ */
