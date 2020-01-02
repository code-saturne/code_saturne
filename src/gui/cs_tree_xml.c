/*============================================================================
 * Tree structure and XML file mapping.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_tree.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_tree_xml.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_tree_xml.c
        Tree structure and XML file mapping.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type and structure definitions
 *============================================================================*/

typedef struct {

  const char   *buffer_name;    /*!< pointer to name of buffer */

  char         *buf;            /*!< pointer to working buffer
                                  (length size+1, NULL-terminated) */

  size_t        size;           /*!< buffer size */
  size_t        byte;           /*!< current byte in buffer */
  size_t        line;           /*!< current line in buffer */

  char          s_char;         /*!< current starting character */

  int           depth;          /*!< depth */
  bool          have_attrs;     /*!< have attributes for tags */
  bool          first;          /*!< first node (descriptor) */

  cs_tree_node_t  *node;        /*!< current node */
  cs_tree_node_t  *parent;      /*!< parent node */

} cs_xml_t;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle element for a given XML buffer node.
 *
 * \param[in, out]  n        node to which XML contents are added
 * \param[in, out]  doc      XML parser structure
 * \param[in]       is_attr  if true, key-value attribute
 * \param[in]       name     node name
 * \param[in]       value    associated value if leaf, or NULL
 *
 * \return pointer to extracted string
 */
/*----------------------------------------------------------------------------*/

static void
_handle_element(cs_xml_t    *doc,
                bool         is_attr,
                const char  *name,
                const char  *value)
{
  if (is_attr) {
    assert(doc->node != NULL);
    cs_tree_node_set_tag(doc->node, name, value);
  }
  else {
    if (doc->node == NULL) {
      if (doc->first) {
        if (doc->parent->parent != NULL)
          doc->parent = doc->parent->parent;
        doc->first = false;
      }
      doc->node = cs_tree_add_child_str(doc->parent,
                                        name,
                                        value);
    }
    else if (doc->parent == NULL) {
      /* We choose here to place the top (descriptor) and associated
         tags in a child node, to avoid renaming the root node;
         its sub-nodes will be raised one level so as to be merged in
         the input tree */
      doc->parent = doc->node;
      doc->node = cs_tree_add_child_str(doc->node, name, value);
    }
    else {
      if (name != NULL)
        cs_tree_node_set_name(doc->node, name);
      if (value != NULL)
        cs_tree_node_set_value_str(doc->node, value);
    }
  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)

  int i = doc->depth;
  if (is_attr)
    i += 1;

  while (i > 0) {
    printf("  ");
    i--;
  }

  if (name != NULL)
    printf("%s:", name);
  else
    printf("  ");
  if (value != NULL) {
    if (doc->have_attrs)
      printf("  %s", value);
    else
      printf("%s", value);
  }

  printf("\n");

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find next entry in XML parser.
 *
 * \param[in, out]  doc      XML parser structure
 */
/*----------------------------------------------------------------------------*/

static void
_next(cs_xml_t  *doc)
{
  size_t i = doc->byte;

  /* find start if required */

  while (i < doc->size && isspace(doc->buf[i])) {
    if (doc->buf[i] == '\n')
      doc->line += 1;
    i++;
  }

  if (doc->buf[i] == '<' || doc->buf[i] == '>') {
    doc->s_char = doc->buf[i];
    doc->buf[i] = '\0';
    i += 1;
  }
  else if (i > doc->byte)
    doc->s_char = '\0';

  doc->byte = i;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check for comment and find following entry if present.
 *
 * \param[in, out]  doc     XML parser structure
 */
/*----------------------------------------------------------------------------*/

static void
_check_and_skip_comment(cs_xml_t  *doc)
{
  if (doc->byte < doc->size - 3) {

    assert(doc->s_char == '<');

    /* Skip comments */

    size_t i = doc->byte;

    if (doc->buf[i] == '!' && doc->buf[i+1] == '-' && doc->buf[i+2] == '-') {
      printf("%c%c%c\n", doc->buf[i], doc->buf[i+1], doc->buf[i+2]);
      bool closed = false;
      while (i < doc->size && closed == false) {
        while (i < doc->size && doc->buf[i] != '>')
          i++;
        if (i > 1) {
          if (doc->buf[i-1] == '-' && doc->buf[i-2] == '-')
            closed = true;
          else
            i++;
        }
      }
      if (i < doc->size)
        i++;
      _next(doc);
      if (doc->s_char != '<')
        _next(doc);
    }

    doc->byte = i;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read attribute string in XML parser.
 *
 * \param[in, out]  doc      XML parser structure
 *
 * \return pointer to next attribute key, or NULL
 */
/*----------------------------------------------------------------------------*/

static const char *
_read_attr_key(cs_xml_t  *doc)
{
  /* find start */

  _next(doc);

  size_t i = doc->byte;

  if (doc->s_char == '>')
    return NULL;

  assert(! isspace(doc->buf[i]));

  /* find '=' sign */

  doc->byte = i;
  size_t si = i;

  while (doc->buf[i] != '=' && i < doc->size) {
    if (doc->buf[i] == '\n')
      doc->line += 1;
    else if (isspace(doc->buf[i])) /* remove possible trailing blanks ? */
      doc->buf[i] = '\0';
    i++;
  }

  if (i >= doc->size)
    bft_error(__FILE__, __LINE__, 0,
              _("In XML data \"%s\", line %d"
                "malformed or unhandled key: %s ..."),
              doc->buffer_name, (int)(doc->line), doc->buf + doc->byte);

  assert(doc->buf[i] == '=');
  doc->s_char = doc->buf[i];
  doc->buf[i] = '\0';
  doc->byte = i+1;

  return doc->buf + si;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read string in XML file.
 *
 * \param[in, out]  doc       XML parser structure
 *
 * \return pointer to extracted string
 */
/*----------------------------------------------------------------------------*/

static const char *
_read_string(cs_xml_t  *doc)
{
  /* find start */

  _next(doc);

  size_t _si = doc->byte, _ei = 0;

  /* quoted string */

  if (doc->buf[doc->byte] == '"') {
    size_t i = doc->byte + 1;
    _si = i;
    while (i < doc->size) {
      if (doc->buf[i] == '"') {
        _ei = i;
        doc->buf[i] = '\0';
        doc->byte = i+1;
        doc->s_char = '\0';
        break;
      }
      else if (doc->buf[i] == '>' || doc->buf[i] == '<') {
        doc->buf[i+1] = '\0';
        bft_error(__FILE__, __LINE__, 0,
                  _("In XML data (%s, line %d)\n"
                    "malformed string: %s"),
                  doc->buffer_name, (int)doc->line, doc->buf + doc->byte);
        return NULL;
      }
      else {
        if (doc->buf[i] == '\n')
          doc->line += 1;
        i++;
      }
    }
    doc->buf[_si - 1] = '\0';
  }

  /* unquoted string */

  else {
    size_t i = doc->byte;
    while (i < doc->size) {
      doc->s_char = doc->buf[i];
      if (doc->buf[i] == '<') {
        _ei = i;
        doc->buf[i] = '\0';
        doc->byte = i+1;
        break;
      }
      else if (doc->buf[i] == '"') {
        doc->buf[i+1] = '\0';
        bft_error(__FILE__, __LINE__, 0,
                  _("In XML data (%s, line %d)\n"
                    "malformed string: %s"),
                  doc->buffer_name, (int)doc->line, doc->buf + doc->byte);
        return NULL;
      }
      else {
        if (doc->buf[i] == '\n')
          doc->line += 1;
        i++;
      }
    }
  }

  /* restore escaped characters */

  if (_ei > _si  + 1) {

    size_t j = _si;
    for (size_t i = _si; i < _ei; i++) {
      if (doc->buf[i] == '&') {
        if (strncmp(doc->buf+i, "&quot;", 6) == 0) {
          doc->buf[j++] = '"';
          i += 5;
        }
        else if (strncmp(doc->buf+i, "&apos;", 6) == 0) {
          doc->buf[j++] = '\'';
          i += 5;
        }
        else if (strncmp(doc->buf+i, "&amp;", 5) == 0) {
          doc->buf[j++] = '&';
          i += 4;
        }
        else if (strncmp(doc->buf+i, "&lt;", 4) == 0) {
          doc->buf[j++] = '<';
          i += 3;
        }
        else if (strncmp(doc->buf+i, "&gt;", 4) == 0) {
          doc->buf[j++] = '>';
          i += 3;
        }
        else {
          assert(0);
        }
      }
      else
        doc->buf[j++] = doc->buf[i];
    }
    assert(j <= _ei);
    doc->buf[j] = '\0';
  }

  return (doc->buf + _si);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recursively read next XML buffer tag.
 *
 * \param[in, out]  doc     XML parser structure
 * \param[out]      closed  true if this tag ends with closing '/'
 *
 * \return pointer to extracted string
 */
/*----------------------------------------------------------------------------*/

static const char *
_read_tag(cs_xml_t  *doc,
          bool      *closed)
{
  const char *tag = NULL;
  *closed = false;

  /* Search for next start tag */

  if (doc->s_char != '<')
    _next(doc);

  _check_and_skip_comment(doc);

  size_t i = doc->byte;

  if (doc->s_char != '<' && i < doc->size) {
    if (i < doc->size + 64)
      doc->buf[i+64] = '\0';
    bft_error(__FILE__, __LINE__, 0,
              _("In XML data (%s, line %d)\n"
                "expected tag instead of: %s"),
              doc->buffer_name, (int)(doc->line), doc->buf + i);
  }

  /* Start parsing tag */

  size_t s_id = i;

  while (i < doc->size) {
    doc->s_char = doc->buf[i];
    if (doc->buf[i] == '>' || isspace(doc->buf[i])) {
      if (doc->buf[i] == '>' && doc->buf[i-1] == '/')
        *closed = true;
      doc->buf[i] = '\0';
      tag = doc->buf + s_id;
      if (tag[0] == '/')
        *closed = true;
      if (*closed == false)
        _handle_element(doc, false, tag, NULL);
      i++;
      break;
    }
    i++;
  }
  doc->byte = i;

  doc->have_attrs = false;

  /* Possible attributes */

  while (doc->s_char != '>' && doc->byte < doc->size) {
    _next(doc);
    if (doc->s_char == '>')
      break;
    i = doc->byte;
    if (doc->buf[i] == '/' && doc->buf[i+1] == '>') {
      *closed = true;
      doc->s_char = '>';
      doc->byte = i+2;
    }
    else {
      const char *key = NULL, *val = NULL;
      key = _read_attr_key(doc);
      if (key != NULL) {
        val = _read_string(doc);
        if (val != NULL) {
          doc->have_attrs = true;
          _handle_element(doc, true, key, val);
        }
      }
    }
  }
  if (   doc->s_char != '>'
      && (doc->depth > 0 || doc->byte < doc->size))
    bft_error(__FILE__, __LINE__, 0,
              _("In XML data (%s, line %d)\n"
                "malformed tag: %s"),
              doc->buffer_name, (int)(doc->line), tag);

  return tag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Recursively read next XML buffer node.
 *
 * \param[in, out]  doc  XML parser structure
 *
 * \return pointer to extracted string
 */
/*----------------------------------------------------------------------------*/

static const char *
_read_element(cs_xml_t  *doc)
{
  bool tag_closed;

  const char *tag = _read_tag(doc, &tag_closed);
  const char *tag_c = NULL;

  if (tag_closed)
    return tag;

  /* Search for next start tag or value */

  while (tag_c == NULL && doc->byte < doc->size) {
    _next(doc);
    if (doc->s_char == '<') {
      if (doc->buf[doc->byte] == '/')
        tag_c = _read_tag(doc, &tag_closed);
      else {
        doc->depth++;
        doc->parent = doc->node;
        doc->node = NULL;
        (void)_read_element(doc);
        doc->node = doc->parent;
        doc->parent = doc->node->parent;
        doc->depth--;
      }
    }
    else {
      const char *val = _read_string(doc);
      if (val != NULL) {
        tag_c = _read_tag(doc, &tag_closed);
        _handle_element(doc, false, NULL, val);
      }
    }
  }

  /* Check closing tag match */

  if (tag_c != NULL) {
    if (tag_c[0] == '/') {
      if (strcmp(tag_c + 1, tag) != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("In XML data (%s, line %d)\n"
                    "closing tag <%s> does not match opening tag <%s>"),
                  doc->buffer_name, (int)doc->line, tag_c, tag);
    }
  }

  return tag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read and check XML header.
 *
 * \param[in, out]  doc  XML parser structure
 */
/*----------------------------------------------------------------------------*/

static void
_read_header(cs_xml_t  *doc)
{
  size_t i = doc->byte;
  bool found = false;

  /* find start if required */

  while (doc->byte < doc->size && found == false) {

    while (i < doc->size && doc->buf[i] != '<') {
      if (doc->buf[i] == '\n')
        doc->line += 1;
      i++;
    }

    if (strncmp(doc->buf + doc->byte, "<?xml", 5) == 0) {
      found = true;
      doc->byte += 5;
      while (doc->byte < doc->size) {
        _next(doc);
        if (strncmp(doc->buf + doc->byte, "?>", 2) == 0) {
          doc->byte += 1;
          doc->s_char = '\0';
          doc->buf[doc->byte] = '\0';
          doc->byte += 1;
          break;
        }
        doc->s_char = '\0';
        const char *attr_key = _read_attr_key(doc);
        const char *attr_val = NULL;
        if (attr_key != NULL) {
          attr_val = _read_string(doc);
          if (strcmp(attr_key, "version") == 0) {
            if (strcmp(attr_val, "1.0") != 0)
              bft_printf(_("XML (%s) %s %s unexpected\n"),
                         doc->buffer_name, attr_key, attr_val);
          }
          else if (strcmp(attr_key, "encoding") == 0) {
            if (strcmp(attr_val, "utf-8") != 0)
              bft_printf(_("XML (%s) %s %s unexpected\n"),
                         doc->buffer_name, attr_key, attr_val);
          }
        }
        else
          break;
      }
    }

    /* Variant with no XML version info */

    else if (strncmp(doc->buf + doc->byte, "<", 1) == 0) {
      found = true;
    }

  }

  assert(found);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read and parse XML file to tree.
 *
 * \param[in, out]   r    root node to which XML contents are read
 * \param[in]  path  path to XML file
 */
/*----------------------------------------------------------------------------*/

void
cs_tree_xml_read(cs_tree_node_t  *r,
                 const char       path[])
{
  /* Read buffer */

  cs_xml_t *doc = NULL;
  BFT_MALLOC(doc, 1, cs_xml_t);

  cs_gnum_t f_size;
  if (cs_glob_rank_id < 1)
    f_size = cs_file_size(path);
  cs_parall_bcast(0, 1, CS_GNUM_TYPE, &f_size);

  if (f_size <= 0)
    bft_error(__FILE__, __LINE__, 0,
              _("File \"%s\" seems empty."), path);

  doc->size = f_size;
  BFT_MALLOC(doc->buf, doc->size + 1, char);
  doc->buffer_name = path;
  doc->byte = 0;
  doc->line = 1;
  doc->s_char = '\0';
  doc->depth = 0;
  doc->have_attrs = false;
  doc->first = true;
  doc->node = r;
  doc->parent = NULL;

  cs_file_t *f = cs_file_open_serial(path, CS_FILE_MODE_READ);
  cs_file_read_global(f, doc->buf, 1, f_size);
  f = cs_file_free(f);

  doc->buf[doc->size] = '\0';

  /* Now parse buffer */

  {
    /* Read XML header */
    _read_header(doc);

    /* Now parse tree */

    const char *s = NULL;
    do {
      s = _read_element(doc);
    } while (s != NULL);
  }

  BFT_FREE(doc->buf);
  BFT_FREE(doc);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_tree_dump(CS_LOG_DEFAULT, 0, r);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
