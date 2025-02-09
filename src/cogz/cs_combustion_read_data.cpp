/*============================================================================
 * Gas combustion model: read thermochemical data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_assert.h"
#include "base/cs_field.h"
#include "base/cs_file.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "gui/cs_gui_specific_physics.h"
#include "pprt/cs_physical_model.h"
#include "pprt/cs_combustion_model.h"
#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_bsh.h"
#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_read_data.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_read_data.cpp

  \brief Combustion combustion model: setup data from input.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Point to next line in text buffer
 *
 * \param[in]       name       name of associated file (for error logging)
 * \param[in, out]  line_num   line number counter
 * \param[in, out]  start      any position in current line before
 *                             line break in, pointer to start of line out
 * \param[in]       last       last possible position in buffer
 * \param[in]       allow_end  if true, return null if last position in buffer
 *                             is reached; otherwise, this is considered a
 *                              premature end
 *
 * \return  pointer to start of next line, or null if end reached
 */
/*----------------------------------------------------------------------------*/

static void
_next_line(const char*  name,
           int         &line_num,
           char*       &s,
           char*        last,
           bool         allow_end)
{
  line_num++;

  char *n = s;
  while (n < last) {
    if (*n == '\n') {
      n++;
      break;
    }
    else
      n++;
  }
  if (n == last) {
    if (allow_end)
      n = nullptr;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Premature end of \"%s\" at line %d."), name, line_num);
  }

  s = n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract i-th token from buffer.
 *
 * A null-termination character is added to to end of the token unless
 * the end of line has been reached, so the buffer values are modified.
 * The pointer to the next portion of the string is set just after this
 * null character.
 *
 * \param[in]       buffer_name      buffer name (for error logging)
 * \param[in]       token_name       token type name (for error logging)
 * \param[in, out]  next             pointer to next part of buffer
 * \param[in]       last             pointer to last byte of buffer
 * \param[in]       line_num         line number (for error logging)
 * \param[in]       idx              token index in line (for error logging)
 */
/*----------------------------------------------------------------------------*/

static char *
_extract_token(const char*   buffer_name,
               const char*   token_name,
               char*        &next,
               char*         last,
               int           line_num,
               int           idx)
{
  char *tok = nullptr;

  char *ss = next;
  while (ss < last && (*ss == ' ' || *ss == '\t'))
    ss++;
  tok = ss;
  char *se = ss;
  while (se < last && (*se != ' ' && *se != '\t' && *se != '\n' && *se != '\r'))
    se++;
  if (se == ss)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: error reading %s %d in %s, line %d."),
              __func__, token_name, idx+1, buffer_name, line_num);
  else if (isblank(*se)) {
    *se = '\0';
    se++;
  }

  next = se;
  return tok;
}

/*----------------------------------------------------------------------------
 * Solve Ax = B for a 4x4 system.
 *
 * Here, a is stored in column major due to caller constraints.
 *
 * parameters:
 *   n                  <-- row and column size
 *   a                  <-- matrix (column major)
 *   b                  <-- right hand side
 *   x                  --> solution of Ax = b
 *----------------------------------------------------------------------------*/

static void
_solve_ax_b_gauss(int               n,
                  double            a[],
                  double  *restrict b,
                  double  *restrict x)
{
  double factor;

  double _a[  CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES
            * CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES];
  double _b[CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES];
  double _a_swap[CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES];
  double _b_swap;

  const double _epsilon = 1.e-24;

  /* Copy array and RHS (avoid modifying a and b,
     and switch to row major) */

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      _a[j*n + i] = a[i*n + j];
    }
    _b[i] = b[i];
  }

  /* Forward elimination */

  for (int i = 0; i < n-1; i++) {

    /* Seek best pivot */

    int k_pivot = i;
    double abs_pivot = fabs(_a[i*n + i]);

    for (int k = i+1; k < n; k++) {
      double abs_a_ki = fabs(_a[k*n + i]);
      if (abs_a_ki > abs_pivot) {
        abs_pivot = abs_a_ki;
        k_pivot = k;
      }
    }

    /* Swap if necessary */

    if (k_pivot != i) {
      for (int j = 0; j < n; j++) {
        _a_swap[j] = _a[i*n + j];
        _a[i*n + j] = _a[k_pivot*n + j];
        _a[k_pivot*n + j] = _a_swap[j];
      }
      _b_swap = _b[i];
      _b[i] = _b[k_pivot];
      _b[k_pivot] = _b_swap;
    }

    if (abs_pivot < _epsilon) {
      bft_error(__FILE__, __LINE__, 0,
                _("%s: n non-zero pivot for Gaussian elimination."), __func__);
    }

    /* Eliminate values */

    for (int k = i+1; k < n; k++) {
      factor = _a[k*n + i] / _a[i*n + i];
      _a[k*n + i] = 0.0;
      for (int j = i+1; j < n; j++) {
        _a[k*n + j] -= _a[i*n + j]*factor;
      }
      _b[k] -= _b[i]*factor;
    }
  }

  /* Solve system */

  for (int k = n-1; k > -1; k--) {
    x[k] =  _b[k];
    for (int j = k+1; j < n; j++)
      x[k] -= _a[k*n + j]*x[j];
    x[k] /= _a[k*n + k];
  }
}

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read gas combustion thermochemistry data.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_read_data(void)
{
  /* Initializations and checks
     -------------------------- */

  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

  const int ngasem = CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS;
  const int ngasgm = CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES;
  const int npot = CS_COMBUSTION_GAS_MAX_TABULATION_POINTS;
  const int nrgasm = CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS;
  int cm_type = (int)cm->type;

  // rank of fuel in the r-th reaction
  int igfuel[CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS];
  // rank of oxydizer in the r-th reaction
  int igoxy[CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS];
  // rank of products in the r-th reaction
  int igprod[CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS];

  int    iereac[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
  double wmolce[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  double cpgase[CS_COMBUSTION_GAS_MAX_TABULATION_POINTS]
               [CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
  double ehgase[CS_COMBUSTION_GAS_MAX_TABULATION_POINTS]
               [CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

  // Stoichiometry in reaction global species.  Negative for the reactants,
  // and positive for the products.
  double stoeg[CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS]
              [CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  // stoichiometric coefficient of global species
  double nreact[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];

  char nomcoe[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS][13];
  char nomcog[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES][151];

  cs_real_t epsi = 1.e-9;

  // Initialializations to 0 (should not be absolutely necessary,
  // but make for nicer printing under a debugger.

  for (int it = 0; it < npot; it++) {
    for (int ige = 0; ige < ngasem; ige++) {
      cpgase[it][ige] = 0.;
      ehgase[it][ige] = 0.;
    }
  }

  for (int ige = 0; ige < ngasem; ige++) {
    iereac[ige] = 0;
    wmolce[ige] = 0;
  }

  for (int ir = 0; ir < nrgasm; ir++) {
    for (int igg = 0; igg < ngasgm; igg++)
      stoeg[ir][igg] = 0.;
    igfuel[ir] = 0;
    igoxy[ir] = 0;
    igprod[ir] = 0;
  }

  for (int igg = 0; igg < ngasgm; igg++) {
    nreact[igg] = 0;
  }

  memset(nomcoe, 0, sizeof(nomcoe));
  memset(nomcog, 0, sizeof(nomcog));

  // FIXME: check if we should really call this function here,
  // as it is already called before, and this call could overwrite
  // user-defined settings.
  cs_gui_combustion_gas_model_temperatures();

  const char *path = cm->data_file_name;

  /* Read buffer */

  char *buf = nullptr;
  cs_gnum_t f_size;
  if (cs_glob_rank_id < 1)
    f_size = cs_file_size(path);
  cs_parall_bcast(0, 1, CS_GNUM_TYPE, &f_size);

  if (f_size <= 0)
    bft_error(__FILE__, __LINE__, 0,
              _("File \"%s\" seems empty."), path);

  CS_MALLOC(buf, f_size + 1, char);

  cs_file_t *f = cs_file_open_serial(path, CS_FILE_MODE_READ);
  cs_file_read_global(f, buf, 1, f_size);
  buf[f_size] = '\0';
  f = cs_file_free(f);

  int line_num = 1;
  char *line = buf, *last = &buf[f_size], *ss = nullptr;  // pointers in buffer

  /* Using JANAF with local thermochemistry reaction description
     ----------------------------------------------------------- */

  if (cm->use_janaf) {

    // number of elementary gas constituents
    int ns = sscanf(line, "%d", &(cm->n_gas_el_comp));
    if (ns !=1)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error reading number of gases in %s, line %d"),
                __func__, path, line_num);
    else if (cm->n_gas_el_comp > CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: the number species read in %s is %d, but it must not be\n"
           "greater than CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS = %d."),
         __func__, path, cm->n_gas_el_comp,
         CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS);

    // number of enthalpy-temperature tabulation points

    _next_line(path, line_num, line, last, false);
    ns = sscanf(line, "%d", &(cm->n_tab_points));
    if (ns !=1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error reading number of tabulation points in %s, line %d."),
         __func__, path, line_num);
    else if (cm->n_tab_points > CS_COMBUSTION_GAS_MAX_TABULATION_POINTS)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: the number of tabulation points read in %s is %d,\n"
           "but it must not be greater than\n"
           "CS_COMBUSTION_GAS_MAX_TABULATION_POINTS = %d."),
         __func__, path, cm->n_tab_points,
         CS_COMBUSTION_GAS_MAX_TABULATION_POINTS);

    // Min and max temperature

    double tmin, tmax;
    _next_line(path, line_num, line, last, false);
    ns = sscanf(line, "%lg", &tmin);
    if (ns < 1)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error reading tmin in %s, line %d."),
         __func__, path, line_num);
    _next_line(path, line_num, line, last, false);
    ns = sscanf(line, "%lg", &tmax);
    if (ns < 1)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error reading tmax in %s, line %d."),
         __func__, path, line_num);

    // Names of elementary gas components

    _next_line(path, line_num, line, last, false); // skip line

    _next_line(path, line_num, line, last, false); ss = line;
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      char *tok = _extract_token(path, "component name",
                                 ss, last, line_num, ige);
      int l = 0;
      while (isalnum(tok[l])) {
        if (l < 13)  // truncate string if needed
          nomcoe[ige][l] = tok[l];
        l++;
      }
    }

    // Absorption coefficient of current species

    double kabse[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

    _next_line(path, line_num, line, last, false); ss = line;
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      char *tok = _extract_token(path, "absorption coefficient",
                                 ss, last, line_num, ige);
      ns = sscanf(tok, "%lg", &(kabse[ige]));
      if (ns < 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: error parsing absorption coeffient in %s, line %d."),
                  __func__, path, line_num);
    }

    // Number of atomic species

    _next_line(path, line_num, line, last, false);
    ns = sscanf(line, "%d", &(cm->n_atomic_species));
    if (ns !=1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error reading number of atomic species in %s, line %d."),
         __func__, path, line_num);
    else if (cm->n_atomic_species > CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: the number of tabulation points read in %s is %d,\n"
           "but it must not be greater than\n"
           "CS_COMBUSTION_GAS_MAX_TABULATION_POINTS = %d."),
         __func__, path, cm->n_atomic_species,
         CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES);

    // Molar mass of atomic species
    // Composition of current species relative to elementary species

    double atgase[CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES]
                 [CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

    for (int iat = 0; iat < cm->n_atomic_species; iat++) {
      _next_line(path, line_num, line, last, false); ss = line;
      char *tok = _extract_token(path, "molar mass",
                                 ss, last, line_num, 0);
      ns = sscanf(tok, "%lg", &(cm->wmolat[iat]));
      if (ns < 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error parsing molar mass in %s, line %d."),
           __func__, path, line_num);

      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        tok = _extract_token(path, "composition",
                             ss, last, line_num, ige+1);
        ns = sscanf(tok, "%lg", &(atgase[iat][ige]));
        if (ns < 1)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error parsing composition (column %d) in %s, line %d."),
             __func__, ige+2, path, line_num);
      }
    }

    // Number of global species

    _next_line(path, line_num, line, last, false);
    ns = sscanf(line, "%d", &(cm->n_gas_species));
    if (ns !=1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error reading number of global species in %s, line %d."),
         __func__, path, line_num);

    // We consider a single global reaction
    // so n_gas_species = 3 (F, O, P)

    // Either the user balances the reaction, in which case 3 global species
    // are read, or we compute a balance and the user must provide the 2
    // global reactive species.

    if (cm->n_gas_species < 2 || cm->n_gas_species > 3)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: the number of global species can only be equal to\n"
           "  2 (automatic reaction balance) or\n"
           "  3 (composition of known products)\n"
           " but the number of global species read in %s is %d."),
         __func__, path, cm->n_gas_species);

    /* a) If the 3 global species are known, the user provides the
     *    stoichiometry in moles of global species. */

    if (cm->n_gas_species == 3) {

      // Composition of global species based on current species
      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        _next_line(path, line_num, line, last, false); ss = line;
        for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
          char *tok = _extract_token(path, "composition",
                                     ss, last, line_num, ige);
          ns = sscanf(tok, "%lg", &(cm->compog[igg][ige]));
          if (ns < 1)
            bft_error
              (__FILE__, __LINE__, 0,
               _("%s: error parsing composition (column %d) in %s, line %d."),
               __func__, ige+1, path, line_num);
        }
      }

      // Number of global reactions

      _next_line(path, line_num, line, last, false);
      ns = sscanf(line, "%d", &(cm->n_reactions));
      if (ns !=1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error reading number of global reactions in %s, line %d."),
           __func__, path, line_num);

      // We consider a single global reaction

      if (cm->n_reactions != 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: the number of global reactions must be 1\n"
             "but the value given in %s is %d."),
           __func__, path, cm->n_reactions);

      // Numbers of species in stoichiometry,
      // given in global reaction species.
      // igfuel[ir] : global fuel species of reaction ir
      // igoxy[ir]  : global oxydant of reaction ir
      // igprod[ir] : global product of reaction ir

      for (int ir = 0; ir < cm->n_reactions; ir++) {
        _next_line(path, line_num, line, last, false); ss = line;
        char *tok = _extract_token(path, "igfuel",
                                   ss, last, line_num, ir);
        ns = sscanf(tok, "%d", &(igfuel[ir]));
        if (ns < 1)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error parsing igfuel (column %d) in %s, line %d."),
             __func__, 1, path, line_num);
        igfuel[ir] -= 1;

        tok = _extract_token(path, "igoxy",
                             ss, last, line_num, ir);
        ns = sscanf(tok, "%d", &(igoxy[ir]));
        if (ns < 1)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error parsing igoxy (column %d) in %s, line %d."),
             __func__, 2, path, line_num);
        igoxy[ir] -= 1;

        for (int igg = 0; igg < cm->n_gas_species; igg++) {
          tok = _extract_token(path, "coefficient",
                                     ss, last, line_num, igg+2);
          ns = sscanf(tok, "%lg", &(stoeg[ir][igg]));
          if (ns < 1)
            bft_error
              (__FILE__, __LINE__, 0,
               _("%s: error parsing composition (column %d) in %s, line %d."),
               __func__, igg+3, path, line_num);
        }
      }

      // Product
      igprod[0] = 2;

    }

    /* b) If the user provides only the reactive species, the reaction
     *    is balanced automatially. */

    else if (cm->n_gas_species == 2) {

      // NComposition of global species based on current species.
      // igfuel[ir]  : global fuel species of reaction ir
      // igoxy[ir]   : global oxydant of reaction ir
      // igprod[ir]  : global product of reaction ir
      // iereac[igg] : reactive elementary species of global species igg

      //  Fuel
      igfuel[0] = 0;
      iereac[igfuel[0]] = 0;

      _next_line(path, line_num, line, last, false); ss = line;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        char *tok = _extract_token(path, "coefficient",
                                   ss, last, line_num, ige);
        ns = sscanf(tok, "%lg", &(cm->compog[igfuel[0]][ige]));
        if (ns < 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("%s: error parsing coeffient in %s, line %d."),
                  __func__, path, line_num);
      }

      if (cm->compog[igfuel[0]][iereac[igfuel[0]]] <= 0) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: the reactive elementary species %d is not present in the\n"
             "global species %d.\n"
             "Its ratio is %g in file %s."),
           __func__, iereac[igfuel[0]], igfuel[0],
           cm->compog[igfuel[0]][iereac[igfuel[0]]], path);
      }

      // Oxydant
      igoxy[0] = 1;
      iereac[igoxy[0]] = 1;

      _next_line(path, line_num, line, last, false); ss = line;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        char *tok = _extract_token(path, "coefficient",
                                   ss, last, line_num, ige);
        ns = sscanf(tok, "%lg", &(cm->compog[igoxy[0]][ige]));
        if (ns < 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("%s: error parsing coeffient in %s, line %d."),
                  __func__, path, line_num);
      }

      if (cm->compog[igoxy[0]][iereac[igoxy[0]]] <= 0) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: the reactive elementary species %d is not present in the\n"
             "global species %d.\n"
             "Its ratio is %g in file %s."),
           __func__, iereac[igoxy[0]], igoxy[0],
           cm->compog[igoxy[0]][iereac[igoxy[0]]], path);
      }

      // Product
      igprod[0] = 2;
      iereac[igprod[0]] = -1;

      // The composition is computed
      // If we have more than 3 elementary species in the products,
      // we need additional data (part of product based on combustible).

      _next_line(path, line_num, line, last, false); ss = line;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        char *tok = _extract_token(path, "coefficient",
                                   ss, last, line_num, ige);
        ns = sscanf(tok, "%lg", &(cm->compog[igprod[0]][ige]));
        if (ns < 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("%s: error parsing coeffient in %s, line %d."),
                  __func__, path, line_num);
      }

      // These species must be placed at the end to avoid 0 pivots when
      // inverting the matrx (we could maybe move them...)

      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        if (cm->compog[igprod[0]][ige] > 0. && ige < 5)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: the species whose part is specified must be positionned\n"
               "last among the products.\n"
               "Elementary species %s is at position %d in the compositions\n"
               "of the global species product in file %s.n"),
             __func__, nomcoe[ige], ige+1, path);
      }

      //! Number of global reactions: we can only balance one reaction.
      cm->n_reactions = 1;

    }

    // Effective Heat of Combustion (J/kg)

    cm->pcigas = 0;
    _next_line(path, line_num, line, last, true);
    if (line != nullptr) {
      char *tok = strstr(line, "EHC");
      if (tok != nullptr) {
        ns = sscanf(tok + 3, "%lg", &(cm->pcigas));
        if (ns < 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("%s: error parsing pcigas in %s, line %d."),
                    __func__, path, line_num);
      }
    }

    /* Compute additional data
       ----------------------- */

    // Compute molar masses of current species

    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      cm->wmole[ige] = 0.;
      for (int iat = 0; iat < cm->n_atomic_species; iat++)
        cm->wmole[ige] += atgase[iat][ige]*cm->wmolat[iat];
    }

    // Compute compostion of products and stoichiometric coefficients

    int igf = igfuel[0];
    int igo = igoxy[0];
    int igp = igprod[0];

    // Automatic equilibrium of the reaction
    // -------------------------------------

    if (cm->n_gas_species == 2) {

      // Add the product global species
      cm->n_gas_species = 3;

      // Composition of known species in products, by unit of combustible mass.
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
        cm->compog[igp][ige] *= cm->wmole[igf] / cm->wmole[ige];
      }

      // Compute stoichiometric coefficients
      double aa[  CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES
                * CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      double bb[  CS_COMBUSTION_GAS_MAX_ATOMIC_SPECIES];
      double xx[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

      const int nato = cm->n_atomic_species;
      const int ngase = cm->n_gas_el_comp;
      const int n_gas_species = cm->n_gas_species;

      for (int ige = 0; ige < ngase; ige++) {
        for (int iat = 0; iat < nato; iat++)
          aa[ige*nato + iat] = 0.;
        xx[ige] = 0.;
      }
      for (int iat = 0; iat < nato; iat++) {
        bb[iat] = 0.;
      }

      for (int iat = 0; iat < nato; iat++) {
        for (int ige = 0; ige < ngase; ige++) {
          // fuel
          bb[iat] = bb[iat] + atgase[iat][ige]*cm->compog[igf][ige];
          // oxydizer
          aa[iereac[igo]*nato+iat] -= atgase[iat][ige]*cm->compog[igo][ige];
          // products
          bb[iat] -= atgase[iat][ige]*cm->compog[igp][ige];
          if ((  ige != iereac[igf] && ige != iereac[igo])
              && fabs(cm->compog[igp][ige]) <= 0.) {
            aa[ige*nato + iat] += atgase[iat][ige];
          }
        }
      }

      _solve_ax_b_gauss(nato, aa+nato, bb, xx+1);

      // we now know the stoichiometric coefficients of global species.
      nreact[igf] = -1.;
      nreact[igo] = -xx[iereac[igo]];
      nreact[igp] = 1.;

      for (int igg = 0; igg < n_gas_species; igg++) {
        // reactive global species are already known
        if (igg != igf && igg != igo) {
          for (int ige = 0; ige < ngase; ige++) {
            // so are reactive elementary species
            if (ige != iereac[igf] && ige != iereac[igo]
                && fabs(cm->compog[igg][ige]) <= 0.) {
              cm->compog[igg][ige] = fabs(xx[ige]);
            }
          }
        }
      }

      // Stoichiometry in global reaction species
      for (int ir= 0; ir < n_gas_species; ir++) {
        for (int igg = 0; igg < n_gas_species; igg++) {
          stoeg[ir][igg] = 0.;
          for (int ige = 0; ige < ngase; ige++) {
            stoeg[ir][igg] += cm->compog[igg][ige]*nreact[igg];
          }
        }
      }

    }

    // Equilibrium defined by the user
    // -------------------------------

    else {

      // global stoechiometric coefficient, to write the reaction

      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        double nmolg = 0.;
        for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
          nmolg += cm->compog[igg][ige];
        if (fabs(nmolg) <= 0) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error in input data (\"%s\")\n"
               "The number of moles in the global species %d is %g."),
             __func__, path, igg+1, nmolg);
        }
        nreact[igg] = stoeg[0][igg] / nmolg;
      }
    }

    // Reaction logging
    // ----------------

    // fuel
    strcpy(nomcog[igf], nomcoe[igf]);

    // oxydizer
    nomcog[igo][0] = '\0';
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      if (fabs(cm->compog[igo][ige] -1.) < 1e-16) {
        char sbuf[151];
        snprintf(sbuf, 150, "%s %s +", nomcog[igo], nomcoe[ige]);
        strncpy(nomcog[igo], sbuf, 150); nomcog[igo][150] = '\0';
      }
      else if (cm->compog[igo][ige] > 0) {
        char sbuf[256];
        snprintf(sbuf, 256, "%s %6.3f %s +",
                 nomcog[igo], cm->compog[igo][ige], nomcoe[ige]);
        strncpy(nomcog[igo], sbuf, 150); nomcog[igo][150] = '\0';
      }
    }

    // product
    nomcog[igp][0] = '\0';
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      if (cm->compog[igp][ige] > 0) {
        char sbuf[256];
        snprintf(sbuf, 256, "%s %6.3f %s +",
                 nomcog[igp], cm->compog[igp][ige], nomcoe[ige]);
        strncpy(nomcog[igp], sbuf, 150); nomcog[igp][150] = '\0';
      }
    }

    // Temperature discretization

    const int npo = cm->n_tab_points;
    for (int it = 0; it < npo; it++) {
      cm->th[it] = (double)(it)*(tmax-tmin)/(double)(npo-1) + tmin;
    }

    // Enthalpy calculation of elementary species

    // If user specifies a EHC, fuel enthalpy is not computed from pptbht
    // but set to zero waiting its calculation from EHC.
    // Otherwise, it is computed.

    int ncoel = cm->n_gas_el_comp, icoel = 0;

    if (cm->pcigas > 0) {
      ncoel = cm->n_gas_el_comp - 1;
      icoel = 1;
      for (int it = 0; it < npo; it++)
        ehgase[it][0] = 0.;
      // Fuel must be placed in first position in the elementary species.
      if (fabs(cm->compog[igfuel[0]][0] - 1.) > 1e-16) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error in input data (\"%s\")\n"
             "Fuel must be placed at first place in the list of\n"
             "elementary species. First elementary specie is now %s."),
           __func__, path, nomcoe[0]);
      }
    }

    for (int ige = icoel; ige < cm->n_gas_el_comp; ige++)
      wmolce[ige] = cm->wmole[ige];

    cs_combustion_enthalpy_and_cp_from_janaf
      (ncoel, ngasem, npo, &(nomcoe[icoel]),
       &(ehgase[0][icoel]), &(cpgase[0][icoel]), &(wmolce[icoel]),
       cm->th);

    // Masses of global species (becoming molar masses below).

    for (int igg = 0; igg < cm->n_gas_species; igg++) {
      cm->wmolg[igg] = 0.;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
        cm->wmolg[igg] += cm->compog[igg][ige] * cm->wmole[ige];
    }

    // Absorption coefficients of global species in case of
    // radiation calculation.
    // TODO check why this is not needed or used with SLFM model.

    if (cm_type%2 == 1 && cm_type/100 != 2) {
      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        cm->ckabsg[igg] = 0.;
        for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
          cm->ckabsg[igg] += cm->compog[igg][ige]*kabse[ige]*cm->wmole[ige];
        cm->ckabsg[igg] /= cm->wmolg[igg];
      }
    }

    // Enthalpies and mass heat capacity of global species

    for (int igg = icoel; igg < cm->n_gas_species; igg++) {
      for (int it = 0; it < npo; it++) {
        cm->eh_gas_g[it][igg] = 0.;
        cm->cp_gas_g[it][igg] = 0.;
        for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
          cm->eh_gas_g[it][igg]
            += cm->compog[igg][ige]*cm->wmole[ige]*ehgase[it][ige];
          cm->cp_gas_g[it][igg]
            += cm->compog[igg][ige]*cm->wmole[ige]*cpgase[it][ige];
        }
        cm->eh_gas_g[it][igg] /= cm->wmolg[igg];
        cm->cp_gas_g[it][igg] /= cm->wmolg[igg];
      }
    }

    // Molar masses of global species

    for (int igg = 0; igg < cm->n_gas_species; igg++) {
      double nmolg = 0.;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
        nmolg += cm->compog[igg][ige];
      if (fabs(nmolg) <= 0) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error in input data (\"%s\")\n"
             "The number of moles in the global species %d is %g."),
           __func__, path, igg+1, nmolg);
      }
      cm->wmolg[igg] /= nmolg;
      for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
        cm->compog[igg][ige] /= nmolg;
    }

    // Estimation of fuel enthalpy in case of user EHC

    if (cm->pcigas > 0) {
      for (int it = 0; it < npo; it++) {
        cm->eh_gas_g[it][0] = 0.;
        for (int igg = icoel; igg < cm->n_gas_species; igg++)
          cm->eh_gas_g[it][0] += stoeg[0][igg]*cm->wmolg[igg]
                                              *cm->eh_gas_g[it][igg];
        cm->eh_gas_g[it][0] =   cm->pcigas
                              - cm->eh_gas_g[it][0]/(cm->wmolg[0]*stoeg[0][0]);
      }
    }

    // Molar coefficients XCO2 , XH2O

    int iio2  = -1;
    int iic   = -1;
    int iico2 = -1;
    int iih2o = -1;

    const char *gas_name = nullptr;

    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      gas_name = nomcoe[ige];
      if (strcmp(gas_name, "C(S)") == 0) iic = ige;
      // if (strcmp(gas_name, "CO") == 0) iico = ige;
      if (strcmp(gas_name, "O2") == 0) iio2 = ige;
      if (strcmp(gas_name, "CO2") == 0) iico2 = ige;
      if (strcmp(gas_name, "H2O") == 0) iih2o = ige;
    }

    cm->iic = iic+1;
    cm->iio2 = iio2+1;
    cm->iico2 = iico2+1;

    // Not used in model

    cm->xco2 = cm->compog[2][iico2];
    cm->xh2o = cm->compog[2][iih2o];

    // Balance verification and stoichiometric ratio for each reaction

    double mfuel = 0, moxyd = 0, mreac = 0;

    for (int ir = 0; ir < cm->n_reactions; ir++) {
      for (int iat = 0; iat < cm->n_atomic_species; iat++) {
        double balance = 0.;
        for (int igg = 0; igg < cm->n_gas_species; igg++) {
          for (int ige = 0; ige < cm->n_gas_el_comp; ige++)
            balance += stoeg[ir][igg]*cm->compog[igg][ige]*atgase[iat][ige];
        }
        if (fabs(balance) > epsi) {
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error in input data (\"%s\")\n"
               "Conservation problem encountered in reaction %d\n"
               "for element %d\n."
               "The molar balance is %g mol, whereas it should be 0."),
             __func__, path, ir+1, iat+1, balance);
        }
      }
      mfuel = stoeg[ir][igfuel[ir]]*cm->wmolg[igfuel[ir]];
      moxyd = stoeg[ir][igoxy[ir]]*cm->wmolg[igoxy[ir]];
      mreac = mfuel + moxyd;
      cm->fs[ir] = mfuel/mreac;
    }

    // Compute coefficients of the oxydant mass fraction for pdflwc.

    cm->lw.coeff1 = cm->compog[1][1]*cm->wmole[1]/cm->wmolg[1];
    cm->lw.coeff3 = moxyd / mfuel;
    cm->lw.coeff2 = cm->lw.coeff3 * cm->lw.coeff1;

  }

  /* Direct use of an enthalpy-temperature tabulation
     ------------------------------------------------ */

  else {

    // number of enthalpy-temperature tabulation points

    int ns = sscanf(line, "%d", &(cm->n_tab_points));
    if (ns !=1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error reading number of tabulation points in %s, line %d."),
         __func__, path, line_num);
    else if (cm->n_tab_points > CS_COMBUSTION_GAS_MAX_TABULATION_POINTS)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: the number of tabulation points read in %s is %d,\n"
           "but it must not be greater than\n"
           "CS_COMBUSTION_GAS_MAX_TABULATION_POINTS = %d."),
         __func__, path, cm->n_tab_points,
         CS_COMBUSTION_GAS_MAX_TABULATION_POINTS);

    // enthalpy-temperature for global species

    for (int ip = 0; ip < cm->n_tab_points; ip++) {
      _next_line(path, line_num, line, last, false); ss = line;
      char *tok = _extract_token(path, "temperature",
                                 ss, last, line_num, 0);
      ns = sscanf(tok, "%lg", &(cm->th[ip]));
      if (ns < 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error parsing temperature mass in %s, line %d."),
           __func__, path, line_num);

      for (int j = 0; j < 3; j++) {
        tok = _extract_token(path, "enthalpy",
                             ss, last, line_num, j+1);
        ns = sscanf(tok, "%lg", &(cm->eh_gas_g[ip][j]));
        if (ns < 1)
          bft_error
            (__FILE__, __LINE__, 0,
             _("%s: error parsing enthalpy (column %d) in %s, line %d."),
             __func__, j+2, path, line_num);
      }
    }

    // We consider a single global reaction
    // so we have 3 global species (fuel, oxydizer, product)

    cm->n_gas_species = 3;
    cm->n_reactions = 1;

    // Molar masses for global species

    _next_line(path, line_num, line, last, false); ss = line;
    for (int j = 0; j < 3; j++) {
      char *tok = _extract_token(path, "wmolg",
                                 ss, last, line_num, j);
      ns = sscanf(tok, "%lg", &(cm->wmolg[j]));
      if (ns < 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error parsing wmolg (column %d) in %s, line %d."),
           __func__, j+1, path, line_num);
    }

    // Mixture fraction at stoichiometry

    _next_line(path, line_num, line, last, false); ss = line;
    ns = sscanf(line, "%lg", &(cm->fs[0]));
    if (ns < 1)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error parsing fs in %s, line %d."),
         __func__, path, line_num);

    // Absorption coefficients for global species

    _next_line(path, line_num, line, last, false); ss = line;
    for (int j = 0; j < 3; j++) {
      char *tok = _extract_token(path, "ckabs",
                                 ss, last, line_num, j);
      ns = sscanf(tok, "%lg", &(cm->ckabsg[j]));
      if (ns < 1)
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: error parsing ckabsg (column %d) in %s, line %d."),
           __func__, j+1, path, line_num);
    }

    // Molar CO2 and H2O coefficients in products (for radiation)

    _next_line(path, line_num, line, last, false); ss = line;
    ns = sscanf(line, "%lg %lg", &(cm->xco2), &(cm->xh2o));
    if (ns < 2)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: error parsing xco2 and xh2o in %s, line %d."),
         __func__, path, line_num);

    // Compute oxygen mass fraction coefficients
    // We consider that the oxydant is a mix of O2 and N2

    cm->lw.coeff1 = ((cm->wmolg[1]-0.028)/(0.032-0.028)) * (0.032/cm->wmolg[1]);
    cm->lw.coeff3 = (1-cm->fs[0]) / cm->fs[0];
    cm->lw.coeff2 = cm->lw.coeff3 * cm->lw.coeff1;

    // FIXME: we need som additional coefficient computations.
    // Otherwise, we have a division by 0 later, due to having
    // at least m->n_gas_el_comp = 0, and cm->coefeg = 0.

    // This mode is broken since at least 2014-11-20, so
    // if nobody has noticed this in 10 years, perhaps we can simply
    // remove this "non-JANAF" mode, which seems to be unused.
  }

  CS_FREE(buf);

  // Convert coefficients from global species to elementary species

  for (int igg = 0; igg < cm->n_gas_species; igg++) {
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      cm->coefeg[igg][ige]
        = cm->compog[igg][ige]*cm->wmole[ige]/cm->wmolg[igg];
    }
  }

  // gas name storage
  const char *gas_name = nomcoe[0];

  // PCI calculation

  if (cm->use_janaf) {

    cm->pcigas = 0.;

    cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
    for (int i = 0; i < ngasem; i++)
      coefg[i] = 0;

    for (int ir = 0; ir < cm->n_reactions; ir++) {
      for (int igg = 0; igg < cm->n_gas_species; igg++) {

        // formation enthalpies
        coefg[0]  = 0.0;
        coefg[1]  = 0.0;
        coefg[2]  = 0.0;
        coefg[igg] = 1.0;
        cs_real_t tgas = 300.0;
        cs_real_t efgas = cs_gas_combustion_t_to_h(coefg, tgas);

        cm->pcigas += stoeg[ir][igg]*cm->wmolg[igg]*efgas;
      }

      // PCI dimension is J/kg of combustible
      cm->pcigas /= (stoeg[ir][0]*cm->wmolg[0]);
    }
  }

  /* Burke Schumann Relationships
     ---------------------------- */

  if (   cm->type == CS_COMBUSTION_BSH_ADIABATIC
      || cm->type == CS_COMBUSTION_BSH_PERMEATIC)
    cs_burke_schumann();

  /* Logging
     ------- */

  // TODO: move this to standard setup_log calls once C migration is
  // finished, so required fields may be saved in combustion model
  // structure without complex mappings.

  if (cm->use_janaf && cm_type/100 == 1) { // 3-point chemistry

    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    " ** SPECIFIC PHYSICS:\n"
                    "    ----------------\n\n"));

    cs_log_printf(CS_LOG_SETUP,
                  _(" --- Diffusion Flame: 3 Point Chemistry\n"
                    "       Option = %d\n\n"),
                  (cm->type)%100);

    cs_log_printf(CS_LOG_SETUP,
                  _(" --- Combustible characteristics\n"
                    "       Combustible :  %s\n"
                    "       PCI =          %14.5e J/kg\n\n"),
                  gas_name, cm->pcigas);

    cs_log_printf(CS_LOG_SETUP,
                  _(" --- Chemical reaction: \n"
                    "       %s + %6.3f (%s) --> %s\n\n"),
                  nomcog[igfuel[0]], nreact[igoxy[0]],
                  nomcog[igoxy[0]], nomcog[igprod[0]]);

    cs_log_printf
      (CS_LOG_SETUP,
       _("    Mass composition           Fuel       Oxydizer       Products\n"
         "    ----------------           ----       --------       --------\n"));
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      cs_log_printf(CS_LOG_SETUP, "    %16s ", nomcoe[ige]);
      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        cs_log_printf(CS_LOG_SETUP, " %14.5f", cm->coefeg[igg][ige]);
      }
      cs_log_printf(CS_LOG_SETUP, "\n");
    }
    cs_log_printf(CS_LOG_SETUP, "\n");

    cs_log_printf
      (CS_LOG_SETUP,
       _("   Molar composition           Fuel       Oxydizer       Products\n"
         "   -----------------           ----       --------       --------\n"));
    for (int ige = 0; ige < cm->n_gas_el_comp; ige++) {
      cs_log_printf(CS_LOG_SETUP, "    %16s ", nomcoe[ige]);
      for (int igg = 0; igg < cm->n_gas_species; igg++) {
        cs_log_printf(CS_LOG_SETUP, " %14.5f", cm->compog[igg][ige]);
      }
      cs_log_printf(CS_LOG_SETUP, "\n");
    }
    cs_log_printf(CS_LOG_SETUP, "\n");
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
