/*============================================================================
 * Management of the GUI parameters file: mobile mesh
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
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "mei_evaluate.h"

#include "cs_ale.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_zone.h"
#include "cs_convection_diffusion.h"
#include "cs_field_pointer.h"
#include "cs_gui.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_mesh.h"
#include "cs_timer.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_mobile_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Static variables
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Possible values for boundary nature
 *----------------------------------------------------------------------------*/

enum ale_boundary_nature
{
  ale_boundary_nature_none,
  ale_boundary_nature_fixed_wall,
  ale_boundary_nature_sliding_wall,
  ale_boundary_nature_internal_coupling,
  ale_boundary_nature_external_coupling,
  ale_boundary_nature_fixed_velocity,
  ale_boundary_nature_fixed_displacement
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return value for ALE method
 *
 * parameters:
 *   tn_ale    <-- tree node associated with ALE
 *
 * return:
 *   0 for isotropic, 1 for orthotropic
 *----------------------------------------------------------------------------*/

static int
_ale_visc_type(cs_tree_node_t  *tn_ale)
{
  int mvisc_type = 0;

  cs_tree_node_t *tn_mv = cs_tree_get_node(tn_ale, "mesh_viscosity");

  const char *type = cs_tree_node_get_tag(tn_mv, "type");
  if (type != NULL) {
    if (strcmp(type, "isotrop") != 0) {
      if (strcmp(type, "orthotrop") == 0)
        mvisc_type = 1;
      else
        bft_error(__FILE__, __LINE__, 0,
                  "invalid mesh viscosity type: %s", type);
    }
  }

  return mvisc_type;
}

/*-----------------------------------------------------------------------------
 * Initialize mei tree and check for symbols existence
 *
 * parameters:
 *   formula        <-- mei formula
 *   symbols        <-- array of symbol to check
 *   symbol_nbr     <-- number of symbol in symbols
 *   variables      <-- variables required in the formula
 *   variable_nbr   <-- number of variable in variables
 *   dtref          <-- time step
 *   ttcabs         <-- current time
 *   ntcabs         <-- current iteration number
 *----------------------------------------------------------------------------*/

static mei_tree_t *
_init_mei_tree(const char    *formula,
               const char   **symbols,
               int            symbol_nbr,
               const char   **variables,
               const double  *variables_value,
               int            variable_nbr,
               const double   dtref,
               const double   ttcabs,
               const int      ntcabs)
{
  int i = 0;

  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* Insert variables into mei_tree */
  for (i = 0; i < variable_nbr; ++i) {
    double value = 0;

    /* Read value from variables_value if it is not null 0 otherwise */
    if (variables_value)
      value = variables_value[i];
    mei_tree_insert(tree, variables[i], value);
  }

  /* Add commun variables: dt, t, nbIter */
  mei_tree_insert(tree, "dt",   dtref);
  mei_tree_insert(tree, "t",    ttcabs);
  mei_tree_insert(tree, "iter", ntcabs);

  /* add variable from notebook */
  cs_gui_add_notebook_variables(tree);

  /* try to build the interpreter */
  if (mei_tree_builder(tree))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not interpret expression: %s\n"), tree->string);

  /* Check for symbols */
  for (i = 0; i < symbol_nbr; ++i) {
    const char* symbol = symbols[i];

    if (mei_tree_find_symbol(tree, symbol)) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not find the required symbol: %s\n"), symbol);
    }
  }

  return tree;
}

/*-----------------------------------------------------------------------------
 * Get the ale boundary formula
 *
 * parameters:
 *   tn_w    <-- pointer to tree node for a given wall BC
 *   choice  <-- nature: "fixed_velocity" or "fixed_displacement"
 *----------------------------------------------------------------------------*/

static const char*
_get_ale_boundary_formula(cs_tree_node_t  *tn_w,
                          const char      *choice)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_w, "ale");
  tn = cs_tree_node_get_sibling_with_tag(tn, "choice", choice);

  const char *formula = cs_tree_node_get_child_value_str(tn, "formula");

  return formula;
}

/*-----------------------------------------------------------------------------
 * Get uialcl data for fixed displacement
 *
 * parameters:
 *   tn_w           <-- pointer to tree node for a given wall BC
 *   begin          <-- begin index for nodfbr
 *   end            <-- end index for nodfbr
 *   nnod           <-> number of nodes
 *   b_face_vtx_lst <-- nodfbr
 *   impale         --> impale
 *   disale         --> disale
 *   dtref          <-- time step
 *   ttcabs         <-- current time
 *   ntcabs         <-- current iteration number
 *----------------------------------------------------------------------------*/

static void
_uialcl_fixed_displacement(cs_tree_node_t   *tn_w,
                           const cs_lnum_t   begin,
                           const cs_lnum_t   end,
                           const cs_lnum_t   b_face_vtx_lst[],
                           int              *impale,
                           cs_real_3_t      *disale,
                           const double      dtref,
                           const double      ttcabs,
                           const int         ntcabs)
{
  const char*  variables[3] = {"mesh_x", "mesh_y", "mesh_z"};
  int variable_nbr = 3;

  /* Get formula */
  const char *formula = _get_ale_boundary_formula(tn_w, "fixed_displacement");

  if (!formula)
    bft_error(__FILE__, __LINE__, 0,
              _("Boundary nature formula is null for %s."),
              cs_gui_node_get_tag(tn_w, "label"));

  /* Init mei */
  mei_tree_t *ev = _init_mei_tree(formula, variables, variable_nbr,
                                  0, 0, 0, dtref, ttcabs, ntcabs);

  mei_evaluate(ev);

  /* Get mei results */
  cs_real_t x_mesh = mei_tree_lookup(ev, "mesh_x");
  cs_real_t y_mesh = mei_tree_lookup(ev, "mesh_y");
  cs_real_t z_mesh = mei_tree_lookup(ev, "mesh_z");

  mei_tree_destroy(ev);

  /* Set disale and impale */
  for (cs_lnum_t ii = begin; ii < end; ++ii) {
    cs_lnum_t inod = b_face_vtx_lst[ii];
    if (impale[inod] == 0) {
      disale[inod][0] = x_mesh;
      disale[inod][1] = y_mesh;
      disale[inod][2] = z_mesh;
      impale[inod] = 1;
    }
  }
}

/*-----------------------------------------------------------------------------
 * Get uialcl data for fixed velocity
 *
 * parameters:
 *   tn_w         <-- pointer to tree node for a given wall BC
 *   iuma         <-- iuma
 *   ivma         <-- ivma
 *   iwma         <-- iwma
 *   nfabor       <-- Number of boundary faces
 *   n_faces      <-- number of selected faces
 *   faces_list   <-- listof selected faces
 *   rcodcl       --> rcodcl
 *   dtref        <-- time step
 *   ttcabs       <-- current time
 *   ntcabs       <-- current iteration number
 *----------------------------------------------------------------------------*/

static void
_uialcl_fixed_velocity(cs_tree_node_t  *tn_w,
                       int              iuma,
                       int              ivma,
                       int              iwma,
                       int              ivimpo,
                       cs_lnum_t        nfabor,
                       cs_lnum_t        n_faces,
                       const cs_lnum_t  faces_list[],
                       int              ialtyb[],
                       double          *rcodcl,
                       double           dtref,
                       double           ttcabs,
                       int              ntcabs)
{
  const char*  variables[3] = {"mesh_velocity_U",
                               "mesh_velocity_V",
                               "mesh_velocity_W" };

  /* Get formula */
  const char *formula = _get_ale_boundary_formula(tn_w, "fixed_velocity");

  if (!formula)
    bft_error(__FILE__, __LINE__, 0,
              _("Boundary nature formula is null for %s."),
              cs_gui_node_get_tag(tn_w, "label"));

  /* Init MEI */
  mei_tree_t *ev = _init_mei_tree(formula, variables, 3, 0, 0, 0,
                                  dtref, ttcabs, ntcabs);

  for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {

    cs_lnum_t ifbr = faces_list[ifac];

    mei_evaluate(ev);

    /* Fill  rcodcl */
    rcodcl[(iuma-1) * nfabor + ifbr] = mei_tree_lookup(ev, "mesh_velocity_U");
    rcodcl[(ivma-1) * nfabor + ifbr] = mei_tree_lookup(ev, "mesh_velocity_V");
    rcodcl[(iwma-1) * nfabor + ifbr] = mei_tree_lookup(ev, "mesh_velocity_W");

    ialtyb[ifbr]  = ivimpo;

  }

  mei_tree_destroy(ev);
}

/*-----------------------------------------------------------------------------
 * Return the ale boundary nature of a given boundary condition
 *
 * parameters:
 *   tn_w <-- pointer to tree node for a given wall BC
 *
 * return:
 *   associated nature
 *----------------------------------------------------------------------------*/

static enum  ale_boundary_nature
_get_ale_boundary_nature(cs_tree_node_t  *tn_w)
{
  enum ale_boundary_nature nature = ale_boundary_nature_none;

  cs_tree_node_t *tn = cs_tree_get_node(tn_w, "ale/choice");
  const char *nat = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(nat, "fixed_boundary"))
    nature = ale_boundary_nature_fixed_wall;
  if (cs_gui_strcmp(nat, "sliding_boundary"))
    nature = ale_boundary_nature_sliding_wall;
  else if (cs_gui_strcmp(nat, "internal_coupling"))
    nature = ale_boundary_nature_internal_coupling;
  else if (cs_gui_strcmp(nat, "external_coupling"))
    nature = ale_boundary_nature_external_coupling;
  else if (cs_gui_strcmp(nat, "fixed_velocity"))
    nature = ale_boundary_nature_fixed_velocity;
  else if (cs_gui_strcmp(nat, "fixed_displacement"))
    nature = ale_boundary_nature_fixed_displacement;

  return nature;
}

/*-----------------------------------------------------------------------------
 * Return the ale boundary type of a given boundary condition
 *
 * parameters:
 *   tn_w <-- pointer to tree node for a given wall BC
 *
 * return:
 *   associated boundary type
 *----------------------------------------------------------------------------*/

static cs_boundary_type_t
_get_ale_boundary_type(cs_tree_node_t  *tn_w)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_w, "ale/choice");
  const char *nat = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(nat, "fixed_boundary"))
    return CS_BOUNDARY_ALE_FIXED;
  else if (cs_gui_strcmp(nat, "sliding_boundary"))
    return CS_BOUNDARY_ALE_SLIDING;
  else if (cs_gui_strcmp(nat, "internal_coupling"))
    return CS_BOUNDARY_ALE_INTERNAL_COUPLING;
  else if (cs_gui_strcmp(nat, "external_coupling"))
    return CS_BOUNDARY_ALE_EXTERNAL_COUPLING;
  else if (cs_gui_strcmp(nat, "fixed_velocity"))
    return CS_BOUNDARY_ALE_IMPOSED_VEL;
  else if (cs_gui_strcmp(nat, "fixed_displacement"))
    return CS_BOUNDARY_ALE_IMPOSED_DISP;
  else
    return CS_BOUNDARY_N_TYPES;
}

/*-----------------------------------------------------------------------------
 * Retrieve internal coupling x, y and z values
 *
 * parameters:
 *   tn_ic <-- tree node for a given BC's internal coupling definitions
 *   name  <-- node name
 *   xyz   --> result matrix
 *----------------------------------------------------------------------------*/

static void
_get_internal_coupling_xyz_values(cs_tree_node_t  *tn_ic,
                                  const char      *name,
                                  double           xyz[3])
{
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_ic, name);

  cs_gui_node_get_child_real(tn, "X", xyz);
  cs_gui_node_get_child_real(tn, "Y", xyz+1);
  cs_gui_node_get_child_real(tn, "Z", xyz+2);
}

/*-----------------------------------------------------------------------------
 * Retrieve the internal coupling matrices
 *
 * parameters:
 *   tn_ic            <-- node for a given BC's internal coupling definitions
 *   name             <-- matrix name
 *   symbols          <-- see _init_mei_tree
 *   symbol_nbr       <-- see _init_mei_tree
 *   variables        <-- see _init_mei_tree
 *   variables_value  <-- see _init_mei_tree
 *   variable_nbr     <-- see _init_mei_tree
 *   output_matrix,   --> result matrix
 *   dtref            <-- time step
 *   ttcabs           <-- current time
 *   ntcabs           <-- current iteration number
 *----------------------------------------------------------------------------*/

static void
_get_internal_coupling_matrix(cs_tree_node_t  *tn_ic,
                              const char      *name,
                              const char      *symbols[],
                              int              symbol_nbr,
                              const char     **variables,
                              const double    *variables_value,
                              int              variable_nbr,
                              double          *output_matrix,
                              double           dtref,
                              double           ttcabs,
                              int              ntcabs)
{
  /* Get the formula */
  mei_tree_t *tree;

  int i = 0;
  cs_tree_node_t  *tn = cs_tree_node_get_child(tn_ic, name);

  const char *matrix = cs_tree_node_get_child_value_str(tn, "formula");

  if (!matrix)
    bft_error(__FILE__, __LINE__, 0,
              _("Formula is null for %s %s"), tn_ic->name, name);

  /* Initialize mei */
  tree = _init_mei_tree(matrix, symbols, symbol_nbr,
                        variables, variables_value, variable_nbr,
                        dtref, ttcabs, ntcabs);
  mei_evaluate(tree);

  /* Read matrix values */
  for (i = 0; i < symbol_nbr; ++i) {
    const char *symbol = symbols[i];
    output_matrix[i] = mei_tree_lookup(tree, symbol);
  }

  mei_tree_destroy(tree);
}

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling for a specific boundary
 *
 * parameters:
 *   label    <-- boundary label
 *   xmstru   --> Mass matrix
 *   xcstr    --> Damping matrix
 *   xkstru   --> Stiffness matrix
 *   forstr   --> Fluid force matrix
 *   istruc   <-- internal coupling boundary index
 *   dtref    <-- time step
 *   ttcabs   <-- current time
 *   ntcabs   <-- current iteration number
 *----------------------------------------------------------------------------*/

static void
_get_uistr2_data(cs_tree_node_t  *tn_ic,
                 double          *xmstru,
                 double          *xcstru,
                 double          *xkstru,
                 double          *forstr,
                 int              istruc,
                 double           dtref,
                 double           ttcabs,
                 int              ntcabs)
{
  const char  *m_symbols[] = {"m11", "m12", "m13",
                              "m21", "m22", "m23",
                              "m31", "m32", "m33"};
  const char  *c_symbols[] = {"c11", "c12", "c13",
                              "c21", "c22", "c23",
                              "c31", "c32", "c33"};
  const char  *k_symbols[] = {"k11", "k12", "k13",
                              "k21", "k22", "k23",
                              "k31", "k32", "k33"};

  int symbol_nbr = sizeof(m_symbols) / sizeof(m_symbols[0]);

  const char   *force_symbols[] = {"fx", "fy", "fz"};
  int force_symbol_nbr = sizeof(force_symbols) / sizeof(force_symbols[0]);

  const int  variable_nbr = 3;
  const char *variables[3] = {"fluid_fx", "fluid_fy", "fluid_fz"};
  double variable_values[3];

  /* Get mass matrix, damping matrix and stiffness matrix */

  _get_internal_coupling_matrix(tn_ic, "mass_matrix", m_symbols,
                                symbol_nbr, 0, 0, 0,
                                &xmstru[istruc * symbol_nbr],
                                dtref, ttcabs, ntcabs);

  _get_internal_coupling_matrix(tn_ic, "damping_matrix", c_symbols,
                                symbol_nbr, 0, 0, 0,
                                &xcstru[istruc * symbol_nbr],
                                dtref, ttcabs, ntcabs);

  _get_internal_coupling_matrix(tn_ic, "stiffness_matrix", k_symbols,
                                symbol_nbr, 0, 0, 0,
                                &xkstru[istruc * symbol_nbr],
                                dtref, ttcabs, ntcabs);

  /* Set variable for fluid force matrix */
  variable_values[0] = forstr[istruc * force_symbol_nbr + 0];
  variable_values[1] = forstr[istruc * force_symbol_nbr + 1];
  variable_values[2] = forstr[istruc * force_symbol_nbr + 2];

  /* Get fluid force matrix */
  _get_internal_coupling_matrix(tn_ic, "fluid_force_matrix",
                                force_symbols, force_symbol_nbr,
                                variables, variable_values, variable_nbr,
                                &forstr[istruc * force_symbol_nbr],
                                dtref, ttcabs, ntcabs);
}

/*-----------------------------------------------------------------------------
 * Return the external coupling dof ("DDL") value
 *
 *  <boundary_conditions>
 *      <wall label=label_argument">
 *          <ale choice="external_coupling">
 *              <node_name_argument choice="off"/>
 *
 * parameters:
 *   label     <-- boundary label
 *   node_name <--  Node name: DDLX, DDLY or DDLZ.
 *----------------------------------------------------------------------------*/

static int
_get_external_coupling_dof(cs_tree_node_t  *tn_ec,
                           const char      *name)
{
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_ec, name);
  const char *choice = cs_tree_node_get_child_value_str(tn, "choice");

  int is_on = cs_gui_strcmp(choice, "on");

  return is_on;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * ALE related keywords
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALIN
 * *****************
 *
 * nalinf  <->  number of sub iteration of initialization
 *              of fluid
 * nalimx  <->  max number of iterations of implicitation of
 *              the displacement of the structures
 * epalim  <->  realtive precision of implicitation of
 *              the displacement of the structures
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialin, UIALIN) (int     *nalinf,
                                int     *nalimx,
                                double  *epalim)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/ale_method");

  cs_gui_node_get_status_int(tn, &cs_glob_ale);

  if (cs_glob_ale) {
    cs_gui_node_get_child_int(tn, "fluid_initialization_sub_iterations",
                              nalinf);
    cs_gui_node_get_child_int(tn, "max_iterations_implicitation",
                              nalimx);
    cs_gui_node_get_child_real(tn, "implicitation_precision",
                               epalim);
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--cs_glob_ale = %i\n", cs_glob_ale);
  if (cs_glob_ale > 0) {
    bft_printf("--nalinf = %i\n", *nalinf);
    bft_printf("--nalimx = %i\n", *nalimx);
    bft_printf("--epalim = %g\n", *epalim);
  }
#endif
}

/*----------------------------------------------------------------------------
 * ALE diffusion type
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALVM
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialvm, UIALVM) ()
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/ale_method");

  int iortvm = _ale_visc_type(tn);

  cs_var_cal_opt_t vcopt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_field_t *f_mesh_u = cs_field_by_name("mesh_velocity");
  cs_field_get_key_struct(f_mesh_u, key_cal_opt_id, &vcopt);

  if (iortvm == 1) { /* orthotropic viscosity */
    vcopt.idften = CS_ANISOTROPIC_LEFT_DIFFUSION;
  } else { /* isotropic viscosity */
    vcopt.idften = CS_ISOTROPIC_DIFFUSION;
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iortvm = %i\n",  iortvm);
  }
#endif
}

/*-----------------------------------------------------------------------------
 * uialcl
 *
 * parameters:
 *   ialtyb       <-- ialtyb
 *   impale       <-- uialcl_fixed_displacement
 *   disale       <-- See uialcl_fixed_displacement
 *   iuma         <-- See uialcl_fixed_velocity
 *   ivma         <-- See uialcl_fixed_velocity
 *   iwma         <-- See uialcl_fixed_velocity
 *   rcodcl       <-- See uialcl_fixed_velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialcl, UIALCL) (const int *const    ibfixe,
                                const int *const    igliss,
                                const int *const    ivimpo,
                                const int *const    ifresf,
                                int       *const    ialtyb,
                                int       *const    impale,
                                cs_real_3_t        *disale,
                                const int *const    iuma,
                                const int *const    ivma,
                                const int *const    iwma,
                                double    *const    rcodcl)
{
  double t0;

  const cs_mesh_t *m = cs_glob_mesh;

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  /* Loop on wall boundary zones */

  for (cs_tree_node_t *tn_w = cs_tree_node_get_child(tn_b0, "wall");
       tn_w != NULL;
       tn_w = cs_tree_node_get_next_of_name(tn_w)) {

    const char *label = cs_tree_node_get_tag(tn_w, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
    if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
      continue;

    cs_lnum_t n_faces = z->n_elts;
    const cs_lnum_t *faces_list = z->elt_ids;

    enum ale_boundary_nature nature = _get_ale_boundary_nature(tn_w);

    if (nature == ale_boundary_nature_none)
      continue;

    if (nature ==  ale_boundary_nature_fixed_wall) {
      for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
        cs_lnum_t ifbr = faces_list[ifac];
        ialtyb[ifbr] = *ibfixe;
      }
    }
    else if (nature ==  ale_boundary_nature_sliding_wall) {
      for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
        cs_lnum_t ifbr = faces_list[ifac];
        ialtyb[ifbr] = *igliss;
      }
    }
    else if (nature == ale_boundary_nature_fixed_displacement) {
      t0 = cs_timer_wtime();
      for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
        cs_lnum_t ifbr = faces_list[ifac];
        _uialcl_fixed_displacement(tn_w,
                                   m->b_face_vtx_idx[ifbr],
                                   m->b_face_vtx_idx[ifbr+1],
                                   m->b_face_vtx_lst, impale, disale,
                                   cs_glob_time_step->dt_ref,
                                   cs_glob_time_step->t_cur,
                                   cs_glob_time_step->nt_cur);
      }
      cs_gui_add_mei_time(cs_timer_wtime() - t0);
    }
    else if (nature == ale_boundary_nature_fixed_velocity) {
      t0 = cs_timer_wtime();
      _uialcl_fixed_velocity(tn_w, *iuma, *ivma, *iwma, *ivimpo,
                             m->n_b_faces, n_faces, faces_list,
                             ialtyb, rcodcl,
                             cs_glob_time_step->dt_ref,
                             cs_glob_time_step->t_cur,
                             cs_glob_time_step->nt_cur);
      cs_gui_add_mei_time(cs_timer_wtime() - t0);
    }
  }

  /* Loop on free surface boundary zones */

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_b0, "free_surface");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
    if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
      continue;

    cs_lnum_t n_faces = z->n_elts;
    const cs_lnum_t *faces_list = z->elt_ids;

    for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
      cs_lnum_t ifbr = faces_list[ifac];
      ialtyb[ifbr]  = *ifresf;
    }

  }
}

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called once at initialization
 *
 * parameters:
 *   idfstr   --> Structure definition
 *   mbstru   <-- number of previous structures (-999 or by restart)
 *   aexxst   --> Displacement prediction alpha
 *   bexxst   --> Displacement prediction beta
 *   cfopre   --> Stress prediction alpha
 *   ihistr   --> Monitor point synchronisation
 *   xstr0    <-> Values of the initial displacement
 *   xstreq   <-> Values of the equilibrium displacement
 *   vstr0    <-> Values of the initial velocity
 *----------------------------------------------------------------------------*/

void CS_PROCF (uistr1, UISTR1) (cs_lnum_t        *idfstr,
                                const int        *mbstru,
                                double           *aexxst,
                                double           *bexxst,
                                double           *cfopre,
                                int              *ihistr,
                                double           *xstr0,
                                double           *xstreq,
                                double           *vstr0)
{
  int  istruct = 0;

  cs_tree_node_t *tn0
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/ale_method");

  /* Get advanced data */
  cs_gui_node_get_child_real(tn0, "displacement_prediction_alpha", aexxst);
  cs_gui_node_get_child_real(tn0, "displacement_prediction_beta", bexxst);
  cs_gui_node_get_child_real(tn0, "stress_prediction_alpha", cfopre);
  cs_gui_node_get_child_status_int(tn0, "monitor_point_synchronisation", ihistr);

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  cs_tree_node_t *tn_b0 = cs_tree_node_get_child(tn, "boundary");
  cs_tree_node_t *tn_w0 = cs_tree_node_get_child(tn, "wall");

  /* At each time-step, loop on boundary faces */

  int izone = 0;

  for (tn = tn_b0;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), izone++) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    cs_tree_node_t *tn_w
      = cs_tree_node_get_sibling_with_tag(tn_w0, "label", label);

    /* Keep only internal coupling */
    if (  _get_ale_boundary_nature(tn_w)
        == ale_boundary_nature_internal_coupling) {

      if (istruct+1 > *mbstru) { /* Do not overwrite restart data */
        /* Read initial_displacement, equilibrium_displacement
           and initial_velocity */
        cs_tree_node_t *tn_ic = cs_tree_get_node(tn_w, "ale");
        tn_ic = cs_tree_node_get_sibling_with_tag(tn_ic,
                                                  "choice",
                                                  "internal_coupling");
        _get_internal_coupling_xyz_values(tn_ic, "initial_displacement",
                                          &xstr0[3 * istruct]);
        _get_internal_coupling_xyz_values(tn_ic, "equilibrium_displacement",
                                          &xstreq[3 * istruct]);
        _get_internal_coupling_xyz_values(tn_ic, "initial_velocity",
                                          &vstr0[3 * istruct]);
      }

      const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
      if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
        continue;

      cs_lnum_t n_faces = z->n_elts;
      const cs_lnum_t *faces_list = z->elt_ids;

      /* Set idfstr to positive index starting at 1 */
      for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
        cs_lnum_t ifbr = faces_list[ifac];
        idfstr[ifbr] = istruct + 1;
      }
      ++istruct;
    }
  }
}

/*-----------------------------------------------------------------------------
 * Retrieve data for internal coupling. Called at each step
 *
 * parameters:
 *   xmstru       --> Mass matrix
 *   xcstr        --> Damping matrix
 *   xkstru       --> Stiffness matrix
 *   forstr       --> Fluid force matrix
 *   dtref        <--  time step
 *   ttcabs       <-- current time
 *   ntcabs       <-- current iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uistr2, UISTR2) (double *const  xmstru,
                                double *const  xcstru,
                                double *const  xkstru,
                                double *const  forstr,
                                double *const  dtref,
                                double *const  ttcabs,
                                int    *const  ntcabs)
{
  int istru   = 0;

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  cs_tree_node_t *tn_b0 = cs_tree_node_get_child(tn, "boundary");
  cs_tree_node_t *tn_w0 = cs_tree_node_get_child(tn, "wall");

  /* At each time-step, loop on boundary faces */

  for (tn = tn_b0; tn != NULL; tn = cs_tree_node_get_next_of_name(tn)) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    cs_tree_node_t *tn_w
      = cs_tree_node_get_sibling_with_tag(tn_w0, "label", label);

    enum ale_boundary_nature nature =_get_ale_boundary_nature(tn_w);

    /* Keep only internal coupling */
    if (nature == ale_boundary_nature_internal_coupling) {

      cs_tree_node_t *tn_ic = cs_tree_get_node(tn_w, "ale");
      tn_ic = cs_tree_node_get_sibling_with_tag(tn_ic,
                                                "choice",
                                                "internal_coupling");

      _get_uistr2_data(tn_ic,
                       xmstru,
                       xcstru,
                       xkstru,
                       forstr,
                       istru,
                       *dtref,
                       *ttcabs,
                       *ntcabs);
      ++istru;
    }
  }
}

/*-----------------------------------------------------------------------------
 * Retrieve data for external coupling
 *
 * parameters:
 *   idfstr    <-- Structure definition
 *   asddlf    --> Block of the DDL forces
 *----------------------------------------------------------------------------*/

void
CS_PROCF(uiaste, UIASTE)(int       *idfstr,
                         cs_int_t  *asddlf)
{
  if (!cs_gui_file_is_loaded())
    return;

  int istruct     = 0;

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  cs_tree_node_t *tn_b0 = cs_tree_node_get_child(tn, "boundary");
  cs_tree_node_t *tn_w0 = cs_tree_node_get_child(tn, "wall");

  /* Loop on boundary faces */

  for (tn = tn_b0; tn != NULL; tn = cs_tree_node_get_next_of_name(tn)) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    cs_tree_node_t *tn_w
      = cs_tree_node_get_sibling_with_tag(tn_w0, "label", label);

    enum ale_boundary_nature nature =_get_ale_boundary_nature(tn_w);

    if (nature == ale_boundary_nature_external_coupling) {

      const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
      if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
        continue;

      cs_lnum_t n_faces = z->n_elts;
      const cs_lnum_t *faces_list = z->elt_ids;

      cs_tree_node_t *tn_ec = cs_tree_get_node(tn_w, "ale");
      tn_ec = cs_tree_node_get_sibling_with_tag(tn_ec,
                                                "choice",
                                                "external_coupling");

      /* Get DDLX, DDLY and DDLZ values */
      asddlf[istruct*3 + 0] = _get_external_coupling_dof(tn_ec, "DDLX") ? 0 : 1;
      asddlf[istruct*3 + 1] = _get_external_coupling_dof(tn_ec, "DDLY") ? 0 : 1;
      asddlf[istruct*3 + 2] = _get_external_coupling_dof(tn_ec, "DDLZ") ? 0 : 1;

      /* Set idfstr with negative value starting from -1 */
      for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) {
        cs_lnum_t ifbr = faces_list[ifac];
        idfstr[ifbr]  = -istruct - 1;
      }
      istruct++;
    }

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the viscosity's type of ALE method
 *
 * parameters:
 *   type --> type of viscosity's type
 *----------------------------------------------------------------------------*/

void
cs_gui_get_ale_viscosity_type(int  *type)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/ale_method");

  *type = _ale_visc_type(tn);
}

/*----------------------------------------------------------------------------
 * Mesh viscosity setting.
 *----------------------------------------------------------------------------*/

void
cs_gui_mesh_viscosity(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  cs_tree_node_t *tn0
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/ale_method");

  /* Get formula */
  const char *mvisc_expr = cs_tree_node_get_child_value_str(tn0, "formula");

  if (mvisc_expr == NULL)
    return;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;
  const cs_lnum_t n_cells   = cs_glob_mesh->n_cells;
  const char * symbols[3]   = {"x", "y", "z" };
  const char * variables[3] = {"mesh_viscosity_1",
                               "mesh_viscosity_2",
                               "mesh_viscosity_3" };

  int orthotropic = _ale_visc_type(tn0);

  cs_lnum_t nd = (orthotropic) ? 3 : 1;

  /* Init mei */
  mei_tree_t *ev = _init_mei_tree(mvisc_expr, variables,
                                  nd, symbols, 0, 3,
                                  cs_glob_time_step->dt_ref,
                                  cs_glob_time_step->t_cur,
                                  cs_glob_time_step->nt_cur);

  /* for each cell, update the value of the table of symbols for each scalar
     (including the thermal scalar), and evaluate the interpreter */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    /* insert symbols */
    mei_tree_insert(ev, "x", cell_cen[c_id][0]);
    mei_tree_insert(ev, "y", cell_cen[c_id][1]);
    mei_tree_insert(ev, "z", cell_cen[c_id][2]);

    mei_evaluate(ev);

    /* Set mesh_u components */
    CS_F_(vism)->val[nd*c_id] = mei_tree_lookup(ev, "mesh_viscosity_1");
    if (orthotropic) {
      CS_F_(vism)->val[nd*c_id + 1] = mei_tree_lookup(ev, "mesh_viscosity_2");
      CS_F_(vism)->val[nd*c_id + 2] = mei_tree_lookup(ev, "mesh_viscosity_3");
    }
  }
  mei_tree_destroy(ev);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Translate the user settings for the domain boundaries into a
 *         structure storing the ALE boundaries (New mechanism used in CDO)
 *
 * \param[in, out]  domain   pointer to a \ref cs_domain_t structure
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_get_boundaries(cs_domain_t     *domain)
{
  assert(domain != NULL);

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  /* Loop on wall boundary zones */

  for (cs_tree_node_t *tn_w = cs_tree_node_get_child(tn_b0, "wall");
       tn_w != NULL;
       tn_w = cs_tree_node_get_next_of_name(tn_w)) {

    const char *label = cs_tree_node_get_tag(tn_w, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
    if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
      continue;

    cs_boundary_type_t  ale_bdy = _get_ale_boundary_type(tn_w);

    if (ale_bdy == CS_BOUNDARY_N_TYPES)
      continue;

    cs_boundary_add(domain->ale_boundaries,
                    ale_bdy,
                    z->name);

  } /* Loop on wall boundary zones */

    /* TODO */
    /* else if (nature == ale_boundary_nature_fixed_displacement) { */
    /*   t0 = cs_timer_wtime(); */
    /*   for (cs_lnum_t ifac = 0; ifac < n_faces; ifac++) { */
    /*     cs_lnum_t ifbr = faces_list[ifac]; */
    /*     _uialcl_fixed_displacement(tn_w, */
    /*                                m->b_face_vtx_idx[ifbr], */
    /*                                m->b_face_vtx_idx[ifbr+1], */
    /*                                m->b_face_vtx_lst, impale, disale, */
    /*                                cs_glob_time_step->dt_ref, */
    /*                                cs_glob_time_step->t_cur, */
    /*                                cs_glob_time_step->nt_cur); */
    /*   } */
    /*   cs_gui_add_mei_time(cs_timer_wtime() - t0); */
    /* } */
    /* else if (nature == ale_boundary_nature_fixed_velocity) { */
    /*   t0 = cs_timer_wtime(); */
    /*   _uialcl_fixed_velocity(tn_w, *iuma, *ivma, *iwma, *ivimpo, */
    /*                          m->n_b_faces, n_faces, faces_list, */
    /*                          ialtyb, rcodcl, */
    /*                          cs_glob_time_step->dt_ref, */
    /*                          cs_glob_time_step->t_cur, */
    /*                          cs_glob_time_step->nt_cur); */
    /*   cs_gui_add_mei_time(cs_timer_wtime() - t0); */
    /* } */

  /* Loop on free surface boundary zones */

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_b0, "free_surface");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name_try(label);
    if (z == NULL)  /* possible in case of old XML file with "dead" nodes */
      continue;

    cs_boundary_add(domain->ale_boundaries,
                    CS_BOUNDARY_ALE_FREE_SURFACE,
                    z->name);

  } /* Loop on free surface boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the fixed velocity for a boundary
 *
 * \param[in]  label boundary condition label
 * \param[out] vel   imposed mesh velocity
 *----------------------------------------------------------------------------*/

void
cs_gui_mobile_mesh_get_fixed_velocity(const char*    label,
                                      cs_real_3_t    vel)
{
  cs_real_t dtref = cs_glob_time_step->dt_ref;
  cs_real_t ttcabs = cs_glob_time_step->t_cur;
  int ntcabs = cs_glob_time_step->nt_cur;

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree, "boundary_conditions");

  /* Loop on wall boundary zones */

  for (cs_tree_node_t *tn_w = cs_tree_node_get_child(tn_b0, "wall");
       tn_w != NULL;
       tn_w = cs_tree_node_get_next_of_name(tn_w)) {

    const char *_node_label = cs_tree_node_get_tag(tn_w, "label");

    if (strcmp(_node_label, label) == 0) {

      const char*  variables[3] = {"mesh_velocity_U",
                                   "mesh_velocity_V",
                                   "mesh_velocity_W" };

      /* Get formula */
      const char *formula = _get_ale_boundary_formula(tn_w, "fixed_velocity");

      if (!formula)
        bft_error(__FILE__, __LINE__, 0,
                  _("Boundary nature formula is null for %s."),
                  cs_gui_node_get_tag(tn_w, "label"));

      /* Init MEI */
      mei_tree_t *ev = _init_mei_tree(formula, variables, 3, 0, 0, 0,
                                      dtref, ttcabs, ntcabs);

      mei_evaluate(ev);

      vel[0] = mei_tree_lookup(ev, "mesh_velocity_U");
      vel[1] = mei_tree_lookup(ev, "mesh_velocity_V");
      vel[2] = mei_tree_lookup(ev, "mesh_velocity_W");

      mei_tree_destroy(ev);

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
