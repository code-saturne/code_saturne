/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CCMIO files
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

/*----------------------------------------------------------------------------*/

#if defined(HAVE_CCM)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * CCM library headers
 *----------------------------------------------------------------------------*/

#include <libccmio/ccmio.h>
#include <libccmio/ccmioversion.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_base.h"
#include "cs_file.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_ccm.h"

/*----------------------------------------------------------------------------
 *  Constants
 *----------------------------------------------------------------------------*/

#define MESH_TIME 0
#define FIELD_TIME 1

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Definitions missing in older CCMIO versions */

#if !defined(CCMIOSIZEC)
  #define CCMIOSIZEC(x)  (x)
#endif

#if !defined(CCMIOINDEXC)
  #define CCMIOINDEXC(x) (x)
#endif

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

#if (kCCMIOVersion == 20601)
typedef int CCMIOSize_t;
typedef int CCMIOIndex_t;
#endif

typedef int cs_ccm_num_t;          /* CCM integer for connectivity */

/*----------------------------------------------------------------------------
 * FVM nodal to writer section translation list
 *----------------------------------------------------------------------------*/

typedef struct _ccm_writer_section_t {

  struct _ccm_writer_section_t  *next;  /* Pointer to next element
                                           in list (NULL at end) */

  cs_lnum_t   n_elts;                   /* number of asociated elements
                                           (or vertices) */

  cs_lnum_t   num_shift;                /* Element number shift when no
                                           parent lists are used */

  const cs_lnum_t    *parent_elt_num;   /* pointer to parent list */

} ccm_writer_section_t;

/*----------------------------------------------------------------------------
 * Time set CCM mesh structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int           n_time_values;   /* Number of time step values */
  int           last_time_step;  /* Last (current) time step number */

  double       *time_value;      /* Time step values */

} fvm_to_ccm_time_t;

/*----------------------------------------------------------------------------
 * CCM writer structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char         *name;                /* Writer name */
  char         *mesh_filename;       /* associated CCMIO geometry file name */
  char         *mesh_basename;       /* base CCMIO geometry file name */
  char         *solution_filename;   /* associated CCMIO solution file name */

  CCMIOID       root_id;             /* Id of the root_node */
  CCMIOID       vertices_id;         /* Id of the vertices node */
  CCMIOID       topology_id;         /* Id of the topology node */
  CCMIOID       state_id;            /* Id of the master state */
  CCMIOID       processor_id;        /* Id of the processor node */
  CCMIOID       solution_id;         /* Id of the solution node */
  CCMIOID       cell_map_id;         /* Id of the cell map */
  CCMIOID       b_face_map_id;       /* Id of the boundary faces map */
  CCMIOID       i_face_map_id;       /* Id of the internal faces map */
  CCMIOID       vtx_map_id;          /* Id of the cell map */

  int           time_step;           /* Current time step number */
  double        time_value;          /* Current time value  */

  bool          is_open;             /* True if CCM file is open */

  int           rank;                /* Local process rank in communicator */
  int           n_ranks;             /* Number of ranks in communicator */

  unsigned long state_counter;       /* Next state number */

  char         *path;                /* Path to the ccmg/ccmp files */

  cs_gnum_t     n_g_perio_faces;     /* Associated number of
                                        periodic faces */

  const fvm_nodal_t  *v_mesh;        /* Reference volume mesh */
  const fvm_nodal_t  *b_mesh;        /* Reference boundary mesh */

  fvm_writer_time_dep_t   time_dependency;  /* Time dependency */
  fvm_to_ccm_time_t       mesh_time;        /* Mesh time structure */
  fvm_to_ccm_time_t       field_time;       /* Field time structure */

  int                   n_time_fields[3]; /* Number of fields for a
                                             given time for cells,
                                             boundary faces, and vertices */

#if defined(HAVE_MPI)
  MPI_Comm     comm;                 /* Associated MPI communicator */
#endif

} fvm_to_ccm_writer_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_datatype_t _ccm_num_datatype = CS_INT32;
static char _ccm_version_string[32] = {"CCM"};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/* prototype for use by following function */

void ADF_Database_Version(const double Root_ID,
                          char *version,
                          char *creation_date,
                          char *modification_date,
                          int *error_return);

/*----------------------------------------------------------------------------
 * Function used only to ensure link with adf library.
 *
 * For shared library builds on versions of Linux recent enough to use
 * the gold linker, linker commands to use the ADF library seem to be ignored,
 * as the libccmio.so library does not include dependency info to libadf
 * (from LibCCMIO versions 2.6.1 to 2.06.023 at least).
 *
 * parameters:
 *   do_something <-- if true, call ADF function (should be called with false)
 *----------------------------------------------------------------------------*/

static void
_force_adf_link(bool do_something)
{
  if (do_something) {
    char *version = NULL, *creation_date = NULL, *modification_date = NULL;
    int error_return;
    double root_id = 0;

    ADF_Database_Version(root_id,
                         version,
                         creation_date,
                         modification_date,
                         &error_return);
  }
}

/*----------------------------------------------------------------------------
 * Get the global number of entities associated to a mesh.
 *
 * This function assumes non-duplicated entities, such as cells or boundary
 * faces.
 *
 * parameters:
 *   mesh     <-- pointer to nodal mesh structure
 *   elt_dim  <-- dimension of the entities to consider
 *
 * returns:
 *   number of associated entities
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_n_g_mesh_elts(const fvm_nodal_t  *mesh,
               int                 ent_dim)
{
  int i;
  cs_gnum_t retval = 0;

  for (i = 0; i < mesh->n_sections; i++) {
    const fvm_nodal_section_t  *const  section = mesh->sections[i];
    if (section->entity_dim == ent_dim) {
      if (section->global_element_num != NULL)
        retval += fvm_io_num_get_global_count(section->global_element_num);
      else {
        retval += section->n_elements;
      }
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Build cell global numbering array in order of output (as defined by
 * nodal mesh sections)
 *
 * parameters:
 *   mesh     <-- pointer to nodal mesh structure
 *   elt_dim  <-- dimension of the entities to consider
 *   elt_gnum --> associated global numbering array
 *----------------------------------------------------------------------------*/

static void
_build_ordered_elt_gnum(const fvm_nodal_t  *mesh,
                        int                 ent_dim,
                        cs_gnum_t          *elt_gnum)
{
  int i;
  cs_lnum_t j;
  cs_gnum_t num_shift = 0;

  for (i = 0; i < mesh->n_sections; i++) {
    const fvm_nodal_section_t  *const  section = mesh->sections[i];
    if (section->entity_dim == ent_dim) {
      if (section->global_element_num != NULL) {
        const cs_gnum_t *g_num
          = fvm_io_num_get_global_num(section->global_element_num);
        if (section->parent_element_num != NULL) {
          const cs_lnum_t *p_num = section->parent_element_num;
          for (j = 0; j < section->n_elements; j++)
            elt_gnum[p_num[j]-1] = g_num[j] + num_shift;
        }
        else {
          for (j = 0; j < section->n_elements; j++)
            elt_gnum[j] = g_num[j] + num_shift;
        }
        num_shift += fvm_io_num_get_global_count(section->global_element_num);
      }
      else {
        if (section->parent_element_num != NULL) {
          const cs_lnum_t *p_num = section->parent_element_num;
          for (j = 0; j < section->n_elements; j++)
            elt_gnum[p_num[j]-1] = j+1 + num_shift;
        }
        else {
          for (j = 0; j < section->n_elements; j++)
            elt_gnum[j] = j+1 + num_shift;
        }
        num_shift += section->n_elements;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Build cell global numbering array in order of output (as defined by
 * nodal mesh sections)
 *
 * parameters:
 *   b_mesh           <-- pointer to base mesh structure
 *   mesh             <-- pointer to nodal mesh structure
 *
 * returns:
 *   array of global cell numbers
 *----------------------------------------------------------------------------*/

static cs_gnum_t *
_build_ordered_cell_gnum(const cs_mesh_t    *b_mesh,
                         const fvm_nodal_t  *mesh)
{
  cs_gnum_t *cell_gnum = NULL;

  /* Allocate array */

  BFT_MALLOC(cell_gnum, b_mesh->n_cells_with_ghosts, cs_gnum_t);

  /* Build global numbering */

  _build_ordered_elt_gnum(mesh, 3, cell_gnum);

  /* Synchronize halo, blanking periodicity */

  if (b_mesh->halo != NULL) {

    const cs_halo_t *halo = b_mesh->halo;

    cs_lnum_t i;
    cs_lnum_t  rank_id, t_id, shift;
    cs_lnum_t  start = 0, end = 0;

    const cs_int_t  n_transforms = halo->n_transforms;
    const cs_int_t  n_elts = halo->n_local_elts;

    cs_halo_sync_untyped(b_mesh->halo,
                         CS_HALO_EXTENDED,
                         sizeof(cs_gnum_t),
                         cell_gnum);

    for (t_id = 0; t_id < n_transforms; t_id++) {

      shift = 4 * halo->n_c_domains * t_id;

      for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

        start = halo->perio_lst[shift + 4*rank_id];
        end = start + halo->perio_lst[shift + 4*rank_id + 1];
        for (i = start; i < end; i++)
          cell_gnum[n_elts+i] = 0;

        start = halo->perio_lst[shift + 4*rank_id + 2];
        end = start + halo->perio_lst[shift + 4*rank_id + 3];
        for (i = start; i < end; i++)
          cell_gnum[n_elts+i] = 0;

      } /* End of loop on ranks */

    } /* End of loop on transformation */

  }

  return cell_gnum;
}

/*----------------------------------------------------------------------------
 * Build face global numbering array in order of output (as defined by
 * nodal mesh sections)
 *
 * parameters:
 *   b_mesh <-- pointer to base mesh structure
 *
 * returns:
 *   array of global boundary face numbers
 *----------------------------------------------------------------------------*/

static cs_gnum_t *
_build_ordered_b_face_gnum(const cs_mesh_t  *b_mesh)
{
  cs_gnum_t *face_gnum = NULL;

  /* Allocate array */

  BFT_MALLOC(face_gnum, b_mesh->n_b_faces, cs_gnum_t);

  /* As the nodal mesh used to build the CCM mesh is the (full) volume
     mesh, we need to build a local nodal boundary mesh to ensure
     numberings are consistent */

  fvm_nodal_t *mesh = cs_mesh_connect_faces_to_nodal(b_mesh,
                                                     NULL,
                                                     false,
                                                     0,
                                                     b_mesh->n_b_faces,
                                                     NULL,
                                                     NULL);
  fvm_nodal_reduce(mesh, 0);

  /* Build global numbering */

  _build_ordered_elt_gnum(mesh, 2, face_gnum);

  fvm_nodal_destroy(mesh);

  return face_gnum;
}

/*----------------------------------------------------------------------------
 * Build an array of element's ordering based on their global numbers.
 *
 * parameters:
 *   n_elts   <-- number of elements
 *   elt_gnum <-- array of associated global numbers (or NULL)
 *
 * returns:
 *   element ordering array
 *----------------------------------------------------------------------------*/

static cs_lnum_t *
_build_order_by_gnum(cs_lnum_t         n_elts,
                     const cs_gnum_t  *elt_gnum)
{
  cs_lnum_t i;
  cs_lnum_t *order = NULL;

  /* Allocate array */

  BFT_MALLOC(order, n_elts, cs_lnum_t);

  /* Build global numbering */

  if (elt_gnum != NULL) {
    for (i = 0; i < n_elts; i++)
      order[elt_gnum[i] - 1] = i;
  }
  else {
    for (i = 0; i < n_elts; i++)
      order[i] = i;
  }

  return order;
}

/*----------------------------------------------------------------------------
 * Build list of sections to output
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   export_dim       <-- minimum dimension of sections to export
 *
 * returns:
 *   array of section translations (must be freed by caller),
 *   or NULL if section list is completely empty
 *----------------------------------------------------------------------------*/

static ccm_writer_section_t *
_build_export_list(const fvm_nodal_t  *mesh,
                   int                 export_dim)
{
  int  i;
  int  n_sections = 0;
  ccm_writer_section_t *export_list = NULL;

  cs_lnum_t num_shift = 0;

  /* Initial count and allocation */

  n_sections = 0;

  if (export_dim == 0)
    n_sections = 1;
  else if (export_dim > 1) {
    for (i = 0; i < mesh->n_sections; i++) {
      const fvm_nodal_section_t  *const  section = mesh->sections[i];
      if (section->entity_dim == export_dim)
        n_sections += 1;
    }
  }

  /* If no sections are present no list is returned */

  if (n_sections == 0)
    return NULL;

  BFT_MALLOC(export_list, n_sections, ccm_writer_section_t);

  /* Build unsorted list */

  if (export_dim == 0) {
    (export_list[0]).n_elts = mesh->n_vertices;
    (export_list[0]).num_shift = 0;
    (export_list[0]).parent_elt_num = mesh->parent_vertex_num;
  }
  else if (export_dim > 1) {
    n_sections = 0;
    for (i = 0; i < mesh->n_sections; i++) {
      const fvm_nodal_section_t  *const  section = mesh->sections[i];
      if (section->entity_dim != export_dim)
        continue;
      (export_list[n_sections]).n_elts = section->n_elements;
      (export_list[n_sections]).num_shift = num_shift;
      (export_list[n_sections]).parent_elt_num = section->parent_element_num;
      n_sections++;
      num_shift += section->n_elements;
    }
  }

  for (i = 0; i < n_sections - 1; i++)
    (export_list[i]).next = &(export_list[i+1]);
  export_list[n_sections - 1].next = NULL;

  return export_list;
}

/*----------------------------------------------------------------------------
 * Create and write a state.
 *
 * parameters:
 *   w  <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_state(fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    CCMIOID state_id;
    char *state_full_name;
    char state_number[20];
    const char *state_name = "State ";

    /* Prepare state name */
    sprintf(state_number, "%lu", w->state_counter);
    BFT_MALLOC(state_full_name,
               strlen(state_name) + strlen(state_number) + 1,
               char);
    strcpy(state_full_name, state_name);
    strcat(state_full_name, state_number);

    /* Write_state */
    CCMIOError error = kCCMIONoErr, *err = &error;
    CCMIONewState(err, w->root_id, state_full_name, NULL, NULL, &state_id);
    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing state."), (int)error);
    w->state_id = state_id;

    BFT_FREE(state_full_name);
  }

  w->state_counter++;
}

/*----------------------------------------------------------------------------
 * Create and write a processor node; clean the node if it is not new.
 *
 * parameters:
 *   w <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_processor(fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {
    CCMIOSize_t i = 0;
    CCMIOError error = kCCMIONoErr, *err = &error;
    CCMIOID processor_id;

    /* Check if the current state already has a Processor node */
    if (   CCMIONextEntity(NULL, w->state_id, kCCMIOProcessor, &i, &processor_id)
        != kCCMIONoErr)
      CCMIONewEntity(err, w->state_id, kCCMIOProcessor, NULL, &processor_id);

    /* Clear the node data in any case */
    CCMIOClearProcessor(err,
                        w->state_id,
                        processor_id,
                        TRUE,
                        TRUE,
                        TRUE,
                        TRUE,
                        TRUE);

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing processor node."), (int)error);

    w->processor_id = processor_id;
  }
}

/*----------------------------------------------------------------------------
 * Finalize a processor node.
 *
 * parameters:
 *   vertices_path  <-> vertices path in file structure
 *   topology_path  <-> topology path in file structure
 *   w              <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_finalize_processor(char                 *vertices_path,
                    char                 *topology_path,
                    fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {
    CCMIOError error = kCCMIONoErr, *err = &error;
    CCMIOWriteProcessor(err,
                        w->processor_id,
                        vertices_path,
                        &w->vertices_id,
                        topology_path,
                        &w->topology_id,
                        NULL,
                        NULL,
                        NULL,
                        &w->solution_id);
    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d finalizing processor node."),
                (int)error);
  }
}

/*----------------------------------------------------------------------------
 * Update time, either for field data or mesh.
 *
 * parameters:
 *   type        <-- type of the time structure {MESH_TIME, FIELD_TIME}
 *   time        <-> pointer to time structure
 *   time_step   <-- current time step
 *   time_value  <-- current time value
 *   w           <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_update_time(int                   type,
             fvm_to_ccm_time_t    *time,
             int                   time_step,
             double                time_value,
             fvm_to_ccm_writer_t  *w)
{
  if (type == MESH_TIME) {

    /* If fixed mesh */
    if (time_step == -1 && w->state_counter == 1) {
      time->n_time_values = 1;
      time->last_time_step = -1;
      BFT_MALLOC(time->time_value, time->n_time_values, double);
      time->time_value[time->n_time_values-1] = time_value;
    }

    /* If variable mesh */
    else if (time_step != -1) {
      time->n_time_values++;
      time->last_time_step = time_step;
      BFT_REALLOC(time->time_value, time->n_time_values, double);
      time->time_value[time->n_time_values-1] = time_value;
    }

  }

  /* Field time */
  else if (type == FIELD_TIME) {
    time->n_time_values++;
    time->last_time_step = time_step;
    BFT_REALLOC(time->time_value, time->n_time_values, double);
    time->time_value[time->n_time_values-1] = time_value;
  }

}

/*----------------------------------------------------------------------------
 *  Create and write a phase if the node does not exist
 *
 *  parameters:
 *    phase_id     --> id of the face that is written
 *    w            <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_phase(CCMIOID              *phase_id,
             fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {
    CCMIOSize_t i = 0;
    CCMIOError error = kCCMIONoErr, *err =  &error;

    /* Check if the current solution node already has a phase node */
    if (   CCMIONextEntity(NULL, w->solution_id, kCCMIOFieldPhase, &i, phase_id)
        != kCCMIONoErr) {
      CCMIONewIndexedEntity(err,
                            w->solution_id,
                            kCCMIOFieldPhase,
                            0,
                            NULL,
                            phase_id);
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d creating new phase."),
                  (int)error);

    }
  }
}

/*----------------------------------------------------------------------------
 * Create and write a problem description if none exists.
 *
 * parameters:
 *   w <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_problem_description(fvm_to_ccm_writer_t  *w)
{
  CCMIOID problem_id, id;
  CCMIOError error = kCCMIONoErr, *err =  &error;
  CCMIOSize_t i = 0;

  if (w->rank < 1) {

    if (CCMIONextEntity(NULL,
                        w->root_id,
                        kCCMIOProblemDescription,
                        &i,
                        &problem_id) != kCCMIONoErr) {
      CCMIONewEntity(err,
                     w->root_id,
                     kCCMIOProblemDescription,
                     NULL,
                     &problem_id);
      CCMIONewIndexedEntity(err,
                            problem_id,
                            kCCMIOCellType,
                            1,
                            "Region_1",
                            &id);

      /* TODO write additional regions here */
    }
    CCMIOWriteState(err, w->state_id, problem_id, NULL);

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing problem description."),
                (int)error);

  }
}

/*----------------------------------------------------------------------------
 * Create and write a solution.
 *
 * parameters:
 *   w <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_solution(fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    CCMIOError error = kCCMIONoErr, *err = &error;
    CCMIOID solution_id;

    CCMIONewEntity(err, w->root_id, kCCMIOFieldSet, "Field set", &solution_id);
    w->solution_id = solution_id;

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing solution."), (int)error);

  }
}

/*----------------------------------------------------------------------------
 * Write restart info.
 *
 * parameters:
 *   time_step     <-- time step to write
 *   start_angle   <-- beginning start angle
 *   time_value    <-- time value to write
 *   w             <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_restart_info(int                   time_step,
                    double                time_value,
                    double                start_angle,
                    fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    /* Prepare solver name */
    char solver_info[128];
    sprintf(solver_info,
            "Code_Saturne "VERSION" with libCCMIO %d",
            kCCMIOVersion);

    /* Write node */
    CCMIOError error = kCCMIONoErr, *err = &error;
    CCMIOID restart_id;
    CCMIONewEntity(err, w->solution_id, kCCMIORestart, NULL, &restart_id);
    CCMIOWriteRestartInfo(err,
                          restart_id,
                          solver_info,
                          time_step,
                          time_value,
                          NULL,
                          start_angle);

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing restart info."), (int)error);

  }

}

/*----------------------------------------------------------------------------
 * Create and write a map.
 *
 * parameters:
 *   name          <-- optional map name, or NULL
 *   n_g_elts      <-- global number of elements associated to this entity
 *   bi            <-- part to block info structure fo this entity
 *   map_num_shift <-- shift associated with this map
 *   map_id        --> CCMIO id for map (for rank 0 only)
 *   w             <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_map(const char             *name,
           cs_ccm_num_t            n_g_elts,
           cs_block_dist_info_t    bi,
           cs_ccm_num_t            map_num_shift,
           CCMIOID                *map_id,
           fvm_to_ccm_writer_t    *w)
{
  if (w->rank < 1) {

    cs_ccm_num_t start_id, end_id, i, j;
    cs_ccm_num_t *map_data = NULL;
    CCMIOError error = kCCMIONoErr, *err = &error;

    CCMIONewEntity(err, w->root_id, kCCMIOMap, name, map_id);

    BFT_MALLOC(map_data, bi.block_size, cs_ccm_num_t);

    for (start_id = 1, end_id = bi.block_size;
         start_id < n_g_elts;
         start_id += bi.block_size, end_id += bi.block_size) {

      CCMIOIndex_t start, end;

      if (start_id > n_g_elts)
        start_id = n_g_elts;
      if (end_id > n_g_elts)
        end_id = n_g_elts;

      for (i = start_id, j = 0; i <= end_id; i++, j++)
        map_data[j] = i + map_num_shift;

      /* The index of the first table of data written should be kCCMIOStart */

      if (start_id == 1)
        start = kCCMIOStart;
      else
        start = CCMIOINDEXC(start_id-1);
      end = CCMIOINDEXC(end_id);
      CCMIOWriteMap(err, *map_id, CCMIOSIZEC(n_g_elts), CCMIOSIZEC(n_g_elts),
                    map_data, start, end);

      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing map."),
                  (int)error);

    }
    BFT_FREE(map_data);
  }
}

/*----------------------------------------------------------------------------
 * Write a cell map for the ccmp file to be able to write the fields.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   w        <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_cells_map(const cs_mesh_t      *mesh,
                 fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    CCMIOID map_id;

    cs_block_dist_info_t cell_bi
      = cs_block_dist_compute_sizes(w->rank,
                                    w->n_ranks,
                                    0,
                                    cs_parall_get_min_coll_buf_size(),
                                    mesh->n_g_cells);

    /* Write map and define entity */

    _write_map("Cell map",
               mesh->n_g_cells,
               cell_bi,
               0,
               &map_id,
               w);

    w->cell_map_id = map_id;

  }
}

/*----------------------------------------------------------------------------
 * Write a boundary faces map for the ccmp file to be able to write the fields.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   w        <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_faces_map(const cs_mesh_t      *mesh,
                 fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    CCMIOID map_id;
    cs_ccm_num_t  map_num_shift = 0;
    cs_ccm_num_t n_g_faces = mesh->n_g_b_faces + w->n_g_perio_faces;

    cs_block_dist_info_t face_bi
      = cs_block_dist_compute_sizes(w->rank,
                                    w->n_ranks,
                                    0,
                                    cs_parall_get_min_coll_buf_size(),
                                    n_g_faces);

    _write_map(NULL, n_g_faces, face_bi, map_num_shift, &map_id, w);

    w->b_face_map_id = map_id;

  }
}

/*----------------------------------------------------------------------------
 * Write a vertices map for the ccmp file to be able to write the fields.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   w        <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_vertices_map(const cs_mesh_t      *mesh,
                    fvm_to_ccm_writer_t  *w)
{
  if (w->rank < 1) {

    CCMIOID map_id;
    cs_ccm_num_t  map_num_shift = 0;
    cs_ccm_num_t n_g_vertices = mesh->n_g_vertices;

    cs_block_dist_info_t vtx_bi
      = cs_block_dist_compute_sizes(w->rank,
                                    w->n_ranks,
                                    0,
                                    cs_parall_get_min_coll_buf_size(),
                                    n_g_vertices);

    _write_map("Vertex map", n_g_vertices, vtx_bi, map_num_shift, &map_id, w);

    w->vtx_map_id = map_id;

  }
}

/*----------------------------------------------------------------------------
 * Count global periodic faces.
 *
 * parameters:
 *   b_mesh        <-- pointer to base mesh structure
 *   cell_gnum     <-- array of global cell numbers, ordered by nodal mesh
 *   w             <-> pointer to writer structure
 *
 * returns:
 *   global number of periodic faces
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_count_faces_perio_g(const cs_mesh_t      *b_mesh,
                     const cs_gnum_t      *cell_gnum,
                     fvm_to_ccm_writer_t  *w)
{
  cs_lnum_t i;
  cs_gnum_t n_g_perio_faces = 0;

  if (b_mesh->periodicity != NULL) {

    const cs_lnum_2_t *face_cells = (const cs_lnum_2_t *)(b_mesh->i_face_cells);

    for (i = 0; i < b_mesh->n_i_faces; i++) {
      if (   cell_gnum[face_cells[i][0]] == 0
          || cell_gnum[face_cells[i][1]] == 0)
        n_g_perio_faces += 1;
    }

    cs_parall_sum(1, CS_GNUM_TYPE, &n_g_perio_faces);

  }

  w->n_g_perio_faces = n_g_perio_faces;

  return n_g_perio_faces;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write vertices in parallel.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   w        <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_vertices_g(const cs_mesh_t      *mesh,
                  fvm_to_ccm_writer_t  *w)
{
  cs_part_to_block_t *d = NULL;
  cs_file_serializer_t *s = NULL;

  void *_vtx_coords = NULL, *_vtx_coords_s = NULL;

  const cs_datatype_t real_type
    = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  CCMIOError error = kCCMIONoErr, *err = &error;

  cs_block_dist_info_t vtx_bi
    = cs_block_dist_compute_sizes(w->rank,
                                  w->n_ranks,
                                  0,
                                  cs_parall_get_min_coll_buf_size(),
                                  mesh->n_g_vertices);

  /* Write map and define entity */

  _write_map("Vertex map",
             mesh->n_g_vertices,
             vtx_bi,
             0,
             &(w->vtx_map_id),
             w);

  if (w->rank < 1) {
    CCMIONewEntity(err,
                   w->root_id,
                   kCCMIOVertices,
                   "Vertices",
                   &(w->vertices_id));
    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing new vertices entity."), (int)error);
  }

  /* Create distribution structure */

  d = cs_part_to_block_create_by_gnum(w->comm,
                                      vtx_bi,
                                      mesh->n_vertices,
                                      mesh->global_vtx_num);

  BFT_MALLOC(_vtx_coords,
             (vtx_bi.gnum_range[1] - vtx_bi.gnum_range[0]) * 3 *sizeof(cs_real_t),
             cs_byte_t);

  cs_part_to_block_copy_array(d,
                              real_type,
                              3,
                              mesh->vtx_coord,
                              _vtx_coords);

  /* Now write vertex coordinates */

  cs_gnum_t range[2] = {vtx_bi.gnum_range[0],
                        vtx_bi.gnum_range[1]};

  s = cs_file_serializer_create(sizeof(cs_real_t),
                                3,
                                range[0],
                                range[1],
                                0,
                                _vtx_coords,
                                w->comm);

  do {

    _vtx_coords_s = cs_file_serializer_advance(s, range);

    if (_vtx_coords_s != NULL) { /* only possible on rank 0 */
      double scale = 1.0;
      if (sizeof(cs_real_t) == 8)
        CCMIOWriteVerticesd(err,
                            w->vertices_id,
                            CCMIOSIZEC(3),
                            scale,
                            w->vtx_map_id,
                            (const double *)_vtx_coords_s,
                            CCMIOINDEXC(range[0]-1),
                            CCMIOINDEXC(range[1]-1));
      else if (sizeof(cs_real_t) == 4)
        CCMIOWriteVerticesf(err,
                            w->vertices_id,
                            CCMIOSIZEC(3),
                            scale,
                            w->vtx_map_id,
                            (const float *)_vtx_coords_s,
                            CCMIOINDEXC(range[0]-1),
                            CCMIOINDEXC(range[1]-1));
    }

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing vertices."), (int)error);

  } while (_vtx_coords_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_vtx_coords);

  cs_part_to_block_destroy(&d);
}

/*----------------------------------------------------------------------------
 * Write cells in parallel.
 *
 * parameters:
 *   b_mesh    <-- pointer to base mesh structure
 *   cell_gnum <-- pointer to global cell number, ordered by nodal mesh
 *   w         <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_cells_g(const cs_mesh_t      *b_mesh,
               const cs_gnum_t      *cell_gnum,
               fvm_to_ccm_writer_t  *w)
{
  cs_part_to_block_t *d = NULL;
  cs_file_serializer_t *s = NULL;

  int *_cell_gc_id = NULL, *_cell_gc_id_s = NULL;

  CCMIOID map_id, topology_id, cells_id;
  CCMIOError error = kCCMIONoErr, *err = &error;

  const cs_datatype_t lnum_type
    = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  cs_block_dist_info_t cell_bi
    = cs_block_dist_compute_sizes(w->rank,
                                  w->n_ranks,
                                  0,
                                  cs_parall_get_min_coll_buf_size(),
                                  b_mesh->n_g_cells);

  /* Write map and define entity */

  _write_map("Cell map",
             b_mesh->n_g_cells,
             cell_bi,
             0,
             &map_id,
             w);

  w->cell_map_id = map_id;

  if (w->rank < 1) {
    CCMIONewEntity(err, w->root_id, kCCMIOTopology, "Topology", &topology_id);
    CCMIONewEntity(err, topology_id, kCCMIOCells, "Cells", &cells_id);
    w->topology_id = topology_id;
  }

  /* Create distribution structure */

  BFT_MALLOC(_cell_gc_id,
             (cell_bi.gnum_range[1] - cell_bi.gnum_range[0]),
             cs_lnum_t);

  d = cs_part_to_block_create_by_gnum(w->comm,
                                      cell_bi,
                                      b_mesh->n_cells,
                                      cell_gnum);

  cs_part_to_block_copy_array(d,
                              lnum_type,
                              1,
                              b_mesh->cell_family,
                              _cell_gc_id);

  /* Now write cell family info */

  cs_gnum_t range[2] = {cell_bi.gnum_range[0],
                        cell_bi.gnum_range[1]};

  s = cs_file_serializer_create(sizeof(cs_ccm_num_t),
                                1,
                                range[0],
                                range[1],
                                0,
                                _cell_gc_id,
                                w->comm);

  do {

    _cell_gc_id_s = cs_file_serializer_advance(s, range);

    if (_cell_gc_id_s != NULL) { /* only possible on rank 0 */

      CCMIOWriteCells(err, cells_id, map_id, _cell_gc_id_s,
                      CCMIOINDEXC(range[0]-1), CCMIOINDEXC(range[1]-1));

      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing cells."), (int)error);

    }

  } while (_cell_gc_id_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_cell_gc_id);

  cs_part_to_block_destroy(&d);
}

/*----------------------------------------------------------------------------
 * Write face -> vertices connectivity in parallel.
 *
 * parameters:
 *   b_mesh        <-- pointer to base mesh structure
 *   face_bi       <-- face part to block info structure
 *   entity        <-- interior or boundary faces
 *   entity_id     <-- CCMIO id for this entity
 *   map_id        <-- CCMIO id for map
 *   d             <-- face part to block distribution helper
 *   w             <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_g(const cs_mesh_t         *b_mesh,
                       cs_block_dist_info_t     face_bi,
                       CCMIOEntity              entity,
                       CCMIOID                  entity_id,
                       CCMIOID                  map_id,
                       cs_part_to_block_t      *d,
                       fvm_to_ccm_writer_t     *w)
{
  cs_lnum_t i, j, k;
  cs_ccm_num_t n_face_vertices;

  cs_lnum_t *face_vtx_idx = NULL, *face_vtx_lst = NULL;
  cs_ccm_num_t *face_connect_idx = NULL, *_face_connect_idx = NULL;
  cs_ccm_num_t *face_connect_g = NULL, *_face_connect_g = NULL;

  cs_lnum_t n_faces = 0, face_connect_size = 0;
  cs_ccm_num_t block_size = 0;

  cs_ccm_num_t *_face_connect_g_s = NULL;

  cs_file_serializer_t *s = NULL;

  const cs_datatype_t ccm_num_type
    = (sizeof(cs_ccm_num_t) == 8) ? CS_INT64 : CS_INT32;

  if (entity == kCCMIOInternalFaces) {
    n_faces = b_mesh->n_i_faces;
    face_connect_size = b_mesh->i_face_vtx_connect_size;
    face_vtx_idx = b_mesh->i_face_vtx_idx;
    face_vtx_lst = b_mesh->i_face_vtx_lst;
  }
  else if (entity == kCCMIOBoundaryFaces) {
    n_faces = b_mesh->n_b_faces;
    face_connect_size = b_mesh->b_face_vtx_connect_size;
    face_vtx_idx = b_mesh->b_face_vtx_idx;
    face_vtx_lst = b_mesh->b_face_vtx_lst;
  }

  block_size = face_bi.gnum_range[1] - face_bi.gnum_range[0];

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  BFT_MALLOC(face_connect_idx, n_faces + 1, cs_ccm_num_t);
  BFT_MALLOC(_face_connect_idx, block_size + 1, cs_ccm_num_t);

  face_connect_idx[0] = 0;
  for (i = 0; i < n_faces; i++) {
    n_face_vertices = face_vtx_idx[i+1] - face_vtx_idx[i];
    face_connect_idx[i+1] = face_connect_idx[i] + n_face_vertices + 1;
  }

  cs_part_to_block_copy_index(d, face_connect_idx, _face_connect_idx);

  /* Build connectivity */

  BFT_MALLOC(face_connect_g,
             n_faces + face_connect_size,
             cs_ccm_num_t);

  k = 0;
  for (i = 0; i < n_faces; i++) {
    face_connect_g[k++] = face_vtx_idx[i+1] - face_vtx_idx[i];
    for (j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++)
      face_connect_g[k++] = b_mesh->global_vtx_num[face_vtx_lst[j]];
  }

  BFT_MALLOC(_face_connect_g, _face_connect_idx[block_size], cs_ccm_num_t);

  cs_part_to_block_copy_indexed(d,
                                ccm_num_type,
                                face_connect_idx,
                                face_connect_g,
                                _face_connect_idx,
                                _face_connect_g);

  BFT_FREE(face_connect_g);
  BFT_FREE(face_connect_idx);

  /* Now write face -> vertices connectivity */

  cs_gnum_t g_connect_size = 0;
  cs_gnum_t range_base = _face_connect_idx[block_size], range_shift = 0;
  MPI_Scan(&range_base, &range_shift, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
  range_shift -= _face_connect_idx[block_size];

  cs_gnum_t range[2] = {range_shift + _face_connect_idx[0] + 1,
                        range_shift + _face_connect_idx[block_size] + 1};

  MPI_Allreduce(&(range[1]), &g_connect_size, 1, CS_MPI_GNUM, MPI_MAX, w->comm);
  if (g_connect_size > 0)
    g_connect_size -= 1;

  BFT_FREE(_face_connect_idx);

  s = cs_file_serializer_create(sizeof(ccm_num_type),
                                1,
                                range[0],
                                range[1],
                                0,
                                _face_connect_g,
                                w->comm);

  do {

    _face_connect_g_s = cs_file_serializer_advance(s, range);

    if (_face_connect_g_s != NULL) { /* only possible on rank 0 */
      CCMIOError error = kCCMIONoErr, *err = &error;
      CCMIOWriteFaces(err, entity_id, entity, map_id,
                      CCMIOSIZEC(g_connect_size), _face_connect_g_s,
                      CCMIOINDEXC(range[0]-1), CCMIOINDEXC(range[1]-1));
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing face -> vertices connectivity."),
                  (int)error);
    }

  } while (_face_connect_g_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_face_connect_g);
}

/*----------------------------------------------------------------------------
 * Write face -> cells connectivity in parallel.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   face_bi     <-- face part to block info structure
 *   entity      <-- interior or boundary faces
 *   entity_id   <-- CCMIO id for this entity
 *   map_id      <-- CCMIO id for map
 *   cell_gnum   <-- array of global cell numbers, ordered by nodal mesh
 *   d           <-- face part to block distribution helper
 *   w           <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_face_cells_g(const cs_mesh_t        *b_mesh,
                    cs_block_dist_info_t    face_bi,
                    CCMIOEntity             entity,
                    CCMIOID                 entity_id,
                    CCMIOID                 map_id,
                    const cs_gnum_t        *cell_gnum,
                    cs_part_to_block_t     *d,
                    fvm_to_ccm_writer_t    *w)
{
  cs_lnum_t i;

  cs_lnum_t n_cells_per_face = 0;
  cs_lnum_t n_faces = 0;

  const cs_lnum_t *face_cells = NULL;
  cs_ccm_num_t *face_cell_g = NULL, *_face_cell_g = NULL;
  cs_ccm_num_t *_face_cell_g_s = NULL;
  cs_file_serializer_t *s = NULL;

  const cs_datatype_t ccm_num_type
    = (sizeof(cs_ccm_num_t) == 8) ? CS_INT64 : CS_INT32;

  if (entity == kCCMIOInternalFaces) {
    n_cells_per_face = 2;
    n_faces = b_mesh->n_i_faces;
    face_cells = (const cs_lnum_t *)(b_mesh->i_face_cells);
  }
  else if (entity == kCCMIOBoundaryFaces) {
    n_cells_per_face = 1;
    n_faces = b_mesh->n_b_faces;
    face_cells = b_mesh->b_face_cells;
  }

  /* Face -> cell connectivity */
  /*---------------------------*/

  BFT_MALLOC(face_cell_g, n_faces*n_cells_per_face, cs_ccm_num_t);

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < n_faces; i++) {
      face_cell_g[i*2]     = cell_gnum[face_cells[i*2]];
      face_cell_g[i*2 + 1] = cell_gnum[face_cells[i*2 + 1]];
    }
  }
  else {
    for (i = 0; i < n_faces; i++)
      face_cell_g[i] = cell_gnum[face_cells[i]];
  }

  /* Distribute to blocks and write */

  BFT_MALLOC(_face_cell_g,
             (face_bi.gnum_range[1] - face_bi.gnum_range[0]) * n_cells_per_face,
             cs_ccm_num_t);

  cs_part_to_block_copy_array(d,
                              ccm_num_type,
                              n_cells_per_face,
                              face_cell_g,
                              _face_cell_g);

  BFT_FREE(face_cell_g);

  /* Now write face -> cells connectivity */

  cs_gnum_t range[2] = {face_bi.gnum_range[0],
                        face_bi.gnum_range[1]};

  s = cs_file_serializer_create(sizeof(ccm_num_type),
                                n_cells_per_face,
                                range[0],
                                range[1],
                                0,
                                _face_cell_g,
                                w->comm);

  do {

    _face_cell_g_s = cs_file_serializer_advance(s, range);

    if (_face_cell_g_s != NULL) { /* only possible on rank 0 */
      CCMIOError error = kCCMIONoErr, *err = &error;
      CCMIOWriteFaceCells(err, entity_id, entity, map_id,
                          _face_cell_g_s,
                          CCMIOINDEXC(range[0]-1), CCMIOINDEXC(range[1]-1));
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing face -> cells connectivity."),
                  (int)error);
    }

  } while (_face_cell_g_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_face_cell_g);
}

/*----------------------------------------------------------------------------
 * Write face -> vertices connectivity in parallel with periodic faces.
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   face_bi       <-- face part to block info structure
 *   entity        <-- interior or boundary faces
 *   entity_id     <-- CCMIO id for this entity
 *   map_id        <-- CCMIO id for map
 *   cell_gnum     <-- array of global cell numbers, ordered by nodal mesh
 *   d             <-- face part to block distribution helper
 *   w             <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_perio_g(const cs_mesh_t        *b_mesh,
                             cs_block_dist_info_t    face_bi,
                             CCMIOEntity             entity,
                             CCMIOID                 entity_id,
                             CCMIOID                 map_id,
                             const cs_gnum_t        *cell_gnum,
                             cs_part_to_block_t     *d,
                             fvm_to_ccm_writer_t    *w)
{
  cs_lnum_t i, j, k;
  cs_ccm_num_t n_face_vertices;

  cs_ccm_num_t *face_connect_idx = NULL, *_face_connect_idx = NULL;
  cs_ccm_num_t *face_connect_g = NULL, *_face_connect_g = NULL;

  cs_lnum_t n_faces = b_mesh->n_i_faces;
  cs_lnum_t face_connect_size = b_mesh->i_face_vtx_connect_size;
  cs_ccm_num_t block_size = 0;

  cs_ccm_num_t *_face_connect_g_s = NULL;

  cs_file_serializer_t *s = NULL;

  const cs_datatype_t ccm_num_type
    = (sizeof(cs_ccm_num_t) == 8) ? CS_INT64 : CS_INT32;

  const cs_lnum_t *face_vtx_idx = b_mesh->i_face_vtx_idx;
  const cs_lnum_t *face_vtx_lst = b_mesh->i_face_vtx_lst;
  const cs_lnum_t *face_cells = (const cs_lnum_t *)(b_mesh->i_face_cells);

  /* Allocate arrays large enough for both periodic boundary + true interior
     faces to avoid counting loop */

  BFT_MALLOC(face_connect_idx, n_faces + 1, cs_ccm_num_t);

  block_size = face_bi.gnum_range[1] - face_bi.gnum_range[0];

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  face_connect_idx[0] = 0;

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < n_faces; i++) {
      if (   cell_gnum[face_cells[2*i]] > 0
          && cell_gnum[face_cells[2*i + 1]] > 0) {
        n_face_vertices = face_vtx_idx[i+1] - face_vtx_idx[i];
        face_connect_idx[i+1] = face_connect_idx[i] + n_face_vertices + 1;
      }
      else
        face_connect_idx[i+1] = face_connect_idx[i];
    }
  }
  else if (entity == kCCMIOBoundaryFaces) {
    for (i = 0; i < n_faces; i++) {
      if (   cell_gnum[face_cells[2*i]] == 0
          || cell_gnum[face_cells[2*i + 1]] == 0) {
        n_face_vertices = face_vtx_idx[i+1] - face_vtx_idx[i];
        face_connect_idx[i+1] = face_connect_idx[i] + n_face_vertices + 1;
      }
      else
        face_connect_idx[i+1] = face_connect_idx[i];
    }
  }

  BFT_MALLOC(_face_connect_idx, block_size + 1, cs_ccm_num_t);

  cs_part_to_block_copy_index(d, face_connect_idx, _face_connect_idx);

  /* Build connectivity */

  BFT_MALLOC(face_connect_g,
             n_faces + face_connect_size,
             cs_ccm_num_t);

  k = 0;

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < n_faces; i++) {
      if (face_vtx_idx[i+1] > face_vtx_idx[i]) {
        face_connect_g[k++] = face_vtx_idx[i+1] - face_vtx_idx[i];
        for (j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++)
          face_connect_g[k++] = b_mesh->global_vtx_num[face_vtx_lst[j]];
      }
    }
  }
  else if (entity == kCCMIOBoundaryFaces) {
    for (i = 0; i < n_faces; i++) {
      if (cell_gnum[face_cells[2*i]] == 0) {
        face_connect_g[k++] = face_vtx_idx[i+1] - face_vtx_idx[i];
        for (j = face_vtx_idx[i+1] - 1; j >= face_vtx_idx[i]; j--)
          face_connect_g[k++] = b_mesh->global_vtx_num[face_vtx_lst[j]];
      }
      else if (cell_gnum[face_cells[2*i + 1]] == 0) {
        face_connect_g[k++] = face_vtx_idx[i+1] - face_vtx_idx[i];
        for (j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++)
          face_connect_g[k++] = b_mesh->global_vtx_num[face_vtx_lst[j]];
      }
    }
  }

  BFT_MALLOC(_face_connect_g, _face_connect_idx[block_size], cs_ccm_num_t);

  cs_part_to_block_copy_indexed(d,
                                ccm_num_type,
                                face_connect_idx,
                                face_connect_g,
                                _face_connect_idx,
                                _face_connect_g);

  BFT_FREE(face_connect_g);
  BFT_FREE(face_connect_idx);

  /* Now write face -> vertices connectivity */

  cs_gnum_t g_connect_size = 0;
  cs_gnum_t range_base = _face_connect_idx[block_size], range_shift = 0;
  MPI_Scan(&range_base, &range_shift, 1, CS_MPI_GNUM, MPI_SUM, w->comm);
  range_shift -= _face_connect_idx[block_size];

  cs_gnum_t range[2] = {range_shift + _face_connect_idx[0] + 1,
                        range_shift + _face_connect_idx[block_size] + 1};

  MPI_Allreduce(&(range[1]), &g_connect_size, 1, CS_MPI_GNUM, MPI_MAX, w->comm);
  if (g_connect_size > 0)
    g_connect_size -= 1;

  BFT_FREE(_face_connect_idx);

  s = cs_file_serializer_create(sizeof(ccm_num_type),
                                1,
                                range[0],
                                range[1],
                                0,
                                _face_connect_g,
                                w->comm);

  do {

    _face_connect_g_s = cs_file_serializer_advance(s, range);

    if (_face_connect_g_s != NULL) { /* only possible on rank 0 */
      CCMIOError error = kCCMIONoErr, *err = &error;
      CCMIOWriteFaces(err, entity_id, entity, map_id,
                      CCMIOSIZEC(g_connect_size), _face_connect_g_s,
                      CCMIOINDEXC(range[0]-1), CCMIOINDEXC(range[1]-1));
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing face -> vertices connectivity."),
                  (int)error);
    }

  } while (_face_connect_g_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_face_connect_g);
}

/*----------------------------------------------------------------------------
 * Write face -> cells connectivity in parallel with periodic faces.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   face_bi     <-- face part to block info structure
 *   entity      <-- interior or boundary faces
 *   entity_id   <-- CCMIO id for this entity
 *   map_id      <-- CCMIO id for map
 *   cell_gnum   <-- array of global cell numbers, ordered by nodal mesh
 *   d           <-- face part to block distribution helper
 *   w           <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_face_cells_perio_g(const cs_mesh_t        *b_mesh,
                          cs_block_dist_info_t    face_bi,
                          CCMIOEntity             entity,
                          CCMIOID                 entity_id,
                          CCMIOID                 map_id,
                          const cs_gnum_t        *cell_gnum,
                          cs_part_to_block_t     *d,
                          fvm_to_ccm_writer_t    *w)
{
  cs_lnum_t i;

  cs_lnum_t n_cells_per_face = 0;
  cs_lnum_t n_faces = b_mesh->n_i_faces;

  cs_ccm_num_t *face_cell_g = NULL, *_face_cell_g = NULL;
  cs_ccm_num_t *_face_cell_g_s = NULL;
  cs_file_serializer_t *s = NULL;

  const cs_datatype_t ccm_num_type
    = (sizeof(cs_ccm_num_t) == 8) ? CS_INT64 : CS_INT32;

  const cs_lnum_t *face_cells = (const cs_lnum_t *)(b_mesh->i_face_cells);

  if (entity == kCCMIOInternalFaces)
    n_cells_per_face = 2;
  else if (entity == kCCMIOBoundaryFaces)
    n_cells_per_face = 1;

  /* Face -> cell connectivity */
  /*---------------------------*/

  BFT_MALLOC(face_cell_g, n_faces*n_cells_per_face, cs_ccm_num_t);

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < n_faces; i++) {
      face_cell_g[i*2]     = cell_gnum[face_cells[i*2]];
      face_cell_g[i*2 + 1] = cell_gnum[face_cells[i*2 + 1]];
    }
  }
  else {
    for (i = 0; i < n_faces; i++) {
      if (cell_gnum[face_cells[i*2] - 1] == 0)
        face_cell_g[i] = cell_gnum[face_cells[i*2 + 1]];
      else if (cell_gnum[face_cells[i*2 + 1] - 1] == 0)
        face_cell_g[i] = cell_gnum[face_cells[i*2]];
      else
        face_cell_g[i] = 0;
    }
  }

  /* Distribute to blocks and write */

  BFT_MALLOC(_face_cell_g,
             (face_bi.gnum_range[1] - face_bi.gnum_range[0]) * n_cells_per_face,
             cs_ccm_num_t);

  cs_part_to_block_copy_array(d,
                              ccm_num_type,
                              n_cells_per_face,
                              face_cell_g,
                              _face_cell_g);

  BFT_FREE(face_cell_g);

  /* Now write face -> cells connectivity */

  cs_gnum_t range[2] = {face_bi.gnum_range[0],
                        face_bi.gnum_range[1]};

  s = cs_file_serializer_create(sizeof(ccm_num_type),
                                n_cells_per_face,
                                range[0],
                                range[1],
                                0,
                                _face_cell_g,
                                w->comm);

  cs_ccm_num_t write_range[2] = {0, 0};

  do {

    _face_cell_g_s = cs_file_serializer_advance(s, range);

    if (_face_cell_g_s != NULL) { /* only possible on rank 0 */

      cs_lnum_t j = 0;
      cs_lnum_t n_elts = range[1] - range[0];
      CCMIOError error = kCCMIONoErr, *err = &error;

      /* Compact array before write */

      if (entity == kCCMIOInternalFaces) {
        for (i = 0; i < n_elts; i++) {
          if (_face_cell_g_s[i*2] > 0 && _face_cell_g_s[i*2+1] > 0) {
            _face_cell_g_s[j*2] = _face_cell_g_s[i*2];
            _face_cell_g_s[j*2+1] = _face_cell_g_s[i*2+1];
            j++;
          }
        }
      }
      else {
        for (i = 0; i < n_elts; i++) {
          if (_face_cell_g_s[i] > 0) {
            _face_cell_g_s[j] = _face_cell_g_s[i];
            j++;
          }
        }
      }

      write_range[0] = write_range[1];
      write_range[1] = write_range[0] + j;

      if (write_range[1] > write_range[0]) {

        CCMIOWriteFaceCells(err, entity_id, entity, map_id,
                            _face_cell_g_s,
                            CCMIOINDEXC(write_range[0]),
                            CCMIOINDEXC(write_range[1]));

        if (error != kCCMIONoErr)
          bft_error(__FILE__, __LINE__, 0,
                    _("CCMIO error %d writing face -> cells connectivity."),
                    (int)error);

      }

    }

  } while (_face_cell_g_s != NULL);

  cs_file_serializer_destroy(&s);

  BFT_FREE(_face_cell_g);
}

/*----------------------------------------------------------------------------
 * Write face information in parallel.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   entity      <-- interior or boundary faces
 *   topology_id <-- CCMIO id for this topology
 *   cell_gnum   <-- array of global cell numbers, ordered by nodal mesh
 *   w           <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_faces_g(const cs_mesh_t       *b_mesh,
               CCMIOEntity            entity,
               CCMIOID                topology_id,
               const cs_gnum_t       *cell_gnum,
               fvm_to_ccm_writer_t   *w)
{
  cs_block_dist_info_t face_bi;

  cs_lnum_t n_faces = 0;
  cs_ccm_num_t n_g_faces = 0;
  cs_ccm_num_t map_num_shift = 0;
  cs_gnum_t n_g_perio_faces = _count_faces_perio_g(b_mesh, cell_gnum, w);

  cs_gnum_t *_face_gnum = NULL;
  const cs_gnum_t *face_gnum = NULL;

  cs_part_to_block_t *d = NULL;

  CCMIOID entity_id, map_id;
  CCMIOError error = kCCMIONoErr, *err = &error;

  if (entity == kCCMIOInternalFaces) {
    n_faces = b_mesh->n_i_faces;
    n_g_faces = b_mesh->n_g_i_faces;
    map_num_shift = b_mesh->n_g_b_faces + n_g_perio_faces;
    face_gnum = b_mesh->global_i_face_num;
  }
  else if (entity == kCCMIOBoundaryFaces) {
    n_faces = b_mesh->n_b_faces;
    n_g_faces = b_mesh->n_g_b_faces;
    _face_gnum = _build_ordered_b_face_gnum(b_mesh);
    face_gnum = _face_gnum;
  }

  /* Build global face part to block distribution structures */

  face_bi = cs_block_dist_compute_sizes(w->rank,
                                        w->n_ranks,
                                        0,
                                        cs_parall_get_min_coll_buf_size(),
                                        n_g_faces);

  d = cs_part_to_block_create_by_gnum(w->comm,
                                      face_bi,
                                      n_faces,
                                      face_gnum);

  /* Create map and entity */

  if (w->rank < 1) {

    cs_gnum_t n_g_map_faces = n_g_faces;
    if (entity == kCCMIOInternalFaces)
      n_g_map_faces -= n_g_perio_faces;

    cs_block_dist_info_t map_face_bi = face_bi;

    if (n_g_perio_faces > 0 && entity == kCCMIOInternalFaces)
      map_face_bi = cs_block_dist_compute_sizes(w->rank,
                                                w->n_ranks,
                                                0,
                                                cs_parall_get_min_coll_buf_size(),
                                                n_g_faces - n_g_perio_faces);

    _write_map(NULL, n_g_map_faces, map_face_bi, map_num_shift, &map_id, w);

    if (entity == kCCMIOInternalFaces)
      CCMIONewEntity(err, topology_id, entity,
                     "Internal faces", &entity_id);

    else if (entity == kCCMIOBoundaryFaces) {
      w->b_face_map_id = map_id;
      CCMIONewIndexedEntity(err, topology_id, entity, 0,
                            "Boundary faces", &entity_id);
    }

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d writing faces entity."), (int)error);

  }

  /* When there is no periodicity or we are handling "true"
     boundary faces, use basic output functions */

  if (n_g_perio_faces == 0 || entity == kCCMIOBoundaryFaces) {

    _write_face_vertices_g(b_mesh,
                           face_bi,
                           entity,
                           entity_id,
                           map_id,
                           d,
                           w);

    _write_face_cells_g(b_mesh,
                        face_bi,
                        entity,
                        entity_id,
                        map_id,
                        cell_gnum,
                        d,
                        w);

  }

  /* In case of periodicity, use special functions, and
     add a second boundary section for interior periodic faces */

  if (n_g_perio_faces > 0) {

    if (entity == kCCMIOBoundaryFaces) {

      /* Rebuild global face part to block distribution structures */

      cs_part_to_block_destroy(&d);

      if (_face_gnum != NULL) {
        BFT_FREE(_face_gnum);
        face_gnum = NULL;
      }

      face_bi = cs_block_dist_compute_sizes(w->rank,
                                            w->n_ranks,
                                            0,
                                            cs_parall_get_min_coll_buf_size(),
                                            b_mesh->n_g_i_faces);

      d = cs_part_to_block_create_by_gnum(w->comm,
                                          face_bi,
                                          b_mesh->n_i_faces,
                                          b_mesh->global_i_face_num);

      if (w->rank < 1) {

        map_num_shift = b_mesh->n_g_b_faces;

        cs_block_dist_info_t map_face_bi
          = cs_block_dist_compute_sizes(w->rank,
                                        w->n_ranks,
                                        0,
                                        cs_parall_get_min_coll_buf_size(),
                                        n_g_perio_faces);

        _write_map(NULL, n_g_perio_faces, map_face_bi, map_num_shift,
                   &map_id, w);

        CCMIONewIndexedEntity(err, topology_id, entity, 1,
                              "Periodic faces", &entity_id);

        if (error != kCCMIONoErr)
          bft_error(__FILE__, __LINE__, 0,
                    _("CCMIO error %d writing faces entity."), (int)error);

      }

    }

    /* Use periodic output function, which filter faces
       to separate true internal from periodic faces */

    _write_face_vertices_perio_g(b_mesh,
                                 face_bi,
                                 entity,
                                 entity_id,
                                 map_id,
                                 cell_gnum,
                                 d,
                                 w);

    _write_face_cells_perio_g(b_mesh,
                              face_bi,
                              entity,
                              entity_id,
                              map_id,
                              cell_gnum,
                              d,
                              w);

  }

  /* Free face part to block distribution structures */

  cs_part_to_block_destroy(&d);

  if (_face_gnum != NULL)
    BFT_FREE(_face_gnum);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write vertices in serial mode.
 *
 * parameters:
 *   mesh     <-- pointer to mesh structure
 *   w        <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_vertices_l(const cs_mesh_t      *mesh,
                  fvm_to_ccm_writer_t  *w)
{
  CCMIOID map_id;
  CCMIOError error = kCCMIONoErr, *err = &error;

  cs_block_dist_info_t vtx_bi
    = cs_block_dist_compute_sizes(w->rank,
                                  w->n_ranks,
                                  0,
                                  cs_parall_get_min_coll_buf_size(),
                                  mesh->n_g_vertices);

  /* Write map and define entity */

  _write_map("Vertex map",
             mesh->n_g_vertices,
             vtx_bi,
             0,
             &map_id,
             w);

  CCMIONewEntity(err, w->root_id, kCCMIOVertices, "Vertices", &(w->vertices_id));

  /* Now write vertex coordinates */

  cs_gnum_t range[2] = {vtx_bi.gnum_range[0],
                        vtx_bi.gnum_range[1]};

  double scale = 1.0;
  if (sizeof(cs_real_t) == 8)
    CCMIOWriteVerticesd(err,
                        w->vertices_id,
                        CCMIOSIZEC(3),
                        scale,
                        map_id,
                        (const double *)(mesh->vtx_coord),
                        CCMIOINDEXC(range[0]-1),
                        CCMIOINDEXC(range[1]-1));
  else if (sizeof(cs_real_t) == 4)
    CCMIOWriteVerticesf(err,
                        w->vertices_id,
                        CCMIOSIZEC(3),
                        scale,
                        map_id,
                        (const float *)(mesh->vtx_coord),
                        CCMIOINDEXC(range[0]-1),
                        CCMIOINDEXC(range[1]-1));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing vertices."), (int)error);
}

/*----------------------------------------------------------------------------
 * Write cells in serial mode.
 *
 * parameters;
 *   b_mesh    <-- pointer to nodal mesh structure
 *   cell_gnum <-- pointer to global cell number, ordered by nodal mesh
 *   w         <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_cells_l(const cs_mesh_t      *b_mesh,
               const cs_gnum_t      *cell_gnum,
               fvm_to_ccm_writer_t  *w)
{
  cs_lnum_t  i;
  int *cell_family = NULL;
  CCMIOID map_id, topology_id, cells_id;
  CCMIOError error = kCCMIONoErr, *err = &error;

  cs_block_dist_info_t cell_bi
    = cs_block_dist_compute_sizes(w->rank,
                                  w->n_ranks,
                                  0,
                                  0,
                                  b_mesh->n_g_cells);

  /* Write map and define entity */

  _write_map("Cell map",
             b_mesh->n_g_cells,
             cell_bi,
             0,
             &map_id,
             w);

  w->cell_map_id = map_id;

  CCMIONewEntity(err, w->root_id, kCCMIOTopology, "Topology", &topology_id);
  CCMIONewEntity(err, topology_id, kCCMIOCells, "Cells", &cells_id);

  w->topology_id = topology_id;

  BFT_MALLOC(cell_family, b_mesh->n_cells, int);

  for (i = 0; i < b_mesh->n_cells; i++)
    cell_family[cell_gnum[i] - 1] = b_mesh->cell_family[i];

  /* Now write cell family info */

  cs_gnum_t range[2] = {cell_bi.gnum_range[0],
                        cell_bi.gnum_range[1]};

  CCMIOWriteCells(err, cells_id, map_id, cell_family,
                  CCMIOINDEXC(range[0]-1), CCMIOINDEXC(range[1]-1));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing cells."), (int)error);

  BFT_FREE(cell_family);
}

/*----------------------------------------------------------------------------
 * Write face -> vertices connectivity in serial mode.
 *
 * parameters:
 *   b_mesh        <-- pointer to base mesh structure
 *   entity        <-- interior or boundary faces
 *   entity_id     <-- CCMIO id for this entity
 *   map_id        <-- CCMIO id for map
 *   face_order    <-- face ordering array
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_l(const cs_mesh_t         *b_mesh,
                       CCMIOEntity              entity,
                       CCMIOID                  entity_id,
                       CCMIOID                  map_id,
                       const cs_lnum_t         *face_order)
{
  cs_lnum_t i, j, k;
  cs_lnum_t *face_vtx_idx = NULL, *face_vtx_lst = NULL;
  cs_ccm_num_t *face_connect = NULL;

  cs_lnum_t n_faces = 0;
  size_t face_connect_size = 0;

  if (entity == kCCMIOInternalFaces) {
    n_faces = b_mesh->n_i_faces;
    face_vtx_idx = b_mesh->i_face_vtx_idx;
    face_vtx_lst = b_mesh->i_face_vtx_lst;
  }
  else if (entity == kCCMIOBoundaryFaces) {
    n_faces = b_mesh->n_b_faces;
    face_vtx_idx = b_mesh->b_face_vtx_idx;
    face_vtx_lst = b_mesh->b_face_vtx_lst;
  }

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  face_connect_size = n_faces;
  for (i = 0; i < n_faces; i++)
    face_connect_size += face_vtx_idx[i+1] - face_vtx_idx[i];

  /* Build connectivity */

  BFT_MALLOC(face_connect, face_connect_size, cs_ccm_num_t);

  k = 0;
  for (i = 0; i < n_faces; i++) {
    cs_lnum_t face_id = face_order[i];
    face_connect[k++] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
    for (j = face_vtx_idx[face_id]; j < face_vtx_idx[face_id+1]; j++)
      face_connect[k++] = face_vtx_lst[j] + 1;
  }

  /* Now write face -> vertices connectivity */

  CCMIOError error = kCCMIONoErr, *err = &error;

  CCMIOWriteFaces(err, entity_id, entity, map_id,
                  CCMIOSIZEC(face_connect_size), face_connect,
                  CCMIOINDEXC(0), CCMIOINDEXC(face_connect_size));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing face -> vertices connectivity."),
              (int)error);

  BFT_FREE(face_connect);
}

/*----------------------------------------------------------------------------
 * Write face -> cells connectivity in serial mode.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   entity      <-- interior or boundary faces
 *   entity_id   <-- CCMIO id for this entity
 *   map_id      <-- CCMIO id for map
 *   face_order  <-- face ordering array
 *   cell_gnum   <-- pointer to global cell number, ordered by nodal mesh
 *----------------------------------------------------------------------------*/

static void
_write_face_cells_l(const cs_mesh_t        *b_mesh,
                    CCMIOEntity             entity,
                    CCMIOID                 entity_id,
                    CCMIOID                 map_id,
                    const cs_lnum_t        *face_order,
                    const cs_gnum_t        *cell_gnum)
{
  cs_lnum_t i;
  cs_lnum_t *face_cells = NULL;

  /* Face -> cell connectivity */
  /*---------------------------*/

  if (entity == kCCMIOInternalFaces) {
    BFT_MALLOC(face_cells, b_mesh->n_i_faces * 2, cs_ccm_num_t);
    for (i = 0; i < b_mesh->n_i_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      face_cells[i*2]     = cell_gnum[b_mesh->i_face_cells[face_id][0]];
      face_cells[i*2 + 1] = cell_gnum[b_mesh->i_face_cells[face_id][1]];
    }
  }
  else if (entity == kCCMIOBoundaryFaces) {
    BFT_MALLOC(face_cells, b_mesh->n_b_faces, cs_ccm_num_t);
    for (i = 0; i < b_mesh->n_b_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      face_cells[i] = cell_gnum[b_mesh->b_face_cells[face_id]];
    }
  }

  /* Now write face -> cells connectivity */

  CCMIOError error = kCCMIONoErr, *err = &error;

  CCMIOWriteFaceCells(err, entity_id, entity, map_id,
                      face_cells,
                      CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing face -> cells connectivity."),
              (int)error);

  BFT_FREE(face_cells);
}

/*----------------------------------------------------------------------------
 * Write face -> vertices connectivity in serial mode with periodic faces.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   entity      <-- interior or boundary faces
 *   entity_id   <-- CCMIO id for this entity
 *   map_id      <-- CCMIO id for map
 *   face_order  <-- face ordering array
 *   cell_gnum   <-- pointer to global cell number, ordered by nodal mesh
 *----------------------------------------------------------------------------*/

static void
_write_face_vertices_perio_l(const cs_mesh_t        *b_mesh,
                             CCMIOEntity             entity,
                             CCMIOID                 entity_id,
                             CCMIOID                 map_id,
                             const cs_lnum_t        *face_order,
                             const cs_gnum_t        *cell_gnum)
{
  cs_lnum_t i, j, k;

  cs_lnum_t n_faces = 0;

  cs_ccm_num_t *face_connect = NULL;

  const cs_lnum_t *face_vtx_idx = b_mesh->i_face_vtx_idx;
  const cs_lnum_t *face_vtx_lst = b_mesh->i_face_vtx_lst;
  const cs_lnum_2_t *face_cells = (const cs_lnum_2_t *)(b_mesh->i_face_cells);

  /* Allocate array large enough for both periodic boundary + true interior
     faces to avoid counting loop */

  const cs_lnum_t n_max_faces
    = CS_MAX(b_mesh->i_face_vtx_connect_size + b_mesh->n_i_faces,
             b_mesh->b_face_vtx_connect_size + b_mesh->n_b_faces);

  BFT_MALLOC(face_connect, n_max_faces, cs_ccm_num_t);

  /* Face -> vertex connectivity */
  /*-----------------------------*/

  k = 0;

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < b_mesh->n_i_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      if (   cell_gnum[face_cells[face_id][0]] > 0
          && cell_gnum[face_cells[face_id][1]] > 0) {
        n_faces += 1;
        face_connect[k++] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
        for (j = face_vtx_idx[face_id]; j < face_vtx_idx[face_id+1]; j++)
          face_connect[k++] = face_vtx_lst[j] + 1;
      }
    }
  }
  else if (entity == kCCMIOBoundaryFaces) {
    for (i = 0; i < b_mesh->n_i_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      if (cell_gnum[face_cells[i][0]] == 0) {
        n_faces += 1;
        face_connect[k++] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
        for (j = face_vtx_idx[face_id+1] - 1; j >= face_vtx_idx[face_id]; j--)
          face_connect[k++] = face_vtx_lst[j] + 1;
      }
      else if (cell_gnum[face_cells[i][1]] == 0) {
        n_faces += 1;
        face_connect[k++] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
        for (j = face_vtx_idx[face_id]; j < face_vtx_idx[face_id+1]; j++)
          face_connect[k++] = face_vtx_lst[j] + 1;
      }
    }
  }

  CCMIOError error = kCCMIONoErr, *err = &error;

  CCMIOWriteFaces(err, entity_id, entity, map_id,
                  CCMIOSIZEC(k), face_connect,
                  CCMIOINDEXC(0), CCMIOINDEXC(k));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing face -> vertices connectivity."),
              (int)error);

  BFT_FREE(face_connect);
}

/*----------------------------------------------------------------------------
 * Write face -> cells connectivity in serial mode with periodic faces.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   entity      <-- interior or boundary faces
 *   entity_id   <-- CCMIO id for this entity
 *   map_id      <-- CCMIO id for map
 *   face_order  <-- face ordering array
 *   cell_gnum   <-- array of global cell numbers, ordered by nodal mesh
 *----------------------------------------------------------------------------*/

static void
_write_face_cells_perio_l(const cs_mesh_t        *b_mesh,
                          CCMIOEntity             entity,
                          CCMIOID                 entity_id,
                          CCMIOID                 map_id,
                          const cs_lnum_t        *face_order,
                          const cs_gnum_t        *cell_gnum)
{
  cs_lnum_t i, j;

  cs_lnum_t n_cells_per_face = 0;

  cs_ccm_num_t *face_cells = NULL;

  if (entity == kCCMIOInternalFaces)
    n_cells_per_face = 2;
  else if (entity == kCCMIOBoundaryFaces)
    n_cells_per_face = 1;

  /* Face -> cell connectivity */
  /*---------------------------*/

  BFT_MALLOC(face_cells, b_mesh->n_i_faces*n_cells_per_face, cs_ccm_num_t);

  j = 0;

  if (entity == kCCMIOInternalFaces) {
    for (i = 0; i < b_mesh->n_i_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      if (   cell_gnum[b_mesh->i_face_cells[face_id][0]] > 0
          && cell_gnum[b_mesh->i_face_cells[face_id][1]] > 0) {
        face_cells[j*2]     = cell_gnum[b_mesh->i_face_cells[face_id][0]];
        face_cells[j*2 + 1] = cell_gnum[b_mesh->i_face_cells[face_id][1]];
        j += 1;
      }
    }
  }
  else {
    for (i = 0; i < b_mesh->n_i_faces; i++) {
      cs_lnum_t face_id = face_order[i];
      if (cell_gnum[b_mesh->i_face_cells[face_id][0]] == 0) {
        face_cells[j] = cell_gnum[b_mesh->i_face_cells[face_id][1]];
        j += 1;
      }
      else if (cell_gnum[b_mesh->i_face_cells[face_id][1]] == 0) {
        face_cells[j] = cell_gnum[b_mesh->i_face_cells[face_id][0]];
        j += 1;
      }
    }
  }

  /* Write connectivity */

  CCMIOError error = kCCMIONoErr, *err = &error;

  CCMIOWriteFaceCells(err, entity_id, entity, map_id,
                      face_cells,
                      CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing face -> cells connectivity."),
              (int)error);

  BFT_FREE(face_cells);
}

/*----------------------------------------------------------------------------
 * Write face information in serial mode.
 *
 * parameters:
 *   b_mesh      <-- pointer to base mesh structure
 *   mesh        <-- pointer to nodal mesh structure
 *   entity      <-- interior or boundary faces
 *   topology_id <-- CCMIO id for this topology
 *   cell_gnum   <-- array of global cell numbers, ordered by nodal mesh
 *   w           <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_faces_l(const cs_mesh_t       *b_mesh,
               CCMIOEntity            entity,
               CCMIOID                topology_id,
               const cs_gnum_t       *cell_gnum,
               fvm_to_ccm_writer_t   *w)
{
  cs_block_dist_info_t face_bi;

  cs_lnum_t n_faces = 0;
  cs_ccm_num_t n_g_faces = 0;
  cs_ccm_num_t  map_num_shift = 0;
  cs_gnum_t n_g_perio_faces = _count_faces_perio_g(b_mesh, cell_gnum, w);

  cs_lnum_t *face_order = NULL;

  CCMIOID entity_id, map_id;
  CCMIOError error = kCCMIONoErr, *err = &error;

  if (entity == kCCMIOInternalFaces) {
    n_faces = b_mesh->n_i_faces;
    n_g_faces = b_mesh->n_g_i_faces;
    map_num_shift = b_mesh->n_g_b_faces + n_g_perio_faces;
    face_order = _build_order_by_gnum(n_faces, b_mesh->global_i_face_num);
  }
  else if (entity == kCCMIOBoundaryFaces) {
    cs_gnum_t *face_gnum = _build_ordered_b_face_gnum(b_mesh);
    n_faces = b_mesh->n_b_faces;
    n_g_faces = b_mesh->n_g_b_faces;
    face_order = _build_order_by_gnum(n_faces, face_gnum);
    BFT_FREE(face_gnum);
  }

  if (n_g_perio_faces > 0 && entity == kCCMIOInternalFaces)
    n_g_faces -= n_g_perio_faces;

  face_bi = cs_block_dist_compute_sizes(w->rank,
                                        w->n_ranks,
                                        0,
                                        0,
                                        n_g_faces);

  /* Create map and entity */

  _write_map(NULL, n_g_faces, face_bi, map_num_shift, &map_id, w);

  if (entity == kCCMIOInternalFaces)
    CCMIONewEntity(err, topology_id, entity,
                   "Internal faces", &entity_id);

  else if (entity == kCCMIOBoundaryFaces) {
    w->b_face_map_id = map_id;
    CCMIONewIndexedEntity(err, topology_id, entity, 0,
                          "Boundary faces", &entity_id);
  }

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing faces entity."), (int)error);

  /* When there is no periodicity or we are handling "true"
     boundary faces, use basic output functions */

  if (n_g_perio_faces == 0 || entity == kCCMIOBoundaryFaces) {

    _write_face_vertices_l(b_mesh,
                           entity,
                           entity_id,
                           map_id,
                           face_order);

    _write_face_cells_l(b_mesh,
                        entity,
                        entity_id,
                        map_id,
                        face_order,
                        cell_gnum);

  }

  /* In case of periodicity, use special functions, and
     add a second boundary section for interior periodic faces */

  if (n_g_perio_faces > 0) {

    if (entity == kCCMIOBoundaryFaces) {

      /* Rebuild global face part to block distribution structures */

      BFT_FREE(face_order);
      face_order = _build_order_by_gnum(b_mesh->n_i_faces,
                                        b_mesh->global_i_face_num);

      face_bi = cs_block_dist_compute_sizes(w->rank,
                                            w->n_ranks,
                                            0,
                                            0,
                                            n_g_perio_faces);

      map_num_shift = b_mesh->n_g_b_faces;

      _write_map(NULL, n_g_perio_faces, face_bi, map_num_shift, &map_id, w);

      CCMIONewIndexedEntity(err, topology_id, entity, 1,
                            "Periodic faces", &entity_id);

      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing faces entity."), (int)error);

    }

    /* Use periodic output function, which filter faces
       to separate true internal from periodic faces */

    _write_face_vertices_perio_l(b_mesh,
                                 entity,
                                 entity_id,
                                 map_id,
                                 face_order,
                                 cell_gnum);

    _write_face_cells_perio_l(b_mesh,
                              entity,
                              entity_id,
                              map_id,
                              face_order,
                              cell_gnum);

  }

  BFT_FREE(face_order);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build cell global numbering array in order of output (as defined by
 * nodal mesh sections)
 *
 * parameters:
 *   mesh     <-- pointer to nodal mesh structure
 *   n_elts   <-- number of associated elements
 *   elt_dim  <-- dimension of the entities to consider
 *
 * returns:
 *   pointer to global numbering array
 *----------------------------------------------------------------------------*/

static cs_gnum_t *
_build_buffer_elt_gnum(const fvm_nodal_t  *mesh,
                       cs_lnum_t           n_elts,
                       int                 ent_dim)
{
  int i;
  cs_lnum_t j, k;
  cs_gnum_t num_shift = 0;
  cs_gnum_t *elt_gnum = NULL;

  BFT_MALLOC(elt_gnum, n_elts, cs_gnum_t);

  k = 0;

  for (i = 0; i < mesh->n_sections; i++) {
    const fvm_nodal_section_t  *const  section = mesh->sections[i];
    if (section->entity_dim == ent_dim) {
      const cs_gnum_t *g_num
        = fvm_io_num_get_global_num(section->global_element_num);
      for (j = 0; j < section->n_elements; j++, k++)
        elt_gnum[k] = g_num[j] + num_shift;
      num_shift += fvm_io_num_get_global_count(section->global_element_num);
    }
  }

  return elt_gnum;
}

/*----------------------------------------------------------------------------
 * Write field data in parallel mode.
 *
 * parameters:
 *   data_id          <-- entity to write the field data
 *   data_location    <-- field location (cells, faces)
 *   map_id           <-- map of the entity
 *   interlace        <-- indicates if variable in memory is interlaced
 *   dim_shift        <-- used to indicate the component in interlaced data
 *   datatype         <-- indicates the data type of (source) field values
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   field_values     <-- table of pointers on the field values
 *   b_mesh           <-- pointer to base mesh structure
 *   mesh             <-- pointer to nodal mesh structure
 *   w                <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_field_data_g(CCMIOID                 data_id,
                    CCMIODataLocation       data_location,
                    CCMIOID                 map_id,
                    int                     interlace,
                    int                     dim_shift,
                    int                     dimension,
                    cs_datatype_t           datatype,
                    int                     n_parent_lists,
                    const cs_lnum_t         parent_num_shift[],
                    const void       *const field_values[],
                    const cs_mesh_t        *b_mesh,
                    const fvm_nodal_t      *mesh,
                    fvm_to_ccm_writer_t    *w)
{
  cs_part_to_block_t *d = NULL;
  cs_file_serializer_t *s = NULL;
  void *_field_values_s = NULL;

  int ent_dim = -1;
  cs_lnum_t n_elts = 0;
  cs_ccm_num_t n_g_elts = 0;
  cs_gnum_t *_elt_gnum = NULL;
  const cs_gnum_t *elt_gnum = NULL;

  /* Choose if we have to write cell data or face data */
  switch (data_location) {
  case kCCMIOCell:
    ent_dim = 3;
    n_g_elts = b_mesh->n_g_cells;
    n_elts = b_mesh->n_cells;
    _elt_gnum = _build_buffer_elt_gnum(mesh, n_elts, 3);
    elt_gnum = _elt_gnum;
    break;
  case kCCMIOFace:
    ent_dim = 2;
    n_g_elts = b_mesh->n_g_b_faces;
    n_elts = b_mesh->n_b_faces;
    _elt_gnum = _build_buffer_elt_gnum(mesh, n_elts, 2);
    elt_gnum = _elt_gnum;
    break;
  case kCCMIOVertex:
    ent_dim = 0;
    n_g_elts = b_mesh->n_g_vertices;
    n_elts = b_mesh->n_vertices;
    elt_gnum = b_mesh->global_vtx_num;
    break;
  }

  /* Prepare part to block distribution */

  cs_block_dist_info_t elt_bi
    = cs_block_dist_compute_sizes(w->rank,
                                  w->n_ranks,
                                  0,
                                  cs_parall_get_min_coll_buf_size(),
                                  n_g_elts);

  d = cs_part_to_block_create_by_gnum(w->comm,
                                      elt_bi,
                                      n_elts,
                                      elt_gnum);

  unsigned char *_field_values_p = NULL, *_field_values_b = NULL;

  cs_lnum_t buffer_size = elt_bi.gnum_range[1] - elt_bi.gnum_range[0];
  cs_lnum_t part_size = cs_part_to_block_get_n_part_ents(d);

  cs_gnum_t range[2] = {elt_bi.gnum_range[0],
                        elt_bi.gnum_range[1]};

  /* Determine output data type */

  cs_datatype_t dst_datatype = CS_DATATYPE_NULL;

  if (   datatype == CS_INT32 || datatype == CS_INT64
      || datatype == CS_UINT32 || datatype == CS_UINT64)
    dst_datatype = _ccm_num_datatype;
  else
    dst_datatype = datatype;

  BFT_MALLOC(_field_values_b,
             buffer_size  * cs_datatype_size[dst_datatype],
             unsigned char);
  BFT_MALLOC(_field_values_p,
             part_size * cs_datatype_size[dst_datatype],
             unsigned char);

  /* Initialize the serializer with the appropriate type.

     We have to use the fvm_convert_array function because data
     may need to be de-interleaved, as shown in the following example for
     the velocity vector field:
     _________________________________________________________________________
     |                       |                        |                      |
     |  X1, X2,... Xn        |   Y1, Y2,... Yn        |  Z1, Z2,... Zn       |
     |_______________________|________________________|______________________|

     |__________|_________|_________|.............................|__________|
      rank 1     rank 2     rank 3                                 rank n

  */

  ccm_writer_section_t  *export_list = _build_export_list(mesh, ent_dim);

  for (ccm_writer_section_t *section = export_list;
       section != NULL;
       section = section->next) {

    cs_lnum_t src_shift = (n_parent_lists == 0) ? section->num_shift : 0;
    size_t dest_shift = section->num_shift * cs_datatype_size[dst_datatype];

    fvm_convert_array(dimension,
                      dim_shift,
                      1,
                      src_shift,
                      section->n_elts + src_shift,
                      interlace,
                      datatype,
                      dst_datatype,
                      n_parent_lists,
                      parent_num_shift,
                      section->parent_elt_num,
                      field_values,
                      _field_values_p + dest_shift);

  }

  /* TODO
     handle case where the nodal mesh is a subset of the base mesh */
  assert(fvm_nodal_get_n_entities(mesh, ent_dim) == part_size);

  BFT_FREE(export_list);

  /* Switch from partition to block distribution */

  cs_part_to_block_copy_array(d,
                              dst_datatype,
                              1,
                              _field_values_p,
                              _field_values_b);

  BFT_FREE(_field_values_p);

  /* Prepare Serializer */

  s = cs_file_serializer_create(cs_datatype_size[dst_datatype],
                                1,
                                range[0],
                                range[1],
                                0,
                                _field_values_b,
                                w->comm);

  /* Write the data in the file */

  do {

    _field_values_s = cs_file_serializer_advance(s, range);

    if (_field_values_s != NULL) { /* Only possible on rank 0 */

      CCMIOError    error = kCCMIONoErr, *err = &error;
      CCMIOIndex_t  start = CCMIOINDEXC(range[0]-1);
      CCMIOIndex_t  end = CCMIOINDEXC(range[1]-1);

      if (start == 0)
        start = kCCMIOStart;

      /* Write integer data */

      if (   datatype == CS_INT32 || datatype == CS_INT64
          || datatype == CS_UINT32 || datatype == CS_UINT64)
        CCMIOWriteFieldDatai(err,
                             data_id,
                             map_id,
                             data_location,
                             _field_values_s,
                             start,
                             end);

      /* Write float data */

      else if (datatype == CS_FLOAT)
        CCMIOWriteFieldDataf(err,
                             data_id,
                             map_id,
                             data_location,
                             _field_values_s,
                             start,
                             end);

      /* Write double data */

      else if (datatype == CS_DOUBLE)
        CCMIOWriteFieldDatad(err,
                             data_id,
                             map_id,
                             data_location,
                             _field_values_s,
                             start,
                             end);


      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d writing field data."), (int)error);

    }

  } while (_field_values_s != NULL);

  /* Free allocated memory */

  cs_file_serializer_destroy(&s);

  BFT_FREE(_field_values_b);

  cs_part_to_block_destroy(&d);

  if (_elt_gnum != NULL)
    BFT_FREE(_elt_gnum);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Write field data in serial mode
 *
 * parameters:
 *   data_id          <-- entity to write the field data
 *   data_location    <-- field location (cells, faces)
 *   map_id           <-- map of the entity
 *   interlace        <-- indicates if variable in memory is interlaced
 *   dim_shift        <-- used to indicate the component in interlaced data
 *   dimension        <-- dimension of the field
 *   datatype         <-- indicates the data type of (source) field values
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   field_values     <-- table of pointers on the field values
 *   b_mesh           <-- pointer to base mesh structure
 *   mesh             <-- pointer to nodal mesh structure
 *----------------------------------------------------------------------------*/

static void
_write_field_data_l(CCMIOID                     data_id,
                    CCMIODataLocation           data_location,
                    CCMIOID                     map_id,
                    int                         interlace,
                    int                         dim_shift,
                    int                         dimension,
                    cs_datatype_t               datatype,
                    int                         n_parent_lists,
                    const cs_lnum_t             parent_num_shift[],
                    const void           *const field_values[],
                    const cs_mesh_t            *b_mesh,
                    const fvm_nodal_t          *mesh)
{
  int ent_dim = -1;
  CCMIOError error = kCCMIONoErr, *err = &error;
  cs_ccm_num_t start = kCCMIOStart;
  cs_ccm_num_t end = 0;
  cs_datatype_t dst_datatype = CS_DATATYPE_NULL;

  unsigned char *_field_values = NULL;

  /* Choose if we have to write cell data or face data */

  switch (data_location) {
  case kCCMIOCell:
    end = b_mesh->n_cells;
    ent_dim = 3;
    break;
  case kCCMIOFace:
    end = b_mesh->n_b_faces;
    ent_dim = 2;
    break;
  case kCCMIOVertex:
    end = b_mesh->n_vertices;
    ent_dim = 0;
    break;
  }

  if (   datatype == CS_INT32 || datatype == CS_INT64
      || datatype == CS_UINT32 || datatype == CS_UINT64)
    dst_datatype = _ccm_num_datatype;
  else
    dst_datatype = datatype;

  BFT_MALLOC(_field_values,
             end * cs_datatype_size[dst_datatype],
             unsigned char);

  ccm_writer_section_t  *export_list = _build_export_list(mesh, ent_dim);

  for (ccm_writer_section_t *section = export_list;
       section != NULL;
       section = section->next) {

    cs_lnum_t src_shift = (n_parent_lists == 0) ? section->num_shift : 0;
    size_t dest_shift = section->num_shift * cs_datatype_size[dst_datatype];

    fvm_convert_array(dimension,
                      dim_shift,
                      1,
                      src_shift,
                      section->n_elts + src_shift,
                      interlace,
                      datatype,
                      dst_datatype,
                      n_parent_lists,
                      parent_num_shift,
                      section->parent_elt_num,
                      field_values,
                      _field_values + dest_shift);

  }

  /* TODO
     handle case where the nodal mesh is a subset of the base mesh */
  assert(fvm_nodal_get_n_entities(mesh, ent_dim) == end);

  BFT_FREE(export_list);

  /* Write integer data */

  if (   datatype == CS_INT32 || datatype == CS_INT64
      || datatype == CS_UINT32 || datatype == CS_UINT64)
    CCMIOWriteFieldDatai(err,
                         data_id,
                         map_id,
                         data_location,
                         (void *)_field_values,
                         start,
                         end);

  /* Write float data */

  else if (datatype == CS_FLOAT)
    CCMIOWriteFieldDataf(err,
                         data_id,
                         map_id,
                         data_location,
                         (void *)_field_values,
                         start,
                         end);

  /* Write double data */

  else if (datatype == CS_DOUBLE)
    CCMIOWriteFieldDatad(err,
                         data_id,
                         map_id,
                         data_location,
                         (void *)_field_values,
                         start,
                         end);

  BFT_FREE(_field_values);

  if (error != kCCMIONoErr)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO error %d writing field data."), (int)error);
}

/*----------------------------------------------------------------------------
 * Write a new field
 *
 * parameters:
 *   name           <-- the field's name
 *   short_name     <-- a shorter name for the field
 *   phase_id       <-- phase to write the field data
 *   field_id       <-- id of the parent field
 *   dimension      <-- variable dimension
 *   datatype       <-- indicates the data type of (source) field values
 *   interlace      <-- indicates if variable in memory is interlaced
 *   field_values   <-- table of pointers to field values
 *   data_location  <-- datalocation kCCMIOCell, kCCMIOFace...
 *   map_id         <-- id of the map relative to the datalocation
 *   b_mesh         <-- pointer to base mesh structure
 *   mesh           <-- pointer to nodal mesh structure
 *   w              <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_multidimensional_field_data(const char              *name,
                                   char                    *short_name,
                                   CCMIOID                  phase_id,
                                   CCMIOID                  field_id,
                                   int                      dimension,
                                   cs_datatype_t            datatype,
                                   CCMIODataLocation        data_location,
                                   int                      interlace,
                                   int                      n_parent_lists,
                                   const cs_lnum_t          parent_num_shift[],
                                   const void        *const field_values[],
                                   CCMIOID                  map_id,
                                   const cs_mesh_t         *b_mesh,
                                   const fvm_nodal_t       *mesh,
                                   fvm_to_ccm_writer_t     *w)
{
  /* Set up the variable names */

  char *full_name = NULL;
  BFT_MALLOC(full_name, strlen("CS_") + strlen(name) + strlen("XX")  + 1, char);
  full_name[0] = '\0';

  char vect_short_name[7]="";
  char vector_char1[2] = "X";
  char vector_char2[2] = "X";
  int i;

  CCMIOID data_id;
  CCMIOComponent component;

  CCMIOError error = kCCMIONoErr, *err = &error;

  /* Each iteration writes one component of the multidimensional field */

  for (i = 0; i < dimension; i++) {

    CCMIOID child_field_id;

    /* Set up the component to be written and the name of this component */
    switch (dimension) {

    case 3: /* Vector */
      switch (i) {
      case 0:
        component = kCCMIOVectorX;
        vector_char1[0] = 'X';
        break;
      case 1:
        component = kCCMIOVectorY;
        vector_char1[0] = 'Y';
        break;
      case 2:
        component = kCCMIOVectorZ;
        vector_char1[0] = 'Z';
        break;
      }
      break;

    case 6: /* Symmetrical tensor */
      switch (i) {
      case 0:
        component = kCCMIOTensorXX;
        vector_char1[0] = 'X';
        vector_char2[0] = 'X';
            break;
      case 1:
        component = kCCMIOTensorXY;
        vector_char1[0] = 'Y';
        vector_char2[0] = 'X';
        break;
      case 2:
        component = kCCMIOTensorXZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'X';
        break;
      case 3:
        component = kCCMIOTensorYY;
        vector_char1[0] = 'Y';
        vector_char2[0] = 'Y';
        break;
      case 4:
        component = kCCMIOTensorYZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'Y';
        break;
      case 5:
        component = kCCMIOTensorZZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'Z';
        break;
      }
      break;

    case 9: /* Unsymmetrical tensor */
      switch (i) {
      case 0:
        component = kCCMIOTensorXX;
        vector_char1[0] = 'X';
        vector_char2[0] = 'X';
        break;
      case 1:
        component = kCCMIOTensorXY;
        vector_char1[0] = 'Y';
        vector_char2[0] = 'X';
        break;
      case 2:
        component = kCCMIOTensorXZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'X';
        break;
      case 3:
        component = kCCMIOTensorYX;
        vector_char1[0] = 'X';
        vector_char2[0] = 'Y';
        break;
      case 4:
        component = kCCMIOTensorYY;
        vector_char1[0] = 'Y';
        vector_char2[0] = 'Y';
        break;
      case 5:
        component = kCCMIOTensorYZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'Y';
        break;
      case 6:
        component = kCCMIOTensorZX;
        vector_char1[0] = 'X';
        vector_char2[0] = 'Z';
        break;
      case 7:
        component = kCCMIOTensorZY;
        vector_char1[0] = 'Y';
        vector_char2[0] = 'Z';
        break;
      case 8:
        component = kCCMIOTensorZZ;
        vector_char1[0] = 'Z';
        vector_char2[0] = 'Z';
        break;
      }
      break;

    default: /* Unknown */
      bft_error(__FILE__, __LINE__, 0,
                  _("Incorrect multidimensional field data format"));
      break;
    }

    /* Set up field name */
    strcpy(full_name, "CS_");
    strcat(full_name, name);

    /* Add the second dimension if we write a tensor */
    if (i >= 3 || dimension >= 6)
      strcat(full_name, vector_char2);
    strcat(full_name, vector_char1);

    /* Star-ccm+ will import properly the vector fields only if the short
       names of the components are SU, SV and SW */
    if (dimension == 3) {
      if (tolower(name[0]) == 'e')
        strcpy(vect_short_name, "E");      /* efforts */
      else if (tolower(name[0]) == 'm')
        strcpy(vect_short_name, "M");      /* mesh velocity */
      else
        strcpy(vect_short_name, "S");      /* velocity */

      switch (i) {
      case(0):
        strcat(vect_short_name, "U");
        break;
      case(1):
        strcat(vect_short_name, "V");
        break;
      case(2):
        strcat(vect_short_name, "W");
        break;
      }
    }

    /* For the other field dimensions it's OK */
    else {
      strcpy(vect_short_name, short_name);
      if (i >= 3 || dimension >= 6)
        strcat(full_name, vector_char2);
      strcat(vect_short_name, vector_char1);
    }

    /* Write the component as a scalar field and the data entity */
    if (w->rank < 1) {
      CCMIONewField(err,
                    phase_id,
                    full_name,
                    vect_short_name,
                    kCCMIOScalar,
                    &child_field_id);
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d creating field: %s."),
                  (int)error, full_name);
      CCMIONewEntity(err,
                     child_field_id,
                     kCCMIOFieldData,
                     NULL,
                     &data_id);
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d creating new entity."), (int)error);
    }

    /* Write data in parallel */

#if defined(HAVE_MPI)

    if (w->n_ranks > 1)
      _write_field_data_g(data_id,
                          data_location,
                          map_id,
                          interlace,
                          i,
                          dimension,
                          datatype,
                          n_parent_lists,
                          parent_num_shift,
                          field_values,
                          b_mesh,
                          mesh,
                          w);

#endif

    /* Write data in serial mode */
    if (w->n_ranks == 1)
      _write_field_data_l(data_id,
                          data_location,
                          map_id,
                          interlace,
                          i,
                          dimension,
                          datatype,
                          n_parent_lists,
                          parent_num_shift,
                          field_values,
                          b_mesh,
                          mesh);

    /* Link the data node with component of the parent field node */
    if (w->rank < 1) {
      CCMIOWriteMultiDimensionalFieldData(err,
                                          field_id,
                                          component,
                                          child_field_id);
      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d creating multidimensional field data."),
                  (int)error);
    }
  }

  BFT_FREE(full_name);
}

/*----------------------------------------------------------------------------
 * Write a new field

 * parameters:
 *   name             <-- name of the field to be written
 *   phase_id         <-- phase to write the field data
 *   datatype         <-- indicates the data type of (source) field values
 *   dimension        <-- variable dimension
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   field_values     <-- table of pointers to field values
 *   entity_dim       <-- entity type, i.e 3: cell, 2:face
 *   b_mesh           <-- pointer to base mesh structure
 *   mesh             <-- pointer to nodal mesh structure
 *   w                <-> pointer to writer structure
 *----------------------------------------------------------------------------*/

static void
_write_field(const char                 *name,
             CCMIOID                     phase_id,
             int                         dimension,
             cs_datatype_t               datatype,
             int                         interlace,
             int                         n_parent_lists,
             const cs_lnum_t             parent_num_shift[],
             const void           *const field_values[],
             int                         entity_dim,
             const cs_mesh_t            *b_mesh,
             const fvm_nodal_t          *mesh,
             fvm_to_ccm_writer_t        *w)
{
  CCMIOError error = kCCMIONoErr, *err = &error;
  CCMIOID field_id, map_id, data_id;
  CCMIODataLocation data_location;

  char short_name[15];
  strncpy(short_name, name, 4);
  short_name[4] = '\0';

  CCMIODimensionality dimensionality;

  switch (dimension) {

  case 0: /* Constant */
    dimensionality = kCCMIODimNull;
    break;

  case 1: /* Scalar */
    dimensionality = kCCMIOScalar;
    break;

  case 3: /* Vector */
    dimensionality = kCCMIOVector;
    break;

  case 6: /* Tensor */
  case 9:
    dimensionality = kCCMIOTensor;
    break;

  default: /* Unknown */
    bft_error(__FILE__, __LINE__, 0,
              _("Unhandled field data format"));
    dimensionality = kCCMIODimNull;
    break;

  }

  switch (entity_dim) {

  case 0: /* Vertex */
    map_id = w->vtx_map_id;
    data_location = kCCMIOVertex;
    break;

  case 2: /* Face */
    map_id = w->b_face_map_id;
    data_location = kCCMIOFace;
    break;

  case 3: /* Cell */
    map_id = w->cell_map_id;
    data_location = kCCMIOCell;
    break;

  default: /* Unknown */
    bft_error(__FILE__, __LINE__, 0,
              _("Incorrect entity type to store field data"));
    break;

  }

  char *full_name = NULL;
  BFT_MALLOC(full_name, strlen(name) + strlen("CS_") + 1, char);
  strcpy(full_name, "CS_");
  strcat(full_name, name);

  /* Write the field */

  if (w->rank < 1) {

    CCMIONewField(err,
                  phase_id,
                  full_name,
                  short_name,
                  dimensionality,
                  &field_id);

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d creating field: %s."),
                  (int)error, full_name);

  }

  /* Write the data now */

  switch (dimensionality) {

  case kCCMIODimNull: /* Constant (not handled) */
    bft_error(__FILE__, __LINE__, 0,
              _("Type of post data not handled at the moment"));
    break;

  case kCCMIOScalar:

    if (w->rank < 1)
      CCMIONewEntity(err, field_id, kCCMIOFieldData, NULL, &data_id);

    /* Write data in parallel */

#if defined(HAVE_MPI)

    if (w->n_ranks > 1)
      _write_field_data_g(data_id,
                          data_location,
                          map_id,
                          interlace,
                          0,
                          dimension,
                          datatype,
                          n_parent_lists,
                          parent_num_shift,
                          field_values,
                          b_mesh,
                          mesh,
                          w);

#endif

    /* Write data in serial mode */

    if (w->n_ranks == 1)
      _write_field_data_l(data_id,
                          data_location,
                          map_id,
                          interlace,
                          0,
                          dimension,
                          datatype,
                          n_parent_lists,
                          parent_num_shift,
                          field_values,
                          b_mesh,
                          mesh);

    break;

  case kCCMIOVector:
  case kCCMIOTensor:
    _write_multidimensional_field_data(name,
                                       short_name,
                                       phase_id,
                                       field_id,
                                       dimension,
                                       datatype,
                                       data_location,
                                       interlace,
                                       n_parent_lists,
                                       parent_num_shift,
                                       field_values,
                                       map_id,
                                       b_mesh,
                                       mesh,
                                       w);
    break;

  }

  BFT_FREE(full_name);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Returns number of library version strings associated with the CCMIO format.
 *
 * returns:
 *   number of library version strings associated with the CCMIO format.
 *----------------------------------------------------------------------------*/

int
fvm_to_ccm_n_version_strings(void)
{
  return 1;
}

/*----------------------------------------------------------------------------
 * Returns a library version string associated with the CCMIO format.
 *
 * In certain cases, when using dynamic libraries, fvm may be compiled
 * with one library version, and linked with another. If both run-time
 * and compile-time version information is available, this function
 * will return the run-time version string by default.
 *
 * Setting the compile_time flag to 1, the compile-time version string
 * will be returned if this is different from the run-time version.
 * If the version is the same, or only one of the 2 version strings are
 * available, a NULL character string will be returned with this flag set.
 *
 * parameters:
 *   string_index <-- index in format's version string list (0 to n-1)
 *   compile_time <-- 0 by default, 1 if we want the compile-time version
 *                    string, if different from the run-time version.
 *
 * returns:
 *   pointer to constant string containing the library's version.
 *----------------------------------------------------------------------------*/

const char *
fvm_to_ccm_version_string(int string_index,
                          int compile_time_version)
{
  CS_UNUSED(string_index);
  CS_UNUSED(compile_time_version);

#if    defined(kCCMIOMajorVersion) && defined(kCCMIOMinorVersion) \
    && defined(kCCMIORevision)
  snprintf(_ccm_version_string, 31, "%s %d.%d.%d", "CCMIO",
           kCCMIOMajorVersion, kCCMIOMinorVersion, kCCMIORevision);
#else
  sprintf(_ccm_version_string, "%s", "CCMIO");
#endif

  return _ccm_version_string;
}

/*----------------------------------------------------------------------------
 * Initialize FVM to CCMIO file writer.
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque CCMIO writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_ccm_init_writer(const char             *name,
                       const char             *path,
                       const char             *options,
                       fvm_writer_time_dep_t   time_dependency,
                       MPI_Comm                comm)
#else
void *
fvm_to_ccm_init_writer(const char             *name,
                       const char             *path,
                       const char             *options,
                       fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(options);

  int  i;
  int  mesh_filename_length, mesh_basename_length, name_length, path_length;

  int writer_index = 0;
  fvm_to_ccm_writer_t  *writer = NULL;

  /* Initialize writer */

  BFT_MALLOC(writer, 1, fvm_to_ccm_writer_t);

  /* Mesh metadata */

  writer->v_mesh = NULL;
  writer->b_mesh = NULL;

  writer->n_g_perio_faces = 0;

  /* Mesh time dependency */

  if (time_dependency != FVM_WRITER_FIXED_MESH)
    bft_error(__FILE__, __LINE__, 0,
              _("CCMIO output can currently handle only "
                "non-time-dependent meshes."));

  writer->time_dependency = time_dependency;

  writer->state_counter = 1;
  writer->mesh_time.n_time_values = 0;
  writer->mesh_time.last_time_step = -1;
  writer->mesh_time.time_value = NULL;

  /* Field time dependency */

  writer->field_time.n_time_values = 0;
  writer->field_time.last_time_step = 0;
  writer->field_time.time_value = NULL;

  writer->n_time_fields[0] = 0;
  writer->n_time_fields[1] = 0;
  writer->n_time_fields[2] = 0;

  /* Writer name */

  name_length = strlen(name);
  if (name_length == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty CCMIO filename."));
  BFT_MALLOC(writer->name, name_length + 1, char);
  strcpy(writer->name, name);

  for (i = 0; i < name_length; i++) {
    if (writer->name[i] == ' ' || writer->name[i] == '\t')
      writer->name[i] = '_';
  }

  /* Writer's associated filename */

  if (path != NULL)
    path_length = strlen(path);
  else
    path_length = 0;
  mesh_filename_length = path_length + name_length + strlen(".ccmg") + 1;
  BFT_MALLOC(writer->mesh_filename, mesh_filename_length, char);

  writer->solution_filename = NULL;

  if (path != NULL)
    strcpy(writer->mesh_filename, path);
  else
    writer->mesh_filename[0] = '\0';

  strcat(writer->mesh_filename, writer->name);
  strcat(writer->mesh_filename, ".ccmg");

  mesh_basename_length = name_length + strlen(".ccmg") + 1;
  BFT_MALLOC(writer->mesh_basename, mesh_basename_length, char);
  strcpy(writer->mesh_basename, writer->name);
  strcat(writer->mesh_basename, ".ccmg");

  BFT_MALLOC(writer->path, strlen(path)+1, char);
  strcpy(writer->path, path);

  /* CCMIO Base structure */

  /* Other variables */

  writer->rank = 0;
  writer->n_ranks = 1;

  /* Open CCMIO file */

  writer->is_open = false;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag && comm != MPI_COMM_NULL) {
      writer->comm = comm;
      MPI_Comm_rank(writer->comm, &rank);
      MPI_Comm_size(writer->comm, &n_ranks);
      writer->rank = rank;
      writer->n_ranks = n_ranks;
    }
    else
      writer->comm = MPI_COMM_NULL;
  }
#endif /* defined(HAVE_MPI) */

#if defined(HAVE_MPI)
  if (writer->n_ranks > 1)
    MPI_Bcast(&writer_index, 1, MPI_INT, 0, writer->comm);
#endif

#if 0
  cs_parall_set_min_coll_buf_size(0); /* for testing */
#endif

  /* Artificially force link of ADF library using gold linker */
  _force_adf_link(false);

  return writer;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to CCMIO file writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque CCMIO writer structure.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

void *
fvm_to_ccm_finalize_writer(void  *this_writer_p)
{
  fvm_to_ccm_writer_t *w = this_writer_p;

  w->is_open = false;

  /* Free memory */

  BFT_FREE(w->path);
  BFT_FREE(w->solution_filename);
  BFT_FREE(w->mesh_basename);
  BFT_FREE(w->mesh_filename);
  BFT_FREE(w->name);
  BFT_FREE(w->mesh_time.time_value);
  BFT_FREE(w->field_time.time_value);
  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time value with a writer structure if necessary.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_ccm_set_mesh_time(void     *this_writer_p,
                         int       time_step,
                         double    time_value)
{
  fvm_to_ccm_writer_t *w = this_writer_p;

 /* Mark meshes as unset to allow re-export */

  if (time_step != w->mesh_time.last_time_step) {
    w->v_mesh = NULL;
    w->b_mesh = NULL;
  }

  /* Update the current mesh time */

  _update_time(MESH_TIME, &w->mesh_time, time_step, time_value, w);

}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a CCMIO file
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer.
 *   mesh          <-- pointer to nodal mesh structure that should be written.
 *----------------------------------------------------------------------------*/

void
fvm_to_ccm_export_nodal(void               *this_writer_p,
                        const fvm_nodal_t  *mesh)
{
  fvm_to_ccm_writer_t *w = this_writer_p;

  CCMIOError error = kCCMIONoErr, *err = &error;

  bool allow_export = false;

  const cs_mesh_t *b_mesh = cs_glob_mesh;

  /* Only export complete volume mesh, with usable boundary.
     Boundary is a special case: it may only be mapped to boundary
     of complete volume mesh, after that mesh has been exported.
     Also, the mesh builder is tested for, as it is a simple
     means of checking we are finished building the main mesh. */

  int entity_dim = fvm_nodal_get_max_entity_dim(mesh);

  if (entity_dim == 3) {
    if (  w->v_mesh == NULL
        && _n_g_mesh_elts(mesh, entity_dim) == b_mesh->n_g_cells
        && cs_glob_mesh_builder == NULL) {
      allow_export = true;
      w->v_mesh = mesh;
    }
  }
  else if (entity_dim == 2) {
    if (   w->b_mesh == NULL
        && _n_g_mesh_elts(mesh, entity_dim) == b_mesh->n_g_b_faces)
      w->b_mesh = mesh;
  }

  if (allow_export == false)
    return;

  if (w->rank < 1) {

    CCMIOID root;

    /* Generate a new ccmg filename */

    if (   (w->mesh_time.last_time_step == -1 && w->state_counter == 1)
        ||  w->mesh_time.last_time_step != -1) {

      if (w->time_dependency != FVM_WRITER_FIXED_MESH) {
        char s_time_step[16] ="";
        sprintf(s_time_step, "%d", w->mesh_time.last_time_step);
        int path_length = strlen(w->path) + strlen(w->name) + 1
                          + strlen(s_time_step) + strlen(".ccmg") + 1;
        BFT_REALLOC(w->mesh_filename, path_length, char);
        sprintf(w->mesh_filename, "%s%s-%s.ccmg",
                w->path, w->name, s_time_step);
        path_length =   strlen(w->name) + 1
                      + strlen(s_time_step) + strlen(".ccmg") + 1;
        BFT_REALLOC(w->mesh_basename, path_length, char);
        sprintf(w->mesh_basename, "%s-%s.ccmg",
                w->name, s_time_step);
      }
      else {
        int path_length = strlen(w->path) + strlen(w->name) + strlen(".ccmg") + 1;
        BFT_REALLOC(w->mesh_filename, path_length, char);
        sprintf(w->mesh_filename, "%s%s.ccmg", w->path, w->name);
        path_length = strlen(w->name) + strlen(".ccmg") + 1;
        BFT_REALLOC(w->mesh_basename, path_length, char);
        sprintf(w->mesh_basename, "%s.ccmg", w->name);
      }
    }

    /* Open file for output */

    cs_file_remove(w->mesh_filename);
    CCMIOOpenFile(err, w->mesh_filename, kCCMIOWrite, &root);
    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIOOpenFile() failed to open file \"%s\"\n"
                  "CCMIO error %d."),
                w->mesh_filename, (int)error);
    w->root_id = root;
    w->is_open = true;
  }

  /* Write mesh data */
  /*-----------------*/

#if defined(HAVE_MPI)

  if (w->n_ranks > 1) {

    /* Export mesh only one time if it is static */

    if (   (w->mesh_time.last_time_step == -1 && w->state_counter == 1)
        || w->mesh_time.last_time_step != -1) {

      /* Build global cell numbering including parallel halos,
         except for periodic values */

      cs_gnum_t  *cell_gnum = _build_ordered_cell_gnum(b_mesh, mesh);

      _write_state(w);
      _write_processor(w);
      _write_vertices_g(b_mesh, w);
      _write_cells_g(b_mesh, cell_gnum, w);
      _write_faces_g(b_mesh,
                     kCCMIOInternalFaces,
                     w->topology_id,
                     cell_gnum,
                     w);
      _write_faces_g(b_mesh,
                     kCCMIOBoundaryFaces,
                     w->topology_id,
                     cell_gnum,
                     w);
      _write_problem_description(w);
      _write_solution(w);
      _write_restart_info(w->mesh_time.last_time_step,
                          w->mesh_time.time_value[w->state_counter-2],
                          0.0,
                          w);

      _finalize_processor(NULL, NULL, w);

      BFT_FREE(cell_gnum);
    }
  }

#endif

  if (w->n_ranks == 1) {

    /* Export mesh only one time if it is static */
    if (  (w->mesh_time.last_time_step == -1 && w->state_counter == 1)
        || w->mesh_time.last_time_step != -1) {

      /* Build global cell numbering including parallel halos,
         except for periodic values */

      cs_gnum_t  *cell_gnum = _build_ordered_cell_gnum(b_mesh, mesh);

      _write_state(w);
      _write_processor(w);
      _write_vertices_l(b_mesh, w);
      _write_cells_l(b_mesh, cell_gnum, w);
      _write_faces_l(b_mesh,
                     kCCMIOInternalFaces,
                     w->topology_id,
                     cell_gnum,
                     w);
      _write_faces_l(b_mesh,
                     kCCMIOBoundaryFaces,
                     w->topology_id,
                     cell_gnum,
                     w);
      _write_problem_description(w);
      _write_solution(w);
      _write_restart_info(w->mesh_time.last_time_step,
                          w->mesh_time.time_value[w->state_counter-2],
                          0.0,
                          w);

      _finalize_processor(NULL, NULL, w);

      BFT_FREE(cell_gnum);
    }

    if (w->rank < 1) {

      CCMIOCloseFile(err, w->root_id);

      if (error != kCCMIONoErr)
        bft_error(__FILE__, __LINE__, 0,
                  _("CCMIO error %d closing file."), (int)error);

    }

    w->is_open = false;

  }

}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a CCMIO file.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- fvm grid location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_ccm_export_field(void                   *this_writer_p,
                        const fvm_nodal_t      *mesh,
                        const char             *name,
                        fvm_writer_var_loc_t    location,
                        int                     dimension,
                        cs_interlace_t          interlace,
                        int                     n_parent_lists,
                        const cs_lnum_t         parent_num_shift[],
                        cs_datatype_t           datatype,
                        int                     time_step,
                        double                  time_value,
                        const void       *const field_values[])
{
  fvm_to_ccm_writer_t *w = this_writer_p;

  const cs_mesh_t *b_mesh = cs_glob_mesh;

  CCMIOID root;
  CCMIOError error = kCCMIONoErr, *err = &error;
  CCMIOID phase_id;

  /* Get entity dimension (0: vertex; 2: face, 3: cell) */

  int entity_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Ensure the matching mesh has been exported */

  if (   (entity_dim != 3 || mesh != w->v_mesh)
      && (entity_dim != 2 || mesh != w->b_mesh))
    return;

  if (location == FVM_WRITER_PER_NODE && entity_dim == 2)
    return;

  if (location == FVM_WRITER_PER_NODE)
    entity_dim = 0;

  if (w->rank < 1) {

    bool new_time_value = false;

    if (    w->field_time.n_time_values == 0
        ||  w->field_time.time_value[w->field_time.n_time_values-1]
            < time_value)
      new_time_value = true;

    /* Prepare new ccmp filename */

    if (time_step > -1) {
      char s_time_step[16] ="";
      sprintf(s_time_step, "%d", time_step);
      int path_length =   strlen(w->path) + strlen(w->name) + 1
                        + strlen(s_time_step) + strlen(".ccmp") + 1;
      BFT_REALLOC(w->solution_filename, path_length, char);
      sprintf(w->solution_filename, "%s%s-%s.ccmp",
              w->path, w->name, s_time_step);
    }
    else {
      int path_length = strlen(w->path) + strlen(w->name) + strlen(".ccmp") + 1;
      BFT_REALLOC(w->solution_filename, path_length, char);
      sprintf(w->solution_filename, "%s%s.ccmp",
              w->path, w->name);
    }

    /* Open file */

    if (new_time_value)
      cs_file_remove(w->solution_filename);

    CCMIOOpenFile(err, w->solution_filename, kCCMIOWrite, &root);
    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIOOpenFile() failed to open file \"%s\"\n"
                  "CCMIO error %d."),
                w->solution_filename, (int)error);
    w->root_id = root;
    w->is_open = true;

    /* Write new field info */

    if (new_time_value) {

      _write_state(w);
      _write_processor(w);
      _write_solution(w);

      w->n_time_fields[0] = 0;
      w->n_time_fields[1] = 0;
      w->n_time_fields[2] = 0;

      _write_problem_description(w);
      _write_restart_info(time_step,
                          time_value,
                          0.0,
                          w);

      _update_time(FIELD_TIME, &w->field_time, time_step, time_value, w);

      /* Prepare mesh_filename */

      /* Link solution file to geometry file */

      _finalize_processor(w->mesh_basename,
                          w->mesh_basename,
                          w);
    }

    if (location == FVM_WRITER_PER_NODE) {
      if (w->n_time_fields[2] == 0)
        _write_vertices_map(b_mesh, w);
      w->n_time_fields[2] += 1;
    }
    else {
      if (entity_dim == 2) {
        if (w->n_time_fields[1] == 0)
          _write_faces_map(b_mesh, w);
        w->n_time_fields[1] += 1;
      }
      else if (entity_dim == 3) {
        if (w->n_time_fields[0] == 0)
          _write_cells_map(b_mesh, w);
        w->n_time_fields[0] += 1;
      }
    }

  }

  /* Write phase if not already written */

  _write_phase(&phase_id, w);

  /* Write new field */

  _write_field(name,
               phase_id,
               dimension,
               datatype,
               interlace,
               n_parent_lists,
               parent_num_shift,
               field_values,
               entity_dim,
               b_mesh,
               mesh,
               w);

  /* Close file */

  if (w->rank < 1) {

    CCMIOCloseFile(err, w->root_id);

    if (error != kCCMIONoErr)
      bft_error(__FILE__, __LINE__, 0,
                _("CCMIO error %d closing file."), (int)error);

  }

  w->is_open = false;
}

/*----------------------------------------------------------------------------*/

#endif /* defined(HAVE_CCM) */

END_C_DECLS
