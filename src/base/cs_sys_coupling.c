/*============================================================================
 * SYSTEM Scale code coupling (0D/1D equations)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>
#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#if defined(HAVE_MPI)
#include "cs_coupling.h"
#endif

#include "cs_base.h"
#include "cs_field_pointer.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sys_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

const int  cs_sys_coupling_tag =
  'C'+'S'+'_'+'S'+'Y'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'I'+'N'+'G';

/*============================================================================
 *  Global variables
 *============================================================================*/

static int            _sys_n_couplings = 0;
static cs_sys_cpl_t **_sys_couplings   = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and return a new CFD<-->SYS coupled BC structre.
 *
 * \param[in] type  type of boundary condition.
 *
 * \return pointer to new cs_cfd_sys_cplbc_t structre.
 */
/*----------------------------------------------------------------------------*/

static cs_cfd_sys_cplbc_t *
_create_cfd_sys_cplbc(cs_syscpl_bc_type_t type)
{
  cs_cfd_sys_cplbc_t *retval = NULL;

  BFT_MALLOC(retval, 1, cs_cfd_sys_cplbc_t);

  retval->type = type;

  retval->input_zone_id = -1;
  retval->selection_criteria_output = NULL;

  retval->bnd_dir = 1;
  retval->surf_coeff = 1.;

  retval->n_send_fields  = 0;
  retval->send_field_ids = NULL;

  retval->n_recv_fields  = 0;
  retval->recv_field_ids = NULL;

  retval->n_sys_elts = 0;
  retval->im         = NULL;

  retval->element_name = NULL;
  retval->sys_elt_idx[0] = -1;
  retval->sys_elt_idx[1] = -1;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and return new CFD<-->SYS coupling structure.
 *
 * \param[in] sys_name      Name of the coupled SYS instance
 * \param[in] n_cpl_phases  Number of coupled fluid phases.
 *
 * \return pointer to newly create cs_sys_cpl_t structure.
 */
/*----------------------------------------------------------------------------*/

static cs_sys_cpl_t *
_create_cs_sys_coupling(const char *sys_name,
                        const int   n_cpl_phases)
{
  assert(sys_name != NULL);

  if (n_cpl_phases < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("A CFD-SYS coupling requires at least 1 coupled phase.\n"));

  cs_sys_cpl_t *cpl = NULL;

  BFT_MALLOC(cpl, 1, cs_sys_cpl_t);

  cpl->sys_name = NULL;
  BFT_MALLOC(cpl->sys_name, strlen(sys_name) + 1, char);
  strcpy(cpl->sys_name, sys_name);

  cpl->n_cpl_phases = n_cpl_phases;

  cpl->n_cpl_bcs = 0;
  cpl->cplbc     = NULL;

  cpl->n_send_vals = 0;
  cpl->send_vals   = NULL;

  cpl->n_recv_vals = 0;
  cpl->recv_vals   = NULL;

#if defined(HAVE_MPI)
  cpl->comm = MPI_COMM_NULL;
#else
  bft_error(__FILE__, __LINE__, 0,
            "Error: CFD/System scale coupling requires MPI.\n");
#endif
  cpl->cfd_root    = -1;
  cpl->sys_root    = -1;
  cpl->sys_n_ranks = 0;

  return cpl;
}

/*----------------------------------------------------------------------------
 * Initialize communicator for system coupling.
 *
 * parameters:
 *   sys_coupling  <-> System coupling structure
 *   coupling_id   <-- id of this coupling (for log messages)
 *
 *----------------------------------------------------------------------------*/

static void
_init_comm(cs_sys_cpl_t *sys_coupling,
           int           coupling_id)
{
#if defined(HAVE_MPI)

  int mpi_flag = 0;
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  MPI_Initialized(&mpi_flag);

  if (!mpi_flag)
    return;

  bft_printf("Initializing MPI coupling \"%d\"with system code \"%s\".",
             coupling_id, sys_coupling->sys_name);
  bft_printf_flush();

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  const MPI_Comm _set_comm = ple_coupling_mpi_set_get_base_comm(mpi_apps);

  ple_coupling_mpi_intracomm_create(_set_comm,
                                    cs_glob_mpi_comm,
                                    sys_coupling->sys_root,
                                    &(sys_coupling->comm),
                                    local_range,
                                    distant_range);

  bft_printf(_("[ok]\n"));
  bft_printf(_("  Local ranks = [%d..%d], distant ranks = [%d..%d].\n\n"),
             local_range[0], local_range[1] - 1,
             distant_range[0], distant_range[1] - 1);
  bft_printf_flush();

  sys_coupling->sys_n_ranks = distant_range[1] - distant_range[0];
  sys_coupling->sys_root    = distant_range[0];

  sys_coupling->cfd_root    = local_range[0];
#endif
}

/*----------------------------------------------------------------------------
 * Initialize communicator for CFD<->SYSTEM coupling
 *
 * parameters:
 *   sys_coupling  <-> System coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *   sys_root_rank <-- System code root rank
 *   n_sys_ranks   <-- Number of ranks associated with System code
 *---------------------------------------------------------------------------*/

static void
_sys_coupling_init_comm(cs_sys_cpl_t *sys_coupling,
                        int           coupling_id,
                        int           sys_root_rank,
                        int           n_sys_ranks)
{
#if defined(HAVE_MPI)

  sys_coupling->sys_root = sys_root_rank;
  sys_coupling->sys_n_ranks = n_sys_ranks;

  _init_comm(sys_coupling, coupling_id);

  /* Send data to CATHARE instance */

  /* Share number of coupled phases and zones */
  int _buff_glob[2] = {sys_coupling->n_cpl_phases,
                       sys_coupling->n_cpl_bcs};
  MPI_Send(_buff_glob, 2, MPI_INT,
           sys_coupling->sys_root, cs_sys_coupling_tag,
           sys_coupling->comm);

  /* For each coupled zone, share CATHARE element name, and:
   * type of coupled zone, number of exchanged fields, and
   * ids of inner/outer (0D) or first/last (1D) scalar cells.
   */
  int n_cpl_bcs = sys_coupling->n_cpl_bcs;
  for (int i = 0; i < n_cpl_bcs; i++) {
    cs_cfd_sys_cplbc_t *cplbc = sys_coupling->cplbc[i];

    char elt_name_buff[256] = "";
    strncpy(elt_name_buff, cplbc->element_name, 255);
    MPI_Send(elt_name_buff, 255, MPI_CHAR,
             sys_coupling->sys_root, cs_sys_coupling_tag,
             sys_coupling->comm);

    int _buff_bc[5] = {cplbc->type,
                       cplbc->n_send_fields,
                       cplbc->n_recv_fields,
                       cplbc->sys_elt_idx[0],
                       cplbc->sys_elt_idx[1]};
    MPI_Send(_buff_bc, 5, MPI_INT,
             sys_coupling->sys_root, cs_sys_coupling_tag,
             sys_coupling->comm);
  }

#endif
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI CFD<-->SYSTEM couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   n_unmatched    <->  pointer to number of unmatched couplings
 *   unmatched_ids  <->  pointer to array of unmatched couplings
 *---------------------------------------------------------------------------*/

static void
_init_all_mpi_sys(int  *n_unmatched,
                  int **unmatched_ids)
{
  int _n_unmatched = *n_unmatched;
  int *_unmatched_ids = *unmatched_ids;

  const int n_couplings = _sys_n_couplings;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* Loop on applications */
  for (int i = 0; i < n_apps; i++) {
    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);

    if (strncmp(ai.app_type, "CATHARE", 7) == 0) {
      int match_queue_id = -1;
      int coupling_id    = -1;

      if (n_apps == 2 && n_couplings == 1 && _n_unmatched == 1) {
        match_queue_id = 0;
        coupling_id = 0;
      }
      else if (ai.app_name != NULL) {
        for (int j = 0; j < _n_unmatched; j++) {
          int k = _unmatched_ids[j];
          cs_sys_cpl_t *cpl = _sys_couplings[k];
          if (strcmp(ai.app_name, cpl->sys_name) == 0) {
            coupling_id = k;
            match_queue_id = j;
            break;
          }
        }
      }

      if (coupling_id > -1) {
        _n_unmatched -= 1;
        for (int l = match_queue_id; l < _n_unmatched; l++)
          _unmatched_ids[l] = _unmatched_ids[l+1];

        if (_n_unmatched == 0)
          BFT_FREE(_unmatched_ids);

        /* Set communicator */
        _sys_coupling_init_comm(_sys_couplings[coupling_id],
                                coupling_id,
                                ai.root_rank,
                                ai.n_ranks);

        /* Print matching info */
        const char *sys_version  = cs_empty_string;
        const char *local_name   = cs_empty_string;
        const char *distant_name = cs_empty_string;

        if (ai.app_name != NULL)
          local_name = ai.app_name;
        if (ai.app_type != NULL)
          sys_version = ai.app_type;
        if (ai.app_name != NULL)
          distant_name = ai.app_name;

        bft_printf(_(" CATHARE coupling           :\n"
                     "   coupling id              : \"%d\"\n"
                     "   version                  : \"%s\"\n"
                     "   local name               : \"%s\"\n"
                     "   distant application name : \"%s\"\n"
                     "   MPI application id       : %d\n"
                     "   MPI root rank            : %d\n"
                     "   number of MPI ranks      : %d\n\n"),
                   coupling_id, sys_version, local_name, distant_name,
                   i, ai.root_rank, ai.n_ranks);

      }
    }
  } /* End loop on applications */

  bft_printf_flush();

  /* Set return values */
  *n_unmatched   = _n_unmatched;
  *unmatched_ids = _unmatched_ids;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Set coupled fields for a coupled condition depending on number of phases.
 *
 * parameters:
 *   cplbc        <-- pointer to coupled condition
 *   n_cpl_phases <-- number of coupled phases
 *---------------------------------------------------------------------------*/

static void
_sys_coupling_set_fields(cs_cfd_sys_cplbc_t *cplbc,
                         const int           n_cpl_phases)
{
  if (n_cpl_phases == 1) {
    cs_sys_cplbc_add_field_to_send(cplbc, cs_thermal_model_field()->id);
    /* CATHARE uses total pressure.
     * Use correct field according to energy model.
     * Same treatment as in _physical_property_thermal_law (cs_gui.c)
     */
    cs_field_t *field_p_tot = cs_field_by_name_try("total_pressure");
    if (field_p_tot != NULL)
      cs_sys_cplbc_add_field_to_send(cplbc, field_p_tot->id);
    else
      cs_sys_cplbc_add_field_to_send(cplbc, CS_F_(p)->id);


    // Recv
    cs_sys_cplbc_add_field_to_recv(cplbc, CS_F_(vel)->id);
    cs_sys_cplbc_add_field_to_recv(cplbc, CS_F_(vel)->id);
    cs_sys_cplbc_add_field_to_recv(cplbc, cs_thermal_model_field()->id);
    cs_sys_cplbc_add_field_to_recv(cplbc, CS_F_(p)->id);
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              "Error: Only single phase is deployed.\n");
  }
}

/*----------------------------------------------------------------------------
 * Allocate arrays used for data exchange.
 *
 * parameters:
 *   cpl <-- pointer to cfd<-->sys coupling structure
 *---------------------------------------------------------------------------*/

static void
_sys_coupling_finish_initialization(cs_sys_cpl_t *cpl)
{

  assert(cpl != NULL);


  for (int i = 0; i < cpl->n_cpl_bcs; i++) {
    cs_cfd_sys_cplbc_t *cplbc = cpl->cplbc[i];
    _sys_coupling_set_fields(cplbc, cpl->n_cpl_phases);

    cpl->n_send_vals += cplbc->n_send_fields * cplbc->n_sys_elts;
    cpl->n_recv_vals += cplbc->n_recv_fields * cplbc->n_sys_elts;
  }

  BFT_MALLOC(cpl->send_vals, cpl->n_send_vals, cs_real_t);
  memset(cpl->send_vals, 0, cpl->n_send_vals*sizeof(cs_real_t));

  BFT_MALLOC(cpl->recv_vals, cpl->n_recv_vals, cs_real_t);
  memset(cpl->recv_vals, 0, cpl->n_recv_vals*sizeof(cs_real_t));

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a cfd<-->sys coupling structure by its id
 *
 * \param[in] cpl_id  id of the requested coupling
 *
 * \return pointer to coupling structure if found, raises an error otherwise.
 */
/*----------------------------------------------------------------------------*/

cs_sys_cpl_t *
cs_sys_coupling_by_id(const int cpl_id)
{

  cs_sys_cpl_t *cpl = NULL;

  if (cpl_id < 0 || cpl_id >= _sys_n_couplings)
    bft_error(__FILE__, __LINE__, 0,
              _("Requested id \"%d\" is out of bounds [0, %d]\n"),
              cpl_id, _sys_n_couplings-1);

  cpl = _sys_couplings[cpl_id];

  return cpl;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Try getting a cfd<-->sys coupling structure by its name
 *
 * \param[in] sys_name  name of the requested coupling
 *
 * \return pointer to coupling structure if found, NULL if not found.
 */
/*----------------------------------------------------------------------------*/

cs_sys_cpl_t *
cs_sys_coupling_by_name_try(const char *sys_name)
{

  cs_sys_cpl_t *cpl = NULL;

  if (sys_name != NULL) {
    for (int cpl_id = 0; cpl_id < _sys_n_couplings; cpl_id++) {
      if (strcmp(_sys_couplings[cpl_id]->sys_name, sys_name) == 0) {
        cpl = _sys_couplings[cpl_id];
        break;
      }
    }
  }

  return cpl;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a cfd<-->sys coupling structure by its name
 *
 * \param[in] sys_name  name of the requested coupling
 *
 * \return pointer to coupling structure if found, raises and error if not found.
 */
/*----------------------------------------------------------------------------*/

cs_sys_cpl_t *
cs_sys_coupling_by_name(const char *sys_name)
{
  cs_sys_cpl_t *retval = cs_sys_coupling_by_name_try(sys_name);

  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "Error: CFD/System scale coupling \"%s\" does not exits.\n",
              sys_name);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a field to send during coupling to a given coupled BC
 *
 * \param[in] cplbc     pointer to coupled condition
 * \param[in] field_id  id of the field to send
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_cplbc_add_field_to_send(cs_cfd_sys_cplbc_t *cplbc,
                               const int           field_id)
{

  int new_id = cplbc->n_send_fields;

  cplbc->n_send_fields += 1;

  BFT_REALLOC(cplbc->send_field_ids, cplbc->n_send_fields, int);
  cplbc->send_field_ids[new_id] = field_id;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a field to recieve during coupling to a given coupled BC
 *
 * \param[in] cplbc     pointer to coupled condition
 * \param[in] field_id  id of the field to recieve
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_cplbc_add_field_to_recv(cs_cfd_sys_cplbc_t *cplbc,
                               const int           field_id)
{

  int new_id = cplbc->n_recv_fields;

  cplbc->n_recv_fields += 1;

  BFT_REALLOC(cplbc->recv_field_ids, cplbc->n_recv_fields, int);
  cplbc->recv_field_ids[new_id] = field_id;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a surface coefficient to a given coupled BC
 *
 * \param[in] cplbc  pointer to coupled condition
 * \param[in] coeff  surface coefficient to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_cplbc_define_surf_coeff(cs_cfd_sys_cplbc_t *cplbc,
                               const cs_real_t     coeff)
{
  assert(cplbc != NULL);

  cplbc->surf_coeff = coeff;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a flowrate inversion between CFD and System codes if signs
 * are inversed for a given coupled BC
 *
 * \param[in] cplbc     pointer to coupled condition
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_cplbc_inverse_bnd_dir(cs_cfd_sys_cplbc_t *cplbc)
{
  assert(cplbc != NULL);

  cplbc->bnd_dir = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a field to send/recv during coupling to a given coupled BC
 *
 * \param[in] cplbc     pointer to coupled condition
 * \param[in] dir       0 send; 1 recv
 * \param[in] field_id  id of the field to exchange
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_cplbc_add_exchanged_field(cs_cfd_sys_cplbc_t *cplbc,
                                 const int           dir,
                                 const int           field_id)
{

  assert(cplbc != NULL);

  switch(dir) {
  case 0:
    {
      int fs_id = cplbc->n_send_fields;
      cplbc->n_send_fields += 1;
      BFT_REALLOC(cplbc->send_field_ids,
                  cplbc->n_send_fields,
                  int);
      cplbc->send_field_ids[fs_id] = field_id;
      break;
  }
  case 1:
    {
      int fr_id = cplbc->n_recv_fields;
      cplbc->n_recv_fields += 1;
      BFT_REALLOC(cplbc->recv_field_ids,
                  cplbc->n_recv_fields,
                  int);
      cplbc->recv_field_ids[fr_id] = field_id;
      break;
    }
  default:
    {
      bft_error(__FILE__, __LINE__, 0,
                "Error: direction value \"%d\" is neither 0 nor 1.\n",
                dir);
      break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a coupled condition to a cfd<-->sys coupling
 *
 * \param[in] sys_coupling         pointer to cfd<->sys coupling
 * \param[in] type                 type of coupled condition
 * \param[in] z_input              coupled zone (boundary or volume)
 * \param[in] sel_criteria_output  selection criteria for cfd->sys data selection
 * \param[in] element_name         name of coupled sys element
 * \param[in] c0                   first sys cell index
 * \param[in] c1                   second sys cell index
 * \param[in] n_sys_elts           number of coupled cells in the system code
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_coupling_add_cplbc(cs_sys_cpl_t        *sys_coupling,
                          cs_syscpl_bc_type_t  type,
                          const cs_zone_t     *z_input,
                          const char          *sel_criteria_output,
                          const char          *element_name,
                          const int            c0,
                          const int            c1,
                          const int            n_sys_elts)
{

  assert(sys_coupling != NULL && element_name != NULL);

  int cpl_id = sys_coupling->n_cpl_bcs;

  sys_coupling->n_cpl_bcs += 1;

  BFT_REALLOC(sys_coupling->cplbc,
              sys_coupling->n_cpl_bcs,
              cs_cfd_sys_cplbc_t *);

  cs_cfd_sys_cplbc_t *cplbc = _create_cfd_sys_cplbc(type);

  cplbc->input_zone_id = z_input->id;
  if (sel_criteria_output != NULL) {
    size_t _l = strlen(sel_criteria_output);
    BFT_MALLOC(cplbc->selection_criteria_output, _l + 1, char);
    strcpy(cplbc->selection_criteria_output, sel_criteria_output);
  }

  if (element_name == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "Error: element name is NULL.\n");
  }
  else {
    size_t _l = strlen(element_name);
    BFT_MALLOC(cplbc->element_name, _l + 1, char);
    strncpy(cplbc->element_name, element_name, _l);
    cplbc->element_name[_l] = '\0';
  }

  cplbc->sys_elt_idx[0] = c0;
  cplbc->sys_elt_idx[1] = c1;
  cplbc->n_sys_elts = n_sys_elts;

  sys_coupling->cplbc[cpl_id] = cplbc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a cfd<->sys coupling
 *
 * \param[in] sys_name      name of the new coupling
 * \param[in] n_cpl_phases  number of phases to coupled
 *
 * \return id of the newly created coupling
 */
/*----------------------------------------------------------------------------*/

int
cs_sys_coupling_add(const char *sys_name,
                    const int   n_cpl_phases)
{

  // Check that coupling doesn't allready exist
  cs_sys_cpl_t *cpl = cs_sys_coupling_by_name_try(sys_name);
  if (cpl != NULL)
    bft_error(__FILE__, __LINE__, 0,
              "Error: coupling \"%s\" allready exists.\n",
              sys_name);

  // Create new coupling structure
  int new_id = _sys_n_couplings;
  cpl = _create_cs_sys_coupling(sys_name, n_cpl_phases);

  // Reallocate arrays
  _sys_n_couplings += 1;

  BFT_REALLOC(_sys_couplings, _sys_n_couplings, cs_sys_cpl_t *);
  _sys_couplings[new_id] = cpl;

  return new_id;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief send data to system code
 *
 * \param[in] cpl pointer to coupling structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_coupling_send_data(cs_sys_cpl_t *cpl)
{

#if defined(HAVE_MPI)
  if (cs_glob_rank_id <= 0)
    MPI_Send(cpl->send_vals,
             cpl->n_send_vals,
             CS_MPI_REAL,
             cpl->sys_root,
             cs_sys_coupling_tag,
             cpl->comm);
#endif

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief recieve data from system code
 *
 * \param[in] cpl pointer to coupling structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_coupling_recv_data(cs_sys_cpl_t *cpl)
{

#if defined(HAVE_MPI)

  MPI_Status status;
  if (cs_glob_rank_id <= 0)
    MPI_Recv(cpl->recv_vals,
             cpl->n_recv_vals,
             CS_MPI_REAL,
             cpl->sys_root,
             cs_sys_coupling_tag,
             cpl->comm,
             &status);

  // If multi CPU broadcast to other threads
  if (cs_glob_rank_id >= 0)
    MPI_Bcast(cpl->recv_vals,
              cpl->n_recv_vals,
              CS_MPI_REAL,
              0,
              cs_glob_mpi_comm);
#endif

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize cfd<->system coupling once all couplings are defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_coupling_all_init(void)
{

  // For the moment only user definition.
  // TODO: Add GUI definition
  cs_user_cathare_coupling();

  // Finalize creation
  for (int i = 0; i < _sys_n_couplings; i++) {
    cs_sys_cpl_t *cpl = _sys_couplings[i];
    _sys_coupling_finish_initialization(cpl);
  }

  /* MPI */
  int n_unmatched = _sys_n_couplings;

  int *unmatched_ids = NULL;
  BFT_MALLOC(unmatched_ids, n_unmatched, int);
  for (int i = 0; i < n_unmatched; i++)
    unmatched_ids[i] = i;

#if defined(HAVE_MPI)
  if (n_unmatched > 0)
    _init_all_mpi_sys(&n_unmatched, &unmatched_ids);
#endif

  if (n_unmatched > 0) {
    bft_printf("Unmatched SYSTEM couplings:\n"
               "---------------------------\n\n");

    BFT_FREE(unmatched_ids);
    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 SYSTEM coupling was defined for which\n"
                "no communication with a SYSTEM scale instance is possible.\n"));
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize all cfd<->sys couplings
 */
/*----------------------------------------------------------------------------*/

void
cs_sys_coupling_all_finalize(void)
{

  return;
}

END_C_DECLS
