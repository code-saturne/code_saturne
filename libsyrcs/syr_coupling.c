/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

/*===========================================================================
 * Main API functions for coupling between Syrthes 3 and Code_Saturne
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "syr_defs.h"
#include "syr_comm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "syr_coupling.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Auxiliary structure to determine how data is distributed across
   distant ranks */

typedef struct {

  int       n_ranks;    /* Number of ranks containing data */
  int       n_elts_g;   /* Global number of elements */
  int       n_elts_max; /* Maximum local number of elements */

  int      *proc_id;    /* Proc id (global rank) for each local rank */
  int      *index;      /* Start index for each rank (n_ranks + 1) */
  int      *elt_num;    /* Global number (1 to n_elts_g) associated with each
                           element */

} syr_coupling_dist_t;

struct _syr_coupling_t {

  int           comm_echo;    /* Optional echo to standard output */

  syr_comm_t   *comm;         /* Communicator */

  int          *cs_rank;      /* Array of ranks in Code_Saturne */

  syr_coupling_dist_t  dist;  /* Distribution structure for variables */

};

/*===========================================================================
 * Static global variables
 *===========================================================================*/

static char  syr_glob_coupling_sec_name_error[]
               = "Erreur dans la communication : \"%s\".\n"
                 "Nom de rubrique <%s> recue du rang <%d>\n"
                 "est different de <%s> recue du rang <1>";

/*===========================================================================
 * Private function definitions
 *===========================================================================*/

/*---------------------------------------------------------------------------
 * Finalize syr_coupling_dist_t structure
 *---------------------------------------------------------------------------*/

static void
_syr_coupling_dist_finalize(syr_coupling_dist_t  *dist)
{
  if (dist != NULL) {

    if (dist->n_elts_g != 0) {
      dist->n_elts_g = 0;
      dist->n_elts_max = 0;
      PLE_FREE(dist->proc_id);
      PLE_FREE(dist->index);
      PLE_FREE(dist->elt_num);
    }

  }
}

/*---------------------------------------------------------------------------
 * Gather variable values from a communication buffer
 *---------------------------------------------------------------------------*/

static void
_syr_coupling_gather_var(const syr_coupling_t  *coupling,
                         int                    proc_id,
                         const double           buffer[],
                         double                 var[])
{
  /* Gather values */

  int ii;
  const syr_coupling_dist_t *dist = &(coupling->dist);
  int n_elts = dist->index[proc_id+1] - dist->index[proc_id];

  if (dist->elt_num != NULL) {

    const int *elt_num = dist->elt_num + dist->index[proc_id];
    for (ii = 0; ii < n_elts; ii++)
      var[elt_num[ii] - 1] = buffer[ii];

  }
  else {

    for (ii = 0; ii < n_elts; ii++)
      var[dist->index[proc_id] + ii] = buffer[ii];

  }
}

/*---------------------------------------------------------------------------
 * Scatter variable values to a communication buffer
 *---------------------------------------------------------------------------*/

static void
_syr_coupling_scatter_var(const syr_coupling_t  *coupling,
                          int                    proc_id,
                          double                 buffer[],
                          const double           var[])
{
  /* Scatter values */

  int ii;
  const syr_coupling_dist_t *dist = &(coupling->dist);
  int n_elts = dist->index[proc_id+1] - dist->index[proc_id];

  if (dist->elt_num != NULL) {

    const int *elt_num = dist->elt_num + dist->index[proc_id];
    for (ii = 0; ii < n_elts; ii++)
      buffer[ii] = var[elt_num[ii] - 1];

  }
  else {

    for (ii = 0; ii < n_elts; ii++)
      buffer[ii] = var[dist->index[proc_id] + ii];

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize syr_coupling_t structure
 *
 * arguments:
 *   coupling_id <-- Id of Syrthes coupling (0 to n-1)
 *   cs_app_name <-- Application name of Code_Saturne MPI process
 *   comm_echo   <-- Optional echo to standard output
 *----------------------------------------------------------------------------*/

syr_coupling_t  *
syr_coupling_initialize(int               coupling_id,
                        const char       *cs_app_name,
                        int               comm_echo)
{
  int root_rank = -1, n_ranks = -1;

  syr_coupling_t *coupling = NULL;

  /* Allocate syr_coupling_t structure */

  PLE_MALLOC(coupling, 1, syr_coupling_t);

  /* Default initialization */

  coupling->comm_echo = comm_echo;

  coupling->cs_rank = NULL;

  syr_mpi_appinfo(cs_app_name, &root_rank, &n_ranks);

  /* Initialize communicator */

  coupling->comm = NULL;

  coupling->comm = syr_comm_initialize(coupling_id + 1,
                                       root_rank,
                                       n_ranks,
                                       comm_echo);

  /* Variable distribution will be defined later */

  coupling->dist.n_ranks = 0;
  coupling->dist.n_elts_g = 0;
  coupling->dist.n_elts_max = 0;
  coupling->dist.proc_id = NULL;
  coupling->dist.index = NULL;
  coupling->dist.elt_num = NULL;

  return coupling;
}

/*---------------------------------------------------------------------------
 * Finalize syr_coupling_t structure
 *---------------------------------------------------------------------------*/

syr_coupling_t *
syr_coupling_finalize(syr_coupling_t  *coupling)
{
  /* Finalize send and receive communicators */

  coupling->comm = syr_comm_finalize(coupling->comm);

  if (coupling->cs_rank != NULL)
    PLE_FREE(coupling->cs_rank);

  /* Finalize distribution */

  _syr_coupling_dist_finalize(&(coupling->dist));

  PLE_FREE(coupling);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Receive Code_Saturne's coupled surface mesh
 *
 * the elt_connect and vtx_coords arrays are allocated here.
 *
 * arguments:
 *   coupling    <-- Associated coupling object
 *   sp_dim      --> Spatial dimension for fluid
 *   n_vtx,      --> Number of vertices
 *   n_elts,     --> Number of segments or triangles
 *   vtx_coords  --> Vertex coordinates
 *   elt_connect --> Segment or triangle connectivity
 *----------------------------------------------------------------------------*/

void
syr_coupling_receive_bc_mesh(syr_coupling_t  *coupling,
                             int             *sp_dim,
                             int             *n_vtx,
                             int             *n_elts,
                             double          *vtx_coords[],
                             int             *elt_connect[])
{
  char sec_name[SYR_COMM_L_SEC_NAME + 1];
  char sec_name_cmp[SYR_COMM_L_SEC_NAME + 1];
  int  ii, jj, kk;
  int  proc_id;
  syr_type_t elt_type;

  int n_ranks = 0;
  int stride = 0;
  int max_msg_size = 0;
  int msg_size = 0;
  syr_coupling_dist_t *dist = NULL;
  syr_coupling_dist_t  elt_dist;

  const syr_comm_t  *comm  = coupling->comm;

  /* Temporary distribution for elements */

  elt_dist.n_ranks = 0;
  elt_dist.n_elts_g = 0;
  elt_dist.n_elts_max = 0;
  elt_dist.proc_id = NULL;
  elt_dist.index = NULL;
  elt_dist.elt_num = NULL;

  /* Number of ranks of coupled Code_Saturne processes */

  n_ranks = syr_comm_get_n_procs(coupling->comm);

  /* Default initializations */

  *sp_dim = 3;

  *n_elts = 0;

  *vtx_coords = NULL;
  *elt_connect = NULL;

  /* Loop on messages */

  while (1) {

    /* Read message header */
    /*---------------------*/

    syr_comm_read_header(sec_name,
                         &msg_size,
                         &elt_type,
                         comm,
                         0);

    /* Partial decoding so as to handle messages by rank 0 only */
    /*----------------------------------------------------------*/

    /* End of initialization data (try this first) */

    if (!strncmp("coupl:b:start", sec_name, strlen("coupl:b:start")))
      break;

    /* Spatial dimension */

    else if (!strncmp("coupl:b:ndim_", sec_name, strlen("coupl:b:ndim_"))) {
      syr_comm_read_body(1,
                         sp_dim,
                         elt_type,
                         comm,
                         0);
      continue;
    }

    /* Other messages should be received from all ranks, so read headers */
    /*-------------------------------------------------------------------*/

    strncpy(sec_name_cmp, sec_name, SYR_COMM_L_SEC_NAME);
    sec_name_cmp[SYR_COMM_L_SEC_NAME] = '\0';
    max_msg_size = msg_size;

    for (proc_id = 1; proc_id < n_ranks; proc_id++) {

      syr_comm_read_header(sec_name,
                           &msg_size,
                           &elt_type,
                           comm,
                           proc_id);

      /* Sanity test */

      if (strncmp(sec_name_cmp, sec_name, SYR_COMM_L_SEC_NAME))
        ple_error(__FILE__, __LINE__, 0,
                  syr_glob_coupling_sec_name_error,
                  syr_comm_get_name(comm), sec_name,
                  proc_id + 1, sec_name_cmp);

      /* Number of elements in body */

      max_msg_size = SYR_MAX(msg_size, max_msg_size);

    }

    /* Decode headers received from all ranks */
    /*-----------------------------------------*/

    /* Read the local number of vertices or elements */
    /*-----------------------------------------------*/

    if (   !strncmp("coupl:b:npoinf", sec_name, strlen("coupl:b:npoinf"))
        || !strncmp("coupl:b:nelebf", sec_name, strlen("coupl:b:nelebf"))) {

      int rank_count = 0;
      int *_n_elts = NULL;

      PLE_MALLOC(_n_elts, n_ranks, int);

      /* Set distribution based on vertices or elements */

      if (!strncmp("coupl:b:npoinf", sec_name, strlen("coupl:b:npoinf")))
        dist = &(coupling->dist);
      else
        dist = &elt_dist;

      /* Receive data */

      for (proc_id = 0; proc_id < n_ranks ; proc_id++) {

        syr_comm_read_body(1,
                           &(_n_elts[proc_id]),
                           elt_type,
                           comm,
                           proc_id);

      }

      /* Reduce n_ranks to those ranks owning at least one element */

      for (proc_id = 0; proc_id < n_ranks ; proc_id++) {
        if (_n_elts[proc_id] > 0)
          rank_count += 1;
      }

      dist->n_ranks = rank_count;

      /* Build the list of coupled ranks, and associated index */

      if (dist->proc_id != NULL) {/* in case of re-definition */
        PLE_FREE(dist->proc_id);
        PLE_FREE(dist->index);
      }

      PLE_MALLOC(dist->proc_id, dist->n_ranks, int);
      PLE_MALLOC(dist->index, dist->n_ranks + 1, int);

      dist->index[0] = 0;
      rank_count = 0;

      for (proc_id = 0; proc_id < n_ranks ; proc_id++) {

        if (_n_elts[proc_id] > 0) {

          if (_n_elts[proc_id] > dist->n_elts_max)
            dist->n_elts_max = _n_elts[proc_id];

          dist->proc_id[rank_count] = proc_id;
          dist->index[rank_count + 1] = (  dist->index[rank_count]
                                         + _n_elts[proc_id]);
          rank_count += 1;

        }

      }

      /* If there is only one distant rank, the local and
         global number of vertices are the same, and the latter
         might not be re-sent later; for elements, we can
         always compute the global number of elements, as
         they do not overlap on processor boundaries */

      if (n_ranks == 1)
        dist->n_elts_g = _n_elts[0];

      if (dist == &elt_dist)
        dist->n_elts_g = dist->index[dist->n_ranks];

      PLE_FREE(_n_elts);

    } /* End reading the local number of entities */

    /* Read the global number of vertices */
    /*------------------------------------*/

    else if (!strncmp("coupl:b:g:npoinf",
                      sec_name,
                      strlen("coupl:b:g:npoinf"))) {

      for (proc_id = 0; proc_id < n_ranks; proc_id++) {

        syr_comm_read_body(1,
                           &(coupling->dist.n_elts_g),
                           elt_type,
                           comm,
                           proc_id);

      }

    }

    /* Read the global numbering of vertices or elements */
    /*---------------------------------------------------*/

    else if (   !strncmp("coupl:b:g:vtxnum",
                         sec_name,
                         strlen("coupl:b:g:vtxnum"))
             || !strncmp("coupl:b:g:eltnum",
                         sec_name,
                         strlen("coupl:b:g:eltnum"))) {

      /* Set distribution based on vertices or elements */

      if (!strncmp("coupl:b:g:vtxnum", sec_name, strlen("coupl:b:g:vtxnum")))
        dist = &(coupling->dist);
      else
        dist = &elt_dist;

      /* Receive data */

      if (dist->elt_num != NULL) /* in case of re-definition */
        PLE_FREE(dist->elt_num);

      PLE_MALLOC(dist->elt_num, dist->index[dist->n_ranks], int);

      for (ii = 0; ii < dist->n_ranks; ii++) {

        syr_comm_read_body((dist->index[ii+1] - dist->index[ii]),
                           (dist->elt_num + dist->index[ii]),
                           elt_type,
                           comm,
                           dist->proc_id[ii]);

      }

    }

    /* Read the vertices's coordinates (non-interlaced) */
    /*--------------------------------------------------*/

    else if (!strncmp("coupl:b:xyzf", sec_name, strlen("coupl:b:xyzf") )) {

      double *buffer = NULL;

      dist = &(coupling->dist);
      stride = *sp_dim;

      /* Allocate necessary memory then receive data */

      PLE_MALLOC(*vtx_coords, dist->n_elts_g * stride, double);

      if (dist->n_ranks > 1)
        PLE_MALLOC(buffer, max_msg_size * stride, double);
      else
        buffer = *vtx_coords;

      for (ii = 0; ii < dist->n_ranks; ii++) {

        int glob_id;
        int n_elts_loc = dist->index[ii+1] - dist->index[ii];

        syr_comm_read_body(n_elts_loc * stride,
                           buffer,
                           elt_type,
                           comm,
                           dist->proc_id[ii]);

        /* Assemble coordinates if necessary */

        if (buffer != *vtx_coords) {
          for (jj = 0; jj < n_elts_loc; jj++) {
            if (dist->elt_num != NULL)
              glob_id = dist->elt_num[dist->index[ii] + jj] - 1;
            else
              glob_id = dist->index[ii] + jj;

            for (kk = 0; kk < stride; kk++) {
              (*vtx_coords)[(dist->n_elts_g * kk) + glob_id]
                = buffer[(n_elts_loc * kk) + jj];
            }
          }
        }

      } /* End of loop on distant ranks */

      if (buffer != *vtx_coords) {
        PLE_FREE(buffer);
      }

      *n_vtx = dist->n_elts_g;

    }

    /* Read segment or triangle connectivity (non-interlaced) */
    /*--------------------------------------------------------*/

    else if (!strncmp("coupl:b:nodebf", sec_name, strlen("coupl:b:nodebf"))) {

      int *buffer = NULL;
      syr_coupling_dist_t *vtx_dist = &(coupling->dist);

      stride = *sp_dim; /* 2 for 2d segments, 3 for 3d triangles */
      dist = &elt_dist;

      /* Allocate necessary memory then receive data */

      PLE_MALLOC(*elt_connect, dist->n_elts_g * stride, int);

      if (dist->n_ranks > 1)
        PLE_MALLOC(buffer, max_msg_size * stride, int);
      else
        buffer = *elt_connect;

      for (ii = 0; ii < dist->n_ranks; ii++) {

        int elt_id, vtx_id, global_vtx_num;
        int n_elts_loc = dist->index[ii+1] - dist->index[ii];

        syr_comm_read_body(n_elts_loc * stride,
                           buffer,
                           elt_type,
                           comm,
                           dist->proc_id[ii]);

        /* Assemble connectivity if necessary */

        if (buffer != *elt_connect) {
          for (jj = 0; jj < n_elts_loc; jj++) {
            if (dist->elt_num != NULL)
              elt_id = dist->elt_num[dist->index[ii] + jj] - 1;
            else
              elt_id = dist->index[ii] + jj;

            assert( elt_id != -1 && elt_id <= dist->n_elts_g );

            for (kk = 0; kk < stride; kk++) {

              vtx_id = buffer[(n_elts_loc * kk) + jj] - 1;

              assert(vtx_id != -1);

              global_vtx_num = vtx_dist->elt_num[vtx_dist->index[ii] + vtx_id];

              assert(global_vtx_num != 0);
              assert(global_vtx_num <= vtx_dist->n_elts_g);

              (*elt_connect)[(dist->n_elts_g * kk) + elt_id] = global_vtx_num;


            } /* End of loop on stride */

          } /* End of loop on local elements */

        } /* End if buffer != *elt_connect */

      } /* End of loop on distant ranks */

      if (buffer != (*elt_connect)) {
        PLE_FREE(buffer);
      }

      /* Set return values */

      *n_elts = dist->n_elts_g;

    }

    /* Otherwise unknown or unexpected message */
    /*-----------------------------------------*/

    else
      ple_error(__FILE__, __LINE__, 0,
                "Message \"%s\" inconnu ou inattendu a cette etape :\n"
                "--> abandon.",
                sec_name);

  } /* End of loop on messages */

  /* Free values that will be needed no more */

  _syr_coupling_dist_finalize(&elt_dist);

  /* Check all data was received */

  if (*n_vtx == 0 || *vtx_coords == NULL)
    ple_error(__FILE__, __LINE__, 0,
              "Aucune donnee sur les sommets n'a ete recue.");

  if (*n_elts == 0 || *elt_connect == NULL)
    ple_error(__FILE__, __LINE__, 0,
              "Aucune donnee sur les elements n'a ete recue.");

}

/*----------------------------------------------------------------------------
 * Exchange of synchronization (supervision) messages
 *
 * parameters:
 *  coupling <-- Associated coupling object
 *  is_last  --> Last time step or iteration indicator
 *  is_end   --> Calculation stop indicator
 *----------------------------------------------------------------------------*/

void
syr_coupling_supervise(syr_coupling_t  *coupling,
                       int             *is_last,
                       int             *is_end)
{
  char  sec_name[SYR_COMM_L_SEC_NAME + 1];

  syr_type_t elt_type;

  int msg_size = 0;

  const int comm_echo = coupling->comm_echo;
  const syr_comm_t *comm = coupling->comm;

  /* Send a command message to Code_Saturne */
  /*----------------------------------------*/

  if (*is_end == 1)
    syr_comm_write_section("cmd:stop",
                           0,
                           NULL,
                           SYR_TYPE_void,
                           comm,
                           0);

  else /* Ready to start a new iteration */
    syr_comm_write_section("cmd:iter:start",
                           0,
                           NULL,
                           SYR_TYPE_void,
                           comm,
                           0);

  /* Receive a command message from Code_Saturne */
  /*---------------------------------------------*/

  syr_comm_read_header(sec_name,
                       &msg_size,
                       &elt_type,
                       comm,
                       0);

  /* If message indicates calculation stop */
  /*---------------------------------------*/

  if (   !strncmp("EOF", sec_name, strlen("EOF"))
      || !strncmp("cmd:stop", sec_name, strlen("cmd:stop"))) {

    printf("\txxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
           "\tx  Couplage arrete par Code_Saturne  x\n"
           "\txxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    fflush(stdout);

    *is_end = 1;

  }

  /* If message indicates last iteration */
  /*-------------------------------------*/

  else if (!strncmp("cmd:iter:start:last", sec_name,
                    strlen("cmd:iter:start:last"))) {

    if (comm_echo >= 0) {
      printf("\t*** Code_Saturne indique une derniere iteration.\n\n");
      fflush(stdout);
    }

    *is_last = 1;

  }

  /* If message indicates a new iteration */
  /*--------------------------------------*/

  else if (!strncmp("cmd:iter:start", sec_name, strlen("cmd:iter:start"))) {

    if (comm_echo >= 0) {
      printf("\t*** Code_Saturne indique une nouvelle iteration.\n\n");
      fflush(stdout);
    }

  }

  /* Otherwise, message is unknown */
  /*-------------------------------*/

  else
    ple_error(__FILE__, __LINE__, 0,
              "Message \"%s\" inconnu ou inattendu a cette etape :\n"
              "--> abandon.",
              sec_name);

  assert(msg_size == 0);
}

/*----------------------------------------------------------------------------
 * Data exchange prior to iteration
 *
 * Send wall temperature
 * Receive fluid temperature and pseudo-exchange coefficient
 * Possibly receive CFD code time step
 *
 * parameters:
 *   coupling <-- Associated coupling object
 *   tpf      <-> Wall Temperature in, Fluid temperature out
 *   hht      --> Pseudo-exchange coefficient
 *   dtfluid  --> Time step set by CFD code if received, -1 otherwise
 *----------------------------------------------------------------------------*/

void
syr_coupling_exchange_var(syr_coupling_t  *coupling,
                          double          *tpf,
                          double          *hht,
                          double          *dtfluid)
{
  char sec_name[SYR_COMM_L_SEC_NAME + 1];
  char sec_name_cmp[SYR_COMM_L_SEC_NAME + 1];

  int ii;
  int r_tf = 0, r_hht = 0;
  int n_elts = 0;
  syr_type_t elt_type;
  double *buffer = NULL, *recv_var = NULL;
  syr_coupling_dist_t *dist = NULL;

  const syr_comm_t *comm = coupling->comm;

  /* Initialize dtfluid in case the corresponding information
     is not set by the CFD code */

  *dtfluid = -1;

  /* Send wall temperature calculated by Syrthes */
  /*---------------------------------------------*/

  /* Define section name */

  sprintf(sec_name, "coupl:b:tparoi");

  /* Loop on communicating ranks */

  dist = &(coupling->dist);

  PLE_MALLOC(buffer, dist->n_elts_max, double);

  for (ii = 0; ii < dist->n_ranks; ii++) {

    /* Scatter to local buffer then write */

    _syr_coupling_scatter_var(coupling,
                              ii,
                              buffer,
                              tpf);

    syr_comm_write_section(sec_name,
                           (dist->index[ii + 1] - dist->index[ii]),
                           (void *)buffer,
                           SYR_TYPE_double,
                           comm,
                           dist->proc_id[ii]);

  } /* End of loop on communicating ranks */

  /*-------------------------------------------------------------------------*/
  /* Receive fluid temperature and exchange coefficient computed by the      */
  /* Code_Saturne kernel(s)                                                  */
  /*-------------------------------------------------------------------------*/

  while (r_tf == 0 || r_hht == 0) {

    n_elts = 0;
    recv_var = NULL;

    /* Read message header */
    /*-------------------- */

    for (ii = 0; ii < dist->n_ranks; ii++) {

      syr_comm_read_header(sec_name,
                           &n_elts,
                           &elt_type,
                           comm,
                           dist->proc_id[ii]);

      /* Optional time step message may have been inserted first */

      if (!strncmp("coupl:dtfluid:",
                   sec_name,
                   strlen("coupl:dtfluid:"))) {

        assert(n_elts == 1);
        assert(elt_type == SYR_TYPE_double);

        syr_comm_read_body(1,
                           dtfluid,
                           elt_type,
                           comm,
                           dist->proc_id[ii]);

        /* Read next header */

        syr_comm_read_header(sec_name,
                             &n_elts,
                             &elt_type,
                             comm,
                             dist->proc_id[ii]);
      }

      if (ii == 0) {

        if (!strncmp("coupl:b:tfluid",
                     sec_name,
                     strlen("coupl:b:tfluid"))) {
          recv_var = tpf;
          r_tf = 1;
        }
        else if (!strncmp("coupl:b:hparoi",
                          sec_name,
                          strlen("coupl:b:hparoi"))) {
          recv_var = hht;
          r_hht = 1;
        }
        else
          ple_error(__FILE__, __LINE__, 0,
                    "Message \"%s\" inconnu ou inattendu a cette etape",
                    sec_name);

        strncpy(sec_name_cmp, sec_name, SYR_COMM_L_SEC_NAME);
        sec_name_cmp[SYR_COMM_L_SEC_NAME] = '\0';

      }

      /* Sanity test */

      else if (   ii > 0
               && strncmp(sec_name_cmp, sec_name, SYR_COMM_L_SEC_NAME))
        ple_error(__FILE__, __LINE__, 0,
                  syr_glob_coupling_sec_name_error,
                  syr_comm_get_name(comm), sec_name, ii + 1,
                  sec_name_cmp);

      /* Receive fluid temperature or exchange coefficient */
      /*-------------------------------------------------- */

      if (recv_var != NULL) {

        assert(n_elts <= dist->n_elts_max);

        syr_comm_read_body(n_elts,
                           buffer,
                           elt_type,
                           comm,
                           dist->proc_id[ii]);

        _syr_coupling_gather_var(coupling,
                                 ii,
                                 buffer,
                                 recv_var);

      }

    }

  }

  PLE_FREE(buffer);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
