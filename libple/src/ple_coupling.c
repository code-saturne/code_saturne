/*============================================================================
 * Set up communication with coupled codes.
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2010  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"
#include "ple_config_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "ple_coupling.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure used manage information about coupling with MPI_COMM_WORLD */

#if defined(PLE_HAVE_MPI)

struct _ple_coupling_mpi_world_t {

  int    n_apps;       /* Number of distinct applications */
  int    app_id;       /* Id of the local application in the application info */
  int    app_names_l;  /* Length of application names array */

  int   *app_info;     /* For each application, 5 integers: application number,
                          associated root, n_ranks, and indexes in app_names */
  char  *app_names;    /* Array of application type names and instance names */

};

#endif /* defined(PLE_HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of an array of strings.
 *
 * parameters:
 *   name   <-- pointer to array of names that should be ordered.
 *   level  <-- level of the binary tree to descend
 *   n_ents <-- number of entities in the binary tree to descend
 *   order  <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_names_descend_tree(const char  *name[],
                          int          level,
                          const int    n_ents,
                          int          order[])
{
  int i_save, i1, i2, lv_cur;

  i_save = order[level];

  while (level <= (n_ents/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_ents - 1) {

      i1 = order[lv_cur+1];
      i2 = order[lv_cur];

      if (strcmp(name[i1], name[i2]) > 0) lv_cur++;
    }

    if (lv_cur >= n_ents) break;

    i1 = i_save;
    i2 = order[lv_cur];

    if (strcmp(name[i1], name[i2]) >= 0) break;

    order[level] = order[lv_cur];
    level = lv_cur;
  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of strings.
 *
 * parameters:
 *   name   <-- pointer to array of names that should be ordered.
 *   order  <-- pre-allocated ordering table
 *   n_ents <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_names(const char  *name[],
             int          order[],
             const int    n_ents)
{
  int i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < n_ents ; i++)
    order[i] = i;

  if (n_ents < 2)
    return;

  /* Create binary tree */

  i = (n_ents / 2) ;
  do {
    i--;
    _order_names_descend_tree(name, i, n_ents, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n_ents - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_names_descend_tree(name, 0, i, order);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * This function builds a group id within a communicator based on its name.
 *
 * If multiple groups are present, ids are number from 0 to n_groups - 1,
 * based on the odering of group names. If all processes have the same
 * group name, the returned value is -1.
 *
 * The returned id may typically be used as a "color" argument for
 * MPI_Comm_split().
 *
 * As this function requires communication between applications, it
 * is a collective function in comm.
 *
 * parameters:
 *   comm       <-- MPI communicator.
 *   group_name <-- name associated with current group
 *
 * returns:
 *   id associated with local name.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_name_to_id(MPI_Comm     comm,
                            const char  *group_name)
{
  int i, eq_prev, eq_all;
  MPI_Status status;

  int l = 0, l_prev = 0;
  int rank_prev = MPI_PROC_NULL, rank_next = MPI_PROC_NULL;
  int rank_id = 0, n_ranks = 1, tag = 1;

  char *_group_name = NULL;
  char *buf = NULL, *names_buf = NULL;
  int *recv_count = NULL, *recv_displ = NULL, *app_id = NULL;

  int retval = -1;

  /* Initialization */

  if (_group_name == NULL) {
    l = strlen(group_name);
    PLE_MALLOC(_group_name, l + 1, char);
    strcpy(_group_name, group_name);
  }
  else {
    PLE_MALLOC(_group_name, 1, char);
    _group_name[0] = '\0';
  }

  if (comm != MPI_COMM_NULL) {
    MPI_Comm_rank(comm, &rank_id);
    MPI_Comm_size(comm, &n_ranks);
    if (rank_id > 0)
      rank_prev = rank_id -1;
    if (rank_id + 1 < n_ranks)
      rank_next = rank_id + 1;
  }

  /* Check if multiple names are present using "light" algorithm */
  /*-------------------------------------------------------------*/

  if (rank_id %2 == 0) {
    MPI_Send(&l, 1, MPI_INT, rank_next, tag, comm);
    MPI_Recv(&l_prev, 1, MPI_INT, rank_prev, tag, comm, &status);
  }
  else {
    MPI_Recv(&l_prev, 1, MPI_INT, rank_prev, tag, comm, &status);
    MPI_Send(&l, 1, MPI_INT, rank_next, tag, comm);
  }

  PLE_MALLOC(buf, l_prev + 1, char);

  if (rank_id %2 == 0) {
    MPI_Send(_group_name, l, MPI_CHAR, rank_next, tag, comm);
    MPI_Recv(buf, l_prev, MPI_CHAR, rank_prev, tag, comm, &status);
  }
  else {
    MPI_Recv(buf, l_prev, MPI_CHAR, rank_prev, tag, comm, &status);
    MPI_Send(_group_name, l, MPI_CHAR, rank_next, tag, comm);
  }

  eq_prev = 1;
  if (rank_id > 0) {
    buf[l_prev] = '\0';
    if (strcmp(_group_name, buf))
      eq_prev = 0;
  }
  MPI_Allreduce(&eq_prev, &eq_all, 1, MPI_INT, MPI_MIN, comm);

  PLE_FREE(buf);

  if (eq_all == 1) {
    PLE_FREE(_group_name);
    return -1;
  }

  /* Now gather to rank 0 for full algorithm */
  /*-----------------------------------------*/

  if (rank_id == 0) {
    PLE_MALLOC(recv_count, n_ranks, int);
    PLE_MALLOC(recv_displ, n_ranks, int);
  }

  MPI_Gather(&l, 1, MPI_INT, recv_count, 1, MPI_INT, 0, comm);

  if (rank_id == 0) {

    recv_displ[0] = 0;
    for (i = 1; i < n_ranks; i++)
      recv_displ[i] = recv_displ[i-1] + recv_count[i-1] + 1;

    PLE_MALLOC(names_buf,
               recv_displ[n_ranks - 1] + recv_count[n_ranks - 1] + 1,
               char);
  }

  MPI_Gatherv(_group_name, l, MPI_CHAR,
              names_buf, recv_count, recv_displ, MPI_CHAR,
              0, comm);

  PLE_FREE(_group_name);

  /* Order groups for rank 0 */

  if (rank_id == 0) {

    int n_apps = 1;
    int *order = NULL;
    char *name_prev = NULL;
    char **names_ptr = NULL;

    PLE_MALLOC(names_ptr, n_ranks, char *);
    for (i = 0; i < n_ranks; i++) {
      names_ptr[i] = names_buf + recv_displ[i];
      (names_ptr[i])[recv_count[i]] = '\0';
      recv_count[i] = -1;
    }

    /* Re-use arrays */
    order = recv_displ; recv_displ = NULL;
    app_id = recv_count; recv_count = NULL;

    _order_names((const char **)names_ptr, order, n_ranks);

    name_prev = names_ptr[order[0]];
    app_id[order[0]] = n_apps - 1;
    for (i = 1; i < n_ranks; i++) {
      if (strcmp(names_ptr[order[i]], name_prev)) {
        n_apps += 1;
        name_prev = names_ptr[order[i]];
      }
      app_id[order[i]] = n_apps - 1;
    }

    PLE_FREE(names_ptr);
    PLE_FREE(names_buf);
    PLE_FREE(order);
  }

  /* Now send app_id to all ranks */

  MPI_Scatter(app_id, 1, MPI_INT, &retval, 1, MPI_INT, 0, comm);

  if (rank_id == 0)
    PLE_FREE(app_id);

  return retval;
}

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI_COMM_WORLD.
 *
 * The application communicator app_comm is usually obtained from
 * MPI_COMM_WORLD using MPI_Comm_split, with app_num corresponding to
 * the "color" argument in that function.
 *
 * As this function requires communication between applications, it
 * is a collective function in MPI_COMM_WORLD.
 *
 * parameters:
 *   app_num   <-- application number in MPI_COMM_WORLD (nonnegative).
 *   app_name  <-- name of current application.
 *   case_name <-- name of current case, or NULL.
 *   app_comm  <-- communicator associated with local application.
 *
 * returns:
 *   PLE coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_world_t *
ple_coupling_mpi_world_create(int          app_num,
                              const char  *app_type,
                              const char  *app_name,
                              MPI_Comm     app_comm)
{
  int i;
  MPI_Status status;

  int world_rank = -1, app_rank = -1, n_app_ranks = 0;
  int root_marker = 0;
  int info_tag = 1, name_tag = 2;

  int counter[2] = {0, 0};
  int l_rank_info[5] = {-1, -1, -1, 1, 1};
  int *rank_info = NULL;
  char *app_names = NULL;

  ple_coupling_mpi_world_t *w = NULL;

  /* Initialization */

  PLE_MALLOC(w, 1, ple_coupling_mpi_world_t);

  w->n_apps = 0;
  w->app_id = -1;
  w->app_names_l = 0;
  w->app_info = NULL;
  w->app_names = NULL;

  /* Initialization */

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (app_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(app_comm, &app_rank);
    MPI_Comm_size(app_comm, &n_app_ranks);
  }
  else {
    app_rank = 0;
    n_app_ranks = 1;
  }

  l_rank_info[0] = app_num;
  l_rank_info[1] = world_rank;
  l_rank_info[2] = n_app_ranks;
  if (app_type != NULL)
    l_rank_info[3] = strlen(app_type) + 1;
  if (app_name != NULL)
    l_rank_info[4] = strlen(app_name) + 1;

  if (app_rank == 0)
    root_marker = 1;

  /* Root rank of MPI_COMM_WORLD counts applications and receives info */

  MPI_Reduce(&root_marker, &(counter[0]), 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);

  /* Root of MPI_COMM_WORLD collects all info */

  if (world_rank == 0) {

    int start = 0;

    PLE_MALLOC(rank_info, counter[0]*5, int);

    if (app_rank == 0) {
      for (i = 0; i < 5; i++)
        rank_info[i] = l_rank_info[i];
      start = 1;
    }

    /* Use of different tags for info and strings is important
       here as we use MPI_ANY_SOURCE and messages could be mixed */

    for (i = start; i < counter[0]; i++)
      MPI_Recv(rank_info + i*5, 5, MPI_INT, MPI_ANY_SOURCE, info_tag,
               MPI_COMM_WORLD, &status);

    /* Convert rank_info count to index values */

    for (i = 0; i < counter[0]; i++)
      counter[1] += (rank_info[i*5 + 3] + rank_info[i*5 + 4]);

    PLE_MALLOC(app_names, counter[1], char);
    memset(app_names, 0, counter[1]);

    counter[1] = 0;

    if (app_rank == 0) {
      strcpy(app_names, app_type);
      if (app_name != NULL)
        strcpy(app_names + rank_info[3], app_name);
      else
        app_names[rank_info[3]] = '\0';
      counter[1] += (rank_info[3] + rank_info[4]);
      rank_info[4] = rank_info[3];
      rank_info[3] = 0;
    }

    for (i = start; i < counter[0]; i++) {
      int app_type_size = rank_info[i*5 + 3];
      int app_name_size = rank_info[i*5 + 4];
      int msg_len = app_type_size + app_name_size;
      rank_info[i*5 + 3] = counter[1];
      rank_info[i*5 + 4] = counter[1] + app_type_size;
      MPI_Recv(app_names + counter[1], msg_len, MPI_CHAR, rank_info[i*5 +1],
               name_tag, MPI_COMM_WORLD, &status);
      counter[1] += msg_len;
    }

  }

  /* Other root ranks send info to root */

  else if (app_rank == 0) { /* world_rank != 0 */

    char *sendbuf = NULL;
    int   sendbuf_l = l_rank_info[3] + l_rank_info[4];

    PLE_MALLOC(sendbuf, sendbuf_l, char);

    if (app_type != NULL)
      strcpy(sendbuf, app_type);
    else
      sendbuf[0] = '\0';
    if (app_name != NULL)
      strcpy(sendbuf + l_rank_info[3], app_name);
    else
      sendbuf[l_rank_info[3]] = '\0';

    MPI_Send(l_rank_info, 5, MPI_INT, 0, info_tag, MPI_COMM_WORLD);
    MPI_Send(sendbuf, sendbuf_l, MPI_CHAR, 0, name_tag, MPI_COMM_WORLD);

    PLE_FREE(sendbuf);
  }

  /* Now root broadcasts application info */

  MPI_Bcast(counter, 2, MPI_INT, 0, MPI_COMM_WORLD);

  if (world_rank != 0) {
    PLE_MALLOC(rank_info, counter[0]*5, int);
    PLE_MALLOC(app_names, counter[1], char);
  }

  MPI_Bcast(rank_info, counter[0]*5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(app_names, counter[1], MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Set global values */

  w->n_apps = counter[0];
  w->app_names_l = counter[1];
  w->app_info = rank_info;
  w->app_names = app_names;

  for (i = 0; i < w->n_apps && w->app_id < 0; i++) {
    if (w->app_info[i*5] == app_num)
      w->app_id = i;
  }

  return w;
}

/*----------------------------------------------------------------------------
 * Free an PLE coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-> pointer to structure that should be freed.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_world_destroy(ple_coupling_mpi_world_t **w)
{
  ple_coupling_mpi_world_t *_w = *w;

  if (_w != NULL) {
    PLE_FREE(_w->app_info);
    PLE_FREE(_w->app_names);
    PLE_FREE(*w);
  }
}

/*----------------------------------------------------------------------------
 * Return the number of applications in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   number of application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_world_n_apps(const ple_coupling_mpi_world_t  *w)
{
  int retval = 0;

  if (w != NULL)
    retval = w->n_apps;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the id of the local application in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   id of the local application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_world_get_app_id(const ple_coupling_mpi_world_t  *w)
{
  int retval = -1;

  if (w != NULL)
    retval = w->app_id;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return application information in MPI_COMM_WORLD.
 *
 * parameters:
 *   w      <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *   app_id <-- application id
 *
 * returns:
 *   application information structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_world_info_t
ple_coupling_mpi_world_get_info(const ple_coupling_mpi_world_t  *w,
                                int                              app_id)
{
  ple_coupling_mpi_world_info_t  retval;

  retval.app_num = -1;
  retval.root_rank = -1;
  retval.n_ranks = 0;
  retval.app_type = NULL;
  retval.app_name = NULL;

  if (w != NULL) {
    if (app_id >= 0 && app_id < w->n_apps) {
      retval.app_num = w->app_info[app_id*5];
      retval.root_rank = w->app_info[app_id*5 + 1];
      retval.n_ranks = w->app_info[app_id*5 + 2];
      retval.app_type = w->app_names + w->app_info[app_id*5 + 3];
      retval.app_name = w->app_names + w->app_info[app_id*5 + 4];
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Create an intracommunicator from a local and distant communicator
 * within MPI_COMM_WORLD.
 *
 * parameters:
 *   app_comm      <-- communicator associated with local application
 *   distant_root  <-- rank of distant group leader in MPI_COMM_WORLD
 *   new_comm      --> pointer to new communicator
 *   local_range   --> first and past-the last ranks of local application
 *                     in new communicator
 *   distant_range --> first and past-the last ranks of distant application
 *                     in new communicator
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_intracomm_create(MPI_Comm   app_comm,
                                  int        distant_root,
                                  MPI_Comm  *new_comm,
                                  int        local_range[2],
                                  int        distant_range[2])
{
  int coupling_tag = ('P'+'L'+'E'+'_'+'C'+'O'+'U'+'P'+'L'+'I'+'N'+'G') % 512;
  int  mpi_flag = 0;
  int  n_dist_ranks = 0;
  int  n_loc_ranks, r_glob, r_loc_max;
  MPI_Comm  intercomm_tmp;
  int  r_coupl, r_coupl_min;
  int  high = 1;

  /* Initialization */

  *new_comm = MPI_COMM_NULL;

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    return;

  MPI_Comm_rank(MPI_COMM_WORLD, &r_glob);

  MPI_Allreduce(&r_glob, &r_loc_max, 1, MPI_INT, MPI_MAX, app_comm);

  if (distant_root > r_loc_max)
    high = 0;

  MPI_Comm_size(app_comm, &n_loc_ranks);

  /* Create a reserved communicator */

  MPI_Intercomm_create(app_comm, 0, MPI_COMM_WORLD,
                       distant_root, coupling_tag, &intercomm_tmp);

  MPI_Intercomm_merge(intercomm_tmp, high, new_comm);

  MPI_Comm_free(&intercomm_tmp);

  /* Compute number of distant ranks and first distant rank */

  MPI_Comm_size(*new_comm, &n_dist_ranks);
  n_dist_ranks -= n_loc_ranks;

  /* Check rank in new communicator (should not be necessary with correctly
     set "high" value, but seems to be with Open MPI 1.0.1) */

  MPI_Comm_rank(*new_comm, &r_coupl);
  MPI_Allreduce(&r_coupl, &r_coupl_min, 1, MPI_INT, MPI_MIN, app_comm);
  high = (r_coupl_min == 0) ? 0 : 1;

  /* Deduce the position of the first distant rank in the new communicator */

  if (high == 0) {
    local_range[0] = 0;
    distant_range[0] = n_loc_ranks;
  }
  else {
    local_range[0] = n_dist_ranks;
    distant_range[0] = 0;
  }

  local_range[1] = local_range[0] + n_loc_ranks;
  distant_range[1] = distant_range[0] + n_dist_ranks;
}

/*----------------------------------------------------------------------------
 * Dump printout of an PLE coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_world_dump(const ple_coupling_mpi_world_t  *w)
{
  int i;

  if (w == NULL) {
    ple_printf("  Coupling MPI_COMM_WORLD info: nil\n");
    return;
  }

  ple_printf("  Coupling MPI_COMM_WORLD info: %p\n"
             "    number of applications:     %d\n"
             "    local application id:       %d\n"
             "    app_names_size:             %d\n\n",
             w, w->n_apps, w->app_id, w->app_names_l);

  for (i = 0; i < w->n_apps; i++)
    ple_printf("    Application number:  %d\n"
               "      root_rank:         %d\n"
               "      n_ranks:           %d\n"
               "      app_type:          \"%s\"\n"
               "      app_name:          \"%s\"\n\n",
               w->app_info[i*5], w->app_info[i*5 + 1], w->app_info[i*5 + 2],
               w->app_names + w->app_info[i*5 + 3],
               w->app_names + w->app_info[i*5 + 4]);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

