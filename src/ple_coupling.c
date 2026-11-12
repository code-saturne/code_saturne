/*============================================================================
 * Set up communication with coupled codes.
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2024  EDF S.A.

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

/*!
 * \file ple_coupling.c
 *
 * \brief Set up communication with coupled codes.
 */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdint.h>
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

/* Structure used to manage information about coupling */

#if defined(PLE_HAVE_MPI)

struct _ple_coupling_mpi_set_t {

  int    n_apps;       /* Number of distinct applications */
  int    app_id;       /* Id of the local application in the application info */
  int    app_names_l;  /* Length of application names array */

  int   *app_info;     /* For each application, 4 integers: associated root in
                          base_comm, n_ranks, and indexes in app_names */
  char  *app_names;    /* Array of application type names and instance names */

  int      *app_status;   /* Synchronization status for each application */
  double   *app_timestep; /* Current time step for each application */

  MPI_Comm  base_comm;    /* Handle to base communicator */
  MPI_Comm  app_comm;     /* Handle to local application communicator */

};

/* Structure for MPI_DOUBLE_INT */

typedef struct {
  double  d;
  int     i;
} _mpi_double_int_t;

#endif /* defined(PLE_HAVE_MPI) */

/*============================================================================
 * Static global variables
 *============================================================================*/

static int _coupling_tag
= ('P'+'L'+'E'+'_'+'C'+'O'+'U'+'P'+'L'+'I'+'N'+'G') % 512;

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

/*----------------------------------------------------------------------------
 * Synchronize root applications in a set.
 *
 * This function should only be called by application roots.
 *
 * parameters:
 *   s         <-- pointer to PLE coupling MPI set info structure.
 *   sync_flag <-- synchronization info for current application.
 *   time_step <-- time step for current application.
 *   glob_vals <-> temporary synchronization values.
 *----------------------------------------------------------------------------*/

static void
_coupling_mpi_set_synchronize_roots(ple_coupling_mpi_set_t  *s,
                                    int                      sync_flag,
                                    double                   time_step,
                                    _mpi_double_int_t        glob_vals[])
{
  int i;
  MPI_Status status;
  int app_rank;

  int sync_root = -1;

  MPI_Comm_rank(s->app_comm, &app_rank);

  /* Return immediately if not application root */

  if (app_rank != 0 || (s->app_status[s->app_id] & PLE_COUPLING_NO_SYNC))
    return;

  /* First, sync data to root */

  for (i = 0; i < s->n_apps; i++) {
    if (! (s->app_status[i] & PLE_COUPLING_NO_SYNC)) {
      sync_root = i;
      break;
    }
  }

  if (sync_root == s->app_id) {

    for (i = 0; i < s->n_apps; i++) {
      if (s->app_status[i] & PLE_COUPLING_NO_SYNC) { /* Keep previous values */
        glob_vals[i].i = s->app_status[i];
        glob_vals[i].d = s->app_timestep[i];
      }
      else {
        if (i != sync_root)
          MPI_Recv(&(glob_vals[i]), 1, MPI_DOUBLE_INT, s->app_info[i*4],
                   _coupling_tag, s->base_comm, &status);
        else {
          glob_vals[i].i = sync_flag;
          glob_vals[i].d = time_step;
        }
      }
    }
  }
  else if (! (s->app_status[s->app_id] & PLE_COUPLING_NO_SYNC)) {
    _mpi_double_int_t send_vals;
    send_vals.i = sync_flag;
    send_vals.d = time_step;
    MPI_Send(&send_vals, 1, MPI_DOUBLE_INT, s->app_info[sync_root],
             _coupling_tag, s->base_comm);
  }

  /* Now, root sends data to all */

  if (sync_root == s->app_id) {
    for (i = 0; i < s->n_apps; i++) {
      if (i != sync_root && ! (s->app_status[i] & PLE_COUPLING_NO_SYNC))
        MPI_Send(glob_vals, s->n_apps, MPI_DOUBLE_INT, s->app_info[i*4],
                 _coupling_tag, s->base_comm);
    }
  }
  else
    MPI_Recv(glob_vals, s->n_apps, MPI_DOUBLE_INT, s->app_info[sync_root],
             _coupling_tag, s->base_comm, &status);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a group id within a communicator based on its name.
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
 * \param  [in]  comm       MPI communicator.
 * \param  [in]  group_name name associated with current group.
 *
 * \return id associated with local name.
 */
/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discover other applications in a set with a common communicator.
 *
 * In most cases, the base communicator is MPI_COMM_WORLD, and the local
 * application communicator app_comm is usually obtained from it using
 * MPI_Comm_split, but other combinations may be possible using MPI-2
 * process management functions.
 *
 * As this function requires communication between applications, it
 * is a collective function in base_comm.
 *
 * \param[in] sync_flag 1 if application is to be synchronized at each
 *                      time step, 0 if independent from others.
 * \param[in] app_type  name of current application type (software name).
 * \param[in] app_name  name of current application (data/case name).
 * \param[in] base_comm communicator associated with all applications.
 * \param[in] app_comm  communicator associated with local application.
 *
 * \return PLE coupling MPI set info structure.
 */
/*----------------------------------------------------------------------------*/

ple_coupling_mpi_set_t *
ple_coupling_mpi_set_create(int          sync_flag,
                            const char  *app_type,
                            const char  *app_name,
                            MPI_Comm     base_comm,
                            MPI_Comm     app_comm)
{
  int i, j;
  MPI_Status status;

  int set_rank = -1, app_rank = -1, n_app_ranks = 0;
  int root_marker = 0;
  int info_tag = 1, name_tag = 2;

  int counter[2] = {0, 0};
  int l_rank_info[5] = {-1, -1, -1, 1, 1}; /* Status (1) + rank info (4) */
  int *rank_info = NULL;
  char *app_names = NULL;

  ple_coupling_mpi_set_t *s = NULL;

  /* Initialization */

  PLE_MALLOC(s, 1, ple_coupling_mpi_set_t);

  s->n_apps = 0;
  s->app_id = -1;
  s->app_names_l = 0;
  s->app_info = NULL;
  s->app_names = NULL;

  s->base_comm = base_comm;
  s->app_comm = app_comm;

  /* Initialization */

  MPI_Comm_rank(base_comm, &set_rank);

  if (app_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(app_comm, &app_rank);
    MPI_Comm_size(app_comm, &n_app_ranks);
  }
  else {
    app_rank = 0;
    n_app_ranks = 1;
  }

  l_rank_info[0] = sync_flag | PLE_COUPLING_INIT;
  l_rank_info[1] = set_rank;
  l_rank_info[2] = n_app_ranks;
  if (app_type != NULL)
    l_rank_info[3] = strlen(app_type) + 1;
  if (app_name != NULL)
    l_rank_info[4] = strlen(app_name) + 1;

  if (app_rank == 0)
    root_marker = 1;

  /* Root rank of base_comm counts applications and receives info */

  MPI_Reduce(&root_marker, &(counter[0]), 1, MPI_INT, MPI_SUM, 0,
             base_comm);

  /* Root of base_comm collects all info */

  if (set_rank == 0) {

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
               base_comm, &status);

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
               name_tag, base_comm, &status);
      counter[1] += msg_len;
    }

  }

  /* Other root ranks send info to root */

  else if (app_rank == 0) { /* set_rank != 0 */

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

    MPI_Send(l_rank_info, 5, MPI_INT, 0, info_tag, base_comm);
    MPI_Send(sendbuf, sendbuf_l, MPI_CHAR, 0, name_tag, base_comm);

    PLE_FREE(sendbuf);
  }

  /* Now root broadcasts application info */

  MPI_Bcast(counter, 2, MPI_INT, 0, base_comm);

  if (set_rank != 0) {
    PLE_MALLOC(rank_info, counter[0]*5, int);
    PLE_MALLOC(app_names, counter[1], char);
  }

  MPI_Bcast(rank_info, counter[0]*5, MPI_INT, 0, base_comm);
  MPI_Bcast(app_names, counter[1], MPI_CHAR, 0, base_comm);

  /* Set global values */

  s->n_apps = counter[0];
  s->app_names_l = counter[1];
  s->app_names = app_names;

  PLE_MALLOC(s->app_info, s->n_apps*4, int);
  PLE_MALLOC(s->app_status, s->n_apps, int);
  PLE_MALLOC(s->app_timestep, s->n_apps, double);

  for (i = 0; i < s->n_apps && s->app_id < 0; i++) {
    for (j = 0; j < 4; j++)
      s->app_info[i*4 + j] = rank_info[i*5 + j+1];
    s->app_status[i] = rank_info[i*5];
    s->app_timestep[i] = 0.;
  }

  PLE_FREE(rank_info);

  /* Set rank set to that of the application root for matching */

  MPI_Bcast(&set_rank, 1, MPI_INT, 0, app_comm);

  for (i = 0; i < s->n_apps && s->app_id < 0; i++) {
    if (s->app_info[i*4] == set_rank)
      s->app_id = i;
  }

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free an PLE coupling MPI set info structure.
 *
 * \param[in, out] s pointer to structure that should be freed.
 */
/*----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_destroy(ple_coupling_mpi_set_t **s)
{
  ple_coupling_mpi_set_t *_s = *s;

  if (_s != NULL) {
    PLE_FREE(_s->app_info);
    PLE_FREE(_s->app_names);
    PLE_FREE(_s->app_status);
    PLE_FREE(_s->app_timestep);
    PLE_FREE(*s);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of applications in a coupled set.
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 *
 * \return number of application in set's common communicator.
 */
/*----------------------------------------------------------------------------*/

int
ple_coupling_mpi_set_n_apps(const ple_coupling_mpi_set_t  *s)
{
  int retval = 0;

  if (s != NULL)
    retval = s->n_apps;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of the local application in a coupled set.
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 *
 * \return id of the local application in set's common communicator.
 */
/*----------------------------------------------------------------------------*/

int
ple_coupling_mpi_set_get_app_id(const ple_coupling_mpi_set_t  *s)
{
  int retval = -1;

  if (s != NULL)
    retval = s->app_id;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return application information in set's common communicator.
 *
 * \param[in] s       pointer to PLE coupling MPI set info structure.
 * \param[in] app_id  application id
 *
 * \return application information structure.
 */
/*----------------------------------------------------------------------------*/

ple_coupling_mpi_set_info_t
ple_coupling_mpi_set_get_info(const ple_coupling_mpi_set_t  *s,
                              int                            app_id)
{
  ple_coupling_mpi_set_info_t  retval;

  retval.status = 0;
  retval.root_rank = -1;
  retval.n_ranks = 0;
  retval.app_type = NULL;
  retval.app_name = NULL;

  if (s != NULL) {
    if (app_id >= 0 && app_id < s->n_apps) {
      retval.status = s->app_status[app_id];
      retval.root_rank = s->app_info[app_id*4];
      retval.n_ranks = s->app_info[app_id*4 + 1];
      retval.app_type = s->app_names + s->app_info[app_id*4 + 2];
      retval.app_name = s->app_names + s->app_info[app_id*4 + 3];
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize applications in a set.
 *
 * Note that if a member of the set has used a PLE_COUPLING_STOP or
 * PLE_COUPLING_LAST flag when calling ple_coupling_mpi_set_create() or
 * or at the previous call to this function, it will not be synchronized
 * anymore (i.e. the PLE_COUPLING_NO_SYNC flag will be added).
 *
 * param[in] s         pointer to PLE coupling MPI set info structure.
 * param[in] sync_flag synchronization info for current application.
 * param[in] time_step time step for current application.
 */
/*----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_synchronize(ple_coupling_mpi_set_t  *s,
                                 int                      sync_flag,
                                 double                   time_step)
{
  int i;

  int last_sync_mask = (  PLE_COUPLING_NO_SYNC
                        | PLE_COUPLING_STOP
                        | PLE_COUPLING_LAST);

  _mpi_double_int_t *glob_vals = NULL;

  /* Update synchronization flag first. */

  if (s->app_status[s->app_id] & last_sync_mask)
    sync_flag = sync_flag | PLE_COUPLING_NO_SYNC;
  if (sync_flag & PLE_COUPLING_INIT)
    sync_flag -= PLE_COUPLING_INIT;

  for (i = 0; i < s->n_apps; i++) {
    if (s->app_status[i] & last_sync_mask)
      s->app_status[i] = (s->app_status[i] | PLE_COUPLING_NO_SYNC);
    if (s->app_status[i] & PLE_COUPLING_INIT)
      s->app_status[i] -= PLE_COUPLING_INIT;
  }

  /* Return immediately if not synchronized */

  if (s->app_status[s->app_id] & PLE_COUPLING_NO_SYNC)
    return;

  /* Synchronize roots, then broad cast locally */

  PLE_MALLOC(glob_vals, s->n_apps, _mpi_double_int_t);

  _coupling_mpi_set_synchronize_roots(s, sync_flag, time_step, glob_vals);

  MPI_Bcast(glob_vals, s->n_apps, MPI_DOUBLE_INT, 0, s->app_comm);

  /* Save values */

  for (i = 0; i < s->n_apps; i++) {
    s->app_status[i] = glob_vals[i].i;
    s->app_timestep[i] = glob_vals[i].d;
  }

  PLE_FREE(glob_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get status of applications in a set.
 *
 * This function allows access to the status flag of each synchronized
 * application in the set. It may be used immediately after
 * ple_coupling_mpi_set_create(), and flags are updated after each
 * call to ple_coupling_mpi_set_synchronize().
 *
 * param[in] s pointer to PLE coupling MPI set info structure.
 *
 * \return a pointer to the set's array of status flags
 */
/*----------------------------------------------------------------------------*/

const int *
ple_coupling_mpi_set_get_status(const ple_coupling_mpi_set_t  *s)
{
  const int *retval;

  if (s != NULL)
    return s->app_status;
  else retval = NULL;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get time steps in a set.
 *
 * This function may be called after ple_coupling_mpi_set_synchronize()
 * to query the time step values of each synchronized application in the set.
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 *
 * \return a pointer to the set's array of time steps
 */
/*----------------------------------------------------------------------------*/

const double *
ple_coupling_mpi_set_get_timestep(const ple_coupling_mpi_set_t  *s)
{
  const double *retval;

  if (s != NULL)
    retval = s->app_timestep;
  else
    retval = NULL;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute recommended time step for the current application based on
 * provided flags and values of applications in a set.
 *
 * The flags and values used to compute this recommended time step value
 * are update at each call to ple_coupling_mpi_set_synchronize().
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 *
 * \return computed application time step
 */
/*----------------------------------------------------------------------------*/

double
ple_coupling_mpi_set_compute_timestep(const ple_coupling_mpi_set_t  *s)
{
  if (s == NULL)
    return -1;

  double retval = s->app_timestep[s->app_id];
  int self_status = s->app_status[s->app_id];

  if (   (self_status & PLE_COUPLING_NO_SYNC)
      || (self_status & PLE_COUPLING_TS_INDEPENDENT))
    return retval;

  int leader_id = -1;
  double ts_min = -1.;

  /* Check for leader and minimum time step;
     unsynced and follower apps do not influence the result. */

  const int ignore_bits = PLE_COUPLING_NO_SYNC | PLE_COUPLING_TS_FOLLOWER;

  for (int i = 0; i < s->n_apps; i++) {

    if (s->app_status[i] & ignore_bits)
      continue;

    /* Handle leader */

    if (s->app_status[i] & PLE_COUPLING_TS_LEADER) {
      if (leader_id > -1) {
        ple_coupling_mpi_set_info_t ai_prev
          = ple_coupling_mpi_set_get_info(s, leader_id);
         ple_coupling_mpi_set_info_t ai
           = ple_coupling_mpi_set_get_info(s, i);
        ple_error
          (__FILE__, __LINE__, 0,
           _("\nApplication \"%s\" (%s) tried to set the group time step, but\n"
             "application \"%s\" (%s) has already done so."),
           ai.app_name, ai.app_type, ai_prev.app_name, ai_prev.app_type);
      }
      else {
        leader_id = i;
        retval = s->app_timestep[i];
      }
    }

    /* Update minimum time in all coupled and non-follower cases. */

    if (ts_min > 0) {
      if (s->app_timestep[i] < ts_min)
        ts_min = s->app_timestep[i];
    }
    else
      ts_min = s->app_timestep[i];
  }

  /* If aligning with minimal time step, update return value */

  if (self_status & PLE_COUPLING_TS_MIN) {
    if (ts_min > 0)
      retval = ts_min;
    else if (s->n_apps > 1 && !(self_status & PLE_COUPLING_STOP))
      ple_error
        (__FILE__, __LINE__, 0,
         _("\nCoupling parameters require following minimal time step,\n"
           "but all other synchronized applications are also followers."));
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an intracommunicator from a local and distant communicator
 * within a base communicator.
 *
 * \param[in] base_comm     communicator associated with both applications
 * \param[in] app_comm      communicator associated with local application
 * \param[in] distant_root  rank of distant group leader in base_comm
 * \param[in] new_comm      pointer to new communicator
 * \param[in] local_range   first and past-the last ranks of local application
 *                          in new communicator
 * \param[in] distant_range first and past-the last ranks of distant
 *                          application in new communicator
 */
/*----------------------------------------------------------------------------*/

void
ple_coupling_mpi_intracomm_create(MPI_Comm   base_comm,
                                  MPI_Comm   app_comm,
                                  int        distant_root,
                                  MPI_Comm  *new_comm,
                                  int        local_range[2],
                                  int        distant_range[2])
{
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

  MPI_Comm_rank(base_comm, &r_glob);

  MPI_Allreduce(&r_glob, &r_loc_max, 1, MPI_INT, MPI_MAX, app_comm);

  if (distant_root > r_loc_max)
    high = 0;

  MPI_Comm_size(app_comm, &n_loc_ranks);

  /* Create a reserved communicator */

  MPI_Intercomm_create(app_comm, 0, base_comm,
                       distant_root, _coupling_tag, &intercomm_tmp);

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get base communicator of an PLE coupling MPI set
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 */
/*----------------------------------------------------------------------------*/

MPI_Comm
ple_coupling_mpi_set_get_base_comm(const ple_coupling_mpi_set_t *s)
{
  assert(s != NULL);

  if (s == NULL) {
    ple_printf("  Coupling MPI set info: nil\n");
    return MPI_COMM_NULL;
  }

  return s->base_comm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump printout of an PLE coupling MPI set info structure.
 *
 * \param[in] s pointer to PLE coupling MPI set info structure.
 */
/*----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_dump(const ple_coupling_mpi_set_t  *s)
{
  int i;

  if (s == NULL) {
    ple_printf("  Coupling MPI set info: nil\n");
    return;
  }

  ple_printf("  Coupling MPI set info:        %p\n"
             "    number of applications:     %d\n"
             "    local application id:       %d\n"
             "    app_names_size:             %d\n\n",
             s, s->n_apps, s->app_id, s->app_names_l);

  for (i = 0; i < s->n_apps; i++)
    ple_printf("    Application id:      %d\n"
               "      root_rank:         %d\n"
               "      n_ranks:           %d\n"
               "      app_type:          \"%s\"\n"
               "      app_name:          \"%s\"\n"
               "      status:            %d\n"
               "      time step:         %f\n\n",
               i, s->app_info[i*4], s->app_info[i*4 + 1],
               s->app_names + s->app_info[i*4 + 2],
               s->app_names + s->app_info[i*4 + 3],
               s->app_status[i], s->app_timestep[i]);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

