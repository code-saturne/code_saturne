/*============================================================================
 * Numbering information for vectorization or multithreading
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_numbering.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_numbering.c
        Numbering information for vectorization or multithreading.
*/

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names for numbering types */

const char  *cs_numbering_type_name[] = {N_("default"),
                                         N_("vectorization"),
                                         N_("threads")};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return number of elements in a given thread exclusion group
 *
 * parameters:
 *   n     <-- pointer to numbering considered
 *   g_id  <-- associated group id
 *
 * returns: number of elements in given group
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_n_group_elts(const cs_numbering_t  *n,
              int                    g_id)
{
  cs_lnum_t n_elts = 0;
  for (int t_id = 0; t_id < n->n_threads; t_id++) {
    if (n->group_index[(t_id*n->n_groups + g_id)*2 + 1] > 0)
      n_elts +=   n->group_index[(t_id*n->n_groups + g_id)*2 + 1]
                - n->group_index[(t_id*n->n_groups + g_id)*2];
  }

  return n_elts;
}

/*----------------------------------------------------------------------------
 * Return number of elements in a numbering
 *
 * parameters:
 *   n     <-- pointer to numbering considered
 *
 * returns: number of elements
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_n_total_elts(const cs_numbering_t  *n)
{
  cs_lnum_t n_elts = 0;
  for (int g_id = 0; g_id < n->n_groups; g_id++)
    n_elts += _n_group_elts(n, g_id);

  return n_elts;
}

/*----------------------------------------------------------------------------
 * Log statistics for default numberings in serial mode.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *----------------------------------------------------------------------------*/

static void
_log_default_info_l(cs_log_t               log,
                    const cs_numbering_t  *numbering)
{
  if (numbering->n_no_adj_halo_elts > 0)
    cs_log_printf
      (log,
       _("  number of halo-independent elements: %7u\n"),
       (unsigned)(numbering->n_no_adj_halo_elts));

  cs_lnum_t n_elts = _n_total_elts(numbering);
  if (n_elts >= numbering->n_no_adj_halo_elts)
    cs_log_printf
      (log,
       _("  number of halo-adjacent elements:  %9u\n"),
       (unsigned)(n_elts - numbering->n_no_adj_halo_elts));
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log statistics for default numberings.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_log_default_info(cs_log_t               log,
                  const cs_numbering_t  *numbering,
                  MPI_Comm               comm)
{
  int n_domains = 1;

  if (comm != MPI_COMM_NULL)
    MPI_Comm_size(comm, &n_domains);

  if (n_domains == 1) {
    _log_default_info_l(log, numbering);
    return;
  }

  cs_gnum_t nnh_elts[2] = {numbering->n_no_adj_halo_elts, 0};
  cs_gnum_t nnh_tot[2], nnh_min[2], nnh_max[2];

  cs_lnum_t n_elts = _n_total_elts(numbering);
  if (n_elts >= numbering->n_no_adj_halo_elts)
    nnh_elts[1] = n_elts - numbering->n_no_adj_halo_elts;

  MPI_Allreduce(nnh_elts, nnh_tot, 2, CS_MPI_GNUM, MPI_SUM, comm);

  if (nnh_tot[0] > 0) {

    MPI_Allreduce(nnh_elts, nnh_min, 2, CS_MPI_GNUM, MPI_MIN, comm);
    MPI_Allreduce(nnh_elts, nnh_max, 2, CS_MPI_GNUM, MPI_MAX, comm);

    cs_log_printf
      (log,
       _("                                       minimum   maximum      mean\n"
         "  number of halo-independent elements: %7u %9u %9u\n"),
       (unsigned)nnh_min[0], (unsigned)nnh_max[0],
       (unsigned)(nnh_tot[0]/(cs_gnum_t)n_domains));

    cs_log_printf
      (log,
       _("  number of halo-adjacent elements:  %9u %9u %9u\n"),
       (unsigned)nnh_min[1], (unsigned)nnh_max[1],
       (unsigned)(nnh_tot[1]/(cs_gnum_t)n_domains));

  }
}

#endif /* have_MPI */

/*----------------------------------------------------------------------------
 * Log statistics for vectorization in serial mode.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *----------------------------------------------------------------------------*/

static void
_log_vector_info_l(cs_log_t               log,
                   const cs_numbering_t  *numbering)
{
  cs_log_printf(log,
                _("  vector size:                             %3d\n"),
                numbering->vector_size);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log statistics for vectorization.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_log_vector_info(cs_log_t               log,
                 const cs_numbering_t  *numbering,
                 MPI_Comm               comm)
{
  int n_domains = 1;

  if (comm != MPI_COMM_NULL)
    MPI_Comm_size(comm, &n_domains);

  if (n_domains == 1) {
    _log_vector_info_l(log, numbering);
    return;
  }

  int vs = numbering->vector_size;

  int vs_min = vs, vs_max = vs;

  /* plan for cases where vector sizes migh differ across nodes
     (heterogeneous architectures) */

  MPI_Allreduce(&vs, &vs_max, 1, MPI_INT, MPI_MAX, comm);
  if (numbering->type != CS_NUMBERING_VECTORIZE)
    vs_min = vs_max;
  MPI_Allreduce(&vs, &vs_min, 1, MPI_INT, MPI_MIN, comm);

  if (vs_min != vs_max) {
    int vs_tot;
    MPI_Allreduce(&vs, &vs_tot, 1, MPI_INT, MPI_SUM, comm);
    vs = vs_tot / n_domains;
  }

  if (vs_min != vs_max)
    cs_log_printf
      (log,
       _("                                       minimum   maximum      mean\n"
         "  vector size:                             %3d       %3d       %3d\n"),
       vs_min, vs_max, vs);
  else
    cs_log_printf(log,
                  _("  vector size:                             %3d\n"),
                  vs);
}

#endif /* have_MPI */

/*----------------------------------------------------------------------------
 * Estimate unbalance between threads of a given group.
 *
 * Test local operations related to renumbering.
 *
 * Unbalance is considered to be: (max/mean - 1)
 *
 * parameters:
 *   numbering <-- pointer to numbering structure
 *
 * returns:
 *   estimated unbalance for this group
 *----------------------------------------------------------------------------*/

static double
_estimate_imbalance(const cs_numbering_t  *numbering)
{
  double t_imbalance_tot = 0.0;

  if (numbering == NULL)
    return 0;

  if (numbering->type == CS_NUMBERING_THREADS) {

    int g_id;

    cs_lnum_t n_elts = 0;

    const int n_threads = numbering->n_threads;
    const int n_groups = numbering->n_groups;
    const cs_lnum_t *group_index = numbering->group_index;

    for (g_id = 0; g_id < n_groups; g_id++) {

      int t_id;
      double n_t_elts_mean, imbalance;

      cs_lnum_t n_t_elts_sum = 0;
      cs_lnum_t n_t_elts_max = 0;

      for (t_id = 0; t_id < n_threads; t_id++) {
        cs_lnum_t n_t_elts =   group_index[(t_id*n_groups + g_id)*2 + 1]
                             - group_index[(t_id*n_groups + g_id)*2];
        n_t_elts = CS_MAX(n_t_elts, 0);
        n_t_elts_sum += n_t_elts;
        n_t_elts_max = CS_MAX(n_t_elts, n_t_elts_max);
      }

      n_elts += n_t_elts_sum;

      n_t_elts_mean = (double)n_t_elts_sum / n_threads;

      imbalance = (n_t_elts_max / n_t_elts_mean) - 1.0;
      t_imbalance_tot += imbalance*n_t_elts_sum;

    }

    t_imbalance_tot /= n_elts;

  }

  return t_imbalance_tot;
}

/*----------------------------------------------------------------------------
 * Log statistics for threads and groups in serial mode.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *----------------------------------------------------------------------------*/

static void
_log_threading_info_l(cs_log_t               log,
                      const cs_numbering_t  *numbering)
{
  double imbalance = _estimate_imbalance(numbering);

  const int n_groups = numbering->n_groups;
  const int n_threads = numbering->n_threads;

  cs_log_printf
    (log,
     _("  number of threads:                       %3d\n"
       "  number of exclusive groups:              %3d\n"),
     n_threads, n_groups);

  for (int g_id = 0; g_id < n_groups; g_id++) {
    cs_lnum_t n_elts = _n_group_elts(numbering, g_id);
    cs_log_printf
      (log,
       _("   number of elements in group %2d:   %9u\n"),
       g_id, (unsigned)n_elts);
  }

  cs_log_printf
    (log,
     _("  estimated thread imbalance:            %5.3f\n"),
     imbalance);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Log statistics for threads and groups.
 *
 * parameters:
 *   log           <-- log type
 *   numbering     <-- pointer to numbering considered
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_log_threading_info(cs_log_t               log,
                    const cs_numbering_t  *numbering,
                    MPI_Comm               comm)
{
  int n_domains = 1;

  if (comm != MPI_COMM_NULL)
    MPI_Comm_size(comm, &n_domains);

  if (n_domains == 1) {
    _log_threading_info_l(log, numbering);
    return;
  }

  double imbalance = _estimate_imbalance(numbering);

  const int n_groups = numbering->n_groups;
  const int n_threads = numbering->n_threads;

  int count_l[3] = {n_threads,
                    n_groups,
                    numbering->n_no_adj_halo_groups};
  int count_tot[3], count_min[3], count_max[3];
  double imb_tot, imb_min, imb_max;

  MPI_Allreduce(count_l, count_tot, 3, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(count_l, count_min, 3, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(count_l, count_max, 3, MPI_INT, MPI_MAX, comm);

  cs_log_printf
    (log,
     _("                                       minimum   maximum      mean\n"
       "  number of threads:                       %3d       %3d       %3d\n"
       "  number of exclusive groups:              %3d       %3d       %3d\n"),
     count_min[0], count_max[0], (int)(count_tot[0]/n_domains),
     count_min[1], count_max[1], (int)(count_tot[1]/n_domains));

  if (count_tot[2] > 0)
    cs_log_printf
      (log,
       _("  number of halo-independent groups:       %3d       %3d       %3d\n"),
       count_min[2], count_max[2], (int)(count_tot[2]/n_domains));

#if defined(HAVE_MPI_IN_PLACE)

  cs_gnum_t *group_sum;
  cs_lnum_t *group_min;
  cs_lnum_t *group_max;

  BFT_MALLOC(group_sum, count_max[1]*2, cs_gnum_t);
  BFT_MALLOC(group_min, count_max[1], cs_lnum_t);
  BFT_MALLOC(group_max, count_max[1], cs_lnum_t);

  for (int g_id = 0; g_id < n_groups; g_id++) {
    cs_lnum_t n_elts = _n_group_elts(numbering, g_id);
    group_sum[g_id*2]   = n_elts;
    group_sum[g_id*2+1] = 1;
    group_min[g_id]     = n_elts;
    group_max[g_id]     = n_elts;
  }
  for (int g_id = n_groups; g_id < count_max[1]; g_id++) {
    group_sum[g_id*2] = 0;
    group_sum[g_id*2+1] = 0;
    group_max[g_id] = 0;
  }

  MPI_Allreduce(MPI_IN_PLACE, group_sum, count_max[1]*2, CS_MPI_GNUM,
                MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, group_max, count_max[1], CS_MPI_LNUM,
                MPI_MAX, comm);

  for (int i = n_groups; i < count_max[1]; i++)
    group_min[i] = group_max[i]; /* only groups locally present contribute */
  MPI_Allreduce(MPI_IN_PLACE, group_min, count_max[1], CS_MPI_LNUM,
                MPI_MIN, comm);

  for (int i = 0; i < count_max[1]; i++) {
    cs_lnum_t group_mean = 0;
    if (group_sum[i*2+1] > 0)
      group_mean = group_sum[i*2] / group_sum[i*2+1];
    cs_log_printf
      (log,
       _("   number of elements in group %2d:   %9u %9u %9u\n"),
       i,
       (unsigned)group_min[i], (unsigned)group_max[i], (unsigned)group_mean);
  }

  BFT_FREE(group_sum);
  BFT_FREE(group_min);
  BFT_FREE(group_max);

#endif

  MPI_Allreduce(&imbalance, &imb_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(&imbalance, &imb_min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&imbalance, &imb_max, 1, MPI_DOUBLE, MPI_MAX, comm);

  cs_log_printf
    (log,
     _("  estimated thread imbalance:            %5.3f     %5.3f     %5.3f\n"),
     imb_min, imb_max, (double)(imb_tot/n_domains));
}

#endif /* have_MPI */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure.
 *
 * \param[in]  n_elts  number of associated elements
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_default(cs_lnum_t  n_elts)
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_DEFAULT;

  numbering->vector_size = 1;

  numbering->n_threads = 1;
  numbering->n_groups = 1;

  numbering->n_no_adj_halo_groups = 0;
  numbering->n_no_adj_halo_elts = 0;

  BFT_MALLOC(numbering->group_index, 2, cs_lnum_t);
  numbering->group_index[0] = 0;
  numbering->group_index[1] = n_elts;

  return numbering;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure
 *        in case of vectorization.
 *
 * \param[in]  n_elts       number of associated elements
 * \param[in]  vector_size  vector size used for this vectorization
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_vectorized(cs_lnum_t  n_elts,
                               int        vector_size)
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_VECTORIZE;

  numbering->vector_size = vector_size;

  numbering->n_threads = 1;
  numbering->n_groups = 1;

  numbering->n_no_adj_halo_groups = 0;
  numbering->n_no_adj_halo_elts = 0;

  BFT_MALLOC(numbering->group_index, 2, cs_lnum_t);
  numbering->group_index[0] = 0;
  numbering->group_index[1] = n_elts;

  return numbering;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure
 *        in case of threading.
 *
 * \param[in]  n_elts       number of associated elements
 * \param[in]  n_groups     number of groups
 * \param[in]  group_index  group_index[thread_id*group_id*2 + group_id*2] and
 *                          group_index[thread_id*group_id*2 + group_id*2 +1]
 *                          define the start and end ids for entities in a
 *                          given group and thread;
 *                          (size: n_groups *2 * n_threads)
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_threaded(int        n_threads,
                             int        n_groups,
                             cs_lnum_t  group_index[])
{
  cs_numbering_t  *numbering = NULL;

  BFT_MALLOC(numbering, 1, cs_numbering_t);

  numbering->type = CS_NUMBERING_THREADS;

  numbering->vector_size = 1;

  numbering->n_threads = n_threads;
  numbering->n_groups = n_groups;

  numbering->n_no_adj_halo_groups = 0;
  numbering->n_no_adj_halo_elts = 0;

  BFT_MALLOC(numbering->group_index, n_threads*2*n_groups, cs_lnum_t);

  memcpy(numbering->group_index,
         group_index,
         (n_threads*2*n_groups) * sizeof(cs_lnum_t));

  return numbering;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a numbering information structure.
 *
 * \param[in, out]  numbering  pointer to cs_numbering_t structure pointer
 *                             (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_destroy(cs_numbering_t  **numbering)
{
  if (*numbering != NULL) {

    cs_numbering_t  *_n = *numbering;

    BFT_FREE(_n->group_index);

    BFT_FREE(*numbering);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log information relative to a cs_numbering_t structure.
 *
 * \param[in]  log          log type
 * \param[in]  description  description of numbering type
 * \param[in]  numbering    pointer to cs_numbering_t structure (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_log_info(cs_log_t               log,
                      const char            *description,
                      const cs_numbering_t  *numbering)
{
  if (numbering == NULL)
    return;

  cs_log_printf(log, _("\n Numbering for %s:\n"), description);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int type_counts[3];
    int type_counts_l[3] = {0, 0, 0};
    type_counts_l[numbering->type] = 1;

    MPI_Allreduce(type_counts_l, type_counts, 3, MPI_INT,
                  MPI_SUM, cs_glob_mpi_comm);

    MPI_Comm comm = cs_glob_mpi_comm;
    for (cs_numbering_type_t i = 0; i <= CS_NUMBERING_THREADS; i++) {

      if (type_counts[i] > 0) {

        if (type_counts[i] < cs_glob_n_ranks) {
          if (comm == cs_glob_mpi_comm)
            MPI_Comm_split(cs_glob_mpi_comm, numbering->type, cs_glob_rank_id,
                           &comm);
          cs_log_printf(log, _("\n Numbering type: %s for %d rank(s)\n"),
                        _(cs_numbering_type_name[i]), type_counts[i]);
        }
        else
          cs_log_printf(log, _("\n type: %s\n"),
                        _(cs_numbering_type_name[i]));

        switch(i) {
        case CS_NUMBERING_DEFAULT:
          _log_default_info(log, numbering, comm);
          break;
        case CS_NUMBERING_VECTORIZE:
          _log_vector_info(log, numbering, comm);
          break;
        case CS_NUMBERING_THREADS:
          _log_threading_info(log, numbering, comm);
          break;
        }

      }

    }

    if (comm != cs_glob_mpi_comm)
      MPI_Comm_free(&comm);

  }

#endif /* if !defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {

    cs_log_printf(log, _("\n type: %s\n"),
                  _(cs_numbering_type_name[numbering->type]));

    switch(numbering->type) {
    case CS_NUMBERING_DEFAULT:
      _log_default_info_l(log, numbering);
      break;
    case CS_NUMBERING_VECTORIZE:
      _log_vector_info_l(log, numbering);
      break;
    case CS_NUMBERING_THREADS:
      _log_threading_info_l(log, numbering);
      break;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_numbering_t structure.
 *
 * \param[in]  numbering    pointer to cs_numbering_t structure (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_dump(const cs_numbering_t  *numbering)
{
  int  i, j;

  if (numbering == NULL) {
    bft_printf("\n  Numbering: nil (default)\n");
    return;
  }

  bft_printf("\n  Numbering:           %p\n"
             "  type:                  %s\n"
             "  vector_size:           %d\n"
             "  n_threads:             %d\n"
             "  n_groups:              %d\n"
             "  n_no_adj_halo_groups:  %d\n"
             "  n_no_adj_halo_elts:    %ld\n",
             (const void *)numbering, cs_numbering_type_name[numbering->type],
             numbering->vector_size,
             numbering->n_threads, numbering->n_groups,
             numbering->n_no_adj_halo_groups,
             (long)(numbering->n_no_adj_halo_elts));

  if (numbering->group_index != NULL) {

    bft_printf("\n  group start index:\n"
               "\n    group_id thread_id (id) start_index\n");

    for (i = 0; i < numbering->n_groups; i++) {
      for (j = 0; j < numbering->n_threads; j++) {
        int k = (j*numbering->n_groups + i);
        bft_printf("      %2d       %2d      %3d   %d\n",
                   i, j, k, (int)(numbering->group_index[k*2]));
      }
      int l = (numbering->n_threads-1)*numbering->n_groups + i;
      bft_printf("      %2d                     %d\n",
                 i, (int)(numbering->group_index[l*2+1]));
    }
  }

  bft_printf("\n\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
