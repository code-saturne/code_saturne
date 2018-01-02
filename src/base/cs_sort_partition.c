/*============================================================================
 * Data partitioning for parallel sort or order operations.
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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sort_partition.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  _distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int _distrib_n_iter_max = 5;

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Evaluate a distribution array.
 *
 * parameters:
 *   n_ranges     <-- Number of ranges in the distribution
 *   distribution <-- Number of elements associated to each range of
 *                    the distribution
 *   optim        <-- Optimal count in each range
 *
 * returns:
 *   a fit associated to the distribution. If fit = 0,
 *   distribution is perfect.
 *----------------------------------------------------------------------------*/

static double
_evaluate_distribution(cs_lnum_t    n_ranges,
                       cs_gnum_t   *distribution,
                       double       optim)
{
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (cs_lnum_t i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = CS_MAX(d_up, distribution[i] - optim);
    else
      d_low = CS_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    bft_printf("<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}

/*----------------------------------------------------------------------------
 * Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 * parameters:
 *   n_ranks         <-- number of ranks (= number of ranges)
 *   sampling_factor <-- number of samples per rank
 *   elt_size        <-- size associated with each element
 *   gsum_weight     <-- global sum of all weightings
 *   n_elts          <-- local number of elements
 *   elts            <-- local list of elements to distribute
 *   weight          <-- optional weight related to each element, or NULL
 *   order           <-- ordering array
 *   sampling        <-- sampling array
 *   s_to_elt        <-- pointer to coordinate to element conversion function
 *   compare         <-- pointer to comparison function
 *   f_input         <-- optional input to s_to_elt and compare, or NULL
 *   c_freq          <-> pointer to the cumulative frequency array
 *   g_distrib       <-> pointer to a distribution array
 *   comm            <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_define_rank_distrib(int                           n_ranks,
                     cs_lnum_t                     sampling_factor,
                     cs_gnum_t                     gsum_weight,
                     size_t                        elt_size,
                     cs_lnum_t                     n_elts,
                     const void                   *elts,
                     const cs_lnum_t              *weight,
                     const cs_lnum_t               order[],
                     const double                  sampling[],
                     cs_sort_partition_s_to_elt_t  s_to_elt,
                     cs_sort_partition_compare_t   compare,
                     const void                   *f_input,
                     double                        cfreq[],
                     cs_gnum_t                     g_distrib[],
                     MPI_Comm                      comm)
{
  unsigned char  _sample_buffer[1024];
  unsigned char  *sample_code = _sample_buffer;
  const unsigned char  *restrict _elts = elts;
  int  bucket_id = 1;
  cs_gnum_t   *l_distrib = NULL;

  const cs_lnum_t  n_samples = sampling_factor * n_ranks;

  if (elt_size > 1024)
    BFT_MALLOC(sample_code, elt_size, unsigned char);

  /* Initialization */

  BFT_MALLOC(l_distrib, n_samples, cs_gnum_t);

  for (cs_lnum_t id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  /* Elements are assumed to be pre- ordered */

  s_to_elt(sampling[bucket_id], sample_code, f_input);

  if (weight != NULL) {

    for (size_t i = 0; i < (size_t)n_elts; i++) {

      size_t o_id = order[i];

      if (compare(sample_code, _elts + o_id*elt_size, f_input) >= 0)
        l_distrib[bucket_id - 1] += weight[o_id];

      else {

        while (compare(_elts + o_id*elt_size, sample_code, f_input) > 0) {

          bucket_id++;
          assert(bucket_id < n_samples + 1);

          s_to_elt(sampling[bucket_id], sample_code, f_input);
        }

        l_distrib[bucket_id - 1] += weight[o_id];

      }

    } /* End of loop on elements */

  }
  else { /* weight == NULL */

    for (size_t i = 0; i < (size_t)n_elts; i++) {

      size_t o_id = order[i];

      if (compare(sample_code, _elts + o_id*elt_size, f_input) >= 0)
        l_distrib[bucket_id - 1] += 1;

      else {

        while (compare(_elts + o_id*elt_size, sample_code, f_input) > 0) {

          bucket_id++;
          assert(bucket_id < n_samples + 1);

          s_to_elt(sampling[bucket_id], sample_code, f_input);
        }

        l_distrib[bucket_id - 1] += 1;

      }

    } /* End of loop on elements */

  }

  /* Define the global distribution */

  MPI_Allreduce(l_distrib, g_distrib, n_samples, CS_MPI_GNUM, MPI_SUM,
                comm);

  BFT_FREE(l_distrib);

  /* Define the cumulative frequency related to g_distribution */

  cfreq[0] = 0.;
  for (cs_lnum_t id = 0; id < n_samples; id++)
    cfreq[id+1] = cfreq[id] + (double)g_distrib[id]/(double)gsum_weight;
  cfreq[n_samples] = 1.0;

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  if (cs_glob_rank_id <= 0) {

    FILE  *dbg_file = NULL;
    char  *rfilename = NULL;
    int  len;
    static int  loop_id1 = 0;

    len = strlen("DistribOutput_l.dat")+1+2;
    BFT_MALLOC(rfilename, len, char);
    sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1);

    loop_id1++;

    dbg_file = fopen(rfilename, "w");

    fprintf(dbg_file,
            "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |"
            "Global Distrib\n");
    for (i = 0; i < n_samples; i++)
      fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
              i, (double)i/(double)n_samples, cfreq[i],
              sampling[i], g_distrib[i]);
    fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
            i, 1.0, 1.0, 1.0, 0);

    fclose(dbg_file);
    BFT_FREE(rfilename);

  }

#endif /* debugging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {

    cs_gnum_t   sum = 0;
    cs_lnum_t   shift = rank_id * sampling_factor;

    for (cs_lnum_t id = 0; id < sampling_factor; id++)
      sum += g_distrib[shift + id];
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check in debug */
  {
    cs_gnum_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight)
      bft_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
  }
#endif /* sanity check */

  if (sample_code != _sample_buffer)
    BFT_FREE(sample_code);
}

/*----------------------------------------------------------------------------
 * Update a distribution associated to sampling to assume a well-balanced
 * distribution of the leaves of the tree.
 *
 * parameters:
 *   n_ranks         <-- number of ranks (= number of ranges)
 *   sampling_factor <-- number of samples per rank
 *   sampling        <-> pointer to pointer to a sampling array
 *   comm            <-- mpi communicator
 *----------------------------------------------------------------------------*/

static void
_update_sampling(int        n_ranks,
                 cs_lnum_t  sampling_factor,
                 double     c_freq[],
                 double    *sampling[])
{
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  double  *new_sampling = NULL, *_sampling = *sampling;

  const cs_lnum_t  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute new_sampling */

  BFT_MALLOC(new_sampling, n_samples + 1, double);

  new_sampling[0] = _sampling[0];

  cs_lnum_t next_id = 1;

  for (cs_lnum_t i = 0; i < n_samples; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (cs_lnum_t j = next_id; j < n_samples + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = s_low + delta;
    }
    else /* f_high = f_low */
      new_sampling[i+1] = s_low + 0.5 * (s_low + s_high);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf(" <_update_distrib> (rank: %d) delta: %g, target: %g,"
               " next_id: %d, f_low: %g, f_high: %g, s_low: %g, s_high: %g\n"
               "\t => new_sampling: %g\n",
               cs_glob_rank_id, delta, target_freq, next_id,
               f_low, f_high, s_low, s_high, new_sampling[i+1]);
#endif

  } /* End of loop on samples */

  new_sampling[n_samples] = 1.0;

  BFT_FREE(_sampling);

  /* Return pointers */

  *sampling = new_sampling;
}

/*----------------------------------------------------------------------------
 * Compute a sampling array which assumes a well-balanced distribution of
 * leaves of the tree among the ranks.
 *
 * parameters:
 *   n_ranks         <--  number of ranks
 *   sampling_factor <-- number of samples per rank
 *   elt_size        <-- size associated with each element
 *   n_elts          <-- local number of elements
 *   elts            <-- local list of elements to distribute
 *   weight          <-- optional weight for each element, or NULL
 *   order           <-- ordering array
 *   s_to_elt        <-- pointer to coordinate to element conversion function
 *   compare         <-- pointer to comparison function
 *   f_input         <-- optional input to s_to_elt and compare, or NULL
 *   sampling        <-> pointer to pointer to a sampling array
 *   comm            <-- mpi communicator
 *
 * returns:
 *   fit associated to the returned sampling array
 *----------------------------------------------------------------------------*/

static double
_bucket_sampling(int                            n_ranks,
                 cs_lnum_t                      sampling_factor,
                 size_t                         elt_size,
                 cs_lnum_t                      n_elts,
                 const void                    *elts,
                 const cs_lnum_t               *weight,
                 const cs_lnum_t                order[],
                 cs_sort_partition_s_to_elt_t   s_to_elt,
                 cs_sort_partition_compare_t    compare,
                 const void                    *f_input,
                 double                        *sampling[],
                 MPI_Comm                       comm)
{
  double  fit, best_fit;

  cs_gnum_t   *distrib = NULL;
  double  *cfreq = NULL, *best_sampling = NULL;
  double  *_sampling = *sampling;

  const cs_lnum_t  n_samples = sampling_factor * n_ranks;
  const double  unit = 1/(double)n_samples;

  /* Compute the global number of elements and the optimal number of elements
     on each rank */

  cs_gnum_t  lsum_weight = 0, gsum_weight = 0;

  if (weight != NULL) {
    for (cs_lnum_t j = 0; j < n_elts; j++)
      lsum_weight += weight[j];
  }
  else
    lsum_weight = n_elts;

  MPI_Allreduce(&lsum_weight, &gsum_weight, 1, CS_MPI_GNUM, MPI_SUM, comm);

  double optim = (double)gsum_weight / (double)n_ranks;

  /* Define a naive sampling (uniform distribution) */

  for (cs_lnum_t i = 0; i < n_samples + 1; i++)
    _sampling[i] = i*unit;

  /* Define the distribution associated to the current sampling array */

  BFT_MALLOC(distrib, n_samples, cs_gnum_t);
  BFT_MALLOC(cfreq, n_samples + 1, double);

  _define_rank_distrib(n_ranks,
                       sampling_factor,
                       gsum_weight,
                       elt_size,
                       n_elts,
                       elts,
                       weight,
                       order,
                       _sampling,
                       s_to_elt,
                       compare,
                       f_input,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  fit = _evaluate_distribution(n_ranks, distrib, optim);
  best_fit = fit;

  BFT_MALLOC(best_sampling, n_samples + 1, double);

  for (cs_lnum_t i = 0; i < n_samples + 1; i++)
    best_sampling[i] = _sampling[i];

  /* Loop to get a better sampling array */

  for (int n_iters = 0;
       n_iters < _distrib_n_iter_max && fit > _distrib_tol;
       n_iters++) {

    _update_sampling(n_ranks, sampling_factor, cfreq, &_sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(n_ranks,
                         sampling_factor,
                         gsum_weight,
                         elt_size,
                         n_elts,
                         elts,
                         weight,
                         order,
                         _sampling,
                         s_to_elt,
                         compare,
                         f_input,
                         cfreq,
                         distrib,
                         comm);

    fit = _evaluate_distribution(n_ranks, distrib, optim);

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (cs_lnum_t i = 0; i < n_samples + 1; i++)
        best_sampling[i] = _sampling[i];

    }

  } /* End of while */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    bft_printf("\n  <_bucket_sampling> n_iter: %d, opt: %g, best_fit: %g\n",
               n_iters, optim, best_fit);
#endif

  /* Free memory */

  BFT_FREE(cfreq);
  BFT_FREE(distrib);
  BFT_FREE(_sampling);

  *sampling = best_sampling;

  return best_fit;
}

/*----------------------------------------------------------------------------
 * Build a real values rank index.
 *
 * The rank_index[i] contains the first element assigned to rank [i].
 *
 * parameters:
 *   sampling_factor <-- number of samples per rank
 *   elt_size        <-- size associated with each element
 *   n_elts          <-- number of elements to be indexed
 *   elt             <-- array of elements
 *   weight          <-- optional weight related to each element, or NULL
 *   order           <-- ordering array
 *   s_to_elt        <-- pointer to coordinate to element conversion function
 *   compare         <-- pointer to comparison function
 *   f_input         <-- optional input to s_to_elt and compare, or NULL
 *   rank_index      <-> pointer to the global values rank index
 *   comm            <-- MPI communicator on which we build the global index
 *
 * returns:
 *   the fit related to the element distribution (lower is better).
 *----------------------------------------------------------------------------*/

static double
_build_rank_index(cs_lnum_t                      sampling_factor,
                  size_t                         elt_size,
                  cs_lnum_t                      n_elts,
                  const void                    *elts,
                  const cs_lnum_t               *weight,
                  const cs_lnum_t                order[],
                  void                          *rank_index,
                  cs_sort_partition_s_to_elt_t   s_to_elt,
                  cs_sort_partition_compare_t    compare,
                  const void                    *f_input,
                  MPI_Comm                       comm)
{
  int     n_ranks;

  double  *sampling = NULL;

  /* Allocations and Initialization */

  MPI_Comm_size(comm, &n_ranks);

  cs_lnum_t n_samples = sampling_factor * n_ranks;

  BFT_MALLOC(sampling, n_samples + 1, double);

  for (cs_lnum_t i = 0; i < n_samples + 1; i++)
    sampling[i] = 0.0;

  double best_fit = _bucket_sampling(n_ranks,
                                     sampling_factor,
                                     elt_size,
                                     n_elts,
                                     elts,
                                     weight,
                                     order,
                                     s_to_elt,
                                     compare,
                                     f_input,
                                     &sampling,
                                     comm);

  /* Define rank index */

  unsigned char *_rank_index = rank_index;

  for (size_t rank_id = 0; rank_id < (size_t)n_ranks + 1; rank_id++) {
    cs_lnum_t id = rank_id * sampling_factor;
    s_to_elt(sampling[id], _rank_index + rank_id*elt_size, f_input);
  }

  /* Free memory */

  BFT_FREE(sampling);

  return best_fit;
}

/*----------------------------------------------------------------------------
 * Get the quantile associated to a Morton code using a binary search.
 *
 * No check is done to ensure that the code is present in the quantiles.
 *
 * parameters:
 *   n_quantiles  <-- number of quantiles
 *   elt_size     <-- size associated with each element
 *   elt          <-- pointer to data element
 *   rank_index   <-- pointer to the global values rank index
 *   compare      <-- pointer to comparison function
 *   f_input      <-- optional input to s_to_elt and compare, or NULL
 *
 * returns:
 *   id associated to the given code in the codes array.
 *----------------------------------------------------------------------------*/

static size_t
_quantile_search(size_t                        n_quantiles,
                 size_t                        elt_size,
                 const void                   *elt,
                 const void                   *rank_index,
                 cs_sort_partition_compare_t   compare,
                 const void                   *f_input)
{
  size_t mid_id = 0;
  size_t start_id = 0;
  size_t end_id = n_quantiles;

  const unsigned char  *_rank_index = rank_index;

  /* use binary search */

  while (start_id + 1 < end_id) {
    mid_id = start_id + ((end_id -start_id) / 2);
    if (compare(_rank_index + mid_id*elt_size, elt, f_input) > 0)
      end_id = mid_id;
    else
      start_id = mid_id;
  }

  /* We may have stopped short of the required value,
     or have multiple occurences of a quantile start
     (in case of empty quantiles), of which we want to
     find the find highest one */

  while (   start_id < n_quantiles - 1
         && (compare(elt, _rank_index + (start_id+1)*elt_size, f_input) >= 0))
    start_id++;

  return start_id;
}

#endif /* HAVE_MPI */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine to which rank data elements should be sent for parallel
 *        sorting or ordering.
 *
 * \param[in]   sampling_factor  number of samples per rank
 * \param[in]   elt_size         size associated with each element
 * \param[in]   n_elts           number of elements to be indexed
 * \param[in]   elts             array of elements
 * \param[in]   weight           optional weight of each element, or NULL
 * \param[in]   order            ordering array
 * \param[out]  dest_rank_id     destination rank id (size: n_elts)
 * \param[in]   s_to_elt         coordinate to element conversion function
 * \param[in]   compare          comparison function
 * \param[in]   f_input          optional input to s_to_elt and compare, or NULL
 * \param[in]   comm             MPI communicator on which we build the
 *                               global index
 */
/*----------------------------------------------------------------------------*/

void
cs_sort_partition_dest_rank_id(cs_lnum_t                      sampling_factor,
                               size_t                         elt_size,
                               cs_lnum_t                      n_elts,
                               const void                    *elts,
                               const cs_lnum_t               *weight,
                               const cs_lnum_t                order[],
                               int                            dest_rank_id[],
                               cs_sort_partition_s_to_elt_t   s_to_elt,
                               cs_sort_partition_compare_t    compare,
                               const void                    *f_input,
                               MPI_Comm                       comm)
{
  int n_ranks;
  MPI_Comm_size(comm, &n_ranks);

  unsigned char *_rank_index;
  const unsigned char *_elts = elts;

  BFT_MALLOC(_rank_index, (size_t)(n_ranks+1)*elt_size, unsigned char);

  _build_rank_index(sampling_factor,
                    elt_size,
                    n_elts,
                    elts,
                    weight,
                    order,
                    _rank_index,
                    s_to_elt,
                    compare,
                    f_input,
                    comm);

# pragma omp parallel for  if(n_elts > CS_THR_MIN)
  for (size_t i = 0; i < (size_t)n_elts; i++) {
    dest_rank_id[i] = _quantile_search(n_ranks,
                                       elt_size,
                                       _elts + i*elt_size,
                                       _rank_index,
                                       compare,
                                       f_input);
  }

  BFT_FREE(_rank_index);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/

END_C_DECLS
