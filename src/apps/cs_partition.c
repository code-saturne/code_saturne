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

/*============================================================================
 * Main program
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * METIS library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HAVE_PARMETIS)

#include <parmetis.h>

#endif

#if defined(HAVE_METIS_H)

#include <metis.h>

#else

/* Caution:
   --------
   File included by <metis.h> for METIS 4.0 lead to compilation
   errors on Linux in C99 mode due to redeclaration of functions from
   <stdlib.h>, so we include our own prototypes here.
   Also, ParMETIS's parmetis.h does not include prototypes for
   direct calling of METIS functions, so we also include prototypes
   in this case. */

#if !defined(IDXTYPE_INT)
typedef int idxtype;
#endif

void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                              int *, int *, int *, int *, int *, idxtype *);
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                         int *, int *, int *, int *, int *, idxtype *);

#endif /* defined(HAVE_METIS_H) */

#ifdef __cplusplus
}
#endif

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */


/*----------------------------------------------------------------------------
 * SCOTCH library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

#include <scotch.h>

#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_block_to_part.h>
#include <fvm_file.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_io.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert an argument to an integer and check its validity
 *
 * parameters:
 *   arg_id  <-- index of argument in argv
 *   argc    <-- number of command line arguments
 *   argv    <-- array of command line arguments
 *   argerr  --> error indicator
 *
 * returns:
 *   integer value
 *----------------------------------------------------------------------------*/

static int
_arg_to_int(int    arg_id,
            int    argc,
            char  *argv[],
            int   *argerr)
{
  char  *start = NULL;
  char  *end = NULL;
  int  retval = 0;

  *argerr = 0;

  if (arg_id < argc) {
    start = argv[arg_id];
    end = start + strlen(start);
    retval = strtol(start, &end, 0);
    if (end != start + strlen(start)) *argerr = 1;
  }
  else {
    *argerr = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Print command line help
 *
 * parameters:
 *   name  <-- name of executable program
 *----------------------------------------------------------------------------*/

static void
_arg_env_help(const char  *name)
{
  if (cs_glob_rank_id >= 1)
    return;

  bft_printf
    (_("Usage: %s [options] <n_ranks> [<n_ranks_2> ... <n_ranks_n>]\n\n"
       "Compute domain partition(s) of Code_Saturne Preprocessor output.\n\n"
       "<n_ranks> number of ranks of destination partition (> 1).\n"),
     name);

  bft_printf (_("\nCommand line options:\n\n"));

  bft_printf
    (_(" --no-perio        ignore periodicity information if present.\n"
       " --no-write        do not write output.\n"));

#if defined(HAVE_PARMETIS)
  bft_printf
    (_(" --metis           use ParMETIS for partitioning (default).\n"));
#elif defined(HAVE_METIS)
  bft_printf
    (_(" --metis           use METIS for partitioning (default).\n"));
#endif
#if defined(HAVE_PTSCOTCH)
  bft_printf
    (_(" --scotch          use PT-SCOTCH for partitioning.\n"));
#elif defined(HAVE_SCOTCH)
  bft_printf
    (_(" --scotch          use SCOTCH for partitioning.\n"));
#endif
  bft_printf
    (_(" --mpi             use MPI for parallelism or coupling\n"));
  bft_printf
    (_(" --mpi-io          <mode> set parallel I/O behavior\n"
       "                     off: do not use MPI-IO\n"
       "                     eo:  MPI-IO with explicit offsets\n"
       "                          (default if available)\n"
       "                     ip:  MPI-IO with individual file pointers\n"));
  bft_printf
    (_(" --log             output redirection for rank -1 or 0:\n"
       "                     0: standard output\n"
       "                     1: output in \"partition.log\" (default)\n"));
  bft_printf
    (_(" --logp            output redirection for rank > 0:\n"
       "                    -1: remove output (default)\n"
       "                     0: no redirection (if independant\n"
       "                        terminals, debugger type)\n"
       "                     1: output in \"partition_n<rank>.log\"\n"));
  bft_printf
    (_(" -v, --version     print version number\n"));
  bft_printf
    (_(" -h, --help        this help message\n\n"));
}

/*----------------------------------------------------------------------------
 * Print version number
 *----------------------------------------------------------------------------*/

static void
_print_version(void)
{
  if (cs_glob_rank_id >= 1)
    return;

  printf(_("%s version %s\n"), CS_APP_NAME, CS_APP_VERSION);
}

/*----------------------------------------------------------------------------
 * Define options and call some associated initializations
 * based on command line arguments
 *
 * parameters:
 *   argc     <-- number of command line arguments
 *   argv     <-- array of command line arguments
 *   alg_opt  --> 1 for METIS, 2 for SCOTCH
 *   no_write --> do not write output if 1
 *   no_perio --> ignore periodicity information if 1
 *   n_parts  --> number of partitionings required
 *   n_ranks  --> array of required ranks per partitioning
 *----------------------------------------------------------------------------*/

static void
_define_options(int     argc,
                char   *argv[],
                int    *alg_opt,
                int    *no_write,
                int    *no_perio,
                int    *n_parts,
                int   **n_ranks)
{
  /* Local variables */

  int _n_ranks;
  const char *s;
  int arg_id = 0, argerr = 0;

  /* Default initialization */

  int ilisr0 = 1;
  int ilisrp = 2;

  int mpi_io_mode = -1;

  /* Initialize return arguments */

  *no_write = 0;
  *no_perio  = 0;
  *n_parts   = 0;
  *n_ranks   = NULL;

  *alg_opt = 0;

  if (cs_glob_n_ranks == 1) {
#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  *alg_opt = 1;
#endif
#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)
  if (*alg_opt == 0)
    *alg_opt = 2;
#endif
  }
  else { /* if cs_glob_n_ranks > 1) */
#if defined(HAVE_PARMETIS)
    *alg_opt = 1;
#endif
#if defined(HAVE_PTSCOTCH)
    if (*alg_opt == 0)
      *alg_opt = 2;
#endif
  }

  /* Parse command line arguments */

  while (++arg_id < argc && argerr == 0) {

    s = argv[arg_id];

    /* Version number */

    if (strcmp(s, "-v") == 0 || strcmp(s, "--version") == 0)
      argerr = 3;

    /* Usage */

    else if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0)
      argerr = 2;

#if defined(HAVE_MPI)

    else if (strcmp(s, "--mpi") == 0) {
    }

    else if (strcmp(s, "--mpi-io") == 0) {
      if (arg_id + 1 < argc) {
        const char *s_n = argv[arg_id + 1];
        if (strcmp(s_n, "off") == 0)
          mpi_io_mode = 0;
        else if (strcmp(s_n, "eo") == 0)
          mpi_io_mode = 1;
        else if (strcmp(s_n, "ip") == 0)
          mpi_io_mode = 2;
        else
          argerr = 1;
        if (argerr == 0)
          arg_id++;
      }
      else
        argerr = 1;
    }

#else /* !defined(HAVE_MPI) */

    else if (strcmp(s, "--mpi") == 0 || strcmp(s, "--mpi-io") == 0) {
      fprintf(stderr, _("%s was built without MPI support,\n"
                        "so option \"%s\" may not be used.\n"),
              argv[0], s);
      cs_exit(EXIT_FAILURE);
    }

#endif /* defined(HAVE_MPI) */

    else if (strcmp(s, "--log") == 0) {
      int n1 = 0;
      n1 = _arg_to_int(++arg_id, argc, argv, &argerr);
      if (n1 == 0)
        ilisr0 = 0;
      else if (n1 == 1)
        ilisr0 = 1;
      else
        argerr = 1;
    }

    else if (strcmp(s, "--logp") == 0) {
      int n1 = 0;
      n1 = _arg_to_int(++arg_id, argc, argv, &argerr);
      if (n1 == -1)
        ilisrp = 2;
      else if (n1 == 0)
        ilisrp = 0;
      else if (n1 == 1)
        ilisrp = 1;
      else
        argerr = 1;
    }

    else if (strcmp(argv[arg_id], "--no-write") == 0)
      *no_write = 1;

    else if (strcmp(argv[arg_id], "--no-perio") == 0)
      *no_perio = 1;

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

    else if (strcmp(argv[arg_id], "--metis") == 0)
      *alg_opt = 1;

#else /* !defined(HAVE_METIS) && ! defined(HAVE_PARMETIS) */

    else if (strcmp(s, "--metis") == 0) {
      fprintf(stderr, _("%s was built without METIS or ParMETIS support,\n"
                        "so option \"%s\" may not be used.\n"),
              argv[0], s);
      cs_exit(EXIT_FAILURE);
    }

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

    else if (strcmp(argv[arg_id], "--scotch") == 0)
      *alg_opt = 2;

#else /* !defined(HAVE_SCOTCH) && ! defined(HAVE_PTSCOTCH) */

    else if (strcmp(s, "--scotch") == 0) {
      fprintf(stderr, _("%s was built without SCOTCH or PT-SCOTCH support,\n"
                        "so option \"%s\" may not be used.\n"),
              argv[0], s);
      cs_exit(EXIT_FAILURE);
    }

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

    else {
      _n_ranks = atoi(argv[arg_id]);

      if (_n_ranks <= 1)
        argerr = 1;

      else {
        BFT_REALLOC(*n_ranks, *n_parts + 1, int);
        (*n_ranks)[*n_parts] = _n_ranks;
        *n_parts += 1;
      }

    }
  } /* End parsing command line */

  /* Print version/help and exit if required or in case of command line error */

  if (argerr != 0) {
    if (cs_glob_rank_id <= 0) {
      switch (argerr) {
      case 1:
      case 2:
        cs_base_logfile_head(argc, argv);
        _arg_env_help(argv[0]);
        break;
      case 3:
        _print_version();
        break;
      default:
        break;
      }
    }
    if (argerr == 1)
      cs_exit(EXIT_FAILURE);
    else
      cs_exit(EXIT_SUCCESS);
  }

  /* Define logging and MPI-IO options */

  cs_base_bft_printf_set("partition.log", ilisr0, ilisrp);

  cs_io_set_defaults(mpi_io_mode);
}

/*----------------------------------------------------------------------------
 * Display the distribution of cells per partition in serial mode
 *
 * parameters:
 *   n_cells  <-- number of cells
 *   n_parts  <-- number of partitions
 *   part     <-- cell partition number
 *----------------------------------------------------------------------------*/

static void
_cell_part_histogram(size_t      n_cells,
                     int         n_parts,
                     const   int part[])
{
  int i, k;
  size_t j;
  double step;

  fvm_lnum_t *n_part_cells;
  fvm_lnum_t n_min, n_max;
  size_t count[10];
  int n_steps = 10;

  if (n_parts <= 1) /* Should never happen */
    return;

  bft_printf(_("  Number of cells per domain (histogramm):\n"));

  BFT_MALLOC(n_part_cells, n_parts, fvm_lnum_t);

  for (i = 0; i < n_parts; i++)
    n_part_cells[i] = 0;

  for (j = 0; j < n_cells; j++)
    n_part_cells[part[j] - 1] += 1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    fvm_lnum_t *n_part_cells_sum;
    BFT_MALLOC(n_part_cells_sum, n_parts, fvm_lnum_t);
    MPI_Allreduce(n_part_cells, n_part_cells_sum, n_parts,
                  FVM_MPI_LNUM, MPI_SUM, cs_glob_mpi_comm);
    BFT_FREE(n_part_cells);
    n_part_cells = n_part_cells_sum;
    n_part_cells_sum = NULL;
  }

#endif /* defined(HAVE_MPI) */

  /* Compute min and max */

  n_min = n_part_cells[0];
  n_max = n_part_cells[0];

  for (i = 1; i < n_parts; i++) {
    if (n_part_cells[i] > n_max)
      n_max = n_part_cells[i];
    else if (n_part_cells[i] < n_min)
      n_min = n_part_cells[i];
  }

  /* Define axis subdivisions */

  for (i = 0; i < n_steps; i++)
    count[i] = 0;

  if (n_max - n_min > 0) {

    if (n_max-n_min < n_steps)
      n_steps = n_max-n_min > 0 ? n_max-n_min : 1;

    step = (double)(n_max - n_min) / n_steps;

    /* Loop on partitions */

    for (i = 0; i < n_parts; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (n_part_cells[i] < n_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    [ %10d ; %10d [ = %10d\n",
                 (int)(n_min + i*step),
                 (int)(n_min + j*step),
                 (int)(count[i]));

    bft_printf("    [ %10d ; %10d ] = %10d\n",
               (int)(n_min + (n_steps - 1)*step),
               (int)n_max,
               (int)(count[n_steps - 1]));

  }

  else { /* if (n_max == n_min) */
    bft_printf("    [ %10d ; %10d ] = %10d\n",
               (int)(n_min), (int)n_max, (int)n_parts);
  }

  BFT_FREE(n_part_cells);
}

/*----------------------------------------------------------------------------
 * Read a global count type section from the input mesh.
 *
 * parameters:
 *   inp    <-- input file pointer
 *   rec_id <-- record id
 *   header <-> record header
 *
 * returns: value associated with section.
 *----------------------------------------------------------------------------*/

static fvm_gnum_t
_read_global_count(cs_io_t             *inp,
                   int                  rec_id,
                   cs_io_sec_header_t  *header)
{
  fvm_gnum_t retval;

  assert(inp != NULL);

  /* Check data consistency */

  if (header->n_location_vals != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Section: \"%s\" of file \"%s\" has  %d values per\n"
                "location where 0 were expected."),
              header->sec_name, cs_io_get_name(inp), header->n_location_vals);

  /* Ensure type of value matches */

  if (header->elt_type == FVM_UINT32 || header->elt_type == FVM_UINT64)
    cs_io_set_fvm_gnum(header, inp);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Section: \"%s\" of file \"%s\" is not\n"
                "of unsigned integer type."),
              header->sec_name, cs_io_get_name(inp));

  /* Now set position in file to read data */

  cs_io_set_indexed_position(inp, header, rec_id);

  cs_io_read_global(header, &retval, inp);

  return retval;
}

/*----------------------------------------------------------------------------
 * Read and global array section from the input mesh.
 *
 * Each rank may read a different portion of the array, base on the
 * given block info structure.
 *
 * parameters:
 *   inp       <-- input file pointer
 *   rank_id   <-- MPI rank id, or 0
 *   n_ranks   <-- number of MPI ranks, or 1
 *   rec_id    <-- record id
 *   bi        <-- associated block information
 *   header    <-> record header
 *   adjacency --> adjacency array (pre-allocated)
 *
 * returns: pointer to array allocated and read.
 *----------------------------------------------------------------------------*/

static void
_read_adjacency_array(cs_io_t                         *inp,
                      int                              rec_id,
                      const fvm_block_to_part_info_t  *bi,
                      cs_io_sec_header_t              *header,
                      fvm_gnum_t                      *adjacency)
{
  assert(inp != NULL);

  /* Check data consistency */

  if (header->n_location_vals != 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Section: \"%s\" of file \"%s\" has  %d values per\n"
                "location where 2 were expected."),
              header->sec_name, cs_io_get_name(inp), header->n_location_vals);

  /* Ensure type of value matches */

  if (   header->elt_type == FVM_INT32 || header->elt_type == FVM_INT64
      || header->elt_type == FVM_UINT32 || header->elt_type == FVM_UINT64)
    cs_io_set_fvm_gnum(header, inp);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Section: \"%s\" of file \"%s\" is not\n"
                "of integer type."),
              header->sec_name, cs_io_get_name(inp));

  /* Now set position in file and read data */

  cs_io_set_indexed_position(inp, header, rec_id);
  cs_io_set_fvm_gnum(header, inp);

  cs_io_read_block(header,
                   bi->gnum_range[0],
                   bi->gnum_range[1],
                   adjacency,
                   inp);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Update global face -> cells adjacency for periodicity in parallel mode.
 *
 * If face Fi is periodic with face Fj, and face Fi is adjacent to cell Ci,
 * while face Fj is adjacent to face Cj, we add Cj to Fi's adjacent cells,
 * and Ci to Fj's adjacent cells, just as if Fi and Fj were regular interior
 * faces (this ignores the geometric transformation, but is suffient to
 * build the cells -> cells graph for domain partitioning).
 *
 * This function should be called when faces are distributed by blocks,
 * that is prior to redistribution based on face -> cells adjacency.
 *
 * arguments:
 *   bi                 <-- block size and range info for faces
 *   g_face_cells       <-> global face->cells adjacency to update
 *   n_periodic_couples <-- number of periodic couples
 *   periodic_couples   <-- array indicating periodic couples (interlaced,
 *                          using global numberings)
 *   comm               <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_add_perio_to_face_cells_g(fvm_block_to_part_info_t  bi,
                           fvm_gnum_t                g_face_cells[],
                           fvm_lnum_t                n_periodic_couples,
                           const fvm_gnum_t          periodic_couples[],
                           MPI_Comm                  comm)
{
  int i;
  fvm_lnum_t j;
  size_t k;

  int n_ranks;
  int *send_count = NULL, *recv_count = NULL;
  int *send_displ = NULL, *recv_displ = NULL;
  fvm_gnum_t *send_buf = NULL, *recv_buf = NULL, *send_adj = NULL;

  size_t n_send = 0, n_recv = 0;

  const int rank_step = bi.rank_step;
  const fvm_lnum_t block_size = bi.block_size;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  BFT_MALLOC(send_count, n_ranks, int);
  BFT_MALLOC(recv_count, n_ranks, int);
  BFT_MALLOC(send_displ, n_ranks, int);
  BFT_MALLOC(recv_displ, n_ranks, int);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  /* Count number of values to send to each rank */

  for (j = 0; j < n_periodic_couples; j++) {
    int rank_0 = ((periodic_couples[j*2] -1)/block_size) * rank_step;
    int rank_1 = ((periodic_couples[j*2+1] -1)/block_size) * rank_step;
    send_count[rank_0] += 1;
    send_count[rank_1] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_displ[0] = 0;
  recv_displ[0] = 0;
  for (i = 1; i < n_ranks; i++) {
    send_displ[i] = send_displ[i-1] + send_count[i-1];
    recv_displ[i] = recv_displ[i-1] + recv_count[i-1];
  }
  n_send = send_displ[n_ranks - 1] + send_count[n_ranks - 1];
  n_recv = recv_displ[n_ranks - 1] + recv_count[n_ranks - 1];

  /* Now prepare data for request to matching cell */

  BFT_MALLOC(send_buf, n_send, fvm_gnum_t);
  BFT_MALLOC(recv_buf, n_recv, fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (j = 0; j < n_periodic_couples; j++) {
    fvm_gnum_t f_id_0 = periodic_couples[j*2] - 1;
    fvm_gnum_t f_id_1 = periodic_couples[j*2 + 1] - 1;
    int rank_0 = (f_id_0/block_size) * rank_step;
    int rank_1 = (f_id_1/block_size) * rank_step;
    send_buf[send_displ[rank_0] + send_count[rank_0]] = f_id_0 + 1;
    send_buf[send_displ[rank_1] + send_count[rank_1]] = f_id_1 + 1;
    send_count[rank_0] += 1;
    send_count[rank_1] += 1;
  }

  MPI_Alltoallv(send_buf, send_count, send_displ, FVM_MPI_GNUM,
                recv_buf, recv_count, recv_displ, FVM_MPI_GNUM,
                comm);

  /* Receive buffer contains global cell face whose cell adjacency
     is defined on the local rank, and we replace the received value
     by the adjacent cell number, for return exchange */

  for (k = 0; k < n_recv; k++) {
    fvm_gnum_t g_face_id = recv_buf[k] - bi. gnum_range[0];
    fvm_gnum_t c_num_0 = g_face_cells[g_face_id*2];
    const fvm_gnum_t c_num_1 = g_face_cells[g_face_id*2 + 1];
    assert(c_num_0 == 0 || c_num_1 == 0);
    /* assign c_num_0 or c_num_1 depending on which is nonzero
       (or 0 if both are 0, which should not happen) */
    recv_buf[k] = c_num_0 + c_num_1;
  }

  /* Return message (send and receive buffers inverted) */

  MPI_Alltoallv(recv_buf, recv_count, recv_displ, FVM_MPI_GNUM,
                send_buf, send_count, send_displ, FVM_MPI_GNUM,
                comm);

  /* Now send_buf contains the global cell number matching a given face */

  BFT_MALLOC(send_adj, n_send*2, fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (j = 0; j < n_periodic_couples; j++) {

    fvm_gnum_t f_id_0 = periodic_couples[j*2] - 1;
    fvm_gnum_t f_id_1 = periodic_couples[j*2 + 1] - 1;
    int rank_0 = (f_id_0/block_size) * rank_step;
    int rank_1 = (f_id_1/block_size) * rank_step;

    /* Retrieve adjacent cell number from previous exchange */
    fvm_gnum_t c_num_0 = send_buf[send_displ[rank_0] + send_count[rank_0]];
    fvm_gnum_t c_num_1 = send_buf[send_displ[rank_1] + send_count[rank_1]];

    /* Send global face number and cell number adjacent with its
       periodic face to rank defining the adjacency for this face */

    send_adj[(send_displ[rank_0] + send_count[rank_0])*2] = f_id_0 + 1;
    send_adj[(send_displ[rank_0] + send_count[rank_0])*2+1] = c_num_1;

    send_adj[(send_displ[rank_1] + send_count[rank_1])*2] = f_id_1 + 1;
    send_adj[(send_displ[rank_1] + send_count[rank_1])*2+1] = c_num_0;

    send_count[rank_0] += 1;
    send_count[rank_1] += 1;
  }

  /* Adjust counts and disps */

  for (i = 1; i < n_ranks; i++) {
    send_count[i] *= 2;
    recv_count[i] *= 2;
    send_displ[i] *= 2;
    recv_displ[i] *= 2;
  }

  BFT_FREE(send_buf);
  BFT_REALLOC(recv_buf, n_recv*2, fvm_gnum_t);

  MPI_Alltoallv(send_adj, send_count, send_displ, FVM_MPI_GNUM,
                recv_buf, recv_count, recv_displ, FVM_MPI_GNUM,
                comm);

  BFT_FREE(send_adj);
  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_displ);
  BFT_FREE(recv_displ);

  /* Update face -> cell connectivity */

  for (k = 0; k < n_recv; k++) {
    const fvm_gnum_t g_face_id = recv_buf[k*2] - bi. gnum_range[0];
    const fvm_gnum_t g_cell_num = recv_buf[k*2 + 1];
    if (g_face_cells[g_face_id*2] == 0)
      g_face_cells[g_face_id*2] = g_cell_num;
    else {
      assert(g_face_cells[g_face_id*2 + 1] == 0);
      g_face_cells[g_face_id*2 + 1] = g_cell_num;
    }
  }

  BFT_FREE(recv_buf);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Update global face -> cells adjacency for periodicity in local mode.
 *
 * If face Fi is periodic with face Fj, and face Fi is adjacent to cell Ci,
 * while face Fj is adjacent to face Cj, we add Cj to Fi's adjacent cells,
 * and Ci to Fj's adjacent cells, just as if Fi and Fj were regular interior
 * faces (this ignores the geometric transformation, but is suffient to
 * build the cells -> cells graph for domain partitioning).
 *
 * This function should be called when faces are distributed by blocks,
 * that is prior to redistribution based on face -> cells adjacency.
 *
 * arguments:
 *   bi                 <-- block size and range info for faces
 *   g_face_cells       <-> global face->cells adjacency to update
 *   n_periodic_couples <-- number of periodic couples
 *   periodic_couples   <-- array indicating periodic couples (interlaced,
 *                          using global numberings)
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

static void
_add_perio_to_face_cells_l(fvm_gnum_t        g_face_cells[],
                           fvm_lnum_t        n_periodic_couples,
                           const fvm_gnum_t  periodic_couples[])
{
  fvm_lnum_t j;

  /* Update face -> cell connectivity */

  for (j = 0; j < n_periodic_couples; j++) {

    fvm_gnum_t f_id_0 = periodic_couples[j*2] - 1;
    fvm_gnum_t f_id_1 = periodic_couples[j*2 + 1] - 1;
    fvm_gnum_t c_num_00 = g_face_cells[f_id_0*2];
    fvm_gnum_t c_num_01 = g_face_cells[f_id_0*2 + 1];
    fvm_gnum_t c_num_10 = g_face_cells[f_id_1*2];
    fvm_gnum_t c_num_11 = g_face_cells[f_id_1*2 + 1];

    assert(c_num_00 == 0 || c_num_01 == 0);
    assert(c_num_10 == 0 || c_num_11 == 0);

    /* assign c_num_i0 or c_num_i1 depending on which is nonzero
       (or 0 if both are 0, which should not happen) */

    if (g_face_cells[f_id_0*2] == 0)
      g_face_cells[f_id_0*2] = c_num_10 + c_num_11;
    else
      g_face_cells[f_id_0*2 + 1] = c_num_00 + c_num_01;

    if (g_face_cells[f_id_1*2] == 0)
      g_face_cells[f_id_1*2] = c_num_00 + c_num_01;
    else
      g_face_cells[f_id_1*2 + 1] = c_num_10 + c_num_11;
  }
}

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

/*----------------------------------------------------------------------------
 * Build cell -> cell connectivity
 *
 * parameters:
 *   n_cells        <-- number of cells in mesh
 *   n_faces        <-- number of cells in mesh
 *   start_cell     <-- number of first cell for the curent rank
 *   face_cells     <-- face->cells connectivity
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *----------------------------------------------------------------------------*/

static void
_metis_cell_cells(size_t        n_cells,
                  size_t        n_faces,
                  fvm_gnum_t    start_cell,
                  fvm_gnum_t   *face_cells,
                  idxtype     **cell_idx,
                  idxtype     **cell_neighbors)
{
  size_t i, id_0, id_1;

  fvm_gnum_t  c_num[2];

  idxtype  *n_neighbors;
  idxtype  *_cell_idx;
  idxtype  *_cell_neighbors;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, idxtype);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    assert(   c_num[0] - start_cell < n_cells
           || c_num[1] - start_cell < n_cells);

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells)
      n_neighbors[c_num[0] - start_cell] += 1;
    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells)
      n_neighbors[c_num[1] - start_cell] += 1;
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, idxtype);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++)
    _cell_idx[i + 1] = _cell_idx[i] + n_neighbors[i];

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_cells], idxtype);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells) {
      id_0 = c_num[0] - start_cell;
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = c_num[1] - 1;
      n_neighbors[id_0] += 1;
    }

    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells) {
      id_1 = c_num[1] - start_cell;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = c_num[0] - 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_FREE(n_neighbors);

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
}

/*----------------------------------------------------------------------------
 * Compute partition using METIS
 *
 * parameters:
 *   n_cells       <-- number of cells in mesh
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *----------------------------------------------------------------------------*/

static void
_part_metis(size_t    n_cells,
            int       n_parts,
            idxtype  *cell_idx,
            idxtype  *cell_neighbors,
            int      *cell_part)
{
  size_t i;
  double  start_time[2], end_time[2];

  int     wgtflag    = 0; /* No weighting for faces or cells */
  int     numflag    = 0; /* 0 to n-1 numbering (C type) */
  int     options[5] = {0, 3, 1, 1, 0}; /* By default if options[0] = 0 */
  int     edgecut    = 0; /* <-- Number of faces on partition */

  int       _n_cells = n_cells;
  idxtype  *_cell_part = NULL;

  start_time[0] = bft_timer_wtime();
  start_time[1] = bft_timer_cpu_time();

  if (sizeof(idxtype) == sizeof(int))
    _cell_part = cell_part;

  else
    BFT_MALLOC(_cell_part, n_cells, idxtype);

  if (n_parts < 8) {

    bft_printf(_("\n"
                 "Partitioning %d cells to %d domains"
                 " (METIS_PartGraphRecursive).\n"),
               (int)n_cells, n_parts);

    METIS_PartGraphRecursive(&_n_cells,
                             cell_idx,
                             cell_neighbors,
                             NULL,       /* vwgt:   cell weights */
                             NULL,       /* adjwgt: face weights */
                             &wgtflag,
                             &numflag,
                             &n_parts,
                             options,
                             &edgecut,
                             _cell_part);

  }

  else {

    bft_printf(_("\n"
                 "Partitioning %d cells to %d domains"
                 " (METIS_PartGraphKway).\n"),
               (int)n_cells, n_parts);

    METIS_PartGraphKway(&_n_cells,
                        cell_idx,
                        cell_neighbors,
                        NULL,       /* vwgt:   cell weights */
                        NULL,       /* adjwgt: face weights */
                        &wgtflag,
                        &numflag,
                        &n_parts,
                        options,
                        &edgecut,
                        _cell_part);

  }

  end_time[0] = bft_timer_wtime();
  end_time[1] = bft_timer_cpu_time();

  bft_printf(_("\n"
               "  Total number of faces on parallel boundaries: %llu\n"
               "  wall-clock time: %f s; CPU time: %f s\n\n"),
             (unsigned long long)edgecut,
             (double)(end_time[0] - start_time[0]),
             (double)(end_time[1] - start_time[1]));

  if (sizeof(idxtype) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
  for (i = 0; i < n_cells; i++)
    cell_part[i] += 1;

  _cell_part_histogram(n_cells, n_parts, cell_part);
}

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_PARMETIS)

/*----------------------------------------------------------------------------
 * Compute partition using ParMETIS
 *
 * parameters:
 *   n_g_cells     <-- global number of cells
 *   cell_range    <-- first and past-the-last cell numbers for this rank
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_part_parmetis(fvm_gnum_t  n_g_cells,
               fvm_gnum_t  cell_range[2],
               int         n_parts,
               idxtype    *cell_idx,
               idxtype    *cell_neighbors,
               int        *cell_part,
               MPI_Comm    comm)
{
  size_t i;
  double  start_time, end_time;

  int     wgtflag    = 0; /* No weighting for faces or cells */
  int     numflag    = 0; /* 0 to n-1 numbering (C type) */
  int     ncon       = 0; /* number of weights for each vertex */
  int     options[3] = {0, 1, 15}; /* By default if options[0] = 0 */
  int     edgecut    = 0; /* <-- Number of faces on partition */

  int       n_ranks;
  size_t    n_cells = cell_range[1] - cell_range[0];
  idxtype   vtxstart = cell_range[0] - 1;
  idxtype   vtxend = cell_range[1] - 1;
  idxtype  *vtxdist = NULL;
  idxtype  *_cell_part = NULL;
  MPI_Datatype mpi_idxtype = MPI_DATATYPE_NULL;

  start_time = bft_timer_wtime();

  MPI_Comm_size(comm, &n_ranks);

  /* Adjust mpi_idxtype if necessary */

  if (sizeof(idxtype) == sizeof(short))
    mpi_idxtype = MPI_SHORT;
  else if (sizeof(idxtype) == sizeof(int))
    mpi_idxtype = MPI_INT;
  else if (sizeof(idxtype) == sizeof(long))
    mpi_idxtype = MPI_LONG; /* standard ParMETIS 3.1.1 only short or int */
  else {
    assert(0); /* porting error, possible with future or modified ParMETIS */
  }

  if (sizeof(idxtype) == sizeof(int))
    _cell_part = cell_part;

  else
    BFT_MALLOC(_cell_part, n_cells, idxtype);

  bft_printf(_("\n"
               "Partitioning %llu cells to %d domains"
               " (ParMETIS_V3_PartKway).\n"),
             (unsigned long long)n_g_cells, n_parts);

  /* Build vtxdist */

  BFT_MALLOC(vtxdist, n_ranks + 1, idxtype);

  MPI_Allgather(&vtxstart, 1, mpi_idxtype,
                vtxdist, 1, mpi_idxtype, comm);
  MPI_Allreduce(&vtxend, vtxdist + n_ranks, 1, mpi_idxtype, MPI_MAX, comm);

  /* Call ParMETIS partitioning */

  ParMETIS_V3_PartKway
    (vtxdist,
     cell_idx,
     cell_neighbors,
     NULL,       /* vwgt:   cell weights */
     NULL,       /* adjwgt: face weights */
     &wgtflag,
     &numflag,
     &ncon,
     &n_parts,
     NULL,       /* tpwgts: size ncon, vtx weight fraction */
     NULL,       /* ubvec: size ncon, vtx imbalance */
     options,
     &edgecut,
     _cell_part,
     &comm);

  end_time = bft_timer_wtime();

  BFT_FREE(vtxdist);

  bft_printf(_("\n"
               "  Total number of faces on parallel boundaries: %llu\n"
               "  wall-clock time: %f s\n\n"),
             (unsigned long long)edgecut,
             (double)(end_time - start_time));

  if (sizeof(idxtype) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
  for (i = 0; i < n_cells; i++)
    cell_part[i] += 1;

  _cell_part_histogram(n_cells, n_parts, cell_part);
}

#endif /* defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Print an error message and exit with an EXIT_FAILURE code.
 *
 * An implementation of this function is required by libScotch.
 *
 * parameters:
 *   errstr <-- format string, as printf() and family.
 *   ...    <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

void
SCOTCH_errorPrint(const char  *errstr,
                  ...)
{
  if (cs_glob_rank_id < 1) {

    va_list  errlist;

    fflush(stdout);

    fprintf(stderr, "\n");

    fprintf(stderr, _("\nFatal error encountered by libScotch.\n\n"));

    va_start(errlist, errstr);
    vfprintf(stderr, errstr, errlist);
    va_end(errlist);
    fprintf(stderr, "\n\n");
    fflush(stderr);
  }

  assert(0);

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Abort(cs_glob_mpi_comm, EXIT_FAILURE);
  }
#endif /* HAVE_MPI */

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Print a warning message.
 *
 * An implementation of this function is required by libScotch.
 *
 * parameters:
 *   errstr <-- format string, as printf() and family.
 *   ...    <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

void
SCOTCH_errorPrintW (const char *errstr,
                    ...)
{
  if (cs_glob_rank_id < 1) {

    va_list  errlist;

    fflush(stdout);

    fprintf(stdout, "\n");

    fprintf(stdout, _("\nWarning (libScotch):\n\n"));

    va_start(errlist, errstr);
    vfprintf(stdout, errstr, errlist);
    va_end(errlist);
    fprintf(stdout, "\n\n");
    fflush(stdout);
  }
}

/*----------------------------------------------------------------------------
 * Build cell -> cell connectivity
 *
 * parameters:
 *   n_cells        <-- number of cells in mesh
 *   n_faces        <-- number of cells in mesh
 *   start_cell     <-- number of first cell for the curent rank
 *   face_cells     <-- face->cells connectivity
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *----------------------------------------------------------------------------*/

static void
_scotch_cell_cells(size_t        n_cells,
                   size_t        n_faces,
                   fvm_gnum_t    start_cell,
                   fvm_gnum_t   *face_cells,
                   SCOTCH_Num  **cell_idx,
                   SCOTCH_Num  **cell_neighbors)
{
  size_t i, id_0, id_1;

  fvm_gnum_t  c_num[2];

  SCOTCH_Num  *n_neighbors;
  SCOTCH_Num  *_cell_idx;
  SCOTCH_Num  *_cell_neighbors;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, SCOTCH_Num);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    assert(   c_num[0] - start_cell < n_cells
           || c_num[1] - start_cell < n_cells);

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells)
      n_neighbors[c_num[0] - start_cell] += 1;
    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells)
      n_neighbors[c_num[1] - start_cell] += 1;
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, SCOTCH_Num);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++)
    _cell_idx[i + 1] = _cell_idx[i] + n_neighbors[i];

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_cells], SCOTCH_Num);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells) {
      id_0 = c_num[0] - start_cell;
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = c_num[1] - 1;
      n_neighbors[id_0] += 1;
    }

    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells) {
      id_1 = c_num[1] - start_cell;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = c_num[0] - 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_FREE(n_neighbors);

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
}

/*----------------------------------------------------------------------------
 * Compute partition using SCOTCH
 *
 * parameters:
 *   n_cells       <-- number of cells in mesh
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *----------------------------------------------------------------------------*/

static void
_part_scotch(SCOTCH_Num   n_cells,
             int          n_parts,
             SCOTCH_Num  *cell_idx,
             SCOTCH_Num  *cell_neighbors,
             int         *cell_part)
{
  SCOTCH_Num  i;
  SCOTCH_Graph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;

  double  start_time[2], end_time[2];

  int     retval = 0;

  SCOTCH_Num    edgecut = 0; /* <-- Number of faces on partition */
  SCOTCH_Num  *_cell_part = NULL;

  /* Initialization */

  start_time[0] = bft_timer_wtime();
  start_time[1] = bft_timer_cpu_time();

  if (sizeof(SCOTCH_Num) == sizeof(int))
    _cell_part = cell_part;
  else
    BFT_MALLOC(_cell_part, n_cells, SCOTCH_Num);

  bft_printf(_("\n"
               "Partitioning %d cells to %d domains"
               " (SCOTCH_graphPart).\n"),
             (int)n_cells, n_parts);

  /* Partition using libScotch */

  SCOTCH_graphInit(&grafdat);

  retval
    = SCOTCH_graphBuild(&grafdat,
                        0,                  /* baseval; 0 to n -1 numbering */
                        n_cells,            /* vertnbr */
                        cell_idx,           /* verttab */
                        NULL,               /* vendtab: verttab + 1 or NULL */
                        NULL,               /* velotab: vertex weights */
                        NULL,               /* vlbltab; vertex labels */
                        cell_idx[n_cells],  /* edgenbr */
                        cell_neighbors,     /* edgetab */
                        NULL);              /* edlotab */

  if (retval == 0) {

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_graphCheck(&grafdat) == 0)
      retval = SCOTCH_graphPart(&grafdat, n_parts, &stradat, _cell_part);

    SCOTCH_stratExit(&stradat);
  }

  SCOTCH_graphExit(&grafdat);

  /* Shift cell_part values to 1 to n numbering and free possible temporary */

  if (sizeof(SCOTCH_Num) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
  for (i = 0; i < n_cells; i++)
    cell_part[i] += 1;

  /* Compute edge cut */

  if (retval == 0) {

    SCOTCH_Num cell_id, edgenum, commcut;

    commcut = 0;

    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      SCOTCH_Num  edgennd,  partval;
      partval = cell_part[cell_id];
      for (edgenum = cell_idx[cell_id], edgennd = cell_idx[cell_id + 1];
           edgenum < edgennd;
           edgenum++) {
        if (cell_part[cell_neighbors[edgenum]] != partval)
          commcut++;
      }
    }

    edgecut = commcut / 2;
  }

  /* Finalization */

  end_time[0] = bft_timer_wtime();
  end_time[1] = bft_timer_cpu_time();

  bft_printf(_("\n"
               "  Total number of faces on parallel boundaries: %llu\n"
               "  wall-clock time: %f s; CPU time: %f s\n\n"),
             (unsigned long long)edgecut,
             (double)(end_time[0] - start_time[0]),
             (double)(end_time[1] - start_time[1]));

  _cell_part_histogram(n_cells, n_parts, cell_part);
}

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

#if defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Compute partition using PT-SCOTCH
 *
 * parameters:
 *   n_g_cells     <-- global number of cells
 *   cell_range    <-- first and past-the-last cell numbers for this rank
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_part_ptscotch(fvm_gnum_t   n_g_cells,
               fvm_gnum_t   cell_range[2],
               int          n_parts,
               SCOTCH_Num  *cell_idx,
               SCOTCH_Num  *cell_neighbors,
               int         *cell_part,
               MPI_Comm     comm)
{
  SCOTCH_Num  i;
  SCOTCH_Dgraph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;

  double  start_time, end_time;

  int     retval = 0;

  SCOTCH_Num    n_cells = cell_range[1] - cell_range[0];
  SCOTCH_Num  *_cell_part = NULL;

  /* Initialization */

  start_time = bft_timer_wtime();

  if (sizeof(SCOTCH_Num) == sizeof(int))
    _cell_part = cell_part;
  else
    BFT_MALLOC(_cell_part, n_cells, SCOTCH_Num);

  bft_printf(_("\n"
               "Partitioning %llu cells to %d domains"
               " (SCOTCH_dgraphPart).\n"),
             (unsigned long long)n_g_cells, n_parts);

  /* Partition using libScotch */

  retval = SCOTCH_dgraphInit(&grafdat, comm);

  if (retval == 0) {
    retval = SCOTCH_dgraphBuild
               (&grafdat,
                0,                  /* baseval; 0 to n -1 numbering */
                n_cells,            /* vertlocnbr */
                n_cells,            /* vertlocmax (= vertlocnbr) */
                cell_idx,           /* vertloctab */
                NULL,               /* vendloctab: vertloctab + 1 or NULL */
                NULL,               /* veloloctab: vertex weights */
                NULL,               /* vlblloctab; vertex labels */
                cell_idx[n_cells],  /* edgelocnbr */
                cell_idx[n_cells],  /* edgelocsiz */
                cell_neighbors,     /* edgeloctab */
                NULL,               /* edgegstab */
                NULL);              /* edloloctab */
  }

  if (retval == 0) {

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_dgraphCheck(&grafdat) == 0)
      retval = SCOTCH_dgraphPart(&grafdat, n_parts, &stradat, _cell_part);

    SCOTCH_stratExit(&stradat);
  }

  SCOTCH_dgraphExit(&grafdat);

  /* Shift cell_part values to 1 to n numbering and free possible temporary */

  if (sizeof(SCOTCH_Num) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
  for (i = 0; i < n_cells; i++)
    cell_part[i] += 1;

  /* Finalization */

  end_time = bft_timer_wtime();

  bft_printf(_("\n"
               "  wall-clock time: %f s;\n\n"),
             (double)(end_time - start_time));

  _cell_part_histogram(n_cells, n_parts, cell_part);
}

#endif /* defined(HAVE_PTSCOTCH) */

/*----------------------------------------------------------------------------
 * Read data from a mesh file.
 *
 * parameters:
 *   path       <-- path to file
 *   cell_range <-- first and past-the-last cell numbers for this rank
 *   no_perio   <-- ignore periodicity information if 1
 *----------------------------------------------------------------------------*/

static void
_read_input(const char   *path,
            int           no_perio,
            fvm_gnum_t   *n_g_cells,
            fvm_gnum_t    cell_range[2],
            fvm_lnum_t   *n_faces,
            fvm_gnum_t  **g_face_cells)
{
  fvm_block_to_part_info_t cell_bi;
  fvm_block_to_part_info_t face_bi;

  int rank_id = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;

  int min_rank_step = 0;
  int min_block_size = 0;

  size_t index_size = 0;
  size_t rec_id;

  fvm_gnum_t n_g_faces = 0;
  fvm_gnum_t *_g_face_cells = NULL;

  fvm_lnum_t n_per_face_couples = 0;
  fvm_gnum_t n_g_per_face_couples = 0;
  fvm_gnum_t *g_per_face_couples = NULL;

#if defined(HAVE_MPI)
  fvm_block_to_part_t *d = NULL;
#endif

  cs_io_t *inp = NULL;

  const char *read_type_name[2] = {N_("Read:   "),
                                   N_("Ignored:")};

  const char *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                  "unexpected or of incorrect size.");

  /* Open mesh file */

#if defined(HAVE_MPI)
  inp = cs_io_initialize_with_index(path,
                                    "Face-based mesh definition, R0",
                                    cs_glob_io_hints,
                                    CS_IO_ECHO_OPEN_CLOSE,
                                    cs_glob_mpi_comm);
#else
  inp = cs_io_initialize_with_index(path,
                                    "Face-based mesh definition, R0",
                                    cs_glob_io_hints,
                                    CS_IO_ECHO_OPEN_CLOSE);
#endif

  bft_printf("\n");

  /* Loop over records in the index */

  index_size = cs_io_get_index_size(inp);

  for (rec_id = 0; rec_id < index_size; rec_id++) {

    int read_type = 0;
    cs_io_sec_header_t header = cs_io_get_indexed_sec_header(inp, rec_id);

    /* Number of cells or faces */

    if (strcmp("n_cells", header.sec_name) == 0) {
      *n_g_cells = _read_global_count(inp, rec_id, &header);
      cell_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                n_ranks,
                                                0,
                                                0,
                                                *n_g_cells);

    }
    else if (strcmp("n_faces", header.sec_name) == 0)
      n_g_faces = _read_global_count(inp, rec_id, &header);

    /* Face -> cells connectivity */

    else if (strcmp("face_cells", header.sec_name) == 0) {

      fvm_gnum_t n_elts = 0;

      if (header.n_vals != (fvm_file_off_t)(n_g_faces*2))
        bft_error(__FILE__, __LINE__, 0, unexpected_msg,
                  header.sec_name, cs_io_get_name(inp));
      face_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                n_ranks,
                                                min_rank_step,
                                                min_block_size,
                                                n_g_faces);

      n_elts = (face_bi.gnum_range[1] - face_bi.gnum_range[0])*2;
      BFT_MALLOC(_g_face_cells, n_elts, fvm_gnum_t);

      _read_adjacency_array(inp,
                            rec_id,
                            &face_bi,
                            &header,
                            _g_face_cells);
    }

    /* Additional connectivity in case of periodicity */

    else if (!no_perio && !strncmp("periodicity_faces_",
                                   header.sec_name,
                                   strlen("periodicity_faces_"))) {

      fvm_block_to_part_info_t per_face_bi;
      fvm_gnum_t n_per_face_couples_prev = n_per_face_couples;

      if (n_g_faces == 0)
        bft_error(__FILE__, __LINE__, 0, unexpected_msg,
                  header.sec_name, cs_io_get_name(inp));

      per_face_bi = fvm_block_to_part_compute_sizes(rank_id,
                                                    n_ranks,
                                                    min_rank_step,
                                                    min_block_size,
                                                    header.n_vals/2);

      n_per_face_couples  += (  per_face_bi.gnum_range[1]
                              - per_face_bi.gnum_range[0]);
      BFT_REALLOC(g_per_face_couples, n_per_face_couples*2, fvm_gnum_t);

      _read_adjacency_array(inp,
                            rec_id,
                            &per_face_bi,
                            &header,
                            g_per_face_couples + n_per_face_couples_prev*2);

      n_g_per_face_couples += header.n_vals/2;
    }

    else
      read_type = 1;

    /* Print info */

    if (header.n_vals > 0) {
      bft_printf(_("  %s \"%-32s\"; Type: %-6s; Size: %llu\n"),
                 _(read_type_name[read_type]),
                 header.sec_name,
                 fvm_datatype_name[header.type_read],
                 (unsigned long long)(header.n_vals));
    }
    else {
      bft_printf(_("  %s \"%-32s\"\n"),
                 _(read_type_name[read_type]), header.sec_name);
    }

  }

  /* Close mesh file */

  bft_printf("\n");
  cs_io_finalize(&inp);

  /* In case of periodicity, update global face -> cells connectivity:
     if face Fi is periodic with face Fj, and face Fi is connected
     to cell Ci, while face Fj is connected to face Cj, we add Cj
     to the cells connected to Fi, and Ci to the cells connected to Fj,
     just as if Fi and Fj were regular interior faces (this ignores
     the geometric transformation, but is suffient to build the
     cells -> cells graph for domain partitioning). */

  if (n_g_per_face_couples > 0) {

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      _add_perio_to_face_cells_g(face_bi,
                                 _g_face_cells,
                                 n_per_face_couples,
                                 g_per_face_couples,
                                 cs_glob_mpi_comm);
#endif

    if (n_ranks == 1)
      _add_perio_to_face_cells_l(_g_face_cells,
                                 n_per_face_couples,
                                 g_per_face_couples);

    BFT_FREE(g_per_face_couples);
    n_per_face_couples = 0;
  }

  /* Distribute faces if necessary */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    fvm_datatype_t gnum_type
      = (sizeof(fvm_gnum_t) == 8) ? FVM_UINT64 : FVM_UINT32;
    fvm_gnum_t *_g_face_cells_tmp = NULL;

    d = fvm_block_to_part_create_by_adj_s(cs_glob_mpi_comm,
                                          face_bi,
                                          cell_bi,
                                          2,
                                          _g_face_cells,
                                          NULL);

    *n_faces = fvm_block_to_part_get_n_part_ents(d);

    BFT_MALLOC(_g_face_cells_tmp, (*n_faces)*2, fvm_gnum_t);

    /* Face -> cell connectivity */

    fvm_block_to_part_copy_array(d,
                                 gnum_type,
                                 2,
                                 _g_face_cells,
                                 _g_face_cells_tmp);

    BFT_FREE(_g_face_cells);
    _g_face_cells = _g_face_cells_tmp;
    _g_face_cells_tmp = NULL;

    fvm_block_to_part_destroy(&d);
  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1)
    *n_faces = n_g_faces;

  /* Prepare return values */

  *g_face_cells = _g_face_cells;
  cell_range[0] = cell_bi.gnum_range[0];
  cell_range[1] = cell_bi.gnum_range[1];
}

/*----------------------------------------------------------------------------
 * Write output file.
 *
 * parameters:
 *   n_g_cells  <-- global number of cells
 *   cell_range <-- first and past-the-last cell numbers for this rank
 *   n_ranks    <-- number of ranks corresonding to output file
 *   domain_num <-- domain number array (size: n_cells)
 *----------------------------------------------------------------------------*/

static void
_write_output(fvm_gnum_t  n_g_cells,
              fvm_gnum_t  cell_range[2],
              int         n_ranks,
              int         domain_num[])
{
  size_t i;
  int n_ranks_size;
  char *filename = NULL;
  cs_io_t *fh = NULL;
  fvm_datatype_t datatype_gnum = FVM_DATATYPE_NULL;
  fvm_datatype_t datatype_int = FVM_DATATYPE_NULL;
  const char magic_string[] = "Domain partitioning, R0";

  if (sizeof(int) == 4)
    datatype_int = FVM_INT32;
  else if (sizeof(int) == 8)
    datatype_int = FVM_INT64;
  else {
    assert(0);
  }

  if (sizeof(fvm_gnum_t) == 4)
    datatype_gnum = FVM_UINT32;
  else if (sizeof(fvm_gnum_t) == 8)
    datatype_gnum = FVM_UINT64;
  else {
    assert(0);
  }

  /* Open file */

  for (i = n_ranks, n_ranks_size = 1;
       i >= 10;
       i /= 10, n_ranks_size += 1);

  BFT_MALLOC(filename,
             strlen("domain_number_") + n_ranks_size + 1,
             char);

  sprintf(filename, "domain_number_%d", n_ranks);

#if defined(HAVE_MPI)
  fh = cs_io_initialize(filename,
                        magic_string,
                        CS_IO_MODE_WRITE,
                        cs_glob_io_hints,
                        CS_IO_ECHO_OPEN_CLOSE,
                        cs_glob_mpi_comm);
#else
  fh = cs_io_initialize(filename,
                        magic_string,
                        CS_IO_MODE_WRITE,
                        cs_glob_io_hints,
                        CS_IO_ECHO_OPEN_CLOSE);
#endif

  BFT_FREE(filename);

  /* Write headers */

  cs_io_write_global("n_cells",
                     1,
                     1,
                     0,
                     1,
                     datatype_gnum,
                     &n_g_cells,
                     fh);

  cs_io_write_global("n_ranks",
                     1,
                     1,
                     0,
                     1,
                     datatype_int,
                     &n_ranks,
                     fh);

  cs_io_write_block_buffer("cell:domain number",
                           n_g_cells,
                           cell_range[0],
                           cell_range[1],
                           1,
                           0,
                           1,
                           datatype_int,
                           domain_num,
                           fh);

  cs_io_finalize(&fh);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*============================================================================
 * Main program
 *============================================================================*/

int
main(int    argc,
     char  *argv[])
{
  int  n_ext_libs = 0;
  int  app_num = -1;
  int  alg_opt = 0;
  int  no_write = 0;
  int  no_perio = 0;
  int  n_parts  = 0;
  int  *n_ranks = NULL;
  int  *cell_part = NULL;

  fvm_gnum_t  n_g_cells = 0;
  fvm_gnum_t  cell_range[2] = {0, 0};
  fvm_lnum_t  n_cells = 0;
  fvm_lnum_t  n_faces = 0;
  fvm_gnum_t  *face_cells = NULL;

  /* First analysis of the command line to determine if MPI is required,
     and MPI initialization if it is. */

#if defined(HAVE_MPI)
  app_num = cs_base_mpi_init(&argc, &argv);
#endif

  /* Default initialization */

#if defined(_CS_ARCH_Linux)

  if (getenv("LANG") != NULL)
    setlocale(LC_ALL,"");
  else
    setlocale(LC_ALL, "C");
  setlocale(LC_NUMERIC, "C");

#endif

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);
#endif

  (void)bft_timer_wtime();

  /* Initialize memory management and signals */

  cs_base_mem_init();
  cs_base_error_init();

  /* Parse command line and define I/O options */

  _define_options(argc,
                  argv,
                  &alg_opt,
                  &no_write,
                  &no_perio,
                  &n_parts,
                  &n_ranks);

  /* Log-file header and command line arguments recap */

  cs_base_logfile_head(argc, argv);

  /* Print header */

  bft_printf(_(" .------------------------------.\n"
               " |   Code_Saturne Partitioner   |\n"
               " `------------------------------'\n"));

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  n_ext_libs++;
#endif
#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)
  n_ext_libs++;
#endif

  if (n_ext_libs == 1)
    bft_printf(_("\n  External library:\n"));
  else if (n_ext_libs > 1)
    bft_printf(_("\n  External libraries:\n"));

#if defined(HAVE_METIS)
#if defined(HAVE_METIS_H) && defined(METIS_VER_MAJOR)
  bft_printf("    METIS %d.%d.%d\n",
             METIS_VER_MAJOR, METIS_VER_MINOR, METIS_VER_SUBMINOR);
#else
  bft_printf("    METIS\n");
#endif
#elif defined(HAVE_PARMETIS)
#if defined(HAVE_METIS_H) && defined(METIS_VER_MAJOR)
  bft_printf("    ParMETIS %d.%d.%d\n",
             METIS_VER_MAJOR, METIS_VER_MINOR, METIS_VER_SUBMINOR);
#else
  bft_printf("    ParMETIS\n");
#endif
#endif

#if defined(HAVE_SCOTCH)
  bft_printf("    SCOTCH\n");
#elif defined(HAVE_PTSCOTCH)
  bft_printf("    PT-SCOTCH\n");
#endif

  /* System information */

  cs_base_system_info();
  cs_io_defaults_info();

  /* Read selected data */

  _read_input("preprocessor_output",
              no_perio,
              &n_g_cells,
              cell_range,
              &n_faces,
              &face_cells);

  n_cells = cell_range[1] - cell_range[0];

  /* Build and partition graph */

  BFT_MALLOC(cell_part, n_cells, int);

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

  if (alg_opt == 1) {

    int  i;
    idxtype  *cell_idx = NULL, *cell_neighbors = NULL;

    _metis_cell_cells(n_cells,
                      n_faces,
                      cell_range[0],
                      face_cells,
                      &cell_idx,
                      &cell_neighbors);

    BFT_FREE(face_cells);

#if defined(HAVE_PARMETIS)

    if (cs_glob_n_ranks > 1) {

      for (i = 0; i < n_parts; i++) {

        _part_parmetis(n_g_cells,
                       cell_range,
                       n_ranks[i],
                       cell_idx,
                       cell_neighbors,
                       cell_part,
                       cs_glob_mpi_comm);

        if (!no_write)
          _write_output(n_g_cells, cell_range, n_ranks[i], cell_part);
      }
    }

#endif

    if (cs_glob_n_ranks == 1) {

      for (i = 0; i < n_parts; i++) {

        _part_metis(n_cells,
                    n_ranks[i],
                    cell_idx,
                    cell_neighbors,
                    cell_part);

        if (!no_write)
          _write_output(n_g_cells, cell_range, n_ranks[i], cell_part);
      }
    }

    BFT_FREE(cell_idx);
    BFT_FREE(cell_neighbors);
  }

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

  if (alg_opt == 2) {

    int  i;
    SCOTCH_Num  *cell_idx = NULL, *cell_neighbors = NULL;

    _scotch_cell_cells(n_cells,
                       n_faces,
                       cell_range[0],
                       face_cells,
                       &cell_idx,
                       &cell_neighbors);

    BFT_FREE(face_cells);

#if defined(HAVE_PTSCOTCH)

    if (cs_glob_n_ranks > 1) {

      for (i = 0; i < n_parts; i++) {

        _part_ptscotch(n_g_cells,
                       cell_range,
                       n_ranks[i],
                       cell_idx,
                       cell_neighbors,
                       cell_part,
                       cs_glob_mpi_comm);

        if (!no_write)
          _write_output(n_g_cells, cell_range, n_ranks[i], cell_part);
      }
    }

#endif

    if (cs_glob_n_ranks == 1) {

      for (i = 0; i < n_parts; i++) {

        _part_scotch(n_cells,
                     n_ranks[i],
                     cell_idx,
                     cell_neighbors,
                     cell_part);

        if (!no_write)
          _write_output(n_g_cells, cell_range, n_ranks[i], cell_part);
      }
    }

    BFT_FREE(cell_idx);
    BFT_FREE(cell_neighbors);
  }

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

  BFT_FREE(n_ranks);
  BFT_FREE(cell_part);

  /* Call main calculation initialization function or help */

  cs_base_time_summary();
  cs_base_mem_finalize();

  bft_printf(_("\n"
               " .------------------------.\n"
               " |   Partitioner finish   |\n"
               " `------------------------'\n\n"));

  /* Return */

  cs_exit(EXIT_SUCCESS);

  /* Never called, but avoids compiler warning */
  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
