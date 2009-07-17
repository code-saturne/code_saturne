/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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
 * Functions associated with code coupling.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_coupling.h>
#include <fvm_locator.h>
#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sat_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_sat_coupling_t {

  fvm_locator_t   *localis_cel;  /* Locator associated with cells */
  fvm_locator_t   *localis_fbr;  /* Locator associated with boundary faces */

  cs_int_t         nbr_cel_sup;  /* Number of associated cell locations */
  cs_int_t         nbr_fbr_sup;  /* Number of associated face locations */
  fvm_nodal_t     *cells_sup;    /* Local cells at which distant values are
                                    interpolated*/
  fvm_nodal_t     *faces_sup;    /* Local faces at which distant values are
                                    interpolated*/

  cs_real_t       *distant_dist_fbr; /* Distant vectors (distance JJ') */
  cs_real_t       *distant_of;
  cs_real_t       *local_of;
  cs_real_t       *distant_pond_fbr; /* Distant weighting coefficient */
  cs_real_t       *local_pond_fbr;   /* Local weighting coefficient */

#if defined(HAVE_MPI)

  MPI_Comm         comm;         /* Associated MPI communicator */

  cs_int_t         n_dist_ranks;    /* Number of associated distant ranks */
  cs_int_t         dist_root_rank;  /* First associated distant rank */

#endif

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Array of couplings */

static int                  cs_glob_nbr_couplages = 0;
static int                  cs_glob_nbr_couplages_max = 0;
static cs_sat_coupling_t  **cs_glob_couplages = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a coupling.
 *
 * Couplings are allowed either with process totally distinct from the
 * application communicator (cs_glob_mpi_comm), or within this same
 * communicator.
 *
 * parameters:
 *   root_rank <-- root rank of distant process leader in MPI_COMM_WORLD
 *
 * returns:
 *   pointer to new coupling structure
 *----------------------------------------------------------------------------*/

static cs_sat_coupling_t *
_sat_coupling_create(cs_int_t  root_rank)
{
  int  mpi_flag = 0;
  int  n_dist_ranks = 0;
  int  dist_root_rank = 0;
  cs_sat_coupling_t  *couplage = NULL;

  const double  tolerance = 0.1;

  /* Create associated structure and MPI communicator */

  BFT_MALLOC(couplage, 1, cs_sat_coupling_t);

#if defined(HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    couplage->comm = MPI_COMM_NULL;

  else {

    int  n_loc_ranks, n_glob_ranks, r_glob, r_loc_min, r_loc_max;

    /* Check that coupled processes overlap exactly or not at all */

    MPI_Comm_rank(MPI_COMM_WORLD, &r_glob);
    MPI_Comm_size(MPI_COMM_WORLD, &n_glob_ranks);
    MPI_Comm_size(cs_glob_mpi_comm, &n_loc_ranks);

    MPI_Allreduce(&r_glob, &r_loc_min, 1, MPI_INT, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&r_glob, &r_loc_max, 1, MPI_INT, MPI_MAX,
                  cs_glob_mpi_comm);

    if (root_rank > r_loc_min && root_rank <= r_loc_max)
      bft_error(__FILE__, __LINE__, 0,
                _("Coupling definition is impossible: a distant root rank equal to\n"
                  "%d is required, whereas the local group corresponds to\n"
                  "rank %d to %d\n"),
                (int)root_rank, r_loc_min, r_loc_max);

    else if (root_rank < 0 || root_rank >= n_glob_ranks)
      bft_error(__FILE__, __LINE__, 0,
                _("Coupling definition is impossible: a distant root rank equal to\n"
                  "%d is required, whereas the global ranks (MPI_COMM_WORLD)\n"
                  "range from to 0 to %d\n"),
                (int)root_rank, n_glob_ranks - 1);

    /* Case for a coupling internal to the process group */

    if (root_rank == r_loc_min) {
      if (n_loc_ranks == 1)
        couplage->comm = MPI_COMM_NULL;
      else
        couplage->comm = cs_glob_mpi_comm;
      n_dist_ranks = n_loc_ranks;
    }

    /* Case for a coupling external to the process group */

    else {

      int local_range[2] = {-1, -1};
      int distant_range[2] = {-1, -1};

      fvm_coupling_mpi_intracomm_create(cs_glob_mpi_comm,
                                        root_rank,
                                        &(couplage->comm),
                                        local_range,
                                        distant_range);

      bft_printf(_("coupling: local_ranks = [%d..%d], distant ranks = [%d..%d]\n"),
                 local_range[0], local_range[1] - 1,
                 distant_range[0], distant_range[1] - 1);

      n_dist_ranks = distant_range[1] - distant_range[0];
      dist_root_rank = distant_range[0];
    }

  }

  couplage->n_dist_ranks = n_dist_ranks;
  couplage->dist_root_rank = dist_root_rank;

#endif

  /* Creation of the localization structures */

#if defined(FVM_HAVE_MPI)

  couplage->localis_cel = fvm_locator_create(tolerance,
                                             couplage->comm,
                                             n_dist_ranks,
                                             dist_root_rank);

  couplage->localis_fbr = fvm_locator_create(tolerance,
                                             couplage->comm,
                                             n_dist_ranks,
                                             dist_root_rank);

#else

  couplage->localis_cel = fvm_locator_create(tolerance);
  couplage->localis_fbr = fvm_locator_create(tolerance);

#endif

  couplage->nbr_cel_sup = 0;
  couplage->nbr_fbr_sup = 0;
  couplage->cells_sup = NULL;
  couplage->faces_sup = NULL;

  couplage->distant_dist_fbr = NULL;
  couplage->distant_of       = NULL;
  couplage->local_of         = NULL;
  couplage->distant_pond_fbr = NULL;
  couplage->local_pond_fbr   = NULL;

  return couplage;
}

/*----------------------------------------------------------------------------
 * Destroy a coupling structure
 *
 * parameters:
 *   couplage <-> pointer to coupling structure to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static cs_sat_coupling_t *
_sat_coupling_destroy(cs_sat_coupling_t  *couplage)
{
  fvm_locator_destroy(couplage->localis_cel);
  fvm_locator_destroy(couplage->localis_fbr);

  if (couplage->cells_sup != NULL)
    fvm_nodal_destroy(couplage->cells_sup);
  if (couplage->faces_sup != NULL)
    fvm_nodal_destroy(couplage->faces_sup);

  BFT_FREE(couplage->distant_dist_fbr);
  BFT_FREE(couplage->distant_of);
  BFT_FREE(couplage->local_of);
  BFT_FREE(couplage->distant_pond_fbr);
  BFT_FREE(couplage->local_pond_fbr);

#if defined(HAVE_MPI)
  if (   couplage->comm != MPI_COMM_WORLD
      && couplage->comm != cs_glob_mpi_comm)
    MPI_Comm_free(&(couplage->comm));
#endif

  BFT_FREE(couplage);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Computed some quantities needed for a centred-like interpolation
 *  - distance JJ' for distant boundary faces
 *  - local weighting coefficients
 *----------------------------------------------------------------------------*/

static void
_sat_coupling_interpolate(cs_sat_coupling_t  *couplage)
{
  int    icoo;
  int    reverse;

  cs_int_t    ind;
  cs_int_t    iel;
  cs_int_t    ifac;

  cs_int_t    n_fbr_loc  = 0;
  cs_int_t    n_fbr_dist = 0;

  cs_real_t   pdt_scal;
  cs_real_t   surface;

  cs_real_t   distance_fbr_cel;
  cs_real_t   distance_cel_cel;

  cs_real_t   dist_cel_fbr[3];
  cs_real_t   vect_surf_norm[3];

  cs_real_t  *local_surf     = NULL;
  cs_real_t  *local_xyzcen   = NULL;
  cs_real_t  *distant_surf   = NULL;
  cs_real_t  *distant_xyzcen = NULL;

  const fvm_lnum_t   *lstfbr        = NULL;
  const fvm_lnum_t   *element       = NULL;
  const fvm_coord_t  *distant_coord = NULL;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;


  /* Removing the connectivity and localization informations in case of
     coupling update */

  if (couplage->distant_dist_fbr != NULL)
    BFT_FREE(couplage->distant_dist_fbr);
  if (couplage->distant_of != NULL)
    BFT_FREE(couplage->distant_of);
  if (couplage->local_of != NULL)
    BFT_FREE(couplage->local_of);
  if (couplage->distant_pond_fbr != NULL)
    BFT_FREE(couplage->local_pond_fbr);
  if (couplage->local_pond_fbr != NULL)
    BFT_FREE(couplage->local_pond_fbr);


  /* Interpolation structure */

  n_fbr_loc  = fvm_locator_get_n_interior(couplage->localis_fbr);
  lstfbr     = fvm_locator_get_interior_list(couplage->localis_fbr);

  n_fbr_dist    = fvm_locator_get_n_dist_points(couplage->localis_fbr);
  element       = fvm_locator_get_dist_locations(couplage->localis_fbr);
  distant_coord = fvm_locator_get_dist_coords(couplage->localis_fbr);


  /* Calculation of the distance DJJPB defining the distance from */
  /* the local cell centre to the distant boundary face norm      */
  /*--------------------------------------------------------------*/

  BFT_MALLOC(couplage->distant_dist_fbr, 3*n_fbr_dist, cs_real_t);

  /* Store the local surface vector of the coupled boundary faces */

  BFT_MALLOC(local_surf, 3*n_fbr_loc, cs_real_t);

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    ifac = lstfbr[ind] - 1;

    for (icoo = 0 ; icoo < 3 ; icoo++)
      local_surf[ind*3 + icoo] = mesh_quantities->b_face_normal[ifac*3 + icoo];

 }

  /* Get the distant faces surface vector (reverse = 1) */

  reverse = 1;

  BFT_MALLOC(distant_surf, 3*n_fbr_dist, cs_real_t);

  fvm_locator_exchange_point_var(couplage->localis_fbr,
                                 distant_surf,
                                 local_surf,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(local_surf);

  /* Calculation of the JJ' vectors */

  BFT_MALLOC(distant_xyzcen, 3*n_fbr_dist, cs_real_t);

  for (ind = 0; ind < n_fbr_dist; ind++) {

    iel = element[ind] - 1;

    surface = 0.;
    for (icoo = 0; icoo < 3; icoo++)
      surface += distant_surf[ind*3 + icoo]*distant_surf[ind*3 + icoo];
    surface = sqrt(surface);

    pdt_scal = 0.;
    for (icoo = 0; icoo < 3; icoo++) {

      dist_cel_fbr[icoo] =
        distant_coord[ind*3 + icoo] - mesh_quantities->cell_cen[iel*3 + icoo];

      /* Store the distant coordinates to compute the weighting coefficients */
      distant_xyzcen[ind*3 + icoo] = mesh_quantities->cell_cen[iel*3 + icoo];

      vect_surf_norm[icoo] =
        distant_surf[ind*3 + icoo] / surface;

      pdt_scal += dist_cel_fbr[icoo]*vect_surf_norm[icoo];

    }

    for (icoo = 0; icoo < 3; icoo++)
      couplage->distant_dist_fbr[ind*3 + icoo] =
        dist_cel_fbr[icoo] - pdt_scal*vect_surf_norm[icoo];

  }

  BFT_FREE(distant_surf);


  /* Calculation of the local weighting coefficient */
  /*------------------------------------------------*/

  BFT_MALLOC(couplage->distant_pond_fbr, n_fbr_dist, cs_real_t);
  BFT_MALLOC(couplage->local_pond_fbr, n_fbr_loc, cs_real_t);

  /* Get the cell centres coordinates (reverse = 0) */

  reverse = 0;

  BFT_MALLOC(local_xyzcen, 3*n_fbr_loc, cs_real_t);

  fvm_locator_exchange_point_var(couplage->localis_fbr,
                                 distant_xyzcen,
                                 local_xyzcen,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(distant_xyzcen);

  /* Calculation of the local weighting coefficients */

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    ifac = lstfbr[ind] - 1;
    iel  = mesh->b_face_cells[ifac] - 1;

    surface = 0.;

    distance_fbr_cel = 0.;
    distance_cel_cel = 0.;

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      surface += mesh_quantities->b_face_normal[ifac*3 + icoo]
        * mesh_quantities->b_face_normal[ifac*3 + icoo];

      distance_fbr_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3 + icoo] - mesh_quantities->b_face_cog[ifac*3 + icoo]);

      distance_cel_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3 + icoo] - mesh_quantities->cell_cen[iel*3 + icoo]);

    }

    surface = sqrt(surface);

    distance_fbr_cel /= surface;
    distance_cel_cel /= surface;

    if (fabs(distance_cel_cel) > 1.e-12)
      couplage->local_pond_fbr[ind] = distance_fbr_cel / distance_cel_cel;
    else
      couplage->local_pond_fbr[ind] = 0.5;

  }


  /* Get the distant weighting coefficients (reverse = 1) */

  reverse = 1;

  fvm_locator_exchange_point_var(couplage->localis_fbr,
                                 couplage->distant_pond_fbr,
                                 couplage->local_pond_fbr,
                                 NULL,
                                 sizeof(cs_real_t),
                                 1,
                                 reverse);



  /* Calculation of the OF distance */
  /*--------------------------------*/

  BFT_MALLOC(couplage->distant_of, 3*n_fbr_dist, cs_real_t);
  BFT_MALLOC(couplage->local_of, 3*n_fbr_loc, cs_real_t);

  for (ind = 0 ; ind < n_fbr_loc ; ind++) {

    cs_real_t pond = couplage->local_pond_fbr[ind];

    ifac = lstfbr[ind] - 1;
    iel  = mesh->b_face_cells[ifac] - 1;

    surface = 0.;

    distance_fbr_cel = 0.;
    distance_cel_cel = 0.;

    for (icoo = 0 ; icoo < 3 ; icoo++) {

      surface += mesh_quantities->b_face_normal[ifac*3 + icoo]
        * mesh_quantities->b_face_normal[ifac*3 + icoo];

      distance_fbr_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3+icoo] - mesh_quantities->b_face_cog[ifac*3+icoo]);

      distance_cel_cel += mesh_quantities->b_face_normal[ifac*3 + icoo] *
        (local_xyzcen[ind*3+icoo] - mesh_quantities->cell_cen[iel*3+icoo]);

    }

    surface = sqrt(surface);

    distance_fbr_cel /= surface;
    distance_cel_cel /= surface;

    for (icoo = 0 ; icoo < 3 ; icoo++)

      couplage->local_of[ind*3 + icoo] =
        mesh_quantities->b_face_cog[ifac*3 + icoo]
        -  (mesh_quantities->b_face_cog[ifac*3 + icoo] /*  O'  */
            + mesh_quantities->b_face_normal[ifac*3 + icoo] *distance_fbr_cel/surface   /*J'=F+n*FJ'*/
            - 0.5*mesh_quantities->b_face_normal[ifac*3 + icoo] *distance_cel_cel/surface );  /*-n*I'J'/2*/

  }

  reverse = 1;

  fvm_locator_exchange_point_var(couplage->localis_fbr,
                                 couplage->distant_of,
                                 couplage->local_of,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 reverse);

  BFT_FREE(local_xyzcen);


}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of code coupling
 *
 * Fortran interface:
 *
 * SUBROUTINE NBCCPL
 * *****************
 *
 * INTEGER          NBRCPL         : <-- : number of code couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbccpl, NBCCPL)
(
 cs_int_t  *nbrcpl
)
{
  *nbrcpl = cs_glob_nbr_couplages;
}


/*----------------------------------------------------------------------------
 * Set the list of cells and boundary faces associated to a coupling
 * and a cloud of point.
 *
 * The local "support" cells and boundary faces are used to localize
 * the values in the distant "coupled" cells and faces.
 * Depending on the role of sender and/or receiver of the current process
 * in the coupling, some of these sets can be empty or not.
 *
 * The cell values are always localized and interpolated on the distant
 * "cells" support. The face values are localized and interpolated on
 * the distant "face" support if present, or on the distant "cell" support
 * if not.
 *
 * If the input arrays LCESUP and LFBSUP are not ordered, they will be
 * orderd in output.
 *
 * Fortran interface:
 *
 * SUBROUTINE DEFCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCESUP         : --> : number of "support" cells
 * INTEGER          NFBSUP         : --> : number of "support" boundary faces
 * INTEGER          NCECPL         : --> : number of coupled cells
 * INTEGER          NFBCPL         : --> : number of coupled boundary faces
 * INTEGER          LCESUP(NCESUP) : <-> : list of "support" cells
 * INTEGER          LFBSUP(NFBSUP) : <-> : list of "support" boundary faces
 * INTEGER          LCECPL(NCECPL) : --> : list of coupled cells
 * INTEGER          LFBCPL(NFBCPL) : --> : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (defcpl, DEFCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *ncesup,
 const cs_int_t  *nfbsup,
 const cs_int_t  *ncecpl,
 const cs_int_t  *nfbcpl,
       cs_int_t   lcesup[],
       cs_int_t   lfbsup[],
 const cs_int_t   lcecpl[],
 const cs_int_t   lfbcpl[]
)
{
  cs_int_t  ind;

  int  indic_glob[2] = {0, 0};
  int  indic_loc[2] = {0, 0};

  cs_sat_coupling_t  *coupl = NULL;
  fvm_nodal_t  *support_fbr = NULL;
  cs_mesh_quantities_t  *mesh_quantities = cs_glob_mesh_quantities;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  /* Removing the connectivity and localization informations in case of
     coupling update */

  if (coupl->cells_sup != NULL)
    fvm_nodal_destroy(coupl->cells_sup);
  if (coupl->faces_sup != NULL)
    fvm_nodal_destroy(coupl->faces_sup);

  /* Create the local lists */

  coupl->nbr_cel_sup = *ncesup;
  coupl->nbr_fbr_sup = *nfbsup;

  /* Create the corresponding FVM structures */

  if (*ncesup > 0)
    indic_loc[0] = 1;
  if (*nfbsup > 0)
    indic_loc[1] = 1;

  for (ind = 0 ; ind < 2 ; ind++)
    indic_glob[ind] = indic_loc[ind];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allreduce (indic_loc, indic_glob, 2, MPI_INT, MPI_MAX,
                   cs_glob_mpi_comm);
#endif

  if (indic_glob[0] > 0)
    coupl->cells_sup = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                      "coupled_cells",
                                                      *ncesup,
                                                      lcesup);
  if (indic_glob[1] > 0)
    coupl->faces_sup = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                      "coupled_boundary_faces",
                                                      0,
                                                      *nfbsup,
                                                      NULL,
                                                      lfbsup);

  /* Initialization of the distant point localization */

  fvm_locator_set_nodal(coupl->localis_cel,
                        coupl->cells_sup,
                        1,
                        3,
                        *ncecpl,
                        lcecpl,
                        mesh_quantities->cell_cen);

  if (indic_glob[1] > 0)
    support_fbr = coupl->faces_sup;
  else
    support_fbr = coupl->cells_sup;

  fvm_locator_set_nodal(coupl->localis_fbr,
                        support_fbr,
                        1,
                        3,
                        *nfbcpl,
                        lfbcpl,
                        mesh_quantities->b_face_cog);


  /* Computed some quantities needed for a centred-like interpolation */

  if (coupl->localis_fbr != NULL)
    _sat_coupling_interpolate(coupl);


#if 0
  /* TODO: associate the FVM meshes to the post-processing,
     with a fonction giving a pointer to the associated FVM structures,
     and another enabling its compacting or removing */
  {
    fvm_writer_t *w = fvm_writer_init("coupled_mesh",
                                      NULL,
                                      "EnSight Gold",
                                      "binary",
                                      FVM_WRITER_FIXED_MESH);

    fvm_writer_export_nodal(w, coupl->cells_sup);
    fvm_writer_finalize(w);

  }
#endif

  /* Compacting the interpolation support (could be removed) */

  if (coupl->cells_sup != NULL)
    fvm_nodal_reduce(coupl->cells_sup, 1);
  if (coupl->faces_sup != NULL)
    fvm_nodal_reduce(coupl->faces_sup, 1);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  fvm_locator_dump(coupl->localis_cel);
  fvm_locator_dump(coupl->localis_fbr);
#endif

}

/*----------------------------------------------------------------------------
 * Get the number of cells and boundary faces, "support", coupled and not
 * localized associated to a given coupling
 *
 * Fortran interface:
 *
 * SUBROUTINE NBECPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCESUP         : <-- : number of "support" cells
 * INTEGER          NFBSUP         : <-- : number of "support" boundary faces
 * INTEGER          NCECPL         : <-- : number of coupled cells
 * INTEGER          NFBCPL         : <-- : number of coupled boundary faces
 * INTEGER          NCENCP         : <-- : number of not coupled cells
 *                                 :     : (since not localized)
 * INTEGER          NFBNCP         : <-- : number of not coupled boundary faces
 *                                 :     : (since not localized)
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbecpl, NBECPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *ncesup,
       cs_int_t  *nfbsup,
       cs_int_t  *ncecpl,
       cs_int_t  *nfbcpl,
       cs_int_t  *ncencp,
       cs_int_t  *nfbncp
)
{
  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  *ncesup = coupl->nbr_cel_sup;
  *nfbsup = coupl->nbr_fbr_sup;

  *ncecpl = 0;
  *nfbcpl = 0;

  *ncencp = 0;
  *nfbncp = 0;

  if (coupl->localis_cel != NULL) {
    *ncecpl = fvm_locator_get_n_interior(coupl->localis_cel);
    *ncencp = fvm_locator_get_n_exterior(coupl->localis_cel);
  }

  if (coupl->localis_fbr != NULL) {
    *nfbcpl = fvm_locator_get_n_interior(coupl->localis_fbr);
    *nfbncp = fvm_locator_get_n_exterior(coupl->localis_fbr);
  }

}

/*----------------------------------------------------------------------------
 * Get the lists of coupled cells and boundary faces (i.e. receiving)
 * associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LELCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCECPL         : --> : number of coupled cells
 * INTEGER          NFBCPL         : --> : number of coupled boundary faces
 * INTEGER          LCECPL(*)      : <-- : list of coupled cells
 * INTEGER          LFBCPL(*)      : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lelcpl, LELCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *ncecpl,
 const cs_int_t  *nfbcpl,
       cs_int_t  *lcecpl,
       cs_int_t  *lfbcpl
)
{
  cs_int_t  ind;

  cs_int_t  _ncecpl = 0;
  cs_int_t  _nfbcpl = 0;

  cs_sat_coupling_t  *coupl = NULL;

  const cs_int_t  *lst = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (coupl->localis_cel != NULL)
    _ncecpl = fvm_locator_get_n_interior(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    _nfbcpl = fvm_locator_get_n_interior(coupl->localis_fbr);

  if (*ncecpl != _ncecpl || *nfbcpl != _nfbcpl)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for LELCPL()\n"
                "NCECPL = %d and NFBCPL = %d are indicated.\n"
                "The values for this coupling should be %d and %d."),
              *numcpl, (int)(*ncecpl), (int)(*nfbcpl),
              (int)_ncecpl, (int)_nfbcpl);

  /* Copy lists (would be useless with a pure C API) */

  if (_ncecpl > 0) {
    lst = fvm_locator_get_interior_list(coupl->localis_cel);
    for (ind = 0 ; ind < _ncecpl ; ind++)
      lcecpl[ind] = lst[ind];
  }

  if (_nfbcpl > 0) {
    lst = fvm_locator_get_interior_list(coupl->localis_fbr);
    for (ind = 0 ; ind < _nfbcpl ; ind++)
      lfbcpl[ind] = lst[ind];
  }
}

/*----------------------------------------------------------------------------
 * Get the lists of not coupled cells and boundary faces (i.e. receiving but
 * not localized) associated to a given coupling
 *
 * The number of cells and boundary faces, got with NBECPL(), are used
 * for arguments coherency checks.
 *
 * Fortran interface:
 *
 * SUBROUTINE LENCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCENCP         : --> : number of not coupled cells
 * INTEGER          NFBNCP         : --> : number of not coupled boundary faces
 * INTEGER          LCENCP(*)      : <-- : list of not coupled cells
 * INTEGER          LFBNCP(*)      : <-- : list of not coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (lencpl, LENCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *ncencp,
 const cs_int_t  *nfbncp,
       cs_int_t  *lcencp,
       cs_int_t  *lfbncp
)
{
  cs_int_t  ind;

  cs_int_t  _ncencp = 0;
  cs_int_t  _nfbncp = 0;
  cs_sat_coupling_t  *coupl = NULL;

  const cs_int_t  *lst = NULL;


  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (coupl->localis_cel != NULL)
    _ncencp = fvm_locator_get_n_exterior(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    _nfbncp = fvm_locator_get_n_exterior(coupl->localis_fbr);

  if (*ncencp != _ncencp || *nfbncp != _nfbncp)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for LELNCP()\n"
                "NCENCP = %d and NFBNCP = %d are indicated.\n"
                "The values for this coupling should be %d and %d."),
              *numcpl, (int)(*ncencp), (int)(*nfbncp),
              (int)_ncencp, (int)_nfbncp);

  /* Copy lists (would be useless with a pure C API) */

  if (_ncencp > 0) {
    lst = fvm_locator_get_exterior_list(coupl->localis_cel);
    for (ind = 0 ; ind < _ncencp ; ind++)
      lcencp[ind] = lst[ind];
  }

  if (_nfbncp > 0) {
    lst = fvm_locator_get_exterior_list(coupl->localis_fbr);
    for (ind = 0 ; ind < _nfbncp ; ind++)
      lfbncp[ind] = lst[ind];
  }
}

/*----------------------------------------------------------------------------
 * Get the number of distant point associated to a given coupling
 * and localized on the local domain
 *
 * Fortran interface:
 *
 * SUBROUTINE NPDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NCEDIS         : <-- : number of distant cells
 * INTEGER          NFBDIS         : <-- : numbre de distant boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (npdcpl, NPDCPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *ncedis,
       cs_int_t  *nfbdis
)
{
  cs_sat_coupling_t  *coupl = NULL;

  /* Verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  /* Get the number of points */

  *ncedis = 0;
  *nfbdis = 0;

  if (coupl->localis_cel != NULL)
    *ncedis = fvm_locator_get_n_dist_points(coupl->localis_cel);

  if (coupl->localis_fbr != NULL)
    *nfbdis = fvm_locator_get_n_dist_points(coupl->localis_fbr);

}

/*----------------------------------------------------------------------------
 * Get the distant points coordinates associated to a given coupling
 * and a list of points, and the elements number and type (cell or face)
 * "containing" this points.
 *
 * The number of distant points NBRPTS must be equal to one the arguments
 * NCEDIS or NFBDIS given by NPDCPL(), and is given here for coherency checks
 * between the arguments NUMCPL and ITYSUP.
 *
 * Fortran interface:
 *
 * SUBROUTINE COOCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRPTS         : --> : number of distant points
 * INTEGER          ITYDIS         : --> : 1 : access to the points associated
 *                                 :     :     to the distant cells
 *                                 :     : 2 : access to the points associated
 *                                 :     :     to the distant boundary faces
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * INTEGER          LOCPTS(*)      : <-- : "containing" number associated to
 *                                 :     :   each point
 * DOUBLE PRECISION COOPTS(3,*)    : <-- : distant point coordinates
 * DOUBLE PRECISION DJPPTS(3,*)    : <-- : distant vectors to the coupled face
 * DOUBLE PRECISION PNDPTS(*)      : <-- : distant weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (coocpl, COOCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrpts,
 const cs_int_t  *itydis,
       cs_int_t  *ityloc,
       cs_int_t  *locpts,
       cs_real_t *coopts,
       cs_real_t *djppts,
       cs_real_t *dofpts,
       cs_real_t *pndpts
)
{
  cs_int_t  ind, icoo;

  cs_int_t  n_pts_dist = 0;
  cs_sat_coupling_t  *coupl = NULL;
  fvm_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  *ityloc = 0;

  if (*itydis == 1) {
    localis = coupl->localis_cel;
    *ityloc = 1;
  }
  else if (*itydis == 2) {
    localis = coupl->localis_fbr;
    if (coupl->nbr_fbr_sup > 0)
      *ityloc = 2;
    else
      *ityloc = 1;
  }

  if (localis != NULL)
    n_pts_dist = fvm_locator_get_n_dist_points(localis);

  if (*nbrpts != n_pts_dist)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for COOCPL()\n"
                "ITYDIS = %d and NBRPTS = %d are indicated.\n"
                "The value for NBRPTS should be %d."),
              *numcpl, (int)(*itydis), (int)(*nbrpts), (int)n_pts_dist);

  /* Creation the local lists */

  if (localis != NULL) {

    n_pts_dist = fvm_locator_get_n_dist_points(localis);

    if (n_pts_dist > 0) {

      const fvm_lnum_t   *element;
      const fvm_coord_t  *coord;

      element = fvm_locator_get_dist_locations(localis);
      coord   = fvm_locator_get_dist_coords(localis);

      for (ind = 0 ; ind < n_pts_dist ; ind++) {
        locpts[ind] = element[ind];
        for (icoo = 0 ; icoo < 3 ; icoo++)
          coopts[ind*3 + icoo] = coord[ind*3 + icoo];
      }

      if (*itydis == 2)
        for (ind = 0 ; ind < n_pts_dist ; ind++)
          for (icoo = 0 ; icoo < 3 ; icoo++) {
            djppts[ind*3 + icoo] = coupl->distant_dist_fbr[ind*3 + icoo];
            dofpts[ind*3 + icoo] = coupl->distant_of[ind*3 + icoo];
            pndpts[ind] = coupl->distant_pond_fbr[ind];
          }

    }

  }

}

/*----------------------------------------------------------------------------
 * Get the weighting coefficient needed for a centred-like interpolation
 * in the case of a coupling on boundary faces.
 *
 * Fortran interface:
 *
 * SUBROUTINE PNDCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRCPL         : --> : number of distant points
 * INTEGER          ITYLOC         : <-- : 1 : localization on the local cells
 *                                 :     : 2 : localization on the local faces
 * DOUBLE PRECISION PONDCP(*)      : <-- : weighting coefficients
 *----------------------------------------------------------------------------*/

void CS_PROCF (pndcpl, PNDCPL)
(
 const cs_int_t  *const numcpl,
 const cs_int_t  *const nbrpts,
       cs_int_t  *const ityloc,
       cs_real_t *const pondcp,
       cs_real_t *const distof
)
{
  int             icoo;
  cs_int_t        ind;
  cs_int_t        nfbcpl = 0;
  cs_sat_coupling_t  *coupl = NULL;
  fvm_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (*ityloc == 1)
    bft_error(__FILE__, __LINE__, 0,
              _("The centred interpolation scheme is not available\n"
                "when coupling cells"));
  else if (*ityloc == 2)
    localis = coupl->localis_fbr;


  if (localis != NULL)
    nfbcpl = fvm_locator_get_n_interior(localis);

  if (*nbrpts != nfbcpl)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for PNDCPL().\n"
                "ITYLOC = %d and NBRPTS = %d are indicated.\n"
                "NBRPTS should be %d."),
              *numcpl, (int)(*ityloc), (int)(*nbrpts), (int)nfbcpl);

  /* Creation of the local lists */

  if (localis != NULL) {

    if (nfbcpl > 0) {

      for (ind = 0 ; ind < nfbcpl ; ind++) {
        pondcp[ind] = coupl->local_pond_fbr[ind];
        for (icoo = 0 ; icoo < 3 ; icoo++)
          distof[ind*3 + icoo] = coupl->local_of[ind*3 + icoo];
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Exchange a variable associated to a set of point and a coupling.
 *
 * Fortran interface:
 *
 * SUBROUTINE VARCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          ITYVAR         : --> : 1 : variables defined at cells
 *                                 :     : 2 : variables defined at faces
 * DOUBLE PRECISION VARDIS(*)      : --> : distant variable(to send)
 * DOUBLE PRECISION VARLOC(*)      : <-- : local variable (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (varcpl, VARCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
 const cs_int_t  *ityvar,
       cs_real_t *vardis,
       cs_real_t *varloc
)
{
  cs_int_t  n_val_dist_ref = 0;
  cs_int_t  n_val_loc_ref = 0;
  cs_real_t  *val_dist = NULL;
  cs_real_t  *val_loc = NULL;
  cs_sat_coupling_t  *coupl = NULL;
  fvm_locator_t  *localis = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (*ityvar == 1)
    localis = coupl->localis_cel;
  else if (*ityvar == 2)
    localis = coupl->localis_fbr;

  if (localis != NULL) {
    n_val_dist_ref = fvm_locator_get_n_dist_points(localis);
    n_val_loc_ref  = fvm_locator_get_n_interior(localis);
  }

  if (*nbrdis > 0 && *nbrdis != n_val_dist_ref)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for VARCPL()\n"
                "ITYVAR = %d and NBRDIS = %d are indicated.\n"
                "NBRDIS should be 0 or %d."),
              *numcpl, (int)(*ityvar), (int)(*nbrdis), (int)n_val_dist_ref);

  if (*nbrloc > 0 && *nbrloc != n_val_loc_ref)
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling %d: inconsistent arguments for VARCPL()\n"
                "ITYVAR = %d and NBRLOC = %d are indicated.\n"
                "NBRLOC should be 0 or %d."),
              *numcpl, (int)(*ityvar), (int)(*nbrloc), (int)n_val_loc_ref);

  /* Create the local lists */

  if (localis != NULL) {

    if (*nbrdis > 0)
      val_dist = vardis;
    if (*nbrloc > 0)
      val_loc = varloc;

    fvm_locator_exchange_point_var(localis,
                                   val_dist,
                                   val_loc,
                                   NULL,
                                   sizeof(cs_real_t),
                                   1,
                                   0);

  }

}

/*----------------------------------------------------------------------------
 * Array of integers exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * INTEGER          TABDIS(*)      : --> : distant values (to send)
 * INTEGER          TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbicpl, TBICPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
       cs_int_t  *vardis,
       cs_int_t  *varloc
)
{
  cs_int_t  ind;
  cs_int_t  nbr = 0;
  cs_bool_t  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(vardis, *nbrdis, CS_MPI_INT, coupl->dist_root_rank, 0,
                   varloc, *nbrloc, CS_MPI_INT, coupl->dist_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast (varloc, *nbrloc, CS_MPI_INT, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    nbr = CS_MIN(*nbrdis, *nbrloc);

    for (ind = 0; ind < nbr; ind++)
      varloc[ind] = vardis[ind];

  }
}

/*----------------------------------------------------------------------------
 * Array of reals exchange, associated to a given coupling.
 *
 * It is assumed that the arrays have the same size and the same values on
 * each group of processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE TBRCPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          NBRDIS         : --> : number of values to send
 * INTEGER          NBRLOC         : --> : number of values to receive
 * DOUBLE PRECISION TABDIS(*)      : --> : distant values (to send)
 * DOUBLE PRECISION TABLOC(*)      : <-- : local values (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tbrcpl, TBRCPL)
(
 const cs_int_t  *numcpl,
 const cs_int_t  *nbrdis,
 const cs_int_t  *nbrloc,
       cs_real_t *vardis,
       cs_real_t *varloc
)
{
  cs_int_t  ind;
  cs_int_t  nbr = 0;
  cs_bool_t  distant = false;

#if defined(HAVE_MPI)

  MPI_Status  status;

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    /* Exchange between the groups master node */

    if (cs_glob_rank_id < 1)
      MPI_Sendrecv(vardis, *nbrdis, CS_MPI_REAL, coupl->dist_root_rank, 0,
                   varloc, *nbrloc, CS_MPI_REAL, coupl->dist_root_rank, 0,
                   coupl->comm, &status);

    /* Synchronization inside a group */

    if (cs_glob_n_ranks > 1)
      MPI_Bcast(varloc, *nbrloc, CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif /* defined(HAVE_MPI) */

  if (distant == false) {

    nbr = CS_MIN(*nbrdis, *nbrloc);

    for (ind = 0; ind < nbr; ind++)
      varloc[ind] = vardis[ind];

  }
}

/*----------------------------------------------------------------------------
 * Compute the maximum value of an integer variable associated to a coupling.
 *
 * It is assumed that the integer value is the same for each group of
 * processus (local and distant).
 *
 * Fortran interface:
 *
 * SUBROUTINE MXICPL
 * *****************
 *
 * INTEGER          NUMCPL         : --> : coupling number
 * INTEGER          VALDIS         : --> : distant value (to send)
 * INTEGER          VALMAX         : <-- : local maximum (to receive)
 *----------------------------------------------------------------------------*/

void CS_PROCF (mxicpl, MXICPL)
(
 const cs_int_t  *numcpl,
       cs_int_t  *vardis,
       cs_int_t  *varmax
)
{
  cs_bool_t  distant = false;

#if defined(_CS_HAVE_MPI)

  cs_sat_coupling_t  *coupl = NULL;

  /* Initializations and verifications */

  if (*numcpl < 1 || *numcpl > cs_glob_nbr_couplages)
    bft_error(__FILE__, __LINE__, 0,
              _("Impossible coupling number %d; there are %d couplings"),
              *numcpl, cs_glob_nbr_couplages);
  else
    coupl = cs_glob_couplages[*numcpl - 1];

  if (coupl->comm != MPI_COMM_NULL) {

    distant = true;

    MPI_Allreduce(vardis, varmax, 1, CS_MPI_INT, MPI_MAX, coupl->comm);

  }

#endif /* defined(_CS_HAVE_MPI) */

  if (distant == false) {

    *varmax = *vardis;

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add a coupling.
 *
 * Couplings are allowed either with process totally distinct from the
 * application communicator (cs_glob_mpi_comm), or within this same
 * communicator.
 *
 * parameters:
 *   rang_deb <-- root rank of distant process leader in MPI_COMM_WORLD
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_add(cs_int_t   rang_deb)
{
  cs_sat_coupling_t  *couplage = NULL;

  /* Create the associated structure */

  couplage = _sat_coupling_create(rang_deb);

  /* Increase the couplings global array if necessary */

  if (cs_glob_nbr_couplages == cs_glob_nbr_couplages_max) {

    if (cs_glob_nbr_couplages_max == 0)
      cs_glob_nbr_couplages_max = 2;
    else
      cs_glob_nbr_couplages_max *= 2;

    BFT_REALLOC(cs_glob_couplages,
                cs_glob_nbr_couplages_max,
                cs_sat_coupling_t *);

  }

  /* Associate the new coupling to the structure */

  cs_glob_couplages[cs_glob_nbr_couplages] = couplage;

  cs_glob_nbr_couplages += 1;

  return;

}

/*----------------------------------------------------------------------------
 * Destroy all couplings
 *----------------------------------------------------------------------------*/

void
cs_sat_coupling_all_finalize(void)
{
  int  i;

  for (i = 0 ; i < cs_glob_nbr_couplages ; i++)
    _sat_coupling_destroy(cs_glob_couplages[i]);

  BFT_FREE(cs_glob_couplages);

  cs_glob_nbr_couplages = 0;
  cs_glob_nbr_couplages_max = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
