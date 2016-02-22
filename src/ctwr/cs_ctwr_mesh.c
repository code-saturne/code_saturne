/*============================================================================
 * Definitions, Global variables variables, and functions associated with the
 * exchange zones
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"
#include "bft_mem.h"

#include "fvm_nodal_extract.h"
#include "fvm_nodal_extrude.h"
#include "fvm_nodal_project.h"

#include "cs_base.h"
#include "cs_coupling.h"
#include "cs_ctwr.h"
#include "cs_ctwr_air_props.h"
#include "cs_ctwr_halo.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh_connect.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_CT_MPI_TAG    (int)('C'+'S'+'Z'+'E') /* MPI tag for FVM operations */

/*============================================================================
 * Local variables
 *============================================================================*/

static double _epsilon_denom = 1.e-14; /* Minimum denominator */

static cs_lnum_t            cs_ctwr_nmaxvoi  = 50;

#if defined(HAVE_MPI)
MPI_Status status;
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOct
 * *****************
 *
 * INTEGER          n_ct     : <-- : number of exchange area
 *----------------------------------------------------------------------------*/

void CS_PROCF(geoct, GEOCT) (void)
{
  /* construction du maillage eau*/
  cs_ctwr_maille(cs_glob_mesh, cs_glob_mesh_quantities);

  /* chainage des ct*/
  cs_ctwr_stacking();
  /* construction de l'interpolation  AIR -> EAU    */
  cs_ctwr_adeau(cs_glob_mesh, cs_glob_mesh_quantities);
  /* construction de l'interpolation  EAU -> AIR   */
  cs_ctwr_adair();
}

/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*---------------------------------------------------------------------------
 * Solve the equation "matrix.x = b" with Cramer's rule.
 *
 * parameters:
 *   m[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   1 if matrix is singular, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if (CS_ABS(det) < _epsilon_denom)
    return 1;
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Test of coplanarity
 *----------------------------------------------------------------------------*/

static int
_is_coplanar(const cs_real_t  *coord,
             const cs_lnum_t   nvoi[cs_ctwr_nmaxvoi],
             const cs_lnum_t   nbvoi,
             cs_real_t         mat[4][4],
             cs_real_t         vectBase[3][3],
             const cs_lnum_t   decF,
             const cs_real_t   dh)
{
  cs_real_t det,norme1, norme2, min;
  cs_real_t tmpRes[3];
  cs_lnum_t  i,ii,iii,ind,i0,i1,i2;
  cs_lnum_t numi, numii, numiii;

  det = 0.0;
  i0 = 0;
  i1 = 0;
  i2 = 0;
  min = 1000.0;

  for (i = 0; i < 4; i++) {
    ind = 0;
    for (ii = 0; ii < 4; ii++) {
      if (ii != i) {
        vectBase[0][ind] = mat[ii][1];
        vectBase[1][ind] = mat[ii][2];
        vectBase[2][ind] = mat[ii][3];
        ind++;
      }
    }
    cs_math_3_cross_product(vectBase[0],vectBase[1], tmpRes);
    det += pow(-1,i) * mat[i][0] * cs_math_3_dot_product(vectBase[2], tmpRes);
  }

  if (CS_ABS(det) <= (0.000001*dh)) {
    /*2D*/
    for (i=0; i< nbvoi; i++) {
      for (ii=i+1; ii< nbvoi; ii++) {
        for (iii=ii+1; iii< nbvoi; iii++) {
          numi   = nvoi[i]   + decF;
          numii  = nvoi[ii]  + decF;
          numiii = nvoi[iii] + decF;
          vectBase[0][0] = coord[3*numii   ] - coord[3*numi  ];
          vectBase[0][1] = coord[3*numii +1] - coord[3*numi+1];
          vectBase[0][2] = coord[3*numii +2] - coord[3*numi+2];
          vectBase[1][0] = coord[3*numiii  ] - coord[3*numi  ];
          vectBase[1][1] = coord[3*numiii+1] - coord[3*numi+1];
          vectBase[1][2] = coord[3*numiii+2] - coord[3*numi+2];
          cs_math_3_cross_product(vectBase[0],vectBase[1], tmpRes);
          if (   (cs_math_3_norm(tmpRes) > (0.000001*dh))
              && (CS_ABS(cs_math_3_dot_product(vectBase[0],
                                               vectBase[1])) < min)) {
            i0 = i;
            i1 = ii;
            i2 = iii;
          }
        }
      }
    }
  }
  else{
    /*3D*/
    vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
    vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
    vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;
    return 3;
  }

  if (i0 == 0 && i1 == 0 && i2 == 0) {
    /*1D*/
    vectBase[0][0] = coord[3*(nvoi[1]+ decF)  ] - coord[3*(nvoi[0]+ decF)  ];
    vectBase[0][1] = coord[3*(nvoi[1]+ decF)+1] - coord[3*(nvoi[0]+ decF)+1];
    vectBase[0][2] = coord[3*(nvoi[1]+ decF)+2] - coord[3*(nvoi[0]+ decF)+2];
    return 1;
  }
  else{
    vectBase[0][0] = coord[3*(nvoi[i1]+ decF)  ] - coord[3*(nvoi[i0]+ decF)  ];
    vectBase[0][1] = coord[3*(nvoi[i1]+ decF)+1] - coord[3*(nvoi[i0]+ decF)+1];
    vectBase[0][2] = coord[3*(nvoi[i1]+ decF)+2] - coord[3*(nvoi[i0]+ decF)+2];
    vectBase[1][0] = coord[3*(nvoi[i2]+ decF)  ] - coord[3*(nvoi[i0]+ decF)  ];
    vectBase[1][1] = coord[3*(nvoi[i2]+ decF)+1] - coord[3*(nvoi[i0]+ decF)+1];
    vectBase[1][2] = coord[3*(nvoi[i2]+ decF)+2] - coord[3*(nvoi[i0]+ decF)+2];
    cs_math_3_cross_product(vectBase[0],vectBase[1], tmpRes);
    cs_math_3_cross_product(tmpRes, vectBase[0], vectBase[1]);
    norme1= cs_math_3_norm(vectBase[0]);
    norme2= cs_math_3_norm(vectBase[1]);
    vectBase[0][0] /= norme1; vectBase[1][0] /= norme2;
    vectBase[0][1] /= norme1; vectBase[1][1] /= norme2;
    vectBase[0][2] /= norme1; vectBase[1][2] /= norme2;
    vectBase[2][0] = 0.0; vectBase[2][1] = 0.0; vectBase[2][2] = 0.0;
    return 2;
  }

}

/*----------------------------------------------------------------------------
 * Inversion matrice (methode de Jordan)
 *----------------------------------------------------------------------------*/

static int
_invmat(cs_real_t mat[4][4],
        cs_real_t matInv[4][4],
        int       idim)
{
  cs_real_t aux;
  int i, j, k, err;

  err=1;

  for (i = 0; i < idim + 1; i++)
    for (j = 0; j < idim + 1; j++)
      matInv[i][j] = mat[i][j];


  i = 0;
  while (err == 1 && (i < (idim+1))) {

    if (CS_ABS(matInv[i][i]) > 1.e-15) {

      aux = 1.0/matInv[i][i];

      for (j = 0; j < idim + 1; j++)
         matInv[i][j] *=  aux;

      matInv[i][i]= aux;

      for (k = 0; k < i; k++) {
        aux = matInv[k][i];
        for (j = 0; j < idim + 1; j++) {
           matInv[k][j] -= aux * matInv[i][j];
        }
        matInv[k][i] = -aux *matInv[i][i];
      }

      for (k = i + 1; k < idim + 1; k++) {
        aux = matInv[k][i];
        for (j = 0; j < idim + 1; j++) {
           matInv[k][j] -= aux * matInv[i][j];
        }
        matInv[k][i] = -aux *matInv[i][i];
      }
      i++;
    }
    else{
      err =0;
    }
  }

  return err;
}

/*----------------------------------------------------------------------------
 * Calculation of the ponderation function
 *----------------------------------------------------------------------------*/

static cs_real_t
_weighting(const cs_real_t  dx,
           const cs_real_t  dy,
           const cs_real_t  dz,
           const cs_real_t  ouv,
           const cs_lnum_t  lf,
           const cs_real_t  epgauss,
           const cs_real_t  cx,
           const cs_real_t  cy,
           const cs_real_t  cz)
{
  cs_real_t pi, lambda;
  cs_real_t poids = 0.0;
  cs_real_t xy = sqrt(pow(dx/cx,2.) + pow(dy/cy,2.) + pow(dz/cz,2.));

  if (xy < ouv)
  {
    switch(lf) {
    case 1:
      poids=1.-xy/ouv;
      break;
    case 2:
      pi = acos(-1.);
      poids = 0.5*(1.+cos(pi*xy/ouv));
      break;
    case 3:
      lambda = ouv/sqrt(epgauss*log(10.));
      poids = exp(-pow((xy/lambda),2.));
      break;
    default:
      assert(lf == 1 || lf == 2 || lf == 3);
    }
  }

  return poids;
}

/*----------------------------------------------------------------------------
 * Scalar product between the face normal vector and the gravity vector
 *----------------------------------------------------------------------------*/

static cs_real_t
_dot_product_ng(const cs_lnum_t  ifac,
                const int        dim,
                const cs_real_t *surf_f,
                const cs_real_t  gravite[3],
                const cs_real_t  direction)
{
  cs_real_t n_sortant[3], g_uni[3];
  int idim;
  cs_real_t aux1,aux2;

  aux1 = cs_math_3_norm(gravite);

  for (idim = 0; idim < 3; idim++) {
    n_sortant[idim] = direction * surf_f[ifac*3+idim];
    g_uni[idim]= gravite[idim];
  }

  aux2 = cs_math_3_norm(n_sortant);
  for (idim = 0; idim < dim; idim++) {
    n_sortant[idim] /= aux2;
    g_uni[idim] /= aux1;

  }

  return cs_math_3_dot_product(n_sortant, g_uni);
}


/*---------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/

static void
_search_height(cs_ctwr_zone_t   *ct,
               const cs_real_t  *gravite,
               cs_real_t        *hmin,
               cs_real_t        *hmax)
{
  cs_lnum_t    i, ifac, nb,nb_dist, axe, idx, ii, jj;
  cs_real_t   *lst_xyz_sup; /* coord des sommets de la face Sup */
  cs_real_t   *lst_xyz_inf; /* coord des sommets de la face Inf*/
  cs_real_t   *lst_xyz_fi;  /* coord des sommets proj de la face Inf*/
  cs_real_t   *lst_xyz_fs;  /* coord des sommets proj de la face Inf*/
  cs_real_t   *hmin_dist;
  const cs_coord_t *lst_xyz_dist = NULL;

  cs_real_t   aux;
  const cs_lnum_t  *location_fac = NULL;

  cs_lnum_t  *faces_vtx_idx   = NULL;
  cs_lnum_t  *faces_vtx_lst   = NULL;

  const double tolerance = 0.1;
  fvm_nodal_t *fs_tmp_mesh = NULL;
  fvm_nodal_t *fi_tmp_mesh = NULL;

  double coeff[3], v_aux[3], vertex_coords[2];
  double  v_x, v_y;
  double  v_f_x = 0., v_f_y = 0., v_f_z = 0.;
  double a[3][3] = {{0., 0., 0.},
                    {0., 0., 0.},
                    {0., 0., 0.} };

  double b_x[3] = {0., 0., 0. };
  double b_y[3] = {0., 0., 0. };
  double b_z[3] = {0., 0., 0. };

  double matrice[6] = {0., 0., 0., 0., 0., 0. };

  ple_locator_t   *locator = NULL;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 1;

  nb = (cs_lnum_t) fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);
  BFT_MALLOC(lst_xyz_sup, nb*3, cs_coord_t);
  fvm_nodal_get_vertex_coords(ct->face_sup_mesh,
                              CS_INTERLACE,
                              lst_xyz_sup);

  nb = (cs_lnum_t) fvm_nodal_get_n_entities(ct->face_inf_mesh, 0);
  BFT_MALLOC(lst_xyz_inf, nb*3, cs_coord_t);
  fvm_nodal_get_vertex_coords(ct->face_inf_mesh,
                              CS_INTERLACE,
                              lst_xyz_inf);

  fs_tmp_mesh = fvm_nodal_copy(ct->face_sup_mesh);
  fi_tmp_mesh = fvm_nodal_copy(ct->face_inf_mesh);

  aux = 0.;

  for (i = 0; i < 3; i++)
    if (CS_ABS(gravite[i]) > aux) {
      axe = i;
      aux = CS_ABS (gravite [i]);
    }

  if (axe == 0) {
    matrice[1] = 1;
    matrice[5] = 1;
  }
  else {
    matrice[0] = 1;
    if (axe == 1)
      matrice[5] = 1;
    else
      matrice[4] = 1;
  }

  fvm_nodal_project_coords(fs_tmp_mesh, matrice);
  fvm_nodal_project_coords(fi_tmp_mesh, matrice);

  nb = (cs_lnum_t) fvm_nodal_get_n_entities(fs_tmp_mesh, 0);

  BFT_MALLOC(lst_xyz_fs, nb*2, cs_coord_t);

  fvm_nodal_get_vertex_coords(fs_tmp_mesh, CS_INTERLACE, lst_xyz_fs);


  nb = (cs_lnum_t) fvm_nodal_get_n_entities(fi_tmp_mesh, 0);

  BFT_MALLOC(lst_xyz_fi, nb*2, cs_coord_t);

  fvm_nodal_get_vertex_coords(fi_tmp_mesh, CS_INTERLACE, lst_xyz_fi);

  /* Create locator on the proj surf  */

#if defined(PLE_HAVE_MPI)
  locator = ple_locator_create(cs_glob_mpi_comm,
                               cs_glob_n_ranks,
                               0);
#else
  locator = ple_locator_create();
#endif

  nb = (cs_lnum_t) fvm_nodal_get_n_entities(fs_tmp_mesh, 0);

  ple_locator_set_mesh(locator,
                       fi_tmp_mesh,
                       locator_options,
                       0,
                       tolerance,
                       2,
                       nb,
                       NULL,
                       NULL,
                       lst_xyz_fs,
                       NULL,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh);

  nb_dist = ple_locator_get_n_dist_points(locator);

  /* Construction de la connectivite Face->sommet du projete  */

  BFT_MALLOC(hmin_dist, nb_dist, cs_coord_t);

  cs_reverse_vtx_faces_connect(fi_tmp_mesh,
                                &(faces_vtx_idx),
                                &(faces_vtx_lst));

  location_fac = ple_locator_get_dist_locations(locator);
  lst_xyz_dist = ple_locator_get_dist_coords(locator);

  for (i = 0; i < nb_dist; i++) {

    ifac = location_fac[i] - 1;

    vertex_coords [0] = lst_xyz_dist [i*2    ];
    vertex_coords [1] = lst_xyz_dist [i*2 + 1];

    for (ii = 0; ii < 3; ii++) {
      b_x[ii] = 0.;
      b_y[ii] = 0.;
      b_z[ii] = 0.;
      for (jj = 0; jj < 3; jj++)
        a[ii][jj] = 0.;

    }


    for (idx = faces_vtx_idx[ifac   ];
         idx < faces_vtx_idx[ifac +1]; idx++) {

      v_x = lst_xyz_fi[faces_vtx_lst[idx]* 2    ];
      v_y = lst_xyz_fi[faces_vtx_lst[idx]* 2 + 1];

      v_f_x = lst_xyz_inf[faces_vtx_lst[idx]* 3    ];
      v_f_y = lst_xyz_inf[faces_vtx_lst[idx]* 3 + 1];
      v_f_z = lst_xyz_inf[faces_vtx_lst[idx]* 3 + 2];

      a[0][0] += v_x * v_x;
      a[0][1] += v_x * v_y;
      a[0][2] += v_x;

      a[1][1] += v_y * v_y;
      a[1][2] += v_y;


      a[2][2] += 1.;

      b_x[0] += v_x * v_f_x;
      b_x[1] += v_y * v_f_x;
      b_x[2] += v_f_x;

      b_y[0] += v_x * v_f_y;
      b_y[1] += v_y * v_f_y;
      b_y[2] += v_f_y;

      b_z[0] += v_x * v_f_z;
      b_z[1] += v_y * v_f_z;
      b_z[2] += v_f_z;


    }

    /* Matrix is symmetric */

    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];

    if (_inverse_3x3(a, b_x, coeff) == 0) {

      v_aux[0] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
       v_aux[0] = -v_f_x;

    if (_inverse_3x3(a, b_y, coeff) == 0) {

      v_aux[1] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
      v_aux[1] = -v_f_y;

    if (_inverse_3x3(a, b_z, coeff) == 0) {

      v_aux[2] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
      v_aux[2] = -v_f_z;

    hmin_dist[i] =  cs_math_3_dot_product(v_aux, gravite)
                  / cs_math_3_norm(gravite);
  }


  ple_locator_exchange_point_var(locator,
                                 hmin_dist, hmin, NULL, sizeof(cs_real_t),1,0);

  for (i = 0; i < nb; i++) {

      v_aux[0] = -lst_xyz_sup[i*3    ];/* Opposite Vector to g */
      v_aux[1] = -lst_xyz_sup[i*3 + 1];
      v_aux[2] = -lst_xyz_sup[i*3 + 2];

      aux = cs_math_3_dot_product(v_aux, gravite);
      hmax[i] = aux / cs_math_3_norm(gravite); /* project on "g" axis */

  }

   BFT_FREE(lst_xyz_inf);
   BFT_FREE(lst_xyz_sup);
   BFT_FREE(lst_xyz_fi);
   BFT_FREE(lst_xyz_fs);
   BFT_FREE(hmin_dist);

   locator = ple_locator_destroy(locator);
   fs_tmp_mesh = fvm_nodal_destroy(fs_tmp_mesh);
   fi_tmp_mesh = fvm_nodal_destroy(fi_tmp_mesh);
}


/*----------------------------------------------------------------------------
 * Function cs_ctwr_maille
 * Construction du maillage eau
 *----------------------------------------------------------------------------*/

void cs_ctwr_maille(const cs_mesh_t             *mesh,
                    const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t   icel_1, icel_2, ii, length, nb, rank,
             dist_rank, res_loc, res_dist;
  cs_lnum_t   ifac, ict, icpt, icpti, icptla, icptfac,
             iaux, i, j;
  cs_real_t  aux, gravite[3], v_aux[3], alpha;
  cs_coord_t *extrusion_vectors, *lst_xyz_cel, *lst_xyz;
  cs_lnum_t   *lst_par_fac_sup;
  cs_gnum_t   *fsup_gb_vt_num = NULL;
  cs_real_t   *hmin_vect;
  cs_real_t   *hmax_vect;

  char  *mesh_name       = NULL;
  char  *export_name     = NULL;
  const double tolerance = 0.1;

  cs_lnum_t   n_vertices;

  cs_lnum_t   *face_sup;      /* liste des faces internes superieures de la ct
                                 de taille  (nbfac_sct) */
  cs_lnum_t   *fbr_sup;       /* liste des faces de bord superieures de la ct
                                 de taille  (nbfbr_sct) */
  cs_lnum_t   *face_inf;      /* liste des faces internes inferieures de la ct
                                 de taille  (nbfac_ict) */
  cs_lnum_t   *fbr_inf;       /* liste des faces de bord  inferieures de la ct
                                 de taille  (nbfac_ict) */
  cs_lnum_t   *face_lat;      /* liste des faces internes laterales de la ct
                                 de taille  (nbfac_lct) */
  cs_lnum_t   *fbr_lat;       /* liste des faces de bord laterales de la ct
                                 de taille  (nbfbr_lct) */
  cs_lnum_t   *face_ct;       /* liste des faces interne de la ct
                                 de taille  (nbfac_ct) */

  const cs_lnum_2_t  *i_face_cells  = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const cs_lnum_t *b_face_cells  = mesh->b_face_cells;
  const cs_real_t *i_face_normal = mesh_quantities->i_face_normal;
  const cs_real_t *b_face_normal = mesh_quantities->b_face_normal;

  cs_interface_set_t  *interface_set = NULL;
  cs_ctwr_zone_t  *ct;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 1;

  iaux = 0;
  alpha = 0.875;

  /* Vecteur gravite */
  gravite[0] = ct_prop->gravx;
  gravite[1] = ct_prop->gravy;
  gravite[2] = ct_prop->gravz;

  /*--------------------------------------------*/
  /* List of air nodes for each Exchange Area   */
  /*--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    icpt = 0;
    ct = cs_glob_ct_tab[ict];
    length = strlen("cell_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "cell_mesh_ct_%d", ict);

    ct->cell_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                                   mesh_name,
                                                   false,
                                                   ct->nbevct,
                                                   ct->ze_cell_list);

    BFT_MALLOC(ct->mark_ze,mesh->n_cells, cs_lnum_t);

    /*----------------------------------------------------------*
     * Begin identification of air nodes for each Exchange Area *
     *----------------------------------------------------------*/

    for (i = 0; i < mesh->n_cells; i++)
          for (j = 0; j < ct->nbevct; j++) {
                  if ((ct->ze_cell_list[j]) == i+1) {
                  ct->mark_ze[i]=1;
                  break;
        }
        else
                  ct->mark_ze[i]=0;
        }
     }

  /*---------------------------------------------*
   * End list of air nodes for each Exchange Area*
   *---------------------------------------------*/

  /*--------------------------------------------------------*
   * Calcul du nombre de noeuds eau des faces superieures   *
   * des zones d'echanges et du nombre de faces superieures *
   * et inferieures                                         *
   *--------------------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {
    ct = cs_glob_ct_tab[ict];
    /* Contribution faces internes */
    for (ifac = 0; ifac < mesh->n_i_faces; ifac++) {
      assert((ifac * 2 + 1) < (2*mesh->n_i_faces));
      icel_1 = i_face_cells[ifac][0];
      icel_2 = i_face_cells[ifac][1];
      /* Comparaison  des couleurs des cellules 1 et 2 */
      if ((ct->mark_ze[icel_1] == 1) ||
          (ct->mark_ze[icel_2] == 1)) {
        if  (ct->mark_ze[icel_1] != ct->mark_ze[icel_2]) {
          if (ct->mark_ze[icel_1] == 1) {
            aux = _dot_product_ng(ifac, ct->idimct, i_face_normal, gravite, 1);
          }
          if (ct->mark_ze[icel_2] == 1) {
            aux = _dot_product_ng(ifac, ct->idimct, i_face_normal, gravite, -1);
          }

          if (aux < (-alpha)) {
            ct->nnpsct++;
            ct->nbfac_sct++;
          }else{
            if (aux > alpha) {
              ct->nbfac_ict++;
            }else{
              ct->nbfac_lct++;
            }
          }
        }else{
          ct->nbfac_ct++;
        }
      }

    }  /* fin contribution faces internes */

    /* Contribution faces externes */
    for (ifac = 0; ifac < mesh->n_b_faces; ifac++) {
      icel_1 = b_face_cells[ifac]; /* indice de la cellule  */
      if (ct->mark_ze[icel_1] == 1) {

        aux = _dot_product_ng(ifac,ct->idimct, b_face_normal, gravite, 1);

        if (aux < (-alpha)) {
          ct->nnpsct++;
          ct->nbfbr_sct++;
        }else{
          if (aux > alpha) {
            ct->nbfbr_ict++;
          }else{
            ct->nbfbr_lct++;
          }
        }
      }
    }/* fin contribution faces externes */


    /* allocation memoire pour la liste des faces superieures et inferieures
    * des ct */
    BFT_MALLOC(face_sup,ct->nbfac_sct, cs_lnum_t);
    BFT_MALLOC(face_inf,ct->nbfac_ict, cs_lnum_t);
    BFT_MALLOC(face_lat,ct->nbfac_lct, cs_lnum_t);
    BFT_MALLOC(fbr_sup, ct->nbfbr_sct, cs_lnum_t);
    BFT_MALLOC(fbr_inf, ct->nbfbr_ict, cs_lnum_t);
    BFT_MALLOC(fbr_lat, ct->nbfbr_lct, cs_lnum_t);
    BFT_MALLOC(face_ct, ct->nbfac_ct, cs_lnum_t);


  /* --------------------------------------------------------*
   * Fin Calcul du nombre de noeuds eau des faces superieures*
   * des zones d'echanges et du nombre de faces superieures  *
   * et inferieures                                          *
   *---------------------------------------------------------*/


  /*-----------------------------------------------------------------*
   * Liste des faces superieures et inferieures des zones d echanges *
   * et liste des noeuds eau des faces sup de la ct sans ct amont    *
   *-----------------------------------------------------------------*/

    /* Contribution faces internes */
    icpt   = 0; /*indice tableau des faces  sup */
    icpti  = 0; /*indice tableau des faces  inf */
    icptla = 0; /*indice tableau des faces  laterales */
    icptfac  = 0; /*indice tableau des noeuds sup ct */
    /* Boucle sur les faces internes du domaine */
    for (ifac = 0; ifac < mesh->n_i_faces; ifac++) {
      icel_1 = i_face_cells[ifac][0];
      icel_2 = i_face_cells[ifac][1];
      /* Comparaison  couleur de la ct et couleur des cellules 1 et 2 */
      if ((ct->mark_ze[icel_1] == 1) ||
          (ct->mark_ze[icel_2] ==1)) {
        if  (ct->mark_ze[icel_1] != ct->mark_ze[icel_2]) {
          if (ct->mark_ze[icel_1] ==1) {
            aux = _dot_product_ng(ifac, ct->idimct, i_face_normal, gravite, 1);
          }
          if (ct->mark_ze[icel_2] == 1) {
            aux = _dot_product_ng(ifac, ct->idimct, i_face_normal, gravite, -1);
          }

          if (aux < (-alpha)) {
            /*ajout d'une face sup de la ct*/
            face_sup[icpt] = ifac + 1;
            icpt ++;
          }else{
            if (aux > alpha) {
            /*ajout d'un face inf de la ct*/
            face_inf[icpti]  = ifac + 1;
            icpti ++;

            }else{
              assert(icptla < ct->nbfac_lct+ct->nbfbr_lct);
              face_lat[icptla] = ifac  + 1;
              icptla ++;
            }
          }
        }else{
          face_ct[icptfac] = ifac  + 1;
          icptfac ++;
        }
      }
    }/* fin contribution faces internes */

    /* Contribution faces de bords */
    /* initialisation des indices */
    icpt   = 0; /*indice tableau des faces  sup */
    icpti  = 0; /*indice tableau des faces  inf */
    icptla = 0; /*indice tableau des faces  laterales */

    for (ifac = 0; ifac < mesh->n_b_faces; ifac++) {

      icel_1 = b_face_cells[ifac];/* indice de la cellule  */
      if (ct->mark_ze[icel_1]== 1) {

        aux = _dot_product_ng(ifac, ct->idimct, b_face_normal, gravite, 1);

        if (aux < (-alpha)) {
          /* ajout d'une face sup de la ct */
          fbr_sup[icpt]= ifac + 1;
          icpt ++;
        }else{
          if (aux > alpha) {
            /*ajout d'un face inf de la ct*/
            fbr_inf[icpti]= ifac + 1;
            icpti ++;
          }else{
            fbr_lat[icptla]= ifac + 1;
            icptla ++;
          }
        }
      }
    } /* fin contribution faces externes */

    /*---------------------------------------------------------*
    * Creation des maillages surfacique en connectivite nodale*
    *---------------------------------------------------------*/

    /* mesh for superiors faces */

    BFT_FREE(mesh_name);
    length = strlen("face_sup_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_sup_mesh_ct_%d", ict);


    ct->face_sup_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       false,
                                                       ct->nbfac_sct,
                                                       ct->nbfbr_sct,
                                                       face_sup,
                                                       fbr_sup);

    /* mesh for inferiors faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_inf_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_inf_mesh_ct_%d", ict);


    ct->face_inf_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       false,
                                                       ct->nbfac_ict,
                                                       ct->nbfbr_ict,
                                                       face_inf,
                                                       fbr_inf);
    /* mesh for laterals faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_lat_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_lat_mesh_ct_%d", ict);


    ct->face_lat_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       false,
                                                       ct->nbfac_lct,
                                                       ct->nbfbr_lct,
                                                       face_lat,
                                                       fbr_lat);
    /* mesh for laterals faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_mesh_ct_%d", ict);


    ct->fac_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                  mesh_name,
                                                  false,
                                                  ct->nbfac_ct,
                                                  0,
                                                  face_ct,
                                                  NULL);
    /* water mesh*/
    BFT_FREE(mesh_name);
    length = strlen("water_mesh_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "water_mesh_%d", ict);

    ct->water_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                    mesh_name,
                                                    false,
                                                    ct->nbfac_sct,
                                                    ct->nbfbr_sct,
                                                    face_sup,
                                                    fbr_sup);

    /*--------------------------------------------------------------*
    *  Fin creation des maillages surfacique en connectivite nodale*
    *--------------------------------------------------------------*/

    /*--------------------------------------------------------------*
    * Construct cs_array_rank                                  *
    *--------------------------------------------------------------*/

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {

      nb   = cs_glob_n_ranks;
      rank = cs_glob_rank_id;
      BFT_MALLOC(ct->cs_array_rank, nb, cs_lnum_t);


      ct->cs_array_rank[rank] = res_loc = ct->nbevct;


      for (dist_rank = 0; dist_rank <  nb; dist_rank++)
        if (dist_rank != rank) {
          MPI_Sendrecv(&res_loc, 1, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                        &res_dist, 1, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                        cs_glob_mpi_comm, &status);

          ct->cs_array_rank[dist_rank] = res_dist;

        }
    }
#endif

    /*--------------------------------------------------------------*
    * End of Construction cs_array_rank                         *
    *--------------------------------------------------------------*/

    /*--------------------------------------------------------------*
    * Reseach of  hmin and hmax                                    *
    *--------------------------------------------------------------*/

    /* loop on the superior faces for hmax */
    nb = (cs_lnum_t) fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);

    BFT_MALLOC(hmax_vect, nb, cs_coord_t);
    BFT_MALLOC(hmin_vect, nb, cs_coord_t);

    _search_height(ct,
                   gravite,
                   hmin_vect,
                   hmax_vect);

    for (i = 0; i < nb; i++) {

      aux = hmax_vect[i];
      if (aux >= ct->hmax)
        ct->hmax = aux;

      aux = hmin_vect[i];
      if (aux <= ct->hmin)
        ct->hmin = aux;
    }

    /* loop on the sup faces for surface_in and surface_out */
    BFT_MALLOC(lst_par_fac_sup, ct->nnpsct, cs_lnum_t);

    fvm_nodal_get_parent_num(ct->face_sup_mesh, 2, lst_par_fac_sup);

    BFT_MALLOC(ct->surf_fac_sup, ct->nnpsct, cs_real_t);

    for (ifac = 0; ifac < ct->nnpsct; ifac++) {
      if (ifac< ct->nbfbr_sct) {
        for (ii = 0; ii < 3; ii++)
          v_aux[ii] = b_face_normal[3 * (lst_par_fac_sup[ifac] -1) + ii];
      }
      else{
        for (ii = 0; ii < 3; ii++)
          v_aux[ii] = i_face_normal[3 * (  lst_par_fac_sup[ifac]
                                         - mesh->n_b_faces - 1) + ii];
      }
      aux = cs_math_3_norm(v_aux);
      ct->surface_in += aux;
      ct->surface_out += aux;
      ct->surf_fac_sup[ifac] = aux;

    }



#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {

      nb   = cs_glob_n_ranks;
      rank = cs_glob_rank_id;

      /* TODO : changer ce bordel !!!!!!!!!!!!!!!!!
         sans doute equivalent a MPI_Allreduce(ct-hmax, ..., MPI_MAX) */

      if (ct->cs_array_rank[rank] != 0) {
        for (dist_rank = 0; dist_rank < nb; dist_rank++) {
          if (dist_rank != rank) {
            if (ct->cs_array_rank [dist_rank] != 0) {

              MPI_Sendrecv(&ct->hmax, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           &aux, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           cs_glob_mpi_comm, &status);

              if (aux > ct->hmax) ct->hmax = aux;

              MPI_Sendrecv(&ct->hmin, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           &aux, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           cs_glob_mpi_comm, &status);

              if (aux < ct->hmin) ct->hmin = aux;

            }
          }
        }
      }

      MPI_Allreduce (&ct->surface_in, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->surface_in = aux;

      MPI_Allreduce (&ct->surface_out, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->surface_out = aux;
    }
#endif

   /* -------------------------------------------------------------*
    * End of Reseach  hmin et hmax                                  *
    *---------------------------------------------------------------*/

    nb = fvm_nodal_get_n_entities(ct->water_mesh, 0);

    BFT_MALLOC(extrusion_vectors, (nb*3), cs_coord_t);

    for (i=0; i < nb; i++) {

      aux =CS_ABS(hmax_vect[i] -  hmin_vect[i])/cs_math_3_norm(gravite);

      extrusion_vectors[i*3]     =  gravite[0] * aux;
      extrusion_vectors[i*3 + 1] =  gravite[1] * aux;
      extrusion_vectors[i*3 + 2] =  gravite[2] * aux;
    }

    fvm_nodal_extrude(ct->water_mesh,
                      ct->nelect,
                      extrusion_vectors,
                      NULL);

    BFT_FREE(extrusion_vectors);


    /* Set halo structure for the water mesh */

    n_vertices = fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);

    BFT_MALLOC(fsup_gb_vt_num, n_vertices, cs_gnum_t);

    fvm_nodal_get_global_vertex_num(ct->face_sup_mesh, fsup_gb_vt_num);

    interface_set = cs_interface_set_create(n_vertices, NULL, fsup_gb_vt_num,
                                            NULL, 0, NULL, NULL, NULL);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    cs_interface_set_dump(interface_set);
#endif

    /* Creation of the cs_mesh_halo_t structure. */

    bft_printf(_(" Creating halos\n"));
    bft_printf_flush();

    ct->water_halo = cs_halo_create(interface_set);

    bft_printf(_(" Defining halos\n"));
    bft_printf_flush();

    cs_ctwr_halo_define(ct, interface_set);

    cs_interface_set_destroy(&interface_set);

    /* Create locator for interpolate */

#if defined(PLE_HAVE_MPI)
    ct->locat_water_air = ple_locator_create(cs_glob_mpi_comm,
                                             cs_glob_n_ranks,
                                             0);
#else
    ct->locat_water_air = ple_locator_create();
#endif

    BFT_MALLOC(lst_xyz_cel, ct->nbevct*3, cs_coord_t);

    fvm_nodal_get_element_centers(ct->cell_mesh, CS_INTERLACE, 3, lst_xyz_cel);

    ple_locator_set_mesh(ct->locat_water_air,
                         ct->water_mesh,
                         locator_options,
                         0,
                         tolerance,
                         3,
                         ct->nbevct,
                         NULL,
                         NULL,
                         lst_xyz_cel,
                         NULL,
                         cs_coupling_mesh_extents,
                         cs_coupling_point_in_mesh);


#if defined(PLE_HAVE_MPI)
    ct->locat_air_water = ple_locator_create(cs_glob_mpi_comm,
                                             cs_glob_n_ranks,
                                             0);
#else
    ct->locat_air_water = ple_locator_create();
#endif

    BFT_MALLOC(lst_xyz, ct->nnpsct*ct->nelect*3, cs_coord_t);

    fvm_nodal_get_element_centers(ct->water_mesh, CS_INTERLACE, 3, lst_xyz);

    ple_locator_set_mesh(ct->locat_air_water,
                         ct->cell_mesh,
                         locator_options,
                         3,
                         0,
                         tolerance,
                         ct->nnpsct*ct->nelect,
                         NULL,
                         NULL,
                         lst_xyz,
                         NULL,
                         cs_coupling_mesh_extents,
                         cs_coupling_point_in_mesh_p);

    BFT_FREE(mesh_name);
    BFT_FREE(export_name);
    BFT_FREE(face_sup);
    BFT_FREE(face_inf);
    BFT_FREE(face_lat);
    BFT_FREE(fbr_sup);
    BFT_FREE(lst_par_fac_sup);
    BFT_FREE(fbr_inf);
    BFT_FREE(fbr_lat);
    BFT_FREE(face_ct);
    BFT_FREE(lst_xyz);
    BFT_FREE(lst_xyz_cel);
    BFT_FREE(fsup_gb_vt_num);

  }
  /*--------------------------------------------*
   * Fin Liste des faces superieures et inferieures des zones d echanges
   * et liste des noeuds eau des faces sup de la ct*
   *--------------------------------------------*/

  /*--------------------------------------------*
   * Initialization of the water variables      *
   *--------------------------------------------*/

  for (ict=0; ict < cs_glob_ct_nbr; ict++) {
    ct = cs_glob_ct_tab[ict];
    /* Te */
    BFT_MALLOC(ct->teau, (ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);
    /* Fe */
    BFT_MALLOC(ct->fem, (ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);
    /* vg */
    BFT_MALLOC(ct->vgoutte,(ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);

    /* initialisation*/
    for (iaux = 0; iaux < (ct->nnpsct_with_ghosts*ct->nelect); iaux++) {
      /* temperature de l eau*/
      ct->teau[iaux]    = ct->cl_teau;
      /* debit massique par unite de surface */
      ct->fem[iaux]     = ct->cl_fem/ct->surface;
      /* vitesse des gouttes */
      ct->vgoutte[iaux] = 0.0;
    }
    /* Initialisation en tenant compte de l'ecart de temperature impose*/
    aux = ct->deltat / ((cs_real_t) (ct->nelect - 1));
    for (i = 0; i < ct->nnpsct_with_ghosts; i++) {
      for (j = 1; j < ct->nelect; j++) {
          ii = i*ct->nelect + j;
          ct->teau[ii] =  ct->teau[ii - 1] - aux;
        }
    }


  }/*  fin de la boucle sur les zones d'echanges */


  /*----------------------------------------------*
   * End of fnitialization of the water variables *
   *----------------------------------------------*/


  /*--------------------------------------------*
   * Initialisation des tableaux d interpolation*
   *--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[ict];
    /* Liste des voisins eau des cellules air */
    nb = (int) ple_locator_get_n_dist_points(ct->locat_air_water);
    BFT_MALLOC(ct->voiseau,(nb * cs_ctwr_nmaxvoi), cs_lnum_t);
    /* Coefficients d interpolation eau pour l air*/
    BFT_MALLOC(ct->coefeau, (nb * cs_ctwr_nmaxvoi), cs_real_t);
    /* Positions dans la liste des voisins eau */
    BFT_MALLOC(ct->pvoiseau, (nb + 1), cs_lnum_t);
    /* Liste des voisins air des noeuds eau */

    ct->pvoiseau[0] = 0;

    for (iaux = 0; iaux < (nb *cs_ctwr_nmaxvoi); iaux++) {
      ct->voiseau[iaux] = -1;
    }

    nb = (int) ple_locator_get_n_dist_points(ct->locat_water_air);

    BFT_MALLOC(ct->voisair,(nb * cs_ctwr_nmaxvoi), cs_lnum_t);
    /* Positions dans la liste voisins air */
    BFT_MALLOC(ct->pvoisair, (nb + 1), cs_lnum_t);
    /* Coefficients d interpolation air pour l eau */
    BFT_MALLOC(ct->coefair,(nb * cs_ctwr_nmaxvoi), cs_real_t);

    ct->pvoisair[0] = 0;

    for (iaux = 0; iaux < (nb* cs_ctwr_nmaxvoi); iaux++) {
      ct->voisair[iaux] = -1;
    }
  }
  /*------------------------------------------------------*
   * Fin de l initialisation des tableaux d interpolation *
   *------------------------------------------------------*/
}

/*----------------------------------------------------------------------------
 * Function cs_ctwr_adeau
 * Interpolation AIR -> EAU
 *----------------------------------------------------------------------------*/

void cs_ctwr_adeau
(
  const cs_mesh_t             *mesh,
  const cs_mesh_quantities_t  *mesh_quantities
)
{
  /* Coordonnees des centres des cellules  */
  const cs_real_t *coo_cel        = mesh_quantities->cell_cen;
  const cs_lnum_2_t  *i_face_cells = (const cs_lnum_2_t  *)(mesh->i_face_cells);

  cs_lnum_t   ict, iwat,nb_node_water, ii, jj, iair, nbvois,
             nbn, ifac, icel_1, icel_2, lf, indice, dim;
  cs_lnum_t   nvois[cs_ctwr_nmaxvoi];
  cs_real_t  dhi;
  cs_real_t  xwat, ywat, zwat, dx, dy, dz, dxx, dyy, dzz, ouv, aux;
  cs_real_t  coeff[cs_ctwr_nmaxvoi];
  cs_real_t  vectBase[3][3];
  cs_real_t  cx, cy, cz, epgauss, w,
             pp[4][4],ppInv[4][4];

  const cs_coord_t *lst_xyz_water = NULL;
  cs_coord_t *lst_xyz_cel ;
  cs_lnum_t   *lst_par_fac ;
  cs_lnum_t   *lst_par_cel ;

  const cs_lnum_t  *location_cel = NULL;

  /*--------------------------------------------*
   * parametres et initialisation               *
   *--------------------------------------------*/
  cs_ctwr_zone_t  *ct;

  /*--------------------------------------------*
   * fin parametres et initialisation           *
   *--------------------------------------------*/

  /* Make sure with have extended neighborhood */

#if 0 /* Is it no more needed ? */
  assert(cell_cells_idx != NULL);
#endif

  /*---------------------------------------------*
   * Construction des coefficient d'interpolation*
   * sur chaque zone d echange ict               *
   *---------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[ict];

    nbn = 3;
    if (ct->idimct==3) {
      nbn = 4;
    }

    /* Calcul de dh */
    dhi = (ct->hmax-ct->hmin)/(ct->nelect-1);

    /* Copy element centers of the water mesh to an array.*/

    nb_node_water = (int) ple_locator_get_n_dist_points(ct->locat_air_water);

    /* */
    location_cel = ple_locator_get_dist_locations(ct->locat_air_water);

    /* */
    lst_xyz_water   = ple_locator_get_dist_coords(ct->locat_air_water);

    BFT_MALLOC(lst_xyz_cel, ct->nbevct*3, cs_coord_t);
    fvm_nodal_get_element_centers(ct->cell_mesh, CS_INTERLACE, 3, lst_xyz_cel);

    BFT_MALLOC(lst_par_cel, ct->nbevct, cs_lnum_t);
    fvm_nodal_get_parent_num(ct->cell_mesh, 3, lst_par_cel);

    BFT_MALLOC(lst_par_fac, ct->nbfac_ct, cs_lnum_t);
    fvm_nodal_get_parent_num(ct->fac_mesh, 2, lst_par_fac);


    /* boucle sur les noeuds eau */
    for (iwat = 0; iwat < nb_node_water; iwat++) {

      /*--------------------------------------------*
       * Calcul des coord. du noeud a partir de     *
       *  celles du noeud de la face sup            *
       *   Noeud  = NoeudSup + iloc * dh *g / ||g|| *
       *--------------------------------------------*/
      xwat = (cs_real_t) lst_xyz_water[iwat*3    ];
      ywat = (cs_real_t) lst_xyz_water[iwat*3 + 1];
      zwat = (cs_real_t) lst_xyz_water[iwat*3 + 2];

      /*--------------------------------------------*
       * boucle sur les cellules appartenant a la ct*
       * recherche du noeud air le plus proche      *
       *--------------------------------------------*/
      iair = location_cel[iwat] -1;

      /*--------------------------------------------*
       * initialiation matrice d interpolation et   *
       *  tableau nvois[] et coeff[]                *
       *--------------------------------------------*/

      for (ii=0; ii<4; ii++) {
        for (jj=0; jj<4; jj++) {
          pp[jj][ii]   = 0.0;
          ppInv[jj][ii]= 0.0;
        }
      }

      for (jj=0;jj< cs_ctwr_nmaxvoi;jj++) {
        coeff[jj]= -1.0;
        nvois[jj]= -1 ;
      }
      /* fin initialisation */

      /*-------------------------------------------------*
       * Recherche des voisins du noeuds air le + proche *
       *  boucle sur les faces internes du maillage      *
       *-------------------------------------------------*/

      nbvois = 1;
      nvois[0] = iair;

      for (ifac = 0; ifac < ct->nbfac_ct; ifac++) {
        icel_1 = i_face_cells[lst_par_fac[ifac]- mesh->n_b_faces - 1][0];
        icel_2 = i_face_cells[lst_par_fac[ifac]- mesh->n_b_faces - 1][1];
        if (icel_1==iair) {
          nvois[nbvois]=icel_2;
          nbvois += 1;
        }
        else if (icel_2==iair) {
          nvois[nbvois]=icel_1;
          nbvois += 1;
        }
      }

#if 0 /* Is it no more needed ? */
      for (icel = cell_cells_idx[iair    ];
           icel < cell_cells_idx[iair + 1];
           icel++) {

         indice = cell_cells_lst[icel];
         if (ct->mark_ze[indice+1]==1) {
            nvois[nbvois]= indice;
            nbvois += 1;
         }
      }
#endif

      /* fin Recherche */

      /*--------------------------------------------*
       *nombre de voisins insuffisant               *
       *--------------------------------------------*/

      if (nbvois<nbn) {
        nbvois = 1;
        nvois[0] = iair;
        coeff[0] = 1.0;
        goto enregistre;
      }

      dim = ct->idimct;
      vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
      vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
      vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;
      passage2D :;


      /*--------------------------------------------*/
      /* Calcul de l'ouverture de la fonction  de ponderation*/
      /*  egale au max de la distance entre le noeud eau et les voisins air */
      /*--------------------------------------------*/
      ouv = 0.;
      for (ii = 0; ii < nbvois; ii++) {
        iair = nvois[ii];
        dxx  =  (cs_real_t) (coo_cel[iair*3+0] - xwat);
        dyy  =  (cs_real_t) (coo_cel[iair*3+1] - ywat);
        dzz  =  (cs_real_t) (coo_cel[iair*3+2] - zwat);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        aux =  pow(dx,2.)+pow(dy,2.)+pow(dz,2.);
        if (ouv<aux) {
          ouv = aux;
        }
      }
      ouv = sqrt(ouv)*1.1;
      /* fin calcul de l'ouverture */

      /*--------------------------------------------*/
      /*Construction de la matrice A               */
      /*--------------------------------------------*/
      for (ii = 0; ii < nbvois; ii++) {

        indice = nvois[ii];
        dxx = (cs_real_t) (coo_cel[indice*3+0] - xwat);
        dyy = (cs_real_t) (coo_cel[indice*3+1] - ywat);
        dzz = (cs_real_t) (coo_cel[indice*3+2] - zwat);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /* parametre de la fonction de ponderation*/
        cx = 1.0;
        cy = 1.0;
        cz = 1.e10;
        lf = 3;
        epgauss = 5.0;
        if (dim==3) {
          cz = 1.0;
        }
        if (dim==1) {
          cy = 1.e10;
          cz = 1.e10;
        }
        /*fonction de ponderation*/
        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz);
        if (dim == 1) { /* 1D */
          pp[0][0] = w      + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
        }
        else if (dim == 2) { /* 2D */
          pp[0][0] = w       + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[0][2] = w*dy    + pp[0][2];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
          pp[1][2] = w*dx*dy + pp[1][2];
          pp[2][0] = w*dy    + pp[2][0];
          pp[2][1] = w*dx*dy + pp[2][1];
          pp[2][2] = w*dy*dy + pp[2][2];
        }
        else if (dim == 3) {/* 3D */
          pp[0][0] = w       + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[0][2] = w*dy    + pp[0][2];
          pp[0][3] = w*dz    + pp[0][3];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
          pp[1][2] = w*dy*dx + pp[1][2];
          pp[1][3] = w*dz*dx + pp[1][3];
          pp[2][0] = w*dy    + pp[2][0];
          pp[2][1] = w*dx*dy + pp[2][1];
          pp[2][2] = w*dy*dy + pp[2][2];
          pp[2][3] = w*dz*dy + pp[2][3];
          pp[3][0] = w*dz    + pp[3][0];
          pp[3][1] = w*dx*dz + pp[3][1];
          pp[3][2] = w*dy*dz + pp[3][2];
          pp[3][3] = w*dz*dz + pp[3][3];
        }
      }
      /*Fin Construction de la matrice A */
      if (ct->idimct == 3 && dim == 3) {
        dim  = _is_coplanar(coo_cel, nvois, nbvois, pp, vectBase, 0, dhi);

        if (dim!= 3) {
          for (ii=0; ii<3; ii++)
            for (jj=0; jj<3; jj++)
              pp[jj][ii]=0.0;
          goto passage2D;
        }

      }

      /*--------------------------------------------*/
      /* inversion de la matrice par la methode de  */
      /* jordan                                     */
      /*--------------------------------------------*/

      if (_invmat(pp, ppInv, dim) == 0) cs_exit(EXIT_FAILURE);

      /*--------------------------------------------*/
      /* Calcul des coefficients                    */
      /*--------------------------------------------*/
      for (ii = 0; ii < nbvois; ii++) {

        indice = nvois[ii];
        dxx = (cs_real_t) (coo_cel[indice*3+0] - xwat);
        dyy = (cs_real_t) (coo_cel[indice*3+1] - ywat);
        dzz = (cs_real_t) (coo_cel[indice*3+2] - zwat);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /* parametre de la fonction de ponderation*/
        cx = 1.0;
        cy = 1.0;
        cz = 1.e10;
        lf = 3;
        epgauss = 5.0;
        if (dim==3) {
          cz = 1.0;
        }
        if (dim==1) {
          cy = 1.e10;
          cz = 1.e10;
        }

        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz);

        if (dim == 1) {
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx);
        }
        else if (dim ==2) {
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx+ppInv[0][2]*dy);
        }
        else if (dim ==3) {
          coeff[ii] = w*( ppInv[0][0]
                          +ppInv[0][1]*dx
                          +ppInv[0][2]*dy
                          +ppInv[0][3]*dz);
        }

      }
      /* Fin Calcul des coefficients */

      enregistre :;

      /*--------------------------------------------*/
      /* boucle while sur pvoiseau pour trouver le  */
      /* dernier indice                             */
      /*--------------------------------------------*/
      indice = 0;
      while (ct->voiseau[indice]!=-1) {
        indice += 1;
      }
      /*--------------------------------------------*
       * Ajout des voisins et des coefficients      *
       *--------------------------------------------*/
      for (ii = 0; ii < nbvois; ii++) {
        ct->voiseau[indice+ii]  = nvois[ii];
        ct->coefeau[indice+ii]  = coeff[ii];
      }

      ct->pvoiseau[iwat+1] = ct->pvoiseau[iwat] + nbvois;

    }/* fin boucle sur iseg */
    BFT_FREE(lst_par_fac);
    BFT_FREE(lst_xyz_cel);
    BFT_FREE(lst_par_cel);

  } /* fin boucle sur les zones d echanges ict */

}

/*-----------------------------------------------------------------------------*
 * Function cs_ctwr_adair                                                        *
 * Interpolation EAU -> AIR                                                    *
 *-----------------------------------------------------------------------------*/
void cs_ctwr_adair (void)
{
  /* Coordonnees des centres des cellules  */

  const cs_coord_t  *lst_xyz_cel   = NULL;
  cs_coord_t  *lst_xyz_water = NULL;
  const cs_lnum_t  *location_cel  = NULL;
  cs_lnum_t   ict,icol, ilig,ieau,ii,jj,iair,nbvois,nbn,lf,indice;
  cs_lnum_t   nvois[cs_ctwr_nmaxvoi];
  cs_lnum_t   dim, nb_air_node;
  cs_real_t  dhi,dmin,dist,ouv,aux;
  cs_real_t  coeff[cs_ctwr_nmaxvoi];
  cs_real_t  dx,dy,dz,dxx,dyy,dzz;
  cs_real_t  cx,cy,cz,epgauss,w;
  cs_real_t  pp[4][4], ppInv[4][4];
  cs_real_t  vectBase[3][3];
  cs_lnum_t loca_cel;
  cs_ctwr_zone_t  *ct;

   ieau = -1;

  /*---------------------------------------------*
   * Construction des coefficient d'interpolation*
   * sur chaque zone d echange ict               *
   *---------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[ict];

    nbn = 3;
    if (ct->idimct==3) {
      nbn = 4;
    }


    /* Calcul de dh */
    dhi = (ct->hmax-ct->hmin)/(ct->nelect-1);

    /* Copy element centers of the water mesh to an array.*/


    nb_air_node = (int)ple_locator_get_n_dist_points(ct->locat_water_air);

    lst_xyz_cel = ple_locator_get_dist_coords(ct->locat_water_air);

    location_cel = ple_locator_get_dist_locations(ct->locat_water_air);


    BFT_MALLOC(lst_xyz_water, (3*ct->nelect*ct->nnpsct), cs_coord_t);
    fvm_nodal_get_element_centers(ct->water_mesh, CS_INTERLACE,3, lst_xyz_water);

    if (ct->water_halo != NULL) {
      BFT_REALLOC(lst_xyz_water, (3*ct->nelect*ct->nnpsct_with_ghosts), cs_coord_t);
      cs_halo_sync_var_strided(ct->water_halo, ct->halo_type, lst_xyz_water, 3);

    }
    /*--------------------------------------------*
     * Loops on the air nodes of teh exchange area*
     *--------------------------------------------*/

    for (iair = 0; iair < nb_air_node; iair++)  {

       loca_cel = location_cel[iair] -1;
      /*--------------------------------------------*
       * initialiation matrice d interpolation et   *
       * tableau nvois[] et coeff[]                 *
       *--------------------------------------------*/

      for (ii=0; ii<4; ii++) {
        for (jj=0; jj<4; jj++) {
          pp[jj][ii]=0.0;
          ppInv[jj][ii]=0.0;
        }
      }

      for (jj=0; jj< cs_ctwr_nmaxvoi; jj++) {
        coeff[jj]=0.0;
        nvois[jj]= -1;
      }/* fin initialisation */

      /*--------------------------------------------*
       * Traitement particulier pour les noeuds air *
       * en bordure inferieure ou superieure des ct *
       *--------------------------------------------*/

      /* indice du noeud air dans le maillage */
      dmin   = 1000.;

        if ((loca_cel%(ct->nelect) == 0) ||
                  (loca_cel% (ct->nelect) == (ct->nelect-1)))  {
          for (jj = 0; jj < (ct->nelect*ct->nnpsct_with_ghosts); jj++) {
            dx = (cs_real_t) (lst_xyz_water[3*jj  ] - lst_xyz_cel[iair*3  ]);
            dy = (cs_real_t) (lst_xyz_water[3*jj+1] - lst_xyz_cel[iair*3+1]);
            dz = (cs_real_t) (lst_xyz_water[3*jj+2] - lst_xyz_cel[iair*3+2]);
            dist = (pow(dx,2.)+pow(dy,2.)+pow(dz,2.));
            if (dmin>dist) {
              dmin = dist;
              ieau = jj;
            }
          }
        }

      if (dmin<1000.) {
        /* Cellule air est en bordure inf ou sup, on saute l'etape suivante */
        nbvois   = 1;
        nvois[0] = ieau;
        coeff[0] = 1.0;
        goto enregistre;
      }
      /*------------------------------------------------*
       * Fin Traitement particulier pour les noeuds air *
       * en bordure inferieure ou superieure des ct     *
       *------------------------------------------------*/

      /*--------------------------------------------*
       * On continue avec les noeuds air qui ne sont*
       * pas en bordure                             *
       * recherche des cellules air voisines pour   *
       * les noeuds air qui ne sont                 *
       *  pas en bordure inferieure ou superieure   *
       *--------------------------------------------*/
      nbvois = 1;
      nvois[0] = loca_cel;
      /*---------------------------------------------*
       * Recherche du nombre de voisins du noeuds air*
       * boucle sur les faces internes du maillage   *
       *---------------------------------------------*/


      indice=1;
      nvois[indice++] = loca_cel + 1;
      nvois[indice++] = loca_cel - 1;
      nbvois+=2;
      icol = loca_cel/(ct->nelect);
      ilig = loca_cel%(ct->nelect);
      for (ii = ct->fac_sup_connect_idx[icol];
           ii < ct->fac_sup_connect_idx[icol + 1];
           ii++) {

           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ii] + ilig-1;
           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ii] + ilig ;
           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ii] + ilig+1;
           nbvois += 3;

      }

      /*---------------------------------------------*
       * nombre de voisin eau insuffisants           *
       * meme que noeuds en bordure                  *
       *-------------------------------------------- */

      if (nbvois<nbn) {
        nbvois   = 1;
        coeff[0] = 1.0;
        goto enregistre;
      }

       dim = ct->idimct;
       vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
       vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
       vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;

       passage2D :;

      /*--------------------------------------------*
       * Calcul de l'ouverture de la fonction  de   *
       * ponderation egale au max de la distance    *
       *  entre le noeud air et les voisins eau     *
       *--------------------------------------------*/
      ouv = 0.;
      for (ii = 0; ii < nbvois; ii++) {
        ieau = nvois[ii];
        dxx = (cs_real_t) (lst_xyz_water[3*ieau]   - lst_xyz_cel[iair*3+0]);
        dyy = (cs_real_t) (lst_xyz_water[3*ieau+1] - lst_xyz_cel[iair*3+1]);
        dzz = (cs_real_t) (lst_xyz_water[3*ieau+2] - lst_xyz_cel[iair*3+2]);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        aux =  pow(dx,2.)+pow(dy,2.)+pow(dz,2.);
        if (ouv<aux) {
          ouv = aux;
        }
      }
      ouv = sqrt(ouv)*1.1;
      /* Fin de calcul de l'ouverture */

      /*--------------------------------------------*
       *Construction de la matrice A                *
       *--------------------------------------------*/

      for (ii = 0; ii < nbvois; ii++) {

        ieau = nvois[ii];
        dxx = (cs_real_t) (lst_xyz_water[3*ieau +0] - lst_xyz_cel[iair*3+0]);
        dyy = (cs_real_t) (lst_xyz_water[3*ieau +1] - lst_xyz_cel[iair*3+1]);
        dzz = (cs_real_t) (lst_xyz_water[3*ieau +2] - lst_xyz_cel[iair*3+2]);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];
        /* parametre de la fonction de ponderation*/
        cx = 1.0;
        cy = 1.0;
        cz = 1.e10;
        lf = 3;
        epgauss = 5.0;
        if (dim==3) {
          cz = 1.0;
        }
        if (dim==1) {
          cy = 1.e10;
          cz = 1.e10;
        }
        /*fonction de ponderation*/
        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz);

        if (dim == 1) { /* 1D */
          pp[0][0] = w       + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
        }
        else if (dim == 2) { /* 2D */
          pp[0][0] = w     + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[0][2] = w*dy    + pp[0][2];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
          pp[1][2] = w*dx*dy + pp[1][2];
          pp[2][0] = w*dy    + pp[2][0];
          pp[2][1] = w*dx*dy + pp[2][1];
          pp[2][2] = w*dy*dy + pp[2][2];
        }
        else if (dim == 3) { /* 3D */
          pp[0][0] = w       + pp[0][0];
          pp[0][1] = w*dx    + pp[0][1];
          pp[0][2] = w*dy    + pp[0][2];
          pp[0][3] = w*dz    + pp[0][3];
          pp[1][0] = w*dx    + pp[1][0];
          pp[1][1] = w*dx*dx + pp[1][1];
          pp[1][2] = w*dy*dx + pp[1][2];
          pp[1][3] = w*dz*dx + pp[1][3];
          pp[2][0] = w*dy    + pp[2][0];
          pp[2][1] = w*dx*dy + pp[2][1];
          pp[2][2] = w*dy*dy + pp[2][2];
          pp[2][3] = w*dz*dy + pp[2][3];
          pp[3][0] = w*dz    + pp[3][0];
          pp[3][1] = w*dx*dz + pp[3][1];
          pp[3][2] = w*dy*dz + pp[3][2];
          pp[3][3] = w*dz*dz + pp[3][3];
        }

      } /* Fin de construction de la matrice A */
      if (ct->idimct == 3 && dim == 3) {
        dim  = _is_coplanar(lst_xyz_water, nvois, nbvois, pp, vectBase,0,dhi);
        if (dim!= 3) {
          for (ii=0; ii<3; ii++)
            for (jj=0; jj<3; jj++)
              pp[jj][ii]=0.;
          goto passage2D;
        }
      }

      /*--------------------------------------------*
       * inversion de la matrice par la methode de  *
       * jordan                                     *
       *--------------------------------------------*/
      if (_invmat(pp, ppInv, dim) == 0) cs_exit(EXIT_FAILURE);

      /*--------------------------------------------*
       * Calcul des coefficients                    *
       *--------------------------------------------*/
      for (ii = 0; ii < nbvois; ii++) {
        ieau = nvois[ii];
        dxx = (cs_real_t) (lst_xyz_water[3*ieau   ] - lst_xyz_cel[iair*3+0]);
        dyy = (cs_real_t) (lst_xyz_water[3*ieau +1] - lst_xyz_cel[iair*3+1]);
        dzz = (cs_real_t) (lst_xyz_water[3*ieau +2] - lst_xyz_cel[iair*3+2]);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /*parametre de la fonction de ponderation*/
        cx = 1.0;
        cy = 1.0;
        cz = 1.e10;
        lf = 3;
        epgauss = 5.0;
        if (dim==3) {
          cz = 1.0;
        }
        if (dim==1) {
          cy = 1.e10;
          cz = 1.e10;
        }

        w = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz);



        if (dim == 1) {
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx);
        }
        else if (dim == 2) {
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx+ppInv[0][2]*dy);
        }
        else if (dim == 3) {
          coeff[ii] = w*(ppInv[0][0]
                         +ppInv[0][1]*dx
                         +ppInv[0][2]*dy
                         +ppInv[0][3]*dz);
        }

      }
      /* Fin de calcul des coefficients */

      /*---------------------------------------------*
       * note :Reprise pour les noeuds air en bordure*
       * ou avec un nbre de voisin insuffisant       *
       *---------------------------------------------*/
      enregistre :;

      /*--------------------------------------------*
       * trouver le dernier indice sur pvoisair     *
       *--------------------------------------------*/
      indice = 0;
      while (ct->voisair[indice]!=-1) {
        indice += 1;
      }

      /*--------------------------------------------*
       * Ajout des voisins et des coefficients      *
       *--------------------------------------------*/
      for (icol = 0; icol < nbvois; icol++)
        {
          ct->voisair[indice + icol]  = nvois[icol];
          ct->coefair[indice + icol]  = coeff[icol];
         }


      ct->pvoisair[iair+1] = ct->pvoisair[iair] + nbvois;
    }
    /*---------------------------------------------*
     * fin de la boucle sur les noeuds air de la ct*
     *---------------------------------------------*/

    BFT_FREE(lst_xyz_water);


  }/* fin de la boucle sur les ct */

}


/*----------------------------------------------------------------------------*
 * Chaining of the exchange area                                              *
 *----------------------------------------------------------------------------*/

void
cs_ctwr_stacking(void)
{
  cs_lnum_t i, j, rank, dist_rank, nb, nb_ct, itmp, ict, ict_uw;
  cs_lnum_t * aux;
  cs_ctwr_zone_t  *ct, *ct_upw;
  cs_real_t tmp;
  cs_real_t gravite[3];
  cs_coord_t * lst_xyz;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;
  const double tolerance = 0.1;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 1;

  nb = cs_glob_ct_nbr  * cs_glob_ct_nbr;

  BFT_MALLOC(cs_stack_ct, nb, cs_lnum_t);
  BFT_MALLOC(cs_chain_ct, cs_glob_ct_nbr, cs_lnum_t);

  gravite[0]= ct_prop->gravx;
  gravite[1]= ct_prop->gravy;
  gravite[2]= ct_prop->gravz;

  for (i=0; i < cs_glob_ct_nbr; i++)
    for (j=0; j < cs_glob_ct_nbr; j++)
      cs_stack_ct[i*cs_glob_ct_nbr + j]=0;

  for (i=0; i < cs_glob_ct_nbr; i++)
    for (j=0; j < cs_glob_ct_nbr; j++)
      if (CS_ABS(cs_glob_ct_tab[i]->hmax - cs_glob_ct_tab[j]->hmin)< 1.e-6)
        cs_stack_ct[i*cs_glob_ct_nbr + j] =1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    BFT_MALLOC(aux, nb, cs_lnum_t);
    rank = cs_glob_rank_id;

    for (dist_rank = 0; dist_rank < cs_glob_n_ranks; dist_rank++)
      if (dist_rank != rank) {

        MPI_Sendrecv(cs_stack_ct, nb, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                     aux, nb, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                     cs_glob_mpi_comm, &status);
        for (i=0; i < cs_glob_ct_nbr; i++)
          for (j=0; j < cs_glob_ct_nbr; j++) {
            if (aux[i*cs_glob_ct_nbr + j] > cs_stack_ct[i*cs_glob_ct_nbr + j])
              cs_stack_ct[i*cs_glob_ct_nbr + j] = aux[i*cs_glob_ct_nbr + j];
          }
      }

    BFT_FREE(aux);

  }
#endif

  /* to order the exchange area */
    /*Init the chaining array */
  for (i = 0; i < cs_glob_ct_nbr; i++)
    cs_chain_ct[i] = i;

  for (i = 0; i < cs_glob_ct_nbr; i++)
    for (j = i+1; j < cs_glob_ct_nbr; j++)
      if (cs_stack_ct[cs_chain_ct[i]*cs_glob_ct_nbr + cs_chain_ct[j]] == 1) {
        itmp = cs_chain_ct [i];
        cs_chain_ct [i] = cs_chain_ct [j];
        cs_chain_ct [j] = itmp;
      }

  for (ict = 0; ict< cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    nb_ct = 0;

    for (ict_uw = 0; ict_uw < cs_glob_ct_nbr; ict_uw++)
      if (   cs_stack_ct[cs_chain_ct[ict]*cs_glob_ct_nbr + cs_chain_ct[ict_uw]]
          == 1) {

        nb_ct++;
        ct_upw = cs_glob_ct_tab[cs_chain_ct[ict_uw]];

        BFT_MALLOC(lst_xyz,
                   3*(ct_upw->nbfac_ict+ct_upw->nbfbr_ict), cs_coord_t);

        fvm_nodal_get_element_centers
          (ct_upw->face_inf_mesh, CS_INTERLACE,2,lst_xyz);

        tmp  = CS_ABS(ct_upw->hmax - ct_upw->hmin)/(ct_upw->nelect-1);
        tmp /= cs_math_3_norm(gravite);

        for (i=0; i < (ct_upw->nbfac_ict+ct_upw->nbfbr_ict); i++) {
          lst_xyz[3*i + 0] -= tmp * gravite[0];
          lst_xyz[3*i + 1] -= tmp * gravite[1];
          lst_xyz[3*i + 2] -= tmp * gravite[2];
        }

        BFT_REALLOC(ct->locat_cell_ct_upwind, nb_ct, ple_locator_t *);

#if defined(PLE_HAVE_MPI)
        ct->locat_cell_ct_upwind[nb_ct-1] =
          ple_locator_create(cs_glob_mpi_comm,
                             cs_glob_n_ranks,
                             0);
#else
        ct->locat_cell_ct_upwind[nb_ct-1] = ple_locator_create();
#endif

        ple_locator_set_mesh(ct->locat_cell_ct_upwind[nb_ct-1],
                             ct_upw->water_mesh,
                             locator_options,
                             0,
                             tolerance,
                             3,
                             ct_upw->nbfac_ict+ct_upw->nbfbr_ict,
                             NULL,
                             NULL,
                             lst_xyz,
                             NULL,
                             cs_coupling_mesh_extents,
                             cs_coupling_point_in_mesh);
        BFT_FREE(lst_xyz);

      }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
