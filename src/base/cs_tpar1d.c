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
 * Modelling the thermal wall with 1D approach
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_suite.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_tpar1d.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure definitions
 *============================================================================*/

struct par1d
{
  cs_int_t    n;     /* Number of discretization points on the coupled face */
  cs_real_t  *z;     /* Discretization points coordinates                   */
  cs_real_t   e;     /* Thickness associated to the coupled face            */
  cs_real_t  *t;     /* Temperature at each point of discretization         */
};


/*============================================================================
 * Static global variable
 *============================================================================*/

static struct par1d *cs_glob_par1d = NULL;
static cs_suite_t   *cs_glob_tpar1d_suite = NULL;


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate the cs_glob_par1d structure
 *
 * parameters:
 *   nfpt1d <- number of coupled boundary faces
 *   nppt1d <- number of discretization points on each coupled faces
 *----------------------------------------------------------------------------*/

static void
cs_loc_tpar1d_cree(cs_int_t         nfpt1d,
                   const cs_int_t  *nppt1d)
{
  cs_int_t   nb_pts_tot;
  cs_int_t   i;

  if (nfpt1d == 0)
    return;

  /* Allocate the cs_glob_par1d structure */
  BFT_MALLOC(cs_glob_par1d, nfpt1d, struct par1d);

  /* Initialization of the number of discretization points in each structure
     Computation of the toatl number of discretization points */
  nb_pts_tot = 0;

  for (i = 0; i < nfpt1d; i++) {
    cs_glob_par1d[i].n = nppt1d[i];
    nb_pts_tot += nppt1d[i];
  }

  /* Allocate the "t" arrays: Temperature in each point of discretization
          and the "z" arrays: Coordonnates of each point of discretization */

  BFT_MALLOC(cs_glob_par1d->z, 2 * nb_pts_tot, cs_real_t);
  cs_glob_par1d->t = cs_glob_par1d->z + nb_pts_tot;

  for (i = 1; i < nfpt1d; i++) {
    cs_glob_par1d[i].z = cs_glob_par1d[i-1].z + nppt1d[i-1];
    cs_glob_par1d[i].t = cs_glob_par1d[i-1].t + nppt1d[i-1];
  }

}

/*----------------------------------------------------------------------------
 * Open the restart file associated to cs_tpar1d
 * Allocate cs_glob_tpar1d_suite
 *
 * parameters:
 *   nomsui <- name of the restart file
 *   lngnom <- name length
 *   ireawr <- 1 for reading, 2 for writing
 *----------------------------------------------------------------------------*/

static void
cs_loc_tpar1d_opnsuite(const char            *nomsui,
                       const cs_int_t        *lngnom,
                       const cs_suite_mode_t  ireawr)
{
  char            *nombuf;

  /* Name treatment for the C API */
  nombuf = cs_base_string_f_to_c_create(nomsui, *lngnom);

  cs_glob_tpar1d_suite = cs_suite_cree(nombuf, ireawr);

  /* Free the memory if necessary */
  cs_base_string_f_to_c_free(&nombuf);
}

/*============================================================================
 *  Public functions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create the 1D mesh for each face and initialize the temperature
 *
 * Fortran interface:
 *
 * SUBROUTINE  MAIT1D
 * ******************
 *
 * INTEGER          NFPT1D         : <-  : number of coupled faces
 * INTEGER          NPPT1D(NFPT1D) : <-  : number of mesh points for each face
 * DOUBLE PRECISION EPPT1D(NFPT1D) : <-  : wall thickness for each face
 * DOUBLE PRECISION RGPT1D(NFPT1D) : <-  : mesh geometric ratio for each face
 * DOUBLE PRECISION TPPT1D(NFPT1D) : <-  : temperature initizalition value
 *----------------------------------------------------------------------------*/

void CS_PROCF (mait1d,MAIT1D)
(
 cs_int_t   *nf,
 cs_int_t    n[],
 cs_real_t   e[],
 cs_real_t   r[],
 cs_real_t   tp[]
)
{
  cs_int_t i, k;
  cs_real_t m, rr;
  cs_real_t *zz;

  /* Allocate the global structure: cs_glob_par1d and the number of
     discretization points on each face */
  cs_loc_tpar1d_cree(*nf, n);

  /* Initialization of each coupled face thickness */
  for (i = 0; i < *nf; i++) {
    cs_glob_par1d[i].e = e[i];
  }

  for (i = 0; i < *nf; i++) {

    /* Initialization of the Temperature */
    for (k = 0; k<n[i]; k++) {
      (cs_glob_par1d[i].t)[k] = tp[i];
    }

    /* Mesh */
    zz = cs_glob_par1d[i].z;
    rr = r[i];

    /* Regular */
    if (fabs(rr-1.0) <= 1.0e-6) {
      zz[0] = e[i]/n[i]/2.;
      for (k = 1; k < n[i]; k++) {
        zz[k]=zz[k-1]+e[i]/n[i];
      }
    }

    /* Geometric */
    else {
      m = e[i]*(1.0-rr)/(1.0-pow(rr,n[i]));
      *zz = m/2.;
      for (k = 1; k< n[i]; k++) {
        zz[k] = zz[k-1]+m/2.;
        m = m*rr;
        zz[k] = zz[k]+m/2.;
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Solve the 1D equation for a given face
 *
 * Fortran interface:
 *
 * SUBROUTINE  TPAR1D
 * ******************
 *
 * INTEGER          II     : <-  : face number
 * INTEGER          ICLT1D : <-  : type of exterior boundary condition
 * DOUBLE PRECISION TBORD  : <-  : fluid temperature at the boundary
 * DOUBLE PRECISION HBORD  : <-  : exchange coefficient for the fluid
 *                         :     : at the boundary
 * DOUBLE PRECISION TET1D  : <-  : temperature on the exterior boundary
 *                         :     : (Dirichlet boundary condition)
 * DOUBLE PRECISION HET1D  : <-  : exchange coefficient on the exterior wall
 * DOUBLE PRECISION FET1D  : <-  : flux on the exterior wall
 *                         :     : (Neumann boundary condition)
 * DOUBLE PRECISION LAMT1D : <-  : conductivity (lambda)
 * DOUBLE PRECISION RCPT1D : <-  : rho*Cp product
 * DOUBLE PRECISION DTPT1D : <-> : time-step for the solid resolution
 * DOUBLE PRECISION TPPT1D : <-> : physical temperature at the fluid/solid
 *                         :     : interface
 *----------------------------------------------------------------------------*/

void CS_PROCF (tpar1d,TPAR1D)
(
 cs_int_t *ii,
 cs_int_t *icdcle,
 cs_real_t *tf,
 cs_real_t *hf,
 cs_real_t *te,
 cs_real_t *he,
 cs_real_t *fe,
 cs_real_t *lb,
 cs_real_t *rocp,
 cs_real_t *dtf,
 cs_real_t *tp
)
{
  cs_int_t k;

  cs_real_t a1; /* extrapolation coefficient for temperature1 */
  cs_real_t h2; /* thermal exchange coefficient on T(1) */
  cs_real_t f3; /* thermal flux on Tfluide */
  cs_real_t a4; /* extrapolation coefficient for temperature4 */
  cs_real_t h5; /* thermal exchange coefficient on T(n) */
  cs_real_t f6; /* thermal flux on Text */

  cs_real_t m;

  cs_real_t *al, *bl, *cl, *dl;
  cs_real_t *zz;
  cs_int_t n;

  n = cs_glob_par1d[*ii].n;

  BFT_MALLOC(al, 4*n, cs_real_t);
  bl = al+n;
  cl = bl+n;
  dl = cl+n;

  zz = cs_glob_par1d[*ii].z;

  /* Build the tri-diagonal matrix */

  /* Boundary conditions on the fluid side: flux conservation */
  /*   flux in the fluid = flux in the solid = f3 + h2*T1 */
  a1 = 1./(*hf)+zz[0]/(*lb);
  h2 = -1./a1;
  f3 = -h2*(*tf);

  /* Boundary conditions on the exterior */
  /*   flux in the fluid = flux in the solid = f6 + h5*T(n-1) */

  /* Dirichlet condition */
  if (*icdcle == 1) {
    a4 = 1./(*he)+(cs_glob_par1d[*ii].e - zz[n-1])/(*lb);
    h5 = -1./a4;
    f6 = -h5*(*te);
  }
  /* Forced flux condition */
  else if (*icdcle == 3) {
    h5 = 0.;
    f6 = *fe;
  }

  /* Mesh interior points */
  for (k=1; k <= n-1; k++) {
    al[k] = -(*lb)/(zz[k]-zz[k-1]);
  }

  m = 2*zz[0];
  for (k=1; k <= n-2; k++) {
    m = 2*(zz[k]-zz[k-1])-m;
    bl[k] = (*rocp)/(*dtf)*m +(*lb)/(zz[k+1]-zz[k]) +(*lb)/(zz[k]-zz[k-1]);
  }

  for (k=0; k <= n-2; k++) {
    cl[k] =  -(*lb)/(zz[k+1]-zz[k]);
  }

  m = 2*zz[0];
  dl[0] = (*rocp)/(*dtf)*m*(cs_glob_par1d[*ii].t)[0];

  for (k=1; k <= n-1; k++) {
    m = 2*(zz[k]-zz[k-1])-m;
    dl[k] = (*rocp)/(*dtf)*m*(cs_glob_par1d[*ii].t)[k];
  }

  /* Boundary points */
  /* bl[0] and bl[n-1] are initialized here and set later,
     in the case where 0 = n-1 */
  bl[0] = 0.;
  bl[n-1] = 0.;
  al[0] = 0.;
  bl[0] = bl[0] + (*rocp)/(*dtf)*2*zz[0] + (*lb)/(zz[1]-zz[0]) - h2;
  cl[0] = cl[0];
  dl[0] = dl[0] +f3;
  al[n-1] = al[n-1];
  bl[n-1] =   bl[n-1] + (*rocp)/(*dtf)*2*(cs_glob_par1d[*ii].e-zz[n-1])
            + (*lb)/(zz[n-1]-zz[n-2]) -h5;
  cl[n-1] = 0.;
  dl[n-1] = dl[n-1] +f6;

  /* System resolution by a Cholesky method ("dual-scan") */
  for (k=1; k<=n-1; k++) {
    bl[k] = bl[k] -al[k]*cl[k-1]/bl[k-1];
    dl[k] = dl[k] -al[k]*dl[k-1]/bl[k-1];
  }

  cs_glob_par1d[*ii].t[n-1] = dl[n-1]/bl[n-1];

  for (k=n-2; k>=0; k--) {
    cs_glob_par1d[*ii].t[k] = (dl[k] -cl[k]*cs_glob_par1d[*ii].t[k+1])/bl[k];
  }


  /* Compute the new value of tp */
  *tp = (*hf)+(*lb)/zz[0];
  *tp = 1/(*tp)*((*lb)*cs_glob_par1d[*ii].t[0]/zz[0]+(*hf)*(*tf));

  BFT_FREE(al);
}

/*----------------------------------------------------------------------------
 * Read the restart file of the 1D-wall thermal module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECT1D
 * *****************
 *
 * CHARACTER        NOMSUI : <-- : Name of the restart file
 * INTEGER          LNGNOM : <-- : Name length
 * INTEGER          NFPT1D : <-- : Number of coupled faces
 * INTEGER          NFPT1T : <-- : Total number of coupled faces
 * INTEGER          NMXT1D : <-- : Max number of points on the 1D meshes
 * INTEGER          NFABOR : <-- : Number of boundary faces
 * INTEGER          NPPT1D : <-- : Number of points of each face 1D-mesh
 * INTEGER          IFPT1D : <-- : Indirection array for 1D-module faces
 * DOUBLE PRECISION EPPT1D : <-- : Wall thickness of each face
 * DOUBLE PRECISION RGPT1D : <-- : Geometric reason associated to faces
 * DOUBLE PRECISION TPPT1D : --> : Wall temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (lect1d,LECT1D)
(
 const char       *const nomsui,
 const cs_int_t   *const lngnom,
 const cs_int_t   *const nfpt1d,
 const cs_int_t   *const nfpt1t,
 const cs_int_t   *const nmxt1d,
 const cs_int_t   *const nfabor,
 const cs_int_t   *const nppt1d,
 const cs_int_t   *const ifpt1d,
 const cs_real_t  *const eppt1d,
 const cs_real_t  *const rgpt1d,
       cs_real_t  *const tppt1d
 CS_ARGF_SUPP_CHAINE
)
{
  cs_bool_t           corresp_cel, corresp_fac, corresp_fbr, corresp_som;
  cs_int_t            nbvent;
  cs_int_t            i, j, ifac, indfac, ierror;
  cs_int_t            version;    /* Not used at the moment */

  cs_suite_t          *suite;
  cs_suite_mode_t     suite_mode;
  cs_suite_support_t  support;
  cs_type_t           typ_val;


  ierror = CS_SUITE_SUCCES;
  suite_mode = CS_SUITE_MODE_LECTURE;

  /* Open the restart file */
  cs_loc_tpar1d_opnsuite(nomsui,
                         lngnom,
                         suite_mode);

  if (cs_glob_tpar1d_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the 1D-wall thermal module restart file "
                "in read mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              *nomsui);


  /* Pointer to the global restart structure */
  suite = cs_glob_tpar1d_suite;

  /* Verification of the associated "support" to the restart file */
  cs_suite_verif_support_base(suite, &corresp_cel, &corresp_fac,
                              &corresp_fbr, &corresp_som);

  /* Only boundary faces are of interest */
  indfac = (corresp_fbr == true ? 1 : 0);
  if (indfac == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while reading the 1D-wall thermal module restart file.\n"
                "The number of boundary faces has been modified\n"
                "Verify that the restart file corresponds to "
                "the present study.\n"));


  { /* Read the header */
    char       nomrub[] = "version_fichier_suite_module_1d";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_SCAL;
    typ_val = CS_TYPE_cs_int_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       INCORRECT FILE TYPE\n"
                  "\n"
                  "The file %s does not seem to be a restart file\n"
                  "for the 1D-wall thermal module.\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to a\n"
                  "restart file for the 1D-wall thermal module.\n"),
                *nomsui);

    version = *tabvar;

    BFT_FREE(tabvar);
  }

  { /* Read the number of discretiaztion points and coherency checks
       with the input data of USPT1D */
    char       nomrub[] = "nb_pts_discretis";
    cs_int_t   *tabvar;
    cs_int_t   mfpt1d, mfpt1t;
    cs_int_t   iok;

    BFT_MALLOC(tabvar, *nfabor, cs_int_t);

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_int_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Coherency checks between the read NFPT1T and the one from USPT1D */
    mfpt1d = 0;
    for (ifac = 0; ifac < *nfabor; ifac++) {
      if (tabvar[ifac] > 0) mfpt1d++;
    }
    mfpt1t = mfpt1d;
    /* if necessary, sum over all the processors */
#if defined(_CS_HAVE_MPI)
    if (cs_glob_base_nbr > 1)
      MPI_Allreduce (&mfpt1d, &mfpt1t, 1, CS_MPI_INT, MPI_SUM,
                     cs_glob_base_mpi_comm);
#endif
    if (mfpt1t != *nfpt1t)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "The number of faces with 1D thermal module has\n"
                  "been modified.\n"
                  "PREVIOUS: %d boundary faces (total)\n"
                  "CURRENT:  %d boundary faces (total)\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to a\n"
                  "restart file for the 1D-wall thermal module.\n"
                  "Verify uspt1d.\n"), mfpt1t, *nfpt1t);

    /* Coherency check between read NFPT1D/IFPT1D and the ones from USPT1D
       One already knows that the number of faces are equal, it is then
       enough to check that all selected faces in USPT1D were also
       selected in the previous calculation */
    iok = 0;
    for (i = 0; i < *nfpt1d; i++){
        ifac = ifpt1d[i]-1;
        if (tabvar[ifac] != nppt1d[i]) iok++;
    }
    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "IFPT1D or NPPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file correspond to\n"
                  "the present study.\n"
                  "Verify uspt1d\n"
                  "(refer to the user manual for the specificities\n"
                  "of the test on IFPT1D)"));

    /* Allocate the cs_glob_par1d structure */

    cs_loc_tpar1d_cree (*nfpt1d, nppt1d);

    BFT_FREE(tabvar);
  }

  { /* Read the wall thickness and check the coherency with USPT1D*/
    char        nomrub[] = "epaisseur_paroi";
    cs_real_t   *tabvar;
    cs_int_t    iok;

    BFT_MALLOC(tabvar, *nfabor, cs_real_t);

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Coherence check between the read EPPT1D and the one from USPT1D */
    iok = 0;
    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i]-1;
      if (fabs(tabvar[ifac]-eppt1d[i])/eppt1d[i] > 1.e-10) iok++;
    }
    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                  "\n"
                  "The parameter EPPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                  "Verify that the restart file corresponds to\n"
                  "the present study.\n"
                  "Verify uspt1d\n"));

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;
      cs_glob_par1d[i].e = tabvar[ifac];
    }

    BFT_FREE(tabvar);
  }

  { /* Read the interior boundary temperature */
    char       nomrub[] = "temperature_bord_int";
    cs_real_t  *tabvar;

    BFT_MALLOC(tabvar, *nfabor, cs_real_t);

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;
      tppt1d[i] = tabvar[ifac];
    }

    BFT_FREE(tabvar);
  }

  { /* Read the 1D-mesh coordinates */
    char        nomrub[] = "coords_maillages_1d";
    cs_int_t    nptmx;
    cs_int_t    iok;
    cs_real_t   *tabvar;
    cs_real_t   zz1, zz2, rrgpt1;

    nptmx = (*nfabor) * (*nmxt1d);
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    nbvent  = *nmxt1d;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES)
      bft_error(__FILE__, __LINE__, 0,
                _("Problem while reading section in the restart file\n"
                  "for the 1D-wall thermal module:\n"
                  "<%s>\n"
                  "The calculation will not be run.\n"), nomrub);

    /* Now one have the cell centers, RGPT1D can be tested */
    iok = 0;
    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i]-1;
      if (nppt1d[i] > 1) {
        zz1 = tabvar[0 + (*nmxt1d)*ifac];
        zz2 = tabvar[1 + (*nmxt1d)*ifac];
        rrgpt1 = (zz2-2.*zz1)/zz1;
        if (fabs(rrgpt1-rgpt1d[i])/rgpt1d[i] > 1.e-10) iok++;
      }
    }

    if (iok > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                  "********               1D-WALL THERMAL MODULE\n"
                  "       CURRENT AND OLD DATA ARE DIFFERENT\n"
                  "\n"
                  "The parameter RGPT1D has been modified with respect\n"
                  "to the restart file on at least on face with\n"
                  "1D thermal module\n"
                  "\n"
                  "The calculation will not be run.\n"
                  "\n"
                   "Verify that the restart file correspond to\n"
                  "the present study\n"
                  "Verify uspt1d\n"));

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i]-1;
      /* The array is filled until the number of discretization points of
         the given face is reached */
      for (j = 0; j < cs_glob_par1d[i].n; j++)
        cs_glob_par1d[i].z[j] = tabvar[j + (*nmxt1d)*ifac];
    }

    BFT_FREE(tabvar);
  }

  { /* Read the wall temperature */
    char        nomrub[] = "temperature_interne";
    cs_int_t    nptmx;
    cs_real_t   *tabvar;

    nptmx = (*nfabor) * (*nmxt1d);
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    nbvent  = *nmxt1d;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    ierror = cs_suite_lit_rub(suite,
                              nomrub,
                              support,
                              nbvent,
                              typ_val,
                              tabvar);

    if (ierror < CS_SUITE_SUCCES) {
      cs_base_warn(__FILE__,__LINE__);
      bft_printf(_("Problem while reading the section in the restart file\n"
                   "for the 1D-wall thermal module:\n"
                   "<%s>\n"), nomrub);
    }

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached */
      for (j = 0; j < cs_glob_par1d[i].n; j++)
        cs_glob_par1d[i].t[j] = tabvar[j + (*nmxt1d)*ifac];

    }

    BFT_FREE(tabvar);
  }

  /* Close the restart file and free structures */
  cs_suite_detruit(cs_glob_tpar1d_suite);
  cs_glob_tpar1d_suite = NULL;
}

/*----------------------------------------------------------------------------
 * Write the restart file of the 1D-wall thermal module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECT1D
 * *****************
 *
 * CHARACTER        NOMSUI : <-- : Name of the restart file
 * INTEGER          LNGNOM : <-- : Name length
 * INTEGER          NFPT1D : <-- : Number of coupled faces
 * INTEGER          NMXT1D : <-- : Max number of points on the 1D meshes
 * INTEGER          NFABOR : <-- : Number of boundary faces
 * DOUBLE PRECISION TPPT1D : --> : Wall temperature
 * INTEGER          IFPT1D : <-- : Indirection array for 1D-module faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrt1d,ECRT1D)
(
 const char       *const nomsui,
 const cs_int_t   *const lngnom,
 const cs_int_t   *const nfpt1d,
 const cs_int_t   *const nmxt1d,
 const cs_int_t   *const nfabor,
 const cs_real_t  *const tppt1d,
 const cs_int_t   *const ifpt1d
 CS_ARGF_SUPP_CHAINE
)
{
  cs_int_t            nbvent, ierror;
  cs_int_t            i, j, ifac;

  cs_suite_t          *suite;
  cs_suite_support_t  support;
  cs_suite_mode_t     suite_mode;
  cs_type_t           typ_val;


  ierror = CS_SUITE_SUCCES;
  suite_mode = CS_SUITE_MODE_ECRITURE;

  /* Open the restart file */
  cs_loc_tpar1d_opnsuite(nomsui,
                         lngnom,
                         suite_mode);

  if (cs_glob_tpar1d_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the 1D-wall thermal module restart "
                "file in write mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              *nomsui);


  /* Pointer to the global restart structure */
  suite = cs_glob_tpar1d_suite;

  { /* Write the header */
    char       nomrub[] = "version_fichier_suite_module_1d";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, 1, cs_int_t);

    *tabvar = 120;

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_SCAL;
    typ_val = CS_TYPE_cs_int_t;

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the number of discretization points */
    char       nomrub[] = "nb_pts_discretis";
    cs_int_t   *tabvar;

    BFT_MALLOC(tabvar, *nfabor, cs_int_t);

    for (i = 0; i < *nfabor; i++)
      tabvar[i] = 0;

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_int_t;

    for (i = 0; i < *nfpt1d; i++) {
            ifac = ifpt1d[i] - 1;
            tabvar[ifac] = cs_glob_par1d[i].n;
          }

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the wall thickness */
    char        nomrub[] = "epaisseur_paroi";
    cs_real_t   *tabvar;

    BFT_MALLOC(tabvar, *nfabor, cs_real_t);

    for (i = 0; i < *nfabor; i++)
      tabvar[i] = 0.0;

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    for (i = 0; i < *nfpt1d; i++) {
            ifac = ifpt1d[i] - 1;
            tabvar[ifac] = cs_glob_par1d[i].e;
          }

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the internal wall-boundary temperature */
    char       nomrub[] = "temperature_bord_int";
    cs_real_t  *tabvar;

    BFT_MALLOC(tabvar, *nfabor, cs_real_t);

    for (i = 0; i < *nfabor; i++)
      tabvar[i] = 0.0;

    nbvent  = 1;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;
      tabvar[ifac] = tppt1d[i];
    }

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the 1D-mesh coordinates */
    char        nomrub[] = "coords_maillages_1d";
    cs_int_t    nptmx;
    cs_real_t   *tabvar;

    nptmx = (*nfabor) * (*nmxt1d);
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    for (i = 0; i < nptmx; i++)
      tabvar[i] = 0.0;

    nbvent  = *nmxt1d;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached (the following cases, up to nmxt1d, are
         already set to 0 during the initalization of tabvar */
      for (j = 0; j < cs_glob_par1d[i].n; j++)
        tabvar[j + (*nmxt1d)*ifac] = cs_glob_par1d[i].z[j];
    }

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  { /* Write the wall-interior temperature */
    char        nomrub[] = "temperature_interne";
    cs_int_t    nptmx;
    cs_real_t   *tabvar;

    nptmx = (*nfabor) * (*nmxt1d);
    BFT_MALLOC(tabvar, nptmx, cs_real_t);

    for (i = 0; i < nptmx; i++)
      tabvar[i] = 0.0;

    nbvent  = *nmxt1d;
    support = CS_SUITE_SUPPORT_FAC_BRD;
    typ_val = CS_TYPE_cs_real_t;

    for (i = 0; i < *nfpt1d; i++) {
      ifac = ifpt1d[i] - 1;

      /* The array is filled until the number of discretization points of
         the given face is reached (the following cases, up to nmxt1d, are
         already set to 0 during the initalization of tabvar */
      for (j = 0; j < cs_glob_par1d[i].n; j++)
        tabvar[j + (*nmxt1d)*ifac] = cs_glob_par1d[i].t[j];

    }

    cs_suite_ecr_rub(suite,
                     nomrub,
                     support,
                     nbvent,
                     typ_val,
                     tabvar);

    BFT_FREE(tabvar);
  }

  /* Close the restart file and free structures */
  cs_suite_detruit(cs_glob_tpar1d_suite);
  cs_glob_tpar1d_suite = NULL;

}


/*----------------------------------------------------------------------------
 * Free memory
 *
 * Interface Fortran :
 *
 * SUBROUTINE  LBRT1D ()
 * ******************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (lbrt1d,LBRT1D)(void)
{
  BFT_FREE(cs_glob_par1d->z);
  BFT_FREE(cs_glob_par1d);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
