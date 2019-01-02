/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_math.h"
#include "cs_sdm.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Static global variables
 *============================================================================*/

static FILE  *sdm = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform several basic tests/operations on cs_sdm_t structures
 */
/*----------------------------------------------------------------------------*/

static void
_test_sdm(FILE  *out)
{
  cs_real_33_t  mpty = {{1.0, 0.5, 0.0},
                        {0.5, 1.0, 0.5},
                        {0.0, 0.5, 1.0}};
  cs_real_t  eigen_ratio, eigen_max;

  /* Useful for a weak enforcement of the BC */
  cs_math_33_eigen((const cs_real_t (*)[3])mpty,
                   &eigen_ratio, &eigen_max);

  fprintf(out, " Matrix property eig: ratio % .4e max. % .4e\n",
          eigen_ratio, eigen_max);

  {
    fprintf(out, "\n Matrix factorization\n");

    const int  max_size = 6;
    cs_sdm_t  *m = cs_sdm_square_create(max_size);

    cs_sdm_square_init(3, m);
    m->val[0] = 2, m->val[1] = -1, m->val[2] = 0;
    m->val[3] =-1, m->val[4] =  2, m->val[5] =-1;
    m->val[6] = 0, m->val[7] = -1, m->val[8] = 1;

    cs_real_6_t  b = {1, 0, 0, 0, 0, 0};

    /* Compute the L.D.L^T decomposition and then solve */
    cs_real_6_t  sol, tmp;
    cs_real_t  facto[21];

    cs_sdm_33_ldlt_compute(m, facto);
    cs_sdm_33_ldlt_solve(facto, b, sol);

    fprintf(out, "\n3x3 matrix\n");
    cs_sdm_fprintf(out, NULL, cs_math_zero_threshold, m);

    fprintf(out, " Solution l.d.l^T 33: % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2]);

    cs_sdm_ldlt_compute(m, facto, tmp);
    cs_sdm_ldlt_solve(3, facto, b, sol);

    fprintf(out, " Solution l.d.l^T:    % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2]);

    cs_sdm_square_init(4, m);
    m->val[ 0] = 2, m->val[ 1] = -1, m->val[ 2] = 0, m->val[ 3] = 0;
    m->val[ 4] =-1, m->val[ 5] =  2, m->val[ 6] =-1, m->val[ 7] = 0;
    m->val[ 8] = 0, m->val[ 9] = -1, m->val[10] = 2, m->val[11] =-1;
    m->val[12] = 0, m->val[13] =  0, m->val[14] =-1, m->val[15] = 1;

    cs_sdm_44_ldlt_compute(m, facto);
    cs_sdm_44_ldlt_solve(facto, b, sol);

    fprintf(out, "\n4x4 matrix\n");
    cs_sdm_fprintf(out, NULL, cs_math_zero_threshold, m);

    fprintf(out, " Solution l.d.l^T 44: % .4e % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2], sol[3]);

    cs_sdm_ldlt_compute(m, facto, tmp);
    cs_sdm_ldlt_solve(4, facto, b, sol);

    fprintf(out, " Solution l.d.l^T   : % .4e % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2], sol[3]);

    cs_sdm_square_init(6, m);
    cs_real_t *a = m->val;
    a[ 0] = 2, a[ 1] = -1, a[ 2] = 0, a[ 3] = 0, a[ 4] =  0, a[ 5] =  0;
    a[ 6] =-1, a[ 7] =  2, a[ 8] =-1, a[ 9] = 0, a[10] =  0, a[11] =  0;
    a[12] = 0, a[13] = -1, a[14] = 2, a[15] =-1, a[16] =  0, a[17] =  0;
    a[18] = 0, a[19] =  0, a[20] =-1, a[21] = 2, a[22] = -1, a[23] =  0;
    a[24] = 0, a[25] =  0, a[26] = 0, a[27] =-1, a[28] =  2, a[29] = -1;
    a[30] = 0, a[31] =  0, a[32] = 0, a[33] = 0, a[34] = -1, a[35] =  1;

    cs_sdm_66_ldlt_compute(m, facto);
    cs_sdm_66_ldlt_solve(facto, b, sol);

    fprintf(out, "\n6x6 matrix\n");
    cs_sdm_fprintf(out, NULL, cs_math_zero_threshold, m);

    fprintf(out, " Solution l.d.l^T 66: % .4e % .4e % .4e % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2], sol[3], sol[4], sol[5]);

    cs_sdm_ldlt_compute(m, facto, tmp);
    cs_sdm_ldlt_solve(6, facto, b, sol);

    fprintf(out, " Solution l.d.l^T   : % .4e % .4e % .4e % .4e % .4e % .4e\n",
            sol[0], sol[1], sol[2], sol[3], sol[4], sol[5]);

    m = cs_sdm_free(m);
  }

  { /* Test symmetry */
    const int  max_size = 6;
    cs_sdm_t  *m = cs_sdm_square_create(max_size);

    fprintf(out, "\n Test symmetry (non-block version)\n");
    cs_sdm_square_init(3, m);
    m->val[0] = 2, m->val[1] = -1, m->val[2] = 0;
    m->val[3] =-1, m->val[4] =  2, m->val[5] =-1;
    m->val[6] = 0, m->val[7] = -1, m->val[8] = 1;

    cs_sdm_fprintf(out, NULL, cs_math_zero_threshold, m);

    double  eval_sym = cs_sdm_test_symmetry(m);
    fprintf(out, " symmetry evaluation = %g\n", eval_sym);

    cs_sdm_square_init(3, m);
    m->val[0] = 2, m->val[1] = -1, m->val[2] = 0.25;
    m->val[3] =-3, m->val[4] =  2, m->val[5] =-1;
    m->val[6] = 0, m->val[7] = -0.5, m->val[8] = 1;

    cs_sdm_fprintf(out, NULL, cs_math_zero_threshold, m);

    double  eval_sym2 = cs_sdm_test_symmetry(m);
    fprintf(out, " symmetry evaluation = %g\n", eval_sym2);

    m = cs_sdm_free(m);
  }


  /* Use block-matrix */
  {
    fprintf(out, "\n Test copy (block version)\n");

    short int  bsize[3] = {2, 2, 3};

    cs_sdm_t  *mb = cs_sdm_block_create(3, 3, bsize, bsize);

    /* Set values */
    cs_sdm_block_init(mb, 3, 3, bsize, bsize);

    cs_sdm_t  *b11 = cs_sdm_get_block(mb, 1, 1);
    b11->val[0] = b11->val[1] = b11->val[2] = 1;
    b11->val[3] = 2;

    /* 3 rows and 2 columns */
    cs_sdm_t  *b21 = cs_sdm_get_block(mb, 2, 1);
    b21->val[0] = 0.5, b21->val[1] = 0.25;
    b21->val[2] = 1, b21->val[3] = 0.75;
    b21->val[4] = 2, b21->val[5] = 0.1;

    cs_sdm_t  *b12 = cs_sdm_get_block(mb, 1, 2);
    /* 2 rows and 3 columns */
    b12->val[0] = 0.5, b12->val[1] = 0.25, b12->val[2] = 1;
    b12->val[3] = 0.75, b12->val[4] = 2, b12->val[5] = 0.1;

    /* cs_sdm_block_dump(0, mb); */
    fprintf(out, " Reference matrix\n");
    cs_sdm_block_fprintf(out, NULL, cs_math_zero_threshold, mb);

    cs_sdm_t  *cpy = cs_sdm_block_create_copy(mb);

    fprintf(out, " Copy of the previous matrix\n");
    cs_sdm_block_fprintf(out, NULL, cs_math_zero_threshold, cpy);

    fprintf(out, "\n Test symmetry (block version)\n");
    fprintf(out, " Symmetry evaluation for the reference = %g\n",
            cs_sdm_test_symmetry(mb));

    b12 = cs_sdm_get_block(cpy, 1, 2);    /* 2 rows and 3 columns */
    b12->val[0] = 0.5, b12->val[1] = 1, b12->val[2] = 2;
    b12->val[3] = 0.25, b12->val[4] = 0.75, b12->val[5] = 0.1;

    fprintf(out, " Symmetric matrix defined by block\n");
    cs_sdm_block_fprintf(out, NULL, cs_math_zero_threshold, cpy);
    fprintf(out, " Symmetry evaluation for a symmetric matrix = %g\n",
            cs_sdm_test_symmetry(cpy));

    mb = cs_sdm_free(mb);
    cpy = cs_sdm_free(cpy);
  }

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main program to check CDO/HHO algorithms
 *
 * \param[in]    argc
 * \param[in]    argv
 */
/*----------------------------------------------------------------------------*/

int
main(int    argc,
     char  *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  {
    int t_id;
#pragma omp parallel private(t_id)
    {
      t_id = omp_get_thread_num();
      if (t_id == 0)
        cs_glob_n_threads = omp_get_max_threads();
    }
  }
#endif

  sdm = fopen("SDM_tests.log", "w");

  /* =======================================
   * TEST on Small Dense Matrices operations
   * ======================================= */

  _test_sdm(sdm);

  fclose(sdm);

  printf("\n\n -->> SDM Tests (Done)\n");
  exit (EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
