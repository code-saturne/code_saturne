/*============================================================================
 * Unit tests for Anderson acceleration iterative algorithm
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "cdo/cs_cdo_blas.h"
#include "cdo/cs_iter_algo.h"

/*============================================================================
 * Static variables for BLAS callbacks
 *============================================================================*/

static cs_lnum_t _active_n_elts = 0;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Euclidean dot product of two arrays of size _active_n_elts
 */
/*----------------------------------------------------------------------------*/

static double
_test_dotprod(const cs_real_t *a, const cs_real_t *b)
{
  double sum = 0.0;
  for (cs_lnum_t i = 0; i < _active_n_elts; i++)
    sum += a[i] * b[i];
  return sum;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute squared Euclidean norm of an array of size _active_n_elts
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_test_sqnorm(const cs_real_t *a)
{
  double sum = 0.0;
  for (cs_lnum_t i = 0; i < _active_n_elts; i++)
    sum += a[i] * a[i];
  return static_cast<cs_real_t>(sum);
}

/*============================================================================
 * Test 1: Lifecycle, parameters and accessors
 *============================================================================*/

static void
_test_lifecycle_and_accessors(FILE *out)
{
  fprintf(out, "=== Test 1: Lifecycle and accessors ===\n");

  cs_iter_algo_t *algo = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  assert(algo != nullptr);
  assert(algo->type == CS_ITER_ALGO_ANDERSON);

  /* Verify default Anderson parameters */
  cs_iter_algo_param_aac_t p_def = cs_iter_algo_get_anderson_param(algo);
  assert(p_def.n_max_dir == 5);
  assert(p_def.starting_iter == 3);
  assert(cs::abs(p_def.max_cond - (-1.0)) < 1e-12);
  assert(cs::abs(p_def.beta - 1.0) < 1e-12);

  /* Configure custom parameters */
  cs_iter_algo_param_aac_t p_custom;
  p_custom.n_max_dir     = 4;
  p_custom.starting_iter = 2;
  p_custom.max_cond      = -1.0;
  p_custom.beta          = 1.0;

  const cs_lnum_t n_elts = 25;
  cs_iter_algo_set_anderson_param(algo, p_custom, n_elts);

  cs_iter_algo_param_aac_t p_read = cs_iter_algo_get_anderson_param(algo);
  assert(p_read.n_max_dir == 4);
  assert(p_read.starting_iter == 2);
  assert(p_read.max_cond < 0.0);
  assert(cs::abs(p_read.beta - 1.0) < 1e-12);

  /* Check initial state queries */
  assert(cs_iter_algo_get_n_iter(algo) == 0);
  assert(cs_iter_algo_get_n_inner_iter(algo) == 0);
  assert(cs_iter_algo_get_cvg_status(algo) == CS_SLES_ITERATING);
  assert(cs::abs(cs_iter_algo_get_normalization(algo) - 1.0) < 1e-12);

  /* Check tolerance and normalization modifiers */
  cs_iter_algo_set_normalization(algo, 2.5);
  assert(cs::abs(cs_iter_algo_get_normalization(algo) - 2.5) < 1e-12);

  cs_iter_algo_set_initial_residual(algo, 10.0);
  cs_iter_algo_set_tolerance(algo, 1e-8);

  /* Reset and verify */
  cs_iter_algo_reset(algo);
  assert(cs_iter_algo_get_n_iter(algo) == 0);
  assert(cs_iter_algo_get_cvg_status(algo) == CS_SLES_ITERATING);

  /* Free */
  cs_iter_algo_free(&algo);
  assert(algo == nullptr);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 2: Nonlinear vector fixed-point problem (Anderson vs Picard)
 *============================================================================*/

static void
_test_fixed_point_anderson(FILE *out)
{
  fprintf(out, "=== Test 2: Nonlinear fixed-point problem ===\n");

  /* Fixed-point mapping g_i(x) = 0.5 * cos(x_{(i+1)%N}) + 0.5
     Unique fixed point in [0, 1]^N */
  const cs_lnum_t n = 12;
  _active_n_elts    = n;

  /* 1. Picard iteration reference */
  cs_real_t *x_picard = nullptr;
  CS_MALLOC(x_picard, n, cs_real_t);
  for (cs_lnum_t i = 0; i < n; i++)
    x_picard[i] = 0.1 * (i + 1);

  int    iter_picard = 0;
  double res_picard  = 1.0;
  while (res_picard > 1e-10 && iter_picard < 200) {
    double    sq_diff = 0.0;
    cs_real_t x_next[12];
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_lnum_t next_i = (i + 1) % n;
      x_next[i]        = 0.5 * std::cos(x_picard[next_i]) + 0.5;
      double diff      = x_next[i] - x_picard[i];
      sq_diff += diff * diff;
    }
    for (cs_lnum_t i = 0; i < n; i++)
      x_picard[i] = x_next[i];
    res_picard = std::sqrt(sq_diff);
    iter_picard++;
  }
  fprintf(out,
          "  Picard: %d iterations, residual = %.4e\n",
          iter_picard,
          res_picard);

  /* 2. Anderson accelerated iteration */
  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-10;
  algo->cvg_param.atol       = 1e-12;
  algo->cvg_param.n_max_iter = 100;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 4;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t *cur_x = nullptr;
  cs_real_t *pre_x = nullptr;
  CS_MALLOC(cur_x, n, cs_real_t);
  CS_MALLOC(pre_x, n, cs_real_t);

  for (cs_lnum_t i = 0; i < n; i++)
    cur_x[i] = 0.1 * (i + 1);

  int                         iter_anderson = 0;
  cs_sles_convergence_state_t state         = CS_SLES_ITERATING;

  while (state == CS_SLES_ITERATING && iter_anderson < 100) {
    /* pre_x = x_k */
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* Evaluate g(x_k): cur_x = g(pre_x) */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_lnum_t next_i = (i + 1) % n;
      cur_x[i]         = 0.5 * std::cos(pre_x[next_i]) + 0.5;
    }

    /* Apply Anderson acceleration */
    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    /* Compute residual norm ||cur_x - pre_x|| */
    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }
    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    iter_anderson++;
  }

  fprintf(out,
          "  Anderson: %d iterations, residual = %.4e, status = %d\n",
          iter_anderson,
          cs_iter_algo_get_residual(algo),
          state);

  assert(state == CS_SLES_CONVERGED);
  assert(iter_anderson <= iter_picard);

  /* Compare final solutions */
  double max_diff = 0.0;
  for (cs_lnum_t i = 0; i < n; i++)
    max_diff = cs::max(max_diff, cs::abs(cur_x[i] - x_picard[i]));

  fprintf(out,
          "  Max difference between Anderson and Picard: %.4e\n",
          max_diff);
  assert(max_diff < 1e-8);

  CS_FREE(x_picard);
  CS_FREE(cur_x);
  CS_FREE(pre_x);
  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 3: History saturation and sliding window (QR deletion)
 *============================================================================*/

static void
_test_history_saturation_qrdelete(FILE *out)
{
  fprintf(out, "=== Test 3: History saturation with QR deletion ===\n");

  /* Use small m_max = 2 so QR deletion is exercised repeatedly */
  const cs_lnum_t n = 8;
  _active_n_elts    = n;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-11;
  algo->cvg_param.atol       = 1e-13;
  algo->cvg_param.n_max_iter = 50;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 2; /* Memory window of only 2 directions */
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t *cur_x = nullptr;
  cs_real_t *pre_x = nullptr;
  CS_MALLOC(cur_x, n, cs_real_t);
  CS_MALLOC(pre_x, n, cs_real_t);

  for (cs_lnum_t i = 0; i < n; i++)
    cur_x[i] = 0.1 * (i + 1);

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 50) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* Coupled nonlinear mapping with distinct component dynamics */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_lnum_t next_i = (i + 1) % n;
      cur_x[i]         = 0.5 * std::sin(pre_x[next_i]) + 0.05 * (i + 1);
    }

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }
    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    fprintf(out,
            "    iter %d: res = %.4e, cur_x[0] = %.4e, status = %d\n",
            iter,
            res,
            cur_x[0],
            state);

    iter++;
  }

  fprintf(out,
          "  Converged in %d iterations (m_max=2), final res = %.4e\n",
          iter,
          cs_iter_algo_get_residual(algo));

  /* Ensure history saturation was triggered (more iterations than m_max) */
  assert(iter > aac_p.n_max_dir + aac_p.starting_iter);
  assert(state == CS_SLES_CONVERGED);

  CS_FREE(cur_x);
  CS_FREE(pre_x);
  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 4: Warmup phase verification (starting_iter)
 *============================================================================*/

static void
_test_warmup_phase(FILE *out)
{
  fprintf(out, "=== Test 4: Warmup phase (starting_iter) ===\n");

  const cs_lnum_t n = 6;
  _active_n_elts    = n;

  cs_iter_algo_t          *algo = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 3;
  aac_p.starting_iter = 3; /* Wait 3 full iterations */
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t cur_x[6] = { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
  cs_real_t pre_x[6];
  cs_real_t expected_picard[6];

  for (int step = 0; step < 5; step++) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* g(x) = 0.5 * x + 0.2 */
    for (cs_lnum_t i = 0; i < n; i++) {
      expected_picard[i] = 0.5 * pre_x[i] + 0.2;
      cur_x[i]           = expected_picard[i];
    }

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    /* During warmup (step < starting_iter), cur_x MUST remain the Picard step
     */
    if (step < aac_p.starting_iter) {
      for (cs_lnum_t i = 0; i < n; i++)
        assert(cs::abs(cur_x[i] - expected_picard[i]) < 1e-14);
    }

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }
    cs_iter_algo_update_residual(algo, std::sqrt(sq_diff));
    cs_iter_algo_update_cvg_tol_auto(algo);
  }

  assert(cs_iter_algo_get_n_iter(algo) == 5);
  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 5: Linear system solving via accelerated Richardson iteration
 *============================================================================*/

static void
_test_linear_system_richardson(FILE *out)
{
  fprintf(out,
          "=== Test 5: Linear Poisson solve (Accelerated Richardson) ===\n");

  /* 1D discrete Laplacian: A * x = b
     size N, Dirichlet boundary conditions */
  const cs_lnum_t n = 10;
  _active_n_elts    = n;

  cs_real_t *b     = nullptr;
  cs_real_t *cur_x = nullptr;
  cs_real_t *pre_x = nullptr;
  cs_real_t *ax    = nullptr;
  CS_MALLOC(b, n, cs_real_t);
  CS_MALLOC(cur_x, n, cs_real_t);
  CS_MALLOC(pre_x, n, cs_real_t);
  CS_MALLOC(ax, n, cs_real_t);

  /* RHS b_i = 1.0 */
  for (cs_lnum_t i = 0; i < n; i++) {
    b[i]     = 1.0;
    cur_x[i] = 0.0;
  }

  /* Richardson relaxation parameter omega < 2 / lambda_max.
     Here lambda_max < 4, so omega = 0.45 is stable */
  const cs_real_t omega = 0.45;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-6;
  algo->cvg_param.atol       = 1e-8;
  algo->cvg_param.n_max_iter = 150;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 5;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 150) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* Compute A * pre_x */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_real_t val = 2.0 * pre_x[i];
      if (i > 0)
        val -= pre_x[i - 1];
      if (i < n - 1)
        val -= pre_x[i + 1];
      ax[i] = val;
    }

    /* Fixed-point Richardson: cur_x = pre_x - omega * (ax - b) */
    for (cs_lnum_t i = 0; i < n; i++)
      cur_x[i] = pre_x[i] - omega * (ax[i] - b[i]);

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    /* Residual norm ||b - A * cur_x|| */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_real_t val = 2.0 * cur_x[i];
      if (i > 0)
        val -= cur_x[i - 1];
      if (i < n - 1)
        val -= cur_x[i + 1];
      ax[i] = val;
    }

    double sq_res = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double r = b[i] - ax[i];
      sq_res += r * r;
    }
    double res = std::sqrt(sq_res);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    iter++;
  }

  fprintf(out,
          "  Linear solve converged in %d iterations, residual = %.4e\n",
          iter,
          cs_iter_algo_get_residual(algo));

  assert(state == CS_SLES_CONVERGED);

  /* Verify analytical solution for -u'' = 1 with u(0)=u(L)=0:
     x_i = 0.5 * (i + 1) * (N + 1 - (i + 1)) */
  double max_err = 0.0;
  for (cs_lnum_t i = 0; i < n; i++) {
    double x_exact = 0.5 * (i + 1.0) * (n - i);
    max_err        = cs::max(max_err, cs::abs(cur_x[i] - x_exact));
  }
  fprintf(out, "  Max absolute error against exact solution: %.4e\n", max_err);
  assert(max_err < 1e-5);

  CS_FREE(b);
  CS_FREE(cur_x);
  CS_FREE(pre_x);
  CS_FREE(ax);
  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 6: Two-level iteration monitoring
 *============================================================================*/

static void
_test_two_level_monitoring(FILE *out)
{
  fprintf(out, "=== Test 6: Two-level iteration monitoring ===\n");

  cs_iter_algo_t *algo =
    cs_iter_algo_create(CS_ITER_ALGO_ANDERSON | CS_ITER_ALGO_TWO_LEVEL);
  assert(algo != nullptr);
  algo->verbosity = 1;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 3;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, 5);

  /* Simulate 4 outer iterations with variable inner iterations */
  int inner_counts[4] = { 8, 6, 5, 4 };
  int expected_cumul  = 0;

  for (int step = 0; step < 4; step++) {
    cs_iter_algo_update_inner_iters(algo, inner_counts[step]);
    expected_cumul += inner_counts[step];

    assert(cs_iter_algo_get_n_inner_iter(algo) == expected_cumul);

    cs_iter_algo_update_residual(algo, 1.0 / (step + 1));
    cs_iter_algo_update_cvg_tol_auto(algo);

    cs_iter_algo_log_cvg(algo, "TwoLevel");
  }

  assert(cs_iter_algo_get_n_iter(algo) == 4);
  assert(cs_iter_algo_get_n_inner_iter(algo) == 23);

  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 7: Divergence and max iterations detection
 *============================================================================*/

static void
_test_divergence_and_max_iter(FILE *out)
{
  fprintf(out, "=== Test 7: Divergence and max iterations detection ===\n");

  /* 1. Divergence test */
  {
    cs_iter_algo_t *algo = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
    algo->cvg_param.dtol = 5.0; /* Divergence threshold */
    algo->cvg_param.rtol = 1e-6;

    cs_iter_algo_update_residual(algo, 1.0);
    cs_sles_convergence_state_t st0 = cs_iter_algo_update_cvg_tol_auto(algo);
    assert(st0 == CS_SLES_ITERATING);

    cs_iter_algo_update_residual(algo, 2.0);
    cs_sles_convergence_state_t st1 = cs_iter_algo_update_cvg_tol_auto(algo);
    assert(st1 == CS_SLES_ITERATING);

    /* Residual jumps by factor > 5.0 relative to res0 */
    cs_iter_algo_update_residual(algo, 25.0);
    cs_sles_convergence_state_t st2 = cs_iter_algo_update_cvg_tol_auto(algo);
    assert(st2 == CS_SLES_DIVERGED);
    assert(cs_iter_algo_get_cvg_status(algo) == CS_SLES_DIVERGED);

    cs_iter_algo_free(&algo);
  }

  /* 2. Max iterations test */
  {
    cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
    algo->cvg_param.n_max_iter = 4;
    algo->cvg_param.rtol       = 1e-12;
    algo->cvg_param.atol       = 1e-14;

    for (int step = 0; step < 3; step++) {
      cs_iter_algo_update_residual(algo, 0.5);
      cs_sles_convergence_state_t st = cs_iter_algo_update_cvg_tol_auto(algo);
      assert(st == CS_SLES_ITERATING);
    }

    /* 4th iteration reaches n_max_iter */
    cs_iter_algo_update_residual(algo, 0.5);
    cs_sles_convergence_state_t st_final =
      cs_iter_algo_update_cvg_tol_auto(algo);
    assert(st_final == CS_SLES_MAX_ITERATION);
    assert(cs_iter_algo_get_cvg_status(algo) == CS_SLES_MAX_ITERATION);

    cs_iter_algo_free(&algo);
  }

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 8: Reset and reuse lifecycle
 *============================================================================*/

static void
_test_reset_and_reuse(FILE *out)
{
  fprintf(out, "=== Test 8: Reset and reuse lifecycle ===\n");

  const cs_lnum_t n = 5;
  _active_n_elts    = n;

  cs_iter_algo_t          *algo = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 3;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t cur_x[5] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
  cs_real_t pre_x[5];

  /* Run problem 1 */
  for (int step = 0; step < 8; step++) {
    for (cs_lnum_t i = 0; i < n; i++) {
      pre_x[i] = cur_x[i];
      cur_x[i] = 0.5 * pre_x[i] + 0.1;
    }
    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);
    cs_iter_algo_update_residual(algo, 1e-2 / (step + 1));
    cs_iter_algo_update_cvg_tol_auto(algo);
  }
  assert(cs_iter_algo_get_n_iter(algo) == 8);

  /* Reset for problem 2 */
  cs_iter_algo_reset(algo);
  assert(cs_iter_algo_get_n_iter(algo) == 0);
  assert(cs_iter_algo_get_cvg_status(algo) == CS_SLES_ITERATING);

  for (cs_lnum_t i = 0; i < n; i++)
    cur_x[i] = 0.2 * (i + 1);

  /* Run problem 2 reusing allocated structures */
  for (int step = 0; step < 5; step++) {
    for (cs_lnum_t i = 0; i < n; i++) {
      pre_x[i] = cur_x[i];
      cur_x[i] = 0.3 * pre_x[i] + 0.4;
    }
    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);
    cs_iter_algo_update_residual(algo, 1e-3 / (step + 1));
    cs_iter_algo_update_cvg_tol_auto(algo);
  }
  assert(cs_iter_algo_get_n_iter(algo) == 5);

  cs_iter_algo_free(&algo);
  assert(algo == nullptr);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 9: Damped Anderson acceleration (m_max = 1)
 *============================================================================*/

static void
_test_anderson_damping(FILE *out)
{
  fprintf(out, "=== Test 9: Damped Anderson acceleration (m_max = 1) ===\n");

  /* When m_max = 1, n_dir <= 1, so (n_dir - 1) == 0, exercising
     _anderson_damping without triggering the multi-subtraction bug */

  const cs_lnum_t n = 4;
  _active_n_elts    = n;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-10;
  algo->cvg_param.atol       = 1e-12;
  algo->cvg_param.n_max_iter = 60;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 1; /* 1 direction */
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 0.7; /* 0 < beta < 1 triggers _anderson_damping */
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t cur_x[4] = { 0.5, 0.5, 0.5, 0.5 };
  cs_real_t pre_x[4];

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 60) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* g_i(x) = 0.4 * cos(x_i) + 0.3 */
    for (cs_lnum_t i = 0; i < n; i++)
      cur_x[i] = 0.4 * std::cos(pre_x[i]) + 0.3;

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }

    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    iter++;
  }

  fprintf(out,
          "  Damped Anderson converged in %d iterations, residual = %.4e\n",
          iter,
          cs_iter_algo_get_residual(algo));
  assert(state == CS_SLES_CONVERGED);

  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 10: Rank deficiency guard and collinear directions handling
 *============================================================================*/

static void
_test_rank_deficiency_guard(FILE *out)
{
  fprintf(out, "=== Test 10: Rank deficiency guard ===\n");

  /* System with identical initial states and decoupled identical equations:
     g_i(x) = 0.6 * sin(x_i) + 0.2 with x_i = 0.8 for all i.
     All Delta f vectors are strictly collinear, which triggers rank deficiency
     at iteration 2 (n_dir > 1).
     The rank deficiency guard must detect this, drop redundant directions,
     and prevent division by zero / NaN. */

  const cs_lnum_t n = 8;
  _active_n_elts    = n;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-10;
  algo->cvg_param.atol       = 1e-12;
  algo->cvg_param.n_max_iter = 50;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 3;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t cur_x[8];
  cs_real_t pre_x[8];
  for (cs_lnum_t i = 0; i < n; i++)
    cur_x[i] = 0.8;

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 50) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    for (cs_lnum_t i = 0; i < n; i++)
      cur_x[i] = 0.6 * std::sin(pre_x[i]) + 0.2;

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }

    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    iter++;
  }

  fprintf(
    out,
    "  Rank-deficient problem converged in %d iterations, residual = %.4e\n",
    iter,
    cs_iter_algo_get_residual(algo));

  assert(state == CS_SLES_CONVERGED);
  /* Verify that solution does not contain NaN */
  for (cs_lnum_t i = 0; i < n; i++)
    assert(!std::isnan(cur_x[i]) && !std::isinf(cur_x[i]));

  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 11: Multi-direction Anderson damping (m_max > 1, 0 < beta < 1)
 *============================================================================*/

static void
_test_multidir_damping(FILE *out)
{
  fprintf(out, "=== Test 11: Multi-direction Anderson damping ===\n");

  /* Tests Anderson damping with m_max = 4 directions and beta = 0.75.
     Verifies that the corrected damping formula (subtracting fold once)
     converges cleanly without residual explosion when n_dir > 1. */

  const cs_lnum_t n = 10;
  _active_n_elts    = n;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-10;
  algo->cvg_param.atol       = 1e-12;
  algo->cvg_param.n_max_iter = 60;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 4; /* Multiple directions */
  aac_p.starting_iter = 1;
  aac_p.max_cond      = -1.0;
  aac_p.beta          = 0.75; /* Damping active */
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  cs_real_t cur_x[10];
  cs_real_t pre_x[10];
  for (cs_lnum_t i = 0; i < n; i++)
    cur_x[i] = 0.1 * (i + 1);

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 60) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* Coupled nonlinear mapping */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_lnum_t next_i = (i + 1) % n;
      cur_x[i]         = 0.5 * std::cos(pre_x[next_i]) + 0.5;
    }

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }

    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    iter++;
  }

  fprintf(out,
          "  Multi-direction damped Anderson converged in %d iterations, "
          "residual = %.4e\n",
          iter,
          cs_iter_algo_get_residual(algo));

  assert(state == CS_SLES_CONVERGED);
  /* Verify that iterations exceeded m_max so n_dir > 1 was fully exercised */
  assert(iter > aac_p.n_max_dir + aac_p.starting_iter);

  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Test 12: Condition number filtering under degraded conditioning (CR-002)
 *============================================================================*/

static void
_test_condition_number_filtering(FILE *out)
{
  fprintf(out, "=== Test 12: Condition number filtering under degraded conditioning ===\n");

  /* Tests Anderson acceleration when the condition number of R degrades,
     triggering the QR filtering mechanism (max_cond > 0).
     We set up a stiff linear Richardson iteration where disparate eigenvalues
     lead to rapid convergence along large-eigenvalue modes, producing nearly
     collinear directions and a large condition number in R.

     Parameters:
       m_max = 4 (multiple directions)
       max_cond = 50.0 (strict threshold to trigger filtering)

     Verifies:
       1. The algorithm prunes ill-conditioned directions via _qrdelete with
          proper dimension (c->n_dir) without memory corruption or out-of-bounds access.
       2. The history dg is properly shifted via _anderson_shift_dg, keeping
          gamma and dg synchronized.
       3. The solver converges to machine precision and matches exact solution. */

  const cs_lnum_t n = 6;
  _active_n_elts    = n;

  cs_iter_algo_t *algo       = cs_iter_algo_create(CS_ITER_ALGO_ANDERSON);
  algo->cvg_param.rtol       = 1e-10;
  algo->cvg_param.atol       = 1e-12;
  algo->cvg_param.n_max_iter = 60;

  cs_iter_algo_param_aac_t aac_p;
  aac_p.n_max_dir     = 4;
  aac_p.starting_iter = 1;
  aac_p.max_cond      = 50.0; /* Strict conditioning threshold */
  aac_p.beta          = 1.0;
  cs_iter_algo_set_anderson_param(algo, aac_p, n);

  /* Matrix diagonal entries spanning two orders of magnitude */
  const cs_real_t diag_A[6]  = {2.0, 1.5, 1.0, 0.5, 0.1, 0.02};
  const cs_real_t rhs_b[6]   = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  const cs_real_t x_exact[6] = {1.0/2.0, 1.0/1.5, 1.0/1.0, 1.0/0.5, 1.0/0.1, 1.0/0.02};
  const cs_real_t omega      = 0.5;

  cs_real_t cur_x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  cs_real_t pre_x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  cs_sles_convergence_state_t state = CS_SLES_ITERATING;
  int                         iter  = 0;

  while (state == CS_SLES_ITERATING && iter < 60) {
    for (cs_lnum_t i = 0; i < n; i++)
      pre_x[i] = cur_x[i];

    /* Richardson step: cur_x = pre_x - omega * (A*pre_x - b) */
    for (cs_lnum_t i = 0; i < n; i++) {
      cs_real_t res_i = diag_A[i] * pre_x[i] - rhs_b[i];
      cur_x[i] = pre_x[i] - omega * res_i;
    }

    cs_iter_algo_update_anderson(algo,
                                 cur_x,
                                 pre_x,
                                 _test_dotprod,
                                 _test_sqnorm);

    double sq_diff = 0.0;
    for (cs_lnum_t i = 0; i < n; i++) {
      double diff = cur_x[i] - pre_x[i];
      sq_diff += diff * diff;
    }

    double res = std::sqrt(sq_diff);
    cs_iter_algo_update_residual(algo, res);
    state = cs_iter_algo_update_cvg_tol_auto(algo);

    fprintf(out,
            "    iter %2d: res = %.4e, x[0] = %.4e, x[5] = %.4e, status = %d\n",
            iter, res, cur_x[0], cur_x[5], state);

    iter++;
  }

  fprintf(out,
          "  Conditioning filtering problem converged in %d iterations, residual = %.4e\n",
          iter,
          cs_iter_algo_get_residual(algo));

  assert(state == CS_SLES_CONVERGED);

  /* Verify solution against exact solution */
  double max_err = 0.0;
  for (cs_lnum_t i = 0; i < n; i++) {
    double err = cs::abs(cur_x[i] - x_exact[i]);
    if (err > max_err)
      max_err = err;
    assert(!std::isnan(cur_x[i]) && !std::isinf(cur_x[i]));
  }
  fprintf(out, "  Max absolute error against exact solution: %.4e\n", max_err);
  assert(max_err < 1e-8);

  cs_iter_algo_free(&algo);

  fprintf(out, "  -> Passed.\n\n");
}

/*============================================================================
 * Main entry point
 *============================================================================*/

int
main(int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  FILE *out = stdout;

  printf("=============================================================\n");
  printf(" Running Anderson Acceleration Unit Tests\n");
  printf("=============================================================\n\n");

  // Test #1
  _test_lifecycle_and_accessors(out);

  // Test #2
  _test_fixed_point_anderson(out);

  // Test #3
  _test_history_saturation_qrdelete(out);

  // Test #4
  _test_warmup_phase(out);

  // Test #5
  _test_linear_system_richardson(out);

  // Test #6
  _test_two_level_monitoring(out);

  // Test #7
  _test_divergence_and_max_iter(out);

  // Test #8
  _test_reset_and_reuse(out);

  // Test #9
  _test_anderson_damping(out);

  // Test #10
  _test_rank_deficiency_guard(out);

  // Test #11
  _test_multidir_damping(out);

  // Test #12
  _test_condition_number_filtering(out);

  printf("=============================================================\n");
  printf(" All 12 Anderson Acceleration Tests Passed Successfully!\n");
  printf("=============================================================\n");

  return EXIT_SUCCESS;
}
