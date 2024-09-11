#ifndef __CS_PARAM_SADDLE_H__
#define __CS_PARAM_SADDLE_H__

/*============================================================================
 * Routines to handle the SLES settings
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_param_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 * \file cs_param_saddle.h

 * \brief Handle the settings of saddle-point systems.
 *        These systems arise from the monolithic coupling of the Navier-Stokes
 *        equations or in mixed formulation of scalar-valued equations.
 */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_param_saddle_precond_t
 *
 * \brief Type of preconditioner used to solve a saddle-point system. Up to
 * now, this happens only in two cases: (1) with CDO cell-based schemes and (2)
 * with the Stokes or Navier-Stokes equations with CDO face-based schemes and a
 * monolithic approach (fully coupled) velocity-pressure coupling.
 *
 * \var CS_PARAM_SADDLE_PRECOND_NONE
 * No preconditioner to apply.
 *
 * \var CS_PARAM_SADDLE_PRECOND_DIAG
 * A block-diagonal preconditioner is used. The (1,1)-block is an approximation
 * of the (1,1)-block in the saddle-point system based on a cheap
 * resolution. The parameter settings for this resolution relies on the
 * structure \ref cs_param_sles_t. The (2,2)-block is approximated by the Schur
 * approximation given (by default the identity).
 *
 * \var CS_PARAM_SADDLE_PRECOND_LOWER
 * A 2x2 block matrix is used as preconditioner with the (1,2)-block fills with
 * zero. The (1,1)-block is an approximation of the (1,1)-block in the
 * saddle-point system based on a cheap resolution. The parameter settings for
 * this resolution relies on the structure \ref cs_param_sles_t. The
 * (2,2)-block is approximated by the Schur approximation given (by default the
 * identity).
 *
 * \var CS_PARAM_SADDLE_PRECOND_SGS
 * A symmetric Gauss-Seidel 2x2 block preconditioner is used. The (1,1)-block
 * is an approximation of the (1,1)-block in the saddle-point system based on a
 * cheap resolution. The parameter settings for this resolution relies on the
 * structure \ref cs_param_sles_t The (2,2)-block is approximated by the Schur
 * approximation given (by default the identity).
 *
 * \var CS_PARAM_SADDLE_PRECOND_UPPER
 * A 2x2 block matrix is used as preconditioner with the (2,1)-block fills with
 * zero. The (1,1)-block is an approximation of the (1,1)-block in the
 * saddle-point system based on a cheap resolution. The parameter settings for
 * this resolution relies on the structure \ref cs_param_sles_t The (2,2)-block
 * is approximated by the Schur approximation given (by default the identity).
 *
 * \var CS_PARAM_SADDLE_PRECOND_UZAWA
 * An Uzawa-like 2x2 block preconditioner. One needs a Schur approximation.
 */

typedef enum {

  CS_PARAM_SADDLE_PRECOND_NONE,
  CS_PARAM_SADDLE_PRECOND_DIAG,
  CS_PARAM_SADDLE_PRECOND_LOWER,
  CS_PARAM_SADDLE_PRECOND_SGS,
  CS_PARAM_SADDLE_PRECOND_UPPER,
  CS_PARAM_SADDLE_PRECOND_UZAWA,

  CS_PARAM_SADDLE_N_PRECOND

} cs_param_saddle_precond_t;

/*!
 * \enum cs_param_saddle_solver_t
 *
 * \brief Type of solver used to solve a saddle-point system. Up to now, this
 * happens only with CDO cell-based schemes or when solving the fully coupled
 * Navier-Stokes system (monolithic) with CDO face-based schemes. A
 * saddle-point system is an indefinite system.
 *
 * \var CS_PARAM_SADDLE_SOLVER_NONE
 * No solver defined. No saddle-point to solve.
 *
 * \var CS_PARAM_SADDLE_SOLVER_ALU
 * In-house solver based on the Uzawa algorithm with an Augmented Lagrangian
 * acceleration. In this case, one has to specify the scaling coefficient for
 * the augmented system (in the current case, a grad-div operator)
 *
 * \var CS_PARAM_SADDLE_SOLVER_FGMRES
 * Flexible variant of the GMRES Krylov iterative solver. This solver can be
 * applied to a general indefinite systems. Up to now, this solver can be
 * called only through the PETSc library.
 *
 * \var CS_PARAM_SADDLE_SOLVER_GCR
 * GCR = Genreralized Conjugate Residual. This iterative solver can be applied
 * to a general indefinite systems. It relies on a specific storage of the
 * saddle-point system: (1,1)-block is assembled and the (2,1)-block is
 * unassembled. The (1,2) is not stored since this is the transposition of the
 * (2,1)-block. This solver is an in-house version of the Notay's variant for
 * the GCR. One can specify the solver for the (1,1)-block preconditioner and
 * also the type of Schur approximation for the (2,2)-block.
 *
 * \var CS_PARAM_SADDLE_SOLVER_GKB
 * GKB = Golub-Kahan bi-diagonalization. This solver can be applied to a
 * symmetrice saddle-point system (or very slightly unsymmetric system). The
 * stopping crietrion is different from the other solvers since it relies on an
 * approximation of the residual in the energy norm (evaluation on-the-fly).
 * This solver is available throught PETSc or directly with
 * code_saturne. Moreover, one can add an augmentation term as in the ALU
 * algorithm.
 *
 * \var CS_PARAM_SADDLE_SOLVER_MINRES
 * MinRes = Minimal Residual algorithm. This iterative solver can be applied to
 * an indefinite symmetric systems. This is like a CG solver for symmetric
 * saddle-point systems. This solver relies on a specific storage of the
 * saddle-point system: (1,1)-block is assembled and the (2,1)-block is
 * unassembled. The (1,2) is not stored since this the transposition of the
 * (2,1)-block. One can specify the solver for the (1,1)-block preconditioner
 * and also the type of Schur approximation for the (2,2)-block.
 *
 * \var CS_PARAM_SADDLE_SOLVER_MUMPS
 * Sparse direct solver based on the external library called MUMPS. This
 * iterative solver can be applied to a general indefinite systems. Be careful
 * that according to the boundary conditions. The pressure can be defined up to
 * a constant. In this situation, MUMPS can breakdown since there is a non
 * empty kernel.
 *
 * \var CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM
 * Transform the saddle-point problem into an equivalent system without a zero
 * block for the (2,2) block. This transformation allows one to consider
 * standard preconditionner and iterative Krylov solver. This is an
 * experimental feature (only tested in sequential run)
 *
 * \var CS_PARAM_SADDLE_SOLVER_UZAWA_CG
 * Resolution using an uzawa algorithm optimized using a conjugate gradient
 * reformulation. Two systems are solved at each iteration (one related to the
 * velocity block, one related to the Schur complement approximation - size of
 * the pressure space). The first one is specified using the cs_sles_param_t
 * structure related to the (1,1)-block and the second system is specified
 * using the cs_sles_param_t structure related to the Schur approximation.
 */

typedef enum {

  CS_PARAM_SADDLE_SOLVER_NONE,

  CS_PARAM_SADDLE_SOLVER_ALU,
  CS_PARAM_SADDLE_SOLVER_FGMRES,
  CS_PARAM_SADDLE_SOLVER_GCR,
  CS_PARAM_SADDLE_SOLVER_GKB,
  CS_PARAM_SADDLE_SOLVER_MINRES,
  CS_PARAM_SADDLE_SOLVER_MUMPS,
  CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM,
  CS_PARAM_SADDLE_SOLVER_SIMPLE,
  CS_PARAM_SADDLE_SOLVER_UZAWA_CG,

  CS_PARAM_SADDLE_N_SOLVERS

} cs_param_saddle_solver_t;

/*! \enum cs_param_saddle_schur_approx_t
 *
 *  \brief Strategy to build the Schur complement approximation. This appears
 *         in block preconditioning or Uzawa algorithms when a fully coupled
 *         (also called monolithic) approach) is used.
 *  \f[ \begin{bmatrix} A&B^t\\ B&O  \end{bmatrix}\f]
 *  The exact Schur complement is then
 * \f[ S = -B \cdot A^{-1} \cdot B \f]
 *
 *  \var CS_PARAM_SADDLE_SCHUR_NONE
 *  There is no Schur complement approximation.
 *
 *  \var CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE
 *  The Schur complement approximation is defined as
 *  \f[ S \approx -B \cdot diag(A)^{-1} \cdot B^t \f]
 *
 *  \var CS_PARAM_SADDLE_SCHUR_IDENTITY
 *  The Schur complement approximation is simply the identity matrix
 *
 *  \var CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE
 *  The Schur complement approximation is defined as
 * \f[ B \cdot lumped(A^{-1}) \cdot B^t \f]
 *  where \f$x=lumped(A^{-1})\f$ results from \f$A.x = \bf{1}\f$ (\f$\bf{1}\f$
 *  is the array fills with 1 in each entry)
 *
 *  \var CS_PARAM_SADDLE_SCHUR_MASS_SCALED
 *  The Schur complement approximation is simply a scaled **diagonal** mass
 *  matrix related to the (2,2)-block
 *
 *  \var CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE
 *  The Schur complement approximation is defined as
 *  \f[ S \approx \alpha M_{22} + \frac{1}{dt} B.diag(A)^{-1}.B^t \f]
 *  where \f$M_{22}\f$ is the mass matrix related to the (2,2)-block
 *
 *  \var CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE
 *  The Schur complement approximation is defined as
 *  \f[ S \approx \alpha \cdot M_{22} + \frac{1}{dt} B\cdot lumped(A^{-1})
 *  \cdot B^t \f]
 *  where \f$ M_{22} \f$ is the mass matrix related to the (2,2) block and where
 *  \f$x=lumped(A^{-1})\f$ results from \f$A.x = \bf{1}\f$ (\f$\bf{1}\f$
 *  is the array fills with 1 in each entry)
 */

typedef enum {

  CS_PARAM_SADDLE_SCHUR_NONE,

  CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE,
  CS_PARAM_SADDLE_SCHUR_IDENTITY,
  CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE,
  CS_PARAM_SADDLE_SCHUR_MASS_SCALED,
  CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE,
  CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE,

  CS_PARAM_SADDLE_N_SCHUR_APPROX

} cs_param_saddle_schur_approx_t;

/*!
 * \struct cs_param_saddle_t
 * \brief Structure storing all metadata related to the resolution of a
 *        saddle-point linear system. A saddle-point system is depicted as
 *
 * |  A  :  Bt |  where A is the (1,1) block
 * |-----:-----|        Bt is the (1,2) block which the transposed operator
 * |  B  :  0  |        w.r.t. B, the (2,1)-block
 *
 *
 * According t the choice of the solver, the saddle-point system is totally or
 * partially assembled.
 */

typedef struct {

  /*!
    * \var verbosity
    * verbosity (level of information displayed)
    */

  int                         verbosity;

  /*!
    * \var name
    * name of the saddle system or NULL (not mandatory). If NULL, the name
    * of the (1,1)-block is used.
    */

  char                       *name;

  /*!
   * @name Main parameter settings
   *
   * Set of parameters to drive the resolution of a saddle-point system
   */

  /*! @{ */

   /*!
    * \var solver_class
    * class of SLES to consider (advanced usage)
    */

  cs_param_solver_class_t     solver_class;

  /*!
   * \var solver
   *  Type of solver to solve the saddle-point system
   *  If solver is set to CS_PARAM_SADDLE_N_SOLVERS, then there is no need to
   *  solve a saddle-point system.
   */

  cs_param_saddle_solver_t    solver;

  /*! \var precond
   *  Type of preconditioner for the saddle-point system which is viewed as a
   *  2x2 block matrix.
   */

  cs_param_saddle_precond_t   precond;

  /*! \var cvg_param
   *  Structure storing the parameters to know if an iterative process has to
   *  stop (convergence or divergence). These criteria are related the
   *  iterative algorithm used to solve the saddle-point system. This is the
   *  case for instance with an Uzawa or GKB algorithm.
   */

  cs_param_convergence_t      cvg_param;

  /*! @} */

  /*! \var block11_sles_param
   * Set of parameters used to solve (1,1)-block i.e. the A matrix. This
   * is shared with the \ref cs_equation_param_t structure.
   */

  const cs_param_sles_t      *block11_sles_param;

  /*!
   * @name Schur complement approximation
   *
   * Set of parameters to drive the resolution of the pressure-related
   * block. This is often a Schur complement approximation to B.A^-1.Bt
   */
  /*! @{ */

  /*!
   * \var schur_approx
   *
   * \var schur_approx
   * Choice of the way of preconditioning the schur approximation
   */

  cs_param_saddle_schur_approx_t  schur_approx;

  /*! \var schur_sles_param
   * Set of parameters used to solve the Schur complement if needed. This
   * depends on the type of Schur approximation which has been chosen.
   */

  cs_param_sles_t                *schur_sles_param;

  /*! @} */

  /* Additional parameters which are advanced/developper settings */

  void                           *context;

} cs_param_saddle_t;


/* Set of advanced settings according to the type of saddle-point solver */
/* ===================================================================== */

/* Augmented Lagrangian Uzawa algorithm */
/* ------------------------------------ */

typedef struct {

  /*! \var augmentation_scaling
   *  Value of the scaling coefficient in front of the augmented system.
   *  This is only useful when a ALU algorithm is used.
   */

  double            augmentation_scaling;

  /* \var dedicated_init_sles
   * Define an additional SLES to perform the resolution associated to the
   * transformation of the right-hand side. By default, this is false in order
   * to not compute two setup steps in one call. */

  bool              dedicated_init_sles;

  /*! \var init_sles_param
   * The initial transformation of the linear system requires a more accurate
   * resolution. Thus, one adds a dedicated \ref cs_param_sles_t structure for
   * this purpose)
   */

  cs_param_sles_t  *init_sles_param;

} cs_param_saddle_context_alu_t;


/* Block preconditioner algorithm with a Krylov solver */
/* --------------------------------------------------- */

typedef struct {

  /*! \var n_stored_directions
   *  Number of iterations to perform before restarting the solver. This
   *  quantity is useful when a GCR or FMGRES is used.
   */

  int               n_stored_directions;

  /*! \var xtra_sles_param
   * Set of parameters only used in some situations such as the need to solve
   * approximately the A.x = 1 linear system (A is the (1,1)-block
   * matrix). This is a complementary step in the approximation of the Schur
   * complement.
   *
   * By default, this is a copy of the \ref block11_sles_param with less
   * restrictive convergence criteria
   */

  cs_param_sles_t  *xtra_sles_param;

} cs_param_saddle_context_block_krylov_t;


/* GKB algorithm */
/* ------------- */

typedef struct {

  /*! \var augmentation_scaling
   *  Value of the scaling coefficient in front of the augmented system.
   *  This is only useful when a GKB algorithm is used. By default, there is
   *  no augmentation.
   */

  double            augmentation_scaling;

  /*! \var truncation_threshold
   *  Default number of values used to estimation the residual in the energy
   *  norm.
   */

  int               truncation_threshold;

  /* \var dedicated_init_sles
   * Define an additional SLES to perform the resolution associated to the
   * transformation of the right-hand side. By default, this is false in order
   * to not compute two setup steps in one call. */

  bool              dedicated_init_sles;

  /*! \var init_sles_param
   * The initial transformation of the linear system requires a more accurate
   * resolution. Thus, one adds a dedicated \ref cs_param_sles_t structure for
   * this purpose)
   */

  cs_param_sles_t  *init_sles_param;

} cs_param_saddle_context_gkb_t;


/* Notay's algebraic transformation */
/* -------------------------------- */

typedef struct {

  /* \var scaling_coef
   * This coefficient is used in Notay's transformation devised in "Algebraic
   * multigrid for Stokes equations" SIAM J. Sci. Comput. Vol. 39 (5), 2017
   */

  double  scaling_coef;

} cs_param_saddle_context_notay_t;

/* Uzawa-CG algorithm */
/* ------------------ */

typedef struct {

  /* \var dedicated_init_sles
   * Define an additional SLES to perform the initial resolution. By default,
   * this is false in order to not compute two setup steps in one call. */

  bool              dedicated_init_sles;

  /*! \var init_sles_param
   * The initial linear system requires a more accurate resolution. Thus, one
   * adds a dedicated \ref cs_param_sles_t structure for this purpose)
   */

  cs_param_sles_t  *init_sles_param;

  /*! \var xtra_sles_param
   * Set of parameters only used in some situations such as the need to solve
   * approximately the A.x = 1 linear system (A is the (1,1)-block
   * matrix). This is a complementary step in the approximation of the Schur
   * complement.
   *
   * By default, this is a copy of the \ref block11_sles_param with less
   * restrictive convergence criteria
   */

  cs_param_sles_t  *xtra_sles_param;

} cs_param_saddle_context_uzacg_t;

/* SIMPLE-like algorithm */
/* --------------------- */

typedef struct {

  /* \var dedicated_init_sles
   * Define an additional SLES to perform the initial resolution. By default,
   * this is false in order to not compute two setup steps in one call. */

  bool              dedicated_init_sles;

  /*! \var init_sles_param
   * The initial linear system requires a more accurate resolution. Thus, one
   * adds a dedicated \ref cs_param_sles_t structure for this purpose)
   */

  cs_param_sles_t  *init_sles_param;

  /*! \var xtra_sles_param
   * Set of parameters only used in some situations such as the need to solve
   * approximately the A.x = 1 linear system (A is the (1,1)-block
   * matrix). This is a complementary step in the approximation of the Schur
   * complement.
   *
   * By default, this is a copy of the \ref block11_sles_param with less
   * restrictive convergence criteria
   */

  cs_param_sles_t  *xtra_sles_param;

} cs_param_saddle_context_simple_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the number of iterations to store before starting a Krylov solver
 *
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 * \param[in]      restart_range  number of directions
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_restart_range(cs_param_saddle_t  *saddlep,
                                  int                 restart_range);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling coefficient used in the Notay's transformation
 *        devised in
 *        "Algebraic multigrid for Stokes equations" SIAM J. Sci. Comput.
 *        Vol. 39 (5), 2017
 *        In this article, this scaling is denoted by alpha
 *
 * \param[in, out] saddlep       set of parameters for solving a saddle-point
 * \param[in]      scaling_coef  value of the scaling coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_notay_scaling(cs_param_saddle_t  *saddlep,
                                  double              scaling_coef);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling in front of the augmentation term when an ALU or a
 *        GKB algorithm is considered
 *
 * \param[in, out] saddlep  set of parameters for solving a saddle-point
 * \param[in]      coef     value of the scaling coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_augmentation_coef(cs_param_saddle_t  *saddlep,
                                      double              coef);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the scaling coefficient in front of the augmentation term when an
 *        ALU or a GKB algorithm is considered.
 *
 * \param[in] saddlep  set of parameters for solving a saddle-point
 *
 * \return 0 if not relevant or the value of the augmentation coefficient
 */
/*----------------------------------------------------------------------------*/

double
cs_param_saddle_get_augmentation_coef(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the name of the type of saddle-point solver
 *
 * \param[in] type  type of saddle-point solver
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_saddle_get_type_name(cs_param_saddle_solver_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a cs_param_saddle_t structure
 *        No solver is set by default.
 *
 * \return a pointer to the new cs_param_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_saddle_t *
cs_param_saddle_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_free(cs_param_saddle_t  **p_saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the name of the saddle-point solver
 *
 * \param[in] saddlep  pointer to a set of saddle-point parameters
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_saddle_get_name(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the name of the saddle-point system.
 *
 * \param[in]      name     name associated to this saddle-point system
 * \param[in, out] saddlep  pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_name(const char         *name,
                         cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the \ref cs_param_sles_t structure (shared) related to the
 *        (1,1)-block to the structure managing the resolution of the
 *        saddle-point problems.
 *
 * \param[in, out] saddlep         pointer to the structure to update
 * \param[in]      block11_slesp   set of parameters for the (1,1) block
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_block11_sles_param(cs_param_saddle_t      *saddlep,
                                       const cs_param_sles_t  *block11_slesp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of preconditioning to apply for this saddle-point system
 *
 * \param[in]      keyval     value of the key for the preconditioner
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_precond(const char          *keyval,
                            cs_param_saddle_t   *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of Schur approximation to apply to this saddle-point
 *        system
 *
 * \param[in]      keyval     value of the key for the schur approx.
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_schur_approx(const char          *keyval,
                                 cs_param_saddle_t   *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the class of solver to apply for this saddle-point system
 *
 * \param[in]      keyval     value of the key for the preconditioner
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_solver_class(const char          *keyval,
                                 cs_param_saddle_t   *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of solver to apply for this saddle-point system
 *
 * \param[in]      keyval     value of the key for the preconditioner
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_solver(const char          *keyval,
                           cs_param_saddle_t   *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a \ref cs_param_sles_t structure for the Schur
 *        approximation nested inside a \ref cs_param_saddle_t structure. By
 *        default, this member is not allocated. Do nothing if the related
 *        structure is already allocated.
 *
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_try_init_schur_sles_param(cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to the approximation of the Schur complement. This is
 *        only useful for solving a saddle-point problem relying on an
 *        elaborated approximation of the Schur complement.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_schur_sles_param(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to an extra-operation specific to a saddle-point solver
 *        It returns a non NULL pointer only for some sadlle-point solver
 *        relying on a more elaborated Schur complement approximation.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_xtra_sles_param(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to the initial saddle-point problem. It returns a non
 *        NULL pointer only for some sadlle-point solver.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_init_sles_param(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_param_saddle_t structure from ref to dest
 *
 * \param[in]      ref     reference structure to be copied
 * \param[in, out] dest    destination structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_copy(const cs_param_saddle_t  *ref,
                     cs_param_saddle_t        *dest);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup information for the given cs_param_saddle_t structure
 *
 * \param[in] saddlep     pointer to the structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_log(const cs_param_saddle_t  *saddlep);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PARAM_SADDLE_H__ */
