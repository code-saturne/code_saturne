#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for user-defined functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_base.h"
#include "base/cs_field.h"
#include "base/cs_mobile_structures.h"
#include "base/cs_probe.h"
#include "base/cs_time_control.h"
#include "base/cs_volume_zone.h"
#include "cdo/cs_domain.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_bad_cells.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  User function prototypes
 *============================================================================*/

/*! \file cs_user_1d_wall_thermal.cpp */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Data Entry of the 1D wall thermal module.
 *
 * \param[in]   iappel   Call number:
 *                       - 1: first call during the initialization (called once)
 *                       Setting the number of cells nfpt1d.
 *                       - 2: second call during the initialization (called once)
 *                       Filling ifpt1d, nppt1d, eppt1d and rgpt1d arrays.
 *                       - 3: called at each time step
 *                       Filling iclt1d, xlmbt1, rcpt1d and dtpt1d arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_1d_wall_thermal(int iappel);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_atmo.cpp
 *
 * \brief User-defined functions specific to amospheric flow models.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill in vertical profiles of atmospheric properties prior to solving
 *        1D radiative transfers.
 *
 * \param[in, out] preray        pressure vertical profile
 * \param[in, out] temray        real temperature vertical profile
 * \param[in, out] romray        density vertical profile
 * \param[in, out] qvray         water vapor content vertical profile
 * \param[in, out] qlray         water liquid content vertical profile
 * \param[in, out] ncray         droplets density vertical profile
 * \param[in, out] aeroso        aerosol concentration vertical profile
 */
/*----------------------------------------------------------------------------*/

void
cs_user_atmo_1d_rad_prf(cs_real_t   preray[],
                        cs_real_t   temray[],
                        cs_real_t   romray[],
                        cs_real_t   qvray[],
                        cs_real_t   qlray[],
                        cs_real_t   ncray[],
                        cs_real_t   aeroso[]);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.cpp
 *
 * \brief User functions for boundary condition definitions.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_setup(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 *
 * This function is called at each times step.
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", access
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]\n\n
 *
 * The icodcl values are common to all components, so icodcl values are
 * defined as for a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions for ALE.
 *
 * See \ref sec_ug_bc_ale_legacy for additional details.
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in, out]  bc_type      boundary face types
 * \param[in, out]  ale_bc_type  boundary face types for mesh velocity
 *                               (see \ref cs_boundary_ale_subtype_bits_t)
 * \param[in]       impale       indicator for prescribed node displacement
 *                               (0 or 1)
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", access
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]\n\n
 *
 * The icodcl values are common to all components, so icodcl values are
 * defined as for a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_ale(cs_domain_t  *domain,
                                int           bc_type[],
                                int           ale_bc_type[],
                                int           impale[]);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_coupling.cpp
 *
 * \brief Code couplings definition with SYRTHES and code_saturne.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with other instances of code_saturne.
 *
 * This is done by calling the \ref cs_sat_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_saturne_coupling(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with the SYRTHES code.
 *
 * This is done by calling the \ref cs_syr_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a volume exchange coefficient for SYRTHES couplings.
 *
 * \param[in]   coupling_id   Syrthes coupling id
 * \param[in]   syrthes_name  name of associated Syrthes instance
 * \param[in]   n_elts        number of associated cells
 * \param[in]   elt_ids       associated cell ids
 * \param[out]  h_vol         associated exchange coefficient (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling_volume_h(int               coupling_id,
                                  const char       *syrthes_name,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t   elt_ids[],
                                  cs_real_t         h_vol[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with CATHARE code.
 *
 * This is done by calling the \ref cs_sys_coupling_add function for each
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cathare_coupling(void);

/*! \file cs_user_electric_scaling.cpp */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rescale all electro-magnetic physical fields
 *        (electric potential, current density and Joule effect).
 *
 * \param[in]      mesh             pointer to a cs_mesh_t structure
 * \param[in,out]  mesh_quantities  pointer to a cs_mesh_quantities_t structure
 * \param[in]      dt               time step per cell
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_scaling_elec(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     cs_real_t                   *dt);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_electric_scaling.cpp
 *
 * \brief Define scaling parameter for electric model.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.cpp
 *
 * \brief The `cs_user_extra_operations` series of function have a very
 * broad purpose (i.e. handing anything that does not have another dedicated
 * user function).
 *
 * In most cases, `cs_user_extra_operations` is used for data extraction
 * and postprocessing, but may also be used to steer a computation,
 * for example.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables or setup extra operations.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_initialize(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User operations called at the end of each time step.
 *
 * This function has a very general purpose, although it is recommended to
 * handle mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User operations called at the end of the calculation.
 *
 * This function has a very general purpose, although it is recommended to
 * handle mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_finalize(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_fluid_structure_interaction.cpp
 *
 * \brief User-defined functions dedicated to Fluid-Structure interaction
 *        modeling.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of internal mobile structures and corresponding initial
 *        conditions (initial displacement and velocity ).
 *
 * \param[in]       is_restart         indicate if computation is restarted
 * \param[in]       n_structs          number of mobile structures
 * \param[in, out]  plot;              monitoring format mask
 *                                       0: no plot
 *                                       1: plot to text (.dat) format
 *                                       2: plot to .csv format
 *                                       3: plot to both formats
 * \param[in, out]  plot_time_control  plot time output frequency control
 * \param[in, out]  aexxst             coefficient for predicted displacement
 * \param[in, out]  bexxst             coefficient for predicted displacement
 * \param[in, out]  cfopre             coefficient for predicted force
 * \param[in, out]  xstr0              initial displacement per structure
 * \param[in, out]  vstr0              initial velocity per structure
 * \param[in, out]  xstreq             displacement of initial mesh relative to
 *                                     structures position at equilibrium
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_define(int                 is_restart,
                             int                 n_structs,
                             int                *plot,
                             cs_time_control_t  *plot_time_control,
                             cs_real_t          *aexxst,
                             cs_real_t          *bexxst,
                             cs_real_t          *cfopre,
                             cs_real_t           xstr0[][3],
                             cs_real_t           vstr0[][3],
                             cs_real_t           xstreq[][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Time-based settings for internal mobile structures.
 *
 * \param[in]       n_structs  number of mobile structures
 * \param[in]       ts         time step structure
 * \param[in]       xstreq     displacement of initial mesh rel. to equilibrium
 * \param[in]       xstr       structural displacement
 * \param[in]       vstr       structural velocity
 * \param[in, out]  xmstru     matrix of structural mass
 * \param[in, out]  xcstru     matrix of structural friction
 * \param[in, out]  xkstru     matrix of structural stiffness
 * \param[in, out]  forstr     forces acting on structures (take forces)
 * \param[in, out]  dtstr      structural time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_values(int                    n_structs,
                             const cs_time_step_t  *ts,
                             const cs_real_t        xstreq[][3],
                             const cs_real_t        xstr[][3],
                             const cs_real_t        vstr[][3],
                             cs_real_t              xmstru[][3][3],
                             cs_real_t              xcstru[][3][3],
                             cs_real_t              xkstru[][3][3],
                             cs_real_t              forstr[][3],
                             cs_real_t              dtstr[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define structure numbers for faces associated with internal
 *        or external (code_aster) structures.
 *
 * Structure numbers associated to a given face have the following values:
 * - -i where coupled to  i-th (1-to n) external (code_aster) structure.
 * - 0 where not coupled with an internal or external structure.
 * - i  where coupled to  i-th (1-to n) internal (mass-spring) structure.
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in, out]  structure_num  structure id associated to each face
 * \param[in, out]  structure_typ  structure type associated to each face
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_structure_num(cs_domain_t                *domain,
                          int                         structure_num[],
                          cs_mobile_structure_type_t  structure_typ[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute displacement fluid boundary for the external structures
 *        (code_aster excluded).
 *
 * \param[in]       domain    pointer to a cs_domain_t structure
 * \param[in]       b_stress  pointer to boundary stress array
 * \param[in, out]  disaple   pointer to mesh_displacement array
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_external_displacement(const cs_domain_t  *domain,
                                  const cs_real_3_t  *b_stress,
                                  cs_real_3_t        *disaple);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute convergence state of the external structures
 *        (code_aster excluded).
 *
 * Compute converge status and residual of the external coupling
 * Convergence status: 0 - external coupling has not converged
 *                     1 - external coupling has converged
 *
 * \param[in]   domain        pointer to a cs_domain_t structure
 * \param[in]   epsilon       convergence criterion
 * \param[out]  cvg_status    convergence status
 * \param[out]  residual      value of the residual
 */
/*----------------------------------------------------------------------------*/

void
cs_user_fsi_external_cvg(const cs_domain_t *domain,
                         const cs_real_t    epsilon,
                         int               *cvg_status,
                         cs_real_t         *residual);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_head_losses.cpp
 *
 * \brief User head loss definitions.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * ck11, ck22, ck33, ck12, ck23, ck13.
 *
 * Coefficients are set to zero (then computed based on definitions provided
 * through the GUI if this is the case) before calling this function, so
 * setting values to zero is usually not necessary, unless we want to fully
 * overwrite a GUI-based definition.
 *
 * Diagonal coefficients must be positive; the calculation may diverge
 * if this is not the case.
 *
 * \param[in]       zone  pointer to zone structure
 * \param[in, out]  cku   head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_user_head_losses(const cs_zone_t  *zone,
                    cs_real_t         cku[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_hgn.cpp
 *
 * \brief Define user properties for two-phase homogeneous compressible model.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the relaxation time-scale to equilibrium in the
 *        frame of the homogeneous two-phase model..
 *
 * \param[in]  mesh       pointer to mesh
 * \param[in]  alpha_eq   equilibrium volume fraction
 * \param[in]  y_eq       equilibrium mass fraction
 * \param[in]  z_eq       equilibrium energy fraction
 * \param[in]  ei         specific internal energy
 * \param[in]  v          specific volume
 * \param[in]  relax_tau  relaxation time scale towards equilibrium
 */
/*----------------------------------------------------------------------------*/

void
cs_user_hgn_thermo_relax_time(const cs_mesh_t  *mesh,
                              const cs_real_t  *alpha_eq,
                              const cs_real_t  *y_eq,
                              const cs_real_t  *z_eq,
                              const cs_real_t  *ei,
                              const cs_real_t  *v,
                              cs_real_t        *relax_tau);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.cpp
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define initial conditions for variables.
 *
 * This function is not yet called by code_saturne, use
 * `cs_user_initialization` for the time being
 * This function is called before the beginning of the computation
 * allowing an overload of the GUI defined initialization (called just
 * after cs_gui_initial_conditions).
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initial_conditions(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.cpp
 *
 * \brief Definition and modification of the calculation mesh.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Force preprocessing behavior in case of restart.
 *
 * By default, in case of restart, if a "restart/mesh_input.csm" file
 * is present, it will be read and proprocessing will be skipped.
 *
 * This behavior may be changed in the GUI (in the "Mesh" section, unchecking
 * "Use unmodified checkpoint mesh in case of restart"), or by calling
 * \ref cs_preprocessor_data_set_restart_mode in this function.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_restart_mode(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh files to read and optional associated transformations.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_input(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cartesian mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_cartesian_define(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh joinings.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_join(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define periodic faces.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_periodicity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set options for cutting of warped faces.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_warping(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundaries into a mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_boundary(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mesh smoothing.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_smoothe(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable or disable mesh saving.
 *
 * By default, mesh is saved when modified.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_save(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tag bad cells within the mesh based on user-defined geometric criteria.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_bad_cells_tag(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply partial modifications to the mesh after the preprocessing
 *        and initial postprocessing mesh building stage.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify_partial(cs_mesh_t             *mesh,
                            cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_paramedmem_coupling.cpp
 *
 * \brief User functions for input of ParaMEDMEM coupling parameters
 *
 * \brief User functions for input of calculation parameters.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define ParaMEDMEM coupling(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_couplings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled meshes
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_meshes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define fields to couple with ParaMEDMEM
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.cpp
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *
 * Available native iterative linear solvers include conjugate gradient,
 * Jacobi, BiCGStab, BiCGStab2, and GMRES. For symmetric linear systems,
 * an algebraic multigrid solver is available (and recommended).
 *
 * External solvers may also be setup using this function, the cs_sles_t
 * mechanism allowing such through user-define functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time moments.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined, and before fine control of field output options
 * is defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_moments(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define internal coupling options.
 *
 * Options are usually defined using cs_internal_coupling_add_entity.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_performance_tuning.cpp
 *
 * \brief Definition of advanced options relative to parallelism.
 *
 * See \ref cs_user_performance_tuning for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced mesh numbering options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_numbering(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define advanced partitioning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_partition(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define parallel IO settings.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parallel_io(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define sparse matrix tuning options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_matrix_tuning(void);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.cpp
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of enthalpy to temperature conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        enthalpy values
 * \param[in, out]  t        temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_h_to_t(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   h[],
                                   cs_real_t         t[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of temperature to enthalpy conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        temperature values
 * \param[in, out]  t        enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_t_to_h(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   t[],
                                   cs_real_t         h[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of the turbulence viscosity.
 *
 * Turbulent viscosity \f$ \mu_T \f$ (kg/(m s)) can be modified.
 * You can access the field by its name.
 *
 * \param[in, out]   domain      pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_turb_viscosity(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function to define a custom law for the thermodynamic pressure.
 *
 * Allows to define a custom law for the constant uniform thermodynamic
 * pressure (whenn \ref cs_velocity_pressure_model_t::idilat = 3 or
 * \ref cs_fluid_properties_t::ipthrm = 1).
 *
 * The density is then updated
 * (in \ref cs_compute_thermo_pressure_density.c) as:
 * \f[\rho^{n+1} =\rho^{n} \cdot \frac{P_{th}^{n+1}}{P_{th}^{n}}\f].
 *
 * \param[in, out]  td_p  Updated value of the thermodynamic pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_td_pressure(cs_real_t  *td_p);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_porosity.cpp
 *
 * \brief User definitions of porous media.
 *
 * See \ref cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the porosity (volume factor \f$ \epsilon \f$
 *        when the porosity model is activated.
 *        (\ref cs_glob_porous_model > 0).
 *
 * This function is called at the beginning of the simulation only.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_porosity(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_postprocess.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, null otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param[in]  nt_max_abs  maximum time step number
 * \param[in]  nt_cur_abs  current time step number
 * \param[in]  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs);

/*----------------------------------------------------------------------------*/
/*! \file cs_user_radiative_transfer.cpp
 *
 * \brief User function for input of radiative transfer parameters:
 *        absorption coefficient and net radiation flux.
 *
 *  See \ref cs_user_radiative_transfer for examples.
 */
/*----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Absorption coefficient for radiative module
 *
 * It is necessary to define the value of the fluid's absorption coefficient Ck.
 *
 * This value is defined automatically for specific physical models, such
 * as gas and coal combustion, so this function is not used by these models.
 *
 * For a transparent medium, the coefficient should be set to 0.
 *
 * In the case of the P-1 model, we check that the optical length is at
 * least of the order of 1.
 *
 * \param[in]   bc_type  boundary face types
 * \param[out]  ck       medium's absorption coefficient (zero if transparent)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_absorption(const int  bc_type[],
                                cs_real_t  ck[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the net radiation flux.
 *
 * The density of net radiation flux must be calculated
 * consistently with the boundary conditions of the intensity.
 * The density of net flux is the balance between the radiative
 * emiting part of a boundary face (and not the reflecting one)
 * and the radiative absorbing part.
 *
 * \param[in]   bc_type   boundary face types
 * \param[in]   twall     inside current wall temperature (K)
 * \param[in]   qincid    radiative incident flux  (W/m2)
 * \param[in]   xlam      conductivity (W/m/K)
 * \param[in]   epa       thickness (m)
 * \param[in]   eps       emissivity (>0)
 * \param[in]   ck        absorption coefficient
 * \param[out]  net_flux  net flux (W/m2)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_net_flux(const int        itypfb[],
                              const cs_real_t  twall[],
                              const cs_real_t  qincid[],
                              const cs_real_t  xlam[],
                              const cs_real_t  epa[],
                              const cs_real_t  eps[],
                              const cs_real_t  ck[],
                              cs_real_t        net_flux[]);

/*----------------------------------------------------------------------------*/
/*
 *! \file cs_user_radiative_transfer_bcs.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of radiative transfer boundary conditions.
 *
 * See \ref cs_user_radiative_transfer for examples.
 *
 * \warning the temperature unit here is the Kelvin
 *
 * For each boundary face face_id, a specific output (logging and
 * postprocessing) class id may be assigned. This allows realizing balance
 * sheets by treating them separately for each zone. By default, the
 * output class id is set to the general (input) zone id associated to a face.
 *
 * To access output class ids (both for reading and modifying), use the
 * \ref cs_boundary_zone_face_class_id function.
 * The zone id values are arbitrarily chosen by the user, but must be
 * positive integers; very high numbers may also lead to higher memory
 * consumption.
 *
 * \section cs_user_radiative_transfer_bcs_wall  Wall characteristics
 *
 * The following face characteristics must be set:
 *  - isothp(face_id) boundary face type
 *    * CS_BOUNDARY_RAD_WALL_GRAY:
 *      Gray wall with temperature based on fluid BCs
 *    * CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T:
 *      Gray wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T:
 *      Reflecting wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX:
 *      Gray wall with fixed conduction flux
 *    * CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX:
 *      Reflecting wall with fixed conduction flux
 *
 * Depending on the value of isothp, other values may also need to be set:
 *  - rcodcl = conduction flux
 *  - epsp   = emissivity
 *  - xlamp  = conductivity (W/m/K)
 *  - epap   = thickness (m)
 *  - textp  = outside temperature (K)
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       bc_type       boundary face types
 * \param[in]       isothp        boundary face type for radiative transfer
 * \param[out]      tmin          min allowed value of the wall temperature
 * \param[out]      tmax          max allowed value of the wall temperature
 * \param[in]       tx            relaxation coefficient (0 < tx < 1)
 * \param[in]       dt            time step (per cell)
 * \param[in]       thwall        inside current wall temperature (K)
 * \param[in]       qincid        radiative incident flux  (W/m2)
 * \param[in]       hfcnvp        convective exchange coefficient (W/m2/K)
 * \param[in]       flcnvp        convective flux (W/m2)
 * \param[out]      xlamp         conductivity (W/m/K)
 * \param[out]      epap          thickness (m)
 * \param[out]      epsp          emissivity (>0)
 * \param[out]      textp         outside temperature (K)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_radiative_transfer_bcs(cs_domain_t      *domain,
                               const int         bc_type[],
                               int               isothp[],
                               cs_real_t        *tmin,
                               cs_real_t        *tmax,
                               cs_real_t        *tx,
                               const cs_real_t   dt[],
                               const cs_real_t   thwall[],
                               const cs_real_t   qincid[],
                               cs_real_t         hfcnvp[],
                               cs_real_t         flcnvp[],
                               cs_real_t         xlamp[],
                               cs_real_t         epap[],
                               cs_real_t         epsp[],
                               cs_real_t         textp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_solver.cpp
 *
 * \brief User solver setting and implementation
 *
 * See \ref user_solver for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set user solver
 *
 * \return  1 if user solver is called, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_user_solver_set(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Main call to user solver
 *
 * \param[in] mesh pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_solver(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms.cpp
 *
 * \brief Additional source terms for variable equations.
 *
 * See \ref user_source_terms for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Additional user-defined source terms for variable equations
 *        (momentum, scalars, turbulence...).
 *
 * This function is called at each time step, for each relevant field.
 * It is therefore necessary to
 * test the value of the field id or name to separate
 * the treatments of the different variables.
 *
 * The additional source term is decomposed into an explicit part (st_exp) and
 * an implicit part (st_imp) that must be provided here.
 * The resulting equation solved by the code for a scalar f is:
 *
 *   \f[ \rho*volume*\frac{df}{dt} + .... = st\_imp*f + st\_exp \f]
 *
 * Note that st_exp and st_imp are defined after the Finite Volume integration
 * over the cells, so they include the "volume" term. More precisely:
 *   - st_exp is expressed in kg.[var]/s, where [var] is the unit of the
 *     variable.
 *     Its dimension is the one of the variable (3 for vectors)
 *   - st_imp is expressed in kg/s.
 *     Its dimension is 1 for scalars, 3x3 for vectors.
 *
 * The st_exp and st_imp arrays are already initialized to 0 (or a value
 * defined through the GUI or defined by a model) before entering
 * the function. It is generally not useful to reset them here.
 *
 * For stability reasons, code_saturne will not add -st_imp directly to the
 * diagonal of the matrix, but Max(-st_imp,0). This way, the st_imp term is
 * treated implicitely only if it strengthens the diagonal of the matrix.
 * However, when using the second-order in time scheme, this limitation cannot
 * be done anymore and -st_imp is added directly. The user should therefore
 * check for the negativity of st_imp.
 *
 * When using the second-order in time scheme, one should supply:
 *   - st_exp at time n
 *   - st_imp at time n+1/2
 *
 * \warning
 * \parblock
 *
 * If the variable is a temperature, the resulting equation solved is:
 *
 *   rho*Cp*volume*dT/dt + .... = st_imp*T + st_exp
 *
 * \endparblock
 *
 * Note that st_exp and st_imp are defined after the Finite Volume integration
 * over the cells, so they include the "volume" term. More precisely:
 *   - st_exp is expressed in W
 *   - st_imp is expressed in W/K
 *
 * \par Steep source terms
 * \parblock
 *
 * In case of a complex, non-linear source term, say F(f), for variable f, the
 * easiest method is to implement the source term explicitly.
 *
 *   df/dt = .... + F(f(n))
 *   where f(n) is the value of f at time tn, the beginning of the time step.
 *
 * This yields:
 *   st_exp = volume*F(f(n))
 *   st_imp = 0
 *
 * However, if the source term is potentially steep, this fully explicit
 * method will probably generate instabilities. It is therefore wiser to
 * partially implicit the term by writing:
 *
 *   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
 *
 * This yields:
 *   st_exp = volume*( F(f(n)) - dF/df*f(n) )
 *   st_imp = volume*dF/df
 *
 * \endparblock
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_time_table.cpp
 *
 * \brief User definitions of time tables
 *
 * See \ref cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define time tables using C API.
 * This function is called at the begin of the simulation only.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_time_table(void);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_turbomachinery.cpp
 *
 * \brief Definition of turbomachinery related options.
 *
 * See \ref turbomachinery for examples.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotor/stator model.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotor axes, associated cells, and rotor/stator faces.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_rotor(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotation velocity of rotor.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_set_rotation_velocity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_wall_condensation.cpp
 *
 * \brief Source terms associated at the boundary faces and the neighboring
 *      cells with surface condensation.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Source terms associated at the boundary faces and the neighboring
 *        cells with surface condensation.
 *
 * This function fills the condensation source terms for each variable at
 * the cell center associated to the boundary faces identifed in the mesh.
 * The fluid exchange coefficient is computed with a empirical law to be
 * imposed at the boundary face where the condensation phenomenon occurs.
 *
 * This user function is called which allows the setting of
 * \f$ \gamma_{\mbox{cond}} \f$ the condensation source term.
 *
 * This function fills the condensation source term array gamma_cond adding
 * to the following equations:
 *
 * - The equation for mass conservation:
 * \f[ D\frac{rho}{dt} + divs \left( \rho \vect{u}^n\right) = \Gamma _{cond}
 * \f]
 *
 * - The equation for a variable \f$\Phi \f$:
 * \f[ D\frac{\phi}{dt} = ... + \Gamma _{cond}*(\Phi _i - \Phi)
 * \f]
 *
 * discretized as below:
 *
 * \f[ \rho*\dfrac{\Phi^{n+1}-\Phi^{n}}/dt = ...
 *                            + \Gamma _{cond}*(\Phi _i - \Phi^{n+1})
 * \f]
 *
 * \remarks
 *  - \f$ \Phi _i \f$ is the value of \f$ \Phi \f$ associated to the
 *    injected condensation rate.
 *
 *    With 2 options are available:
 *       - the condensation rate is injected with the local value
 *         of variable \f$ \Phi = \Phi ^{n+1}\f$
 *         in this case the \f$ \Phi \f$ variable is not modified.
 *
 *       - the condensation rate is injected with a specific value
 *         for \f$ \Phi = \Phi _i \f$ the specified value given by the
 *         user.
 *
 * \section use Usage
 *
 * The three stages in the code where this User subroutine
 * is called (with \code iappel = 1, 2 and 3\endcode)
 *
 * \code iappel = 1 \endcode
 *  - Calculation of the number of cells where a mass source term is
 *    imposed: ncesmp
 *    Called once at the beginning of the calculation
 *
 * \code iappel = 2 \endcode
 *   - Identification of the cells where a mass source term is imposed:
 *     array icesmp(ncesmp)
 *     Called once at the beginning of the calculation
 *
 * \code iappel = 3 \endcode
 *   - Calculation of the values of the mass source term
 *     Called at each time step
 *
 * \section the specific variables to define with is user subroutine
 *
 *  - ifbpcd(ieltcd): identification of the faces where a condensation
 *                    source term is imposed.
 *
 *  - itypcd(ieltcd,ivar): type of treatment for variable ivar in the
 *                       ieltcd cell with condensation source term.
 *                     - itypcd = 0 -- * injection of ivar at local value
 *                     - itypcd = 1 -- * injection of ivar at user
 *                                      specified value.
 *
 *  - spcond(ielscd,ipr): value of the injection condensation rate
 *                       gamma_cond (kg/m3/s) in the ieltcd cell
 *                       with condensation source term.
 *
 *  - spcond(ieltcd,ivar): specified value for variable ivar associated
 *                        to the injected condensation in the ieltcd
 *                        cell with a condensation source term except
 *                        for ivar=ipr.
 *
 * \remarks
 *  - For each face where a condensation source term is imposed ielscd
 *    in [1;nfbpcd]), ifbpcd(ielscd) is the global index number of the
 *    corresponding face (ifbpcd(ieltcd) in [1;ncel]).
 *  - if itypcd(ieltcd,ivar)=0, spcond(ielpcd,ivar) is not used.
 *  - if spcond(ieltcd,ipr)<0, mass is removed from the system,
 *     therefore Code_Saturna automatically considers f_i=f^(n+1),
 *     whatever the values of itypcd or smacel specified by the user
 *
 *   \par Examples of settings for boundary condensation mass source terms
 *        Examples are available
 *        \ref condens_h_boundary "here".
 *
 * \param[in]  iappel     indicates which at which stage the routine is
 */
/*----------------------------------------------------------------------------*/

void
cs_user_wall_condensation(int  iappel);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_zones.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volume and surface zones.
 *
 * See \ref sec_selection_criteria for details on selection criteria.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_zones(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
