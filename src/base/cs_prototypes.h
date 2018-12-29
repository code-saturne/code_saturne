#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_domain.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_probe.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Fortran function/subroutine prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Main Fortran subroutine
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (caltri, CALTRI)
(
 void
);

/*----------------------------------------------------------------------------
 * Convert gas temperature to and from enthalpy based on concentrations
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (cpthp1, CPTHP1)
(
 const cs_int_t   *mode,    /* <-- 1: h to t, 2: t to h */
 cs_real_t        *eh,      /* <-> enthalpy of gas mix */
 cs_real_t        *xesp,    /* <-- mas fraction of species */
 cs_real_t        *f1mc,    /* <-- mean f1 */
 cs_real_t        *f2mc,    /* <-- mean f2 */
 cs_real_t        *tp       /* <-- gas temperature (K) */
);

/*----------------------------------------------------------------------------
 * Initialize Fortran base common block values
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csinit, CSINIT)
(
 const cs_int_t  *irgpar,  /* <-- MPI Rank in parallel, -1 otherwise */
 const cs_int_t  *nrgpar   /* <-- Number of MPI processes, or 1 */
);

/*----------------------------------------------------------------------------
 * Compute distance to wall by solving a 3d diffusion equation.
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (distpr, DISTPR)
(
 const cs_int_t  *itypfb,    /* <-- boudnary face types */
 cs_real_t       *distpa     /* <-- wall distance */
);

/*----------------------------------------------------------------------------
 * Developer function for output of variables on a post-processing mesh
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (dvvpst, DVVPST)
(
 const cs_int_t  *nummai,    /* <-- number or post-processing mesh */
 const cs_int_t  *numtyp,    /* <-- number or post-processing type
                              *     (-1 as volume, -2 as boundary, or nummai) */
 const cs_int_t  *nvar,      /* <-- number of variables */
 const cs_int_t  *ncelps,    /* <-- number of post-processed cells */
 const cs_int_t  *nfbrps,    /* <-- number of post processed boundary faces */
 const cs_int_t   lstcel[],  /* <-- list of post-processed cells */
 const cs_int_t   lstfbr[],  /* <-- list of post-processed boundary faces */
 cs_real_t        tracel[],  /* --- work array for output cells */
 cs_real_t        trafbr[]   /* --- work array for output boundary faces */
);

/*----------------------------------------------------------------------------
 * Find the nearest cell's center from a node
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (findpt, FINDPT)
(
 const cs_int_t   *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t   *ncel,     /* <-- number of cells */
 const cs_real_t  *xyzcen,   /* <-- cell centers */
 const cs_real_t  *xx,       /* <-- node coordinate X */
 const cs_real_t  *yy,       /* <-- node coordinate Y */
 const cs_real_t  *zz,       /* <-- node coordinate Z */
       cs_int_t   *node,     /* --> node we are looking for, zero if error */
       cs_int_t   *ndrang    /* --> rank of associated process */
);

/*----------------------------------------------------------------------------
 * Check necessity of extended mesh from FORTRAN options.
 *
 * Interface Fortran :
 *
 * SUBROUTINE HALTYP (IVOSET)
 * *****************
 *
 * INTEGER          IVOSET      : <-- : Indicator of necessity of extended mesh
 *----------------------------------------------------------------------------*/

extern void
CS_PROCF (haltyp, HALTYP)(const cs_int_t   *ivoset);

/*----------------------------------------------------------------------------
 * Main Fortran options initialization
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (initi1, INITI1)
(
 void
);

/*----------------------------------------------------------------------------
 * Set the CDO mode in the FORTRAN part
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (set_cdo_mode, SET_CDO_MODE)
(
 const cs_int_t   *mode     /* <-- -1: no CDO, 1: with CDO, 2: CDO only */
);

/*----------------------------------------------------------------------------
 * User function for enthalpy <-> temperature conversion
 *----------------------------------------------------------------------------*/

void CS_PROCF (usthht, USTHHT)
(
 const cs_int_t  *mode,      /* <-- -1 : t -> h ; 1 : h -> t */
 cs_real_t       *enthal,    /* <-- enthalpy */
 cs_real_t       *temper     /* <-- temperature */
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute GUI-defined head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * cku11, cku22, cku33, cku12, cku13, cku23.
 *
 * \param[in]       zone  pointer to zone structure
 * \param[in, out]  cku   head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_user_head_losses(const cs_zone_t  *zone,
                    cs_real_t         cku[][6]);

/*----------------------------------------------------------------------------
 * Absorption coefficient for radiative module
 *----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_absorption(const int         bc_type[],
                                const cs_real_t   dt[],
                                cs_real_t         ck[]);

/*----------------------------------------------------------------------------
 * Compute the net radiation flux
 *----------------------------------------------------------------------------*/

void
cs_user_rad_transfer_net_flux(const int        itypfb[],
                              const cs_real_t  dt[],
                              const cs_real_t  coefap[],
                              const cs_real_t  coefbp[],
                              const cs_real_t  cofafp[],
                              const cs_real_t  cofbfp[],
                              const cs_real_t  twall[],
                              const cs_real_t  qincid[],
                              const cs_real_t  xlam[],
                              const cs_real_t  epa[],
                              const cs_real_t  eps[],
                              const cs_real_t  ck[],
                              cs_real_t        net_flux[]);

/*----------------------------------------------------------------------------
 * Convert temperature to enthalpy at boundary
 *----------------------------------------------------------------------------*/

void CS_PROCF (b_t_to_h, b_t_to_h)
(
 const cs_lnum_t *nlst,          /* --> number of faces in list */
 const cs_lnum_t *lstfac,        /* --> list of boundary faces at which
                                    conversion is requested */
 const cs_real_t *t_b,           /* --> temperature at boundary */
 cs_real_t       *h_b            /* --> enthalpy at boundary */
);

/*----------------------------------------------------------------------------
 * Convert enthalpy to temperature at cells
 *----------------------------------------------------------------------------*/

void CS_PROCF (c_h_to_t, c_h_to_t)
(
 const cs_real_t *h,           /* --> enthalpy */
 cs_real_t       *t            /* --> temperature */
);

/*----------------------------------------------------------------------------
 * Add field indexes associated with a new non-user solved variable,
 * with default options
 *
 * parameters:
 *   f_id <--   field id
 *
 * returns:
 *   scalar number for defined field
 *----------------------------------------------------------------------------*/

int
cs_add_model_field_indexes(int f_id);

/*----------------------------------------------------------------------------
 * Return Lagrangian model status.
 *
 * parameters:
 *   model_flag   --> 0 without Lagrangian, 1 or 2 with Lagrangian
 *   restart_flag --> 1 for Lagrangian restart, 0 otherwise
 *   frozen_flag  --> 1 for frozen Eulerian flow, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_lagr_status(int  *model_flag,
               int  *restart_flag,
               int  *frozen_flag);

/*============================================================================
 *  User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Data Entry of the 1D wall thermal module.
 *----------------------------------------------------------------------------*/

void
cs_user_1d_wall_thermal(int iappel,
                        int isuit1);

/*----------------------------------------------------------------------------
 * Define global options for couplings.
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 *----------------------------------------------------------------------------*/

void
cs_user_coupling(void);

/*----------------------------------------------------------------------------
 * This function is called at each time step for boundary conditions.
 *----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         icodcl[],
                            int         bc_type[],
                            cs_real_t   rcodcl[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables.
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
cs_user_extra_operations_initialize(cs_domain_t     *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of the calculation.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.

 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations_finalize(cs_domain_t     *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called one time step to initialize problem.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain);

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
 * \brief Define volumes as internal coupling zones.
 *
 * These zones will be separated from the rest of the domain using automatically
 * defined thin walls.
 *
 * \param[in, out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_add_volumes(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volumesi from separated meshes as internal coupling zones.
 *
 * These zones must be disjoint and the face selection criteria must be
 * specified.
 *
 * \param[in, out]  mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_internal_coupling_from_disjoint_meshes(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * This function is called each time step to define physical properties.
 *----------------------------------------------------------------------------*/

void
cs_user_physical_properties(const cs_mesh_t             *mesh,
                            const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the porosity (volume factor \f$ \epsilon \f$
 *        when the porosity model is activated
 *        (iporos greater than 1 in cs_user_parameters.f90).
 *
 * This function is called at the begin of the simulation only.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_porosity(void);

/*----------------------------------------------------------------------------
 * Define mesh joinings.
 *----------------------------------------------------------------------------*/

void
cs_user_join(void);

/*----------------------------------------------------------------------------
 * Define linear solver options.
 *
 * This function is called at the setup stage, once user and most model-based
 * fields are defined.
 *----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_output(void);

/*----------------------------------------------------------------------------
 * Tag bad cells within the mesh based on geometric criteria.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_bad_cells_tag(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Define mesh files to read and optional associated transformations.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_input(void);

/*----------------------------------------------------------------------------
 * Modifiy geometry and mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Insert boundary wall into a mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_boundary(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Mesh smoothing.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_smoothe(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Enable or disable mesh saving.
 *
 * By default, mesh is saved when modified.
 *
 * parameters:
 *   mesh <-> pointer to mesh structure
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_save(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Set options for cutting of warped faces
 *
 * parameters:
 *   mesh <-> pointer to mesh structure to smoothe
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_warping(void);

/*----------------------------------------------------------------------------
 * Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup.
 *----------------------------------------------------------------------------*/

void
cs_user_model(void);

/*----------------------------------------------------------------------------
 * Define advanced mesh numbering options.
 *----------------------------------------------------------------------------*/

void
cs_user_numbering(void);

/*----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_user_parallel_io(void);

/*----------------------------------------------------------------------------
 * Define advanced partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_user_partition(void);

/*----------------------------------------------------------------------------
 * Define sparse matrix tuning options.
 *----------------------------------------------------------------------------*/

void
cs_user_matrix_tuning(void);

/*----------------------------------------------------------------------------
 * Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so speciic settings related to those
 * fields may be set here.
 *----------------------------------------------------------------------------*/

void
cs_user_parameters(void);

/*----------------------------------------------------------------------------
 * User function for input of radiative transfer module options.
 *----------------------------------------------------------------------------*/

void
cs_user_radiative_transfer_parameters(void);

/*-----------------------------------------------------------------------------
 * User subroutine for input of radiative transfer boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_user_radiative_transfer_bcs(int               nvar,
                               const int         bc_type[],
                               int               icodcl[],
                               int               isothp[],
                               cs_real_t        *tmin,
                               cs_real_t        *tmax,
                               cs_real_t        *tx,
                               const cs_real_t   dt[],
                               cs_real_t         rcodcl[],
                               const cs_real_t   thwall[],
                               const cs_real_t   qincid[],
                               cs_real_t         hfcnvp[],
                               cs_real_t         flcnvp[],
                               cs_real_t         xlamp[],
                               cs_real_t         epap[],
                               cs_real_t         epsp[],
                               cs_real_t         textp[],
                               cs_real_t         tintp[]);

/*----------------------------------------------------------------------------
 * Define periodic faces.
 *----------------------------------------------------------------------------*/

void
cs_user_periodicity(void);

/*----------------------------------------------------------------------------
 * Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_writers(void);

/*-----------------------------------------------------------------------------
 * Define monitoring probes and profiles. A profile is seen as a set of probes.
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void);

/*----------------------------------------------------------------------------
 * Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_meshes(void);

/*----------------------------------------------------------------------------
 * User function for output of values on a post-processing mesh.
 *----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * parameters:
 *   nt_max_abs <-- maximum time step number
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs);

/*----------------------------------------------------------------------------
 * Define couplings with other instances of Code_Saturne.
 *----------------------------------------------------------------------------*/

void
cs_user_saturne_coupling(void);

/*----------------------------------------------------------------------------
 * Set user solver.
 *----------------------------------------------------------------------------*/

int
cs_user_solver_set(void);

/*----------------------------------------------------------------------------
 * Main call to user solver.
 *----------------------------------------------------------------------------*/

void
cs_user_solver(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Define couplings with SYRTHES code.
 *----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void);

/*----------------------------------------------------------------------------
 * Define time moments.
 *----------------------------------------------------------------------------*/

void
cs_user_time_moments(void);

/*----------------------------------------------------------------------------
 * Define rotor/stator model.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery(void);

/*----------------------------------------------------------------------------
 * Define rotor axes, associated cells, and rotor/stator faces.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_rotor(void);

/*----------------------------------------------------------------------------
 * Define rotation velocity of rotor.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_set_rotation_velocity(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volume and surface zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_zones(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define scaling parameter for electric model
*/
/*----------------------------------------------------------------------------*/

void
cs_user_scaling_elec(const cs_mesh_t             *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities,
                     cs_real_t                   *dt);

/*----------------------------------------------------------------------------
 * Computation of the relaxation time-scale to equilibrium in the frame of
 * the homogeneous two-phase model.
 *----------------------------------------------------------------------------*/

void
cs_user_hgn_thermo_relax_time(const cs_mesh_t *mesh,
                              const cs_real_t *alpha_eq,
                              const cs_real_t *y_eq,
                              const cs_real_t *z_eq,
                              const cs_real_t *ei,
                              const cs_real_t *v,
                              cs_real_t       *relax_tau);

/*============================================================================
 *  CDO User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  After the first step: cs_user_cdo_init_setup(), this second step
 *         concludes the setup of properties, equations, source terms...
 *         At this step, mesh quantities and connectivities are build as well
 *         as the field arrays.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_finalize_setup(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional user-defined operations on results provided by the CDO
 *         kernel. Define advanced post-processing and analysis for example.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(cs_domain_t     *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Final step for user-defined operations on results provided by the
 *         CDO kernel.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_end_extra_op(cs_domain_t     *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for each soil and tracer how is defined each term of the
 *         the tracer equation. Soils and tracer equations have to be added
 *         previously
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_gwf_setup(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
