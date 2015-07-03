#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"

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
 * Initialize Fortran base common block values
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (csinit, CSINIT)
(
 const cs_int_t  *irgpar,  /* <-- MPI Rank in parallel, -1 otherwise */
 const cs_int_t  *nrgpar   /* <-- Number of MPI processes, or 1 */
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
 const cs_int_t  *nscal,     /* <-- number of scalars */
 const cs_int_t  *nvlsta,    /* <-- number of statistical variables (lagr) */
 const cs_int_t  *nvisbr,    /* <-- number of boundary stat. variables (lagr) */
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
 * Generator for distribution function of p's
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (fische, FISCHE)
(
 const cs_int_t   *n,
 const cs_real_t  *mu,
       cs_int_t    p[]);

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
 * Free Fortran allocated memory
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (memfin, MEMFIN) (void);

/*----------------------------------------------------------------------------
 * User function for enthalpy <-> temperature conversion
 *----------------------------------------------------------------------------*/

void CS_PROCF (usthht, USTHHT)
(
 const cs_int_t  *mode,      /* <-- -1 : t -> h ; 1 : h -> t */
 cs_real_t       *enthal,    /* <-- enthalpy */
 cs_real_t       *temper     /* <-- temperature */
);

/*----------------------------------------------------------------------------
 * User function for output of variables on a post-processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (usvpst, USVPST)
(
 const cs_int_t  *nummai,    /* <-- number or post-processing mesh */
 const cs_int_t  *nvar,      /* <-- number of variables */
 const cs_int_t  *nscal,     /* <-- number of scalars */
 const cs_int_t  *nvlsta,    /* <-- number of statistical variables (lagr) */
 const cs_int_t  *ncelps,    /* <-- number of post-processed cells */
 const cs_int_t  *nfacps,    /* <-- number of post processed interior faces */
 const cs_int_t  *nfbrps,    /* <-- number of post processed boundary faces */
 const cs_int_t   itypps[3], /* <-- flag (0 or 1) for presence of cells, */
                             /*     interior faces, and boundary faces */
 const cs_int_t   lstcel[],  /* <-- list of post-processed cells */
 const cs_int_t   lstfac[],  /* <-- list of post-processed interior faces */
 const cs_int_t   lstfbr[]   /* <-- list of post-processed boundary faces */
);

/*----------------------------------------------------------------------------
 * Uniform random number generator
 *----------------------------------------------------------------------------*/

void CS_PROCF (zufall, zufall)
(
 const cs_int_t   *n,             /* --> size of the vector */
 const cs_real_t  *a              /* <-- generated random number vector */
);

/*----------------------------------------------------------------------------
 * Gaussian random number generator
 *----------------------------------------------------------------------------*/

void CS_PROCF (normalen, normalen)
(
 const cs_int_t   *n,             /* --> size of the vector */
 const cs_real_t  *x              /* <-- generated random number vector */
);

/*----------------------------------------------------------------------------
 * Initialize Lagrangian module parameters for a given zone and class
 *
 * parameters:
 *   i_cz_params <-- integer parameters for this class and zone
 *   r_cz_params <-- real parameters for this class and zone
 *----------------------------------------------------------------------------*/

void
cs_lagr_init_zone_class_param(const cs_int_t   i_cs_params[],
                              const cs_real_t  r_cs_params[]);

/*----------------------------------------------------------------------------
 * Define Lagrangian module parameters for a given zone and class
 *
 * parameters:
 *   class_id    <-- id of given particle class
 *   zone_id     <-- id of given boundary zone
 *   i_cz_params <-- integer parameters for this class and zone
 *   r_cz_params <-- real parameters for this class and zone
 *----------------------------------------------------------------------------*/

void
cs_lagr_define_zone_class_param(cs_int_t         class_id,
                                cs_int_t         zone_id,
                                const cs_int_t   i_cs_params[],
                                const cs_real_t  r_cs_params[]);

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
 * Define global options for couplings.
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 *----------------------------------------------------------------------------*/

void
cs_user_coupling(void);

/*----------------------------------------------------------------------------
 * This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *----------------------------------------------------------------------------*/

void
cs_user_extra_operations(void);

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
 * Insert thin wall into a mesh.
 *----------------------------------------------------------------------------*/

void
cs_user_mesh_thinwall(cs_mesh_t  *mesh);

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


/*============================================================================
 *  CDO User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the specifications for the current simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the way geometric quantities
 *         are built
 *
 * \return the type of computation to evaluate the cell center
 */
/*----------------------------------------------------------------------------*/

cs_cdo_cc_algo_t
cs_user_cdo_geometric_settings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the numerical parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_numeric_settings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the discrete Hodge operators
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_hodge_settings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional operations on results produced by CDO schemes.
 *         Define advanced post-processing and/or analysis for instance.
 *
 * \param[in]  m        pointer to a cs_mesh_t struct.
 * \param[in]  mq       pointer to a cs_quantities_t struct.
 * \param[in]  connect  pointer to a cs_cdo_connect_t structure
 * \param[in]  quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]  tcur     current physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_extra_op(const cs_mesh_t             *m,
                     const cs_mesh_quantities_t  *mq,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     double                       tcur);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
