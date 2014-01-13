#ifndef __CS_PROTOTYPES_H__
#define __CS_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 const cs_real_t  rtp[],     /* <-- cell variables */
 const cs_real_t  propce[],  /* <-- cell physical properties */
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
 * Mesh renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (numvec, NUMVEC)
(
 const cs_int_t  *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t  *ncel,     /* <-- number of local cells */
 const cs_int_t  *nfac,     /* <-- number of interior faces */
 const cs_int_t  *nfabor,   /* <-- number of boundary faces */
 const cs_int_t  *lregis,   /* <-- vector registor length */
 cs_int_t        *irveci,   /* <-> interior face vectorization indic. */
 cs_int_t        *irvecb,   /* <-> boundary face vectorization indic. */
 cs_int_t         ifacel[], /* <-- interior face->cell connectivity */
 cs_int_t         ifabor[], /* <-- boundary face->cell connectivity */
 cs_int_t         inumfi[], /* <-> interior faces renumbering (size: nfac) */
 cs_int_t         inumfb[], /* <-> boundary faces renumbering (size: nfabor) */
 cs_int_t         iworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_int_t         ismbs[]   /* --- work array, size: ncelet */
);

/*----------------------------------------------------------------------------
 * Test renumbering for vector processors
 *----------------------------------------------------------------------------*/

extern void CS_PROCF (tstvec, TSTVEC)
(
 const cs_int_t  *ncelet,   /* <-- number of extended (real + ghost) cells */
 const cs_int_t  *ncel,     /* <-- number of local cells */
 const cs_int_t  *nfac,     /* <-- number of interior faces */
 const cs_int_t  *nfabor,   /* <-- number of boundary faces */
 const cs_int_t   ifacel[], /* <-- interior face->cell connectivity */
 const cs_int_t   ifabor[], /* <-- boundary face->cell connectivity */
 cs_int_t         iworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_int_t         ismbs[],  /* --- work array, size: ncelet */
 cs_int_t         ismbv[],  /* --- work array, size: ncelet */
 cs_real_t        rworkf[], /* --- work array, size: max(nfac, nfabor) */
 cs_real_t        rsmbs[],  /* --- work array, size: ncelet */
 cs_real_t        rsmbv[]   /* --- work array, size: ncelet */
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
 const cs_int_t   lstfbr[],  /* <-- list of post-processed boundary faces */
 const cs_real_t  dt[],      /* <-- local time step */
 const cs_real_t  rtpa[],    /* <-- cell variables at previous time step */
 const cs_real_t  rtp[],     /* <-- cell variables */
 const cs_real_t  propce[]   /* <-- cell physical properties */
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
 * Define mesh joinings.
 *----------------------------------------------------------------------------*/

void
cs_user_join(void);

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
 * Define rotor/stator model.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery(void);

/*----------------------------------------------------------------------------
 * Define rotor axes, associated cells, and rotor/stator faces.
 *----------------------------------------------------------------------------*/

void
cs_user_turbomachinery_rotor(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROTOTYPES_H__ */
