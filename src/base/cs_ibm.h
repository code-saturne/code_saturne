#ifndef __CS_IBM_H__
#define __CS_IBM_H__

/*============================================================================
 * Time and space immersed boundaries model.
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

#include "base/cs_defs.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_medcoupling_intersector.h"
#include "mesh/cs_stl.h"
#include "cdo/cs_xdef.h"

/*----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

BEGIN_C_DECLS

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*---------------------------------------------------------------------------
 * Choice for porosity model
 *---------------------------------------------------------------------------*/

typedef enum
{
  CS_IBM_OFF          =  0,  /* No immersed boundary model    */
  CS_IBM_FIXED_SOLID  =  1,  /* Fix solid porosity model      */

} cs_ibm_model_type_t;

/*---------------------------------------------------------------------------
 * Detection of problem dimension type for immersed boundaries
 *---------------------------------------------------------------------------*/

typedef enum
{
  CS_IBM_3D    = -1,  /* 3D computation                              */
  CS_IBM_2D_X  =  0,  /* 2D computation with symmetry in X direction */
  CS_IBM_2D_Y  =  1,  /* 2D computation with symmetry in Y direction */
  CS_IBM_2D_Z  =  2   /* 2D computation with symmetry in Z direction */

} cs_ibm_prob_dim_type_t;

/*---------------------------------------------------------------------------
 * Algorithm choice for porosity computation
 *---------------------------------------------------------------------------*/

typedef enum
{
  CS_IBM_ALGO_NONE       ,  /* Default init before knowing the chosen algo */
  CS_IBM_ALGO_CUT_CELLS  ,  /* Cut-cells: optimised cutting     */
  CS_IBM_ALGO_MEDCOUPLING,  /* MEDCoupling sequential volume interpolation */
  CS_IBM_ALGO_STL        ,  /* Computation from STL ASCII file */

  CS_N_VAR_PORO_ALGO_TYPES  /* Number of algo choices */

} cs_ibm_algo_type_t;

/*---------------------------------------------------------------------------
 * Choice for B.C. at walls
 *---------------------------------------------------------------------------*/

typedef enum
{
  CS_IBM_SLIP_WALL_CONDITION    ,  /* Slip wall b.c.     */
  CS_IBM_NO_SLIP_WALL_CONDITION ,  /* No-Slip wall b.c.  */
  CS_IBM_WALL_LAW_WALL_CONDITION   /* Wall law wall b.c. */

} cs_ibm_wall_cond_type_t;

/*---------------------------------------------------------------------------
 * Objects mechanical or thermophysical properties
 *---------------------------------------------------------------------------*/

typedef enum
{

  CS_IBM_OBJ_PROP_DENSITY,
  NC_IBM_OBJ_PROP_MASS,
  NC_IBM_OBJ_PROP_INERTIA_MATRIX,
  CS_IBM_OBJ_PROP_CP,
  CS_IBM_OBJ_PROP_LAMBDA,
  CS_IBM_OBJ_PROP_STIFFNESS,
  CS_IBM_OBJ_PROP_DAMPING,
  CS_IBM_OBJ_PROP_YOUNG_MODULE,
  CS_IBM_OBJ_PROP_INERTIA_MOM,
  CS_IBM_OBJ_PROP_CROSS_SECTION,
  CS_IBM_OBJ_PROP_RAYLEIGH_DAMP_A,
  CS_IBM_OBJ_PROP_RAYLEIGH_DAMP_B,

  CS_N_IBM_OBJ_PROP_TYPES

} cs_ibm_object_property_type_t;

typedef enum
{
  CS_IBM_OBJ_INIT_COG_EQ,
  CS_IBM_OBJ_INIT_COG,
  CS_IBM_OBJ_INIT_ANGLE,
  CS_IBM_OBJ_INIT_VELOCITY,
  CS_IBM_OBJ_INIT_ACCELERATION,
  CS_IBM_OBJ_INIT_ANGULAR_VEL,
  CS_IBM_OBJ_INIT_FLUID_FORCE,

  CS_N_IBM_OBJ_INIT_TYPES

} cs_ibm_object_init_param_t;

/*---------------------------------------------------------------------------
 * Structure associated to objects management
 *---------------------------------------------------------------------------*/

typedef int
(cs_cutcell_func_t)(const cs_lnum_t    c_id,
                    const cs_real_3_t  xyz,
                    const cs_real_t    t,
                    const int          num_object);

/*---------------------------------------------------------------------------
 * Structure associated to objects management
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Porosity immersed boundaries model options descriptor
 *----------------------------------------------------------------------------*/

typedef struct {
  int  porosity_mode;
} cs_porosity_ibm_opt_t;

/* Pointer to options structure */
extern cs_porosity_ibm_opt_t *cs_glob_porosity_ibm_opt;

typedef struct
{

  /* --------------------------------------- */
  /* Name of the object */
  /* --------------------------------------- */
  char                          *name;

  /* --------------------------------------- */
  /* Computation method */
  /* --------------------------------------- */
  cs_ibm_algo_type_t             method;

  /* --------------------------------------- */
  /* Method specific pointers */
  cs_cutcell_func_t             *cutcell_func;
  cs_stl_mesh_t                 *stl;
  cs_medcoupling_intersector_t  *mi;
  /* --------------------------------------- */

  /* --------------------------------------- */
  /* Object thermophysical properties */
  cs_xdef_t *property_defs[CS_N_IBM_OBJ_PROP_TYPES];
  /* --------------------------------------- */

  /* --------------------------------------- */
  /* Initial parameters */
  cs_xdef_t *init_vals_defs[CS_N_IBM_OBJ_INIT_TYPES];
  /* --------------------------------------- */

} cs_ibm_object_t;

/*---------------------------------------------------------------------------
 * Structure associated to immersed boundaries management
 *---------------------------------------------------------------------------*/

typedef struct
{
  /* ---------------------------------------------- */
  /* Objects defined using stl or medcoupling. */
  int                           n_objects;
  cs_ibm_object_t             **objects;
  /* ---------------------------------------------- */

  /* Detection of prob dim type for immersed boundaries   */
  cs_ibm_prob_dim_type_t   prob_dim;
  /* Algorithm choice for porosity computation            */
  cs_ibm_algo_type_t       algo_choice;
  /* Choice for velocity B.C. at walls                    */
  cs_ibm_wall_cond_type_t  wall_condition;
  /* Number of sub-cut for cells using the Cut-Cells algo */
  int           nb_cut_cells;
  /* Number of sub-cut for faces using the Cut-Cells algo */
  int           nb_cut_faces;
  /* Porosity dynamic modification                        */
  bool          porosity_user_source_term_modification;
  /* Limitation area for porosity calculation             */
  cs_real_3_t   xyzmin_moving_porosity;
  /* Limitation area for porosity calculation             */
  cs_real_3_t   xyzmax_moving_porosity;
  /* Solid internal porosity                              */
  cs_real_t     *solid_porosity;
  /* Fluid volume at the first time step                  */
  cs_real_t     isovol;
  /* Keep same volume for porous object at each iteration */
  bool          ensure_isovol;
  /* Cell porosity based on nodes porosity (smoothing)    */
  bool          porosity_from_nodes;

} cs_ibm_t;


/*============================================================================
 *  Global variables
 *============================================================================*/

/* Pointer to the cs_ibm structure for various arrays */

extern cs_ibm_t  *cs_ibm;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its id.
 *
 * \param[in] obj_id  id of the object
 *
 * \return pointer to object structure.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_id(int obj_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Try to get an object based on its name. Returns NULL if not found
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure, NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name_try(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its name. Error if not found
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure, NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name(const char *name);

/*----------------------------------------------------------------------------
 * Create an empty cs_ibm structure for various arrays
 *
 * returns:
 *   pointer to created cs_ibm structure
 *----------------------------------------------------------------------------*/

cs_ibm_t *
cs_ibm_create(void);

/*----------------------------------------------------------------------------
 * Destroy a cs_ibm structure
 *
 * cs_ibm <-- pointer to a cs_ibm_t structure
 *
 * returns:
 *   NULL
 *----------------------------------------------------------------------------*/

void
cs_ibm_finalize(void);

/*----------------------------------------------------------------------------
 * Define immersed boundaries in time and space (solid(s) interior part).
 *----------------------------------------------------------------------------*/

void cs_immersed_boundaries(const cs_mesh_t *mesh,
                            const cs_mesh_quantities_t *mesh_quantities);

/*----------------------------------------------------------------------------
 * Define space immersed boundaries on a set of zones defined by the user in the
 * GUI
 *----------------------------------------------------------------------------*/

void cs_volumic_zone_porosity(const cs_mesh_quantities_t *mesh_quantities);

/*----------------------------------------------------------------------------*
 * User function. Locally Modify porosity (erosion, fouling effects, ...)
 *----------------------------------------------------------------------------*/

void cs_user_ibm_modify(const cs_mesh_t *mesh,
                        const cs_mesh_quantities_t *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its id.
 *
 * \param[in] obj_id  id of the object
 *
 * \return pointer to object structure.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_id(int obj_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an object based on its name.
 *
 * \param[in] name  name of the object to get
 *
 * \return pointer to object structure.
 */
/*----------------------------------------------------------------------------*/

cs_ibm_object_t *
cs_ibm_object_by_name_try(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a point is solid or fluid based on the cut-cell method.
 *
 * \param[in]  c_id         local cell number
 * \param[out] ipenal       indicator for cut cells algo
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

int
cs_ibm_object_compute_cut_porosity(const cs_lnum_t    c_id,
                                   const cs_real_3_t  xyz,
                                   const cs_real_t    t,
                                   const int          num_object);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an object from a file using STL or MED formats
 *
 * \param[in] name          name of the object
 * \param[in] method        Porosity computation method
 * \param[in] file_name     file name
 * \param[in] solve_fsi     Is the object used in the FSI resolution ?
 *                          (currently ignored)
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_object_from_file(const char          *name,
                            cs_ibm_algo_type_t   method,
                            const char          *file_name,
                            bool                 solve_fsi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define an object from a function used in the cutcell algorithm
 *
 * \param[in] name          name of the object
 * \param[in] cutcell_func  pointer to the cutcell function of the object
 * \param[in] solve_fsi     Is the object used in the FSI resolution ?
 *                          (currently ignored)
 * \param[in] n_nodes       Number of nodes if the object is deformable
 *                          (currently ignored)
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_add_object_from_func(const char        *name,
                            cs_cutcell_func_t *cutcell_func,
                            bool               solve_fsi,
                            int                n_nodes);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define exterior points for an STL object.
 *
 * \param[in] name          name of the object
 * \param[in] n_pts         number of points
 * \param[in] pts_coords    coordinates of the points
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_stl_define_ext_points(const char      *name,
                             const int        n_pts,
                             cs_real_t       *pts_coords);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotate an object based on the STL or MED algorithms
 *
 * \param[in] name          name of the object
 * \param[in] angle         angle of rotation
 * \param[in] axis          axis of rotation
 * \param[in] center        center of rotation
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_rotate(const char *name,
                     cs_real_t   angle,
                     cs_real_t   axis[3],
                     cs_real_t   center[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a new constant property definition for an object.
 *
 * \param[in] obj       pointer to object
 * \param[in] ppty_id   property id (si enum for list)
 * \param[in] val       property constant value
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_set_property_const(cs_ibm_object_t               *obj,
                                 cs_ibm_object_property_type_t  ppty_id,
                                 cs_real_t                      val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Translate an object based on the STL or MED algorithms
 *
 * \param[in] name          name of the object
 * \param[in] vector        translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_translate(const char *name,
                        cs_real_t   vector[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scale an object based on a factor
 *
 * \param[in] name          name of the object
 * \param[in] factor        scaling factor
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_scale(const char *name,
                    cs_real_t   factor);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply user parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_user_parameters(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Init writers for STL or MED objects.
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_init_writer(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transform an object from its initial state using a transformation
 *         matrix.
 *
 * \param[in] obj     pointer to object structure
 * \param[in] matrix  transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_transform_from_init(cs_ibm_object_t *obj,
                                  cs_real_34_t     matrix);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the volume fraction of an object over all cells.
 *
 * \param[in]  obj            pointer to object structure
 * \param[in]  m              pointer to mesh structure
 * \param[in]  cell_vol       pointer to cell volume array
 * \param[out] obj_frac_tot   array containing the total vol fraction of solids
 * \param[in]  indic          indicator array
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_object_compute_intersect_vol(cs_ibm_object_t            *obj,
                                    const cs_mesh_t            *m,
                                    const cs_real_t            *cell_vol,
                                    cs_real_t                  *obj_frac_tot,
                                    int                        *indic);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function in which the user defines the objects to model.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_define_objects(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function to set global parameters for the immersed boundaries
 *         module.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_parameters(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function where to apply predefined transformations to med/stl
 *         based objects.
 *
 * \param[in]  t            time value for the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_object_transformations(const cs_real_t time);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function which allows the definition of a 'porous' object.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_ibm_solid_por(const cs_lnum_t    c_id,
                      const cs_real_3_t  xyz,
                      const cs_real_t    t,
                      const int          num_object);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define space immersed boundaries on a set of zones defined by the user
 *         in the GUI.
 *
 * \param[in]  mesh_quantities  pointer to associated mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ibm_volumic_zone(const cs_mesh_quantities_t *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_IBM_H__ */
