#ifndef __CS_PROPERTY_H__
#define __CS_PROPERTY_H__

/*============================================================================
 * Manage the definition/setting of properties
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

#include "cs_field.h"
#include "cs_flag.h"
#include "cs_param_types.h"
#include "cs_xdef.h"
#include "cs_xdef_cw_eval.h"
#include "cs_xdef_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Common property names (property which is shared between different module) */

#define CS_PROPERTY_MASS_DENSITY   "mass_density"

/*!
 * @defgroup cdo_property_flags Flags specifying metadata related to the
 *  post-processing for a property
 * @{
 */

/*!  1: Perform the computation and post-processing of the Fourier number */

#define CS_PROPERTY_POST_FOURIER  (1 << 0)

/*! @} */

/*!
 * @defgroup cdo_property_type Flags specifying metadata related to the
 *  the type of property
 * @{
 */

/*! \var CS_PROPERTY_ISO
 *  1: Isotropic behavior (one real number is sufficient to describe the
 *  property) */

#define CS_PROPERTY_ISO                      (1 << 0)

/*! \var CS_PROPERTY_ORTHO
 *  2: Orthotropic behavior (three real numbers describe the behavior assuming
 *  that the different behavior is aligned with Cartesian axis) */

#define CS_PROPERTY_ORTHO                    (1 << 1)

/*! \var CS_PROPERTY_ANISO
 *  4: Anisotropic behavior (a 3x3 tensor describe the behavior). This tensor
 *  should be symmetric positive definite (i.e 6 real numbers describe the
 *  behavior) but by default a 3x3 tensor is used. */

#define CS_PROPERTY_ANISO                    (1 << 2)

/*! \var CS_PROPERTY_ANISO_SYM
 *  8: Anisotropic behavior. This tensor is represented with 6 real numbers
 *  since the tensor is symmetric */

#define CS_PROPERTY_ANISO_SYM                (1 << 3)

/*! \var CS_PROPERTY_BY_PRODUCT
 *  16: The property is defined as the product of two other properties
 */

#define CS_PROPERTY_BY_PRODUCT               (1 << 4)

/*! \var CS_PROPERTY_SCALED
 *  32: The property is defined up to a scaling factor
 */

#define CS_PROPERTY_SCALED                   (1 << 5)

/*! \var CS_PROPERTY_SUBCELL_DEFINITION

 *  64: The property is defined such that one wants to evaluate the definition
 *  on entities which are a partition of a cell. By default, one perfoms only
 *  one evaluation in each cell
 */

#define CS_PROPERTY_SUBCELL_DEFINITION       (1 << 6)

/*! @} */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef cs_flag_t cs_property_type_t;

/*! \enum cs_property_key_t
 *  \brief List of available keys for setting options on a property
 *
 * \var CS_PTYKEY_POST_FOURIER
 * Perform the computation (and post-processing) of the Fourier number
 */

typedef enum {

  CS_PTYKEY_POST_FOURIER,
  CS_PTYKEY_N_KEYS

} cs_property_key_t;

/* ======================================== */
/* Set of parameters attached to a property */
/* ======================================== */

/*!
 * \struct cs_property_t
 * \brief Structure associated to the definition of a property relying on the
 * \ref cs_xdef_t structure
 */

typedef struct _cs_property_t cs_property_t;

struct _cs_property_t {

  char *restrict       name;
  int                  id;
  cs_flag_t            state_flag;
  cs_flag_t            process_flag;
  cs_property_type_t   type;

  /* Reference value wich is used as default when nothing else is set. This
   * value can also be used to renormalized quantities related to this property
   * By default, this is set to 1
   */

  cs_real_t            ref_value;

  /* Value of the scaling factor to apply when a property is defined up to a
   * scaling factor. By default the value is 1.0 and the property is not
   * scaled.
   */

  cs_real_t            scaling_factor;

  /* Property is up to now only defined on the whole domain (volume) */

  int                  n_definitions;  /* Current number of definitions used */
  cs_xdef_t          **defs;           /* List of definitions */

  /* Store the definition id for each cell, null if there is only one
     definition set */

  short int           *def_ids;

  /* Function pointers to handle generic tasks related to a property. There
     is one function related to each definition. Some functions may not be
     allocated according to the kind of property */

  /* Retrieve the evaluation of the property at the cell center for each
     definition */

  cs_xdef_eval_t        **get_eval_at_cell;

  /* Same thing as the previous one but now with the usage of cellwise algo.
     relying on a cs_cell_mesh_t structure */

  cs_xdef_cw_eval_t     **get_eval_at_cell_cw;

  /* For properties relying on other properties for their definition, one
   * stores the pointers to these related properties */

  int                     n_related_properties;
  const cs_property_t   **related_properties;

  /* Optional: Definition(s) of the property at the boundary. If nothing is
     set, then one uses the definition on the volume to get information on the
     boundary. This is the default behavior. */

  int                     n_b_definitions;
  cs_xdef_t             **b_defs;

  /* Store the definition id for each boundary face, null if there is only one
     boundary definition set */

  short int              *b_def_ids;

};

/*!
 * \struct cs_property_data_t
 * \brief Structure storing the evaluation of a property and its related
 *        data
 */

typedef struct {

  const cs_property_t   *property; /* shared pointer */

  bool                   is_iso;   /* Detect if this an easier case */
  bool                   is_unity; /* Detect if this the simplest case */

  /* Way to evaluate the property and buffer to store the evaluation done
     on-the-fly */

  bool                   need_tensor;
  cs_real_t              tensor[3][3];
  cs_real_t              value;

  /* Additional metadata */

  bool                   need_eigen;
  cs_real_t              eigen_max;
  cs_real_t              eigen_ratio;

} cs_property_data_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns true if the property is steady and uniform, otherwise false
 *
 * \param[in] pty    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_property_is_constant(const cs_property_t   *pty)
{
  if (pty == NULL)
    return true; /* Treated as the "unity" property */

  if (pty->state_flag & CS_FLAG_STATE_STEADY) {
    if (pty->state_flag & CS_FLAG_STATE_UNIFORM)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns true if the property is steady, otherwise false
 *
 * \param[in] pty    pointer to a property to test
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_property_is_steady(const cs_property_t   *pty)
{
  if (pty == NULL)
    return true; /* Treated as the "unity" property */

  if (pty->state_flag & CS_FLAG_STATE_STEADY)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns true if the property is uniform, otherwise false
 *
 * \param[in] pty    pointer to a property to test
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_property_is_uniform(const cs_property_t   *pty)
{
  if (pty == NULL)
    return true; /* Treated as the "unity" property */

  if (pty->state_flag & CS_FLAG_STATE_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns true if the property is isotropic, otherwise false
 *
 * \param[in] pty    pointer to a property to test
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_property_is_isotropic(const cs_property_t   *pty)
{
  if (pty == NULL)
    return false;

  if (pty->type & CS_PROPERTY_ISO)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns true if the property is allowed to be evaluated on a
 *        sub-partition of a cell
 *
 * \param[in] pty    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_property_is_subcell(const cs_property_t   *pty)
{
  if (pty == NULL)
    return false;

  if (pty->type & CS_PROPERTY_SUBCELL_DEFINITION)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the name of a property
 *
 * \param[in] pty    pointer to a property
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

static inline const char *
cs_property_get_name(const cs_property_t   *pty)
{
  if (pty == NULL)
    return NULL;

  return pty->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the type of a property
 *
 * \param[in] pty    pointer to a property
 *
 * \return the type of the related property
 */
/*----------------------------------------------------------------------------*/

static inline cs_property_type_t
cs_property_get_type(const cs_property_t   *pty)
{
  if (pty == NULL)
    return 0; /* means undefined */

  return pty->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the dimension of the property
 *
 * \param[in] pty    pointer to a property
 *
 * \return the value of the dimension
 */
/*----------------------------------------------------------------------------*/

static inline int
cs_property_get_dim(const cs_property_t   *pty)
{
  if (pty == NULL)
    return 0; /* means undefined */

  if (pty->type & CS_PROPERTY_ISO)
    return 1;
  else if (pty->type & CS_PROPERTY_ORTHO)
    return 3;
  else if (pty->type & CS_PROPERTY_ANISO_SYM)
    return 6;
  else if (pty->type & CS_PROPERTY_ANISO)
    return 9;
  else
    return 0; /* means undefined */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set shared pointers to main domain members
 *
 * \param[in] mesh        mesh structure shared between FV and CDO
 * \param[in] quant       additional mesh quantities struct.
 * \param[in] connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_property_init_sharing(const cs_mesh_t              *mesh,
                         const cs_cdo_quantities_t    *quant,
                         const cs_cdo_connect_t       *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the number of properties
 *
 * \return the number of properties
 */
/*----------------------------------------------------------------------------*/

int
cs_property_get_n_properties(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a new property structure
 *
 * \param[in] name        name of the property
 * \param[in] type        type of property
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add(const char            *name,
                cs_property_type_t     type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a new property structure with an evaluation
 *        which can be called on a sub-partition of a cell.
 *        This kind of property is not available for all numerical scheme.
 *        By default, only one evaluation is performed in each cell.
 *
 * \param[in] name          name of the property
 * \param[in] type          type of property
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_subcell_add(const char            *name,
                        cs_property_type_t     type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to the product of two
 *        properties
 *        The type is infered from that of the related properties
 *        The value of the property is given as value_ab = value_a * value_b
 *
 * \param[in] name      name of the property
 * \param[in] pty_a     pointer to a cs_property_t structure
 * \param[in] pty_b     pointer to a cs_property_t structure
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_add_as_product(const char             *name,
                           const cs_property_t    *pty_a,
                           const cs_property_t    *pty_b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the related property definition from its name
 *
 * \param[in] name    name of the property to find
 *
 * \return null if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_name(const char   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find the related property definition from its id
 *
 * \param[in] id      id of the property to find
 *
 * \return null if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_by_id(int         id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set optional parameters related to a cs_property_t structure
 *
 * \param[in, out] pty       pointer to a cs_property_t structure
 * \param[in]      key       key related to a setting option
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_option(cs_property_t       *pty,
                       cs_property_key_t    key);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the reference value associated to a \ref cs_property_t structure
 *        This is a real number even whatever the type of property is.
 *
 * \param[in, out] pty      pointer to a cs_property_t structure
 * \param[in]      refval   value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_reference_value(cs_property_t    *pty,
                                double            refval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling factor associated to a \ref cs_property_t structure
 *        This is a real number whatever the type of property is. If the
 *        property was not defined as CS_PROPERTY_SCALED, then this tag is
 *        added.
 *
 * \param[in, out] pty    pointer to a cs_property_t structure
 * \param[in]      val    value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_scaling_factor(cs_property_t    *pty,
                               double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling factor associated to a \ref cs_property_t structure
 *        This is a real number whatever the type of property is. If the
 *        property was not defined as CS_PROPERTY_SCALED, then this tag is
 *        added.
 *
 * \param[in, out] pty    pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_unscale(cs_property_t    *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all cs_property_t structures and the array storing all the
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last stage of the definition of a property based on several
 *        definitions (i.e. definition by subdomains)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_finalize_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the array associated to the volume definition for the given
 *        property.
 *        Available only if there is one definition by array for the volume.
 *
 * \param[in] pty     pointer to the property structure
 *
 * \return a pointer to the array or null
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_property_get_array(const cs_property_t     *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a \ref cs_property_data_t structure (not a pointer to this
 *         structure). If property is null then one considers that this is a
 *         unitary property
 *
 * \param[in]   need_tensor  true if one needs a tensor-valued evaluation
 * \param[in]   need_eigen   true if one needs an evaluation of eigen values
 * \param[in]   property     pointer to the \ref cs_property_t structure
 *
 * \return an initialized structure
 */
/*----------------------------------------------------------------------------*/

cs_property_data_t
cs_property_data_define(bool                     need_tensor,
                        bool                     need_eigen,
                        const cs_property_t     *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_property_data_t structure. If property is null
 *         then one considers that this is a unitary property
 *
 * \param[in]      need_tensor  true if one needs a tensor-valued evaluation
 * \param[in]      need_eigen   true if one needs an evaluation of eigen values
 * \param[in]      property     pointer to the \ref cs_property_t structure
 * \param[in, out] data         structure to initialize (already allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_data_init(bool                     need_tensor,
                      bool                     need_eigen,
                      const cs_property_t     *property,
                      cs_property_data_t      *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a single uniform and steady isotropic definition for the
 *         given cs_property_t structure.
 *         This is a specialized variant of \ref cs_property_def_iso_by_value
 *         since several assumptions are satisfied.
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_constant_value(cs_property_t    *pty,
                               double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_iso_by_value(cs_property_t    *pty,
                             const char       *zname,
                             double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value at the boundary for the given (isotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       val      value to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_iso_by_value(cs_property_t    *pty,
                                      const char       *zname,
                                      double            val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an orthotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]       val      values to set (vector of size 3)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_ortho_by_value(cs_property_t    *pty,
                               const char       *zname,
                               double            val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value at the boundary for the given (orthotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       vals     values to set
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_ortho_by_value(cs_property_t    *pty,
                                        const char       *zname,
                                        double            vals[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         related to a volume zone
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]       tens     values to set (3x3 tensor)
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_by_value(cs_property_t    *pty,
                               const char       *zname,
                               cs_real_t         tens[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value at the boundary for the given (anisotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       tens     values to set given as a 3x3 tensor
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_aniso_by_value(cs_property_t    *pty,
                                        const char       *zname,
                                        double            tens[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure thanks to a time
 *        function for all cells associated to the zone named zname.
 *        Optimized case with a symmetric storage.
 *
 * \param[in, out] pty       pointer to a cs_property_t structure
 * \param[in]      zname     name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]      symtens   symmetric tensor given as an array of 6 values
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_aniso_sym_by_value(cs_property_t    *pty,
                                   const char       *zname,
                                   cs_real_t         symtens[6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value at the boundary for the given (anisotropic)
 *        property. This value is uniform and steady for all the boundary
 *        faces associated to the boundary zone named zname
 *        Optimized case with a symmetric storage.
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       tens     tensor to set given as an array of 6 values
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_aniso_sym_by_value(cs_property_t    *pty,
                                            const char       *zname,
                                            double            tens[6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to a time function for all
 *        cells associated to the zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_time_func_t function
 * \param[in]       input    null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_time_func(cs_property_t      *pty,
                             const char         *zname,
                             cs_time_func_t     *func,
                             void               *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure at the boundary thanks
 *        to a time function in a subdomain attached to a zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       func     pointer to a cs_time_func_t function
 * \param[in]       input    null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_time_func(cs_property_t      *pty,
                                      const char         *zname,
                                      cs_time_func_t     *func,
                                      void               *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_property_t structure thanks to an analytic
 *         function for all cells associated to the zone named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           cells are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *zname,
                            cs_analytic_func_t   *func,
                            void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure at the boundary thanks
 *        to a time function for all boundary faces associated to the zone
 *        named zname
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       zname    name of the associated zone (if null or "" all
 *                           boundary faces are considered)
 * \param[in]       func     pointer to a cs_analytic_func_t function
 * \param[in]       input    null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_analytic(cs_property_t        *pty,
                                     const char           *zname,
                                     cs_analytic_func_t   *func,
                                     void                 *input);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of a cs_property_t structure thanks to low-level
 *        functions specifying how to evaluate the property in each cell (with
 *        cell-wise structures or not). This definition applies to all cells
 *        associated to the zone named zname
 *
 * \param[in, out] pty           pointer to a cs_property_t structure
 * \param[in]      zname         name of the zone (if null or "" then all cells
 *                               are selected)
 * \param[in]      context       pointer to a structure (may be null)
 * \param[in]      get_eval_at_cell      pointer to a function
 * \param[in]      get_eval_at_cell_cw   pointer to a function
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_func(cs_property_t         *pty,
                        const char            *zname,
                        void                  *context,
                        cs_xdef_eval_t        *get_eval_at_cell,
                        cs_xdef_cw_eval_t     *get_eval_at_cell_cw);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_property_t structure thanks to an array of values. If an
 *        advanced usage of the definition by array is needed, then call \ref
 *        cs_xdef_array_set_adjacency and/or \ref cs_xdef_array_set_sublist
 *
 * \param[in, out] pty           pointer to a cs_property_t structure
 * \param[in]      zname         name of the zone (if null or "" then all cells
 *                               are selected)
 * \param[in]      val_location  information to know where are located values
 * \param[in]      array         pointer to an array
 * \param[in]      is_owner      transfer the lifecycle to the cs_xdef_t struc.
 *                               (true or false)
 * \param[in]      full_length   if true, the size of "array" is allocated to
 *                               the total numbers of entities related to the
 *                               given location. If false, a new list is
 *                               allocated and filled with the related subset
 *                               indirection.
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_array(cs_property_t      *pty,
                         const char         *zname,
                         cs_flag_t           val_location,
                         cs_real_t          *array,
                         bool                is_owner,
                         bool                full_length);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the values of a property at the boundary thanks to an array.
 *        If an advanced usage of the definition by array is needed, then call
 *        \ref cs_xdef_array_set_adjacency and/or \ref
 *        cs_xdef_array_set_sublist
 *
 * \param[in, out] pty          pointer to a cs_property_t structure
 * \param[in]      zname        null or name of the boundary zone
 * \param[in]      val_loc      information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, the size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_array(cs_property_t      *pty,
                                  const char         *zname,
                                  cs_flag_t           val_loc,
                                  cs_real_t          *array,
                                  bool                is_owner,
                                  bool                full_length);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a field structure. One
 *         assumes that all cells are defined using this array.
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_def_by_field(cs_property_t    *pty,
                         cs_field_t       *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the values of a property at the boundary thanks to a field
 *        structure. One assumes that all boundary faces are defined using
 *        this array.
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 *
 * \return a pointer to the resulting cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_property_boundary_def_by_field(cs_property_t    *pty,
                                  cs_field_t       *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the values of the property at cells from the given
 *        definition. According to the parameter "dense_ouput", the "eval"
 *        array should be allocated with a size equal to pty->dim*n_cells
 *        (where "dim" depends on the type of property to handle) when no dense
 *        ouput is requested. Otherwise, an allocation size equal to pty->dim *
 *        the number of cells associated to the definition "def" is enough.
 *
 *        No scaling is applied to the value. This should be done with a higher
 *        level function like \ref cs_property_eval_at_cells or
 *        \ref cs_property_get_cell_tensor
 *
 * \param[in]      pty           pointer to a property structure
 * \param[in]      def_id        id associated to the definition
 * \param[in]      dense_output  true/false
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in, out] eval          array storing the result of the evaluations
 */
/*----------------------------------------------------------------------------*/

void
cs_property_evaluate_def(const cs_property_t    *pty,
                         int                     def_id,
                         bool                    dense_output,
                         double                  t_eval,
                         cs_real_t              *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the values of a property at boundary faces from the given
 *        boundary definition. If a dense ouput is not requested, then the size
 *        of the resulting array should be allocated at least to pty->dim *
 *        n_b_faces. Otherwise, n_b_faces can be replaced by the number of
 *        boundary faces associated to the current definition.
 *
 *        No scaling is applied to the value. This should be done with a higher
 *        level function like \ref cs_property_eval_at_boundary_faces
 *
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      def_id        id associated to the definition
 * \param[in]      dense_output  true/false
 * \param[in]      t_eval        time at which one performs the evaluation
 * \param[in, out] array         array storing the result of the evaluation(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_evaluate_boundary_def(const cs_property_t  *pty,
                                  int                   def_id,
                                  bool                  dense_output,
                                  double                t_eval,
                                  cs_real_t            *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each cell. Store the evaluation
 *        in the given array.
 *
 * \param[in]      t_eval       physical time at which one evaluates the term
 * \param[in]      pty          pointer to a cs_property_t structure
 * \param[out]     pty_stride   = 0 if uniform, =1 otherwise
 * \param[in, out] p_pty_vals   pointer to an array of values. Allocated if not
 *                              The size of the allocation depends on the value
 *                              of the pty_stride
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_get_cell_values(cs_real_t               t_eval,
                                const cs_property_t    *pty,
                                int                    *pty_stride,
                                cs_real_t             **p_pty_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each cell. Store the
 *        evaluation in the given array.
 *
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      pty      pointer to a cs_property_t structure
 * \param[in, out] array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_cells(cs_real_t               t_eval,
                          const cs_property_t    *pty,
                          cs_real_t              *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the value of the property at each boundary face. Store the
 *        result of the evaluation in the given array.
 *
 * \param[in]      t_eval   physical time at which one evaluates the term
 * \param[in]      pty      pointer to a cs_property_t structure
 * \param[in, out] array    pointer to an array of values (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_eval_at_boundary_faces(cs_real_t               t_eval,
                                   const cs_property_t    *pty,
                                   cs_real_t              *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of the tensor attached to a property at the cell
 *        center
 *
 * \param[in]      c_id          id of the current cell
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      do_inversion  true or false
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_cell_tensor(cs_lnum_t               c_id,
                            cs_real_t               t_eval,
                            const cs_property_t    *pty,
                            bool                    do_inversion,
                            cs_real_t               tensor[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the cell center
 *
 * \param[in] c_id     id of the current cell
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_cell_value(cs_lnum_t              c_id,
                           cs_real_t              t_eval,
                           const cs_property_t   *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of the tensor attached to a property at the cell
 *        center
 *        Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm            pointer to a cs_cell_mesh_t structure
 * \param[in]      pty           pointer to a cs_property_t structure
 * \param[in]      t_eval        physical time at which one evaluates the term
 * \param[in]      do_inversion  true or false
 * \param[in, out] tensor        3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_tensor_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_property_t    *pty,
                           cs_real_t               t_eval,
                           bool                    do_inversion,
                           cs_real_t               tensor[3][3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the cell center
 *        Version using a cs_cell_mesh_t structure
 *
 * \param[in] cm        pointer to a cs_cell_mesh_t structure
 * \param[in] pty       pointer to a cs_property_t structure
 * \param[in] t_eval    physical time at which one evaluates the term
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_value_in_cell(const cs_cell_mesh_t   *cm,
                          const cs_property_t    *pty,
                          cs_real_t               t_eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of an isotropic property in each portion of dual
 *        cell in a (primal) cell. This relies on the c2v connectivity.
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty       pointer to a cs_property_t structure
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in, out] eval      array of values storing the evaluation
 */
/*----------------------------------------------------------------------------*/

void
cs_property_c2v_values(const cs_cell_mesh_t   *cm,
                       const cs_property_t    *pty,
                       cs_real_t               t_eval,
                       cs_real_t              *eval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the value of a property at the face center
 *
 * \param[in] f_id     id of the current face
 * \param[in] t_eval   physical time at which one evaluates the term
 * \param[in] pty      pointer to a cs_property_t structure
 *
 * \return the value of the property for the given face
 */
/*----------------------------------------------------------------------------*/

cs_real_t cs_property_get_face_value(cs_lnum_t            f_id,
                                     cs_real_t            t_eval,
                                     const cs_property_t *pty);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Fourier number in each cell
 *
 * \param[in]      pty       pointer to the diffusive property struct.
 * \param[in]      t_eval    physical time at which one evaluates the term
 * \param[in]      dt        value of the current time step
 * \param[in, out] fourier   pointer to an array storing Fourier numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_fourier(const cs_property_t    *pty,
                        cs_real_t               t_eval,
                        double                  dt,
                        cs_real_t               fourier[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print a summary of the settings for all defined cs_property_t
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_property_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PROPERTY_H__ */
