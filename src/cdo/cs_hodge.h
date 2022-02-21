#ifndef __CS_HODGE_H__
#define __CS_HODGE_H__

/*============================================================================
 * Manage discrete Hodge operators and closely related operators
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_param_cdo.h"
#include "cs_property.h"
#include "cs_sdm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/*! \enum cs_hodge_type_t
 *  \brief Type of discrete Hodge operators
 *
 * \var CS_HODGE_TYPE_VPCD
 * Hodge operator from primal vertices to dual cells. This operator is also
 * equivalent to a mass matrix for vertex-based schemes
 *
 * \var CS_HODGE_TYPE_EPFD
 * Hodge operator from primal edges to dual faces
 *
 * \var CS_HODGE_TYPE_FPED
 * Hodge operator from primal faces to dual edges
 *
 * \var CS_HODGE_TYPE_EDFP
 * Hodge operator from dual edges to primal faces
 *
 * \var CS_HODGE_TYPE_CPVD
 * Hodge operator from primal cells to dual vertices
 *
 * \var CS_HODGE_TYPE_FB
 * Hodge operator playing the role of mass matrix for face-based schemes.
 * This operator deals with a hybrid space of degrees of freedom (faces and
 * cells)
 *
 * \var CS_HODGE_TYPE_VC
 * Hodge operator playing the role of mass matrix for vertex+cell-based schemes.
 * This operator deals with a hybrid space of degrees of freedom (vertices and
 * cells)
 */

typedef enum {

  CS_HODGE_TYPE_VPCD,
  CS_HODGE_TYPE_EPFD,
  CS_HODGE_TYPE_FPED,
  CS_HODGE_TYPE_EDFP,
  CS_HODGE_TYPE_CPVD,
  CS_HODGE_TYPE_FB,
  CS_HODGE_TYPE_VC,

  CS_HODGE_N_TYPES

} cs_hodge_type_t;

/*! \enum cs_hodge_algo_t
 *  \brief Type of algorithm to build a discrete Hodge operator
 *
 * \var CS_HODGE_ALGO_VORONOI
 * This algorithm assumes that an orthogonality condition holds between
 * geometrical entities (e.g. the face normal/dual edge or the primal edge/dual
 * face normal). This algorithm leads to a diagonal discrete Hodge operator
 * (and is related to a two-point flux approximation). Be aware that using this
 * technique on a non-orthogonal mesh leads to a consistency error which does
 * not decrease when refining the mesh.
 *
 * \var CS_HODGE_ALGO_WBS
 * WBS means "Whitney Barycentric Subdivision".
 * This algorithm relies on a subdivision into tetrahedra of each polyhedral
 * cell and then applies algorithms on this subdivision shape functions close
 * to what exists in the Finite Element litterature
 *
 * \var CS_HODGE_ALGO_COST
 * Orthogonal decomposition between the consistency (CO) and the stabilization
 * (ST) parts. Several discretizations share this way to build discrete Hodge
 * operators like SUSHI, DGA and GCR (these discretizations only differ from
 * the scaling coefficient in front of the stabilization term)
 *
 * \var CS_HODGE_ALGO_OCS2
 * This algorithm is close to \ref CS_HODGE_ALGO_COST but relies on a
 * subdivision of each polyhedral cell corresponding to the refinement
 * considered in \ref CS_HODGE_ALGO_COST for the building of the stabilization
 * term
 *
 * \var CS_HODGE_ALGO_BUBBLE
 * This algorithm also relies on an orthogonal decomposition between the
 * consistency (CO) and the stabilization (ST) parts but the stabilization part
 * relies on a bubble function (similar to what can be found in the Finite
 * Element literature)
 *
 * \var CS_HODGE_ALGO_AUTO
 * Automatic switch between the above-mentioned algorithms (not used for the
 * moment).
 */

typedef enum {

  CS_HODGE_ALGO_VORONOI,
  CS_HODGE_ALGO_WBS,
  CS_HODGE_ALGO_COST,
  CS_HODGE_ALGO_OCS2,
  CS_HODGE_ALGO_BUBBLE,
  CS_HODGE_ALGO_AUTO,

  CS_HODGE_N_ALGOS

} cs_hodge_algo_t;

/*!
 * \struct cs_hodge_param_t
 * \brief Structure storing all metadata/parameters related to the usage of a
 *        discrete Hodge operator
 */

typedef struct {

  bool              inv_pty;

  /*!< \brief Inversion of the property evaluation or not
   *
   * Each Hodge operator is associated to a property. Since either the value or
   * its reciprocal is considered, this parameter helps one to know in which
   * situation we are.  If this is set to true, then one needs to consider the
   * reciprocal of the associated property */

  cs_hodge_type_t   type;   /*!< Type of discrete Hodge operator */
  cs_hodge_algo_t   algo;   /*!< Type of algorithm used to build the operator */

  double            coef;

  /*!< \brief Scaling coefficient value
   *
   * Value of the coefficient scaling the stabilization part if the COST or
   * OCS2 algorithm is used. Otherwise the value is set to 0 and ignored.
   */

} cs_hodge_param_t;

/* DISCRETE HODGE OPERATORS */
/* ======================== */

/*!
 * \struct cs_hodge_t
 * \brief Structure associated to a discrete Hodge operator *
 */

typedef struct {

  const cs_hodge_param_t  *param;      /*!< Set of parameters (shared) */

  cs_property_data_t      *pty_data;

  /*!< \brief Data evaluation
   *
   * Each Hodge operator is associated to a property. Pointer to a structure
   * storing the evaluation of the associated property
   */

  cs_sdm_t                *matrix;    /*!< Matrix storing operator values  */

} cs_hodge_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a discrete Hodge operator or a related operator (such as the
 *         stiffmess matrix) for a given cell.
 *         The discrete Hodge operator is stored in hodge->matrix whereas the
 *         associated operator is stored in cb->loc
 *         One checks if something has to be computed. One skips the
 *         computation if the value of the associated property is equal to
 *         zero.
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] hodge     pointer to a cs_hodge_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise.
 */
/*----------------------------------------------------------------------------*/

typedef bool
(cs_hodge_compute_t)(const cs_cell_mesh_t   *cm,
                     cs_hodge_t             *hodge,
                     cs_cell_builder_t      *cb);

/*============================================================================
 * Static inline public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check if two sets of parameters related to how build a discrete
 *          Hodge operator are similar.
 *
 * \param[in]  h1_info     pointer to a first cs_hodge_param_t structure
 * \param[in]  h2_info     pointer to a second cs_hodge_param_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_hodge_param_is_similar(const cs_hodge_param_t    h1_info,
                          const cs_hodge_param_t    h2_info)
{
  if (h1_info.type != h2_info.type)
    return false;
  if (h1_info.algo != h2_info.algo)
    return false;
  if (h1_info.algo == CS_HODGE_ALGO_COST ||
      h1_info.algo == CS_HODGE_ALGO_BUBBLE) {
    if (fabs(h1_info.coef - h2_info.coef) > 0)
      return false;
    else
      return true;
  }
  else
    return true;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a pointer to a cs_hodge_t structure
 *
 * \param[in] connect        pointer to cs_cdo_connect_t structure
 * \param[in] property       pointer to a property structure
 * \param[in] hp             pointer to a cs_hodge_param_t structure
 * \param[in] need_tensor    true if one needs a tensor otherwise false
 * \param[in] need_eigen     true if one needs to compute eigen valuese
 *
 * \return a pointer to the new allocated cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_t *
cs_hodge_create(const cs_cdo_connect_t   *connect,
                const cs_property_t      *property,
                const cs_hodge_param_t   *hp,
                bool                      need_tensor,
                bool                      need_eigen);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize an array of pointers to a cs_hodge_t
 *         structures. This array is of size the number of OpenMP threads.
 *         Only the one associated to the current thread is set.
 *
 * \param[in] connect        pointer to cs_cdo_connect_t structure
 * \param[in] property       pointer to a property structure
 * \param[in] hp             pointer to a cs_hodge_param_t structure
 * \param[in] need_tensor    true if one needs a tensor otherwise false
 * \param[in] need_eigen     true if one needs to compute eigen valuese
 *
 * \return an array of pointers of cs_hodge_t structures
 */
/*----------------------------------------------------------------------------*/

cs_hodge_t **
cs_hodge_init_context(const cs_cdo_connect_t   *connect,
                      const cs_property_t      *property,
                      const cs_hodge_param_t   *hp,
                      bool                      need_tensor,
                      bool                      need_eigen);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_hodge_t structure
 *
 * \param[in, out] p_hodge    double pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_free(cs_hodge_t    **p_hodge);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a set of cs_hodge_t structures
 *
 * \param[in, out] p_hodges    triple pointer to cs_hodge_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_free_context(cs_hodge_t    ***p_hodges);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a function pointer to compute a discrete Hodge operator
 *
 * \param[in] calling_func    name of the calling function
 * \param[in] hp              a cs_hodge_param_t structure
 *
 * \return a pointer to the corresponding function
 */
/*----------------------------------------------------------------------------*/

cs_hodge_compute_t *
cs_hodge_get_func(const char               *calling_func,
                  const cs_hodge_param_t    hp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the consistency of the settings between terms related to a
 *         mass matrix and define the common algorithm to use.
 *         If a term should not be considered, set the algorithm to
 *         CS_HODGE_N_ALGOS
 *
 * \param[in] eqname     name of the equation to check
 * \param[in] reac_algo  optional algo. used for the reaction term
 * \param[in] time_algo  optional algo. used for the unsteady term
 * \param[in] srct_algo  optional algo. used for the source term
 *
 * \return the common algorithm to use
 */
/*----------------------------------------------------------------------------*/

cs_hodge_algo_t
cs_hodge_set_mass_algo(const char         *eqname,
                       cs_hodge_algo_t     reac_algo,
                       cs_hodge_algo_t     time_algo,
                       cs_hodge_algo_t     srct_algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a cs_hodge_param_t structure
 *
 * \param[in] prefix    optional string
 * \param[in] property  optional pointer to a property structure
 * \param[in] hp        a cs_hodge_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_param_log(const char               *prefix,
                   const cs_property_t      *property,
                   const cs_hodge_param_t    hp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the set of parameters associated to a discrete Hodge operator
 *         to another one
 *
 * \param[in]       h_ref   reference set of parameters
 * \param[in, out]  h_cpy   set of parameters to update
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_copy_parameters(const cs_hodge_param_t   *h_ref,
                         cs_hodge_param_t         *h_cpy);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the property value (scalar- or tensor-valued) related to a
 *         discrete Hodge operator inside a cell and if needed other related
 *         quantities
 *
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_set_property_value(const cs_lnum_t       c_id,
                            const cs_real_t       t_eval,
                            const cs_flag_t       c_flag,
                            cs_hodge_t           *hodge);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the property value (scalar- or tensor-valued) related to a
 *         discrete Hodge operator inside a cell and if needed ohter related
 *         quantities.
 *         Cell-wise variant (usage of cs_cell_mesh_t structure)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_set_property_value_cw(const cs_cell_mesh_t   *cm,
                               const cs_real_t         t_eval,
                               const cs_flag_t         c_flag,
                               cs_hodge_t             *hodge);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_cost_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          with the usage of bubble stabilization.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_bubble_get_stiffness(const cs_cell_mesh_t    *cm,
                                 cs_hodge_t              *hodge,
                                 cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_voro_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and an isotropic property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_iso_stiffness(const cs_cell_mesh_t   *cm,
                                   cs_hodge_t             *hodge,
                                   cs_cell_builder_t      *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and an anistropic property
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge    pointer to a cs_hodge_t structure
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_aniso_stiffness(const cs_cell_mesh_t    *cm,
                                     cs_hodge_t              *hodge,
                                     cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and isotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_bubble_get_iso_stiffness(const cs_cell_mesh_t    *cm,
                                     cs_hodge_t              *hodge,
                                     cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and anisotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_bubble_get_aniso_stiffness(const cs_cell_mesh_t    *cm,
                                       cs_hodge_t              *hodge,
                                       cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case Vb schemes and an anisotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_ocs2_get_aniso_stiffness(const cs_cell_mesh_t     *cm,
                                     cs_hodge_t               *hodge,
                                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_voro_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS)
 *          algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_wbs_get_stiffness(const cs_cell_mesh_t     *cm,
                              cs_hodge_t               *hodge,
                              cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS) algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_get_stiffness(const cs_cell_mesh_t     *cm,
                           cs_hodge_t               *hodge,
                           cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator on a given cell which is equivalent of
 *         a mass matrix. It relies on a CO+ST algo. and is specific to CDO-Fb
 *         schemes.
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_get(const cs_cell_mesh_t     *cm,
                cs_hodge_t               *hodge,
                cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Voronoi
 *          algo. This leads to a diagonal operator.
 *          This function is specific for vertex+cell-based schemes
 *          The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_voro_get(const cs_cell_mesh_t     *cm,
                      cs_hodge_t               *hodge,
                      cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the WBS algo.
 *          This function is specific for vertex+cell-based schemes
 *          The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_wbs_get(const cs_cell_mesh_t     *cm,
                     cs_hodge_t               *hodge,
                     cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator for a given cell using WBS algo.
 *         Hodge op. from primal vertices to dual cells.
 *         This function is specific for vertex-based schemes
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vpcd_wbs_get(const cs_cell_mesh_t    *cm,
                      cs_hodge_t              *hodge,
                      cs_cell_builder_t       *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator for a given cell using VORONOI algo.
 *         Hodge op. from primal vertices to dual cells.
 *         This function is specific for vertex-based schemes
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vpcd_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal edges to dual faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal edges to dual faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          with a bubble stabilization.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          Hodge op. from primal edges to dual faces. This function is
 *          specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_ocs2_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] hodge     pointer to a cs_hodge_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] hodge     pointer to a cs_hodge_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_cost_get_opt(const cs_cell_mesh_t     *cm,
                           cs_hodge_t               *hodge,
                           cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator and multiple it with
 *          a vector
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      hodgep      cs_hodge_param_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      in_vals   vector to multiply with the discrete Hodge op.
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] result    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_matvec(const cs_cdo_connect_t       *connect,
                const cs_cdo_quantities_t    *quant,
                const cs_hodge_param_t        hodgep,
                const cs_property_t          *pty,
                const cs_real_t               in_vals[],
                cs_real_t                     t_eval,
                cs_real_t                     result[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator in order to define
 *          a circulation array from a flux array
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      hodgep      cs_hodge_param_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      flux      vector to multiply with the discrete Hodge op.
 * \param[in, out] circul    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_circulation_from_flux(const cs_cdo_connect_t       *connect,
                               const cs_cdo_quantities_t    *quant,
                               cs_real_t                     t_eval,
                               const cs_hodge_param_t        hodgep,
                               const cs_property_t          *pty,
                               const cs_real_t               flux[],
                               cs_real_t                     circul[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the hodge operator related to a face (i.e. a mass matrix
 *          with unity property) using a Whitney Barycentric Subdivision (WBS)
 *          algorithm
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] hf        pointer to a cs_sdm_t structure to define
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_compute_wbs_surfacic(const cs_face_mesh_t    *fm,
                              cs_sdm_t                *hf);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HODGE_H__ */
