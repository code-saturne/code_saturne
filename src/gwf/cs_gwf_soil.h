#ifndef __CS_GWF_SOIL_H__
#define __CS_GWF_SOIL_H__

/*============================================================================
 * Set of main functions to handle soils in the groundwater flow module
 * when using CDO schemes
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

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_gwf_param.h"
#include "cs_gwf_hydraulic_model.h"
#include "cs_mesh.h"
#include "cs_property.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function pointer prototypes
 *============================================================================*/

typedef struct _gwf_soil_vgm_tpf_param_t  cs_gwf_soil_vgm_tpf_param_t;
typedef struct _gwf_soil_t cs_gwf_soil_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to update the properties related to a hydraulic
 *         model given the soil model. The soil parameters depend on the type
 *         of soil model.
 *
 * \param[in]      t_eval         time at which one performs the evaluation
 * \param[in]      mesh           pointer to a cs_mesh_t structure
 * \param[in]      connect        pointer to a cs_cdo_connect_t structure
 * \param[in]      quant          pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone           pointer to a cs_zone_t
 * \param[in, out] soil           pointer to the soil structure to update
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_update_t) (const cs_real_t              t_eval,
                        const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_zone_t             *zone,
                        cs_gwf_soil_t               *soil);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function to set free the parameter structure associated to
 *         a soil
 *
 * \param[in, out] p_param     double pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_free_param_t)(void         **p_param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of the different properties related to a soil in
 *        the case of a Van Genuchten-Mualem model and a two-phase flow model.
 *
 * \param[in]  sp        set of modelling parameters
 * \param[out] sl        liquid saturation
 * \param[out] dsldpc    liquid capacity
 * \param[out] krl       relative permeability for the liquid phase
 * \param[out] krg       relative permeability for the gas phase
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_gwf_soil_tpf_eval_t)(const cs_gwf_soil_vgm_tpf_param_t    *sp,
                         const double                          pc,
                         double                               *sl,
                         double                               *dsldpc,
                         double                               *krl,
                         double                               *krg);

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_gwf_soil_join_type_t
 *  \brief Kind of joining function used with closure laws
 *
 * \var CS_GWF_SOIL_JOIN_NOTHING
 *      No joining function
 *
 * \var CS_GWF_SOIL_JOIN_C1_HYPERBOLIC
 *      C1 join using a hyperbolic function
 *
 * \var CS_GWF_SOIL_JOIN_C1_EXPONENTIAL
 *      C1 join using an exponential function
 *
 * \var CS_GWF_SOIL_JOIN_C1_POLY_ORDER2
 *      C1 join using a second order polynomial
 *
 * \var CS_GWF_SOIL_N_JOINS
 */

typedef enum {

  CS_GWF_SOIL_JOIN_NOTHING,
  CS_GWF_SOIL_JOIN_C1_HYPERBOLIC,
  CS_GWF_SOIL_JOIN_C1_EXPONENTIAL,
  CS_GWF_SOIL_JOIN_C1_POLY_ORDER2,

  CS_GWF_SOIL_N_JOINS

} cs_gwf_soil_join_type_t;

/*! \enum cs_gwf_soil_state_t
 *  \brief Kind of state in which a cell is
 *
 * \var CS_GWF_SOIL_STATE_SATURATED
 *      All the available space in the porous media is liquid
 *
 * \var CS_GWF_SOIL_STATE_UNSATURATED
 *      Only a part of the porous media is filled with liquid
 *
 * \var CS_GWF_SOIL_STATE_DRY
 *      All the available space in the porous media is gas or void
 *
 * \var CS_GWF_SOIL_N_STATES
 */

typedef enum {

  CS_GWF_SOIL_STATE_SATURATED    = 0,
  CS_GWF_SOIL_STATE_UNSATURATED  = 1,
  CS_GWF_SOIL_STATE_DRY          = 2,

  CS_GWF_SOIL_N_STATES           = 3,

} cs_gwf_soil_state_t;

/*! \struct cs_gwf_soil_vgm_spf_param_t
 *
 * \brief Structure to handle the Van Genuchten-Mualem model of soil in the
 *        case of a single-phase flow in a porous media
 *
 *        This structure stores the parameters defining the evolution laws for
 *        the liquid saturation and the relative permeability.
 */

typedef struct {

  /*!
   * \var residual_moisture
   *      Also called residual liquid saturation
   *
   * \var n
   *      Shape parameter. Should be 1.25 < n < 6
   *
   * \var m
   *      Value depending of that of n
   *      m = 1 - 1/n
   *
   * \var scale
   *      Value of a scaling parameter in [m^-1]
   *
   * \var tortuosity
   *      Tortuosity parameter for saturated hydraulic conductivity
   */

  double             residual_moisture;

  double             n;
  double             m;
  double             scale;
  double             tortuosity;

} cs_gwf_soil_vgm_spf_param_t;

/*! \struct cs_gwf_soil_vgm_tpf_param_t
 *
 * \brief Structure to handle the Van Genuchten-Mualem model of soil in the
 *        case of a two-phase flow in a porous media
 *
 *        This structure stores the parameters defining the evolution laws for
 *        the liquid saturation and the relative permeabilities (in the liquid
 *        and gas phases).
 */

struct _gwf_soil_vgm_tpf_param_t {

  /*!
   * \var n
   *      Shape parameter. This value should be strictly greater than 1.0.
   *
   * \var m
   *      Value depending of that of n
   *      Derived quantity: m = 1 - 1/n
   *
   * \var inv_m
   *      Derived quantity which is the reciprocal of m. This value is given
   *      by the identity 1/m = 1 + 1/(n-1)
   *
   * \var pr_r
   *      Reference (capillarity) pressure
   *
   * \var inv_pr_r
   *      Derived quantity: Reciprocal of pr_r
   *
   * \var sl_r
   *      Residual liquid saturation
   *
   * \var sl_s
   *      Saturated (i.e. maximum) liquid saturation
   *
   * \var sl_range
   *      Derived quantity: sl_s - sl_r
   */

  double       n;
  double       m;
  double       inv_m;
  double       pr_r;
  double       inv_pr_r;
  double       sl_r;
  double       sl_s;
  double       sl_range;

  /*!
   * Parameters to handle a joining function
   *
   * \var sle_jtype
   *      type of joining function to consider for the Sle(Pc) curve
   *
   * \var kr_jtype
   *      type of joining function to consider for the krg(Sl) and krl(Sl)
   *      curves
   *
   * \var sle_thres
   *      Value above which the suction law is replaced with a joining function
   *      (for instance sle_thres = 0.999 is the default value). If the value
   *      is greater or equal than 1.0, there is no polynomial joining.
   *
   * \var eval_properties
   *      function performing the evaluation of the soil laws with/without a
   *      joining
   *
   * \var pc_star
   *      capillarity pressure related to the value of sle_thres
   *
   * \var dsldpc_star
   *      derivative of the liquid saturation with respect to the capillarity
   *      pressure at pc_star
   *
   * \var sle_alpha
   *      optional pre-computed coefficient when a joining function is used
   *      for the effective liquid saturation
   *
   * \var sle_beta
   *      optional pre-computed coefficient when a joining function is used
   *      for the effective liquid saturation
   *
   * \var krg_star
   *      relative permeability in the gas phase for the value sle_thres
   *
   * \var dkrgdsl_star
   *      derivative of the relative permeability in the gas with respect to
   *      the liquid saturation at sle_thres
   *
   * \var krg_alpha
   *      pre-computed coefficient when a joining function is used for the
   *      relative permeability in the gaz
   *
   * \var krl_star
   *      relative permeability in the liquid phase for the value sle_thres
   *
   * \var dkrldsl_star
   *      derivative of the relative permeability in the liquid with respect to
   *      the liquid saturation at sle_thres
   *
   * \var krl_alpha
   *      pre-computed coefficient when a joining function is used for the
   *      relative permeability in the liquid
   */

  cs_gwf_soil_join_type_t    sle_jtype;
  cs_gwf_soil_join_type_t    kr_jtype;
  double                     sle_thres;

  /* Function pointer */

  cs_gwf_soil_tpf_eval_t    *eval_properties;

  /* Derived quantities */

  double                     pc_star;
  double                     dsldpc_star;
  double                     sle_alpha;
  double                     sle_beta;

  double                     krg_star;
  double                     dkrgdsl_star;
  double                     krg_alpha;

  double                     krl_star;
  double                     dkrldsl_star;
  double                     krl_alpha;

};

typedef struct _gwf_soil_vgm_tpf_param_t cs_gwf_soil_vgm_tpf_param_t;

/*! \struct _gwf_soil_t
 *
 * \brief Main structure to handle a soil in the groundwater flow module: its
 *        definition, the way to update its related properties
 *
 *        Store a set of parameters and pointers describing a soil and its
 *        related hydraulic model (this hydraulic model is shared with the main
 *        structure \ref cs_gwf_t)
 */

struct _gwf_soil_t {

  /*!
   * @name Metadata
   * @{

   * \var id
   *      id associated to a soil. Position in the array of soils.
   *
   * \var zone_id
   *      id related to a volumic cs_zone_t structure (based on cells)
   *
   * \var hydraulic_model
   *      Type of model use in the groundwater flow module to describe the
   *      hydraulic (see \ref cs_gwf_model_type_t for more details)
   *
   * \var hydraulic_context
   * Structure cast on-the-fly. This structure contains parameters, arrays,
   * properties and fields describing the hydraulic state. It depends on the
   * type of hydraulic model which is considered.
   *
   * @}
   * @name Soil features (whatever is the soil model)
   * @{
   *
   * \var bulk_density
   * Value of the mass density of the soil
   *
   * \var porosity
   *      Max. portion of volume in a soil where the liquid (or the gas) can
   *      be.  In a single-phase saturated model this corresponds to the
   *      saturated moisture or the max. liquid saturation.
   *
   * \var abs_permeability_dim
   *      =1 if isotropic or =3 if orthotropic or =9 if anisotropic
   *
   * \var abs_permeability
   *      Value of the intrisic permeability in the soil when all the porous
   *      media is filled with water (= saturated permeability)
   *
   * @}
   * @name Soil modelling
   * @{
   *
   * \var model
   *      Type of model describing the hydraulic behaviour of a soil (cf. \ref
   *      cs_gwf_soil_model_t for more details)
   *
   * \var model_param
   *      Pointer to a structure cast on-the-fly (it depends on the type of
   *      soil model). This structure contains the set of parameters describing
   *      a soil. Can be set to null in the case of a saturated single-phase
   *      model.
   *
   * @}
   * @name Function pointers
   * @{
   *
   * \var update_properties
   *      Pointer to a function which manages the update of the properties
   *      describing the porous media or used in associated equations
   *      (diffusion terms for instance). These functions depend on the model
   *      of soil and the type of model used in the groundwater flow
   *      module. May be set to null if there is no need to update soil
   *      properties.
   *
   * \var free_model_param
   *      Pointer to a function which free the param structure if needed. May
   *      be set to null if there is nothing to free inside the param
   *      structure.
   *
   * @}
   */

  int                             id;
  int                             zone_id;

  cs_gwf_model_type_t             hydraulic_model;
  void                           *hydraulic_context;


  double                          bulk_density;
  double                          porosity;
  int                             abs_permeability_dim;
  double                          abs_permeability[3][3];

  cs_gwf_soil_model_t             model;
  void                           *model_param;

  /* Pointers to functions */

  cs_gwf_soil_update_t           *update_properties;
  cs_gwf_soil_free_param_t       *free_model_param;

};

/*============================================================================
 * User-defined function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the number of allocated soils
 *
 * \return the number of allocated soils
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_get_n_soils(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_gwf_soil_t structure and add it to the array of
 *        soils. An initialization by default of all members is performed.
 *
 * \param[in] zone                pointer to a volume zone structure
 * \param[in] hydraulic_model     main hydraulic model for the module
 * \param[in] model               type of model for the soil behavior
 * \param[in] perm_type           type of permeability (iso/anisotropic)
 * \param[in] k_abs               absolute (intrisic) permeability
 * \param[in] porosity            porosity or max. moisture content
 * \param[in] bulk_density        value of the mass density
 * \param[in] hydraulic_context   pointer to the context structure
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_create(const cs_zone_t                 *zone,
                   cs_gwf_model_type_t              hydraulic_model,
                   cs_gwf_soil_model_t              model,
                   cs_property_type_t               perm_type,
                   double                           k_abs[3][3],
                   double                           porosity,
                   double                           bulk_density,
                   void                            *hydraulic_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its id
 *
 * \param[in]  id      id to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_id(int   id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a soil structure from its name
 *
 * \param[in]  name      name to look for
 *
 * \return a pointer to a cs_gwf_soil_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_soil_by_name(const char    *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve a zone associated to a soil from its id
 *
 * \param[in] soil_id      id to look for
 *
 * \return a pointer to a zone structure or null
 */
/*----------------------------------------------------------------------------*/

const cs_zone_t *
cs_gwf_soil_get_zone(int   soil_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_free_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last initialization step for the soil structures/parameters
 *
 * \param[in] gwf_model      modelling used for the GWF module
 * \param[in] post_flag      which post-processing to do
 * \param[in] n_cells        number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_finalize_setup(cs_gwf_model_type_t    gwf_model,
                           cs_flag_t              post_flag,
                           cs_lnum_t              n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the settings related to all cs_gwf_soil_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the soil properties
 *
 * \param[in] time_eval      time at which one evaluates properties
 * \param[in] mesh           pointer to the mesh structure
 * \param[in] connect        pointer to the cdo connectivity
 * \param[in] cdoq           pointer to the cdo quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_update(cs_real_t                     time_eval,
                   const cs_mesh_t              *mesh,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the soil state associated to each cell w.r.t. the given
 *        liquid saturation
 *
 * \param[in] n_cells      number of mesh cells
 * \param[in] sliq         values of the liquid saturation in each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_update_soil_state(cs_lnum_t            n_cells,
                              const cs_real_t     *sliq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the soil porosity and absolute permeability
 *        (which are properties always defined in the GWF module). One relies
 *        on the definition of these properties in each soil.
 *
 * \param[in, out] abs_permeability    pointer to a cs_property_t structure
 * \param[in, out] soil_porosity       pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_define_shared_properties(cs_property_t   *abs_permeability,
                                     cs_property_t   *soil_porosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the definition of the soil porosity and absolute porosity (which
 *        are properties always defined). This relies on the definition of
 *        each soil.
 *
 * \param[in, out] moisture_content   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_define_sspf_property(cs_property_t   *moisture_content);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build an array storing the dual volume associated to each vertex
 *        taking into account the porosity of the soil
 *        The computed quantity is stored as a static array. Use the function
 *        cs_gwf_soil_get_dual_vol_l()
 *
 * \param[in] cdoq     additional geometrical quantities for CDO schemes
 * \param[in] connect  additional connectivities for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_build_dual_porous_volume(const cs_cdo_quantities_t    *cdoq,
                                     const cs_cdo_connect_t       *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the dual volume weighted by the soil porosity
 *        Array of size n_vertices
 *
 * \return a pointer to the requested array
 */
/*----------------------------------------------------------------------------*/

const double *
cs_gwf_soil_get_dual_porous_volume(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the associated soil for each cell
 *
 * \return a pointer to the array
 */
/*----------------------------------------------------------------------------*/

const short int *
cs_gwf_soil_get_cell2soil(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the array storing the soil state associated to each cell
 *
 * \return a pointer to the array (may be null)
 */
/*----------------------------------------------------------------------------*/

const int *
cs_gwf_soil_get_soil_state(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the porosity value for the given soil id
 *
 * \param[in] soil_id      id of the requested soil
 *
 * \return the value of the soil porosity
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_porosity(int   soil_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the saturated moisture for the given soil id
 *
 * \param[in]  soil_id     id of the requested soil
 *
 * \return the value of the saturated moisture
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_gwf_soil_get_saturated_moisture(int   soil_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the max dim (aniso=9; iso=1) for the absolute permeability
 *         associated to each soil
 *
 * \return the associated max. dimension
 */
/*----------------------------------------------------------------------------*/

int
cs_gwf_soil_get_permeability_max_dim(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a soil defined by a Van Genuchten-Mualem model in the case of
 *        single-phase flow in an (unsaturated) porous media
 *
 *        The (effective) liquid saturation (also called moisture content)
 *        follows the identity
 *        S_l,eff = (S_l - theta_r)/(theta_s - theta_r)
 *                = (1 + |alpha . h|^n)^(-m)
 *
 *        The isotropic relative permeability is defined as:
 *        k_r = S_l,eff^L * (1 - (1 - S_l,eff^(1/m))^m))^2
 *        where m = 1 -  1/n
 *
 * \param[in, out] soil       pointer to a cs_gwf_soil_t structure
 * \param[in]      theta_r    residual moisture
 * \param[in]      alpha      scale parameter (in m^-1)
 * \param[in]      n          shape parameter
 * \param[in]      L          tortuosity parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_spf_param(cs_gwf_soil_t         *soil,
                              double                 theta_r,
                              double                 alpha,
                              double                 n,
                              double                 L);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters related to a Van Genuchten-Mualen model to defined
 *        the behavior of a soil in the case of two-phase flow in an porous
 *        media
 *
 *        The (effective) liquid saturation follows the identity
 *        sl_eff = (sl - sl_r)/(sl_s - sl_r)
 *                = (1 + |Pc/Pr_r|^n)^(-m)
 *        where m = 1 -  1/n
 *
 *        The isotropic relative permeability in the liquid and gaz are defined
 *        as:
 *        krl = sl_eff^(1/2) * (1 - (1 - sl_eff^(1/m))^m))^2
 *        krg = (1 - sl_eff)^(1/2) * (1 - sl_eff^(1/m))^(2m)
 *
 * \param[in, out] soil         pointer to a cs_gwf_soil_t structure
 * \param[in]      n            shape parameter
 * \param[in]      pr_r         reference (capillarity) pressure
 * \param[in]      sl_r         residual liquid saturation
 * \param[in]      sl_s         saturated (max.) liquid saturation
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_tpf_param(cs_gwf_soil_t         *soil,
                              double                 n,
                              double                 pr_r,
                              double                 sl_r,
                              double                 sl_s);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set advanced parameter settings related to a Van Genuchten-Mualen
 *        soil model
 *
 * \param[in, out] soil        pointer to a cs_gwf_soil_t structure
 * \param[in]      sle_jtype   type of joining function for the effective Sl
 * \param[in]      kr_jtype    type of joining function for krg and krl
 * \param[in]      sle_thres   value of the effective liquid saturation above
 *                             which a joining function is used
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_vgm_tpf_advanced_param(cs_gwf_soil_t             *soil,
                                       cs_gwf_soil_join_type_t    sle_jtype,
                                       cs_gwf_soil_join_type_t    kr_jtype,
                                       double                     sle_thres);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a soil defined by a user-defined model
 *
 * \param[in, out] soil              pointer to a cs_gwf_soil_t structure
 * \param[in]      param             pointer to a structure cast on-the-fly
 * \param[in]      update_func       function pointer to update propoerties
 * \param[in]      free_param_func   function pointer to free the param struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_user_model_param(cs_gwf_soil_t               *soil,
                                 void                        *param,
                                 cs_gwf_soil_update_t        *update_func,
                                 cs_gwf_soil_free_param_t    *free_param_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the value of the polynomial order considered when regularizing
 *        the Van Genuchten-Mualen soil law near the saturated regime. Advanced
 *        usage. This function has to be called before calling the function
 *        \ref cs_gwf_soil_set_vgm_tpf_param
 *        Default: 2. Available values: 2 or 3.
 *
 * \param[in] order       value of the polynomial order (2 or 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_soil_set_joining_poly_order(int    order);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_SOIL_H__ */
