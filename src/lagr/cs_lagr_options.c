/*============================================================================
 * Lagrangian module options setting
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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Functions dealing with Lagrangian module options
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_gui_particles.h"
#include "cs_gui_util.h"
#include "cs_mesh_location.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_physical_model.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_post.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_options.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create source term fields for lagrangian module
 *
 * \param[in]  name  source term field name
 * \param[in]  name  source term field dimension
 */
/*----------------------------------------------------------------------------*/

static void
_define_st_field(const char  *name,
                 int          dim)
{
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  int location_id = CS_MESH_LOCATION_CELLS;

  cs_field_create(name,
                  field_type,
                  location_id,
                  dim,
                  false);
}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the boundary variable names array
 *
 * parameters:
 *   ipp     <-- index from the fortran array associated to varname
 *   varname <-- name or label of the variable/scalar/property
 *----------------------------------------------------------------------------*/

static void
_copy_boundary_varname(int          ipp,
                       const char  *varname)
{
  size_t  l;
  assert(ipp >= 0);

  int nvplmx = 50+4*cs_glob_lagr_const_dim->nlayer;

  if (cs_glob_lagr_boundary_interactions->nombrd == NULL) {

    BFT_MALLOC(cs_glob_lagr_boundary_interactions->nombrd,
               nvplmx,
               char *);
    for (int i = 0; i < nvplmx; i++)
      cs_glob_lagr_boundary_interactions->nombrd[i] = NULL;
  }

  l = strlen(varname);

  BFT_REALLOC(cs_glob_lagr_boundary_interactions->nombrd[ipp], l + 1, char);

  strcpy(cs_glob_lagr_boundary_interactions->nombrd[ipp], varname);
}

/*----------------------------------------------------------------------------
 * Initialize Encrustation pointers.
 *----------------------------------------------------------------------------*/

static void
_init_lagr_encrustation_pointers(void)
{
  if (cs_glob_lagr_encrustation->enc1 == NULL)
    BFT_MALLOC(cs_glob_lagr_encrustation->enc1,
               cs_glob_lagr_const_dim->ncharm2,
               cs_real_t);
  if (cs_glob_lagr_encrustation->enc2 == NULL)
    BFT_MALLOC(cs_glob_lagr_encrustation->enc2,
               cs_glob_lagr_const_dim->ncharm2,
               cs_real_t);
  if (cs_glob_lagr_encrustation->tprenc == NULL)
    BFT_MALLOC(cs_glob_lagr_encrustation->tprenc,
               cs_glob_lagr_const_dim->ncharm2,
               cs_real_t);
  if (cs_glob_lagr_encrustation->visref == NULL)
    BFT_MALLOC(cs_glob_lagr_encrustation->visref,
               cs_glob_lagr_const_dim->ncharm2,
               cs_real_t);

  for (int icha = 0; icha < cs_glob_lagr_const_dim->ncharm2; icha++) {

    cs_glob_lagr_encrustation->tprenc[icha] = -999.0;
    cs_glob_lagr_encrustation->visref[icha] = -999.0;
    cs_glob_lagr_encrustation->enc1[icha] = -999.0;
    cs_glob_lagr_encrustation->enc2[icha] = -999.0;

  }
}

/*----------------------------------------------------------------------------
 * Free encrustation pointers.
 *----------------------------------------------------------------------------*/

static void
_free_lagr_encrustation_pointers(void)
{
  BFT_FREE(cs_glob_lagr_encrustation->enc1);
  BFT_FREE(cs_glob_lagr_encrustation->enc2);
  BFT_FREE(cs_glob_lagr_encrustation->tprenc);
  BFT_FREE(cs_glob_lagr_encrustation->visref);
}

/*----------------------------------------------------------------------------
 * Initialize boundary interaction pointers.
 *----------------------------------------------------------------------------*/

static void
_init_lagr_boundary_interaction_pointers(void)
{
  if (cs_glob_lagr_boundary_interactions->iusb == NULL)
    BFT_MALLOC(cs_glob_lagr_boundary_interactions->iusb,
               cs_glob_lagr_boundary_interactions->nusbor, int);

  if (cs_glob_lagr_boundary_interactions->imoybr == NULL)
    BFT_MALLOC(cs_glob_lagr_boundary_interactions->imoybr,
               cs_glob_lagr_const_dim->nusbrd + 10,
               int);
}

/*----------------------------------------------------------------------------
 * Free boundary interaction pointers.
 *----------------------------------------------------------------------------*/

static void
_free_lagr_boundary_interaction_pointers(void)
{
  BFT_FREE(cs_glob_lagr_boundary_interactions->iusb);
  BFT_FREE(cs_glob_lagr_boundary_interactions->imoybr);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/* ---------------------------------------------------------------------- */
/*!
 * \brief Lagrangian module options definition.
 *
 * - default initialization
 * - read user settings
 * - check settings coherency
 * - initialize some structures relative to Lagrangian module
 *
 * \param[in]  isuite
 * \param[in]  iccvfg
 * \param[in]  iscalt
 * \param[in]  dtref
 */
/* ---------------------------------------------------------------------- */

void
CS_PROCF (lagopt, LAGOPT) (cs_int_t   *isuite,
                           cs_int_t   *iccvfg,
                           cs_int_t   *iscalt,
                           cs_real_t  *dtref)
{
  cs_lagr_option_definition(isuite, iccvfg, iscalt, dtref);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian module options definition.
 *
 * - default initialization
 * - read user settings
 * - check settings coherency
 * - initialize some structures relative to Lagrangian module
 *
 * \param[in]  isuite
 * \param[in]  iccvfg
 * \param[in]  iscalt
 * \param[in]  dtref
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_option_definition(cs_int_t   *isuite,
                          cs_int_t   *iccvfg,
                          cs_int_t   *iscalt,
                          cs_real_t  *dtref)
{
  /* Short-name and write access pointers to global variables */

  const cs_lagr_const_dim_t *const_dim = cs_glob_lagr_const_dim;

  cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  cs_lagr_time_scheme_t * lagr_time_scheme = cs_glob_lagr_time_scheme;
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_dim_t *lagdim = cs_glob_lagr_dim;

  /* Default initializations for Lagrangian module. */

  lagr_time_scheme->iilagr = 0;
  lagr_time_scheme->isuila = 0;

  cs_glob_lagr_stat_options->isuist = 0;

  lagr_model->physical_model = 0;

  cs_glob_lagr_specific_physics->idpvar = 0;

  cs_glob_lagr_specific_physics->itpvar = 0;

  cs_glob_lagr_specific_physics->impvar = 0;

  cs_glob_lagr_specific_physics->tpart = -999.0;

  cs_glob_lagr_specific_physics->cppart = -999.0;

  lagr_model->fouling = 0;

  /* Initializations for physical models */
  _init_lagr_encrustation_pointers();
  _init_lagr_boundary_interaction_pointers();

  lagr_time_scheme->isttio = 0;

  cs_glob_lagr_source_terms->nstits = 1;
  cs_glob_lagr_source_terms->ltsdyn = 0;
  cs_glob_lagr_source_terms->ltsmas = 0;
  cs_glob_lagr_source_terms->ltsthe = 0;

  cs_glob_lagr_stat_options->idstnt = 1;
  cs_glob_lagr_stat_options->nstist = 1;

  cs_glob_lagr_boundary_interactions->nombrd = NULL;

  lagr_time_scheme->t_order = 2;
  lagr_time_scheme->idistu = 1;
  lagr_time_scheme->idiffl = 0;
  lagr_time_scheme->modcpl = 0;
  lagr_time_scheme->idirla = 0;
  lagr_time_scheme->ilapoi = 0;
  lagr_time_scheme->iadded_mass = 0;
  lagr_time_scheme->added_mass_const = 1.0;

  cs_glob_lagr_boundary_interactions->inbrbd = 0;
  cs_glob_lagr_boundary_interactions->iflmbd = 0;
  cs_glob_lagr_boundary_interactions->iangbd = 0;
  cs_glob_lagr_boundary_interactions->ivitbd = 0;
  cs_glob_lagr_boundary_interactions->iencnbbd = 0;
  cs_glob_lagr_boundary_interactions->iencmabd = 0;
  cs_glob_lagr_boundary_interactions->iencdibd = 0;
  cs_glob_lagr_boundary_interactions->iencckbd = 0;
  cs_glob_lagr_boundary_interactions->nusbor = 0;

  for (int ii = 0; ii < cs_glob_lagr_const_dim->nusbrd + 10; ii++)
    cs_glob_lagr_boundary_interactions->imoybr[ii] = 0;

  /* User setup
     ---------- */

  if (cs_gui_file_is_loaded())
    cs_gui_particles_model();

  cs_user_lagr_model();

  if (lagr_time_scheme->iilagr == 0) {

    _free_lagr_encrustation_pointers();
    _free_lagr_boundary_interaction_pointers();

    BFT_FREE(cs_glob_lagr_source_terms->itsmv1);
    BFT_FREE(cs_glob_lagr_source_terms->itsmv2);

    cs_lagr_finalize_zone_conditions();

    return;
  }

  /* Check user initializations of Lagrangian module
     ----------------------------------------------- */

  int iok = 0;

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->iilagr",
                                lagr_time_scheme->iilagr,
                                0, 4);

  /* Restart needed if computation on frozen field.
     Note that for the Lagrangian module, frozen field also includes scalars. */

  if (lagr_time_scheme->iilagr == 3 && *isuite != 1)
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("in Lagrangian module"),
       _("The specified Lagrangian time scheme requires frozen fields\n"
         "(cs_glob_lagr_time_scheme->iilagr == %d)\n"
         "but the background Eulerian computation is not a restart.\n"),
       lagr_time_scheme->iilagr);

  if (lagr_time_scheme->iilagr == 3)
    *iccvfg = 1;

  if (   lagr_time_scheme->iilagr != 2
      && cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 1)
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("in Lagrangian module"),
       _("The pulverized coal coupled with Lagrangian particle transport\n"
         "is activated, but the return coupling of the dispersed phase\n"
         "on the continuous phase is not activated:\n"
         "  cs_glob_lagr_time_scheme->iilagr = %d\n"
         "The return coupling must be acivated for this model:\n"
         "  cs_glob_lagr_time_scheme->iilagr = 2\n"),
       lagr_time_scheme->iilagr);

  if (lagr_time_scheme->iilagr > 0
      && (   cs_glob_time_step->is_local
          || cs_glob_time_step->is_variable)) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L''INDICATEUR SUR LE MODULE LAGRANGIEN IILAGR ET\n"
               "@      LE CHOIX DU TYPE DE PAS DE TEMPS IDTVAR\n"
               "@      ONT DES VALEURS INCOMPATIBLES (LAGOPT).\n"
               "@\n"
               "@       IILAGR = %d\n"
               "@\n"
               "@  Le module lagrangien ne peut pas etre active avec un pas\n"
               "@   de temps variable en temps et en espace. Seuls les pas\n"
               "@   de temps uniforme et constant, et variable en temps et\n"
               "@   uniforme en espace sont possibles.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier les valeurs de IILAGR et IDTVAR.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->iilagr);
    iok++;

  }

  /* ISUILA ISUIST   */
  if (lagr_time_scheme->isuila < 0 ||
      lagr_time_scheme->isuila > 1) {
    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L'INDICATEUR DE SUITE DU MODULE LAGRANGIEN A UNE\n"
               "@       VALEUR NON PERMISE (LAGOPT).\n"
               "@\n"
               "@    ISUILA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
               "@       IL VAUT ICI ISUILA = %d\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IILAGR.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->isuila);
    iok++;

  }

  if (lagr_time_scheme->isuila == 1 &&
      *isuite == 0) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ALERTE A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@                                                  (LAGOPT).\n"
               "@\n"
               "@  Le module lagrangien est active en suite de calcul,\n"
               "@   alors que le calcul de la phase continue n'est pas\n"
               "@   une suite.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de ISUILA.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n");
    iok++;

  }

  if (lagr_time_scheme->isuila == 1) {

    if (cs_glob_lagr_stat_options->isuist < 0 ||
        cs_glob_lagr_stat_options->isuist > 1)  {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR DE SUITE DE CALCUL SUR LES STATISTIQUES\n"
                 "@       VOLUMIQUE ET AUX FRONTIERES, AINSI QUE SUR LES\n"
                 "@       TERMES SOURCES DE COUPLAGES RETOUR\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@    ISUIST DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI ISUIST = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de ISUIST.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
        cs_glob_lagr_stat_options->isuist);
      iok++;

    }
  }
  else
    cs_glob_lagr_stat_options->isuist = 0;

  /* IPHYLA     */
  if (lagr_model->physical_model < 0 ||
      lagr_model->physical_model > 2) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L'INDICATEUR DES MODELES PHYSIQUES LIES AUX PARTICULES\n"
               "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
               "@\n"
               "@    IPHYLA DEVRAIT ETRE UN ENTIER EGAL A 0 1 OU 2\n"
               "@       IL VAUT ICI IPHYLA = %d\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IPHYLA.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
      lagr_model->physical_model);
    iok++;

  }

  if (iok != 0)
    cs_exit(1);

  /* IDPVAR ITPVAR IMPVAR
   * Couplage-retour uniquement vers la phase continue  */
  if (lagr_model->physical_model == 1) {

    if (cs_glob_lagr_specific_physics->idpvar < 0 ||
        cs_glob_lagr_specific_physics->idpvar > 1) {
      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DU DIAMETRE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IDPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IDPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IDPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->idpvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->itpvar < 0 ||
        cs_glob_lagr_specific_physics->itpvar > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA TEMPERATURE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     ITPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI ITPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de ITPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->itpvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->impvar < 0 ||
        cs_glob_lagr_specific_physics->impvar > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA MASSE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IMPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IMPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IMPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->impvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->itpvar == 1 &&
        *iscalt ==  -1) { //TODO module param

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA TEMPERATURE DES\n"
                 "@       PARTICULES EST ACTIVE (ITPVAR = %d)\n"
                 "@       ALORS QU'AUCUN SCALAIRE THERMIQUE N'EST DISPONIBLE\n"
                 "@\n"
                 "@     ISCALT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL 1\n"
                 "@       IL VAUT ICI ISCALT = %d\n"
                 "@\n"
                 "@  La valeur de ISCALT est renseignee automatiquement\n"
                 "@    si une physique particuliere est activee dans USPPMO.\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
        cs_glob_lagr_specific_physics->itpvar,
        *iscalt);
      iok++;

    }

  }
  else {

    cs_glob_lagr_specific_physics->itpvar = 0;
    cs_glob_lagr_specific_physics->impvar = 0;
    cs_glob_lagr_specific_physics->idpvar = 0;

  }

  if (lagr_time_scheme->isuila == 1 &&
      lagr_model->physical_model == 1 &&
      cs_glob_lagr_specific_physics->itpvar == 1) {

    if (cs_glob_lagr_specific_physics->cppart < 0) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LA CHALEUR MASSIQUE D'INITIALISATION DES PARTICULES\n"
                 "@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     CPPART DEVRAIT ETRE UN REEL STRICTEMENT POSITIF\n"
                 "@       IL VAUT ICI CPPART = %14.5E\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de CPPART.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
        cs_glob_lagr_specific_physics->cppart);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->tpart <  -273.15) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LA TEMPERATURE D'INITIALISATION DES PARTICULES\n"
                 "@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     TPART DEVRAIT ETRE UN REEL SUPERIEUR A %14.5E\n"
                 "@       (EN DEGRES CELSIUS)\n"
                 "@       IL VAUT ICI TPART = %14.5E@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de TPART.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 -273.15,
                 cs_glob_lagr_specific_physics->tpart);
      iok++;

    }

  }

  if (iok != 0)
    cs_exit(1);

  /* IENCRA TPRENC VISREF */
  if (lagr_model->physical_model == 2) {

    if (lagr_time_scheme->t_order == 2) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
                 "@      EST ACTIVE (LAGOPT) AVEC UN SCHEMA D'INTEGRATION\n"
                 "@      AU SECOND ORDRE\n"
                 "@\n"
                 "@       IPHYLA = %d\n"
                 "@       NORDRE = %d\n"
                 "@\n"
                 "@  Le transport Lagrangien de particule de charbon ne peut\n"
                 "@   etre resolu au second ordre. Il faudrait mettre\n"
                 "@   à jour les équations\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de NORDRE\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->physical_model,
                 lagr_time_scheme->t_order);
      iok++;

    }

    if (cs_glob_lagr_source_terms->ltsthe == 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
                 "@      EST ACTIVE (LAGOPT) AVEC COUPLAGE RETOUR THERMIQUE\n"
                 "@\n"
                 "@       IPHYLA = %d\n"
                 "@       LTSTHE = %d\n"
                 "@\n"
                 "@  Le transport Lagrangien de particule de charbon ne peut\n"
                 "@   etre couple avec la phase Eulerienne. Il faudrait mettre\n"
                 "@   à jour les équations\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de LTSTHE\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->physical_model,
                 cs_glob_lagr_source_terms->ltsthe);
      iok++;

    }

    if (lagr_model->fouling < 0 ||
        lagr_model->fouling > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'ENCRASSEMENT DES PARTICULES\n"
                 "@       DE CHARBON A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IENCRA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IENCRA = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IENCRA.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->fouling);
      iok++;

    }

    for (int icha = 0; icha < extra->ncharb; icha++) {

      if (lagr_model->fouling == 1 &&
          cs_glob_lagr_encrustation->visref[icha] < 0) {
        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES\n"
                   "@       DE CHARBON EST ACTIVE (IENCRA = %d)\n"
                   "@       AVEC UNE VALEUR DE VISCOSITE CRITIQUE\n"
                   "@       NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@     VISREF DEVRAIT ETRE UN REEL STRICTEMENT POSITIF (Pa.s)\n"
                   "@       IL VAUT ICI VISREF = %14.5E\n"
                   "@       POUR LE CHARBON :%d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute.\n"
                   "@\n"
                   "@  Verifier la valeur de VISREF.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   lagr_model->fouling,
                   cs_glob_lagr_encrustation->visref[icha],
                   icha);
        iok++;

      }

      if (lagr_model->fouling == 1 &&
          cs_glob_lagr_encrustation->tprenc[icha] < 150.0) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR SUR L'ENCRASSEMENT DES PARTICULES\n"
                   "@       DE CHARBON EST ACTIVE (IENCRA = %d)\n"
                   "@       AVEC UNE VALEUR DE TEMPERATURE SEUIL\n"
                   "@       NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@     TPRENC DEVRAIT ETRE UN REEL SUPERIEUR A %14.5E\n"
                   "@       (EN DEGRES CELSIUS)\n"
                   "@       IL VAUT ICI TPRENC = %14.5E\n"
                   "@       POUR LE CHARBON :%d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute. Risque de division par\n"
                   "@  zero lors du calcul de la viscosite du charbon dans\n"
                   "@  cs_lagr_tracking\n"
                   "@\n"
                   "@  Verifier la valeur de TPRENC.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   lagr_model->fouling,
                   150.0e0,
                   cs_glob_lagr_encrustation->tprenc[icha],
                   icha);
        iok++;

      }

    }

  }
  else
    lagr_model->fouling = 0;

  if (lagr_model->physical_model != 2 &&
      cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE\n"
               "@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES\n"
               "@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE L'OPTION\n"
               "@      TRANSPORT DE PARTICULE DE CHARBON\n"
               "@      N'EST PAS ENCLENCHEE (LAGOPT).\n"
               "@\n"
               "@       IPHYLA = %d\n"
               "@       IPPMOD(ICPL3C) = %d\n"
               "@\n"
               "@  Le module lagrangien doit etre active en mode transport\n"
               "@   de particules de charbon pour etre couple avec la\n"
               "@   combustion d'une flamme de charbon pulverise en phase\n"
               "@   continue.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IPHYLA et\n"
               "@  verifier la valeur de IPPMOD dans la subroutine USPPMO.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model,
               cs_glob_physical_model_flag[CS_COMBUSTION_PCLC]);
    iok++;

  }
  if (lagr_model->physical_model == 2 &&
      (cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] < 0 &&
       cs_glob_physical_model_flag[CS_COMBUSTION_COAL] < 0)) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE\n"
               "@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES\n"
               "@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE LE COUPLAGE\n"
               "@      RETOUR DE LA PHASE DISPERSEE SUR LE PHASE CONTINUE\n"
               "@      N''EST PAS ENCLENCHE (LAGOPT).\n"
               "@\n"
               "@       IILAGR = %d\n"
               "@       IPPMOD(ICPL3C) = %d\n"
               "@\n"
               "@  Le module lagrangien doit etre active en mode couplage\n"
               "@   retour pour etre couple avec la combustion d'une\n"
               "@   flamme de charbon pulverise en phase continue.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IILAGR et\n"
               "@  verifier la valeur de IPPMOD dans la subroutine USPPMO.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->iilagr,
               cs_glob_physical_model_flag[CS_COMBUSTION_PCLC]);
    iok++;

  }

  /* ISUILA ISUIST   */
  if (lagr_time_scheme->isuila < 0 ||
      lagr_time_scheme->isuila > 1) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L'INDICATEUR DE SUITE DU MODULE LAGRANGIEN A UNE\n"
               "@       VALEUR NON PERMISE (LAGOPT).\n"
               "@\n"
               "@    ISUILA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
               "@       IL VAUT ICI ISUILA = %d\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IILAGR.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->isuila);
    iok++;

  }

  if (lagr_time_scheme->isuila == 1
      && *isuite == 0) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ALERTE A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@                                                  (LAGOPT).\n"
               "@\n"
               "@  Le module lagrangien est active en suite de calcul,\n"
               "@   alors que le calcul de la phase continue n'est pas\n"
               "@   une suite.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de ISUILA.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n");
    iok++;

  }

  if (lagr_time_scheme->isuila == 1) {

    if (   cs_glob_lagr_stat_options->isuist < 0
           || cs_glob_lagr_stat_options->isuist > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR DE SUITE DE CALCUL SUR LES STATISTIQUES\n"
                 "@       VOLUMIQUE ET AUX FRONTIERES, AINSI QUE SUR LES\n"
                 "@       TERMES SOURCES DE COUPLAGES RETOUR\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@    ISUIST DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI ISUIST = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de ISUIST.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_stat_options->isuist);
      iok++;

    }

  }
  else
    cs_glob_lagr_stat_options->isuist = 0;

  /* IPHYLA     */
  if (   lagr_model->physical_model < 0
         || lagr_model->physical_model > 2) {
    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L'INDICATEUR DES MODELES PHYSIQUES LIES AUX PARTICULES\n"
               "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
               "@\n"
               "@    IPHYLA DEVRAIT ETRE UN ENTIER EGAL A 0 1 OU 2\n"
               "@       IL VAUT ICI IPHYLA = %d\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IPHYLA.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model);
    iok++;

  }

  if (iok != 0)
    cs_exit(1);

  /* IDPVAR ITPVAR IMPVAR */
  /* Couplage-retour uniquement vers la phase continue  */
  if (lagr_model->physical_model == 1) {

    if (cs_glob_lagr_specific_physics->idpvar < 0 ||
        cs_glob_lagr_specific_physics->idpvar > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DU DIAMETRE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IDPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IDPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IDPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->idpvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->itpvar < 0 ||
        cs_glob_lagr_specific_physics->itpvar > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA TEMPERATURE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     ITPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI ITPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de ITPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->itpvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->impvar < 0 ||
        cs_glob_lagr_specific_physics->impvar > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA MASSE DES\n"
                 "@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IMPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IMPVAR = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IMPVAR.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "",
                 cs_glob_lagr_specific_physics->impvar);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->itpvar == 1 &&
        *iscalt == -1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'EQUATION DE LA TEMPERATURE DES\n"
                 "@       PARTICULES EST ACTIVE (ITPVAR = %d)\n"
                 "@       ALORS QU'AUCUN SCALAIRE THERMIQUE N'EST DISPONIBLE\n"
                 "@\n"
                 "@     ISCALT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL 1\n"
                 "@       IL VAUT ICI ISCALT = %d\n"
                 "@\n"
                 "@  La valeur de ISCALT est renseignee automatiquement\n"
                 "@    si une physique particuliere est activee dans USPPMO.\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->itpvar,
                 *iscalt);
      iok++;

    }

  }
  else{

    cs_glob_lagr_specific_physics->itpvar = 0;
    cs_glob_lagr_specific_physics->impvar = 0;
    cs_glob_lagr_specific_physics->idpvar = 0;

  }

  if (   lagr_time_scheme->isuila == 1
      && lagr_model->physical_model == 1
      && cs_glob_lagr_specific_physics->itpvar == 1) {

    if (cs_glob_lagr_specific_physics->cppart < 0) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LA CHALEUR MASSIQUE D'INITIALISATION DES PARTICULES\n"
                 "@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     CPPART DEVRAIT ETRE UN REEL STRICTEMENT POSITIF\n"
                 "@       IL VAUT ICI CPPART = %14.5E\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de CPPART.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_specific_physics->cppart);
      iok++;

    }

    if (cs_glob_lagr_specific_physics->tpart <  -273.15) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LA TEMPERATURE D'INITIALISATION DES PARTICULES\n"
                 "@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     TPART DEVRAIT ETRE UN REEL SUPERIEUR A %14.5E\n"
                 "@       (EN DEGRES CELSIUS)\n"
                 "@       IL VAUT ICI TPART = %14.5E@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de TPART.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 -273.15,
                 cs_glob_lagr_specific_physics->tpart);
      iok++;

    }

  }


  if (iok != 0)
    cs_exit(1);


  /* IENCRA TPRENC VISREF */
  if (lagr_model->physical_model == 2) {

    if (lagr_time_scheme->t_order == 2) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
                 "@      EST ACTIVE (LAGOPT) AVEC UN SCHEMA D'INTEGRATION\n"
                 "@      AU SECOND ORDRE\n"
                 "@\n"
                 "@       IPHYLA = %d\n"
                 "@       NORDRE = %d\n"
                 "@\n"
                 "@  Le transport Lagrangien de particule de charbon ne peut\n"
                 "@   etre resolu au second ordre. Il faudrait mettre\n"
                 "@   à jour les équations\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de NORDRE\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->physical_model,
                 lagr_time_scheme->t_order);
      iok++;

    }

    if (cs_glob_lagr_source_terms->ltsthe == 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
                 "@      EST ACTIVE (LAGOPT) AVEC COUPLAGE RETOUR THERMIQUE\n"
                 "@\n"
                 "@       IPHYLA = %d\n"
                 "@       LTSTHE = %d\n"
                 "@\n"
                 "@  Le transport Lagrangien de particule de charbon ne peut\n"
                 "@   etre couple avec la phase Eulerienne. Il faudrait mettre\n"
                 "@   à jour les équations\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de LTSTHE\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->physical_model,
                 cs_glob_lagr_source_terms->ltsthe);
      iok++;

    }

    if (lagr_model->fouling < 0 ||
        lagr_model->fouling > 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR L'ENCRASSEMENT DES PARTICULES\n"
                 "@       DE CHARBON A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@     IENCRA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                 "@       IL VAUT ICI IENCRA = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IENCRA.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_model->fouling);
      iok++;

    }

    for (int icha = 0; icha < extra->ncharb; icha++) {

      if (lagr_model->fouling == 1 &&
          cs_glob_lagr_encrustation->visref[icha] < 0) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR SUR L'ENCRASSEMENT DES PARTICULES\n"
                   "@       DE CHARBON EST ACTIVE (IENCRA = %d)\n"
                   "@       AVEC UNE VALEUR DE VISCOSITE CRITIQUE\n"
                   "@       NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@     VISREF DEVRAIT ETRE UN REEL STRICTEMENT POSITIF (Pa.s)\n"
                   "@       IL VAUT ICI VISREF = %14.5E\n"
                   "@       POUR LE CHARBON : %d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute.\n"
                   "@\n"
                   "@  Verifier la valeur de VISREF.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   lagr_model->fouling,
                   cs_glob_lagr_encrustation->visref[icha],
                   icha);
        iok++;

      }

      if (lagr_model->fouling == 1 &&
          cs_glob_lagr_encrustation->tprenc[icha] < 150.0) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR SUR L'ENCRASSEMENT DES PARTICULES\n"
                   "@       DE CHARBON EST ACTIVE (IENCRA = %d)\n"
                   "@       AVEC UNE VALEUR DE TEMPERATURE SEUIL\n"
                   "@       NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@     TPRENC DEVRAIT ETRE UN REEL SUPERIEUR A %14.5E\n"
                   "@       (EN DEGRES CELSIUS)\n"
                   "@       IL VAUT ICI TPRENC = %14.5E\n"
                   "@       POUR LE CHARBON :%d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute. Risque de division par\n"
                   "@  zero lors du calcul de la viscosite du charbon dans\n"
                   "@  cs_lagr_tracking\n"
                   "@\n"
                   "@  Verifier la valeur de TPRENC.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   lagr_model->fouling,
                   150.0e0,
                   cs_glob_lagr_encrustation->tprenc[icha],
                   icha);
        iok++;

      }

    }

  }
  else
    lagr_model->fouling = 0;


  if (lagr_model->physical_model != 2 &&
      cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE\n"
               "@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES\n"
               "@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE L'OPTION\n"
               "@      TRANSPORT DE PARTICULE DE CHARBON\n"
               "@      N'EST PAS ENCLENCHEE (LAGOPT).\n"
               "@\n"
               "@       IPHYLA = %d\n"
               "@       IPPMOD(ICPL3C) = %d\n"
               "@\n"
               "@  Le module lagrangien doit etre active en mode transport\n"
               "@   de particules de charbon pour etre couple avec la\n"
               "@   combustion d'une flamme de charbon pulverise en phase\n"
               "@   continue.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IPHYLA et\n"
               "@  verifier la valeur de IPPMOD dans la subroutine USPPMO.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model,
               cs_glob_physical_model_flag[CS_COMBUSTION_PCLC]);
    iok++;

  }

  if (   lagr_model->physical_model == 2
      && (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] < 0
          && cs_glob_physical_model_flag[CS_COMBUSTION_COAL] < 0)) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
               "@      EST ACTIVE (LAGOPT), ALORS QU'AUCUNE PHYSIQUE\n"
               "@      PARTICULIERE SUR LA COMBUSTION DU CHABON PULVERISE\n"
               "@      N'EST PAS ENCLENCHE (USPPMO).\n"
               "@\n"
               "@       IPHYLA = %d\n"
               "@       IPPMOD(ICPL3C) = %d\n"
               "@       IPPMOD(ICP3PL) = %d\n"
               "@\n"
               "@  Le transport lagrangien de particule de charbon doit\n"
               "@   etre couple avec la combustion d'une flamme de charbon\n"
               "@   pulverise en phase continue.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IPHYLA et\n"
               "@  verifier la valeur de IPPMOD dans la subroutine USPPMO.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model,
               cs_glob_physical_model_flag[CS_COMBUSTION_PCLC],
               cs_glob_physical_model_flag[CS_COMBUSTION_COAL]);
    iok++;

  }

  if (lagr_model->physical_model == 2 &&
      const_dim->nlayer < 1) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
               "@      EST ACTIVE (LAGOPT), ALORS QUE LA PARTICULE N'EST\n"
               "@      DISCRETISEE EN AUCUN COUCHE\n"
               "@\n"
               "@       IPHYLA = %d\n"
               "@       NLAYER = %d\n"
               "@\n"
               "@  Il doit y avoir au moins une couche par particule pour\n"
               "@    que le transport lagrangien de particule de charbon\n"
               "@    soit calculé.\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de NLAYER dans la subroutine LAGPAR\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model,
               const_dim->nlayer);
    iok++;

  }

  if (lagr_model->physical_model == 2 &&
      const_dim->nlayer > 99) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON\n"
               "@      EST ACTIVE (LAGOPT) AVEC CALCUL DES STATISTIQUES\n"
               "@    LA PARTICULE EST DISCRETISEE EN COUCHE:\n"
               "@\n"
               "@       IPHYLA = %d\n"
               "@       NLAYER = %d\n"
               "@\n"
               "@  Il y a trop de couche de discrétisation. nlayer devrait\n"
               "@    etre inferieur à 99. Sinon, il y a un problème au\n"
               "@    niveau du nom des variables (XXXXXX_layer_XX).\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de NLAYER dans la subroutine LAGPAR\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_model->physical_model,
               const_dim->nlayer);
    iok++;

  }

  if (iok != 0)
    cs_exit(1);

  if (iok != 0)
    cs_exit(1);

  /* ISTTIO NSTITS LTSDYN LTSMAS LTSTHE  */
  /* Si champs figes alors forcement en stationnaire    */
  if (lagr_time_scheme->iilagr == 3)
    lagr_time_scheme->isttio = 1;

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->isttio",
                                  lagr_time_scheme->isttio,
                                  0, 2);

  if (lagr_time_scheme->iilagr == 2) {

    if (lagr_time_scheme->isttio == 1 &&
        cs_glob_lagr_source_terms->nstits < 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR LE DECLENCHEMENT DU CALCUL\n"
                 "@       STATIONNAIRE DES STATISTIQUES POUR UN COUPLAGE RETOUR\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@    NSTITS DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1\n"
                 "@       IL VAUT ICI NSTITS = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de NSTITS.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_source_terms->nstits);
      iok++;

    }

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("in Lagrangian module"),
                                  "cs_glob_lagr_source_terms->ltsdyn",
                                  cs_glob_lagr_source_terms->ltsdyn,
                                  0, 2);

    if (     lagr_model->physical_model == 1
        && (  cs_glob_lagr_specific_physics->impvar == 1
            ||cs_glob_lagr_specific_physics->idpvar == 1)) {

      if (cs_glob_lagr_source_terms->ltsmas < 0 ||
          cs_glob_lagr_source_terms->ltsmas > 1) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR SUR LE COUPLAGE RETOUR SUR LA MASSE\n"
                   "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@    LTSMAS DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                   "@       IL VAUT ICI LTSMAS = %d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute.\n"
                   "@\n"
                   "@  Verifier la valeur de LTSMAS.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   cs_glob_lagr_source_terms->ltsmas);
        iok++;

      }

    }
    else
      cs_glob_lagr_source_terms->ltsmas = 0;

    if ((lagr_model->physical_model == 1 && cs_glob_lagr_specific_physics->itpvar == 1)
        || lagr_model->physical_model == 2) {

      if (cs_glob_lagr_source_terms->ltsthe < 0 ||
          cs_glob_lagr_source_terms->ltsthe > 1) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR SUR LE COUPLAGE RETOUR SUR LA THERMIQUE\n"
                   "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@    LTSTHE DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1\n"
                   "@       IL VAUT ICI LTSTHE = %d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute.\n"
                   "@\n"
                   "@  Verifier la valeur de LTSTHE.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   cs_glob_lagr_source_terms->ltsthe);
        iok++;

      }

    }
    else
      cs_glob_lagr_source_terms->ltsthe = 0;

    if (cs_glob_lagr_source_terms->ltsdyn == 1 &&
        *iccvfg == 1) {
      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR LE COUPLAGE RETOUR SUR LA DYNAMIQUE\n"
                 "@       EST ACTIVE (LTSDYN = %d) (LAGOPT)\n"
                 "@       ALORS QUE LA PHASE PORTEUSE EST CALCULEE AVEC\n"
                 "@       L'OPTION CHAMP FIGE  (ICCVFG = %d).\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de LTSDYN et\n"
                 "@  verifier la valeur de ICCVFG.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_source_terms->ltsdyn,
                 *iccvfg);
      iok++;
    }

    if (cs_glob_lagr_source_terms->ltsdyn != 1 &&
        cs_glob_lagr_source_terms->ltsthe != 1 &&
        cs_glob_lagr_source_terms->ltsmas != 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========   (LAGOPT).\n"
                 "@\n"
                 "@    L'INDICATEUR SUR LE COUPLAGE RETOUR EST ACTIVE\n"
                 "@        IILAGR = %d\n"
                 "@      ALORS QU'AUCUN COUPLAGE RETOUR N'EST ENCLENCHE\n"
                 "@        DYNAMIQUE : LTSDYN = %d\n"
                 "@        THERMIQUE : LTSTHE = %d\n"
                 "@        MASSIQUE  : LTSMAS = %d\n"
                 "@\n"
                 "@    LES COUPLAGES RETOUR SUR LA THERMIQUE ET SUR LA MASSE\n"
                 "@      NECESSITENT L'ACTIVATION D'UNE PHYSIQUE ADEQUATE\n"
                 "@      ASSOCIEE AUX PARTICULES.\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IILAGR.\n"
                 "@  Verifier la valeur de IPHYLA.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_time_scheme->iilagr,
                 cs_glob_lagr_source_terms->ltsdyn,
                 cs_glob_lagr_source_terms->ltsthe,
                 cs_glob_lagr_source_terms->ltsmas);
      iok++;

    }

  }
  else{

    cs_glob_lagr_source_terms->ltsdyn = 0;
    cs_glob_lagr_source_terms->ltsmas = 0;
    cs_glob_lagr_source_terms->ltsthe = 0;
  }

  if (iok != 0)
    cs_exit(1);

  {
    if (cs_glob_lagr_stat_options->idstnt < 1) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@    IDSTNT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1\n"
                 "@       IL VAUT ICI IDSTNT = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IDSTNT.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 cs_glob_lagr_stat_options->idstnt);
      iok++;

    }

    if (lagr_time_scheme->isttio == 1) {

      if (cs_glob_lagr_stat_options->nstist < cs_glob_lagr_stat_options->idstnt) {

        bft_printf("@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n"
                   "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                   "@    =========\n"
                   "@    L'INDICATEUR DE CALCUL STATIONNAIRE DES STATISTIQUES\n"
                   "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                   "@\n"
                   "@    NSTIST DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL\n"
                   "@       A IDSTNT = %d\n"
                   "@       IL VAUT ICI NSTIST = %d\n"
                   "@\n"
                   "@  Le calcul ne sera pas execute.\n"
                   "@\n"
                   "@  Verifier la valeur de NSTIST.\n"
                   "@\n"
                   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                   "@\n",
                   cs_glob_lagr_stat_options->idstnt, cs_glob_lagr_stat_options->nstist);
        iok++;

      }

    }

  }


  if (iok != 0)
    cs_exit(1);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->t_order",
                                lagr_time_scheme->t_order,
                                1, 3);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->idistu",
                                lagr_time_scheme->idistu,
                                0, 2);

  if (   lagr_time_scheme->idistu == 1
         && extra->itytur != 2
         && extra->itytur != 3
         && extra->iturb != 50
         && extra->iturb != 60) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE\n"
               "@    DE TURBULENCE SELECTIONNE (LAGOPT).\n"
               "@\n"
               "@   Le module lagrangien a ete active avec IILAGR = %d\n"
               "@     et la dispersion turbulente est prise en compte\n"
               "@                                     avec IDISTU = %d\n"
               "@   Le modele de turbulence\n"
               "@     correspond a ITURB = %d\n"
               "@   Or, les seuls traitements de la turbulence compatibles\n"
               "@     avec le module Lagrangien et la dispersion turbulente\n"
               "@     sont k-epsilon et Rij-epsilon, v2f et k-omega\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de IILAGR, IDISTU, et ITURB.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->iilagr,
               lagr_time_scheme->idistu,
               extra->iturb);
    iok++;

  }
  else if (lagr_time_scheme->idistu == 0 &&
           extra->iturb != 0  &&
           extra->itytur!= 2  &&
           extra->itytur!= 3  &&
           extra->iturb != 50 &&
           extra->iturb != 60) {
    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE\n"
               "@    DE TURBULENCE SELECTIONNE (LAGOPT).\n"
               "@\n"
               "@   Le module lagrangien a ete active avec IILAGR = %d\n"
               "@     et la dispersion turbulente est prise en compte\n"
               "@                                     avec IDISTU = %d\n"
               "@   Le modele de turbulence\n"
               "@     correspond a ITURB = %d\n"
               "@   Or, les seuls traitements de la turbulence compatibles\n"
               "@     avec le module Lagrangien et la dispersion turbulente\n"
               "@     sont k-epsilon et Rij-epsilon, v2f et k-omega\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->iilagr,
               lagr_time_scheme->idistu,
               extra->iturb);
    iok++;

  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->idiffl",
                                lagr_time_scheme->idiffl,
                                0, 2);

  /* MODCPL IDIRLA   */
  if (lagr_time_scheme->modcpl < 0) {

    bft_printf("@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
               "@    =========\n"
               "@    L'INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION\n"
               "@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).\n"
               "@\n"
               "@    MODCPL DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0\n"
               "@       IL VAUT ICI MODCPL = %d\n"
               "@\n"
               "@  Le calcul ne sera pas execute.\n"
               "@\n"
               "@  Verifier la valeur de MODCPL.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               lagr_time_scheme->modcpl);
    iok++;

  }

  if (lagr_time_scheme->modcpl > 0) {

    if (lagr_time_scheme->modcpl < cs_glob_lagr_stat_options->idstnt) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    L'INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION\n"
                 "@       TURBULENTE EST INCOMPATIBLE AVEC CELUI DU CALCUL\n"
                 "@       DES STATISTIQUES (LAGOPT).\n"
                 "@\n"
                 "@    LE MODELE COMPLET DE DISPERSION TURBULENTE EST ACTIVE\n"
                 "@      (MODCPL = %d)\n"
                 "@      AVANT LE DEBUT DU CALCUL DES STATISTIQUES\n"
                 "@      (IDSTNT = %d)\n"
                 "@\n"
                 "@  Il est necessaire d'avoir calcule des statistiques\n"
                 "@  pour declencher le modele de dispersion turbulent complet.\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de MODCPL.\n"
                 "@  Verifier la valeur de IDSTNT.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_time_scheme->modcpl,
                 cs_glob_lagr_stat_options->idstnt);
      iok++;

    }

    /* Velocity statistics are needed for this model */
    cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

    if (lagr_time_scheme->idirla != 1 &&
        lagr_time_scheme->idirla != 2 &&
        lagr_time_scheme->idirla != 3) {

      bft_printf("@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n"
                 "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                 "@    =========\n"
                 "@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET\n"
                 "@       A UNE VALEUR NON PERMISE (LAGOPT).\n"
                 "@\n"
                 "@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1, 2 OU 3\n"
                 "@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L'AXE X,\n"
                 "@        LA VALEUR 2 POUR UN ECOULEMENT SELON L'AXE Y,\n"
                 "@        LA VALEUR 3 POUR UN ECOULEMENT SELON L'AXE Z)\n"
                 "@       IL VAUT ICI IDIRLA = %d\n"
                 "@\n"
                 "@  Le calcul ne sera pas execute.\n"
                 "@\n"
                 "@  Verifier la valeur de IDIRLA.\n"
                 "@\n"
                 "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                 "@\n",
                 lagr_time_scheme->idirla);
      iok++;

    }

  }

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_time_scheme->ilapoi",
                                lagr_time_scheme->ilapoi,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_boundary_interactions->inbrbd",
                                cs_glob_lagr_boundary_interactions->inbrbd,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_boundary_interactions->iflmbd",
                                cs_glob_lagr_boundary_interactions->iflmbd,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_boundary_interactions->iangbd",
                                cs_glob_lagr_boundary_interactions->iangbd,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_boundary_interactions->ivitbd",
                                cs_glob_lagr_boundary_interactions->ivitbd,
                                0, 2);

  cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                _("in Lagrangian module"),
                                "cs_glob_lagr_boundary_interactions->nusbor",
                                cs_glob_lagr_boundary_interactions->nusbor,
                                0, cs_glob_lagr_const_dim->nusbrd + 1);

  if (lagr_model->physical_model == 2 &&
      lagr_model->fouling == 1) {

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("in Lagrangian module"),
                                  "cs_glob_lagr_boundary_interactions->iencnbbd",
                                  cs_glob_lagr_boundary_interactions->iencnbbd,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("in Lagrangian module"),
                                  "cs_glob_lagr_boundary_interactions->iencmabd",
                                  cs_glob_lagr_boundary_interactions->iencmabd,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("in Lagrangian module"),
                                  "cs_glob_lagr_boundary_interactions->iencdibd",
                                  cs_glob_lagr_boundary_interactions->iencdibd,
                                  0, 2);

    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("in Lagrangian module"),
                                  "cs_glob_lagr_boundary_interactions->iencckbd",
                                  cs_glob_lagr_boundary_interactions->iencckbd,
                                  0, 2);


  }
  else{

    cs_glob_lagr_boundary_interactions->iencnbbd = 0;
    cs_glob_lagr_boundary_interactions->iencmabd = 0;
    cs_glob_lagr_boundary_interactions->iencdibd = 0;
    cs_glob_lagr_boundary_interactions->iencckbd = 0;

  }

  if (iok != 0)
    cs_exit (1);

  cs_parameters_error_barrier();

  /* ============================================================================== */
  /* 3. INITIALISATIONS DES VARIABLES EN COMMON    */
  /* ATTENTION :     */
  /* ^^^^^^^^^^^     */
  /* CES INITIALISATIONS NE DOIVENT PAS ETRE MODIFIEES PAR L'UTILISATEUR   */
  /* ============================================================================== */

  /* 3.1 GENERALITES (D'AUTRES INITIALISATIONS SONT FAITES DANS LAGLEC) */

  /* PAS DE TEMPS LAGRANGIEN (LAGUNE) : Par defaut le pas de temps     */
  /* de reference de la phase continue   */
  cs_glob_lagr_time_step->dtp = *dtref;

  /* TEMPS COURANT PHYSIQUE LAGRANGIEN   */
  cs_glob_lagr_time_step->ttclag = 0.0;

  /* STATISTIQUES AUX FRONTIERES    */
  /* Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES    */
  /* des stats aux frontieres  */
  cs_glob_lagr_boundary_interactions->npstf = 0;

  /* Nombre de pas de temps total des stats aux frontieres   */
  /* depuis le debut du calcul, partie instationnaire comprise    */
  cs_glob_lagr_boundary_interactions->npstft = 0;

  /* Temps physique des stats aux frontieres  */
  cs_glob_lagr_boundary_interactions->tstatp = 0.0;

  /* COUPLAGE RETOUR */
  /* Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES    */
  /* des termes sources pour le couplage retour    */
  cs_glob_lagr_source_terms->npts = 0;

  /* Initialisation du sous-pas     */
  cs_glob_lagr_time_step->nor = 0;

  /* 3.7 DEFINITION DES POINTEURS LIES AUX STATISTIQUES AUX FRONTIERES
   * INBRBD : NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES
   * IFLMBD : FLUX DE MASSE PARTICULAIRE
   * IANGBD : ANGLE VITESSE
   * IVITBD : VITESSE DE LA PARTICULE

   * IENCNBBD : NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES AVEC ENCRASSEMENT CHARBON
   * IENCMABD : MASSE DE GRAINS DE CHARBON ENCRASSES
   * IENCDIBD : DIAMETRE DES GRAINS DE CHARBON ENCRASSES
   * IENCCKBD : FRACTION DE COKE DES GRAINS DE CHARBON ENCRASSES
   * NUSBOR : INFORMATIONS UTILISATEUR SUPPLEMENTAIRES
   * NVISBR : NOMBRE TOTAL D'INTERACTIONS A ENREGISTRER */

  int irf = -1;

  if (cs_glob_lagr_boundary_interactions->inbrbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->inbr = irf;
    _copy_boundary_varname(irf, "Part_impact_number");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;

  }

  if (cs_glob_lagr_boundary_interactions->iflmbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iflm = irf;
    _copy_boundary_varname(irf, "Part_bndy_mass_flux");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 1;

  }

  if (cs_glob_lagr_boundary_interactions->iangbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iang = irf;
    _copy_boundary_varname(irf, "Part_impact_angle");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 2;

  }

  if (cs_glob_lagr_boundary_interactions->ivitbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->ivit = irf;
    _copy_boundary_varname(irf, "Part_impact_velocity");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 2;

  }

  if (lagr_model->resuspension == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->ires = irf;
    _copy_boundary_varname(irf, "Part_resusp_number");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->iflres = irf;
    _copy_boundary_varname(irf, "Part_resusp_mass_flux");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 1;

  }

  if (lagr_model->clogging == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->inclg = irf;
    _copy_boundary_varname(irf, "Part_deposited_number");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->inclgt = irf;
    _copy_boundary_varname(irf, "Part_deposited_part");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->iclogt = irf;
    _copy_boundary_varname(irf, "Part_deposited_time");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->iclogh = irf;
    _copy_boundary_varname(irf, "Part_consolidation_height");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->iscovc = irf;
    _copy_boundary_varname(irf, "Part_surf_coverage");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->ihdepm = irf;
    _copy_boundary_varname(irf, "Part_dep_height_mean");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->ihdiam = irf;
    _copy_boundary_varname(irf, "Part_dep_diameter_mean");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->ihsum = irf;
    _copy_boundary_varname(irf, "Part_dep_diameter_sum");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;
    irf++;
    cs_glob_lagr_boundary_interactions->ihdepv = irf;
    _copy_boundary_varname(irf, "Part_dep_height_variance");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;

  }

  if (lagr_model->physical_model == 2 &&
      lagr_model->fouling == 1 &&
      cs_glob_lagr_boundary_interactions->iencnbbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iencnb = irf;
    _copy_boundary_varname(irf, "Part_fouled_impact_number");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;

  }

  if (lagr_model->physical_model == 2 &&
      lagr_model->fouling == 1 &&
      cs_glob_lagr_boundary_interactions->iencmabd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iencma = irf;
    _copy_boundary_varname(irf, "Part_fouled_mass_flux");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 1;

  }

  if (lagr_model->physical_model == 2 &&
      lagr_model->fouling == 1 &&
      cs_glob_lagr_boundary_interactions->iencdibd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iencdi = irf;
    _copy_boundary_varname(irf, "Part_fouled_diam");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 3;

  }

  if (lagr_model->physical_model == 2 &&
      lagr_model->fouling == 1 &&
      cs_glob_lagr_boundary_interactions->iencckbd == 1) {

    irf++;
    cs_glob_lagr_boundary_interactions->iencck = irf;
    _copy_boundary_varname(irf, "Part_fouled_Xck");
    cs_glob_lagr_boundary_interactions->imoybr[irf] = 3;

  }

  if (cs_glob_lagr_boundary_interactions->nusbor > 0) {

    for (int ii = 0; ii<cs_glob_lagr_boundary_interactions->nusbor; ii++) {

      irf++;
      char buf[64];
      cs_glob_lagr_boundary_interactions->iusb[ii] = irf;
      snprintf(buf, 64, "addRec%d", ii);
      _copy_boundary_varname(irf, buf);
      cs_glob_lagr_boundary_interactions->imoybr[irf] = 0;

    }

  }

  lagdim->nvisbr = irf + 1;

  /* 3.8 DEFINITION DES POINTEURS LIES AUX TERMES SOURCES LAGRANGIEN
   * POUR COUPLAGE RETOUR
   * Nombre de termes sources de couplage-retour   */

  irf = 0;

  lagdim->ntersl = 0;
  cs_glob_lagr_source_terms->itsli  = 0;
  cs_glob_lagr_source_terms->itske  = 0;
  cs_glob_lagr_source_terms->itsmas = 0;
  cs_glob_lagr_source_terms->itste  = 0;
  cs_glob_lagr_source_terms->itsti  = 0;

  /* if we don't use fortran initialization (NEPTUNE_CFD)
   * we need allocate memory */
  if (cs_glob_lagr_source_terms->itsmv1 == NULL)
    BFT_MALLOC(cs_glob_lagr_source_terms->itsmv1,
               cs_glob_lagr_const_dim->ncharm2, int);

  if (cs_glob_lagr_source_terms->itsmv2 == NULL)
    BFT_MALLOC(cs_glob_lagr_source_terms->itsmv2,
               cs_glob_lagr_const_dim->ncharm2, int);

  for (int icha = 0; icha < const_dim->ncharm2; icha++) {

    cs_glob_lagr_source_terms->itsmv1[icha] = 0;
    cs_glob_lagr_source_terms->itsmv2[icha] = 0;

  }

  cs_glob_lagr_source_terms->itsco = 0;
  cs_glob_lagr_source_terms->itsfp4 = 0;

  /* Dynamique : Vitesse + Turbulence    */
  if (cs_glob_lagr_source_terms->ltsdyn == 1) {

    _define_st_field("velocity_st_lagr", 3);

    lagdim->ntersl += 1;
    cs_glob_lagr_source_terms->itsli = ++irf;

    if (   extra->itytur == 2
        || extra->iturb == 50
        || extra->iturb == 60) {

      /* K-eps, v2f et k-omega     */
      lagdim->ntersl += 1;
      cs_glob_lagr_source_terms->itske = ++irf;

    }
    else if (extra->itytur == 3) {

      /* RIJ   */

      _define_st_field("rij_st_lagr", 6);

    }
    else {

      bft_error(__FILE__, __LINE__, 0,
                "@\n"
                "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                "@\n"
                "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
                "@    =========\n"
                "@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE\n"
                "@    DE TURBULENCE SELECTIONNE (LAGOPT).\n"
                "@\n"
                "@   Le module lagrangien a ete active avec IILAGR = %d\n"
                "@     et le couplage inverse sur la dynamique est pris en\n"
                "@                              compte avec LTSDYN = %d\n"
                "@   Le modele de turbulence\n"
                "@     correspond a ITURB = %d\n"
                "@   Or, les seuls traitements de la turbulence compatibles\n"
                "@     avec le module Lagrangien et le couplage inverse sur\n"
                "@     la dynamique sont k-epsilon, Rij-epsilon, v2f\n"
                "@     et k-omega\n"
                "@\n"
                "@  Le calcul ne sera pas execute.\n"
                "@\n"
                "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                "@\n",
                lagr_time_scheme->iilagr,
                cs_glob_lagr_source_terms->ltsdyn,
                extra->iturb);

    }

  }

  /* Modele de depot */
  if (   lagr_model->deposition == 1
      && lagr_time_scheme->t_order == 2) {

    bft_error(__FILE__, __LINE__, 0,
              "@\n"
              "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              "@\n"
              "@ @@ ATTENTION : ARRET A L'EXECUTION DU MODULE LAGRANGIEN\n"
              "@    =========\n"
              "@    LE MODELE SPECIFIQUE DE DEPOT (Guingo & Minier, 2008)\n"
              "@     EST ACTIVABLE UNIQUEMENT AVEC  UN SCHEMA D'ORDRE  1\n"
              "@\n"
              "@\n"
              "@  Le calcul ne sera pas execute.\n"
              "@\n"
              "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              "@\n");

  }

  /* Masse */
  if (cs_glob_lagr_source_terms->ltsmas == 1) {

    lagdim->ntersl += 1;
    cs_glob_lagr_source_terms->itsmas = irf + 1;
    irf = cs_glob_lagr_source_terms->itsmas;

  }

  /* Thermique  */
  if (cs_glob_lagr_source_terms->ltsthe == 1) {

    if (lagr_model->physical_model == 1) {

      /* Temperature     */
      if (cs_glob_lagr_specific_physics->itpvar == 1) {

        lagdim->ntersl += 2;
        cs_glob_lagr_source_terms->itste = irf + 1;
        cs_glob_lagr_source_terms->itsti = cs_glob_lagr_source_terms->itste + 1;
        irf = cs_glob_lagr_source_terms->itsti;

      }

    }

    /* Charbon    */
    else if (lagr_model->physical_model == 2) {

      lagdim->ntersl += 4 + 2 * extra->ncharb;
      cs_glob_lagr_source_terms->itste = irf + 1;
      cs_glob_lagr_source_terms->itsti
        = cs_glob_lagr_source_terms->itste + 1;

      for (int icha = 0; icha < extra->ncharb; icha++)
        cs_glob_lagr_source_terms->itsmv1[icha]
          = cs_glob_lagr_source_terms->itsti + icha;

      for (int icha = 0; icha < extra->ncharb; icha++)
        cs_glob_lagr_source_terms->itsmv2[icha]
          = cs_glob_lagr_source_terms->itsmv1[extra->ncharb] + icha;

      cs_glob_lagr_source_terms->itsco
        = cs_glob_lagr_source_terms->itsmv2[extra->ncharb] + 1;
      cs_glob_lagr_source_terms->itsfp4
        = cs_glob_lagr_source_terms->itsco + 1;

      irf = cs_glob_lagr_source_terms->itsfp4;

    }

  }

  /* Now define particle map */
  cs_lagr_particle_attr_initialize();

  if (lagr_model->deposition > 0)
    cs_field_find_or_create("ustar",
                            CS_FIELD_PROPERTY | CS_FIELD_PROPERTY,
                            CS_MESH_LOCATION_BOUNDARY_FACES,
                            1);

  /* Now activate basic statistics */
  cs_lagr_stat_initialize();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
