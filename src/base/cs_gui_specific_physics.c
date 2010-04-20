/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include "fvm_selector.h"

/*----------------------------------------------------------------------------
 * MEI library headers
 *----------------------------------------------------------------------------*/

#include "mei_evaluate.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_mesh.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_specific_physics.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the activated specific physics scalar number
 *----------------------------------------------------------------------------*/

static int
_scalar_number(const char* model)
{
  char *path = NULL;
  int   nb;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");

  nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return nb;
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Predefined physics indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPPMO
 * *****************
 *
 * INTEGER          IPPMOD <--  specific physics indicator array
 * INTEGER          ICOD3P  --> diffusion flame in fast complete chemistry
 * INTEGER          ICODEQ  --> diffusion flame in fast chemistry towards balance
 * INTEGER          ICOEBU  --> Eddy Break Up premixing flame
 * INTEGER          ICOBML  --> Bray - Moss - Libby premixing flame
 * INTEGER          ICOLWC  --> Libby Williams premixing flame
 * INTEGER          ICP3PL  --> Coal combustion. Combustible moyen local
 * INTEGER          ICPL3C  --> Coal combustion coupled with lagrangien approach
 * INTEGER          ICFUEL  --> Fuel combustion
 * INTEGER          IELJOU  --> Joule effect
 * INTEGER          IELARC  --> electrical arc
 * INTEGER          IELION  --> ionique mobility
 * INTEGER          ICOMPF  --> compressible without shock
 * INTEGER          IATMOS  --> atmospheric flows
 * INTEGER          IAEROS  --> cooling tower
 * INTEGER          INDJON  --> INDJON=1: a JANAF enthalpy-temperature
 *                              tabulation is used. INDJON=1: users tabulation
 * INTEGER          IEQCO2  --> CO2 massic fraction transport
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uippmo, UIPPMO)(int *const ippmod,
                               int *const icod3p,
                               int *const icodeq,
                               int *const icoebu,
                               int *const icobml,
                               int *const icolwc,
                               int *const icp3pl,
                               int *const icpl3c,
                               int *const icfuel,
                               int *const ieljou,
                               int *const ielarc,
                               int *const ielion,
                               int *const icompf,
                               int *const iatmos,
                               int *const iaeros,
                               int *const indjon,
                               int *const ieqco2)
{
    int isactiv = 0;
    int nscapp = 0;

    cs_var_t  *vars = cs_glob_var;

    ippmod[*icod3p - 1] = -1;
    ippmod[*icodeq - 1] = -1;
    ippmod[*icoebu - 1] = -1;
    ippmod[*icobml - 1] = -1;
    ippmod[*icolwc - 1] = -1;
    ippmod[*icp3pl - 1] = -1;
    ippmod[*icpl3c - 1] = -1;
    ippmod[*icfuel - 1] = -1;
    ippmod[*ieljou - 1] = -1;
    ippmod[*ielarc - 1] = -1;
    ippmod[*ielion - 1] = -1;
    ippmod[*icompf - 1] = -1;
    ippmod[*iatmos - 1] = -1;
    ippmod[*iaeros - 1] = -1;

    *indjon = 1;
    *ieqco2 = 0;

    /* Look for the active specific physics and give the value of the associated
       model attribute */
    isactiv = cs_gui_get_activ_thermophysical_model();

    if (isactiv)
    {
        if (cs_gui_strcmp(vars->model, "pulverized_coal"))
        {
            if (cs_gui_strcmp(vars->model_value, "coal_homo"))
                ippmod[*icp3pl - 1] = 0;
            else if (cs_gui_strcmp(vars->model_value, "coal_homo2"))
                ippmod[*icp3pl - 1] = 1;
            else
                bft_error(__FILE__, __LINE__, 0,
                          _("Invalid coal model: %s.\n"), vars->model_value);
        }
        else if  (cs_gui_strcmp(vars->model, "atmospheric_flows"))
        {
            if (cs_gui_strcmp(vars->model_value, "constant"))
                ippmod[*iatmos - 1] = 0;
            else if (cs_gui_strcmp(vars->model_value, "dry"))
                ippmod[*iatmos - 1] = 1;
            else if (cs_gui_strcmp(vars->model_value, "humid"))
                ippmod[*iatmos - 1] = 2;
            else
                bft_error(__FILE__, __LINE__, 0,
                          _("Invalid atmospheric flow model: %s.\n"),
                          vars->model_value);
        }

        /* If the model is active, one only takes the specific physics scalars */
        nscapp = _scalar_number(vars->model);
    }

    vars->nscapp = nscapp;

#if _XML_DEBUG_
    bft_printf("==>UIPPMO\n");
    if (isactiv)
    {
        bft_printf("--thermophysical model: %s\n", vars->model);
        bft_printf("--thermophysical value: %s\n", vars->model_value);
        bft_printf("--model scalars number: %i\n", vars->nscapp);
    }
#endif

}

/*----------------------------------------------------------------------------
 * Density under relaxation
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI1 (SRROM)
 * *****************
 * DOUBLE PRECISION SRROM   <--   density relaxation
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi1, UICPI1) (double *const srrom)
{
  cs_gui_numerical_double_parameters("density_relaxation", srrom);

#if _XML_DEBUG_
  bft_printf("==>UICPI1\n");
  bft_printf("--srrom  = %f\n", *srrom);
#endif
}

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicppr, UICPPR) (const int *const nclass,
                                const int *const nsalpp,
                                const int *const nsalto,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ihtco2,
                                const int *const itemp1,
                                const int *const irom1,
                                const int *const iym1,
                                const int *const immel,
                                const int *const itemp2,
                                const int *const ix2,
                                const int *const irom2,
                                const int *const idiam2,
                                const int *const igmdch,
                                const int *const igmdv1,
                                const int *const igmdv2,
                                const int *const igmhet,
                                const int *const ighco2,
                                const int *const igmsec)
{
  int i = 0;
  int n;
  char *name = NULL;
  char *snumpp = NULL;

  cs_var_t  *vars = cs_glob_var;

  n = vars->nprop;
  vars->nprop  += *nsalpp;
  vars->nsalpp  = *nsalpp;

  BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
  BFT_REALLOC(vars->propce,          vars->nprop, int);
  BFT_REALLOC(vars->properties_name, vars->nprop, char*);

 /* ITEMP1 */
  vars->properties_ipp[n] = ipppro[ ipproc[ *itemp1 -1 ]-1 ];
  vars->propce[n] = *itemp1;
  BFT_MALLOC(vars->properties_name[n], strlen("Temp_GAZ")+1, char);
  strcpy(vars->properties_name[n++], "Temp_GAZ");

 /* IROM1 */
  vars->properties_ipp[n] = ipppro[ ipproc[ *irom1 -1 ]-1 ];
  vars->propce[n] = *irom1;
  BFT_MALLOC(vars->properties_name[n], strlen("ROM_GAZ")+1, char);
  strcpy(vars->properties_name[n++], "ROM_GAZ");

 /*  YM_CHX1M */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[0] -1 ]-1 ];
  vars->propce[n] = iym1[0];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx1m")+1, char);
  strcpy(vars->properties_name[n++], "YM_CHx1m");

 /*  YM_CHX2M */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[1] -1 ]-1 ];
  vars->propce[n] = iym1[1];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx2m")+1, char);
  strcpy(vars->properties_name[n++], "YM_CHx2m");

 /*  YM_CO */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[2] -1 ]-1 ];
  vars->propce[n] = iym1[2];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO")+1, char);
  strcpy(vars->properties_name[n++], "YM_CO");

 /*  YM_O2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[3] -1 ]-1 ];
  vars->propce[n] = iym1[3];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_O2")+1, char);
  strcpy(vars->properties_name[n++], "YM_O2");

 /*  YM_CO2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[4] -1 ]-1 ];
  vars->propce[n] = iym1[4];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO2")+1, char);
  strcpy(vars->properties_name[n++], "YM_CO2");

 /*  YM_H2O */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[5] -1 ]-1 ];
  vars->propce[n] = iym1[5];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_H2O")+1, char);
  strcpy(vars->properties_name[n++], "YM_H2O");

 /*  YM_N2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[6] -1 ]-1 ];
  vars->propce[n] = iym1[6];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_N2")+1, char);
  strcpy(vars->properties_name[n++], "YM_N2");

 /* IMEL */
  vars->properties_ipp[n] = ipppro[ ipproc[ *immel -1 ]-1 ];
  vars->propce[n] = *immel;
  BFT_MALLOC(vars->properties_name[n], strlen("XM")+1, char);
  strcpy(vars->properties_name[n++], "XM");

  /* ITEMP2 loop on classes */
  BFT_MALLOC(name, strlen("Temp_CP")+1 + 2, char);
  BFT_MALLOC(snumpp, 1 + 2, char);
  strcpy(name, "Temp_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ itemp2[i] -1 ]-1 ];
    vars->propce[n] = itemp2[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Temp_CP");
  }

 /* IX2 loop on classes */
  BFT_REALLOC(name, strlen("Frm_CP")+1 + 2, char);
  strcpy(name, "Frm_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ ix2[i] -1 ]-1 ];
    vars->propce[n] = ix2[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Frm_CP");
  }

 /* IROM2 loop on classes */
  BFT_REALLOC(name, strlen("Rho_CP")+1 + 2, char);
  strcpy(name, "Rho_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ irom2[i] -1 ]-1 ];
    vars->propce[n] = irom2[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Rho_CP");
  }

 /* IDIAM2 loop on classes */
  BFT_REALLOC(name, strlen("Dia_CK")+1 + 2, char);
  strcpy(name, "Dia_CK");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ idiam2[i] -1 ]-1 ];
    vars->propce[n] = idiam2[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Dia_CK");
  }

 /* IGMDCH loop on classes */
  BFT_REALLOC(name, strlen("Ga_DCH")+1 + 2, char);
  strcpy(name, "Ga_DCH");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdch[i] -1 ]-1 ];
    vars->propce[n] = igmdch[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DCH");
  }

 /* IGMDV1 loop on classes */
  BFT_REALLOC(name, strlen("Ga_DV1")+1 + 2, char);
  strcpy(name, "Ga_DV1");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdv1[i] -1 ]-1 ];
    vars->propce[n] = igmdv1[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV1");
  }

 /* IGMDV2 loop on classes */
  BFT_REALLOC(name, strlen("Ga_DV2")+1 + 2, char);
  strcpy(name, "Ga_DV2");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdv2[i] -1 ]-1 ];
    vars->propce[n] = igmdv2[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV2");
  }

 /* IGMHET loop on classes */
  BFT_REALLOC(name, strlen("Ga_HET_O2")+1 + 2, char);
  strcpy(name, "Ga_HET_O2");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmhet[i] -1 ]-1 ];
    vars->propce[n] = igmhet[i];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_HET_O2");
  }

    if (*ihtco2 == 1)
    {
        /* IGHCO2 loop on classes */
        BFT_REALLOC(name, strlen("Ga_HET_CO2")+1 + 2, char);
        strcpy(name, "Ga_HET_CO2");
        for (i = 0; i < *nclass; i++)
        {
            sprintf(snumpp, "%2.2i", i+1);
            strcat(name, snumpp);

            vars->properties_ipp[n] = ipppro[ ipproc[ ighco2[i] -1 ]-1 ];
            vars->propce[n] = ighco2[i];
            BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
            strcpy(vars->properties_name[n++], name);

           strcpy(name, "Ga_HET_CO2");
        }
    }

    if (ippmod[*icp3pl -1] == 1)
    {
        /* IGMSEC loop on classes */
        BFT_REALLOC(name, strlen("Ga_SEC")+1 + 2, char);
        strcpy(name, "Ga_SEC");
        for (i = 0; i < *nclass; i++)
        {
            sprintf(snumpp, "%2.2i", i+1);
            strcat(name, snumpp);

            vars->properties_ipp[n] = ipppro[ ipproc[ igmsec[i] -1 ]-1 ];
            vars->propce[n] = igmsec[i];
            BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
            strcpy(vars->properties_name[n++], name);

            strcpy(name, "Ga_SEC");
        }
    }

  BFT_FREE(name);
  BFT_FREE(snumpp);

  if (n != vars->nsalpp)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
                n, vars->nsalpp);

#if _XML_DEBUG_
  bft_printf("==>UICPPR\n");
  bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
  for (i=0 ; i<vars->nprop ; i++)
    bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
                 "properties_name[%i]: %s\n",
                 i, vars->properties_ipp[i],
                 i, vars->propce[i],
                 i, vars->properties_name[i]);
#endif
}

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const noxyd,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ieqco2,
                                const int *const ihtco2,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if3m,
                                const int *const if3mc2,
                                const int *const if4p2m,
                                const int *const if5m,
                                const int *const if6m,
                                const int *const if7m,
                                const int *const iyco2)
{
  int i;
  char *name = NULL;
  char *snumsca = NULL;

  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0) {
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  } else {
    BFT_MALLOC(vars->label, vars->nscapp, char*);
  }

  /* IHM */
  BFT_MALLOC(vars->label[*ihm -1], strlen("Enthalpy")+1, char);
  strcpy(vars->label[*ihm -1], "Enthalpy");

  /* Loop on classes IH2, INP, IXCH, IXCK */
  BFT_MALLOC(snumsca, 1 + 2, char);

  /* IH2 */
  BFT_MALLOC(name, strlen("ENT_CP")+1 + 2, char);
  strcpy(name, "ENT_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ih2[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ih2[i] -1], name);

    strcpy(name, "ENT_CP");
  }

  /* INP */
  BFT_REALLOC(name, strlen("NP_CP")+1 + 2, char);
  strcpy(name, "NP_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[inp[i] -1], strlen(name)+1, char);
    strcpy(vars->label[inp[i] -1], name);

    strcpy(name, "NP_CP");
  }

  /* IXCH */
  BFT_REALLOC(name, strlen("XCH_CP")+1 + 2, char);
  strcpy(name, "XCH_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ixch[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ixch[i] -1], name);

    strcpy(name, "XCH_CP");
  }

  /* IXCK */
  BFT_REALLOC(name, strlen("XCK_CP")+1 + 2, char);
  strcpy(name, "XCK_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ixck[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ixck[i] -1], name);

    strcpy(name, "XCK_CP");
  }

  /* Loop on coals IFM1 IFM2 */

  BFT_REALLOC(name, strlen("Fr_MV1")+1 + 2, char);
  strcpy(name, "Fr_MV1");
  for (i = 0; i < *ncharb; i++) {
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[if1m[i] -1], strlen(name)+1, char);
    strcpy(vars->label[if1m[i] -1], name);

    strcpy(name, "Fr_MV1");
  }

  BFT_REALLOC(name, strlen("Fr_MV2")+1 + 2, char);
  strcpy(name, "Fr_MV2");
  for (i = 0; i < *ncharb; i++) {
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[if2m[i] -1], strlen(name)+1, char);
    strcpy(vars->label[if2m[i] -1], name);

    strcpy(name, "Fr_MV2");
  }

    /* IF3M */
    BFT_MALLOC(vars->label[*if3m -1], strlen("Fr_HET_O2")+1, char);
    strcpy(vars->label[*if3m -1], "Fr_HET_O2");

    if (*ihtco2 == 1)
    {
        /* IF3MC2 */
        BFT_MALLOC(vars->label[*if3mc2 -1], strlen("Fr_HET_CO2")+1, char);
        strcpy(vars->label[*if3mc2 -1], "Fr_HET_CO2");
    }

    /* IF4P2M */
    BFT_MALLOC(vars->label[*if4p2m -1], strlen("Var_AIR")+1, char);
    strcpy(vars->label[*if4p2m -1], "Var_AIR");

    if (ippmod[*icp3pl -1] == 1)
    {
        /* IXWT */
        BFT_MALLOC(name, strlen("XWT_CP")+1 + 2, char);
        strcpy(name, "XWT_CP");
        for (i = 0; i < *nclass; i++)
        {
            sprintf(snumsca,"%2.2i", i+1);
            strcat(name, snumsca);

            BFT_MALLOC(vars->label[ixwt[i] -1], strlen(name)+1, char);
            strcpy(vars->label[ixwt[i] -1], name);
            strcpy(name, "XWT_CP");
        }

        /* IF5M */
        BFT_MALLOC(vars->label[*if5m -1], strlen("FR_H20")+1, char);
        strcpy(vars->label[*if5m -1], "FR_H20");
    }


    if (*noxyd >= 2)
    {
        /* IF6M */
        BFT_MALLOC(vars->label[*if6m -1], strlen("FR_OXYD2")+1, char);
        strcpy(vars->label[*if6m -1], "FR_OXYD2");
    }

    if (*noxyd == 3)
    {
        /* IF7M */
        BFT_MALLOC(vars->label[*if7m -1], strlen("FR_OXYD3")+1, char);
        strcpy(vars->label[*if7m -1], "FR_OXYD3");
    }

    if (*ieqco2 == 1)
    {
        /* IYCO2 */
        BFT_MALLOC(vars->label[*iyco2 -1], strlen("FR_CO2")+1, char);
        strcpy(vars->label[*iyco2 -1], "FR_CO2");
    }

  BFT_FREE(name);
  BFT_FREE(snumsca);

#if _XML_DEBUG_
  bft_printf("==>UICPSC\n");
  for (i=0; i< vars->nscaus+vars->nscapp; i++)
      bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
#endif

}

/*----------------------------------------------------------------------------
 * Atmospheric flows: read of meteorological file of data
 *
 * Fortran Interface:
 *
 * subroutine uiati1
 * *****************
 * integer         imeteo   <--   on/off index
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiati1, UIATI1) (int *const imeteo)
{
    char *path   = NULL;
    int   status = 0;

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                    "atmospheric_flows",
                                    "read_meteo_data");

    cs_xpath_add_attribute(&path, "status");
    if (cs_gui_get_status(path, &status))
        *imeteo = status;
    BFT_FREE(path);

#if _XML_DEBUG_
    bft_printf("==>UIATI1\n");
    bft_printf("--imeteo  = %i\n", *imeteo);
#endif
}

/*----------------------------------------------------------------------------
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for physical properties.
 *
 * Fortran Interface:
 *
 * subroutine uiatpr
 * *****************
 * integer         nsalpp   -->
 * integer         nsalto   -->
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         ipppro   -->
 * integer         ipproc   -->
 * integer         itempc   -->   index for real temperature
 * integer         iliqwt   -->   index for liquid water
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatpr, UIATPR) (const int *const nsalpp,
                                const int *const nsalto,
                                const int *const ippmod,
                                const int *const iatmos,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const itempc,
                                const int *const iliqwt)
{
    int n;
    cs_var_t *vars = cs_glob_var;

    n = vars->nprop;
    vars->nprop  += *nsalpp;
    vars->nsalpp  = *nsalpp;

    BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
    BFT_REALLOC(vars->propce,          vars->nprop, int);
    BFT_REALLOC(vars->properties_name, vars->nprop, char*);

    /* itempc */
    vars->properties_ipp[n] = ipppro[ ipproc[ *itempc -1 ]-1 ];
    vars->propce[n] = *itempc;
    BFT_MALLOC(vars->properties_name[n], strlen("real_temperature")+1, char);
    strcpy(vars->properties_name[n++], "real_temperature");

    if (ippmod[*iatmos -1] == 2)
    {
        /* iliqwt */
        vars->properties_ipp[n] = ipppro[ ipproc[ *iliqwt -1 ]-1 ];
        vars->propce[n] = *iliqwt;
        BFT_MALLOC(vars->properties_name[n], strlen("liquid_water")+1, char);
        strcpy(vars->properties_name[n++], "liquid_water");
    }
#if _XML_DEBUG_
    {
        int i;
        bft_printf("==>UIATPR\n");
        bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
        for (i=0 ; i<vars->nprop ; i++)
            bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
                       "properties_name[%i]: %s\n",
                       i, vars->properties_ipp[i],
                       i, vars->propce[i],
                       i, vars->properties_name[i]);
    }
#endif
}

/*----------------------------------------------------------------------------
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for models scalars.
 *
 * Fortran Interface:
 *
 * subroutine uiatsc
 * *****************
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         itempp   -->   index for potential temperature
 * integer         itempl   -->   index for liquid potential temperature
 * integer         itotwt   -->   index for total water content
 * integer         intdrp   -->   index for total number of droplets
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatsc, UIATSC) (const int *const ippmod,
                                const int *const iatmos,
                                const int *const itempp,
                                const int *const itempl,
                                const int *const itotwt,
                                const int *const intdrp)
{
    cs_var_t  *vars = cs_glob_var;

    if (vars->nscaus > 0)
    {
        BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
    }
    else
    {
        BFT_MALLOC(vars->label, vars->nscapp, char*);
    }

    if (ippmod[*iatmos -1] == 1)
    {
        /* itempp */
        BFT_MALLOC(vars->label[*itempp -1], strlen("potential_temperature")+1, char);
        strcpy(vars->label[*itempp -1], "potential_temperature");
    }
    else if (ippmod[*iatmos -1] == 2)
    {
        /* itempl */
        BFT_MALLOC(vars->label[*itempl -1], strlen("liquid_potential_temperature")+1, char);
        strcpy(vars->label[*itempl -1], "liquid_potential_temperature");

        /* itotwt */
        BFT_MALLOC(vars->label[*itotwt -1], strlen("total_water")+1, char);
        strcpy(vars->label[*itotwt -1], "total_water");

        /* intdrp */
        BFT_MALLOC(vars->label[*intdrp -1], strlen("number_of_droplets")+1, char);
        strcpy(vars->label[*intdrp -1], "number_of_droplets");
    }
#if _XML_DEBUG_
    {
        int i;
        bft_printf("==>UIATSC\n");
        for (i=0; i< vars->nscaus+vars->nscapp; i++)
            bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
    }
#endif
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo          -->  thermophysical model
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_thermophysical_model(const char *const model_thermo)
{
  char *model = NULL;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", model_thermo);
  cs_xpath_add_attribute(&path, "model");

  model = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return model;
}

/*-----------------------------------------------------------------------------
 * Return 1 if a specific physics model is activated. Store in the global
 * structure vars:
 *   vars->model         <= thermophysical model
 *   vars->model_value   <= related model name
 *----------------------------------------------------------------------------*/

int
cs_gui_get_activ_thermophysical_model(void)
{
    int i, isactiv = 0;
    char *value = NULL;

    cs_var_t  *vars = cs_glob_var;

    const char *name[] = { "pulverized_coal",
                           "gas_combustion",
                           "joule_effect",
                           "atmospheric_flows" };
    int name_nbr = sizeof(name) / sizeof(name[0]);

    if (vars->model != NULL && vars->model_value != NULL)
    {
        isactiv = 1;
        return isactiv;
    }
    else
    {
        vars->model = NULL;
        vars->model_value = NULL;
    }

    for (i = 0; i < name_nbr; i++)
    {
        value = cs_gui_get_thermophysical_model(name[i]);

        if (value && !cs_gui_strcmp(value, "off"))
        {
            BFT_MALLOC(vars->model, strlen(name[i])+1, char);
            strcpy(vars->model, name[i]);

            BFT_MALLOC(vars->model_value, strlen(value)+1, char);
            strcpy(vars->model_value, value);

            isactiv = 1;
            break;
        }
    }

    BFT_FREE(value);

    return isactiv;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
