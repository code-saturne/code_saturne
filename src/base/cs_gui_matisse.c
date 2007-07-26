/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
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
 * Reader of the parameters file: matisse
 *============================================================================*/


#if defined(_CS_HAVE_XML)


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
 * Local headers
 *----------------------------------------------------------------------------*/


#include "cs_base.h"
#include "cs_gui_util.h"


/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/


#include "cs_gui_matisse.h"


/*----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/


/* debugging switch */
#define _XML_DEBUG_ 0


/*============================================================================
 *  Static global variables
 *============================================================================*/


static const char *const cs_matisse_map_type[4]=
{ "inlet_range",
  "outlet_range",
  "network",
  "thermal_capacity" };

static const char *const cs_matisse_map_axis[3]=
{ "line",
  "row",
  "height" };


/*============================================================================
 * Private functions prototypes
 *============================================================================*/


double cs_gui_data_matisse_double(const char *const markup1,
                                  const char *const markup2,
                                  const char *const data);

int cs_gui_data_matisse_int(const char *const markup,
                            const char *const data);

int cs_gui_data_matisse_att_status(const char *const data);

int cs_gui_warehousing_type(void);


/*============================================================================
 * Public functions API Fortran
 *============================================================================*/


/*----------------------------------------------------------------------------
 * Traitement des parametres geometriques de type entier de Matisse
 *----------------------------------------------------------------------------*/


void CS_PROCF (csgein, CSGEIN) (int *const nptran,
                                int *const nplgrs,
                                int *const nelgrs,
                                int *const nchest,
                                int *const netran,
                                int *const itypen)
{
  *nptran = cs_gui_data_matisse_int("compute", "nptran");
  *nplgrs = cs_gui_data_matisse_int("compute", "nplgrs");
  *nelgrs = cs_gui_data_matisse_int("compute", "nelgrs");
  *nchest = cs_gui_data_matisse_int("compute", "nchest");
  *netran = cs_gui_data_matisse_int("compute", "netran");
  *itypen = cs_gui_warehousing_type();

#if _XML_DEBUG_
  bft_printf(_("==>CSGEIN\n"));
  bft_printf(_("--nptran = %d\n"), *nptran);
  bft_printf(_("--nplgrs = %d\n"), *nplgrs);
  bft_printf(_("--nelgrs = %d\n"), *nelgrs);
  bft_printf(_("--nchest = %d\n"), *nchest);
  bft_printf(_("--netran = %d\n"), *netran);
  bft_printf(_("--itypen = %d\n"), *itypen);
#endif
}


/*----------------------------------------------------------------------------
 * Traitement des parametres geometriques de type réel de Matisse
 * non stockés dans les COMMON
 *----------------------------------------------------------------------------*/


void CS_PROCF (csmhdb, CSMHDB) (double *const jeuchr,
                                double *const jeurcl,
                                double *const jeuclr,
                                double *const jeurch,
                                int    *const nechrg,
                                int    *const nergrs,
                                int    *const neclrg,
                                int    *const nergch,
                                double *const hbdtoi,
                                int    *const neciel)
{
  *nechrg = cs_gui_data_matisse_int("mesh", "nechrg");
  *nergrs = cs_gui_data_matisse_int("mesh", "nergrs");
  *neclrg = cs_gui_data_matisse_int("mesh", "neclrg");
  *nergch = cs_gui_data_matisse_int("mesh", "nergch");
  *neciel = cs_gui_data_matisse_int("mesh", "neciel");

  *jeuchr = cs_gui_data_matisse_double("mesh", "geometry", "jeuchr");
  *jeurcl = cs_gui_data_matisse_double("mesh", "geometry", "jeurcl");
  *jeuclr = cs_gui_data_matisse_double("mesh", "geometry", "jeuclr");
  *jeurch = cs_gui_data_matisse_double("mesh", "geometry", "jeurch");
  *hbdtoi = cs_gui_data_matisse_double("mesh", "geometry", "hbdtoi");

#if _XML_DEBUG_
  bft_printf(_("==>CSMHDB\n"));
  bft_printf(_("--nechrg = %i\n"), *nechrg);
  bft_printf(_("--nergrs = %i\n"), *nergrs);
  bft_printf(_("--neclrg = %i\n"), *neclrg);
  bft_printf(_("--nergch = %i\n"), *nergch);
  bft_printf(_("--jeuchr = %f\n"), *jeuchr);
  bft_printf(_("--jeurcl = %f\n"), *jeurcl);
  bft_printf(_("--jeuclr = %f\n"), *jeuclr);
  bft_printf(_("--jeurch = %f\n"), *jeurch);
#endif
}


/*----------------------------------------------------------------------------
 * Traitement des parametres geometriques de type réel de Matisse
 *----------------------------------------------------------------------------*/


void CS_PROCF (csgedb, CSGEDB) (double *const epregi,
                                double *const epchem,
                                double *const hconve,
                                double *const rconve,
                                double *const hchali,
                                double *const hcheva,
                                double *const hfttoi,
                                double *const ptrres,
                                double *const frdtra,
                                double *const plgres,
                                double *const epchel,
                                double *const dmcont)
{
  *epregi = cs_gui_data_matisse_double("compute", "geometry", "epregi");
  *epchem = cs_gui_data_matisse_double("compute", "geometry", "epchem");
  *hconve = cs_gui_data_matisse_double("compute", "geometry", "hconve");
  *rconve = cs_gui_data_matisse_double("compute", "geometry", "rconve");
  *hchali = cs_gui_data_matisse_double("compute", "geometry", "hchali");
  *hcheva = cs_gui_data_matisse_double("compute", "geometry", "hcheva");
  *hfttoi = cs_gui_data_matisse_double("compute", "geometry", "hfttoi");
  *ptrres = cs_gui_data_matisse_double("compute", "geometry", "ptrres");
  *frdtra = cs_gui_data_matisse_double("compute", "geometry", "frdtra");
  *plgres = cs_gui_data_matisse_double("compute", "geometry", "plgres");
  *epchel = cs_gui_data_matisse_double("compute", "geometry", "epchel");
  *dmcont = cs_gui_data_matisse_double("compute", "geometry", "dmcont");

#if _XML_DEBUG_
  bft_printf(_("==>CSGEDB\n"));
  bft_printf(_("--epregi = %f\n"), *epregi);
  bft_printf(_("--epchem = %f\n"), *epchem);
  bft_printf(_("--hconve = %f\n"), *hconve);
  bft_printf(_("--rconve = %f\n"), *rconve);
  bft_printf(_("--hchali = %f\n"), *hchali);
  bft_printf(_("--hcheva = %f\n"), *hcheva);
  bft_printf(_("--hfttoi = %f\n"), *hfttoi);
  bft_printf(_("--ptrres = %f\n"), *ptrres);
  bft_printf(_("--frdtra = %f\n"), *frdtra);
  bft_printf(_("--plgres = %f\n"), *plgres);
  bft_printf(_("--epchel = %f\n"), *epchel);
  bft_printf(_("--dmcont = %f\n"), *dmcont);
#endif
}


/*----------------------------------------------------------------------------
 * Traitement des parametres physiques de type double precision de Matisse
 *----------------------------------------------------------------------------*/


void CS_PROCF (csphdb, CSPHDB) (double *const dtdtmx,
                                double *const puicon,
                                double *const tinit,
                                double *const tcrit,
                                double *const emicon,
                                double *const emimur,
                                double *const hepcnt,
                                double *const dhpcnt,
                                double *const debmas,
                                double *const pdccha,
                                double *const pdcfch,
                                double *const dhchea,
                                double *const sdchea,
                                double *const pdcche,
                                double *const pdccch,
                                double *const dhches,
                                double *const sdches,
                                double *const pdcalg,
                                double *const pdcatv,
                                double *const argamt,
                                double *const pdcslg,
                                double *const pdcstv,
                                double *const argavl,
                                double *const amppdc,
                                double *const dhalve,
                                double *const hreso,
                                double *const hplen,
                                double *const dpvent)
{
  *dtdtmx = cs_gui_data_matisse_double("compute", "physical_model", "dtdtmx");
  *puicon = cs_gui_data_matisse_double("compute", "physical_model", "puicon");
  *tinit  = cs_gui_data_matisse_double("compute", "physical_model", "tinit");
  *tcrit  = cs_gui_data_matisse_double("compute", "physical_model", "tcrit");
  *emicon = cs_gui_data_matisse_double("compute", "physical_model", "emicon");
  *emimur = cs_gui_data_matisse_double("compute", "physical_model", "emimur");
  *hepcnt = cs_gui_data_matisse_double("compute", "physical_model", "hepcnt");
  *dhpcnt = cs_gui_data_matisse_double("compute", "physical_model", "dhpcnt");
  *debmas = cs_gui_data_matisse_double("compute", "physical_model", "debmas");
  *pdccha = cs_gui_data_matisse_double("compute", "physical_model", "pdccha");
  *pdcfch = cs_gui_data_matisse_double("compute", "physical_model", "pdcfch");
  *dhchea = cs_gui_data_matisse_double("compute", "physical_model", "dhchea");
  *sdchea = cs_gui_data_matisse_double("compute", "physical_model", "sdchea");
  *pdcche = cs_gui_data_matisse_double("compute", "physical_model", "pdcche");
  *pdccch = cs_gui_data_matisse_double("compute", "physical_model", "pdccch");
  *dhches = cs_gui_data_matisse_double("compute", "physical_model", "dhches");
  *sdches = cs_gui_data_matisse_double("compute", "physical_model", "sdches");
  *pdcalg = cs_gui_data_matisse_double("compute", "physical_model", "pdcalg");
  *pdcatv = cs_gui_data_matisse_double("compute", "physical_model", "pdcatv");
  *argamt = cs_gui_data_matisse_double("compute", "physical_model", "argamt");
  *pdcslg = cs_gui_data_matisse_double("compute", "physical_model", "pdcslg");
  *pdcstv = cs_gui_data_matisse_double("compute", "physical_model", "pdcstv");
  *argavl = cs_gui_data_matisse_double("compute", "physical_model", "argavl");
  *amppdc = cs_gui_data_matisse_double("compute", "physical_model", "amppdc");
  *dhalve = cs_gui_data_matisse_double("compute", "physical_model", "dhalve");
  *hreso  = cs_gui_data_matisse_double("compute", "physical_model", "hreso");
  *hplen  = cs_gui_data_matisse_double("compute", "physical_model", "hplen");
  *dpvent = cs_gui_data_matisse_double("compute", "physical_model", "dpvent");

#if _XML_DEBUG_
  bft_printf(_("==>CSPHDB\n"));
  bft_printf(_("--dtdtmx = %f\n"), *dtdtmx);
  bft_printf(_("--puicon = %f\n"), *puicon);
  bft_printf(_("--tinit  = %f\n"), *tinit);
  bft_printf(_("--tcrit  = %f\n"), *tcrit);
  bft_printf(_("--emicon = %f\n"), *emicon);
  bft_printf(_("--emimur = %f\n"), *emimur);
  bft_printf(_("--hepcnt = %f\n"), *hepcnt);
  bft_printf(_("--dhpcnt = %f\n"), *dhpcnt);
  bft_printf(_("--debmas = %f\n"), *debmas);
  bft_printf(_("--pdccha = %f\n"), *pdccha);
  bft_printf(_("--pdcfch = %f\n"), *pdcfch);
  bft_printf(_("--dhchea = %f\n"), *dhchea);
  bft_printf(_("--sdchea = %f\n"), *sdchea);
  bft_printf(_("--pdcche = %f\n"), *pdcche);
  bft_printf(_("--pdccch = %f\n"), *pdccch);
  bft_printf(_("--dhches = %f\n"), *dhches);
  bft_printf(_("--sdches = %f\n"), *sdches);
  bft_printf(_("--pdcalg = %f\n"), *pdcalg);
  bft_printf(_("--pdcatv = %f\n"), *pdcatv);
  bft_printf(_("--argamt = %f\n"), *argamt);
  bft_printf(_("--pdcslg = %f\n"), *pdcslg);
  bft_printf(_("--pdcstv = %f\n"), *pdcstv);
  bft_printf(_("--argavl = %f\n"), *argavl);
  bft_printf(_("--amppdc = %f\n"), *amppdc);
  bft_printf(_("--dhalve = %f\n"), *dhalve);
  bft_printf(_("--hreso  = %f\n"), *hreso);
  bft_printf(_("--hplen  = %f\n"), *hplen);
  bft_printf(_("--dpvent = %f\n"), *dpvent);
#endif
}


/*----------------------------------------------------------------------------
 * Traitement des parametres physiques de type attribut (sens XML) de Matisse
 *----------------------------------------------------------------------------*/


void CS_PROCF (csphat, CSPHAT)(int *const imdcnt,
                               int *const icofor,
                               int *const iconlg,
                               int *const ialveo)
{
  *imdcnt = cs_gui_data_matisse_att_status("imdcnt");
  *icofor = cs_gui_data_matisse_att_status("icofor");
  *iconlg = cs_gui_data_matisse_att_status("iconlg");
  *ialveo = cs_gui_data_matisse_att_status("ialveo");

#if _XML_DEBUG_
  bft_printf(_("==>CSPHAT\n"));
  bft_printf(_("--imdcnt = %d\n"), *imdcnt);
  bft_printf(_("--icofor = %d\n"), *icofor);
  bft_printf(_("--iconlg = %d\n"), *iconlg);
  bft_printf(_("--ialveo = %d\n"), *ialveo);
#endif
}


/*----------------------------------------------------------------------------
 * Test si la balise Matisse se trouve dans le fichier XML
 *----------------------------------------------------------------------------*/


void CS_PROCF(csmtpr,CSMTPR)(int *imatis)
{
 char *path;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "matisse");

  if (cs_gui_get_nb_element(path) > 0 )
    *imatis = 1;
  else
    *imatis = 0;

#if _XML_DEBUG_
  bft_printf(_("==>CSMTPR\n"));
  bft_printf(_("--imatis = %d\n"), *imatis);
#endif

  BFT_FREE(path);
}


/*----------------------------------------------------------------------------
 * Calcul le nombre de zones d'une carte et d'une direction donnée
 *----------------------------------------------------------------------------*/


void CS_PROCF(csnbmp,CSNBMP) (int *const direction,
                              int *const carte,
                              int *const nb)
{
  char *path;
  int icarte = (*carte)-1;
  int idirec = (*direction)-1;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "matisse", "compute", "map");
  if (!cs_gui_strcmp(cs_matisse_map_type[icarte], "thermal_capacity"))
    cs_xpath_add_element(&path, "headloss");
  cs_xpath_add_element(&path, cs_matisse_map_type[icarte]);
  cs_xpath_add_element(&path, cs_matisse_map_axis[idirec]);
  cs_xpath_add_element(&path, "area");

  *nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf(_("==>CSNBMP\n"));
  bft_printf(_("--Zones number for the map %i in the direction %i: %i\n"),
                *carte, *direction, *nb);
#endif
}


/*----------------------------------------------------------------------------
 * Rempli les cartes 2D et 3D de pertes de charge et de puissance thermique
 *----------------------------------------------------------------------------*/


void CS_PROCF(csdfmp,CSDFMP) (   int *const zone,
                                 int *const direction,
                                 int *const carte,
                              double *const min,
                              double *const max,
                              double *const value)
{
  char *path;
  char *pathtmp;
  int icarte = (*carte)-1;
  int idirec = (*direction)-1;
  int izone  = (*zone)-1;

  /* Construction de la requete */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "matisse", "compute", "map");
  if (!cs_gui_strcmp(cs_matisse_map_type[icarte], "thermal_capacity"))
    cs_xpath_add_element(&path, "headloss");
  cs_xpath_add_element(&path, cs_matisse_map_type[icarte]);
  cs_xpath_add_element(&path, cs_matisse_map_axis[idirec]);
  cs_xpath_add_element_num(&path, "area" , izone+1);


  /* Détermination de min */
  BFT_MALLOC(pathtmp, strlen(path)+1, char);
  strcpy(pathtmp, path);
  cs_xpath_add_element(&path, "min");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, min))
    bft_error(__FILE__, __LINE__, 0,
              _("Missing 'min' markup for xpath : %s\n"), path);


  /* Détermination de max */
  strcpy(path, pathtmp);
  cs_xpath_add_element(&path, "max");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, max))
    bft_error (__FILE__, __LINE__, 0,
              _("Missing 'max' markup for xpath : %s\n"), path);

  /* Détermination de value */
  if (cs_gui_strcmp(cs_matisse_map_type[icarte], "thermal_capacity")) {
    strcpy(path, pathtmp);
    cs_xpath_add_element(&path, "value");
    cs_xpath_add_function_text(&path);

    if (!cs_gui_get_double(path, value))
      bft_error(__FILE__, __LINE__, 0,
                _("Missing 'value' markup for xpath : %s\n"), path);

  }

  BFT_FREE(path);
  BFT_FREE(pathtmp);

#if _XML_DEBUG_
  bft_printf(_("==>CSDFMP\n"));
  if (cs_gui_strcmp(cs_matisse_map_type[icarte], "thermal_capacity"))
    bft_printf(_("--Zone %i Direction %i Map %i = %f %f %f \n"),
                  *zone, *direction, *carte, *min, *max, *value);
  else
    bft_printf(_("--Zone %i Direction %i Map %i = %f %f \n"),
                  *zone, *direction, *carte, *min, *max );
#endif
}


/*-----------------------------------------------------------------------------
 * Retourne une donnee matisse de type double
 *----------------------------------------------------------------------------*/


double cs_gui_data_matisse_double(const char *const markup1,
                                  const char *const markup2,
                                  const char *const data)
{
  char   *path;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "matisse", markup1, markup2, data);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}


/*-----------------------------------------------------------------------------
 * Retourne une donnee matisse de type entier
 *----------------------------------------------------------------------------*/


int cs_gui_data_matisse_int(const char *const markup,
                            const char *const data)
{
  char *path;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "matisse", markup, "geometry", data);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_int(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}


/*-----------------------------------------------------------------------------
 * Retourne une donnee matisse de type entier (attribut status XML)
 *----------------------------------------------------------------------------*/


int cs_gui_data_matisse_att_status(const char *const data)
{
  char *path;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "matisse", "compute", "physical_model", data);
  cs_xpath_add_attribute(&path, "status");

  if (!cs_gui_get_status(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}


/*-----------------------------------------------------------------------------
 * Retourne Le type d'entreposage (1 pour Emm, 0 pour Vault)
 *----------------------------------------------------------------------------*/


int cs_gui_warehousing_type(void)
{
  char *path;
  char *value;
  int   intval;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "matisse", "compute", "geometry", "typent");
  cs_xpath_add_attribute(&path, "label");

  value = cs_gui_get_attribute_value(path);

  if (cs_gui_strcmp(value, "vault"))
    intval = 0;

  else if (cs_gui_strcmp(value,"emm"))
    intval = 1;

  else
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  BFT_FREE(value);

  return intval;
}


/*----------------------------------------------------------------------------*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_HAVE_XML */

