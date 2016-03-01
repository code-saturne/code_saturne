#ifndef __CS_GUI_BOUNDARY_CONDITION_H__
#define __CS_GUI_BOUNDARY_CONDITION_H__

/*============================================================================
 * Management of the GUI parameters file: boundary conditions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * MEI library headers
 *----------------------------------------------------------------------------*/

#include "mei_evaluate.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structures associated to boundary conditions definition
 *----------------------------------------------------------------------------*/

typedef struct {
  double val1;             /* fortran array RCODCL(.,.,1) mapping             */
  double val2;             /* fortran array RCODCL(.,.,2) mapping             */
  double val3;             /* fortran array RCODCL(.,.,3) mapping             */
} cs_val_t;

typedef struct {
  int        read_data;    /* 1 if profile is calculated from data            */
  int        automatic;    /* 1 if nature of the boundary is automatic        */
} cs_meteo_t;

typedef struct {
  char      **label;       /* label for each boundary zone                    */
  char      **nature;      /* nature for each boundary zone                   */
  int        *iqimp;       /* 1 if a flow rate is applied                     */
  int        *ientfu;      /* 1  for a fuel flow inlet (gas combustion - D3P) */
  int        *ientox;      /* 1 for an air flow inlet (gas combustion - D3P)  */
  int        *ientgb;      /* 1 for burned gas inlet (gas combustion)         */
  int        *ientgf;      /* 1 for unburned gas inlet (gas combustion)       */
  int        *ientat;      /* 1 if inlet for oxydant (coal combustion)        */
  int        *ientcp;      /* 1 if inlet for oxydant+coal (coal combustion)   */
  int        *icalke;      /* automatic boundaries for turbulent variables    */
  double     *qimp;        /* oxydant flow rate (coal combustion)             */
  int        *inmoxy;      /* oxydant number (coal combustion)                */
  double     *timpat;      /* inlet temperature of oxydant (coal combustion)  */
  double     *tkent;       /* inlet temperature (gas combustion)              */
  double    **qimpcp;      /* inlet coal flow rate (coal combustion)          */
  double    **timpcp;      /* inlet coal temperature (coal combustion)        */
  double     *fment;       /* Mean Mixture Fraction at Inlet (gas combustion) */
  int        *itype;       /* type of inlet/outlet (compressible model)       */
  double     *prein;       /* inlet pressure (compressible model)             */
  double     *rhoin;       /* inlet density  (compressible model)             */
  double     *tempin;      /* inlet temperature (compressible model)          */
  double     *entin;       /* inlet total energy (compressible model)         */
  double     *preout;      /* outlet pressure for subsonic(compressible model)*/
  double     *dh;          /* inlet hydraulic diameter                        */
  double     *xintur;      /* inlet turbulent intensity                       */
  int       **type_code;   /* type of boundary for each variables             */
  cs_val_t  **values;      /* fortran array RCODCL mapping                    */
  double   ***distch;      /* ratio for each coal                             */
  double     *rough;       /* roughness size                                  */
  double     *norm;        /* norm of velocity vector                         */
  double     *dirx;        /* directions x inlet velocity                     */
  double     *diry;        /* directions y inlet velocity                     */
  double     *dirz;        /* directions z inlet velocity                     */
  mei_tree_t **velocity;   /* formula for norm or mass flow rate of velocity  */
  mei_tree_t **direction;  /* formula for direction of velocity               */
  cs_meteo_t  *meteo;      /* inlet or outlet info for atmospheric flow       */
  mei_tree_t ***scalar;    /* formula for scalar (neumann, dirichlet or
                              exchange coefficient)*/
  mei_tree_t **headLoss;   /* formula for head loss (free inlet/outlet)       */
  mei_tree_t **groundwat;  /* formula for hydraulic head (groundwater)        */
} cs_boundary_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main boundaries structure */

extern cs_boundary_t *boundaries;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/
/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Remember: rdoccl[k][j][i] = rcodcl[ k * dim1 *dim2 + j *dim1 + i]
 *
 * Fortran Interface:
 *
 * subroutine uiclim
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const int  *ntcabs,
                               const int  *nfabor,
                               const int  *idarcy,
                               const int  *nozppm,
                               const int  *ncharm,
                               const int  *ncharb,
                               const int  *nclpch,
                               const int  *iindef,
                               const int  *ientre,
                               const int  *iesicf,
                               const int  *isspcf,
                               const int  *iephcf,
                               const int  *isopcf,
                               const int  *iparoi,
                               const int  *iparug,
                               const int  *isymet,
                               const int  *isolib,
                               const int  *ifrent,
                               const int  *ifresf,
                               int        *iqimp,
                               int        *icalke,
                               int        *ientat,
                               int        *ientcp,
                               int        *inmoxy,
                               int        *ientox,
                               int        *ientfu,
                               int        *ientgf,
                               int        *ientgb,
                               int        *iprofm,
                               int        *itypfb,
                               int        *izfppp,
                               int        *icodcl,
                               double     *dtref,
                               double     *ttcabs,
                               double     *surfbo,
                               double     *cdgfbo,
                               double     *qimp,
                               double     *qimpat,
                               double     *qimpcp,
                               double     *dh,
                               double     *xintur,
                               double     *timpat,
                               double     *timpcp,
                               double     *tkent,
                               double     *fment,
                               double     *distch,
                               int        *nvarcl,
                               double     *rcodcl);

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * integer          nfabor  <-- number of boundary faces
 * integer          nozppm  <-- max number of boundary conditions zone
 * integer          iindef  <-- type of boundary: not defined
 * integer          ientre  <-- type of boundary: inlet
 * integer          iesicf  --> type of boundary: imposed inlet (compressible)
 * integer          isspcf  --> type of boundary: supersonic outlet (compressible)
 * integer          iephcf  --> type of boundary: subsonic inlet at imposed
 *                              total pressure and enthalpy (compressible)
 * integer          isopcf  --> type of boundary: subsonic outlet (compressible)
 * integer          iparoi  <-- type of boundary: smooth wall
 * integer          iparug  <-- type of boundary: rough wall
 * integer          isymet  <-- type of boundary: symetry
 * integer          isolib  <-- type of boundary: outlet
 * integer          ifrent  <-- type of boundary: free inlet outlet
 * integer          ifresf  <-- type of boundary: free surface
 * integer          iale    <-- ale module activated
 * integer          itypfb  <-- type of boundary for each face
 * integer          izfppp  <-- zone number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int  *nfabor,
                               const int  *nozppm,
                               const int  *iindef,
                               const int  *ientre,
                               const int  *iesicf,
                               const int  *iephcf,
                               const int  *isspcf,
                               const int  *isopcf,
                               const int  *iparoi,
                               const int  *iparug,
                               const int  *isymet,
                               const int  *isolib,
                               const int  *ifrent,
                               const int  *ifresf,
                               const int  *iale,
                               int        *itypfb,
                               int        *izfppp);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return number of boundary regions definition
 *----------------------------------------------------------------------------*/

int
cs_gui_boundary_zones_number(void);

/*-----------------------------------------------------------------------------
 * Return the nature of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *
cs_gui_boundary_zone_nature(int  ith_zone);

/*-----------------------------------------------------------------------------
 * Return the label of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *
cs_gui_boundary_zone_label(int  ith_zone);

/*-----------------------------------------------------------------------------
 * Return the zone number of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

int
cs_gui_boundary_zone_number(int  ith_zone);

/*-----------------------------------------------------------------------------
 * Return the description of a boundary zone
 *
 * parameters:
 *   label <--  label of boundary zone
 *----------------------------------------------------------------------------*/

char *
cs_gui_boundary_zone_localization(const char *label);

/*-----------------------------------------------------------------------------
 * Helper to get the face list for the izone
 *
 * parameters:
 *   izone   <--  zone index
 *   label   <--  boundary label
 *   nfabor  <--  number of boundary faces
 *   nozppm  <--  max number of boundary zone for preefined physics
 *   faces   -->  number of faces
 *----------------------------------------------------------------------------*/

cs_lnum_t*
cs_gui_get_faces_list(int          izone,
                      const char  *label,
                      cs_lnum_t    nfabor,
                      int          nozppm,
                      cs_lnum_t   *faces);

/*----------------------------------------------------------------------------
 * Free boundary conditions structures
 *
 * parameters:
 *   ncharb  <-- number of coals
 *----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(const int  *ncharb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GUI_BOUNDARY_CONDITION_H__ */
