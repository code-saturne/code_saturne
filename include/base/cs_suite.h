/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_SUITE_H__
#define __CS_SUITE_H__

/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Error codes */

#define CS_SUITE_SUCCES          0 /* Success */
#define CS_SUITE_ERR_NUM_FIC    -1 /* No restart file for the given number */
#define CS_SUITE_ERR_SUPPORT    -2 /* Undefined location / incorrect size */
#define CS_SUITE_ERR_TYPE_VAL   -3 /* Unknown or unexpected value type */
#define CS_SUITE_ERR_NBR_VAL    -4 /* Number of values does not match */
#define CS_SUITE_ERR_MODE       -5 /* Incompatible access mode */
#define CS_SUITE_ERR_EXISTE     -6 /* Section not available */

/*============================================================================
 * Local type definitions
 *============================================================================*/

/* Read or write mode */

typedef enum {

  CS_SUITE_MODE_LECTURE,         /* Read mode */
  CS_SUITE_MODE_ECRITURE         /* Write mode */

} cs_suite_mode_t;

/* Predefined location types for a given section */

typedef enum {

  CS_SUITE_SUPPORT_SCAL,         /* Scalare (no location) */
  CS_SUITE_SUPPORT_CEL,          /* Values defined at cells */
  CS_SUITE_SUPPORT_FAC_INT,      /* Values defined at interior faces */
  CS_SUITE_SUPPORT_FAC_BRD,      /* Values defined at boundary faces */
  CS_SUITE_SUPPORT_SOM           /* Values defined at vertices */

} cs_suite_support_t;

/*
  Pointeur associated with a restart file structure. The structure itself
  is defined in "cs_suite.c", and is opaque outside that unit.
*/

typedef struct _cs_suite_t cs_suite_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Open a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE OPNSUI (NOMSUI, LNGNOM, IREAWR, NUMSUI, IERROR)
 * *****************
 *
 * CHARACTER*       NOMSUI      : --> : Restart file name
 * INTEGER          LNGNOM      : --> : Restart file name length
 * INTEGER          IREAWR      : --> : 1: read; 2: write
 * INTEGER          NUMSUI      : <-- : Number of opened restart file
 * INTEGER          IERROR      : <-- : 0: success; < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (opnsui, OPNSUI)
(
 const char       *nomsui,
 const cs_int_t   *lngnom,
 const cs_int_t   *ireawr,
       cs_int_t   *numsui,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Close a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE CLSSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : <-> : number of restart file to close
 * INTEGER          IERROR      : <-- : 0: success; < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (clssui, CLSSUI)
(
 const cs_int_t   *numsui,
       cs_int_t   *ierror
);

/*----------------------------------------------------------------------------
 * Check the locations associated with a restart file.
 *
 * For each type of entity, return 1 if the associated number of entities
 * matches the current value (and so that we consider the mesh locations are
 * the same), 0 otherwise.
 *
 * Fortran interface
 *
 * SUBROUTINE TSTSUI (NUMSUI, INDCEL, INDFAC, INDFBR, INDSOM)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * INTEGER          INDCEL      : <-- : Matching cells flag
 * INTEGER          INDFAC      : <-- : Matching interior faces flag
 * INTEGER          INDFBR      : <-- : Matching boundary faces flag
 * INTEGER          INDSOM      : <-- : Matching vertices flag
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstsui, TSTSUI)
(
 const cs_int_t  *numsui,
       cs_int_t  *indcel,
       cs_int_t  *indfac,
       cs_int_t  *indfbr,
       cs_int_t  *indsom
);

/*----------------------------------------------------------------------------
 * Print index associated with a restart file in read mode
 *
 * Fortran interface
 *
 * SUBROUTINE INFSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 *----------------------------------------------------------------------------*/

void CS_PROCF (infsui, INFSUI)
(
 const cs_int_t  *numsui
);

/*----------------------------------------------------------------------------
 * Read a section from a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE LECSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * CHARACTER*       NOMRUB      : --> : Section name
 * INTEGER          LNGNOM      : --> : Section name length
 * INTEGER          ITYSUP      : --> : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * INTEGER          NBVENT      : --> : N. values per location entity
 * INTEGER          IRTYPE      : --> : 1 for integers, 2 for double precision
 * (?)              TABVAR      : <-> : Array of values to read
 * INTEGER          IERROR      : <-- : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecsui, LECSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *nbvent,
 const cs_int_t   *irtype,
       void       *tabvar,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Write a section to a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE ECRSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * CHARACTER*       NOMRUB      : --> : Section name
 * INTEGER          LNGNOM      : --> : Section name length
 * INTEGER          ITYSUP      : --> : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * INTEGER          NBVENT      : --> : N. values per location entity
 * INTEGER          IRTYPE      : --> : 1 for integers, 2 for double precision
 * (?)              TABVAR      : --> : Array of values to write
 * INTEGER          IERROR      : <-- : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrsui, ECRSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *nbvent,
 const cs_int_t   *irtype,
 const void       *tabvar,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a restart file
 *
 * parameters:
 *   nom  <-- file name
 *   mode <-- read or write
 *
 * returns:
 *   pointer to initialized restart file structure
 *----------------------------------------------------------------------------*/

cs_suite_t *
cs_suite_cree(const char             *nom,
              cs_suite_mode_t         mode);

/*----------------------------------------------------------------------------
 * Destroy structure associated with a restart file (and close the file).
 *
 * parameters:
 *   suite <-- pointer to restart file structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_suite_t *
cs_suite_detruit(cs_suite_t  *suite);

/*----------------------------------------------------------------------------
 * Check the locations associated with a restart file.
 *
 * For each type of entity, the correspondinf flag is set to true if the
 * associated number of entities matches the current value (and so that we
 * consider the mesh locations are the same), false otherwise.
 *
 * parameters:
 *   suite        <-- associated restart file pointer
 *   corresp_cell <-- matching cells flag
 *   corresp_fac  <-- matching interior faces flag
 *   corresp_fbr  <-- matching boundary faces flag
 *   corresp_som  <-- matching vertices flag
 *----------------------------------------------------------------------------*/

void
cs_suite_verif_support_base(const cs_suite_t  *suite,
                            cs_bool_t         *corresp_cel,
                            cs_bool_t         *corresp_fac,
                            cs_bool_t         *corresp_fbr,
                            cs_bool_t         *corresp_som);

/*----------------------------------------------------------------------------
 * Add a location definition.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   location_name   <-- name associated with the location
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers, or NULL
 *
 * returns:
 *   the location id assigned to the location, or -1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_suite_ajoute_support(cs_suite_t        *suite,
                        const char        *location_name,
                        fvm_gnum_t         n_glob_ents,
                        fvm_lnum_t         n_ents,
                        const fvm_gnum_t  *ent_global_num);

/*----------------------------------------------------------------------------
 * Print the index associated with a restart file in read mode
 *
 * parameters:
 *   suite <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_suite_affiche_index(const cs_suite_t  *suite);

/*----------------------------------------------------------------------------
 * Read a section from a restart file.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   nom_rub         <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   typ_val         <-- value type
 *   val             --> array of values
 *
 * returns: 0 (CS_SUITE_SUCCES) in case of success,
 *          or error code (CS_SUITE_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_suite_lit_rub(cs_suite_t  *suite,
                 const char  *nom_rub,
                 int          location_id,
                 cs_int_t     n_location_vals,
                 cs_type_t    typ_val,
                 void        *val);

/*----------------------------------------------------------------------------
 * Write a section to a restart file.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   nom_rub         <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   typ_val         <-- value type
 *   val             <-- array of values
 *----------------------------------------------------------------------------*/

void
cs_suite_ecr_rub(cs_suite_t   *suite,
                 const char   *nom_rub,
                 int           location_id,
                 cs_int_t      n_location_vals,
                 cs_type_t     typ_val,
                 const void   *val);

/*----------------------------------------------------------------------------
 * Initialize the restart file Fortran API
 *----------------------------------------------------------------------------*/

void
cs_suite_f77_api_init(void);

/*----------------------------------------------------------------------------
 * Finalize the restart file Fortran API
 *----------------------------------------------------------------------------*/

void
cs_suite_f77_api_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SUITE_H__ */
