#ifndef __CS_RESTART_H__
#define __CS_RESTART_H__

/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Error codes */

#define CS_RESTART_SUCCES         0 /* Success */
#define CS_RESTART_ERR_FILE_NUM  -1 /* No restart file for the given number */
#define CS_RESTART_ERR_LOCATION  -2 /* Undefined location / incorrect size */
#define CS_RESTART_ERR_VAL_TYPE  -3 /* Unknown or unexpected value type */
#define CS_RESTART_ERR_N_VALS    -4 /* Number of values does not match */
#define CS_RESTART_ERR_MODE      -5 /* Incompatible access mode */
#define CS_RESTART_ERR_EXISTS    -6 /* Section not available */

/*============================================================================
 * Local type definitions
 *============================================================================*/

/* Read or write mode */

typedef enum {

  CS_RESTART_MODE_READ,         /* Read mode */
  CS_RESTART_MODE_WRITE         /* Write mode */

} cs_restart_mode_t;

/* Predefined location types for a given section */

typedef enum {

  CS_RESTART_LOCATION_NONE,        /* Global (no location) */
  CS_RESTART_LOCATION_CELL,        /* Values defined at cells */
  CS_RESTART_LOCATION_I_FACE,      /* Values defined at interior faces */
  CS_RESTART_LOCATION_B_FACE,      /* Values defined at boundary faces */
  CS_RESTART_LOCATION_VERTEX       /* Values defined at vertices */

} cs_restart_location_t;

/*
  Pointeur associated with a restart file structure. The structure itself
  is defined in "cs_restart.c", and is opaque outside that unit.
*/

typedef struct _cs_restart_t cs_restart_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Indicate if a restart directory is present.
 *
 * Fortran interface
 *
 * subroutine indsui (isuite)
 * *****************
 *
 * integer          isuite      : --> : 1 for restart, 0 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF (indsui, INDSUI)
(
 cs_int_t   *isuite
);

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
 * Check if we have a restart directory.
 *
 * returns:
 *   1 if a restart directory is present, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_restart_present(void);

/*----------------------------------------------------------------------------
 * Initialize a restart file
 *
 * parameters:
 *   name <-- file name
 *   mode <-- read or write
 *
 * returns:
 *   pointer to initialized restart file structure
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_create(const char         *name,
                  cs_restart_mode_t   mode);

/*----------------------------------------------------------------------------
 * Destroy structure associated with a restart file (and close the file).
 *
 * parameters:
 *   restart <-- pointer to restart file structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_destroy(cs_restart_t  *restart);

/*----------------------------------------------------------------------------
 * Check the locations associated with a restart file.
 *
 * For each type of entity, the corresponding flag is set to true if the
 * associated number of entities matches the current value (and so that we
 * consider the mesh locations are the same), false otherwise.
 *
 * parameters:
 *   restart      <-- associated restart file pointer
 *   match_cell   <-- matching cells flag
 *   match_i_face <-- matching interior faces flag
 *   match_b_face <-- matching boundary faces flag
 *   match_vertex <-- matching vertices flag
 *----------------------------------------------------------------------------*/

void
cs_restart_check_base_location(const cs_restart_t  *restart,
                               cs_bool_t           *match_cell,
                               cs_bool_t           *match_i_face,
                               cs_bool_t           *match_b_face,
                               cs_bool_t           *match_vertex);

/*----------------------------------------------------------------------------
 * Add a location definition.
 *
 * parameters:
 *   restart        <-- associated restart file pointer
 *   location_name  <-- name associated with the location
 *   n_glob_ents    <-- global number of entities
 *   n_ents         <-- local number of entities
 *   ent_global_num <-- global entity numbers, or NULL
 *
 * returns:
 *   the location id assigned, or -1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_add_location(cs_restart_t      *restart,
                        const char        *location_name,
                        fvm_gnum_t         n_glob_ents,
                        fvm_lnum_t         n_ents,
                        const fvm_gnum_t  *ent_global_num);

/*----------------------------------------------------------------------------
 * Print the index associated with a restart file in read mode
 *
 * parameters:
 *   restart <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_dump_index(const cs_restart_t  *restart);

/*----------------------------------------------------------------------------
 * Read a section from a restart file.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *   val             --> array of values
 *
 * returns: 0 (CS_RESTART_SUCCES) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_section(cs_restart_t  *restart,
                        const char    *sec_name,
                        int            location_id,
                        cs_int_t       n_location_vals,
                        cs_type_t      val_type,
                        void          *val);

/*----------------------------------------------------------------------------
 * Write a section to a restart file.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *   val             <-- array of values
 *----------------------------------------------------------------------------*/

void
cs_restart_write_section(cs_restart_t  *restart,
                         const char    *sec_name,
                         int            location_id,
                         cs_int_t       n_location_vals,
                         cs_type_t      val_type,
                         const void    *val);

/*----------------------------------------------------------------------------
 * Print statistics associated with restart files
 *----------------------------------------------------------------------------*/

void
cs_restart_print_stats(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_H__ */
