#ifndef __CS_RESTART_H__
#define __CS_RESTART_H__

/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Error codes */

#define CS_RESTART_SUCCESS        0 /* Success */
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

/* Datatype enumeration to transmit a data's type to a function */

typedef enum {
  CS_TYPE_char,
  CS_TYPE_cs_int_t,
  CS_TYPE_cs_gnum_t,
  CS_TYPE_cs_real_t,
} cs_restart_val_type_t;

/*
  Pointer associated with a restart file structure. The structure itself
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
 * subroutine dflsui (ntsuit, ttsuit, wtsuit)
 * *****************
 *
 * integer          ntsuit      : <-- : > 0: checkpoint time step interval
 *                              :     : 0: default interval
 *                              :     : -1: checkpoint at end
 *                              :     : -2: no checkpoint
 * double precision ttsuit      : <-- : if> 0, checkpoint time interval
 * double precision wtsuit      : <-- : if> 0, checkpoint wall time interval
 *----------------------------------------------------------------------------*/

void CS_PROCF (dflsui, DFLSUI)
(
 cs_int_t   *ntsuit,
 cs_real_t  *ttsuit,
 cs_real_t  *wtsuit
);

/*----------------------------------------------------------------------------
 * Check if checkpointing is recommended at a given time.
 *
 * Fortran interface
 *
 * subroutine reqsui (iisuit)
 * *****************
 *
 * integer          iisuit      : --> : 0 if no restart required, 1 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF (reqsui, RESSUI)
(
 cs_int_t   *iisuit
);

/*----------------------------------------------------------------------------
 * Indicate checkpointing has been done at a given time.
 *
 * This updates the status for future checks to determine
 * if checkpointing is recommended at a given time.
 *
 * Fortran interface
 *
 * subroutine stusui
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (stusui, STUSUI)
(
 void
);

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
 * subroutine opnsui (nomsui, lngnom, ireawr, numsui, ierror)
 * *****************
 *
 * character*       nomsui      : <-- : Restart file name
 * integer          lngnom      : <-- : Restart file name length
 * integer          ireawr      : <-- : 1: read; 2: write
 * integer          numsui      : --> : Number of opened restart file
 * integer          ierror      : --> : 0: success; < 0: error code
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
 * subroutine clssui (numsui)
 * *****************
 *
 * integer          numsui      : <-> : Number of restart file to close
 * integer          ierror      : --> : 0: success; < 0: error code
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
 * subroutine tstsui (numsui, indcel, indfac, indfbr, indsom)
 * *****************
 *
 * integer          numsui      : <-- : Restart file number
 * integer          indcel      : --> : Matching cells flag
 * integer          indfac      : --> : Matching interior faces flag
 * integer          indfbr      : --> : Matching boundary faces flag
 * integer          indsom      : --> : Matching vertices flag
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
 * INTEGER          NUMSUI      : <-- : Restart file number
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
 * subroutine lecsui (numsui, nomrub, lngnom, itysup, nbvent, irtype, tabvar)
 * *****************
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Section name
 * integer          lngnom      : <-- : Section name length
 * integer          itysup      : <-- : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          nbvent      : <-- : N. values per location entity
 * integer          irtype      : <-- : 1 for integers, 2 for double precision
 * (?)              tabvar      : <-> : Array of values to read
 * integer          ierror      : --> : 0: success, < 0: error code
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
 * subroutine ecrsui (numsui, nomrub, lngnom, itysup, nbvent, irtype, tabvar)
 * *****************
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Section name
 * integer          lngnom      : <-- : Section name length
 * integer          itysup      : <-- : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          nbvent      : <-- : N. values per location entity
 * integer          irtype      : <-- : 1 for integers, 2 for double precision
 * (?)              tabvar      : <-- : Array of values to write
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrsui, ECRSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *nbvent,
 const cs_int_t   *irtype,
 const void       *tabvar
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Read basic particles information from a restart file.
 *
 * Fortran interface
 *
 * subroutine lipsui (numsui, nomrub, lngnom, itysup, nbvent, irtype, tabvar)
 * *****************
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Particles location name
 * integer          lngnom      : <-- : Particles location name length
 * integer          nbpart      : --> : Number of particles
 * integer          itysup      : --> : Particles location id,
 *                                      or -1 in case of error
 *----------------------------------------------------------------------------*/

void CS_PROCF (lipsui, LIPSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
       cs_int_t   *nbpart,
       cs_int_t   *itysup
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Read basic particles information from a restart file.
 *
 * Fortran interface
 *
 * subroutine lepsui (numsui, nomrub, lngnom, inmcoo, nbpart, ipcell,
 * *****************
 *                    coopar, itysup, ierror)
 *
 * integer          numsui      : <-- : Restart file number
 * integer          ipcell      : --> : Particle -> cell number
 * double precision coopar      : --> : Particle coordinate
 * integer          ipsup       : <-- : Particles location id
 * integer          ierror      : --> : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (lepsui, LEPSUI)
(
 const cs_int_t   *numsui,
       cs_int_t   *ipcell,
       cs_real_t  *coopar,
 const cs_int_t   *itysup,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Write basic particles information to a restart file.
 *
 * This includes defining a matching location and associated global numbering,
 * then writing particle coordinates and cell ids.
 *
 * Fortran interface
 *
 * subroutine ecpsui (numsui, nomrub, lngnom, inmcoo, nbpart, ipcell,
 * *****************
 *                    coopar, itysup, ierror)
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Particles location name
 * integer          lngnom      : <-- : Particles location name length
 * integer          inmcoo      : <-- : Number by coords
 * integer          nbpart      : <-- : Number of particles
 * integer          ipcell      : <-- : Particle -> cell number
 * double precision coopar      : <-- : Particle coordinates
 * integer          ipsup       : --> : Particles location id
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecpsui, ECPSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *inmcoo,
 const cs_int_t   *nbpart,
 const cs_int_t   *ipcell,
 const cs_real_t  *coopar,
       cs_int_t   *itysup
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Read a referenced location id section from a restart file.
 *
 * The section read from file contains the global ids matching the local
 * element ids of a given location. Global id's are transformed to local
 * ids by this function.
 *
 * In case global ids read do not match those of local elements,
 * id_base - 1 is assigned to the corresponding local ids.
 *
 * Fortran interface
 *
 * subroutine leisui (numsui, nomrub, lngnom, itysup, irfsup, idbase, tabid, &
 * *****************
 *                    ierror)
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Section name
 * integer          lngnom      : <-- : Section name length
 * integer          itysup      : <-- : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          irfsup      : <-- : Referenced location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          idbase      : <-- : base of referenced entity id numbers
 *                              :     : (usually 0 or 1)
 * (?)              tabid       : <-> : Array of ids to read
 * integer          ierror      : --> : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (leisui, LEISUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *irfsup,
 const cs_int_t   *idbase,
       void       *tabid,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
 );

/*----------------------------------------------------------------------------
 * Write a referenced location id section to a restart file.
 *
 * The section written to file contains the global ids matching the
 * local element ids of a given location.
 *
 * Fortran interface
 *
 * subroutine ecisui (numsui, nomrub, lngnom, itysup, irfsup, idbase, tabid, &
 * *****************
 *                    ierror)
 *
 * integer          numsui      : <-- : Restart file number
 * character*       nomrub      : <-- : Section name
 * integer          lngnom      : <-- : Section name length
 * integer          itysup      : <-- : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          irfsup      : <-- : Referenced location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * integer          idbase      : <-- : base of referenced entity id numbers
 *                              :     : (usually 0 or 1)
 * (?)              tabid       : <-- : Array of ids to write
 * integer          ierror      : --> : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecisui, ECISUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *irfsup,
 const cs_int_t   *idbase,
 const cs_int_t   *tabid
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define default checkpoint interval.
 *
 * parameters
 *   nt_interval <-- if > 0 time step interval for checkpoint
 *                   if 0, default of 4 checkpoints per run
 *                   if -1, checkpoint at end
 *                   if -2, no checkpointing
 *   t_interval  <-- if > 0, time value interval for checkpoint
 *   wt_interval <-- if > 0, wall-clock interval for checkpoints
 *----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_defaults(int     nt_interval,
                                   double  t_interval,
                                   double  wt_interval);

/*----------------------------------------------------------------------------
 * Define next forced checkpoint time step
 *
 * parameters
 *   nt_next <-- next time step for forced checkpoint
 *----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_ts(int  nt_next);

/*----------------------------------------------------------------------------
 * Define next forced checkpoint time value
 *
 * parameters
 *   t_next <-- next time value for forced checkpoint
 *----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_tv(double  t_next);

/*----------------------------------------------------------------------------
 * Define next forced checkpoint wall-clock time value
 *
 * parameters
 *   wt_next <-- next wall-clock time value for forced checkpoint
 *----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_wt(double  wt_next);

/*----------------------------------------------------------------------------
 * Check if checkpointing is recommended at a given time.
 *
 * parameters
 *   ts <-- time step status structure
 *
 * returns:
 *   true if checkpointing is recommended, 0 otherwise
 *----------------------------------------------------------------------------*/

bool
cs_restart_checkpoint_required(const cs_time_step_t  *ts);

/*----------------------------------------------------------------------------
 * Indicate checkpointing has been done at a given time.
 *
 * This updates the status for future checks to determine
 * if checkpointing is recommended at a given time.
 *
 * parameters
 *   ts <-- time step status structure
 *----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_done(const cs_time_step_t  *ts);

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
 *   path <-- optional directory name for output
 *            (directory automatically created if necessary)
 *   mode <-- read or write
 *
 * returns:
 *   pointer to initialized restart file structure
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_create(const char         *name,
                  const char         *path,
                  cs_restart_mode_t   mode);

/*----------------------------------------------------------------------------
 * Destroy structure associated with a restart file (and close the file).
 *
 * parameters:
 *   restart <-- pointer to restart file structure pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_destroy(cs_restart_t  **restart);

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
                               bool                *match_cell,
                               bool                *match_i_face,
                               bool                *match_b_face,
                               bool                *match_vertex);

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
                        cs_gnum_t          n_glob_ents,
                        cs_lnum_t          n_ents,
                        const cs_gnum_t   *ent_global_num);

/*----------------------------------------------------------------------------
 * Print the index associated with a restart file in read mode
 *
 * parameters:
 *   restart <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_restart_dump_index(const cs_restart_t  *restart);

/*----------------------------------------------------------------------------
 * Check the presence of a given section in a restart file.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_check_section(cs_restart_t           *restart,
                         const char             *sec_name,
                         int                     location_id,
                         int                     n_location_vals,
                         cs_restart_val_type_t   val_type);

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
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_section(cs_restart_t           *restart,
                        const char             *sec_name,
                        int                     location_id,
                        int                     n_location_vals,
                        cs_restart_val_type_t   val_type,
                        void                   *val);

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
cs_restart_write_section(cs_restart_t           *restart,
                         const char             *sec_name,
                         int                     location_id,
                         int                     n_location_vals,
                         cs_restart_val_type_t   val_type,
                         const void             *val);

/*----------------------------------------------------------------------------
 * Read basic particles information from a restart file.
 *
 * This includes building a matching location and associated global numbering.
 *
 * parameters:
 *   restart      <-- associated restart file pointer
 *   name         <-- name of particles set
 *   n_particles  --> number of particles, or NULL
 *
 * returns:
 *   the location id assigned to the particles, or -1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_particles_info(cs_restart_t  *restart,
                               const char    *name,
                               cs_lnum_t     *n_particles);

/*----------------------------------------------------------------------------
 * Read basic particles information from a restart file.
 *
 * parameters:
 *   restart               <-- associated restart file pointer
 *   particles_location_id <-- location id of particles set
 *   particle_cell_id      --> local cell id to which particles belong
 *   particle_coords       --> local particle coordinates (interleaved)
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_particles(cs_restart_t  *restart,
                          int            particles_location_id,
                          cs_lnum_t     *particle_cell_id,
                          cs_real_t     *particle_coords);

/*----------------------------------------------------------------------------
 * Write basic particles information to a restart file.
 *
 * This includes defining a matching location and associated global numbering,
 * then writing particle coordinates and cell ids.
 *
 * parameters:
 *   restart           <-- associated restart file pointer
 *   name              <-- name of particles set
 *   number_by_coords  <-- if true, numbering is based on current coordinates;
 *                         otherwise, it is simply based on local numbers,
 *                         plus the sum of particles on lower MPI ranks
 *   n_particles       <-- local number of particles
 *   particle_cell_num <-- local cell number (1 to n) to which particles
 *                         belong; 0 for untracked particles
 *   particle_coords   <-- local particle coordinates (interleaved)
 *
 * returns:
 *   the location id assigned to the particles
 *----------------------------------------------------------------------------*/

int
cs_restart_write_particles(cs_restart_t     *restart,
                           const char       *name,
                           bool              number_by_coords,
                           cs_lnum_t         n_particles,
                           const cs_lnum_t  *particle_cell_num,
                           const cs_real_t  *particle_coords);

/*----------------------------------------------------------------------------
 * Read a referenced location id section from a restart file.
 *
 * The section read from file contains the global ids matching the local
 * element ids of a given location. Global id's are transformed to local
 * ids by this function.
 *
 * In case global referenced ids read do not match those of local elements,
 * id_base - 1 is assigned to the corresponding local ids.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of location on which id_ref is defined
 *   ref_location_id <-- id of referenced location
 *   ref_id_base     <-- base of location entity id numbers (usually 0 or 1)
 *   ref_id          --> array of location entity ids
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_ids(cs_restart_t     *restart,
                    const char       *sec_name,
                    int               location_id,
                    int               ref_location_id,
                    cs_lnum_t         ref_id_base,
                    cs_lnum_t        *ref_id);

/*----------------------------------------------------------------------------
 * Write a referenced location id section to a restart file.
 *
 * The section written to file contains the global ids matching the
 * local element ids of a given location.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of location on which id_ref is defined
 *   ref_location_id <-- id of referenced location
 *   ref_id_base     <-- base of location entity id numbers (usually 0 or 1)
 *   ref_id          <-- array of location entity ids
 *----------------------------------------------------------------------------*/

void
cs_restart_write_ids(cs_restart_t           *restart,
                     const char             *sec_name,
                     int                     location_id,
                     int                     ref_location_id,
                     cs_lnum_t               ref_id_base,
                     const cs_lnum_t        *ref_id);

/*----------------------------------------------------------------------------
 * Read a section from a restart file, when that section may have used a
 * different name in a previous version.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *   val             --> array of values
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_section_compat(cs_restart_t           *restart,
                               const char             *sec_name,
                               const char             *old_name,
                               int                     location_id,
                               int                     n_location_vals,
                               cs_restart_val_type_t   val_type,
                               void                   *val);

/*----------------------------------------------------------------------------
 * Read a cs_real_t section from a restart file, when that section may
 * have used a different name in a previous version.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val             --> array of values
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_real_t_compat(cs_restart_t  *restart,
                              const char    *sec_name,
                              const char    *old_name,
                              int            location_id,
                              int            n_location_vals,
                              cs_real_t     *val);

/*----------------------------------------------------------------------------
 * Read a cs_real_3_t vector section from a restart file, when that
 * section may have used a different name and been non-interleaved
 * in a previous version.
 *
 * This file assumes a mesh-base location (i.e. location_id > 0)
 *
 * parameters:
 *   restart     <-- associated restart file pointer
 *   sec_name    <-- section name
 *   old_name_x  <-- old name, x component
 *   old_name_y  <-- old name, y component
 *   old_name_y  <-- old name, z component
 *   location_id <-- id of corresponding location
 *   val         --> array of values
 *
 * returns: 0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_real_3_t_compat(cs_restart_t  *restart,
                                const char    *sec_name,
                                const char    *old_name_x,
                                const char    *old_name_y,
                                const char    *old_name_z,
                                int            location_id,
                                cs_real_3_t   *val);

/*----------------------------------------------------------------------------
 * Print statistics associated with restart files
 *----------------------------------------------------------------------------*/

void
cs_restart_print_stats(void);

/*----------------------------------------------------------------------------
 * Return pointer to restart file based on Fortran id
 *
 * parameters:
 *   r_num  <-- associated fortran restart number
 *
 * returns:
 *   pointer to restart file, or NULL
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_by_fortran_id(int  r_num);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RESTART_H__ */
