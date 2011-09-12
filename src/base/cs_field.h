#ifndef __CS_FIELD_H__
#define __CS_FIELD_H__

/*============================================================================
 * Field management.
 *============================================================================*/

/*
  This file is part of the Code_Saturne Kernel, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1998-2011 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Kernel is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Kernel is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Kernel; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Field property type
 */

#define CS_FIELD_INTENSIVE           (1 << 0)
#define CS_FIELD_EXTENSIVE           (1 << 1)

/* Field category */

#define CS_FIELD_VARIABLE            (1 << 2)
#define CS_FIELD_PROPERTY            (1 << 3)
#define CS_FIELD_POSTPROCESS         (1 << 4)

#define CS_FIELD_USER                (1 << 5)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Field handling error types */

typedef enum {

  CS_FIELD_OK,
  CS_FIELD_INVALID_KEY_NAME,
  CS_FIELD_INVALID_KEY_ID,
  CS_FIELD_INVALID_CATEGORY,
  CS_FIELD_INVALID_TYPE

} cs_field_error_type_t;

/* Field location types */

typedef enum {

  CS_FIELD_CELLS,
  CS_FIELD_INTERIOR_FACES,
  CS_FIELD_BOUNDARY_FACES,
  CS_FIELD_VERTICES,
  CS_FIELD_OTHER
  // CS_FIELD_PARTICLES

} cs_field_location_type_t;

/* Field descriptor */

typedef struct {

  const char        *name;         /* Canonical name */

  char               label[64];    /* Field label */

  int                id;           /* Field id */
  int                type;         /* Field type flag */

  int                dim;          /* Field dimension */
  bool               interleaved;  /* is field interleaved ? */

  int                location_id;  /* Id of matching location */

  cs_real_t         *val;          /* For each active location, pointer
                                      to matching values array */

  cs_real_t         *val_pre;      /* For each active location, pointer
                                      to matching previous values array
                                      (for CS_FIELD_VARIABLE and
                                      CS_FIELD_PROPERTY categories only) */

  bool               is_owner;     /* Ownership flag for values */

} cs_field_t;

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Map an existing field.
 *
 * Fortran interface
 *
 * subroutine fldmap (name,   lname,  ireawr, numsui, ierror)
 * *****************
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          iexten      : <-- : 0: intensive; 1: extensive
 * integer          itycat      : <-- : Field category
 *                              :     :  0: variable
 *                              :     :  1: property
 * integer          ityloc      : <-- : Location type
 *                              :     :  0: cells
 *                              :     :  1: interior faces
 *                              :     :  2: interior faces
 *                              :     :  3: vertices
 * integer          idim        : <-- : Field dimension
 * integer          ilved       : <-- : 0: not intereaved; 1: interleaved
 * integer          iprev       : <-- : 0: no previous values, 1: previous
 * cs_real_t*       val         : <-- : Pointer to field values array
 * cs_real_t*       valp        : <-- : Pointer to values at previous
 *                              :     : time step if iprev = 1
 * integer          idfld       : --> : id of mapped field
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldmap, FLDMAP)
(
 const char       *name,
 const cs_int_t   *lname,
 const cs_int_t   *iexten,
 const cs_int_t   *itycat,
 const cs_int_t   *ityloc,
 const cs_int_t   *idim,
 const cs_int_t   *ilved,
 const cs_int_t   *iprev,
 cs_real_t        *val,
 cs_real_t        *valp,
 cs_int_t         *idfld
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Return an id associated with a given field name if present.
 *
 * If the field has not been defined previously, -1 is returned.
 *
 * Fortran interface
 *
 * subroutine fldfid (name,   lname,  ifield)
 * *****************
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          ifield      : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldfid, FLDFID)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ifield
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Return an id associated with a given key name if present.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * Fortran interface
 *
 * subroutine fldkid (name,   lname,  ikeyid)
 * *****************
 *
 * character*       name        : <-- : Key name
 * integer          lname       : <-- : Key name length
 * integer          ikeyid      : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldkid, FLDKID)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ikeyid
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Assign an integer value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldski (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * integer          value       : <-- : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldski, FLDSKI)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_int_t         *value
);

/*----------------------------------------------------------------------------
 * Return a integer value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldgki (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * integer          value       : --> : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgki, FLDGKI)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_int_t         *value
);

/*----------------------------------------------------------------------------
 * Assign a floating point value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldskd (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * double precision value       : <-- : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldskd, FLDSKD)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_real_t        *value
);

/*----------------------------------------------------------------------------
 * Return a floating point value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldgkd (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * double precision value       : --> : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgkd, FLDGKD)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_real_t        *value
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the number of defined fields.
 *
 * returns:
 *   number of defined fields.
 *----------------------------------------------------------------------------*/

int
cs_field_n_fields(void);

/*----------------------------------------------------------------------------
 * Create a field descriptor.
 *
 * For fields with a dimension greater than 1, components are interleaved.
 *
 * parameters:
 *   name         <-- field name
 *   type_flag    <-- mask of field property and category values
 *   location_id  <-- id of associated location
 *   dim          <-- field dimension (number of components)
 *   has_previous <-- maintain values at the previous time step ?
 *
 * returns:
 *   pointer to new field.
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_field_create(const char   *name,
                int           type_flag,
                int           location_id,
                int           dim,
                bool          has_previous);

/*----------------------------------------------------------------------------
 * Create a field descriptor for an existing array of values.
 *
 * parameters:
 *   name        <-- field name
 *   type_flag   <-- mask of field property and category values
 *   location_id <-- id of associated location
 *   dim         <-- field dimension (number of components)
 *   interleaved <-- true if field components are interleaved
 *   val         <-- pointer to array of values
 *   val_pre     <-- pointer to array of previous values, or NULL
 *
 * returns:
 *   pointer to new field.
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_field_create_mapped(const char   *name,
                       int           type_flag,
                       int           location_id,
                       int           dim,
                       bool          interleaved,
                       cs_real_t    *val,
                       cs_real_t    *val_pre);

/*----------------------------------------------------------------------------
 * Destroy all defined fields.
 *----------------------------------------------------------------------------*/

void
cs_field_destroy_all(void);

/*----------------------------------------------------------------------------
 * Return a pointer to a field based on its id.
 *
 * This function requires that a field of the given id is defined.
 *
 * parameters:
 *   id <-- field id
 *
 * returns:
 *   pointer to the field structure
 *----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_id(int  id);

/*----------------------------------------------------------------------------
 * Return a pointer to a field based on its name.
 *
 * This function requires that a field of the given name is defined.
 *
 * parameters:
 *   name <-- field name
 *
 * returns:
 *   pointer to the field structure
 *----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_name(const char *name);

/*----------------------------------------------------------------------------
 * Return a pointer to a field based on its name if present.
 *
 * If no field of the given name is defined, NULL is returned.
 *
 * parameters:
 *   name <-- field name
 *
 * returns:
 *   pointer to the field structure, or NULL
 *----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_name_try(const char *name);

/*----------------------------------------------------------------------------
 * Return an id associated with a given key name.
 *
 * The key must have been defined previously.
 *
 * parameters:
 *   name <-- key name
 *
 * returns:
 *   id associated with key
 *----------------------------------------------------------------------------*/

int
cs_field_key_id(const char  *name);

/*----------------------------------------------------------------------------
 * Return an id associated with a given key name if present.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * parameters:
 *   name <-- key name
 *
 * returns:
 *   id associated with key, or -1
 *----------------------------------------------------------------------------*/

int
cs_field_key_id_try(const char  *name);

/*----------------------------------------------------------------------------
 * Define a key for an integer value by its name and return an associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * parameters:
 *   name          <-- key name
 *   default_value <-- default value associated with key
 *   type flag     <-- mask associated with field types with which the
 *                     key may be associated, or 0
 *
 * returns:
 *   id associated with key
 *----------------------------------------------------------------------------*/

int
cs_field_define_key_int(const char  *name,
                        int          default_value,
                        int          type_flag);

/*----------------------------------------------------------------------------
 * Define a key for an floating point value by its name and return an
 * associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * parameters:
 *   name          <-- key name
 *   default_value <-- default value associated with key
 *   type flag     <-- mask associated with field types with which the
 *                     key may be associated, or 0
 *
 * returns:
 *   id associated with key
 *----------------------------------------------------------------------------*/

int
cs_field_define_key_double(const char  *name,
                           double       default_value,
                           int          type_flag);

/*----------------------------------------------------------------------------
 * Destroy all defined field keys and associated values.
 *----------------------------------------------------------------------------*/

void
cs_field_destroy_all_keys(void);

/*----------------------------------------------------------------------------
 * Get the type flag associated with a given key id.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * parameters:
 *   key_id  <-- id of associated key
 *
 * returns:
 *   type flag associated with key, or -1
 *----------------------------------------------------------------------------*/

int
cs_field_key_flag(int key_id);

/*----------------------------------------------------------------------------
 * Assign a integer value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * parameters:
 *   f             <-- pointer to field structure
 *   key_id        <-- id of associated key
 *   value         <-- value associated with key
 *
 * returns:
 *   0 in case of success, > 1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_field_set_key_int(cs_field_t  *f,
                     int          key_id,
                     int          value);

/*----------------------------------------------------------------------------
 * Return a integer value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * parameters:
 *   f             <-- pointer to field structure
 *   key_id        <-- id of associated key
 *   value         <-- value associated with key
 *
 * returns:
 *   integer value associated with the key id for this field
 *----------------------------------------------------------------------------*/

int
cs_field_get_key_int(const cs_field_t  *f,
                     int                key_id);

/*----------------------------------------------------------------------------
 * Assign a floating point value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * parameters:
 *   f             <-- pointer to field structure
 *   key_id        <-- id of associated key
 *   value         <-- value associated with key
 *
 * returns:
 *   0 in case of success, > 1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_field_set_key_double(cs_field_t  *f,
                        int          key_id,
                        double       value);

/*----------------------------------------------------------------------------
 * Return a floating point value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * parameters:
 *   f             <-- pointer to field structure
 *   key_id        <-- id of associated key
 *   value         <-- value associated with key
 *
 * returns:
 *   floating point value associated with the key id for this field
 *----------------------------------------------------------------------------*/

double
cs_field_get_key_double(const cs_field_t  *f,
                        int                key_id);

/*----------------------------------------------------------------------------
 * Print info relative to a given field to log file.
 *
 * parameters:
 *   f <-- pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_log_info(const cs_field_t  *f);

/*----------------------------------------------------------------------------
 * Print info relative to all defined fields to log file.
 *----------------------------------------------------------------------------*/

void
cs_field_log_fields(void);

/*----------------------------------------------------------------------------
 * Print info relative to all defined field keys to log file.
 *----------------------------------------------------------------------------*/

void
cs_field_log_keys(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_H__ */
