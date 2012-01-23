#ifndef __CS_FIELD_H__
#define __CS_FIELD_H__

/*============================================================================
 * Field management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
#define CS_FIELD_ACCUMULATOR         (1 << 5)

#define CS_FIELD_USER                (1 << 6)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Field handling error types */
/*----------------------------*/

typedef enum {

  CS_FIELD_OK,
  CS_FIELD_INVALID_KEY_NAME,
  CS_FIELD_INVALID_KEY_ID,
  CS_FIELD_INVALID_CATEGORY,
  CS_FIELD_INVALID_TYPE

} cs_field_error_type_t;

/* Field boundary condition descriptor (for variables) */
/*-----------------------------------------------------*/

typedef struct {

  int                location_id;  /* Id of matching location */

  cs_real_t         *a;            /* Explicit coefficient */
  cs_real_t         *b;            /* Implicit coefficient */
  cs_real_t         *af;           /* Explicit coefficient for flux */
  cs_real_t         *bf;           /* Implicit coefficient for flux */

} cs_field_bc_coeffs_t;

/* Field descriptor */
/*------------------*/

typedef struct {

  const char             *name;         /* Canonical name */

  int                     id;           /* Field id */
  int                     type;         /* Field type flag */

  int                     dim;          /* Field dimension */
  bool                    interleaved;  /* is field interleaved ? */

  int                     location_id;  /* Id of matching location */

  int                     n_time_vals;  /* Number of time values (1 or 2) */

  cs_real_t              *val;          /* For each active location, pointer
                                           to matching values array */

  cs_real_t              *val_pre;      /* For each active location, pointer
                                           to matching previous values array
                                           (if n_time_vals == 2) */

  cs_field_bc_coeffs_t   *bc_coeffs;    /* Boundary condition coefficients,
                                           for variable type fields */

  bool                    is_owner;     /* Ownership flag for values
                                           and boundary coefficients */

} cs_field_t;

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a field.
 *
 * Fortran interface; use flddef (see field.f90)
 *
 * subroutine fldde1 (name, lname, iexten, itycat, ityloc, idim, ilved,
 * *****************
 *                    iprev, idfld)
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          iexten      : <-- : 1: intensive; 2: extensive
 * integer          itycat      : <-- : Field category (may be added)
 *                              :     :   4: variable
 *                              :     :   8: property
 *                              :     :  16: postprocess
 *                              :     :  32: accumulator
 *                              :     :  64: user
 * integer          ityloc      : <-- : Location type
 *                              :     :  0: none
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: interior faces
 *                              :     :  4: vertices
 * integer          idim        : <-- : Field dimension
 * integer          ilved       : <-- : 0: not intereaved; 1: interleaved
 * integer          iprev       : <-- : 0: no previous values, 1: previous
 * integer          idfld       : --> : id of defined field
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldde1, FLDDE1)
(
 const char       *name,
 const cs_int_t   *lname,
 const cs_int_t   *iexten,
 const cs_int_t   *itycat,
 const cs_int_t   *ityloc,
 const cs_int_t   *idim,
 const cs_int_t   *ilved,
 const cs_int_t   *iprev,
 cs_int_t         *idfld
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Allocate field values
 *
 * Fortran interface
 *
 * subroutine fldalo (ifield)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldalo, FLDALO)
(
 const cs_int_t   *ifield
);

/*----------------------------------------------------------------------------
 * Map values to a field.
 *
 * Fortran interface
 *
 * subroutine fldmap (ifield, val, valp)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * cs_real_t*       val         : <-- : Pointer to field values array
 * cs_real_t*       valp        : <-- : Pointer to values at previous
 *                              :     : time step if field was defined
 *                              :     : with iprev = 1
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldmap, FLDMAP)
(
 const cs_int_t   *ifield,
 cs_real_t        *val,
 cs_real_t        *valp
);

/*----------------------------------------------------------------------------
 * Map field boundary coefficient arrays.
 *
 * Fortran interface
 *
 * subroutine fldbcm (ifield, icpled, a, b, af, bf)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * cs_real_t*       a           : <-- : explicit BC coefficients array
 * cs_real_t*       b           : <-- : implicit BC coefficients array
 * cs_real_t*       af          : <-- : explicit flux BC coefficients array,
 *                              :     : or a (or NULL)
 * cs_real_t*       bf          : <-- : implicit flux BC coefficients array,
 *                              :     : or a (or NULL)
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldbcm, FLDBCM)
(
 const cs_int_t   *ifield,
 cs_real_t        *a,
 cs_real_t        *b,
 cs_real_t        *af,
 cs_real_t        *bf
);

/*----------------------------------------------------------------------------
 * Allocate arrays for all defined fields based on their location.
 *
 * Location sized must thus be known.
 *
 * Fields that do not own their data should all have been mapped at this
 * stage, and are checked.
 *
 * Fortran interface
 *
 * subroutine fldama
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldama, FLDAMA)
(
 void
);

/*----------------------------------------------------------------------------
 * Retrieve field value pointer for scalar field
 *
 * Fortran interface; use fldpts
 *
 * function fldps1 (ifield, iprev)
 * ***************
 *
 * integer          ifield      : <-- : Field id
 * integer          iprev       : <-- : if 1, pointer to previous values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldps1, FLDPS1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *iprev
);

/*----------------------------------------------------------------------------
 * Retrieve field value pointer for vector field
 *
 * Fortran interface; use fldptv
 *
 * function fldpv1 (ifield, iprev)
 * ***************
 *
 * integer          ifield      : <-- : Field id
 * integer          iprev       : <-- : if 1, pointer to previous values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldpv1, FLDPV1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *iprev
);

/*----------------------------------------------------------------------------
 * Return an id associated with a given field name if present.
 *
 * If the field has not been defined previously, -1 is returned.
 *
 * Fortran interface; use fldfid (see field.f90)
 *
 * subroutine fldfi1 (name,   lname,  ifield)
 * *****************
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          ifield      : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldfi1, FLDFI1)
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
 * Fortran interface; use fldkid (see field.f90)
 *
 * subroutine fldki1 (name,   lname,  ikeyid)
 * *****************
 *
 * character*       name        : <-- : Key name
 * integer          lname       : <-- : Key name length
 * integer          ikey        : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldki1, FLDKI1)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ikey
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

/*----------------------------------------------------------------------------
 * Assign a character string for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * Fortran interface; use fldsks (see field.f90)
 *
 * subroutine fldsk1 (ifield, ikey, str, lstr)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * character*       str         : <-- : Associated string
 * integer          lstr        : <-- : Associated string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldsk1, FLDSK1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 const char       *str,
 const cs_int_t   *lstr
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Return a character string for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * Fortran interface; use fldgk1 (see field.f90)
 *
 * subroutine fldgk1 (ifield, ikey, str, lstr)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * character*       str         : --> : Associated string
 * integer          lstr        : <-- : Associated string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgk1, FLDGK1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 char             *str,
 const cs_int_t   *lstr
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
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
 *   interleaved  <-- indicate if values ar interleaved
 *                    (ignored if number of components < 2)
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
                bool          interleaved,
                bool          has_previous);

/*----------------------------------------------------------------------------
 * Allocate arrays for field values.
 *
 * parameters:
 *   f <-- pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_allocate_values(cs_field_t  *f);

/*----------------------------------------------------------------------------
 * Map existing values to field descriptor.
 *
 * parameters:
 *   f           <-> pointer to field structure
 *   val         <-- pointer to array of values
 *   val_pre     <-- pointer to array of previous values, or NULL
 *----------------------------------------------------------------------------*/

void
cs_field_map_values(cs_field_t   *f,
                    cs_real_t    *val,
                    cs_real_t    *val_pre);

/*----------------------------------------------------------------------------
 * Allocate boundary condition coefficient arrays.
 *
 * For fields on location CS_MESH_LOCATION_CELLS, boundary conditions
 * are located on CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * Boundary condition coefficients are not currently supported for other
 * locations (though support could be added by mapping a boundary->location
 * indirection array in the cs_mesh_location_t structure).
 *
 * For multidimensional fields, arrays are assumed to have the same
 * interleaving behavior as the field, unless components are coupled.
 *
 * For multidimensional fields with coupled components, interleaving
 * is the norm, and implicit coefficients arrays are arrays of block matrices,
 * not vectors, so the number of entris for each boundary face is
 * dim*dim instead of dim.
 *
 * parameters:
 *   f            <-- pointer to field structure
 *   have_flux_bc <-- if true, flux BC coefficients (af and bf) are added
 *----------------------------------------------------------------------------*/

void
cs_field_allocate_bc_coeffs(cs_field_t  *f,
                            bool         have_flux_bc);

/*----------------------------------------------------------------------------
 * Map existing field boundary condition coefficient arrays.
 *
 * For fields on location CS_MESH_LOCATION_CELLS, boundary conditions
 * are located on CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * Boundary condition coefficients are not currently supported for other
 * locations (though support could be added by mapping a boundary->location
 * indirection array in the cs_mesh_location_t structure).
 *
 * For multidimensional fields, arrays are assumed to have the same
 * interleaving behavior as the field, unless components are coupled.
 *
 * For multidimensional fields with coupled components, interleaving
 * is the norm, and implicit coefficients arrays are arrays of block matrices,
 * not vectors, so the number of entris for each boundary face is
 * dim*dim instead of dim.
 *
 * parameters:
 *   f  <-> pointer to field structure
 *   a  <-- explicit BC coefficients array
 *   b  <-- implicit BC coefficients array
 *   af <-- explicit flux BC coefficients array, or NULL
 *   bf <-- implicit flux BC coefficients array, or NULL
 *----------------------------------------------------------------------------*/

void
cs_field_map_bc_coeffs(cs_field_t  *f,
                       cs_real_t   *a,
                       cs_real_t   *b,
                       cs_real_t   *af,
                       cs_real_t   *bf);

/*----------------------------------------------------------------------------
 * Destroy all defined fields.
 *----------------------------------------------------------------------------*/

void
cs_field_destroy_all(void);

/*----------------------------------------------------------------------------
 * Allocate arrays for all defined fields based on their location.
 *
 * Location sized must thus be known.
 *
 * Fields that do not own their data should all have been mapped at this
 * stage, and are checked.
 *----------------------------------------------------------------------------*/

void
cs_field_allocate_or_map_all(void);

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
 * Define a key for an string point value by its name and return an
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
cs_field_define_key_str(const char  *name,
                        const char  *default_value,
                        int          type_flag);

/*----------------------------------------------------------------------------
 * Define a sub key.
 *
 * The sub key is the same type as the parent key.
 *
 * For a given field, when querying a sub key's value and that value has not
 * been set, the query will return the value of the parent key.
 *
 * parameters:
 *   name      <-- key name
 *   parent_id <-- parent key id
 *
 * returns:
 *   id associated with key
 *----------------------------------------------------------------------------*/

int
cs_field_define_sub_key(const char  *name,
                        int          parent_id);

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
 *   key_id <-- id of associated key
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
 *   f      <-- pointer to field structure
 *   key_id <-- id of associated key
 *   value  <-- value associated with key
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
 *   f      <-- pointer to field structure
 *   key_id <-- id of associated key
 *   value  <-- value associated with key
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
 *
 * returns:
 *   floating point value associated with the key id for this field
 *----------------------------------------------------------------------------*/

double
cs_field_get_key_double(const cs_field_t  *f,
                        int                key_id);

/*----------------------------------------------------------------------------
 * Assign a character string value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * parameters:
 *   f      <-- pointer to field structure
 *   key_id <-- id of associated key
 *   str    <-- string associated with key
 *
 * returns:
 *   0 in case of success, > 1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_field_set_key_str(cs_field_t  *f,
                     int          key_id,
                     const char  *str);

/*----------------------------------------------------------------------------
 * Return a string for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * parameters:
 *   f             <-- pointer to field structure
 *   key_id        <-- id of associated key
 *
 * returns:
 *   pointer to character string associated with the key id for this field
 *----------------------------------------------------------------------------*/

const char *
cs_field_get_key_str(const cs_field_t  *f,
                     int                key_id);

/*----------------------------------------------------------------------------
 * Print info relative to all field definitions to log file.
 *----------------------------------------------------------------------------*/

void
cs_field_log_defs(void);

/*----------------------------------------------------------------------------
 * Print info relative to a given field to log file.
 *
 * parameters:
 *   f            <-- pointer to field structure
 *   log_keywords <-- log level for keywords (0: do not log,
 *                    1: log non-default values, 2: log all)
 *----------------------------------------------------------------------------*/

void
cs_field_log_info(const cs_field_t  *f,
                  int                log_keywords);

/*----------------------------------------------------------------------------
 * Print info relative to all defined fields to log file.
 *
 * parameters:
 *   log_keywords <-- log level for keywords (0: do not log,
 *                    1: log non-default values, 2: log all)
 *----------------------------------------------------------------------------*/

void
cs_field_log_fields(int  log_keywords);

/*----------------------------------------------------------------------------
 * Print info relative to all key definitions to log file.
 *----------------------------------------------------------------------------*/

void
cs_field_log_key_defs(void);

/*----------------------------------------------------------------------------
 * Print info relative to a given field key to log file.
 *
 * parameters:
 *   int key_id   <-- id of associated key
 *   log_defaults <-- if true, log default field values in addition to
 *                    defined field values
 *----------------------------------------------------------------------------*/

void
cs_field_log_key_vals(int   key_id,
                      bool  log_defaults);

/*----------------------------------------------------------------------------
 * Print info relative to all given field keys to log file.
 *
 * parameters:
 *   log_defaults <-- if true, log default field values in addition to
 *                    defined field values
 *----------------------------------------------------------------------------*/

void
cs_field_log_all_key_vals(bool  log_defaults);

/*----------------------------------------------------------------------------
 * Define base keys.
 *
 * Keys defined by this function are:
 *   "label"     (string)
 *   "post_vis"  (integer)
 *   "coupled"   (integer, restricted to CS_FIELD_VARIABLE)
 *   "moment_dt" (integer, restricted to CS_FIELD_PROPERTY);
 *
 * A recommened practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 *----------------------------------------------------------------------------*/

void
cs_field_define_keys_base(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_H__ */
