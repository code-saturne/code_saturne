#ifndef __CS_INTERFACE_H__
#define __CS_INTERFACE_H__

/*============================================================================
 * Main structure for handling of interfaces associating mesh elements
 * (such as inter-processor or periodic connectivity between cells, faces,
 * or vertices);
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "fvm_defs.h"
#include "fvm_periodicity.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining an I/O numbering scheme
 *----------------------------------------------------------------------------*/

/*
  Pointer to structures representing an interface and a list of interfaces.
  The structures themselves are private, and are defined in cs_interface.c
*/

typedef struct _cs_interface_t     cs_interface_t;
typedef struct _cs_interface_set_t cs_interface_set_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return process rank associated with an interface's distant elements.
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   process rank associated with the interface's distant elements
 *----------------------------------------------------------------------------*/

int
cs_interface_rank(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Return number of local and distant elements defining an interface.
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   number of local and distant elements defining the interface
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_size(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Return pointer to array of local element ids defining an interface.
 *
 * The size of the array may be obtained by cs_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local element ids (0 to n-1) defining the interface
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_elt_ids(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Return pointer to array of matching element ids defining an interface.
 *
 * This array is only available if cs_interface_set_add_match_ids() has
 * been called for the containing interface set.
 *
 * The size of the array may be obtained by cs_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   pointer to array of local element ids (0 to n-1) defining the interface
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_match_ids(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Return size of index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * Its size should thus be equal to 1 + number of periodic transforms + 1,
 * In absence of periodicity, it may be 0, as the index is not needed.
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   transform index size for the interface
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_get_tr_index_size(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Return pointer to index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * In absence of periodicity, it may be NULL, as it is not needed.
 *
 * parameters:
 *   itf <-- pointer to interface structure
 *
 * returns:
 *   pointer to transform index for the interface
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_tr_index(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------
 * Creation of a list of interfaces between elements of a same type.
 *
 * These interfaces may be used to identify equivalent vertices or faces using
 * domain splitting, as well as periodic elements (on the same or on
 * distant ranks).
 *
 * Note that periodicity information will be completed and made consistent
 * based on the input, so that if a periodic couple is defined on a given rank,
 * the reverse couple wil be defined, whether it is also defined on the same
 * or a different rank.
 *
 * In addition, multiple periodicity interfaces will be built automatically
 * if the periodicity structure provides for composed periodicities, so they
 * need not be defined prior to this function.
 *
 * parameters:
 *   n_elts             <-- number of local elements considered
 *                          (size of parent_element_id[]
 *   parent_element_id  <-- pointer to list of selected elements local
 *                          numbers (0 to n-1), or NULL if all first n_elts
 *                          elements are used
 *   global_number      <-- pointer to list of global (i.e. domain splitting
 *                          independent) element numbers
 *   periodicity        <-- periodicity information (NULL if none)
 *   n_periodic_lists   <-- number of periodic lists (may be local)
 *   periodicity_num    <-- periodicity number (1 to n) associated with
 *                          each periodic list (primary periodicities only)
 *   n_periodic_couples <-- number of periodic couples associated with
 *                          each periodic list
 *   periodic_couples   <-- array indicating periodic couples (using
 *                          global numberings) for each list
 *
 * returns:
 *  pointer to list of interfaces (possibly NULL in serial mode)
 *----------------------------------------------------------------------------*/

cs_interface_set_t *
cs_interface_set_create(cs_lnum_t                 n_elts,
                        const cs_lnum_t           parent_element_id[],
                        const cs_gnum_t           global_number[],
                        const fvm_periodicity_t  *periodicity,
                        int                       n_periodic_lists,
                        const int                 periodicity_num[],
                        const cs_lnum_t           n_periodic_couples[],
                        const cs_gnum_t    *const periodic_couples[]);

/*----------------------------------------------------------------------------
 * Destruction of an interface set.
 *
 * parameters:
 *   ifs <-> pointer to pointer to structure to destroy
 *----------------------------------------------------------------------------*/

void
cs_interface_set_destroy(cs_interface_set_t  **ifs);

/*----------------------------------------------------------------------------
 * Return number of interfaces associated with an interface set.
 *
 * parameters:
 *   ifs <-- pointer to interface set structure
 *
 * returns:
 *   number of interfaces in set
 *----------------------------------------------------------------------------*/

int
cs_interface_set_size(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Return total number of elements in interface set.
 *
 * This is equal to the sum of cs_interface_size() on the cs_interface_size()
 * interfaces of a set.
 *
 * parameters:
 *   ifs <-- pointer to interface set structure
 *
 * returns:
 *   number of interfaces in set
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_set_n_elts(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Return pointer to a given interface in an interface set.
 *
 * parameters:
 *   ifs          <-- pointer to interface set structure
 *   interface_id <-- index of interface in set (0 to n-1)
 *
 * returns:
 *   pointer to interface structure
 *----------------------------------------------------------------------------*/

const cs_interface_t *
cs_interface_set_get(const cs_interface_set_t  *ifs,
                     int                        interface_id);

/*----------------------------------------------------------------------------
 * Return pointer to the periocicity structure associated of an interface set.
 *
 * parameters:
 *   ifs <-- pointer to interface set structure
 *
 * returns:
 *   pointer to periodicity structure, or NULL
 *----------------------------------------------------------------------------*/

const fvm_periodicity_t *
cs_interface_set_periodicity(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Apply renumbering of elements referenced by an interface set.
 *
 * For any given element i, a negative old_to_new[i] value means that that
 * element does not appear anymore in the new numbering.
 *
 * parameters:
 *   ifs        <-> pointer to interface set structure
 *   old_to_new <-- renumbering array (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

void
cs_interface_set_renumber(cs_interface_set_t  *ifs,
                          const cs_lnum_t      old_to_new[]);

/*----------------------------------------------------------------------------
 * Add matching element id information to an interface set.
 *
 * This information is required by calls to cs_interface_get_dist_ids(),
 * and may be freed using cs_interface_set_free_match_ids().
 *
 * parameters:
 *   ifs <-> pointer to interface set structure
 *----------------------------------------------------------------------------*/

void
cs_interface_set_add_match_ids(cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Free matching element id information of an interface set.
 *
 * This information is used by calls to cs_interface_get_dist_ids(),
 * and may be defined using cs_interface_set_add_match_ids().
 *
 * parameters:
 *   ifs <-> pointer to interface set structure
 *----------------------------------------------------------------------------*/

void
cs_interface_set_free_match_ids(cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------
 * Copy array from distant or matching interface elements to local elements.
 *
 * Source and destination arrays define values for all elements in the
 * interface set (i.e. all elements listed by cs_interface_get_elt_ids()
 * when looping over interfaces of a set,
 *
 * parameters:
 *   ifs           <-- pointer to interface set structure
 *   datatype      <-- type of data considered
 *   stride        <-- number of values per entity (interlaced)
 *   src_on_parent <-- true if source array is defined on the elements
 *                     defined by ifs->elt_ids, false if source array
 *                     defined directly on cs_interface_set_n_elts(ifs)
 *   src           <-- source array (size: cs_interface_set_n_elts(ifs)*stride
 *                     or parent array size * stride)
 *   dest          <-- destination array
 *                     (size: cs_interface_set_n_elts(ifs)*stride)
 *----------------------------------------------------------------------------*/

void
cs_interface_set_copy_array(const cs_interface_set_t  *ifs,
                            cs_datatype_t              datatype,
                            int                        stride,
                            bool                       src_on_parent,
                            const void                *src,
                            void                      *dest);

/*----------------------------------------------------------------------------
 * Copy indexed array from distant or matching interface elements to
 * local elements.
 *
 * Source and destination arrays define values for all elements in the
 * interface set (i.e. all elements listed by cs_interface_get_elt_ids()
 * when looping over interfaces of a set,
 *
 * Note that when copying the same type of data to all matching elements,
 * the source and destination index may be the same, if src_on_parent is true.
 * To avoid requiring a separate destination index, the dest_index argument
 * may be set to NULL, in which case it is assumed that source and destination
 * are symmetric, and src_index is sufficient to determine sizes (whether
 * src_on_parent is true or not).
 *
 * In some use cases, for example when copying values only in one direction,
 * the copying is not symmetric, so both a source and destination buffer must
 * be provided.
 *
 * parameters:
 *   ifs           <-- pointer to interface set structure
 *   datatype      <-- type of data considered
 *   src_on_parent <-- true if source array is defined on the elements
 *                     defined by ifs->elt_ids, false if source array
 *                     defined directly on cs_interface_set_n_elts(ifs)
 *   src_index     <-- index for source array
 *   dest_index    <-- index for destination array, or NULL
 *   src           <-- source array (size:
 *                     src_index[cs_interface_set_n_elts(ifs)]
 *                     or parent array size * stride)
 *   dest          <-- destination array (size:
 *                     src_index[cs_interface_set_n_elts(ifs)] or
 *                     dest_index[cs_interface_set_n_elts(ifs)])
 *----------------------------------------------------------------------------*/

void
cs_interface_set_copy_indexed(const cs_interface_set_t  *ifs,
                              cs_datatype_t              datatype,
                              bool                       src_on_parent,
                              const cs_lnum_t            src_index[],
                              const cs_lnum_t            dest_index[],
                              const void                *src,
                              void                      *dest);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update values using the bitwise inclusive or operation for elements
 *        associated with an interface set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been processed.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * \param[in]       ifs        pointer to a fvm_interface_set_t structure
 * \param[in]       n_elts     number of elements in var buffer
 * \param[in]       stride     number of values (non interlaced) by entity
 * \param[in]       interlace  true if variable is interlaced (for stride > 1)
 * \param[in]       datatype   type of data considered
 * \param[in, out]  var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_inclusive_or(const cs_interface_set_t  *ifs,
                              cs_lnum_t                  n_elts,
                              cs_lnum_t                  stride,
                              bool                       interlace,
                              cs_datatype_t              datatype,
                              void                      *var);

/*----------------------------------------------------------------------------
 * Update the sum of values for elements associated with an interface set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been added.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * parameters:
 *   ifs       <-- pointer to a fvm_interface_set_t structure
 *   n_elts    <-- number of elements in var buffer
 *   stride    <-- number of values (non interlaced) by entity
 *   interlace <-- true if variable is interlaced (for stride > 1)
 *   datatype  <-- type of data considered
 *   var       <-> variable buffer
 *----------------------------------------------------------------------------*/

void
cs_interface_set_sum(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the sum of values for elements associated with an
 * interface set, allowing control over periodicity of rotation.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been added.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * \param[in]       ifs        pointer to a fvm_interface_set_t structure
 * \param[in]       n_elts     number of elements in var buffer
 * \param[in]       stride     number of values (non interlaced) by entity
 * \param[in]       interlace  true if variable is interlaced (for stride > 1)
 * \param[in]       datatype   type of data considered
 * \param[in]       ignore_rotation  ignore rotation if present ?
 * \param[in, out]  var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_sum_tr(const cs_interface_set_t  *ifs,
                        cs_lnum_t                  n_elts,
                        cs_lnum_t                  stride,
                        bool                       interlace,
                        cs_datatype_t              datatype,
                        bool                       ignore_rotation,
                        void                      *var);

/*----------------------------------------------------------------------------
 * Update to minimum value for elements associated with an interface set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been added.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * parameters:
 *   ifs       <-- pointer to a fvm_interface_set_t structure
 *   n_elts    <-- number of elements in var buffer
 *   stride    <-- number of values (non interlaced) by entity
 *   interlace <-- true if variable is interlaced (for stride > 1)
 *   datatype  <-- type of data considered
 *   var       <-> variable buffer
 *----------------------------------------------------------------------------*/

void
cs_interface_set_min(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------
 * Update to maximum value for elements associated with an interface set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been added.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * parameters:
 *   ifs       <-- pointer to a fvm_interface_set_t structure
 *   n_elts    <-- number of elements in var buffer
 *   stride    <-- number of values (non interlaced) by entity
 *   interlace <-- true if variable is interlaced (for stride > 1)
 *   datatype  <-- type of data considered
 *   var       <-> variable buffer
 *----------------------------------------------------------------------------*/

void
cs_interface_set_max(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------
 * Dump printout of an interface list.
 *
 * parameters:
 *   ifs <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
cs_interface_set_dump(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTERFACE_H__ */
