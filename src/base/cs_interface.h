#ifndef __CS_INTERFACE_H__
#define __CS_INTERFACE_H__

/*============================================================================
 * Main structure for handling of interfaces associating mesh elements
 * (such as inter-processor or periodic connectivity between cells, faces,
 * or vertices);
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/
/*
 * \brief Return process rank associated with an interface's distant elements.
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  process rank associated with the interface's distant elements
 */
/*----------------------------------------------------------------------------*/

int
cs_interface_rank(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return number of local and distant elements defining an interface.
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  number of local and distant elements defining the interface
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_size(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to array of local element ids defining an interface.
 *
 * The size of the array may be obtained by cs_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  pointer to array of local element ids (0 to n-1) defining
 * the interface
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_elt_ids(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to array of matching element ids defining an interface.
 *
 * This array is only available if cs_interface_set_add_match_ids() has
 * been called for the containing interface set.
 *
 * The size of the array may be obtained by cs_interface_size().
 * The array is owned by the interface structure, and is not copied
 * (hence the constant qualifier for the return value).
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  pointer to array of local element ids (0 to n-1) defining
 * the interface
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_match_ids(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return size of index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * Its size should thus be equal to 1 + number of periodic transforms + 1,
 * In absence of periodicity, it may be 0, as the index is not needed.
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  transform index size for the interface
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_get_tr_index_size(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to index of sub-sections for different transformations.
 *
 * The index is applicable to both local_num and distant_num arrays,
 * with purely parallel equivalences appearing at position 0, and
 * equivalences through periodic transform i at position i+1;
 * In absence of periodicity, it may be NULL, as it is not needed.
 *
 * \param[in]  itf  pointer to interface structure
 *
 * \return  pointer to transform index for the interface
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_interface_get_tr_index(const cs_interface_t  *itf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Creation of a list of interfaces between elements of a same type.
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
 * need not be defined prior to this function being called.
 *
 * \param[in]  n_elts                 number of local elements considered
 *                                    (size of parent_element_id[])
 * \param[in]  parent_element_id      pointer to list of selected elements
 *                                    local ids (0 to n-1), or NULL if all
 *                                    first n_elts elements are used
 * \param[in]  global_number          pointer to list of global (i.e. domain
 *                                    splitting independent) element numbers
 * \param[in]  periodicity            periodicity information (NULL if none)
 * \param[in]  n_periodic_lists       number of periodic lists (may be local)
 * \param[in]  periodicity_num        periodicity number (1 to n) associated
 *                                    with each periodic list (primary
 *                                    periodicities only)
 * \param[in]  n_periodic_couples     number of periodic couples associated
 *                                    with each periodic list
 * \param[in]  periodic_couples       array indicating periodic couples
 *                                    (interlaced, using global numberings)
 *                                    for each list
 *
 * \return  pointer to list of interfaces (possibly NULL in serial mode)
 */
/*----------------------------------------------------------------------------*/

cs_interface_set_t *
cs_interface_set_create(cs_lnum_t                 n_elts,
                        const cs_lnum_t           parent_element_id[],
                        const cs_gnum_t           global_number[],
                        const fvm_periodicity_t  *periodicity,
                        int                       n_periodic_lists,
                        const int                 periodicity_num[],
                        const cs_lnum_t           n_periodic_couples[],
                        const cs_gnum_t    *const periodic_couples[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Destruction of an interface set.
 *
 * \param[in, out]  ifs  pointer to pointer to structure to destroy
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_destroy(cs_interface_set_t  **ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Duplicate an interface set, applying an optional constant stride.
 *
 * \param[in, out]  ifs     pointer to interface set structure
 * \param[in]       stride  if > 1, each element subdivided in stride elements
 *
 * \return  pointer to new interface set
 */
/*----------------------------------------------------------------------------*/

cs_interface_set_t  *
cs_interface_set_dup(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  stride);

/*----------------------------------------------------------------------------*/
/*
 * \brief Duplicate an interface set for coupled variable blocks.
 *
 * \param[in, out]  ifs         pointer to interface set structure
 * \param[in]       block_size  local block size (number of elements)
 * \param[in]       n_blocks    number of associated blocks
 *
 * \return  pointer to new interface set
 */
/*----------------------------------------------------------------------------*/

cs_interface_set_t  *
cs_interface_set_dup_blocks(cs_interface_set_t  *ifs,
                            cs_lnum_t            block_size,
                            cs_lnum_t            n_blocks);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return number of interfaces associated with an interface set.
 *
 * \param[in]  ifs  pointer to interface set structure
 *
 * \return  number of interfaces in set
 */
/*----------------------------------------------------------------------------*/

int
cs_interface_set_size(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return total number of elements in interface set.
 *
 * This is equal to the sum of cs_interface_size() on the cs_interface_size()
 * interfaces of a set.
 *
 * \param[in]  ifs  pointer to interface set structure
 *
 * \return  number of interfaces in set
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_interface_set_n_elts(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to a given interface in an interface set.
 *
 * \param[in]  ifs           pointer to interface set structure
 * \param[in]  interface_id  index of interface in set (0 to n-1)
 *
 * \return  pointer to interface structure
 */
/*----------------------------------------------------------------------------*/

const cs_interface_t *
cs_interface_set_get(const cs_interface_set_t  *ifs,
                     int                        interface_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to the periocicity structure associated of an
 * interface set.
 *
 * \param[in]  ifs  pointer to interface set structure
 *
 * \return  pointer to periodicity structure, or NULL
 */
/*----------------------------------------------------------------------------*/

const fvm_periodicity_t *
cs_interface_set_periodicity(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Apply renumbering of elements referenced by an interface set.
 *
 * For any given element i, a negative old_to_new[i] value means that that
 * element does not appear anymore in the new numbering.
 *
 * \param[in, out]  ifs         pointer to interface set structure
 * \param[in]       old_to_new  renumbering array (0 to n-1 numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_renumber(cs_interface_set_t  *ifs,
                          const cs_lnum_t      old_to_new[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy array from distant or matching interface elements to
 * local elements.
 *
 * Source and destination arrays define values for all elements in the
 * interface set (i.e. all elements listed by cs_interface_get_elt_ids()
 * when looping over interfaces of a set.
 *
 * \param[in]   ifs            pointer to interface set structure
 * \param[in]   datatype       type of data considered
 * \param[in]   stride         number of values per entity (interlaced)
 * \param[in]   src_on_parent  true if source array is defined on the elements
 *                             defined by ifs->elt_ids, false if source array
 *                             defined directly on cs_interface_set_n_elts(ifs)
 * \param[in]   src            source array (size:
 *                             cs_interface_set_n_elts(ifs)*stride
 *                             or parent array size * stride)
 * \param[out]  dest           destination array (size:
 *                             cs_interface_set_n_elts(ifs)*stride)
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_copy_array(const cs_interface_set_t  *ifs,
                            cs_datatype_t              datatype,
                            int                        stride,
                            bool                       src_on_parent,
                            const void                *src,
                            void                      *dest);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy indexed array from distant or matching interface elements to
 * local elements.
 *
 * Source and destination arrays define values for all elements in the
 * interface set (i.e. all elements listed by cs_interface_get_elt_ids()
 * when looping over interfaces of a set.
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
 * \param[in]   ifs            pointer to interface set structure
 * \param[in]   datatype       type of data considered
 * \param[in]   src_on_parent  true if source array is defined on the elements
 *                             defined by ifs->elt_ids, false if source array
 *                             defined directly on cs_interface_set_n_elts(ifs)
 * \param[in]   src_index      index for source array
 * \param[in]   dest_index     index for destination array, or NULL
 * \param[in]   src            source array (size:
 *                             src_index[cs_interface_set_n_elts(ifs)]
 *                             or parent array size)
 * \param[out]  dest           destination array (size:
 *                             src_index[cs_interface_set_n_elts(ifs)] or
 *                             dest_index[cs_interface_set_n_elts(ifs)])
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_copy_indexed(const cs_interface_set_t  *ifs,
                              cs_datatype_t              datatype,
                              bool                       src_on_parent,
                              const cs_lnum_t            src_index[],
                              const cs_lnum_t            dest_index[],
                              const void                *src,
                              void                      *dest);

/*----------------------------------------------------------------------------*/
/*
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

/*----------------------------------------------------------------------------*/
/*
 * \brief Update the sum of values for elements associated with an
 * interface set.
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
 * \param[in, out]  var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_sum(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update the sum of values for elements associated with an
 * interface set, allowing control over periodicity.
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
 * \param[in]       tr_ignore  if > 0, ignore periodicity with rotation;
 *                             if > 1, ignore all periodic transforms
 * \param[in, out]  var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_sum_tr(const cs_interface_set_t  *ifs,
                        cs_lnum_t                  n_elts,
                        cs_lnum_t                  stride,
                        bool                       interlace,
                        cs_datatype_t              datatype,
                        int                        tr_ignore,
                        void                      *var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update the minimum value of elements associated with an interface
 *        set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been updated.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * \param[in]      ifs        pointer to a cs_interface_set_t structure
 * \param[in]      n_elts     number of elements in var buffer
 * \param[in]      stride     number of values (non interlaced) by entity
 * \param[in]      interlace  true if variable is interlaced (for stride > 1)
 * \param[in]      datatype   type of data considered
 * \param[in, out] var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_min(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update the maximum value of elements associated with an interface
 *        set.
 *
 * On input, the variable array should contain local contributions. On output,
 * contributions from matching elements on parallel or periodic boundaries
 * have been updated.
 *
 * Only the values of elements belonging to the interfaces are modified.
 *
 * \param[in]      ifs        pointer to a cs_interface_set_t structure
 * \param[in]      n_elts     number of elements in var buffer
 * \param[in]      stride     number of values (non interlaced) by entity
 * \param[in]      interlace  true if variable is interlaced (for stride > 1)
 * \param[in]      datatype   type of data considered
 * \param[in, out] var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_max(const cs_interface_set_t  *ifs,
                     cs_lnum_t                  n_elts,
                     cs_lnum_t                  stride,
                     bool                       interlace,
                     cs_datatype_t              datatype,
                     void                      *var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update the maximum of values for elements associated with an
 * interface set, allowing control over periodicity.
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
 * \param[in]       tr_ignore  if > 0, ignore periodicity with rotation;
 *                             if > 1, ignore all periodic transforms
 * \param[in, out]  var        variable buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_max_tr(const cs_interface_set_t  *ifs,
                        cs_lnum_t                  n_elts,
                        cs_lnum_t                  stride,
                        bool                       interlace,
                        cs_datatype_t              datatype,
                        int                        tr_ignore,
                        void                      *var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Add matching element id information to an interface set.
 *
 * This information is required by calls to cs_interface_get_match_ids(),
 * and may be freed using cs_interface_set_free_match_ids().
 *
 * \param[in]  ifs  pointer to interface set structure
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_add_match_ids(cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Free matching element id information of an interface set.
 *
 * This information is used by calls to cs_interface_get_match_ids(),
 * and may be defined using cs_interface_set_add_match_ids().
 *
 * \param[in]  ifs  pointer to interface set structure
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_free_match_ids(cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/
/*
 * \brief Tag mutiple elements of local interface with a given values.
 *
 * This is effective only on an interface matching the current rank,
 * and when multiple (periodic) instances of a given element appear on that
 * rank, al instances except the first are tagged with the chosen value.
 *
 * \param[in]       itf          pointer to interface structure
 * \param[in]       periodicity  periodicity information (NULL if none)
 * \param[in]       tr_ignore   if > 0, ignore periodicity with rotation;
 *                              if > 1, ignore all periodic transforms
 * \param[in]       tag_value   tag to assign
 * \param[in, out]  tag         global tag array for elements
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_tag_local_matches(const cs_interface_t     *itf,
                               const fvm_periodicity_t  *periodicity,
                               int                       tr_ignore,
                               cs_gnum_t                 tag_value,
                               cs_gnum_t                *tag);
/*----------------------------------------------------------------------------*/
/*
 * \brief Dump printout of an interface list.
 *
 * \param[in]  ifs  pointer to structure that should be dumped
 */
/*----------------------------------------------------------------------------*/

void
cs_interface_set_dump(const cs_interface_set_t  *ifs);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTERFACE_H__ */
