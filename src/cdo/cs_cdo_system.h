#ifndef __CS_CDO_SYSTEM_H__
#define __CS_CDO_SYSTEM_H__

/*============================================================================
 * Structure and functions used to manipulate elementary structures related to
 * the definition/usage of linear systems (objects manipulated with this
 * structure are matrices, matrix strcutures and assemblers, range set and
 * interface set for parallel synchronization/operations)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_assembly.h"
#include "cs_cdo_connect.h"
#include "cs_matrix.h"
#include "cs_matrix_assembler.h"
#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_CDO_SYSTEM_DEFAULT,      /*!> One matrix with one block */
  CS_CDO_SYSTEM_COUPLED,      /*!> System arising of a system of equations  */
  CS_CDO_SYSTEM_SADDLE_POINT, /*!> System arising in coupled (Navier-)Stokes */

} cs_cdo_system_type_t;

/*!
 * \enum cs_cdo_system_block_type_t
 * \brief type of block composing a (block) matrix
 *
 * \var CS_CDO_SYSTEM_BLOCK_DEFAULT
 * Simplest block type. One matrix inside a block.
 *
 * \var CS_CDO_SYSTEM_BLOCK_SPLIT
 * The block is described by several (sub)matrices.
 *
 * \var CS_CDO_SYSTEM_BLOCK_UNASS
 * The block is unassembled so that there is no matrix and related structures
 * to describe this type of block. For instance, this block can be associated
 * to the pressure block in some strategy for solving sadle-point problems.
 *
 * \var CS_CDO_SYSTEM_BLOCK_EXT
 * External way to define a block. This block does not rely on the generic
 * function to assemble and build the associated structures. For instance, it
 * corresponds to the full assembly strategy in case of a saddle-point problem.
 */

typedef enum {

  CS_CDO_SYSTEM_BLOCK_DEFAULT,
  CS_CDO_SYSTEM_BLOCK_SPLIT,
  CS_CDO_SYSTEM_BLOCK_UNASS,
  CS_CDO_SYSTEM_BLOCK_EXT,

  CS_CDO_SYSTEM_N_BLOCK_TYPES

} cs_cdo_system_block_type_t;

/*!
 * \enum cs_cdo_system_matrix_class_t
 * \brief Class of matrices to consider
 *
 * \var CS_CDO_SYSTEM_MATRIX_NONE
 * No matrix is considered (unassembled treatment for instance)
 *
 * \var CS_CDO_SYSTEM_MATRIX_CS
 * Matrix format specific to code_saturne (native)
 *
 * \var CS_CDO_SYSTEM_MATRIX_PETSC
 * Matrix format specific to the library PETSC (external)
 *
 * \var CS_CDO_SYSTEM_MATRIX_HYPRE
 * Matrix format specific to the library HYPRE (external)
 */

typedef enum {

  CS_CDO_SYSTEM_MATRIX_NONE,
  CS_CDO_SYSTEM_MATRIX_CS,
  CS_CDO_SYSTEM_MATRIX_PETSC,
  CS_CDO_SYSTEM_MATRIX_HYPRE,

  CS_CDO_SYSTEM_N_MATRIX_CLASSES

} cs_cdo_system_matrix_class_t;


typedef struct {

  /* If the block is associated to one or several matrices, one keeps
   * the information about the class of the matrix
   */

  cs_cdo_system_matrix_class_t    matrix_class;

  /* An element is for instance a vertex, edge, face or cell */

  cs_flag_t   location;    /*!> DoF and element location */
  cs_lnum_t   n_elements;  /*!> number of elements (for this location) */
  int         stride;      /*!> number of DoFs by element */

  /* Description how is stored the matrix values for this block.
   *
   * When stride > 1, there are several options to store the matrix
   * (all possibilities are not available up to now)
   *
   * unrolled: True means that a vector-valued DoF for instance is treated as
   * a scalar-valued DoF which is three times larger.
   *
   * interlaced: True means that a vector-valued DoF is stored as follows:
   * U0x U0y U0z U1x U1y U1z ... UNx UNy UNz
   * false means that its representation is as follows:
   * U0x U1x ... UNx U0y U1y ... UNy U0z U1z ... UNz
   */

  bool        unrolled;
  bool        interlaced;

} cs_cdo_system_block_info_t;


/*! \struct cs_cdo_system_dblock_t
 *  \brief Structure associated to the default type of block
 *
 * The standard case is to have only one (sub-)matrix which is the system
 * matrix. This type of block is used in scalar-valued and vector-valued
 * problems and corresponds to the simplest usage of block to define a linear
 * system.
 */

typedef struct {

  /*
   * \var matrix
   *      matrix of the system to solve
   *
   * \var matrix_class
   *      class of the matrix
   *
   * \var mav
   *      structure to manage the assembly of matrix values
   *
   * \var assembly_func
   *      function pointer to operate the assembly stage
   *
   * \var slave_assembly_func
   *      function pointer to operate the assembly stage when the system helper
   *      is declared as slave (this is the same for all matrices). Useful for
   *      coupled systems.
   */

  cs_matrix_t                    *matrix;
  cs_matrix_assembler_values_t   *mav;
  cs_cdo_assembly_func_t         *assembly_func;
  cs_cdo_assembly_func_t         *slave_assembly_func;

  /* The following structures can be shared if the same block configuration
     is requested */

  cs_range_set_t                 *range_set;
  cs_interface_set_t             *interface_set;
  cs_matrix_assembler_t          *matrix_assembler;
  cs_matrix_structure_t          *matrix_structure;

} cs_cdo_system_dblock_t;


/*! \struct cs_cdo_system_sblock_t
 *  \brief Structure associated to the split type of block
 *
 * In this case, the block is split into several matrices sharing the same
 * pattern (structure and way to be assembled)
 */

typedef struct {

  /*
   * \var n_matrices
   *      number of matrices to consider
   *
   * \var matrices
   *      matrix of the system to solve (the matrix id in the list of matrices
   *      corresponding to the matrix at row i and column j is i*n_matrices+j)
   *
   * \var mav_array
   *      array of structures to manage the assembly of matrix values
   *
   * \var assembly_func
   *      function pointer to operate the assembly stage (this is the same for
   *      all matrices)
   *
   * \var slave_assembly_func
   *      function pointer to operate the assembly stage when the system helper
   *      is declared as slave (this is the same for all matrices)
   */

  int                             n_matrices;
  cs_matrix_t                   **matrices;
  cs_matrix_assembler_values_t  **mav_array;
  cs_cdo_assembly_func_t         *assembly_func;
  cs_cdo_assembly_func_t         *slave_assembly_func;

  /* The two following structures are always shared. They correspond to the
     scalar-valued version. */

  cs_range_set_t                 *range_set;
  cs_interface_set_t             *interface_set;

  bool                            matrix_struct_ownership;
  cs_matrix_assembler_t          *matrix_assembler;
  cs_matrix_structure_t          *matrix_structure;

} cs_cdo_system_sblock_t;

/*! \struct cs_cdo_system_ublock_t
 *  \brief Structure associated to the unassembled type of block
 *
 * With this type of block, there is no matrix and its associated
 * structures. Only the pointer to a range set and its associated interface set
 * are stored.
 */

typedef struct {

  /*! \var values
   *       Shared. values of the operator. According to the stride value in the
   *       block metadata, one can apply a stride to access to the array of
   *       values.
   *
   * \var _values
   *      Private. The structure has the ownership of the operator
   *      values. According to the stride value in the block metadata, one can
   *      apply a stride to access to the array of values.
   *
   *  \var adjacency
   *       Always shared. The associated adjacency to acces to the values. It
   *       corresponds to the connectivity x2y. n_x is equal to
   *       adjacency->n_elts and n_y is equal to info.n_elts
   */

  cs_real_t               *values;
  cs_real_t               *_values;

  const cs_adjacency_t    *adjacency;

  /* The following structures can be shared if the same block configuration is
     requested. In this case, shared_structures is set to true. By default,
     there is no sharing. */

  bool                     shared_structures;
  cs_range_set_t          *range_set;
  cs_interface_set_t      *interface_set;

} cs_cdo_system_ublock_t;


/*! \struct cs_cdo_system_xblock_t
 *  \brief Structure associated to the extended type of block
 *
 * There is no predefined usage for this type of block. There is no shared
 * structure in this case and no pointer to a generic assembly function.
 */

typedef struct {

  /*
   * \var matrix
   *      matrix of the system to solve
   *
   * \var mav
   *      structure to manage the assembly of matrix values
   */

  cs_matrix_t                    *matrix;
  cs_matrix_assembler_values_t   *mav;

  /* Private structures */

  cs_range_set_t                 *range_set;
  cs_interface_set_t             *interface_set;
  cs_matrix_assembler_t          *matrix_assembler;
  cs_matrix_structure_t          *matrix_structure;

} cs_cdo_system_xblock_t;


/*! \struct cs_cdo_system_block_t
 *          Block structure which defines a part of the globality of the system
 *          matrix
 */

typedef struct {

  /* \var info
   *      Set of metadata to describe the block
   */

  cs_cdo_system_block_info_t   info;

  /*! \var type
   *  Type of block
   *
   * \var owner
   *      If true, the lifecycle is managed by this structure
   *
   * \var id
   *      block id in the list of defined block structures
   */

  cs_cdo_system_block_type_t   type;

  bool                         owner;

  int                          id;

  /*!
   * \var block_pointer
   *      Untyped block. One has to cast the block on-the-fly according to the
   *      type of block
   */

  void                        *block_pointer;

} cs_cdo_system_block_t;

/*! \struct cs_cdo_system_helper_t
 *          Structure to handle linear systems which may be simply described by
 *          one block and one matrix or can be defined by a set of different
 *          blocks in more complex situations. The layout of the matrix is
 *          detailed thanks to block of different nature. The rhs can be also
 *          defined by block or thanks to a full-length array.
 */

typedef struct {

  /*!
   * \var type
   *      type of system. The type predefines some behaviors.
   *
   * \var n_col_blocks
   *      number of blocks in a row
   *
   * \var col_block_sizes
   *      size of each block (size = n_elts * stride). This size corresponds to
   *      the number of DoFs in a "scattered" view
   *
   * \var max_col_block_sizes
   *      max. size of each block. This size takes into account DoFs shared
   *      across ranks. In a sequential run, the values in the
   *      max_col_block_sizes array should be equal to the values in the
   *      col_block_sizes array.
   *
   * \var full_rhs_size
   *      size of the rhs when taking into account all column blocks
   *
   * \var rhs_array
   *      rhs array which is split in as many parts as there are blocks in a
   *      row (size n_col_blocks)
   *
   * \var rhs
   *      If allocated, it means that the system has its own rhs which is
   *      allocated in one array of the maximal size
   *
   * \var n_blocks
   *      number of blocks used to describe the system layout
   *
   * \var blocks
   *      array of pointers to block structures
   */

  cs_cdo_system_type_t     type;

  int                      n_col_blocks;
  cs_lnum_t               *col_block_sizes;
  cs_lnum_t               *max_col_block_sizes;

  cs_lnum_t                full_rhs_size;
  cs_real_t               *rhs;    /* shared */
  cs_real_t               *_rhs;   /* private */
  cs_real_t              **rhs_array;

  int                      n_blocks;
  cs_cdo_system_block_t  **blocks;

} cs_cdo_system_helper_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize matrix-related structures according to
 *         the type of discretization used for this simulation
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_init_sharing(cs_mesh_t           *mesh,
                           cs_cdo_connect_t    *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a system_helper structure from its set of metadata.
 *        n_col_blocks and n_blocks may differ according to the settings. For
 *        instance, for a saddle-point system, n_col_blocks >= n_blocks
 *
 * \param[in] type             type of system to handle
 * \param[in] n_col_blocks     number of blocks in a row
 * \param[in] col_block_sized  number of DoFs in each block of the row
 * \param[in] n_blocks         number of blocks associated to this system
 *
 * \return the pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_helper_t *
cs_cdo_system_helper_create(cs_cdo_system_type_t    type,
                            int                     n_col_blocks,
                            const cs_lnum_t        *col_block_sizes,
                            int                     n_blocks);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an unassembled block definition at position "block_id" in the
 *        helper structure Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      matclass    class of the matrix to handle
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 * \param[in]      interlaced  useful if stride > 1; way to store components
 * \param[in]      unrolled    useful if stride > 1; true=as scalar-valued
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_dblock(cs_cdo_system_helper_t       *sh,
                         int                           block_id,
                         cs_cdo_system_matrix_class_t  matclass,
                         cs_flag_t                     location,
                         cs_lnum_t                     n_elements,
                         int                           stride,
                         bool                          interlaced,
                         bool                          unrolled);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a split block definition at position "block_id" in the helper
 *        structure. Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      matclass    class of the matrix to handle
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_sblock(cs_cdo_system_helper_t       *sh,
                         int                           block_id,
                         cs_cdo_system_matrix_class_t  matclass,
                         cs_flag_t                     location,
                         cs_lnum_t                     n_elements,
                         int                           stride);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an unassembled block definition at position "block_id" in the
 *        helper structure Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      adjacency   shared adjacency structure
 * \param[in]      location    where DoFs are defined
 * \param[in]      n_elements  number of elements (support entities for DoFs)
 * \param[in]      stride      number of DoFs by element
 * \param[in]      interlaced  useful if stride > 1; way to store components
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_ublock(cs_cdo_system_helper_t   *sh,
                         int                       block_id,
                         const cs_adjacency_t     *adjacency,
                         cs_flag_t                 location,
                         cs_lnum_t                 n_elements,
                         int                       stride,
                         bool                      interlaced);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an external block definition at position "block_id" in the helper
 *        structure. Only metadata are set at this stage.
 *
 * \param[in, out] sh          pointer to the system helper to update
 * \param[in]      block_id    id in blocks array in a system helper
 * \param[in]      n_dofs      number of degrees of freedom
 *
 * \return a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_system_block_t *
cs_cdo_system_add_xblock(cs_cdo_system_helper_t   *sh,
                         int                       block_id,
                         cs_lnum_t                 n_dofs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to the given block_id
 *
 * \param[in]  sh         pointer to the system_helper structure to update
 * \param[in]  block_id   id of the block to consider
 *
 * \return a pointer to a range set structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_range_set_t *
cs_cdo_system_get_range_set(const cs_cdo_system_helper_t  *sh,
                            int                            block_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the matrix associated to the given block_id. If the type of
 *        block is either CS_CDO_SYSTEM_BLOCK_DEFAULT or
 *        CS_CDO_SYSTEM_BLOCK_EXT. In other cases, a NULL pointer is
 *        returned. The unassembled block has no matrix and to get a matrix of
 *        a split block, one should use cs_cdo_system_get_sub_matrix(sh,
 *        block_id, sub_id)
 *
 * \param[in]  sh         pointer to the system_helper structure to update
 * \param[in]  block_id   id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_cdo_system_get_matrix(const cs_cdo_system_helper_t  *sh,
                         int                            block_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the (sub-)matrix associated to a split block with id equal
 *        to block_id. sub_id is the id in the list of matrices of size equal
 *        to stride*stride.
 *        If the type of the block is not CS_CDO_SYSTEM_BLOCK_SPLIT, then a
 *        NULL pointer is returned.
 *
 * \param[in, out]  sh         pointer to the system_helper structure to update
 * \param[in]       block_id   id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_cdo_system_get_sub_matrix(cs_cdo_system_helper_t  *sh,
                             int                      block_id,
                             int                      sub_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the associated structures for the given system_helper structure
 *        If a member is already allocated, one keeps the member as it is.
 *
 * \param[in, out]  sh         pointer to the system_helper structure to update
 * \param[in]       block_id   specific block to handle or -1 for all blocks
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_build_block(cs_cdo_system_helper_t  *sh,
                          int                      block_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize the matrix, rhs and the matrix assembler
 *        values
 *
 * \param[in, out] sh       pointer to a system helper structure
 * \param[in, out] p_rhs    double pointer to the RHS array to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_init_system(cs_cdo_system_helper_t    *sh,
                                 cs_real_t                **p_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize the assembly after the cellwise building and assembly
 *
 * \param[in, out]      sh       pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_finalize_assembly(cs_cdo_system_helper_t    *sh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free matrix and rhs after the solve step
 *
 * \param[in, out]      sh       pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_reset(cs_cdo_system_helper_t    *sh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_cdo_system_helper_t structure
 *
 * \param[in, out] p_helper    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_helper_free(cs_cdo_system_helper_t   **p_helper);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize shared assembly structures from the existing helper
 *        structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_allocate_assembly(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all members associated to system helpers
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_system_destroy_all(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_SYSTEM_H__ */
