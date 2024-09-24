#ifndef __CS_EXECUTION_CONTEXT_H__
#define __CS_EXECUTION_CONTEXT_H__

/*============================================================================
 * Class to handle different execution policies (MPI, OpenMP, CUDA, ...)
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

// Valid only for C++

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
#ifdef __cplusplus
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * Class used to define execution context for different methods and algorithms.
 *
 * It handles both the distributed memory context (i.e. MPI communicator),
 * and shared memory context (i.e. OpenMP thread, Cuda stream, SYCL queue).
 */
/*----------------------------------------------------------------------------*/

class cs_execution_context {

/* Public methods
   --------------*/

public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Constructor.
   */
  /*--------------------------------------------------------------------------*/

  cs_execution_context()
  {
#if defined(HAVE_MPI)
    this->_comm = MPI_COMM_NULL;
#endif
    this->_comm_rank    = 0;
    this->_comm_n_ranks = 0;
    this->_thread_id    = 0;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Destructor.
   */
  /*--------------------------------------------------------------------------*/

  virtual ~cs_execution_context()
  {
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for MPI communicator rank id.
   *
   * \return rank id inside the related MPI communicator.
   */
  /*--------------------------------------------------------------------------*/

  int rank() const
  {
    return _comm_rank;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for MPI communicator number of ranks.
   *
   * \return number of ranks inside the related MPI communicator.
   */
  /*--------------------------------------------------------------------------*/

  int n_ranks() const
  {
    return _comm_n_ranks;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Does the execution context uses MPI parallelism ?
   *
   * \return true if MPI is used and n_ranks > 1, false otherwise
   */
  /*--------------------------------------------------------------------------*/

  bool use_mpi() const
  {
    bool retval = (_comm_n_ranks > 1) ? true : false;
    return retval;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Is the current cpu/task the MPI root (rank 0) ?
   *
   * \return true if caller is rank 0, false otherwise
   */
  /*--------------------------------------------------------------------------*/

  bool is_mpi_root() const
  {
    bool retval = (_comm_rank < 1) ? true : false;
    return retval;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for thread id.
   *
   * \return thread id
   */
  /*--------------------------------------------------------------------------*/

  int thread_id() const
  {
    return _thread_id;
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Free the MPI communicator.
   */
  /*--------------------------------------------------------------------------*/

  void
  comm_free()
  {
#if defined(HAVE_MPI)
    if (_comm != MPI_COMM_NULL)
      MPI_Comm_free(&(this->_comm));
    this->_comm = MPI_COMM_NULL;
#endif
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Call MPI_Barrier over the MPI communicator, do nothing if no MPI.
   */
  /*--------------------------------------------------------------------------*/

  int
  barrier()
  {
    int retval = 0;
#if defined(HAVE_MPI)
    if (_comm != MPI_COMM_NULL)
      retval = MPI_Barrier(_comm);
#endif
    return retval;
  }

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
#if defined(HAVE_MPI)
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Set the execution context MPI communicator.
   */
  /*--------------------------------------------------------------------------*/

  void
  set_comm
  (
    MPI_Comm comm /*!<[in] MPI communicator to set */
  )
  {
    int _initialized;
    MPI_Initialized(&_initialized);
    if (_initialized) {
      this->_comm = comm;
      MPI_Comm_rank(this->_comm, &(this->_comm_rank));
      MPI_Comm_size(this->_comm, &(this->_comm_n_ranks));
    }
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Getter function for MPI communicator.
   *
   * \return MPI communicator.
   */
  /*--------------------------------------------------------------------------*/

  MPI_Comm comm() const
  {
    return _comm;
  };

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
#endif // defined(HAVE_MPI)
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/* Private attributes
 -------------------- */

private:

#if defined(HAVE_MPI)
  MPI_Comm _comm; /*!< MPI Communicator */
#endif

  int _comm_rank;    /*!< Rank id inside the MPI communicator */
  int _comm_n_ranks; /*!< Number of ranks inside the MPI communicator */
  int _thread_id;    /*!< Thread id */

#if defined(__NVCC__)
  cudaStream_t _stream; /*!< Cuda stream */

#elif defined(SYCL_LANGUAGE_VERSION)
  sycl::queue _queue; /*!< SYCL queue */
#endif

};

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
#endif /* __cplusplus */
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

#endif /* __CS_EXECUTION_CONTEXT_H__ */
