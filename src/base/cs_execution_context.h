#ifndef __CS_EXECUTION_CONTEXT_H__
#define __CS_EXECUTION_CONTEXT_H__

/*============================================================================
 * Class to handle different execution policies (MPI, OpenMP, CUDA, ...)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include "base/cs_dispatch.h"

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>

/*============================================================================
 * Type definitions
 *============================================================================*/

namespace cs {
namespace execution {
class mpi_wrapper {

public:

CS_F_HOST
mpi_wrapper()
{};

CS_F_HOST
~mpi_wrapper()
{};

CS_F_HOST_DEVICE
inline
int
rank() const
{
  return _comm_rank;
};

CS_F_HOST_DEVICE
inline
int
n_ranks() const
{
  return _comm_n_ranks;
}

CS_F_HOST_DEVICE
inline
bool
active() const
{
  bool retval = (_comm_n_ranks > 1) ? true : false;
  return retval;
}

CS_F_HOST_DEVICE
inline
bool
is_root() const
{
  bool retval = (_comm_rank < 1) ? true : false;
  return retval;
}

CS_F_HOST
void
free()
{
#if defined(HAVE_MPI)
  if (_comm != MPI_COMM_NULL)
    MPI_Comm_free(&(this->_comm));
  this->_comm = MPI_COMM_NULL;
#endif
}

#if defined(HAVE_MPI)
CS_F_HOST
void
set_comm
(
  MPI_Comm comm
)
{
  int _initialized;
  MPI_Initialized(&_initialized);
  if (_initialized) {
    this->_comm = comm;
    MPI_Comm_rank(this->_comm, &(this->_comm_rank));
    MPI_Comm_size(this->_comm, &(this->_comm_n_ranks));
  }
}

CS_F_HOST
inline
MPI_Comm
comm() const
{
  return _comm;
}

#endif

private:

#if defined(HAVE_MPI)
  MPI_Comm _comm = MPI_COMM_NULL; /*!< MPI Communicator */
#endif

  int _comm_rank = -1;    /*!< Rank id inside the MPI communicator */
  int _comm_n_ranks = 0; /*!< Number of ranks inside the MPI communicator */
};

class environment {

public:
  environment()
  {
  }

  ~environment()
  {
  }

  mpi_wrapper         mpi;
  cs_dispatch_context g_ctx;
};

const environment *
default_env(void);

cs_dispatch_context&
default_context(void);

cs_host_context&
default_h_context(void);

mpi_wrapper&
default_mpi(void);

}; // namespace execution
}; // namespace cs

void
cs_execution_default_env_init(void);

void
cs_execution_default_env_finalize();

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

    this->g_ctx = cs_dispatch_context();
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
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Call MPI_Barrier over the MPI communicator, do nothing if no MPI.
   */
  /*--------------------------------------------------------------------------*/

  int
  barrier
  (
    const bool verbosity = false, /*!<[in] Activate verbosity mode */
#if (defined(__GNUC__) || defined(__clang__)) && \
  __has_builtin(__builtin_LINE) && \
  __has_builtin(__builtin_FILE)
    const char *file_name   = __builtin_FILE(), /*!<[in] Caller file (for log) */
    const int   line_number = __builtin_LINE()  /*!<[in] Caller line (for log) */
#else
    const char *file_name   = __FILE__, /*!<[in] Caller file (for log) */
    const int   line_number = __LINE__  /*!<[in] Caller line (for log) */
#endif
  ) const;

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

public:
  cs_dispatch_context g_ctx;
};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the current execution context. For the moment only global
 * context is returned.
 *
 * \return pointer to current execution context.
 */
/*----------------------------------------------------------------------------*/

const cs_execution_context *
cs_execution_context_get(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the global execution context.
 *
 * \return pointer to global execution context.
 */
/*----------------------------------------------------------------------------*/

const cs_execution_context *
cs_execution_context_glob_get(void);

/*--------------------------------------------------------------------------*/
/*!
 * \brief Get the global default dispatch context
 *
 * \return reference of default global dispatch context
 */
/*--------------------------------------------------------------------------*/

cs_dispatch_context&
cs_execution_context_glob_get_ctx(void);

/*--------------------------------------------------------------------------*/
/*!
 * \brief  Get the global default host context
 *
 * \return reference of default global host context
 */
/*--------------------------------------------------------------------------*/

cs_host_context&
cs_execution_context_glob_get_h_ctx(void);

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
#endif /* __cplusplus */
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the global execution context.
 */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the global execution context pointer.
 */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_finalize(void);

/*----------------------------------------------------------------------------*/


#if defined(__cplusplus)

namespace cs {

enum class exec_type {
  device,
  host,
  unknown
};


/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using a dispatch context based on template
 * parameter (enum) "exec_type"
 */
/*--------------------------------------------------------------------------*/

template<exec_type _EXEC_, class F, class... Args>
bool
parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
  auto& _ctx = (_EXEC_ != exec_type::host) ?
    cs::execution::default_context() : cs::execution::default_h_context();
  return _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using the default global dispatch context
 */
/*--------------------------------------------------------------------------*/

template <class F, class... Args>
bool
parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
  auto& _ctx = cs::execution::default_context();
  return _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using a dispatch context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class F, class... Args>
bool
parallel_for(cs_dispatch_context& ctx, cs_lnum_t n, F&& f, Args&&... args) {
  return ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using a host context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class F, class... Args>
bool
parallel_for(cs_host_context& ctx, cs_lnum_t n, F&& f, Args&&... args) {
  return ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

// paralell_for_i_faces
/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_i_faces construct using a dispatch context based on
 * template parameter (enum) "exec_type"
 */
/*--------------------------------------------------------------------------*/

template <exec_type _EXEC_, class M, class F, class... Args>
bool
parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
  const cs_lnum_t n = m->n_i_faces;
  auto& _ctx = (_EXEC_ != exec_type::host) ?
    cs::execution::default_context() : cs::execution::default_h_context();
  return _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_i_faces construct using the default global dispatch
 * context
 */
/*--------------------------------------------------------------------------*/

template <class M, class F, class... Args>
bool
parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
  const cs_lnum_t n = m->n_i_faces;
  auto& _ctx = cs::execution::default_context();
  return _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_i_faces construct using dispatch context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class M, class F, class... Args>
bool
parallel_for_i_faces(cs_dispatch_context& ctx, const M* m, F&& f, Args&&... args) {
  const cs_lnum_t n = m->n_i_faces;
  return ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_i_faces construct using host context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class M, class F, class... Args>
bool
parallel_for_i_faces(cs_host_context& ctx, const M* m, F&& f, Args&&... args) {
  const cs_lnum_t n = m->n_i_faces;
  return ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using a dispatch context based on
 * template parameter (enum) "exec_type"
 */
/*--------------------------------------------------------------------------*/

template <exec_type _EXEC_, class T, class F, class... Args>
bool
parallel_for_reduce_sum(
  cs_lnum_t n,
  T&        sum,
  F&&       f,
  Args&&... args)
{
  auto& _ctx = (_EXEC_ != exec_type::host) ?
    cs::execution::default_context() : cs::execution::default_h_context();
  return _ctx.parallel_for_reduce_sum(n,
                                      sum,
                                      static_cast<F&&>(f),
                                      static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using the default global
 * dispatch context
 */
/*--------------------------------------------------------------------------*/

template <class T, class F, class... Args>
bool
parallel_for_reduce_sum(
  cs_lnum_t n,
  T&        sum,
  F&&       f,
  Args&&... args)
{
  auto& _ctx = cs::execution::default_context();
  return _ctx.parallel_for_reduce_sum(n,
                                      sum,
                                      static_cast<F&&>(f),
                                      static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using the default global
 * dispatch context
 */
/*--------------------------------------------------------------------------*/

template <class T, class F, class... Args>
bool
parallel_for_reduce_sum(
  cs_dispatch_context& ctx,
  cs_lnum_t            n,
  T&                   sum,
  F&&                  f,
  Args&&...            args)
{
  return ctx.parallel_for_reduce_sum(n,
                                     sum,
                                     static_cast<F&&>(f),
                                     static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using the default global
 * host context
 */
/*--------------------------------------------------------------------------*/

template <class T, class F, class... Args>
bool
parallel_for_reduce_sum(
  cs_host_context& ctx,
  cs_lnum_t        n,
  T&               sum,
  F&&              f,
  Args&&...        args)
{
  return ctx.parallel_for_reduce_sum(n,
                                     sum,
                                     static_cast<F&&>(f),
                                     static_cast<Args&&>(args)...);
}

// Parallel reduction with general reducer.
/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce construct using a dispatch context based on
 * template parameter (enum) "exec_type"
 */
/*--------------------------------------------------------------------------*/

template <exec_type _EXEC_, class T, class R, class F, class... Args>
bool
parallel_for_reduce (
  cs_lnum_t n,
  T&        result,
  R&        reducer,
  F&&       f,
  Args&&... args)
{
  auto& _ctx = (_EXEC_ != exec_type::host) ?
    cs::execution::default_context() : cs::execution::default_h_context();
  return _ctx.parallel_for_reduce(n,
                                  result,
                                  reducer,
                                  static_cast<F&&>(f),
                                  static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce construct using the default global
 * dispatch context
 */
/*--------------------------------------------------------------------------*/

template <class T, class R, class F, class... Args>
bool
parallel_for_reduce (
  cs_lnum_t n,
  T&        result,
  R&        reducer,
  F&&       f,
  Args&&... args)
{
  auto& _ctx = cs::execution::default_context();
  return _ctx.parallel_for_reduce(n,
                                  result,
                                  reducer,
                                  static_cast<F&&>(f),
                                  static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using the default global
 * dispatch context
 */
/*--------------------------------------------------------------------------*/

template <class T, class R, class F, class... Args>
bool
parallel_for_reduce (
  cs_dispatch_context& ctx,
  cs_lnum_t            n,
  T&                   result,
  R&                   reducer,
  F&&                  f,
  Args&&...            args)
{
  return ctx.parallel_for_reduce(n,
                                 result,
                                 reducer,
                                 static_cast<F&&>(f),
                                 static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for_reduce_sum construct using the default global
 * host context
 */
/*--------------------------------------------------------------------------*/

template <class T, class R, class F, class... Args>
bool
parallel_for_reduce (
  cs_host_context& ctx,
  cs_lnum_t        n,
  T&               result,
  R&               reducer,
  F&&              f,
  Args&&...        args)
{
  return ctx.parallel_for_reduce(n,
                                 result,
                                 reducer,
                                 static_cast<F&&>(f),
                                 static_cast<Args&&>(args)...);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief wait construct using a dispatch context based on template parameter
 * of exec_type (enum)
 */
/*--------------------------------------------------------------------------*/

template<exec_type _EXEC_, class... Args>
bool
wait(void)
{
  auto& _ctx = (_EXEC_ != exec_type::host) ?
    cs::execution::default_context() : cs::execution::default_h_context();
  return _ctx.wait();
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief wait construct using the default global dispatch context
 */
/*--------------------------------------------------------------------------*/

template<class... Args>
bool
wait
(void)
{
  auto& _ctx = cs::execution::default_context();
  _ctx.wait();
  return true;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief wait construct using a dispatch context as input argument
 */
/*--------------------------------------------------------------------------*/

template<class... Args>
bool
wait
(
  cs_dispatch_context& ctx
)
{
  ctx.wait();
  return true;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief wait construct using a host context as input argument
 */
/*--------------------------------------------------------------------------*/

template<class... Args>
bool
wait
(
  cs_host_context& ctx
)
{
  return ctx.wait();
}

}

using cs_mpi_wrapper = cs::execution::mpi_wrapper;
#endif
#endif /* __CS_EXECUTION_CONTEXT_H__ */
