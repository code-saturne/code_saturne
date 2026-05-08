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
    mpi = nullptr;
    ctx = nullptr;
  }

  ~environment()
  {
    if (mpi != nullptr)
      delete mpi;

    if (ctx != nullptr)
      delete ctx;
  }

  mpi_wrapper         *mpi;
  cs_dispatch_context *ctx;
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
 * \brief Initialize the global execution context.
 */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_init(void);

/*--------------------------------------------------------------------------*/
/*!
 * \brief Initialize context (CPU/GPU)
 */
/*--------------------------------------------------------------------------*/

void
cs_execution_default_env_init_context(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the global execution context pointer.
 */
/*----------------------------------------------------------------------------*/

void
cs_execution_default_env_finalize(void);

/*----------------------------------------------------------------------------*/

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
  _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using a dispatch context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class F, class... Args>
bool
parallel_for(cs_dispatch_context& ctx, cs_lnum_t n, F&& f, Args&&... args) {
  ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief parallel_for construct using a host context as input argument
 */
/*--------------------------------------------------------------------------*/

template <class F, class... Args>
bool
parallel_for(cs_host_context& ctx, cs_lnum_t n, F&& f, Args&&... args) {
  ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for_reduce_sum(n,
                               sum,
                               static_cast<F&&>(f),
                               static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for_reduce_sum(n,
                               sum,
                               static_cast<F&&>(f),
                               static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for_reduce_sum(n,
                              sum,
                              static_cast<F&&>(f),
                              static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for_reduce_sum(n,
                              sum,
                              static_cast<F&&>(f),
                              static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for_reduce(n,
                           result,
                           reducer,
                           static_cast<F&&>(f),
                           static_cast<Args&&>(args)...);
  return true;
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
  _ctx.parallel_for_reduce(n,
                           result,
                           reducer,
                           static_cast<F&&>(f),
                           static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for_reduce(n,
                          result,
                          reducer,
                          static_cast<F&&>(f),
                          static_cast<Args&&>(args)...);
  return true;
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
  ctx.parallel_for_reduce(n,
                          result,
                          reducer,
                          static_cast<F&&>(f),
                          static_cast<Args&&>(args)...);
  return true;
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
  _ctx.wait();
  return true;
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
  ctx.wait();
  return true;
}

} // namespace cs

using cs_mpi_wrapper = cs::execution::mpi_wrapper;

/*--------------------------------------------------------------------------*/

#endif /* __CS_EXECUTION_CONTEXT_H__ */
