#ifndef CS_REDISTRIBUTE_H
#define CS_REDISTRIBUTE_H

/*============================================================================
 * Redistribution of mesh and field data.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_all_to_all.h"
#include "base/cs_algorithm.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

struct cs_distributor_t {
  cs_all_to_all_t *d;
  cs_lnum_t n_send;
  cs_lnum_t n_recv;
  cs_lnum_t n_uniq;
  cs_lnum_t *send_list;
  cs_lnum_t *recv_order;
  cs_lnum_t *ordered_to_local;
};

#endif

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute a buffer using the distributor structure.
 *
 * This function guarantees the uniqueness of the elements post-redistribution.
 *
 * \param[in]       db          pointer to the distributor
 * \param[in]       stride      stride of the array to distribute
 * \param[in,out]   buffer      pointer to the array to redistribute
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_distribute_buffer(const cs_distributor_t   *db,
                     const int                 stride,
                     T                        *buffer[])
{
  if (buffer == nullptr || db == nullptr)
    return;

  T *_buffer = nullptr;

  if (db->send_list) {

    T *send_data = nullptr;
    CS_MALLOC(send_data, stride * db->n_send, T);

    _buffer = *buffer;

    for (cs_lnum_t i = 0; i < db->n_send; i++) {
      cs_lnum_t id = db->send_list[i];
      for (int j = 0; j < stride; j++)
        send_data[stride*i+j] = _buffer[stride*id+j];
    }

    CS_FREE(*buffer);
    *buffer = send_data;

  }

  T *recv_data = cs_all_to_all_copy_array(db->d,
                                          stride,
                                          false,
                                          *buffer);

  CS_REALLOC(*buffer, db->n_uniq * stride, T);
  _buffer = *buffer;

  for (cs_lnum_t i = 0; i < db->n_recv; i++) {
    cs_lnum_t ordered = db->recv_order[i];
    cs_lnum_t local = db->ordered_to_local[ordered];
    if (local == -1) continue;

    for (int j = 0; j < stride; j++)
      _buffer[stride*local+j] = recv_data[stride*ordered+j];
  }

  CS_FREE(recv_data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute a buffer using the distributor structure.
 *
 * Unlike cs_distribute_buffer, the distribution is not in-place, and the
 * recv_buf array must be allocated to the correct size by the caller.
 *
 * This function guarantees the uniqueness of the elements post-redistribution.
 *
 * \param[in]   db          pointer to the distributor
 * \param[in]   stride      stride of the array to distribute
 * \param[in]   send_buf    pointer to the array to redistribute
 * \param[out]  recv_buf    pointer to the redistributed array
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_distribute_buffer_allocated(const cs_distributor_t   *db,
                               const int               stride,
                               const T                 send_buf[],
                               T                       recv_buf[])
{
  if (recv_buf == nullptr || db == nullptr)
    return;

  T *_send_buf = nullptr;

  if (db->send_list) {

    CS_MALLOC(_send_buf, stride * db->n_send, T);

    for (cs_lnum_t i = 0; i < db->n_send; i++) {
      cs_lnum_t id = db->send_list[i];
      for (int j = 0; j < stride; j++)
        _send_buf[stride*i+j] = send_buf[stride*id+j];
    }

    send_buf = _send_buf;

  }

  T *_recv_buf = cs_all_to_all_copy_array(db->d,
                                          stride,
                                          false,
                                          send_buf);

  CS_FREE(_send_buf);

  for (cs_lnum_t i = 0; i < db->n_recv; i++) {
    cs_lnum_t ordered = db->recv_order[i];
    cs_lnum_t local = db->ordered_to_local[ordered];
    if (local == -1) continue;

    for (int j = 0; j < stride; j++)
      recv_buf[stride*local+j] = _recv_buf[stride*ordered+j];
  }

  CS_FREE(_recv_buf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Distribute an indexed array and its corresponding index array using
 * the distributor structure.
 *
 * This function guarantees the uniqueness of the elements post-redistribution.
 *
 * \param[in]        db        pointer to the distributor
 * \param[in,out]    idx       pointer to the index array to redistribute
 * \param[in,out]    buffer    pointer to the indexed array to redistribute
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_distribute_buffer_indexed(const cs_distributor_t *db,
                             cs_lnum_t            *idx[],
                             T                    *buffer[])
{
  T *_buffer = nullptr;

  cs_lnum_t *idx_s = nullptr;
  CS_MALLOC(idx_s, db->n_send+1, cs_lnum_t);

  cs_lnum_t *_idx = *idx;

  for (cs_lnum_t i = 0; i < db->n_send; i++) {
    cs_lnum_t id = db->send_list ? db->send_list[i] : i;
    idx_s[i] = _idx[id+1] - _idx[id];
  }

  cs_lnum_t *idx_r = nullptr;
  CS_MALLOC(idx_r, db->n_recv+1, cs_lnum_t);

  cs_all_to_all_copy_array(db->d,
                           1,
                           false,
                           idx_s,
                           idx_r);

  // Build the buffer to send.

  cs_dispatch_context  ctx;
  ctx.set_use_gpu(false);

  cs::algorithm::count_to_index(ctx, db->n_send, idx_s);

  if (db->send_list) {

    T* send_data = nullptr;
    CS_MALLOC(send_data, idx_s[db->n_send], T);

    T* _send_data = send_data;
    _buffer = *buffer;

    for (cs_lnum_t i = 0; i < db->n_send; i++) {
      cs_lnum_t id = db->send_list ? db->send_list[i] : i;
      for (cs_lnum_t j = _idx[id]; j < _idx[id+1]; j++)
        *_send_data++ = _buffer[j];
    }

    assert(_send_data - send_data == idx_s[db->n_send]);

    CS_FREE(*buffer);
    CS_FREE(*idx);

    *buffer = send_data;
  }

  // Build the correct idx.

  CS_REALLOC(*idx, db->n_uniq+1, cs_lnum_t);
  _idx = *idx;

  for (cs_lnum_t i = 0; i < db->n_recv; i++) {
    cs_lnum_t ordered = db->recv_order[i];
    cs_lnum_t local = db->ordered_to_local[ordered];
    if (local == -1) continue;

    _idx[local] = idx_r[ordered];
  }

  cs::algorithm::count_to_index(ctx, db->n_uniq, _idx);

  // Exchange buffer.

  cs::algorithm::count_to_index(ctx, db->n_recv, idx_r);

  T *_recv_data = cs_all_to_all_copy_indexed(db->d,
                                             false,
                                             idx_s,
                                             *buffer,
                                             idx_r);

  CS_FREE(idx_s);

  // Build the unique buffer.

  CS_REALLOC(*buffer, _idx[db->n_uniq], T);
  _buffer = *buffer;

  for (cs_lnum_t i = 0; i < db->n_recv; i++) {
    cs_lnum_t ordered = db->recv_order[i];
    cs_lnum_t local = db->ordered_to_local[ordered];
    if (local == -1) continue;

    T *dst = _buffer + _idx[local];

    for (cs_lnum_t j = idx_r[ordered]; j < idx_r[ordered+1]; j++)
      *dst++ = _recv_data[j];
  }

  CS_FREE(_recv_data);
  CS_FREE(idx_r);
}

void
cs_distributor_destroy(cs_distributor_t  **db);

void
cs_redistribute(const int           cell_dest_rank[],
                cs_distributor_t  **cd,
                cs_distributor_t  **ifd,
                cs_distributor_t  **bfd,
                cs_distributor_t  **vd);

#endif // defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/

#endif /* CS_REDISTRIBUTE_H */
