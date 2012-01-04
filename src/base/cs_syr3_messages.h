#ifndef __CS_SYR3_MESSAGES_H__
#define __CS_SYR3_MESSAGES_H__

/*============================================================================
 * Manage messages for Syrthes coupling: sending, receiving and interpolation
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if Syrthes coupling continues or if we must finalize communications.
 *
 * parameters:
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs <-> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_test_iter(int   nt_cur_abs,
                           int  *nt_max_abs);

/*----------------------------------------------------------------------------
 * Synchronize new time step
 *
 * parameters:
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs --> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_new_time_step(int  nt_cur_abs,
                               int  nt_max_abs);

/*----------------------------------------------------------------------------
 * Receive coupling variables from Syrthes
 *
 * parameters:
 *   syr_num <-- Syrthes 3 coupling number
 *   twall   --> wall temperature
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_recv_twall(cs_int_t   syr_num,
                            cs_real_t  twall[]);

/*----------------------------------------------------------------------------
 * Send coupling variables to Syrthes
 *
 * parameters:
 *   syr_num <-- Syrthes 3 coupling number
 *   tfluid  <-- wall exchange coefficient
 *   hwall   <-- wall exchange coefficient
 *----------------------------------------------------------------------------*/

void
cs_syr3_messages_send_tf_hwall(cs_int_t   syr_num,
                               cs_real_t  tfluid[],
                               cs_real_t  hwall[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR3_MESSAGES_H__ */
