#ifndef __CS_ASSERT_H__
#define __CS_ASSERT_H__

/*============================================================================
 * File and directory operations, with parallel file I/O
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "bft_error.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_assert.h
        Assertion usable in all build types.
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Abort the program if the given assertion is false.
 *
 * Contrary to the standard C assert macro, this is always defined, so
 * the standard C assert macro should be use for debug build-only checking,
 * whereas this macros should be used mainly for less expensive argument or
 * precondition checks that should be present in all build types.
 *
 * \param[in]  expr  expression to verify
 */
/*----------------------------------------------------------------------------*/

# define cs_assert(expr)                                       \
if (!(expr)) bft_error(__FILE__, __LINE__, 0,                  \
                       "Assertion failed in function %s: %s",  \
                       __func__, # expr)

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FILE_H__ */
