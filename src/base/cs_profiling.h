#pragma once

/*============================================================================
 * Building with profiling annotations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

// Numeric values for different profiling libraries.

/// No profiling enabled.
#define CS_PROFILING_NONE 0

/// Enable NVTX profiling.
#define CS_PROFILING_NVTX 1

// No profiling library is used by default.
#ifndef CS_PROFILING
#define CS_PROFILING CS_PROFILING_NONE
#endif

// Make sure NVTX profiling is used in host code only, not CUDA.
#if defined(__CUDA_ARCH__)
#undef CS_PROFILING
#define CS_PROFILING CS_PROFILING_NONE
#endif

#define CS_COMBINE_DETAIL(x, y) x##y
#define CS_COMBINE(x, y) CS_COMBINE_DETAIL(x, y)

#define CS_STRINGIFY_DETAIL(x) #x
#define CS_STRINGIFY(x) CS_STRINGIFY_DETAIL(x)

#if CS_PROFILING == CS_PROFILING_NONE

/*----------------------------------------------------------------------------
 * Default: no profiling library
 *----------------------------------------------------------------------------*/

/// Annotates a whole function.
#define CS_PROFILE_FUNC_RANGE()

/// Annotates a range delimited by the lifetime of a variable.
#define CS_PROFILE_RANGE(range_name)

/// Adds a mark in a profile that corresponds to the current file and line.
#define CS_PROFILE_MARK_LINE()

/*----------------------------------------------------------------------------
 * Profiling with NVIDIA NVTX3
 *----------------------------------------------------------------------------*/

#elif CS_PROFILING == CS_PROFILING_NVTX

#include <nvtx3/nvtx3.hpp>

/// Annotates a whole function.
#ifndef __CUDA_ARCH__
#define CS_PROFILE_FUNC_RANGE() NVTX3_FUNC_RANGE()
#else
#define CS_PROFILE_FUNC_RANGE()
#endif

/// Annotates a range delimited by the lifetime of a variable.
#define CS_PROFILE_RANGE(range_name)                                           \
  nvtx3::scoped_range CS_COMBINE(__cs__profile_range_, __LINE__){              \
    range_name " (" __FILE__ ":" CS_STRINGIFY(__LINE__) ")"                    \
  };

/// Adds a mark in a profile that corresponds to the current file and line.
#define CS_PROFILE_MARK_LINE() nvtx3::mark(__FILE__ ":" CS_STRINGIFY(__LINE__))

#endif

/*----------------------------------------------------------------------------*/
