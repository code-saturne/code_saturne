//-------------------------------------------------------------------------------
//   This file is part of the "Parallel Location and Exchange" library,
//   intended to provide mesh or particle-based code coupling services.
//
//   Copyright (C) 2005-2018  EDF S.A.
//
//   This library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License as published by the Free Software Foundation; either
//   version 2.1 of the License, or (at your option) any later version.
//
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the GNU Lesser General Public
//   License along with this library; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//-------------------------------------------------------------------------------

// file: swig.i

%module pyple
%{
#include "ple_coupling.h"
#include "ple_coupling.c"
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

typedef struct {
  int          status;
  int          root_rank;
  int          n_ranks;
  const char  *app_type;
  const char  *app_name;
} ple_coupling_mpi_set_info_t;

%include "ple_coupling.h"
%include "ple_coupling.c"
