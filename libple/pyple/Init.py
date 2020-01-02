# ==============================================================================
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2020  EDF S.A.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# ==============================================================================

# ==============================================================================
# MPI UTILITIES FOR EASIER USE BY EXTERNAL CODES
# ==============================================================================

# ==============================================================================
def mpi_initialized():
    """
    This function checks if the python script was launched with MPI.
    """

    try:
        import mpi4py
    except ImportError:
        return False

    try:
        import mpi4py.MPI
    except ImportError:
        return False

    if mpi4py.MPI.COMM_WORLD.size == 1 and mpi4py.MPI.COMM_WORLD.rank == 0:
        return False

    return True
# ==============================================================================

# ==============================================================================
def get_comm_world():
    """
    This function returns the MPI_COMM_WORLD communicator which can then
    be used in Python for the different MPI operations.
    """

    from mpi4py import MPI

    base_comm = MPI.COMM_WORLD

    return base_comm
# ==============================================================================
