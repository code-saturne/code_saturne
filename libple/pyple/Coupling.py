# ==============================================================================
#  This file is part of the "Parallel Location and Exchange" library,
#  intended to provide mesh or particle-based code coupling services.
#
#  Copyright (C) 2005-2019  EDF S.A.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# ==============================================================================

# ==============================================================================
# Python wrapping for ple_coupling_* functions calls.
# Contains the ple_mpi_set python class which provides methods which call
# upon the ple_coupling C functions.
# ==============================================================================

# ==============================================================================
# Import the .so libray created using .c wrapper
import libpyplecoupling as pyplec
# ==============================================================================

# ==============================================================================
# Retrieve the PLE flags masks
PLE_COUPLING_FLAGS = pyplec.coupling_get_masks()
# ==============================================================================

# ==============================================================================
def name_to_id(Comm, name):
    """
    Return the id of the app with name "name" within the communicator
    """

    new_id = pyplec.coupling_mpi_name_to_id(Comm, name)

    return new_id
# ==============================================================================

# ==============================================================================
# Finalize usage of the Coupling module
def finalize():

    pyplec.coupling_mpi_set_destroy_all()

    return

# ==============================================================================

# ==============================================================================
# MPI set class
class ple_mpi_set(object):
    """
    Python class which allows a manipulation of a ple_coupling_mpi_set_t.
    This class also contains methods which correspond to the different
    PLE functions which call upon a ple_coupling_mpi_set.
    """

    def __init__(self, sync_flags, app_type, app_name, base_comm, app_comm):
        """
        Initialization function
        """

        try:
            self.idx = pyplec.coupling_mpi_set_create(sync_flags,
                                                      app_type,
                                                      app_name,
                                                      base_comm,
                                                      app_comm)
        except:
            raise Exception("ple mpi set create failed...")

        self.n_apps = pyplec.coupling_mpi_set_n_apps(self.idx)

        self.app_id = pyplec.coupling_mpi_set_get_app_id(self.idx)

    def get_info(self, app_id):
        """
        Return app information under the form of a dictionnary.
        """

        return pyplec.coupling_mpi_set_get_info(self.idx, app_id)

    def get_n_apps(self):
        """
        Return the number of apps in the MPI set
        """

        return self.n_apps

    def get_app_id(self):
        """
        Return the app id of the code making the call
        """

        return self.app_id

    def synchronize(self, sync_flags, time_step):
        """
        Synchronize the different codes in the MPI set.
        Flag choices are available in Coupling.PLE_FLAGS dictionnary
        """

        pyplec.coupling_mpi_set_synchronize(self.idx, sync_flags, time_step)

    def get_status(self):
        """
        Return a list of status' for all applications within the set
        """

        return pyplec.coupling_mpi_set_get_status(self.idx)

    def get_timesteps(self):
        """
        Return a list containing the time steps of all applications after a
        call to synchronize.
        """

        return pyplec.coupling_mpi_set_get_timestep(self.idx)

    def dump(self):
        """
        Dump MPI set information to stdout.
        """

        pyplec.coupling_mpi_set_dump(self.idx)
# ==============================================================================

