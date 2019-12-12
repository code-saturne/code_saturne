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
# Python utility module, facilitating the use of PLE for the coupling of
# a python code with Code_Saturne.
# ==============================================================================

try:
    from ple import Init as ple_init
    from ple import Coupling as ple_coupling
    _ple_avail = True
except:
    _ple_avail = False

class pyple_coupler():
    """
    Coupling class based on PLE
    """

    # ----------------------------------
    def __init__(self, name, verbosity=1):
        """
        Init method.
        :param name: string
                     name of the app
        :param verbosity: int, optional
                          log verbosity level
                          default value is 1
        """

        self.name      = name
        self.log_name  = "saturne_%s_pycoupler_run.log" % name
        self.verbosity = verbosity
        if self.verbosity < 0:
            self.verbosity = 0

        self.last_cpl_iter = False

        self.init_log()

        self.has_ple = False
        # Check for MPI and PLE presence
        self.check_mpi_status()

        # Launch PLE parallel initialization
        self.coupling_init()
    # ----------------------------------


    # ----------------------------------
    def init_log(self):
        """
        Init the log file
        """
        f = open(self.log_name, "w")
        f.write(" ===================================== \n")
        f.write("  Code_Saturne/%s pycoupler log \n" % self.name)
        f.write(" ===================================== \n\n")
        f.close()
    # ----------------------------------


    # ----------------------------------
    def log(self, msg, verbosity = 0):
        """
        Write a message to the log file.
        :param msg:       string
                          message to be printed
        :param verbosity: int, optional
                          verbosity level of the message
                          default is 0
        """

        if self.verbosity >= verbosity:
            f = open(self.log_name, "a")
            f.write(" ===================================== \n")
            msg_lines = msg.split('\n')
            for line in msg_lines:
                f.write("  "+line+"\n")
            f.write(" ===================================== \n\n")
            f.close()
    # ----------------------------------


    # ----------------------------------
    def check_mpi_status(self):
        """
        Check that MPI and PLE are available and initialized.
        """

        self.log("MPI INIT AND CHECK FOR PLE PRESENCE")

        init_msg = ""

        if _ple_avail:
                self.has_ple = ple_init.mpi_initialized()
                if self.has_ple:
                    init_msg="Code %s was launched with PLE using MPI support" % name
                else:
                    init_msg="Code %s was launched, but without MPI..." % name
        else:
            self.has_ple = False
            init_msg="Code %s was launched, but could not import PLE..." % name

        self.log(init_msg)
        if not self.has_ple:
            raise Exception(init_msg)

        self.base_comm = ple_init.get_comm_world()
    # ----------------------------------


    # ----------------------------------
    def coupling_init(self):
        """
        Start the coupling communicator based on the PLE coupling library
        of Code_Saturne
        """
        # Names to be used by PLE
        self.app_type = 'PYTHON'

        if self.name == None:
            self.name = 'PYCODE'

        # PLE calls to create the subcommunicator
        self.app_num = ple_coupling.name_to_id(self.base_comm, self.name)

        world_rank = self.base_comm.rank

        if self.app_num > -1:
            self.my_comm = self.base_comm.Split(self.app_num, world_rank)
        else:
            self.my_comm = self.base_comm.Dup()

        _sync_flag = 0

        self.ple_set = ple_coupling.ple_mpi_set(_sync_flag,
                                                self.app_type,
                                                self.app_name,
                                                self.base_comm,
                                                self.my_comm)

        self.n_apps = self.ple_set.get_n_apps()
        self.app_id = self.ple_set.get_app_id()

        # Store the root ranks of each app, which is useful for direct
        # send/recv operations
        self.root_ranks = {}
        for app_id in range(self.n_apps):
            ai = self.ple_set.get_info(app_id)

            self.root_ranks[ai['app_name']] = ai['root_rank']

        # Print info concerning the coupling mpi set to log
        self.ple_set.dump()
    # ----------------------------------


    # ----------------------------------
    def get_app_ranks(self, app_name):
        """
        Retrieve the ranks in MPI_COMM_WORLD of a given app
        :param app_name: string
                         name of the app
        """

        app_ranks = []
        for app_id in range(self.n_apps):
            ai = self.ple_set.get_info(app_id)

            if ai['app_name'] == app_name:
                r0 = ai['root_rank']
                r1 = ai['root_rank'] + ai['n_ranks']

                app_ranks = [rank for rank in range(r0, r1)]
                break

        return app_ranks
    # ----------------------------------


    # ----------------------------------
    def sync_coupling_status(self, end_coupling=False):
        """
        Sync with the other codes/apps using PLE.
        Returns exit_status which tells us if the coupling needs to stop.
        0: continue
        1: stop

        :param end_coupling: bool, optional
                             Flag to force the end of the coupled run.
                             defualt is False
        """

        sync_flags = ple_coupling.PLE_COUPLING_FLAGS['NO_SYNC']
        sync_flags = sync_flags | ple_coupling.PLE_COUPLING_FLAGS['TS_MIN']

        if end_coupling:
            sync_flags = sync_flags | ple_coupling.PLE_COUPLING_FLAGS['STOP']

        self.ple_set.synchronize(sync_flags, 1000000000.0)

        exit_status = 0

        sync_msg = ''
        for app_id in range(self.n_apps):

            if app_id != self.app_id:
                ai = self.ple_set.get_info(app_id)
                if (ai['status'] & ple_coupling.PLE_COUPLING_FLAGS['STOP']) != 0:
                    sync_msg = "%s asked to stop the coupling" % (ai.app_name)
                    exit_status = 1
                    break

                if (ai['status'] & ple_coupling.PLE_COUPLING_FLAGS['LAST']) != 0:

                    if self.last_cpl_iter:
                        sync_msg =  "%s is at its last iteration. " % (ai.app_name)
                        sync_msg += "Coupling is to stop"
                        exit_status = 1
                        break

                    else:
                        exit_status = 0
                        self.last_cpl_iter = True
                        break

        if exit_status == 0 and end_coupling:
            exit_status = 1
            sync_msg = "%s asked to stop the coupling" % (self.name)

        if sync_msg != '':
            self.log(sync_msg)

        return exit_status

    # ----------------------------------


    # ----------------------------------
    def finalize(self, force_exit=True):
        """
        Finalize the coupling and exit if necessary
        :param force_exit: bool, optional
                           call sys.exit() to provide a proper exit
                           default value is True
        """

        if self.has_ple:
            self.log("I'm at the MPI_COMM_WORLD exit barrier", 1)
            self.base_comm.Barrier()

        self.log("Exiting, coupling is over")

        if force_exit:
            sys.exit(0)
    # ----------------------------------
