# ==============================================================================
#  This file is part of the "Parallel Location and Exchange" library,
#  intended to provide mesh or particle-based code coupling services.
#
#  Copyright (C) 2005-2020  EDF S.A.
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
    def __init__(self,
                 run_dir=None,
                 verbosity=1,
                 logfile = None,
                 output_mode="master",):
        """
        Init method.
        :param name: string
                     name of the app
        :param verbosity:   int, optional
                            log verbosity level
                            default value is 1
                            If negative value, no output
        :param logfile:     string, optional
                            User defined name for the logfile.
                            Default is None, leading to the default value of
                            init_log method
        :param output_mode: string, optional
                            Log output mode. Options are "all" or "master"
                            default is "master"
        """

        self.name = None

        # Read command line arguments
        import sys
        self._parse_cmd_line(sys.argv)

        # Check if we need to change directory
        if run_dir:
            # user imposed directly
            import os
            os.chdir(run_dir)

        else:
            # imposed from from command line arguments (mpmd like launch)
            if self.cmd_line_options.wdir:
                run_dir = self.cmd_line_options.wdir
                import os
                os.chdir(run_dir)

        self.user_options = {'rundir':run_dir,
                             'verbos':verbosity,
                             'logfile':logfile,
                             'output':output_mode}

        self.verbosity = verbosity
        self.print_info = False

        self.last_cpl_iter = False

        self.has_ple = False
        # Check for MPI and PLE presence
        self.check_mpi_status()
    # ----------------------------------


    # ----------------------------------
    def _parse_cmd_line(self, argv):
        """
        Parse the command line in launch script
        :param argv: list
                     list of command line input arguments
        """
        from argparse import ArgumentParser

        parser = ArgumentParser()

        parser.add_argument("--app-type", "--type", "-t", dest="app_type",
                            metavar="<apptype>",
                            help="Application type")

        parser.add_argument("--app-name", "--name", dest="app_name",
                            metavar="<appname>",
                            help="Name of the application")

        parser.add_argument("--wdir", "-wdir", dest="wdir",
                            metavar="<wdir>",
                            help="Working directory for the application")

        parser.set_defaults(app_type=None)
        parser.set_defaults(app_name=None)
        parser.set_defaults(wdir=None)

        (options, args) = parser.parse_known_args(argv)

        self.cmd_line_options = options
    # ----------------------------------


    # ----------------------------------
    def init_log(self):
        """
        Init the log file
        default logfile name is:
        <app_name>_<rank>_pyple_coupler.log'
        """
        if self.verbosity > 0:
            if self.user_options['output'] == "all":
                self.print_info = True
            elif self.local_rank == 0:
                self.print_info = True

            if self.user_options['logfile'] != None:
                self.log_name = self.user_options['logfile']
            else:
                self.log_name = "%s_r%05d_pyple_coupler.log" % (self.name, self.local_rank)

        if self.print_info:
            f = open(self.log_name, "w")
            f.write(" ===================================== \n")
            f.write("  '%s' pycoupler usage log \n" % self.name)
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

        if self.verbosity >= verbosity and self.print_info:
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

        init_msg = ""

        if _ple_avail:
                self.has_ple = ple_init.mpi_initialized()
                if self.has_ple:
                    init_msg="Code was launched with PLE using MPI support"
                else:
                    init_msg="Code was launched, but without MPI..."
        else:
            self.has_ple = False
            init_msg="Code was launched, but could not import PLE..."

        if not self.has_ple:
            raise Exception(init_msg)

        self.base_comm  = ple_init.get_comm_world()
        self.world_rank = self.base_comm.rank
    # ----------------------------------


    # ----------------------------------
    def init_coupling(self, app_name=None, app_type=None):
        """
        Start the coupling communicator based on the PLE coupling library
        of Code_Saturne
        """

        # Define app name and type to be used by PLE.
        # First we check if the user provided a name and type to this method.
        # If not, we check for command line values.
        # If neither is provided, defaults are used: name=PYCODE and type=PYTHON

        if app_type != None:
            self.app_type = app_type
        elif self.cmd_line_options.app_type != None:
            self.app_type = self.cmd_line_options.app_type
        else:
            self.app_type = "PYTHON"

        if self.name == None:
            if app_name != None:
                self.name = app_name
            elif self.cmd_line_options.app_name != None:
                self.name = self.cmd_line_options.app_name
            else:
                self.name = 'PYCODE'

        # PLE calls to create the subcommunicator
        self.app_num = ple_coupling.name_to_id(self.base_comm, self.name)


        if self.app_num > -1:
            self.my_comm = self.base_comm.Split(self.app_num, self.world_rank)
        else:
            self.my_comm = self.base_comm.Dup()

        self.local_rank = self.my_comm.rank

        self.init_log()

        _sync_flag = 0

        self.ple_set = ple_coupling.ple_mpi_set(_sync_flag,
                                                self.app_type,
                                                self.name,
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
    def sync_coupling_status(self, end_coupling=False, app_dt=100000000.0):
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

        self.ple_set.synchronize(sync_flags, app_dt)

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
