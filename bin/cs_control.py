#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

"""
This module describes the script used to run a study/case for Code_Saturne.

This module defines the following functions:
- process_cmd_line
- main
"""

#===============================================================================
# Import required Python modules
#===============================================================================

import os, sys
import types, string, re, fnmatch

import socket

from optparse import OptionParser

#===============================================================================
# Classes
#===============================================================================

class controller:
    """
    Controller class for running computation.
    """

    #---------------------------------------------------------------------------

    def __init__(self,
                 path = None,      # run directory
                 package = None):  # main package

        # Package specific information

        self.package = package

        # run directory

        self.path = path

        # Initialize connection

        import random

        key = str(random.randrange(0, 2**31))
        hostname = socket.getfqdn()

        self.s = socket.socket()
        self.s.bind(('', 0)) # OS will pick an available port.

        port = self.s.getsockname()[1]

        f_path = 'control_file'
        if self.path != None:
            f_path = os.path.append(self.path, f_path)
        c = open(f_path, 'w')
        c.write('connect ' + hostname+':'+str(port) + ' ' + key + '\n')
        c.close()

        self.s.listen(0)
        (self.conn, self.address) = self.s.accept()

        cmp_key = self.conn.recv(len(key)).decode("utf-8")
        if cmp_key != key:
            print('incorrect key returned: expected ' + key + ', got ' + cmp_key)

        magic_string = 'CFD_control_comm_socket'
        cmp_string = self.conn.recv(len(magic_string)).decode("utf-8")

        if cmp_string != magic_string:
            print('incorrect magic string returned: expected ' + magic_string)
            print('got ' + cmp_string)

        self.conn.send(magic_string.encode("utf-8"))

    #---------------------------------------------------------------------------

    def send_string(self, str):

        self.conn.send((str + "\n").encode("utf-8"))
        # print("sent ", str)

    #---------------------------------------------------------------------------

    def recv_string(self, l=32768):

        str = self.conn.recv(l).decode("utf-8")
        lr = len(str)
        while lr >= l:
            str0 = self.conn.recv(l).decode("utf-8")
            l = len(str0)
            str += str0
        return str

    #---------------------------------------------------------------------------

    def advance(self, n = 1):

        self.send_string("advance " + str(n))
        retcode = self.recv_string()
        # print(retcode)

    #---------------------------------------------------------------------------

    def disconnect(self):

        self.send_string("disconnect ")

        self.conn.shutdown(socket.SHUT_RDWR)
        self.conn.close()

        self.s.shutdown(socket.SHUT_RDWR)
        self.s.close()

    #---------------------------------------------------------------------------

    def __del__(self):

        self.disconnect()

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog run [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("--exec-dir", dest="exec_dir", type="string",
                      metavar="<case>",
                      help="server's exection directory")

    parser.set_defaults(exec_dir=None)

    # Note: we could use args to pass a calculation status file as an argument,
    # which would allow pursuing the later calculation stages.

    (options, args) = parser.parse_args(argv)

    return  (options.exec_dir)

#===============================================================================
# Main function
#===============================================================================

def main(argv, pkg):
    """
    Main function.
    """

    (exec_dir) = process_cmd_line(argv, pkg)

    c = controller(exec_dir, pkg)

    c.advance(3)
    c.advance(3)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Run package

    try:
        from cs_package import package
        pkg = package()
    except Exception:
        pkg = None

    retval = main(sys.argv[1:], pkg)

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

