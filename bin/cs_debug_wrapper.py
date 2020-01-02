#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
This module defines a wrapper to launch an executable under a debugger.
"""

#===============================================================================
# Import required Python modules
#===============================================================================

import sys, os, stat
import subprocess

#-------------------------------------------------------------------------------
# Global state variables
#-------------------------------------------------------------------------------

# MPI rank

rank_id = -1

# List of supported debuggers

debuggers = {"gdb": "GNU gdb debugger",
             "cgdb": "Console front-end to gdb",
             "gdbgui": "gdbgui gdb web browser interface",
             "ddd": "Data Display Debugger",
             "emacs": "Emacs with gdb debugger",
             "emacs23": "Emacs 23 or older with gdb debugger",
             "kdbg": "KDbg",
             "kdevelop": "Kdevelop",
             "gede": "Gede",
             "nemiver": "Nemiver"}

#-------------------------------------------------------------------------------
# Enquote arguments if required
#-------------------------------------------------------------------------------

def enquote_arg(s):
    """
    Add quotes around argument if it contains whitespace, leave it
    unchanged otherwise; if the argument already contains unprotected
    quotes, do not add any more (so for example --option="string 1"
    is unchanged).
    """

    if s:
        if (s.find(' ') > -1):
            protect = False
            for i in range(len(s)):
                if s[i] == '\\':
                    protect = not protect
                if not protect and s[i] == '"':
                    return s
            return '"' + s + '"'
        else:
            return s
    else:
        return s

#-------------------------------------------------------------------------------
# Print a help page.
#-------------------------------------------------------------------------------

def print_help():
    """
    Print a help page.
    """

    help_string = \
"""
This is a debugger launcher wrapper. Usage:

%s [debugger opts] [mpiexec opts] [valgrind opts] --program=<program> [arguments]

This wrapper may be run either under MPI, or include an mpi run command.
If no debugger or Valgrind options are specified, the program is run under
the gdb debugger.

Debugger options:

  --debugger             Indicates the program should be debugged
  --debugger=DEBUGGER    Allows selection of the debugger
  DEBUGGER               Same as above (if debugger in known list)
  --debugger=list        Lists debuggers supported by this script
  --asan-bp              Adds a breakpoint for gcc's Address-Sanitizer
  --back-end=GDB         Path to debugger back-end (for graphical front-ends)
  --breakpoints=LIST     Comma-separated list of breakpoints to insert
  --terminal=TERM        Select terminal type to use for console debugger

  Other, standard options specific to each debugger may also be
  used, as long as they do not conflict with options in this script.

MPI options:

  If it recognizes mpi execution commands, this launcher will re-run
  itself under MPI, using the provided options. This allows simply
  prefixing debug options to a command-line, and may allow for
  support of parallel debuggers.

Valgrind options:

  --valgrind             Indicates the program should be run under Valgrind
  --valgrind=VALGRIND    Allows selection of the valgrind path
  VALGRIND               Same as above

  Other valgrind options may be used. Most importantly, when the
  --vgdb-errors=<num> option is used, the progam is run under a Valgrind
  gdb server, and debugged with the specified gdb debugger interface.
  When no debugger option is provided, or the (deprecated)
  --db-attach=yes option is used, the program is run under Valgrind,
  possibly under separate terminals when run under MPI.

Program options:

  --program=PROGRAM      Specifies the program that should be debugged
  PROGRAM                Same as above

  Program arguments should follow the program name or path.

"""
    print(help_string % sys.argv[0])

#-------------------------------------------------------------------------------
# Process the command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Process the passed command line arguments.
    """

    positions = {"debugger": -1,
                 "valgrind": -1,
                 "mpiexec": -1,
                 "program": -1}

    debugger_options = ("asan-bp", "back-end", "breakpoints", "terminal")

    files = []

    p = ['.'] + os.getenv('PATH').split(':')

    # First loop on options to determine option groups and executables

    idx = 0
    for a in argv:

        ie = a[2:].find("=")
        if ie > -1:
            b = a[2:ie+2]
        else:
            b = a[2:]
        if b in positions:
            if positions[b] < 0 or positions[b] > idx:
                positions[b] = idx
            if a == '--debugger=list':
                for k in debuggers.keys():
                    print(k + ': ' + debuggers[k])
                return None

        if not a[0] == '-':
            # Check for executable files
            rname = None
            ename = os.path.expandvars(os.path.expanduser(a))
            sname = os.path.split(ename)
            if os.path.isabs(ename) or (len(sname) > 1 and sname[0]):
                rname = os.path.realpath(ename)
            else:
                for d in p:
                    absname = os.path.join(d, ename)
                    if os.path.isfile(absname) or os.path.islink(absname):
                        rname = os.path.realpath(absname)
                        break
            if rname:
                if os.path.isfile(rname):
                    files.append((a, idx))

        idx += 1

    # Check if mpiexec or equivalent is inside the command, so
    # as to place it first (we want to debug the application,
    # not mpiexec). Also check in case Valgrind is passed
    # directly as a path, not --valgrind.

    for f in files:
        b = os.path.basename(f[0]).split('.')[0]
        if b in ['mpiexec', 'mpirun', 'srun', 'aprun', 'poe']:
            positions['mpiexec'] = f[1]
        elif b in debuggers.keys():
            positions['debugger'] = f[1]
        elif b == 'valgrind' and positions['valgrind'] == -1:
            positions['valgrind'] = f[1]

    # Now check for executable:

    if positions['program'] == -1:
        s_id = -1
        for k in ['mpiexec', 'debugger', 'valgrind']:
            if positions[k] > s_id:
                s_id = positions[k]
        for f in files:
            if f[1] > s_id and positions['program'] == -1:
                p = os.path.realpath(f[0])
                if os.path.isfile(p) or os.path.islink(p):
                    if os.stat(p).st_mode & stat.S_IXUSR:
                        positions['program'] = f[1]

    # Second check if Valgrind is used, based on the fact that all
    # Valgrind options start with "-".

    if positions['program'] == -1 and positions['valgrind'] > -1:
        i = positions['valgrind'] + 1
        for p in argv[i:]:
            if p[0] != '-':
                positions['program'] = i
                break
            else:
                i += 1

    # Check for error

    if positions['program'] == -1:
        for a in argv:
            if a == '-h' or a == '--help':
                print_help()
                return None
        print("usage: %s [debugger opts] [mpiexec opts] [valgrind opts] "
               % sys.argv[0] + "--program=<program> [arguments]")
        return None

    # Now get all parts

    o_idx = list(positions.values())
    o_idx.sort()

    cmds = {}
    for s in positions.keys():
        if positions[s] > -1:
            s_idx = positions[s]
            e_idx = len(argv)
            for p in o_idx:
                if p > s_idx and p < e_idx:
                    e_idx = p
            # Unify --<tool>=<path> with <path> logic
            tool = argv[s_idx]
            ie = tool[2:].find("=")
            if tool[0] == '-' and ie > -1:
                tool = tool[ie+3:]
            cmds[s] = [tool] + argv[s_idx+1:e_idx]

    # Debugger options

    if positions['debugger'] == -1:
        need_gdb = True
        if 'valgrind' in cmds.keys():
            need_gdb = False
            for k in cmds['valgrind'][1:]:
                if k[0:6] == '--vgdb':
                    need_gdb = True
            for k in cmds['valgrind'][1:]:
                if k == '--vgdb=off':
                    need_gdb = False
        if need_gdb:
            debugger = 'gdb'
            p_min = len(argv)
            for s in positions.keys():
                if positions[s] > -1 and positions[s] < p_min:
                    p_min = positions[s]
            cmds['debugger'] = [debugger] + argv[0:p_min]

    # For for debugger options which might have appeared first
    # check only now to avoid minor risk of similarly-named options
    # for other tools

    idx_s = 0
    idx_e = -1
    for k in ('debugger', 'mpiexec', 'valgrind', 'program'):
        if positions[k] > -1:
            idx_e = positions[k]
            break;

    idx_mpi = positions['mpiexec'] # in case placed first in the future
    if idx_mpi > -1 and idx_mpi < idx_e:
        idx_e = idx_mpi

    idx = 0
    for idx in range(idx_s, idx_e):
        a = argv[idx]
        ie = a[2:].find("=")
        if ie > -1:
            b = a[2:ie+2]
        else:
            b = a[2:]

        if b in debugger_options:
            idx_s = idx
            break

    try:
        if cmds['debugger']:
            for idx in range(idx_s, idx_e):
                cmds['debugger'].insert(idx-idx_s + 1, argv[idx])
    except Exception:
        pass

    return cmds

#-------------------------------------------------------------------------------
# Determine rank id when run under MPI
#-------------------------------------------------------------------------------

def init_rank_id():
    """
    Determine rank id when run under MPI.
    """

    rank_env_vars = ['PMI_RANK', 'OMPI_COMM_WORLD_RANK', 'PMI_ID',
                     'SLURM_PROCID', 'MPI_RANKID', 'MP_CHILD', 'MPIRUN_RANK']

    global rank_id
    for v in rank_env_vars:
        if os.getenv(v):
            rank_id = int(os.getenv(v))
            break

#-------------------------------------------------------------------------------
# Generate a gdb command file.
#-------------------------------------------------------------------------------

def gen_cmd_file(cmds):
    """
    Generate a gdb command file.
    """

    # Build file name

    if rank_id > -1:
        f_name = "./commands_r" + str(rank_id) + ".gdb"
    else:
        f_name = "./commands.gdb"
    f = open(f_name, "w")
    f.write("set breakpoint pending on\n")
    for c in cmds:
        f.write(c + '\n')
    f.close()

    return f_name

#-------------------------------------------------------------------------------
# Run gdb-based or compatible debugger with minimal options
#-------------------------------------------------------------------------------

def run_minimal_debug(path, args=None,
                      debugger='gdb', debugger_opts=None,
                      debugger_ui='terminal'):
    """
    Run gdb-based or compatible debugger with minimal options.
    This is useful for debuggers such as KDevelop and Nemiver, whose
    command lines do not allow passing GDB options or command files.
    """

    # Start building debugger command options

    cmd = [debugger]

    gdb = 'gdb'
    if debugger_opts:
        for o in debugger_opts:
            if o.find('--back-end=') == 0: # Specify back-end
                gdb = o[o.find('=')+1:]

    if debugger_ui == 'kdevelop':
        cmd += ['--debug', gdb]
        cmd += [path]
        if args:
            cmd += args

    elif debugger_ui == 'kdbg':
        if args:
            cmd += ['-a'] + args
        cmd += [path]

    elif debugger_ui == 'gede':
        cmd += ['--args', path]
        if args:
            cmd += args

    elif debugger_ui == 'nemiver':
        if gdb != 'gdb':
            cmd += ['--gdb-binary=' + gdb]
        cmd += [path]
        if args:
            cmd += args

    else:
        cmd += [path]
        if args:
            cmd += args

    return subprocess.call(cmd)

#-------------------------------------------------------------------------------
# Find available port (for gdbgui)
#-------------------------------------------------------------------------------

def find_free_port():
    import socket
    tcp = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    tcp.bind(('', 0))
    addr, port = tcp.getsockname()
    tcp.close()
    return port

#-------------------------------------------------------------------------------
# Run gdb-based or compatible debugger.
#-------------------------------------------------------------------------------

def run_gdb_debug(path, args=None, gdb_cmds=None,
                  debugger='gdb', debugger_opts=None, debugger_ui='terminal'):
    """
    Run gdb-based or compatible debugger.
    """

    # Combine given (gdb_cmds) and needed gdb commands in "cmds" list

    cmds = []

    target = False
    if gdb_cmds:
        if gdb_cmds[0].find('target') == 0:
            target = True
            cmds.append(gdb_cmds[0])

    # Check if we already have a command file

    cmds.append('b main')

    if not target:
        run_cmd = "run"
        if args:
            for a in args:
                if a.find(' ') > -1:
                    run_cmd += ' "' + a + '"'
                else:
                    run_cmd += ' ' + a
        cmds.append(run_cmd)

    if gdb_cmds:
        for c in gdb_cmds:
            if c.find('target') < 0:
                cmds.append(c)

    if target:
        cmds.append('continue')

    cmds.append('b exit')

    # Start building debugger command options

    cmd = []

    # Add debugger options, with a specific handling of
    # existing debugger command files

    gdb = 'gdb'

    debugger_command_file = None

    term = os.getenv('TERM')
    if not term:
        term = 'xterm'
    elif term[:5] == 'xterm':
        term = 'xterm'

    if debugger_opts:
        file_next = False
        for o in debugger_opts:
            if o == '-x':
                file_next = True
            elif file_next:
                debugger_command_file = cmd
                file_next = False
            elif o.find('--back-end=') == 0: # Specify back-end
                gdb = o[o.find('=')+1:]
            elif o.find('--breakpoints=') == 0: # Specify breakpoints
                for bp in o[o.find('=')+1:].split(','):
                    cmds.append('b ' + bp)
            elif o == '--asan-bp': # gcc Adress sanitizer breakpoints
                cmds.append('b __asan::ReportGenericError')
            elif o.find('--terminal=') == 0: # Specify terminal
                term = o[o.find('=')+1:]
            else:
                cmd.append(o)

    f = gen_cmd_file(cmds)

    cmd.append('-x')
    cmd.append(f)

    if debugger_command_file:
        fp = open(debugger_command_file)
        lines = f.readlines()
        fp.close()
        fp = open(f, 'a')
        for l in lines:
            fp.write(l)
        fp.close()

    # Finalize command

    cmd.append(path)

    # Start building command to run

    if debugger_ui in ['terminal', 'cgdb']:
        cmd_string = str(debugger)
        for c in cmd:
            cmd_string += ' ' + enquote_arg(str(c))
        cmd = [term, '-e', cmd_string]

    elif debugger_ui == 'gdbgui':
        cmd.insert(0, debugger)
        if gdb != 'gdb':
            cmd.insert(1, gdb)
            cmd.insert(1, '--gdb')
        p = find_free_port()
        cmd.insert(1, str(p))
        cmd.insert(1, '-p')

    elif debugger_ui == 'ddd':
        cmd.insert(0, debugger)
        if gdb != 'gdb':
            cmd.insert(1, gdb)
            cmd.insert(1, '--debugger')

    elif debugger_ui == 'emacs' or debugger_ui == 'emacs23':
        cmd_string = r'"(gdb \"' + gdb
        if debugger_ui == 'emacs23': # emacs 23 or older
            cmd_string += ' --annotate=3'
        else:
            cmd_string += ' -i=mi'   # emacs 24 and newer
        for c in cmd:
            if cmd != gdb:
                cmd_string += ' ' + enquote_arg(str(c))
        cmd_string += r'\")"'
        cmd = [debugger, '--eval', cmd_string]

    if not cmd[0]:
        cmd[0] = 'gdb'

    if not debugger_ui in ("emacs", "emacs23"):
        p = subprocess.Popen(cmd)
    else:
        cmd_line = cmd[0]
        for c in cmd[1:]:
            cmd_line += ' ' + c
        p = subprocess.Popen(cmd_line, shell=True)

    return p

#-------------------------------------------------------------------------------
# Run Valgrind gdb server and debugger.
#-------------------------------------------------------------------------------

def run_vgdb_debug(path,
                   args = None,
                   valgrind = 'valgrind',
                   valgrind_opts = None,
                   debugger = 'gdb',
                   debugger_opts=None,
                   debugger_ui = 'terminal'):
    """
    Run Valgrind gdb server and debugger.
    """

    kwargs = {}
    kwargs['stderr'] = subprocess.PIPE

    cmds = []

    cmd = [valgrind]
    if valgrind_opts:
       cmd += valgrind_opts
    cmd += [path]
    if args:
       cmd += args
    p0 = subprocess.Popen(cmd, universal_newlines=True, **kwargs)

    cmd = None
    p1 = None

    vgdb_pid = None

    while p0.poll() == None:
        output = p0.stderr.readline()
        if not cmd:
            idx = output.find("target remote")
            if idx > -1:
                cmd = output[idx:].rstrip()
                try:
                    vgdb_pid = int(cmd[cmd.index("--pid")+6:])
                except Exception:
                    pass
                cmds.insert(0, cmd)
                p1 = run_gdb_debug(path, args, cmds,
                                   debugger, debugger_opts,
                                   debugger_ui)
        print(output.rstrip())
        if p1:
            if p1.poll() != None:
                break

    if p1:
        if p1.returncode != None: # Debugger has finished/exited
            p0.poll()
            if p0.returncode == None:
                if vgdb_pid:
                    # make sure vgdb is from the same path as Valgrind
                    if os.path.isabs(valgrind):
                        vgdb = os.path.join(os.path.dirname(valgrind), "vgdb")
                    else:
                        vgdb = "vgdb"
                    subprocess.call([vgdb, "--pid="+str(vgdb_pid), "v.kill"])
                else:
                    p0.kill()

    p0.communicate()

    if p1:
        p1.communicate()

    return p0.returncode

#-------------------------------------------------------------------------------
# Run program under Valgrind.
#-------------------------------------------------------------------------------

def run_valgrind(path, args=None, valgrind='valgrind', valgrind_opts=None,
                 debugger_opts=None):
    """
    Run program under Valgrind.
    """

    # Start building valgrind command options

    cmd = [valgrind]

    # Add valgrind options, with a specific handling of
    # existing debugger command files

    debugger_attach = False

    if valgrind_opts:
        cmd += valgrind_opts
        for o in valgrind_opts:
            if o[0:5] == '--db-':
                debugger_attach = True
                break

    cmd += [path]

    if args:
        cmd += args

    # Start building command to run

    if debugger_attach:
        if rank_id < 0:
            return subprocess.call(cmd)
        else:
            term = os.getenv('TERM')
            if not term:
                term = xterm
            if debugger_opts:
                for o in debugger_opts:
                    if o.find('--terminal=') == 0: # Specify terminal
                        term = o[o.find('=')+1:]

            cmd_line = term + ' -e "' + valgrind
            for c in cmd[1:]:
                cmd_line += ' ' + enquote_arg(str(c))
            cmd_line += ' ; read"'

            return subprocess.call(cmd_line, shell=True)

    else:
        if rank_id > 0:
            vg_f = open('valgrind.out.r' + str(rank_id), 'w')
            returncode = subprocess.call(cmd,
                                         stdout=vg_f,
                                         stderr=subprocess.STDOUT)
            vg_f.close()
            return returncode
        else:
            return subprocess.call(cmd)

#-------------------------------------------------------------------------------
# Run debugger
#-------------------------------------------------------------------------------

def run_debug(cmds):
    """
    Run debugger.
    """

    # Initializations

    cmd = []

    debugger_type = 'gdb'
    vgdb = False
    need_terminal = False

    # Tests for Valgrind

    valgrind = None

    if 'valgrind' in cmds.keys():

        valgrind = cmds['valgrind'][0]
        for cmd in cmds['valgrind'][1:]:
            if cmd[0:6] == '--vgdb':
                vgdb = True
            else:
                if cmd.find("--db-attach") > -1:
                    need_terminal = True
        for cmd in cmds['valgrind'][1:]:
            if cmd == '--vgdb=off':
                vgdb = False

    # Tests for debugger choice

    debugger_ui = 'terminal'

    if 'debugger' in cmds.keys():
        debugger = cmds['debugger'][0]
        dbg_name = os.path.basename(debugger)
        if dbg_name in debuggers.keys() and dbg_name != 'gdb':
            debugger_ui = dbg_name
        elif dbg_name.find("emacs23") > -1:
            debugger_ui = 'emacs23'
        elif dbg_name.find("emacs") > -1:
            debugger_ui = 'emacs'
        commands = None
        for cmd in cmds['debugger'][1:]:
            if debugger_type == 'gdb' and cmd == '--tui':
                need_terminal = True

    if need_terminal:
        debugger_ui = 'terminal'

    # Now run appropriate tool:

    if vgdb:
        return run_vgdb_debug(path = cmds['program'][0],
                              args = cmds['program'][1:],
                              valgrind = cmds['valgrind'][0],
                              valgrind_opts = cmds['valgrind'][1:],
                              debugger = cmds['debugger'][0],
                              debugger_opts = cmds['debugger'][1:],
                              debugger_ui = debugger_ui)

    elif 'debugger' in cmds.keys():
        if debugger in ['kdbg', 'kdevelop', 'gede', 'nemiver']:
            return run_minimal_debug(path = cmds['program'][0],
                                     args = cmds['program'][1:],
                                     debugger = cmds['debugger'][0],
                                     debugger_opts = cmds['debugger'][1:],
                                     debugger_ui = debugger_ui)
        else:
            p = run_gdb_debug(path = cmds['program'][0],
                              args = cmds['program'][1:],
                              gdb_cmds = None,
                              debugger = cmds['debugger'][0],
                              debugger_opts = cmds['debugger'][1:],
                              debugger_ui = debugger_ui)
            p.communicate()
            return p.returncode

    elif 'valgrind' in cmds.keys():
        debugger_opts = None
        if 'debugger' in cmds.keys():
            debugger_opts = cmds['debugger'][1:]
        return run_valgrind(path = cmds['program'][0],
                            args = cmds['program'][1:],
                            valgrind = cmds['valgrind'][0],
                            valgrind_opts = cmds['valgrind'][1:],
                            debugger_opts = debugger_opts)

    # We should not reach this area.

    else:
        return 1

#===============================================================================
# Run the calculation
#===============================================================================

def main(argv, pkg=None):
    """
    Main function.
    """

    cmds = process_cmd_line(argv, pkg)

    if not cmds:
        return 1

    # When command includes mpiexec, re-run script under MPI

    if 'mpiexec' in cmds.keys():
        cmd_mpi = cmds['mpiexec']
        cmd_mpi.append(sys.argv[0])
        for k in ['debugger', 'valgrind', 'program']:
            if k in cmds:
                cmd_mpi += cmds[k]
        return subprocess.call(cmd_mpi)

    else:
        init_rank_id()
        return run_debug(cmds)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    retval = main(sys.argv[1:])

    sys.exit(retval)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
