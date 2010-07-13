#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Solver.
#
#   Copyright (C) 2009  EDF
#
#   The Code_Saturne Preprocessor is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Preprocessor is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import optparse
import os, re
import sys, shutil
import string
import tempfile

import cs_config

from cs_exec_environment import run_command

#-------------------------------------------------------------------------------

def subst_name(s):

    s_up = string.upper(s)

    stdgeomfile = 'chr.geo'
    stdcasefile = 'CHR.case'

    geomfile = s + '.geo'
    casefile = s_up + '.case'

    fd = file(stdcasefile, 'r')
    fdt = file(casefile, 'w')
    kwd = re.compile('chr.geo')
    for line in fd:
        line = re.sub(kwd, s+'.geo', line)
        fdt.write(line)
    fd.close()
    fdt.close()

    os.remove(stdcasefile)
    shutil.move(stdgeomfile, geomfile)


#-------------------------------------------------------------------------------

def run_check(opts):
    """
    Run Code_Saturne preprocessor and solver.
    """
    retval = 0

    cur_dir = os.getcwd()

    cmd = os.path.join(cs_config.dirs.bindir, 'cs_preprocess')

    for o in ('-h', '--help'):
        if o in opts:
            cmd = cmd + " " + str(o)
            retval = run_command(cmd)
            return retval

    for o in opts:
        cmd = cmd + " " + str(o)
    cmd = cmd + " --ensight --case check_mesh"
    retval = run_command(cmd)

    os.chdir('check_mesh.ensight')
    subst_name('preprocessor')
    os.chdir(cur_dir)


    retval = 0

    cmd = os.path.join(cs_config.dirs.bindir, 'cs_solver')
    cmd = cmd + " --quality --log 0"
    retval = run_command(cmd)

    dir_files = os.listdir('chr.ensight')
    for f in dir_files:
        shutil.move(os.path.join('chr.ensight',f), 'check_mesh.ensight')
    os.rmdir('chr.ensight')

    os.chdir('check_mesh.ensight')
    subst_name('quality')
    os.chdir(cur_dir)


    os.remove('preprocessor_output')

    return retval

#-------------------------------------------------------------------------------

def main(argv):
    """
    Main function.
    """

    retcode = run_check(argv)
    sys.exit(retcode)


if __name__ == '__main__':
    main(sys.argv[1:])
