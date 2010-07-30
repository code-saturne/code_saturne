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


"""
This module describes the script used to create a study/case for Code_Saturne.

This module defines the following functions:
- usage
- process_cmd_line
- make_executable
- comments
and the following classes:
- study
- class
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, pwd, shutil, stat
import types, string, re, fnmatch
from optparse import OptionParser
import ConfigParser

import cs_config


#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------


def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-s", "--study", dest="study_name", type="string",
                      metavar="<study>",
                      help="create a new study")

    parser.add_option("-c", "--case", dest="cases_name", type="string",
                      metavar="<case>", action="append",
                      help="create a new case")

    parser.add_option("--copy-from", dest="copy", type="string",
                      metavar="<case>",
                      help="create a case from another one")

    parser.add_option("--nogui", dest="use_gui",
                      action="store_false",
                      help="don't use the GUI")

    parser.add_option("--noref", dest="use_ref",
                      action="store_false",
                      help="don't copy references")

    parser.add_option("-q", "--quiet",
                      action="store_const", const=0, dest="verbose",
                      help="do not output any information")

    parser.add_option("-v", "--verbose",
                      action="store_const", const=2, dest="verbose",
                      help="dump study creation parameters")

    parser.add_option("--nsat", dest="n_sat", type="int",
                      metavar="<nsat>",
                      help="specify the number of Code_Saturne instances")

    parser.add_option("--nsyr", dest="n_syr", type="int",
                      metavar="<nsyr>",
                      help="specify the number of SYRTHES instances")

    parser.add_option("--nsyr4", dest="n_syr4", type="int",
                      metavar="<nsyr4>",
                      help="specify the number of SYRTHES 4 instances")

    parser.set_defaults(use_gui=True)
    parser.set_defaults(use_ref=True)
    parser.set_defaults(study_name=os.path.basename(os.getcwd()))
    parser.set_defaults(cases_name=[])
    parser.set_defaults(copy=None)
    parser.set_defaults(verbose=1)
    parser.set_defaults(n_sat=1)
    parser.set_defaults(n_syr=0)
    parser.set_defaults(n_syr4=0)

    (options, args) = parser.parse_args(argv)

    if options.cases_name == []:
        if len(args) > 0:
            options.cases_name = args
        else:
            options.cases_name = ["CASE1"]

    return Study(options.study_name,
                 options.cases_name,
                 options.copy,
                 options.use_gui,
                 options.use_ref,
                 options.verbose,
                 options.n_sat,
                 options.n_syr,
                 options.n_syr4)


#-------------------------------------------------------------------------------
# Give executable mode (chmod +x) to a file
#-------------------------------------------------------------------------------


def make_executable(filename):
    """
    Give executable permission to a given file.
    Equivalent to `chmod +x` shell function.
    """

    st   = os.stat(filename)
    mode = st[stat.ST_MODE]
    os.chmod(filename, mode | stat.S_IEXEC)

    return


#-------------------------------------------------------------------------------
# Comment or uncomment examples in user files
#-------------------------------------------------------------------------------


def comments(filename, use_gui):
    """
    Comment or uncomment examples in user files.
    """

    fd = file(filename, 'r')
    fdt = file(filename+'.tmp','w')

    kwd_beg = re.compile('EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START')
    kwd_end = re.compile('EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END')


    if use_gui:

        comment_line = 0
        for line in fd:
            if kwd_beg.search(line): comment_line = 1
            if comment_line == 0:
                fdt.write(line)
            else:
                if len(line) > 1:
                    fdt.write('!ex '+line)
                else:
                    fdt.write('!ex'+line)
                if kwd_end.search(line): comment_line = 0

    else:

        for line in fd:
            if not kwd_beg.search(line) and not kwd_end.search(line):
                fdt.write(line)

    fd.close()
    fdt.close()

    shutil.move(filename+'.tmp', filename)

    return


#-------------------------------------------------------------------------------
# Definition of a class for a study
#-------------------------------------------------------------------------------

class Study:


    def __init__(self, name, cases, copy, use_gui, use_ref,
                 verbose, n_sat, n_syr, n_syr4):
        """
        Initialize the structure for a study.
        """

        self.name = name
        self.cases = []
        for c in cases:
            self.cases.append(string.upper(c))
        self.copy = copy
        if self.copy is not None:
            self.copy = os.path.abspath(self.copy)
        self.use_gui = use_gui
        self.use_ref = use_ref
        self.verbose = verbose
        self.n_sat = n_sat
        self.n_syr = n_syr
        self.n_syr4 = n_syr4


    def create(self):
        """
        Create a study.
        """

        if self.name != os.path.basename(os.getcwd()):

            if self.verbose > 0:
                sys.stdout.write("  o Creating study '%s'...\n" % self.name)

            os.mkdir(self.name)
            os.chdir(self.name)
            os.mkdir('MESH')
            os.mkdir('POST')

        # Creating cases
        repbase = os.getcwd()
        for c in self.cases:
            os.chdir(repbase)
            self.createCase(c)


    def createCase(self, casename):
        """
        Create a case for a Code_Saturne study.
        """

        if self.verbose > 0:
            sys.stdout.write("  o Creating case  '%s'...\n" % casename)

        datadir = os.path.join(cs_config.dirs.pkgdatadir)
        data_distpath  = os.path.join(datadir, 'data')
        users_distpath = os.path.join(datadir, 'users')

        try:
            os.mkdir(casename)
        except:
            sys.exit(1)

        os.chdir(casename)

        # Loop for dependency on the number of instances of Code_Saturne

        for i in range(self.n_sat):

            # Data directory

            if self.n_sat == 1:
                data = 'DATA'
            else:
                data = 'DATA.%(inst)d' % { 'inst' : i+1 }

            os.mkdir(data)

            if self.use_ref:

                thch_distpath = os.path.join(data_distpath, 'thch')
                thch          = os.path.join(data, 'THCH')
                os.mkdir(thch)
                for f in ['dp_C3P', 'dp_C3PSJ', 'dp_ELE', 'dp_FCP', 'dp_FUE', 'meteo']:
                    shutil.copy(os.path.join(thch_distpath, f), thch)

            if self.use_gui:

                csguiname = 'SaturneGUI'
                csguiscript = os.path.join(datadir, csguiname)

                shutil.copy(csguiscript, data)
                make_executable(os.path.join(data, csguiname))

            # User source files directory

            if self.n_sat == 1:
                src = 'SRC'
            else:
                src = 'SRC.%(inst)d' % { 'inst' : i+1 }

            os.mkdir(src)

            if self.use_ref:

                users = os.path.join(src, 'REFERENCE')
                shutil.copytree(users_distpath, users)

                for file in ['usini1.f90','usalin.f90']:
                    f = os.path.join(users, 'base', file)
                    comments(f, self.use_gui)

            # Copy data and source files from another case

            if self.copy is not None:

                ref_data = os.path.join(self.copy, data)
                data_files = os.listdir(ref_data)
                for f in data_files:
                    abs_f = os.path.join(ref_data, f)
                    if os.path.isfile(abs_f) and \
                            f not in ['SaturneGUI', 'preprocessor_output']:
                        shutil.copy(os.path.join(ref_data, abs_f), data)


                ref_src = os.path.join(self.copy, src)
                src_files = os.listdir(ref_src)

                c_files = fnmatch.filter(src_files, '*.c')
                h_files = fnmatch.filter(src_files, '*.h')
                f_files = fnmatch.filter(src_files, '*.[fF]90')

                for f in c_files + h_files + f_files:
                    shutil.copy(os.path.join(ref_src, f), src)


        # Loop for dependency on the number of instances of SYRTHES

        for i in range(self.n_syr):

            # Data directory

            if self.n_syr == 1:
                data_syr = 'DATA_SYR'
            else:
                data_syr = 'DATA_SYR.%(inst)d' % { 'inst' : i+1 }

            os.mkdir(data_syr)

            if self.use_ref:
                data_ref_syr = os.path.join(data_syr, 'REFERENCE')
                shutil.copytree(os.path.join(cs_config.dirs.syrthes_prefix,
                                             'data'),
                                data_ref_syr)

            # User source files directory

            if self.n_syr == 1:
                src_syr = 'SRC_SYR'
            else:
                src_syr = 'SRC_SYR.%(inst)d' % { 'inst' : i+1 }

            os.mkdir(src_syr)

            if self.use_ref:
                users_syr = os.path.join(src_syr, 'REFERENCE')
                shutil.copytree(os.path.join(cs_config.dirs.syrthes_prefix,
                                             'usr'),
                                users_syr)

        # Loop for dependency on the number of instances of SYRTHES 4

        if self.n_syr4 > 0:
            config = ConfigParser.ConfigParser()
            config.read([os.path.expanduser('~/.code_saturne.cfg'),
                         os.path.join(cs_config.dirs.sysconfdir,
                                      'code_saturne.cfg')])
            if config.has_option('install', 'syrthes'):
                syr4home = config.get('install', 'syrthes')
            else:
                sys.stderr.write("Cannot localize SYRTHES installation")
                sys.exit(1)

            for i in range(self.n_syr4):
        
                # Data directory

                if self.n_syr == 1:
                    data_syr = 'DATA_SYR'
                else:
                    data_syr = 'DATA_SYR.%(inst)d' % { 'inst' : i+1 }

                os.mkdir(data_syr)

                # Reference files path

                ref_syr_data = os.path.join(syr4home, 'data')

                if self.use_ref:
                    data_ref_syr = os.path.join(data_syr, 'REFERENCE')
                    os.mkdir(data_ref_syr)
                    files = ["syrthes.data", "syrthes.py"]
                    for f in files:
                        src = os.path.join(ref_syr_data, f)
                        dst = os.path.join(data_ref_syr, f)
                        shutil.copy2(src, dst)

                # User source files directory

                if self.n_syr == 1:
                    src_syr = 'SRC_SYR'
                else:
                    src_syr = 'SRC_SYR.%(inst)d' % { 'inst' : i+1 }

                os.mkdir(src_syr)

                # Always copy Makefile

                files = ["Makefile"]
                for f in files:
                    src = os.path.join(syr4home, f)
                    dst = os.path.join(src_syr, f)
                    shutil.copy2(src, dst)

                ref_syr_user = os.path.join(syr4home, 'usr')

                if self.use_ref:
                    users_syr = os.path.join(src_syr, 'REFERENCE')
                    shutil.copytree(ref_syr_user, users_syr)

        # Results directory (only one for all instances)

        resu = 'RESU'
        os.mkdir(resu)

        # Script directory (only one for all instances)

        scripts = 'SCRIPTS'
        os.mkdir(scripts)

        shutil.copy(os.path.join(datadir, 'runcase'), scripts)
        runcase = os.path.join(scripts, 'runcase')

        kwd = re.compile('CASEDIRNAME')

        repbase = os.getcwd()

        runcase_tmp = runcase + '.tmp'

        fd  = open(runcase, 'r')
        fdt = open(runcase_tmp,'w')

        for line in fd:
            line = re.sub(kwd, repbase, line)
            fdt.write(line)

        fd.close()
        fdt.close()

        shutil.move(runcase_tmp, runcase)
        make_executable(runcase)

        # On clusters, also copy the batch card (if defined)

        config = ConfigParser.ConfigParser()
        config.read([os.path.expanduser('~/.code_saturne.cfg'),
                     os.path.join(cs_config.dirs.sysconfdir,
                                  'code_saturne.cfg')])

        if config.has_option('install', 'batch'):

            batchfile = 'batch.' + config.get('install', 'batch')

            shutil.copy(os.path.join(datadir, 'batch', batchfile),
                        os.path.join(scripts, 'batch'))

            kwd = re.compile('nameandcase')

            studycasename = string.lower(self.name) + string.lower(casename)
            # In the cluster, names are limited to 15 caracters
            studycasename = studycasename[:15]

            batchfile = os.path.join(scripts, 'batch')
            batchfile_tmp = batchfile + '.tmp'

            fd  = open(batchfile, 'r')
            fdt = open(batchfile_tmp,'w')

            for line in fd:
                line = re.sub(kwd, repbase, line)
                fdt.write(line)

            fd.close()
            fdt.close()

            shutil.move(batchfile_tmp, batchfile)


    def dump(self):
        """
        Dump the structure of a study.
        """

        print()
        print("Name  of the study:", self.name)
        print("Names of the cases:", self.cases)
        if self.copy is not None:
            print("Copy from case:", self.copy)
        print("Use the GUI:", self.use_gui)
        print("Copy references:", self.use_ref)
        if self.n_sat > 1:
            print("Number of Code_Saturne instances:", self.n_sat)
        if self.n_syr > 0:
            print("Number of SYRTHES instances:", self.n_syr)
        if self.n_syr4 > 0:
            print("Number of SYRTHES 4 instances:", self.n_syr4)
        print()


#-------------------------------------------------------------------------------
# Creation of the study directory
#-------------------------------------------------------------------------------

def main(argv):
    """
    Main function.
    """

    welcome = """\
Code_Saturne %s study/case generation
"""

    study = process_cmd_line(argv)

    if study.verbose > 0:
        sys.stdout.write(welcome % cs_config.package.version)

    study.create()

    if study.verbose > 1:
        study.dump()


if __name__ == "__main__":
    main(sys.argv[1:])


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
