#!/usr/bin/env python
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

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import os, sys
from optparse import OptionParser


#-------------------------------------------------------------------------------
# Process command line
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string the
             string is split to create a list of arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-p", "--param", dest="param", type="string",
                      metavar="<param>",
                      help="give the parameters file")

    parser.add_option("-s", "--source", dest="src_dir", type="string",
                      metavar="<src_dir>",
                      help="choose source file directory")

    parser.add_option("-n", "--nproc", dest="nproc", type="int",
                      metavar="<nproc>",
                      help="specify the number of processors")

    parser.set_defaults(param="")
    parser.set_defaults(src_dir=os.getcwd())
    parser.set_defaults(nproc=1)

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)

    param = options.param
    if param != "":
        param = os.path.expanduser(options.param)
        param = os.path.expandvars(param)
        param = os.path.abspath(param)
        if not os.path.isfile(param):
            sys.stderr.write("Error: cannot access parameter file %s\n" % param)
            sys.exit(1)

    src_dir = os.path.expanduser(options.src_dir)
    src_dir = os.path.expandvars(src_dir)
    src_dir = os.path.abspath(src_dir)
    if not os.path.isdir(src_dir):
        sys.stderr.write("Error: %s is not a directory\n" % src_dir)
        sys.exit(1)

    return options.nproc, param, src_dir


#-------------------------------------------------------------------------------
# Checking user files consistency
#-------------------------------------------------------------------------------

def check_consistency(param, src_dir, n_procs):
    """
    Return 0 if parameter file options and user subroutines are consistent,
    or 1 if they are incompatible.
    """

    # List of the different available modules in Code_Saturne
    modules = ['base', 'lagr', 'rayt', 'cplv', 'fuel', 'c3pt', 'cebu',
               'clwc', 'elec', 'cfbl', 'atmo', 'ctwr']

    # Dictionnary of module's name
    moduleName = {
        'base':'standard',
        'lagr':'lagrangian',
        'rayt':'radiative transfer',
        'cplv':'pulverized coal combustion',
        'fuel':'heavy-fuel oil combustion',
        'c3pt':'3-point chemistry combustion',
        'cebu':'EBU combustion',
        'clwc':'LWC combustion',
        'elec':'electric arcs',
        'cfbl':'compressible',
        'atmo':'atmospheric',
        'ctwr':'cooling towers'}

    # Dictionnary for consistancy check (on boundary conditions definition)
    moduleCheck = {
        'base':False, 'lagr':False, 'rayt':False,
        'cplv':True,  'fuel':True,  'c3pt':True,  'cebu':True,  'clwc':True,
        'elec':True,  'cfbl':True,  'atmo':True,  'ctwr':True}

    # Dictionnary of the correspondance between a module and a specific user file
    moduleFile = {
        'base':'usclim', 'lagr':'uslag2', 'rayt':'usray2', 'cplv':'uscpcl',
        'fuel':'usfucl', 'c3pt':'usd3pc', 'cebu':'usebuc', 'clwc':'uslwcc',
        'elec':'uselcl', 'cfbl':'uscfcl', 'atmo':'usatcl', 'ctwr':'usctcl'}

    # Dictionnary to know if a module is used
    moduleUse = {
        'base':True,  'lagr':False, 'rayt':False, 'cplv':False, 'fuel':False,
        'c3pt':False, 'cebu':False, 'clwc':False, 'elec':False, 'cfbl':False,
        'atmo':False, 'ctwr':False}


    # Function returning a boolean according to the presence of the file given
    # in argument
    def isPresent(file):
        filename = file + '.f90'
        return os.path.isfile(os.path.join(src_dir, filename))


    # Consistancy test, following the definition of the boundary conditions
    errorMsg = """
     --ERROR --
     When %(f1)s is used, %(f2)s must not be
      (%(mod)s module)
      Boundary conditions are defined in %(f1)s.

    """

    if isPresent(moduleFile['base']):
        for mod in modules:
            if moduleCheck[mod] and isPresent(moduleFile[mod]):
                sys.stderr.write(errorMsg % {'f1':moduleFile[mod],
                                             'f2':moduleFile['base'],
                                             'mod':moduleName[mod]})
                return 1


    # Test on the module used (standard module is not considered here)
    moduleMsg = """
     Use the %(mod)s module

    """

    for mod in modules[1:]:
        if isPresent(moduleFile[mod]):
            sys.stdout.write(moduleMsg % {'mod':moduleName[mod]})
            moduleUse[mod] = True


    # Check if current module is incompatible with parallel runs
    errorMsg = """
     The %(mod)s module is incompatible with
     parallel runs as of the current version.

    """

    if n_procs != None:
        if n_procs > 1:
            for mod in ['ctwr','lagr']:
                if moduleUse[mod]:
                    sys.stderr.write(errorMsg % {'mod':moduleName[mod]})
                    return 1

    return 0

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    n_procs, param, src_dir = process_cmd_line(argv, pkg)

    ret_val = check_consistency(param, src_dir, n_procs)
    sys.exit(ret_val)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
