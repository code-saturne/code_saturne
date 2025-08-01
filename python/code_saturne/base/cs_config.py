# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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

import os
import sys
from optparse import OptionParser

#-------------------------------------------------------------------------------

# Prerequisites libraries
#------------------------

class prerequisite:

    def __init__(self, name, key, config_dict, have = None,
                 add_rpath = True,
                 prefix = None,
                 bindir = None, includedir = None, libdir = None,
                 flags = None):

        # Library name
        self.name = name

        if config_dict and key:
            d = config_dict.get(key, {})
        else:
            d = {}

        # Have

        false_t = ('False', 'false', 'no', '0')
        true_t = ('True', 'true', 'yes', '1')

        if have != None:
            self.have = have
        else:
            self.have = d.get('have', have)
            if self.have in false_t:
                self.have = False
            elif self.have in true_t:
                self.have = True

        # Loaded dynamically on demand (linked separately)
        self.dynamic_load = d.get('dynamic_load', False)
        if self.dynamic_load in false_t:
            self.dynamic_load = False
        elif self.dynamic_load in true_t:
            self.dynamic_load = True

        # Add in rpath
        self.add_rpath = d.get('add_rpath', None)
        if self.add_rpath == None:
            self.add_rpath = add_rpath

        # Library variant
        self.variant = d.get('variant', None)

        # Library installation directories
        for k in ('prefix', 'bindir', 'libdir'):
            setattr(self, k, d.get(k, ''))

        # Library build (dictionnary {cppflags, ldflags, libs})
        if flags != None:
            self.flags = flags
        else:
            self.flags = {}
            if self.have != 'no':
                for k in ('cppflags', 'ldflags', 'libs'):
                    self.flags[k] = d.get(k, '')

                ld_add_path = d.get('ld_add_path', None)
                if ld_add_path != None:
                    self.flags['ld_add_path'] = ld_add_path

                pythonpath = d.get('pythonpath', None)
                if pythonpath != None:
                    self.flags['pythonpath'] = pythonpath

    def print_config(self):

        print("Prequisite: " + self.name)


# Configuration info
#-------------------

class config:

    def __init__(self, config_dict):
        """
        constructor for configuration class, using a dictionnary
        read from a configuration file.
        """

        # List of optional features

        self.optfeatures = ['debug', 'relocatable', 'shared',
                            'gui', 'frontend',
                            'mpi', 'openmp', 'socket',
                            'long-gnum', 'host']
        self.features = {}

        # List of mandatory and optionnal libraries to link with
        # The order is important so as to have a coherent link command

        self.deplibs = ['saturne',                      # Code_Saturne
                        'ple',                          # PLE
                        'eos', 'coolprop',              # Equations of state
                        'ccm', 'cgns', 'med', 'hdf5',   # Mesh filters
                        'catalyst', 'melissa',          # co-processing libraries
                        'medcoupling',                  # MED coupling
                        'mumps', 'cudss',               # Sparse direct solver
                        'amgx', 'hypre', 'petsc',       # Linear algebra
                        'metis', 'scotch',              # Partionning libraries
                        'mpi',                          # MPI
                        'cuda',                         # CUDA
                        'blas',                         # BLAS (benchmark use)
                        'system']                       # User & system libraries

        # Compilers, flags and special commands

        d = config_dict.get('compilers', {})

        self.compilers = {}
        for k in ('cc', 'cxx', 'fc', 'nvcc', 'ld'):
            self.compilers[k] = d.get(k, None)

        self.flags = {}
        for k in ('cflags', 'cxxflags', 'fcflags', 'nvccflags', 'nvccflags_cpp',
                  'cflags_shared', 'cxxflags_shared', 'fcflags_shared'):
            self.flags[k] = d.get(k, '')

        self.fcmodinclude = d.get('fcmodinclude', '')
        self.rpath = d.get('rpath', '')
        self.special_user_link = d.get('cs_special_user_link', '')
        self.ld_default_search_path = d.get('ld_default_search_path', '')

        system_flags = {'cppflags': d.get('cppflags', ''),
                        'ldflags': d.get('ldflags', ''),
                        'ldflags_shared': d.get('ldflags_shared', ''),
                        'libs': d.get('libs', '')}

        # Constants for system-dependant file extensions
        if sys.platform.startswith("win"):
            self.cfgext = ".ini"
            self.exeext = ".exe"
            self.shext = ".bat"
        else:
            self.cfgext = ".cfg"
            self.exeext = ""
            self.shext = ""

        # Known module names (code_saturne always present; others may
        # depend on availability, but names are reserved here);
        # Must be all lowercase to simplify tests across scripts.

        self.solver_modules = {}

        self.solver_modules['code_saturne'] \
            = {'solver': 'cs_solver' + self.exeext}

        self.solver_modules['neptune_cfd'] \
            = {'solver': 'nc_solver' + self.exeext}

        # Headers and libraries to add for specific executable/module names

        self.exec_include = {'nc_solver' + self.exeext: "neptune_cfd"}

        self.exec_libs = {'cs_solver' + self.exeext: ["-lcs_solver"],
                          'nc_solver' + self.exeext: ["-lnc_solver",
                                                      "-lneptune"]}

        # Python-related information

        d = config_dict.get('python', {})

        self.python = d.get('python', None)
        self.pyuic = d.get('pyuic5', None)
        self.pyrcc = d.get('pyurcc5', None)

        # Execution environment

        d = config_dict.get('environment', {})

        self.env_modules = d.get('env_modules', '')
        self.env_modulecmd = d.get('env_modulecmd', '')

        self.salome_env = d.get('salome_env', '')

        # Setup the optionnal features

        d = config_dict.get('features', {})

        self.features = {}
        self.features['debug'] = d.get('debug', 'no')
        self.features['relocatable'] = d.get('relocatable', 'no')
        self.features['shared'] = d.get('shared', 'yes')
        self.features['gui'] = d.get('gui', 'yes')
        self.features['frontend'] = d.get('frontend', 'yes')
        self.features['mpi'] = d.get('mpi', 'yes')
        self.features['openmp'] = d.get('openmp', 'yes')
        self.features['cuda'] = d.get('cuda', 'no')
        self.features['long-gnum'] = d.get('long-gnum', 'yes')
        self.features['build_os'] = d.get('build_os', '')

        # Now, one can setup the prerequisites information

        self.libs = {}

        # Setup code_saturne libraries
        # Here, CPPFLAGS and LDFLAGS will be provided by a get_dir method

        self.libs['saturne'] = \
            prerequisite('Code_Saturne', None, None,
                         have = True,
                         flags = {'cppflags': "",
                                  'ldflags': "",
                                  'libs': "-lsaturne"})

        # Setup PLE library
        # Here, the variant (internal or external) will be used to add
        # paths to the command line

        self.libs['ple'] = \
            prerequisite('PLE', 'ple', config_dict, have=True)

        # Setup user and system libraries

        self.libs['system'] = \
            prerequisite('System', None, None,
                         have = True,
                         flags = system_flags,
                         add_rpath = False)

        # Setup the optionnal libraries

        self.libs['blas'] = prerequisite('BLAS', 'blas', config_dict)

        self.libs['ccm']  = prerequisite('CCM', 'ccm', config_dict)
        self.libs['cgns'] = prerequisite('CGNS', 'cgns', config_dict)
        self.libs['hdf5'] = prerequisite('HDF5', 'hdf5', config_dict)
        self.libs['med']  = prerequisite('MED', 'med', config_dict)

        self.libs['catalyst'] = prerequisite('CATALYST', 'catalyst', config_dict)
        self.libs['melissa']  = prerequisite('MELISSA', 'melissa', config_dict)
        self.libs['medcoupling'] = prerequisite('MEDCOUPLING',
                                                'medcoupling', config_dict)

        self.libs['eos']       = prerequisite('EOS', 'eos', config_dict)
        self.libs['coolprop']  = prerequisite('COOLPROP', 'coolprop', config_dict)

        self.libs['mpi'] = prerequisite('MPI', 'mpi', config_dict, add_rpath=False)

        self.libs['scotch'] = prerequisite('SCOTCH', 'scotch', config_dict)
        self.libs['metis']  = prerequisite('METIS', 'metis', config_dict)

        self.libs['cuda'] =  prerequisite('CUDA', 'cuda', config_dict)

        self.libs['petsc'] = prerequisite('PETSc', 'petsc', config_dict)
        self.libs['amgx']  = prerequisite('Amgx', 'amgx', config_dict)
        self.libs['hypre'] = prerequisite('HYPRE', 'hypre', config_dict)
        self.libs['mumps'] = prerequisite('MUMPS', 'mumps', config_dict)
        self.libs['cudss'] = prerequisite('CUDSS', 'cudss', config_dict)

    def __get_search_paths_catalyst__(self):
        """
        return Catalyst library search path, Python search paths,
        and other environment variables if available
        """

        lib_dirs = []
        pythonpath_dirs = []
        env_vars = None

        catalyst_lib_dir = None

        libs = self.libs['catalyst'].flags['libs']
        for l in libs.split('-Wl,'):
            if l.find('libvtkPVPythonCatalyst') > -1:
                catalyst_lib_dir = os.path.dirname(l)
                break

        if catalyst_lib_dir and self.features['relocatable']:
            catalyst_root_dir = os.getenv('CATALYST_ROOT_DIR')
            if catalyst_root_dir:
                subdir_idx = catalyst_lib_dir.rfind('lib')
                if subdir_idx > 1:
                    catalyst_lib_dir = os.path.join(catalyst_root_dir,
                                                    catalyst_lib_dir[subdir_idx:])

        if catalyst_lib_dir:
            if self.libs['catalyst'].dynamic_load:
                lib_dirs = [catalyst_lib_dir]
            sp_dir = os.path.join(catalyst_lib_dir, 'site-packages')
            if os.path.isdir(sp_dir):
                pythonpath_dirs = [sp_dir]

        # Add possible additional Catalyst dependency paths

        catalyst_ld_add_path = os.getenv('CATALYST_LD_ADD_PATH')
        if not catalyst_ld_add_path:
            catalyst_ld_add_path = self.libs['catalyst'].flags.get('ld_add_path',
                                                                   None)

        if catalyst_ld_add_path:
            for d in catalyst_ld_add_path.split(os.pathsep):
                if d: # avoid empty values before first or after last separator
                    lib_dirs.append(d)

        # Add additional environment variables

        catalyst_plugin_path = os.getenv('CATALYST_PLUGIN_PATH')
        if not catalyst_plugin_path:
            catalyst_plugin_path = ''
        env_vars = {'PV_PLUGIN_PATH':catalyst_plugin_path}

        return lib_dirs, pythonpath_dirs, env_vars

    def get_run_environment_dependencies(self):
        """
        return library search path, Python search paths,
        and other environment variables if available or required
        """

        lib_dirs = []
        pythonpath_dirs = []
        env_vars = {}

        for lib in self.deplibs:
            if self.libs[lib].have == True:

                if lib == 'catalyst':
                    catalyst_lib_dirs, catalyst_pythonpath_dirs, catalyst_env_vars \
                        = self.__get_search_paths_catalyst__()

                    for d in catalyst_lib_dirs:
                        lib_dirs.append(d)
                    for p in catalyst_pythonpath_dirs:
                        pythonpath_dirs.append(p)
                    env_vars.update(catalyst_env_vars)

        return lib_dirs, pythonpath_dirs, env_vars

    def __get_dep_libs_path_catalyst__(self):
        """
        return Catalyst dependency path required for compilation
        """

        lib_dirs = []

        # Add possible additional Catalyst dependency paths

        catalyst_ld_add_path = os.getenv('CATALYST_LD_ADD_PATH')
        if catalyst_ld_add_path:
            for d in catalyst_ld_add_path.split(os.pathsep):
                if d: # avoid empty values before first or after last separator
                    lib_dirs.append(d)

        return lib_dirs

    def get_compile_dependency_paths(self):
        """
        return additional library search if available or required
        """

        lib_dirs = []

        for lib in self.deplibs:
            if self.libs[lib].have == True:

                if lib == 'catalyst':
                    catalyst_lib_dirs = self.__get_dep_libs_path_catalyst__()
                    for d in catalyst_lib_dirs:
                        lib_dirs.append(d)

        return lib_dirs

    def print_config(self):
        """
        Print configuration info
        """

        for lib in self.deplibs:
            self.libs[lib].print_config()

#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.

    Input Argument:
      arg -- This can be either a list of arguments as in
             sys.argv[1:] or a string that is similar to the one
             passed on the command line.  If it is a string,
             it is split to create a list of arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("--cc", dest="print_cc",
                      action="store_true",
                      help="C compiler used for build")

    parser.add_option("--cxx", dest="print_cxx",
                      action="store_true",
                      help="C++ compiler used for build")

    parser.add_option("--fc", dest="print_fc",
                      action="store_true",
                      help="Fortran compiler used for build")

    parser.add_option("--cflags", dest="print_cflags",
                      action="store_true",
                      help="C compiler flags")

    parser.add_option("--cxxflags", dest="print_cxxflags",
                      action="store_true",
                      help="C++ compiler flags")

    parser.add_option("--fcflags", dest="print_fcflags",
                      action="store_true",
                      help="Fortran compiler flags")

    parser.add_option("--rpath", dest="print_rpath",
                      action="store_true",
                      help="Linker rpath command line")

    parser.add_option("--python", dest="print_python",
                      action="store_true",
                      help="Python interpreter")

    parser.add_option("--pyuic", dest="print_pyuic",
                      action="store_true",
                      help="pyuic tool for PyQt support")

    parser.add_option("--pyrcc", dest="print_pyrcc",
                      action="store_true",
                      help="pyrcc tool for PyQt support")

    parser.add_option("--have", dest="have", metavar="<lib>",
                      help="supported feature or library")

    parser.add_option("--cppflags", dest="cppflags", metavar="<lib>",
                      help="C preprocessor flags (e.g. -D<macro>, ...)")

    parser.add_option("--ldflags", dest="ldflags", metavar="<lib>",
                      help="linker flags (e.g. -g, -L<path>, ...)")

    parser.add_option("--libs", dest="libs", metavar="<lib>",
                      help="librairies used (e.g. -l<libname>, ...)")

    parser.add_option("--pythondir", dest="print_pythondir",
                      action="store_true",
                      help="directory for the 'site-packages' subdirectory" \
                          " of the standard Python install tree")

    parser.add_option("--datarootdir", dest="print_datarootdir",
                      action="store_true",
                      help="directory where architecture-independent" \
                          + " files are installed (e.g. <prefix>/share)")

    parser.set_defaults(print_cc=False)
    parser.set_defaults(print_cxx=False)
    parser.set_defaults(print_fc=False)

    parser.set_defaults(print_cflags=False)
    parser.set_defaults(print_cxxflags=False)
    parser.set_defaults(print_fcflags=False)

    parser.set_defaults(print_rpath=False)

    parser.set_defaults(print_python=False)
    parser.set_defaults(print_pyrcc=False)
    parser.set_defaults(print_pyuic=False)

    parser.set_defaults(have=None)
    parser.set_defaults(cppflags=None)
    parser.set_defaults(ldflags=None)
    parser.set_defaults(libs=None)

    parser.set_defaults(print_pythondir=False)
    parser.set_defaults(print_datarootdir=False)

    (options, args) = parser.parse_args(argv)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)

    return options

#-------------------------------------------------------------------------------

def get_config(pkg):
    """
    Get the configuration information.
    """
    msg = """\
Compilers and associated options:
  cc = %(cc)s
  cxx = %(cxx)s
  fc = %(fc)s
  cflags = %(cflags)s
  cxxflags = %(cxxflags)s
  fcflags = %(fcflags)s
  rpath = %(rpath)s\
"""

    return msg \
        % { 'cc':pkg.config.compilers['cc'],
            'cxx': pkg.config.compilers['cxx'],
            'fc':pkg.config.compilers['fc'],
            'cflags':pkg.config.flags['cflags'],
            'cxxflags':pkg.config.flags['cxxflags'],
            'fcflags':pkg.config.flags['fcflags'],
            'rpath':pkg.config.rpath }

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main configuration info function.
    """

    opts = process_cmd_line(argv)

    cfg = pkg.config

    if opts.print_cc  == True: print(cfg.compilers['cc'])
    if opts.print_cxx == True: print(cfg.compilers['cxx'])
    if opts.print_fc  == True: print(cfg.compilers['fc'])

    if opts.print_cflags   == True: print(cfg.flags['cflags'])
    if opts.print_cxxflags == True: print(cfg.flags['cxxflags'])
    if opts.print_fcflags  == True: print(cfg.flags['fcflags'])

    if opts.print_rpath == True: print(cfg.rpath)

    if opts.print_python  == True: print(cfg.python)
    if opts.print_pyuic  == True: print(cfg.pyuic)
    if opts.print_pyrcc  == True: print(cfg.pyrcc)

    if opts.have is not None:
        if opts.have in cfg.deplibs: print(cfg.libs[opts.have].have)
        if opts.have in cfg.optfeatures: print(cfg.features[opts.have])

    if opts.cppflags is not None:
        # Specific handling of code_saturne has pkgincludedir has to be
        # correctly expended. Likewise for PLE, if internal version is used
        if opts.cppflags == "saturne":
            print("-I" + pkg.get_dir("pkgincludedir"))
        elif opts.cppflags == "ple":
            if cfg.libs['ple'].variant == "internal":
                print("-I" + pkg.get_dir("includedir"))
            else:
                print(cfg.libs[opts.cppflags].flags['cppflags'])
        else:
            print(cfg.libs[opts.cppflags].flags['cppflags'])

    if opts.ldflags is not None:
        # Specific handling of code_saturne has pkgincludedir has to be
        # correctly expended. Likewise for PLE, if internal version is used
        if opts.ldflags == "saturne":
            print("-L" + pkg.get_dir("libdir"))
        elif opts.ldflags == "ple":
            if cfg.libs['ple'].variant == "internal":
                print("-L" + pkg.get_dir("libdir"))
            else:
                print(cfg.libs[opts.cppflags].flags['ldflags'])
        else:
            if cfg.libs[opts.ldflags].dynamic_load == False:
                print(cfg.libs[opts.ldflags].flags['ldflags'])

    if opts.libs is not None:
        if cfg.libs[opts.libs].dynamic_load == False:
            print(cfg.libs[opts.libs].flags['libs'])

    if opts.print_pythondir: print(pkg.get_dir("pythondir"))
    if opts.print_datarootdir: print(pkg.get_dir("datarootdir"))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    from code_saturne.base import cs_package
    pkg = cs_package.package()
    main(sys.argv[1:], pkg)
