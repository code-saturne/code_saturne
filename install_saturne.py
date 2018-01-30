#!/usr/bin/env python

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys

if sys.version_info[:2] < (2,6):
    sys.stderr.write("This script needs Python 2.6 at least\n")

import platform

if platform.system == 'Windows':
    sys.stderr.write("This script only works on Unix-like platforms\n")

import os, shutil
import string
import subprocess
import types, string, re, fnmatch

#-------------------------------------------------------------------------------
# Global variable
#-------------------------------------------------------------------------------

verbose = 'yes'

#-------------------------------------------------------------------------------
# Global methods
#-------------------------------------------------------------------------------

def run_command(cmd, stage, app, log):
    """
    Run a command via the subprocess module.
    """

    if verbose == 'yes':
        sys.stdout.write("   o " + stage + "...\n")

    p = subprocess.Popen(cmd,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)

    output = p.communicate()
    log.write(output[0])

    if p.returncode != 0:
        sys.stderr.write("Error during " + stage.lower() +
                         " stage of " + app + ".\n")
        sys.stderr.write("See " + log.name + " for more information.\n")
        sys.exit(1)

#-------------------------------------------------------------------------------

def run_test(cmd):
    """
    Run a test for a given command via the subprocess module.
    """

    if verbose == 'yes':
        sys.stdout.write("   o Checking for " + os.path.basename(cmd) + "...  ")

    cmd = "type " + cmd

    p = subprocess.Popen(cmd,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)

    output = p.communicate()
    if verbose == 'yes':
        if p.returncode == 0: log_str = output[0].split()[2]
        else: log_str = "not found"
        sys.stdout.write("%s\n" % log_str)

    return p.returncode

#-------------------------------------------------------------------------------

def check_directory():
    """
    Check we are not in the source directory.
    """

    script_path = os.path.abspath(sys.argv[0])
    top_srcdir, script_name = os.path.split(script_path)
    abscwd = os.path.abspath(os.getcwd())

    if not os.path.relpath(abscwd, top_srcdir)[0:2] == '..':
        message = \
"""
The '%(script_name)s' installer script should not be run from inside the
'%(top_srcdir)s' Code_Saturne source directory,
but from a separate directory.

We recommend running for example:

  cd %(top_srcdir)s/..
  mkdir saturne_build
  cd saturne_build
  ../%(rel_dir_name)s/%(rel_script_path)s

or (using absolute paths):

  mkdir %(build_path)s
  cd %(build_path)s
  %(script_path)s

"""
        rel_script_path = os.path.basename(script_path)
        rel_dir_name = os.path.basename(top_srcdir)
        build_path = top_srcdir + '_build'

        sys.stdout.write(message % {'script_name': script_name,
                                    'top_srcdir': top_srcdir,
                                    'rel_dir_name': rel_dir_name,
                                    'rel_script_path': rel_script_path,
                                    'build_path': build_path,
                                    'script_path': script_path})

        sys.exit(1)

#-------------------------------------------------------------------------------

def find_executable(names, env_var=None):

    """
    Find executable in path using given names.
    """

    # If an associated environment variable id defined,
    # test it first.

    if env_var:
        k = os.environ.get(env_var)
        if k:
            return k
        else:
            for a in sys.argv[1:]:
                if a.find(env_var) == 0:
                    return a.split('=', 1)[1]

    # Otherwise, use standard path.

    p = os.getenv('PATH').split(':')

    for name in names:
        if os.path.isabs(name):
            if os.path.isfile(name):
                return name
        else:
            for d in p:
                absname = os.path.join(d, name)
                if os.path.isfile(absname):
                    if p != '/usr/bin':
                        return absname
                    else:
                        return name
    return None

#-------------------------------------------------------------------------------
# Class definition for a generic package
#-------------------------------------------------------------------------------

class Package:

    def __init__(self, name, description, package, version, archive, url):

        # Package information
        self.name = name
        self.description = description
        self.package = package
        self.version = version
        self.archive = archive
        if self.archive:
            self.url = url % self.archive
        else:
            self.url = None

        # Installation information
        self.shared = True # not modifiable yet
        self.use = 'no'
        self.installation = 'no'
        self.source_dir = None
        self.install_dir = None
        self.config_opts = ''
        self.log_file = sys.stdout
        self.cxx = None
        self.cc = None
        self.fc = None
        self.vpath_support = True
        self.create_install_dirs = False

    #---------------------------------------------------------------------------

    def set_version_from_configure(self, path):
        f = open(path)
        for l in f:
            if not self.version and l[0:15] == 'PACKAGE_VERSION':
                sep = l[16] # quote is usually ', but could be "
                try:
                    self.version = l.split(sep)[1]
                    break
                except Exception:
                    pass
        f.close()

    #---------------------------------------------------------------------------

    def info(self):

        sys.stdout.write("\n"
                         "   %(s_name)s (%(l_name)s)\n"
                         "   version: %(vers)s\n"
                         "   url: %(url)s\n"
                         "   package: %(pack)s\n"
                         "   source_dir: %(src)s\n"
                         "   install_dir: %(inst)s\n"
                         "   config_opts: %(opts)s\n\n"
                         % {'s_name':self.name, 'l_name':self.description,
                            'vers':self.version, 'url':self.url,
                            'pack':self.package,
                            'src':self.source_dir, 'inst':self.install_dir,
                            'opts':self.config_opts})

    #---------------------------------------------------------------------------

    def download(self):

        if sys.version_info[0] < (3):
            import urllib2
            u = urllib2.urlopen(self.url)
            data = u.read()
            f = open(self.archive, 'wb')
            f.write(data)
            f.close()

        else:
            import urllib.request
            urllib.request.urlretrieve(self.url, self.archive)

    #---------------------------------------------------------------------------

    def extract(self):

        if self.archive[-4:] == '.zip':

            import zipfile

            if not zipfile.is_zipfile(self.archive):
                sys.stderr.write("%s is not a zip archive\n" % self.archive)
                sys.exit(1)

            zip = zipfile.ZipFile(self.archive)

            relative_source_dir = zip.namelist()[0].split(os.path.sep)[0]
            self.source_dir = os.path.abspath(relative_source_dir)

            zip.close()

            # Use external unzip command so as to keep file properties

            p = subprocess.Popen('unzip ' + self.archive,
                                 shell=True,
                                 universal_newlines=True,
                                 stdout=sys.stdout,
                                 stderr=sys.stderr)
            output = p.communicate()

            if p.returncode != 0:
                sys.stderr.write("Error unzipping file " + self.archive + ".\n")
                sys.exit(1)

        else:

            import tarfile

            if not tarfile.is_tarfile(self.archive):
                sys.stderr.write("%s is not a tar archive\n" % self.archive)
                sys.exit(1)

            tar = tarfile.open(self.archive)

            first_member = tar.next()
            relative_source_dir = first_member.name.split(os.path.sep)[0]
            self.source_dir = os.path.abspath(relative_source_dir)

            try:
                tar.extractall()
            except AttributeError:
                for tarinfo in tar:
                    tar.extract(tarinfo)

            tar.close()

    #---------------------------------------------------------------------------

    def install(self):

        current_dir = os.getcwd()

        build_dir = self.source_dir + '.build'
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        # Create some install directories in case install script does not work
        if self.create_install_dirs and self.install_dir:
            inc_dir = os.path.join(self.install_dir, 'include')
            lib_dir = os.path.join(self.install_dir, 'lib')
            for dir in [inc_dir, lib_dir]:
                if not os.path.isdir(dir):
                    os.makedirs(dir)

        # Copy source files in build directory if VPATH feature is unsupported
        if self.vpath_support:
            os.makedirs(build_dir)
        else:
            shutil.copytree(self.source_dir, build_dir)
        os.chdir(build_dir)

        configure = os.path.join(self.source_dir, 'configure')
        if os.path.isfile(configure):

            # Set command line for configure pass

            if self.install_dir:
                configure = configure + ' --prefix=' + self.install_dir
            configure = configure + ' ' + self.config_opts

            # Add compilers
            if self.cxx: configure += ' CXX=\"' + self.cxx + '\"'
            if self.cc: configure += ' CC=\"' + self.cc + '\"'
            if self.fc: configure += ' FC=\"' + self.fc + '\"'

            # Install the package and clean build directory
            run_command(configure, "Configure", self.name, self.log_file)
            run_command("make", "Compile", self.name, self.log_file)
            run_command("make install", "Install", self.name, self.log_file)
            run_command("make clean", "Clean", self.name, self.log_file)

        elif os.path.isfile(os.path.join(self.source_dir, 'CMakeLists.txt')):

            # Set command line for CMake pass

            cmake = 'cmake'
            if self.install_dir:
                cmake += ' -DCMAKE_INSTALL_PREFIX=' + self.install_dir
            cmake += ' ' + self.config_opts

            # Add compilers
            if self.cxx: cmake += ' -DCMAKE_CXX_COMPILER=\"' + self.cxx + '\"'
            if self.cc: cmake += ' -DCMAKE_C_COMPILER=\"' + self.cc + '\"'
            if self.fc: cmake += ' -DCMAKE_Fortran_COMPILER=\"' + self.fc + '\"'

            cmake += ' ' + self.source_dir

            # Install the package and clean build directory
            run_command(cmake, "Configure", self.name, self.log_file)
            run_command("make VERBOSE=1", "Compile", self.name, self.log_file)
            run_command("make install VERBOSE=1", "Install", self.name, self.log_file)
            run_command("make clean", "Clean", self.name, self.log_file)

        # End of installation
        os.chdir(current_dir)

    #---------------------------------------------------------------------------

    def install_ptscotch(self):

        current_dir = os.getcwd()

        build_dir = self.source_dir + '.build'
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        # Create some install directories in case install script does not work
        if self.install_dir:
            inc_dir = os.path.join(self.install_dir, 'include')
            lib_dir = os.path.join(self.install_dir, 'lib')
            for dir in [inc_dir, lib_dir]:
                if not os.path.isdir(dir):
                    os.makedirs(dir)

        # Copy source files in build directory as VPATH feature is unsupported
        shutil.copytree(self.source_dir, build_dir)

        os.chdir(os.path.join(build_dir, 'src'))

        # Work around Ubuntu Metis build bug
        ldflags_add = ''
        try:
            p = subprocess.Popen(self.cc + ' -Xlinker --help',
                                 shell=True,
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            output = p.communicate()[0]
            if output.find("--no-as-needed") > -1:
                ldflags_add = ' -Wl,--no-as-needed\n'
        except Exception:
            pass

        if self.shared:
            fdr = open('Make.inc/Makefile.inc.x86-64_pc_linux2.shlib')
        else:
            fdr = open('Make.inc/Makefile.inc.x86-64_pc_linux2')
        fd = open('Makefile.inc','w')

        re_thread = re.compile('-DSCOTCH_PTHREAD')
        re_intsize32 = re.compile('-DINTSIZE32')
        re_intsize64 = re.compile('-DINTSIZE64')
        re_idxsize64 = re.compile('-DIDXSIZE64')

        for line in fdr:

            if line[0:3] in ['CCS', 'CCP', 'CCD']:
                i1 = line.find('=')
                line = line[0:i1] + '= ' + self.cc + '\n'
            line = re.sub(re_thread, '', line)
            line = re.sub(re_intsize32, '', line)
            line = re.sub(re_intsize64, '', line)
            line = re.sub(re_idxsize64, '-DIDXSIZE64 -DINTSIZE64', line)
            if ldflags_add and line[0:7] == 'LDFLAGS':
                line = line[:-1] + ldflags_add

            fd.write(line)

        fdr.close()
        fd.close()

        # Build and install
        for target in ['scotch', 'ptscotch']:
            run_command("make "+target, "Compile", self.name, self.log_file)
            run_command("make install prefix="+self.install_dir,
                        "Install", self.name, self.log_file)
            run_command("make clean", "Clean", self.name, self.log_file)

        # End of installation
        os.chdir(current_dir)

    #---------------------------------------------------------------------------

    def install_parmetis(self):

        current_dir = os.getcwd()

        build_dir = self.source_dir + '.build'
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        # Copy source files in build directory as VPATH feature is unsupported
        shutil.copytree(self.source_dir, build_dir)

        for d in [os.path.join(build_dir, 'metis'), build_dir]:

            os.chdir(d)

            configure = "make config prefix=" + self.install_dir
            configure += " cc=" + self.cc
            if self.cxx:
                configure += " cxx=" + self.cxx
            if self.shared:
                configure += " shared=1 "

            # Install the package and clean build directory
            run_command(configure, "Configure", self.name, self.log_file)
            run_command("make", "Compile", self.name, self.log_file)
            run_command("make install", "Install", self.name, self.log_file)
            run_command("make clean", "Clean", self.name, self.log_file)

        # End of installation
        os.chdir(current_dir)

    #---------------------------------------------------------------------------

    def test_library(self, executables=None, header=None, libname=None):

        libroot = None

        header_found = False
        lib_found = False

        search_dirs = ['/usr/local', '/usr']

        if executables != None:
            e = find_executable(executables)
            if e:
                if os.path.isabs(e):
                    libdir = os.path.split(os.path.dirname(e))[0]
                    if libdir not in search_dirs:
                        search_dirs.insert(0, libdir)

        for d in search_dirs:
            if os.path.isfile(os.path.join(d, 'include', header)):
                header_found = True
                libroot = d
                break

        if header_found:
            d = libroot
            lib_found = False
            if os.path.isfile(os.path.join(d, 'lib',
                                           'lib' + libname + '.so')):
                lib_found = True
            elif self.shared == False:
                if os.path.isfile(os.path.join(d, 'lib',
                                               'lib' + libname + '.a')):
                    lib_found = True

        # If library seems to be found, suggest it

        if header_found and lib_found:
            if libroot in ['/usr']:
                self.use = 'auto'
            else:
                self.use = 'yes'
                self.installation = 'no'
            self.install_dir = libroot

        # If headers found but not library, assume the library
        # is in a system path, so preselect 'auto' mode.

        elif header_found:

            self.use = 'auto'
            self.installation = 'no'
            self.install_dir = None

#-------------------------------------------------------------------------------
# Class definition for Code_Saturne setup
#-------------------------------------------------------------------------------

class Setup:

    def __init__(self):

        # Source directory
        self.top_srcdir = os.path.abspath(os.path.dirname(sys.argv[0]))

        # Optional libraries
        self.optlibs = ['hdf5', 'cgns', 'med', 'scotch', 'parmetis', 'libxml2']

        # Logging file
        self.log_file = sys.stdout

        # Download packages
        self.download = 'yes'

        # Code_Saturne language (may be en/fr)
        self.language = 'en'

        # Code_Saturne installation with debugging symbols
        self.debug = 'no'

        # Installation with shared libraries (not modifiable yet)
        self.shared = True

        # Default compilers
        self.cc = None
        self.fc = None
        self.cxx = None
        self.mpicc = None
        self.mpicxx = None

        # Disable GUI
        self.disable_gui = 'no'

        # Disable frontend
        self.disable_frontend = 'no'

        # Python interpreter path
        self.python = None

        # SALOME libraries (not application) path
        self.salome = None

        # Architecture name
        self.use_arch = 'no'
        self.arch = None

        # Installation prefix (if None, standard directory "/usr/local" will be used)
        self.prefix = None

        # Packages definition
        self.packages = {}

        # Code_Saturne

        self.packages['code_saturne'] = \
            Package(name="Code_Saturne",
                    description="Code_Saturne CFD tool",
                    package="code_saturne",
                    version=None,
                    archive=None,
                    url="http://code-saturne.org")

        p = self.packages['code_saturne']

        p.set_version_from_configure(os.path.join(self.top_srcdir, 'configure'))

        p.use = 'yes'
        p.installation = 'yes'

        # HDF5 library

        self.packages['hdf5'] = \
            Package(name="HDF5",
                    description="Hierarchical Data Format",
                    package="hdf5",
                    version="1.8.17",
                    archive="hdf5-1.8.17.tar.gz",
                    url="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.17/src/%s")

        p = self.packages['hdf5']
        p.config_opts = "--enable-production"

        # CGNS library

        self.packages['cgns'] = \
            Package(name="CGNS",
                    description="CFD General Notation System",
                    package="cgnslib",
                    version="3.2.1",
                    archive="cgnslib_3.2.1.tar.gz",
                    url="http://sourceforge.net/projects/cgns/files/cgnslib_3.2/%s/download")

        p = self.packages['cgns']
        p.config_opts = "-DCGNS_ENABLE_64BIT=ON -DCGNS_ENABLE_SCOPING=ON"

        # MED library

        self.packages['med'] = \
            Package(name="MED",
                    description="Model for Exchange of Data",
                    package="med",
                    version="3.2.1",
                    archive="med-3.2.1.tar.gz",
                    url="http://files.salome-platform.org/Salome/other/%s")

        p = self.packages['med']
        p.config_opts = "--with-med_int=long --disable-fortran --disable-python"

        # Libxml2 library (possible mirror at "ftp://fr.rpmfind.net/pub/libxml/%s")

        self.packages['libxml2'] = \
            Package(name="libxml2",
                    description="XML library",
                    package="libxml2",
                    version="2.9.2",
                    archive="libxml2-sources-2.9.2.tar.gz",
                    url="ftp://xmlsoft.org/libxml2/%s")

        p = self.packages['libxml2']
        p.config_opts = "--with-ftp=no --with-http=no"

        # ParMETIS

        self.packages['parmetis'] = \
            Package(name="parmetis",
                    description="ParMETIS",
                    package="parmetis",
                    version="4.0.3",
                    archive="parmetis-4.0.3.tar.gz",
                    url="http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/%s")

        # SCOTCH

        self.packages['scotch'] = \
            Package(name="scotch",
                    description="PT-Scotch",
                    package="scotch",
                    version="6.0.4",
                    archive="scotch_6.0.4.tar.gz",
                    url="https://gforge.inria.fr/frs/download.php/file/34618/%s")

    #---------------------------------------------------------------------------

    def setup_defaults(self):

        self.cc = find_executable(['cc', 'gcc', 'icc', 'xlc'], 'CC')
        self.fc = find_executable(['f95', 'gfortran', 'ifort'], 'FC')
        self.cxx = find_executable(['c++', 'g++', 'icpc', 'xlc++'], 'CXX')
        self.mpicc = find_executable(['mpicc', 'mpicc.openmpi', 'mpicc.mpich'])
        self.mpicxx = find_executable(['mpicxx', 'mpicxx.openmpi', 'mpicxx.mpich'])
        self.python = find_executable(['python'], 'PYTHON')

        # Architecture name
        self.arch = os.uname()[0] + '_' + os.uname()[4]

        # Installation prefix (if None, standard directory "/usr/local" will be used)
        self.prefix = '/usr/local'

        # Packages definition

        p = self.packages['hdf5']

        p.test_library(executables=['h5cc', 'h5pcc'],
                       header='H5public.h',
                       libname='hdf5')

        # CGNS library

        p = self.packages['cgns']

        p.test_library(executables=['cgnsnames'],
                       header='cgnslib.h',
                       libname='cgns')

        # MED library

        p = self.packages['med']

        p.test_library(executables=['mdump'],
                       header='med.h',
                       libname='medC')

        # Libxml2 library (possible mirror at "ftp://fr.rpmfind.net/pub/libxml/%s")
        p = self.packages['libxml2']

        p.test_library(header='libxml2/libxml/parser.h',
                       libname='xml2')

        # Expand user variables

        p = self.packages['code_saturne']
        from os.path import expanduser
        home = expanduser("~")
        self.prefix = os.path.join(home, 'Code_Saturne', p.version)

    #---------------------------------------------------------------------------

    def check_setup_file(self):

        # If setup file exists, nothing to do here
        if os.path.isfile('setup'):
            return

        # If setup does not exist, define on from template

        message = \
"""

Please edit the 'setup' file in the current directory
to define your Code_Saturne setup options.

You may then re-run '%(script_path)s'
to start the installation.

"""
        sys.stdout.write(message % {'script_path': sys.argv[0]})

        self.setup_defaults()
        self.write_setup()

        sys.exit(0)

    #---------------------------------------------------------------------------

    def read_setup(self):

        #
        # setup file reading
        #
        try:
            setupFile = open('setup', mode='r')
        except IOError:
            sys.stderr.write('Error: opening setup file\n')
            sys.exit(1)

        shutil.copy('setup','setup_ini')

        while 1:

            line = setupFile.readline()
            if line == '': break

            # skip comments
            if line[0] == '#': continue
            line = line.splitlines()
            list = line[0].split()
            # skip blank lines
            if len(list) == 0: continue

            key = list[0]

            if len(list) > 1:
                if key == 'download': self.download = list[1]
                elif key == 'prefix':
                    if not list[1] in ['default', 'auto']:
                        self.prefix = list[1]
                elif key == 'debug': self.debug = list[1]
                elif key == 'language': self.language = list[1]
                elif key == 'use_arch': self.use_arch = list[1]
                elif key == 'arch':
                    self.arch = list[1]
                    if self.arch == 'ignore':
                        self.use_arch = 'no'
                elif key == 'compCxx':
                    if not list[1] in ['default', 'auto']:
                        self.cxx = list[1]
                elif key == 'compC': self.cc = list[1]
                elif key == 'compF': self.fc = list[1]
                elif key == 'mpiCompC':
                    if not list[1] in ['default', 'auto']:
                        self.mpicc = list[1]
                elif key == 'mpiCompCxx':
                    if not list[1] in ['default', 'auto']:
                        self.mpicxx = list[1]
                elif key == 'disable_gui': self.disable_gui = list[1]
                elif key == 'disable_frontend': self.disable_frontend = list[1]
                elif key == 'python':
                    if not list[1] in ['default', 'auto']:
                        self.python = list[1]
                elif key == 'salome':
                    if not list[1] in ['default', 'auto', 'no']:
                        self.salome = list[1]
                else:
                    p = self.packages[key]
                    p.use = list[1]
                    p.installation = list[2]
                    if (p.use != 'no'):
                        if list[3] != 'None':
                            p.install_dir = list[3]

        # Specify architecture name
        if self.use_arch == 'yes' and self.arch is None:
            self.arch = os.uname()[0] + '_' + os.uname()[4]

        # Expand user variables
        if self.prefix:
            self.prefix = os.path.expanduser(self.prefix)
            self.prefix = os.path.expandvars(self.prefix)
            self.prefix = os.path.abspath(self.prefix)

        if self.python:
            self.python = os.path.expanduser(self.python)
            self.python = os.path.expandvars(self.python)
            self.python = os.path.abspath(self.python)

        if self.salome:
            self.salome = os.path.expanduser(self.salome)
            self.salome = os.path.expandvars(self.salome)
            self.salome = os.path.abspath(self.salome)

    #---------------------------------------------------------------------------

    def check_setup(self):

        check = """
Check the setup file and some utilities presence.
"""

        sys.stdout.write(check)
        if verbose == 'yes':
            sys.stdout.write("\n")

        # Testing download option
        if self.download not in ['yes', 'no']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'download\' option in the setup file "
                             "should be \'yes\' or \'no\'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Testing debug option
        if self.debug not in ['yes', 'no']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'debug\' option in the setup file "
                             "should be \'yes\' or \'no\'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Testing GUI option
        if self.disable_gui not in ['yes', 'no']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'disable_gui\' option in the setup file "
                             "should be \'yes\' or \'no\'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Testing frontend option
        if self.disable_frontend not in ['yes', 'no']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'disable_frontend\' option in the setup file "
                             "should be \'yes\' or \'no\'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Testing language option
        if self.language not in ['en', 'fr']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'language\' option in the setup file "
                             "should be \'en\' or \'fr'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Testing prefix directory
        if self.prefix and not os.path.isdir(self.prefix):
            try:
                os.makedirs(self.prefix)
            except Exception:
                pass
        if self.prefix and not os.path.isdir(self.prefix):
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'%s\' prefix directory is provided in the setup "
                             "file but is not a directory.\n"
                             "Please check your setup file.\n\n"
                             % self.prefix)
            sys.exit(1)

        # Testing architecture option
        if self.use_arch not in ['yes', 'no']:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'use_arch\' option in the setup file "
                             "should be \'yes\' or \'no\'.\n"
                             "Please check your setup file.\n\n")
            sys.exit(1)

        # Looking for compilers provided by the user
        for compiler in [self.cc, self.mpicc, self.fc]:
            if compiler:
                ret = run_test(compiler)
                if ret != 0:
                    sys.stderr.write("\n*** Aborting installation:\n"
                                     "\'%s\' compiler is provided in the setup "
                                     "file but cannot be found.\n"
                                     "Please check your setup file.\n\n"
                                     % compiler)
                    sys.exit(1)

        # Looking for Python executable provided by the user
        python = 'python'
        if self.python: python = self.python
        ret = run_test(python)
        if ret != 0:
            if self.python:
                sys.stderr.write("\n*** Aborting installation:\n"
                                 "\'%s\' Python exec is provided in the setup "
                                 "file doesn't not seem to be executable.\n"
                                 "Please check your setup file.\n\n"
                                 % self.python)
            else:
                sys.stderr.write("\n*** Aborting installation:\n"
                                 "Cannot find Python executable.\n"
                                 "Please check your setup file.\n\n")
            sys.exit(1)
        else:
            cmd = python + " -c \'import sys; print(sys.version[:3])\'"
            if verbose == 'yes':
                sys.stdout.write("     Python version is ")
            p = subprocess.Popen(cmd,
                                 shell=True,
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)

            output = p.communicate()
            if verbose == 'yes':
                if p.returncode == 0:
                    sys.stdout.write(output[0])

        # Checking libraries options
        for lib in self.optlibs:
            p = self.packages[lib]
            if p.use not in ['yes', 'no', 'auto']:
                sys.stderr.write("\n*** Aborting installation:\n"
                                 "\'%s\' use option in the setup file "
                                 "should be \'yes\', \'no' or \'auto\'.\n"
                                 "Please check your setup file.\n\n"
                                 % lib)
                sys.exit(1)
            if p.installation not in ['yes', 'no']:
                sys.stderr.write("\n*** Aborting installation:\n"
                                 "\'%s\' install option in the setup file "
                                 "should be \'yes\' or \'no'.\n"
                                 "Please check your setup file.\n\n"
                                 % lib)
                sys.exit(1)
            if p.installation == 'no' and p.use == 'yes':
                if not os.path.isdir(p.install_dir):
                    sys.stderr.write("\n*** Aborting installation:\n"
                                     "\'%(path)s\' path is provided for "
                                     "\'%(lib)s\' in the setup "
                                     "file but is not a directory.\n"
                                     "Please check your setup file.\n\n"
                                     % {'path':p.install_dir, 'lib':lib})
                    sys.exit(1)

        # Looking for SALOME path probided by the user
        if self.salome and not os.path.isdir(self.salome):
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'%s\' SALOME directory is provided in the setup "
                             "file but is not present.\n"
                             "Please check your setup file.\n\n"
                             % self.salome)
            sys.exit(1)


        # Looking for make utility
        ret = run_test("make")
        if ret != 0:
            sys.stderr.write("\n*** Aborting installation:\n"
                             "\'make\' utility is mandatory for Code_Saturne "
                             "compilation.\n"
                             "Please install development tools.\n\n")
            sys.exit(1)

        if verbose == 'yes':
            sys.stdout.write("\n")

    #---------------------------------------------------------------------------

    def update_package_opts(self):

        # Update log file, installation directory and compilers
        for lib in self.optlibs + ['code_saturne']:
            p = self.packages[lib]
            # Update logging file
            p.log_file = self.log_file
            # Installation directory
            if p.installation == 'yes' and not p.install_dir:
                subdir = os.path.join(p.package + '-' + p.version)
                if self.arch:
                    subdir = os.path.join(subdir, 'arch', self.arch)
                p.install_dir = os.path.join(self.prefix, subdir)
            # Compilers
            p.cc = self.cc
            p.cxx = self.cxx
            if lib in ['scotch'] and self.mpicc:
                p.cc = self.mpicc
            elif lib in ['code_saturne', 'parmetis']:
                if self.mpicc:
                    p.cc = self.mpicc
                if self.mpicxx:
                    p.cxx = self.mpicxx
            if lib in ['code_saturne']:
                p.fc = self.fc
            p.shared = self.shared

        # Update configuration options

        config_opts = ''
        if self.debug == 'yes':
            config_opts = config_opts + " --enable-debug"

        hdf5 = self.packages['hdf5']
        cgns = self.packages['cgns']
        med= self.packages['med']
        scotch = self.packages['scotch']
        parmetis = self.packages['parmetis']
        libxml2 = self.packages['libxml2']

        # Disable GUI

        if self.disable_gui == 'yes':
            config_opts = config_opts + " --disable-gui"

        # Disable frontend

        if self.disable_frontend == 'yes':
            config_opts = config_opts + " --disable-frontend"

        # HDF5 (needed for MED and recommended for CGNS)

        if hdf5.use == 'no':
            config_opts = config_opts + " --without-hdf5"
        else:
            cgns.config_opts += " -DCGNS_ENABLE_HDF5=ON"
            if hdf5.install_dir:
                config_opts = config_opts + " --with-hdf5=" + hdf5.install_dir
                med.config_opts += " --with-hdf5=" + hdf5.install_dir
                cgns.config_opts += " -DHDF5_INCLUDE_PATH=" + hdf5.install_dir + "/include" \
                    + " -DHDF5_LIBRARY=" + hdf5.install_dir + "/lib/libhdf5.so"

        # CGNS

        if cgns.use == 'no':
            config_opts = config_opts + " --without-cgns"
        else:
            if cgns.install_dir:
                config_opts = config_opts + " --with-cgns=" + cgns.install_dir

        # MED

        if med.use == 'no':
            config_opts = config_opts + " --without-med"
        else:
            if med.install_dir:
                config_opts = config_opts + " --with-med=" + med.install_dir

        # ParMetis

        if parmetis.use == 'no':
            config_opts = config_opts + " --without-metis"
        else:
            config_opts = config_opts + " --with-metis=" + parmetis.install_dir

        # PT-Scotch

        if scotch.use == 'no':
            config_opts = config_opts + " --without-scotch"
        else:
            config_opts = config_opts + " --with-scotch=" + scotch.install_dir

        # Libxml2

        if libxml2.use == 'no':
            config_opts = config_opts + " --without-libxml2"
        else:
            if libxml2.install_dir:
                config_opts = config_opts + \
                    " --with-libxml2=" + libxml2.install_dir

        # Python

        if self.python:
            config_opts = config_opts + " PYTHON=" + self.python

        # SALOME

        if self.salome:
            config_opts = config_opts + " --with-salome=" + self.salome

        # Language

        if self.language == 'fr':
            config_opts = config_opts + " --enable-french"

        # Build type

        if self.shared:
            config_opts += " --disable-static"
        else:
            config_opts += " --disable-shared"

        self.packages['code_saturne'].config_opts = config_opts

    #---------------------------------------------------------------------------

    def install(self):

        if self.download == 'yes':
            for lib in self.optlibs:
                p = self.packages[lib]
                if p.installation == 'yes':
                    sys.stdout.write("Download of %s\n  (%s)\n" % (p.name, p.url))
                    p.download()
            self.download = 'no'
            self.write_setup()
            sys.stdout.write("\n")

        for lib in self.optlibs:
            p = self.packages[lib]
            p.info()
            if p.installation == 'yes':
                sys.stdout.write("Extract of %s\n" % p.name)
                p.extract()
                sys.stdout.write("Installation of %s\n" % p.name)
                if verbose == 'yes':
                    p.info()
                if lib == 'scotch':
                    p.install_ptscotch()
                elif lib == 'parmetis':
                    p.install_parmetis()
                else:
                    p.install()
                p.installation = 'no'
                self.write_setup()
                if verbose == 'yes':
                    sys.stdout.write("\n")

        p = self.packages['code_saturne']
        p.info()
        if p.installation == 'yes':
            p.source_dir = self.top_srcdir
            sys.stdout.write("Installation of %s\n" % p.name)
            if verbose == 'yes':
                p.info()
            p.install()
            p.installation = 'no'
            self.write_setup()
            if verbose == 'yes':
                sys.stdout.write("\n")

    #---------------------------------------------------------------------------

    def write_setup(self):
        #
        # setup file update
        #
        sf = open(os.path.join(os.getcwd(), "setup"), mode='w')

        setupMain = \
"""#========================================================
# Setup file for Code_Saturne installation
#========================================================
#
#--------------------------------------------------------
# Download packages ?
#--------------------------------------------------------
download  %(download)s
#
#--------------------------------------------------------
# Language
#   default: "en" english
#   others:  "fr" french
#--------------------------------------------------------
language  %(lang)s
#
#--------------------------------------------------------
# Install Code_Saturne with debugging symbols
#--------------------------------------------------------
debug     %(debug)s
#
#--------------------------------------------------------
# Installation directory
#--------------------------------------------------------
prefix    %(prefix)s
#
#--------------------------------------------------------
# Optional architecture Name (installation subdirectory)
#--------------------------------------------------------
use_arch  %(use_arch)s
arch      %(arch)s
#
#--------------------------------------------------------
# C compiler and optional MPI wrapper
#--------------------------------------------------------
compC     %(cc)s
mpiCompC  %(mpicc)s
#
#--------------------------------------------------------
# Fortran compiler
#--------------------------------------------------------
compF    %(fc)s
#
#--------------------------------------------------------
# C++ compiler and MPI wrapper for optional packages
#
# Required only for static builds using the MED library
# or for build of optional modules such as MEDCoupling
# support.
#--------------------------------------------------------
compCxx     %(cxx)s
mpiCompCxx  %(mpicxx)s
#
#--------------------------------------------------------
# Python interpreter.
#--------------------------------------------------------
python    %(python)s
#
#--------------------------------------------------------
# Disable the Graphical user Interface ?
#--------------------------------------------------------
disable_gui  %(disable_gui)s
#
#--------------------------------------------------------
# Disable frontend (also disables GUI) ?
# May be useful for debug builds and HPC cluster builds
# installed side-by side with a full build.
#--------------------------------------------------------
disable_frontend  %(disable_frontend)s
#
#--------------------------------------------------------
# Optional SALOME platform install path.
#
# This is the path for the main SALOME directory,
# not the application directory.
#
# If Code_Saturne is built with SALOME support,
# running "code_saturne salome" will launch the
# associated application, containing the CFDSTUDY module.
#--------------------------------------------------------
salome    %(salome)s
#
#--------------------------------------------------------
# Optional packages:
# ------------------
#
# MED / HDF5  For MED file format support
#             (used by SALOME and by Gmsh)
#
# CGNS / HDF5 For CGNS file support
#             (used by many meshing tools)
#
# Scotch (includes PT-Scotch) and/or ParMetis
# for parallel partitioning
#
#   For Linux workstations, HDF5, CGNS, and even MED
# packages may be available through the package manager.
# HDF5 is also often available on large systems.
# When building with SALOME, the platform distribution's
# packages may be used, by setting 'salome' in the
# matching entry under the "Use" column.
#
# Scotch and Pt-Scotch are available in some Linux
# distributions, but may be built with options
# incompatible with non-threaded Code_Saturne runs.
#
#   To install CGNS or ParMetis, the CMake
# configuration/installation tool is required
# (it is available in most Linux distributions)
#
#   Libxml2 is needed to read XML files output by the
# Graphical User Interface. It is generally available
# through the package manager.
#--------------------------------------------------------
#
#  Name    Use   Install  Path
#
"""
        setupLib= \
"""%(lib)-9s  %(use)-4s  %(install)-3s      %(dir)s
"""
        setupEnd= \
"""#
#========================================================
"""

        prefix = self.prefix
        arch = self.arch
        cc = self.cc
        fc = self.fc
        mpicc = self.mpicc
        cxx = self.cxx
        mpicxx = self.mpicxx
        python = self.python
        salome = self.salome

        # Clean some potentially undefined variables for output
        if not prefix: prefix = 'default'
        if not arch: arch = 'ignore'
        if not cc: cc = 'NEEDS_DEFINITION'
        if not fc: fc = 'NEEDS_DEFINITION'
        if not mpicc: mpicc = 'auto'
        if not cxx: cxx = 'auto'
        if not mpicxx: mpicxx = 'auto'
        if not python: python = 'NEEDS_DEFINITION'
        if not salome: salome = 'no'

        sf.write(setupMain
                 % { 'download':self.download, 'prefix':prefix,
                     'lang':self.language, 'debug':self.debug,
                     'use_arch':self.use_arch, 'arch':arch,
                     'cc':cc, 'mpicc':mpicc,
                     'fc':fc,
                     'cxx':cxx, 'mpicxx':mpicxx,
                     'disable_gui':self.disable_gui,
                     'disable_frontend':self.disable_frontend,
                     'python':self.python, 'salome':salome})

        for lib in self.optlibs:
            p = self.packages[lib]
            sf.write(setupLib % {'lib':lib,
                                 'use':p.use,
                                 'install':p.installation,
                                 'dir':p.install_dir})

        sf.write(setupEnd)
        sf.close()

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # Messages
    # --------
    welcome = \
        """
        Installation of Code_Saturne
        ____________________________

The process will take several minutes.
You can have a look at the log file meanwhile.
"""

    finalize = \
"""
Before using Code_Saturne, please update your environment with:

  cspath=%(cspath)s
  alias code_saturne="$cspath/code_saturne"

The documentation should then be available through the commands:
  code_saturne info -g refcard
  code_saturne info -g user

Do not forget the post-installation steps recommended in the
installation documentation, available using:
  code_saturne info -g install

"""

    thanks = \
"""
Thank you for choosing Code_Saturne!

"""

    # Setup process
    # -------------

    check_directory()

    sys.stdout.write(welcome)

    setup = Setup()

    setup.check_setup_file()

    setup.log_file = open('install_saturne.log', mode='w')

    setup.read_setup()

    setup.check_setup()
    setup.update_package_opts()
    setup.install()

    setup.log_file.close()

    cspath = os.path.join(setup.packages['code_saturne'].install_dir, 'bin')
    sys.stdout.write(finalize % {'cspath':cspath})
    sys.stdout.write(thanks)
