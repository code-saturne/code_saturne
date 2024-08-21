# ==============================================================================
# IMPORTS
# ==============================================================================

from setuptools import setup, find_packages, Command
from setuptools.command.egg_info import egg_info
from setuptools.command.install_egg_info import install_egg_info
from setuptools.command.build_py import build_py as _build_py
import glob
import os
import shutil
import subprocess
import sys

# ==============================================================================
# Get setup.py path needed for autoconf
# ==============================================================================

SRC_PATH = os.path.relpath(os.path.join(os.path.dirname(__file__)))

# ==============================================================================
# Get environment variables needed
# ==============================================================================

_cs_opts = {'enable_gui':False,
            'use_qt5':False,
            'exclude_dirs':[],
            'qtpkg':None,
            'pyuic':None,
            'pyrcc':None,
            'version':None
            }

# version
if '--version' in sys.argv:
    index = sys.argv.index('--version')
    sys.argv.pop(index)
    _cs_opts['version'] = sys.argv.pop(index)
else:
    _cs_opts['version'] = 'master'

# GUI Enabled or not

if '--with-gui' in sys.argv:
    _cs_opts['enable_gui'] = True
    sys.argv.remove('--with-gui')
elif '--without-gui' in sys.argv:
    _cs_opts['enable_gui'] = False
    sys.argv.remove('--without-gui')

if '--use-qt5' in sys.argv:
    _cs_opts['use_qt5'] = True
    sys.argv.remove('--use-qt5')
    _cs_opts['qtpkg'] = 'qt5'
else:
    _cs_opts['use_qt5'] = False
    _cs_opts['qtpkg'] = 'qt4'

if '--pyuic' in sys.argv:
    index = sys.argv.index('--pyuic')
    sys.argv.pop(index)
    _cs_opts['pyuic'] = sys.argv.pop(index)

if '--pyrcc' in sys.argv:
    index = sys.argv.index('--pyrcc')
    sys.argv.pop(index)
    _cs_opts['pyrcc'] = sys.argv.pop(index)

# Sanity checks
if _cs_opts['enable_gui']:
    if not _cs_opts['pyuic']:
        raise Exception("--pyuic must be set to PyQt pyuic")

    if not _cs_opts['pyrcc']:
        raise Exception("--pyrcc must be set to PyQt pyrcc")
else:
    _cs_opts['exclude_dirs'] = ['*gui*']
    # If gui  disabled, remove all options
    _cs_opts['pyuic'] = None
    _cs_opts['pyrcc'] = None
    _cs_opts['qtpkg'] = None
    _cs_opts['use_qt5'] = False


# ==============================================
# Overloading the egg-info command
# ==============================================

# We overload the egg_info command to avoid creating the ".egg-info"
# folder in the sources folder

class _cs_egg_info(egg_info):

    def run(self):
        build_command = self.distribution.command_obj.get('build', None)
        if build_command:
            self.egg_base = build_command.build_base
            self.egg_info = os.path.join(self.egg_base, os.path.basename(self.egg_info))

        super().run()

# We overload the install_egg_info to ensure that the correct "egg_info"
# command is found, otherwise the wrong egg_info path is used...

class _cs_install_egg_info(install_egg_info):

    def finalize_options(self):
        self.run_command('egg_info')
        super().finalize_options()

# ==============================================
# Specific commands to build python files for GUI
# ==============================================

class _cs_build_py_dummy(Command):
    """
    A custom command used to factorize code.
    """

    user_options = []

    # ==============
    # Public methods
    # ==============

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        raise Exception('Method non implemented for dummy class')

class _build_qt_files(_cs_build_py_dummy):
    """
    A custom command to generate the Qt files based on distribution of Qt
    """

    def run(self):
        build_command = self.distribution.command_obj.get('build', None)
        if build_command and _cs_opts['qtpkg']:
            _p1 = os.path.join('code_saturne', 'gui', 'base')
            for _f in glob.glob(os.path.join(SRC_PATH, _p1,
                                             _cs_opts['qtpkg'],
                                             '*.py')):
                _abs_dst = os.path.join(SRC_PATH, _p1, os.path.basename(_f))
                _nf = os.path.join(build_command.build_base,
                                   'lib',
                                   os.path.relpath(_abs_dst, SRC_PATH))

                shutil.copy2(_f, _nf, follow_symlinks=False)

class _build_ui_files(_cs_build_py_dummy):
    """
    A custom command to generate the ui based python files using pyuic
    """

    def run(self):
        build_command = self.distribution.command_obj.get('build', None)
        if build_command and _cs_opts['pyuic']:
            for _f in glob.glob(SRC_PATH+'/**/*.ui', recursive=True):
                _ui_f = os.path.join(build_command.build_base,
                                     'lib',
                                     os.path.relpath(_f.replace('.ui','.py'),
                                                     SRC_PATH)
                                     )
                subprocess.run(args=[_cs_opts['pyuic'], '-o', _ui_f, _f],
                               check=True)

class _build_rc_files(_cs_build_py_dummy):
    """
    A custom command to generate the qrc based python files using pyrcc
    """

    def run(self):
        build_command = self.distribution.command_obj.get('build', None)
        if build_command and _cs_opts['pyrcc']:
            for _f in glob.glob(SRC_PATH+'/**/*.qrc', recursive=True):
                _rc_f = os.path.join(build_command.build_base,
                                     'lib',
                                     os.path.relpath(_f.replace('.qrc','_rc.py'),
                                                     SRC_PATH)
                                     )
                _rc_dir = os.path.dirname(_f)
                subprocess.run(args=[_cs_opts['pyrcc'], '-o', _rc_f, _f],
                               check=True)

# =========================
# Overload build_py command
# =========================

class build_py(_build_py):

    def initialize_options(self):
        super().initialize_options()

    def finalize_options(self):
        super().finalize_options()

    def run(self):
        super().run()
        self.run_command("build_cs_qt_files")
        self.run_command("build_cs_ui_files")
        self.run_command("build_cs_rc_files")

setup(name='code_saturne',
      author='code_saturne dev team',
      author_email='saturne-support@edf.fr',
      description='Python CLI and GUI of code_saturne multiphysics CFD package',
      package_dir={'':SRC_PATH},
      package_data={'code_saturne':['data/icons/*', 'data/icons/22x22/*', 'data/icons/32x32/*'],},
      cmdclass={'build_py':build_py,
                'build_cs_qt_files':_build_qt_files,
                'build_cs_ui_files':_build_ui_files,
                'build_cs_rc_files':_build_rc_files,
                'egg_info':_cs_egg_info,
                'install_egg_info':_cs_install_egg_info},
      version=_cs_opts['version'],
      packages=find_packages(where=SRC_PATH, exclude=_cs_opts['exclude_dirs'])
      )
