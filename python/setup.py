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

# ------------------------------------------------------------------------------
# Simple version string parser to avoid errors
# ------------------------------------------------------------------------------

def _parse_saturne_version_string(version_string):
    """
    Parse a version string and return a PEP440 compliant one.
    """

    retval = ''
    if version_string and len(version_string) > 0:
        _string = version_string.split('-')
        retval = _string[0]
        if len(_string) > 1:
            if _string[1] == "alpha":
                retval += 'a0'
            elif _string[1] == "beta":
                retval += 'b0'
            elif _string[1][:2] == 'rc':
                retval += _string[1]
            elif _string[1] == 'patch':
                retval += '.post0'
            else:
                retval += '.dev0'

    else:
        retval = 'master'

    return retval

# ------------------------------------------------------------------------------

_cs_opts = {'enable_gui':False,
            'use_qt':None,
            'exclude_dirs':[],
            'pyuic':None,
            'uic':None,
            'pyrcc':None,
            'rcc':None,
            'version':None
            }

# version
if '--version' in sys.argv:
    index = sys.argv.index('--version')
    sys.argv.pop(index)
    _cs_opts['version'] = _parse_saturne_version_string(sys.argv.pop(index))
else:
    _cs_opts['version'] = 'master'

# GUI Enabled or not

if '--with-gui' in sys.argv:
    _cs_opts['enable_gui'] = True
    sys.argv.remove('--with-gui')
elif '--without-gui' in sys.argv:
    _cs_opts['enable_gui'] = False
    sys.argv.remove('--without-gui')

if '--use-qt=pyqt5' in sys.argv:
    _cs_opts['use_qt'] = 'pyqt5'
    sys.argv.remove('--use-qt=pyqt5')
elif '--use-qt=pyqt6' in sys.argv:
    _cs_opts['use_qt'] = 'pyqt6'
    sys.argv.remove('--use-qt=pyqt6')
elif '--use-qt=pyside6' in sys.argv:
    _cs_opts['use_qt'] = 'pyside6'
    sys.argv.remove('--use-qt=pyside6')
else:
    _cs_opts['use_qt'] = False

if '--pyuic' in sys.argv:
    index = sys.argv.index('--pyuic')
    sys.argv.pop(index)
    _cs_opts['pyuic'] = sys.argv.pop(index)

if '--pyrcc' in sys.argv:
    index = sys.argv.index('--pyrcc')
    sys.argv.pop(index)
    _cs_opts['pyrcc'] = sys.argv.pop(index)

if '--uic' in sys.argv:
    index = sys.argv.index('--uic')
    sys.argv.pop(index)
    _cs_opts['uic'] = sys.argv.pop(index)

if '--rcc' in sys.argv:
    index = sys.argv.index('--rcc')
    sys.argv.pop(index)
    _cs_opts['rcc'] = sys.argv.pop(index)

# Sanity checks
if _cs_opts['enable_gui']:
    if not (_cs_opts['pyuic'] or _cs_opts['uic']):
        raise Exception("--pyuic must be set to PyQt pyuic or uic set to Qt uic")

    if not (_cs_opts['pyrcc'] or _cs_opts['rcc']):
        raise Exception("--pyrcc must be set to PyQt pyrcc or rcc set to Qt rcc")
else:
    _cs_opts['exclude_dirs'] = ['*gui*']
    # If gui  disabled, remove all options
    _cs_opts['pyuic'] = None
    _cs_opts['uic'] = None
    _cs_opts['pyrcc'] = None
    _cs_opts['rcc'] = None
    _cs_opts['use_qt'] = False

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
        qtpkg = None
        if _cs_opts['use_qt'] == 'pyqt5':
            qtpkg = 'qt5'
        elif _cs_opts['use_qt'] == 'pyqt6':
            qtpkg = 'qt6'
        elif _cs_opts['use_qt'] == 'pyside6':
            qtpkg = 'pyside6'
        if build_command and qtpkg:
            _p1 = os.path.join('code_saturne', 'gui', 'base')
            for _f in glob.glob(os.path.join(SRC_PATH, _p1,
                                             qtpkg,
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
        cmd = None
        cmd_post = None
        if _cs_opts['pyuic']:
            cmd = [_cs_opts['pyuic']]
        elif _cs_opts['uic']:
            cmd = [_cs_opts['uic'], '-g', 'python']
            if _cs_opts['use_qt'] == 'pyqt5':
                cmd_post = "sed -i -e 's|PySide2|PyQt5|'"
            elif _cs_opts['use_qt'] == 'pyqt6':
                cmd_post = "sed -i -e 's|PySide6|PyQt6|'"
        if build_command and cmd:
            for _f in glob.glob(SRC_PATH+'/**/*.ui', recursive=True):
                _ui_f = os.path.join(build_command.build_base,
                                     'lib',
                                     os.path.relpath(_f.replace('.ui','.py'),
                                                     SRC_PATH)
                                     )
                cmd_ui = cmd + ['-o', _ui_f, _f]
                subprocess.run(args=cmd_ui, check=True)
                if cmd_post:
                    cmd_post +=  " " + _rc_f
                    subprocess.run(cmd_post, shell=True, check=True)

class _build_rc_files(_cs_build_py_dummy):
    """
    A custom command to generate the qrc based python files using pyrcc/rcc
    """

    def run(self):
        build_command = self.distribution.command_obj.get('build', None)
        cmd = None
        cmd_post = None
        if _cs_opts['pyrcc']:
            cmd = [_cs_opts['pyrcc']]
        elif _cs_opts['rcc']:
            cmd = [_cs_opts['rcc'], '-g', 'python']
            if _cs_opts['use_qt'] == 'pyqt5':
                cmd_post = "sed -i -e 's|PySide2|PyQt5|'"
            elif _cs_opts['use_qt'] == 'pyqt6':
                cmd_post = "sed -i -e 's|PySide6|PyQt6|'"
        if build_command and cmd:
            for _f in glob.glob(SRC_PATH+'/**/*.qrc', recursive=True):
                _rc_f = os.path.join(build_command.build_base,
                                     'lib',
                                     os.path.relpath(_f.replace('.qrc','_rc.py'),
                                                     SRC_PATH)
                                     )
                _rc_dir = os.path.dirname(_f)
                cmd_rc = cmd + ['-o', _rc_f, _f]
                subprocess.run(args=cmd_rc, check=True)
                if cmd_post:
                    cmd_post +=  " " + _rc_f
                    subprocess.run(cmd_post, shell=True, check=True)

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
