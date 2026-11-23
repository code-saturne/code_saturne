# ==============================================================================
# IMPORTS
# ==============================================================================

from setuptools import setup, find_namespace_packages, Command
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

_setup_dir = os.path.dirname(os.path.abspath(__file__))
SRC_PATH = os.path.relpath(_setup_dir)

# ==============================================================================
# Helpers
# ==============================================================================

def _replace_in_file(filepath, old_text, new_text):
    """
    Portably replace occurrences of old_text with new_text in a text file.
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    content = content.replace(old_text, new_text)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

# ------------------------------------------------------------------------------
# Version string parser (PEP 440 compliant)
# ------------------------------------------------------------------------------

def _parse_saturne_version_string(version_string):
    """
    Parse a version string and return a PEP 440 compliant version.
    """
    if version_string and len(version_string) > 0 and version_string != 'master':
        _string = version_string.split('-')
        retval = _string[0]
        if len(_string) > 1:
            suffix = _string[1]
            if suffix == "alpha":
                retval += 'a0'
            elif suffix == "beta":
                retval += 'b0'
            elif suffix.startswith('rc'):
                retval += suffix
            elif suffix == 'patch':
                retval += '.post0'
            else:
                retval += '.dev0'
        return retval
    return '0.0.0.dev0'

# ------------------------------------------------------------------------------
# Option parsing: environment variables with sys.argv overrides
# ------------------------------------------------------------------------------

_cs_opts = {'enable_gui': os.environ.get('CS_ENABLE_GUI', '').lower() in ('1', 'true', 'yes'),
            'use_qt': os.environ.get('CS_USE_QT', None),
            'exclude_dirs': [],
            'pyuic': os.environ.get('CS_PYUIC', None),
            'uic': os.environ.get('CS_UIC', None),
            'pyrcc': os.environ.get('CS_PYRCC', None),
            'rcc': os.environ.get('CS_RCC', None),
            'version': _parse_saturne_version_string(os.environ.get('CS_VERSION', None)),
            'build_base': os.environ.get('CS_PYBUILD_DIR', None)
            }

def _pop_option(arg_name):
    """Safely extract option and value from sys.argv."""
    if arg_name in sys.argv:
        idx = sys.argv.index(arg_name)
        sys.argv.pop(idx)
        if idx < len(sys.argv):
            return sys.argv.pop(idx)
    return None

# Build directory CLI flag
val = _pop_option('--build-base')
if val:
    _cs_opts['build_base'] = val

# Version CLI flag
val = _pop_option('--version')
if val:
    _cs_opts['version'] = _parse_saturne_version_string(val)
elif not _cs_opts.get('version') or _cs_opts['version'] == '0.0.0.dev0':
    _cs_opts['version'] = '0.0.0.dev0'

# GUI Enabled or not
if '--with-gui' in sys.argv:
    _cs_opts['enable_gui'] = True
    sys.argv.remove('--with-gui')
elif '--without-gui' in sys.argv:
    _cs_opts['enable_gui'] = False
    sys.argv.remove('--without-gui')

# Qt binding
for qt_arg in ('--use-qt=pyqt5', '--use-qt=pyqt6', '--use-qt=pyside6'):
    if qt_arg in sys.argv:
        _cs_opts['use_qt'] = qt_arg.split('=')[1]
        sys.argv.remove(qt_arg)

val = _pop_option('--pyuic')
if val:
    _cs_opts['pyuic'] = val

val = _pop_option('--pyrcc')
if val:
    _cs_opts['pyrcc'] = val

val = _pop_option('--uic')
if val:
    _cs_opts['uic'] = val

val = _pop_option('--rcc')
if val:
    _cs_opts['rcc'] = val

# Sanity checks
if _cs_opts['enable_gui']:
    if not (_cs_opts['pyuic'] or _cs_opts['uic']):
        raise Exception("--pyuic must be set to PyQt pyuic or uic set to Qt uic")

    if not (_cs_opts['pyrcc'] or _cs_opts['rcc']):
        raise Exception("--pyrcc must be set to PyQt pyrcc or rcc set to Qt rcc")
else:
    _cs_opts['exclude_dirs'] = ['*gui*']
    # If gui disabled, reset options
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
        if build_command and getattr(build_command, 'build_base', None):
            self.egg_base = build_command.build_base
            self.egg_info = os.path.join(self.egg_base,
                                         os.path.basename(self.egg_info))

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
        build_py_cmd = self.get_finalized_command('build_py')
        build_lib = getattr(build_py_cmd, 'build_lib', None)
        if not build_lib:
            build_command = self.distribution.command_obj.get('build', None)
            build_lib = os.path.join(build_command.build_base, 'lib') if build_command else None

        qtpkg = None
        if _cs_opts['use_qt'] == 'pyqt5':
            qtpkg = 'qt5'
        elif _cs_opts['use_qt'] == 'pyqt6':
            qtpkg = 'qt6'
        elif _cs_opts['use_qt'] == 'pyside6':
            qtpkg = 'pyside6'

        if build_lib and qtpkg:
            _p1 = os.path.join('code_saturne', 'gui', 'base')
            for _f in glob.glob(os.path.join(SRC_PATH, _p1, qtpkg, '*.py')):
                _abs_dst = os.path.join(SRC_PATH, _p1, os.path.basename(_f))
                _nf = os.path.join(build_lib, os.path.relpath(_abs_dst, SRC_PATH))
                os.makedirs(os.path.dirname(_nf), exist_ok=True)
                shutil.copy2(_f, _nf, follow_symlinks=False)

class _build_ui_files(_cs_build_py_dummy):
    """
    A custom command to generate the ui based python files using pyuic
    """

    def run(self):
        build_py_cmd = self.get_finalized_command('build_py')
        build_lib = getattr(build_py_cmd, 'build_lib', None)
        if not build_lib:
            build_command = self.distribution.command_obj.get('build', None)
            build_lib = os.path.join(build_command.build_base, 'lib') if build_command else None

        cmd = None
        post_replace = None
        if _cs_opts['pyuic']:
            cmd = [_cs_opts['pyuic']]
        elif _cs_opts['uic']:
            cmd = [_cs_opts['uic'], '-g', 'python']
            if _cs_opts['use_qt'] == 'pyqt5':
                post_replace = ('PySide2', 'PyQt5')
            elif _cs_opts['use_qt'] == 'pyqt6':
                post_replace = ('PySide6', 'PyQt6')

        if build_lib and cmd:
            for _f in glob.glob(SRC_PATH + '/**/*.ui', recursive=True):
                _ui_f = os.path.join(build_lib,
                                     os.path.relpath(_f.replace('.ui', '.py'),
                                                     SRC_PATH))
                os.makedirs(os.path.dirname(_ui_f), exist_ok=True)
                cmd_ui = cmd + ['-o', _ui_f, _f]
                subprocess.run(args=cmd_ui, check=True)
                if post_replace:
                    _replace_in_file(_ui_f, post_replace[0], post_replace[1])

class _build_rc_files(_cs_build_py_dummy):
    """
    A custom command to generate the qrc based python files using pyrcc/rcc
    """

    def run(self):
        build_py_cmd = self.get_finalized_command('build_py')
        build_lib = getattr(build_py_cmd, 'build_lib', None)
        if not build_lib:
            build_command = self.distribution.command_obj.get('build', None)
            build_lib = os.path.join(build_command.build_base, 'lib') if build_command else None

        cmd = None
        post_replace = None
        if _cs_opts['pyrcc']:
            cmd = [_cs_opts['pyrcc']]
        elif _cs_opts['rcc']:
            cmd = [_cs_opts['rcc'], '-g', 'python']
            if _cs_opts['use_qt'] == 'pyqt5':
                post_replace = ('PySide2', 'PyQt5')
            elif _cs_opts['use_qt'] == 'pyqt6':
                post_replace = ('PySide6', 'PyQt6')

        if build_lib and cmd:
            for _f in glob.glob(SRC_PATH + '/**/*.qrc', recursive=True):
                _rc_f = os.path.join(build_lib,
                                     os.path.relpath(_f.replace('.qrc', '_rc.py'),
                                                     SRC_PATH))
                os.makedirs(os.path.dirname(_rc_f), exist_ok=True)
                cmd_rc = cmd + ['-o', _rc_f, _f]
                subprocess.run(args=cmd_rc, check=True)
                if post_replace:
                    _replace_in_file(_rc_f, post_replace[0], post_replace[1])

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
        if _cs_opts['enable_gui']:
            self.run_command("build_cs_qt_files")
            self.run_command("build_cs_ui_files")
            self.run_command("build_cs_rc_files")

setup_options = {}
if _cs_opts.get('build_base'):
    setup_options['build'] = {'build_base': _cs_opts['build_base']}

setup(name='code_saturne',
      author='code_saturne dev team',
      author_email='saturne-support@edf.fr',
      description='Python CLI and GUI of code_saturne multiphysics CFD package',
      package_dir={'':SRC_PATH},
      package_data={'code_saturne.data.icons': ['COPYING.Tango', 'code_saturne.ico', '22x22/*', '32x32/*'],},
      options=setup_options,
      cmdclass={'build_py':build_py,
                'build_cs_qt_files':_build_qt_files,
                'build_cs_ui_files':_build_ui_files,
                'build_cs_rc_files':_build_rc_files,
                'egg_info':_cs_egg_info,
                'install_egg_info':_cs_install_egg_info},
      version=_cs_opts['version'],
      packages=find_namespace_packages(where=SRC_PATH, exclude=_cs_opts['exclude_dirs'] + ['build*'])
      )
