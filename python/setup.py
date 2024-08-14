# ==============================================================================
# IMPORTS
# ==============================================================================

from setuptools import setup, find_packages, Command
from setuptools.command.egg_info import egg_info
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
# Specific command to build python files for GUI
# ==============================================

class build_py_from_ui(Command):
    """
    A custom command to build python files from .ui files.
    """

    user_options = []

    # ===============
    # Private methods
    # ===============

    def _copy_qt_files(self, build_base):
        """
        Copy corresponding Qt4 or Qt5 files.
        """

        for _f in glob.glob(os.path.join(SRC_PATH, 'code_saturne', 'gui', 'base', _cs_opts['qtpkg'], '*.py')):
            _abs_dst = os.path.join(SRC_PATH, 'code_saturne', 'gui', 'base', os.path.basename(_f))
            _nf = os.path.join(build_base,
                               'lib',
                               os.path.relpath(_abs_dst, SRC_PATH))

            shutil.copy2(_f, _nf, follow_symlinks=False)


    def _compile_cs_ui_files(self, build_base):
        """
        Compile ".ui" files into ".py" files.
        """
        for _f in glob.glob(SRC_PATH+'/**/*.ui', recursive=True):
            _ui_f = os.path.join(build_base,
                                 'lib',
                                 os.path.relpath(_f.replace('.ui','.py'),
                                                 SRC_PATH)
                                 )
            subprocess.run(args=[_cs_opts['pyuic'], '-o', _ui_f, _f],
                           check=True)

    def _compile_cs_rc_files(self, build_base):
        """
        Compile ".rc" files into ".py" files
        """
        for _f in glob.glob(SRC_PATH+'/**/*.qrc', recursive=True):
            _rc_f = os.path.join(build_base,
                                 'lib',
                                 os.path.relpath(_f.replace('.qrc','_rc.py'),
                                                 SRC_PATH)
                                 )
            _rc_dir = os.path.dirname(_f)
            subprocess.run(args=[_cs_opts['pyrcc'], '-o', _rc_f, _f],
                           check=True)


    # ==============
    # Public methods
    # ==============

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Get ui files
        if 'build' in self.distribution.command_obj:
            build_command = self.distribution.command_obj["build"]

            # Copy Qt files
            if _cs_opts['qtpkg']:
                self._copy_qt_files(build_command.build_base)

            if _cs_opts['pyuic']:
                self._compile_cs_ui_files(build_command.build_base)

            if _cs_opts['pyrcc']:
                self._compile_cs_rc_files(build_command.build_base)

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
        self.run_command("build_py_from_ui")

setup(name='code_saturne',
      author='code_saturne dev team',
      author_email='saturne-support@edf.fr',
      description='Python CLI and GUI of code_saturne multiphysics CFD package',
      package_dir={'':SRC_PATH},
      package_data={'code_saturne':['data/icons/*', 'data/icons/22x22/*', 'data/icons/32x32/*'],},
      cmdclass={'build_py_from_ui':build_py_from_ui,
                'build_py':build_py,},
      version=_cs_opts['version'],
      packages=find_packages(where=SRC_PATH, exclude=_cs_opts['exclude_dirs'])
      )
