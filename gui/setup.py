#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This is the cs_gui setup.py script.

1) for installation:
python setup.py build
python setup.py install
python setup.py install --prefix MY_INSTALL_PATH
python setup.py install --prefix=MY_INSTALL_PATH

2) for distribution:
python setup.py bdist --format=wininst
python setup.py bdist --format=rpm
python setup.py sdist
python setup.py sdist --formats=gztar,zip
python setup.py bdist_rpm
python setup.py bdist_wininst

3) for more options:
python setup.py --help
python setup.py bdist --help-formats
"""

#-------------------------------------------------------------------------------
# Check versions
#-------------------------------------------------------------------------------

import sys, string, os

if not hasattr(sys, 'version_info') or sys.version_info < (2, 4, 0, 'final'):
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires python 2.4 or later."
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui  import *
except ImportError:
    print "\n  Error: Unable to import PyQt4.QtCore or PyQt4.QtGui modules."
    print "  Please check your PyQt4 installation.\n"
    sys.exit(0)


if map(int, string.split(QT_VERSION_STR, ".")) < [4, 3, 0]:
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires Qt 4.3 or later (found %s)." % QT_VERSION_STR


if map(int, string.split(PYQT_VERSION_STR, ".")) < [4, 3, 0]:
    raise SystemExit, "Graphical users interface of Code_Saturne "\
                      "requires PyQt 4.3 or later (found %s)." % PYQT_VERSION_STR

#-------------------------------------------------------------------------------
# Distribution or installation
#-------------------------------------------------------------------------------

from distutils.core import setup
from distutils.sysconfig import get_python_lib

install_dir = os.path.join(get_python_lib(), 'ncs')

for arg in sys.argv[1:]:
    if arg.startswith("--prefix="):
        install_dir = os.path.join(get_python_lib(prefix=arg[9:]), 'ncs')
        break

for i in range(len(sys.argv)):
    if sys.argv[i].startswith("--prefix") and not sys.argv[i].startswith("--prefix="):
        install_dir = os.path.join(get_python_lib(prefix=sys.argv[i+1]), 'ncs')
        break

setup(name='cs_gui',
      url='www.code-saturne.org',
      description='Graphical user interface of Code_Saturne CFD code',
      license='GNU GPL',
      maintainer='Code_Saturne team',
      maintainer_email='saturne-support@edf.fr',
      scripts=['cs_gui'],
      package_dir={ 'ncs' : '.', 'ncs.Base' : 'Base', 'ncs.Pages' : 'Pages'},
      packages=['ncs', 'ncs.Base', 'ncs.Pages'],
      data_files=[(install_dir,
                  ["AUTHORS", "COPYING"]),
                  (os.path.join(install_dir, "Base", "icons", "SplashScreen"),
                  [os.path.join("Base", "icons", "SplashScreen", "logocs.png")])]
      )

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
