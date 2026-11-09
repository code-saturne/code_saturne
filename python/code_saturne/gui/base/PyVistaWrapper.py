# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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

"""
This module defines a wrapper for the PyVista and PyVistaQt library python API.
"""

import logging
import unittest

log = logging.getLogger("PyVistaWrapper")

class pyvistaWrapper(object):
    """
    Wrapper class to check for the presence of PyVista and PyVistaQt
    and provide access to their visualization utilities.
    """

    def __init__(self):
        """
        Constructor.
        """
        self._has_pyvista = True
        self._has_pyvistaqt = True
        self._pv = None
        self._qt_interactor = None

        try:
            import pyvista as pv
            self._pv = pv
        except Exception as e:
            self._has_pyvista = False
            log.debug("PyVista not available: %s", e)

        try:
            from pyvistaqt import QtInteractor
            self._qt_interactor = QtInteractor
        except Exception as e:
            self._has_pyvistaqt = False
            log.debug("PyVistaQt not available: %s", e)

    def isActive(self):
        """
        Indicate if both PyVista and PyVistaQt libraries are available
        for 3D GUI visualization.
        """
        return self._has_pyvista and self._has_pyvistaqt

    def hasPyVista(self):
        """
        Indicate if PyVista library is available.
        """
        return self._has_pyvista

    def hasPyVistaQt(self):
        """
        Indicate if PyVistaQt library is available.
        """
        return self._has_pyvistaqt

    def getVersion(self):
        """
        Get PyVista version string, or None if unavailable.
        """
        if self._has_pyvista and self._pv is not None:
            return getattr(self._pv, '__version__', None)
        return None

    def getPyVista(self):
        """
        Get imported pyvista module if available, or None.
        """
        return self._pv

    def getQtInteractor(self):
        """
        Get QtInteractor class if available, or None.
        """
        return self._qt_interactor


#-------------------------------------------------------------------------------
# Unit tests
#-------------------------------------------------------------------------------

class PyVistaWrapperTestCase(unittest.TestCase):
    def testPyVistaWrapper(self):
        wrapper = pyvistaWrapper()
        self.assertIsInstance(wrapper.isActive(), bool)
        self.assertIsInstance(wrapper.hasPyVista(), bool)
        self.assertIsInstance(wrapper.hasPyVistaQt(), bool)
        if wrapper.isActive():
            self.assertIsNotNone(wrapper.getPyVista())
            self.assertIsNotNone(wrapper.getQtInteractor())
            self.assertIsNotNone(wrapper.getVersion())


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(PyVistaWrapperTestCase))
    return testSuite


def runTest():
    runner = unittest.TextTestRunner()
    runner.run(suite())


if __name__ == '__main__':
    runTest()
