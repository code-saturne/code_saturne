# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module defines a wrapper for the EOS library python API
"""

import unittest

class eosWrapper(object):

    def __init__(self):
        """
        Constructor
        """

        self._has_eos = True
        try:
            import eosAva
        except:
            self._has_eos = False
        else:
            import eosAva

        if self._has_eos:
            self._ava = eosAva.EosAvailable()

    def isActive(self):
        """
        Indicate if EOS library is available.
        """

        return self._has_eos

    def getListOfFluids(self):
        """
        Get list of available fluids.
        """

        retval = []
        if self._has_eos:
            retval = self._ava.whichFluids()

        return retval

    def getFluidMethods(self, fluid):
        """
        Get list of available tables for a given fluid.
        """

        retval = []
        if self._has_eos:
            self._ava.setMethods(fluid)

            retval =  self._ava.whichMethods()

        return retval

    def getFluidReferences(self, fluid, method):
        """
        Get list of available reference methods in a given table/method
        for a fluid liquid phase.
        """

        retval = []
        if self._has_eos:
            self._ava.setMethods(fluid)
            self._ava.setReferences(fluid, method)

            retval =  self._ava.whichReferences()

        return retval

    def getNumberOfFluidReferences(self, fluid, method):
        """
        Get number of available references for a fluid using a given method.
        """

        retval = len(self.getFluidReferences(fluid, method))

        return retval

    def getLiquidReferences(self, fluid, method):
        """
        Get list of available reference methods in a given table/method
        for a fluid liquid phase.
        """

        tmp = self.getFluidReferences(fluid, method)
        retval = [r for r in tmp if 'Liquid' in r]

        return retval

    def getVaporReferences(self, fluid, method):
        """
        Get list of available reference methods in a given table/method
        for a fluid vapor phase.
        """

        tmp = self.getFluidReferences(fluid, method)
        retval = [r for r in tmp if 'Vapor' in r]

        return retval

