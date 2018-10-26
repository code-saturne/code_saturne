# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

import sys, unittest
from code_saturne.Base.XMLvariables import Model
from code_saturne.Base.XMLengine import *
from code_saturne.Base.XMLmodelNeptune import *
from code_saturne.Pages.OutputControlModel import OutputControlModel as SaturneOutputControlModel


class OutputControlModel(SaturneOutputControlModel):
    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        SaturneOutputControlModel.__init__(self, case)


    @Variables.undoGlobal
    def deleteMonitoringPoint(self, num):
        """
        Public method.
        Delete a single monitoring point.
        @type num: C{String}
        @param num: identifier of the monitoring point
        """
        self.isStr(num)
        self.isGreater(int(num), 0)
        self.isLowerOrEqual(int(num), self.getNumberOfMonitoringPoints())

        # delete the node of the monitoring point

        node = self.node_out.xmlGetNode('probe', name=num)
        if node:
            node.xmlRemoveNode()
            self.case.xmlRemoveChild('probe_recording', name=num)

            # renumerotation of all monitoring points

            for p in range(int(num)+1, self.getNumberOfMonitoringPoints()+2):
                probe = self.node_out.xmlGetNode('probe', name=p)
                probe['name'] = p - 1
                for probe_recording in self.case.xmlGetNodeList('probe_recording', name=p):
                    probe_recording['name'] = p - 1

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class OutputControlTestCase(ModelTest):
    """
    """
    def checkOutputControlInstantiation(self):
        """Check whether the OutputControl class could be instantiated"""
        model = None
        model = OutputControlModel(self.case)
        assert model != None, 'Could not instantiate OutputControlModel'

def suite():
    testSuite = unittest.makeSuite(OutputControlTestCase, "check")
    return testSuite

def runTest():
    print("OutputControlCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
