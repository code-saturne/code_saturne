# -*- coding: iso-8859-1 -*-
#
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
This module contains the following classes and function:
- BodyForcesModel
- BodyForcesView
- PageText
- BodyForcesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from math import sqrt
import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Body Force model class
#-------------------------------------------------------------------------------

class BodyForcesModel(Model):

    def __init__(self, case):

        """
        Constuctor.
        """
        self.case = case

        self.node_pp      = self.case.xmlGetNode('physical_properties')
        self.node_gravity = self.node_pp.xmlInitNode('gravity')

        self.nodes = ['gravity_x', 'gravity_y', 'gravity_z']


    def defaultBodyForcesValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['gravity'] = 0.0
        default['status']  = "off"

        return default


    def getGravity(self, var):
        """
        Return value of gravity for var
        """
        self.isInList(var, self.nodes)
        gravity = self.node_gravity.xmlGetDouble(var)
        if gravity == None:
            gravity = self.defaultBodyForcesValues()['gravity']
            self.setGravity(var, gravity)

        return gravity


    def setGravity(self, txml, value):
        """
        Put value of gravity for txml balise
        """
        self.isInList(txml, self.nodes)
        self.isFloat(value)
        self.node_gravity.xmlSetData(txml, value)

        if self.getGravity('gravity_x') == 0.0 and \
           self.getGravity('gravity_y') == 0.0 and \
           self.getGravity('gravity_z') == 0.0:
            from Pages.TimeStepModel import TimeStepModel
            TimeStepModel(self.case).RemoveThermalTimeStepNode()
            del TimeStepModel

    def getHydrostaticPressure(self):
        """
        Return status of hydrostatic pressure :
        'off' if standard, 'on' if improved
        """
        node = self.node_pp.xmlInitNode('hydrostatic_pressure', 'status')
        status = node['status']
        if not status:
            status = self.defaultBodyForcesValues()['status']
            self.setHydrostaticPressure(status)
        return status


    def setHydrostaticPressure(self, var):
        """
        Put status of hydrostatic pressure
        """
        self.isOnOff(var)
        node = self.node_pp.xmlInitNode('hydrostatic_pressure', 'status')
        node['status'] = var


#-------------------------------------------------------------------------------
# BodyForce test case
#-------------------------------------------------------------------------------


class BodyForcesTestCase(ModelTest):
    """
    """
    def checkBodyForcesInstantiation(self):
        """
        Check whether the BodyForcesModel class could be instantiated
        """
        model = None
        model = BodyForcesModel(self.case)
        assert model != None, 'Could not instantiate BodyForcesModel'

    def checkGetandSetGravity(self):
        """
        Check whether the gravity terms could be set and get
        """
        mdl = BodyForcesModel(self.case)
        mdl.setGravity('gravity_x', 1.1)
        mdl.setGravity('gravity_y', 2.2)
        mdl.setGravity('gravity_z', 3.3)
        doc = '''<gravity>
                    <gravity_x>1.1</gravity_x>
                    <gravity_y>2.2</gravity_y>
                    <gravity_z>3.3</gravity_z>
                </gravity>'''
        assert mdl.node_gravity == self.xmlNodeFromString(doc), \
                                        'Could not set gravity'
        assert mdl.getGravity('gravity_y') == 2.2, 'Could not get gravity'

    def checkGetandSetHydrostaticPressure(self):
        """
        Check whether the hydrostatic pressure could be set and get
        """
        mdl = BodyForcesModel(self.case)
        mdl.setHydrostaticPressure('on')
        n = mdl.node_pp.xmlGetNode('hydrostatic_pressure')
        doc = '''<hydrostatic_pressure status="on"/>'''
        assert n == self.xmlNodeFromString(doc), \
                    'Could not set hydrostatic pressure'
        assert mdl.getHydrostaticPressure() == 'on',\
                                'Could not get hydrostatic pressure'


def suite():
    testSuite = unittest.makeSuite(BodyForcesTestCase, "check")
    return testSuite


def runTest():
    print("BodyForcesTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
