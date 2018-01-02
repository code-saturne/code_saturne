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

"""
This module defines the HeadLosses model data management.

This module contains the following classes and function:
- HeadLossesModel
- HeadLossesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Pages.LocalizationModel import LocalizationModel, VolumicLocalizationModel, Zone

#-------------------------------------------------------------------------------
# HeadLosses model class
#-------------------------------------------------------------------------------

class HeadLossesModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about HeadLosses
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_domain  = self.case.xmlGetNode('solution_domain')
        self.node_volzone = self.node_domain.xmlGetNode('volumic_conditions')
        self.node_hloss   = self.node_models.xmlInitNode('head_losses')

        self.coeffNames = ('kxx', 'kyy', 'kzz')
        self.matrix = ('a11', 'a12', 'a13',
                       'a21', 'a22', 'a23',
                       'a31', 'a32', 'a33')
        self.choicevalue = ('choice')

        self.getNameAndLocalizationZone()



    def __defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['kxx']     = 0.0
        default['kyy']     = 0.0
        default['kzz']     = 0.0
        default['a11']     = 1.0
        default['a12']     = 0.0
        default['a13']     = 0.0
        default['a21']     = 0.0
        default['a22']     = 1.0
        default['a23']     = 0.0
        default['a31']     = 0.0
        default['a32']     = 0.0
        default['a33']     = 1.0
        default['choice'] = 'off'
        return default


    @Variables.noUndo
    def getNameAndLocalizationZone(self):
        """
        Return name and localization zone from volume regions definitions.
        """
        zoneDico = {}
        zonesList = LocalizationModel('VolumicZone', self.case).getZones()
        for zone in zonesList:
            if zone.getNature()['head_losses'] == 'on':
                label = zone.getLabel()
                zoneid = zone.getCodeNumber()
                localization = zone.getLocalization()
                zoneDico[label] = (zoneid, localization)
                self.setNameAndLabelZone(zoneid)

        return zoneDico


    @Variables.undoGlobal
    def setNameAndLabelZone(self, zoneid):
        """
        Set name and label zone for head losses markups.
        """
        self.node_hloss.xmlInitChildNode('head_loss', zone_id=zoneid)
        self.getKCoefficients(zoneid)
        self.getMatrix(zoneid)
        self.getMatrixChoice(zoneid,'choice')


    @Variables.noUndo
    def getMatrixChoice(self,zoneid,choice):
        """
        Get the Transfo Matrix choice
        """
        self.isInList(choice, self.choicevalue)
        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        value = node.xmlGetString(choice)
        if value == None:
            value = self.__defaultValues()[choice]
            self.setMatrixChoice(zoneid, choice, value)
        return value


    @Variables.undoLocal
    def setMatrixChoice(self, zoneid, choice, value):
        """
        Set the Transfo Matrix Choice
        """
        self.isInt(int(zoneid))
        self.isInList(choice, self.choicevalue)

        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        node.xmlSetData(choice, value)


    @Variables.noUndo
    def getCoefficient(self, zoneid, k):
        """
        Return value of coefficient k for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        self.isInList(k, self.coeffNames)

        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        value = node.xmlGetDouble(k)
        if value == None:
            value = self.__defaultValues()[k]
            self.setCoefficient(zoneid, k, value)

        return value


    @Variables.noUndo
    def getKCoefficients(self, zoneid):
        """
        Get value of kxx, kyy and kzz from xml file, for the head loss with zone's id.
        """
        self.isInt(int(zoneid))

        kxx = self.getCoefficient(zoneid, 'kxx')
        kyy = self.getCoefficient(zoneid, 'kyy')
        kzz = self.getCoefficient(zoneid, 'kzz')

        return kxx, kyy, kzz


    @Variables.undoLocal
    def setCoefficient(self, zoneid, k, value):
        """
        Set value of coefficient k for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        self.isInList(k, self.coeffNames)
        self.isFloat(value)

        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        node.xmlSetData(k, value)


    @Variables.undoGlobal
    def setKCoefficients(self, zoneid, kxx, kyy, kzz):
        """
        Set value of kxx, kyy and kzz into xml file, for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        self.isFloat(kxx)
        self.isFloat(kyy)
        self.isFloat(kzz)

        self.setCoefficient(zoneid, 'kxx', kxx)
        self.setCoefficient(zoneid, 'kyy', kyy)
        self.setCoefficient(zoneid, 'kzz', kzz)


    @Variables.noUndo
    def getMatrixComposant(self, zoneid, a):
        """
        Get values of one composant of the matrix of the change reference frame,
        for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        self.isInList(a, self.matrix)

        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        value = node.xmlGetDouble(a)
        if value == None:
            value = self.__defaultValues()[a]
            self.setMatrixComposant(zoneid, a, value)

        return value


    @Variables.noUndo
    def getMatrix(self, zoneid):
        """
        Get values of matrix of the change reference frame from xml file,
        for the head loss with zone's id.
        """
        self.isInt(int(zoneid))

        a11 = self.getMatrixComposant(zoneid, 'a11')
        a12 = self.getMatrixComposant(zoneid, 'a12')
        a13 = self.getMatrixComposant(zoneid, 'a13')
        a21 = self.getMatrixComposant(zoneid, 'a21')
        a22 = self.getMatrixComposant(zoneid, 'a22')
        a23 = self.getMatrixComposant(zoneid, 'a23')
        a31 = self.getMatrixComposant(zoneid, 'a31')
        a32 = self.getMatrixComposant(zoneid, 'a32')
        a33 = self.getMatrixComposant(zoneid, 'a33')

        return a11, a12, a13, a21, a22, a23, a31, a32, a33


    @Variables.undoLocal
    def setMatrixComposant(self, zoneid, a, value):
        """
        Set value of composant of matrix of the change reference frame,
        for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        self.isInList(a, self.matrix)
        self.isFloat(value)

        node = self.node_hloss.xmlGetNode('head_loss', zone_id=zoneid)
        node.xmlSetData(a, value)


    @Variables.undoGlobal
    def setMatrix(self, zoneid, a11, a12, a13, a21, a22, a23, a31, a32, a33):
        """
        Set values of the matrix of the change reference frame,
        for the head loss with zone's id.
        """
        self.isInt(int(zoneid))
        for a in (a11, a12, a13, a21, a22, a23, a31, a32, a33):
            self.isFloat(a)

        self.setMatrixComposant(zoneid, 'a11', a11)
        self.setMatrixComposant(zoneid, 'a12', a12)
        self.setMatrixComposant(zoneid, 'a13', a13)
        self.setMatrixComposant(zoneid, 'a21', a21)
        self.setMatrixComposant(zoneid, 'a22', a22)
        self.setMatrixComposant(zoneid, 'a23', a23)
        self.setMatrixComposant(zoneid, 'a31', a31)
        self.setMatrixComposant(zoneid, 'a32', a32)
        self.setMatrixComposant(zoneid, 'a33', a33)


#-------------------------------------------------------------------------------
# HeadLossesModel test case
#-------------------------------------------------------------------------------


class HeadLossesModelTestCase(ModelTest):
    """
    """
    def checkHeadLossesInstantiation(self):
        """Check whether the HeadLossesModel class could be instantiated"""
        model = None
        model = HeadLossesModel(self.case)
        assert model != None, 'Could not instantiate HeadLossesModel'


    def checkSetandGetNameLabelandCoefficients(self):
        """Check whether the head_losses markups could be set and get"""
        # we create a new zone for head losses
        loc = LocalizationModel("VolumicZone", self.case)
        zone = Zone("VolumicZone", label='toto', localization="1 or door", nature="head_losses")
        loc.addZone(zone)

        #we can test
        mdl = HeadLossesModel(self.case)
        mdl.getNameAndLocalizationZone()
        doc = '''<head_losses>
                    <head_loss zone_id="2">
                        <kxx>0</kxx>
                        <kyy>0</kyy>
                        <kzz>0</kzz>
                        <a11>0</a11>
                        <a12>0</a12>
                        <a13>0</a13>
                        <a21>0</a21>
                        <a22>0</a22>
                        <a23>0</a23>
                        <a31>0</a31>
                        <a32>0</a32>
                        <a33>0</a33>
                    </head_loss>
                </head_losses>'''
        assert mdl.node_hloss == self.xmlNodeFromString(doc),\
            'Could not set zone_id and label for head losses'
        assert mdl.getNameAndLocalizationZone() == {'toto': (2, '1 or door')},\
            'Could not get zone_id, label and localization for head losses'

    def checkSetandGetKCoefficients(self):
        """Check whether the head_losses could be set and get kxx, kyy, kzz"""
        # we create a new zone for head losses
        loc = LocalizationModel("VolumicZone", self.case)
        zone = Zone("VolumicZone", label='toto', localization="1 or door", nature="head_losses")
        loc.addZone(zone)

        mdl = HeadLossesModel(self.case)
        mdl.setKCoefficients('2', 10., 100., 1000.)
        doc = '''<head_losses>
                    <head_loss zone_id="2">
                        <kxx>10</kxx>
                        <kyy>100</kyy>
                        <kzz>1000</kzz>
                        <a11>0</a11>
                        <a12>0</a12>
                        <a13>0</a13>
                        <a21>0</a21>
                        <a22>0</a22>
                        <a23>0</a23>
                        <a31>0</a31>
                        <a32>0</a32>
                        <a33>0</a33>
                    </head_loss>
                </head_losses>'''

        assert mdl.node_hloss == self.xmlNodeFromString(doc),\
            'Could not set kxx, kyy, kzz coefficients for head losses'
        assert mdl.getKCoefficients('2') == (10., 100., 1000.),\
            'Could not get kxx, kyy, kzz coefficients for head losses'

        mdl.setCoefficient('2', 'kyy', 555.)
        doc2 = '''<head_losses>
                    <head_loss zone_id="2">
                        <kxx>10</kxx>
                        <kyy>555</kyy>
                        <kzz>1000</kzz>
                        <a11>0</a11>
                        <a12>0</a12>
                        <a13>0</a13>
                        <a21>0</a21>
                        <a22>0</a22>
                        <a23>0</a23>
                        <a31>0</a31>
                        <a32>0</a32>
                        <a33>0</a33>
                    </head_loss>
                </head_losses>'''

        assert mdl.node_hloss == self.xmlNodeFromString(doc2),\
            'Could not set one coefficient for head losses'
        assert mdl.getCoefficient('2', 'kyy') == 555.,\
            'Could not get one coefficient for head losses'

    def checkSetandGetMatrix(self):
        """Check whether the head_losses could be set and get composantes of matrix"""
        # we create a new zone for head losses
        loc = LocalizationModel("VolumicZone", self.case)
        zone = Zone("VolumicZone", label='toto', localization="1 or door", nature="head_losses")
        loc.addZone(zone)

        mdl = HeadLossesModel(self.case)
        mdl.setMatrix('2', 1., 1.2, 1.5, 2., 2.2, 2.5, 3., 3.2, 3.5)
        doc = '''<head_losses>
                    <head_loss zone_id="2">
                        <kxx>0</kxx>
                        <kyy>0</kyy>
                        <kzz>0</kzz>
                        <a11>1</a11>
                        <a12>1.2</a12>
                        <a13>1.5</a13>
                        <a21>2</a21>
                        <a22>2.2</a22>
                        <a23>2.5</a23>
                        <a31>3</a31>
                        <a32>3.2</a32>
                        <a33>3.5</a33>
                    </head_loss>
                </head_losses>'''

        assert mdl.node_hloss == self.xmlNodeFromString(doc),\
            'Could not set matrix for head losses'
        assert mdl.getMatrix('2') == (1., 1.2, 1.5, 2., 2.2, 2.5, 3., 3.2, 3.5),\
            'Could not get matrix for head losses'

        mdl.setMatrixComposant('2', 'a23', 2300.55)
        doc2='''<head_losses>
                    <head_loss zone_id="2">
                        <kxx>0</kxx>
                        <kyy>0</kyy>
                        <kzz>0</kzz>
                        <a11>1</a11>
                        <a12>1.2</a12>
                        <a13>1.5</a13>
                        <a21>2</a21>
                        <a22>2.2</a22>
                        <a23>2300.55</a23>
                        <a31>3</a31>
                        <a32>3.2</a32>
                        <a33>3.5</a33>
                    </head_loss>
                </head_losses>'''
        assert mdl.node_hloss == self.xmlNodeFromString(doc2),\
            'Could not set one composant of the matrix for head losses'
        assert mdl.getMatrixComposant('2', 'a23') == 2300.55,\
            'Could not get one composant of the matrix for head losses'


def suite():
    testSuite = unittest.makeSuite(HeadLossesModelTestCase, "check")
    return testSuite

def runTest():
    print("HeadLossesModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
