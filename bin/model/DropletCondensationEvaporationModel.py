# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.NonCondensableModel import NonCondensableModel
import copy

#-------------------------------------------------------------------------------
# Constructor
#-------------------------------------------------------------------------------


class DropletCondensationEvaporationModel(NonCondensableModel, Variables, Model):
    """
    This class manages the wall tranfer model objects in the XML file
    """


    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        NonCondensableModel.__init__(self, case)
        self.case = case
        self.XMLNodeclosure = self.case.xmlGetNode('closure_modeling')
        self.XMLMassTrans   = self.XMLNodeclosure.xmlInitNode('mass_transfer_model')


    def defaultValues(self):
        default = {}
        default['yplusmodel']               = "diameter"
        default['yplusvalue']               = 250.

        return default


    @Variables.noUndo
    def getYPlusModel(self):
        """
        get Y+ status
        """
        ChildNode = self.XMLMassTrans.xmlGetNode('yplusmodel')
        if ChildNode == None:
           model = self.defaultValues()['yplusmodel']
           self.setYPlusModel(model)
        model = self.XMLMassTrans.xmlGetNode('yplusmodel')['model']
        return model


    @Variables.undoLocal
    def setYPlusModel(self, model):
        """
        put Y+ status
        """
        self.isInList(model, ('center','Yplus_value','diameter'))

        ChildNode = self.XMLMassTrans.xmlInitChildNode('yplusmodel')
        ChildNode.xmlSetAttribute(model = model)

        if model != "Yplus_value" :
           self.XMLMassTrans.xmlRemoveChild('yplusvalue')


    @Variables.noUndo
    def getYPlusValue(self):
        """
        get value for Y Plus
        """
        value = self.XMLMassTrans.xmlGetDouble('yplusvalue')
        if value == None :
           value = self.defaultValues()['yplusvalue']
           self.setYPlusValue(value)
        return value


    @Variables.undoLocal
    def setYPlusValue(self, value):
        """
        put value for Y Plus
        """
        self.isGreater(value, 0.)
        self.XMLMassTrans.xmlSetData('yplusvalue', value)


#-------------------------------------------------------------------------------
# DropletCondensationEvaporation test case
#-------------------------------------------------------------------------------

class DropletCondensationEvaporationTestCase(ModelTest):
    """
    """
    def checkDropletCondensationEvaporationInstantiation(self):
        """Check whether the DropletCondensationEvaporation class could be instantiated"""
        model = None
        model = DropletCondensationEvaporationModel(self.case)
        assert model != None, 'Could not instantiate DropletCondensationEvaporationModel'


    def checkGetandSetYPlusModel(self):
        """Check whether the WallTransferModel class could set and get YPlusModel"""
        mdl = DropletCondensationEvaporationModel(self.case)
        mdl.setYPlusModel('Yplus_value')
        doc = '''<mass_transfer_model>
                         <nucleate_boiling/>
                         <yplusmodel model="Yplus_value"/>
                 </mass_transfer_model>'''
        assert mdl.XMLMassTrans == self.xmlNodeFromString(doc),\
            'Could not set YPlusModel'
        assert mdl.getYPlusModel() == 'Yplus_value',\
            'Could not get YPlusModel'


    def checkGetandSetYPlusValue(self):
        """Check whether the WallTransferModel class could set and get YPlusValue"""
        mdl = DropletCondensationEvaporationModel(self.case)
        mdl.setYPlusModel('Yplus_value')
        mdl.setYPlusValue(260)
        doc = '''<mass_transfer_model>
                         <nucleate_boiling/>
                         <yplusmodel model="Yplus_value"/>
                         <yplusvalue>
                                 260
                         </yplusvalue>
                 </mass_transfer_model>'''
        assert mdl.XMLMassTrans == self.xmlNodeFromString(doc),\
            'Could not set YPlusValue'
        assert mdl.getYPlusValue() == 260,\
            'Could not get YPlusValue'


def suite():
    testSuite = unittest.makeSuite(WallTransferTestCase, "check")
    return testSuite

def runTest():
    print("WallTransferTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
