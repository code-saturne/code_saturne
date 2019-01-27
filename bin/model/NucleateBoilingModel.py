# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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


class NucleateBoilingModel(NonCondensableModel, Variables, Model):
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
        self.XMLnucleate    = self.XMLMassTrans.xmlInitNode('nucleate_boiling')


    def defaultValues(self):
        default = {}
        default['wallmodel']                = "nucleate_boiling"
        default['wallfunction']             = "standard"
        default['heatmodel']                = "extended_kurul-podowski"
        default['yplusmodel']               = "diameter"
        default['yplusvalue']               = 250.
        default['cavities_radius']          = 0.0001
        default['bubbles_diameter']         = 0.003
        default['oversaturate_temperature'] = 1.
        default['thermal_cond']             = 17.
        default['density']                  = 8000.
        default['cp']                       = 531.
        default['thicknessmodel']           = "off"
        default['thicknessvalue']           = 0.0001

        return default


    @Variables.noUndo
    def getHeatTransferModel(self):
        """
        get heat transfer model
        """
        model = self.XMLnucleate['model']
        if model == None :
            model = self.defaultValues()['heatmodel']
            self.setHeatTransferModel(model)
        return model


    @Variables.undoLocal
    def setHeatTransferModel(self, model):
        """
        put heat transfer model
        """
        self.isInList(model, ('standard_kurul-podowski','extended_kurul-podowski'))
        self.XMLnucleate['model'] = model
        if model == "standard_kurul-podowski" :
           self.XMLnucleate.xmlRemoveChild('cavities_radius')
           self.XMLnucleate.xmlRemoveChild('bubbles_diameter')


    @Variables.noUndo
    def getWallFunctionModel(self):
        """
        get wall function model
        """
        model = self.XMLnucleate['wallfunction']
        if model == None :
            model = self.defaultValues()['wallfunction']
            self.setWallFunctionModel(model)
        return model


    @Variables.undoLocal
    def setWallFunctionModel(self, model):
        """
        put wall function model
        """
        self.isInList(model, ('standard','koncar','mimouni'))
        self.XMLnucleate['wallfunction'] = model


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


    @Variables.noUndo
    def getMaxRadius(self):
        """
        get value for maximal cavity radius
        """
        value = self.XMLnucleate.xmlGetDouble('cavities_radius')
        if value == None :
           value = self.defaultValues()['cavities_radius']
           self.setMaxRadius(value)
        return value


    @Variables.undoLocal
    def setMaxRadius(self, value):
        """
        put value for maximal cavity radius
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('cavities_radius', value)


    @Variables.noUndo
    def getMaxDiameter(self):
        """
        get value for maximal bubble diameter
        """
        value = self.XMLnucleate.xmlGetDouble('bubbles_diameter')
        if value == None :
           value = self.defaultValues()['bubbles_diameter']
           self.setMaxDiameter(value)
        return value


    @Variables.undoLocal
    def setMaxDiameter(self, value):
        """
        put value for maximal bubble diameter
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('bubbles_diameter', value)


    @Variables.noUndo
    def getMaxOverSaturation(self):
        """
        get value for maximal over-saturation temperature
        """
        value = self.XMLnucleate.xmlGetDouble('oversaturate_temperature')
        if value == None :
           value = self.defaultValues()['oversaturate_temperature']
           self.setMaxOverSaturation(value)
        return value


    @Variables.undoLocal
    def setMaxOverSaturation(self, value):
        """
        put value for maximal over-saturation temperature
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('oversaturate_temperature', value)


    @Variables.noUndo
    def getThermalConductivity(self):
        """
        get value for conductivity
        """
        value = self.XMLnucleate.xmlGetDouble('thermal_cond')
        if value == None :
           value = self.defaultValues()['thermal_cond']
           self.setThermalConductivity(value)
        return value


    @Variables.undoLocal
    def setThermalConductivity(self, value):
        """
        put value for conductivity
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('thermal_cond', value)


    @Variables.noUndo
    def getDensity(self):
        """
        get value for density
        """
        value = self.XMLnucleate.xmlGetDouble('density')
        if value == None :
           value = self.defaultValues()['density']
           self.setDensity(value)
        return value


    @Variables.undoLocal
    def setDensity(self, value):
        """
        put value for density
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('density', value)


    @Variables.noUndo
    def getSpecificHeat(self):
        """
        get value for specific heat
        """
        value = self.XMLnucleate.xmlGetDouble('cp')
        if value == None :
           value = self.defaultValues()['cp']
           self.setSpecificHeat(value)
        return value


    @Variables.undoLocal
    def setSpecificHeat(self, value):
        """
        put value for specific heat
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('cp', value)


    @Variables.noUndo
    def getThicknessStatus(self):
        """
        get thickness status
        """
        ChildNode = self.XMLnucleate.xmlGetChildNode('thicknessmodel')
        if ChildNode == None:
           status = self.defaultValues()['thicknessmodel']
           self.setThicknessStatus(status)
        status = self.XMLnucleate.xmlGetString('thicknessmodel')
        return status


    @Variables.undoLocal
    def setThicknessStatus(self, status):
        """
        put thickness status
        """
        self.isOnOff(status)
        self.XMLnucleate.xmlSetData('thicknessmodel', status)
        if status == "off" :
           self.XMLnucleate.xmlRemoveChild('thicknessvalue')


    @Variables.noUndo
    def getThicknessValue(self):
        """
        get value for Thickness
        """
        value = self.XMLnucleate.xmlGetDouble('thicknessvalue')
        if value == None :
           value = self.defaultValues()['thicknessvalue']
           self.setThicknessValue(value)
        return value


    @Variables.undoLocal
    def setThicknessValue(self, value):
        """
        put value for Thickness
        """
        self.isGreater(value, 0.)
        self.XMLnucleate.xmlSetData('thicknessvalue', value)


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class WallTransferTestCase(ModelTest):
    """
    """
    def checkWallTransferInstantiation(self):
        """Check whether the WallTransferModel class could be instantiated"""
        model = None
        model = WallTransferModel(self.case)
        assert model != None, 'Could not instantiate WallTransferModel'

    def checkGetandSetWallTransferModel(self):
        """Check whether the WallTransferModel class could set and get WallTransferModel"""
        mdl = WallTransferModel(self.case)
        mdl.setWallTransferModel('nucleate_boiling')
        doc = '''<mass_transfer_model wallmodel="nucleate_boiling">
                         <nucleate_boiling/>
                 </mass_transfer_model>'''
        assert mdl.XMLMassTrans == self.xmlNodeFromString(doc),\
            'Could not set WallTransferModel'
        assert mdl.getWallTransferModel() == 'nucleate_boiling',\
            'Could not get WallTransferModel'


    def checkGetandSetHeatTransferModel(self):
        """Check whether the WallTransferModel class could set and get HeatTransferModel"""
        mdl = WallTransferModel(self.case)
        mdl.setHeatTransferModel('standard_kurul-podowski')
        doc = '''<nucleate_boiling model="standard_kurul-podowski"/>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set HeatTransferModel'
        assert mdl.getHeatTransferModel() == 'standard_kurul-podowski',\
            'Could not get HeatTransferModel'

    def checkGetandSetWallFunctionModel(self):
        """Check whether the WallTransferModel class could set and get WallFunctionModel"""
        mdl = WallTransferModel(self.case)
        mdl.setWallFunctionModel('standard')
        doc = '''<nucleate_boiling wallfunction="standard"/>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set WallFunctionModel'
        assert mdl.getWallFunctionModel() == 'standard',\
            'Could not get WallFunctionModel'

    def checkGetandSetYPlusModel(self):
        """Check whether the WallTransferModel class could set and get YPlusModel"""
        mdl = WallTransferModel(self.case)
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
        mdl = WallTransferModel(self.case)
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

    def checkGetandSetMaxRadius(self):
        """Check whether the WallTransferModel class could set and get MaxRadius"""
        mdl = WallTransferModel(self.case)
        mdl.setMaxRadius(0.01)
        doc = '''<nucleate_boiling>
                         <cavities_radius>
                                 0.01
                         </cavities_radius>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set MaxRadius'
        assert mdl.getMaxRadius() == 0.01,\
            'Could not get MaxRadius'

    def checkGetandSetMaxDiameter(self):
        """Check whether the WallTransferModel class could set and get MaxDiameter"""
        mdl = WallTransferModel(self.case)
        mdl.setMaxDiameter(0.02)
        doc = '''<nucleate_boiling>
                         <bubbles_diameter>
                                 0.02
                         </bubbles_diameter>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set MaxDiameter'
        assert mdl.getMaxDiameter() == 0.02,\
            'Could not get MaxDiameter'

    def checkGetandSetMaxOverSaturation(self):
        """Check whether the WallTransferModel class could set and get MaxOverSaturation"""
        mdl = WallTransferModel(self.case)
        mdl.setMaxOverSaturation(2)
        doc = '''<nucleate_boiling>
                         <oversaturate_temperature>
                                 2
                         </oversaturate_temperature>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set MaxOverSaturation'
        assert mdl.getMaxOverSaturation() == 2,\
            'Could not get MaxOverSaturation'

    def checkGetandSetThermalConductivity(self):
        """Check whether the WallTransferModel class could set and get ThermalConductivity"""
        mdl = WallTransferModel(self.case)
        mdl.setThermalConductivity(15)
        doc = '''<nucleate_boiling>
                         <thermal_cond>
                                 15
                         </thermal_cond>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set ThermalConductivity'
        assert mdl.getThermalConductivity() == 15,\
            'Could not get ThermalConductivity'

    def checkGetandSetDensity(self):
        """Check whether the WallTransferModel class could set and get Density"""
        mdl = WallTransferModel(self.case)
        mdl.setDensity(7000)
        doc = '''<nucleate_boiling>
                         <density>
                                 7000
                         </density>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set Density'
        assert mdl.getDensity() == 7000,\
            'Could not get Density'

    def checkGetandSetSpecificHeat(self):
        """Check whether the WallTransferModel class could set and get SpecificHeat"""
        mdl = WallTransferModel(self.case)
        mdl.setSpecificHeat(542)
        doc = '''<nucleate_boiling>
                         <cp>
                                 542
                         </cp>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set SpecificHeat'
        assert mdl.getSpecificHeat() == 542,\
            'Could not get SpecificHeat'


    def checkGetandSetThicknessStatus(self):
        """Check whether the WallTransferModel class could set and get ThicknessStatus"""
        mdl = WallTransferModel(self.case)
        mdl.setThicknessStatus('on')
        doc = '''<nucleate_boiling>
                         <thicknessmodel>
                                 on
                         </thicknessmodel>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set ThicknessStatus'
        assert mdl.getThicknessStatus() == 'on',\
            'Could not get ThicknessStatus'


    def checkGetandSetThicknessValue(self):
        """Check whether the WallTransferModel class could set and get ThicknessValue"""
        mdl = WallTransferModel(self.case)
        mdl.setThicknessValue(0.01)
        doc = '''<nucleate_boiling>
                         <thicknessvalue>
                                 0.01
                         </thicknessvalue>
                 </nucleate_boiling>'''
        assert mdl.XMLnucleate == self.xmlNodeFromString(doc),\
            'Could not set ThicknessValue'
        assert mdl.getThicknessValue() == 0.01,\
            'Could not get ThicknessValue'


def suite():
    testSuite = unittest.makeSuite(WallTransferTestCase, "check")
    return testSuite

def runTest():
    print("WallTransferTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
