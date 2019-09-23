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
from code_saturne.model.MainFieldsModel import *
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
from code_saturne.model.SpeciesModel import SpeciesModel
from code_saturne.model.NonCondensableModel import NonCondensableModel
from code_saturne.model.NotebookModel import NotebookModel


class MainFieldsInitializationModel(MainFieldsModel, Variables, Model):

    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case                = case
        self.XMLThermo           = self.case.xmlGetNode('thermophysical_models')
        self.XMLvariables        = self.XMLThermo.xmlGetNode('variables')
        self.XMLNonCondvariables = self.XMLThermo.xmlGetNode('non_condensable_list')
        self.XMLUserScalar       = self.case.xmlGetNode('additional_scalars')
        self.XMLScalar           = self.XMLUserScalar.xmlInitNode('scalars')


    def defaultValues(self):
        default = {}

        default['enthalpyModel']  = "temperature"
        return default


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())


    @Variables.undoLocal
    def setFormulaPressure(self, zone, str):
        """
        Gives a formula for pressure
        """
        self.__verifyZone(zone)
        node = self.XMLvariables.xmlGetNode('variable', field_id="none", name="pressure")
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        n.xmlSetData('formula', str)


    @Variables.noUndo
    def getFormulaPressure(self, zone):
        """
        Return a formula for pressure
        """
        self.__verifyZone(zone)
        node = self.XMLvariables.xmlGetNode('variable', field_id="none", name="pressure")
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        return n.xmlGetString('formula')

    @Variables.noUndo
    def getPressureFormulaComponents(self, zone):

        exp = self.getFormulaPressure(zone)
        req = [('pressure', 'pressure')]

        sym = [('x', "X cell's gravity center"),
               ('y', "Y cell's gravity center"),
               ('z', "Z cell's gravity center")]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setFormula(self, zone, fieldId, var_name, str):
        """
        Gives a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(fieldId, self.getFieldIdList())

        node = self.XMLvariables.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        n.xmlSetData('formula', str)


    @Variables.noUndo
    def getFormula(self, zone, fieldId, var_name):
        """
        Return a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLvariables.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        return n.xmlGetString('formula')


    @Variables.noUndo
    def getFormulaComponents(self, zone, fieldId, var_name):
        exp = self.getFormula(zone, fieldId, var_name)

        if var_name == 'velocity':
            req = [('u', 'Velocity along X'),
                   ('v', 'Velocity along Y'),
                   ('w', 'Velocity along Z')]

            sym = [('x', "X cell's gravity center"),
                   ('y', "Y cell's gravity center"),
                   ('z', "Z cell's gravity center")]

        elif var_name == 'volume_fraction':
            req = [('vol_f', 'volume fraction')]

            sym = [('x', "X cell's gravity center"),
                   ('y', "Y cell's gravity center"),
                   ('z', "Z cell's gravity center")]

        elif var_name == 'enthalpy':
            th_sca_label = self.getEnergyModel(zone, fieldId)
            req = [(th_sca_label, str(th_sca_label))]

            sym = [('x', "X cell's gravity center"),
                   ('y', "Y cell's gravity center"),
                   ('z', "Z cell's gravity center")]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setFormulaNonCondensable(self, zone, fieldId, var_name, str):
        """
        Gives a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(fieldId, self.getFieldIdList())

        node = self.XMLNonCondvariables.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        n.xmlSetData('formula', str)


    @Variables.noUndo
    def getFormulaNonCondensable(self, zone, fieldId, var_name):
        """
        Return a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLNonCondvariables.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        return n.xmlGetString('formula')


    @Variables.noUndo
    def getNonCondensableFormulaComponents(self, zone, fieldId, var_name):

        exp   = self.getFormulaNonCondensable(zone, fieldId, var_name)
        label = NonCondensableModel(self.case).getNonCondLabel(var_name)
        req   = [(label, str(label))]

        sym = [('x', "X cell's gravity center"),
               ('y', "Y cell's gravity center"),
               ('z', "Z cell's gravity center")]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym

    @Variables.undoLocal
    def setFormulaScalar(self, zone, fieldId, var_name, str):
        """
        Gives a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(fieldId, self.getFieldIdList())

        node = self.XMLScalar.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        n.xmlSetData('formula', str)


    @Variables.noUndo
    def getFormulaScalar(self, zone, fieldId, var_name):
        """
        Return a formula for initial values
        """
        self.__verifyZone(zone)
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLScalar.xmlGetNode('variable', field_id=fieldId, name=var_name)
        n = node.xmlInitChildNode('initial_value', zone_id=zone)
        return n.xmlGetString('formula')


    @Variables.noUndo
    def getScalarFormulaComponents(self, zone, fieldId, var_name):

        exp = self.getFormulaScalar(zone, fieldId, var_name)
        scalar_label = SpeciesModel(self.case).getScalarLabelByName(var_name)
        req = [(scalar_label, str(scalar_label))]

        sym = [('x', "X cell's gravity center"),
               ('y', "Y cell's gravity center"),
               ('z', "Z cell's gravity center")]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setEnergyModel(self, zone, fieldId, model) :
        """
        Public method.
        Set model for initialization of energy initialization
        """
        self.__verifyZone(zone)
        self.isInList(model, ('enthalpy', 'temperature', 'hsat_P'))
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLvariables.xmlGetNode('variable', field_id=fieldId, name="enthalpy")
        if not node:
            msg = "There is an error: this node " + str(node) + " should be existed"
            raise ValueError(msg)

        n = node.xmlInitChildNode('initial_type', zone_id=zone)
        n.xmlSetTextNode(model)

        if model == "hsat_P" :
            n = node.xmlGetNode('initial_value', zone_id=zone)
            if n :
                n.xmlRemoveNode()


    @Variables.noUndo
    def getEnergyModel(self, zone, fieldId) :
        """
        Public method.
        Return model for initialization of energy initialization
        """
        self.__verifyZone(zone)
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLvariables.xmlGetNode('variable', field_id=fieldId, name="enthalpy")
        if not node:
            msg = "There is an error: this node " + str(node) + " should be existed"
            raise ValueError(msg)

        nodem = node.xmlGetChildNode('initial_type', zone_id=zone)
        if nodem == None:
            model = self.defaultValues()['enthalpyModel']
            if ThermodynamicsModel(self.case).getMaterials(fieldId) == 'user_material' :
                model = 'enthalpy'
            self.setEnergyModel(zone, fieldId, model)
            nodem = node.xmlGetChildNode('initial_type', zone_id=zone)
        model = nodem.xmlGetTextNode()
        return model


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class MainFieldsInitializationTestCase(ModelTest):
    """
    """
    def checkMainFieldsInitializationInstantiation(self):
        """Check whether the MainFieldsInitialization class could be instantiated"""
        model = None
        model = MainFieldsInitializationModel(self.case)
        assert model != None, 'Could not instantiate MainFieldsInitiaziationModel'


    def checkGetandSetInitialVelocity(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get InitialVelocity"""
        MainFieldsModel(self.case).addField()
        mdl = MainFieldsInitializationModel(self.case)
        zone = '1'
        mdl.setInitialVelocity(zone, '1', 'VelocityX', 4.5)
        mdl.setInitialVelocity(zone, '1', 'VelocityY', -2.5)
        mdl.setInitialVelocity(zone, '1', 'VelocityZ', 8.5)
        doc = '''<variables>
                         <variable field_id="none" label="Pressure" name="pressure">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="enthalpy1" name="enthalpy">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="vol_f_1" name="volume_fraction">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="U1" name="VelocityX">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         4.5
                                 </initial_value>
                         </variable>
                         <variable field_id="1" label="V1" name="VelocityY">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         -2.5
                                 </initial_value>
                         </variable>
                         <variable field_id="1" label="W1" name="VelocityZ">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         8.5
                                 </initial_value>
                         </variable>
                 </variables>'''
        assert mdl.XMLvariables == self.xmlNodeFromString(doc),\
            'Could not set InitialVelocity'
        assert mdl.getInitialVelocity('1', '1') == [4.5, -2.5, 8.5],\
            'Could not get InitialVelocity'


    def checkGetandSetInitialFraction(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get InitialFraction"""
        MainFieldsModel(self.case).addField()
        mdl = MainFieldsInitializationModel(self.case)
        mdl.setInitialFraction('1','1',0.1)
        doc = '''<variables>
                         <variable field_id="none" label="Pressure" name="pressure">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="enthalpy1" name="enthalpy">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="vol_f_1" name="volume_fraction">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         0.1
                                 </initial_value>
                         </variable>
                         <variable field_id="1" label="U1" name="VelocityX">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="V1" name="VelocityY">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="W1" name="VelocityZ">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </variables>'''
        assert mdl.XMLvariables == self.xmlNodeFromString(doc),\
            'Could not set InitialFraction'
        assert mdl.getInitialFraction('1','1') == 0.1,\
            'Could not get InitialFraction'


    def checkGetandSetInitialEnergy(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get InitialEnergy"""
        MainFieldsModel(self.case).addField()
        mdl = MainFieldsInitializationModel(self.case)
        mdl.setInitialEnergy('1','1',10)
        doc = '''<variables>
                         <variable field_id="none" label="Pressure" name="pressure">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="enthalpy1" name="enthalpy">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         10
                                 </initial_value>
                         </variable>
                         <variable field_id="1" label="vol_f_1" name="volume_fraction">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="U1" name="VelocityX">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="V1" name="VelocityY">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="W1" name="VelocityZ">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </variables>'''
        assert mdl.XMLvariables == self.xmlNodeFromString(doc),\
            'Could not set InitialEnergy'
        assert mdl.getInitialEnergy('1','1') == 10,\
            'Could not get InitialEnergy'


    def checkGetandSetInitialNonCondensable(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get InitialNonCondensable"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField("2", "field2", 'dispersed', 'gas', 'on', 'on', 'off', 2)
        from code_saturne.model.NonCondensableModel import NonCondensableModel
        NonCondensableModel(self.case).addNonCondensable()
        mdl = MainFieldsInitializationModel(self.case)
        mdl.setInitialNonCondensable('1','2','mass_fraction_non_condensable_gas_0',0.23)
        doc = '''<non_condensable_list>
                         <variable field_id="2" label="XX_2" name="mass_fraction_non_condensable_gas_0">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         0.23
                                 </initial_value>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNonCondvariables == self.xmlNodeFromString(doc),\
            'Could not set InitialNonCondensable'
        assert mdl.getInitialNonCondensable('1','2','mass_fraction_non_condensable_gas_0') == 0.23,\
            'Could not get InitialNonCondensable'


    def checkGetandSetEnergyModel(self):
        """Check whether the MainFieldsInitiaziationModel class could set and get EnergyModel"""
        MainFieldsModel(self.case).addField()
        mdl = MainFieldsInitializationModel(self.case)
        mdl.setEnergyModel('1','1','temperature')
        doc = '''<variables>
                         <variable field_id="none" label="Pressure" name="pressure">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="enthalpy1" name="enthalpy">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_type zone_id="1">
                                         temperature
                                 </initial_type>
                         </variable>
                         <variable field_id="1" label="vol_f_1" name="volume_fraction">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="U1" name="VelocityX">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="V1" name="VelocityY">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                         <variable field_id="1" label="W1" name="VelocityZ">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </variables>'''
        assert mdl.XMLvariables == self.xmlNodeFromString(doc),\
            'Could not set EnergyModel'
        assert mdl.getEnergyModel('1','1') == 'temperature',\
            'Could not get EnergyModel'


    def checkGetandSetInitialValuePressure(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValuePressure"""
        MainFieldsModel(self.case).addField()
        mdl = MainFieldsInitializationModel(self.case)
        zone = '1'
        mdl.setInitialValuePressure(zone, 700000)
        doc = '''<variables>
                         <variable field_id="none" label="Pressure" name="pressure">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value zone_id="1">
                                         700000
                                 </initial_value>
                         </variable>
                 </variables>'''
        assert mdl.XMLvariables == self.xmlNodeFromString(doc),\
            'Could not set Initial Value pressure'
        assert mdl.getInitialValuePressure('1') == 700000,\
            'Could not get Initial Value pressure'



        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValuePressure(1.23)
        doc = '''<thermodynamics>
                         <reference_pressure>
                                 1.23
                         </reference_pressure>
                 </thermodynamics>'''
        assert mdl.getXMLThermo() == self.xmlNodeFromString(doc),\
            'Could not set InitialValuePressure'
        assert mdl.getInitialValuePressure() == 1.23,\
            'Could not get InitialValuePressure'


def suite():
    testSuite = unittest.makeSuite(MainFieldsInitializationTestCase, "check")
    return testSuite


def runTest():
    print("MainFieldsInitializationTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


