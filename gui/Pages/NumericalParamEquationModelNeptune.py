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

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string
import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.XMLvariables import Model
from code_saturne.Base.XMLengine import *
from code_saturne.Base.XMLmodelNeptune import *


class NumericalParamEquatModel(Model):

    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case) :
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case            = case
        self.XMLThermo       = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeVariable = self.XMLThermo.xmlGetNode('variables')
        self.XMLClosure      = self.case.xmlGetNode('closure_modeling')
        self.XMLTurbulence   = self.XMLClosure.xmlInitNode('turbulence')
        self.XMLTurbVariable = self.XMLTurbulence.xmlInitNode('variables')
        self.XMLAddScalar    = self.case.xmlInitNode('additional_scalars')
        self.XMLScalar       = self.XMLAddScalar.xmlInitNode('scalars')
        self.XMLnonCond      = self.XMLThermo.xmlInitNode('non_condensable_list')

        self.imgra_model = ['iterative_handling', \
                            'sweep_linear', \
                            'sweep_linear_ext_neigh', \
                            'least_square_linear', \
                            'least_square_linear_ext_neigh', \
                            'iterative_handling_weighted', \
                            'least_square_linear_weighted', \
                            'least_square_linear_ext_neigh_weighted']


    def defaultValues(self) :
        default = {}

        default['max_iter_number']                  = 10000
        default['solver_precision']                 = 1e-5
        default['solver_precision_alpha']           = 1e-8
        default['solver_precision_enthal']          = 1e-7
        default['slope_test']                       = 'on'
        default['flux_reconstruction']              = 'off'
        default['order_scheme']                     = 'solu'
        default['order_scheme_pressure_turbulence'] = 'upwind'
        default['solver']                           = 'automatic'
        default['gradient_reconstruction']          = 'sweep_linear_ext_neigh'
        return default


    def getVariableList(self) :
        """
        return list of variables
        """
        list = []
        from MainFieldsModel import MainFieldsModel
        for node in self.XMLNodeVariable.xmlGetNodeList('variable') :
            if self._isPressure(node) != 1 :
                # control to add enthalpy only if solved!!!
                if self._isEnthalpy(node) == 1 :
                    field = node['field_id']
                    if MainFieldsModel(self.case).getEnergyResolution(field) == "on":
                        list.append(node['label'])
                else:
                    list.append(node['label'])

        if self.XMLTurbVariable != None:
            for node in self.XMLTurbVariable.xmlGetNodeList('variable') :
                list.append(node['label'])

        for node in self.XMLScalar.xmlGetNodeList('variable') :
            list.append(node['label'])

        for node in self.XMLnonCond.xmlGetNodeList('variable') :
            list.append(node['label'])

        del MainFieldsModel
        return list


    def _isPressure(self, node) :
        """
        Return : 1 if name of node is 'pressure', 0 if not
        """
        if node and node['name'] == 'pressure':
            return 1
        else:
            return 0


    def _isTurbulenceVariable(self, node) :
        """
        Return : 1 if name of node is a turbulence variable, 0 if not
        """
        if node and node['name'] in TurbulenceModelsDescribing.turbulenceVariables['all'] :
            return 1
        else:
            return 0


    def _isAlpha(self, node) :
        """
        Return : 1 if name of node is 'alpha', 0 if not
        """
        if node and node['name'] == 'volume_fraction':
            return 1
        else:
            return 0


    def _isEnthalpy(self, node):
        """
        Return : 1 if name of node is 'enthalpy', 0 if not
        """
        if node and node['name'] == 'enthalpy':
            return 1
        else:
            return 0


    def _getSchemeLabelNode(self, label):
        """ Private method: return node called with label'label' for scheme nodes"""
        for node in self.XMLNodeVariable.xmlGetNodeList('variable') :
            if node['label'] == label:
                if node['label'] == 'Pressure':
                    raise ValueError("This method does not run with pressure")
                else:
                    return node
        for node in self.XMLTurbVariable.xmlGetNodeList('variable') :
            if node['label'] == label:
                return node

        for node in self.XMLScalar.xmlGetNodeList('variable') :
            if node['label'] == label:
                return node

        for node in self.XMLnonCond.xmlGetNodeList('variable') :
            if node['label'] == label:
                return node

        raise ValueError("This label does not exist: " + label)


    @Variables.undoGlobal
    def setSchemeModel(self, label, value) :
        """
        Put value of order scheme for variable or scalar labelled label
        only if it 's different of default value
        """
        self.isInList(value, ('upwind', 'centered', 'solu'))
        node = self._getSchemeLabelNode(label)
        defval = self.defaultValues()['order_scheme']
        if (self._isTurbulenceVariable(node)) :
            defval = self.defaultValues()['order_scheme_pressure_turbulence']

        if value == defval :
            node.xmlRemoveChild('scheme')
        else:
            n = node.xmlInitNode('scheme')
            n['choice'] = value

        if value == "upwind":
            self.setSlopeTestStatus(label, "off")
        else:
            self.setSlopeTestStatus(label, "on")


    @Variables.noUndo
    def getSchemeModel(self, label) :
        """
        Return value of order scheme for variable labelled label
        """
        node = self._getSchemeLabelNode(label)
        value = self.defaultValues()['order_scheme']
        if (self._isTurbulenceVariable(node)) :
            value = self.defaultValues()['order_scheme_pressure_turbulence']
        n = node.xmlGetNode('scheme')
        if n:
            value = n['choice']
        return value


    @Variables.undoLocal
    def setSlopeTestStatus(self, label, status) :
        """
        Put status of slope test for variable labelled label
        """
        self.isOnOff(status)
        node = self._getSchemeLabelNode(label)
        if status == self.defaultValues()['slope_test']:
            node.xmlRemoveChild('slope_test')
        else:
            n = node.xmlInitNode('slope_test')
            n['status'] = status


    @Variables.noUndo
    def getSlopeTestStatus(self, label) :
        """
        Return value of slope test for variable labelled label
        """
        node = self._getSchemeLabelNode(label)
        value = self.defaultValues()['slope_test']
        n = node.xmlGetNode('slope_test')
        if n:
            value = n['status']
        return value


    @Variables.undoLocal
    def setSolverModel(self, label, value) :
        """
        Put value of solver for variable or scalar labelled label
        """
        self.isInList(value, ('automatic', 'jacobi', 'pcg', 'cgstab', 'jacobi_saturne', 'pcg_saturne', 'bicgstab_saturne', 'bicgstab2_saturne', 'gmres_saturne', 'gauss_seidel_saturne', 'sym_gauss_seidel_saturne', 'pcr3_saturne'))
        node = self._getSchemeLabelNode(label)
        default = self.defaultValues()['solver']
        if value == default:
            node.xmlRemoveChild('solver')
        else:
            n = node.xmlInitNode('solver')
            n['choice'] = value


    @Variables.noUndo
    def getSolverModel(self, label) :
        """
        Return value of solver for variable labelled label
        """
        node = self._getSchemeLabelNode(label)
        value = self.defaultValues()['solver']
        n = node.xmlGetNode('solver')
        if n:
            value = n['choice']
        return value


    @Variables.undoLocal
    def setSolverPrecision(self, label, value) :
        """
        Put value of solveur precision for variable labelled label
        """
        self.isPositiveFloat(value)
        node = self._getSchemeLabelNode(label)
        default = self.defaultValues()['solver_precision']
        if self._isAlpha(node):
            default = self.defaultValues()['solver_precision_alpha']
        elif self._isEnthalpy(node):
            default = self.defaultValues()['solver_precision_enthal']

        if value != default:
            node.xmlSetData('solver_precision', value)
        else:
            node.xmlRemoveChild('solver_precision')


    @Variables.noUndo
    def getSolverPrecision(self, label) :
        """
        Return value of solveur precision for variable labelled label
        """
        node = self._getSchemeLabelNode(label)

        default = self.defaultValues()['solver_precision']
        if self._isAlpha(node):
            default = self.defaultValues()['solver_precision_alpha']
        elif self._isEnthalpy(node):
            default = self.defaultValues()['solver_precision_enthal']

        value = node.xmlGetDouble('solver_precision')
        if value == None:
            value = default
        return value


    @Variables.undoLocal
    def setMaximumIteration(self, label, value) :
        """
        Put number of maximum iterations for variable labelled label
        """
        self.isInt(value)
        node = self._getSchemeLabelNode(label)
        default = self.defaultValues()['max_iter_number']
        if value != default :
            node.xmlSetData('max_iter_number', value)
        else:
            node.xmlRemoveChild('max_iter_number')


    @Variables.noUndo
    def getMaximumIteration(self, label) :
        """
        Return number of maximum iterations for variable labelled label
        """
        node = self._getSchemeLabelNode(label)
        value = node.xmlGetInt('max_iter_number')
        if value == None:
            value = self.defaultValues()['max_iter_number']
        return value


    def getPressureNodefields(self):
        for node in self.XMLNodeVariable.xmlGetNodeList('variable') :
            if node['label'] == 'Pressure':
                return node

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class NumericalParamEquatTestCase(ModelTest):
    """
    """
    def checkNumericalParamEquatInstantiation(self):
        """Check whether the NumericalParamEquationModel class could be instantiated"""
        model = None
        model = NumericalParamEquatModel(self.case)
        assert model != None, 'Could not instantiate NumericalParamEquationModel'


    def checkGetVariableList(self):
        """Check whether the NumericalParamEquationModel class could get the VariableList"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        assert mdl.getVariableList() == ['enthalpy1', 'alpha1', 'U1', 'V1', 'W1'],\
            'Could not get VariableList'


    def checkGetandSetSchemeModel(self):
        """Check whether the NumericalParamEquationModel class could be set and get SchemeModel"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        mdl.setSchemeModel('pressure','solu')
        doc = '''<variable field_id="none" label="Pressure" name="pressure">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <scheme choice="solu"/>
                 </variable>'''
        assert mdl.getPressureNodefields() == self.xmlNodeFromString(doc),\
            'Could not set SchemeModel'
        assert mdl.getSchemeModel('Pressure') == 'solu',\
            'Could not get SchemeModel'


    def checkGetandSetSlopeTestStatus(self):
        """Check whether the NumericalParamEquationModel class could be set and get SlopeTestStatus"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        mdl.setSlopeTestStatus('Pressure','off')
        doc = '''<variable field_id="none" label="Pressure" name="pressure">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <slope_test status="off"/>
                 </variable>'''
        assert mdl.getPressureNodefields() == self.xmlNodeFromString(doc),\
            'Could not set SlopeTestStatus'
        assert mdl.getSlopeTestStatus('Pressure') == 'off',\
            'Could not get SlopeTestStatus'


    def checkGetandSetSolverModel(self):
        """Check whether the NumericalParamEquationModel class could be set and get SolverModel"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        mdl.setSolverModel('Pressure','jacobi')
        doc = '''<variable field_id="none" label="Pressure" name="pressure">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <solver choice="jacobi"/>
                 </variable>'''
        assert mdl.getPressureNodefields() == self.xmlNodeFromString(doc),\
            'Could not set SolverModel'
        assert mdl.getSolverModel('Pressure') == 'jacobi',\
            'Could not get SolverModel'


    def checkGetandSetSolverPrecision(self):
        """Check whether the NumericalParamEquationModel class could be set and get SolverPrecision"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        mdl.setSolverPrecision('Pressure',56.23)
        doc = '''<variable field_id="none" label="Pressure" name="pressure">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <solver_precision>
                                 56.23
                         </solver_precision>
                 </variable>'''
        assert mdl.getPressureNodefields() == self.xmlNodeFromString(doc),\
            'Could not set SolverPrecision'
        assert mdl.getSolverPrecision('Pressure') == 56.23,\
            'Could not get SolverPrecision'


    def checkGetandSetMaximumIteration(self):
        """Check whether the NumericalParamEquationModel class could be set and get MaximumIteration"""
        from MainFieldsModel import MainFieldsModel
        MainFieldsModel(self.case).addField()
        mdl = NumericalParamEquatModel(self.case)
        mdl.setMaximumIteration('Pressure',18)
        doc = '''<variable field_id="none" label="Pressure" name="pressure">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <max_iter_number>
                                 18
                         </max_iter_number>
                 </variable>'''
        assert mdl.getPressureNodefields() == self.xmlNodeFromString(doc),\
            'Could not set MaximumIteration'
        assert mdl.getMaximumIteration('Pressure') == 18,\
            'Could not get MaximumIteration'


def suite():
    testSuite = unittest.makeSuite(NumericalParamEquatTestCase, "check")
    return testSuite


def runTest():
    print("NumericalParamEquatTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
