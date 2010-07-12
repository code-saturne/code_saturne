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
This module defines the Variables class which creates the <variable>,
<scalar> and <property> markups.

This module contains the following classes and function:
- Variables
- VariablesTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import unittest


#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from Base.Common import *
from Base import Toolbox


#-------------------------------------------------------------------------------
# Class Model
#-------------------------------------------------------------------------------


class Model:
    """
    This class verifies int and float values mathematically
    """
    def isStr(self, string):
        """This method verifies that the string type is str"""
        if type(string) != str:
            msg = "There is an error : " + str(string) + " is not a string\n"
            raise ValueError(msg)
        return True


    def isList(self, liste):
        """This method verifies that list is not empty"""
        if type(liste) != list:
            msg = "There is an error: " + string.join(liste) + " is not a list\n"
            raise ValueError(msg)
        return True


    def isInt(self, ival):
        """This method verifies that ival is a int value"""
        if type(ival) != int:
            msg = "There is an error: this value " + str(ival) + " is not an integer\n"
            raise ValueError(msg)
        return True


    def isPositiveInt(self, ival):
        """This method verifies that ival is a int value > or = 0"""
        if self.isInt(ival):
            if ival < 0:
                msg = "There is an error: this value " + str(ival) + " must not be negative\n"
                raise ValueError(msg)
        return True


#    def isStrictPositiveInt(self, ival):
#        """This method verifies that ival is a int value > 0"""
#        if self.isInt(ival):
#            if ival <= 0:
#                msg = "There is an error: this value " + str(ival) + "\n"\
#                      "must not be neither negative neither 0\n"
#                raise ValueError, msg
#        return True


    def isIntEqual(self, ival1,  ival2):
        """This method verifies that val1 = val2"""
        if self.isInt(ival1) and self.isInt(ival2):
            if ival1 != ival2:
                msg = "There is an error: this value " + str(ival1) + "\n"\
                      "must be equal to " + str(ival2) + "\n"
                raise ValueError(msg)
        return True

##    def isStrictBetweenInt(self, ival,  imin, imax):
##        """This method verifies that ival is in imin and imax"""
##        if self.isInt(ival):
##            if ival <= imin or ival >= imax:
##                msg = "There is an error: this value " + str(ival) + "\n"\
##                      "must be strictly between " + str(imin) + "and" + str(imax) + "\n"
##                raise ValueError, msg
##        return True
##
##
##    def isBetweenInt(self, ival,  imin, imax):
##        """This method verifies that ival is in imin and imax"""
##        if self.isInt(ival):
##            if ival < imin or ival > imax:
##                msg = "There is an error: this value " + str(ival) + "\n"\
##                      "must be between " + str(imin) + "and" + str(imax) + "\n"
##                raise ValueError, msg
##        return True


    def isIntInList(self, ival, list):
        """This method verifies that ival is in list"""
        if self.isInt(ival) and self.isList(list):
            if ival not in list:
                msg = "There is an error: this value " + str(ival) + "\n"\
                      "is not in list " + list + "\n"
                raise ValueError(msg)
        return True


    def isFloat(self, val):
        """This method verifies that val is a float value > 0"""
        if type(val) != float and type(val) != int:
            msg = "There is an error: this value " + str(val) + " is not a float value\n"
            raise ValueError(msg)
        return True


    def isPositiveFloat(self, val):
        """This method verifies that val is a float value > or = 0"""
        if self.isFloat(val):
            if val < 0:
                msg = "There is an error: this value " + str(val) + " must not be negative\n"
                raise ValueError(msg)
        return True


    def isStrictPositiveFloat(self, val):
        """This method verifies that val is a float value > 0"""
        if self.isFloat(val):
            if val <= 0:
                msg = "There is an error: this value " + str(val) + "\n"\
                      "must not be neither negative neither 0\n"
                raise ValueError(msg)
        return True


    def isFloatEqual(self, val1, val2):
        """This method verifies that val1 = val2"""
        if self.isFloat(val1) and self.isFloat(val2):
            if val1 > (val2 + 1e-6) or val1 < (val2 - 1e-6):
                msg = "There is an error: this value " + str(val1) + "\n"\
                      "must be equal to " + str(val2) + "\n"
                raise ValueError(msg)
        return True


    def isGreater(self, val,  min):
        """This method verifies that val > min"""
        if self.isFloat(val):
            if val <= min:
                msg = "There is an error: this value " + str(val) + "\n"\
                      "must be greater than " + str(min) + "\n"
                raise ValueError(msg)
        return True


    def isGreaterOrEqual(self, val,  min):
        """This method verifies that val >= min"""
        if self.isFloat(val):
            if val < min:
                msg = "There is an error: this value " + str(val) + "\n"\
                      "must be greater or equal than " + str(min) + "\n"
                raise ValueError(msg)
        return True


    def isLower(self, val,  max):
        """This method verifies that val < max"""
        if self.isFloat(val):
            if val >= max:
                msg = "There is an error: this value " + str(val) + "\n"\
                      "must be lower than " + str(max) + "\n"
                raise ValueError(msg)
        return True


    def isLowerOrEqual(self, val,  max):
        """This method verifies that val <= max"""
        if self.isFloat(val):
            if val > max:
                msg = "There is an error: this value " + str(val) + "\n"\
                      "must be lower or equal than " + str(max) + "\n"
                raise ValueError(msg)
        return True


    def isFloatInList(self, val, list):
        """This method verifies that val is in list"""
        if self.isFloat(val):
            if val not in list:
                msg = "There is an error: this float value " + str(val) + "\n"\
                      "is not in list " + str(list) + "\n"
                raise ValueError(msg)
        return True


    def isInList(self, val, list):
        """This method verifies that val is in list"""
        if val not in list:
            msg = "There is an error: this value " + str(val) + "\n"\
                  "is not in list " + str(list) + "\n"
            raise ValueError(msg)
        return True


    def isNotInList(self, val, list):
        """This method verifies that val is not in list"""
        if val in list:
            msg = "There is an error: this value " + str(val) + "\n"\
                  "should be not in list " + str(list) + "\n"
            raise ValueError(msg)
        return True


    def isOnOff(self, status):
        """This method verifies that status is on or off"""
        if status not in ['on', 'off']:
            msg = "There is an error: this status " + str(status) + "\n"\
                  "is not in list " + str(list) + "\n"
            raise ValueError(msg)
        return True


#-------------------------------------------------------------------------------
# class Variables : creates <variable>, <scalar> and <property> markups.
#-------------------------------------------------------------------------------


class Variables:
    """
    This class creates <variable>, <scalar> and <property> markups.
    Each new markup has a 'name' and a 'label' attribute.
    Each new markup has <listing_printing status='on'>,
    <postprocessing_recording status='on'> and several
    <probe_recording name="XX"> (if any <probe> exists) as child markups.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case


##    def defaultInitialValues(self):
##        """
##        Return in a dictionnary which contains default values.
##        """
##        default = {}
##        #Initial values for thermal scalar: 20 deg C air at atmospheric pressure.
##        default['temperature_celsius'] = 20.0
##        default['temperature_kelvin'] = 293.15
##        default['enthalpy'] = 297413.
##
##        #Initial values for properties: 20 degC air at atmospheric pressure.
##        default['density'] = 1.17862
##        default['molecular_viscosity'] = 1.83e-05
##        default['diffusion_coefficient'] = 1.83e-05
##        default['specific_heat'] = 1017.24
##        default['thermal_conductivity'] = 0.02495
##
##        return default


##    def setOutputControl(self, variable):
##        """
##        Update the output markups <probe_recording name="XX">,
##        <postprocessing_recording status='on'> and
##        <listing_printing status='on'> for the new 'variable' markup.
##        """
##        analysis_ctrl = self.case.xmlGetNode('analysis_control')
##        node_output   = analysis_ctrl.xmlInitNode('output')
##        for node in node_output.xmlGetNodeList('probe', 'name'):
##            num = node['name']
##            variable.xmlInitChildNode('probe_recording', name=num)

##        variable.xmlInitChildNode('listing_printing', status='on')
##        variable.xmlInitChildNode('postprocessing_recording', status='on')


    def updateLabel(self, vv):
        """
        """
        if not vv['label']:
            vv['label'] = Toolbox.dicoLabel(vv['name'])


    def setNewTurbulenceVariable(self, node, tag):
        """
        Input a new <variable name="my_variable" label="ma_variable">
        in the xmldoc.
        """
        if not node.xmlGetNode('variable', name=tag):
            n = node.xmlInitNode('variable', name=tag)
##            n.xmlSetData('blending_factor', 0)
##            n.xmlInitNode('order_scheme', choice='upwind')

            self.updateLabel(n)


    def setNewVariable(self, node, tag):
        """
        Input a new <variable name="my_variable" label="ma_variable">
        in the xmldoc.
        """
        if not node.xmlGetNode('variable', name=tag):
            v1 = node.xmlInitNode('variable', name=tag)
            self.updateLabel(v1)


##    def setNewThermalScalar(self, node, tag, zone):
##        """
##        Input a child node
##        <scalar name="my_variable" label="ma_variable" type="thermal">
##        to the argument node. Initial values are for air at
##        atmospheric pressure.
##        """
##        if not node.xmlGetNode('scalar', type='thermal', name=tag):
##            s1 = node.xmlInitNode('scalar', type='thermal', name=tag)
##            s1.xmlInitChildNode('initial_value', zone=zone)
##            if tag == "temperature_celsius":
##                s1.xmlSetData('initial_value', self.defaultInitialValues()[tag])
##            elif tag == "temperature_kelvin":
##                s1.xmlSetData('initial_value', self.defaultInitialValues()[tag])
##            elif tag == "enthalpy":
##                s1.xmlSetData('initial_value', self.defaultInitialValues()[tag])
##            else:
##                print("Error in setNewThermalScalar:")
##                print("the given tag is %s" % (tag))
##                exit(0)
##
##            s1.xmlSetData('min_value', SMGRAND)
##            s1.xmlSetData('max_value', SPGRAND)
##
####            self.setOutputControl(s1)
##            self.updateLabel(s1)


    def setNewModelScalar(self, node, tag):
        """
        Input a new <scalar label="ma_variable" type="model">
        in the xmldoc.
        """
        if not node.xmlGetNodeList('scalar', type="model", name=tag, label=tag):
            s1 = node.xmlInitNode('scalar', type="model", name=tag, label=tag)
##            self.setOutputControl(s1)
            self.updateLabel(s1)
        else:
            s1 = node.xmlGetNode('scalar', type="model", name=tag, label=tag )

        return s1

##
##    def deleteAllModelScalars(self, node):
##        """
##        Input a new <scalar label="ma_variable" type="model">
##        in the xmldoc.
##        """
##        nodeList = node.xmlGetNodeList('scalar')
##        if nodeList != None:
##            for node in nodeList :
##                node.xmlRemoveNode()
##

##    def setWeightMatrixComponents(self, event=None):
##        """
##        Input a new <scalar label="ma_variable" type="user">
##        in the xmldoc.
##        """
##        node_np = self.case.xmlInitNode('numerical_parameters')
##        node_ipucou = node_np.xmlInitNode('velocity_pressure_coupling')
##
##        node_Tx = node_ipucou.xmlInitNode('variable', name='weight_matrix_X')
##        node_Ty = node_ipucou.xmlInitNode('variable', name='weight_matrix_Y')
##        node_Tz = node_ipucou.xmlInitNode('variable', name='weight_matrix_Z')
##
##        for (node, val) in [(node_Tx, 'weight_matrix_X'),
##                            (node_Ty, 'weight_matrix_Y'),
##                            (node_Tz, 'weight_matrix_Z')]:
##            self.setOutputControl(node)
##            if not node['label']: node['label'] = Toolbox.dicoLabel(val)

##    def deleteAllModelProperties(self, modelNode):
##        """
##        Input a new <scalar label="ma_variable" type="model">
##        in the xmldoc.
##        """
##        nodeList = modelNode.xmlGetNodeList('property')
##        if nodeList != None:
##            for node in nodeList :
##                node.xmlRemoveNode()


    def setNewProperty(self, node, tag):
        """
        Input a new <property name="my_property" label="ma_propriete">
        in the xmldoc.
        """
        if not node.xmlGetNode('property', name=tag):
            p1 = node.xmlInitNode('property', name=tag)
##            self.setOutputControl(p1)
            try :
                self.updateLabel(p1)
            except:
                pass
            if not p1['label']: p1['label'] = tag

        else:
            p1 = node.xmlGetNode('property', name=tag)

        return p1

    def setNewFluidProperty(self, node, tag, label=None):
        """
        Input a new <property name="my_property" label="ma_propriete", choice="constant">
        in the xmldoc.
        """
        if not node.xmlGetNode('property', name=tag):
            p1 = node.xmlInitNode('property', name=tag, choice='constant')
##            if tag == "density":
##                value = self.defaultInitialValues()[tag]
##            elif tag == "molecular_viscosity" or tag == "diffusion_coefficient":
##                value = self.defaultInitialValues()[tag]
##            elif tag == "specific_heat":
##                value = self.defaultInitialValues()[tag]
##            elif tag == "thermal_conductivity":
##                value = self.defaultInitialValues()[tag]
##            else:
##                print("Error in setNewFluidProperty:")
##                print("the given tag is %s" % (tag))
##                exit(0)

##            p1.xmlSetData('initial_value', value)
##            self.setOutputControl(p1)
            p1.xmlInitChildNode('listing_printing')['status'] = "off"
            p1.xmlInitChildNode('postprocessing_recording')['status'] = "off"

            if label:
                p1['label'] = label
            else:
                self.updateLabel(p1)

        else:
            p1 = node.xmlGetNode('property', name=tag)

        return p1


#-------------------------------------------------------------------------------
# MODEL test case
#-------------------------------------------------------------------------------


class ModelTestCase(unittest.TestCase):
    def checkIsInt(self):
        """ Check whether the Model class could be verify value is integer """
        ival = 3
        xval = 5.2
        assert Model().isInt(ival),'Should be an integer value'


    def checkIsPositiveInt(self):
        """Check whether the Model class could be verify value is a int value > or = 0"""
        ival = 3
        assert Model().isPositiveInt(ival) == True,\
        'Should be a positive integer value'


#    def checkIsStrictPositiveInt(self):
#        """Check whether the Model class could be verify value is a int value > or = 0"""
#        ival = 3
#        assert Model().isStrictPositiveInt(ival) == True,\
#        'Should be a strict positive integer value'


##    def checkIsIntbetween(self):
##        """Check whether the Model class could be verify value is between min-max"""
##        ival = 3
##        imin = 1
##        imax = 10
##        assert Model().isBetweenInt(ival, imin, imax) == True,\
##        'Could not verify value is between min and max integer values'
##        assert Model().isStrictBetweenInt(ival, imin, imax) == True,\
##        'Could not verify value is strictly between min and max integer values'

    def checkIsIntInList(self):
        """Check whether the Model class could be verify value is in list"""
        ival = 3
        list = ['toto', 3.4, 3, "machin"]
        assert Model().isIntInList(ival, list) == True,\
        'Could not verify value is in list '


    def checkIsFloat(self):
        """ Check whether the Model class could be verify value is float """
        ival = 3
        xval = 5.2
        assert Model().isFloat(xval) == True,'Should be a float value'


    def checkIsPositiveFloat(self):
        """Check whether the Model class could be verify value is a float value > or = 0"""
        val = 3.5
        assert Model().isPositiveFloat(val) == True,\
        'Should be a positive float value'


    def checkIsStrictPositiveFloat(self):
        """Check whether the Model class could be verify value is a float value > or = 0"""
        val = 3.5
        assert Model().isStrictPositiveFloat(val) == True,\
        'Should be a strict positive float value'


    def checkIsFloatInList(self):
        """Check whether the Model class could be verify value is in list"""
        val = 3.4
        list = ['toto', 3.4, 3, "machin"]
        assert Model().isFloatInList(val, list) ==True ,\
        'Could not verify value is in list '


    def checkIsList(self):
        """Check whether the Model class could be verify value is not empty"""
        list = ['toto', 3.4, 3, "machin"]
        assert Model().isList(list) == True, 'Should be a list'


    def checkIsInList(self):
        """Check whether the Model class could be verify value is in list"""
        list = ['toto', 3.4, 3, "machin"]
        assert Model().isInList('toto', list) == True, 'Should be in a list'


def suite1():
    testSuite = unittest.makeSuite(ModelTestCase, "check")
    return testSuite

def runTest():
    print("ModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite1())
##    runner.run(suite2())

#-------------------------------------------------------------------------------
# Variables test case
#-------------------------------------------------------------------------------


##class VariablesTestCase(unittest.TestCase):
##    """
##    """
##    def setUp(self):
##        """This method is executed before all "check" methods."""
##        import XMLengine
##        self.doc = XMLengine.XMLDocument("")
##        self.case = self.doc.parseString('<?xml version="2.0" ?>'\
##                                         '<Code_Saturne_GUI/>')
##        self.case.root().xmlAddChild('analysis_control')
##
##    def tearDown(self):
##        """This method is executed after all "check" methods."""
##        del self.case
##
##    def xmlNodeFromString(self, string):
##        """Private method to return a xml node from string"""
##        return self.doc.parseString(string).root()
##
##    def checkVariablesInstantiation(self):
##        """ Check whether the Variable class could be instantiated """
##        doc = None
##        doc = Variables(self.case)
##        assert doc != None, 'Could not instantiate Variables'
##
##    def checkSetOutputControl(self):
##        """ Check whether the <probe_recording name="XX">,
##        <postprocessing_recording status='on'> and
##        <listing_printing status='on'> markups could be set.  """
##        node = self.case.root().xmlAddChild('variable')
##        Variables(self.case).setOutputControl(node)
##        doc = '<variable/>'
##        assert node == self.xmlNodeFromString(doc),\
##            'Could not set the output control markups'
##
##    def checkSetNewVariable(self):
##        """ Check whether a new
##        <variable name="my_variable" label="ma_variable">
##        in the xmldoc could bet set.  """
##        root = self.case.root()
##        Variables(self.case).setNewVariable(root, 'pressure')
##        node = root.xmlGetChildNode('variable')
##        node['label'] = "toto"
##        doc = '<variable label="toto" name="pressure"/>'
##        assert  node == self.xmlNodeFromString(doc),\
##            'Could not set the variable markup'
##
##    def checkSetNewThermalScalar(self):
##        """ Check whether a new
##        <scalar name="my_variable" label="ma_variable" type="">
##        in the xmldoc could bet set.  """
##        root = self.case.root()
##        Variables(self.case).setNewThermalScalar(root, 'temperature_celsius', "0")
##        node = root.xmlGetChildNode('scalar')
##        node['label'] = "toto"
##        doc = '<scalar label="toto" name="temperature_celsius" type="thermal">'\
##                '<initial_value zone="0">20</initial_value>'\
##                '<min_value>-1e+12</min_value>'\
##                '<max_value>1e+12</max_value>'\
##              '</scalar>'
##        assert node == self.xmlNodeFromString(doc),\
##            'Could not set the thermal scalar markups'
##
##    def checkSetNewUserScalar(self):
##        """
##        Check whether a new <scalar label="ma_variable" type="user">
##        in the xmldoc could bet set.
##        """
##        root = self.case.root()
##        Variables(self.case).setNewUserScalar(root, 'scalar1')
##        node = root.xmlGetChildNode('scalar')
##        doc = '<scalar label="scalar1" type="user"/>'
##        assert node == self.xmlNodeFromString(doc),\
##           'Could not set the user scalar markups'
##
##    def checkSetNewProperty(self):
##        """
##        Check whether a new
##        <property name="my_property" label="ma_propriete">
##        in the xmldoc could bet set.
##        """
##        root = self.case.root()
##        Variables(self.case).setNewProperty(root, 'prop')
##        node = root.xmlGetChildNode('property')
##        truc = node.toString()
##        doc = '<property label="prop" name="prop"/>'
##        assert node == self.xmlNodeFromString(doc),\
##           'Could not set the property markups'
##
##    def checkSetNewFluidProperty(self):
##        """
##        Check whether a new
##        <property name="my_property" label="ma_propriete" choice="constant">
##        in the xmldoc could bet set.
##        """
##        root = self.case.root()
##        Variables(self.case).setNewFluidProperty(root, 'density')
##        node = root.xmlGetChildNode('property')
##        node['label'] = "toto"
##        doc = '<property choice="constant" label="toto" name="density">'\
##                '<listing_printing status="off"/>'\
##                '<postprocessing_recording status="off"/>'\
##                '</property>'
##        assert node == self.xmlNodeFromString(doc),\
##           'Could not set the fluid property markups'
##
##def suite():
##    testSuite = unittest.makeSuite(VariablesTestCase, "check")
##    return testSuite
##
##def runTest():
##    print("VariablesTestCase")
##    runner = unittest.TextTestRunner()
##    runner.run(suite())
##

#-------------------------------------------------------------------------------
# End of XMLvariables
#-------------------------------------------------------------------------------
