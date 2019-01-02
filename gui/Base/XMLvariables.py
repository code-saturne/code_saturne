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

"""
This module defines the Variables class which creates the
<variable> and <property> markups.

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

from code_saturne.Base.Common import *
from code_saturne.Base import Toolbox

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
            msg = "There is an error: " + " ".join(liste) + " is not a list\n"
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


    def isIntEqual(self, ival1,  ival2):
        """This method verifies that val1 = val2"""
        if self.isInt(ival1) and self.isInt(ival2):
            if ival1 != ival2:
                msg = "There is an error: this value " + str(ival1) + "\n"\
                      "must be equal to " + str(ival2) + "\n"
                raise ValueError(msg)
        return True


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
# class Variables : creates <variable> and <property> markups.
#-------------------------------------------------------------------------------


class Variables:
    """
    This class creates <variable> and <property> markups.
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


    @staticmethod
    def undoGlobal(f):
        def _wrapper(self, *c, **d):
            """
            we suspend global record to prevent infinity loop
            use when call another class function
            """
            if self.case.record_global == True:
                self.case.undoGlobal(f, c)
                self.case.undoStop()
                r = f(self, *c, **d)
                self.case.undoStart()
            else:
                r = f(self, *c, **d)
            return r
        return _wrapper


    @staticmethod
    def undoLocal(f):
        def _wrapper2(self, *c, **d):
            self.case.undo(f, c)
            return f(self, *c, **d)
        return _wrapper2


    @staticmethod
    def noUndo(f):
        def _wrapper3(self, *c, **d):
            if self.case.record_global == True:
                self.case.undoStopGlobal()
                r = f(self, *c, **d)
                self.case.undoStartGlobal()
            else:
                r = f(self, *c, **d)
            return r
        return _wrapper3


    def updateLabel(self, vv):
        """
        """
        if not vv['label']:
            vv['label'] = Toolbox.dicoLabel(vv['name'])


    def setNewVariable(self, node, tag, dim=None, tpe=None, label=None):
        """
        Input a new <variable name="my_variable" label="ma_variable">
        in the xmldoc.
        """
        if not node.xmlGetNode('variable', name=tag):
            if dim != None:
                if tpe != None:
                    v1 = node.xmlInitNode('variable', name=tag, dimension=dim, type=tpe)
                else:
                    v1 = node.xmlInitNode('variable', name=tag, dimension=dim)
            else:
                if tpe != None:
                    v1 = node.xmlInitNode('variable', name=tag, type=tpe)
                else:
                    v1 = node.xmlInitNode('variable', name=tag)

            if label != None:
                v1['label'] = label
            self.updateLabel(v1)


    def setNewProperty(self, node, tag, dim=None):
        """
        Input a new <property name="my_property" label="ma_propriete">
        in the xmldoc.
        """
        if not node.xmlGetNode('property', name=tag):
            if dim != None:
                p1 = node.xmlInitNode('property', name=tag, dimension=dim)
            else:
                p1 = node.xmlInitNode('property', name=tag)
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
            p1.xmlInitChildNode('listing_printing')['status'] = "off"
            p1.xmlInitChildNode('postprocessing_recording')['status'] = "off"

            if label:
                p1['label'] = label
            else:
                self.updateLabel(p1)

        else:
            p1 = node.xmlGetNode('property', name=tag)

        return p1


    #Copie des methode de Variables pour Neptune (surcharge):

    def __setOutputControl__(self, variable, post=False):
        """
        Update the output markups <probe_recording name="XX">,
        <postprocessing_recording status='on'> and
        <listing_printing status='on'> for the new 'variable' markup.
        """
        variable.xmlInitChildNode('listing_printing', status='on')
        if post:
            variable.xmlInitChildNode('postprocessing_recording', status='on')
        else:
            variable.xmlInitChildNode('postprocessing_recording', status='off')


    def setNewVariableProperty(self, type, choice, node, num, name, label,
                               dim=None, support=None, post=False):
        id = str(num)
        label = label
        if not node.xmlGetNode(type, field_id=id,  name=name):
            if type != 'property' :
                if dim != None:
                    n = node.xmlInitNode(type, field_id=id, name=name, label=label,
                                         dimension=dim)
                else:
                    n = node.xmlInitNode(type, field_id=id, name=name, label=label,
                                         dimension='1')
            else :
                if support != None:
                    n = node.xmlInitNode(type, field_id=id, choice=choice,
                                         name=name, label=label,
                                         dimension='1', support = support)
                else:
                    n = node.xmlInitNode(type, field_id=id, choice=choice,
                                         name=name, label=label, dimension='1')
            self.__setOutputControl__(n, post=post)
            self.updateLabel(n)


    def removeVariableProperty(self, type, node, id, name):
        """
        """
        try :
            noder = node.xmlGetNode(type, field_id=id,  name=name)
            noder.xmlRemoveNode()
        except :
            pass


    def removeVariablesProperties(self, type, node, id):
        """
        """
        nodesList = node.xmlGetNodeList(type, field_id=str(id))
        for variableNode in nodesList:
            variableNode.xmlRemoveNode()


    def setNewTurbField(self, node, id, model, coupling):
        """
        Input a new node
        <field id="my_id" model="my_model", two_way_coupling="coupling">
        in the xmldoc.
        """
        if not node.xmlGetNode('field', field_id=id, model=model):
            v1 = node.xmlInitNode('field', field_id=str(id),
                                        model=model,
                                        two_way_coupling=coupling)
            self.updateLabel(v1)


    def getVariablesPropertiesList(self, average, constant) :
        """
        return list of variables, properties (and scalar)
        for Output field, profiles and averages
        if constant == yes we take account constant variables
        """
        self.XMLNodethermo  = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeclosure = self.case.xmlGetNode('closure_modeling')
        self.XMLNodeTurb    = self.XMLNodeclosure.xmlInitNode('turbulence')
        self.XMLNodeAna     = self.case.xmlGetNode('analysis_control')
        self.XMLNodeAverage = self.XMLNodeAna.xmlGetNode('time_averages')
        self.XMLUserScalar  = self.case.xmlGetNode('additional_scalars')
        self.XMLScalar      = self.XMLUserScalar.xmlInitNode('scalars')
        self.XMLUsers       = self.XMLUserScalar.xmlGetNode('users')
        list = []
        #TODO for variableType in ('variable', 'property', 'scalar') :
        for node in self.XMLNodethermo.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLNodeTurb.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLScalar.xmlGetNodeList('variable') :
            list.append(node)
        for node in self.XMLNodethermo.xmlGetNodeList('property') :
            choice = node['choice']
            if constant == 'yes' :
                list.append(node)
            elif choice != 'constant' :
                list.append(node)
        for node in self.XMLNodeTurb.xmlGetNodeList('property') :
            choice = node['choice']
            if constant == 'yes' :
                list.append(node)
            elif choice != 'constant' :
                list.append(node)
        if average == 'yes':
            # Warning average node is different
            for node in self.XMLNodeAverage.xmlGetNodeList('time_average'):
                list.append(node)
        if self.XMLUsers:
            for node in self.XMLUsers.xmlGetNodeList('variable'):
                list.append(node)

        return list

    #Copie des methode de Variables pour Neptune (surcharge) FIN


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
##        <variable name="my_variable" label="ma_variable" type="">
##        in the xmldoc could bet set.  """
##        root = self.case.root()
##        Variables(self.case).setNewThermalScalar(root, 'temperature_celsius', "0")
##        node = root.xmlGetChildNode('variable')
##        node['label'] = "toto"
##        doc = '<variable label="toto" name="temperature_celsius" type="thermal">'\
##                '<initial_value zone="0">20</initial_value>'\
##                '<min_value>-1e+12</min_value>'\
##                '<max_value>1e+12</max_value>'\
##              '</variable>'
##        assert node == self.xmlNodeFromString(doc),\
##            'Could not set the thermal variable markups'
##
##    def checkSetNewUserScalar(self):
##        """
##        Check whether a new <variable label="ma_variable" type="user">
##        in the xmldoc could bet set.
##        """
##        root = self.case.root()
##        Variables(self.case).setNewUserScalar(root, 'variable1')
##        node = root.xmlGetChildNode('variable')
##        doc = '<variable label="variable1" type="user"/>'
##        assert node == self.xmlNodeFromString(doc),\
##           'Could not set the user variable markups'
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
