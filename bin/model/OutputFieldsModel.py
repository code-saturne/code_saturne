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
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.TimeAveragesModel import TimeAveragesModel

#-------------------------------------------------------------------------------
# Constructor
#-------------------------------------------------------------------------------

class OutputFieldsModel(MainFieldsModel, Variables, Model):
    """
    This class manages the Field objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case           = case
        self.XMLNodethermo  = self.case.xmlGetNode('thermophysical_models')
        self.analysis_ctrl  = self.case.xmlInitNode('analysis_control')
        self.node_output    = self.analysis_ctrl.xmlInitNode('output')
        self.node_probe     = self.node_output.xmlGetNodeList('probe','name')


    def defaultValues(self):
        default = {}
        default['listing'] = 'on'
        default['writer']  = 'on'
        return default


    def getVariableLabelsList(self) :
        """
        Return list of label for variable
        """
        lst = []
        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            for node in self.case.xmlGetNodeList(variableType) :
                lst.append(node['label'])
        return lst


    def getVariableProbeList(self):
        """
        Return list of node for probes
        Only for the View
        """
        probeList = []
        for node in self.node_probe:
            probeList.append(node['name'])
        return probeList


    @Variables.undoLocal
    def setVariableLabel(self, fieldId, label, oldlabel) :
        """
        set label
        """
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), label = oldlabel)
            if node != None:
                break

        if node != None:
            if variableType != 'time_average':
                node['label'] = label
            else:
                TimeAveragesModel(self.case).setLabel(node['label'], label)

        else :
            msg = "This variable " + oldlabel + " doesn't exist"
            raise ValueError(msg)

        # udpate label in formula if needed
        if oldlabel != '':
            for no in self.case.xmlGetNodeList('formula'):
                txt = no.xmlGetTextNode()
                if txt != None:
                    f = txt.replace(oldlabel, label)
                    no.xmlSetTextNode(f)


    @Variables.noUndo
    def getVariableLabel(self, fieldId, name) :
        """
        return label of name variable for fieldId
        """
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            label = node['label']
            return label
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    @Variables.undoLocal
    def setListingStatus(self, fieldId, label, status) :
        """
        return status for listing output for variable name on fieldId
        """
        self.isOnOff(status)
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), label = label)
            if node != None:
                break

        if node != None:
            n = node.xmlGetNode('listing_printing')
            n['status'] = status
        else :
            msg = "This variable " + label + " doesn't exist"
            raise ValueError(msg)


    @Variables.noUndo
    def getListingStatus(self, fieldId, name) :
        """
        return status for listing output for variable name on fieldId
        """
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            value = self.defaultValues()['listing']
            n = node.xmlGetNode('listing_printing')
            if n :
                value = n['status']
            return value
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    @Variables.undoLocal
    def setPostProcessingStatus(self, fieldId, label, status) :
        """
        return status for post processing for variable name on fieldId
        """
        self.isOnOff(status)
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), label = label)
            if node != None:
                break

        if node != None:
            n = node.xmlGetNode('postprocessing_recording')
            n['status'] = status
        else :
            msg = "This variable " + label + " doesn't exist"
            raise ValueError(msg)


    @Variables.noUndo
    def getPostProcessingStatus(self, fieldId, name) :
        """
        return status for post processing for variable name on fieldId
        """
        lst = self.getFieldIdList()
        lst.append("none")
        self.isInList(fieldId, lst)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break

        if node != None:
            value = self.defaultValues()['writer']
            n = node.xmlGetNode('postprocessing_recording')
            if n :
                value = n['status']
            return value
        else :
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


    @Variables.undoGlobal
    def setProbesList(self, fieldId, label, probes):
        """
        return list of probes for variable name on fieldId
        """
        l = self.getFieldIdList()
        l.append("none")
        self.isInList(fieldId, l)

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), label = label)
            if node != None:
                break
        if not node:
            msg = "This variable " + label + " doesn't exist"
            raise ValueError(msg)

        l1 = string.split(probes)
        l2 = self.getVariableProbeList()
        l1.sort()
        l2.sort()

        if l1 == l2:
            # if all probes are selected: no markup is nedded in the xml file
            node.xmlRemoveChild('probes')
            node.xmlRemoveChild('no_probe')
            return
        else:
            try:
                node.xmlRemoveChild('probes')
                node.xmlRemoveChild('no_probe')
            except:
                pass

            if l1:
                n = node.xmlInitNode('probes')
                l3 = []
                for probe in l1:
                    if probe not in l3 and probe in l2:
                        l3.append(probe)
                        n.xmlInitChildNodeList('probe_recording', name = probe)
            else:
                n = node.xmlInitNode('no_probe')


    @Variables.noUndo
    def getProbesList(self, fieldId, name):
        """
        return list of probes for variable name on fieldId
        """
        l = self.getFieldIdList()
        l.append("none")
        self.isInList(fieldId, l)

        l1 = self.getVariableProbeList()

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            node = self.case.xmlGetNode(variableType, field_id = str(fieldId), name = name)
            if node != None:
                break
        else:
            msg = "This variable " + name + " doesn't exist"
            raise ValueError(msg)


        node_probes = node.xmlGetChildNode('probes')

        if node.xmlGetChildNode('no_probe'):
            if node_probes:
                node_probes.xmlRemoveNode()
            return []

        if node_probes:
            l2 = node_probes.xmlGetChildNodeList('probe_recording')
            l3 = []

            if l2 == []:
                node_probes.xmlRemoveNode()
            else:
                for n in l2:
                    if n['name'] not in l3:
                        l3.append(n['name'])

            return l3
        else:
            return l1


    def getGlobalVariables(self) :
        """
        return list of variables with none field criteria
        """
        llst = {}
        lst = []

        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            for node in self.case.xmlGetNodeList(variableType, field_id = "none"):
                if not node['name'].startswith("User_"):
                    lst.append(node['name'])
                else:
                    idx = node['name'].split('_')[1]
                    tmp = idx.zfill(4)
                    llst[tmp] = node['name']
        keys = list(llst.keys())
        keys.sort()

        for key in keys:
            lst.append(llst[key])

        return lst


    def getFieldVariables(self, fieldId) :
        """
        return list of variables with none field criteria
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        lst = []
        for variableType in ('variable', 'property', 'scalar'):
            for node in self.case.xmlGetNodeList(variableType, field_id = str(fieldId)):
                lst.append(node['name'])

        return lst

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class OutputFieldsTestCase(ModelTest):
    """
    """
    def checkOutputFieldsInstantiation(self):
        """Check whether the OutputFieldsModel class could be instantiated"""
        model = None
        model = OutputFieldsModel(self.case)
        assert model != None, 'Could not instantiate OutputFieldsModel'


    def checkGetVariableLabelsList(self):
        """Check whether the  OutputFieldsModel class could get the VariableLabelsList"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        assert mdl.getVariableLabelsList() == ['Pressure', 'enthalpy1', 'alpha1', 'U1', 'V1', 'W1', 'Temp1', 'density1', 'Lam_vis1', 'Sp_heat1', 'Th_cond1', 'mass_trans1'],\
            'Could not get VariableLabelsList'


#    def checkGetVariableProbeList(self):
#        """Check whether the  NonCondensableModel class could get the VariableLabelsList"""
#        MainFieldsModel(self.case).addField()
#        mdl = OutputFieldsModel(self.case)
#        mdl.setProbesList('1','U1','1')
#        print(mdl.getProbesList('1','VelocityX'))
#        print(mdl.getVariableProbeList())
#        print(mdl.case)
#        assert mdl.getVariableLabelsList() == ,\
#            'Could not get VariableLabelsList'


    def checkGetandSetVariableLabel(self):
        """Check whether the OutputFieldsModel class could set and get VariableLabel"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        mdl.setVariableLabel('1','Vitesse1','U1')
        doc = '''<variable field_id="1" label="Vitesse1" name="velocityX">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                 </variable>'''
        assert mdl.case.xmlGetNode('variable', field_id = str('1'), label = 'Vitesse1') == self.xmlNodeFromString(doc),\
            'Could not set VariableLabel'
        assert mdl.getVariableLabel('1','velocityX') == 'Vitesse1',\
            'Could not get VariableLabel'


    def checkGetandSetListingStatus(self):
        """Check whether the OutputFieldsModel class could set and get ListingStatus"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        mdl.setListingStatus('1','U1','off')
        doc = '''<variable field_id="1" label="U1" name="velocityX">
                         <listing_printing status="off"/>
                         <postprocessing_recording status="on"/>
                 </variable>'''
        assert mdl.case.xmlGetNode('variable', field_id = str('1'), label = 'U1') == self.xmlNodeFromString(doc),\
            'Could not set ListingStatus'
        assert mdl.getListingStatus('1','velocityX') == 'off',\
            'Could not get ListingStatus'


    def checkGetandSetPostProcessingStatus(self):
        """Check whether the OutputFieldsModel class could set and get PostProcessingStatus"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        mdl.setPostProcessingStatus('1','U1','off')
        doc = '''<variable field_id="1" label="U1" name="velocityX">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="off"/>
                 </variable>'''
        assert mdl.case.xmlGetNode('variable', field_id = str('1'), label = 'U1') == self.xmlNodeFromString(doc),\
            'Could not set PostProcessingStatus'
        assert mdl.getPostProcessingStatus('1','velocityX') == 'off',\
            'Could not get PostProcessingStatus'


    def checkGetandSetProbesList(self):
        """Check whether the OutputFieldsModel class could set and get ProbesList"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        mdl.setProbesList('1','U1','1')
        doc = '''<variable field_id="1" label="U1" name="velocityX">
                         <listing_printing status="on"/>
                         <postprocessing_recording status="on"/>
                         <probes>
                                 <probe_recording name="1"/>
                         </probes>
                 </variable>'''
        assert mdl.case.xmlGetNode('variable', field_id = str('1'), label = 'U1') == self.xmlNodeFromString(doc),\
            'Could not set ProbesList'
        assert mdl.getProbesList('1','velocityX') == ['1'],\
            'Could not get ProbesList'


    def checkGetGlobalVariables(self):
        """Check whether the  OutputFieldsModel class could get GlobalVariables"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        assert mdl.getGlobalVariables() == ['Pressure'],\
            'Could not get GlobalVariables'


    def checkGetFieldVariables(self):
        """Check whether the  OutputFieldsModel class could get FieldVariables"""
        MainFieldsModel(self.case).addField()
        mdl = OutputFieldsModel(self.case)
        assert mdl.getFieldVariables('1') == ['Enthalpy', 'VolumeFraction', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'density', 'molecular_viscosity', 'specific_heat', 'thermal_conductivity', 'mass_trans'],\
            'Could not get FieldVariables'


def suite():
    testSuite = unittest.makeSuite(OutputFieldsTestCase, "check")
    return testSuite


def runTest():
    print("OutputFieldsTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
