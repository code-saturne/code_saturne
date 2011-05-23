# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
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
This module manages the differents possible outputs :
- listing printing
- post-processing and relationship with the FVM library
- monitoring points

This module defines the following classes:
- OutputControlModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string, sys, unittest
from types import FloatType

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class OutputControlModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor
        """
        self.case = case
        node_control  = self.case.xmlGetNode('analysis_control')
        self.node_out = node_control.xmlInitNode('output')


    def defaultInitialValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['listing_printing_frequency'] = 1
        default['postprocessing_frequency'] = -1
        default['postprocessing_frequency_time'] = 0.1
        default['probe_recording_frequency'] = 1
        default['probe_recording_frequency_time'] = 0.1
        default['postprocessing_options'] = "binary"
        default['postprocessing_format'] = "EnSight"
        if self.case['salome']:
            default['postprocessing_format'] = "MED"
        default['probe_format'] = "DAT"
        default['fluid_domain'] = "on"
        default['domain_boundary'] = "off"
        default['postprocessing_mesh'] = '0'
        default['coordinate'] = 0.0

        return default


    def __getCoordinates(self, name, coord):
        """
        Private method: return coordinate 'coord' for probe named 'name'
        """
        val = self.node_out.xmlGetNode('probe', name = name).xmlGetDouble(coord)
        if val == None:
            val = self.defaultInitialValues()['coordinate']
            self.__setCoordinates(name, coord, val)
        return val


    def __setCoordinates(self, name, coord, val):
        """
        Private method: put value of coordinate 'coord' for probe named 'name'
        """
        self.node_out.xmlGetNode('probe', name=name).xmlSetData(coord, val)


    def getListingFrequency(self):
        """
        Return the frequency for printing listing
        """
        f = self.node_out.xmlGetInt('listing_printing_frequency')
        if f == None:
            f = self.defaultInitialValues()['listing_printing_frequency']
            self.setListingFrequency(f)
        return f


    def setListingFrequency(self, freq):
        """
        Set the frequency for printing listing
        """
        self.isInt(freq)
        self.node_out.xmlSetData('listing_printing_frequency', freq)


    def getPostprocessingType(self):
        """
        Return the type of output for printing listing
        """
        node = self.node_out.xmlGetNode('postprocessing_frequency_time')
        if node != None :
            return 'Frequency_c_x'
        val = self.getPostprocessingFrequency()
        if val == -1 :
            return 'At the end'
        elif val == 1 :
            return 'At each step'
        else :
            return 'Frequency_c'


    def setPostprocessingType(self, type):
        """
        Set the type of output for printing listing
        """
        self.isInList(type, ['At the end', 'At each step', 'Frequency_c', 'Frequency_c_x'])

        if type == 'Frequency_c_x' :
            childNode = self.node_out.xmlGetNode('postprocessing_frequency')
            if childNode != None :
                childNode.xmlRemoveNode()
        else :
            childNode = self.node_out.xmlGetNode('postprocessing_frequency_time')
            if childNode != None :
                childNode.xmlRemoveNode()


    def getPostprocessingFrequency(self):
        """
        Return the frequency for post processing output
        """
        f = self.node_out.xmlGetInt('postprocessing_frequency')
        if f == None:
            f = self.defaultInitialValues()['postprocessing_frequency']
            self.setPostprocessingFrequency(f)
        return f


    def setPostprocessingFrequency(self, freq):
        """
        Set the frequency for post processing output
        """
        self.isInt(freq)
        self.node_out.xmlSetData('postprocessing_frequency', freq)


    def getPostprocessingFrequencyTime(self):
        """
        Return the frequency for post processing output
        """
        f = self.node_out.xmlGetDouble('postprocessing_frequency_time')
        if f == None:
            f = self.defaultInitialValues()['postprocessing_frequency_time']
            self.setPostprocessingFrequencyTime(f)
        return f


    def setPostprocessingFrequencyTime(self, freq):
        """
        Set the frequency for post processing output
        """
        self.isFloat(freq)
        self.node_out.xmlSetData('postprocessing_frequency_time', freq)


    def getFluidDomainPostProStatus(self):
        """
        Return status for traitment of fluid domain
        """
        nod = self.node_out.xmlInitNode('fluid_domain', 'status')
        status = nod['status']
        if not status:
            status = self.defaultInitialValues()['fluid_domain']
            self.setFluidDomainPostProStatus(status)
        return status


    def setFluidDomainPostProStatus(self, status):
        """
        Set status for traitment of fluid domain
        """
        self.isOnOff(status)
        node = self.node_out.xmlInitNode('fluid_domain', 'status')
        node['status'] = status


    def getDomainBoundaryPostProStatus(self):
        """
        Return status for traitment of domain boundary
        """
        nod = self.node_out.xmlInitNode('domain_boundary', 'status')
        status = nod['status']
        if not status:
            status = self.defaultInitialValues()['domain_boundary']
            self.setDomainBoundaryPostProStatus(status)
        return status


    def setDomainBoundaryPostProStatus(self, status):
        """
        Set status for traitment of domain boundary
        """
        self.isOnOff(status)
        node = self.node_out.xmlInitNode('domain_boundary', 'status')
        node['status'] = status


    def getTypePostMeshes(self):
        """
        Return choice of type of post processing for mesh
        """
        node = self.node_out.xmlInitNode('postprocessing_mesh_options', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.defaultInitialValues()['postprocessing_mesh']
            self.setTypePostMeshes(choice)
        return choice


    def setTypePostMeshes(self, choice):
        """
        Set choice of type of post processing for mesh
        """
        self.isInList(choice, ['0', '1', '2', '10', '11', '12'])
        node = self.node_out.xmlInitNode('postprocessing_mesh_options', 'choice')
        node['choice'] = choice


    def getPostProFormat(self):
        """
        Return choice of format for post processing output file
        """
        node = self.node_out.xmlInitNode('postprocessing_format', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.defaultInitialValues()['postprocessing_format']
            self.setPostProFormat(choice)
        return choice


    def setPostProFormat(self, choice):
        """
        Set choice of format for post processing output file
        """
        self.isInList(choice, ('EnSight', 'MED', 'CGNS'))
        node = self.node_out.xmlInitNode('postprocessing_format', 'choice')
        node['choice'] = choice


    def getPostProOptionsFormat(self):
        """
        Return options for post processing output file
        """
        node = self.node_out.xmlInitNode('postprocessing_options', 'choice')
        line = node['choice']
        if not line:
            line = self.defaultInitialValues()['postprocessing_options']
            self.setPostProOptionsFormat(line)
        return line


    def setPostProOptionsFormat(self, line):
        """
        Set options for post processing output file
        """
        list = string.split(line)
        self.isList(list)
        n = self.node_out.xmlInitNode('postprocessing_options', 'choice')
        n['choice'] = line


    def getMonitoringPointType(self):
        """
        Return the type of output for printing listing
        """
        node = self.node_out.xmlGetNode('probe_recording_frequency_time')
        if node != None :
            return 'Frequency_h_x'
        val = self.getMonitoringPointFrequency()
        if val == -1 :
            return 'At the end'
        elif val == 1 :
            return 'At each step'
        else :
            return 'Frequency_h'


    def setMonitoringPointType(self, type):
        """
        Set the type of output for printing listing
        """
        self.isInList(type, ['None', 'At each step', 'Frequency_h', 'Frequency_h_x'])

        if type == 'Frequency_h_x' :
            childNode = self.node_out.xmlGetNode('probe_recording_frequency')
            if childNode != None :
                childNode.xmlRemoveNode()
        else :
            childNode = self.node_out.xmlGetNode('probe_recording_frequency_time')
            if childNode != None :
                childNode.xmlRemoveNode()


    def getMonitoringPointFrequency(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetInt('probe_recording_frequency')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency']
            self.setMonitoringPointFrequency(f)
        return f


    def setMonitoringPointFrequency(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isInt(freq)
        self.node_out.xmlSetData('probe_recording_frequency', freq)


    def getMonitoringPointFrequencyTime(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetInt('probe_recording_frequency_time')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency_time']
            self.setMonitoringPointFrequencyTime(f)
        return f


    def setMonitoringPointFrequencyTime(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isFloat(freq)
        self.node_out.xmlSetData('probe_recording_frequency_time', freq)


    def getMonitoringPointFormat(self):
        """
        Return choice of format for post processing output file
        """
        node = self.node_out.xmlInitNode('probe_format', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.defaultInitialValues()['probe_format']
            self.setMonitoringPointFormat(choice)
        return choice


    def setMonitoringPointFormat(self, choice):
        """
        Set choice of format for probes
        """
        self.isInList(choice, ('DAT', 'CSV'))
        node = self.node_out.xmlInitNode('probe_format', 'choice')
        node['choice'] = choice


    def addMonitoringPoint(self, x=0.0, y=0.0, z=0.0):
        """
        Public method.
        Add a new monitoring point.
        @type x: C{Float}
        @param x: first coordinate
        @type y: C{Float}
        @param y: second coordinate
        @type z: C{Float}
        @param z: third coordinate
        """
        self.isFloat(x)
        self.isFloat(y)
        self.isFloat(z)
        num = str(self.getNumberOfMonitoringPoints() + 1)
        status="on"
        node = self.node_out.xmlInitNode('probe', name=num, status=status)
        for coord, val in [('probe_x', x), ('probe_y', y), ('probe_z', z)]:
            self.__setCoordinates(num, coord, val)


    def replaceMonitoringPointCoordinates(self, name, x=0.0, y=0.0, z=0.0):
        """
        Public method.
        Change the coordinates of a monitoring point.
        @type name: C{String}
        @param name: identifier of the monitoring point
        @type x: C{Float}
        @param x: first new coordinate
        @type y: C{Float}
        @param y: second new coordinate
        @type z: C{Float}
        @param z: third new coordinate
        """
        self.isFloat(x)
        self.isFloat(y)
        self.isFloat(z)
        self.isStr(name)
        self.isGreater(float(name), 0.0)
        self.isLowerOrEqual(float(name), self.getNumberOfMonitoringPoints())

        for coord, val in [('probe_x', x), ('probe_y', y), ('probe_z', z)]:
            self.__setCoordinates(name, coord, val)


    def deleteMonitoringPoints(self, list):
        """
        Public method.
        Conveniant method for the view. Delete a list of monitoring points.
        @type list: C{List} of C{Int}
        @param list: list of identifier of monitoring points to delete
        """
        list.sort()
        r = len(list)
        for n in range(r):
            name = str(list[n])
            self.deleteMonitoringPoint(name)
            for i in range(n, r):
                list[i] = list[i] - 1


    def deleteMonitoringPoint(self, num):
        """
        Public method.
        Delete a single monitoring point.
        @type num: C{String}
        @param num: identifier of the monitoring point
        """
        self.isStr(num)
        self.isGreater(float(num), 0.0)
        self.isLowerOrEqual(float(num), self.getNumberOfMonitoringPoints())

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

            # update the attribute "choice" of the probes markup for variables

            from Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
            listNodeVolum = OutputVolumicVariablesModel(self.case).listNodeVolum
            del OutputVolumicVariablesModel
            for nodeList in listNodeVolum:
                for node in nodeList:
                    n = node.xmlGetChildNode('probes')
                    if n:
                        nlist = n.xmlGetChildNodeList('probe_recording')
                        if not nlist:
                            n.xmlRemoveNode()
                        else:
                            n['choice']= str(len(nlist))


    def getMonitoringPointCoordinates(self, name):
        """
        Public method.
        @type name: C{String}
        @param name: identifier of the monitoring point
        @return: coordinates X, Y, Z for the monitoring point I{name}
        @rtype: C{List} of C{Float}
        """
        self.isStr(name)
        self.isGreater(float(name), 0.0)
        self.isLowerOrEqual(float(name), self.getNumberOfMonitoringPoints())
        X = self.__getCoordinates(name, 'probe_x')
        Y = self.__getCoordinates(name, 'probe_y')
        Z = self.__getCoordinates(name, 'probe_z')
        return X, Y, Z


    def getNumberOfMonitoringPoints(self):
        """
        Public method.
        @return: number of monitoring points already defined.
        @rtype: C{Int}
        """
        return len(self.node_out.xmlGetNodeList('probe'))


#-------------------------------------------------------------------------------
# OutputControlModel test Class
#-------------------------------------------------------------------------------

class OutputControlModelTestCase(ModelTest):
    """
    """
    def checkOutputControlModeInstantiation(self):
        """Check whether the OutputControlModel class could be instantiated"""
        model = None
        model = OutputControlModel(self.case)
        assert model != None, 'Could not instantiate OutputControlModel'

    def checkSetandGetListingFrequency(self):
        """Check whether the frequency of output listing could be set and get"""
        model = OutputControlModel(self.case)
        model.setListingFrequency(12)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <listing_printing_frequency>12</listing_printing_frequency>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
                    'Could not set frequency for listing output control model'
        assert model.getListingFrequency() == 12, \
                    'Could not get frequency for listing output control model'

    def checkSetandGetPostprocessingFrequency(self):
        """Check whether the frequency of post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setPostprocessingFrequency(13)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_frequency>13</postprocessing_frequency>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
                    'Could not set frequency for post processing output control model'
        assert model.getPostprocessingFrequency() == 13, \
                    'Could not get frequency for post processing output control model'

    def checkSetandGetFluidDomainPostProStatus(self):
        """Check whether the status of post processing for fluid domain could be set and get"""
        model = OutputControlModel(self.case)
        model.setFluidDomainPostProStatus('off')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <fluid_domain status="off"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
            'Could not set status of post processing for fluid domain for output control model'
        assert model.getFluidDomainPostProStatus() == 'off', \
            'Could not get status of post processing for fluid domain for output control model'

    def checkSetandGetDomainBoundaryPostProStatus(self):
        """
        Check whether the status of post processing for domain
        boundary could be set and get
        """
        model = OutputControlModel(self.case)
        model.setDomainBoundaryPostProStatus('on')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <domain_boundary status="on"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set status of post processing for domain boundary for output control model'
        assert model.getDomainBoundaryPostProStatus() == 'on', \
        'Could not get status of post processing for domain boundary for output control model'

    def checkSetandGetTypePostMeshes(self):
        """Check whether the type of mesh's post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setTypePostMeshes('2')
        doc = '''<output>
                    <postprocessing_mesh_options choice="2"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set type of post processing for mesh in output control model'
        assert model.getTypePostMeshes() == '2', \
        'Could not get type of post processing for mesh in output control model'

    def checkSetandGetPostProFormat(self):
        """Check whether the format for post processing could be set and get"""
        model = OutputControlModel(self.case)
        model.setPostProFormat('MED')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_format choice="MED"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set format of post processing in output control model'
        assert model.getPostProFormat() == 'MED', \
        'Could not get format of post processing in output control model'

    def checkSetandGetPostProOptionsFormat(self):
        """
        Check whether the options of format for post processing could be set and get
        """
        model = OutputControlModel(self.case)
        model.setPostProOptionsFormat('big_endian,divide_polyhedra')
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <postprocessing_options choice="big_endian,divide_polyhedra"/>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
        'Could not set the options of format for post processing in output control model'
        assert model.getPostProOptionsFormat() == 'big_endian,divide_polyhedra', \
        'Could not get the options of format for post processing in output control model'

    def checkSetandGetMonitoringPointFrequency(self):
        """
        Check whether the frequency of monitoring point could be set and get
        """
        model = OutputControlModel(self.case)
        model.setMonitoringPointFrequency(15)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe_recording_frequency>15</probe_recording_frequency>
                  </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not set the frequency of monitoring point in output control model'
        assert model.getMonitoringPointFrequency() == 15,\
        'Could not get the frequency of monitoring point in output control model'

    def checkAddMonitoringPoint(self):
        """
        Check whether monitoring point could be added
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not add monitoring point in output control model'
        assert model.getMonitoringPointCoordinates("1") == (11.1, 22.2, 33.3),\
        'Could not get monitoring point in output control model'

    def checkReplaceMonitoringPoint(self):
        """
        Check whether monitoring point could be replaced
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        model.addMonitoringPoint(5, 5.1, 5.21)
        model.replaceMonitoringPointCoordinates("2",5., 6, 7.)
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                    <probe name="2" status="on">
                        <probe_x>5</probe_x>
                        <probe_y>6</probe_y>
                        <probe_z>7</probe_z>
                    </probe>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not replace monitoring point in output control model'


    def checkDeleteMonitoringPoint(self):
        """
        Check whether monitoring point could be deleted
        """
        model = OutputControlModel(self.case)
        model.addMonitoringPoint(11.1, 22.2, 33.3)
        model.addMonitoringPoint(5, 5.1, 5.21)
        model.addMonitoringPoint(9.,8.,7.)
        model.deleteMonitoringPoint("2")
        doc = '''<output>
                    <postprocessing_mesh_options choice="0"/>
                    <probe name="1" status="on">
                        <probe_x>11.1</probe_x>
                        <probe_y>22.2</probe_y>
                        <probe_z>33.3</probe_z>
                    </probe>
                    <probe name="2" status="on">
                        <probe_x>9</probe_x>
                        <probe_y>8</probe_y>
                        <probe_z>7</probe_z>
                    </probe>
                </output>'''
        assert model.node_out== self.xmlNodeFromString(doc),\
        'Could not delete monitoring point in output control model'


def suite():
    testSuite = unittest.makeSuite(OutputControlModelTestCase, "check")
    return testSuite

def runTest():
    print("OutputControlModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
