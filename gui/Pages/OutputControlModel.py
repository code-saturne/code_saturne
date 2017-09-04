# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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
This module manages the differents possible outputs :
- listing printing
- post-processing and relationship with the FVM library
- monitoring points
- writer
- mesh

This module defines the following classes:
- OutputControlModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel import ModelTest

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
        default['probe_recording_frequency'] = 1
        default['probe_recording_frequency_time'] = 0.1
        default['probe_format'] = "DAT"
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


    @Variables.noUndo
    def getListingFrequency(self):
        """
        Return the frequency for printing listing
        """
        f = self.node_out.xmlGetInt('listing_printing_frequency')
        if f == None:
            f = self.defaultInitialValues()['listing_printing_frequency']
            self.setListingFrequency(f)
        return f


    @Variables.undoLocal
    def setListingFrequency(self, freq):
        """
        Set the frequency for printing listing
        """
        self.isInt(freq)
        self.node_out.xmlSetData('listing_printing_frequency', freq)


    @Variables.noUndo
    def getListingFrequencyLagrangian(self):
        """
        Return the frequency for printing listing for lagrangian variables
        """
        node_lagr = self.case.root().xmlInitNode('lagrangian', 'model')
        node_out_lag = node_lagr.xmlInitChildNode('output')
        f = node_out_lag.xmlGetInt('listing_printing_frequency')
        if f == None:
            f = self.defaultInitialValues()['listing_printing_frequency']
            self.setListingFrequencyLagrangian(f)
        return f


    @Variables.undoLocal
    def setListingFrequencyLagrangian(self, freq):
        """
        Set the frequency for printing listing for lagrangian variables
        """
        self.isInt(freq)
        node_lagr = self.case.root().xmlInitNode('lagrangian', 'model')
        node_out_lag = node_lagr.xmlInitChildNode('output')
        node_out_lag.xmlSetData('listing_printing_frequency', freq)


    def defaultWriterValues(self):
        """Return the default values - Method also used by ThermalScalarModel"""
        default = {}
        default['frequency_choice']          = 'none'
        default['frequency']                 = 1
        default['frequency_time']            = 1.
        default['output_at_start']           = 'off'
        default['output_at_end']             = 'on'
        default['format']                    = 'ensight'
        default['time_dependency']           = 'fixed_mesh'
        default['options']                   = ''
        default['directory']                 = 'postprocessing'

        if self.case['salome']:
            default['format'] = "med"

        return default


    def __defaultWriterLabelAndId(self):
        """
        Private method.
        Return a default id and label for a new writer.
        """
        id_table = []
        for l in self.getWriterIdList():
            id_table.append(int(l))
        user_table = []
        for l in id_table:
            if l > 0:
                user_table.append(l)
        if user_table != []:
            next_id = max(user_table) +1
        else:
            next_id = 1
        next_label = 'writer(' + str(next_id) + ')'
        n = next_id
        while next_label in self.getWriterLabelList():
            n = n+1
            next_label = 'writer(' + str(n) + ')'
        return str(next_id), next_label


    def __updateWriterId(self):
        #"""
        #Private method.
        #Update suffix number for writer label.
        #"""
        n = 0
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            if int(node['id']) > 0 :
                n = n + 1
                if node['label'] == 'writer(' + node['id'] + ')':
                    node['label'] = 'writer(' + str(n) + ')'
                node['id'] = str(n)


    def addWriter(self):
        """Public method.
        Input a new user writer
        """

        i, l = self.__defaultWriterLabelAndId()
        if l not in self.getWriterIdList():
            self.node_out.xmlInitNode('writer', id = i,label = l)
        self.getWriterFrequencyChoice(i)
        self.getWriterFormat(i)
        self.getWriterDirectory(i)
        self.getWriterOptions(i)
        self.getWriterTimeDependency(i)

        return i


    def addDefaultWriter(self):
        """Public method.
        Input a new user writer
        """
        list_writer = self.getWriterIdList()
        if list_writer == []:
            node = self.node_out.xmlInitNode('writer', id = "-1", label = 'results')
            node.xmlInitNode('frequency', period = 'none')
            node.xmlInitNode('output_at_end', status = 'on')
            node.xmlInitNode('format', name = 'ensight', options = '')
            node.xmlInitNode('directory', name = 'postprocessing')
            node.xmlInitNode('time_dependency', choice = 'fixed_mesh')


    def addDefaultLagrangianWriter(self):
        """Public method.
        Input new lagrangian writer
        """
        node = self.node_out.xmlGetNode('writer', 'label', id = "-3")
        if node == None:
            nodeL = self.node_out.xmlInitNode('writer', id = "-3", label = 'particles')
            nodeL.xmlInitNode('frequency', period = 'none')
            nodeL.xmlInitNode('output_at_end', status = 'on')
            nodeL.xmlInitNode('format', name = 'ensight', options = '')
            nodeL.xmlInitNode('directory', name = 'postprocessing')
            nodeL.xmlInitNode('time_dependency', choice = 'transient_connectivity')

        node = self.node_out.xmlGetNode('writer', 'label', id = "-4")
        if node == None:
            nodeT = self.node_out.xmlInitNode('writer', id = "-4", label = 'trajectories')
            nodeT.xmlInitNode('frequency', period = 'none')
            nodeT.xmlInitNode('output_at_end', status = 'on')
            nodeT.xmlInitNode('format', name = 'ensight', options = '')
            nodeT.xmlInitNode('directory', name = 'postprocessing')
            nodeT.xmlInitNode('time_dependency', choice = 'fixed_mesh')


    def __deleteWriter(self, writer_id):
        """
        Private method.
        Delete writer.
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlRemoveNode()
        self.__updateWriterId()


    def deleteWriter(self, writer_id):
        """
        Public method.
        Delete writer.
        """
        for w in writer_id:
            self.isInList(w, self.getWriterIdList())

            #delete associate with mesh if needed
            for MeshId in self.getMeshIdList():
                for WriterId in self.getAssociatedWriterIdList(MeshId):
                    if (WriterId == w):
                        self.deleteAssociatedWriter(MeshId, WriterId)

        # Delete all scalars
        for writer in reversed(writer_id):
            self.__deleteWriter(writer)


    @Variables.noUndo
    def getWriterIdList(self, lagrangian = None):
        """
        Return a list of writer id already defined
        """
        writer = []
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            if lagrangian != None:
                idd = int(node['id'])
                if lagrangian == 0:
                    if idd >= -1:
                        writer.append(node['id'])
                else:
                    if idd >= 0 or idd == -3 or idd == -4:
                        writer.append(node['id'])
            else:
                writer.append(node['id'])
        return writer


    @Variables.noUndo
    def getWriterLabelList(self):
        """
        Return a list of writer id already defined
        """
        writer = []
        for node in self.node_out.xmlGetNodeList('writer', 'label'):
            writer.append(node['label'])
        return writer


    @Variables.noUndo
    def getWriterIdFromLabel(self, label):
        """
        Return the label of a writer
        """
        node = self.node_out.xmlGetNodeList('writer', 'label')
        for n in node:
            if n['label'] == label:
                writer_id = n['id']
        return writer_id


    @Variables.noUndo
    def getWriterLabel(self, writer_id):
        """
        Return the label of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        n = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        label = n['label']
        if label == None:
            label = __defaultWriterLabelAndId(id=writer_id)
            self.setWriterLabel(writer_id, label)
        return label


    @Variables.undoLocal
    def setWriterLabel(self, writer_id, label):
        """
        Set the label of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        n = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n['label'] = label


    @Variables.noUndo
    def getWriterFrequencyChoice(self, writer_id):
        """
        Return the choice of frequency output for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('frequency')
        frequency_choice = n['period']
        if frequency_choice == None:
            frequency_choice = self.defaultWriterValues()['frequency_choice']
            self.setWriterFrequencyChoice(writer_id, frequency_choice)
        return frequency_choice


    @Variables.undoLocal
    def setWriterFrequencyChoice(self, writer_id, choice):
        """
        Set the choice of frequency output for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(choice, ('none', 'time_step', 'time_value', 'formula'))
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('frequency')
        old_choice = n['period']
        if old_choice != choice:
            n.xmlRemoveNode()
            n = node.xmlInitNode('frequency')
            default = self.defaultWriterValues()
            if choice == 'time_step':
                node.xmlSetData('frequency', str(default['frequency']))
            elif choice == 'time_value':
                node.xmlSetData('frequency', str(default['frequency_time']))
        n['period'] = choice


    @Variables.noUndo
    def getWriterFrequency(self, writer_id):
        """
        Return the frequency of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        freq = node.xmlGetString('frequency')
        return freq


    @Variables.undoLocal
    def setWriterFrequency(self, writer_id, freq):
        """
        Set the frequency of a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        node.xmlSetData('frequency', str(freq))


    @Variables.noUndo
    def getWriterOutputStartStatus(self, writer_id):
        """
        Return the output_at_start status of a mesh
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlGetNode('output_at_start')
        status = None
        if n:
            status = n['status']
        if status == None:
            status = self.defaultWriterValues()['output_at_start']
        return status


    @Variables.undoLocal
    def setWriterOutputStartStatus(self, writer_id, status):
        """
        Set the output_at_start status of a mesh
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isOnOff(status)
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        if status == self.defaultWriterValues()['output_at_start']:
            node.xmlRemoveChild('output_at_start')
        else:
            n = node.xmlInitNode('output_at_start')
            n['status'] = status


    @Variables.noUndo
    def getWriterOutputEndStatus(self, writer_id):
        """
        Return the output_at_end status of a mesh
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlGetNode('output_at_end')
        status = None
        if n:
            status = n['status']
        if status == None:
            status = self.defaultWriterValues()['output_at_end']
        return status


    @Variables.undoLocal
    def setWriterOutputEndStatus(self, writer_id, status):
        """
        Set the output_at_end status of a mesh
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isOnOff(status)
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        if status == self.defaultWriterValues()['output_at_end']:
            node.xmlRemoveChild('output_at_end')
        else:
            n = node.xmlInitNode('output_at_end')
            n['status'] = status


    @Variables.noUndo
    def getWriterFormat(self, writer_id):
        """
        Return the format for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('format')
        format = n['name']
        if format == None:
            format = self.defaultWriterValues()['format']
            self.setWriterFormat(writer_id, format)
        return format


    @Variables.undoLocal
    def setWriterFormat(self, writer_id, format):
        """
        Set the format for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(format, ('ensight', 'med', 'cgns', 'catalyst', 'ccm'))
        node = self.node_out.xmlInitNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('format')
        n['name'] = format


    @Variables.noUndo
    def getWriterDirectory(self, writer_id):
        """
        Return the directory for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('directory')
        directory = n['name']
        if directory == None:
            directory = self.defaultWriterValues()['directory']
            self.setWriterDirectory(writer_id, directory)
        return directory


    @Variables.undoLocal
    def setWriterDirectory(self, writer_id, directory):
        """
        Set the directory for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlInitNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('directory')
        n['name'] = directory


    @Variables.noUndo
    def getWriterOptions(self, writer_id):
        """
        Return the options for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('format')
        options = n['options']
        if options == None:
            options = self.defaultWriterValues()['options']
            self.setWriterOptions(writer_id, options)
        return options


    @Variables.undoLocal
    def setWriterOptions(self, writer_id, options):
        """
        Set the options for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('format')
        n['options'] = options


    @Variables.noUndo
    def getWriterTimeDependency(self, writer_id):
        """
        Return the type of time dependency for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)

        n = node.xmlInitNode('time_dependency')
        choice = n['choice']
        if choice == None:
            choice = self.defaultWriterValues()['time_dependency']
            self.setWriterTimeDependency(writer_id, choice)
        return choice


    @Variables.undoLocal
    def setWriterTimeDependency(self, writer_id, choice):
        """
        Set the type of time dependency for a writer
        """
        self.isInList(writer_id, self.getWriterIdList())
        self.isInList(choice, ('fixed_mesh', 'transient_coordinates', 'transient_connectivity'))
        node = self.node_out.xmlGetNode('writer', 'label', id = writer_id)
        n = node.xmlInitNode('time_dependency')
        n['choice'] = choice


    def defaultMeshValues(self):
        """Return the default values"""
        default = {}
        default['type'] = "cells"
        default['particles_type'] = "particles"
        default['all_variables'] = "on"
        default['location'] = "all[]"
        default['density'] = 1.

        return default


    def getMeshIdList(self):
        """
        Return a list of mesh id already defined
        """
        mesh = []
        for node in self.node_out.xmlGetNodeList('mesh'):
            mesh.append(node["id"])
        return mesh


    def __defaultMeshLabelAndId(self, basename='mesh'):
        """
        Private method.
        Return a default id and label for a new mesh.
        """
        id_table = []
        for l in self.getMeshIdList():
            id_table.append(int(l))
        user_table = []
        for l in id_table:
            if l > 0:
                user_table.append(l)
        if user_table != []:
            next_id = max(user_table) +1
        else:
            next_id = 1
        next_label = basename + '(' + str(next_id) + ')'
        n = next_id
        while next_label in self.getMeshLabelList():
            n = n+1
            next_label = basename + '(' + str(n) + ')'
        return str(next_id), next_label


    def __updateMeshId(self):
        #"""
        #Private method.
        #Update suffixe number for mesh label.
        #"""
        n = 0
        for node in self.node_out.xmlGetNodeList('mesh'):
            if int(node['id']) > 0 :
                n = n + 1
                if node['label'] == 'mesh('+ node['id'] + ')':
                    node['label'] ='mesh('+ str(n) + ')'
                node['id'] = str(n)


    def addMesh(self):
        """Public method.
        Input a new user mesh
        """

        i, l = self.__defaultMeshLabelAndId()
        if l not in self.getMeshIdList():
            self.node_out.xmlInitNode('mesh', id = i,label = l)
        self.getMeshAllVariablesStatus(i)
        self.getMeshType(i)
        self.getMeshLocation(i)

        return i


    def addLagrangianMesh(self):
        """Public method.
        Input a new user mesh
        """

        i, l = self.__defaultMeshLabelAndId(basename = 'particles')
        if l not in self.getMeshIdList():
            self.node_out.xmlInitNode('mesh', id = i,label = l)
        self.getMeshAllVariablesStatus(i)
        self.getLagrangianMeshType(i)
        self.getMeshDensity(i)
        self.getMeshLocation(i)

        return i


    def addDefaultMesh(self):
        """Public method.
        Input default mesh
        """
        list_mesh = self.getMeshIdList()
        if list_mesh == []:
            node1 = self.node_out.xmlInitNode('mesh', id = "-1",
                                              label = 'Fluid domain',
                                              type = 'cells')
            node1.xmlInitNode('all_variables', status = 'on')
            node1.xmlInitNode('location')
            node1.xmlSetData('location','all[]')
            node1.xmlInitNode('writer', id = '-1')

            node2 = self.node_out.xmlInitNode('mesh',
                                              id = "-2",
                                              label = 'Boundary',
                                              type = 'boundary_faces')
            node2.xmlInitNode('all_variables', status = 'on')
            node2.xmlInitNode('location')
            node2.xmlSetData('location', 'all[]')
            node2.xmlInitNode('writer', id = '-1')


    def addDefaultLagrangianMesh(self):
        """Public method.
        Input default particles mesh
        """
        list_mesh = self.getMeshIdList()
        node1 = self.node_out.xmlInitNode('mesh', id = "-3",
                                          label = 'particles',
                                          type = 'particles')
        node1.xmlInitNode('all_variables', status = 'on')
        node1.xmlInitNode('location')
        node1.xmlSetData('location','all[]')
        node1.xmlInitNode('density')
        node1.xmlSetData('density', 1)
        node1.xmlInitNode('writer', id = '-3')


    def __deleteMesh(self, mesh_id):
        """
        Private method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node.xmlRemoveNode()
        self.__updateMeshId()


    def deleteMesh(self, mesh_id):
        """
        Public method.
        Delete mesh.
        """
        for mesh in mesh_id:
            self.isInList(mesh, self.getMeshIdList())

        # Delete all scalars
        for mesh in reversed(mesh_id):
            self.__deleteMesh(mesh)


    @Variables.noUndo
    def getMeshLabelList(self):
        """
        Return a list of mesh id already defined
        """
        mesh = []
        for node in self.node_out.xmlGetNodeList('mesh'):
            mesh.append(node['label'])
        return mesh


    @Variables.noUndo
    def getMeshLabel(self, mesh_id):
        """
        Return the label of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        label = n['label']
        if label == None:
            label = __defaultMeshLabelAndId(id=mesh_id)
            self.setMeshLabel(mesh_id, label)
        return label


    @Variables.undoLocal
    def setMeshLabel(self, mesh_id, label):
        """
        Set the label of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n['label'] = label


    @Variables.noUndo
    def getMeshType(self, mesh_id):
        """
        Return the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        mesh_type = n['type']
        if mesh_type == None:
            mesh_type = self.defaultMeshValues()['type']
            self.setMeshType(mesh_id, mesh_type)
        return mesh_type


    @Variables.undoLocal
    def setMeshType(self, mesh_id, mesh_type):
        """
        Set the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(mesh_type, ('cells', 'interior_faces', 'boundary_faces'))
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node['type'] = mesh_type


    @Variables.noUndo
    def getLagrangianMeshType(self, mesh_id):
        """
        Return the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        mesh_type = n['type']
        if mesh_type == None:
            mesh_type = self.defaultMeshValues()['particles_type']
            self.setLagrangianMeshType(mesh_id, mesh_type)
        return mesh_type


    @Variables.undoLocal
    def setLagrangianMeshType(self, mesh_id, mesh_type):
        """
        Set the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(mesh_type, ('particles', 'trajectories'))
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node['type'] = mesh_type


    @Variables.noUndo
    def getMeshAllVariablesStatus(self, mesh_id):
        """
        Return the all_variables status of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n = node.xmlInitNode('all_variables')
        status = n['status']
        if status == None:
            status = self.defaultMeshValues()['all_variables']
            self.setMeshAllVariablesStatus(mesh_id, status)
        return status


    @Variables.undoLocal
    def setMeshAllVariablesStatus(self, mesh_id, status):
        """
        Set the all_variables status of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isOnOff(status)
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n = node.xmlInitNode('all_variables')
        n['status'] = status


    @Variables.noUndo
    def getMeshLocation(self, mesh_id):
        """
        Return the location of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        loc = node.xmlGetString('location')
        if loc == '':
            loc = self.defaultMeshValues()['location']
            self.setMeshLocation(mesh_id, loc)
        return loc


    @Variables.undoLocal
    def setMeshLocation(self, mesh_id, location):
        """
        Set the location of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        node.xmlSetData('location', location)


    @Variables.noUndo
    def getMeshDensity(self, mesh_id):
        """
        Return the density of a lagrangian mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        den = node.xmlGetString('density')
        if den == '':
            den = self.defaultMeshValues()['density']
            self.setMeshDensity(mesh_id, den)
        return den


    @Variables.undoLocal
    def setMeshDensity(self, mesh_id, density):
        """
        Set the density of a lagrangian mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        if float(density) < 1.0:
            node.xmlSetData('density', density)
        else:
            childNode = node.xmlGetNode('density')
            if childNode != None :
                childNode.xmlRemoveNode()


    @Variables.noUndo
    def getAssociatedWriterIdList(self, mesh_id):
        """
        Return a list of associated writer to a mesh already defined
        """
        self.isInList(mesh_id, self.getMeshIdList())
        associated_writer = []
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        for n in node.xmlGetNodeList('writer'):
            associated_writer.append(n["id"])
        return associated_writer


    def addAssociatedWriter(self, mesh_id, lagrangian):
        """Public method.
        Input a new user associated writer to a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        n = self.node_out.xmlGetNode('mesh', id = mesh_id)
        writer_list = self.getWriterIdList(lagrangian)
        associated_writer_list = []
        for writer in writer_list:
            if writer not in self.getAssociatedWriterIdList(mesh_id):
                associated_writer_list.append(writer)
        writer_id = None
        if associated_writer_list:
            n.xmlInitNode('writer', id = associated_writer_list[0])
            writer_id = associated_writer_list[0]
        return writer_id


    def __deleteAssociatedWriter(self, mesh_id, writer_id):
        """
        Private method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(writer_id, self.getAssociatedWriterIdList(mesh_id))
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        n = node.xmlGetNode('writer', id = writer_id)
        n.xmlRemoveNode()


    def deleteAssociatedWriter(self, mesh_id, writer_id):
        """
        Public method.
        Delete mesh.
        """
        self.isInList(mesh_id, self.getMeshIdList())
        self.isInList(writer_id, self.getAssociatedWriterIdList(mesh_id))

        self.__deleteAssociatedWriter(mesh_id, writer_id)

    @Variables.undoLocal
    def setAssociatedWriterChoice(self, mesh_id, writer_list):
        """
        Set the type of a mesh
        """
        self.isInList(mesh_id, self.getMeshIdList())
        for w in writer_list :
            self.isInList(w, self.getWriterIdList())
        node = self.node_out.xmlGetNode('mesh', id = mesh_id)
        for w in node.xmlGetNodeList('writer'):
            w.xmlRemoveNode()
        for w in writer_list:
            node.xmlInitNode('writer', id = w)


    @Variables.noUndo
    def getMonitoringPointType(self):
        """
        Return the type of output for printing listing
        """
        node = self.node_out.xmlGetNode('probe_recording_frequency_time')
        if node != None :
            return 'Frequency_h_x'
        val = self.getMonitoringPointFrequency()
        if val == -1 :
            return 'None'
        elif val == 1 :
            return 'At each step'
        else :
            return 'Frequency_h'


    @Variables.undoLocal
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


    @Variables.noUndo
    def getMonitoringPointFrequency(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetInt('probe_recording_frequency')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency']
            self.setMonitoringPointFrequency(f)
        return f


    @Variables.undoLocal
    def setMonitoringPointFrequency(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isInt(freq)
        self.node_out.xmlSetData('probe_recording_frequency', freq)


    @Variables.noUndo
    def getMonitoringPointFrequencyTime(self):
        """
        Return the frequency for recording probes
        """
        f = self.node_out.xmlGetDouble('probe_recording_frequency_time')
        if f == None:
            f = self.defaultInitialValues()['probe_recording_frequency_time']
            self.setMonitoringPointFrequencyTime(f)
        return f


    @Variables.undoLocal
    def setMonitoringPointFrequencyTime(self, freq):
        """
        Set the frequency for recording probes
        """
        self.isFloat(freq)
        self.node_out.xmlSetData('probe_recording_frequency_time', freq)


    @Variables.noUndo
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


    @Variables.undoLocal
    def setMonitoringPointFormat(self, choice):
        """
        Set choice of format for probes
        """
        self.isInList(choice, ('DAT', 'CSV'))
        node = self.node_out.xmlInitNode('probe_format', 'choice')
        node['choice'] = choice


    @Variables.undoLocal
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


    @Variables.undoLocal
    def ImportProbesFromCSV(self, fle):
        """
        Public method.
        Read a csv file to add monitoring probes
        """
        probFile = open(fle, "r")
        lines = probFile.readlines()

        num = 0

        for line in lines:
            tmp = line.split(',')
            if len(tmp) != 3:
                pass
            self.addMonitoringPoint(float(tmp[0]), float(tmp[1]), float(tmp[2]))
            num = num + 1
        return num


        propFile.close()

        self.isFloat(x)
        self.isFloat(y)
        self.isFloat(z)
        num = str(self.getNumberOfMonitoringPoints() + 1)
        status="on"
        node = self.node_out.xmlInitNode('probe', name=num, status=status)
        for coord, val in [('probe_x', x), ('probe_y', y), ('probe_z', z)]:
            self.__setCoordinates(num, coord, val)


    @Variables.undoLocal
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


    def deleteMonitoringPoints(self, lst):
        """
        Public method.
        Conveniant method for the view. Delete a list of monitoring points.
        @type list: C{List} of C{Int}
        @param list: list of identifier of monitoring points to delete
        """
        lst.sort()
        r = len(lst)
        for n in range(r):
            name = str(lst[n])
            self.deleteMonitoringPoint(name)
            for i in range(n, r):
                lst[i] = lst[i] - 1


    @Variables.undoLocal
    def deleteMonitoringPoint(self, num):
        """
        Public method.
        Delete a single monitoring point.
        @type num: C{String}
        @param num: identifier of the monitoring point
        """
        self.isStr(num)
        self.isGreater(int(num), 0)
        self.isLowerOrEqual(int(num), self.getNumberOfMonitoringPoints())

        # delete the node of the monitoring point

        node = self.node_out.xmlGetNode('probe', name=num)
        if node:
            node.xmlRemoveNode()
            self.case.xmlRemoveChild('probe_recording', name=num)

            # renumbering of all monitoring points

            for p in range(int(num)+1, self.getNumberOfMonitoringPoints()+2):
                probe = self.node_out.xmlGetNode('probe', name=p)
                probe['name'] = p - 1
                for probe_recording in self.case.xmlGetNodeList('probe_recording', name=p):
                    probe_recording['name'] = p - 1

    @Variables.noUndo
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


    @Variables.noUndo
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
                    <listing_printing_frequency>12</listing_printing_frequency>
                 </output>'''
        assert model.node_out== self.xmlNodeFromString(doc), \
                    'Could not set frequency for listing output control model'
        assert model.getListingFrequency() == 12, \
                    'Could not get frequency for listing output control model'

    def checkSetandGetMonitoringPointFrequency(self):
        """
        Check whether the frequency of monitoring point could be set and get
        """
        model = OutputControlModel(self.case)
        model.setMonitoringPointFrequency(15)
        doc = '''<output>
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
