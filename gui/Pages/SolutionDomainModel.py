# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
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
This module defines the XML calls for ecs execution
This module contains the following classes and function:
- MeshModel
- SolutionDomainModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, string, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.XMLvariables import Variables, Model
from Base.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Class Mesh Model
#-------------------------------------------------------------------------------

class MeshModel:
    """
    This class manages meshes's extension file and formats
    """
    def __init__(self):
        """
        Constructor.

        Initialize the dictionary file extension => format.
        """
        self.ext = {}
        self.ext['case']  = "ensight"
        self.ext['cgns']  = "cgns"
        self.ext['des']   = "des"
        self.ext['med']   = "med"
        self.ext['msh']   = "gmsh"
        self.ext['neu']   = "gambit"
        self.ext['ccm']   = "ccm"
        self.ext['ngeom'] = "ngeom"
        self.ext['unv']   = "ideas"
        self.ext['mesh']  = "meta"


    def getMeshExtension(self, mesh):
        """
        Public method.

        @return: Extension of the mesh file if it exists.
        @rtype: C{String}
        """
        # first check if the mesh is compressed with gzip
        if mesh.endswith(".gz"):
            mesh = mesh[:-3]

        extension = ""
        last_caracters = (string.split(mesh, ".")[-1:])[0]
        if last_caracters in self.ext.keys():
            extension = last_caracters
        return extension


    def getMeshFormat(self, mesh):
        """
        Public method.

        @return: Format of the mesh, if the extension is given.
        @rtype: C{String}
        """
        format = ""
        extension = self.getMeshExtension(mesh)
        if extension:
            format = self.ext[extension]
        return format


    def getExtensionFileList(self):
        """
        Public method.

        @return: List of all authorized extensions for mesh files.
        @rtype: C{List}
        """
        return self.ext.keys()


    def getBuildFormatList(self):
        """
        Public method.

        @return: List of number, format and description for view of popup.
        @rtype: C{List} of C{3-tuple}
        """
        list = [(0,  'ensight', 'EnSight (6 or Gold) ".case"' ),
                (1,  'cgns',    'CGNS ".cgns"'                ),
                (2,  'des',     'Simail (NOPO) ".des"'        ),
                (3,  'med',     'MED ".med"'                  ),
                (4,  'gmsh',    'Gmsh ".msh"'                 ),
                (5,  'gambit',  'GAMBIT Neutral ".neu"'       ),
                (6,  'ccm',     'STAR-CCM+ ".ccm"'            ),
                (7,  'ngeom',   'pro-STAR/STAR4 ".ngeom"'     ),
                (8,  'ideas',   'I-deas universal ".unv"'     ),
                (9,  'meta',    'Meta-mesh file ".mesh"'      )]

        return list


    def getMeshFormatDescription(self, format):
        """
        Public method.

        @type format: C{String}
        @param format: name of the I{format} (ensight, cgns,...)
        @return: Description of I{format}.
        @rtype: C{String}
        """
        label = ''
        for line in self.getBuildFormatList():
            if format == line[1]:
                label = line[2]
        return label


    def getFileFormatList(self):
        """
        Public method.

        @return: List of format and associated text for research file.
        @rtype: C{List} of C{2-tuple}
        """
        list = [("All files",                 "*"      ),
                ("EnSight (6 or Gold) files", "*.case" ),
                ("CGNS files",                "*.cgns" ),
                ("Simail (NOPO) files",       "*.des"  ),
                ("MED  files",                "*.med"  ),
                ("GMSH files",                "*.msh"  ),
                ("GAMBIT Neutral files",      "*.neu"  ),
                ("STAR-CCM+",                 "*.ccm"  ),
                ("pro-STAR/STAR4 files",      "*.ngeom"),
                ("I-deas universal files",    "*.unv"  ),
                ("Meta-mesh files",           "*.mesh" ) ]

        return list


#-------------------------------------------------------------------------------
# Class SolutionDomainModel
#-------------------------------------------------------------------------------

class SolutionDomainModel(MeshModel, Model):
    """
    This class allow to call function for fill saturne running file (lance)
    """

    def __init__(self, case):
        """
        Simple constructor.
        """
        self.case = case

        self.node_ecs        = self.case.xmlGetNode('solution_domain')
        self.node_meshes     = self.node_ecs.xmlInitNode('meshes_list')
        self.node_cut        = self.node_ecs.xmlInitNode('faces_cutting', "status")
        self.node_orient     = self.node_ecs.xmlInitNode('reorientation', "status")
        self.node_join       = self.node_ecs.xmlInitNode('joining')
        self.node_perio      = self.node_ecs.xmlInitNode('periodicity')
        self.node_standalone = self.node_ecs.xmlInitNode('standalone')


#************************* Private methods *****************************

    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['cutting_status'] = "off"
        defvalue['reorientation']  = "off"
        defvalue['select_status']  = "off"
        defvalue['selector']       = "all[]"
        defvalue['fraction']       = 0.1
        defvalue['plane']          = 25.0
        defvalue['verbosity']      = 1
        defvalue['angle']          = 0.01
        defvalue['syrth_status']   = "off"
        defvalue['syrth_mesh_2d']  = "off"
        defvalue['sim_status']     = "on"
        defvalue['verif_mail']     = "on"
        defvalue['postprocessing_format'] = "EnSight"
        defvalue['postprocessing_options'] = "binary"
        defvalue['dir_cas']        = "default_case"
        defvalue['poly_status']    = "off"
        defvalue['perio_mode']     = "translation"
        defvalue['transfo_val']    = 0.0

        return defvalue


    def _updateBatchScriptFile(self, keyword):
        """
        Update, for keyword, the backup file if it's ready to run.
        """
        self.isInList(keyword,('MESHES',
                               'REORIENT'))
        key = self.case['computer']
        if key:
            if not self.case['batchScript'][key]: return

            from Pages.BatchRunningModel import BatchRunningModel
            batch = BatchRunningModel(self.case)
            batch.initializeBatchScriptFile()
            batch.updateBatchScriptFile(keyword)
            del BatchRunningModel


#To follow : private methods to get or put faces
#================================================

    def _getTagNode(self, tagName):
        """
        Private method: Return node corresponding at item "tag"
        """
        self.isInList(tagName, ('face_joining', 'face_periodicity'))
        if tagName == 'face_joining':
            node = self.node_join
        elif tagName == 'face_periodicity':
            node = self.node_perio
        return node


    def _getJoinNode(self, join_id):
        """
        Get node for a given joining
        """
        node = None
        listNode = self.node_join.xmlGetNodeList('face_joining')
        if join_id < len(listNode):
            node = listNode[join_id]

        return node

    def _updateJoinSelectionNumbers(self):
        """
        Update names of join selection
        """
        listNode = self.node_join.xmlGetNodeList('face_joining')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    def _addJoinSelect(self, node, select):
        """
        Private method: Add faces to node (join, periodic) with dictionary select.
        """
        for sel, txt in [ (select['selector'],  'selector'),
                          (select['fraction'],  'fraction'),
                          (select['plane'],     'plane'),
                          (select['verbosity'], 'verbosity')]:
            if sel:
                node.xmlSetData(txt, sel)
            else:
                node.xmlRemoveChild(txt)


    def _getFaces(self, node):
        """
        Private method: Return values found for joining for a given node
        """
        default = {}
        default['selector'] =""
        default['fraction'] = ""
        default['plane'] = ""
        default['verbosity'] = ""

        if node:
            default['selector']  = node.xmlGetString('selector')
            default['fraction']  = node.xmlGetString('fraction')
            default['plane']     = node.xmlGetString('plane')
            default['verbosity'] = node.xmlGetString('verbosity')

        if     default['selector'] == '' and default['fraction'] == "" \
           and default['plane'] == "" and default['verbosity'] == "":
            default = {}

        return default


    def _removeJoinChildren(self, node):
        """
        Private method: Remove all child nodes of node for one selection
        """
        for tag in ('selector',
                    'fraction',
                    'plane',
                    'verbosity',):
            node.xmlRemoveChild(tag)


#To follow : private methods for periodicity:
#===========================================

    def _getPerioNode(self, perio_id):
        """
        Get node for a given periodicity
        """
        node = None
        listNode = self.node_perio.xmlGetNodeList('face_periodicity')
        if perio_id < len(listNode):
            node = listNode[perio_id]

        return node

    def _updatePerioSelectionNumbers(self):
        """
        Update names of periodicity selections
        """
        listNode = self.node_perio.xmlGetNodeList('face_periodicity')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    def _setPeriodicNewMode(self, perio_id, new_mode):
        """
        Private method: Set node and mode of given periodicity'
        """
        listNode = self.node_perio.xmlGetNodeList('face_periodicity')
        if perio_id < len(listNode):
            node = listNode[perio_id]
            node['mode'] = new_mode


    def _setTranslationDefault(self, perio_id):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self._getPerioNode(perio_id)

        if node:
            nList = node.xmlInitChildNodeList('translation')
            for n in nList:
                n.xmlSetData('translation_x', self.defaultValues()['transfo_val'])
                n.xmlSetData('translation_y', self.defaultValues()['transfo_val'])
                n.xmlSetData('translation_z', self.defaultValues()['transfo_val'])


    def _setRotationDefault(self, perio_id):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self._getPerioNode(perio_id)

        if node:
            nList = node.xmlInitChildNodeList('rotation')
            for n in nList:
                n.xmlSetData('angle', self.defaultValues()['transfo_val'])
                n.xmlSetData('axis_x', self.defaultValues()['transfo_val'])
                n.xmlSetData('axis_y', self.defaultValues()['transfo_val'])
                n.xmlSetData('axis_z', self.defaultValues()['transfo_val'])
                n.xmlSetData('invariant_x', self.defaultValues()['transfo_val'])
                n.xmlSetData('invariant_y', self.defaultValues()['transfo_val'])
                n.xmlSetData('invariant_z', self.defaultValues()['transfo_val'])


    def _setMixedDefault(self, perio_id):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self._getPerioNode(perio_id)

        if node:
            nList = node.xmlInitChildNodeList('mixed')
            for n in nList:
                for txt in ('matrix_12', 'matrix_13', 'matrix_14',
                            'matrix_21', 'matrix_23', 'matrix_24',
                            'matrix_31', 'matrix_32', 'matrix_34'):
                    n.xmlSetData(txt, 0.0)
                for txt in ('matrix_11', 'matrix_22', 'matrix_33'):
                    n.xmlSetData(txt, 1.0)


#************************* Methods callable by users*****************************

# Methods to manage meshes :
#=========================

    def addMesh(self, mesh, format=None):
        """
        Public method. Add mesh name in xml file.
        """
        self.isNotInList(mesh, self.getMeshList())
        if not format:
            format = MeshModel().getMeshFormat(mesh)
        else:
            self.isInList(format, MeshModel().ext.values())

        self.node_meshes.xmlInitNode('mesh', name=mesh, format=format)
        self._updateBatchScriptFile('MESHES')


    def delMesh(self, mesh):
        """
        Public method. Delete node for mesh named "mesh" in xml file
        """
        self.isInList(mesh, self.getMeshList())
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            if node['name'] == mesh:
                node.xmlRemoveNode()
        self._updateBatchScriptFile('MESHES')


    def getMeshList(self):
        """
        Public method. Return the meshes name list already put in the case.
        """
        meshList = []
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            meshList.append(node['name'])
        return meshList


    def getMeshFormat(self, mesh):
        """
        Public method. Return the mesh name format recorded in the case.
        """
        self.isInList(mesh, self.getMeshList())
        return self.node_meshes.xmlGetNode('mesh', name=mesh)['format']


    def getMeshExtendedFormat(self, mesh):
        """
        Public method. Return the mesh extended format.
        """
        fmt = MeshModel().getMeshFormat(mesh)
        for line in MeshModel().getBuildFormatList():
            if fmt == line[1]:
                label = line[2]
        return label


    def setMeshNumber(self, mesh, num):
        """
        Public method. Put the mesh number.
        """
        self.isInList(mesh, self.getMeshList())
        self.isInt(num)
        self.isGreater(num, 0)

        self.node_meshes.xmlGetNode('mesh', name=mesh)['num'] = num
        self._updateBatchScriptFile('MESHES')


    def getMeshNumber(self, mesh):
        """
        Public method. Return the mesh number recorded in the case.
        """
        self.isInList(mesh, self.getMeshList())
        num = self.node_meshes.xmlGetNode('mesh', name=mesh)['num']
        if num:
            num = int(num)
        return num


    def setMeshGroupCells(self, mesh, grp_cel):
        """
        Public method. Put the grp-cel option.
        """
        self.isInList(mesh, self.getMeshList())
        self.isInList(grp_cel, ('off', 'section', 'zone'))

        if grp_cel == "off":
            del self.node_meshes.xmlGetNode('mesh', name=mesh)['grp_cel']
        else:
            self.node_meshes.xmlGetNode('mesh', name=mesh)['grp_cel'] = grp_cel
        self._updateBatchScriptFile('MESHES')


    def getMeshGroupCells(self, mesh):
        """
        Public method. Return the mesh 'grp-cel' sub-option recorded in the case.
        """
        return self.__getMeshGroup(mesh, 'grp_cel')


    def setMeshGroupFaces(self, mesh, grp_fac):
        """
        Public method. Put the 'grp-fac' sub-option.
        """
        self.isInList(mesh, self.getMeshList())
        self.isInList(grp_fac, ('off', 'section', 'zone'))

        if grp_fac == "off":
            del self.node_meshes.xmlGetNode('mesh', name=mesh)['grp_fac']
        else:
            self.node_meshes.xmlGetNode('mesh', name=mesh)['grp_fac'] = grp_fac
        self._updateBatchScriptFile('MESHES')


    def getMeshGroupFaces(self, mesh):
        """
        Public method. Return the mesh 'grp_fac' option recorded in the case.
        """
        return self.__getMeshGroup(mesh, 'grp_fac')


    def __getMeshGroup(self, mesh, group):
        """
        Private method. Return the mesh 'grp_fac' or 'grp_cel' sub-option recorded in the case.
        """
        self.isInList(mesh, self.getMeshList())
        grp = self.node_meshes.xmlGetNode('mesh', name=mesh)[group]
        if grp == None:
            grp = 'off'
        return grp

# Methods to manage status of all main tags :
#==========================================

    def getCutStatus(self):
        """
        Get status on tag "faces_cutting" from xml file
        """
        status = self.node_cut['status']
        if not status:
            status = self.defaultValues()['cutting_status']
            self.setCutStatus(status)
        return status


    def setCutStatus(self, status):
        """
        Put status on tag "faces_cutting" in xml file
        """
        self.isOnOff(status)
        self.node_cut['status'] = status


    def setCutAngle(self, var):
        """
        input '--cut_warped_faces' parameter.
        """
        self.isGreaterOrEqual(var, 0.0)
        if var != self.defaultValues()['angle']:
            self.node_cut.xmlSetData('warp_angle_max', var)
        else:
            self.node_cut.xmlRemoveChild('warp_angle_max')


    def getCutAngle(self):
        """
        get '--cut_warped_faces' parameters.
        """
        angle = self.node_cut.xmlGetDouble('warp_angle_max')
        if angle == None:
            angle = self.defaultValues()['angle']
        return angle


    def getOrientation(self):
        """
        Get status on tag "reorientation" from xml file
        """
        status = self.node_orient['status']
        if not status:
            status = self.defaultValues()['reorientation']
            self.setOrientation(status)
        return status


    def setOrientation(self, status):
        """
        Put status on tag "reorientation" in xml file
        """
        self.isOnOff(status)
        self.node_orient['status'] = status
        self._updateBatchScriptFile('REORIENT')


    def getSimCommStatus(self):
        """
        Get status of tag ''similation_communication' into xml file
        """
        node = self.node_standalone.xmlInitNode('simulation_communication', 'status')
        status = node['status']
        if not status:
            status = self.defaultValues()['sim_status']
            self.setSimCommStatus(status)
        return status


    def setSimCommStatus(self, status):
        """
        Put status of tag ''similation_communication' into xml file
        """
        self.isOnOff(status)
        node = self.node_standalone.xmlInitNode('simulation_communication', 'status')
        node['status'] = status


    def getPostProFormat(self):
        """
        Return choice of format for post processing output file
        """
        node = self.node_standalone.xmlInitNode('postprocessing_format', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.defaultValues()['postprocessing_format']
            self.setPostProFormat(choice)
        return choice


    def setPostProFormat(self, choice):
        """
        Set choice of format for post processing output file
        """
        self.isInList(choice, ('EnSight', 'MED_fichier', 'CGNS'))
        node = self.node_standalone.xmlInitNode('postprocessing_format', 'choice')
        node['choice'] = choice


    def getPostProOptionsFormat(self):
        """
        Return options for post processing output file
        """
        node = self.node_standalone.xmlInitNode('postprocessing_options', 'choice')
        line = node['choice']
        if not line:
            line = self.defaultValues()['postprocessing_options']
            self.setPostProOptionsFormat(line)
        return line


    def setPostProOptionsFormat(self, line):
        """
        Set options for post processing output file
        """
        list = string.split(line)
        self.isList(list)
        n = self.node_standalone.xmlInitNode('postprocessing_options', 'choice')
        n['choice'] = line


# Methods to manage periodicity :
#==============================

    def getPeriodicSelectionsCount(self):
        """
        Public method.

        @return: number of periodic faces selections
        @rtype: C{int}
        """
        return len(self.node_perio.xmlGetNodeList('face_periodicity'))


    def getPeriodicityMode(self, perio_id):
        """
        Public method.

        @type perio_id: C{int}
        @param perio_id: id of the periodic boundary
        @return: mode of transformation of periodic boundary I{perio_id}
        @rtype: C{str}
        """
        node = self._getPerioNode(perio_id)

        mode = node['mode']
        if not mode:
            mode = self.defaultValues()['perio_mode']

        return mode


    def updatePeriodicityMode(self, perio_id, mode):
        """
        Public method.

        Update transformation mode from a periodic boundary

        @type perio_id: C{int}
        @param perio_id: id of the periodic boundary
        @type mode: C{str}
        @param mode: mode of the periodic boundary (i.e.: 'translation', 'rotation', 'mixed')
        """
        node = self._getPerioNode(perio_id)

        self.isInList(mode, ('translation', 'rotation', 'mixed'))

        if node['mode'] != mode:
            node['mode'] = mode

            if mode in ('translation', 'rotation', 'mixed'):
                if not node.xmlGetChildNode(mode):
                  if mode =="translation":
                      self._setTranslationDefault(perio_id)
                  elif mode =="rotation":
                      self._setRotationDefault(perio_id)
                  elif mode =="mixed":
                      self._setMixedDefault(perio_id)


    def deletePeriodicity(self, perio_id):
        """
        Public method.

        Delete a transformation in periodic boundary.

        @type perio_id: C{str}
        @param perio_id: name of the periodic boundary
        """
        node = self._getPerioNode(perio_id)
        node.xmlRemoveNode()
        if perio_id < self.getPeriodicSelectionsCount():
            self._updatePerioSelectionNumbers()


    def getTranslationDirection(self, perio_id):
        """
        Public method.

        Get values of translation for periodic translation

        @type perio_id: C{int}
        @param perio_id: name of the periodic boundary
        @return: values of translation for periodic translation
        @rtype: 3 C{float}
        """
        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('translation')
        dx = n.xmlGetString('translation_x')
        dy = n.xmlGetString('translation_y')
        dz = n.xmlGetString('translation_z')

        return dx, dy, dz


    def setTranslationDirection(self, perio_id, dir, valcoor):
        """
        Put values of translation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(dir, ('translation_x', 'translation_y', 'translation_z'))

        node = self._getPerioNode(perio_id)

        for n in node.xmlGetChildNodeList('translation'):
            n.xmlSetData(dir, valcoor)


    def getRotationDirection(self, perio_id):
        """
        Get values for director vector rotation for periodic translation
        """
        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        rx = n.xmlGetString('axis_x')
        ry = n.xmlGetString('axis_y')
        rz = n.xmlGetString('axis_z')

        return rx, ry, rz


    def setRotationVector(self, perio_id, dir, valcoor):
        """
        Put values for director vector rotation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(dir, ("axis_x", "axis_y", "axis_z"))

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData(dir,valcoor)


    def getRotationAngle(self, perio_id):
        """
        Get angle for rotation for periodic rotation
        """
        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        angle = n.xmlGetString('angle')

        return angle


    def setRotationAngle(self, perio_id, angle):
        """
        Put angle for rotation for periodic rotation
        """
        self.isGreaterOrEqual(angle, 0.0)

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData('angle', angle)


    def getRotationCenter(self, perio_id):
        """
        Get coordinates of center of rotation for periodic transformation
        """
        mode = self.getPeriodicityMode(perio_id)
        self.isInList(mode, ('rotation'))

        node = self._getPerioNode(perio_id)

        if mode == "rotation":
            n = node.xmlGetChildNode('rotation')
        px = n.xmlGetString('invariant_x')
        py = n.xmlGetString('invariant_y')
        pz = n.xmlGetString('invariant_z')

        return px, py, pz


    def setRotationCenter(self, perio_id, pos, val):
        """
        Put coordinates of center of rotation for periodic transformation
        """
        self.isFloat(val)
        self.isInList(pos, ('invariant_x', 'invariant_y', 'invariant_z'))
        mode = self.getPeriodicityMode(perio_id)
        self.isInList(mode, ('rotation'))

        node = self._getPerioNode(perio_id)

        if mode == 'rotation':
            n = node.xmlGetChildNode('rotation')
        n.xmlSetData(pos, val)


    def getTransformationMatrix(self, perio_id):
        """
        Get values of matrix of rotation for periodic transformation
        """
        mode = self.getPeriodicityMode(perio_id)
        self.isInList(mode, ('mixed'))

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('mixed')
        m11 = n.xmlGetString('matrix_11')
        m12 = n.xmlGetString('matrix_12')
        m13 = n.xmlGetString('matrix_13')
        m14 = n.xmlGetString('matrix_14')
        m21 = n.xmlGetString('matrix_21')
        m22 = n.xmlGetString('matrix_22')
        m23 = n.xmlGetString('matrix_23')
        m24 = n.xmlGetString('matrix_24')
        m31 = n.xmlGetString('matrix_31')
        m32 = n.xmlGetString('matrix_32')
        m33 = n.xmlGetString('matrix_33')
        m34 = n.xmlGetString('matrix_34')

        return m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34


    def setTransformationMatrix(self, perio_id, pos, val):
        """
        Put values of matrix of rotation for periodic transformation
        """
        self.isFloat(val)
        self.isInList(pos, ('matrix_11','matrix_12', 'matrix_13','matrix_14',
                            'matrix_21','matrix_22', 'matrix_23','matrix_24',
                            'matrix_31','matrix_32', 'matrix_33','matrix_34'))
        mode = self.getPeriodicityMode(perio_id)
        self.isInList(mode, ('mixed'))

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('mixed')
        n.xmlSetData(pos, val)


    def addPeriodicFaces(self, select):
        """
        Add faces selection for periodic transformation.
        Select is a dictionary with 'selector', 'fraction', 'plane', 'verbosity'
        """
        nb = self.getPeriodicSelectionsCount()
        name = str(nb +1)
        node = self.node_perio.xmlAddChild('face_periodicity', mode="", name=name)
        self._addJoinSelect(node, select)
        self.updatePeriodicityMode(nb, 'translation')


    def getPeriodicFaces(self, perio_id):
        """
        Public method.

        @return: faces selection for given periodic transformation
        @rtype: C{dictionary}
        """
        result = {}

        node = self._getPerioNode(perio_id)
        if node:
            result = self._getFaces(node)

        return result


    def replacePeriodicFaces(self, perio_id, select):
        """
        Replace values of faces selection for periodic transformation, by select
        """

        node = self._getPerioNode(perio_id)

        if node:
            self._removeJoinChildren(node)
            self._addJoinSelect(node, select)


# Methods to manage faces :
#========================

    def getJoinSelectionsCount(self):
        """
        Public method.

        @return: number of join faces selections
        @rtype: C{int}
        """
        return len(self.node_join.xmlGetNodeList('face_joining'))


    def addJoinFaces(self, select):
        """
        Add faces selection for face joining.
        Select is a dictionary with 'selector', 'fraction', 'plane', 'verbosity'
        """
        nb = self.getJoinSelectionsCount()
        name = str(nb +1)
        node = self.node_join.xmlAddChild('face_joining', name=name)
        self._addJoinSelect(node, select)


    def getJoinFaces(self, join_id):
        """
        Return faces selection named 'number' for face joining .
        """
        node = self._getJoinNode(join_id)
        return self._getFaces(node)


    def replaceJoinFaces(self, join_id, select):
        """
        Replace values of faces selection named 'number' for face joining, by select
        """
        node = self._getJoinNode(join_id)
        self._removeJoinChildren(node)
        self._addJoinSelect(node, select)


    def deleteJoinFaces(self, join_id):
        """
        Delete faces selection named 'number' for face joining
        """
        node = self._getJoinNode(join_id)
        node.xmlRemoveNode()
        if join_id < self.getJoinSelectionsCount():
            self._updateJoinSelectionNumbers()


# In following methods we build command for "runcase" file
#=========================================================

    def getMeshCommand(self):
        """
        Get mesh command line for preprocessor execution
        """
        lines = ''
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        mesh_mdl = MeshModel()
        for meshNode in nodeList:
            name   = meshNode['name']
            format = meshNode['format']
            #f = string.split(string.split(name,".")[1]," ")[0]
            #if f not in .getExtensionFileList():
            #         name=string.split(name," ")[0]
            mesh = self.case['mesh_path'] + '/' + name
            lines += " -m " + mesh

            if meshNode['num']:
                lines += " --num " + meshNode['num']
            if meshNode['grp_fac']:
                lines += " --grp-fac " + meshNode['grp_fac']
            if meshNode['grp_cel']:
                lines += " --grp-cel " + meshNode['grp_cel']

        lines += " "
        return lines


    def getReorientCommand(self):
        """
        Get reorient value for preprocessor execution
        """
        line = 'False'
        if self.node_orient and self.node_orient['status'] == 'on':
            line = 'True'

        return line


    def getSimCommCommand(self):
        """
        Get " --sim-comm " command line for preprocessor execution
        """
        lines = ''
        node = self.node_standalone.xmlGetNode('simulation_communication')
        if node and node['status'] == 'on':
            lines = " -sc "
        return lines


    def getPostCommand(self):
        """
        Get "--ensight" "--med" and/or "--cgns" command line for preprocessor execution
        """
        line  = ''
        iok = 0
        if self.getPostProFormat() == "EnSight":
            line = ' --ensight '
        if self.getPostProFormat() == "MED_fichier":
            line = ' --med '
        if self.getPostProFormat() == "CGNS":
            line = ' --cgns '

        options = self.getPostProOptionsFormat()
        options = string.split(options, ',')

        for opt in options:
            if opt in ('binary', 'discard_polygons', 'discard_polyhedra'):
                line = line
            if opt in ('divide_polygons', 'divide_polyhedra'):
                opt = line + ' --divide-poly '
            if opt == "big_endian":
                line = line + ' --big_endian '

        return line


# These following methods are kept only for tkInter View
#=======================================================

    def getListNodes(self, tagName):
        """
        Return node corresponding at the selection (only for view)
        """
        node = self._getTagNode(tagName)
        listNode = node.xmlGetNodeList(tagName)

        return listNode


    def getStatusNode(self, node):
        """
        Return status for node "node"
        """
        if node : status = node['status']
        return status


#-------------------------------------------------------------------------------
# SolutionDomain Model test case
#-------------------------------------------------------------------------------

class SolutionDomainTestCase(ModelTest):
    """
    """
    def checkSolutionDomainInstantiation(self):
        """ Check whether the SolutionDomainModel class could be instantiated """
        model = None
        model = SolutionDomainModel(self.case)
        assert model != None, 'Could not instantiate SolutionDomainModel'

    def checkAddDelMeshandGetMeshList(self):
        """ Check whether the meshes could be added and deleted and list of meshes could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.addMesh('fdc','des')
        mdl.addMesh('pic','des')
        mdl.addMesh('down','des')
        mdl.addMesh('up','des')
        doc = '<meshes_list>'\
                '<mesh format="des" name="fdc"/>'\
                '<mesh format="des" name="pic"/>'\
                '<mesh format="des" name="down"/>'\
                '<mesh format="des" name="up"/>'\
                '</meshes_list>'
        assert mdl.node_meshes == self.xmlNodeFromString(doc), \
            'Could not add meshes in SolutionDomainModel'
        mdl.delMesh('down')
        assert mdl.getMeshList() == ['fdc','pic','up'],\
            'Could not get mesh list'


    def checkSetandGetCutStatusAndAngleValue(self):
        """ Check whether the status of node cut and value of angle could be set and get"""
        mdl = SolutionDomainModel(self.case)
        mdl.setCutStatus('on')
        doc1 = '''<faces_cutting status="on"/>'''

        assert mdl.node_cut == self.xmlNodeFromString(doc1), \
            'Could not set status of faces_cutting'
        assert mdl.getCutStatus() == 'on',\
            'Could not get status of faces_cutting'

        mdl.setCutAngle(90.)
        doc2 = '''<faces_cutting status="on">
                    <warp_angle_max>90</warp_angle_max>
                  </faces_cutting>'''

        assert mdl.node_cut == self.xmlNodeFromString(doc2), \
            'Could not set angle for faces_cutting'
        assert mdl.getCutAngle() == 90, \
            'Could not get angle for faces_cutting'

    def checkSetandGetSimCommStatus(self):
        """ Check whether the status of node simulation_communication could be set and get"""
        mdl = SolutionDomainModel(self.case)
        mdl.setSimCommStatus('on')
        doc = '''<standalone>
                    <simulation_communication status="on"/>
                 </standalone>'''

        assert mdl.node_standalone == self.xmlNodeFromString(doc), \
            'Could not set status of node simulation_communication'
        assert mdl.getSimCommStatus() == 'on', \
            'Could not get status of node simulation_communication'

    def checkgetPeriodicSelectionsCount(self):
        """ Check whether the number of periodicities could be get"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.addPeriodicFaces(select)
        doc ='''<face_periodicity mode="translation" name="1">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                </face_periodicity>
                <face_periodicity mode="translation" name="2">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                </face_periodicity>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set number of periodicities'
        assert mdl.getPeriodicSelectionsCount() == 2,\
            'Could not get number for periodicities'

    def checkSetandgetPeriodicityMode(self):
        """ Check whether the mode of transformation could be set and get """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.addPeriodicFaces(select)
        mdl.updatePeriodicityMode('2', "rotation")
        doc ='''<face_periodicity mode="translation" name="1">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                </face_periodicity mode="rotation" name="2">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                      <rotation>
                            <angle>0</angle>
                            <axis_x>0</axis_x>
                            <axis_y>0</axis_y>
                            <axis_z>0</axis_z>
                            <invariant_x>0</invariant_x>
                            <invariant_y>0</invariant_y>
                            <invariant_z>0</invariant_z>
                      </rotation>
                </face_periodicity>'''


        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set mode of transformation for periodicities'
        assert mdl.getPeriodicityMode('2') == "rotation",\
            'Could not get mode of transformation for periodicities'

    def checkSetandgetTranslationDirection(self):
        """ Check whether the dir values translation mode of periodicity could be set and get"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.setTranslationDirection('1','translation_y',3.0)
        doc ='''<face_periodicity mode="translation" name="1">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>3</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                </face_periodicity>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set one direction values for translation'
        assert mdl.getTranslationDirection('1') == ('0', '3', '0'),\
            'Could not get one direction values for translation'

    def checkSetandgetRotationDirectionandAngleandCenter(self):
        """ Check whether the values for rotation's mode of periodicity could be set and get"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.setTranslationDirection('1','translation_y', 3.0)
        mdl.addPeriodicFaces(select)
        mdl.updatePeriodicityMode('2', "rotation")
        mdl.setRotationAngle('2', 180.)
        mdl.setRotationVector('2', 'axis_x', 0.5)
        mdl.setRotationVector('2', 'axis_z', 2.5)
        mdl.setRotationCenter('2', 'invariant_y', 9.8)
        doc ='''<face_periodicity mode="translation" name="1">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>3</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                </face_periodicity>
                <face_periodicity mode="rotation" name="2">
                      <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                      </translation>
                      <rotation>
                            <angle>180</angle>
                            <axis_x>0.5</axis_x>
                            <axis_y>0.0</axis_y>
                            <axis_z>2.5</axis_z>
                            <invariant_x>0</invariant_x>
                            <invariant_y>9.8</invariant_y>
                            <invariant_z>0</invariant_z>
                      </rotation>
                </face_periodicity>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set values for rotation transformation mode'
        assert mdl.getRotationAngle('2') == '180',\
            'Could not get value of angle for rotation transformation mode'
        assert mdl.getRotationDirection('2') == ('0.5', '0', '2.5'),\
            'Could not get values of direction for rotation transformation mode'
        assert mdl.getRotationCenter('2') == ('0', '9.8', '0'),\
            'Could not get value of center of rotation for rotation transformation mode'

    def checkSetandgetTransformationMatrix(self):
        """ Check whether the matrix of rotation for mixed mode could be set """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.updatePeriodicityMode('1','mixed')
        mdl.setTransformationMatrix('1', 'matrix_31', 31.31)
        doc = '''<face_periodicity mode="mixed" name="1">
                      <mixed>
                            <matrix_11>0.0</matrix_11>
                            <matrix_12>0.0</matrix_12>
                            <matrix_13>0.0</matrix_13>
                            <matrix_14>0.0</matrix_14>
                            <matrix_21>0.0</matrix_21>
                            <matrix_22>0.0</matrix_22>
                            <matrix_23>0.0</matrix_23>
                            <matrix_24>0.0</matrix_24>
                            <matrix_31>31.31</matrix_31>
                            <matrix_32>0.0</matrix_32>
                            <matrix_33>0.0</matrix_33>
                            <matrix_34>0.0</matrix_34>
                      </mixed>
                 </face_periodicity>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set values for matrix of rotation for mixed transformation mode'
        assert mdl.getTransformationMatrix('1') == ('0', '0', '0', '0',
                                                    '0', '0', '0', '0',
                                                    '31.31','0', '0', '0'),\
            'Could not get values for matrix of rotation for mixed transformation mode'

    def checkAddandGetJoinFaces(self):
        """ Check whether faces of face joining could be added and get """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addJoinFaces(select)
        doc = '''<face_joining name="1">
                    <selector>1 or 2 or 3 or toto</selector>
                    <fraction>0.1</fraction>
                    <plane>20</plane>
                    <verbosity>1</verbosity>
                 </face_joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for face joining'
        assert mdl.getJoinFaces('1') == {'selector': '1 or 2 or 3 or toto', 'plane': '20',
                                         'fraction': '0.1', 'verbosity': '1'},\
            'Could not get values of faces join for face joining'

    def checkReplaceandDeleteandSetandGetForJoinFaces(self):
        """
        Check whether faces of face joining could be replaced and deleted
        and status could be set and get
        """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = '1'
        deux = {}
        deux['selector'] = '9 or 8 or 7 or coucou'
        deux['fraction'] = '0.2'
        deux['plane'] = '20'
        deux['verbosity'] = '2'
        mdl = SolutionDomainModel(self.case)
        mdl.addJoinFaces(select)
        mdl.addJoinFaces(deux)
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                    </face_joining>
                    <face_joining name="2">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                    </face_joining>
                 </joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for face joining'
        assert mdl.getJoinFaces('1') == {'selector': '1 or 2 or 3 or toto',
                                        'plane': '20', 'fraction': '0.1',
                                        'verbosity': '1'},\
            'Could not get values of faces join for face joining'

        select['selector'] = 'je vais partir'
        mdl.replaceJoinFaces('1', select)
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                    </face_joining>
                    <face_joining name="2">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                    </face_joining>
                 </joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not replace values of faces join for face joining'

        mdl.deleteJoinFaces('1')
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                    </face_joining>
                 </joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not delete faces join for face joining'

        mdl.addJoinFaces(select)
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                    </face_joining>
                    <face_joining name="2">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                    </face_joining>
                 </joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set faces join for joinings'

    def checkAddandGetPeriodicFaces(self):
        """ Check whether faces of periodicity could be added and get """
        select = {}
        select['selector'] = '5 or 6 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = '2'
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.addPeriodicFaces(select)
        mdl.updatePeriodicityMode('2', 'rotation')
        doc = '''<face_periodicity>
                 <face_periodicity mode="translation" name="1">
                        <selector>5 or 6 or toto</selector>
                        <fraction>0.1</fraction>
                        <plane>0.8</plane>
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                 </face_periodicity>
                 <face_periodicity mode="rotation" name="2">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                        <rotation>
                            <angle>0</angle>
                            <axis_x>0</axis_x>
                            <axis_y>0</axis_y>
                            <axis_z>0</axis_z>
                            <invariant_x>0</invariant_x>
                            <invariant_y>0</invariant_y>
                            <invariant_z>0</invariant_z>
                        </rotation>
                        <face_periodicity status="on">
                            <selector>5 or 6 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>0.8</plane>
                        </face_periodicity>
                 </face_periodicity>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not add values of faces for periodicities'
        assert mdl.getPeriodicFaces('1') == {'selector': '5 or 6 or toto',
                                        'plane': '30', 'fraction': '0.1', 'verbosity': '2'},\
            'Could not get values of faces for periodicities'

    def checkReplaceandDeleteandSetandGetStatusForPeriodicFaces(self):
        """
        Check whether faces of of periodicity could be replaced and deleted
        and status could be set and get
        """
        select = {}
        select['selector'] = '5 or 6 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '25'
        select['verbosity'] = '1'
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicFaces(select)
        mdl.addPeriodicFaces(select)
        mdl.updatePeriodicityMode('2', 'rotation')
        doc = '''<face_periodicity mode="translation" name="1">
                       <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                       </translation>
                       <selector>5 or 6 or toto</selector>
                       <fraction>0.1</fraction>
                       <plane>25</plane>
                       <verbosity>1</verbosity>
                 </face_periodicity>
                 <face_periodic mode="rotation" name="2">
                       <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                       </translation>
                       <rotation>
                            <angle>0.0</angle>
                            <axis_x>0.0</axis_x>
                            <axis_y>0.0</axis_y>
                            <axis_z>0.0</axis_z>
                            <invariant_x>0.0</invariant_x>
                            <invariant_y>0.0</invariant_y>
                            <invariant_z>0.0</invariant_z>
                       </rotation>
                 </face_periodicity>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not delete one selection of faces for periodicities'

        select['selector'] = '147 or 963 or PERIODIC'
        select['fraction'] = '0.1'
        select['plane']  = '20'
        select['verbosity'] = '2'
        mdl.replacePeriodicFaces('1', select)
        doc = '''<face_periodicity mode="translation" name="1">
                     <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                     </translation>
                     <selector>147 or 963 or PERIODIC</selector>
                     <fraction>0.1</fraction>
                     <plane>30</plane>
                     <verbosity>1</verbosity>
                 </face_periodicity>
                 </face_periodicity mode="rotation" name="2">
                     <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                     </translation>
                     <rotation>
                            <angle>0.0</angle>
                            <axis_x>0.0</axis_x>
                            <axis_y>0.0</axis_y>
                            <axis_z>0.0</axis_z>
                            <invariant_x>0.0</invariant_x>
                            <invariant_y>0.0</invariant_y>
                            <invariant_z>0.0</invariant_z>
                     </rotation>
                 </face_periodicity>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not replace values of faces for periodicities'

        doc = '''<face_periodicity mode="translation" name="1">
                      <selector>147 or 963 or PERIODIC</selector>
                      <fraction>0.1</fraction>
                      <plane>30</plane>
                      <verbosity>1</verbosity>
                      <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                      </translation>
                 </face_periodicity>
                 <face_periodicity mode="rotation" name="2">
                      <translation>
                            <translation_x>0.0</translation_x>
                            <translation_y>0.0</translation_y>
                            <translation_z>0.0</translation_z>
                      </translation>
                      <rotation>
                            <angle>0.0</angle>
                            <axis_x>0.0</axis_x>
                            <axis_y>0.0</axis_y>
                            <axis_z>0.0</axis_z>
                            <invariant_x>0.0</invariant_x>
                            <invariant_y>0.0</invariant_y>
                            <invariant_z>0.0</invariant_z>
                     </rotation>
                 </face_periodicity>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set selection of faces for periodicities'


    def checkMeshCommand(self):
        """ Check whether  command for meshes could be set """
        mdl = SolutionDomainModel(self.case)
        mdl.case['mesh_path'] = 'MESH'
        mdl.addMesh('fdc.des','des')
        mdl.addMesh('pic.des','des')
        mdl.addMesh('truc.ngeom','ngeom')
        line = ''' -m MESH/fdc.des -m MESH/pic.des -m MESH/truc.ngeom '''

        assert mdl.getMeshCommand() == line,\
            'Mesh command is not verified in SolutionDomain Model'


    def checkReorientSetAndGetStatusAndCommand(self):
        """ Check whether reorient status could be set and get and command line could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.setOrientation('on')
        doc = '''<reorientation status="on"/>'''

        assert mdl.node_orient == self.xmlNodeFromString(doc),\
            'Could not set reorient status in SolutionDomain Model'
        assert mdl.getOrientation() == "on",\
            'Could not get reorient status in SolutionDomain Model'

        cmd_orient = ' --reorient '
        assert mdl.getReorientCommand() == cmd_orient,\
            'Reorient command is not verified in SolutionDomain Model'


    def checkSimCommAndVerifMaillCommand(self):
        """ Check whether simulation_communication command line could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.setSimCommStatus('on')
        cmd_sim = mdl.getSimCommCommand()
        sim =' -sc '
        assert mdl.getSimCommCommand() == sim,\
            'Simulation_communication command is not verified in SolutionDomain Model'

    def checkPostCommand(self):
        """Check whether output postprocessing format command line could be get"""
        mdl = SolutionDomainModel(self.case)
        mdl.setPostProFormat('CGNS')

        cmd_post = ' --cgns '
        assert mdl.getPostCommand() == cmd_post,\
            'Post command is not verified for postprocessing format in SolutionDomain Model'


#-------------------------------------------------------------------------------
# SolutionDomain Model test case
#-------------------------------------------------------------------------------


class MeshModelTestCase(unittest.TestCase):
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        self.files = [ ("toto.case",     "case") ,
                       ("toto.cgns.gz",  "cgns") ,
                       ("toto.des",      "des")  ,
                       ("toto.ccm.gz",   "ccm")  ,
                       ("toto.med",      "med")  ,
                       ("toto.msh.gz",   "msh")  ,
                       ("toto.neu",      "neu")  ,
                       ("toto.ngeom.gz", "ngeom"),
                       ("toto.unv.gz",   "unv")  ,
                       ("toto",          "")     ,
                       ("toto.gz",       "")     ]

    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        pass

    def checkGetMeshExtension(self):
        """Check whether mesh extension could be get"""
        mdl = MeshModel()

        for f in self.files:
          ext = mdl.getMeshExtension(f[0])
          assert ext == f[1], 'could not get the mesh extension'

    def checkGetMeshFormat(self):
        """Check whether mesh extension could be get"""
        mdl = MeshModel()

        for f in self.files:
          fmt = mdl.getMeshFormat(f[0])
          if fmt:
              assert fmt == mdl.ext[f[1]], 'could not get the mesh format'


def suite1():
    testSuite = unittest.makeSuite(SolutionDomainTestCase, "check")
    return testSuite

def suite2():
    testSuite = unittest.makeSuite(MeshModelTestCase, "check")
    return testSuite

def runTest():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite1())
    runner.run(suite2())





