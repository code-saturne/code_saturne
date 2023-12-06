# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This module defines the XML calls for preprocessor execution
This module contains the following classes and function:
- MeshModel
- SolutionDomainModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
import os, sys, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Utility function
#-------------------------------------------------------------------------------

def RelOrAbsPath(path, case_dir):
    """
    Return a relative filepath in a same study, an absolute path otherwise.
    """

    study_dir = os.path.split(case_dir)[0]

    if path.find(study_dir) == 0:

        if hasattr(os.path, 'relpath'):
            return os.path.relpath(path, case_dir)

        elif path.find(case_dir) == 0:
            return path[len(case_dir)+1:]

        else:
            return os.path.join('..', path[len(study_dir)+1:])

    else:
        return path

#-------------------------------------------------------------------------------
# Utility function: get run type from case
#-------------------------------------------------------------------------------

def getRunType(case):
    """
    This function allows querying a case for its run type
    """

    val = None

    node_mgt = case.xmlGetNode('calculation_management')
    if node_mgt:
        val = node_mgt.xmlGetString('run_type')

    if not val:
        val = 'standard'

    case['run_type'] = val

    return val


#-------------------------------------------------------------------------------
# Utility function: get run type from case
#-------------------------------------------------------------------------------

def getMeshOriginType(case):

    val = None

    node_sold = case.xmlGetNode('solution_domain')
    if node_sold:
        node = node_sold.xmlGetChildNode('mesh_origin')
        if node:
            val = node['choice']

    return val


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
        self.ext = {'case':'ensight',
                    'cgns':'cgns',
                    'csm':'code_saturne',
                    'des':'des',
                    'med':'med',
                    'msh':'gmsh',
                    'neu':'gambit',
                    'ccm':'ccm',
                    'unv':'ideas'}


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
        last_caracters = (mesh.split(".")[-1:])[0]
        if last_caracters in list(self.ext.keys()):
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
        return list(self.ext.keys())


    def getBuildFormatList(self):
        """
        Public method.

        @return: List of number, format and description for view of popup.
        @rtype: C{List} of C{2-tuple}
        """
        lst = [('ensight',      'EnSight',           ' (*.case)' ),
               ('cgns',         'CGNS',              ' (*.cgns)' ),
               ('code_saturne', 'code_saturne mesh', ' (*.csm)'  ),
               ('des',          'Simail/NOPO',       ' (*.des)'  ),
               ('med',          'MED',               ' (*.med)'  ),
               ('gmsh',         'Gmsh',              ' (*.msh)'  ),
               ('gambit',       'GAMBIT Neutral',    ' (*.neu)'  ),
               ('ccm',          'STAR-CCM+',         ' (*.ccm)'  ),
               ('ideas',        'I-deas universal',  ' (*.unv)'  )]

        return lst


    def getFileFormatList(self):
        """
        Public method.

        @return: List of format and associated text for research file.
        @rtype: C{List} of C{2-tuple}
        """
        lst = [("All files",                 "*"      ),
               ("EnSight (6 or Gold) files", "*.case" ),
               ("CGNS files",                "*.cgns" ),
               ("code_saturne mesh files",    "*.csm"  ),
               ("Simail (NOPO) files",       "*.des"  ),
               ("MED files",                 "*.med"  ),
               ("GMSH files",                "*.msh"  ),
               ("GAMBIT Neutral files",      "*.neu"  ),
               ("STAR-CCM+",                 "*.ccm"  ),
               ("I-deas universal files",    "*.unv"  )]

        return lst


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
        self.node_cut        = self.node_ecs.xmlInitNode('faces_cutting')
        self.node_smooth     = self.node_ecs.xmlInitNode('mesh_smoothing')
        self.node_join       = self.node_ecs.xmlInitNode('joining')
        self.node_perio      = self.node_ecs.xmlInitNode('periodicity')
        self.node_thinwall   = self.node_ecs.xmlInitNode('thin_walls')
        self.node_extrude    = self.node_ecs.xmlInitNode('extrusion')

        self.node_mgt = self.case.xmlInitNode('calculation_management')


#************************* Private methods *****************************

    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['cutting_status'] = "off"
        defvalue['smooth_status']  = "off"
        defvalue['select_status']  = "off"
        defvalue['selector']       = "all[]"
        defvalue['fraction']       = 0.1
        defvalue['plane']          = 25.0
        defvalue['verbosity']      = 1
        defvalue['visualization']  = 1
        defvalue['angle']          = 0.01
        defvalue['smooth_angle']   = 25.
        defvalue['syrth_status']   = "off"
        defvalue['syrth_mesh_2d']  = "off"
        defvalue['sim_status']     = "on"
        defvalue['verif_mail']     = "on"
        defvalue['dir_cas']        = "default_case"
        defvalue['poly_status']    = "off"
        defvalue['perio_mode']     = "translation"
        defvalue['transfo_val']    = 0.0
        defvalue['layer']          = 2
        defvalue['thickness']      = 1.0
        defvalue['reason']         = 1.5
        defvalue['mesh_origin']    = "mesh_import"

        return defvalue

    def defaultValuesCartesian(self):
        """
        Return a dictionary with default values
        """
        d = {"ncells": "1",
             "min": "0.0",
             "max": "1.0",
             "prog": "1.0",
             "law": "constant"}

        return d

    def _getMeshNode(self, mesh):
        """
        Public method. Return the node matching a mesh.
        """
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            name = node['name']
            path = node['path']
            if path == '':
                path = None
            if (name, path) == mesh:
                return node

        msg = "There is an error: this value " + str(mesh) + "\n"\
            "is not in list " + str(nodeList) + "\n"
        raise ValueError(msg)


#To follow : private methods to get or put faces
#================================================

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
                          (select['verbosity'], 'verbosity'),
                          (select['visualization'], 'visualization')]:
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
        default['visualization'] = ""

        if node:
            default['selector']  = node.xmlGetString('selector')
            default['fraction']  = node.xmlGetString('fraction')
            default['plane']     = node.xmlGetString('plane')
            default['verbosity'] = node.xmlGetString('verbosity')
            default['visualization'] = node.xmlGetString('visualization')
            if not default['selector']:
                default['selector'] = "all[]"
            if not default['fraction']:
                default['fraction'] = 0.1
            if not default['plane']:
                default['plane'] = 25
            if not default['verbosity']:
                default['verbosity'] = 1
            if not default['visualization']:
                default['visualization'] = 1

        else:
            default = {}

        return default


    def _removeJoinChildren(self, node):
        """
        Private method: Remove all child nodes of node for one selection
        """
        for tag in ('selector',
                    'fraction',
                    'plane',
                    'verbosity',
                    'visualization',):
            node.xmlRemoveChild(tag)


#To follow : private methods to get or put thin wall
#===================================================

    def _getThinWallNode(self, thin_id):
        """
        Get node for a given joining
        """
        node = None
        listNode = self.node_thinwall.xmlGetNodeList('thin_wall')
        if thin_id < len(listNode):
            node = listNode[thin_id]

        return node

    def _updateThinWallSelectionNumbers(self):
        """
        Update names of join selection
        """
        listNode = self.node_thinwall.xmlGetNodeList('thin_wall')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    def _addThinWallSelect(self, node, select):
        """
        Private method: Add faces to node.
        """
        node.xmlSetData('selector', select)


    def _removeThinWallChildren(self, node):
        """
        Private method: Remove all child nodes of node for one selection
        """
        node.xmlRemoveChild('selector')


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

# Methods to manage mesh choice (mesh_input, mesh import or cartesian
#======================================

    @Variables.noUndo
    def getMeshOrigin(self):
        """
        Public method. Return the kind of mesh to use (preprocessed, to import
        or cartesian generated by the code.
        """
        mesh_origin = self.node_ecs.xmlGetNode('mesh_origin')
        if mesh_origin:
            return mesh_origin['choice']
        else:
            mi = self.getMeshInput()
            if mi:
                return 'mesh_input'
            else:
                return 'mesh_import'

    @Variables.undoLocal
    def setMeshOrigin(self, choice):

        node = self.node_ecs.xmlInitNode('mesh_origin', 'choice')
        if choice not in ('mesh_import', 'mesh_input', 'mesh_cartesian'):
            choice = 'mesh_import'

        node['choice'] = choice

        if choice != 'mesh_import':
            for m in self.getMeshList():
                self.delMesh(m)

        if choice != 'mesh_input':
            self.delMeshInput()

        if choice != 'mesh_cartesian':
            n = self.node_ecs.xmlGetNode('mesh_cartesian')
            if n:
                n.xmlRemoveNode()


# Methods to manage the cartesian mesh
#======================================

    @Variables.noUndo
    def getCartesianParam(self, d, p):

        if self.getMeshOrigin() != 'mesh_cartesian':
            return None

        node = self.node_ecs.xmlGetNode('mesh_cartesian')
        val = None
        if node:
            n = node.xmlGetChildNode(d)
            if n:
                if p in n.xmlGetAttributeDictionary():
                    val = n[p]

        if val is None:
            val = self.defaultValuesCartesian()[p]

        return val


    @Variables.undoLocal
    def setCartesianParam(self, d, p, val):

        node = self.node_ecs.xmlGetNode('mesh_cartesian')
        if node is None:
            node = self.node_ecs.xmlInitNode('mesh_cartesian')

        n = node.xmlGetChildNode(d)
        if n is None:
            n = node.xmlInitNode(d)

        n[p] = val


# Methods to manage the mesh_input path
#======================================

    @Variables.noUndo
    def getMeshInput(self):
        """
        Public method. Return the mesh_input file or directory path.
        """
        mesh_input = self.node_ecs.xmlGetNode('mesh_input', 'path')
        if mesh_input:
            return mesh_input['path']
        else:
            return None


    @Variables.undoLocal
    def delMeshInput(self):
        """
        Public method. Delete node for mesh_input named "mesh_input" in xml file
        """
        node = self.node_ecs.xmlGetNode('mesh_input', 'path')
        if node != None:
            node.xmlRemoveNode()


    @Variables.undoLocal
    def setMeshInput(self, mesh_input):
        """
        Public method. Add mesh_input path name in xml file.
        """

        if mesh_input == '':
            mesh_input = None

        node = self.node_ecs.xmlInitNode('mesh_input', 'path')
        if mesh_input:
            node['path'] = mesh_input
            if self.getRunType() == 'none':
                self.setRunType('mesh preprocess')
        else:
            node.xmlRemoveNode()


# Methods to manage meshes :
#=========================

    @Variables.undoLocal
    def addMesh(self, mesh):
        """
        Public method. Add mesh name in xml file.
        """
        self.isNotInList(mesh, self.getMeshList())

        if mesh[1] != None:
            self.node_meshes.xmlInitNode('mesh', name=mesh[0], path=mesh[1])
        else:
            self.node_meshes.xmlInitNode('mesh', name=mesh[0])


    @Variables.undoLocal
    def delMesh(self, mesh):
        """
        Public method. Delete node for mesh named "mesh" in xml file
        """
        node = self._getMeshNode(mesh)
        node.xmlRemoveNode()


    def getMeshList(self):
        """
        Public method. Return the meshes name list already put in the case.
        """
        meshList = []
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            name = node['name']
            path = node['path']
            if path == '':
                path = None
            meshList.append((name, path))
        return meshList


    @Variables.undoLocal
    def setMeshFormat(self, mesh, format):
        """
        Public method. Set the mesh format.
        """
        node = self._getMeshNode(mesh)
        if not format or format == MeshModel().getMeshFormat(mesh[0]):
            del node['format']
        else:
            node['format'] = format


    @Variables.noUndo
    def getMeshFormat(self, mesh):
        """
        Public method. Return the mesh format recorded in the case.
        """
        node = self._getMeshNode(mesh)
        fmt = node['format']
        if not fmt:
            fmt = MeshModel().getMeshFormat(mesh[0])
        return fmt


    @Variables.undoLocal
    def setMeshNumbers(self, mesh, num):
        """
        Public method. Set the mesh number(s).
        """
        node = self._getMeshNode(mesh)
        if not num:
            del node['num']
        else:
            node['num'] = num


    @Variables.noUndo
    def getMeshNumbers(self, mesh):
        """
        Public method. Return the mesh number recorded in the case.
        """
        node = self._getMeshNode(mesh)
        return node['num']


    @Variables.undoLocal
    def setMeshGroupCells(self, mesh, grp_cel):
        """
        Public method. Put the grp-cel option.
        """
        node = self._getMeshNode(mesh)
        self.isInList(grp_cel, ('off', 'section', 'zone'))

        if grp_cel == "off":
            del node['grp_cel']
        else:
            node['grp_cel'] = grp_cel


    @Variables.undoLocal
    def setMeshDiscardBadCells(self, mesh, discard_bad_cells):
        """
        Sets the --discard-bad-cells option
        """
        node = self._getMeshNode(mesh)

        if discard_bad_cells == False:
            del node['discard_bad_cells']
        else:
            node['discard_bad_cells'] = 'on'


    @Variables.noUndo
    def getMeshDiscardBadCells(self, mesh):
        """
        Returns the --discard-bad-cells option
        """
        node = self._getMeshNode(mesh)
        discard_bad_cells = False
        if node['discard_bad_cells'] == 'on':
            discard_bad_cells = True
        return discard_bad_cells


    @Variables.undoLocal
    def setMeshReorient(self, mesh, reorient):
        """
        Public method. Put the grp-cel option.
        """
        node = self._getMeshNode(mesh)

        if reorient == False:
            del node['reorient']
        else:
            node['reorient'] = 'on'


    @Variables.noUndo
    def getMeshReorient(self, mesh):
        """
        Public method. Return the mesh 'grp-cel' sub-option recorded in the case.
        """
        node = self._getMeshNode(mesh)
        reorient = False
        if node['reorient'] == 'on':
            reorient = True
        return reorient


    @Variables.noUndo
    def getMeshGroupCells(self, mesh):
        """
        Public method. Return the mesh 'grp-cel' sub-option recorded in the case.
        """
        return self.__getMeshGroup(mesh, 'grp_cel')


    @Variables.undoLocal
    def setMeshGroupFaces(self, mesh, grp_fac):
        """
        Public method. Put the 'grp-fac' sub-option.
        """
        node = self._getMeshNode(mesh)
        self.isInList(grp_fac, ('off', 'section', 'zone'))

        if grp_fac == "off":
            del node['grp_fac']
        else:
            node['grp_fac'] = grp_fac


    @Variables.noUndo
    def getMeshGroupFaces(self, mesh):
        """
        Public method. Return the mesh 'grp_fac' option recorded in the case.
        """
        return self.__getMeshGroup(mesh, 'grp_fac')


    def __getMeshGroup(self, mesh, group):
        """
        Private method. Return the mesh 'grp_fac' or 'grp_cel' sub-option recorded in the case.
        """
        node = self._getMeshNode(mesh)
        grp = node[group]
        if grp is None:
            grp = 'off'
        return grp


    @Variables.noUndo
    def getMeshDir(self):
        """
        Public method. Return the meshdir directory name.
        """
        meshnode = self.node_meshes.xmlGetNode('meshdir', 'name')
        if meshnode != None:
            meshdir = meshnode['name']
            if not os.path.isabs(meshdir):
                meshdir = os.path.join(self.case['case_path'], meshdir)
            meshdir = os.path.abspath(meshdir)
            return meshdir
        else:
            return None


    @Variables.undoLocal
    def setMeshDir(self, mesh_dir):
        """
        Public method. Add mesh name in xml file.
        """

        case_dir = self.case['case_path']

        if mesh_dir:
            if not os.path.isabs(mesh_dir):
                mesh_dir = os.path.join(case_dir, mesh_dir)
            mesh_dir = os.path.abspath(mesh_dir)

        node = self.node_meshes.xmlInitNode('meshdir', 'name')
        if mesh_dir:
            study_dir = os.path.split(case_dir)[0]
            node['name'] = RelOrAbsPath(mesh_dir, case_dir)
        else:
            node.xmlRemoveNode()

        old_mesh_dir = self.case['mesh_path']
        self.case['mesh_path'] = mesh_dir

        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            name = node['name']
            path = node['path']
            # Rebuild absolute name
            if path != None and path !='':
                name = os.path.join(path, name)
            if not os.path.isabs(name) and old_mesh_dir != None:
                name = os.path.join(old_mesh_dir, name)
            # Split components
            if os.path.isfile(name):
                path = os.path.dirname(name)
                if mesh_dir != None:
                    index = path.find(mesh_dir)
                    if index == 0:
                        path = path[len(mesh_dir)+1:]
                if len(path) > 0:
                    node['path'] = path
                else:
                    del node['path']



# Methods to manage status of all main tags :
#==========================================

    @Variables.noUndo
    def getCutStatus(self):
        """
        Get status on tag "faces_cutting" from xml file
        """
        status = self.node_cut['status']
        if not status:
            status = self.defaultValues()['cutting_status']
            self.setCutStatus(status)
        return status


    @Variables.undoLocal
    def setCutStatus(self, status):
        """
        Put status on tag "faces_cutting" in xml file
        """
        self.isOnOff(status)
        self.node_cut['status'] = status


    @Variables.undoLocal
    def setCutAngle(self, var):
        """
        input '--cut_warped_faces' parameter.
        """
        self.isGreaterOrEqual(var, 0.0)
        if var != self.defaultValues()['angle']:
            self.node_cut.xmlSetData('warp_angle_max', var)
        else:
            self.node_cut.xmlRemoveChild('warp_angle_max')


    @Variables.noUndo
    def getCutAngle(self):
        """
        get '--cut_warped_faces' parameters.
        """
        angle = self.node_cut.xmlGetDouble('warp_angle_max')
        if angle is None:
            angle = self.defaultValues()['angle']
        return angle


    @Variables.noUndo
    def getSmoothingStatus(self):
        """
        Get status on tag "mesh_smoothing" from xml file
        """
        status = self.node_smooth['status']
        if not status:
            status = self.defaultValues()['smooth_status']
            self.setSmoothingStatus(status)
        return status


    @Variables.undoLocal
    def setSmoothingStatus(self, status):
        """
        Put status on tag "mesh_smoothing" in xml file
        """
        self.isOnOff(status)
        self.node_smooth['status'] = status


    @Variables.undoLocal
    def setSmoothAngle(self, var):
        """
        input '--mesh_smoothing' parameter.
        """
        self.isGreaterOrEqual(var, 0.0)
        self.isLowerOrEqual(var, 90.0)
        if var != self.defaultValues()['smooth_angle']:
            self.node_smooth.xmlSetData('smooth_angle', var)
        else:
            self.node_smooth.xmlRemoveChild('smooth_angle')


    @Variables.noUndo
    def getSmoothAngle(self):
        """
        get '--mesh_smoothing' parameters.
        """
        angle = self.node_smooth.xmlGetDouble('smooth_angle')
        if angle is None:
            angle = self.defaultValues()['smooth_angle']
        return angle


# Methods to manage periodicity :
#==============================

    @Variables.noUndo
    def getPeriodicSelectionsCount(self):
        """
        Public method.

        @return: number of periodic faces selections
        @rtype: C{int}
        """
        return len(self.node_perio.xmlGetNodeList('face_periodicity'))


    @Variables.noUndo
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


    @Variables.undoGlobal
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


    @Variables.undoGlobal
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


    @Variables.noUndo
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


    @Variables.undoLocal
    def setTranslationDirection(self, perio_id, dir, valcoor):
        """
        Put values of translation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(dir, ('translation_x', 'translation_y', 'translation_z'))

        node = self._getPerioNode(perio_id)

        for n in node.xmlGetChildNodeList('translation'):
            n.xmlSetData(dir, valcoor)


    @Variables.noUndo
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


    @Variables.undoLocal
    def setRotationVector(self, perio_id, dir, valcoor):
        """
        Put values for director vector rotation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(dir, ("axis_x", "axis_y", "axis_z"))

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData(dir,valcoor)


    @Variables.noUndo
    def getRotationAngle(self, perio_id):
        """
        Get angle for rotation for periodic rotation
        """
        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        angle = n.xmlGetString('angle')

        return angle


    @Variables.undoLocal
    def setRotationAngle(self, perio_id, angle):
        """
        Put angle for rotation for periodic rotation
        """
        self.isGreaterOrEqual(angle, 0.0)

        node = self._getPerioNode(perio_id)

        n = node.xmlGetChildNode('rotation')
        n.xmlSetData('angle', angle)


    @Variables.noUndo
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


    @Variables.undoGlobal
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


    @Variables.noUndo
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


    @Variables.undoGlobal
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


    @Variables.undoGlobal
    def addPeriodicFaces(self, select):
        """
        Add faces selection for periodic transformation.
        Select is a dictionary with 'selector', 'fraction', 'plane', 'verbosity', 'visualization'
        """
        nb = self.getPeriodicSelectionsCount()
        name = str(nb +1)
        node = self.node_perio.xmlAddChild('face_periodicity', mode="", name=name)
        self._addJoinSelect(node, select)
        self.updatePeriodicityMode(nb, 'translation')


    @Variables.noUndo
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


    @Variables.undoGlobal
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

    @Variables.noUndo
    def getJoinSelectionsCount(self):
        """
        Public method.

        @return: number of join faces selections
        @rtype: C{int}
        """
        return len(self.node_join.xmlGetNodeList('face_joining'))


    @Variables.undoGlobal
    def addJoinFaces(self, select):
        """
        Add faces selection for face joining.
        Select is a dictionary with 'selector', 'fraction', 'plane', 'verbosity', 'visualization'
        """
        nb = self.getJoinSelectionsCount()
        name = str(nb +1)
        node = self.node_join.xmlAddChild('face_joining', name=name)
        self._addJoinSelect(node, select)


    @Variables.noUndo
    def getJoinFaces(self, join_id):
        """
        Return faces selection named 'number' for face joining .
        """
        node = self._getJoinNode(join_id)
        return self._getFaces(node)


    @Variables.undoGlobal
    def replaceJoinFaces(self, join_id, select):
        """
        Replace values of faces selection named 'number' for face joining, by select
        """
        node = self._getJoinNode(join_id)
        self._removeJoinChildren(node)
        self._addJoinSelect(node, select)


    @Variables.undoGlobal
    def deleteJoinFaces(self, join_id):
        """
        Delete faces selection named 'number' for face joining
        """
        node = self._getJoinNode(join_id)
        node.xmlRemoveNode()
        if join_id < self.getJoinSelectionsCount():
            self._updateJoinSelectionNumbers()


# Methods to manage thin wall :
#==============================

    @Variables.noUndo
    def getThinWallSelectionsCount(self):
        """
        Public method.

        @return: number of thin wall selections
        @rtype: C{int}
        """
        return len(self.node_thinwall.xmlGetNodeList('thin_wall'))


    @Variables.undoGlobal
    def addThinWall(self):
        """
        Add faces selection for thin wall.
        Select is a dictionary with 'id', 'visualization'
        """
        nb = self.getThinWallSelectionsCount()
        name = str(nb +1)
        node = self.node_thinwall.xmlAddChild('thin_wall', name=name)
        select = self.defaultValues()['selector']
        self._addThinWallSelect(node, select)


    @Variables.noUndo
    def getThinWall(self, thin_id):
        """
        Return faces selection named 'number' for thin wall.
        """
        node = self._getThinWallNode(thin_id)
        return node.xmlGetString('selector')


    @Variables.undoGlobal
    def replaceThinWall(self, thin_id, select):
        """
        Replace values of faces selection named 'number' for thin wall, by select
        """
        node = self._getThinWallNode(thin_id)
        self._removeThinWallChildren(node)
        self._addThinWallSelect(node, select)


    @Variables.undoGlobal
    def deleteThinWall(self, thin_id):
        """
        Delete faces selection named 'number' for thin wall
        """
        node = self._getThinWallNode(thin_id)
        node.xmlRemoveNode()
        if thin_id < self.getThinWallSelectionsCount():
            self._updateThinWallSelectionNumbers()


# Methods to manage extrude :
#==============================
    def _getExtrudeNode(self, ext_id):
        """
        """
        node = None
        listNode = self.node_extrude.xmlGetNodeList('extrude_mesh')
        if ext_id < len(listNode):
            node = listNode[ext_id]
        return node

    def _updateExtrudeNumbers(self):
        """
        Update names of join selection
        """
        listNode = self.node_extrude.xmlGetNodeList('extrude_mesh')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    @Variables.noUndo
    def getExtrudeSelectionsCount(self):
        """
        Public method.

        @return: number of thin wall selections
        @rtype: C{int}
        """
        return len(self.node_extrude.xmlGetNodeList('extrude_mesh'))


    @Variables.undoGlobal
    def addExtrude(self):
        """
        Add faces selection for thin wall.
        Select is a dictionary with 'id', 'visualization'
        """
        nb = self.getExtrudeSelectionsCount()
        name = str(nb +1)
        node = self.node_extrude.xmlAddChild('extrude_mesh', name=name)
        select = self.defaultValues()['selector']
        node.xmlSetData('selector', select)
        select = self.defaultValues()['layer']
        node.xmlSetData('layers_number', select)
        select = self.defaultValues()['thickness']
        node.xmlSetData('thickness', select)
        select = self.defaultValues()['reason']
        node.xmlSetData('reason', select)


    @Variables.noUndo
    def getExtrudeSelector(self, thin_id):
        """
        Return faces selection named 'number' for thin wall.
        """
        node = self._getExtrudeNode(thin_id)
        return node.xmlGetString('selector')


    @Variables.undoGlobal
    def setExtrudeSelector(self, thin_id, select):
        """
        Replace values of faces selection named 'number' for thin wall, by select
        """
        node = self._getExtrudeNode(thin_id)
        node.xmlRemoveChild('selector')
        node.xmlSetData('selector', select)


    @Variables.noUndo
    def getExtrudeLayer(self, thin_id):
        """
        Return faces selection named 'number' for thin wall.
        """
        node = self._getExtrudeNode(thin_id)
        return node.xmlGetString('layers_number')


    @Variables.undoGlobal
    def setExtrudeLayer(self, thin_id, select):
        """
        Replace values of faces selection named 'number' for thin wall, by select
        """
        node = self._getExtrudeNode(thin_id)
        node.xmlSetData('layers_number', select)


    @Variables.noUndo
    def getExtrudeThickness(self, thin_id):
        """
        Return faces selection named 'number' for thin wall.
        """
        node = self._getExtrudeNode(thin_id)
        return node.xmlGetString('thickness')


    @Variables.undoGlobal
    def setExtrudeThickness(self, thin_id, select):
        """
        Replace values of faces selection named 'number' for thin wall, by select
        """
        node = self._getExtrudeNode(thin_id)
        node.xmlSetData('thickness', select)


    @Variables.noUndo
    def getExtrudeReason(self, thin_id):
        """
        Return faces selection named 'number' for thin wall.
        """
        node = self._getExtrudeNode(thin_id)
        return node.xmlGetString('reason')


    @Variables.undoGlobal
    def setExtrudeReason(self, thin_id, select):
        """
        Replace values of faces selection named 'number' for thin wall, by select
        """
        node = self._getExtrudeNode(thin_id)
        node.xmlSetData('reason', select)


    @Variables.undoGlobal
    def deleteExtrude(self, thin_id):
        """
        Delete faces selection named 'number' for thin wall
        """
        node = self._getExtrudeNode(thin_id)
        node.xmlRemoveNode()
        if thin_id < self.getExtrudeSelectionsCount():
            self._updateExtrudeNumbers()


    # Methods to manage restart behavior
    #===================================

    @Variables.noUndo
    def getPreprocessOnRestart(self):
        """
        Get preprocess on restart option.
        """

        val = False
        if self.node_ecs:
            s = self.node_ecs.xmlGetString('preprocess_on_restart')
            if str(s) == 'yes':
                val = True

        return val


    @Variables.undoLocal
    def setPreprocessOnRestart(self, val):
        """
        Set run type.
        """
        if val == True:
            self.node_ecs.xmlSetData('preprocess_on_restart', 'yes')
        else:
            node = self.node_ecs.xmlGetNode('preprocess_on_restart')
            if node:
                node.xmlRemoveNode()


    # Methods to manage run type
    #===========================

    @Variables.noUndo
    def getRunType(self):
        """
        Get run type.
        """
        val = self.node_mgt.xmlGetString('run_type')

        if not val:
            val = 'standard'

        self.case['run_type'] = val

        return val


    @Variables.undoLocal
    def setRunType(self, run):
        """
        Set run type.
        """
        self.isInList(run, ('none', 'mesh preprocess', 'mesh quality', 'standard'))
        if run != 'standard':
            self.node_mgt.xmlSetData('run_type', run)
        else:
            node = self.node_mgt.xmlInitNode('run_type')
            node.xmlRemoveNode()

        self.case['run_type'] = run


    # Methods to manage mesh saving behavior
    #=======================================

    @Variables.noUndo
    def getMeshSaveOnModify(self):
        """
        Get save on preprocess/modify option.
        """

        val = True
        if self.node_ecs:
            s = self.node_ecs.xmlGetString('save_mesh_if_modified')
            if str(s) == 'no':
                val = False

        return val


    @Variables.undoLocal
    def setMeshSaveOnModify(self, val):
        """
        Set run type.
        """
        if val == False:
            self.node_ecs.xmlSetData('save_mesh_if_modified', 'no')
        else:
            node = self.node_ecs.xmlGetNode('save_mesh_if_modified')
            if node:
                node.xmlRemoveNode()


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
        """ Check whether the meshes could be added and deleted and list of meshes could be got """
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
        """ Check whether the status of node cut and value of angle could be set and got"""
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


    def checkgetPeriodicSelectionsCount(self):
        """ Check whether the number of periodicities could be got"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        select['visualization'] = 1
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
        """ Check whether the mode of transformation could be set and got """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        select['visualization'] = 1
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
        """ Check whether the dir values translation mode of periodicity could be set and got"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        select['visualization'] = 1
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
        """ Check whether the values for rotation's mode of periodicity could be set and got"""
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = 1
        select['visualization'] = 1
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
        select['visualization'] = 1
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
        select['visualization'] = 1
        mdl = SolutionDomainModel(self.case)
        mdl.addJoinFaces(select)
        doc = '''<face_joining name="1">
                    <selector>1 or 2 or 3 or toto</selector>
                    <fraction>0.1</fraction>
                    <plane>20</plane>
                    <verbosity>1</verbosity>
                    <visualization>1</visualization>
                 </face_joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for face joining'
        assert mdl.getJoinFaces('1') == {'selector': '1 or 2 or 3 or toto', 'plane': '20',
                                         'fraction': '0.1', 'verbosity': '1',
                                         'visualization': '1'},\
            'Could not get values of faces join for face joining'

    def checkReplaceandDeleteandSetandGetForJoinFaces(self):
        """
        Check whether faces of face joining could be replaced and deleted
        and status could be set and got
        """
        select = {}
        select['selector'] = '1 or 2 or 3 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '20'
        select['verbosity'] = '1'
        select['visualization'] = '1'
        deux = {}
        deux['selector'] = '9 or 8 or 7 or coucou'
        deux['fraction'] = '0.2'
        deux['plane'] = '20'
        deux['verbosity'] = '2'
        deux['visualization'] = '2'
        mdl = SolutionDomainModel(self.case)
        mdl.addJoinFaces(select)
        mdl.addJoinFaces(deux)
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                            <visualization>1</visualization>
                    </face_joining>
                    <face_joining name="2">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                            <visualization>2</visualization>
                    </face_joining>
                 </joining>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for face joining'
        assert mdl.getJoinFaces('1') == {'selector': '1 or 2 or 3 or toto',
                                        'plane': '20', 'fraction': '0.1',
                                        'verbosity': '1', 'visualization': '1'},\
            'Could not get values of faces join for face joining'

        select['selector'] = 'je vais partir'
        mdl.replaceJoinFaces('1', select)
        doc = '''<joining>
                    <face_joining name="1">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                            <visualization>1</visualization>
                    </face_joining>
                    <face_joining name="2">
                            <selector>9 or 8 or 7 or coucou</selector>
                            <fraction>0.2</fraction>
                            <plane>30</plane>
                            <verbosity>2</verbosity>
                            <visualization>2</visualization>
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
                            <visualization>2</visualization>
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
                            <visualization>2</visualization>
                    </face_joining>
                    <face_joining name="2">
                            <selector>1 or 2 or 3 or toto</selector>
                            <fraction>0.1</fraction>
                            <plane>20</plane>
                            <verbosity>1</verbosity>
                            <visualization>1</visualization>
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
        select['visualization'] = '2'
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
                                             'plane': '30', 'fraction': '0.1',
                                             'verbosity': '2', 'visualization': '2'},\
            'Could not get values of faces for periodicities'

    def checkReplaceandDeleteandSetandGetStatusForPeriodicFaces(self):
        """
        Check whether faces of of periodicity could be replaced and deleted
        and status could be set and got
        """
        select = {}
        select['selector'] = '5 or 6 or toto'
        select['fraction'] = '0.1'
        select['plane'] = '25'
        select['verbosity'] = '1'
        select['visualization'] = '1'
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
                       <visualization>1</visualization>
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
        select['visualization'] = '2'
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
                     <visualization>1</visualization>
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
                      <visualization>1</visualization>
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


    def checkReorientSetAndGetStatusAndCommand(self):
        """ Check whether reorient status could be set and got and command line could be got """
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
                       ("toto.csm.gz",   "csm")  ,
                       ("toto.des",      "des")  ,
                       ("toto.ccm.gz",   "ccm")  ,
                       ("toto.med",      "med")  ,
                       ("toto.msh.gz",   "msh")  ,
                       ("toto.neu",      "neu")  ,
                       ("toto.unv.gz",   "unv")  ,
                       ("toto",          "")     ,
                       ("toto.gz",       "")     ]

    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        pass

    def checkGetMeshExtension(self):
        """Check whether mesh extension could be got"""
        mdl = MeshModel()

        for f in self.files:
          ext = mdl.getMeshExtension(f[0])
          assert ext == f[1], 'could not get the mesh extension'

    def checkGetMeshFormat(self):
        """Check whether mesh extension could be got"""
        mdl = MeshModel()

        for f in self.files:
            m = (f[0], "")
            fmt = mdl.getMeshFormat(m)
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
