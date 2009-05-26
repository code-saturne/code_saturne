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
        self.ext['unv']   = "unv"
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
                ("I-DEAS universal files",    "*.unv"  ),
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
        self.node_join       = self.node_ecs.xmlInitNode('join_meshes', "status")
        self.node_cut        = self.node_ecs.xmlInitNode('faces_cutting', "status")
        self.node_orient     = self.node_ecs.xmlInitNode('reorientation', "status")
        self.node_syrthes    = self.node_ecs.xmlInitNode('syrthes_coupling', "status")
        self.node_perio      = self.node_ecs.xmlInitNode('periodic_boundary')
        self.node_standalone = self.node_ecs.xmlInitNode('standalone')
##        self.node_select     = self.node_standalone.xmlInitNode('faces_select', "status")


#************************* Private methods *****************************

    def defaultValues(self):
        """
        Return a dictionary with default values
        """
        defvalue = {}
        defvalue['join_status']    = "off"
        defvalue['cutting_status'] = "off"
        defvalue['reorientation']  = "off"
        defvalue['select_status']  = "off"
        defvalue['color']          = ""
        defvalue['group']          = ""
        defvalue['reverse']        = "off"
        defvalue['fraction']       = 0.1
        defvalue['plan']           = 0.8
        defvalue['semiconf']       = "off"
        defvalue['angle']          = 0.01
        defvalue['syrth_status']   = "off"
        defvalue['syrth_mesh_2d']  = "off"
        defvalue['sim_status']     = "on"
        defvalue['verif_mail']     = "on"
        defvalue['postprocessing_format'] = "EnSight"
        defvalue['postprocessing_options'] = "binary"
        defvalue['dir_cas']        = "cas_defaut"
        defvalue['poly_status']    = "off"
        defvalue['perio_mode']     = "translation"
        defvalue['transfo_val']    = 0.0

        return defvalue


    def _updateBatchScriptFile(self, keyword):
        """
        Update, for keyword, the backup file if it's ready to run.
        """
        self.isInList(keyword,('MESH', 
                               'COMMAND_REORIENT',
                               'COMMAND_JOIN',
                               'COMMAND_CWF',
                               'COMMAND_PERIO',
                               'COMMAND_SYRTHES',
                               'CWF_OFF',
                               'JOIN_OFF',
                               'PERIO_OFF',
                               'SYRTHES_OFF'))
        key = self.case['computer']
        if key:
            if not self.case['batchScript'][key]: return

            from BatchRunningModel import BatchRunningModel
            batch = BatchRunningModel(self.case)
            if keyword in ('CWF_OFF', 'JOIN_OFF', 'PERIO_OFF', 'SYRTHES_OFF'):
                cmd = string.split(keyword,'_')[:1][0]
                keyword = 'COMMAND_' + str(cmd)
                batch.dicoValues[keyword] = ''
            batch.initializeBatchScriptFile()
            batch.updateBatchScriptFile(keyword)
            del BatchRunningModel


#To follow : private methods to get or put faces
#================================================

    def _getTagNode(self, tagName):
        """
        Private method: Return node corresponding at item "tag"
        """
        self.isInList(tagName, ('faces_join', 'faces_syrthes', 
                                'faces_select', 'faces_periodic'))
        if tagName == 'faces_join':
            node = self.node_join
        elif tagName == 'faces_syrthes':
            node = self.node_syrthes
        elif tagName == 'faces_select':
            node = self.node_standalone
        elif tagName == 'faces_periodic':
            node = self.node_perio
        return node


    def _updateJoinSelectionsNumbers(self):
        """
        Update names of join selection
        """
        listNode = self.node_join.xmlGetNodeList('faces_join')
        i = 0
        for node in listNode:
            i = i + 1
            if int(node['name']) > i:
                node['name'] = str(i)


    def getJoinSelectionsNumber(self):
        """
        Public method.
        
        @return: number of join faces selections
        @rtype: C{int}
        """
        return len(self.node_join.xmlGetNodeList('faces_join'))


    def _addFacesSelect(self, node, select):
        """
        Private method: Add faces to node (join, periodic or syrthes) with dictionary select.
        """
        for sel, txt in [ (select['color'],    'faces_color'),
                          (select['group'],    'faces_group'),
                          (select['fraction'], 'faces_fraction'),
                          (select['plan'],     'faces_plan')]:
            if sel:
                node.xmlSetData(txt, sel)
            else:
                node.xmlRemoveChild(txt)

        if select['reverse'] == 'on':
            node_revers = node.xmlInitNode('faces_reverse', status='on')
        else:
            node.xmlRemoveChild('faces_reverse')

        if select['semiconf'] == 'on':
            node_semiconf = node.xmlInitNode('faces_semi_conf', status='on')
        else:
            node.xmlRemoveChild('faces_semi_conf')


    def _getFaces(self, node):
        """
        Private method: Return values found for color, group .. for node "node"
        """
        default = {}
        default['color'] =""
        default['group'] = ""
        default['fraction'] = ""
        default['plan'] = ""
        default['reverse'] = ""
        default['semiconf'] = ""

        if node:
            default['color']    = node.xmlGetString('faces_color')
            default['group']    = node.xmlGetString('faces_group')
            default['fraction'] = node.xmlGetString('faces_fraction')
            default['plan']     = node.xmlGetString('faces_plan')

            n_revers = node.xmlGetNode('faces_reverse', 'status')

            if n_revers:
                default['reverse'] = n_revers['status']
            else:
                default['reverse'] = "off"

            n_semi_conf = node.xmlGetNode('faces_semi_conf', 'status')

            if n_semi_conf:
                default['semiconf'] = n_semi_conf['status']
            else:
                default['semiconf'] = "off"
        if default['color'] == '' and default['group'] == '' and default['fraction'] == ""\
                                  and default['plan'] == "" and default['reverse'] == "off" \
                                  and default['semiconf'] == "off":
            default = {}

        return default


    def _removeChildren(self, node):
        """
        Private method: Remove all child nodes of node for one selection 
        """
        for tag in ('faces_color',
                    'faces_group',
                    'faces_fraction',
                    'faces_plan',
                    'faces_reverse',
                    'faces_semi_conf'):
            node.xmlRemoveChild(tag)


    def _getLineCommand(self, node):
        """
        Private method: Get color group faces and revers satus and fraction 
        and plan datas for ommand line to preprocessor execution
        """
        line = ""
        coul = node.xmlGetString('faces_color')
        grp = node.xmlGetString('faces_group')
        n_revers = node.xmlGetNode('faces_reverse', 'status')
        n_semi_conf = node.xmlGetNode('faces_semi_conf',' status')

        if coul:
            line += " --color " + coul
        if grp:
            line += " --group " + grp
        if n_revers and n_revers['status'] == "on":
            line += " --invsel"
        if n_semi_conf and n_semi_conf['status'] == "on":
            line += " --semi-conf"

        return line


#To follow : private methods for periodicity:
#===========================================


    def _setPeriodicNewMode(self, perio_name, new_mode):
        """
        Private method: Set node and mode of periodicity for transformation named 'perio_name'
        """
        #self.isInList(perio_name, ("0", "1", "2", "3"))
        #self.isLowerOrEqual(int(perio_name), int(self.getPeriodicityNumber()))
        self.isInList(perio_name, self.getPeriodicityListName())
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        node['mode'] = new_mode


    def _setTranslationDefault(self, perio_name):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        if node: 
            nList = node.xmlInitChildNodeList('translation')
            for n in nList:
                n.xmlSetData('translation_x', self.defaultValues()['transfo_val'])
                n.xmlSetData('translation_y', self.defaultValues()['transfo_val'])
                n.xmlSetData('translation_z', self.defaultValues()['transfo_val'])


    def _setRotation1Default(self, perio_name):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        nList = node.xmlInitChildNodeList('rotation1')
        for n in nList:
            n.xmlSetData('rotation_angle', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_x', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_y', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_z', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_center_x', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_center_y', self.defaultValues()['transfo_val'])
            n.xmlSetData('rotation_center_z', self.defaultValues()['transfo_val'])


    def _setRotation2Default(self, perio_name):
        """
        Private method: Put default values of translation for periodic translation
        """
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        nList = node.xmlInitChildNodeList('rotation2')
        for n in nList:
            for txt in ('rotation_matrix_11', 'rotation_matrix_12', 'rotation_matrix_13',
                        'rotation_matrix_21', 'rotation_matrix_22', 'rotation_matrix_23',
                        'rotation_matrix_31', 'rotation_matrix_32', 'rotation_matrix_33',
                        'rotation_center_y', 'rotation_center_z', ):
                n.xmlSetData(txt, self.defaultValues()['transfo_val'])


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
        self._updateBatchScriptFile('MESH')


    def delMesh(self, mesh):
        """
        Public method. Delete node for mesh named "mesh" in xml file
        """
        self.isInList(mesh, self.getMeshList())
        nodeList = self.node_meshes.xmlGetNodeList('mesh', 'name')
        for node in nodeList:
            if node['name'] == mesh: 
                node.xmlRemoveNode()
        self._updateBatchScriptFile('MESH')


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
        self._updateBatchScriptFile('MESH')


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
        self._updateBatchScriptFile('MESH')


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
        self._updateBatchScriptFile('MESH')


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

# Methods to manage status of all main balises : 
#=============================================

    def getJoinMeshesStatus(self):
        """
        Get status on balise "join_meshes" from xml file
        """
        status = self.node_join['status']
        if not status:
            status = self.defaultValues()['join_status']
            self.setJoinMeshesStatus(status)
        return status


    def setJoinMeshesStatus(self, status):
        """
        Put status on balise "join_meshes" in xml file
        """
        self.isOnOff(status)
        if self.node_join:
            self.node_join['status'] = status
        else:
            self.node_join = self.node_ecs.xmlInitNode('join_meshes', status=status)
        if status == 'off':
            self._updateBatchScriptFile('JOIN_OFF')
        else:
            self._updateBatchScriptFile('COMMAND_JOIN')


    def getCutStatus(self):
        """
        Get status on balise "faces_cutting" from xml file
        """
        status = self.node_cut['status']
        if not status:
            status = self.defaultValues()['cutting_status']
            self.setCutStatus(status)
        return status


    def setCutStatus(self, status):
        """
        Put status on balise "faces_cutting" in xml file
        """
        self.isOnOff(status)
        self.node_cut['status'] = status
        if status == 'off':
            self._updateBatchScriptFile('CWF_OFF')
        else:
            self._updateBatchScriptFile('COMMAND_CWF')


    def setCutAngle(self, var):
        """
        input '--cwf' parameter.
        """
        self.isGreaterOrEqual(var, 0.0)
        if var != self.defaultValues()['angle']:
            self.node_cut.xmlSetData('warp_angle_max', var)
        else:
            self.node_cut.xmlRemoveChild('warp_angle_max')
        self._updateBatchScriptFile('COMMAND_CWF')


    def getCutAngle(self):
        """
        get '--cwf' parameters.
        """
        angle = self.node_cut.xmlGetDouble('warp_angle_max')
        if angle == None:
            angle = self.defaultValues()['angle']
        return angle


    def getOrientation(self):
        """
        Get status on balise "reorientation" from xml file
        """
        status = self.node_orient['status']
        if not status:
            status = self.defaultValues()['reorientation']
            self.setOrientation(status)
        return status


    def setOrientation(self, status):
        """
        Put status on balise "reorientation" in xml file
        """
        self.isOnOff(status)
        self.node_orient['status'] = status
        self._updateBatchScriptFile('COMMAND_REORIENT')


    def getSyrthesCouplingStatus(self):
        """
        Get status on balise "syrthes_coupling" in xml file
        """
        status = self.node_syrthes['status']
        if not status:
            status = self.defaultValues()['syrth_status']
            self.setSyrthesCouplingStatus(status)
        return status


    def setSyrthesCouplingStatus(self, status):
        """
        Put status on balise "syrthes_coupling" in xml file
        """
        self.isOnOff(status)
        self.node_syrthes['status'] = status
        if status == 'off':
            self._updateBatchScriptFile('SYRTHES_OFF')
        else:
            self._updateBatchScriptFile('COMMAND_SYRTHES')


    def getSyrthes2dMeshStatus(self):
        """
        Return status of balise "'syrthes_mesh_2d'" from xml file
        """
        node = self.node_syrthes.xmlInitNode('syrthes_mesh_2d', 'status')
        status = node['status']
        if not status : 
            status = self.defaultValues()['syrth_mesh_2d']
            self.setSyrthes2dMeshStatus(status)
        return status


    def setSyrthes2dMeshStatus(self, status):
        """
        Put status of balise 'syrthes_mesh_2d' into xml file
        """
        self.isOnOff(status)
        node = self.node_syrthes.xmlInitNode('syrthes_mesh_2d', 'status')
        if status != self.defaultValues()['syrth_mesh_2d']:
            node['status'] = status
        else:
            if node:
                node.xmlRemoveNode()
        self._updateBatchScriptFile('COMMAND_SYRTHES')


    def getSimCommStatus(self):
        """
        Get status of balise ''similation_communication' into xml file             
        """
        node = self.node_standalone.xmlInitNode('simulation_communication', 'status')
        status = node['status']
        if not status:
            status = self.defaultValues()['sim_status']
            self.setSimCommStatus(status)
        return status


    def setSimCommStatus(self, status):
        """
        Put status of balise ''similation_communication' into xml file 
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

    def getPeriodicityListName(self):
        """
        Public method.

        @return: list of name of defined periodic transformation
        @rype: C{list}
        """
        l = []
        for node in self.node_perio.xmlGetNodeList('transformation'):
            l.append(node['name'])
        return l


    def getPeriodicityNumber(self):
        """
        Public method.

        @return: number of "periodic_boundary" markup in xml file
        @rtype: C{int}
        """
        return len(self.node_perio.xmlGetNodeList('transformation'))


    def getPeriodicityMode(self, perio_name):
        """
        Public method.

        @type perio_name: C{str}
        @param perio_name: name of the periodic boundary
        @return: mode of transformation of periodic boundary I{perio_name}
        @rtype: C{str}
        """
        self.isInList(perio_name, self.getPeriodicityListName())

        node = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        mode = node['mode']
        if not mode:
            mode = self.defaultValues()['perio_mode']
            self.addPeriodicity(perio_name)
        return mode


    def addPeriodicity(self, perio_name):
        """
        Public method.

        Add a new transformation in periodic boundary.

        @type perio_name: C{str}
        @param perio_name: name of the periodic boundary
        """
        self.isNotInList(perio_name, self.getPeriodicityListName())

        m = self.defaultValues()['perio_mode']
        self.node_perio.xmlInitNode('transformation', mode="", name=perio_name)
        self.updatePeriodicityMode(perio_name, m)


    def updatePeriodicityMode(self, perio_name, mode):
        """
        Public method.

        Update transformation mode from a periodic boundary

        @type perio_name: C{str}
        @param perio_name: name of the periodic boundary
        @type mode: C{str}
        @param mode: mode of the periodic boundary (i.e.: 'translation', 'rotation1', 'rotation2', 'tr+rota1', 'tr+rota2')
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        self.isInList(mode, ('translation', 'rotation1', 'rotation2', 'tr+rota1', 'tr+rota2'))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        if node['mode'] != mode:
            node['mode'] = mode

            if mode in ('translation', 'rotation1', 'rotation2'):
                if not node.xmlGetChildNode(mode):
                  if mode =="translation":
                      self._setTranslationDefault(perio_name)
                  elif mode =="rotation1":
                      self._setRotation1Default(perio_name)
                  elif mode =="rotation2":
                      self._setRotation2Default(perio_name)
            else:
                if mode =="tr+rota1":
                    if not node.xmlGetChildNode('translation'):
                        self._setTranslationDefault(perio_name)
                    if not node.xmlGetChildNode('rotation1'):
                        self._setRotation1Default(perio_name)
                elif mode =="tr+rota2":
                    if not node.xmlGetChildNodeList('translation'):
                        self._setTranslationDefault(perio_name)
                    if not node.xmlGetChildNodeList('rotation2'):
                        self._setRotation2Default(perio_name)
    
            self._updateBatchScriptFile('COMMAND_PERIO')


    def deletePeriodicity(self, perio_name):
        """
        Public method.

        Delete a transformation in periodic boundary.

        @type perio_name: C{str}
        @param perio_name: name of the periodic boundary
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        self.node_perio.xmlGetNode('transformation', name=perio_name).xmlRemoveNode()

        if len(self.node_perio.xmlGetNodeList('transformation')) == 0:
            self._updateBatchScriptFile('PERIO_OFF')
        else:
            self._updateBatchScriptFile('COMMAND_PERIO')


    def changePeriodicityName(self, perio_name, new_name):
        """
        Public method.

        Change the label name of a periodic boundary

        @type perio_name: C{str}
        @param perio_name: old name of the periodic boundary
        @type new_name: C{str}
        @param new_name: new name of the periodic boundary
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        node['name'] = new_name


    def getTranslationDirection(self, perio_name):
        """
        Public method.

        Get values of translation for periodic translation

        @type perio_name: C{str}
        @param perio_name: name of the periodic boundary
        @return: values of translation for periodic translation
        @rtype: 3 C{float}
        """
        self.isInList(perio_name, self.getPeriodicityListName())

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('translation')
        dx = n.xmlGetString('translation_x')
        dy = n.xmlGetString('translation_y')
        dz = n.xmlGetString('translation_z')

        return dx, dy, dz


    def setTranslationDirection(self, perio_name, dir, valcoor):
        """
        Put values of translation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(perio_name, self.getPeriodicityListName())
        self.isInList(dir, ("translation_x", "translation_y", "translation_z"))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        for n in node.xmlGetChildNodeList('translation'):
            n.xmlSetData(dir, valcoor)
        self._updateBatchScriptFile('COMMAND_PERIO')


    def getRotationDirection(self, perio_name):
        """
        Get values for director vector rotation for periodic translation
        """
        self.isInList(perio_name, self.getPeriodicityListName())

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('rotation1')
        rx = n.xmlGetString('rotation_x')
        ry = n.xmlGetString('rotation_y')
        rz = n.xmlGetString('rotation_z')

        return rx, ry, rz


    def setRotationVector(self, perio_name, dir, valcoor):
        """
        Put values for director vector rotation for periodic translation
        """
        self.isFloat(valcoor)
        self.isInList(perio_name, self.getPeriodicityListName())
        self.isInList(dir, ("rotation_x", "rotation_y", "rotation_z"))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('rotation1')
        n.xmlSetData(dir,valcoor)
        self._updateBatchScriptFile('COMMAND_PERIO')


    def getRotationAngle(self, perio_name):
        """
        Get angle for rotation for periodic rotation
        """
        self.isInList(perio_name, self.getPeriodicityListName())

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('rotation1')
        angle = n.xmlGetString('rotation_angle')

        return angle


    def setRotationAngle(self, perio_name, angle):
        """
        Put angle for rotation for periodic rotation
        """
        self.isGreaterOrEqual(angle, 0.0)
        self.isInList(perio_name, self.getPeriodicityListName())
        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('rotation1')
        n.xmlSetData('rotation_angle', angle)
        self._updateBatchScriptFile('COMMAND_PERIO')


    def getRotationCenter(self, perio_name):
        """
        Get coordinates of center of rotation for periodic transformation
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        mode = self.getPeriodicityMode(perio_name)
        self.isInList(mode, ('rotation1', 'rotation2', 'tr+rota1', 'tr+rota2'))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name) 

        if mode == "rotation1" or mode == "tr+rota1":
            n = node.xmlGetChildNode('rotation1')
        elif mode == "rotation2" or mode == "tr+rota2":
            n = node.xmlGetChildNode('rotation2')
        px = n.xmlGetString('rotation_center_x')
        py = n.xmlGetString('rotation_center_y')
        pz = n.xmlGetString('rotation_center_z')

        return px, py, pz


    def setRotationCenter(self, perio_name, pos, val):
        """
        Put coordinates of center of rotation for periodic transformation
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        self.isFloat(val)
        self.isInList(pos, ("rotation_center_x", "rotation_center_y", "rotation_center_z"))
        mode = self.getPeriodicityMode(perio_name)
        self.isInList(mode, ('rotation1', 'rotation2', 'tr+rota1', 'tr+rota2'))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)

        if mode == "rotation1" or mode == "tr+rota1":
            n = node.xmlGetChildNode('rotation1')
        elif mode == "rotation2" or mode == "tr+rota2":
            n = node.xmlGetChildNode('rotation2')
        n.xmlSetData(pos, val)
        self._updateBatchScriptFile('COMMAND_PERIO')


    def getRotationMatrix(self, perio_name):
        """
        Get values of matrix of rotation for periodic transformation
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        mode = self.getPeriodicityMode(perio_name)
        self.isInList(mode, ('rotation2', 'tr+rota2'))

        node = self.node_perio.xmlGetNode('transformation', name=perio_name)
        n = node.xmlGetChildNode('rotation2')
        m11 = n.xmlGetString('rotation_matrix_11')
        m12 = n.xmlGetString('rotation_matrix_12')
        m13 = n.xmlGetString('rotation_matrix_13')
        m21 = n.xmlGetString('rotation_matrix_21')
        m22 = n.xmlGetString('rotation_matrix_22')
        m23 = n.xmlGetString('rotation_matrix_23')
        m31 = n.xmlGetString('rotation_matrix_31')
        m32 = n.xmlGetString('rotation_matrix_32')
        m33 = n.xmlGetString('rotation_matrix_33')

        return m11, m12, m13, m21, m22, m23, m31, m32, m33


    def setRotationMatrix(self, perio_name, pos, val):
        """
        Put values of matrix of rotation for periodic transformation
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        self.isFloat(val)
        self.isInList(pos, ('rotation_matrix_11','rotation_matrix_12','rotation_matrix_13',
                            'rotation_matrix_21','rotation_matrix_22','rotation_matrix_23',
                            'rotation_matrix_31','rotation_matrix_32','rotation_matrix_33'))
        mode = self.getPeriodicityMode(perio_name)
        self.isInList(mode, ('rotation2', 'tr+rota2'))

        node = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        n = node.xmlGetChildNode('rotation2')
        n.xmlSetData(pos, val)
        self._updateBatchScriptFile('COMMAND_PERIO')


# Methods to manage faces :
#========================

    def addJoinFaces(self, select):
        """
        Add faces selection for join meshes. 
        Select is a dictionary with 'color', 'group', 'fraction', 'plan' 
        """
        nb = self.getJoinSelectionsNumber()
        name = str(nb +1)
        node = self.node_join.xmlAddChild('faces_join', status="on", name=name)
        self._addFacesSelect(node, select)
        self._updateBatchScriptFile('COMMAND_JOIN')


    def getJoinFaces(self, number):
        """
        Return faces selection named 'number' for join meshes . 
        """
        self.isLowerOrEqual(int(number), int(self.getJoinSelectionsNumber()))

        node = self.node_join.xmlGetNode('faces_join', status="on", name=number)
        return self._getFaces(node)


    def replaceJoinFaces(self, number, select):
        """
        Replace values of faces selection named 'number' for join meshes, by select
        """
        self.isLowerOrEqual(int(number), int(self.getJoinSelectionsNumber()))

        node = self.node_join.xmlGetNode('faces_join', status="on", name=number)
        self._removeChildren(node)
        self._addFacesSelect(node, select)
        self._updateBatchScriptFile('COMMAND_JOIN')


    def deleteJoinFaces(self, number):
        """
        Delete faces selection named 'number' for join meshes
        """
        self.isLowerOrEqual(int(number), int(self.getJoinSelectionsNumber()))

        node = self.node_join.xmlGetNode('faces_join', status="on", name=number)
        node.xmlRemoveNode()
        if int(number) <= int(self.getJoinSelectionsNumber()):
            self._updateJoinSelectionsNumbers()
        self._updateBatchScriptFile('COMMAND_JOIN')


    def setJoinStatus(self, number, status):
        """
        Set status of faces selection named 'number' for join meshes
        """
        self.isOnOff(status)
        self.isLowerOrEqual(int(number), int(self.getJoinSelectionsNumber()))

        node = self.node_join.xmlGetNode('faces_join', 'status', name=number)
        node['status'] = status


    def getJoinStatus(self, number):
        """
        Get status of faces selection named 'number' for join meshes
        """
        self.isLowerOrEqual(int(number), int(self.getJoinSelectionsNumber()))

        return self.node_join.xmlGetNode('faces_join', 'status', name=number)['status']


    def addPeriodicFaces(self, perio_name, select):
        """
        Add faces selection for periodic transformation. 
        Select is a dictionary with 'color', 'group', 'fraction', 'plan' ...
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlAddChild('faces_periodic', status="on")
            self._addFacesSelect(node, select)
            self._updateBatchScriptFile('COMMAND_PERIO')


    def getPeriodicFaces(self, perio_name):
        """
        Public method.
        
        @return: faces selection for periodic transformation named perio_name
        @rtype: C{dictionary}
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        result = {}
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlGetChildNode('faces_periodic', 'status')
            if node and node['status'] == 'on':
                result = self._getFaces(node)

        return result

        
    def replacePeriodicFaces(self, perio_name, select):
        """
        Replace values of faces selection for periodic transformation, by select
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlGetChildNode('faces_periodic', status="on")
            self._removeChildren(node)
            self._addFacesSelect(node, select)
            self._updateBatchScriptFile('COMMAND_PERIO')


    def deletePeriodicFaces(self, perio_name):
        """
        Delete faces selection for periodic transformation named 'perio_name'
        """
        self.isInList(perio_name, self.getPeriodicityListName())
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlGetChildNode('faces_periodic', status="on")
            node.xmlRemoveNode()
            self._updateBatchScriptFile('COMMAND_PERIO')


    def setPeriodicStatus(self, perio_name, status):
        """
        Set status of faces for periodic transformation named 'perio_name'
        """
        self.isOnOff(status)
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlGetChildNode('faces_periodic', 'status')
            if node:
                node['status'] = status
            if status == 'on':
                self._updateBatchScriptFile('COMMAND_PERIO')


    def getPeriodicStatus(self, perio_name):
        """
        Get status of faces for periodic transformation named 'perio_name'
        """
        node_tr = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
        if node_tr:
            node = node_tr.xmlGetChildNode('faces_periodic', 'status')
            if node:
                status = node['status']
            return status
        else:
            raise ValueError, "wrong periodicity"


    def addSyrthesFaces(self, select):
        """
        Add faces selection for syrthes coupling. 
        Select is a dictionary with 'color', 'group', 'fraction', 'plan' ...
        """
        node = self.node_syrthes.xmlAddChild('faces_syrthes', status="on")
        select['fraction'] =""
        select['plan'] = ""
        self._addFacesSelect(node, select)
        self._updateBatchScriptFile('COMMAND_SYRTHES')
        
        
    def getSyrthesFaces(self):
        """
        Return faces selection for syrthes coupling (only one authorized selection)
        """
        result = {}
        node  = self.node_syrthes.xmlGetChildNode('faces_syrthes', status="on")
        if node:
            result = self._getFaces(node)
        
        return result


    def replaceSyrthesFaces(self, select):
        """
        Replace values of faces selection for syrthes coupling (only one authorized selection)
        by select
        """
        node  = self.node_syrthes.xmlGetChildNode('faces_syrthes', status="on")
        self._removeChildren(node)
        self._addFacesSelect(node, select)
        self._updateBatchScriptFile('COMMAND_SYRTHES')


    def deleteSyrthesFaces(self):
        """
        Delete faces selection for syrthes coupling (only one authorized selection)
        """
        node = self.node_syrthes.xmlGetChildNode('faces_syrthes', status="on")
        node.xmlRemoveNode()
        self._updateBatchScriptFile('COMMAND_SYRTHES')


    def setSyrthesStatus(self, status):
        """
        Set status of faces selection for syrthes coupling (only one authorized selection)
        """
        self.isOnOff(status)
        node = self.node_syrthes.xmlGetChildNode('faces_syrthes', 'status')
        if node:
            node['status'] = status


    def getSyrthesStatus(self):
        """
        Get status of faces selection for syrthes coupling (only one authorized selection)
        """
        node = self.node_syrthes.xmlGetChildNode('faces_syrthes', 'status')
        if node:
            status = node['status']
        return status


    def addSelectFaces(self, select):
        """
        Add faces selection for standalone selection. 
        Select is a dictionary with 'color', 'group', 'fraction', 'plan' ...
        """
        node = self.node_standalone.xmlGetChildNode('faces_select', 'status')
        if node:
            node['status'] = "on"
        else:
            node = self.node_standalone.xmlAddChild('faces_select', status="on")
        select['fraction'] =""
        select['plan'] = ""
        select['semiconf'] = ""
        self._addFacesSelect(node, select)


    def getSelectFaces(self):
        """
        Return faces selection for standalone selection (only one authorized selection)
        """
        result = {}
        node  = self.node_standalone.xmlGetChildNode('faces_select', status="on")
        if node:
            result = self._getFaces(node)

        return result


    def replaceSelectFaces(self, select):
        """
        Replace values of faces selection for standalone selection (only one authorized selection)
        by select
        """
        node  = self.node_standalone.xmlGetChildNode('faces_select', status="on")
        self._removeChildren(node)
        select['fraction'] =""
        select['plan'] = ""
        select['semiconf'] = ""
        self._addFacesSelect(node, select)


    def deleteSelectFaces(self):
        """
        Delete faces selection for standalone selection (only one authorized selection)
        """
        node = self.node_standalone.xmlGetChildNode('faces_select', 'status')
        node.xmlRemoveNode()


    def setSelectStatus(self, status):
        """
        Set status of faces selection for standalone selection (only one authorized selection)
        """
        self.isOnOff(status)
        node = self.node_standalone.xmlInitNode('faces_select', 'status')
        if status == "off":
            self.deleteSelectFaces()
        else:
            node['status'] = status


    def getSelectStatus(self):
        """
        Get status of faces selection for standalone selection (only one authorized selection)
        """
        node = self.node_standalone.xmlInitChildNode('faces_select', 'status')
        status = node['status']
        if not status:
            status = self.defaultValues()['select_status']
            self.deleteSelectFaces()

        return status


# In following methods we build command for "lance" file
#=======================================================

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
        Get reorient command line for preprocessor execution
        """
        line = ''
        if self.node_orient and self.node_orient['status'] == 'on':
            line += ' --reorient '
        
        return line


    def getJoinCommand(self):
        """
        Get rc command line for preprocessor execution
        """
        lines = ''
        if self.node_join and self.node_join['status'] == 'on':
            node_face_join_list = self.node_join.xmlGetNodeList('faces_join')
            if node_face_join_list:
                for node_face_join in node_face_join_list:
                    if node_face_join['status'] == 'on': 
                        linecoul = ' -j '
                        line = self._getLineCommand(node_face_join)
                        lines = lines + linecoul + line
                        
                        fraction = node_face_join.xmlGetString('faces_fraction')
                        lines = lines + " --fraction " + fraction +" "
                        plan = node_face_join.xmlGetString('faces_plan')
                        lines = lines + " --plane " + plan +" "
                        lines = lines + ' '
            else:
                lines = ' -j '

        return lines


    def getCutCommand(self):
        """
        Get cwf command line for preprocessor execution
        """
        line = ''
        if self.node_cut and self.node_cut['status'] == 'on':
            line += ' --cwf ' + str(self.getCutAngle())

        return line


    def getPerioCommand(self):
        """
        Get perio command line for preprocessor execution
        """
        line   = ''
        for perio_name in self.getPeriodicityListName():
            node = self.node_perio.xmlGetNode('transformation', 'mode', name=perio_name)
            mode = self.getPeriodicityMode(perio_name)
            node_face_perio = node.xmlGetNode('faces_periodic', 'status')
            if node_face_perio:
                node_face_perio_list = node.xmlGetNodeList('faces_periodic', 'status')
                if node_face_perio_list:
                    for node in node_face_perio_list:
                        if node['status'] == 'on':
                            lineperio = ' --perio '
                            l = self._getLineCommand(node)
                            line = line + lineperio + l
                            fraction = node.xmlGetString('faces_fraction')
                            if fraction: 
                                line = line + " --fraction " + fraction +" "
                            plan = node.xmlGetString('faces_plan')
                            if plan: 
                                line = line + " --plane " + plan +" "
            else: 
                line = line + ' --perio '
            mode = self.getPeriodicityMode(perio_name) 
            if mode == 'translation' or mode == 'tr+rota1' or mode == 'tr+rota2':
                dx, dy, dz = self.getTranslationDirection(perio_name)
                line = line + " --trans " + dx +" " + dy +" " + dz +" "
            if mode == 'rotation1' or mode == 'tr+rota1':
                angle = self.getRotationAngle(perio_name)
                rx, ry, rz = self.getRotationDirection(perio_name)
                px, py, pz = self.getRotationCenter(perio_name)
                line = line + " --rota "  + " --angle " + angle + " --dir " + rx +" " + ry +" " + rz +" " + " --invpt " + px +" " + py +" " + pz +" "
            if mode == 'rotation2' or mode == 'tr+rota2':
                m11, m12, m13, m21, m22, m23, m31, m32, m33 = self.getRotationMatrix(perio_name)
                px, py, pz = self.getRotationCenter(perio_name)
                line = line + " --rota "  + " --matrix " + m11 +" " + m12 +" " + m13 +" " + m21 +" " + m22 +" " + m23 +" " + m31 +" " + m32 +" " + m33 +" " + " --invpt " + px +" " + py +" " + pz +" "
    
        return line


    def getSyrthesCommand(self):
        """
        Get syrthes command line for preprocessor execution
        """
        lines = ''
        if self.node_syrthes and self.node_syrthes['status'] == 'on':
            linesyr = ' --syrthes '
            node_list = self.node_syrthes.xmlGetNodeList('faces_syrthes')
            for node in node_list:
                if node['status'] == 'on':
                    line = self._getLineCommand(node)
                    linesyr = linesyr + line
            node_syrth_2d = self.node_syrthes.xmlGetNode('syrthes_mesh_2d', 'status')
            if node_syrth_2d and node_syrth_2d['status'] == 'on':
                linesyr += " --2d "
            lines += linesyr 

        return lines


    def getSimCommCommand(self):
        """
        Get " --sim-comm " command line for preprocessor execution
        """
        lines = ''
        node = self.node_standalone.xmlGetNode('simulation_communication')
        if node and node['status'] == 'on':
            lines = " -sc "
        return lines


    def getSelectCommand(self):
        """
        Get "--int-face" command line for preprocessor execution
        """
        lines  = ''
        node_list = self.node_standalone.xmlGetNodeList('faces_select')
        for nod in node_list:
            if nod['status'] == 'on':
                line = ' --int-face ' + self._getLineCommand(nod)
                lines += line
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


    def getFacesSelect(self, tagName, n):
        """
        Return default filled with values found for color, group ..
        """
        nodeList = self.getListNodes(tagName)
        node = nodeList[n]
        return self._getFaces(node)


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


    def checkSetandGetJoinMeshesStatus(self):
        """ Check whether the status of join meshes could be set and get"""
        mdl = SolutionDomainModel(self.case)
        mdl.setJoinMeshesStatus('on')
        doc = '''<join_meshes status="on"/>'''
        assert mdl.node_join == self.xmlNodeFromString(doc), \
            'Could not set status in join meshes balise'
        assert mdl.getJoinMeshesStatus() == 'on', \
            'Could not get status from join meshes balise'

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

    def checkSetandGetSyrthesCouplingStatusAndMesh2d(self):
        """ Check whether the status of node syrthes_coupling and mesh 2D could be set and get"""
        mdl = SolutionDomainModel(self.case)       
        mdl.setSyrthesCouplingStatus('on')
        doc1 = '''<syrthes_coupling status="on"/>'''
        assert mdl.node_syrthes == self.xmlNodeFromString(doc1), \
            'Could not set status for syrthes_coupling balise'
        assert mdl.getSyrthesCouplingStatus() == 'on',\
            'Could not get status for syrthes_coupling balise'
            
        mdl.setSyrthes2dMeshStatus('on')
        doc2 = '''<syrthes_coupling status="on">
                    <syrthes_mesh_2d status="on"/>
                  </syrthes_coupling>'''
        assert mdl.node_syrthes == self.xmlNodeFromString(doc2), \
            'Could not set status on syrthes_2dMesh balise'
        assert mdl.getSyrthes2dMeshStatus() == 'on',\
            'Could not get status of syrthes_mesh_2d balise'
            
        mdl.setSyrthes2dMeshStatus('off')
        doc3 = '''<syrthes_coupling status="on"/>'''
        assert mdl.node_syrthes == self.xmlNodeFromString(doc3), \
            'Could not set status OFF on syrthes_2dMesh balise'

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

    def checkGetPeriodicityNumber(self):
        """ Check whether the number of periodicities could be get"""
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.addPeriodicity('2')
        doc ='''<periodic_boundary>
                    <transformation mode="translation" name="1">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                    </transformation>
                    <transformation mode="translation" name="2">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                    </transformation>
                 </periodic_boundary>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set number of periodicities'
        assert mdl.getPeriodicityNumber() == 2,\
            'Could not get number for periodicities'

    def checkSetandgetPeriodicityMode(self):
        """ Check whether the mode of transformation could be set and get """         
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.addPeriodicity('2')
        mdl.updatePeriodicityMode('2', "tr+rota1")
        doc ='''<periodic_boundary>
                    <transformation mode="translation" name="1">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                    </transformation>
                    <transformation mode="tr+rota1" name="2">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                        <rotation1>
                            <rotation_angle>0</rotation_angle>
                            <rotation_x>0</rotation_x>
                            <rotation_y>0</rotation_y>
                            <rotation_z>0</rotation_z>
                            <rotation_center_x>0</rotation_center_x>
                            <rotation_center_y>0</rotation_center_y>
                            <rotation_center_z>0</rotation_center_z>
                        </rotation1>
                    </transformation>
            </periodic_boundary>'''

                
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set mode of transformation for periodicities'
        assert mdl.getPeriodicityMode('2') == "tr+rota1",\
            'Could not get mode of transformation for periodicities'

    def checkSetandgetTranslationDirection(self):
        """ Check whether the dir values translation mode of periodicity could be set and get"""
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.setTranslationDirection('1','translation_y',3.0)
        doc ='''<periodic_boundary>
                    <transformation mode="translation" name="1">
                            <translation>
                                    <translation_x>0</translation_x>
                                    <translation_y>3</translation_y>
                                    <translation_z>0</translation_z>
                            </translation>
                    </transformation>
                </periodic_boundary>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set one direction values for translation'
        assert mdl.getTranslationDirection('1') == ('0', '3', '0'),\
            'Could not get one direction values for translation'

    def checkSetandgetRotationDirectionandAngleandCenter(self):
        """ Check whether the values for rotation's mode of periodicity could be set and get"""
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.setTranslationDirection('1','translation_y', 3.0)
        mdl.addPeriodicity('2')
        mdl.updatePeriodicityMode('2', "tr+rota1")
        mdl.setRotationAngle('2', 180.)
        mdl.setRotationVector('2', 'rotation_x', 0.5)
        mdl.setRotationVector('2', 'rotation_z', 2.5)
        mdl.setRotationCenter('2', 'rotation_center_y', 9.8)
        doc ='''<periodic_boundary>
                    <transformation mode="translation" name="1">
                            <translation>
                                    <translation_x>0</translation_x>
                                    <translation_y>3</translation_y>
                                    <translation_z>0</translation_z>
                            </translation>
                    </transformation>
                    <transformation mode="tr+rota1" name="2">
                            <translation>
                                    <translation_x>0</translation_x>
                                    <translation_y>0</translation_y>
                                    <translation_z>0</translation_z>
                            </translation>
                            <rotation1>
                                    <rotation_angle>180</rotation_angle>
                                    <rotation_x>0.5</rotation_x>
                                    <rotation_y>0.0</rotation_y>
                                    <rotation_z>2.5</rotation_z>
                                    <rotation_center_x>0</rotation_center_x>
                                    <rotation_center_y>9.8</rotation_center_y>
                                    <rotation_center_z>0</rotation_center_z>
                            </rotation1>
                    </transformation>
                </periodic_boundary>'''

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set values for rotation transformation mode'
        assert mdl.getRotationAngle('2') == '180',\
            'Could not get value of angle for rotation transformation mode'
        assert mdl.getRotationDirection('2') == ('0.5', '0', '2.5'),\
            'Could not get values of direction for rotation transformation mode'
        assert mdl.getRotationCenter('2') == ('0', '9.8', '0'),\
            'Could not get value of center of rotation for rotation transformation mode'

    def checkSetandgetRotationMatrix(self):
        """ Check whether the matrix of rotation for rotation2 mode could be set """
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.updatePeriodicityMode('1','rotation2')
        mdl.setRotationMatrix('1', 'rotation_matrix_31', 31.31)
        doc = '''<periodic_boundary>
                    <transformation mode="rotation2" name="1">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <rotation2>
                                    <rotation_matrix_11>0.0</rotation_matrix_11>
                                    <rotation_matrix_12>0.0</rotation_matrix_12>
                                    <rotation_matrix_13>0.0</rotation_matrix_13>
                                    <rotation_matrix_21>0.0</rotation_matrix_21>
                                    <rotation_matrix_22>0.0</rotation_matrix_22>
                                    <rotation_matrix_23>0.0</rotation_matrix_23>
                                    <rotation_matrix_31>31.31</rotation_matrix_31>
                                    <rotation_matrix_32>0.0</rotation_matrix_32>
                                    <rotation_matrix_33>0.0</rotation_matrix_33>
                                    <rotation_center_y>0.0</rotation_center_y>
                                    <rotation_center_z>0.0</rotation_center_z>
                            </rotation2>
                    </transformation>
                 </periodic_boundary>'''                                  

        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set values for matrix of rotation for rotation2 transformation mode'
        assert mdl.getRotationMatrix('1') == ('0', '0', '0', '0', '0', '0', '31.31','0', '0'),\
            'Could not get values for matrix of rotation for rotation2 transformation mode'

    def checkAddandGetJoinFaces(self):
        """ Check whether faces of join meshes could be added and get """
        select = {}
        select['color'] = '1 2 3'
        select['group'] = 'toto'
        select['fraction'] = '0.1'
        select['plan'] = '0.8'
        select['reverse'] = 'off'
        select['semiconf'] = 'on'
        mdl = SolutionDomainModel(self.case)
        mdl.setJoinMeshesStatus('on')
        mdl.addJoinFaces(select)
        doc = '''<join_meshes status="on">
                    <faces_join name="1" status="on">
                            <faces_color>1 2 3</faces_color>
                            <faces_group>toto</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                            <faces_semi_conf status="on"/>
                    </faces_join>
                 </join_meshes>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for join meshes'
        assert mdl.getJoinFaces('1') == {'semiconf': 'on', 'reverse': 'off', 
                                         'color': '1 2 3', 'plan': '0.8', 
                                         'fraction': '0.1', 'group': 'toto'},\
            'Could not get values of faces join for join meshes'

    def checkReplaceandDeleteandSetandGetStatusForJoinFaces(self):
        """ 
        Check whether faces of join meshes could be replaced and deleted 
        and status could be set and get
        """
        select = {}
        select['color'] = '1 2 3'
        select['group'] = 'toto'
        select['fraction'] = '0.1'
        select['plan'] = '0.8'
        select['reverse'] = 'off'
        select['semiconf'] = 'on'
        deux = {}
        deux['color'] = '9 8 7'
        deux['group'] = 'coucou'
        deux['fraction'] = '0.2'
        deux['plan'] = '0.82'
        deux['reverse'] = 'off'
        deux['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.setJoinMeshesStatus('on')
        mdl.addJoinFaces(select)
        mdl.addJoinFaces(deux)
        doc = '''<join_meshes status="on">
                    <faces_join name="1" status="on">
                            <faces_color>1 2 3</faces_color>
                            <faces_group>toto</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                            <faces_semi_conf status="on"/>
                    </faces_join>
                    <faces_join name="2" status="on">
                            <faces_color>9 8 7</faces_color>
                            <faces_group>coucou</faces_group>
                            <faces_fraction>0.2</faces_fraction>
                            <faces_plan>0.82</faces_plan>
                    </faces_join>
                 </join_meshes>'''

        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set values of faces join for join meshes'
        assert mdl.getJoinFaces('1') == {'group': 'toto', 'reverse': 'off', 'color': '1 2 3',
                                        'plan': '0.8', 'fraction': '0.1', 
                                        'semiconf': 'on'},\
            'Could not get values of faces join for join meshes'
        
        select['group'] = 'je vais partir'
        mdl.replaceJoinFaces('1', select)
        doc = '''<join_meshes status="on">
                    <faces_join name="1" status="on">
                            <faces_color>1 2 3</faces_color>
                            <faces_group>je vais partir</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                            <faces_semi_conf status="on"/>
                    </faces_join>
                    <faces_join name="2" status="on">
                            <faces_color>9 8 7</faces_color>
                            <faces_group>coucou</faces_group>
                            <faces_fraction>0.2</faces_fraction>
                            <faces_plan>0.82</faces_plan>
                    </faces_join>
                 </join_meshes>'''
                 
        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not replace values of faces join for join meshes'
            
        mdl.deleteJoinFaces('1')
        doc = '''<join_meshes status="on">
                    <faces_join name="1" status="on">
                            <faces_color>9 8 7</faces_color>
                            <faces_group>coucou</faces_group>
                            <faces_fraction>0.2</faces_fraction>
                            <faces_plan>0.82</faces_plan>
                    </faces_join>
                 </join_meshes>'''
                 
        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not delete faces join for join meshes'
            
        mdl.addJoinFaces(select)
        mdl.setJoinStatus('1', 'off')
        doc = '''<join_meshes status="on">
                    <faces_join name="1" status="off">
                            <faces_color>9 8 7</faces_color>
                            <faces_group>coucou</faces_group>
                            <faces_fraction>0.2</faces_fraction>
                            <faces_plan>0.82</faces_plan>
                    </faces_join>
                    <faces_join name="2" status="on">
                            <faces_color>1 2 3</faces_color>
                            <faces_group>je vais partir</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                            <faces_semi_conf status="on"/>
                    </faces_join>
                 </join_meshes>'''
                 
        assert mdl.node_join == self.xmlNodeFromString(doc),\
            'Could not set status for active or not faces join for join meshes'
        assert mdl.getJoinStatus('1') == 'off',\
            'Could not get status for active or not faces join for join meshes'
        
    def checkAddandGetPeriodicFaces(self):
        """ Check whether faces of periodicity could be added and get """
        select = {}
        select['color'] = '5 6'
        select['group'] = 'toto'
        select['fraction'] = '0.1'
        select['plan'] = '0.8'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.addPeriodicity('2')
        mdl.updatePeriodicityMode('2', 'rotation1')
        mdl.addPeriodicFaces('1', select)
        mdl.addPeriodicFaces('2', select)
        doc = '''<periodic_boundary>
                    <transformation mode="translation" name="1">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                        <faces_periodic status="on">
                            <faces_color>5 6</faces_color>
                            <faces_group>toto</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                        </faces_periodic>
                    </transformation>
                    <transformation mode="rotation1" name="2">
                        <translation>
                            <translation_x>0</translation_x>
                            <translation_y>0</translation_y>
                            <translation_z>0</translation_z>
                        </translation>
                        <rotation1>
                            <rotation_angle>0</rotation_angle>
                            <rotation_x>0</rotation_x>
                            <rotation_y>0</rotation_y>
                            <rotation_z>0</rotation_z>
                            <rotation_center_x>0</rotation_center_x>
                            <rotation_center_y>0</rotation_center_y>
                            <rotation_center_z>0</rotation_center_z>
                        </rotation1>
                        <faces_periodic status="on">
                            <faces_color>5 6</faces_color>
                            <faces_group>toto</faces_group>
                            <faces_fraction>0.1</faces_fraction>
                            <faces_plan>0.8</faces_plan>
                        </faces_periodic>
                    </transformation>
                 </periodic_boundary>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not add values of faces for periodicities'
        assert mdl.getPeriodicFaces('1') == {'group': 'toto', 'reverse': 'off', 'color': '5 6',
                                        'plan': '0.8', 'fraction': '0.1', 'semiconf': 'off'},\
            'Could not get values of faces for periodicities'

    def checkReplaceandDeleteandSetandGetStatusForPeriodicFaces(self):
        """ 
        Check whether faces of of periodicity could be replaced and deleted 
        and status could be set and get
        """
        select = {}
        select['color'] = '5 6'
        select['group'] = 'toto'
        select['fraction'] = '0.1'
        select['plan'] = '0.8'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.addPeriodicity('2')
        mdl.updatePeriodicityMode('2', 'rotation1')
        mdl.addPeriodicFaces('1', select)
        mdl.addPeriodicFaces('2', select)
        mdl.deletePeriodicFaces('2')
        doc = '''<periodic_boundary>
                    <transformation mode="translation" name="1">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <faces_periodic status="on">
                                    <faces_color>5 6</faces_color>
                                    <faces_group>toto</faces_group>
                                    <faces_fraction>0.1</faces_fraction>
                                    <faces_plan>0.8</faces_plan>
                            </faces_periodic>
                    </transformation>
                    <transformation mode="rotation1" name="2">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <rotation1>
                                    <rotation_angle>0.0</rotation_angle>
                                    <rotation_x>0.0</rotation_x>
                                    <rotation_y>0.0</rotation_y>
                                    <rotation_z>0.0</rotation_z>
                                    <rotation_center_x>0.0</rotation_center_x>
                                    <rotation_center_y>0.0</rotation_center_y>
                                    <rotation_center_z>0.0</rotation_center_z>
                            </rotation1>
                    </transformation>
                 </periodic_boundary>'''
                 
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not delete one selection of faces for periodicities'
            
        select['color'] = '147 963'
        select['group'] = 'PERIODIC'
        select['fraction'] = '0.1'
        select['plan']  = '0.77'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl.replacePeriodicFaces('1', select)
        doc = '''<periodic_boundary>
                    <transformation mode="translation" name="1">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <faces_periodic status="on">
                                    <faces_color>147 963</faces_color>
                                    <faces_group>PERIODIC</faces_group>
                                    <faces_fraction>0.1</faces_fraction>
                                    <faces_plan>0.77</faces_plan>
                            </faces_periodic>
                    </transformation>
                    <transformation mode="rotation1" name="2">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <rotation1>
                                    <rotation_angle>0.0</rotation_angle>
                                    <rotation_x>0.0</rotation_x>
                                    <rotation_y>0.0</rotation_y>
                                    <rotation_z>0.0</rotation_z>
                                    <rotation_center_x>0.0</rotation_center_x>
                                    <rotation_center_y>0.0</rotation_center_y>
                                    <rotation_center_z>0.0</rotation_center_z>
                            </rotation1>
                    </transformation>
                 </periodic_boundary>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not replace values of faces for periodicities'

        mdl.setPeriodicStatus('1', 'off')
        doc = '''<periodic_boundary>
                    <transformation mode="translation" name="1">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <faces_periodic status="off">
                                    <faces_color>147 963</faces_color>
                                    <faces_group>PERIODIC</faces_group>
                                    <faces_fraction>0.1</faces_fraction>
                                    <faces_plan>0.77</faces_plan>
                            </faces_periodic>
                    </transformation>
                    <transformation mode="rotation1" name="2">
                            <translation>
                                    <translation_x>0.0</translation_x>
                                    <translation_y>0.0</translation_y>
                                    <translation_z>0.0</translation_z>
                            </translation>
                            <rotation1>
                                    <rotation_angle>0.0</rotation_angle>
                                    <rotation_x>0.0</rotation_x>
                                    <rotation_y>0.0</rotation_y>
                                    <rotation_z>0.0</rotation_z>
                                    <rotation_center_x>0.0</rotation_center_x>
                                    <rotation_center_y>0.0</rotation_center_y>
                                    <rotation_center_z>0.0</rotation_center_z>
                            </rotation1>
                    </transformation>
                 </periodic_boundary>'''
        assert mdl.node_perio == self.xmlNodeFromString(doc),\
            'Could not set status for activate or not selection of faces for periodicities'
        assert mdl.getPeriodicStatus('1') == 'off',\
            'Could not get status for activate or not selection of faces for periodicities'

    def checkAddandGetSyrthesFaces(self):
        """ Check whether faces of syrthes coupling could be added and get """
        select = {}
        select['color'] = '85 2'
        select['group'] = 'SYRT'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.setSyrthesCouplingStatus('on')
        mdl.addSyrthesFaces(select)
        doc = '''<syrthes_coupling status="on">
                    <faces_syrthes status="on">
                            <faces_color>85 2</faces_color>
                            <faces_group>SYRT</faces_group>
                    </faces_syrthes>
                 </syrthes_coupling>'''
        
        assert mdl.node_syrthes == self.xmlNodeFromString(doc),\
            'Could not add values of faces for syrthes coupling'
        assert mdl.getSyrthesFaces() == {'group': 'SYRT', 'reverse': 'off', 'color': '85 2',
                                        'plan': '', 'fraction': '', 'semiconf': 'off'},\
            'Could not get values of faces for syrthes coupling'

    def checkReplaceandDeleteandSetandGetStatusForSyrthesFaces(self):
        """ 
        Check whether faces of syrthes coupling could be replaced and deleted 
        and status could be set and get
        """
        select = {}
        select['color'] = '85 2'
        select['group'] = 'SYRT'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.setSyrthesCouplingStatus('on')
        mdl.addSyrthesFaces(select)
        select['color'] = '6 3'
        select['group'] = 'COUPLING'
        mdl.replaceSyrthesFaces(select)
        doc = '''<syrthes_coupling status="on">
                    <faces_syrthes status="on">
                            <faces_color>6 3</faces_color>
                            <faces_group>COUPLING</faces_group>
                    </faces_syrthes>
                 </syrthes_coupling>'''
                 
        assert mdl.node_syrthes == self.xmlNodeFromString(doc),\
            'Could not replace values of faces for syrthes coupling'
            
        mdl.deleteSyrthesFaces()
        doc = '''<syrthes_coupling status="on"/>'''
        
        assert mdl.node_syrthes == self.xmlNodeFromString(doc),\
            'Could not delete values of faces for syrthes coupling'
        
        select['group'] = 'NOUVEAU'
        mdl.addSyrthesFaces(select)
        mdl.setSyrthesStatus('off')
        doc = '''<syrthes_coupling status="on">
                    <faces_syrthes status="off">
                            <faces_color>6 3</faces_color>
                            <faces_group>NOUVEAU</faces_group>
                    </faces_syrthes>
                 </syrthes_coupling>'''
        
        assert mdl.node_syrthes == self.xmlNodeFromString(doc),\
            'Could not set status for activate selection of faces for syrthes coupling'
        assert mdl.getSyrthesStatus() == 'off',\
            'Could not get status for activate selection of faces for syrthes coupling'

    def checkAddandGetSelectFaces(self):
        """ Check whether faces of standalone could be added and get """
        select = {}
        select['color'] = '8 2'
        select['group'] = 'STAND'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.addSelectFaces(select)
        doc = '''<standalone>
                    <faces_select status="on">
                            <faces_color>8 2</faces_color>
                            <faces_group>STAND</faces_group>
                    </faces_select>
                 </standalone>'''
        
        assert mdl.node_standalone == self.xmlNodeFromString(doc),\
            'Could not add values of faces for standalone selection'
        assert mdl.getSelectFaces() == {'group': 'STAND', 'reverse': 'off', 'color': '8 2',
                                        'plan': '', 'fraction': '', 'semiconf': 'off'},\
            'Could not get values of faces for standalone selection'

    def checkReplaceandDeleteandSetandGetStatusForSelectFaces(self):
        """ 
        Check whether faces of standalone could be replaced and deleted 
        and status could be set and get
        """
        select = {}
        select['color'] = '8 2'
        select['group'] = 'STAND'
        select['reverse'] = 'off'
        select['semiconf'] = 'off'
        mdl = SolutionDomainModel(self.case)
        mdl.addSelectFaces(select)
        new = {}
        new['color'] = '7 8 9'
        new['group'] = 'NEW'
        new['reverse'] = 'off'
        new['semiconf'] = 'off'
        mdl.replaceSelectFaces(new)
        doc = '''<standalone>
                    <faces_select status="on">
                        <faces_color>7 8 9</faces_color>
                        <faces_group>NEW</faces_group>
                    </faces_select>
                 </standalone>'''

        assert mdl.node_standalone == self.xmlNodeFromString(doc),\
            'Could not replace values of faces for standalone selection'
        
        mdl.deleteSelectFaces()
        doc = '''<standalone/>'''
        
        assert mdl.node_standalone == self.xmlNodeFromString(doc),\
            'Could not delete values of faces for standalone selection'
            
        select['group'] = 'NOUVEAU'
        mdl.addSelectFaces(select)
        mdl.setSelectStatus('off')
        doc = '''<standalone/>'''

        assert mdl.node_standalone == self.xmlNodeFromString(doc),\
            'Could not set status for activate selection of faces for standalone selection'
        assert mdl.getSelectStatus() == 'off',\
            'Could not get status for activate selection of faces for standalone selection'

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


    def checkJoinAndCutCommand(self):
        """ Check whether join and cut command lines could be get """
        mdl = SolutionDomainModel(self.case)
        select = {}
        select['color'] = '1 2 3'
        select['group'] = 'toto'
        select['fraction'] = '0.1'
        select['plan'] = '0.8'
        select['reverse'] = 'off'
        select['semiconf'] = 'on'
        mdl.setJoinMeshesStatus('on')
        mdl.addJoinFaces(select)
        cmd_join = ' -j  --color 1 2 3 --group toto --fraction 0.1  --plane 0.8  '

        assert mdl.getJoinCommand() == cmd_join,\
            'Join command is not verified in SolutionDomain Model'
        
        mdl.setCutStatus('on')
        mdl.setCutAngle(0.05)        
        cmd_cut = mdl.getCutCommand()
        cut = ' --cwf 0.05'
        assert mdl.getCutCommand() == cmd_cut,\
            'Cut command is not verified in SolutionDomain Model'

    def checkPerioCommand(self):
        """ Check whether perio command line could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.addPeriodicity('1')
        mdl.updatePeriodicityMode('1','rotation1')
        mdl.setRotationVector('1', 'rotation_x', 9.)
        mdl.setRotationVector('1', 'rotation_y', 8.)
        mdl.setRotationVector('1', 'rotation_z', 7.)
        mdl.setRotationAngle('1', 45.)
        mdl.setRotationCenter('1', 'rotation_center_y', 66.)
        mdl.addPeriodicity('2')
        mdl.updatePeriodicityMode('2','tr+rota2')
        mdl.setTranslationDirection('2','translation_y', 3)
        mdl.setRotationMatrix('2', 'rotation_matrix_31', 31.31)
        cmd_perio = " --perio  --rota  --angle 45 --dir 9 8 7  --invpt 0 66 0  --perio  --trans 0 3 0  --rota  --matrix 0 0 0 0 0 0 31.31 0 0  --invpt  0 0 "

        assert mdl.getPerioCommand() == cmd_perio,\
            'Perio command is not verified in SolutionDomain Model'

    def checkSyrthesCommand(self):
        """ Check whether syrthes command line could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.setSyrthesCouplingStatus('on')
        mdl.setSyrthes2dMeshStatus('on')
        select = {}
        select['color'] = '1 2 3'
        select['group'] = 'toto'
        select['reverse'] = 'off'
        select['semiconf'] = 'on'
        mdl = SolutionDomainModel(self.case)
        mdl.addSyrthesFaces(select)
        cmd_syr = ''' --syrthes  --color 1 2 3 --group toto --2d '''

        assert mdl.getSyrthesCommand() == cmd_syr,\
            'Syrthes command is not verified in SolutionDomain Model'

    def checkSimCommAndVerifMaillCommand(self):
        """ Check whether simulation_communication command line could be get """
        mdl = SolutionDomainModel(self.case)
        mdl.setSimCommStatus('on')
        cmd_sim = mdl.getSimCommCommand()
        sim =' -sc '
        assert mdl.getSimCommCommand() == sim,\
            'Simulation_communication command is not verified in SolutionDomain Model'
            
        #A reprendre :
####    def checkSelectCommand(self):
####        """Check whether standalone selection command line could be get"""
####        select = {}
####        select['color'] = '1 2 3'
####        select['group'] = 'toto'
####        select['reverse'] = 'off'
####        select['semiconf'] = 'on'
####        mdl = SolutionDomainModel(self.case)
####        mdl.addSelectFaces(select)
####        print mdl.getSelectCommand()
######        --int-face  --color 1 2 3 --group toto
####        select =' --color 1 2 3  --group toto  --semi-conf '
######        select =' --int-face  --color 1 2 3  --group toto  --semi-conf '
####        assert mdl.getSelectCommand().split() == select.split(),\
####            'Standalone selection faces command is not verified in SolutionDomain Model'
####        
    def checkPostCommand(self):
        """Check whether output postprocessing format command line could be get"""
        mdl = SolutionDomainModel(self.case)
        mdl.setPostProFormat('CGNS')

        cmd_post = ' --cgns '
        assert mdl.getPostCommand() == cmd_post,\
            'Post command is not verified for postprocessing format in SolutionDomain Model'


    def checkGetStatusNodeAndGetFacesSelect(self):
        """Check whether status of node and selection pf fzces could be get - only for view"""
        mdl = SolutionDomainModel(self.case)
        select = {}
        select['color'] = '1 2 3'
        select['group'] = 'toto'
        select['fraction'] = '0.12'
        select['plan'] = '0.88'
        select['reverse'] = 'off'
        select['semiconf'] = 'on'
        mdl.setJoinMeshesStatus('on')
        mdl.addJoinFaces(select)
        node = mdl.node_join.xmlGetChildNode('faces_join')
        
        assert mdl.getStatusNode(node) == 'on',\
            'Could not get status of node in SolutionDomainModel'
            
        deux = {}
        deux['color'] = '4 5 6'
        deux['group'] = 'coucou'
        deux['fraction'] = '0.1'
        deux['plan'] = '0.8'
        deux['reverse'] = 'off'
        deux['semiconf'] = 'on'
        mdl.addJoinFaces(deux)
        list = {'group': 'coucou', 'reverse': 'off', 'color': '4 5 6', 'plan': '0.8', 'fraction': '0.1', 'semiconf': 'on'}
        
        assert mdl.getFacesSelect('faces_join', 1) == list,\
            'Could not get faces for view in SolutionDomainModel'



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
    print __file__
    runner = unittest.TextTestRunner()
    runner.run(suite1())
    runner.run(suite2())





