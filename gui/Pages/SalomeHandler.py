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
This module defines hooks with the SALOME plate-forme concerning
the graphical selection of the Groups.

This module contains the following classes and function:
- BoundaryGroup
- VolumeGroup
- runSolver
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import os
import os.path
import logging

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI
import CFDSTUDYGUI_ProcessMgr
from CFDSTUDYGUI_DataModel import _getStudy, _getEngine
from CFDSTUDYGUI_Commons import sg, sgPyQt
from salome import lcc
import smesh
import GEOM

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SalomeHandler")
log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------

aStudy = _getStudy()
aSMESH_SO = aStudy.FindComponent("SMESH")
aGEOM_SO = aStudy.FindComponent("GEOM")

#loading IORs
builder = aStudy.NewBuilder()
if aSMESH_SO != None:
    aSMESHEngine = lcc.FindOrLoadComponent("FactoryServer", "SMESH")
    builder.LoadWith(aSMESH_SO, aSMESHEngine)
if aGEOM_SO != None:
    aGEOMEngine = lcc.FindOrLoadComponent("FactoryServer", "GEOM")
    builder.LoadWith(aGEOM_SO, aGEOMEngine)

#-------------------------------------------------------------------------------

#def isSmeshAndGeomActivated():
#    if aSMESH_SO == None and aGEOM_SO == None:
#        return False
#    else:
#        return True


#def isSmeshActivated():
#    if aSMESH_SO == None:
#        tkMessageBox.showwarning ("WARNING", "MESH module is not activated.")
#        return False
#    else:
#        return True


def BoundaryGroup():
    """
    Import groups of faces.
    """
    if aSMESH_SO == None and aGEOM_SO == None:
        raise ValueError,  "Component SMESH and GEOM not found"

    local = ""
    if sg.SelectedCount() > 0:
        for i in range (sg.SelectedCount()):
            entry = sg.getSelected(i)
            if entry != '':
                sobj = aStudy.FindObjectID(entry)
                if sobj != None:
                    anObjectDS = sobj.GetObject()
                    if anObjectDS !=  None:

                        # check for smesh group
                        aSmeshObject = anObjectDS._narrow(smesh.SMESH_GroupBase)
                        #if aSmeshObject == None:
                        #    aSmeshObject = anObjectDS._narrow(smesh.SMESH_Group)
                        #if aSmeshObject == None:
                        #    aSmeshObject = anObjectDS._narrow(smesh.SMESH_GroupOnGeom)

                        if aSmeshObject != None and aSmeshObject.GetType() == smesh.FACE:
                            if not local:
                                local = aSmeshObject.GetName()
                            else:
                                local += ' or ' + aSmeshObject.GetName()

                        # check for geom group of faces
                        aGeomObject = anObjectDS._narrow(GEOM.GEOM_Object)
                        if aGeomObject != None and aGeomObject.GetType() == 37:
                            # check the group
                            # get all possible faces
                            import geompy
                            all_ids = geompy.SubShapeAllIDs(aGeomObject.GetMainShape(), geompy.ShapeType["FACE"])
                            cur_ids = geompy.GetObjectIDs(aGeomObject)
                            isValid = len(cur_ids) > 0 # not include empty list
                            if isValid:
                                for face_id in cur_ids:
                                    if not face_id in all_ids:
                                        #invalid id
                                        isValid = False
                                        break

                            if isValid:
                                if not local:
                                    local = aGeomObject.GetName()
                                else:
                                    local += ' or ' + aGeomObject.GetName()

    log.debug("BoundaryGroup -> %s" % str(local))
    return local


def VolumeGroup():
    """
    Import groups of solid.
    """
    if aSMESH_SO == None and aGEOM_SO == None:
        raise ValueError,  "Component SMESH and GEOM not found"

    local = ""
    if sg.SelectedCount() > 0:
        for i in range (sg.SelectedCount()):
            entry = sg.getSelected(i)
            if entry != '':
                sobj = aStudy.FindObjectID(entry)
                if sobj !=  None:
                    anObjectDS = sobj.GetObject()
                    #check for smesh group
                    if anObjectDS !=  None:
                        #aSmeshObject = anObjectDS._narrow(smesh.SMESH_Group)
                        aSmeshObject = anObjectDS._narrow(smesh.SMESH_GroupBase)
                        if aSmeshObject != None and aSmeshObject.GetType() == smesh.VOLUME:
                            if not local:
                                local = aSmeshObject.GetName()
                            else:
                                local += ' or ' + aSmeshObject.GetName()

                        # check for geom group of volumes
                        aGeomObject = anObjectDS._narrow(GEOM.GEOM_Object)
                        if aGeomObject != None and aGeomObject.GetType() == 37:
                            # check the group
                            # get all possible volumes
                            import geompy
                            all_ids = geompy.SubShapeAllIDs(aGeomObject.GetMainShape(), geompy.ShapeType["SOLID"])
                            cur_ids = geompy.GetObjectIDs(aGeomObject)
                            isValid = len(cur_ids) > 0 # not include empty list
                            if isValid:
                                for face_id in cur_ids:
                                    if not face_id in all_ids:
                                        # invalid id
                                        isValid = False
                                        break

                            if isValid:
                                if not local:
                                    local = aGeomObject.GetName()
                                else:
                                    local += ' or ' + aGeomObject.GetName()

    log.debug("VolumeGroup -> %s" % str(local))
    return local


#def VolumeGroupOld():
#    """
#    Import groups of cells.
#    """
#    if aSMESH_SO == None:
#        return False
#
#    #loading IORs
#    builder = aStudy.NewBuilder()
#    if aSMESH_SO != None:
#        aSMESHEngine = lcc.FindOrLoadComponent("FactoryServer", "SMESH")
#        builder.LoadWith(aSMESH_SO, aSMESHEngine)
#
#    local = ""
#    if sg.SelectedCount() > 0:
#        for i in range (sg.SelectedCount()):
#            entry = sg.getSelected(i)
#            if entry != '':
#                sobj = aStudy.FindObjectID(entry)
#                if sobj !=  None:
#                    anObjectDS = sobj.GetObject()
#                    #check for smesh group
#                    if anObjectDS !=  None:
#                        #aSmeshObject = anObjectDS._narrow(smesh.SMESH_Group)
#                        aSmeshObject = anObjectDS._narrow(smesh.SMESH_GroupBase)
#                        if aSmeshObject != None and aSmeshObject.GetType() == smesh.VOLUME:
#                            if not local:
#                                local = aSmeshObject.GetName()
#                            else:
#                                local = local + ' or ' + aSmeshObject.GetName()
#
#    return local


def runSolver(case, cmd, mdl, log_file):
    #_CommandMgr = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_CommandMgr()
    #_CommandMgr.runCommand( cmd, os.path.abspath(os.path.curdir), "", log_file)

    tmp_path = ''
    if mdl.dicoValues['CS_TMP_PREFIX']:
        tmp_path = mdl.dicoValues['CS_TMP_PREFIX']
    elif os.environ["NOM_ARCH"] == "OSF1":
        tmp_path = os.environ["TMPDIR"]
    elif os.environ["NOM_ARCH"] == "Linux_CCRT":
        tmp_path = os.environ["SCRATCHDIR"]
    elif os.environ["NOM_ARCH"] == "Linux_IA64":
        tmp_path = os.environ["SCRATCHDIR"]
    elif os.environ["NOM_ARCH"] == "Linux_SCALI":
        tmp_path = os.environ["HOME"]
    elif "TMPDIR" in os.environ and os.environ["TMPDIR"] != '' and os.environ["TMPDIR"] != '/tmp':
        tmp_path = os.environ["TMPDIR"]
    else:
        tmp_path = os.environ["HOME"]

    study_dir, case_name = os.path.split(case['case_path'])
    study_name = os.path.split(study_dir)[-1]
    tmp_path = os.path.join(tmp_path, "tmp_Saturne")
    tmp_path = os.path.join(tmp_path, study_name) + "." + case_name

    CFDSTUDYGUI.processMgr().addProcess(cmd, \
                                        case['case_path'], \
                                        case['xmlfile'], \
                                        tmp_path, \
                                        mdl.dicoValues['ARG_CS_OUTPUT'])

#-------------------------------------------------------------------------------
