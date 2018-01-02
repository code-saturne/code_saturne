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

"""
This module defines hooks with the SALOME plate-forme concerning
the graphical selection of the Groups.

This module contains the following classes and function:
- BoundaryGroup
- VolumeGroup
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

from code_saturne.Base.Toolbox import GuiParam

import CFDSTUDYGUI
from CFDSTUDYGUI_DataModel import _getStudy, _getEngine
from CFDSTUDYGUI_Commons import sg, sgPyQt

import salome
from salome import *

import SMESH
from salome.smesh import smeshBuilder

import GEOM
from salome.geom import geomBuilder

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SalomeHandler")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------

#loading IORs
aStudy = _getStudy()
builder = aStudy.NewBuilder()

sMeshComponent = salome.myStudy.FindComponent("SMESH")
if sMeshComponent != None:
    aSMESHEngine = lcc.FindOrLoadComponent("FactoryServer", "SMESH")
    builder.LoadWith(sMeshComponent, aSMESHEngine)

sGeomComponent = salome.myStudy.FindComponent("GEOM")
if sGeomComponent != None:
    aGEOMEngine = lcc.FindOrLoadComponent("FactoryServer", "GEOM")
    builder.LoadWith(sGeomComponent, aGEOMEngine)

#-------------------------------------------------------------------------------

def BoundaryGroup():
    """
    Import groups of faces.
    """
    if sMeshComponent == None and sGeomComponent == None:
        raise ValueError("Component SMESH and GEOM not found")

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
                        aSmeshObject = anObjectDS._narrow(SMESH.SMESH_GroupBase)
                        #if aSmeshObject == None:
                        #    aSmeshObject = anObjectDS._narrow(SMESH.SMESH_Group)
                        #if aSmeshObject == None:
                        #    aSmeshObject = anObjectDS._narrow(SMESH.SMESH_GroupOnGeom)
                        if aSmeshObject != None and aSmeshObject.GetType() == SMESH.FACE:
                            if not local:
                                local = aSmeshObject.GetName()
                            else:
                                local += ' or ' + aSmeshObject.GetName()

                        # check for geom group of faces
                        aGeomObject = anObjectDS._narrow(GEOM.GEOM_Object)
                        if aGeomObject != None and aGeomObject.GetType() == 37:
                            # check the group
                            # get all possible faces
                            all_ids = geomBuilder.SubShapeAllIDs(aGeomObject.GetMainShape(), geomBuilder.ShapeType["FACE"])
                            cur_ids = geomBuilder.GetObjectIDs(aGeomObject)
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
    if sMeshComponent == None and sGeomComponent == None:
        raise ValueError("Component SMESH and GEOM not found")

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
                        #aSmeshObject = anObjectDS._narrow(SMESH.SMESH_Group)
                        aSmeshObject = anObjectDS._narrow(SMESH.SMESH_GroupBase)
                        if aSmeshObject != None and aSmeshObject.GetType() == SMESH.VOLUME:
                            if not local:
                                local = aSmeshObject.GetName()
                            else:
                                local += ' or ' + aSmeshObject.GetName()

                        # check for geom group of volumes
                        aGeomObject = anObjectDS._narrow(GEOM.GEOM_Object)
                        if aGeomObject != None and aGeomObject.GetType() == 37:
                            # check the group
                            # get all possible volumes
                            all_ids = geomBuilder.SubShapeAllIDs(aGeomObject.GetMainShape(), geomBuilder.ShapeType["SOLID"])
                            cur_ids = geomBuilder.GetObjectIDs(aGeomObject)
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

#-------------------------------------------------------------------------------
